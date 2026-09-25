// The kappa = 1/3 finite volume MUSCL scheme, which is third order accurate, limited by
// Cada & Torrilhon's compact third-order limiter so that it stays third order at smooth
// extrema instead of clipping them. See validation/MusclRoeM1D.md
// https://doi.org/10.1016/j.jcp.2021.110640
// https://doi.org/10.1016/j.jcp.2009.02.020
//
// The flux is the same RoeM flux as RoeM1D, applied to the reconstructed left/right face
// states instead of to neighbouring cell averages.

use crate::pipes::{BoundaryPair, InteriorSolver, PipeState};
use nalgebra::{Matrix1xX, Matrix3, Matrix3x1, Matrix3xX};
use std::ops::AddAssign;

///MUSCL blend parameter. The limiter's bounds are derived for 1/3 and only this value
/// is third order accurate, so it is fixed rather than a knob.
const KAPPA: f64 = 1.0 / 3.0;

///Radius of the asymptotic region, dimensionless because the indicator below divides by
/// each row's own acoustic scale. Measured, not fitted: see validation/MusclRoeM1D.md
const R_SMOOTH: f64 = 1.0;

const RECONSTRUCT_PRIMITIVE: bool = false;
// false = reconstruct conserved (verified 3rd order);
// true = reconstruct primitive (more robust, expect ~2nd order on nonlinear problems)

///Bundles everything `decode_state` + `euler_flux` produce, so the RoeM face loop can
/// be called on either cell blocks or face blocks without a different function for each.
/// All fields have the same width.
struct Decoded {
    rho: Matrix1xX<f64>, //density
    u: Matrix1xX<f64>,   //velocity
    e: Matrix1xX<f64>,   // specific total energy
    p: Matrix1xX<f64>,   // pressure
    h: Matrix1xX<f64>,   // specific total enthalpy
    f: Matrix3xX<f64>,   // Euler flux F(q)
}

impl Decoded {
    fn zeros(width: usize) -> Self {
        Decoded {
            rho: Matrix1xX::zeros(width),
            u: Matrix1xX::zeros(width),
            e: Matrix1xX::zeros(width),
            p: Matrix1xX::zeros(width),
            h: Matrix1xX::zeros(width),
            f: Matrix3xX::zeros(width),
        }
    }
}

///Scratch buffers for one residual evaluation, sized once outside the time loop
/// (same preallocation style as the rest of the solver).
struct Workspace {
    dq: Matrix3xX<f64>,    // n_total - 1
    scale: Matrix3xX<f64>, // n_total - what the smoothness indicator divides by
    q_l: Matrix3xX<f64>,   // n_faces
    q_r: Matrix3xX<f64>,   // n_faces
    wl: Decoded,           // n_faces - primitives/flux decoded from q_l
    wr: Decoded,           // n_faces - primitives/flux decoded from q_r

    //primitives section
    prim: Matrix3xX<f64>, // n_total -- only touched when RECONSTRUCT_PRIMITIVE
    w_l: Matrix3xX<f64>,  // n_faces -- primitive-form face states, pre-conversion
    w_r: Matrix3xX<f64>,  // n_faces
}

impl Workspace {
    fn new(n_total: usize, n_faces: usize) -> Self {
        Workspace {
            dq: Matrix3xX::zeros(n_total - 1),
            scale: Matrix3xX::zeros(n_total),
            q_l: Matrix3xX::zeros(n_faces),
            q_r: Matrix3xX::zeros(n_faces),
            wl: Decoded::zeros(n_faces),
            wr: Decoded::zeros(n_faces),
            prim: Matrix3xX::zeros(n_total),
            w_l: Matrix3xX::zeros(n_faces),
            w_r: Matrix3xX::zeros(n_faces),
        }
    }
}

///The unlimited reconstruction written as a limiter, which is what makes it comparable
/// with the shock-capturing branch below. At kappa = 1/3 this is Cada's (2 + theta)/3.
fn phi_3(theta: f64) -> f64 {
    0.5 * ((1.0 - KAPPA) * theta + (1.0 + KAPPA))
}

///Cada's shock-capturing branch: the parabola, bounded so it cannot introduce variation
/// where the data is not monotone.
fn phi_hat(theta: f64) -> f64 {
    let p3 = phi_3(theta);
    let inner = (2.0 * theta).min(p3).min(1.6);
    p3.min((-0.5 * theta).max(inner)).max(0.0)
}

///Cada's limiter: the parabola inside the asymptotic region, the bounded branch outside
/// it, blended across a machine-width band so the numerical flux stays Lipschitz.
/// `ref_step` is R*dx scaled by the row's acoustic size, so eta is a pure number.
fn phi_limited(d_far: f64, d_near: f64, ref_step: f64) -> f64 {
    const EPS: f64 = 1.0e-12;
    let eta = (d_far * d_far + d_near * d_near) / (ref_step * ref_step);
    let theta = d_far / d_near;

    if eta <= 1.0 - EPS {
        phi_3(theta)
    } else if eta >= 1.0 + EPS {
        phi_hat(theta)
    } else {
        let w = (eta - 1.0) / EPS;
        0.5 * ((1.0 - w) * phi_3(theta) + (1.0 + w) * phi_hat(theta))
    }
}

///Half-step from a cell average to one of its faces. A zero near-difference means the
/// face sits at the cell value, and short-circuiting it keeps theta away from 0/0.
fn face_step(d_near: f64, d_far: f64, ref_step: f64) -> f64 {
    if d_near == 0.0 {
        return 0.0;
    }
    0.5 * phi_limited(d_far, d_near, ref_step) * d_near
}

///Reference step the smoothness indicator measures a difference against, per cell and
/// per row: R*dx times the row's acoustic size. Dividing by it is what makes a threshold
/// of one mean the same thing in a row of kg/m3 and a row of J/m3. Conserved rows scale
/// as (rho, rho a, rho a^2), primitive rows as (rho, a, rho a^2).
fn fill_scales(
    q: &Matrix3xX<f64>,
    scale: &mut Matrix3xX<f64>,
    dx: f64,
    gamma: f64,
    primitive: bool,
) {
    for col in 0..q.ncols() {
        let rho = q[(0, col)];
        let p = if primitive {
            q[(2, col)]
        } else {
            (gamma - 1.0) * (q[(2, col)] - 0.5 * q[(1, col)] * q[(1, col)] / rho)
        };
        let a = (gamma * p / rho).sqrt();
        let middle = if primitive { a } else { rho * a };

        scale[(0, col)] = R_SMOOTH * dx * rho;
        scale[(1, col)] = R_SMOOTH * dx * middle;
        scale[(2, col)] = R_SMOOTH * dx * rho * a * a;
    }
}

///decodes the state vector into primitives of the conserved variables:
/// rho (density), u (velocity), e (specific total energy), p (pressure)
/// q = [rho, rho*u, rho*E] where rho is the density, u is the velocity, and E is the total specific energy.
fn decode_state(q: &Matrix3xX<f64>, w: &mut Decoded, gamma: f64) {
    w.rho.copy_from(&q.row(0)); //density

    w.u.copy_from(&q.row(1)); //velocity
    w.u.component_div_assign(&w.rho); // u = (rho*u)/rho, in place

    w.e.copy_from(&q.row(2)); // specific total energy, NOT specific internal energy
    w.e.component_div_assign(&w.rho); // e = (rho*E)/rho, in place

    // pressure from equation of state
    //done step by step to avoid allocation
    w.p.copy_from(&w.u);
    w.p.component_mul_assign(&w.u); // p = u*u
    w.p *= -0.5; // p = -0.5*u*u
    w.p.add_assign(&w.e); // p = e - 0.5*u*u
    w.p.component_mul_assign(&w.rho); // p = rho*(e - 0.5*u*u)
    w.p *= gamma - 1.0; // p = (γ-1)*rho*(e - 0.5*u*u)

    // specific total enthalpy
    // computed in steps to avoid extra allocation
    w.h.copy_from(&w.p);
    w.h.component_div_assign(&w.rho); // h = p/rho
    w.h.add_assign(&w.e); // h = e + p/rho
}

///Calculates the Euler flux (F) for every column given the decoded primitives.
/// F = [rho*u, rho*u^2 + p, u*(rho*E + p)] where p is the pressure calculated from the equation of state.
fn euler_flux(w: &mut Decoded) {
    for i in 0..3 {
        for j in 0..w.f.ncols() {
            w.f[(i, j)] = match i {
                0 => w.rho[j] * w.u[j],                     // mass flux
                1 => w.rho[j] * w.u[j] * w.u[j] + w.p[j],   // momentum flux
                2 => w.u[j] * (w.rho[j] * w.e[j] + w.p[j]), // energy flux
                _ => panic!("Invalid index for flux calculation"),
            }
        }
    }
}

///Cell-to-cell differences over the padded (ghosts included) array: d_k = q_{k+1} - q_k.
/// One column shorter than q.
fn cell_differences(q: &Matrix3xX<f64>, dq: &mut Matrix3xX<f64>) {
    let n_diff = dq.ncols();
    q.columns(1, n_diff).sub_to(&q.columns(0, n_diff), dq);
}

///Reconstructs the left/right face states from cell averages using the kappa-family
/// MUSCL blend (kappa=1/3 by default), using a 4 wide stencil.
/// The `first`-relative offsets reduce to the reference's literals when first == 2.
fn reconstruct(
    q: &Matrix3xX<f64>,
    dq: &Matrix3xX<f64>,
    scale: &Matrix3xX<f64>,
    q_l: &mut Matrix3xX<f64>,
    q_r: &mut Matrix3xX<f64>,
    first: usize,
) {
    let n_faces = q_l.ncols();
    debug_assert!(
        first >= 2,
        "MUSCL reconstruction needs at least 2 ghost cells"
    );
    debug_assert!(
        first + n_faces <= dq.ncols(),
        "reconstruction stencil runs past the padded array"
    );

    for k in 0..n_faces {
        // the left state extrapolates forward out of cell first-1+k, the right state
        // backward out of cell first+k, each limited against its own cell's scale
        let (left, right) = (first - 1 + k, first + k);

        for row in 0..3 {
            let back = dq[(row, first - 2 + k)]; // backward difference of the left cell
            let mid = dq[(row, first - 1 + k)]; // shared by both cells
            let fwd = dq[(row, first + k)]; // forward difference of the right cell

            q_l[(row, k)] = q[(row, left)] + face_step(mid, back, scale[(row, left)]);
            q_r[(row, k)] = q[(row, right)] - face_step(mid, fwd, scale[(row, right)]);
        }
    }
}

///Converts an ENTIRE (3, W) block of conserved variables into primitives (rho, u, p).
fn conserved_block_to_primitive(q: &Matrix3xX<f64>, w: &mut Matrix3xX<f64>, gamma: f64) {
    for col in 0..q.ncols() {
        let rho = q[(0, col)];
        let u = q[(1, col)] / rho;
        let e = q[(2, col)] / rho;
        let p = (gamma - 1.0) * rho * (e - 0.5 * u * u);
        w[(0, col)] = rho;
        w[(1, col)] = u;
        w[(2, col)] = p;
    }
}

///Inverse of the above: primitive block -> conserved block. Needed because roe_flux
/// (and the positivity check below) only ever want to see conserved variables --
/// whichever variable set got reconstructed, q_l/q_r must come out in conserved form.
fn primitive_block_to_conserved(w: &Matrix3xX<f64>, q: &mut Matrix3xX<f64>, gamma: f64) {
    for col in 0..w.ncols() {
        let rho = w[(0, col)];
        let u = w[(1, col)];
        let p = w[(2, col)];
        let e = p / ((gamma - 1.0) * rho) + 0.5 * u * u;
        q[(0, col)] = rho;
        q[(1, col)] = rho * u;
        q[(2, col)] = rho * e;
    }
}

///Reads rho and p out of a single conserved-variable column.
fn rho_p_of(qcol: &Matrix3x1<f64>, gamma: f64) -> (f64, f64) {
    let rho = qcol[0];
    let u = qcol[1] / rho;
    let e = qcol[2] / rho;
    let p = (gamma - 1.0) * rho * (e - 0.5 * u * u);
    (rho, p)
}

///Per-face, per-side positivity fallback. If a reconstructed state has non-positive
/// density or pressure, that SIDE of that FACE drops back to its own cell average --
/// a local, first-order-only correction at exactly the offending point, not a panic
/// and not a change to any other face. q is the ORIGINAL padded cell-average array
/// (always conserved, regardless of which path built q_l/q_r), used as the fallback
/// source. Must run after reconstruction AND after any primitive->conserved
/// conversion, since it only knows how to read conserved columns.
fn enforce_positivity(
    q: &Matrix3xX<f64>,
    q_l: &mut Matrix3xX<f64>,
    q_r: &mut Matrix3xX<f64>,
    first: usize,
    gamma: f64,
) {
    for k in 0..q_l.ncols() {
        let (rho_l, p_l) = rho_p_of(&q_l.column(k).into_owned(), gamma);
        if rho_l <= 0.0 || p_l <= 0.0 {
            q_l.set_column(k, &(q.column(first - 1 + k).into_owned()));
        }

        let (rho_r, p_r) = rho_p_of(&q_r.column(k).into_owned(), gamma);
        if rho_r <= 0.0 || p_r <= 0.0 {
            q_r.set_column(k, &(q.column(first + k).into_owned()));
        }
    }
}

///Calculates the RoeM flux at every face given the LEFT and RIGHT reconstructed
/// states. This is face flux F, NOT the flux difference.
fn roe_flux(
    q_l: &Matrix3xX<f64>,
    q_r: &Matrix3xX<f64>,
    wl: &mut Decoded,
    wr: &mut Decoded,
    phi: &mut Matrix3xX<f64>,
    gamma: f64,
) {
    //state and decoded variables from each side
    decode_state(q_l, wl, gamma);
    decode_state(q_r, wr, gamma);

    //Nan Check
    for k in 0..phi.ncols() {
        if wl.rho[k] <= 0.0 || wl.p[k] <= 0.0 || wr.rho[k] <= 0.0 || wr.p[k] <= 0.0 {
            panic!(
                "positivity fallback did not resolve a bad state at face {} -- \
                cell average itself may be non-physical",
                k
            );
        }
    }

    //flux from each side
    euler_flux(wl);
    euler_flux(wr);

    //loop over each real cell face
    phi.column_iter_mut().enumerate().for_each(|(k, mut col)| {
        //intermediate quanties
        let r = (wr.rho[k] / wl.rho[k]).sqrt();

        //Roe averages (built directly from left/right states)
        let roe_rho = r * wl.rho[k]; // Roe average density
        let roe_u = (r * wr.u[k] + wl.u[k]) / (r + 1.0); // Roe average velocity
        let half_roe_u_squared = 0.5 * roe_u * roe_u; //intermediate quantity
        let roe_h = (r * wr.h[k] + wl.h[k]) / (r + 1.0); // Roe average specific total enthalpy
        let roe_a = ((gamma - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound

        //Eigenvalues
        let lambda: [f64; 3] = [roe_u - roe_a, roe_u, roe_u + roe_a];

        //now difference across faces, not cells
        let dq: Matrix3x1<f64> = q_r.column(k) - q_l.column(k);

        let u_l = wl.u[k]; // left veloctity
        let u_r = wr.u[k]; // right velocity

        //Eq 33: the signal velocities take the common speed of sound, which is what
        //lets a contact be captured exactly whichever side is the hotter
        let b1 = lambda[2].max((u_r + roe_a).max(0.0));
        let b2 = lambda[0].min((u_l - roe_a).min(0.0));
        let b3 = b1 + b2;
        let b4 = 2.0 * b1 * b2;
        let b5 = 1.0 / (b1 - b2);
        if !b5.is_normal() {
            panic!(
                "degenerate HLLE bounds at face {}: b1={}, b2={}, u_l={}, u_r={}",
                k, b1, b2, u_l, u_r
            );
        }

        // entropy-wave correction B∆Q,
        let m_hat = roe_u / roe_a; // Roe average Mach number
        let hlle_coeff = b1 * b2 * b5;
        let bdq_coeff = hlle_coeff / (1.0 + m_hat.abs()); // full prefactor

        let dp = wr.p[k] - wl.p[k]; // Δp

        let b_dq_0 = dq[0] - dp / (roe_a * roe_a); // Δρ - Δp/â²
        let b_dq: Matrix3x1<f64> = Matrix3x1::new(
            b_dq_0,
            b_dq_0 * roe_u,
            b_dq_0 * roe_h + roe_rho * (wr.h[k] - wl.h[k]),
        ); // ΔH

        let correction: Matrix3x1<f64> = bdq_coeff * b_dq;

        // accounts for swapping ΔQ -> ΔQ* = Δ(ρ,ρu,ρH) inside the HLLE base term
        let enthalpy_shift: Matrix3x1<f64> = hlle_coeff * Matrix3x1::new(0.0, 0.0, dp);

        //Eigenvector matrix P
        let p_matrix: Matrix3<f64> = Matrix3::new(
            1.0,
            1.0,
            1.0,
            lambda[0],
            lambda[1],
            lambda[2],
            roe_h - roe_u * roe_a,
            half_roe_u_squared,
            roe_h + roe_u * roe_a,
        );

        //more intermediate quantities
        let alpha2 = (gamma - 1.0) / (roe_a * roe_a);
        let alpha1 = alpha2 * half_roe_u_squared;
        let alpha3 = 0.5 / roe_a;

        let p_matrix_inv: Matrix3<f64> = Matrix3::new(
            alpha3 * roe_u + 0.5 * alpha1,
            -alpha3 - alpha2 * roe_u * 0.5,
            alpha2 * 0.5,
            1.0 - alpha1,
            alpha2 * roe_u,
            -alpha2,
            -alpha3 * roe_u + 0.5 * alpha1,
            alpha3 - alpha2 * roe_u * 0.5,
            alpha2 * 0.5,
        );

        //skip building eignvalue matrix since its sparse, construct manually instead
        // 1st, project dq into characteristic (wave) space: w = P^-1 * dq
        let mut w: Matrix3x1<f64> = p_matrix_inv * dq;

        //technically multiplying by diagonal wave speed matrix
        //Modified for RoeM
        w[0] *= (b3 * lambda[0] - b4) * b5;
        w[1] *= (b3 * lambda[1] - b4) * b5;
        w[2] *= (b3 * lambda[2] - b4) * b5;

        // Transform back to physical space: upwind dissipation = P * |Lambda| * P^-1 * dq
        let dissipation: Matrix3x1<f64> = p_matrix * w;

        // Roe flux at the interface: average of the two physical fluxes, minus half the upwind dissipation plus RoeM correction
        col.copy_from(
            &(0.5 * (wl.f.column(k) + wr.f.column(k)) - 0.5 * dissipation + enthalpy_shift
                - correction),
        );
    });
}

///Ties reconstruction and the Riemann solve together, filling phi at every face.
/// q must already have its ghosts filled by apply_bc.
fn fill_phi(
    q: &Matrix3xX<f64>,
    ws: &mut Workspace,
    phi: &mut Matrix3xX<f64>,
    first: usize,
    dx: f64,
    gamma: f64,
) {
    if RECONSTRUCT_PRIMITIVE {
        conserved_block_to_primitive(q, &mut ws.prim, gamma);
        cell_differences(&ws.prim, &mut ws.dq);
        fill_scales(&ws.prim, &mut ws.scale, dx, gamma, true);
        reconstruct(&ws.prim, &ws.dq, &ws.scale, &mut ws.w_l, &mut ws.w_r, first);
        primitive_block_to_conserved(&ws.w_l, &mut ws.q_l, gamma);
        primitive_block_to_conserved(&ws.w_r, &mut ws.q_r, gamma);
    } else {
        cell_differences(q, &mut ws.dq); //calculates dq
        fill_scales(q, &mut ws.scale, dx, gamma, false);
        //reconstructs q_l and q_r, limited so smooth extrema keep third order
        reconstruct(q, &ws.dq, &ws.scale, &mut ws.q_l, &mut ws.q_r, first);
    }

    enforce_positivity(q, &mut ws.q_l, &mut ws.q_r, first, gamma); // last line of defense before the flux

    roe_flux(&ws.q_l, &ws.q_r, &mut ws.wl, &mut ws.wr, phi, gamma); //calculates the roe flux through each face
}

///Third order MUSCL reconstruction around the RoeM flux.
pub struct MusclRoeM1D {
    shared: PipeState,
    ws: Workspace,
}

impl MusclRoeM1D {
    pub(crate) fn new(state: PipeState) -> Self {
        let ws = Workspace::new(state.n_total, state.n_faces);
        Self { shared: state, ws }
    }
}

impl InteriorSolver for MusclRoeM1D {
    fn state(&self) -> &PipeState {
        &self.shared
    }
    fn state_mut(&mut self) -> &mut PipeState {
        &mut self.shared
    }

    ///One residual evaluation for the state q, in the shared df convention.
    /// Decodes its own face states, so it never touches the shared primitives.
    fn residual(&mut self, q: &Matrix3xX<f64>, bc: &BoundaryPair) {
        let (first, dx, gamma) = (self.shared.first, self.shared.dx, self.shared.gamma);
        let Self { shared, ws, .. } = self;
        fill_phi(q, ws, &mut shared.phi, first, dx, gamma);
        shared.difference_flux(bc);
    }
}
