// The MUSCL limiter system is based on the kappa = 1/3, finite volume MUSCL scheme
// from the following papers
// https://doi.org/10.1016/j.jcp.2021.110640
// https://doi.org/10.22055/jacm.2020.32845.2088
//
// The flux is the same RoeM2 scheme used by RoeM1D, applied to the reconstructed
// left/right face states instead of to neighbouring cell averages.

use crate::pipes::{InteriorSolver, PipeState, apply_bc};
use nalgebra::{Matrix1xX, Matrix3, Matrix3x1, Matrix3xX};
use std::ops::AddAssign;

const KAPPA: f64 = 1.0 / 3.0; // MUSCL blend parameter
const C_M: f64 = (1.0 - KAPPA) / 4.0; // weight on the "backward" difference
const C_P: f64 = (1.0 + KAPPA) / 4.0; // weight on the "forward" difference

const RECONSTRUCT_PRIMITIVE: bool = false;
// false = reconstruct conserved (verified 3rd order);
// true = reconstruct primitive (more robust, expect ~2nd order on nonlinear problems)

///Bundles everything `decode_state` + `euler_flux` produce, so the RoeM2 face loop can
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
    dq: Matrix3xX<f64>,   // n_total - 1
    qmin: Matrix3xX<f64>, // n_total - 2
    qmax: Matrix3xX<f64>, // n_total - 2
    q_l: Matrix3xX<f64>,  // n_faces
    q_r: Matrix3xX<f64>,  // n_faces
    wl: Decoded,          // n_faces - primitives/flux decoded from q_l
    wr: Decoded,          // n_faces - primitives/flux decoded from q_r
    phi: Matrix3xX<f64>,  // n_faces - numerical flux at each face

    //primitives section
    prim: Matrix3xX<f64>, // n_total -- only touched when RECONSTRUCT_PRIMITIVE
    w_l: Matrix3xX<f64>,  // n_faces -- primitive-form face states, pre-conversion
    w_r: Matrix3xX<f64>,  // n_faces
}

impl Workspace {
    fn new(n_total: usize, n_faces: usize) -> Self {
        Workspace {
            dq: Matrix3xX::zeros(n_total - 1),
            qmin: Matrix3xX::zeros(n_total - 2),
            qmax: Matrix3xX::zeros(n_total - 2),
            q_l: Matrix3xX::zeros(n_faces),
            q_r: Matrix3xX::zeros(n_faces),
            wl: Decoded::zeros(n_faces),
            wr: Decoded::zeros(n_faces),
            phi: Matrix3xX::zeros(n_faces),
            prim: Matrix3xX::zeros(n_total),
            w_l: Matrix3xX::zeros(n_faces),
            w_r: Matrix3xX::zeros(n_faces),
        }
    }
}

///Componentwise min/max of q over the 3-cell window {i-1, i, i+1}
/// This is the Barth-Jespersen-style monotonicity bound: it defines how far a face
/// value can be pushed away from the cell average before it would create information
/// that didn't exist in any of the three cells feeding the reconstruction.
fn compute_windows(q: &Matrix3xX<f64>, qmin: &mut Matrix3xX<f64>, qmax: &mut Matrix3xX<f64>) {
    let n_slope = qmin.ncols();
    debug_assert_eq!(qmax.ncols(), n_slope);
    debug_assert_eq!(q.ncols(), n_slope + 2);

    for row in 0..3 {
        for i in 0..n_slope {
            let left = q[(row, i)];
            let center = q[(row, i + 1)];
            let right = q[(row, i + 2)];
            qmin[(row, i)] = left.min(center).min(right);
            qmax[(row, i)] = left.max(center).max(right);
        }
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
    q_l: &mut Matrix3xX<f64>,
    q_r: &mut Matrix3xX<f64>,
    qmin: &Matrix3xX<f64>,
    qmax: &Matrix3xX<f64>,
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
    debug_assert!(
        first - 1 + n_faces <= qmin.ncols(),
        "window stencil runs past the padded array"
    );

    // Left state at face k: extrapolated FORWARD from cell first-1+k, using that cell's
    // own backward difference dq[first-2+k] and forward difference dq[first-1+k]
    q_l.copy_from(&q.columns(first - 1, n_faces));
    q_l.zip_zip_apply(
        &dq.columns(first - 2, n_faces), // backward difference of the left cell
        &dq.columns(first - 1, n_faces), // forward difference of the left cell
        |ql, d_minus, d_plus| *ql += C_M * d_minus + C_P * d_plus,
    );

    // Right state at face k: extrapolated BACKWARD from cell first+k, using that cell's
    // own backward difference dq[first-1+k] and forward difference dq[first+k]
    q_r.copy_from(&q.columns(first, n_faces));
    q_r.zip_zip_apply(
        &dq.columns(first - 1, n_faces), // backward difference of the right cell
        &dq.columns(first, n_faces),     // forward difference of the right cell
        |qr, d_minus, d_plus| *qr -= C_P * d_minus + C_M * d_plus,
    );

    // clamp the generated left and right states based on the limits defined by
    // compute_windows in order to stop nonphysical information spread.
    // qmin/qmax are indexed by (cell index - 1).
    for row in 0..3 {
        for k in 0..n_faces {
            let wl = first - 2 + k; // window of the left cell  (first-1+k)
            let wr = first - 1 + k; // window of the right cell (first+k)
            q_l[(row, k)] = q_l[(row, k)].clamp(qmin[(row, wl)], qmax[(row, wl)]);
            q_r[(row, k)] = q_r[(row, k)].clamp(qmin[(row, wr)], qmax[(row, wr)]);
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

///Calculates the RoeM2 flux at every face given the LEFT and RIGHT reconstructed
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
        let a_l = (gamma * wl.p[k] / wl.rho[k]).sqrt(); // left speed of sound
        let u_r = wr.u[k]; // right velocity
        let a_r = (gamma * wr.p[k] / wr.rho[k]).sqrt(); // right speed of sound

        //intermediates
        let b1 = lambda[2].max((u_r + a_r).max(0.0));
        let b2 = lambda[0].min((u_l - a_l).min(0.0));
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

///Ties reconstruction + the Riemann solve + the flux divergence together for every
/// real cell. q must already have its ghosts filled by apply_bc.
///
/// Writes the RAW flux difference phi[k+1] - phi[k] into df, matching the convention
/// every other interior method uses -- `advance` (or the RK3 stages below) supply the
/// -(dt/dx) scaling.
fn residual(
    q: &Matrix3xX<f64>,
    ws: &mut Workspace,
    df: &mut Matrix3xX<f64>,
    first: usize,
    n_real: usize,
    gamma: f64,
) {
    if RECONSTRUCT_PRIMITIVE {
        conserved_block_to_primitive(q, &mut ws.prim, gamma);
        cell_differences(&ws.prim, &mut ws.dq);
        compute_windows(&ws.prim, &mut ws.qmin, &mut ws.qmax);
        reconstruct(
            &ws.prim,
            &ws.dq,
            &mut ws.w_l,
            &mut ws.w_r,
            &ws.qmin,
            &ws.qmax,
            first,
        );
        primitive_block_to_conserved(&ws.w_l, &mut ws.q_l, gamma);
        primitive_block_to_conserved(&ws.w_r, &mut ws.q_r, gamma);
    } else {
        cell_differences(q, &mut ws.dq); //calculates dq
        compute_windows(q, &mut ws.qmin, &mut ws.qmax); //compute bounds for left and right states
        reconstruct(
            q,
            &ws.dq,
            &mut ws.q_l,
            &mut ws.q_r,
            &ws.qmin,
            &ws.qmax,
            first,
        ); //reconstructs q_l and q_r using peicewise linear scheme
    }

    enforce_positivity(q, &mut ws.q_l, &mut ws.q_r, first, gamma); // last line of defense before the flux

    roe_flux(&ws.q_l, &ws.q_r, &mut ws.wl, &mut ws.wr, &mut ws.phi, gamma); //calculates the roe flux through each face

    // df_j = phi_{j+1/2} - phi_{j-1/2} for every real cell j
    ws.phi
        .columns(1, n_real)
        .sub_to(&ws.phi.columns(0, n_real), df);
}

///One SSP-RK3 stage, in place: dst[real] = a*q_n[real] + b*(src[real] + dt*R),
/// where R = -df/dx. Doing the blend in place keeps the time loop allocation free,
/// which the reference's expression-form stages did not.
#[allow(clippy::too_many_arguments)]
fn rk_stage(
    dst: &mut Matrix3xX<f64>,
    q_n: &Matrix3xX<f64>,
    src: &Matrix3xX<f64>,
    df: &Matrix3xX<f64>,
    a: f64,
    b: f64,
    dt_over_dx: f64,
    first: usize,
    n_real: usize,
) {
    let mut d = dst.columns_mut(first, n_real);

    // d = src + dt*R
    d.copy_from(&src.columns(first, n_real));
    d.zip_apply(df, |x, dfv| *x -= dt_over_dx * dfv);

    // d = a*q_n + b*d   (stage 1 is a=0, b=1, so skip the blend entirely)
    if a != 0.0 {
        d *= b;
        d.zip_apply(&q_n.columns(first, n_real), |x, qn| *x += a * qn);
    }
}

///Third order MUSCL reconstruction with SSP-RK3 time integration around the RoeM2 flux.
pub struct MusclRoeM1D {
    shared: PipeState,
    ws: Workspace,
    q_stage1: Matrix3xX<f64>, // n_total
    q_stage2: Matrix3xX<f64>, // n_total
}

impl MusclRoeM1D {
    pub(crate) fn new(state: PipeState) -> Self {
        let ws = Workspace::new(state.n_total, state.n_faces);
        let q_stage1 = state.q1.clone();
        let q_stage2 = state.q1.clone();
        Self {
            shared: state,
            ws,
            q_stage1,
            q_stage2,
        }
    }
}

impl InteriorSolver for MusclRoeM1D {
    fn state(&self) -> &PipeState {
        &self.shared
    }
    fn state_mut(&mut self) -> &mut PipeState {
        &mut self.shared
    }

    ///One residual evaluation from the current state, in the shared df convention.
    /// The RK3 `update` below drives its own stages instead of going through here,
    /// but keeping the convention means the default forward-euler update would still
    /// be correct.
    fn flux_divergence(&mut self) {
        let (first, n_real, gamma) = (self.shared.first, self.shared.n_real, self.shared.gamma);
        let Self { shared, ws, .. } = self;
        residual(&shared.q1, ws, &mut shared.df, first, n_real, gamma);
    }

    ///Advances the interior state by one full SSP-RK3 step of size dt.
    ///
    /// This method decodes its own face states inside `residual` and never reads the
    /// shared primitives, so the "re-decode between stages" rule on the trait does
    /// not apply here. `save_step` is skipped too - nothing reads q0, since `advance`
    /// is unused.
    fn update(&mut self, dt: f64) {
        let (first, n_real, n_ghost, gamma, dx) = (
            self.shared.first,
            self.shared.n_real,
            self.shared.n_ghost,
            self.shared.gamma,
            self.shared.dx,
        );
        let (left_bc, right_bc) = (self.shared.left_bc, self.shared.right_bc);
        let c = dt / dx;

        let Self {
            shared,
            ws,
            q_stage1,
            q_stage2,
        } = self;

        // Stage 1: s1 = q^n + dt*R(q^n)   -- an ordinary forward-Euler step
        residual(&shared.q1, ws, &mut shared.df, first, n_real, gamma);
        rk_stage(
            q_stage1, &shared.q1, &shared.q1, &shared.df, 0.0, 1.0, c, first, n_real,
        );
        apply_bc(q_stage1, first, n_real, n_ghost, left_bc, right_bc);

        // Stage 2: s2 = 3/4 q^n + 1/4 (s1 + dt*R(s1))
        residual(q_stage1, ws, &mut shared.df, first, n_real, gamma);
        rk_stage(
            q_stage2, &shared.q1, q_stage1, &shared.df, 0.75, 0.25, c, first, n_real,
        );
        apply_bc(q_stage2, first, n_real, n_ghost, left_bc, right_bc);

        // Stage 3: q^{n+1} = 1/3 q^n + 2/3 (s2 + dt*R(s2))
        //staged through q_stage1 first, since q1 is both source and destination here
        residual(q_stage2, ws, &mut shared.df, first, n_real, gamma);
        rk_stage(
            q_stage1,
            &shared.q1,
            q_stage2,
            &shared.df,
            1.0 / 3.0,
            2.0 / 3.0,
            c,
            first,
            n_real,
        );
        shared
            .q1
            .columns_mut(first, n_real)
            .copy_from(&q_stage1.columns(first, n_real));
        apply_bc(&mut shared.q1, first, n_real, n_ghost, left_bc, right_bc);
    }
}
