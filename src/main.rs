/*
This code was originally based on P. S. Volpiani's Roe example (07_cfd.py)
https://github.com/psvolpiani/YouTube-CFD-101/tree/main/07_1D_Euler_equations_Roe

The main interior solver scheme is a 1D version of the RoeM2 scheme from the following paper,
witht the exceptions of the f and g functions which only give benefits in higher dimensions
https://doi.org/10.1016/S0021-9991(02)00037-2

The MUSCL limter system is based on the kappa = 1/3, finite volume MUSCL scheme from the following papers
https://doi.org/10.1016/j.jcp.2021.110640
https://doi.org/10.22055/jacm.2020.32845.2088
*/
use nalgebra::{Matrix1xX, Matrix3, Matrix3x1, Matrix3xX};
use std::ops::AddAssign;
use textplots::{Chart, Plot, Shape};
use std::collections::HashMap;

mod junctions;
mod testing;
use junctions::*;

//Input Parameters
const COURANT: f64 = 0.9; //CFL courant number
const GAMMA: f64 = 1.4; //ratio of specific heats
const T_END: f64 = 3.0; //how much virtual time to run the simulation
const KAPPA: f64 = 1.0 / 3.0; // MUSCL blend parameter
const N_CELLS: usize = 2048; //number of REAL (non ghost) cells
const N_GHOST: usize = 2; //ghost cells of each side
const DOMAIN_LENGTH: f64 = 1.0; //basically how long the pipe is

//temporary for junction testing
const N_PIPES: usize = 2;
const JUNCTION_ID: usize = N_PIPES + 1; //next free id

const RECONSTRUCT_PRIMITIVE: bool = false;
// false = reconstruct conserved (verified 3rd order);
// true = reconstruct primitive (more robust, expect ~2nd order on nonlinear problems)

//Calculated Parameters
const C_M: f64 = (1.0 - KAPPA) / 4.0; // weight on the "backward" difference
const C_P: f64 = (1.0 + KAPPA) / 4.0; // weight on the "forward" difference
const N_TOTAL: usize = N_CELLS + 2 * N_GHOST; //total cells stored, ghosts included
const N_DIFF: usize = N_TOTAL - 1; //number of cell-to-cell differences (one less than N_TOTAL)
const N_FACES: usize = N_CELLS + 1; //faces bounding the real cells, domain ends included
const FIRST: usize = N_GHOST; //index of the first REAL cell inside the padded array
const N_SLOPE: usize = N_TOTAL - 2; //One limited slope per cell that has both a backward AND forward difference (everything except the outermost ghost cells)
const DX: f64 = DOMAIN_LENGTH / N_CELLS as f64; //step size

///Componentwise min/max of q over the 3-cell window {i-1, i, i+1}
/// This is the Barth-Jespersen-style monotonicity bound: it defines how far a face
/// value can be pushed away from the cell average before it would create information
/// that didn't exist in any of the three cells feeding the reconstruction.
fn compute_windows(q: &Matrix3xX<f64>, qmin: &mut Matrix3xX<f64>, qmax: &mut Matrix3xX<f64>) {
    debug_assert_eq!(qmin.ncols(), N_SLOPE);
    debug_assert_eq!(qmax.ncols(), N_SLOPE);

    for row in 0..3 {
        for i in 0..N_SLOPE {
            let left = q[(row, i)];
            let center = q[(row, i + 1)];
            let right = q[(row, i + 2)];
            qmin[(row, i)] = left.min(center).min(right);
            qmax[(row, i)] = left.max(center).max(right);
        }
    }
}

#[derive(Clone, Copy)]
enum BoundaryCondition {
    /// zero-gradient: waves pass through the ends undisturbed
    Transmissive,
}

///Fills the ghost cells surrounding the real domain so every real cell — including the
/// two nearest each edge — can be updated with the exact same 4-cell stencil.
/// left/right are independent so e.g. a closed-end shock tube (wall + transmissive)
/// can be built from the same function.
fn apply_bc(q: &mut Matrix3xX<f64>, left: Option<BoundaryCondition>, right: Option<BoundaryCondition>) {
    match left {
        Some(BoundaryCondition::Transmissive) => {
            for g in 0..N_GHOST {
                let mirror = q.column(FIRST).into_owned();
                q.set_column(g, &mirror);
            }
        },
        _ => {}
    }

    match right {
        Some(BoundaryCondition::Transmissive) => {
            for g in 0..N_GHOST {
                let mirror = q.column(FIRST + N_CELLS - 1).into_owned();
                q.set_column(N_TOTAL - 1 - g, &mirror);
            }
        },
        _ => {}
    }
}

///Bundles everything `decode_state` + `euler_flux` produce, so the RoeM2 face loop can be
/// called on either cell blocks or face blocks without a different function for each.
/// All fields have the same width — N_CELLS for cell data, N_FACES for face data.
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

///decodes the state vector into primitives of the conserved variables:
/// rho (density), u (velocity), e (specific total energy), p (pressure)
/// They are vectors with one value per finite-volume cell
/// q = [rho, rho*u, rho*E] where rho is the density, u is the velocity, and E is the total specific energy.
fn decode_state(q: &Matrix3xX<f64>, w: &mut Decoded) {
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
    w.p *= GAMMA - 1.0; // p = (γ-1)*rho*(e - 0.5*u*u)

    // specific total enthalpy
    // computed in steps to avoid extra allocation
    w.h.copy_from(&w.p);
    w.h.component_div_assign(&w.rho); // h = p/rho
    w.h.add_assign(&w.e); // h = e + p/rho 
}

///Calculates the Euler flux (F) for every cell given state vector q and specific heat ratio gamma.
/// q = [rho, rho*u, rho*E] where rho is the density, u is the velocity, and E is the total specific energy.
/// F = [rho*u, rho*u^2 + p, u*(rho*E + p)] where p is the pressure calculated from the equation of state.
fn euler_flux(w: &mut Decoded) {
    for i in 0..3 {
        for j in 0..N_FACES {
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
/// One column shorter than q. Needed by the reconstruction in Phase 3; computed here as
/// its own step so reconstruct() never touches q's raw indices directly.
fn cell_differences(q: &Matrix3xX<f64>, dq: &mut Matrix3xX<f64>) {
    q.columns(1, N_DIFF).sub_to(&q.columns(0, N_DIFF), dq);
}

///Reconstructs the left/right face states from cell averages using the kappa-family
/// MUSCL blend (kappa=1/3 by default), using a 4 wide stencil
fn reconstruct(
    q: &Matrix3xX<f64>,
    dq: &Matrix3xX<f64>,
    q_l: &mut Matrix3xX<f64>,
    q_r: &mut Matrix3xX<f64>,
    qmin: &Matrix3xX<f64>,
    qmax: &Matrix3xX<f64>,
) {
    debug_assert!(
        2 + N_FACES <= N_DIFF,
        "reconstruction stencil runs past the padded array"
    );
    debug_assert!(
        1 + N_FACES <= N_SLOPE,
        "window stencil runs past the padded array"
    );
    debug_assert_eq!(dq.ncols(), N_DIFF);
    debug_assert!(
        2 + N_FACES <= N_DIFF,
        "reconstruction stencil runs past the padded array"
    );

    // Left state at face k: extrapolated FORWARD from cell k+1, using that cell's
    // own backward difference d[k] = q[k+1]-q[k] and forward difference d[k+1] = q[k+2]-q[k+1]
    q_l.copy_from(&q.columns(1, N_FACES));
    q_l.zip_zip_apply(
        &dq.columns(0, N_FACES), // backward difference of the left cell
        &dq.columns(1, N_FACES), // forward difference of the left cell
        |ql, d_minus, d_plus| *ql += C_M * d_minus + C_P * d_plus,
    );

    // Right state at face k: extrapolated BACKWARD from cell k+2, using that cell's
    // own backward difference d[k+1] and forward difference d[k+2]
    q_r.copy_from(&q.columns(2, N_FACES));
    q_r.zip_zip_apply(
        &dq.columns(1, N_FACES), // backward difference of the right cell
        &dq.columns(2, N_FACES), // forward difference of the right cell
        |qr, d_minus, d_plus| *qr -= C_P * d_minus + C_M * d_plus,
    );

    // clamp the generated left and right states based on the limits defined by compute_windows
    // in order to stop nonsphysical information spread
    for row in 0..3 {
        for k in 0..N_FACES {
            q_l[(row, k)] = q_l[(row, k)].clamp(qmin[(row, k)], qmax[(row, k)]);
            q_r[(row, k)] = q_r[(row, k)].clamp(qmin[(row, k + 1)], qmax[(row, k + 1)]);
        }
    }
}

//all the information for a single pipe
struct Pipe {
    workspace: Workspace,
    q1: Matrix3xX<f64>,
    q_stage1: Matrix3xX<f64>,
    q_stage2: Matrix3xX<f64>,
    dq_dt: Matrix3xX<f64>,
    left_bc: Option<BoundaryCondition>,
    right_bc: Option<BoundaryCondition>,
    cells: Decoded,
    x: Vec<f64>,
}

///Scratch buffers for one residual evaluation, sized once outside the time loop
/// (same preallocation style as the rest of the solver).
struct Workspace {
    dq: Matrix3xX<f64>,   // N_DIFF
    qmin: Matrix3xX<f64>, // N_SLOPE
    qmax: Matrix3xX<f64>, // N_SLOPE
    q_l: Matrix3xX<f64>,  // N_FACES
    q_r: Matrix3xX<f64>,  // N_FACES
    wl: Decoded,          // N_FACES — primitives/flux decoded from q_l
    wr: Decoded,          // N_FACES — primitives/flux decoded from q_r
    phi: Matrix3xX<f64>,  // N_FACES — numerical flux at each face

    //primitives section
    prim: Matrix3xX<f64>, // N_TOTAL -- only touched when RECONSTRUCT_PRIMITIVE
    w_l: Matrix3xX<f64>,  // N_FACES -- primitive-form face states, pre-conversion
    w_r: Matrix3xX<f64>,  // N_FACES
}

impl Workspace {
    fn new() -> Self {
        Workspace {
            dq: Matrix3xX::zeros(N_DIFF),
            qmin: Matrix3xX::zeros(N_SLOPE),
            qmax: Matrix3xX::zeros(N_SLOPE),
            q_l: Matrix3xX::zeros(N_FACES),
            q_r: Matrix3xX::zeros(N_FACES),
            wl: Decoded::zeros(N_FACES),
            wr: Decoded::zeros(N_FACES),
            phi: Matrix3xX::zeros(N_FACES),
            prim: Matrix3xX::zeros(N_TOTAL),
            w_l: Matrix3xX::zeros(N_FACES),
            w_r: Matrix3xX::zeros(N_FACES),
            // fallback_count: 0,
        }
    }
}

///Calculates the RoeM2 flux at every face given the LEFT and RIGHT reconstructed
/// This is face flux F, NOT the previously calculated ΔF
fn roe_flux(
    q_l: &Matrix3xX<f64>,
    q_r: &Matrix3xX<f64>,
    wl: &mut Decoded,
    wr: &mut Decoded,
    phi: &mut Matrix3xX<f64>,
) {
    //state and decoded variables from each side
    decode_state(q_l, wl);
    decode_state(q_r, wr);

    //Nan Check
    for k in 0..N_FACES {
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

        //Roe averages (now built directly from left/right states)
        let roe_rho = r * wl.rho[k]; // Roe average density
        let roe_u = (r * wr.u[k] + wl.u[k]) / (r + 1.0); // Roe average velocity
        let half_roe_u_squared = 0.5 * roe_u * roe_u; //intermediate quantity
        let roe_h = (r * wr.h[k] + wl.h[k]) / (r + 1.0); // Roe average specific total enthalpy
        let roe_a = ((GAMMA - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound

        //Eigenvalues
        let lambda: [f64; 3] = [roe_u - roe_a, roe_u, roe_u + roe_a];

        //now difference across faces, not cells
        let dq: Matrix3x1<f64> = q_r.column(k) - q_l.column(k);

        let u_l = wl.u[k]; // left veloctity
        let a_l = (GAMMA * wl.p[k] / wl.rho[k]).sqrt(); // left speed of sound
        let u_r = wr.u[k]; // right velocity
        let a_r = (GAMMA * wr.p[k] / wr.rho[k]).sqrt(); // right speed of sound

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
        let alpha2 = (GAMMA - 1.0) / (roe_a * roe_a);
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

///Ties reconstruction + the Riemann solve + the flux divergence together into a single
/// dq/dt for every real cell. q must already have its ghosts filled by apply_bc.
fn residual(q: &Matrix3xX<f64>, ws: &mut Workspace, dq_dt: &mut Matrix3xX<f64>) {
    if RECONSTRUCT_PRIMITIVE {
        conserved_block_to_primitive(q, &mut ws.prim);
        cell_differences(&ws.prim, &mut ws.dq);
        compute_windows(&ws.prim, &mut ws.qmin, &mut ws.qmax);
        reconstruct(
            &ws.prim,
            &ws.dq,
            &mut ws.w_l,
            &mut ws.w_r,
            &ws.qmin,
            &ws.qmax,
        );
        primitive_block_to_conserved(&ws.w_l, &mut ws.q_l);
        primitive_block_to_conserved(&ws.w_r, &mut ws.q_r);
    } else {
        cell_differences(q, &mut ws.dq); //calculates dq
        compute_windows(q, &mut ws.qmin, &mut ws.qmax); //compute bounds for left and right states
        reconstruct(q, &ws.dq, &mut ws.q_l, &mut ws.q_r, &ws.qmin, &ws.qmax); //reconstructs q_l and q_r using peicewise linear scheme
    }

    enforce_positivity(q, &mut ws.q_l, &mut ws.q_r); // last line of defense before the flux

    roe_flux(&ws.q_l, &ws.q_r, &mut ws.wl, &mut ws.wr, &mut ws.phi); //calculates the roe flux through each face

    // res_j = -(phi_{j+1/2} - phi_{j-1/2}) / dx for every real cell j
    ws.phi
        .columns(1, N_CELLS)
        .sub_to(&ws.phi.columns(0, N_CELLS), dq_dt);
    *dq_dt *= -1.0 / DX; //this is dq/dt, the state update
}

///Advances q by one full SSP-RK3 step of size dt.
fn ssprk3_step(
    q: &mut Matrix3xX<f64>,
    q1: &mut Matrix3xX<f64>,
    q2: &mut Matrix3xX<f64>,
    dq_dt: &mut Matrix3xX<f64>,
    ws: &mut Workspace,
    dt: f64,
    left_bc: Option<BoundaryCondition>,
    right_bc: Option<BoundaryCondition>,
) {
    // Stage 1: q1 = q^n + dt * R(q^n)   -- an ordinary forward-Euler step
    residual(q, ws, dq_dt);
    q1.columns_mut(FIRST, N_CELLS)
        .copy_from(&(q.columns(FIRST, N_CELLS) + dt * &*dq_dt));
    apply_bc(q1, left_bc, right_bc);

    // Stage 2: q2 = 3/4 q^n + 1/4 (q1 + dt * R(q1))
    residual(q1, ws, dq_dt);
    q2.columns_mut(FIRST, N_CELLS).copy_from(
        &(0.75 * q.columns(FIRST, N_CELLS) + 0.25 * (q1.columns(FIRST, N_CELLS) + dt * &*dq_dt)),
    );
    apply_bc(q1, left_bc, right_bc);

    // Stage 3: q^{n+1} = 1/3 q^n + 2/3 (q2 + dt * R(q2))
    //copy into q1 first to avoid borrowing q as mutable and immutable simultaneously
    residual(q2, ws, dq_dt);
    q1.columns_mut(FIRST, N_CELLS).copy_from(
        &((1.0 / 3.0) * q.columns(FIRST, N_CELLS)
            + (2.0 / 3.0) * (q2.columns(FIRST, N_CELLS) + dt * &*dq_dt)),
    );
    q.columns_mut(FIRST, N_CELLS)
        .copy_from(&q1.columns(FIRST, N_CELLS));

    apply_bc(q1, left_bc, right_bc);
}

///Converts an ENTIRE (3, W) block of conserved variables into primitives (rho, u, p) --
/// the vectorized counterpart to the single-column primitive_to_conserved already used
/// by the Enforced boundary condition, run in reverse and over every column at once.
fn conserved_block_to_primitive(q: &Matrix3xX<f64>, w: &mut Matrix3xX<f64>) {
    for col in 0..q.ncols() {
        let rho = q[(0, col)];
        let u = q[(1, col)] / rho;
        let e = q[(2, col)] / rho;
        let p = (GAMMA - 1.0) * rho * (e - 0.5 * u * u);
        w[(0, col)] = rho;
        w[(1, col)] = u;
        w[(2, col)] = p;
    }
}

///Inverse of the above: primitive block -> conserved block. Needed because roe_flux
/// (and the positivity check below) only ever want to see conserved variables --
/// whichever variable set got reconstructed, q_l/q_r must come out in conserved form.
fn primitive_block_to_conserved(w: &Matrix3xX<f64>, q: &mut Matrix3xX<f64>) {
    for col in 0..w.ncols() {
        let rho = w[(0, col)];
        let u = w[(1, col)];
        let p = w[(2, col)];
        let e = p / ((GAMMA - 1.0) * rho) + 0.5 * u * u;
        q[(0, col)] = rho;
        q[(1, col)] = rho * u;
        q[(2, col)] = rho * e;
    }
}

///Reads rho and p out of a single conserved-variable column -- the positivity check's
/// counterpart to pressure_check-style helpers used elsewhere in the solver.
fn rho_p_of(qcol: &Matrix3x1<f64>) -> (f64, f64) {
    let rho = qcol[0];
    let u = qcol[1] / rho;
    let e = qcol[2] / rho;
    let p = (GAMMA - 1.0) * rho * (e - 0.5 * u * u);
    (rho, p)
}

///Per-face, per-side positivity fallback. If a reconstructed state has non-positive
/// density or pressure, that SIDE of that FACE drops back to its own cell average --
/// a local, first-order-only correction at exactly the offending point, not a panic
/// and not a change to any other face. q is the ORIGINAL padded cell-average array
/// (always conserved, regardless of which path built q_l/q_r), used as the fallback
/// source. Must run after reconstruction AND after any primitive->conserved
/// conversion, since it only knows how to read conserved columns.
fn enforce_positivity(q: &Matrix3xX<f64>, q_l: &mut Matrix3xX<f64>, q_r: &mut Matrix3xX<f64>) {
    //let mut count = 0;
    for k in 0..N_FACES {
        let (rho_l, p_l) = rho_p_of(&q_l.column(k).into_owned());
        if rho_l <= 0.0 || p_l <= 0.0 {
            q_l.set_column(k, &(q.column(k + 1).into_owned()));
            //count += 1;
        }

        let (rho_r, p_r) = rho_p_of(&q_r.column(k).into_owned());
        if rho_r <= 0.0 || p_r <= 0.0 {
            q_r.set_column(k, &(q.column(k + 2).into_owned()));
            //count += 1;
        }
    }
    //count
}

//compute max wave speed for real cells only
fn max_wave_speed(u: &Matrix1xX<f64>, a: &Matrix1xX<f64>) -> f64 {
    u.columns(FIRST, N_CELLS)
        .iter()
        .zip(a.columns(FIRST, N_CELLS).iter())
        .fold(0.0_f64, |speed, (&ui, &ai)| speed.max(ui.abs() + ai))
}

fn main() {
    run_pipes();

}

fn run_pipes() {
    /*
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }
    */

    let mut pipes: HashMap<usize, Pipe> = HashMap::new();
    let mut dt;
    let mut next_dt = f64::INFINITY;
    let mut t: f64 = 0.0;
    let mut it: u32 = 0;
    let mut a: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);

    //initialize all pipes
    for id in 0..N_PIPES {
        //for now, do this manually
        let left_bc = Some(BoundaryCondition::Transmissive);
        //junction
        //let right_bc = None;
        let right_bc = Some(BoundaryCondition::Transmissive);


        let (rho0, u0, p0) = match id {
            0 => testing::density_pulse_test(),
            _ => testing::at_rest(),
        };

        //initial total energy
        let e_tot0 = p0.component_div(&((GAMMA - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);

        //initial speed of sound
        let a: Matrix1xX<f64> = (GAMMA * p0.component_div(&rho0)).map(|x| x.sqrt());

        //prepare next timestep
        next_dt = next_dt.min(COURANT * DX / max_wave_speed(&u0, &a));

        //construct conservative state vector, now padded with ghost cells (past copy)
        let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_TOTAL);
        q0.set_row(0, &rho0);
        q0.set_row(1, &(&rho0.component_mul(&u0)));
        q0.set_row(2, &(&rho0.component_mul(&e_tot0)));
        apply_bc(&mut q0, left_bc, right_bc); // initialize ghost cells

        //working copy
        let q1: Matrix3xX<f64> = q0.to_owned();

        // Cell-center coordinates for the REAL cells only: x_j = (j + 1/2)*dx.
        let x: Vec<f64> = (0..N_CELLS).map(|j| (j as f64 + 0.5) * DX).collect();

        let workspace = Workspace::new(); // face-width scratch for residual()
        let dq_dt: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS); // dq/dt, real cells only
        let cells = Decoded::zeros(N_TOTAL); // full-domain primitives, for CFL + display
        let q_stage1: Matrix3xX<f64> = q1.clone();
        let q_stage2: Matrix3xX<f64> = q1.clone();

        pipes.insert(id,Pipe {
            workspace,
            q1,
            q_stage1,
            q_stage2,
            dq_dt,
            left_bc,
            right_bc,
            cells,
            x
        });
    }

    /*
    let left_bc = BoundaryCondition::Periodic;
    let right_bc = BoundaryCondition::Periodic;

    let (rho0, u0, p0) = testing::density_pulse_test();

    //initial total energy
    let e_tot0 = p0.component_div(&((GAMMA - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);

    //initial speed of sound
    let mut a: Matrix1xX<f64> = (GAMMA * p0.component_div(&rho0)).map(|x| x.sqrt());

    //time step
    let mut dt = COURANT * DX / max_wave_speed(&u0, &a);

    //construct conservative state vector, now padded with ghost cells (past copy)
    let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_TOTAL);
    q0.set_row(0, &rho0);
    q0.set_row(1, &(&rho0.component_mul(&u0)));
    q0.set_row(2, &(&rho0.component_mul(&e_tot0)));
    apply_bc(&mut q0, left_bc, right_bc); // initialize ghost cells

    //working copy
    let mut q1: Matrix3xX<f64> = q0.to_owned();

    let mut t: f64 = 0.0;
    let mut it: u32 = 0;

    // Cell-center coordinates for the REAL cells only: x_j = (j + 1/2)*dx.
    let x: Vec<f64> = (0..N_CELLS).map(|j| (j as f64 + 0.5) * DX).collect();

    //assigned only once to avoid repeated allocation
    let mut ws = Workspace::new(); // face-width scratch for residual()
    let mut dq_dt: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS); // dq/dt, real cells only
    let mut cells = Decoded::zeros(N_TOTAL); // full-domain primitives, for CFL + display
    let mut q_stage1: Matrix3xX<f64> = q1.clone();
    let mut q_stage2: Matrix3xX<f64> = q1.clone();
    */

    println!("Beginning Simulation:");

    while t < T_END {
        dt = next_dt;
        next_dt = f64::INFINITY;
        t += dt;
        it += 1;
        println!("Iteration {}", it);
        for (id,pipe) in pipes.iter_mut() {
            //full q update occurs within this function
            ssprk3_step(
                &mut pipe.q1,
                &mut pipe.q_stage1,
                &mut pipe.q_stage2,
                &mut pipe.dq_dt,
                &mut pipe.workspace,
                dt,
                pipe.left_bc,
                pipe.right_bc
            );

            //decode state variables for operations below
            decode_state(&pipe.q1, &mut pipe.cells);

            //reuse shared "a" buffer to calculate local speed of sound for this pipe
            a.copy_from(&pipe.cells.p);
            a.component_div_assign(&pipe.cells.rho); // a = p/rho
            a *= GAMMA; // a = gamma*p/rho
            for x in a.iter_mut() {
                *x = x.sqrt(); // a = sqrt(gamma*p/rho)
            }

            //update timestep
            next_dt = next_dt
                .min(COURANT * DX / max_wave_speed(&pipe.cells.u, &a))
                .min(T_END - t);

            //NaN check for REAL cells only
            if pipe
                .cells
                .rho
                .columns(FIRST, N_CELLS)
                .iter()
                .any(|&x| x < 0.0)
                || pipe
                    .cells
                    .p
                    .columns(FIRST, N_CELLS)
                    .iter()
                    .any(|&x| x < 0.0)
            {
                panic!("NANANNANANNANAN in pipe {}", id);
            }

            //display pressure along the pipe
            println!("Pipe {}", id);
            let points: Vec<(f32, f32)> = pipe
                .x
                .iter()
                .copied()
                .map(|y| y as f32)
                .zip(pipe.cells.rho.iter().copied().map(|y| y as f32))
                .collect();
            Chart::new_with_y_range(50, 50, 0.0, 1.0, 1.0, 2.0)
                .lineplot(&Shape::Points(points.as_slice()))
                .display();
        }
    }
    println!("Done");
}
