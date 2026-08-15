/*
Testing/verification-only architecture for the solver in main.rs: the Riemann/IC test
battery, the Phase 8 order-of-accuracy + shape-preservation apparatus, and the
ENABLE_CLIP toggle those tests need. Nothing in here is exercised by a normal
(non-testing) run of the solver -- see main.rs's run_default() for that.
*/
use crate::{
    BoundaryCondition, COURANT, DX, Decoded, FIRST, GAMMA, N_CELLS, N_TOTAL, T_END, Workspace,
    apply_bc, decode_state, max_wave_speed, ssprk3_step,
};
use nalgebra::{Matrix1xX, Matrix3xX};
use std::f64::consts::PI;

/// false reproduces unlimited kappa=1/3 exactly -- needed for the order test, which
/// must measure the UNLIMITED scheme to see slope 3.
pub const ENABLE_CLIP: bool = true;

///TESTS=============================================
/// Left pressure = 1, right pressure = 0.1.
/// Left density = 1, right density = 0.125.
/// Velocity is 0 everywhere
fn _sods_problem() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let u0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let half = N_TOTAL / 2;

    //left
    rho0.columns_mut(0, half).fill(1.0);
    p0.columns_mut(0, half).fill(1.0);

    //right
    rho0.columns_mut(half, N_TOTAL - half).fill(0.125);
    p0.columns_mut(half, N_TOTAL - half).fill(0.1);

    (rho0, u0, p0)
}

// The following are tests from section 5.1 of the reference paper
fn _shock_tube() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut u0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let half = N_TOTAL / 2;

    //left
    rho0.columns_mut(0, half).fill(3.0);
    u0.columns_mut(0, half).fill(0.9);
    p0.columns_mut(0, half).fill(3.0);

    //right
    rho0.columns_mut(half, N_TOTAL - half).fill(1.0);
    u0.columns_mut(half, N_TOTAL - half).fill(0.9);
    p0.columns_mut(half, N_TOTAL - half).fill(1.0);

    (rho0, u0, p0)
}

fn _contact_discontinuity() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut u0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let half = N_TOTAL / 2;

    //left
    rho0.columns_mut(0, half).fill(10.0);
    u0.columns_mut(0, half).fill(0.1125);
    p0.columns_mut(0, half).fill(1.0);

    //right
    rho0.columns_mut(half, N_TOTAL - half).fill(0.125);
    u0.columns_mut(half, N_TOTAL - half).fill(0.1125);
    p0.columns_mut(half, N_TOTAL - half).fill(1.0);

    (rho0, u0, p0)
}

fn _supersonic_expansion_test() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut u0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let half = N_TOTAL / 2;

    //left
    rho0.columns_mut(0, half).fill(1.0);
    u0.columns_mut(0, half).fill(-2.0);
    p0.columns_mut(0, half).fill(3.0);

    //right
    rho0.columns_mut(half, N_TOTAL - half).fill(1.0);
    u0.columns_mut(half, N_TOTAL - half).fill(2.0);
    p0.columns_mut(half, N_TOTAL - half).fill(3.0);

    (rho0, u0, p0)
}

///Custom test to confirm MUSCL + Limiter functionality
/// Pulse should stay thin and sharp and remain the same shape and size for the entire simulation.
/// Basic RoeM smears this horizontally, and the conservation of area under the curve causes
/// height to decrease as well. This is INCORRECT behavior
pub fn density_pulse_test() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration: density pulse advection test.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::from_element(N_TOTAL, 1.0); // rho_bg
    let u0: Matrix1xX<f64> = Matrix1xX::from_element(N_TOTAL, 0.5);
    let p0: Matrix1xX<f64> = Matrix1xX::from_element(N_TOTAL, 1.0);

    // top-hat: 2% of domain, starting near the inlet
    let pulse_start = (0.3 * N_TOTAL as f64) as usize;
    let pulse_end = (0.32 * N_TOTAL as f64) as usize;
    rho0.columns_mut(pulse_start, pulse_end - pulse_start)
        .fill(2.0);

    (rho0, u0, p0)
}

// Phase 8: verification apparatus. Two new test ICs (smooth periodic entropy wave,
// Gaussian pulse) plus the error/order machinery needed to turn "looks about right"
// into an actual measured number.

/// Exact CELL-AVERAGE of rho0(x) = 1 + 0.5*sin(2*pi*x) over [x_left, x_right], via the
/// closed-form antiderivative -- NOT rho0 evaluated at the cell center. Using the point
/// value here instead of the true average is exactly Pitfall 8/9 from van Leer &
/// Nishikawa: it silently caps the measured order at 2 no matter what the scheme does.
fn exact_density_cell_average(x_left: f64, x_right: f64) -> f64 {
    let h = x_right - x_left;
    let integral = (x_right - x_left)
        - 0.5 / (2.0 * PI) * ((2.0 * PI * x_right).cos() - (2.0 * PI * x_left).cos());
    integral / h
}

/// Periodic entropy-wave IC: rho varies, u and p are exactly constant. Because
/// lambda_2 = u is the only active wave speed (u is spatially uniform), rho0(x)
/// advects RIGIDLY at speed u with no shape change -- an exact solution exists at
/// every t, which is what makes this THE test for order of accuracy. Sod-style tests
/// have no smooth exact solution to measure a convergence rate against.
fn entropy_wave_ic() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    let mut rho = Matrix1xX::zeros(N_CELLS);
    for j in 0..N_CELLS {
        let x_left = j as f64 * DX;
        rho[j] = exact_density_cell_average(x_left, x_left + DX);
    }
    let u = Matrix1xX::from_element(N_CELLS, 1.0);
    let p = Matrix1xX::from_element(N_CELLS, 1.0);
    (rho, u, p)
}

/// Exact cell-averaged density at time t: same closed form, shifted by u*t = t
/// (wrapped for the periodic domain of length 1.0) before integrating -- rho0(x - t).
fn exact_density_at_time(t: f64) -> Matrix1xX<f64> {
    let shift = t.rem_euclid(1.0);
    let mut rho = Matrix1xX::zeros(N_CELLS);
    for j in 0..N_CELLS {
        let x_left = j as f64 * DX - shift;
        rho[j] = exact_density_cell_average(x_left, x_left + DX);
    }
    rho
}

/// L1 error between the solver's cell averages and the exact cell averages -- L1
/// rather than max-norm because it's far less sensitive to a single noisy cell, which
/// matters once the limiter (Phase 4/5) is on and clipping activity varies cell-to-cell
/// in a way L_inf would exaggerate.
fn l1_density_error(numerical: &Matrix1xX<f64>, exact: &Matrix1xX<f64>) -> f64 {
    (numerical - exact).iter().map(|e| e.abs()).sum::<f64>() / N_CELLS as f64
}

/// Observed order of accuracy from two runs at different resolutions. Turns "the error
/// went down when I doubled the cell count" into an actual exponent, so "order ~3" is
/// a printed number rather than something eyeballed off a log-log plot.
fn observed_order(error_coarse: f64, error_fine: f64, refinement_ratio: f64) -> f64 {
    (error_coarse / error_fine).log(refinement_ratio)
}

/// Smooth companion to density_pulse_test()'s top-hat. The top-hat is a genuine
/// discontinuity, so even a perfect scheme is only 1st-order there -- it can only
/// measure SMEARING, not accuracy. A Gaussian is smooth: it SHOULD hold its shape
/// under MUSCL and visibly smear under 1st-order, which is the visual complement to
/// the order-of-accuracy number above.
fn gaussian_pulse_test() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    let (center, width) = (0.25, 0.05);
    let mut rho = Matrix1xX::zeros(N_CELLS);
    for j in 0..N_CELLS {
        let x = (j as f64 + 0.5) * DX;
        rho[j] = 1.0 + 0.5 * (-((x - center) / width).powi(2)).exp();
    }
    let u = Matrix1xX::from_element(N_CELLS, 1.5);
    let p = Matrix1xX::from_element(N_CELLS, 1.0);
    (rho, u, p)
}

///Builds a padded (N_TOTAL-wide) conserved state from an IC given as REAL-cell-only
/// (N_CELLS-wide) primitives -- entropy_wave_ic and gaussian_pulse_test return this
/// shape (unlike the Riemann-test ICs below, which already come back N_TOTAL-wide).
/// Ghosts are left at zero here; the caller fills them with apply_bc.
fn pack_real_cells_into_padded(
    rho0: &Matrix1xX<f64>,
    u0: &Matrix1xX<f64>,
    p0: &Matrix1xX<f64>,
) -> Matrix3xX<f64> {
    let e_tot0 = p0.component_div(&((GAMMA - 1.0) * rho0)) + 0.5 * u0.component_mul(u0);
    let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_TOTAL);
    for j in 0..N_CELLS {
        q0[(0, FIRST + j)] = rho0[j];
        q0[(1, FIRST + j)] = rho0[j] * u0[j];
        q0[(2, FIRST + j)] = rho0[j] * e_tot0[j];
    }
    q0
}

///Recomputes dt from the current padded state, exactly mirroring the CFL update
/// inside the main time loop -- pulled out so the Phase 8 drivers below don't each
/// re-derive it.
fn cfl_dt(
    q: &Matrix3xX<f64>,
    cells: &mut Decoded,
    a: &mut Matrix1xX<f64>,
    t_end: f64,
    t: f64,
) -> f64 {
    decode_state(q, cells);
    a.copy_from(&cells.p);
    a.component_div_assign(&cells.rho);
    *a *= GAMMA;
    for x in a.iter_mut() {
        *x = x.sqrt();
    }
    (COURANT * DX / max_wave_speed(&cells.u, a)).min(t_end - t)
}

/// Phase 8, step 1: order-of-accuracy sweep driver. Runs the periodic entropy wave to
/// T_END=0.1 (short -- only the spatial order is being measured) at whatever N_CELLS /
/// ENABLE_CLIP the binary was built with, and prints a single grep-able RESULT line.
/// No chart, no per-iteration output -- this is meant to run unattended across the
/// resolution sweep in the shell loop below.
pub fn run_order_test() {
    let left_bc = BoundaryCondition::Periodic;
    let right_bc = BoundaryCondition::Periodic;
    let t_end_local = 0.1;

    let (rho0, u0, p0) = entropy_wave_ic();
    let mut q1 = pack_real_cells_into_padded(&rho0, &u0, &p0);
    apply_bc(&mut q1, left_bc, right_bc);

    let mut t: f64 = 0.0;
    let mut ws = Workspace::new();
    let mut dq_dt: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
    let mut cells = Decoded::zeros(N_TOTAL);
    let mut a: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut q_stage1: Matrix3xX<f64> = q1.clone();
    let mut q_stage2: Matrix3xX<f64> = q1.clone();

    let mut dt = cfl_dt(&q1, &mut cells, &mut a, t_end_local, t);

    while t < t_end_local {
        ssprk3_step(
            &mut q1,
            &mut q_stage1,
            &mut q_stage2,
            &mut dq_dt,
            &mut ws,
            dt,
            left_bc,
            right_bc,
        );
        t += dt;
        dt = cfl_dt(&q1, &mut cells, &mut a, t_end_local, t);
    }

    decode_state(&q1, &mut cells);
    let numerical_rho = cells.rho.columns(FIRST, N_CELLS).into_owned();
    let exact_rho = exact_density_at_time(t);
    let err = l1_density_error(&numerical_rho, &exact_rho);

    println!(
        "RESULT N_CELLS={} CLIP={} T_END={} L1={:.10e}",
        N_CELLS, ENABLE_CLIP, t_end_local, err
    );
}

/// Phase 8, step 2: Gaussian pulse shape-preservation driver. Transmissive both ends
/// (per the Phase 8 plan) so the pulse doesn't self-interact through a periodic wrap.
/// NOTE: at u=1.5 in a domain of length 1.0, the pulse (centered at x=0.25) has fully
/// exited by t=~0.5, long before the global T_END=3.0. Past that point the density
/// perturbation driving peak/mean_x/width collapses to ~0 and the width computation
/// divides by ~0 -- expected, not a bug; only t up to ~0.5 is a meaningful shape
/// measurement, and this is called out in the verification report.
pub fn run_gaussian_pulse() {
    let left_bc = BoundaryCondition::Transmissive;
    let right_bc = BoundaryCondition::Transmissive;

    let (rho0, u0, p0) = gaussian_pulse_test();
    let mut q1 = pack_real_cells_into_padded(&rho0, &u0, &p0);
    apply_bc(&mut q1, left_bc, right_bc);

    let mut t: f64 = 0.0;
    let mut it: u32 = 0;
    let x: Vec<f64> = (0..N_CELLS).map(|j| (j as f64 + 0.5) * DX).collect();

    let mut ws = Workspace::new();
    let mut dq_dt: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
    let mut cells = Decoded::zeros(N_TOTAL);
    let mut a: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut q_stage1: Matrix3xX<f64> = q1.clone();
    let mut q_stage2: Matrix3xX<f64> = q1.clone();

    let mut dt = cfl_dt(&q1, &mut cells, &mut a, T_END, t);

    println!("Beginning Gaussian pulse test:");

    while t < T_END {
        ssprk3_step(
            &mut q1,
            &mut q_stage1,
            &mut q_stage2,
            &mut dq_dt,
            &mut ws,
            dt,
            left_bc,
            right_bc,
        );
        t += dt;
        it += 1;
        dt = cfl_dt(&q1, &mut cells, &mut a, T_END, t);

        if cells.rho.columns(FIRST, N_CELLS).iter().any(|&v| v < 0.0)
            || cells.p.columns(FIRST, N_CELLS).iter().any(|&v| v < 0.0)
        {
            panic!("NANANNANANNANAN");
        }

        let numerical_rho = cells.rho.columns(FIRST, N_CELLS).into_owned();
        let peak = numerical_rho.iter().cloned().fold(f64::MIN, f64::max) - 1.0;
        let mean_x: f64 = x
            .iter()
            .zip(numerical_rho.iter())
            .map(|(&xi, &r)| xi * (r - 1.0))
            .sum::<f64>()
            / numerical_rho.iter().map(|&r| r - 1.0).sum::<f64>();
        let variance: f64 = x
            .iter()
            .zip(numerical_rho.iter())
            .map(|(&xi, &r)| (xi - mean_x).powi(2) * (r - 1.0))
            .sum::<f64>()
            / numerical_rho.iter().map(|&r| r - 1.0).sum::<f64>();
        let width = variance.sqrt();
        println!("it={} t={:.6} peak={:.6} width={:.6}", it, t, peak, width);
    }

    println!("RESULT gaussian_pulse final_it={} t={:.6}", it, t);
}

/// Shared driver for the Phase 8 regression suite: builds the IC, runs to the global
/// T_END under the given BCs, and reports the diagnostics from Phase 8 Sec 3.1.
fn run_regression_case(
    name: &str,
    ic: fn() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>),
    left_bc: BoundaryCondition,
    right_bc: BoundaryCondition,
) {
    let (rho0, u0, p0) = ic();

    let e_tot0 = p0.component_div(&((GAMMA - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);
    let a0: Matrix1xX<f64> = (GAMMA * p0.component_div(&rho0)).map(|x| x.sqrt());
    let mut dt = COURANT * DX / max_wave_speed(&u0, &a0);

    let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_TOTAL);
    q0.set_row(0, &rho0);
    q0.set_row(1, &(&rho0.component_mul(&u0)));
    q0.set_row(2, &(&rho0.component_mul(&e_tot0)));
    apply_bc(&mut q0, left_bc, right_bc);

    let mut q1: Matrix3xX<f64> = q0.to_owned();
    let mut t: f64 = 0.0;
    let mut it: u32 = 0;

    let mut ws = Workspace::new();
    let mut dq_dt: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
    let mut cells = Decoded::zeros(N_TOTAL);
    let mut a: Matrix1xX<f64> = Matrix1xX::zeros(N_TOTAL);
    let mut q_stage1: Matrix3xX<f64> = q1.clone();
    let mut q_stage2: Matrix3xX<f64> = q1.clone();

    let mut min_rho_seen = f64::MAX;
    let mut min_p_seen = f64::MAX;

    while t < T_END {
        ssprk3_step(
            &mut q1,
            &mut q_stage1,
            &mut q_stage2,
            &mut dq_dt,
            &mut ws,
            dt,
            left_bc,
            right_bc,
        );
        t += dt;
        it += 1;

        decode_state(&q1, &mut cells);
        let rho_real = cells.rho.columns(FIRST, N_CELLS);
        let p_real = cells.p.columns(FIRST, N_CELLS);

        min_rho_seen = rho_real.iter().fold(min_rho_seen, |m, &v| m.min(v));
        min_p_seen = p_real.iter().fold(min_p_seen, |m, &v| m.min(v));

        if rho_real.iter().any(|&v| v.is_nan()) || p_real.iter().any(|&v| v.is_nan()) {
            panic!("{}: NaN encountered at it={}", name, it);
        }

        dt = cfl_dt(&q1, &mut cells, &mut a, T_END, t);
    }

    // fallback_total reporting disabled -- Workspace.fallback_count is commented out
    // in main.rs (see enforce_positivity's call site in residual()); re-enable both
    // sides together if this metric is needed again.
    println!(
        "test={} final_it={} min_rho={:.6} min_p={:.6}",
        name, it, min_rho_seen, min_p_seen
    );
}

/// Phase 8, step 3: regression suite over the pre-MUSCL Riemann test battery.
/// BC choice: this repo has no committed pre-Phase-3 baseline to read a BC convention
/// off of (the ghost-cell/BoundaryCondition system itself was added as part of this
/// same MUSCL work), so Transmissive/Transmissive is used for all four -- the physically
/// standard choice for an open shock tube and consistent with what the Phase 8 plan
/// specifies for the Gaussian pulse test.
pub fn run_regression_suite() {
    run_regression_case(
        "_sods_problem",
        _sods_problem,
        BoundaryCondition::Transmissive,
        BoundaryCondition::Transmissive,
    );
    run_regression_case(
        "_shock_tube",
        _shock_tube,
        BoundaryCondition::Transmissive,
        BoundaryCondition::Transmissive,
    );
    run_regression_case(
        "_contact_discontinuity",
        _contact_discontinuity,
        BoundaryCondition::Transmissive,
        BoundaryCondition::Transmissive,
    );
    run_regression_case(
        "_supersonic_expansion_test",
        _supersonic_expansion_test,
        BoundaryCondition::Transmissive,
        BoundaryCondition::Transmissive,
    );
}
