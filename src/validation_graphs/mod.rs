// Runs a named benchmark problem instead of a network and graphs it against the analytic
// solution. The problems are the classical Riemann tests every compressible solver is
// expected to reproduce, so this is the visual counterpart to the grading in validation/.

pub mod exact;
mod window;

use crate::boundaries::BoundaryCondition;
use crate::driver::Driver;
use crate::helpers::{Initial, two_state};
use crate::junctions::Junction;
use crate::pipes::{InteriorMethod, MethodKind};
use exact::{Riemann, State};
use nalgebra::Matrix3xX;
use std::collections::BTreeMap;

///Which benchmark to run in place of the network, or None to run the network as usual.
#[derive(Clone, Copy, PartialEq)]
#[allow(dead_code)] //variants are selected by editing main.rs
pub enum Problem {
    None,
    Sod,
}

impl Problem {
    pub fn name(&self) -> &'static str {
        match self {
            Problem::None => "network",
            Problem::Sod => "Sod's shock tube",
        }
    }

    ///Length of the pipe the problem runs in, in meters.
    fn domain(&self) -> f64 {
        match self {
            Problem::None => 1.0,
            Problem::Sod => 1.0,
        }
    }

    ///Where the two initial states meet, which must land on a cell face.
    fn diaphragm(&self) -> f64 {
        self.domain() * 0.5
    }

    ///The two states either side of the diaphragm, which are the whole problem definition.
    /// The initial condition and the exact solution both come from here, so they cannot drift.
    fn states(&self) -> (State, State) {
        match self {
            //Sod 1978, section 3: the standard shock tube, both sides at rest
            Problem::None | Problem::Sod => {
                (State::new(1.0, 0.0, 1.0), State::new(0.125, 0.0, 0.1))
            }
        }
    }

    ///Density, velocity and pressure across the grid at t = 0.
    fn initial(&self, n_cells: usize) -> Initial {
        let (l, r) = self.states();
        two_state(n_cells, [l.rho, l.u, l.p], [r.rho, r.u, r.p])
    }

    ///The analytic solution of this problem, resolved once and then sampled anywhere.
    pub fn riemann(&self, gamma: f64) -> Riemann {
        let (l, r) = self.states();
        Riemann::new(l, r, self.diaphragm(), gamma)
    }
}

///Everything one benchmark run produced, which is both what gets graphed and what gets graded.
pub struct Solved {
    ///cell centres, in meters
    pub x: Vec<f64>,
    ///density, velocity and pressure over the real cells
    pub numerical: [Vec<f64>; 3],
    ///relative L1 error of each, in the same order
    pub errors: [f64; 3],
    pub steps: u32,
}

///Names of the three graphed fields, in the order everything here holds them.
pub const FIELDS: [&str; 3] = ["Density", "Velocity", "Pressure"];
pub const UNITS: [&str; 3] = ["kg/m3", "m/s", "Pa"];

///Runs the problem to exactly `t_end` and measures it against the analytic solution.
/// Opens no window, so a cell count sweep can call this in a loop.
pub fn measure(
    problem: Problem,
    t_end: f64,
    method: MethodKind,
    n_cells: usize,
    gamma: f64,
    courant_scale: f64,
) -> Solved {
    assert!(
        n_cells % 2 == 0,
        "an odd cell count puts the diaphragm inside a cell rather than on a face"
    );

    let dx = problem.domain() / n_cells as f64;
    let (numerical, steps) = solve(problem, t_end, method, n_cells, dx, gamma, courant_scale);
    let exact = problem.riemann(gamma).cell_averages(n_cells, dx, t_end);

    let errors = [
        relative_l1(&numerical[0], &exact[0]),
        relative_l1(&numerical[1], &exact[1]),
        relative_l1(&numerical[2], &exact[2]),
    ];

    Solved {
        x: (0..n_cells).map(|j| (j as f64 + 0.5) * dx).collect(),
        numerical,
        errors,
        steps,
    }
}

///Measures the problem, reports it, and then graphs it.
pub fn run(
    problem: Problem,
    t_end: f64,
    method: MethodKind,
    n_cells: usize,
    gamma: f64,
    courant_scale: f64,
) {
    reference(problem, t_end, gamma);

    let solved = measure(problem, t_end, method, n_cells, gamma, courant_scale);

    println!(
        "{} at t = {t_end}, {n_cells} cells, {} steps",
        problem.name(),
        solved.steps
    );
    for (k, field) in FIELDS.iter().enumerate() {
        println!(
            "  {field:<9} relative L1 error {:.4}%",
            100.0 * solved.errors[k]
        );
    }

    window::show(problem, t_end, n_cells, gamma, solved);
}

///Reports the analytic landmarks of the problem: where each wave has reached and what state
/// it leaves behind. This is what the exact solution itself is checked against, so that a
/// disagreement with the solver can never be blamed on an unverified reference.
pub fn reference(problem: Problem, t_end: f64, gamma: f64) {
    let riemann = problem.riemann(gamma);
    let x0 = problem.diaphragm();
    let speeds = riemann.wave_speeds();

    println!(
        "Exact solution: p* = {:.5}, u* = {:.5}",
        riemann.p_star(),
        riemann.u_star()
    );

    const EDGES: [&str; 5] = [
        "left wave head",
        "left wave tail",
        "contact",
        "right wave tail",
        "right wave head",
    ];
    print!("  waves at t = {t_end} reach x =");
    for (name, speed) in EDGES.iter().zip(speeds.iter()) {
        print!("  {name} {:.5}", x0 + speed * t_end);
    }
    println!();

    //sampling midway between neighbouring waves lands in the middle of each constant region
    let bounds: Vec<f64> = std::iter::once(0.0)
        .chain(speeds.iter().map(|s| x0 + s * t_end))
        .chain(std::iter::once(problem.domain()))
        .collect();

    for pair in bounds.windows(2) {
        //a shock leaves its two edges coincident, so that window holds no region at all
        if pair[1] - pair[0] < 1e-9 {
            continue;
        }
        let x = 0.5 * (pair[0] + pair[1]);
        let (rho, u, p) = riemann.at(x, t_end);
        println!("  x = {x:.4}:  rho {rho:.5}   u {u:.5}   p {p:.5}");
    }
}

///Advances one pipe holding the problem to exactly `t_end`, returning its final primitives.
fn solve(
    problem: Problem,
    t_end: f64,
    method: MethodKind,
    n_cells: usize,
    dx: f64,
    gamma: f64,
    courant_scale: f64,
) -> ([Vec<f64>; 3], u32) {
    let (rho0, u0, p0) = problem.initial(n_cells);
    let e_tot0 = p0.component_div(&((gamma - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);

    let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(n_cells);
    q0.set_row(0, &rho0);
    q0.set_row(1, &rho0.component_mul(&u0));
    q0.set_row(2, &rho0.component_mul(&e_tot0));

    //a single pipe has no junction to bind it, and both ends stay in their initial state for
    //as long as no wave reaches them, so the radius never enters the answer
    let mut pipes: BTreeMap<usize, InteriorMethod> = BTreeMap::new();
    pipes.insert(
        0,
        InteriorMethod::new(
            method,
            q0,
            gamma,
            method.courant_from(courant_scale),
            dx,
            0,
            30.0,
            Some(BoundaryCondition::NonReflecting),
            Some(BoundaryCondition::NonReflecting),
        ),
    );

    let mut junctions: BTreeMap<usize, Junction> = BTreeMap::new();
    let mut driver = Driver::new(&pipes, &junctions);

    let mut t = 0.0;
    let mut steps = 0;

    while t < t_end {
        driver.prepare(&mut pipes, &mut junctions);

        //clipping the last step is what lands the answer on the requested time exactly
        let dt = pipes
            .values_mut()
            .fold(f64::INFINITY, |dt, pipe| dt.min(pipe.get_timestep()))
            .min(t_end - t);

        for pipe in pipes.values() {
            pipe.unreal_check();
        }

        driver.step(dt, &mut pipes, &mut junctions);
        t += dt;
        steps += 1;
    }

    //the step left the primitives stale, since q1 moved underneath them
    let pipe = pipes.get_mut(&0).expect("the benchmark pipe");
    pipe.decode();

    (
        [pipe.rho().to_vec(), pipe.u().to_vec(), pipe.p().to_vec()],
        steps,
    )
}

///Error relative to the size of the exact solution, which is scale free and, unlike a max
/// norm, still converges when the solution contains a shock or a contact.
fn relative_l1(numerical: &[f64], exact: &[f64]) -> f64 {
    let (mut error, mut norm) = (0.0, 0.0);

    for (&n, &e) in numerical.iter().zip(exact.iter()) {
        error += (n - e).abs();
        norm += e.abs();
    }

    error / norm
}
