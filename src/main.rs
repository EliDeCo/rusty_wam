use nalgebra::{Matrix1xX, Matrix3xX, Vector3};
use std::{collections::BTreeMap, env};

mod boundaries;
mod driver;
mod helpers;
mod junctions;
mod pipe_methods;
mod pipes;
use boundaries::BoundaryCondition;
use driver::Driver;
use helpers::*;
use junctions::*;
use pipes::*;

//Input parameters
const COURANT: f64 = 0.9; //CFL courant number
const GAMMA: f64 = 1.4; //ratio of specific heats
const T_END: f64 = 3.0; //how much virtual time to run the simulation
const N_CELLS: usize = 2048; //how many real cells there are
const DOMAIN_LENGTH: f64 = 1.0; //basically how long the pipe is in meters
const N_PIPES: usize = 2; //number of pipes in the simulation
const N_JUNCTIONS: usize = 1; //number of junctions in the simulation
const PIPE_RADIUS: f64 = 30.0; //Pipe radius in mm
const METHOD: MethodKind = MethodKind::RoeM1D; //interior method every pipe uses
const RK_ORDER: TimeIntegrator = TimeIntegrator::Euler; //explicit SSP scheme, any method
const RUN_MODE: RunMode = RunMode::Transient; //how the simulation decides it is finished

//calculated parameters
const DX: f64 = DOMAIN_LENGTH / N_CELLS as f64; //step size

///When the simulation stops: at a fixed end time, or once it stops changing.
#[derive(Clone, Copy)]
#[allow(dead_code)] //selected by editing RUN_MODE
enum RunMode {
    ///advance to T_END
    Transient,
    ///advance until every residual has fallen by `tol`, or `max_it` iterations pass
    Steady { tol: f64, max_it: u32 },
}

///Largest residual magnitude held by each pipe and then each junction.
/// Kept per object so a pipe's flux difference is never compared against a junction's.
fn residuals(
    pipes: &BTreeMap<usize, InteriorMethod>,
    junctions: &BTreeMap<usize, Junction>,
) -> Vec<f64> {
    let from_pipes = pipes.values().map(|pipe| pipe.solver().state().df.amax());
    let from_junctions = junctions.values().map(|junction| junction.df.amax());

    from_pipes.chain(from_junctions).collect()
}

/// Left pressure = 1, right pressure = 0.1.
/// Left density = 1, right density = 0.125.
/// Velocity is 0 everywhere
fn _sods_problem() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let u0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let half = N_CELLS / 2;

    //left
    rho0.columns_mut(0, half).fill(1.0);
    p0.columns_mut(0, half).fill(1.0);

    //right
    rho0.columns_mut(half, N_CELLS - half).fill(0.125);
    p0.columns_mut(half, N_CELLS - half).fill(0.1);

    (rho0, u0, p0)
}

// The following are tests from section 5.1 of the reference paper
fn _shock_tube() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut u0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let half = N_CELLS / 2;

    //left
    rho0.columns_mut(0, half).fill(3.0);
    u0.columns_mut(0, half).fill(0.9);
    p0.columns_mut(0, half).fill(3.0);

    //right
    rho0.columns_mut(half, N_CELLS - half).fill(1.0);
    u0.columns_mut(half, N_CELLS - half).fill(0.9);
    p0.columns_mut(half, N_CELLS - half).fill(1.0);

    (rho0, u0, p0)
}

fn _contact_discontinuity() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut u0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let half = N_CELLS / 2;

    //left
    rho0.columns_mut(0, half).fill(10.0);
    u0.columns_mut(0, half).fill(0.1125);
    p0.columns_mut(0, half).fill(1.0);

    //right
    rho0.columns_mut(half, N_CELLS - half).fill(0.125);
    u0.columns_mut(half, N_CELLS - half).fill(0.1125);
    p0.columns_mut(half, N_CELLS - half).fill(1.0);

    (rho0, u0, p0)
}

fn _supersonic_expansion_test() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut u0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut p0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let half = N_CELLS / 2;

    //left
    rho0.columns_mut(0, half).fill(1.0);
    u0.columns_mut(0, half).fill(-2.0);
    p0.columns_mut(0, half).fill(3.0);

    //right
    rho0.columns_mut(half, N_CELLS - half).fill(1.0);
    u0.columns_mut(half, N_CELLS - half).fill(2.0);
    p0.columns_mut(half, N_CELLS - half).fill(3.0);

    (rho0, u0, p0)
}

///Custom test to confirm MUSCL + Limiter functionality
/// Pulse should stay thin and sharp and remain the same shape and size for the entire simulation.
/// Basic RoeM smears this horizontally, and the conservation of area under the curve causes
/// height to decrease as well. This is INCORRECT behavior
pub fn density_pulse_test() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Configuration: density pulse advection test.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 1.0); // rho_bg
    let u0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 0.5);
    let p0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 1.0);

    // top-hat: 2% of domain, starting near the inlet
    let pulse_start = (0.3 * N_CELLS as f64) as usize;
    let pulse_end = (0.32 * N_CELLS as f64) as usize;
    rho0.columns_mut(pulse_start, pulse_end - pulse_start)
        .fill(2.0);

    (rho0, u0, p0)
}

/// At rest
pub fn at_rest() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Pipe at rest.");

    let rho0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 1.0);
    let u0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let p0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 1.0);

    (rho0, u0, p0)
}

///Which pipe ends attach to which junctions, as a straight-through pair along x.
/// Normals point out of the junction along each pipe axis.
fn topology() -> Vec<Link> {
    vec![
        //pipe 0 runs up to the junction from -x, so it extends back along -x
        Link {
            junction: 0,
            pipe: 0,
            left: false,
            normal: Vector3::new(-1.0, 0.0, 0.0),
        },
        //pipe 1 carries the flow on, extending along +x
        Link {
            junction: 0,
            pipe: 1,
            left: true,
            normal: Vector3::new(1.0, 0.0, 0.0),
        },
    ]
}

fn main() {
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }

    let mut pipes: BTreeMap<usize, InteriorMethod> = BTreeMap::new();
    let mut junctions: BTreeMap<usize, Junction> = BTreeMap::new();

    //the network is described once and then drives both the pipe ends and the junctions
    let links = topology();

    //other variables
    let mut dt;
    let mut t = 0.0;
    let mut it = 0;

    for id in 0..N_PIPES {
        let (rho0, u0, p0) = match id {
            0 => density_pulse_test(),
            _ => at_rest(),
        };

        //initial total energy
        let e_tot0 = p0.component_div(&((GAMMA - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);

        //construct conservative state vector (past copy)
        let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
        q0.set_row(0, &rho0);
        q0.set_row(1, &rho0.component_mul(&u0));
        q0.set_row(2, &rho0.component_mul(&e_tot0));

        //an end named by the topology meets a junction, anything else passes waves out
        let end_bc = |left: bool| {
            links
                .iter()
                .find(|link| link.pipe == id && link.left == left)
                .map_or(BoundaryCondition::NonReflecting, |link| {
                    BoundaryCondition::Junction(link.junction)
                })
        };

        let pipe = InteriorMethod::new(
            METHOD,
            q0,
            GAMMA,
            COURANT,
            DX,
            id,
            PIPE_RADIUS,
            Some(end_bc(true)),
            Some(end_bc(false)),
        );

        pipes.insert(id, pipe);
    }

    //ideally the density pulse should travel through the junciton and continue into pipe 1
    for id in 0..N_JUNCTIONS {
        junctions.insert(id, Junction::new(GAMMA, COURANT, id));
    }

    for link in &links {
        let area = pipes[&link.pipe].solver().state().pipe_area();
        junctions
            .get_mut(&link.junction)
            .expect("link names a junction that does not exist")
            .add_pipe(link.pipe, link.left, link.normal, area);
    }

    for junction in junctions.values_mut() {
        junction.initialize(&pipes);
    }

    let chart = ChartDetails {
        width: 50,
        height: 50,
        x_min: 0.0,
        x_max: 1.0,
        y_min: 0.0,
        y_max: 2.0,
        x: (0..N_CELLS).map(|j| (j as f64 + 0.5) * DX).collect(),
    };

    let mut driver = Driver::new(RK_ORDER, &pipes, &junctions);

    println!("Beginning Simulation:");

    //a steady run has no end time to land on, so it never clips the last step
    let steady = matches!(RUN_MODE, RunMode::Steady { .. });
    let mut peak: Vec<f64> = Vec::new();

    loop {
        dt = pipes
            .values_mut()
            .fold(f64::INFINITY, |dt, pipe| dt.min(pipe.get_timestep()))
            .min(
                junctions
                    .values_mut()
                    .fold(f64::INFINITY, |dt, j| dt.min(j.get_timestep())),
            );
        if !steady {
            dt = dt.min(T_END - t);
        }

        for pipe in pipes.values() {
            //check Nan
            pipe.unreal_check();

            //temp display
            if it % 200 == 0 {
                plot(pipe.rho(), it, pipe.id(), &chart);
            }
        }

        for junction in junctions.values_mut() {
            junction.decode();
            junction.unreal_check();
        }

        //update cell states
        driver.step(dt, &mut pipes, &mut junctions);

        t += dt;
        it += 1;

        match RUN_MODE {
            RunMode::Transient => {
                if t >= T_END {
                    break;
                }
            }
            RunMode::Steady { tol, max_it } => {
                let now = residuals(&pipes, &junctions);
                if peak.is_empty() {
                    peak = now.clone();
                }

                //residuals climb before they fall, so convergence is measured against the
                //worst each object has reached rather than whatever the first step gave
                let worst = now
                    .iter()
                    .zip(peak.iter_mut())
                    .map(|(r, p)| {
                        *p = p.max(*r);
                        r / *p
                    })
                    .fold(0.0_f64, f64::max);

                if worst < tol || it >= max_it {
                    println!("Steady after {it} iterations, residual {worst:.3e}");
                    break;
                }
            }
        }
    }

    println!("Done");
}
