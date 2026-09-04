use nalgebra::{Matrix1xX, Matrix3xX};
use std::{collections::BTreeMap, env};

mod helpers;
mod pipe_methods;
mod pipes;
use helpers::*;
use pipes::*;

//Input parameters
const COURANT: f64 = 0.9; //CFL courant number
const GAMMA: f64 = 1.4; //ratio of specific heats
const T_END: f64 = 1.0; //how much virtual time to run the simulation
const N_CELLS: usize = 2048;
const DOMAIN_LENGTH: f64 = 1.0; //basically how long the pipe is
const N_PIPES: usize = 2; //number of pipes in the simulation

//calculated parameters
const DX: f64 = DOMAIN_LENGTH / N_CELLS as f64; //step size

/// Left pressure = 1, right pressure = 0.1.
/// Left density = 1, right density = 0.125.
/// Velocity is 0 everywhere
fn sods_problem() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
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
pub fn _at_rest() -> (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>) {
    println!("Pipe at rest.");

    let rho0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 1.0);
    let u0: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let p0: Matrix1xX<f64> = Matrix1xX::from_element(N_CELLS, 1.0);

    (rho0, u0, p0)
}

fn main() {
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }

    let mut pipes: BTreeMap<usize, InteriorMethod> = BTreeMap::new();

    //other variables
    let mut dt;
    let mut t = 0.0;
    let mut it = 0;

    for id in 0..N_PIPES {
        let (rho0, u0, p0) = match id {
            0 => density_pulse_test(),
            _ => sods_problem(),
        };

        //initial total energy
        let e_tot0 = p0.component_div(&((GAMMA - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);

        //construct conservative state vector (past copy)
        let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
        q0.set_row(0, &rho0);
        q0.set_row(1, &rho0.component_mul(&u0));
        q0.set_row(2, &rho0.component_mul(&e_tot0));

        let pipe = InteriorMethod::new(
            MethodKind::RoeM1D,
            q0,
            GAMMA,
            COURANT,
            DX,
            N_CELLS,
            id,
            Some(BoundaryCondition::Transmissive),
            Some(BoundaryCondition::Transmissive),
        );

        pipes.insert(id, pipe);
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

    println!("Beginning Simulation:");

    while t < T_END {
        dt = pipes
            .values_mut()
            .fold(f64::INFINITY, |dt, pipe| dt.min(pipe.get_timestep()))
            .min(T_END - t);

        for pipe in pipes.values_mut() {
            //check Nan
            pipe.nan_check();

            //temp display
            if it % 200 == 0 {
                plot(pipe.rho(), it, pipe.id(), &chart);
            }

            //update cell states
            pipe.update(dt);
        }

        t += dt;
        it += 1;
    }
    println!("Done");
}
