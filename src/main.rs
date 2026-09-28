use nalgebra::{Matrix3xX, Vector3};
use std::{collections::BTreeMap /*env*/};

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
//multiplier on the Courant number each method runs best at, so 1.0 is the efficient
//default whichever method is selected. See time_efficiency.md
const COURANT_SCALE: f64 = 1.0;
const GAMMA: f64 = 1.4; //ratio of specific heats
const T_END: f64 = 3.0; //how much virtual time to run the simulation
const N_CELLS: usize = 2048; //how many real cells there are
const DOMAIN_LENGTH: f64 = 1.0; //basically how long the pipe is in meters
const N_PIPES: usize = 2; //number of pipes in the simulation
const N_JUNCTIONS: usize = 1; //number of junctions in the simulation
const PIPE_RADIUS: f64 = 30.0; //Pipe radius in mm
const METHOD: MethodKind = MethodKind::RoeM1D; //interior method every pipe uses
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

///Writes the largest residual magnitude of each pipe and then each junction into `out`.
/// Kept per object so a pipe's flux difference is never compared against a junction's.
fn residuals(
    pipes: &BTreeMap<usize, InteriorMethod>,
    junctions: &BTreeMap<usize, Junction>,
    out: &mut [f64],
) {
    let from_pipes = pipes.values().map(|pipe| pipe.solver().state().df.amax());
    let from_junctions = junctions.values().map(|junction| junction.df.amax());

    for (slot, residual) in out.iter_mut().zip(from_pipes.chain(from_junctions)) {
        *slot = residual;
    }
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
    /*
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }
    */

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
            0 => density_pulse_test(N_CELLS),
            _ => at_rest(N_CELLS),
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
            METHOD.courant_from(COURANT_SCALE),
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
        junctions.insert(id, Junction::new(GAMMA, Junction::courant_from(COURANT_SCALE), id));
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

    let mut chart = ChartDetails {
        width: 50,
        height: 50,
        x_min: 0.0,
        x_max: 1.0,
        y_min: 0.0,
        y_max: 2.0,
        points: (0..N_CELLS)
            .map(|j| (((j as f64 + 0.5) * DX) as f32, 0.0))
            .collect(),
    };

    let mut driver = Driver::new(&pipes, &junctions);

    println!("Beginning Simulation:");

    //a steady run has no end time to land on, so it never clips the last step
    let steady = matches!(RUN_MODE, RunMode::Steady { .. });

    //one residual per pipe and junction, sized here so the loop itself never allocates
    let mut now: Vec<f64> = vec![0.0; pipes.len() + junctions.len()];
    let mut peak: Vec<f64> = now.clone();

    loop {
        //resolves the boundaries and fills the ghosts of the state the step starts from,
        //which the timestep below then decodes once for the whole step
        driver.prepare(&mut pipes, &mut junctions);

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
                plot(pipe.rho(), it, pipe.id(), &mut chart);
            }
        }

        for junction in junctions.values() {
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
                residuals(&pipes, &junctions, &mut now);

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
