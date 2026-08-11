use nalgebra::{Matrix3,Matrix3x1};
use ndarray::{Array1, Array2, Axis, linspace, parallel::prelude::*, s};
use std::env;
use textplots::{Chart, Plot, Shape};

//Parameters
const COURANT: f64 = 1.0; //CFL courant number
const GAMMA: f64 = 1.4; //ratio if specific heats
const T_END: f64 = 1.0; //how much virtual time to run the simulation
const N_CELLS: f64 = 400.0;
const DOMAIN_LENGTH: f64 = 1.0; //basically how long the pipe is
const DX: f64 = DOMAIN_LENGTH / N_CELLS; //step size
const N_BOUNDARIES: usize = N_CELLS as usize + 1;

///decodes the state vector into primitives of the conserved variables:
/// rho (density), u (velocity), e (specific total energy), p (pressure)
/// They are vectors of length nx, where nx is the number of grid points in the simulation domain
fn decode_state(
    q: &Array2<f64>,
    gamma: f64,
) -> (
    Array1<f64>,
    Array1<f64>,
    Array1<f64>,
    Array1<f64>,
    Array1<f64>,
) {
    let rho: Array1<f64> = q.row(0).to_owned();
    let rho_inv: Array1<f64> = 1.0 / &rho;
    let u: Array1<f64> = &q.row(1) * &rho_inv;
    let e: Array1<f64> = &q.row(2) * &rho_inv; // specific total energy, NOT specific internal energy
    let p: Array1<f64> = (gamma - 1.0) * &rho * (&e - 0.5 * &u * &u); // pressure from equation of state

    (rho, rho_inv, u, e, p)
}

///Calculates the Euler flux (F) for ever cell given state vector q and specific heat ratio gamma.
/// q = [rho, rho*u, rho*E] where rho is the density, u is the velocity, and E is the total specific energy.
/// F = [rho*u, rho*u^2 + p, u*(rho*E + p)] where p is the pressure calculated from the equation of state.
fn euler_flux(q: &Array2<f64>, gamma: f64) -> Array2<f64> {
    let (rho, _, u, e, p) = decode_state(q, gamma);

    let flux: Array2<f64> = Array2::from_shape_fn((3, q.ncols()), |(i, j)| {
        match i {
            0 => rho[j] * u[j],                 // mass flux
            1 => rho[j] * u[j] * u[j] + p[j],   // momentum flux
            2 => u[j] * (rho[j] * e[j] + p[j]), // energy flux
            _ => panic!("Invalid index for flux calculation"),
        }
    });

    flux
}
//TODO: Remove extra allocations
///Calculates the Roe flux for every cell given the state vector q and specific heat ratio gamma.
fn roe_flux(q: &Array2<f64>, gamma: f64, n: usize) -> Array2<f64> {
    let (rho, rho_inv, u, e, p) = decode_state(q, gamma);
    let h: Array1<f64> = e + p * &rho_inv; // specific total enthalpy

    //initialize roe flux
    let mut phi: Array2<f64> = Array2::zeros((3, n - 1));

    //loop over each cell boundary (column in phi)
    phi.axis_iter_mut(Axis(1))
        .into_par_iter()
        .enumerate()
        .for_each(|(i, mut col)| {
            //intermediate quanties
            let r = (rho[i + 1] * rho_inv[i]).sqrt();

            //Roe averages
            //let roe_rho = r * rho[i]; // Roe average density, unused in flux calculation
            let roe_u = (r * u[i + 1] + u[i]) / (r + 1.0); // Roe average velocity
            let half_roe_u_squared = 0.5 * roe_u * roe_u; //intermediate quantity
            let roe_h = (r * h[i + 1] + h[i]) / (r + 1.0); // Roe average specific total enthalpy
            let roe_a = ((gamma - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound

            //difference between neighboring cell states
            let dq: Matrix3x1<f64> = Matrix3x1::new(
                q[[0, i + 1]] - q[[0, i]], // d(rho)
                q[[1, i + 1]] - q[[1, i]], // d(rho*u)
                q[[2, i + 1]] - q[[2, i]], // d(rho*E)
            );

            //Eigenvalues
            let lambda: [f64; 3] = [roe_u - roe_a, roe_u, roe_u + roe_a];

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

            //Eigenvalue matrix
            let lambda_matrix: Matrix3<f64> = Matrix3::new(
                    lambda[0].abs(),
                    0.0,
                    0.0,
                    0.0,
                    lambda[1].abs(),
                    0.0,
                    0.0,
                    0.0,
                    lambda[2].abs(),
            );

            //Roe flux jacobian (the 3x3 matrix, constructed from its diagonalization)
            let a: Matrix3<f64> = p_matrix * lambda_matrix * p_matrix_inv;

            let product: Array1<f64> = Array1::from_vec((a * dq).as_slice().to_vec());

            //update roe flux
            col.assign(&product);
        });

    //get physical (euler) flux
    let f: Array2<f64> = euler_flux(q, gamma);

    //Roe numerical flux (phi is currently just 2 * the upwind dissipation, need to divide by 2 and add central flux)
    phi = 0.5 * (&f.slice(s![.., 0..(n - 1)]) + &f.slice(s![.., 1..n]) - phi);

    //Final flux differerence
    let df: Array2<f64> = &phi.slice(s![.., 1..(-1)]) - &phi.slice(s![.., 0..(-2)]);

    df
}

/// Left pressure = 1, right pressure = 0.1.
/// Left density = 1, right density = 0.125.
/// Velocity is 0 everywhere
fn sods_problem() -> (Array1<f64>, Array1<f64>, Array1<f64>) {
    println!("Configuration 1: Sod's problem.");

    let mut rho0: Array1<f64> = Array1::zeros(N_BOUNDARIES);
    let u0: Array1<f64> = Array1::zeros(N_BOUNDARIES);
    let mut p0: Array1<f64> = Array1::zeros(N_BOUNDARIES);
    let half = (0.5 * N_CELLS) as usize;

    //left
    rho0.slice_mut(s![0..half]).fill(1.0);
    p0.slice_mut(s![0..half]).fill(1.0);

    //right
    rho0.slice_mut(s![half..]).fill(0.125);
    p0.slice_mut(s![half..]).fill(0.1);

    (rho0, u0, p0)
}

fn main() {
    unsafe {
        env::set_var("RUST_BACKTRACE", "1");
    }

    let (rho0, u0, p0) = sods_problem();

    //initial total energy
    let e_tot0: Array1<f64> = &p0 / ((GAMMA - 1.0) * &rho0) + 0.5 * &u0 * &u0;

    //initial speed of sound
    let a0: Array1<f64> = (GAMMA * &p0 / &rho0).sqrt();

    //time step
    let mut dt = COURANT * DX / ((&u0).abs() + &a0).iter().fold(0.0_f64, |a, &b| a.max(b));

    //construct conservative state vector
    let mut q: Array2<f64> = Array2::zeros((3, N_BOUNDARIES));
    q.row_mut(0).assign(&rho0);
    q.row_mut(1).assign(&(&rho0 * u0));
    q.row_mut(2).assign(&(rho0 * e_tot0));

    let mut t: f64 = 0.0;
    let mut it: u32 = 0;

    let x: Vec<f64> = linspace(DX / 2., 1.0, N_BOUNDARIES).into_iter().collect();

    println!("Beginning Simulation:");

    while t < T_END {
        //println!("Iteration: {}, t = {}", it, t);

        //copy old solution
        let q0: Array2<f64> = q.clone();

        //calculate change in flux (flux divergence)
        let df: Array2<f64> = roe_flux(&q0, GAMMA, N_BOUNDARIES);

        //finite volume update (not that the boundaries (0 and -1) are unchanged)
        q.slice_mut(s![.., 1..(-2)])
            .assign(&(q0.slice(s![.., 1..(-2)]).to_owned() - (dt / DX) * df));

        //update timestep
        let (_, rho_inv, u, _, p) = decode_state(&q, GAMMA);
        let a: Array1<f64> = (GAMMA * &p * &rho_inv).sqrt();
        dt = COURANT * DX / ((&u).abs() + &a).iter().fold(0.0_f64, |a, &b| a.max(b));

        t += dt;
        it += 1;

        //display pressure along the pipe
        println!("Iteration {}", it);
        let points: Vec<(f32, f32)> = x
            .iter()
            .copied()
            .map(|y| y as f32)
            .zip(p.iter().copied().map(|y| y as f32))
            .collect();
        Chart::new_with_y_range(100, 100, 0.0, 1.0, 0.0, 1.0)
            .lineplot(&Shape::Points(points.as_slice()))
            .display();
    }

    println!("Done");
}
