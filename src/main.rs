use nalgebra::{Matrix1xX, Matrix3, Matrix3x1, Matrix3xX};
use ndarray::linspace;
use std::{env, ops::AddAssign};
use textplots::{Chart, Plot, Shape};

//Parameters
const COURANT: f64 = 0.9; //CFL courant number
const GAMMA: f64 = 1.4; //ratio if specific heats
const T_END: f64 = 1.0; //how much virtual time to run the simulation
const N_CELLS: usize = 2048;
const N_INTERFACES: usize = N_CELLS - 1; //number of interfaces between cells (boundaries don't count)
const N_INTERIOR: usize = N_CELLS - 2; //number of interior (non boundary) cells
const DOMAIN_LENGTH: f64 = 1.0; //basically how long the pipe is
const DX: f64 = DOMAIN_LENGTH / N_CELLS as f64; //step size

///decodes the state vector into primitives of the conserved variables:
/// rho (density), u (velocity), e (specific total energy), p (pressure)
/// They are vectors with one value per finite-volume cell
fn decode_state(
    q: &Matrix3xX<f64>,
    rho: &mut Matrix1xX<f64>,
    u: &mut Matrix1xX<f64>,
    e: &mut Matrix1xX<f64>,
    p: &mut Matrix1xX<f64>,
) {
    rho.copy_from(&q.row(0));

    u.copy_from(&q.row(1)); //velocity
    u.component_div_assign(&rho); // u = (rho*u)/rho, in place

    e.copy_from(&q.row(2)); // specific total energy, NOT specific internal energy
    e.component_div_assign(&rho); // e = (rho*E)/rho, in place

    // pressure from equation of state
    //done step by step to avoid allocation
    p.copy_from(u);
    p.component_mul_assign(u); // p = u*u
    *p *= -0.5; // p = -0.5*u*u
    p.add_assign(&*e); // p = e - 0.5*u*u   (reborrow e as &Matrix1xX)
    p.component_mul_assign(rho); // p = rho*(e - 0.5*u*u)
    *p *= GAMMA - 1.0; // p = (γ-1)*rho*(e - 0.5*u*u)
}

///Calculates the Euler flux (F) for every cell given state vector q and specific heat ratio gamma.
/// q = [rho, rho*u, rho*E] where rho is the density, u is the velocity, and E is the total specific energy.
/// F = [rho*u, rho*u^2 + p, u*(rho*E + p)] where p is the pressure calculated from the equation of state.
fn euler_flux(
    rho: &mut Matrix1xX<f64>,
    u: &mut Matrix1xX<f64>,
    e: &mut Matrix1xX<f64>,
    p: &mut Matrix1xX<f64>,
    f: &mut Matrix3xX<f64>,
) {

    for i in 0..3 {
        for j in 0..N_CELLS {
            f[(i, j)] = match i {
                0 => rho[j] * u[j],                 // mass flux
                1 => rho[j] * u[j] * u[j] + p[j],   // momentum flux
                2 => u[j] * (rho[j] * e[j] + p[j]), // energy flux
                _ => panic!("Invalid index for flux calculation"),
            }
        }
    }
}
///Calculates the Roe flux for every cell given the state vector q and specific heat ratio gamma.
fn roe_flux(
    q: &Matrix3xX<f64>,
    phi: &mut Matrix3xX<f64>,
    f: &mut Matrix3xX<f64>,
    rho: &mut Matrix1xX<f64>,
    u: &mut Matrix1xX<f64>,
    e: &mut Matrix1xX<f64>,
    p: &mut Matrix1xX<f64>,
    df: &mut Matrix3xX<f64>,
    h: &mut Matrix1xX<f64>,
) {
    decode_state(q, rho, u, e, p);
    // specific total enthalpy
    // computed in steps to avoid extra allocation
    h.copy_from(p);
    h.component_div_assign(rho); // h = p/rho
    h.add_assign(&*e); // h = e + p/rho 

    //get physical (euler) flux
    euler_flux(rho, u, e, p, f);

    //loop over each cell interface (column in phi)
    phi.column_iter_mut().enumerate().for_each(|(i, mut col)| {
        //intermediate quanties
        let r = (rho[i + 1] / rho[i]).sqrt();

        //Roe averages
        //let roe_rho = r * rho[i]; // Roe average density, unused in flux calculation
        let roe_u = (r * u[i + 1] + u[i]) / (r + 1.0); // Roe average velocity
        let half_roe_u_squared = 0.5 * roe_u * roe_u; //intermediate quantity
        let roe_h = (r * h[i + 1] + h[i]) / (r + 1.0); // Roe average specific total enthalpy
        let roe_a = ((GAMMA - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound

        //difference between neighboring cell states
        let dq: Matrix3x1<f64> = Matrix3x1::new(
            q[(0, i + 1)] - q[(0, i)], // d(rho)
            q[(1, i + 1)] - q[(1, i)], // d(rho*u)
            q[(2, i + 1)] - q[(2, i)], // d(rho*E)
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
        // Project dq into characteristic (wave) space: w = P^-1 * dq
        let mut w: Matrix3x1<f64> = p_matrix_inv * dq;

        //technically multiplying by diagonal wave speed matrix
        w[0] *= lambda[0].abs();
        w[1] *= lambda[1].abs();
        w[2] *= lambda[2].abs();

        // Transform back to physical space: upwind dissipation = P * |Lambda| * P^-1 * dq
        let dissipation: Matrix3x1<f64> = p_matrix * w;

        // Roe flux at the interface: average of the two physical fluxes, minus half the upwind dissipation
        col.copy_from(&(0.5 * (f.column(i) + f.column(i + 1)) - 0.5 * dissipation));
    });

    // Flux divergence per interior cell: phi_{i+1/2} - phi_{i-1/2}
    phi.columns(1, N_INTERIOR).sub_to(&phi.columns(0, N_INTERIOR), df);
}


fn max_wave_speed(u: &Matrix1xX<f64>, a: &Matrix1xX<f64>) -> f64 {
    u.iter()
        .zip(a.iter())
        .fold(0.0_f64, |speed, (&ui, &ai)| speed.max(ui.abs() + ai))
}

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

fn main() {
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }

    let (rho0, u0, p0) = sods_problem();

    //initial total energy
    let e_tot0 = p0.component_div(&((GAMMA - 1.0) * &rho0)) + 0.5 * u0.component_mul(&u0);

    //initial speed of sound
    let a0: Matrix1xX<f64> = (GAMMA * p0.component_div(&rho0)).map(|x| x.sqrt());

    //time step
    let mut dt = COURANT * DX / max_wave_speed(&u0, &a0);

    //construct conservative state vector (past copy)
    let mut q0: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
    q0.set_row(0, &rho0);
    q0.set_row(1, &(&rho0.component_mul(&u0)));
    q0.set_row(2, &(&rho0.component_mul(&e_tot0)));

    //working copy
    let mut q1: Matrix3xX<f64> = q0.to_owned();

    let mut t: f64 = 0.0;
    let mut it: u32 = 0;

    let x: Vec<f64> = linspace(DX / 2., 1.0, N_CELLS).into_iter().collect();

    //assigned only once to avoid repeated allocation
    let mut df: Matrix3xX<f64> = Matrix3xX::zeros(N_INTERIOR);
    let mut f: Matrix3xX<f64> = Matrix3xX::zeros(N_CELLS);
    let mut a: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut phi: Matrix3xX<f64> = Matrix3xX::zeros(N_INTERFACES);
    let mut rho: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut u: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut e: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut p: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);
    let mut h: Matrix1xX<f64> = Matrix1xX::zeros(N_CELLS);

    println!("Beginning Simulation:");

    while t < T_END {
        //println!("Iteration: {}, t = {}", it, t);

        //move current solution to previous solution buffer
        q0.copy_from(&q1);

        //calculate change in flux (flux divergence)
        roe_flux(
            &q0, &mut phi, &mut f, &mut rho, &mut u, &mut e, &mut p, &mut df, &mut h,
        );

        //finite volume update (the first and last cells are fixed boundary cells)
        df *= -(dt / DX);
        q0.columns(1, N_INTERIOR)
            .add_to(&df, &mut q1.columns_mut(1, N_INTERIOR));

        t += dt;
        it += 1;

        //update timestep
        decode_state(&q1, &mut rho, &mut u, &mut e, &mut p);

        // a = sqrt(gamma * p / rho), computed in place to avoid allocation
        a.copy_from(&p);
        a.component_div_assign(&rho); // a = p/rho
        a *= GAMMA;                   // a = gamma*p/rho
        for x in a.iter_mut() {
            *x = x.sqrt();              // a = sqrt(gamma*p/rho)
        }

        dt = COURANT * DX / max_wave_speed(&u, &a);
        dt = dt.min(T_END - t);

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
