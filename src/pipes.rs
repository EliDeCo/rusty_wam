use nalgebra::{Matrix1xX, Matrix3, Matrix3x1, Matrix3xX};
use std::ops::AddAssign;

pub enum InteriorMethod {
    RoeM1D {
        df: Matrix3xX<f64>,
        f: Matrix3xX<f64>,
        a: Matrix1xX<f64>,
        phi: Matrix3xX<f64>,
        rho: Matrix1xX<f64>,
        u: Matrix1xX<f64>,
        e: Matrix1xX<f64>,
        p: Matrix1xX<f64>,
        h: Matrix1xX<f64>,
        q0: Matrix3xX<f64>,
        q1: Matrix3xX<f64>,
        courant: f64,
        gamma: f64,
        n_cells: usize,
        n_interior: usize,
        dx: f64,
        id: usize,
    },
}
impl InteriorMethod {
    ///Returns a new InteriorMethod instance from the conservative state vector and other parameters
    pub fn new(
        method: &str,
        state: Matrix3xX<f64>,
        gamma: f64,
        courant: f64,
        dx: f64,
        n_cells: usize,
        id: usize,
    ) -> Self {
        return match method {
            "RoeM1D" => {
                let placeholder = Matrix1xX::zeros(n_cells);
                let n_interior = n_cells - 2; //number of interior (non boundary) cells

                Self::RoeM1D {
                    df: Matrix3xX::zeros(n_interior),
                    f: Matrix3xX::zeros(n_cells),
                    a: placeholder.clone(),
                    phi: Matrix3xX::zeros(n_cells - 1), //number of interfaces between cells (boundaries don't count)
                    rho: placeholder.clone(),
                    u: placeholder.clone(),
                    e: placeholder.clone(),
                    p: placeholder.clone(),
                    h: placeholder,
                    q1: state.clone(),
                    q0: state,
                    courant,
                    n_cells,
                    n_interior,
                    dx,
                    id,
                    gamma,
                }
            }
            _ => panic!("Invalid interior method: {}", method),
        };
    }
    ///Returns the dt for this pipe
    pub fn get_timestep(&mut self) -> f64 {
        match self {
            InteriorMethod::RoeM1D {
                a,
                rho,
                u,
                e,
                p,
                h,
                q1,
                courant,
                gamma,
                dx,
                ..
            } => {
                decode_state(&q1, rho, u, e, p, h, *gamma);

                // a = sqrt(gamma * p / rho), computed in place to avoid allocation
                a.copy_from(&p);
                a.component_div_assign(&rho); // a = p/rho
                *a *= *gamma; // a = gamma*p/rho
                for x in a.iter_mut() {
                    *x = x.sqrt(); // a = sqrt(gamma*p/rho)
                }

                *courant * *dx / max_wave_speed(&u, &a)
            }
        }
    }
    pub fn nan_check(&self) {
        match self {
            InteriorMethod::RoeM1D { rho, p, id, .. } => {
                if rho.iter().any(|&x| x < 0.0) || p.iter().any(|&x| x < 0.0) {
                    panic!("Nan in pipe: {}", id);
                }
            }
        }
    }
    ///Updates the interior state forward in time
    pub fn update(&mut self, dt: f64) {
        match self {
            InteriorMethod::RoeM1D {
                df,
                f,
                phi,
                rho,
                u,
                e,
                p,
                h,
                q0,
                q1,
                gamma,
                n_cells,
                n_interior,
                dx,
                ..
            } => {
                //move current solution to previous solution buffer
                q0.copy_from(&q1);

                //change in flux at each interface (df)
                roem_flux(q0, phi, f, rho, u, e, p, df, h, *n_cells, *gamma);

                //finite volume update (the first and last cells are fixed boundary cells)
                *df *= -(dt / *dx);
                q0.columns(1, *n_interior)
                    .add_to(&df, &mut q1.columns_mut(1, *n_interior));
            }
        }
    }
    pub fn rho(&self) -> &Matrix1xX<f64> {
        match self {
            InteriorMethod::RoeM1D { rho, .. } => rho,
        }
    }
    pub fn _u(&self) -> &Matrix1xX<f64> {
        match self {
            InteriorMethod::RoeM1D { u, .. } => u,
        }
    }
    pub fn _p(&self) -> &Matrix1xX<f64> {
        match self {
            InteriorMethod::RoeM1D { p, .. } => p,
        }
    }
}

fn max_wave_speed(u: &Matrix1xX<f64>, a: &Matrix1xX<f64>) -> f64 {
    u.iter()
        .zip(a.iter())
        .fold(0.0_f64, |speed, (&ui, &ai)| speed.max(ui.abs() + ai))
}

///decodes the state vector into primitives of the conserved variables:
/// rho (density), u (velocity), e (specific total energy), p (pressure)
/// They are vectors with one value per finite-volume cell
fn decode_state(
    q: &Matrix3xX<f64>,
    rho: &mut Matrix1xX<f64>,
    u: &mut Matrix1xX<f64>,
    e: &mut Matrix1xX<f64>,
    p: &mut Matrix1xX<f64>,
    h: &mut Matrix1xX<f64>,
    gamma: f64,
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
    *p *= gamma - 1.0; // p = (γ-1)*rho*(e - 0.5*u*u)

    // specific total enthalpy
    // computed in steps to avoid extra allocation
    h.copy_from(p);
    h.component_div_assign(rho); // h = p/rho
    h.add_assign(&*e); // h = e + p/rho 
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
    n_cells: usize,
) {
    for i in 0..3 {
        for j in 0..n_cells {
            f[(i, j)] = match i {
                0 => rho[j] * u[j],                 // mass flux
                1 => rho[j] * u[j] * u[j] + p[j],   // momentum flux
                2 => u[j] * (rho[j] * e[j] + p[j]), // energy flux
                _ => panic!("Invalid index for flux calculation"),
            }
        }
    }
}
///Calculates the RoeM flux for every cell given the state vector q and specific heat ratio gamma.
fn roem_flux(
    q: &Matrix3xX<f64>,
    phi: &mut Matrix3xX<f64>,
    f: &mut Matrix3xX<f64>,
    rho: &mut Matrix1xX<f64>,
    u: &mut Matrix1xX<f64>,
    e: &mut Matrix1xX<f64>,
    p: &mut Matrix1xX<f64>,
    df: &mut Matrix3xX<f64>,
    h: &mut Matrix1xX<f64>,
    n_cells: usize,
    gamma: f64,
) {
    decode_state(q, rho, u, e, p, h, gamma);

    //get physical (euler) flux
    euler_flux(rho, u, e, p, f, n_cells);

    //loop over each cell interface (column in phi)
    phi.column_iter_mut().enumerate().for_each(|(i, mut col)| {
        //intermediate quanties
        let r = (rho[i + 1] / rho[i]).sqrt();

        //Roe averages
        let roe_rho = r * rho[i]; // Roe average density
        let roe_u = (r * u[i + 1] + u[i]) / (r + 1.0); // Roe average velocity
        let half_roe_u_squared = 0.5 * roe_u * roe_u; //intermediate quantity
        let roe_h = (r * h[i + 1] + h[i]) / (r + 1.0); // Roe average specific total enthalpy
        let roe_a = ((gamma - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound

        //Eigenvalues
        let lambda: [f64; 3] = [roe_u - roe_a, roe_u, roe_u + roe_a];

        //difference between neighboring cell states
        let dq: Matrix3x1<f64> = Matrix3x1::new(
            q[(0, i + 1)] - q[(0, i)], // d(rho)
            q[(1, i + 1)] - q[(1, i)], // d(rho*u)
            q[(2, i + 1)] - q[(2, i)], // d(rho*E)
        );

        //RoeM Changes ==================================================
        let u_l = u[i]; // left veloctity
        let a_l = (gamma * p[i] / rho[i]).sqrt(); // left speed of sound
        let u_r = u[i + 1]; // right velocity
        let a_r = (gamma * p[i + 1] / rho[i + 1]).sqrt(); // right speed of sound

        //intermediates
        let b1 = lambda[2].max((u_r + a_r).max(0.0));
        let b2 = lambda[0].min((u_l - a_l).min(0.0));
        let b3 = b1 + b2;
        let b4 = 2.0 * b1 * b2;
        let b5 = 1.0 / (b1 - b2);

        //other quantities
        let m_hat = roe_u / roe_a; // Roe average Mach number
        let p_l = p[i]; // left pressure
        let p_r = p[i + 1]; // right pressure

        // entropy-wave correction B∆Q,
        let hlle_coeff = b1 * b2 * b5;
        let bdq_coeff = hlle_coeff / (1.0 + m_hat.abs()); // full prefactor

        let dp = p_r - p_l; // Δp

        let dh = h[i + 1] - h[i]; // ΔH

        let b_dq_0 = dq[0] - dp / (roe_a * roe_a); // Δρ - Δp/â²
        let b_dq: Matrix3x1<f64> =
            Matrix3x1::new(b_dq_0, b_dq_0 * roe_u, b_dq_0 * roe_h + roe_rho * dh);

        let correction: Matrix3x1<f64> = bdq_coeff * b_dq;

        // accounts for swapping ΔQ -> ΔQ* = Δ(ρ,ρu,ρH) inside the HLLE base term
        let enthalpy_shift: Matrix3x1<f64> = hlle_coeff * Matrix3x1::new(0.0, 0.0, dp);

        // ==============================================================

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
            &(0.5 * (f.column(i) + f.column(i + 1)) - 0.5 * dissipation + enthalpy_shift
                - correction),
        );
    });

    // Flux divergence per interior cell: phi_{i+1/2} - phi_{i-1/2}
    phi.columns(1, n_cells - 2)
        .sub_to(&phi.columns(0, n_cells - 2), df);
}
