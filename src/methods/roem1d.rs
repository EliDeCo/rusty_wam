// This impliments a 1D version of the RoeM2 scheme from the following paper,
// with the exception of the f and g functions which only give benefits in higher dimensions
// https://doi.org/10.1016/S0021-9991(02)00037-2

use crate::pipes::{InteriorSolver, PipeState};
use nalgebra::{Matrix3, Matrix3x1};

pub struct RoeM1D(PipeState);

impl RoeM1D {
    pub(crate) fn new(state: PipeState) -> Self {
        Self(state)
    }
}

impl InteriorSolver for RoeM1D {
    fn state(&self) -> &PipeState {
        &self.0
    }
    fn state_mut(&mut self) -> &mut PipeState {
        &mut self.0
    }

    ///Calculates the RoeM flux at every interface, then differences it into df.
    fn flux_divergence(&mut self) {
        self.0.euler_flux();

        //copy the scalars out first so the buffers below can be split-borrowed
        let gamma = self.0.gamma;
        let n_cells = self.0.n_cells;
        let PipeState {
            q0,
            phi,
            f,
            df,
            rho,
            u,
            p,
            h,
            ..
        } = &mut self.0;

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
                q0[(0, i + 1)] - q0[(0, i)], // d(rho)
                q0[(1, i + 1)] - q0[(1, i)], // d(rho*u)
                q0[(2, i + 1)] - q0[(2, i)], // d(rho*E)
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
}
