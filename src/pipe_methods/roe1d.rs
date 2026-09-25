// This implimentation is based on P.S. Volpani's "07_1D_Euler_equations_Roe" example
//https://github.com/psvolpiani/YouTube-CFD-101

use crate::pipes::{BoundaryPair, InteriorSolver, PipeState};
use nalgebra::{Matrix3, Matrix3x1, Matrix3xX};

pub struct Roe1D(PipeState);

impl Roe1D {
    pub(crate) fn new(state: PipeState) -> Self {
        Self(state)
    }
}

impl InteriorSolver for Roe1D {
    fn state(&self) -> &PipeState {
        &self.0
    }
    fn state_mut(&mut self) -> &mut PipeState {
        &mut self.0
    }

    ///Calculates the Roe1D flux at every interface, then differences it into df.
    fn residual(&mut self, q: &Matrix3xX<f64>, bc: &BoundaryPair) {
        self.0.decode_from(q);
        self.0.euler_flux();

        //copy the scalars out first so the buffers below can be split-borrowed
        let gamma = self.0.gamma;
        let first = self.0.first;
        let PipeState {
            phi,
            f,
            rho,
            u,
            h,
            ..
        } = &mut self.0;

        //loop over each cell interface (column in phi)
        phi.column_iter_mut().enumerate().for_each(|(i, mut col)| {
            //face i sits between cells il and ir
            let il = first - 1 + i;
            let ir = first + i;

            //intermediate quanties
            let r = (rho[ir] / rho[il]).sqrt();

            //Roe averages
            //let roe_rho = r * rho[il]; // Roe average density
            let roe_u = (r * u[ir] + u[il]) / (r + 1.0); // Roe average velocity
            let half_roe_u_squared = 0.5 * roe_u * roe_u; //intermediate quantity
            let roe_h = (r * h[ir] + h[il]) / (r + 1.0); // Roe average specific total enthalpy
            let roe_a = ((gamma - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound

            //Eigenvalues
            let lambda: [f64; 3] = [roe_u - roe_a, roe_u, roe_u + roe_a];

            //difference between neighboring cell states
            let dq: Matrix3x1<f64> = Matrix3x1::new(
                q[(0, ir)] - q[(0, il)], // d(rho)
                q[(1, ir)] - q[(1, il)], // d(rho*u)
                q[(2, ir)] - q[(2, il)], // d(rho*E)
            );

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
            w[0] *= lambda[0].abs();
            w[1] *= lambda[1].abs();
            w[2] *= lambda[2].abs();

            // Transform back to physical space: upwind dissipation = P * |Lambda| * P^-1 * dq
            let dissipation: Matrix3x1<f64> = p_matrix * w;

            // Roe flux at the interface: average of the two physical fluxes minus half the upwind dissipation
            col.copy_from(&(0.5 * (f.column(il) + f.column(ir)) - 0.5 * dissipation));
        });

        self.0.difference_flux(bc);
    }
}
