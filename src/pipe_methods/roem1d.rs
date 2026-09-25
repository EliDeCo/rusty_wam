// This impliments a 1D RoeM flux from the following paper: Eq 27a-c with the signal
// velocities of Eq 33. The f and g of Eq 17b and 20b are omitted, since testing against
// exact solutions showed they only act on captured discontinuities in 1D, and the shock
// instability they exist to cure is a multi-dimensional one. See validation/RoeM1D.md
// https://doi.org/10.1016/S0021-9991(02)00037-2

use crate::pipes::{BoundaryPair, InteriorSolver, PipeState};
use nalgebra::{Matrix3, Matrix3x1, Matrix3xX};

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
            p,
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
            let roe_rho = r * rho[il]; // Roe average density
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

            //RoeM Changes ==================================================
            let u_l = u[il]; // left veloctity
            let u_r = u[ir]; // right velocity

            //Eq 33: the signal velocities take the common speed of sound, which is what
            //lets a contact be captured exactly whichever side is the hotter
            let b1 = lambda[2].max((u_r + roe_a).max(0.0));
            let b2 = lambda[0].min((u_l - roe_a).min(0.0));
            let b3 = b1 + b2;
            let b4 = 2.0 * b1 * b2;
            let b5 = 1.0 / (b1 - b2);

            //other quantities
            let m_hat = roe_u / roe_a; // Roe average Mach number
            let p_l = p[il]; // left pressure
            let p_r = p[ir]; // right pressure

            // entropy-wave correction B∆Q,
            let hlle_coeff = b1 * b2 * b5;
            let bdq_coeff = hlle_coeff / (1.0 + m_hat.abs()); // full prefactor

            let dp = p_r - p_l; // Δp

            let dh = h[ir] - h[il]; // ΔH

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
                &(0.5 * (f.column(il) + f.column(ir)) - 0.5 * dissipation + enthalpy_shift
                    - correction),
            );
        });

        self.0.difference_flux(bc);
    }
}
