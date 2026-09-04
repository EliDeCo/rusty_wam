use crate::methods::roem1d::RoeM1D;
use nalgebra::{Matrix1xX, Matrix3xX};
use std::ops::AddAssign;

///Buffers and parameters shared by every interior method.
/// Every buffer is allocated once here and written in place afterwards, so the
/// simulation loop performs no heap allocation.
pub struct PipeState {
    //conservative state, 3 values per cell
    pub(crate) q0: Matrix3xX<f64>, //n_cells, frozen base for the current step
    pub(crate) q1: Matrix3xX<f64>, //n_cells, working (current) state

    //decoded primitives of q1, one value per cell.
    //`decode` refreshes these; `advance` invalidates them by moving q1 underneath.
    pub(crate) rho: Matrix1xX<f64>,
    pub(crate) u: Matrix1xX<f64>,
    pub(crate) e: Matrix1xX<f64>,
    pub(crate) p: Matrix1xX<f64>,
    pub(crate) h: Matrix1xX<f64>,
    pub(crate) a: Matrix1xX<f64>, //speed of sound, only needed for the CFL condition

    //flux workspace
    pub(crate) f: Matrix3xX<f64>, //physical (euler) flux, n_cells
    pub(crate) phi: Matrix3xX<f64>, //numerical flux at each interface, n_cells - 1
    pub(crate) df: Matrix3xX<f64>, //flux divergence per interior cell, n_interior

    //parameters
    pub(crate) gamma: f64,
    pub(crate) courant: f64,
    pub(crate) dx: f64,
    pub(crate) n_cells: usize,
    pub(crate) n_interior: usize,
    pub(crate) id: usize,
}

impl PipeState {
    ///Allocates every shared buffer from the initial conservative state vector
    pub(crate) fn new(
        state: Matrix3xX<f64>,
        gamma: f64,
        courant: f64,
        dx: f64,
        n_cells: usize,
        id: usize,
    ) -> Self {
        let placeholder = Matrix1xX::zeros(n_cells);
        let n_interior = n_cells - 2; //number of interior (non boundary) cells

        Self {
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

    ///decodes the current state vector (q1) into primitives of the conserved variables:
    /// rho (density), u (velocity), e (specific total energy), p (pressure),
    /// h (specific total enthalpy).
    /// They are vectors with one value per finite-volume cell
    pub(crate) fn decode(&mut self) {
        let Self {
            q1,
            rho,
            u,
            e,
            p,
            h,
            gamma,
            ..
        } = self;
        let gamma = *gamma;

        rho.copy_from(&q1.row(0));

        u.copy_from(&q1.row(1)); //velocity
        u.component_div_assign(rho); // u = (rho*u)/rho, in place

        e.copy_from(&q1.row(2)); // specific total energy, NOT specific internal energy
        e.component_div_assign(rho); // e = (rho*E)/rho, in place

        // pressure from equation of state
        //done step by step to avoid allocation
        p.copy_from(u);
        p.component_mul_assign(u); // p = u*u
        *p *= -0.5; // p = -0.5*u*u
        p.add_assign(&*e); // p = e - 0.5*u*u   
        p.component_mul_assign(rho); // p = rho*(e - 0.5*u*u)
        *p *= gamma - 1.0; // p = (γ-1)*rho*(e - 0.5*u*u)

        // specific total enthalpy
        // computed in steps to avoid extra allocation
        h.copy_from(p);
        h.component_div_assign(rho); // h = p/rho
        h.add_assign(&*e); // h = e + p/rho
    }

    ///Calculates the Euler flux (F) for every cell from the decoded primitives.
    /// q = [rho, rho*u, rho*E] where rho is the density, u is the velocity, and E is the total specific energy.
    /// F = [rho*u, rho*u^2 + p, u*(rho*E + p)] where p is the pressure calculated from the equation of state.
    pub(crate) fn euler_flux(&mut self) {
        let Self {
            rho,
            u,
            e,
            p,
            f,
            n_cells,
            ..
        } = self;

        for i in 0..3 {
            for j in 0..*n_cells {
                f[(i, j)] = match i {
                    0 => rho[j] * u[j],                 // mass flux
                    1 => rho[j] * u[j] * u[j] + p[j],   // momentum flux
                    2 => u[j] * (rho[j] * e[j] + p[j]), // energy flux
                    _ => panic!("Invalid index for flux calculation"),
                }
            }
        }
    }

    ///Returns the dt for this pipe.
    pub fn get_timestep(&mut self) -> f64 {
        self.decode();

        let Self {
            a,
            rho,
            u,
            p,
            courant,
            gamma,
            dx,
            ..
        } = self;

        // a = sqrt(gamma * p / rho), computed in place to avoid allocation
        a.copy_from(p);
        a.component_div_assign(rho); // a = p/rho
        *a *= *gamma; // a = gamma*p/rho
        for x in a.iter_mut() {
            *x = x.sqrt(); // a = sqrt(gamma*p/rho)
        }

        *courant * *dx / max_wave_speed(u, a)
    }

    ///Checks for negative density or pressure, which indicate a numerical blowup.
    pub fn nan_check(&self) {
        if self.rho.iter().any(|&x| x < 0.0) || self.p.iter().any(|&x| x < 0.0) {
            panic!("Nan in pipe: {}", self.id);
        }
    }

    ///Moves current solution to previous solution buffer
    pub(crate) fn save_step(&mut self) {
        self.q0.copy_from(&self.q1);
    }

    ///Finite volume update: q1 = q0 - (dt/dx)*df.
    /// The first and last cells are fixed boundary cells and are left untouched.
    /// Leaves the decoded primitives stale, since q1 has moved.
    pub(crate) fn advance(&mut self, dt: f64) {
        let Self {
            df,
            q0,
            q1,
            n_interior,
            dx,
            ..
        } = self;

        *df *= -(dt / *dx);
        q0.columns(1, *n_interior)
            .add_to(&*df, &mut q1.columns_mut(1, *n_interior));
    }
}

///One interior method. Implementors supply only the spatial discretization;
/// everything method-independent lives on PipeState.
pub trait InteriorSolver {
    ///Shared buffers and parameters
    fn state(&self) -> &PipeState;
    fn state_mut(&mut self) -> &mut PipeState;

    ///Fills state.df from the current state.
    /// This is mainly the only part that differs between interior methods.
    ///
    /// Assumes the decoded primitives are already current for q1, so implementors
    /// do not decode by default. An override calling this more than once per step 
    /// (like RK3) MUST `state_mut().decode()` after each `advance`
    fn flux_divergence(&mut self);

    ///Advances the interior state forward in time by dt.
    /// Override this only when the *time* integration differs from a
    /// single forward euler stage
    fn update(&mut self, dt: f64) {
        self.state_mut().save_step();
        self.flux_divergence();
        self.state_mut().advance(dt);
    }
}

///Selects which interior method InteriorMethod::new constructs
pub enum MethodKind {
    RoeM1D,
}

///Owns one interior method and dispatches to it
pub enum InteriorMethod {
    RoeM1D(RoeM1D),
}

impl InteriorMethod {
    ///Returns a new InteriorMethod instance from the conservative state vector and other parameters
    pub fn new(
        kind: MethodKind,
        state: Matrix3xX<f64>,
        gamma: f64,
        courant: f64,
        dx: f64,
        n_cells: usize,
        id: usize,
    ) -> Self {
        let shared = PipeState::new(state, gamma, courant, dx, n_cells, id);

        match kind {
            MethodKind::RoeM1D => Self::RoeM1D(RoeM1D::new(shared)),
        }
    }

    //The only per-variant match arms exist here
    fn solver(&self) -> &dyn InteriorSolver {
        match self {
            Self::RoeM1D(s) => s,
        }
    }
    fn solver_mut(&mut self) -> &mut dyn InteriorSolver {
        match self {
            Self::RoeM1D(s) => s,
        }
    }

    ///Updates the interior state forward in time
    pub fn update(&mut self, dt: f64) {
        self.solver_mut().update(dt);
    }
    
    ///Returns the dt for this pipe
    pub fn get_timestep(&mut self) -> f64 {
        self.solver_mut().state_mut().get_timestep()
    }
    pub fn nan_check(&self) {
        self.solver().state().nan_check();
    }
    pub fn rho(&self) -> &Matrix1xX<f64> {
        &self.solver().state().rho
    }
    #[allow(dead_code)]
    pub fn u(&self) -> &Matrix1xX<f64> {
        &self.solver().state().u
    }
    #[allow(dead_code)]
    pub fn p(&self) -> &Matrix1xX<f64> {
        &self.solver().state().p
    }
}

fn max_wave_speed(u: &Matrix1xX<f64>, a: &Matrix1xX<f64>) -> f64 {
    u.iter()
        .zip(a.iter())
        .fold(0.0_f64, |speed, (&ui, &ai)| speed.max(ui.abs() + ai))
}
