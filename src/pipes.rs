use crate::pipe_methods::{muscl_roem1d::MusclRoeM1D, roe1d::Roe1D, roem1d::RoeM1D};
use nalgebra::{Matrix1xX, Matrix3xX};
use std::ops::AddAssign;

///How the ghost cells on one end of a pipe are refilled every time the state moves.
#[derive(Clone, Copy)]
pub enum BoundaryCondition {
    ///zero-gradient: waves pass through the end undisturbed
    Transmissive,
}

///Fills the ghost cells surrounding the real domain so every real cell - including the
/// ones nearest each edge - can be updated with the exact same stencil.
/// left/right are independent so e.g. a closed-end shock tube (wall + transmissive)
/// can be built from the same function. `None` leaves that end's ghosts untouched,
/// which reproduces a fixed (frozen) boundary.
///
/// Free function rather than a method because SSP-RK3 has to apply it to its stage
/// buffers, which a `&mut self` method could not reach without a borrow conflict.
pub(crate) fn apply_bc(
    q: &mut Matrix3xX<f64>,
    first: usize,
    n_real: usize,
    n_ghost: usize,
    left: Option<BoundaryCondition>,
    right: Option<BoundaryCondition>,
) {
    let n_total = q.ncols();

    if let Some(BoundaryCondition::Transmissive) = left {
        let mirror = q.column(first).into_owned();
        for g in 0..n_ghost {
            q.set_column(g, &mirror);
        }
    }

    if let Some(BoundaryCondition::Transmissive) = right {
        let mirror = q.column(first + n_real - 1).into_owned();
        for g in 0..n_ghost {
            q.set_column(n_total - 1 - g, &mirror);
        }
    }
}

///Buffers and parameters shared by every interior method.
/// Every buffer is allocated once here and written in place afterwards, so the
/// simulation loop performs no heap allocation.
///
/// The grid is padded: `n_ghost` ghost cells sit on each end of the `n_real` real
/// cells, so the real cells occupy `first .. first + n_real` in every n_total-wide
/// buffer. How many ghosts a method needs depends on its stencil width, so `n_ghost`
/// comes from `MethodKind::n_ghost`.
pub struct PipeState {
    //conservative state, 3 values per cell
    pub(crate) q0: Matrix3xX<f64>, //n_total, frozen base for the current step
    pub(crate) q1: Matrix3xX<f64>, //n_total, working (current) state

    //decoded primitives of q1, one value per cell (ghosts included).
    //`decode` refreshes these; `advance` invalidates them by moving q1 underneath.
    pub(crate) rho: Matrix1xX<f64>,
    pub(crate) u: Matrix1xX<f64>,
    pub(crate) e: Matrix1xX<f64>,
    pub(crate) p: Matrix1xX<f64>,
    pub(crate) h: Matrix1xX<f64>,
    pub(crate) a: Matrix1xX<f64>, //speed of sound, only needed for the CFL condition

    //flux workspace
    pub(crate) f: Matrix3xX<f64>,   //physical (euler) flux, n_total
    pub(crate) phi: Matrix3xX<f64>, //numerical flux at each face, n_faces
    pub(crate) df: Matrix3xX<f64>,  //flux divergence per real cell, n_real

    //geometry
    pub(crate) n_real: usize,  //cells actually advanced in time
    pub(crate) n_ghost: usize, //ghost cells on EACH end
    pub(crate) first: usize,   //index of the first real cell (== n_ghost)
    pub(crate) n_total: usize, //n_real + 2*n_ghost
    pub(crate) n_faces: usize, //faces bounding the real cells (== n_real + 1)

    //parameters
    pub(crate) gamma: f64,
    pub(crate) courant: f64,
    pub(crate) dx: f64,
    pub(crate) id: usize,
    pub(crate) left_bc: Option<BoundaryCondition>,
    pub(crate) right_bc: Option<BoundaryCondition>,
}

impl PipeState {
    ///Allocates every shared buffer, padding the (real-cell-wide) initial conservative
    /// state vector with ghost cells and filling them from the boundary conditions.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        state: Matrix3xX<f64>, //n_real wide: real cells only
        n_ghost: usize,
        gamma: f64,
        courant: f64,
        dx: f64,
        n_real: usize,
        id: usize,
        left_bc: Option<BoundaryCondition>,
        right_bc: Option<BoundaryCondition>,
    ) -> Self {
        let first = n_ghost;
        let n_total = n_real + 2 * n_ghost;
        let n_faces = n_real + 1;

        //pad the initial condition into the full array, then fill the ghosts
        let mut q1: Matrix3xX<f64> = Matrix3xX::zeros(n_total);
        q1.columns_mut(first, n_real).copy_from(&state);
        apply_bc(&mut q1, first, n_real, n_ghost, left_bc, right_bc);

        let placeholder = Matrix1xX::zeros(n_total);

        Self {
            df: Matrix3xX::zeros(n_real),
            f: Matrix3xX::zeros(n_total),
            a: placeholder.clone(),
            phi: Matrix3xX::zeros(n_faces),
            rho: placeholder.clone(),
            u: placeholder.clone(),
            e: placeholder.clone(),
            p: placeholder.clone(),
            h: placeholder,
            q0: q1.clone(),
            q1,
            n_real,
            n_ghost,
            first,
            n_total,
            n_faces,
            courant,
            dx,
            id,
            gamma,
            left_bc,
            right_bc,
        }
    }

    ///decodes the current state vector (q1) into primitives of the conserved variables:
    /// rho (density), u (velocity), e (specific total energy), p (pressure),
    /// h (specific total enthalpy).
    /// They are vectors with one value per finite-volume cell, ghosts included -
    /// the face loops need primitives in the ghosts too.
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
            n_total,
            ..
        } = self;

        for i in 0..3 {
            for j in 0..*n_total {
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
            first,
            n_real,
            ..
        } = self;

        // a = sqrt(gamma * p / rho), computed in place to avoid allocation
        a.copy_from(p);
        a.component_div_assign(rho); // a = p/rho
        *a *= *gamma; // a = gamma*p/rho
        for x in a.iter_mut() {
            *x = x.sqrt(); // a = sqrt(gamma*p/rho)
        }

        *courant * *dx / max_wave_speed(u, a, *first, *n_real)
    }

    ///Checks for negative density or pressure in the real cells, which indicate a
    /// numerical blowup.
    pub fn nan_check(&self) {
        let real = self.first..self.first + self.n_real;
        if self.rho.as_slice()[real.clone()].iter().any(|&x| x < 0.0)
            || self.p.as_slice()[real].iter().any(|&x| x < 0.0)
        {
            panic!("Nan in pipe: {}", self.id);
        }
    }

    ///Moves current solution to previous solution buffer
    pub(crate) fn save_step(&mut self) {
        self.q0.copy_from(&self.q1);
    }

    ///Finite volume update: q1 = q0 - (dt/dx)*df over the real cells, then refills
    /// the ghosts from the boundary conditions.
    /// Leaves the decoded primitives stale, since q1 has moved.
    pub(crate) fn advance(&mut self, dt: f64) {
        let Self {
            df,
            q0,
            q1,
            first,
            n_real,
            dx,
            ..
        } = self;

        *df *= -(dt / *dx);
        q0.columns(*first, *n_real)
            .add_to(&*df, &mut q1.columns_mut(*first, *n_real));

        apply_bc(
            &mut self.q1,
            self.first,
            self.n_real,
            self.n_ghost,
            self.left_bc,
            self.right_bc,
        );
    }

    ///The real-cell window of a full-width per-cell buffer.
    /// Matrix1xX is a single contiguous row, so this is a plain slice.
    pub(crate) fn real<'a>(&self, m: &'a Matrix1xX<f64>) -> &'a [f64] {
        &m.as_slice()[self.first..self.first + self.n_real]
    }
}

///One interior method. Implementors supply only the spatial discretization;
/// everything method-independent lives on PipeState.
pub trait InteriorSolver {
    ///Shared buffers and parameters
    fn state(&self) -> &PipeState;
    fn state_mut(&mut self) -> &mut PipeState;

    ///Fills state.df with the raw flux difference phi[k+1] - phi[k] over the real
    /// cells. `advance` supplies the -(dt/dx) scaling, so every method writes df in
    /// the same convention.
    /// This is mainly the only part that differs between interior methods.
    ///
    /// Assumes the decoded primitives are already current for q1, so implementors
    /// do not decode by default. An override calling this more than once per step
    /// (like RK3) MUST `state_mut().decode()` after each `advance` - unless, like
    /// MusclRoeM1D, it decodes its own face states and never reads the shared ones.
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
#[derive(Clone, Copy)]
#[allow(dead_code)] //variants are selected by editing main.rs
pub enum MethodKind {
    RoeM1D,
    Roe1D,
    MusclRoeM1D,
}

impl MethodKind {
    ///Ghost cells this method needs on each end, set by its stencil width.
    /// The first-order schemes reach one cell past each face; MUSCL reconstruction
    /// reaches two.
    pub(crate) fn n_ghost(&self) -> usize {
        match self {
            MethodKind::RoeM1D | MethodKind::Roe1D => 1,
            MethodKind::MusclRoeM1D => 2,
        }
    }
}

///Owns one interior method and dispatches to it
#[allow(clippy::large_enum_variant)]
pub enum InteriorMethod {
    RoeM1D(RoeM1D),
    Roe1D(Roe1D),
    MusclRoeM1D(MusclRoeM1D),
}

impl InteriorMethod {
    ///Returns a new InteriorMethod instance from the conservative state vector
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        kind: MethodKind,
        state: Matrix3xX<f64>,
        gamma: f64,
        courant: f64,
        dx: f64,
        n_real: usize,
        id: usize,
        left_bc: Option<BoundaryCondition>,
        right_bc: Option<BoundaryCondition>,
    ) -> Self {
        let shared = PipeState::new(
            state,
            kind.n_ghost(),
            gamma,
            courant,
            dx,
            n_real,
            id,
            left_bc,
            right_bc,
        );

        match kind {
            MethodKind::RoeM1D => Self::RoeM1D(RoeM1D::new(shared)),
            MethodKind::Roe1D => Self::Roe1D(Roe1D::new(shared)),
            MethodKind::MusclRoeM1D => Self::MusclRoeM1D(MusclRoeM1D::new(shared)),
        }
    }

    //The only per-variant match arms exist here
    fn solver(&self) -> &dyn InteriorSolver {
        match self {
            Self::RoeM1D(s) => s,
            Self::Roe1D(s) => s,
            Self::MusclRoeM1D(s) => s,
        }
    }
    fn solver_mut(&mut self) -> &mut dyn InteriorSolver {
        match self {
            Self::RoeM1D(s) => s,
            Self::Roe1D(s) => s,
            Self::MusclRoeM1D(s) => s,
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
    pub fn rho(&self) -> &[f64] {
        let s = self.solver().state();
        s.real(&s.rho)
    }
    #[allow(dead_code)]
    pub fn u(&self) -> &[f64] {
        let s = self.solver().state();
        s.real(&s.u)
    }
    #[allow(dead_code)]
    pub fn p(&self) -> &[f64] {
        let s = self.solver().state();
        s.real(&s.p)
    }
    pub fn id(&self) -> usize {
        self.solver().state().id
    }
}

///Largest signal speed over the real cells only.
fn max_wave_speed(u: &Matrix1xX<f64>, a: &Matrix1xX<f64>, first: usize, n_real: usize) -> f64 {
    let real = first..first + n_real;
    u.as_slice()[real.clone()]
        .iter()
        .zip(a.as_slice()[real].iter())
        .fold(0.0_f64, |speed, (&ui, &ai)| speed.max(ui.abs() + ai))
}
