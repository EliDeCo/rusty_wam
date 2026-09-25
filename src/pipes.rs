use crate::boundaries::BoundaryCondition;
use crate::helpers::unphysical;
use crate::pipe_methods::{muscl_roem1d::MusclRoeM1D, roe1d::Roe1D, roem1d::RoeM1D};
use nalgebra::{Matrix1xX, Matrix3xX, Vector3};
use std::ops::AddAssign;

///What a junction supplies to one end of a pipe for a single stage.
/// The driver resolves both from the junction's stage state, so nothing here knows of junctions.
#[derive(Clone, Copy, Default)]
pub struct BoundaryData {
    ///numerical flux forced onto that end's face
    pub flux: Option<Vector3<f64>>,
    ///conservative state that end's ghost cells are filled with
    pub ghost: Option<Vector3<f64>>,
}

///The boundary data for both ends of one pipe.
#[derive(Clone, Copy, Default)]
pub struct BoundaryPair {
    pub left: BoundaryData,
    pub right: BoundaryData,
}

///Refills the ghost cells on each end with the state the driver resolved there.
/// An end with nothing resolved is left untouched, which holds it fixed.
pub(crate) fn apply_bc(q: &mut Matrix3xX<f64>, n_ghost: usize, bc: &BoundaryPair) {
    let n_total = q.ncols();

    if let Some(fill) = bc.left.ghost {
        for g in 0..n_ghost {
            q.set_column(g, &fill);
        }
    }

    if let Some(fill) = bc.right.ghost {
        for g in 0..n_ghost {
            q.set_column(n_total - 1 - g, &fill);
        }
    }
}

///Buffers and parameters shared by every interior method, allocated once each.
/// The grid is padded, so the real cells occupy `first .. first + n_real`.
pub struct PipeState {
    //conservative state, 3 values per cell
    pub(crate) q1: Matrix3xX<f64>, //n_total, working (current) state

    //decoded primitives of q1, one value per cell (ghosts included).
    //`decode` refreshes these; `advance` invalidates them by moving q1 underneath.
    ///kg/m3
    pub(crate) rho: Matrix1xX<f64>,
    ///m/s
    pub(crate) u: Matrix1xX<f64>,
    ///J/kg
    pub(crate) e: Matrix1xX<f64>,
    ///Pa
    pub(crate) p: Matrix1xX<f64>,
    ///J/kg
    pub(crate) h: Matrix1xX<f64>,
    ///m/s
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
    ///radius of the pipe in meters
    pub(crate) r: f64,
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
        id: usize,
        radius: f64,
        left_bc: Option<BoundaryCondition>,
        right_bc: Option<BoundaryCondition>,
    ) -> Self {
        let first = n_ghost;
        let n_real = state.ncols();
        let n_total = n_real + 2 * n_ghost;
        let n_faces = n_real + 1;

        //pad the initial condition into the full array; the driver fills the ghosts
        //before the first residual, so nothing here has to guess at them
        let mut q1: Matrix3xX<f64> = Matrix3xX::zeros(n_total);
        q1.columns_mut(first, n_real).copy_from(&state);

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
            q1,
            n_real,
            n_ghost,
            first,
            n_total,
            n_faces,
            courant,
            dx,
            r: radius,
            id,
            gamma,
            left_bc,
            right_bc,
        }
    }

    ///Decodes the pipe's own current state into the shared primitive buffers.
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
        decode_into(q1, rho, u, e, p, h, *gamma);
    }

    ///Forces any supplied flux onto the two boundary faces, then differences phi into df.
    /// Every interior method fills phi and then ends here, so injection has one home.
    pub(crate) fn difference_flux(&mut self, bc: &BoundaryPair) {
        if let Some(flux) = bc.left.flux {
            self.phi.set_column(0, &flux);
        }
        if let Some(flux) = bc.right.flux {
            self.phi.set_column(self.n_real, &flux);
        }

        let Self {
            phi, df, n_real, ..
        } = self;
        phi.columns(1, *n_real).sub_to(&phi.columns(0, *n_real), df);
    }

    ///Decodes an externally held state into the shared primitive buffers.
    pub(crate) fn decode_from(&mut self, q: &Matrix3xX<f64>) {
        let Self {
            rho,
            u,
            e,
            p,
            h,
            gamma,
            ..
        } = self;
        decode_into(q, rho, u, e, p, h, *gamma);
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

    ///Returns the minimum dt for this pipe.
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

    ///Checks for a density or pressure in the real cells that is not physically
    /// usable, which indicates a numerical blowup.
    pub fn unreal_check(&self) {
        let real = self.first..self.first + self.n_real;
        if self.rho.as_slice()[real.clone()]
            .iter()
            .copied()
            .any(unphysical)
            || self.p.as_slice()[real].iter().copied().any(unphysical)
        {
            panic!("Unphysical state in pipe {}", self.id);
        }
    }

    ///The real-cell window of a full-width per-cell buffer.
    /// Matrix1xX is a single contiguous row, so this is a plain slice.
    pub(crate) fn real<'a>(&self, m: &'a Matrix1xX<f64>) -> &'a [f64] {
        &m.as_slice()[self.first..self.first + self.n_real]
    }

    ///Returns the volume of a single cell in this pipe
    pub fn cell_volume(&self) -> f64 {
        self.pipe_area() * self.dx
    }

    ///Returns the surface area of the pipe cross-sections
    pub fn pipe_area(&self) -> f64 {
        self.r * self.r * std::f64::consts::PI
    }
}

///One interior method. Implementors supply only the spatial discretization;
/// everything method-independent lives on PipeState.
pub trait InteriorSolver {
    ///Shared buffers and parameters
    fn state(&self) -> &PipeState;
    fn state_mut(&mut self) -> &mut PipeState;

    ///Fills state.df with the raw flux difference phi[k+1] - phi[k] for the state q.
    /// The driver supplies the -(dt/dx) scaling, so every method writes df alike.
    fn residual(&mut self, q: &Matrix3xX<f64>, bc: &BoundaryPair);
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

///Which explicit SSP Runge-Kutta scheme the driver advances everything with.
/// Every variant has SSP coefficient 1, so one CFL number covers all of them.
#[derive(Clone, Copy)]
#[allow(dead_code)] //variants are selected by editing main.rs
pub enum TimeIntegrator {
    Euler,
    Ssp2,
    Ssp3,
}

impl TimeIntegrator {
    ///Blend weights (a, b) for each stage in turn, where a + b is always 1.
    pub(crate) fn stages(&self) -> &'static [(f64, f64)] {
        match self {
            TimeIntegrator::Euler => &[(0.0, 1.0)],
            TimeIntegrator::Ssp2 => &[(0.0, 1.0), (0.5, 0.5)],
            TimeIntegrator::Ssp3 => &[(0.0, 1.0), (0.75, 0.25), (1.0 / 3.0, 2.0 / 3.0)],
        }
    }

    ///Scratch registers the driver holds per pipe, on top of the q^n register.
    pub(crate) fn n_registers(&self) -> usize {
        self.stages().len() - 1
    }
}

///One integrator stage in place: dst = a*q_n + b*(src + dt*R), where R = -df/dx.
/// Doing the blend in place keeps the time loop allocation free.
#[allow(clippy::too_many_arguments)]
pub(crate) fn rk_stage(
    dst: &mut Matrix3xX<f64>,
    q_n: &Matrix3xX<f64>,
    src: &Matrix3xX<f64>,
    df: &Matrix3xX<f64>,
    a: f64,
    b: f64,
    dt_over_dx: f64,
    first: usize,
    n_real: usize,
) {
    debug_assert!((a + b - 1.0).abs() < 1e-15, "stage weights must sum to 1");

    let mut d = dst.columns_mut(first, n_real);

    // d = src + dt*R
    d.copy_from(&src.columns(first, n_real));
    d.zip_apply(df, |x, dfv| *x -= dt_over_dx * dfv);

    // d = a*q_n + b*d   (stage 1 is a=0, b=1, so skip the blend entirely)
    if a != 0.0 {
        d *= b;
        d.zip_apply(&q_n.columns(first, n_real), |x, qn| *x += a * qn);
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
    /// # Units
    /// * `state[0]` - mass density in kg/m^3
    /// * `state[1]` - momentum desnity in kg/(m^2-s)
    /// * `state[2]` - energy density in J/m^3
    /// * `radius` - pipe radius in mm
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        kind: MethodKind,
        state: Matrix3xX<f64>,
        gamma: f64,
        courant: f64,
        dx: f64,
        id: usize,
        radius: f64,
        left_bc: Option<BoundaryCondition>,
        right_bc: Option<BoundaryCondition>,
    ) -> Self {
        let shared = PipeState::new(
            state,
            kind.n_ghost(),
            gamma,
            courant,
            dx,
            id,
            radius * 0.001, //mm to m
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
    pub fn solver(&self) -> &dyn InteriorSolver {
        match self {
            Self::RoeM1D(s) => s,
            Self::Roe1D(s) => s,
            Self::MusclRoeM1D(s) => s,
        }
    }
    pub(crate) fn solver_mut(&mut self) -> &mut dyn InteriorSolver {
        match self {
            Self::RoeM1D(s) => s,
            Self::Roe1D(s) => s,
            Self::MusclRoeM1D(s) => s,
        }
    }

    ///Returns the dt for this pipe
    pub fn get_timestep(&mut self) -> f64 {
        self.solver_mut().state_mut().get_timestep()
    }
    pub fn unreal_check(&self) {
        self.solver().state().unreal_check();
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

///Decodes a conservative state into density, velocity, energy, pressure and enthalpy.
/// One value per cell, ghosts included, since the face loops read primitives there too.
#[allow(clippy::too_many_arguments)]
fn decode_into(
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
    u.component_div_assign(rho); // u = (rho*u)/rho, in place

    e.copy_from(&q.row(2)); // specific total energy, NOT specific internal energy
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
