// This impliments the ghost junction method from Hong & Kim, 2011
//https://doi.org/10.1002/fld.2212,

use crate::helpers::unphysical;
use crate::pipes::InteriorMethod;
use nalgebra::{Matrix3xX, Vector3, Vector5};
use std::collections::BTreeMap;

pub struct Junction {
    //state is now 1 vector (only 1 cell) of 5 values, which is the same
    //as the original vector but with momentum into x, y, z directions
    pub q1: Vector5<f64>, //current state

    ///outward flux sum over the attached pipes, which the driver scales by -(dt/volume)
    pub df: Vector5<f64>,

    ///volume of the junction in m^3
    pub volume: f64,

    ///sum of N*s over the attached pipes, which sets the representative wall force
    pub normal_area_sum: Vector3<f64>,

    //decoded primitives
    pub rho: f64,
    ///since there are 3 dimensions, velocity is a vector
    pub u: Vector3<f64>,
    pub p: f64,
    pub h: f64,

    //parameters
    pub gamma: f64,
    pub courant: f64,
    pub id: usize,

    //connections. A junction joins a handful of pipes at most, so a linear scan
    //beats a map lookup and nothing here needs to index by pipe id anyway.
    pub pipes: Vec<InletData>,
}

impl Junction {
    ///Decodes the junction's own current state into primitives.
    pub fn decode(&mut self) {
        self.decode_from(&self.q1.clone());
    }

    ///Decodes a conservative state into density, velocity, pressure and enthalpy.
    pub fn decode_from(&mut self, q: &Vector5<f64>) {
        self.rho = q[0];

        let without_rho = q / self.rho;

        // u = (rho*u)/rho,
        self.u = without_rho.fixed_rows::<3>(1).into_owned();

        // specific total energy, NOT specific internal energy
        // e = (rho*E)/rho
        let e = without_rho[4];

        // pressure from equation of state
        // p = (γ-1)*rho*(e - 0.5*(u_x^2 + u_y^2 + u_z^2))
        self.p = (self.gamma - 1.0) * self.rho * (e - 0.5 * self.u.dot(&self.u));

        // specific total enthalpy
        // h = e + p/rho
        self.h = e + self.p / self.rho;
    }
    ///Checks for a density or pressure that is not physically usable, which indicates
    /// a numerical blowup.
    pub fn unreal_check(&self) {
        if unphysical(self.rho) || unphysical(self.p) {
            panic!("Unphysical state in junction {}", self.id);
        }
    }

    ///Returns the minimum dt for this junction.
    pub fn get_timestep(&mut self) -> f64 {
        self.decode();

        // a = sqrt(gamma * p / rho)
        let a = (self.gamma * self.p / self.rho).sqrt();

        //the normal velocity is recomputed from the freshly decoded u rather than read
        //from InletData::u_dot, which euler_fluxes only refreshes later in the step
        let wave_sum: f64 = self
            .pipes
            .iter()
            .map(|pipe| pipe.area * (self.u.dot(&pipe.normal).abs() + a))
            .sum();

        self.courant * 2.0 * self.volume / wave_sum
    }

    ///Fills df with the outward interface flux sum plus the representative wall force.
    /// The driver stores each interface flux before this runs.
    pub fn residual(&mut self, q: &Vector5<f64>) {
        self.decode_from(q);

        self.df = self
            .pipes
            .iter()
            .fold(Vector5::zeros(), |sum, pipe| sum + pipe.f * pipe.area);

        //Eq 11: the walls left over from the pipe openings push back with the cell pressure
        let momentum = self.df.fixed_rows::<3>(1) - self.normal_area_sum * self.p;
        self.df.fixed_rows_mut::<3>(1).copy_from(&momentum);
    }

    pub fn new(gamma: f64, courant: f64, id: usize) -> Self {
        Self {
            q1: Vector5::zeros(),
            df: Vector5::zeros(),
            volume: 0.0,
            normal_area_sum: Vector3::zeros(),
            //these will get updated by the decode() call in get_timestep()
            rho: 0.0,
            u: Vector3::zeros(),
            p: 0.0,
            h: 0.0,
            gamma,
            courant,
            id,
            pipes: Vec::new(),
        }
    }

    ///Attaches a pipe to this junction based on its id, orientation, and which end
    /// the supplied normal is expected to point in, which is out of the junction
    /// center and along the axis of the pipe.
    pub fn add_pipe(&mut self, pipe_id: usize, left: bool, normal: Vector3<f64>, area: f64) {
        debug_assert!(
            !self.pipes.iter().any(|pipe| pipe.pipe_id == pipe_id),
            "Junction {}: pipe {} attached twice",
            self.id,
            pipe_id
        );

        self.pipes.push(InletData {
            pipe_id,
            normal,
            area,
            left,
            ..Default::default()
        });
    }

    ///Updates the volume and initial condition of the junction based on the attached pipes
    pub fn initialize(&mut self, all_pipes: &BTreeMap<usize, InteriorMethod>) {
        assert!(
            !self.pipes.is_empty(),
            "Junction {}: no pipes are attached",
            self.id
        );

        let n_pipes = self.pipes.len();
        let mut volumes: Vec<f64> = Vec::with_capacity(n_pipes);
        let mut states: Matrix3xX<f64> = Matrix3xX::zeros(n_pipes);
        let mut velocities: Vec<Vector3<f64>> = Vec::with_capacity(n_pipes);
        let mut pressures: Vec<f64> = Vec::with_capacity(n_pipes);
        let mut areas: Vec<f64> = Vec::with_capacity(n_pipes);

        //gather pipe data
        for (i, inlet) in self.pipes.iter().enumerate() {
            let pipe = &all_pipes[&inlet.pipe_id];
            let pipe_state = pipe.solver().state();

            //add the volume of a cell from this pipe
            volumes.push(pipe_state.cell_volume());

            //index of the real cell adjacent to this junction, in the PADDED array -
            //the ghosts shift it over by `first`. Normals point out of the junction,
            //so a pipe met at its left end runs along +normal and one met at its
            //right end runs along -normal.
            let (index, multiplier) = match inlet.left {
                true => (pipe_state.first, 1.),
                false => (pipe_state.first + pipe_state.n_real - 1, -1.),
            };

            //add the initial condition of the adjacent cell from this pipe
            states.column_mut(i).copy_from(&pipe_state.q1.column(index));

            //velocity is read straight back out of that state, so this does not
            //depend on the pipe having been decoded yet
            let u = states[(1, i)] / states[(0, i)];

            //add velocity vector in normal direction scaled by pipe's velocity
            velocities.push(inlet.normal * u * multiplier);
            pressures.push((self.gamma - 1.0) * (states[(2, i)] - 0.5 * states[(0, i)] * u * u));
            areas.push(inlet.area);
        }

        let total_vol = volumes.iter().sum::<f64>();

        //set volume to mean of adjacent cell volumes
        self.volume = total_vol / volumes.len() as f64;

        //fixed geometry, so the wall force only has to scale it by pressure each step
        self.normal_area_sum = self
            .pipes
            .iter()
            .fold(Vector3::zeros(), |sum, pipe| sum + pipe.normal * pipe.area);

        //set mass density and pressure by weighted average
        let mean = |vals: &[f64]| {
            volumes
                .iter()
                .zip(vals)
                .fold(0.0, |sum, (weight, val)| sum + val * weight)
                / total_vol
        };
        let rho_mean = mean(states.row(0).clone_owned().as_slice());
        let p_mean = mean(&pressures);

        //set velocity as a weighted average of vectors by pipe inlet surface area
        let total_surface_area = areas.iter().sum::<f64>();
        let u_mean = velocities
            .iter()
            .zip(areas)
            .fold(Vector3::zeros(), |sum, (val, weight)| sum + val * weight)
            / total_surface_area;

        //energy is rebuilt from that pressure and velocity rather than averaged straight
        //from the branches, which would leave the junction holding more heat than they do
        self.q1[0] = rho_mean;
        self.q1[4] = p_mean / (self.gamma - 1.0) + 0.5 * rho_mean * u_mean.dot(&u_mean);

        //rows 1..4 hold MOMENTUM density, so that averaged velocity is scaled by the
        //density set above before being stored
        self.q1
            .fixed_rows_mut::<3>(1)
            .copy_from(&(u_mean * rho_mean));
    }
}

///Conservative 1D ghost state a junction hands to one attached pipe end.
/// Density, axial velocity and stagnation enthalpy carry across unchanged.
pub fn ghost_state(
    q: &Vector5<f64>,
    gamma: f64,
    normal: &Vector3<f64>,
    left: bool,
) -> Vector3<f64> {
    let rho = q[0];
    let u = q.fixed_rows::<3>(1) / rho;
    let e = q[4] / rho;
    let p = (gamma - 1.0) * rho * (e - 0.5 * u.dot(&u));
    let h0 = e + p / rho;

    //the pipe runs along +normal when met at its left end and along -normal otherwise
    let u_axial = match left {
        true => u.dot(normal),
        false => -u.dot(normal),
    };

    // h0 = gamma*p/((gamma-1)*rho) + 0.5*u^2, solved for p
    let p_g = (gamma - 1.0) / gamma * rho * (h0 - 0.5 * u_axial * u_axial);
    let e_g = p_g / ((gamma - 1.0) * rho) + 0.5 * u_axial * u_axial;

    Vector3::new(rho, rho * u_axial, rho * e_g)
}

// The junction interface uses the same RoeM flux as roem1d.rs, widened to the three
// momentum components a ghost junction cell carries.
// Assembled in the paper's compact form (Eq. 40a) rather than through eigenvectors, which
// is the same algebra without a 5x5 inverse. f_c and g_c are omitted because they suppress
// the carbuncle, which needs a shock spanning several cells across the front, and a
// junction is one control volume with independent two-state interfaces.
// https://doi.org/10.1002/fld.2212

///Density, velocity, pressure and stagnation enthalpy of one conservative state.
fn decode5(q: &Vector5<f64>, gamma: f64) -> (f64, Vector3<f64>, f64, f64) {
    let rho = q[0];
    let u = Vector3::new(q[1], q[2], q[3]) / rho;
    let e = q[4] / rho;
    let p = (gamma - 1.0) * rho * (e - 0.5 * u.dot(&u));

    (rho, u, p, e + p / rho)
}

///Packs density, velocity and specific total energy back into a conservative state.
fn pack5(rho: f64, u: &Vector3<f64>, e: f64) -> Vector5<f64> {
    Vector5::new(rho, rho * u[0], rho * u[1], rho * u[2], rho * e)
}

///Euler flux of one state across a face with the given normal.
fn euler_flux5(q: &Vector5<f64>, p: f64, un: f64, n: &Vector3<f64>) -> Vector5<f64> {
    let momentum = Vector3::new(q[1], q[2], q[3]) * un + n * p;

    Vector5::new(
        q[0] * un,
        momentum[0],
        momentum[1],
        momentum[2],
        (q[4] + p) * un,
    )
}

///Lifts a branch cell into the junction frame, along that pipe's own axis.
fn lift_branch(q: &Vector3<f64>, n: &Vector3<f64>, left: bool) -> Vector5<f64> {
    //the pipe runs along +normal when met at its left end and along -normal otherwise
    let axis = match left {
        true => *n,
        false => -n,
    };

    pack5(q[0], &(axis * (q[1] / q[0])), q[2] / q[0])
}

///Calculates the RoeM flux across one junction interface.
fn roem5(q_l: &Vector5<f64>, q_r: &Vector5<f64>, n: &Vector3<f64>, gamma: f64) -> Vector5<f64> {
    let (rho_l, vel_l, p_l, h_l) = decode5(q_l, gamma);
    let (rho_r, vel_r, p_r, h_r) = decode5(q_r, gamma);

    let un_l = vel_l.dot(n);
    let un_r = vel_r.dot(n);

    //intermediate quanties
    let r = (rho_r / rho_l).sqrt();

    //Roe averages
    let roe_rho = r * rho_l; // Roe average density
    let roe_vel = (vel_r * r + vel_l) / (r + 1.0); // Roe average velocity vector
    let half_roe_u_squared = 0.5 * roe_vel.dot(&roe_vel); //uses the full magnitude, not just U_n
    let roe_h = (r * h_r + h_l) / (r + 1.0); // Roe average specific total enthalpy
    let roe_a = ((gamma - 1.0) * (roe_h - half_roe_u_squared)).sqrt(); // Roe average speed of sound
    let roe_un = roe_vel.dot(n); // Roe average normal velocity

    //RoeM Changes ==================================================
    //Eq 33: the signal velocities take the common speed of sound, which is what
    //lets a contact be captured exactly whichever side is the hotter
    let b1 = (roe_un + roe_a).max((un_r + roe_a).max(0.0));
    let b2 = (roe_un - roe_a).min((un_l - roe_a).min(0.0));
    let b5 = 1.0 / (b1 - b2);

    //other quantities
    let m_hat = roe_un / roe_a; // Roe average Mach number

    // entropy-wave correction B∆Q,
    let hlle_coeff = b1 * b2 * b5;
    let bdq_coeff = hlle_coeff / (1.0 + m_hat.abs()); // full prefactor

    let dp = p_r - p_l; // Δp
    let dh = h_r - h_l; // ΔH
    let d_un = un_r - un_l; // ΔU

    let b_dq_0 = (rho_r - rho_l) - dp / (roe_a * roe_a); // Δρ - Δp/â²

    //the tangential rows vanish identically in 1D, which is why roem1d.rs has only two
    let d_ut = (vel_r - vel_l) - n * d_un;
    let b_dq_mom = roe_vel * b_dq_0 + d_ut * roe_rho;

    let b_dq: Vector5<f64> = Vector5::new(
        b_dq_0,
        b_dq_mom[0],
        b_dq_mom[1],
        b_dq_mom[2],
        b_dq_0 * roe_h + roe_rho * dh,
    );

    let correction: Vector5<f64> = b_dq * bdq_coeff;

    //ΔQ* swaps ρE for ρh, which is where roem1d.rs adds a separate enthalpy_shift
    let mut dq_star: Vector5<f64> = q_r - q_l;
    dq_star[4] = rho_r * h_r - rho_l * h_l;

    // ==============================================================

    let f_l = euler_flux5(q_l, p_l, un_l, n);
    let f_r = euler_flux5(q_r, p_r, un_r, n);

    // RoeM flux at the interface: HLLE base on the physical fluxes, minus the correction
    (f_l * b1 - f_r * b2) * b5 + dq_star * hlle_coeff - correction
}

///Projects a junction interface flux into the attached pipe's 1D coordinates.
/// Mass and energy flip when the pipe runs opposite the outward normal, momentum does not.
fn project_flux(f: &Vector5<f64>, n: &Vector3<f64>, left: bool) -> Vector3<f64> {
    let sign = match left {
        true => 1.0,
        false => -1.0,
    };
    let momentum = Vector3::new(f[1], f[2], f[3]).dot(n);

    Vector3::new(sign * f[0], momentum, sign * f[4])
}

///Interface flux between a junction and one branch, in 3D and projected to that branch.
/// The junction side is reconstructed by the scaling function G before the Riemann solve.
pub fn interface_flux(
    q_j: &Vector5<f64>,
    q_b: &Vector3<f64>,
    n: &Vector3<f64>,
    left: bool,
    gamma: f64,
) -> (Vector5<f64>, Vector3<f64>) {
    let q_r = lift_branch(q_b, n, left);

    //G is measured against the Roe average of the pair BEFORE reconstruction,
    //since reconstruction is what it goes on to change
    let (rho_j, vel_j, p_j, h_j) = decode5(q_j, gamma);
    let (rho_r, vel_r, _, h_r) = decode5(&q_r, gamma);

    let r = (rho_r / rho_j).sqrt();
    let roe_vel = (vel_r * r + vel_j) / (r + 1.0);
    let roe_h = (r * h_r + h_j) / (r + 1.0);
    let roe_a = ((gamma - 1.0) * (roe_h - 0.5 * roe_vel.dot(&roe_vel))).sqrt();

    let un_j = vel_j.dot(n);
    let un_r = vel_r.dot(n);
    let d_un = un_r - un_j;

    //Eq 35: normals point out, so a branch moving along its normal is leaving
    let delta = match un_r > 0.0 {
        true => 1.0,
        false => 0.0,
    };
    let g = delta * (d_un.abs() / roe_a).min(1.0);

    //Eq 37: only the normal component is scaled, density and pressure stay the junction's
    let vel_t = vel_j - n * un_j;
    let vel_l = n * (un_r - g * d_un) + vel_t;
    let e_l = p_j / ((gamma - 1.0) * rho_j) + 0.5 * vel_l.dot(&vel_l);
    let q_l = pack5(rho_j, &vel_l, e_l);

    let f = roem5(&q_l, &q_r, n, gamma);

    (f, project_flux(&f, n, left))
}

///One pipe end attached to one junction, held at runtime so the network can be rewired.
#[derive(Clone, Copy)]
pub struct Link {
    pub junction: usize,
    pub pipe: usize,
    ///Whether the junction sits at the left (index 0) end of the pipe
    pub left: bool,
    ///unit normal pointing out of the junction along the axis of the pipe
    pub normal: Vector3<f64>,
}

//Data clarifying how a pipe interacts with the junction
#[derive(Default)]
pub struct InletData {
    pub pipe_id: usize,
    ///unit normal at this interface, pointing out of the junction along the pipe axis
    pub normal: Vector3<f64>,
    ///surface area of pipe cross section
    pub area: f64,
    ///interface flux at this face, filled by the driver before the junction residual
    pub f: Vector5<f64>,
    ///Whether this junction is on the left (index 0) end of the pipe
    pub left: bool,
}
