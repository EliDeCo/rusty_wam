use crate::helpers::unphysical;
use crate::pipes::InteriorMethod;
use nalgebra::{Matrix3xX, Vector3, Vector5};
use std::collections::BTreeMap;

pub struct Junction {
    //state is now 1 vector (only 1 cell) of 5 values, which is the same
    //as the original vector but with momentum into x, y, z directions
    pub q0: Vector5<f64>, //previous state
    pub q1: Vector5<f64>, //current state

    ///volume of the junction in m^3
    pub volume: f64,

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
    ///decodes the current state vector (q1) into primitives of the conserved variables
    pub fn decode(&mut self) {
        self.rho = self.q1[0];

        let without_rho = self.q1 / self.rho;

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
    ///Calculates the Euler fluxes (F) at every pipe's face normal from the decoded
    /// primitives. Normals point out of the junction, so these are the outward fluxes
    /// the control-volume divergence sum expects and u_dot is positive for outflow.
    pub fn euler_fluxes(&mut self) {
        //f = [rho * (u . n), rho * u * (u . n) + p * n, (rho * E + p) * (u . n)]

        for pipe in self.pipes.iter_mut() {
            let n = pipe.normal;
            pipe.u_dot = self.u.dot(&n);

            pipe.f[0] = self.rho * pipe.u_dot; //mass flux

            //momentum flux
            pipe.f
                .fixed_rows_mut::<3>(1)
                .copy_from(&(self.q1.fixed_rows::<3>(1) * pipe.u_dot + self.p * n));

            pipe.f[4] = (self.q1[4] + self.p) * pipe.u_dot; //energy flux
        }
    }

    ///Checks for a density or pressure that is not physically usable, which indicates
    /// a numerical blowup.
    pub fn unreal_check(&self) {
        if unphysical(self.rho) || unphysical(self.p) {
            panic!("Unphysical state in junction {}", self.id);
        }
    }

    ///Moves current solution to previous solution buffer
    pub fn save_step(&mut self) {
        self.q0 = self.q1;
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

    //Finite volume update, leaves primitives stale
    pub fn advance(&mut self, dt: f64) {
        if self.volume == 0.0 {
            panic!(
                "Junction {}: Volume has not been initialized or no pipes are attached.",
                self.id
            );
        }

        //TODO: Impliment the finite volume update step
    }

    pub fn new(gamma: f64, courant: f64, id: usize) -> Self {
        Self {
            q0: Vector5::zeros(),
            q1: Vector5::zeros(),
            volume: 0.0,
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
        let n_pipes = self.pipes.len();
        let mut volumes: Vec<f64> = Vec::with_capacity(n_pipes);
        let mut states: Matrix3xX<f64> = Matrix3xX::zeros(n_pipes);
        let mut velocities: Vec<Vector3<f64>> = Vec::with_capacity(n_pipes);
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
            areas.push(inlet.area);
        }

        let total_vol = volumes.iter().sum::<f64>();

        //set volume to mean of adjacent cell volumes
        self.volume = total_vol / volumes.len() as f64;

        //set mass density and energy density by weighted average
        self.q0[0] = volumes
            .iter()
            .zip(states.row(0))
            .fold(0.0, |sum, (weight, val)| sum + val * weight)
            / total_vol;
        self.q0[4] = volumes
            .iter()
            .zip(states.row(2))
            .fold(0.0, |sum, (weight, val)| sum + val * weight)
            / total_vol;

        //set velocity as a weighted average of vectors by pipe inlet surface area
        let total_surface_area = areas.iter().sum::<f64>();
        let u_mean = velocities
            .iter()
            .zip(areas)
            .fold(Vector3::zeros(), |sum, (val, weight)| sum + val * weight)
            / total_surface_area;

        //rows 1..4 hold MOMENTUM density, so that averaged velocity is scaled by the
        //density set above before being stored
        let rho_mean = self.q0[0];
        self.q0
            .fixed_rows_mut::<3>(1)
            .copy_from(&(u_mean * rho_mean));

        //q1 is what decode() and every flux loop read, so it has to start from q0
        self.q1 = self.q0;
    }
}

//Data clarifying how a pipe interacts with the junction
#[derive(Default)]
pub struct InletData {
    pub pipe_id: usize,
    ///unit normal at this interface, pointing out of the junction along the pipe axis
    pub normal: Vector3<f64>,
    ///surface area of pipe cross section
    pub area: f64,
    ///the component of velocity in the direction of the pipe's normal,
    /// so positive means flow leaving the junction
    pub u_dot: f64,
    ///physical (euler) flux at the pipe interface
    pub f: Vector5<f64>,
    ///Whether this junction is on the left (index 0) end of the pipe
    pub left: bool,
}
