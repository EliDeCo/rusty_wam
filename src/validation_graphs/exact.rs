// The analytic solution of a Riemann problem, which is what every benchmark here is graded
// against. Wave curves and the sampler follow Toro, Riemann Solvers and Numerical Methods
// for Fluid Dynamics, chapter 4, which is reference [28] of the Sod survey.

use crate::boundaries::root;

///Density, velocity and pressure on one side of the initial discontinuity.
#[derive(Clone, Copy)]
pub struct State {
    pub rho: f64,
    pub u: f64,
    pub p: f64,
}

impl State {
    pub const fn new(rho: f64, u: f64, p: f64) -> Self {
        Self { rho, u, p }
    }

    ///Speed of sound in this state.
    fn a(&self, gamma: f64) -> f64 {
        (gamma * self.p / self.rho).sqrt()
    }
}

///A resolved Riemann problem, which can then answer any position and time.
/// Resolving means finding the pressure and velocity of the contact the two waves leave behind.
pub struct Riemann {
    left: State,
    right: State,
    a_l: f64,
    a_r: f64,
    p_star: f64,
    u_star: f64,
    x0: f64,
    gamma: f64,
}

impl Riemann {
    ///Solves for the star region by matching the two wave curves at a common pressure.
    pub fn new(left: State, right: State, x0: f64, gamma: f64) -> Self {
        let (a_l, a_r) = (left.a(gamma), right.a(gamma));
        let du = right.u - left.u;

        //both curves rise with pressure, so their sum crosses zero once and doubling brackets it
        let total = |p: f64| {
            let (f_l, df_l) = wave(&left, a_l, p, gamma);
            let (f_r, df_r) = wave(&right, a_r, p, gamma);
            (f_l + f_r + du, df_l + df_r)
        };

        let mut hi = left.p.max(right.p) * 2.0;
        while total(hi).0 < 0.0 {
            hi *= 2.0;
        }

        //root expects a residual that falls with pressure, so the rising sum is negated
        let p_star = root(0.0, hi, |p| {
            let (f, df) = total(p);
            (-f, -df)
        });

        let u_star = 0.5 * (left.u + right.u)
            + 0.5 * (wave(&right, a_r, p_star, gamma).0 - wave(&left, a_l, p_star, gamma).0);

        Self {
            left,
            right,
            a_l,
            a_r,
            p_star,
            u_star,
            x0,
            gamma,
        }
    }

    ///The exact state at one position and time, which the whole solution is self-similar in.
    pub fn at(&self, x: f64, t: f64) -> (f64, f64, f64) {
        //before the diaphragm breaks there is no fan to sample, only the two initial states
        if t <= 0.0 {
            let s = if x < self.x0 { self.left } else { self.right };
            return (s.rho, s.u, s.p);
        }

        let xi = (x - self.x0) / t;
        let g = self.gamma;

        //the contact carries the flow, so which side of it the sample sits on picks the wave
        if xi < self.u_star {
            self.side(&self.left, self.a_l, xi, -1.0, g)
        } else {
            self.side(&self.right, self.a_r, xi, 1.0, g)
        }
    }

    ///Samples one side of the contact. `s` is +1 for the right wave and -1 for the left, which
    /// is the only thing that differs between the two mirrored branches.
    fn side(&self, k: &State, a_k: f64, xi: f64, s: f64, g: f64) -> (f64, f64, f64) {
        let ratio = self.p_star / k.p;

        //outside the wave entirely, the initial state still stands
        let past = |speed: f64| s * (xi - speed) > 0.0;

        if ratio > 1.0 {
            //shock, whose speed comes from the Rankine-Hugoniot jump
            let speed =
                k.u + s * a_k * ((g + 1.0) / (2.0 * g) * ratio + (g - 1.0) / (2.0 * g)).sqrt();
            if past(speed) {
                return (k.rho, k.u, k.p);
            }
            let mu2 = (g - 1.0) / (g + 1.0);
            return (
                k.rho * (ratio + mu2) / (mu2 * ratio + 1.0),
                self.u_star,
                self.p_star,
            );
        }

        //rarefaction, bounded by its head at the initial sound speed and its tail at the star one
        let a_star = a_k * ratio.powf((g - 1.0) / (2.0 * g));
        if past(k.u + s * a_k) {
            return (k.rho, k.u, k.p);
        }
        if !past(self.u_star + s * a_star) {
            return (k.rho * ratio.powf(1.0 / g), self.u_star, self.p_star);
        }

        //inside the fan, where the characteristic through this point sets the state
        let w = 2.0 / (g + 1.0) - s * (g - 1.0) / ((g + 1.0) * a_k) * (k.u - xi);
        (
            k.rho * w.powf(2.0 / (g - 1.0)),
            2.0 / (g + 1.0) * (-s * a_k + (g - 1.0) / 2.0 * k.u + xi),
            k.p * w.powf(2.0 * g / (g - 1.0)),
        )
    }

    ///Cell averages of the exact solution on a uniform grid, which is what a finite volume
    /// state actually holds. Sampling cell centres instead would misplace every jump by up to
    /// half a cell, an error of the same order as the scheme's own.
    pub fn cell_averages(&self, n_cells: usize, dx: f64, t: f64) -> [Vec<f64>; 3] {
        const SUB: usize = 64;

        let mut out = [vec![0.0; n_cells], vec![0.0; n_cells], vec![0.0; n_cells]];

        for j in 0..n_cells {
            let (mut rho, mut u, mut p) = (0.0, 0.0, 0.0);
            for k in 0..SUB {
                let x = (j as f64 + (k as f64 + 0.5) / SUB as f64) * dx;
                let (r, v, q) = self.at(x, t);
                rho += r;
                u += v;
                p += q;
            }
            out[0][j] = rho / SUB as f64;
            out[1][j] = u / SUB as f64;
            out[2][j] = p / SUB as f64;
        }

        out
    }

    pub fn p_star(&self) -> f64 {
        self.p_star
    }
    pub fn u_star(&self) -> f64 {
        self.u_star
    }

    ///Speeds of every wave the problem carries, left to right: the leading and trailing edge
    /// of the left wave, the contact, then the trailing and leading edge of the right one.
    /// A shock has no width, so its two edges coincide.
    pub fn wave_speeds(&self) -> [f64; 5] {
        let (l_head, l_tail) = self.edges(&self.left, self.a_l, -1.0);
        let (r_head, r_tail) = self.edges(&self.right, self.a_r, 1.0);

        [l_head, l_tail, self.u_star, r_tail, r_head]
    }

    ///Points the analytic curve should be drawn from: dense through a rarefaction fan, and
    /// only the two ends of every constant region. A jump then stays a jump on the page
    /// instead of becoming a slope, at a fraction of the geometry a uniform sampling costs.
    pub fn curve(&self, length: f64, t: f64, fan: usize) -> [Vec<[f64; 2]>; 3] {
        let mut bounds = vec![0.0];
        bounds.extend(
            self.wave_speeds()
                .iter()
                .map(|s| (self.x0 + s * t).clamp(0.0, length)),
        );
        bounds.push(length);

        let mut out = [Vec::new(), Vec::new(), Vec::new()];

        for (i, pair) in bounds.windows(2).enumerate() {
            let (lo, hi) = (pair[0], pair[1]);
            if hi <= lo {
                continue;
            }

            //only the two waves can spread into a fan; everything else is a constant state
            let n = if i == 1 || i == 4 { fan } else { 1 };

            //nudging off each end keeps the sample on this region's side of the jump
            let edge = (hi - lo) * 1e-9;
            for k in 0..=n {
                let x = lo + edge + (hi - lo - 2.0 * edge) * k as f64 / n as f64;
                let (rho, u, p) = self.at(x, t);
                out[0].push([x, rho]);
                out[1].push([x, u]);
                out[2].push([x, p]);
            }
        }

        out
    }

    ///Leading and trailing edge of the wave on one side of the contact.
    fn edges(&self, k: &State, a_k: f64, s: f64) -> (f64, f64) {
        let g = self.gamma;
        let ratio = self.p_star / k.p;

        if ratio > 1.0 {
            let speed =
                k.u + s * a_k * ((g + 1.0) / (2.0 * g) * ratio + (g - 1.0) / (2.0 * g)).sqrt();
            return (speed, speed);
        }

        let a_star = a_k * ratio.powf((g - 1.0) / (2.0 * g));
        (k.u + s * a_k, self.u_star + s * a_star)
    }
}

///Velocity change across the wave that takes one side to pressure `p`, and its slope.
/// Shock above that side's pressure, rarefaction below, meeting smoothly at it.
fn wave(k: &State, a_k: f64, p: f64, gamma: f64) -> (f64, f64) {
    let g = gamma;

    if p > k.p {
        let a = 2.0 / ((g + 1.0) * k.rho);
        let b = (g - 1.0) / (g + 1.0) * k.p;
        let q = (a / (b + p)).sqrt();

        ((p - k.p) * q, q * (1.0 - 0.5 * (p - k.p) / (b + p)))
    } else {
        let e = (g - 1.0) / (2.0 * g);

        (
            2.0 * a_k / (g - 1.0) * ((p / k.p).powf(e) - 1.0),
            1.0 / (k.rho * a_k) * (p / k.p).powf(-(g + 1.0) / (2.0 * g)),
        )
    }
}
