// Boundaries follow the partial Riemann problem of Dubois: each end is a manifold whose
// codimension is the number of characteristics entering the domain, and the boundary state
// is where the interior's outgoing wave meets it. The face then takes the plain Euler flux
// of that state, so the prescribed quantity is satisfied exactly rather than relaxed toward.
// https://arxiv.org/abs/1101.2752

use nalgebra::Vector3;

///How one end of a pipe is driven.
#[derive(Clone, Copy)]
#[allow(dead_code)] //variants are selected by the topology in main.rs
pub enum BoundaryCondition {
    ///Connected to a junction with the given id
    Junction(usize),
    ///mass flow into this end in kg/s and the stagnation enthalpy it carries, in J/kg
    MassInflow { mdot: f64, h0: f64 },
    ///mass flow out of this end in kg/s
    MassOutflow { mdot: f64 },
    ///static pressure held at this end in Pa
    Pressure(f64),
    ///open end meeting still air at this pressure in Pa and density in kg/m3
    Atmosphere { p: f64, rho: f64 },
    ///lets waves leave without returning, against the state this end started in
    NonReflecting,
}

///Density, velocity and pressure of one conservative state.
fn decode(q: &Vector3<f64>, gamma: f64) -> (f64, f64, f64) {
    let rho = q[0];
    let u = q[1] / rho;

    (rho, u, (gamma - 1.0) * rho * (q[2] / rho - 0.5 * u * u))
}

///Packs density, velocity and pressure back into a conservative state.
pub(crate) fn pack(rho: f64, u: f64, p: f64, gamma: f64) -> Vector3<f64> {
    Vector3::new(rho, rho * u, rho * (p / ((gamma - 1.0) * rho) + 0.5 * u * u))
}

///Euler flux of one conservative state, which is what a resolved end hands its face.
pub fn euler_flux(q: &Vector3<f64>, gamma: f64) -> Vector3<f64> {
    let (rho, u, p) = decode(q, gamma);

    Vector3::new(rho * u, rho * u * u + p, u * (q[2] + p))
}

///State on a pipe's boundary face, where the interior's outgoing wave meets the manifold
/// the condition describes. `reference` is the state this end started in.
pub fn boundary_state(
    interior: &Vector3<f64>,
    reference: &Vector3<f64>,
    bc: &BoundaryCondition,
    area: f64,
    gamma: f64,
    left: bool,
) -> Vector3<f64> {
    //measuring velocity along the outward normal lets one wave curve serve both ends
    let s = match left {
        true => -1.0,
        false => 1.0,
    };
    let (rho_i, u_i, p_i) = decode(interior, gamma);
    let w = Outward::new(rho_i, s * u_i, p_i, gamma);

    //a supersonic outflow admits no datum, so the interior passes straight through
    if w.v >= w.a {
        return *interior;
    }

    let (rho, v, p) = match *bc {
        BoundaryCondition::Pressure(target) => w.at_pressure(target),
        BoundaryCondition::MassOutflow { mdot } => w.at_mass_flux(mdot / area),
        BoundaryCondition::MassInflow { mdot, h0 } => w.at_jet(-mdot / area, h0),
        BoundaryCondition::Atmosphere { p, rho } => w.at_atmosphere(p, rho),
        BoundaryCondition::NonReflecting => {
            let (rho_r, u_r, p_r) = decode(reference, gamma);
            w.at_invariant(rho_r, s * u_r, p_r)
        }
        BoundaryCondition::Junction(_) => unreachable!("a junction end is resolved by the driver"),
    };

    pack(rho, s * v, p, gamma)
}

///The interior cell as the boundary sees it, with velocity along the outward normal so
/// that positive always means leaving the pipe.
struct Outward {
    rho: f64,
    v: f64,
    p: f64,
    a: f64,
    gamma: f64,
}

impl Outward {
    fn new(rho: f64, v: f64, p: f64, gamma: f64) -> Self {
        Self {
            rho,
            v,
            p,
            a: (gamma * p / rho).sqrt(),
            gamma,
        }
    }

    ///Velocity drop across the wave that reaches the boundary, from Eq 1.3.12 and 1.4.16.
    /// Rarefaction below the interior pressure, shock above, meeting smoothly between.
    fn f(&self, p: f64) -> f64 {
        let g = self.gamma;
        match p <= self.p {
            true => 2.0 * self.a / (g - 1.0) * ((p / self.p).powf((g - 1.0) / (2.0 * g)) - 1.0),
            false => {
                (p - self.p) * (2.0 / (self.rho * ((g + 1.0) * p + (g - 1.0) * self.p))).sqrt()
            }
        }
    }

    ///Slope of that curve, which both branches give as 1/(rho*a) at the interior pressure.
    fn df(&self, p: f64) -> f64 {
        let g = self.gamma;
        match p <= self.p {
            true => self.a / (g * self.p) * (p / self.p).powf(-(g + 1.0) / (2.0 * g)),
            false => {
                let d = (g + 1.0) * p + (g - 1.0) * self.p;
                (2.0 / (self.rho * d)).sqrt() * (1.0 - (g + 1.0) * (p - self.p) / (2.0 * d))
            }
        }
    }

    ///Outward velocity the boundary takes if it settles at the given pressure.
    fn v_at(&self, p: f64) -> f64 {
        self.v - self.f(p)
    }

    ///Density the wave leaves behind at the given pressure, from Eq 1.6.12 and 1.6.13:
    /// isentropic through a rarefaction, Hugoniot through a shock.
    fn rho_at(&self, p: f64) -> f64 {
        let g = self.gamma;
        let mu2 = (g - 1.0) / (g + 1.0);
        match p <= self.p {
            true => self.rho * (p / self.p).powf(1.0 / g),
            false => self.rho * (p + mu2 * self.p) / (self.p + mu2 * p),
        }
    }

    ///Slope of that density curve.
    fn drho(&self, p: f64) -> f64 {
        let g = self.gamma;
        let mu2 = (g - 1.0) / (g + 1.0);
        match p <= self.p {
            true => self.rho_at(p) / (g * p),
            false => {
                let d = self.p + mu2 * p;
                self.rho * self.p * (1.0 - mu2 * mu2) / (d * d)
            }
        }
    }

    ///Boundary state at a prescribed static pressure, which needs no iteration.
    fn at_pressure(&self, p: f64) -> (f64, f64, f64) {
        (self.rho_at(p), self.v_at(p), p)
    }

    ///State where the interior's rarefaction reaches Mach 1, which is as much as this end
    /// can ever discharge.
    fn sonic(&self) -> (f64, f64, f64) {
        let g = self.gamma;
        let a = ((g - 1.0) * self.v + 2.0 * self.a) / (g + 1.0);
        let ratio = a / self.a;

        (
            self.rho * ratio.powf(2.0 / (g - 1.0)),
            a,
            self.p * ratio.powf(2.0 * g / (g - 1.0)),
        )
    }

    ///Boundary state carrying a prescribed outward mass flux, saturating at the sonic
    /// state once more is asked for than this end can pass.
    fn at_mass_flux(&self, flux: f64) -> (f64, f64, f64) {
        let (rho_s, v_s, p_s) = self.sonic();
        if flux >= rho_s * v_s {
            return (rho_s, v_s, p_s);
        }

        //mass flux falls with pressure on the subsonic branch, so doubling always brackets
        let mut hi = self.p.max(p_s) * 2.0;
        while self.rho_at(hi) * self.v_at(hi) > flux {
            hi *= 2.0;
        }

        let p = root(p_s, hi, |p| {
            (
                self.rho_at(p) * self.v_at(p) - flux,
                self.drho(p) * self.v_at(p) - self.rho_at(p) * self.df(p),
            )
        });

        (self.rho_at(p), self.v_at(p), p)
    }

    ///Boundary state drawing a prescribed inward mass flux that carries the given
    /// stagnation enthalpy, which is the jet manifold of Eq 3.4.2.
    fn at_jet(&self, flux: f64, h0: f64) -> (f64, f64, f64) {
        let g = self.gamma;

        //nothing coming in is a closed end, which the zero mass flux solve already gives
        if flux == 0.0 {
            return self.at_mass_flux(0.0);
        }

        let scale = g / ((g - 1.0) * flux.abs());

        //eliminating density leaves a quadratic in velocity; the inward root is the one
        let inward = |p: f64| {
            let b = scale * p;
            let s = (b * b + 2.0 * h0).sqrt();
            (b - s, (1.0 - b / s) * scale)
        };

        let mut hi = self.p * 2.0;
        while self.v_at(hi) > inward(hi).0 {
            hi *= 2.0;
        }

        let p = root(0.0, hi, |p| {
            let (v_m, dv_m) = inward(p);
            (self.v_at(p) - v_m, -self.df(p) - dv_m)
        });

        let v = self.v_at(p);
        let rho = flux / v;
        assert!(
            v > -(g * p / rho).sqrt(),
            "supersonic inflow needs a third condition, which this end does not carry"
        );

        (rho, v, p)
    }

    ///Boundary state at an open end meeting still air, discharging against ambient
    /// pressure or drawing in through a nozzle, and choking either way when it must.
    fn at_atmosphere(&self, p_a: f64, rho_a: f64) -> (f64, f64, f64) {
        let g = self.gamma;
        let h0 = g * p_a / ((g - 1.0) * rho_a);

        //an open end sees ambient pressure until its own exit plane goes sonic
        let (rho, v, p) = self.at_pressure(p_a);
        if v >= 0.0 {
            return match v < (g * p / rho).sqrt() {
                true => (rho, v, p),
                false => self.sonic(),
            };
        }

        //drawing in instead, so the gas arrives from still air through what is in effect a
        //convergent nozzle, holding stagnation enthalpy and entropy (Eq 3.4.10)
        let critical = (2.0 / (g + 1.0)).powf(g / (g - 1.0));
        let p_crit = p_a * critical;
        let a_star = (2.0 * (g - 1.0) * h0 / (g + 1.0)).sqrt();
        if self.v_at(p_crit) <= -a_star {
            return (
                rho_a * (2.0 / (g + 1.0)).powf(1.0 / (g - 1.0)),
                -a_star,
                p_crit,
            );
        }

        let p = root(p_crit, p_a, |p| {
            let z = (p / p_a).powf((g - 1.0) / g);
            let s = (2.0 * h0 * (1.0 - z)).sqrt();
            (
                self.v_at(p) + s,
                -self.df(p) - h0 * (g - 1.0) / g * z / (p * s),
            )
        });

        (rho_a * (p / p_a).powf(1.0 / g), self.v_at(p), p)
    }

    ///Boundary state holding the incoming Riemann invariant at the value this end started
    /// with, so a wave reaching it leaves without sending one back.
    fn at_invariant(&self, rho_r: f64, v_r: f64, p_r: f64) -> (f64, f64, f64) {
        let g = self.gamma;
        let a_r = (g * p_r / rho_r).sqrt();

        //the wave curve already carries the outgoing invariant, so the pair is linear
        let incoming = v_r - 2.0 * a_r / (g - 1.0);
        let outgoing = self.v + 2.0 * self.a / (g - 1.0);
        let v = 0.5 * (outgoing + incoming);
        let a = 0.25 * (g - 1.0) * (outgoing - incoming);
        assert!(a > 0.0, "non-reflecting end would cavitate");

        //entropy rides the flow, so it comes from whichever side is upstream
        let (rho_u, a_u, p_u) = match v >= 0.0 {
            true => (self.rho, self.a, self.p),
            false => (rho_r, a_r, p_r),
        };
        let ratio = a / a_u;

        (
            rho_u * ratio.powf(2.0 / (g - 1.0)),
            v,
            p_u * ratio.powf(2.0 * g / (g - 1.0)),
        )
    }
}

///Where a residual that falls with pressure crosses zero, taking Newton steps but
/// bisecting whenever one would leave the bracket.
fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> (f64, f64)) -> f64 {
    let mut p = 0.5 * (lo + hi);

    for _ in 0..60 {
        let (r, dr) = residual(p);
        match r > 0.0 {
            true => lo = p,
            false => hi = p,
        }

        //a NaN step fails both comparisons, which drops it back to bisection
        let step = p - r / dr;
        let next = match step > lo && step < hi {
            true => step,
            false => 0.5 * (lo + hi),
        };

        if (next - p).abs() <= 1e-14 * p.abs() {
            return next;
        }
        p = next;
    }

    p
}
