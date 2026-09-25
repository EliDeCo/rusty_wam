# Boundary validation

What the solver reproduces from Dubois, *Partial Riemann problem, boundary conditions, and
gas dynamics* (arXiv:1101.2752), and how closely.

Air at 101325 Pa and 1.225 kg/m3. Every condition is resolved against 27 interior states,
spanning -250 to +250 m/s over three thermodynamic states, at both ends of a pipe; the wall
is swept further, to Mach 11 outward and to within 0.01 per cent of the cavitation limit
inward. The wave curves the resolved states are measured against are the paper's own psi
(Eq 1.3.19) and phi (Eq 1.4.20), written from the paper rather than taken from the solver.
Network results use 120 to 600 cells, MUSCL-RoeM with third-order SSP Runge-Kutta.

## The construction, common to every condition

| Source | Statement | Worst error |
|---|---|---|
| Theorem 2, Eq 1.6.6 | the resolved state lies on the interior's wave curve | 7.2e-15 |
| Eq 1.6.12, 1.6.13 | its density is on the isentrope or the Hugoniot | 3.4e-16 |
| Eq 3.3.4 | a state already on the manifold comes back unchanged | 8.4e-17 |
| Section 3.2 (iv), 4.3 | a supersonic outflow admits no datum | bit-identical |

The density row covers the codimension-1 conditions. The unchanged-state row covers the
closed-form conditions; the ones needing a root find return it to 1.9e-12.

## Each condition against its manifold

| Condition | Manifold | Statement | Worst error |
|---|---|---|---|
| `Pressure` | Eq 3.4.11 | static pressure held | 2.2e-16 |
| `Pressure` | Fig 3.7 | "an inflow is absolutely compatible" | confirmed |
| `MassOutflow` | Section 3.2 (iii) | mass flux delivered | 1.6e-13 |
| `MassOutflow` | Prop 2, Eq 2.5.3 | saturates on the sonic point, never past it | 4.4e-16 |
| `MassInflow` | Section 3.2 (ii) | mass flux delivered | 2.2e-16 |
| `MassInflow` | Section 3.2 (ii) | stagnation enthalpy delivered | 6.3e-13 |
| `MassInflow` | Theorem 2 | one crossing of wave curve and manifold | confirmed |
| `Atmosphere` | Eq 3.4.11 | ambient pressure held while discharging | 2.2e-16 |
| `Atmosphere` | Eq 3.4.10 | on the nozzle manifold while drawing in | 3.3e-15 |
| `Atmosphere` | Prop 2, Eq 2.5.3 | choked ends are exactly sonic | 3.3e-16 |
| `Atmosphere` | Eq 1.3.7 | choked discharge is isentropic | 6.7e-16 |
| `Atmosphere` | Eq 3.4.10 | choked inlet at `(2/(g+1))^(g/(g-1))` | exact |
| `NonReflecting` | Eq 3.2.4, `Sigma = 0` | incoming invariant pinned to the reference | 1.6e-15 |
| `NonReflecting` | Eq 1.3.8, 1.3.15 | outgoing invariant carried from the interior | 2.2e-15 |
| `NonReflecting` | Eq 1.2.12 | entropy taken from whichever side is upstream | 2.2e-16 |
| `Wall` | Prop 3, Eq 3.5.4 | stationary velocity, so zero mass and energy flux | exact |
| `Wall` | Eq 3.5.11 | the pressure it settles at | 1.1e-14 |
| `Wall` | Eq 3.5.8 | agrees with the mirror Riemann problem | 9.9e-15 |
| `Wall` | Eq 3.5.12 | `p*(V) = p*(0) - rho c V + O(V^2)` | 1.8e-6 at V = 1 mm/s |
| `Wall` | Eq 1.6.10 | closed form approaching cavitation | 7.8e-15 |

`Wall` is exact rather than merely accurate in Prop 3: the resolved momentum is zero to the
bit, so the face flux is exactly `(0, p*, 0)`.

Two of these are independent rather than self-consistent. Eq 3.5.8 solves a full two-wave
Riemann problem against the mirror state `(rho, -u, p)` with its own Newton, seeded from
Eq 1.6.9, and *derives* the zero velocity that the wall imposes; it recovers it to 4.1e-16.
Eq 3.5.12 is confirmed both by value and by order, the relative error halving with V to
within 0.09 per cent over four bisections.

## Choking

The open end chokes on runtime conditions rather than on a prescribed ratio. Starting from
rest, the thresholds found by bisection match closed forms derived from the same wave
curves.

| Direction | Threshold | Closed form | Error |
|---|---|---|---|
| Discharging | `p_i/p_a` = 3.5832 | `((g+1)/2)^(2g/(g-1))` | 1.2e-12 |
| Drawing in | `p_a/p_i` = 8.4125 | shock branch reaching sonic at the throat | 6.5e-12 |

Neither is the steady critical ratio 1.8929. That ratio is where a steady nozzle chokes
against a fixed back pressure; here the back pressure is set by an unsteady wave in the
pipe, so more pressure difference is needed before the face goes sonic. Once choked the
*state* is the steady sonic one: `p_a/p_b` is 1.8929 and the density ratio
`(2/(g+1))^(1/(g-1))` to the bit. Across all three branch switches the resolved state moves
by less than 1.4e-6 of itself over a 1e-6 change in pressure ratio.

## Reflection

Measured as the returned pressure extremum over the incident, for a 1 per cent Gaussian
pulse on 600 cells after a 1.15 m round trip.

| Condition | Dubois | Expected | Measured |
|---|---|---|---|
| `Wall` | Eq 3.2.12, `Sigma = (0, 1)` | +1 | +1.0001 |
| `Pressure` | Eq 3.2.9, `Sigma = (0, -1)` | -1 | -0.9915 |
| `Atmosphere` | Eq 3.4.11, pressure release | -1 | -0.9915 |
| `NonReflecting` | `Sigma = 0` | 0 | -0.0022 |

These carry the interior scheme's dissipation over the round trip as well as the boundary's
error, so the one per cent departures are an upper bound on the boundary alone. They are an
upper bound in a measurable sense: the pulse is a Gaussian, so it carries a smooth extremum,
and when MUSCL-RoeM stopped clipping that extremum the wall figure moved from +0.9995 to
+1.0001 and the pressure figure from -0.9907 to -0.9915 without the boundaries changing at
all. What the table bounds is the pair, not the boundary.

## Network level

| Statement | Worst error |
|---|---|
| Eq 4.2.3: a walled pipe conserves mass and total energy, both solvers | 3.5e-14 |
| Section 3.3: a uniform flow its two ends agree with stays put | 6.8e-14 |

## Note: two conditions are not among the paper's worked examples

Dubois works six manifolds: given state (Eq 3.4.1), jet (3.4.2), nozzle (3.4.8), static
pressure (3.4.11), supersonic outflow (3.4.12) and wall or moving boundary (3.5.3, 3.5.9).
`Pressure`, `Atmosphere` and `Wall` are his; `NonReflecting` is his Eq 3.2.4 with the
reflection operator set to zero.

`MassOutflow` and `MassInflow` are not. They prescribe mass flux alone, and mass flux with
stagnation enthalpy, where his jet prescribes mass flux with static temperature. Both sit
inside his framework at the codimension Section 3.2 requires - one datum for a subsonic
outflow, two for a subsonic inflow, where his list of admissible pairs is explicitly open -
and both are validated as manifolds rather than against a worked example: the resolved
state lies on the interior's wave curve, the prescribed quantities are delivered exactly,
and the crossing is unique.

## Note: codimension-2 ends do not keep the interior's entropy

Eq 2.4.2 gives one simple wave per unit of codimension. A codimension-1 end is joined to
the interior by a single wave, so its density follows Eq 1.6.12 or 1.6.13. A codimension-2
end - `MassInflow`, and `NonReflecting` while gas is entering - is joined by a wave *and* a
contact discontinuity, so its density is set by the incoming side instead. Pressure and
velocity still lie on the wave curve, to 7.2e-15. This is the construction rather than a
departure from it.

## Note: the floor approaching cavitation

The root find bisects at most sixty times over a bracket scaled to the interior pressure,
so relative accuracy in the resolved pressure falls away as that pressure collapses toward
vacuum: 1.3e-15 at `p*/p_i` of 7.8e-3, 9.3e-10 at 1e-14, and below roughly 2^-60 of the
interior pressure it stops resolving and clamps. The clamped value stays positive and
finite. At the reference state that floor is about 1e-16 Pa, far below anything the solver
would still accept as physical.

## Not covered

- The moving boundary of Eq 3.5.9. The velocity manifold takes a target, so the machinery
  is present, but no boundary condition exposes a non-zero one.
- The given-state manifold of Eq 3.4.1, equivalently a supersonic inflow, which prescribes
  the entire state.
- Dubois's own second-order boundary treatment is followed rather than measured: Eq 4.4.18
  sets the boundary face to the cell value instead of a reconstructed one, which is what
  the injected flux does.
