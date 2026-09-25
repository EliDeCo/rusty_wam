# RoeM1D validation

What the interior solver reproduces from Kim, Kim, Rho and Hong 2003
(DOI 10.1016/S0021-9991(02)00037-2), and how closely.

Air at 101325 Pa and 1.225 kg/m3, gamma 1.4, first-order RoeM1D with forward Euler at a
Courant number of 0.5. The exact solutions it is graded against need no Riemann solver: an
isentropic simple wave, which the method of characteristics gives exactly until the
characteristics cross, and shocks, which Rankine-Hugoniot gives algebraically.

## Against the paper's exact flux identities

Taken off the shipped flux directly, by building a piecewise constant pipe and reading the
face flux at the jump.

| Source | Statement | Worst error |
|---|---|---|
| Eq 30 | a stationary contact leaves the flux at `(0, p, 0)` | 2.0e-15 |
| Eq 31 | a moving contact with `T_l > T_r` dissipates at the contact speed | 7.4e-16 |
| Eq 32, 35 | the same with `T_l < T_r`, which Eq 27c cannot do | 1.2e-16 |
| Eq 23 | equal total enthalpy gives `F_E = H F_rho` | 2.0e-16 |
| Section 4.4 | a supersonic jump takes the upwind flux, no intermediate cell | 1.8e-16 |

Eq 34a was also transcribed separately in the paper's own compact form, using no
eigenvectors at all. It agrees with the shipped flux to 6.9e-16 over five jumps, so the Roe
matrix rearrangement, the swap to `dQ*` and the one-dimensional reduction of `B dQ` are
confirmed together.

## Against exact solutions

Isentropic simple wave, run to 60 and 76 per cent of the time its characteristics cross.

| Amplitude | L1 density at 800 cells | Observed order | Entropy drift |
|---|---|---|---|
| Mach 0.088 | 7.3e-4 | 0.86 | 5.5e-5 |
| Mach 0.294 | 1.4e-3 | 0.87 | 6.8e-4 |

Steadily propagating shock at Mach 1.5, 3 and 6, in three frames, on 400 cells. Noise is
the worst departure from the exact plateau more than ten cells clear of the front, as a
fraction of the shock's own pressure jump.

| Frame | Noise behind | Noise ahead | Front position |
|---|---|---|---|
| Stationary | 1.0e-13 | 1.4e-16 | within 0.06 cell |
| Creeping at 0.1a | 6.1e-3 | 9.8e-17 | within 0.05 cell |
| Running at the shock speed | 2.2e-2 | 9.5e-6 | within 1.0 cell |

A stationary shock is captured with no noise at all, to machine precision at every strength
tested.

Shock reflecting off a wall, exact by Rankine-Hugoniot twice, at reflected pressure ratios
to 51.7.

| Quantity | Worst error |
|---|---|
| Reflected plateau pressure | 3.3e-4 |
| Noise behind the reflected front | 2.6e-3 |

## Note: why f and g are omitted

Eq 17b and Eq 20b scale the pressure term in the numerical mass flux and the anti-diffusion
that carries it. Both are written as `abs(M_hat)` raised to `1 - min(p_j/p_j+1, p_j+1/p_j)`,
so on smooth data resolved across several cells the adjacent pressure ratio approaches one,
the exponent approaches zero, and both approach one. They can only act where a pressure
jump spans a single cell, which is to say at a captured discontinuity.

Switched on, over every case above, they change nothing qualitatively:

| Measure | Without | With f | With f and g |
|---|---|---|---|
| L1 density, smooth wave, 100 cells | 7.56e-3 | 7.25e-3 | 7.52e-3 |
| L1 density, smooth wave, 800 cells | 1.41e-3 | 1.39e-3 | 1.40e-3 |
| Worst shock noise, all frames | 2.16e-2 | 2.09e-2 | 2.08e-2 |
| Best shock noise improvement | - | 21 per cent | 13 per cent |
| Reflected plateau pressure | 3.3e-4 | 3.1e-4 | 2.9e-4 |
| Run time | - | +0.0 per cent | +0.9 per cent |

The smooth-wave gap halves between 100 and 800 cells, as the exponent argument predicts.
Stationary shocks are bit identical with and without. Everywhere else `f` is worth between
3 and 21 per cent off a noise figure already at or below 2 per cent of a shock's pressure
jump, which is well inside the roughly one cell of smearing first-order capturing gives
anyway. Adding `g` on top of `f` is worse than `f` alone in six of the eight cases where the
two differ.

They are therefore omitted. The paper's own reason supports it: Section 2.4 states that the
shock instability the two functions exist to cure does not occur in one dimension.

## The same flux in three places

`roem1d.rs`, `muscl_roem1d.rs` and `junctions.rs` all carry this flux: the second applied
to reconstructed face states rather than cell averages, the third widened to the three
momentum components a ghost junction cell holds. All three take the signal velocities of
Eq 33. Handed the two states directly, so that reconstruction and the ghost construction
are out of the way, each satisfies all five identities.

| Identity | roem1d.rs | muscl_roem1d.rs | junctions.rs |
|---|---|---|---|
| Eq 30 | 2.0e-15 | 2.0e-15 | 2.2e-15 |
| Eq 31 | 7.4e-16 | 7.4e-16 | 7.4e-16 |
| Eq 32, 35 | 1.2e-16 | 1.2e-16 | 1.2e-16 |
| Eq 23 | 2.0e-16 | 2.0e-16 | 4.0e-16 |
| Section 4.4 | 1.8e-16 | 1.8e-16 | 1.8e-16 |

The first two columns agree to the bit, which is what sharing the algebra means; the third
differs in the last places only because it carries five components rather than three.

## Not covered

- A Riemann problem with a shock, a contact and a rarefaction interacting at once,
  including a sonic point, which is the paper's Figure 8. Grading it needs an iterative
  exact Riemann solver rather than the closed forms used here.
- Sections 2, 5.2 and 5.3 in full: the odd-even decoupling, the kinked Mach stem, the
  carbuncle, the supersonic corner and every viscous case are two-dimensional or need the
  Navier-Stokes terms and a turbulence model.
- Table 1's 1.2x cost ratio, which is a timing on the authors' own code and mesh.
