# MusclRoeM1D validation

What the reconstruction reproduces from van Leer and Nishikawa 2021
(DOI 10.1016/j.jcp.2021.110640), which sets the order of accuracy, and from Cada and
Torrilhon 2009 (DOI 10.1016/j.jcp.2009.02.020), whose limiter keeps it, and how closely.

The scheme is the kappa = 1/3 finite volume MUSCL reconstruction of van Leer and
Nishikawa's Eq 22, limited by Cada and Torrilhon's Eq 4.41, around the RoeM flux of
`RoeM1D.md`. Order is graded on their Section 4.2 Burgers case, which is the one test in
either paper whose time integrator is the solver's own. Monotonicity is graded on Sod's
classical initial data, and the radius of the smoothness indicator on air at 101325 Pa and
1.225 kg/m3. Gamma is 1.4 throughout.

## Order of accuracy

Their Section 4.2: Burgers with `u(x,0) = 1.5 + sin(2 pi x)` on `[0,1]`, periodic,
three-stage SSP Runge-Kutta, `dt = 1e-4` for 1000 steps to `t = 0.1`. Initial values are
cell averaged by their Eq 85 and graded against the cell-averaged exact solution of their
Eq 91 in the max norm of their Eq 93. The Courant number this produces on the finest grid
is 0.512, against the 0.512 they quote.

| Cells | Unlimited | Order | As shipped | Order | With the indicator off | Order |
|---|---|---|---|---|---|---|
| 127 | 5.0545e-4 | - | 5.0545e-4 | - | 3.2256e-3 | - |
| 255 | 6.5155e-5 | 2.96 | 6.5155e-5 | 2.96 | 1.1385e-3 | 1.50 |
| 511 | 8.2525e-6 | 2.98 | 8.2525e-6 | 2.98 | 3.9531e-4 | 1.53 |
| 1023 | 1.0326e-6 | 3.00 | 1.0326e-6 | 3.00 | 1.3450e-4 | 1.56 |
| 2047 | 1.2769e-7 | 3.02 | 1.2769e-7 | 3.02 | 4.3166e-5 | 1.64 |

The shipped column is bit identical to the unlimited one at every grid, the largest
difference between the two solutions being exactly zero. That is their Fig 2(a) and Table 1
result for MUSCL(1/3), reached with the limiter in place rather than with it removed; the
paper reaches it only limiter-free, stating that limiters "obscure the accuracy of the
underlying scheme". Bit identity rather than mere agreement follows from the shape of their
Eq 4.23, which is a nest of maxima and minima over `(2 + theta)/3` and its bounds: where the
data is smooth the nest returns that argument itself, unaltered, so no arithmetic happens.

The third column is the control. It is the same code with Cada and Torrilhon's smoothness
indicator switched off, leaving their Eq 4.23 alone, and it collapses into the same 1.5 to
1.7 band their Fig 7.3 reports in the max norm for that configuration. The indicator is the
load-bearing part.

Grading the same shipped runs against the *pointwise* exact solution of their Eq 92
instead gives orders of 2.37, 2.24, 2.14 and 2.07, converging on two rather than three.
That is their Pitfall 8 and 9 reproduced: a finite volume scheme compared against point
values reads second order however correct it is.

## Against the limiter's defining identities

Taken off the shipped reconstruction, by handing it a four-cell stencil whose difference
ratio is a chosen `theta` and reading the limiter back out of the face value.

| Source | Statement | Worst error |
|---|---|---|
| Eq 3.34 | the unlimited branch is `(2 + theta)/3` | 2.2e-16 |
| Eq 3.28 | `phi(1) = 1`, the condition for second order | 0.0 |
| Eq 3.28 | `phi'(1) = 1/3`, the condition for third order | 4.7e-11 |
| Eq 4.23 | `phi_hat(1) = 1` | 0.0 |
| Eq 4.23 | `phi_hat(-1) = 1/3`, which is third order at a symmetric extremum | 2.8e-16 |
| Eq 4.23 | `phi_hat(0) = 0`, the cutoff that resolves a discontinuity | 0.0 |
| Eq 4.23 | `phi_hat` reaches the bound `gamma = 1.6` as `theta` grows | 0.0 |
| Eq 4.23 | `phi_hat` vanishes as `theta` falls | 0.0 |
| Section 4.3 | `alpha = 0.5` holds third order across `theta` in `[-2, -0.8]` | 1.4e-16 |
| vLN Eq 22 | the unlimited branch equals the kappa = 1/3 weights | 2.2e-16 |

The `phi'(1)` figure is the central difference used to measure it, not the limiter.

Eq 3.38 is the separating check. A second-order TVD limiter satisfies
`phi(1/theta) = phi(theta)/theta`, which lets it reconstruct both faces from one
evaluation; a third-order one must not. Measured at `theta = 3`, `phi(1/3) = 0.777778`
against `phi(3)/3 = 0.555556`, apart by 2.2e-1, so the right face is built from the
inverse argument as their Eq 4.12 requires.

## The smoothness indicator

Their Eq 4.34 keyed to each conserved row, checked against the rates their Eq 4.35 to
4.39 derive.

| Source | Statement | Observed |
|---|---|---|
| Eq 4.35, 4.36 | `eta` vanishes as `dx^2` at a smooth extremum | order 2.00 |
| Eq 4.38, 4.39 | `eta` grows as `1/dx^2` at a jump | order -2.00 |

The reference step the indicator divides by is `R dx` times the row's acoustic size,
`(rho, rho a, rho a^2)`, reproduced exactly in all three rows.

## Monotonicity

Sod's problem, `(1, 0, 1)` against `(0.125, 0, 0.1)`, Courant 0.5, to `t = 0.15`. The exact
density and pressure are both monotone decreasing, so their total variation is the
end-to-end drop exactly, 0.875 and 0.900; anything above that, and any excursion outside
the initial range, is spurious. Neither measure needs an exact Riemann solve.

| Cells | Excess variation | Order | Worst overshoot | Order |
|---|---|---|---|---|
| 100 | 3.518e-2 | - | 1.518e-3 | - |
| 200 | 2.632e-2 | 0.42 | 1.303e-3 | 0.22 |
| 400 | 1.665e-2 | 0.66 | 1.022e-3 | 0.35 |
| 800 | 1.162e-2 | 0.52 | 7.431e-4 | 0.46 |

The limiter leaves Harten's TVD region for `theta < 0` by construction - that is what
`alpha = 0.5` in their Eq 4.19 buys, and it is why the scheme keeps third order at
extrema. Small overshoots are therefore expected rather than excluded, and what their
Section 7.3 claims for them is that they diminish under refinement. Both columns do,
monotonically, with the worst density excursion below 0.16 per cent at every resolution.

## Note: why the radius is dimensionless here

Their Eq 4.34 is `eta = (d- ^2 + d+ ^2) / (r dx)^2`, which carries the units of the
reconstructed variable, so `r` carries them too. Every case in the thesis is
non-dimensional with `rho` and `p` both near one, and a single `r` near one serves every
row. In SI it cannot: measured on a resolved one per cent acoustic pulse in air, the
largest difference each conserved row sees across the same wave is

| Row | Raw | Divided by the row's acoustic size |
|---|---|---|
| `rho` | 6.81e-2 | 5.56e-2 |
| `rho u` | 2.33e1 | 5.58e-2 |
| `rho E` | 1.98e4 | 1.39e-1 |

a spread of 290000 to 1 raw, and 2.5 to 1 once divided. The indicator therefore divides
by `(rho, rho a, rho a^2)` and the radius is the pure number `R = 1`.

That value is the one the Burgers study above lands on, and the two SI cases bracket it:
the largest reading anywhere on smooth data is 0.139, and the smallest at a genuine
discontinuity, measured on a five-to-one shock tube at equal temperature, is 3.37. `R = 1`
sits inside that window by a factor of seven on one side and 3.4 on the other.

## Note: the flux is RoeM1D's

The flux this reconstruction feeds is the one validated in `RoeM1D.md`, with the signal
velocities of that paper's Eq 33, and it satisfies all five of its exact flux identities to
the same bit. What the identities characterise is the flux function; composing it with a
reconstruction does not carry them over, because the smoothness indicator judges each
conserved row separately and across a contact the energy row is genuinely smooth while the
density row jumps. The boundary results in `Boundaries.md` were re-measured against this
reconstruction and are recorded there.

## Not covered

- Neelan and Nair 2022 (DOI 10.22055/jacm.2020.32845.2088), whose limiters were the other
  candidate. Their tables integrate with HRK42 and use Roe with a Harten entropy fix,
  neither of which is available here, so their columns are not reproducible and nothing is
  graded against them.
- Cada and Torrilhon's Section 6.1 stability headroom. Their von Neumann analysis puts
  this reconstruction with three-stage SSP Runge-Kutta stable to a Courant number near
  1.63, and their Euler cases run at 1.5; the solver assumes 1.0 and that has not been
  tested.
- Their Chapter 10 two-dimensional cases and Chapter 11 stiff relaxation systems, which
  need dimensions and source terms this solver does not have.
- The radius as a tuning parameter. Their Remark 4.5.1 recommends choosing it per problem,
  and their Euler cases span 0.01 to 10; here it is fixed at one and only the two cases
  above bound it.
