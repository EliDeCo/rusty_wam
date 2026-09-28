# Muscl3RoeM1D validation

What the reconstruction reproduces from van Leer and Nishikawa 2021
(DOI 10.1016/j.jcp.2021.110640), which sets the order of accuracy, and from Cada and
Torrilhon 2009 (DOI 10.1016/j.jcp.2009.02.020), whose limiter keeps it, and how closely.

The scheme is the kappa = 1/3 finite volume MUSCL reconstruction of van Leer and
Nishikawa's Eq 22, limited by Cada and Torrilhon's Eq 4.41, around the RoeM flux of
`RoeM1D.md`. Order is graded on their Section 4.2 Burgers case, which is the one test in
either paper whose time integrator is the solver's own. Monotonicity is graded on Sod's
classical initial data, and the radius of the smoothness indicator on air at 101325 Pa and
1.225 kg/m3. The Courant number is graded against their Section 6.1 von Neumann analysis,
which `TimeIntegrator::Ssp43` was added to reach. Gamma is 1.4 throughout.

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

## Order on the Euler equations

The table above grades the reconstruction through a scalar Burgers harness. This one grades the
whole shipped path — `Driver`, the RoeM flux, the reconstruction and `Ssp3` together — on a
problem that is smooth and exactly solvable.

A pure contact wave: air at 101325 Pa and 1.225 kg/m3 with uniform velocity 100 m/s and
uniform pressure, carrying a Gaussian density bump of 10 per cent amplitude and
`sigma = 0.08 m`, on a 1 m pipe with both ends non-reflecting. Every term of the Euler
equations cancels except `rho_t + u rho_x = 0`, so `rho(x,t) = rho_0(x - u t)` holds exactly
for all time while `u` and `p` stay constant. The pulse is advected 0.3 m, which keeps four
sigma of clearance at both ends so no wave ever reaches a boundary. Courant 0.5, graded in L1
against the cell-averaged exact solution and reported relative to the perturbation amplitude.

| Cells | L1 error | Observed order |
|---|---|---|
| 100 | 9.3723e-4 | - |
| 200 | 1.2060e-4 | 2.96 |
| 400 | 1.5140e-5 | 2.99 |
| 800 | 1.8961e-6 | 3.00 |
| 1600 | 2.3842e-7 | 2.99 |

Third order on the system, not only on the scalar case. The max norm over the same runs fits
2.98, so unlike the second order variant this scheme loses nothing at the extremum, which is
what the smoothness indicator is for.

Velocity and pressure hold their initial values to 2e-9 and 8e-10 relative on every grid, so the
wave stays the pure contact the exact solution assumes. That is the contact preservation the
RoeM flux exists for, measured end to end rather than off the flux function.

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

## Courant number

Their Section 6.1: the semi-discrete symbol of the reconstruction in closed form, and the
Courant numbers at which each integrator stays stable on it. Linear advection
`u_t + u_x = 0`, Fourier mode `exp(i j k)`, advanced one step by `R(L(nu,k))` where `L` is
the symbol of Eq 6.21 and `R` the stability polynomial of Eq 6.5 to 6.7. Stability is
`max|R| <= 1` over 4000 wavenumbers in `(0, pi]`, bisected on `nu` to 1e-12. `Ssp3` is
their three-stage SSP33 and `Ssp2` their Heun. What these limits are worth in compute time,
and which Courant number to actually run, is in `time_efficiency.md`.

| Reconstruction | Scheme | Ours | Theirs | Source |
|---|---|---|---|---|
| kappa = 1/3, as shipped | Ssp3 | 1.626 | 1.63 | Fig 6.2 left |
| kappa = 1/3, as shipped | Ssp43 | 2.037 | 2.04 | Fig 6.2 right |
| kappa = 1/3, as shipped | Ssp2 | 0.874 | 0.83 | Section 6.1 |
| kappa = 1/3, as shipped | Euler | unstable at every nu | not given | - |
| kappa = 0, the second order variant | Ssp2 | 1.000 | 1.0 | Fig 6.1 |
| kappa = 0, the second order variant | Ssp3 | 1.176 | not given | - |
| kappa = 0, the second order variant | Ssp43 | 1.601 | not given | - |
| kappa = 0, the second order variant | Euler | unstable at every nu | not given | - |
| phi = theta | Ssp2 | 0.500 | 0.5 | Fig 6.3 |
| phi = theta | Ssp3 | 0.628 | 0.63 | Fig 6.3 |
| phi = 2 theta | Ssp2, Ssp3 | unstable at every nu | absolutely unstable | Section 6.1 |

The last four rows are the controls. They are quoted for reconstructions this solver does
not use, and two of them land on exact round numbers, which is what establishes the symbol
is implemented correctly before the kappa = 1/3 rows are trusted. The Ssp43 row matches the
2.04 printed inside their Fig 6.2 rather than the "about 2.0" of its caption.

First-order upwind is not in their analysis and is added here because it is what `RoeM1D`
and `Roe1D` reduce to on this problem, which the solver table below is graded against:
1.000 with `Euler`, 1.000 with `Ssp2`, 1.256 with `Ssp3`, 2.000 with `Ssp43`.

The `Ssp2` row is the one that misses, 0.874 against 0.83, and the reason is that the
boundary there is soft rather than sharp. `max|R|` is 1.000000000 at `nu = 0.87` and only
1.000157 at `nu = 0.9`, so where the limit is declared depends on the tolerance chosen.
Both numbers describe the same curve.

### Accumulated error against Courant number

Their Eq 6.29 gives the per-step amplification error as `1 - (1/24) nu (2 + nu^3) k^4` for
`Ssp3` and `1 - (1/48) nu (4 + nu^3) k^4` for `Ssp43`. Over a fixed end time the step count
is `T/(nu dx)`, so the accumulated error carries `(2 + nu^3)` and `(4 + nu^3)/2`
respectively. Graded on linear advection, which is the problem Eq 6.29 was derived for:
`u(x,0) = sin(2 pi x)` on `[0,1]`, periodic, to `t = 2`, two full periods, max norm against
the exact cell average, each row normalised by its own `nu = 0.1` error on the 512-cell grid.

| nu | Ssp3 order | Ssp3 e/e(0.1) | Eq 6.29 | Ssp43 order | Ssp43 e/e(0.1) | Eq 6.29 |
|---|---|---|---|---|---|---|
| 0.10 | 3.00 | 1.00 | 1.00 | 3.00 | 1.00 | 1.00 |
| 0.50 | 3.00 | 1.06 | 1.06 | 3.00 | 1.03 | 1.03 |
| 0.90 | 3.00 | 1.36 | 1.36 | 3.00 | 1.18 | 1.18 |
| 1.26 | 3.00 | 2.00 | 2.00 | 3.00 | 1.50 | 1.50 |
| 1.30 | 3.00 | 2.10 | 2.10 | 3.00 | 1.55 | 1.55 |
| 1.59 | 3.00 | 3.01 | 3.01 | 3.00 | 2.00 | 2.00 |
| 1.60 | 3.00 | 3.05 | 3.05 | 3.00 | 2.02 | 2.02 |
| 1.70 | unstable | - | 3.45 | 3.00 | 2.23 | 2.23 |
| 2.10 | unstable | - | 5.63 | unstable | - | 3.31 |

Every stable entry matches Eq 6.29 to the three figures printed, and order stays 3.00
throughout. The two unstable rows bracket the ceilings of the table above from the other
side: `Ssp3` converges cleanly at 1.60 and diverges at 1.70, against the predicted 1.626,
and `Ssp43` converges at 1.70 and diverges at 2.10, against 2.037.

The `nu = 1.26` and `nu = 1.59` rows are there because `2^(1/3)` and `4^(1/3)` are where
each law doubles its `nu = 0` value, and both measure exactly 2.00. What that is worth when
choosing a Courant number to run is in `time_efficiency.md`.

Repeating the sweep on the nonlinear Section 4.2 Burgers case of the order table gives 2.9
to 3.0 at every stable `nu` as well, with ratios of 1.00, 0.98, 1.02, 1.16, 1.46 and 1.64
for `Ssp3` at `nu` of 0.1, 0.5, 0.9, 1.3, 1.6 and 1.7. Those run below the linear law
because the max norm there sits on a steepening profile where spatial truncation carries
more of the error, diluting the temporal share Eq 6.29 describes. The order is the claim
being graded; the ratios are reported for shape, not for agreement.

### Where the solver itself stops

The rows above are scalar. What ships is the Euler equations through the RoeM flux and the
real `Driver`, so its limit is measured separately: a 1e-3 relative standing acoustic mode
on 256 cells, `Wall` both ends, ten acoustic transits, bisected on the Courant number seven
times. At that amplitude the smoothness indicator holds the limiter on `phi_3`, so the
operator is the unlimited kappa = 1/3 one the analysis covers.

| Method | Scheme | Last stable | First unstable | Linear |
|---|---|---|---|---|
| Muscl3RoeM1D | Ssp3 | 1.622 | 1.644 | 1.626 |
| Muscl3RoeM1D | Ssp43 | 2.081 | 2.103 | 2.037 |
| Muscl3RoeM1D | Ssp2 | 1.141 | 1.163 | 0.874 |
| Muscl3RoeM1D | Euler | none | 0.200 | unstable |
| Muscl2RoeM1D | Ssp2 | 1.009 | 1.031 | 1.000 |
| Muscl2RoeM1D | Ssp3 | 1.272 | 1.294 | 1.176 |
| Muscl2RoeM1D | Ssp43 | 1.819 | 1.841 | 1.601 |
| Muscl2RoeM1D | Euler | none | 0.200 | unstable |
| RoeM1D | Euler | 0.988 | 1.009 | 1.000 |
| RoeM1D | Ssp2 | 1.009 | 1.031 | 1.000 |
| RoeM1D | Ssp3 | 1.250 | 1.272 | 1.256 |
| RoeM1D | Ssp43 | 2.016 | 2.038 | 2.000 |

Every bracket either contains the linear prediction or sits above it, never below, which is
the direction a limiter can move a ceiling: the analysis is of the unlimited symbol, and the
dissipation a limiter adds is not in it. The two that sit furthest above, `kappa = 1/3` with
`Ssp2` at 1.14 against 0.874 and `kappa = 0` with `Ssp43` at 1.82 against 1.601, are also
the two whose boundary is softest, so ten transits do not grow enough to trip the test.

`Muscl3RoeM1D` with `Euler` fails at every Courant number tried, down to 0.2, confirming the
analysis: forward Euler's stability region touches the imaginary axis only at the origin,
and this reconstruction has no dissipative part to pull the symbol off it. The failure is
gradual rather than immediate, which is what makes it worth stating -- a run can proceed for
thousands of steps before diverging.

With a shock the ceiling barely moves. Sod on 400 cells to `t = 0.15` completes at every
Courant number up to 1.60 with `Ssp3` and 1.80 with `Ssp43`, and fails at 1.80 and 2.00
respectively. The positivity fallback of `enforce_positivity` fires zero times at every one
of those settings, so on this problem it is a guard that is never reached rather than a
correction the scheme leans on.

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

## Note: how the Ssp43 coefficients were obtained

Their printed Eq 6.3 is not usable as set: its third-stage coefficients 2/3 and 1/2 sum to
7/6, so the stage is not a convex combination and the scheme it describes is not
conservative. The coefficients shipped were recovered from Eq 6.7 instead, which gives the
stability polynomial `1 + z + z^2/2 + z^3/6 + z^4/48`. Expanding the standard Kraaijevanger
four-stage third-order scheme gives `(2a + a^4)/3` with `a = 1 + z/2`, which is that
polynomial exactly, fixing the stages as

```
u1 = u^n                 + (1/2) dt L(u^n)
u2 = u1                  + (1/2) dt L(u1)
u3 = (2/3)u^n + (1/3)u2  + (1/6) dt L(u2)
u4 = u3                  + (1/2) dt L(u3)
```

Driving the shipped `TimeIntegrator::stages()` weights through the same scalar symbol the
table above uses reproduces Eq 6.5 to 6.7 for all four integrators, with the worst deviation
2.24e-16 for `Ssp43` and 1.24e-16 for the others, so the code path and the printed
polynomials agree to rounding. Order on the Burgers and linear cases is third, independently
of the thesis.

## Note: the junction Courant number is not covered by any of this

`Junction::get_timestep` is `courant * 2 V / sum(A (|u.n| + a))`, a control volume bound over
the interfaces rather than a stencil bound, and neither this thesis nor the ghost junction
paper analyses it. The same `courant` constant feeds it, so raising the Courant number on the
strength of the tables above raises it at every junction too, where nothing here says what the
limit is. Everything above was measured on single pipes.

## Not covered

- Neelan and Nair 2022 (DOI 10.22055/jacm.2020.32845.2088), whose limiters were the other
  candidate. Their tables integrate with HRK42 and use Roe with a Harten entropy fix,
  neither of which is available here, so their columns are not reproducible and nothing is
  graded against them.
- Their Chapter 10 two-dimensional cases and Chapter 11 stiff relaxation systems, which
  need dimensions and source terms this solver does not have.
- The radius as a tuning parameter. Their Remark 4.5.1 recommends choosing it per problem,
  and their Euler cases span 0.01 to 10; here it is fixed at one and only the two cases
  above bound it.
