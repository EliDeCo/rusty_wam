# Junction validation

What the solver reproduces from Hong & Kim 2011 (DOI 10.1002/fld.2212), and how closely.

Air at 101325 Pa and 1.225 kg/m3, inlet Mach 0.05, 64 cells per branch, first-order RoeM
with third-order SSP Runge-Kutta, run to a relative residual of 1e-9. Results are
grid-independent to five decimals over 64 to 256 cells. Figure values were read off the
paper's plots as fitted curves.

## Against Table I

At a 90 degree branch angle, where the paper's own results sit on the correlations.

| Source | Coefficient | Swept | Worst error |
|---|---|---|---|
| Type 1 | K13, counter-combining | split, angle, area ratio | 0.3% |
| Type 2 | K31, counter-dividing | split, area ratio | 0.5% |
| Type 3 | K12, straight-branch combining | split, angle | 0.1% |
| Type 3 | K32, straight-branch combining | split | 0.6% (see note) |
| Type 4 | K13, straight-branch dividing | split | 0.2% |
| Figure 12b | K21, sudden expansion | area ratio 0.1 to 0.9 | 0.3% |

## Against the figures

At 45 and 135 degrees, and for Type 4's K12, the paper's own results depart from Table I by
up to 0.55 in K. Sections 5.3 and 5.4.1 describe this. The solver follows the paper rather
than the correlation.

| Figure | Coefficient | Worst departure from the paper's plotted points |
|---|---|---|
| 14a | K31, counter-dividing, 45 and 135 degrees | 0.023 |
| 18a | K12, straight-branch dividing, all angles | 0.030 |
| 18b | K13, straight-branch dividing, 45 and 135 degrees | 0.041 |

Across those 45 sample points the mean departure is 0.013 in K, against plot axes spanning
-1 to 3.5. Figure 13a is matched to within 0.23% across the whole split range.

## Against the scaling function

| Source | Quantity | Result |
|---|---|---|
| Figure 3b | chi with G, Mach 0.01 to 0.5 | 0.985 to 0.999 |
| Section 2.3 | chi without G at Mach 0.01 | 80.6, against "more than eighty times" |

## Courant number

`Junction::get_timestep` is a control volume bound, not a stencil bound, and the paper gives
it only implicitly. Their LU-SGS junction diagonal, Eq 46b, is
`V_j/dt + (1/2) sum_i s_ji rho(A_j)` with `rho(A) = C(|U| + a)` from Eq 43b; balancing the two
terms gives `dt = 2V / sum(A(|u.n| + a))`. Eq 44b does the same for a pipe and collapses to
`dt = dx/(|u| + a)`. So the junction bound is the direct analogue of the pipe bound under
their own convention, and that is where the factor of two comes from.

Swept on `N` branches merging into one outlet held at Mach 0.05, 64 cells per branch,
first-order RoeM with three-stage SSP Runge-Kutta, which is their Eq 42. Stability is judged
on whether the run stays physical, not on whether it reaches the tolerance, because a run
that merely exhausts its iteration budget is slow rather than unstable.

### The junction has its own limit, and it is not the pipes'

Holding the junction bound out of the way and sweeping the **pipe** number finds where the
network stops. The last column is the junction's own Courant number at that point, in the
Eq 46b units above, which is `dt sum(A(|u.n|+a)) / 2V`:

| Branches | Last stable | First unstable | Junction nu there | Limited by |
|---|---|---|---|---|
| 3 | 1.397 | 1.418 | 1.83 | pipe |
| 4 | 1.376 | 1.397 | 2.48 | pipe |
| 7 | 1.186 | 1.207 | 3.32 | junction |
| 9 | 0.785 | 0.806 | 3.17 | junction |
| 11 | 0.616 | 0.637 | 2.92 | junction |
| 13 | 0.490 | 0.511 | 2.67 | junction |
| 17 | 0.363 | 0.384 | 2.79 | junction |

A single pipe with the same method and integrator holds to 1.250 and fails at 1.272. Three
and four branches beat that, so there the junction is not what stops the network. **From
seven branches on it is**, and the network ceiling falls away steeply — by seventeen branches
it is a quarter of what one pipe alone manages.

The last column is what makes that readable: whatever the branch count, the junction gives
out between **2.7 and 3.3**. That is a property of the junction cell rather than of the
network, and it is the number the shipped bound is set against. `Junction::courant_from`
takes 2.0, which is under the lowest observed failure with margin.

### Why it needs a number of its own rather than the pipes'

Their Section 2.1.2 sets the junction volume to the mean of the neighbouring cell volumes,
which is what `initialize` does. For `N` branches of equal area that makes the junction volume
`A dx`, so at the *same* Courant number as the pipes the two bounds sit in a fixed ratio of
`2/N` whatever the flow is doing:

| Branches | dt_junction / dt_pipe | 2/N |
|---|---|---|
| 2 | 1.024 | 1.000 |
| 3 | 0.682 | 0.667 |
| 4 | 0.512 | 0.500 |
| 5 | 0.408 | 0.400 |

Handed the pipes' number, a three-branch junction therefore takes two thirds of their step
while sitting at a Courant number of only 1.9 against a limit of 2.7 — it gives up a third of
the step for nothing. Giving it 2.0 of its own fixes both ends at once. With the pipes at
their optimum of 1.26:

| Branches | dt_junction / dt_pipe at nu_j = 2.0 |
|---|---|
| 3 | 1.107 |
| 4 | 0.832 |
| 7 | 0.476 |
| 9 | 0.370 |
| 13 | 0.256 |
| 17 | 0.196 |

Up to three branches the junction now clears the pipe bound and stops costing anything; from
four on it binds, which is exactly where the table above says it must.

### Their Section 2.1.2 quasi-steady claim

They state that "the volume size and time derivative of the ghost junction cell do not have a
significant impact on the flow around the junction if the flow quickly reaches an equilibrium
state", which is why they discard the `-(dV/dt)Q` term of Eq 7a. Both halves hold, on the
counter-dividing T of Figure 2b at split 0.5, area ratio 1, 90 degrees, Mach 0.05:

| Varied | Range | Worst change in K31 | Effect on iterations |
|---|---|---|---|
| Junction Courant | 0.1 to 3.6 | 3.2e-7 | 275147 down to 21542 |
| Junction volume | 0.25x to 4x | 2.7e-7 | 123483 down to 22206 |

K31 is 1.2433 in every one of those runs. Both knobs move only the timestep, and once either
is large enough for the pipe bound to take over they saturate together at 21542 iterations.

### The same claim on a transient

Their Figures 5(b) and 7(b) show third-order TVD Runge-Kutta giving the same answer at CFL
1.0 and 0.1, but those two figures are a single tube. Carried across to a junction: one branch
of the T closed so the flow leaving along it slams shut at `t = 0`, ten acoustic transits,
total mechanical energy of their Eq 41b sampled at forty times.

| Comparison | Worst relative difference |
|---|---|
| Junction Courant 0.1 against a junction that never binds | 1.1e-5 |

Their two curves overlap on the page, so a per-cent statement would have been the honest bar;
these agree three orders of magnitude inside it. What the junction Courant buys is therefore
speed and nothing else, up to the limit above.

### What this is worth

Every graded coefficient in this document is unchanged by the switch — all 712 printed values
across the Table I and figure comparisons are identical — while the runs reach the same
residual in fewer steps:

| | Iterations, mean over 49 runs | Speedup |
|---|---|---|
| Junction handed the pipes' Courant number | 55291 | - |
| Junction on its own 2.0 | 39215 | 1.40x mean, 1.49x worst case |

The worst case is the 1.5 that `2/N` predicts for the three-branch configurations this
document is built on, and the runs showing no change are the ones where the pipe bound was
already the smaller of the two.

## Note: the Type 3 K32 entry in Table I

Table I prints `K32'' = -1 + 4q - (xi^2 - 2 xi cos(theta) - 2) q^2`. The solver misses that
by 3 to 48 per cent, but matches the same expression with a **plus** on the quadratic term
to within 0.6 per cent. The plus form is taken as correct on three independent grounds:

- A one-dimensional momentum balance across the junction yields the plus form, and that
  same balance reproduces the neighbouring K12 entry character for character.
- In Figure 16b the 45 and 90 degree curves are concave down and 135 is very slightly
  concave up. The printed form gives the opposite curvature on all three.
- In Figure 16b the 45 degree curve peaks near a split of 0.83. The plus form places a
  turning point there; the printed form has no turning point in any curve over the domain.

The paper's own plotted data therefore agrees with the solver and not with its printed
equation, so the discrepancy is a misprint rather than a calibration error.

## Not covered

- Branch counts above seventeen, and junctions whose branches differ in area or whose volume
  departs from the mean-of-neighbours rule. The limit above is flat enough in branch count to
  be worth trusting between three and seventeen, but it was measured on equal-area branches.
- The junction limit under the second and third order interior methods. It was measured with
  first-order RoeM and Ssp3, which is the pairing their Eq 42 and Chapter 5 use; the pipes
  pair differently now, and whether a reconstructing branch changes what the junction can take
  is untested.

