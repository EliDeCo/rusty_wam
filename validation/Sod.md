# Sod's shock tube validation

What the solver reproduces on the benchmark of Sod 1978
(DOI 10.1016/0021-9991(78)90023-2), graded against the analytic solution of the same
Riemann problem rather than against any of the schemes he surveyed.

His Section 3 test problem: `(rho, u, p) = (1, 0, 1)` left of a diaphragm at `x = 0.5` and
`(0.125, 0, 0.1)` right of it, gamma 1.4, run to `t = 0.2`. The pipe is 1 m of uniform area
with both ends non-reflecting, and the run is `Muscl3RoeM1D` with `Ssp3` at a Courant number
of 1.26, which is the pairing `time_efficiency.md` selects at scale 1.0. The last timestep is
clipped so the state is read at exactly `t = 0.2` rather than at the first step past it.

Three fields are graphed against distance — density, velocity and pressure — which is the
first four panels of his Figs 4 to 15 less the energy one. Computed cells are drawn as
markers up to a thousand of them and as a line above that, where they would otherwise merge
into a band; the analytic curve is drawn dense through the rarefaction and from the two ends
of every constant region, so each jump stays vertical on the page.

## The analytic solution, graded first

The exact solution is the reference everything else here is measured against, so it is
checked against the textbook values before any solver output is compared to it. The star
region comes from matching the two wave curves at a common pressure; the rest is sampling.

| Quantity | Computed | Standard value |
|---|---|---|
| `p*`, contact pressure | 0.30313 | 0.30313 |
| `u*`, contact velocity | 0.92745 | 0.92745 |
| `rho_3`, left of the contact | 0.42632 | 0.42632 |
| `rho_4`, right of the contact | 0.26557 | 0.26557 |
| shock speed | 1.75215 | 1.75216 |

At `t = 0.2` the five wave edges sit at `x =` 0.26336, 0.48595, 0.68549, 0.85043, 0.85043.
The last two coincide because a shock has no width, which is the sampler reporting the
degenerate case correctly rather than a duplicated number. The rarefaction head travels at
`-a_1 = -1.18322` and its tail at `u_3 - a_3 = -0.07025`, both recovered to five figures.

Nothing reaches either end of the pipe before `t = 0.2`, so the boundary treatment does not
enter this result.

## What one per cent means here

Error is **relative L1**, per field:

```
sum_j |q_j - qbar_j|  /  sum_j |qbar_j|
```

against the exact solution **cell averaged** over each cell, not sampled at its centre. A
finite volume state holds cell averages, and centre sampling would misplace every jump by up
to half a cell — an error of the same order as the scheme's own, and worth up to a seventh of
the one per cent budget at these resolutions.

A max norm is not used and cannot be. Sod's solution contains a contact and a shock, so
pointwise error at a discontinuity is floored at the size of the jump — 38 per cent at the
contact — however fine the grid. Relative L1 converges; a max norm never falls below that
floor.

## Convergence

| Cells | Density | Order | Velocity | Order | Pressure | Order |
|---|---|---|---|---|---|---|
| 50 | 2.0393% | - | 5.1484% | - | 2.1748% | - |
| 100 | 1.1932% | 0.77 | 2.9553% | 0.80 | 1.2219% | 0.83 |
| 200 | 0.6636% | 0.85 | 1.5240% | 0.96 | 0.6194% | 0.98 |
| 400 | 0.3446% | 0.95 | 0.7597% | 1.00 | 0.3099% | 1.00 |
| 800 | 0.1772% | 0.96 | 0.3559% | 1.09 | 0.1530% | 1.02 |
| 1600 | 0.1016% | 0.80 | 0.1949% | 0.87 | 0.0810% | 0.92 |
| 3200 | 0.0604% | 0.75 | 0.1159% | 0.75 | 0.0462% | 0.81 |

**The order is about one, not three, and that is correct.** Every scheme is first order at a
discontinuity whatever its order on smooth data, and two of this solution's four features are
discontinuous. The third order reconstruction is graded on smooth data in
`Muscl3RoeM1D.md`, where it reads 2.96 to 3.02; what it buys here is a smaller error
constant, not a steeper slope. The drift below 1.0 at the top of the table is the contact,
which has no self-steepening and keeps spreading as the run proceeds.

## Cells to reach one per cent

**316 cells**, at which the three errors are 0.3988%, 0.7558% and 0.3657% in 110 timesteps.

That figure is the smallest count above which **every** grid tested up to 420 stays under one
per cent. It is not the first count that passes, which is 262, and the difference matters:

| | Cells |
|---|---|
| first count that happens to pass | 262 |
| smallest count above which every tested grid passes | 316 |

**Relative L1 is not monotone in the cell count.** Scanning every even count from 240 to 420
shows velocity wandering between 0.72% and 1.28% with no clean crossing — 274 cells is worse
than 250, and 286 is worse than 284. The scatter is about ±0.15 percentage points, which is
the same size as the distance from the threshold in this range. A bisection assumes a
monotone error and would report whichever grid it happened to land on; the scan is what makes
the answer defensible. The cause is that the shock's position relative to the cell faces at
`t = 0.2` changes with the grid, and so does how much of its smearing the norm sees.

## Velocity binds, not density

Velocity carries roughly twice the relative L1 error of the other two at every resolution, so
it is the field that sets the cell count.

That is the opposite of what the feature count suggests. Pressure and velocity are continuous
across the contact — only density jumps there — so velocity sees three features to density's
four. The jump sizes are what decides it instead: **at the shock, velocity falls by the whole
0.92745 while density falls by 0.141**, six times less. Velocity's single discontinuity
carries far more of the norm than density's two together.

The norm's denominator compounds it. Velocity is exactly zero over the 41 per cent of the
domain the waves have not reached, so `sum |qbar|` is smaller for velocity than for density
and the same absolute error reads as a larger relative one. Relative L1 is therefore a
harsher bar for a field that is zero over much of its domain, which is worth knowing before
comparing this number against one quoted elsewhere.

## Not covered

- **Every interior method except `Muscl3RoeM1D`.** `Muscl2RoeM1D`, `RoeM1D` and `Roe1D` run
  the same problem through the same path but none has been graded on it, so no cell count is
  claimed for them.
- **Every time except `t = 0.2`.** The solution is analytic in time and the runner takes any
  end time, but only 0.2 is graded here.
- **The three other Riemann problems** already defined in `helpers.rs` — the moving shock
  tube, the pure contact and the supersonic expansion of the RoeM paper's Section 5.1. The
  exact solver handles all three; none is wired to a `Problem` variant yet.
- **The boundary conditions.** No wave reaches either end inside `t = 0.2`, so this grades the
  interior scheme alone. `Boundaries.md` is where the ends are graded.
- **Saving the graphs as image files.** They are drawn in a window and not written to disk.
