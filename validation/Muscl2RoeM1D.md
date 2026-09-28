# Muscl2RoeM1D validation

What the second order reconstruction reproduces from van Leer and Nishikawa 2021
(DOI 10.1016/j.jcp.2021.110640), and the properties it holds in its own right.

The scheme is the `kappa = 0` finite volume MUSCL reconstruction of their Eq 22 — Fromm's
scheme, and the choice the original MUSCL code was built on — limited by the monotonized
central slope limiter, around the RoeM flux of `RoeM1D.md`. It shares its stencil, its
positivity fallback and its flux with `Muscl3RoeM1D.md`; only the rule that limits a cell
average to its faces differs. Gamma is 1.4 throughout.

One claim in either paper is gradeable here, and it is the order: their Section 1 and
Section 2.4 state that every scheme of the kappa family is second order except `kappa = 1/3`.
They test limiter-free throughout, saying limiters "obscure the accuracy of the underlying
scheme", so that is how the claim is graded below. Everything else in this document is a
property checked against its own definition, not against a published number. See
**Not covered** for why the one published second-order limiter study cannot be used.

## Order of accuracy

Their Section 4.2: Burgers with `u(x,0) = 1.5 + sin(2 pi x)` on `[0,1]`, periodic,
three-stage SSP Runge-Kutta, `dt = 1e-4` for 1000 steps to `t = 0.1`. Initial values are
cell averaged by their Eq 85 and graded against the cell-averaged exact solution in the max
norm of their Eq 93, which is what their Pitfalls 8 and 9 require. The same integrator and
the same case as the third order table, so the two differ only in the reconstruction.

| Cells | Unlimited, Linf | Order | As shipped, Linf | Order | As shipped, L1 | Order |
|---|---|---|---|---|---|---|
| 127 | 2.5549e-3 | - | 2.7189e-3 | - | 3.3797e-4 | - |
| 255 | 6.5643e-4 | 1.96 | 1.0812e-3 | 1.33 | 8.6306e-5 | 1.97 |
| 511 | 1.6663e-4 | 1.98 | 4.6325e-4 | 1.22 | 2.1961e-5 | 1.97 |
| 1023 | 4.1975e-5 | 1.99 | 1.9679e-4 | 1.24 | 5.5057e-6 | 2.00 |
| 2047 | 1.0544e-5 | 1.99 | 7.9262e-5 | 1.31 | 1.3866e-6 | 1.99 |

The unlimited column is the graded one and converges on two, which is their claim for
`kappa = 0`. The third order variant on the same case and the same grids gives 2.96 to 3.02,
so the two columns reproduce both halves of their statement about the kappa family.

The shipped column reads 1.2 to 1.3 in the max norm and 2.0 in L1, and that gap is the
scheme working as designed rather than failing. A TVD limiter clips smooth extrema; this
case is a sine, so its two extrema are clipped to first order while the rest of the profile
stays second order. The max norm sits on the clipped extremum and reports its order, L1
averages over a domain where almost every cell is still second order. The third order
variant exists precisely to avoid this, which is what its smoothness indicator buys, and the
`phi-hat alone` column of `Muscl3RoeM1D.md` shows the same collapse when that indicator is
switched off.

## The limiter against the properties a TVD limiter must have

Taken off the shipped `limited_slope` by handing it a near difference of one and a far
difference of `r`, so that what comes back is `phi(r)` itself. Swept over `r` in `[-4, 4]`
at a spacing of 0.01, plus the reciprocal of every positive value so the symmetry arm is
exercised from both sides.

| Property | Observed | Verdict |
|---|---|---|
| `phi(r) = 0` for `r <= 0` | holds at every r | exact |
| `0 <= phi(r) <= min(2r, 2)`, Harten's TVD region | holds at every r | exact |
| Eq 3.38 symmetry, `phi(1/r) = phi(r)/r` | 0.0e0 | exact |
| `phi(1) = 1`, second order on smooth data | 1.000000 | exact |

`phi(0.5) = 0.75`, `phi(1) = 1`, `phi(2) = 1.5`, `phi(-1) = 0`.

The symmetry is exact rather than approximate because the limiter is written as a function
of the two differences that is symmetric in them, which is the same statement as Eq 3.38
without ever forming `r`. That also means one evaluation is valid for both faces of a cell.
Unlike the third order limiter, this one stays inside Harten's region everywhere, including
`r < 0`, which is the whole point of the variant.

## Monotonicity

Sod's problem, `(1, 0, 1)` against `(0.125, 0, 0.1)`, Courant 0.5, to `t = 0.15`. The exact
density is monotone decreasing, so its total variation is the end-to-end drop exactly,
0.875; anything above that, and any excursion outside the initial range, is spurious.

| Cells | Excess variation | Order | Worst overshoot |
|---|---|---|---|
| 100 | 3.840e-2 | - | 0.000e0 |
| 200 | 2.332e-2 | 0.72 | 0.000e0 |
| 400 | 1.706e-2 | 0.45 | 0.000e0 |
| 800 | 1.239e-2 | 0.46 | 0.000e0 |

**The overshoot is exactly zero at every resolution**, which is what the variant is for and
what the third order scheme cannot offer: on the same case and the same harness the third
order scheme overshoots by 1.559e-3 down to 7.379e-4, because its limiter leaves Harten's
region for `theta < 0` by construction. Excess variation is comparable between the two,
1.239e-2 against 1.203e-2 on the finest grid, so what this scheme buys is the bound, not a
smoother answer. First order `RoeM1D` also never overshoots, as it must not.

## Courant number

`Muscl2RoeM1D` pairs with `Ssp2`, which `MethodKind::integrator` selects for it. The
`kappa = 0` rows of the Section 6.1 analysis in `Muscl3RoeM1D.md` give the linear ceilings,
and the same standing acoustic mode measures what the solver actually does with them:

| Scheme | Linear ceiling | Last stable | First unstable |
|---|---|---|---|
| `Ssp2` (paired) | 1.000 | 1.009 | 1.031 |
| `Ssp3` | 1.176 | 1.272 | 1.294 |
| `Ssp43` | 1.601 | 1.819 | 1.841 |
| `Euler` | unstable at every nu | none | 0.200 |

The paired row lands on its linear ceiling to within 1 per cent. The rows above it are
measured higher than the analysis predicts because the analysis is of the *unlimited*
`kappa = 0` symbol, and the limiter adds dissipation the symbol does not carry.

`Euler` is unstable at every Courant number with this reconstruction exactly as it is with
`kappa = 1/3`, which is why `Driver::with_integrator` refuses the pairing for both.

## Not covered

- Neelan and Nair 2022 (DOI 10.22055/jacm.2020.32845.2088) is the only published study of a
  second order limiter in either reference, and none of it is graded here. Their MMF1 of
  Eq 11 is compared against minmod, superbee and van Albada in their Tables 2, 3 and 4, but
  every run is integrated with HRK42, their own hyperbolic Runge-Kutta, and uses Roe with a
  Harten entropy fix. Neither is implemented here, so their columns cannot be reproduced and
  nothing is measured against them. Their Eq 14 `minmod_s2` is not a candidate either: it
  limits the second derivative of their *third* order scheme, not a slope.
- The choice of `kappa`. Any value in `(-1, 1)` except 1/3 is second order, and `kappa = 0`
  is taken because it is Fromm's scheme and the original MUSCL choice. Their Table 4 notes
  that `kappa = 1/2` reaches third order in the *steady state*, in point values rather than
  cell averages; that is untested here.
- The limiter as a tuning parameter. Monotonized central was selected by measuring the cost
  of reaching one per cent against minmod and van Leer, recorded in `time_efficiency.md`.
  Only the winner ships, and the other two are not in the code.
