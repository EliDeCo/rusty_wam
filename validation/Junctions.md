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
