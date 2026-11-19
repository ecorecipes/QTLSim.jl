# Why the emulator's elasticity in ω looked three times too strong at narrow ω

Exploratory diagnosis (2026-11-18/19), **not a verdict study**: no pass/fail criteria. Each
hypothesis and prediction was written in `validate/narrow_omega.log` before its test, every test
run is reported there, seeds are fixed, and the two design changes made along the way are marked.
Script `validate/narrow_omega.jl`; caches in `validate/narrow_omega_raw/` (git-ignored). No change
to `src/`, `test/` or defaults. 5.44 CPU-hours (`clock()`), six threads.

## The puzzle

Phase B of the elasticity study (`validate/elasticity_b.log`) compared the elasticity of the
half-survival rate k50 in ω = `InitFitnessSD`, cell by cell. In the lowest tercile of ω the
two-stage screen (`loci = :spectrum, spectrum = (; ld = :dynamic, solver = :beta)`, by
`k50_gradient`) gave cell medians of −0.63, −0.55, −0.43, −0.47 (L = 100, 500, 1000, 5000) against
the table's tercile fits of −0.17, −0.19, −0.18, −0.14: sign right, size FAIL in all four cells.

## Answer

**The table summary is what is off, for most of the gap; the emulator carries a smaller, real
excess at L ≥ 500.**

1. The table's tercile slope is not a local elasticity. Inside the "low tercile" (ω from 0.003
   to 2.5–3.3) the table's own e_ω has two branches: positive below ω ≈ 0.2–0.3 and about −0.4 to
   −0.5 above. The tercile fit averages them, with 65–69% of its leverage on the positive branch.
2. On the table's own summary the emulator agrees with the table: −0.21, −0.20, −0.20 against
   −0.17, −0.19, −0.18 (L = 100, 500, 1000).
3. The kernel's own local elasticity, measured directly on 12 inputs with ω from 0.35 to 2.1, is
   −0.20 to −0.61, mean −0.40. It matches the table's fits restricted to ω ≥ 0.3 and is nowhere
   near −0.17.
4. Against that direct measurement the screen is right at L = 100 (difference 0.00) and too
   negative by 0.10 at L = 500 and 0.14 at L = 1000, concentrated in five inputs with L ≥ 500 and
   ω of 0.65–1.8 (−0.14 to −0.33). Of that, 0.07–0.10 in two inputs is the beta solver.

So phase B's "size FAIL" at low ω is mostly a mismatch of summaries (a gap of about 0.3–0.45
shrinks to 0.0–0.14), plus a genuine emulator excess of about 0.1–0.15 at L ≥ 500 on the negative
branch. L = 5000 was not measured with the kernel; there the table's ω ≥ 0.3 fit (−0.20) is still
0.27 from phase B's median (−0.47).

## Established

### The table alone

- **No pre-move deaths.** No row has T ≤ 0. For ω < 0.15, 813 of 1776 rows die within 5
  generations of the move (273 of 2056 at 0.15–0.3, none above 0.5), never at k < 0.01: there the
  first steps of the optimum kill, and a wider ω helps.
- **The tercile fit on sub-ranges of ω** (same model; delta-method SE 0.02–0.06):

  | L | whole tercile | ω < 0.3 | 0.3 ≤ ω | 0.5 ≤ ω | 1 ≤ ω | 0.5 ≤ ω < 1.5 |
  |---|---|---|---|---|---|---|
  | 100 | −0.174 | +0.290 | −0.511 | −0.518 | −0.503 | −0.499 |
  | 500 | −0.191 | +0.201 | −0.422 | −0.454 | −0.427 | −0.512 |
  | 1000 | −0.179 | +0.172 | −0.369 | −0.396 | −0.335 | −0.426 |
  | 5000 | −0.138 | +0.130 | −0.201 | −0.206 | −0.145 | −0.294 |

  Rows with ω < 0.3 are 13% of the tercile's rows and 14–20% of its information, but carry 65–69%
  of the leverage of a slope in log ω.
- **Finer bins** (15 fixed bins, a piecewise-linear spline in log ω per locus number): +0.30 to
  +0.42 below 0.15, then negative; single-bin slopes have SE 0.15–0.4 and are too noisy to use one
  by one. The sub-range fits above are the usable summary.

### The emulator on the table's terms

- Phase B's 126 screen elasticities re-binned follow the same two branches: median +0.28 below
  ω = 0.15, about 0 at 0.15–0.3, −0.36 at 0.3–0.5, −0.46 to −0.60 from 0.5 to 2.2.
- The screen's k50(ω) at reference inputs (V_P 0.0223, h² 0.315, m 0.004, u = 3) has a maximum at
  ω ≈ 0.15–0.2, falls with e of −0.49 to −0.63 between 0.75 and 1.7, and turns up again near
  ω = 4–4.5. With the table's rows and weights, its tercile-style slope is −0.213 / −0.200 /
  −0.202 (table −0.174 / −0.191 / −0.179); over ω ≥ 0.3 it is −0.56 / −0.56 / −0.54 (table −0.51 /
  −0.42 / −0.37).
- **Level.** Screen k50 over the table model's at the reference inputs: 0.81–0.95 at L = 100 for
  ω ≥ 0.15; at L = 500 1.16–1.26 for ω of 0.2–0.55; at L = 1000 1.39–1.51 for ω of 0.2–0.55; back
  to 0.92–1.16 from ω ≈ 1.

### The kernel's local elasticity (12 inputs, arms ω·e^{±0.2})

Kernel: equilibrium start, 1000 generations of burn-in, 21 rates × 10 runs per arm, same seeds in
both arms; Spearman–Kärber estimate, paired bootstrap 95% CI. Emulators: finite difference of
their k50 over the same arms.

| # | input | L | ω | kernel e_ω [95% CI] | screen | grid + LD | no LD | screen/kernel k50 (ω−, ω+) |
|---|---|---|---|---|---|---|---|---|
| 1 | QT0100-5/07137 | 100 | 0.36 | −0.43 [−0.49, −0.38] | −0.31 | | | 0.80, 0.84 |
| 2 | QT0100-5/12231 | 100 | 0.93 | −0.53 [−0.58, −0.49] | −0.54 | | | 0.81, 0.81 |
| 3 | QT0100-3/03404 | 100 | 1.37 | −0.60 [−0.66, −0.54] | −0.62 | | | 0.89, 0.89 |
| 4 | QT0100-5/14789 | 100 | 2.14 | −0.61 [−0.65, −0.56] | −0.68 | | | 0.77, 0.74 |
| 5 | QT0500-3/08249 | 500 | 0.54 | −0.43 [−0.47, −0.40] | −0.46 | −0.46 | −0.72 | 1.20, 1.19 |
| 6 | QT0500-3/04007 | 500 | 0.69 | −0.39 [−0.42, −0.35] | −0.38 | −0.38 | | 1.20, 1.20 |
| 7 | QT0500-3/09773 | 500 | 1.37 | −0.39 [−0.43, −0.35] | −0.61 | −0.55 | −0.74 | 1.07, 0.98 |
| 8 | QT0500-5/08615 | 500 | 1.82 | −0.43 [−0.46, −0.38] | −0.57 | | | 0.94, 0.89 |
| 9 | QT1000-5/11094 | 1000 | 0.35 | −0.26 [−0.29, −0.23] | −0.18 | −0.18 | −0.23 | 1.29, 1.34 |
| 10 | QT1000-5/14211 | 1000 | 0.65 | −0.35 [−0.38, −0.31] | −0.68 | −0.58 | | 1.14, 1.00 |
| 11 | QT1000-5/13607 | 1000 | 1.47 | −0.43 [−0.46, −0.40] | −0.61 | −0.61 | −0.91 | 1.22, 1.13 |
| 12 | QT1000-5/09914 | 1000 | 1.55 | −0.20 [−0.25, −0.16] | −0.38 | | | 0.99, 0.93 |

- Kernel: inverse-variance mean −0.40; by L −0.54, −0.41, −0.32 (table's 0.3 ≤ ω fits: −0.51,
  −0.42, −0.37). Screen − kernel: mean −0.08, median −0.05; by L 0.00, −0.10, −0.14; 8 of 12 screen
  values lie outside the kernel's CI.
- The kernel's survival is nearly a step in the rate per input (logistic width 0.008–0.031 in
  log k).
- `:rare_alleles` has the wrong sign or size (+0.31 to −0.49). Without dynamic LD the spectrum
  closure is more sensitive still (−0.72 to −0.91 in three of four inputs).
- The beta solver against the grid: identical finite differences in inputs 5, 6, 9, 11; 0.07 and
  0.10 more negative in 7 and 10, where the screen's local value also differs from its finite
  difference by 0.13–0.17 (k50(ω) of the screen is not smooth over ±20% there).

### Which links are right

- **State at the move.** Screen V_G over the table's own "Genotype Var at 10000" on 660 rows:
  1.02–1.05 for every ω bin from 0.3 to 2.2 (1.12 at 0.15–0.3, 1.23 below 0.15); the slope of
  log V_G on log ω over 0.3–2.2 is +0.80 (0.03) for the screen and +0.82 (0.04) for the table.
  Against the kernel on the 12 inputs: 0.98–1.04 in ten, 1.11 and 1.17 at ω ≈ 0.3; the sensitivity
  of V_G at the move to ω agrees (median difference −0.01). Mean fitness before the move agrees to
  0.002, except one arm (0.738 vs 0.755).
- **Response and mean fitness.** In the kernel under the move, realised response / (V_G·λ/S) is
  0.994–1.004 and mean fitness / Gaussian formula 0.999–1.000 at both arms of all 12 inputs, also
  where ω² is only 3–10× V_P.
- **V_G under the move (level).** At L ≥ 500 the screen's V_G is right for the first 100–150
  generations (0.99–1.16 of the kernel's) and then too high: 1.18–1.63 at 250–350 and 1.20–1.80
  at 900–1100, in all eight inputs. At L = 100 it is within about ±25% and too low early in two
  inputs; the screen's k50 is 11–26% too low there.

## Suggestive, not established

- **The link carrying the excess sensitivity is V_G under the move, by elimination.** The state at
  the move, the response, the mean fitness and (apparently) the extinction link are right, so what
  is left is the course of V_G during the move, where the screen is known to be wrong in level
  (the decline of `validate/genic_decline.md`).
- The direct test of that did **not** support it as I had stated it: at one common rate the ratio
  screen/kernel V_G is the same at ω− and ω+ in two of the five excess inputs, rises in two and
  falls in one, and the three inputs without an excess show the same level error. The design was
  weak for this question: a common rate puts the two arms at different distances from their own
  thresholds, and near the threshold V_G and the lag are not stationary within 2000 generations.
- The emulator may keep a negative e_ω up to a larger ω than the table (phase B's bin medians of
  −0.34 to −0.45 at ω of 2.2–5 against the table's middle tercile of −0.10 to −0.01). At the
  reference inputs this is weaker (−0.36 to −0.08 at 2.8–3.6, positive at 5). Seen post hoc.

## Untested

- L = 5000 with the kernel (the largest remaining table–emulator gap on the negative branch).
- The block-lineage model over the same ±δ (not run: budget).
- A trajectory comparison at equal relative distance from each arm's own k50, which is what would
  isolate the ω-dependence of the error in V_G under the move.
- The extinction link beyond the observation that deaths are a step in the rate.
- The published protocol's 10,000 generations from zero: the kernel runs used an equilibrium start
  with 1000 generations of burn-in. The kernel's V_G at the move is within 1.3 SD of the table's
  own row in all 12 inputs.
- ω below 0.3 with the kernel (the positive branch rests on the table and the emulator only).

## Proposals (no change made here)

1. **Compare like with like in the elasticity study.** For ω, replace the low tercile cell by two
   cells split at ω = 0.3 (or fit the table on 0.3 ≤ ω ≤ cut), or compare the emulator's
   tercile-style slope (the table's rows and weights applied to the emulator's k50(ω)) with the
   table's slope. Either removes the summary mismatch. The same check is worth making for any
   input whose elasticity changes sign within a cell.
2. **The residual at L ≥ 500 is the open problem already on record**: the kernel's decline of V_G
   under a move, which the spectrum closure lacks. This study adds no new term; it shows that the
   state at the move and the Gaussian links are not where the ω-sensitivity goes wrong.
3. **For finite differences over a wide interval use the grid solver**, or check the beta solver's
   fallback generation at both ends: two of six inputs differed by 0.07–0.10.

## Disclosures

- The main kernel grid was changed after the pilot and before any main run (8 rates × 30 runs at
  steps of 0.08 became 21 rates × 10 at steps of 0.024), because the pilot showed survival to be a
  step in the rate.
- The grid closures were run on a subset of inputs (six with dynamic LD, four without), decided
  before any was run, for cost.
- The reference-input curve compares one parameter set with table fits over all inputs.
- The kernel mode ran across hours of machine sleep; its CPU figure (2.6 h) may include wake-up
  overhead.
