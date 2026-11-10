# Does vignette 11 depend on viability selection? A fecundity-selection rerun

QTLSim `Params(Selection = :fecundity)` against QTSim's viability selection on the same
published inputs and seeds. First run 2026-11-07 at R = 2; extended 2026-11-10 to
R = 5 and 1.5 and to a 5000-locus stratum. The full tables are in
`validate/fecundity_compare.log` (`validate/fecundity_compare.jl 5,2,1.5 2000`).

## What was run

- **Inputs, NoQTL ≤ 1000:** `validate/demography_sweep.jl 2,5,1.5 0.01 2000 1000 fecundity`,
  i.e. scale 0.01 of each stratum with NoQTL ≤ 1000: **900 published inputs**, 150 in each
  of six strata (NoQTL 100, 500, 1000 × u 3, 5). Fixed size and R = 5, 2 and 1.5.
- **Inputs, NoQTL = 5000:** `validate/demography_sweep.jl 2 0.00334 2000 5000 fecundity`,
  a smaller sample for cost: **100 inputs**, 50 per u, the first 50 of each 5000-locus
  stratum's 150 in the 1200-input sweep (the samples are nested). Fixed size and R = 2
  only.
- 2000 generations after the move, as in vignette 11. Strata count equally.
- **Viability** runs are the cached runs of the 1200-input sweep
  (`demography_sweep_summary.md`), restricted to these inputs: the same code and seeds
  (`2026110700000 + row`). Reruns of 8 randomly chosen cached rows per growth rate (500 and
  1000 loci) with the current code reproduced the cached R = 5, 2 and 1.5 rows bit for
  bit. The **fixed-size viability** cache was written by the dense kernel, before the
  block kernel: in 8 checked rows every discrete outcome (time to extinction, 50%
  mortality, dying generation, sizes) is identical and 2 of the 8 differ in the last digit
  of the continuous measures (about 1e-16 relative, the documented rounding difference
  of block genotype sums). It was reused, not rerun, for cost.
- **Fecundity** runs: the same inputs and seeds. The two modes draw random numbers
  differently, so the realisations differ from generation 0. The pairing is by input and
  seed, but there are no common random numbers within a run. Rerunning the R = 2 set on
  2026-11-09 reproduced the 2026-11-07 tables exactly.
- **Intervals:** 95% percentile intervals from 2000 bootstrap resamples of inputs within
  strata.
- **Cost:** the 900-input fixed-size set and R = 2 took 23 + 19 minutes on 3 threads
  (2026-11-07). The R = 1.5 set took 6.3 minutes on 8 threads. The other sets ran on a
  machine that was asleep or heavily loaded for most of their wall time, so their wall
  times are not informative; by CPU time each 900-input growth rate is about 1 CPU-hour,
  the fixed size about 1.5, and the 5000-locus sets about 0.5.

### The definition being compared

Under fecundity selection everyone survives, and each parent is drawn with probability
`W/ΣW` (two distinct parents). Under Beverton–Holt the offspring number is
`Poisson(R S/(1 + (R−1) S/K))` with `S = ΣW = N·W̄`. That is viability's formula with the
number of survivors replaced by its expectation. So the two modes have the same expected
selection differential (exactly at a fixed size, to O(1/N) for viability's ratio
estimator) and the same expected population growth. Viability's growth is lower by a
Jensen gap of at most `R(R−1)N/(4K)`, about 0.5 offspring at R = 2. What differs is the
variance: who survives is random under viability, and viability's parents are fewer
(`N W̄` against `N W̄²/E[W²]`), so its drift is stronger.

## Results

### Beverton–Holt, by growth rate

Time ratios are medians over inputs where the runs compared died after the move.
"Fixed" is QTSim's fixed-size (viability) run of the same input.

| group | R | inputs | alive, viability | alive, fecundity | discordant (V only / F only) | T fecundity / viability | T viability / fixed | T fecundity / fixed |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| NoQTL ≤ 1000 | 5 | 900 | 5.9% (4.3–7.4) | 6.8% (5.2–8.4) | 0 / 8 | 1.06 (1.06–1.07) | 0.60 (0.59–0.60) | 0.64 (0.63–0.64) |
| NoQTL ≤ 1000 | 2 | 900 | 5.2% (3.8–6.8) | 5.3% (3.9–6.9) | 0 / 1 | 1.05 (1.04–1.05) | 0.43 (0.43–0.44) | 0.45 (0.45–0.46) |
| NoQTL ≤ 1000 | 1.5 | 900 | 4.3% (3.1–5.8) | 4.3% (3.1–5.8) | 1 / 1 | 1.03 (1.03–1.04) | 0.37 (0.36–0.37) | 0.38 (0.37–0.39) |
| NoQTL = 5000 | 2 | 100 | 8.0% (3.0–14.0) | 8.0% (3.0–14.0) | 0 / 0 | 1.04 (1.02–1.06) | 0.44 (0.43–0.46) | 0.45 (0.44–0.47) |

The inputs where the time ratios are taken: T fecundity / viability over 838, 848, 848
and 92 inputs; the ratios against fixed size over 831, 833, 825 and 88 inputs where all
three runs died after the move (rows in the order of the table).

- **R = 5:** fecundity selection keeps a few more populations alive (8 inputs alive only
  under fecundity, none the other way; McNemar exact p = 0.008). The shares alive, 6.8%
  against 5.9%, have overlapping intervals.
- **R = 2:** unchanged from the first run. 1 discordant input in 900.
- **R = 1.5:** the modes are indistinguishable: 4.3% alive under both, one discordant
  input each way.
- **5000 loci, R = 2:** the same picture as with fewer loci: 8 of 100 alive under both
  modes, no discordant input, a paired ratio of 1.04, and 0.44 against 0.45 of the
  fixed-size times.
- Populations dying before the move under density regulation (strong stabilising
  selection, mean fitness at the optimum below 1/R): 1, 4 and 12 under viability and 1, 4
  and 11 under fecundity, for R = 5, 2 and 1.5 (of 900).

### Fixed size (QTSim's demography)

| | viability, NoQTL ≤ 1000 | fecundity, NoQTL ≤ 1000 | viability, 5000 | fecundity, 5000 |
|:--|--:|--:|--:|--:|
| inputs | 900 | 900 | 100 | 100 |
| alive 2000 generations after the move | 7.2% (5.7–8.9) | 49.1% (45.7–52.3)* | 12.0% (6.0–19.0) | 48.0% (39.0–58.0)* |
| median generations to "extinction" (those that died) | 124 (115–133) | 975 (875–1045)* | 119 (92–181) | 828 (708–930)* |
| median generations to 50% mortality (W̄ ≈ 0.5 for fecundity) | 45 (41–48) | 45 (42–49) | 50 (33–80) | 53 (32–81) |
| median standardised lag | 0.349 (0.322–0.387) | 0.299 (0.268–0.339) | 0.341 (0.246–0.443) | 0.364 (0.252–0.469) |
| N·W̄ never below 2 within 2000 generations | – | 14.0% (11.8–16.2) | – | 15.0% (8.0–22.0) |
| generations until N·W̄ < 2 (median, those that got there) | – | 128 (120–136) | – | 123 (93–181) |

\* At a fixed size fecundity selection is soft. Nobody dies of low fitness, and the
population is refilled to N from its fittest members. It "dies" only when fitness
underflows to zero in floating point, roughly when the lag exceeds 38 ω. That is a
numerical event, not a biological one, and **the fixed-size fecundity model has no
extinction at all**. So its survival and extinction-time rows, and the paired ratios of
12.0 (11.9–12.2, 458 inputs) and 11.2 (10.8–11.5, 52 inputs at 5000 loci), measure
underflow, not persistence.

The comparable fixed-size quantity is the time until the expected number of survivors,
N·W̄, falls below two: when a viability population of that fitness would most likely have
died. Paired by input, that time divided by viability's time to extinction is
**1.11 (1.09–1.12, 774 inputs)**, and 1.09 (1.07–1.11, 85 inputs) at 5000 loci. The
population's fitness collapses on nearly the same schedule under either mode; fecundity
is about 10% slower, since its lag is a little smaller.

## Answer

**Vignette 11's conclusion does not depend on viability selection, provided the model is
density-regulated, at every growth rate tested and with 5000 loci.** Fecundity
populations take 3–6% longer to die than viability populations of the same input (paired
ratios 1.03–1.06, intervals within 1.02–1.07), the effect being largest at R = 5. Against
QTSim's fixed-size times, density regulation shortens persistence to 0.64, 0.45 and 0.38
under fecundity for R = 5, 2 and 1.5, against 0.60, 0.43 and 0.37 under viability. So
the direction and size of vignette 11's effect stand; fecundity selection makes it at
most 4 percentage points weaker. The share alive agrees within the intervals at every R;
only at R = 5 are the discordant inputs one-sided (8 to 0).

What does depend on the choice is QTSim's **fixed-size** model itself. Its extinctions
exist only because selection is on viability: with a fixed size, fecundity selection can
never empty the population. The paper's extinction rates are therefore a property of hard
(viability) selection with refilling, and vignette 11's point, that a real population
also shrinks when W̄ is low, is what makes extinction meaningful under either form of
selection. Measured by when N·W̄ first drops below two, the fixed-size fitness collapse is
again on nearly the same schedule under the two modes (ratio 1.11, and 1.09 at 5000 loci).

## Caveats

- The 5000-locus check is 100 inputs at R = 2 only, so its share alive has wide
  intervals (3–14%); its time ratios are as precise as those of the other strata because
  they are paired.
- The fixed-size viability runs come from the dense kernel and agree with the current
  kernel to rounding, not bit for bit (see "What was run").
- The intervals cover the sampling of inputs and the runs' own randomness. They do not
  cover model choices: the Beverton–Holt form, and the mean-matching definition above.
- The two modes use independent random streams. A paired ratio's precision comes from
  the shared input, not from shared noise.
- The kernel's fecundity selection agrees with the generated SLiM model in
  `validate/slim/compare_selection.jl` (see `valid_selection.log`).
