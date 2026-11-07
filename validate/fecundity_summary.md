# Does vignette 11 depend on viability selection? A fecundity-selection rerun

2026-11-07. QTLSim `Params(Selection = :fecundity)` against QTSim's viability selection
on the same published inputs and seeds.

## What was run

- **Inputs:** `validate/demography_sweep.jl 2 0.01 2000 1000 fecundity`: scale 0.01 of
  each stratum with NoQTL ≤ 1000, so **900 published inputs**, 150 in each of six strata
  (NoQTL 100, 500, 1000 × u 3, 5). Strata count equally. 2000 generations after the move,
  as in vignette 11. The 5000-locus strata are left out for cost.
- **Viability** runs are the cached runs of the 1200-input sweep
  (`demography_sweep_summary.md`), restricted to these 900 inputs: the same code and seeds
  (`2026110700000 + row`). The default path is bit-identical after this change.
- **Fecundity** runs: the same inputs and seeds, at a fixed size and at R = 2. The two
  modes draw random numbers differently, so the realisations differ from generation 0.
  The pairing is by input and seed, but there are no common random numbers within a run.
- **Tables:** `validate/fecundity_compare.jl 2.0 2000`. Intervals are 95% percentile
  intervals from 2000 bootstrap resamples of inputs within strata.
- **Cost:** 23 + 19 minutes on 3 threads.

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

### Fixed size (QTSim's demography)

| | viability | fecundity |
|:--|--:|--:|
| alive 2000 generations after the move | 7.2% (5.7–8.9) | 49.1% (45.7–52.3)* |
| median generations to "extinction" (those that died) | 124 (115–133) | 975 (875–1045)* |
| median generations to 50% mortality (W̄ ≈ 0.5 for fecundity) | 45 (41–48) | 45 (42–49) |
| median standardised lag | 0.349 (0.322–0.387) | 0.299 (0.268–0.339) |
| generations until N·W̄ < 2 (median, fecundity) | – | 128 (120–136), reached by 86% |

\* At a fixed size fecundity selection is soft. Nobody dies of low fitness, and the
population is refilled to N from its fittest members. It "dies" only when fitness
underflows to zero in floating point, roughly when the lag exceeds 38 ω. That is a
numerical event, not a biological one, and **the fixed-size fecundity model has no
extinction at all**. So its survival and extinction-time rows, and the paired ratio of
12.0 (11.9–12.2, 458 inputs), measure underflow, not persistence.

The comparable fixed-size quantity is the time until the expected number of survivors,
N·W̄, falls below two: when a viability population of that fitness would most likely have
died. Paired by input, that time divided by viability's time to extinction is
**1.11 (1.09–1.12, 774 inputs)**. The population's fitness collapses on nearly the same
schedule under either mode; fecundity is about 10% slower, since its lag is a little
smaller.

### Beverton–Holt, R = 2 (vignette 11's middle growth rate)

| | viability | fecundity |
|:--|--:|--:|
| died before the move (of 900) | 4 | 4 |
| alive 2000 generations after the move | 5.2% (3.8–6.8) | 5.3% (3.9–6.9) |
| median generations to extinction (those that died) | 55 (52–58) | 58 (54–61) |
| median generations to 50% mortality | 48 (45–51) | 44 (41–48) |
| median standardised lag | 0.444 (0.384–0.482) | 0.405 (0.361–0.438) |

Paired by input:

- Time to extinction, fecundity / viability: **1.05 (1.04–1.05)**, over 848 inputs where
  both died after the move.
- Discordant survival: 0 inputs survive under viability but die under fecundity, and 1
  the other way (McNemar exact p = 1).
- Against QTSim's fixed-size viability runs of the same inputs, the median paired ratio
  of times to extinction is **0.43** for viability at R = 2 (vignette 11 reported 0.44 on
  its own sample; the 1200-input sweep 0.43) and **0.45** for fecundity at R = 2 (833
  inputs where all three died).

## Answer

**Vignette 11's conclusion does not depend on viability selection, provided the model is
density-regulated.** At R = 2, fecundity and viability selection give the same share
alive (5.3% against 5.2%, 95% intervals overlapping almost entirely, 1 discordant input
in 900). Fecundity populations take 5% longer to die (paired ratio 1.05, CI 1.04–1.05).
So density regulation shortens persistence to 0.45 of QTSim's fixed-size times under
fecundity, against 0.43 under viability. The direction and size of vignette 11's effect
stand; fecundity selection makes it about 2 percentage points weaker.

What does depend on the choice is QTSim's **fixed-size** model itself. Its extinctions
exist only because selection is on viability: with a fixed size, fecundity selection can
never empty the population. The paper's extinction rates are therefore a property of hard
(viability) selection with refilling, and vignette 11's point, that a real population
also shrinks when W̄ is low, is what makes extinction meaningful under either form of
selection. Measured by when N·W̄ first drops below two, the fixed-size fitness collapse is
again on nearly the same schedule under the two modes (ratio 1.11).

## Caveats

- 900 inputs; the 5000-locus strata are not included. One growth rate (R = 2) only.
- The intervals cover the sampling of inputs and the runs' own randomness. They do not
  cover model choices: the Beverton–Holt form, and the mean-matching definition above.
- The two modes use independent random streams. A paired ratio's precision comes from
  the shared input, not from shared noise.
- The kernel's fecundity selection agrees with the generated SLiM model in
  `validate/slim/compare_selection.jl` (see `valid_selection.log`).
