# QTSimFast's fast modes against the exact program: pre-specified criteria

Written and committed 2026-11-20 before any comparison run. Script
`validate/qtsimfast_fast.jl` (written after this file), log `validate/qtsimfast_fast.log`.

## What is tested

`QTSIM_FAST=1` and `QTSIM_FAST=2` (`reference/QTSimFast`, `FAST.md`) give up bit parity
with QTSim:

- **Level 1**: genotypes are sums of cached block sums instead of locus-by-locus sums. The
  draws are QTSim's, in QTSim's order; only the rounding of the genotypes differs (about
  1e-16 relative). Runs can diverge only where a rounding difference changes a comparison
  (a survival draw U < W, a comparison of phenotypes), which should be rare.
- **Level 2**: level 1, plus another generator (xoshiro256++) and other sampling methods
  (ziggurat normals, inversion Poisson, Lemire integers). Runs do not match QTSim's run for
  run; the claim is equality in distribution.

**Reference.** QTSimFast in its exact mode (`QTSIM_FAST` unset), which is byte-identical to
the original QTSim in all 13,915 comparisons of the parity matrix and much faster. It stands
in for QTSim here; no conclusion depends on the original binary itself.

## Scenarios

Every scenario: PopSize as stated, diploid unless stated, `Sample 1`, the other options as
listed. Fresh seeds: exact runs use seeds 710,000,001 + i, level 2 runs 720,000,001 + i
(i = 0 .. R−1; disjoint sets, never used before in this project). Level 1 uses the exact runs'
seeds (it shares their draws). R = 1000 per program, except D (R = 500, cost).

| | design | generations | move |
|:--|:--|:--|:--|
| A | single-locus protocol: `examples/single_locus.txt` (NoQTL 1, PhenotypeVar 0.01, h² 0.5, m 0.001) | 6500 | phenotype jump of 2 SD at 5000, width 0.17 from 5000 |
| B | 100 loci, Poisson(3) crossovers, PhenotypeVar 0.02, h² 0.6, m 0.001, ω 1.0, μ 1e-4 per gamete | 4000 | none |
| C | as B with 1000 loci, Poisson(5) | 4000 | none |
| D | as B with 5000 loci, Poisson(3) | 2500 | none |
| E | haploid: `examples/haploid_0257_loci.txt` (N 200, 257 loci, μ 0.05, Poisson(12)) | 600 | fitness move 0.3 at 300 |
| F | as B with 500 loci, Poisson(5) | 3000 | continuous 0.05 per generation from 2000 (all die) |
| G | as B (100 loci, Poisson(3)) | 3000 | continuous 0.009 per generation from 2000 (near the survival threshold) |

**Pilot (disclosed).** The rate of G was chosen by a pilot with the exact program only
(seeds 990,001–990,024 and 990,101–990,124, not reused): alive at generation 3000 in 24, 24,
22, 20, 9, 1, 1 and 0 of 24 runs at rates 0.004, 0.006, 0.008, 0.0085, 0.009, 0.0095, 0.010
and 0.013. Nothing else was run before this file.

## Statistics (per run, from the output lines; generation numbers as printed)

"Window mean" = the mean of a column over the printed lines of the window; a run that dies
before the end of a window has no value for it. V_G = column "Genotype Var".

- A: (a1) window mean of V_G over 4000–4999; (a2) window mean of "P Survived" over
  4000–4999; (a3) "Phenotype Mean" at 5049, runs alive then; (a4) alive at the end
  (proportion); (a5) the last printed generation of runs that died (time of death).
- B, C, D (window W = the last 1000 generations of B and C, the last 500 of D): (1) window
  mean of V_G; (2) V_G at the last generation; (3) window mean of "Phenotype Mean"; (4)
  window mean of "P Survived"; (5) window mean of "Surviving Slope".
- E: (e1) window mean of V_G over 200–299; (e2) window mean of "P Survived" over 200–299;
  (e3) alive at the end (proportion); (e4) "Phenotype Mean" at the last generation, runs
  alive.
- F: (f1) window mean of V_G over 1500–1999; (f2) time of death (last printed generation);
  (f3) alive at the end (proportion).
- G: (g1) window mean of V_G over 1500–1999; (g2) alive at the end (proportion); (g3) time
  of death of runs that died; (g4) the lag at 2199, "Fitness Function Mean" − "Phenotype
  Mean", runs alive then.

31 primary statistics in all.

## Tests and verdicts

- **Level 2, primary.** Each continuous statistic: two-sample Kolmogorov–Smirnov test of the
  exact runs' values against level 2's (asymptotic p-value; exact runs and level 2 runs are
  independent). Each proportion: Fisher's exact test (two-sided). A statistic with fewer than
  20 values in either sample is not tested and is reported as such (it is not a pass).
  Holm's procedure over all primary tests actually performed, family-wise α = 0.05.
  **Level 2 PASSES** if no test is rejected by Holm; FAILS otherwise, naming the statistics.
- **Level 2, secondary (reported, no verdict).** The k-sample Anderson–Darling test for
  each continuous statistic; the difference of means with a 95% Welch interval, relative to
  the exact runs' mean; the KS distance D; the empirical CDFs, plotted.
- **Power (stated in advance).** With 1000 runs per sample the KS test at α = 0.05 rejects
  for D ≥ about 0.061; at Holm's smallest level (0.05/31) for D ≥ about 0.084; for D with
  500 per sample about 0.12. Differences in distribution smaller than that are not excluded.
- **Level 1.** Over the exact runs' seeds, the share of runs whose whole output is
  byte-identical to the exact run's, per scenario, and for the others the first generation
  that differs. **Level 1 PASSES** if at least 95% of runs are byte-identical in every
  scenario. If it does not, its runs are also compared with the exact runs by the level-2
  tests above (conservative, since the two share their draws), and the result is reported.
- Exit status 1 if level 2 or level 1 fails.

Timings are not part of the verdict; they are measured separately on a quiet machine.

## Addendum 1 (2026-11-20, written and committed before run 2)

Run 1 (`validate/qtsimfast_fast.log`) passed both levels. Level 2 has since been changed
(commit ed6ca28): every individual's draws in a generation now come from a stream of its
own (survival from a key (Seed, generation, 1, individual); offspring i from a xoshiro256++
stream seeded from (Seed, generation, 2, i)), the draws are made on the pool's threads, and
the output line's statistics are computed in one pass. So level 2 is compared again:

- **Run 2**: the same scenarios, statistics, tests, thresholds and verdict rules as above,
  on fresh seeds: exact runs 730,000,001 + i, level 2 runs 740,000,001 + i (never used
  before). Level 1 is checked again on the exact runs' seeds (it is unchanged; its quick
  parity matrix passed after the change). Every level-2 run uses `QTSIM_THREADS=4`.
- **New criterion T (thread invariance).** For the first 100 runs of every scenario, level 2
  is also run with `QTSIM_THREADS=1` on the same seed; **T PASSES** if all 650 pairs of
  outputs are byte-identical.
- Log `validate/qtsimfast_fast_run2.log`, data `validate/qtsimfast_fast_runs_run2.csv`,
  plots `validate/qtsimfast_fast_figs_run2/`. Run 1's files are kept. The verdict on the
  current level 2 is run 2's. Exit status 1 if level 1, level 2 or T fails.
