# Pre-specified criteria: the two-stage emulator

Written 2026-11-16, after development (validate/two_stage_dev.jl, validate/two_stage_dev.log, on
development inputs only) and before any comparison on fresh data. Committed before
validate/two_stage.jl is run on fresh inputs (its `select` mode, which only draws the samples and
computes the pool's `:rare_alleles` k50s, may run first; no outcome is looked at). Output:
validate/two_stage.log. No default changes whatever the outcome.

## The model under test

`emulate(p; loci = :two_stage)` (src/two_stage.jl) with its defaults, burn-in and horizon 2000:

1. **Screen (stage 1)** on every input: the deterministic emulator `emulate(p; loci = :spectrum,
   spectrum = (; ld = :dynamic, solver = :beta))` at the input's rate k; its call is P2000 > 1/2.
2. **Band.** The input goes to stage 2 iff k is within a factor exp(±w) of the screen's own
   half-survival rate k50_screen, w = **0.25**. Decided without finding k50_screen: one more screen
   evaluation at k·e^{+w} if the screen calls the input alive at k, at k·e^{−w} if dead; inside iff
   the two calls differ (P2000 decreases in k). A screen trajectory that dies before the move is
   outside (the burn-in does not depend on k).
3. **Stage 2** inside the band: the block model, `emulate(p; loci = :blocks, blocks = 160,
   within = :gaussian, replicates = 40, seed)` (the configuration validated in validate/blocks.log),
   replicate r from `Xoshiro(job_seed(seed, r))`.
4. **Combination.** The result is the screen's outside the band and the block model's inside
   (`stage` 1 or 2); the screen's values are recorded beside it.

## Disclosure: how the rule was set (development data only)

Development inputs (validate/two_stage_dev.log): blocks.jl's H (120 near-threshold inputs; block
model = 40 `run_blocks` runs each, cached in wt-blocks-rep/validate/blocks_raw/blocks_traj.csv) and
blocks.jl's F (320 random table inputs; block model = 10 replicates, blocks_emulate.csv). The
kernel's runs on H (blocks_kernel.csv) were printed for information and set nothing. The rule
(distance from the screen's own k50 in log k) was chosen before any distance was computed; the
width is the smallest candidate (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7) at which the screen
disagrees with the block model on < 1% of the development inputs outside the band, on H and on F
separately; the solver by the rule written in the dev log before the sweep. Results (two_stage_dev.log): with the beta solver w* = 0.25 (outside the band: 0 of 90
near-threshold and 0 of 313 table inputs disagree with the block model; at 0.2 one of 96
near-threshold inputs did); stage-2 fraction 0.25 near the threshold and 0.022 (95% CI 0.006–0.038)
of the table sample. With the grid solver w* = 0.4 (Hd 0.38, Fd 0.069). The calls at k of the
two solvers agree on all 440 inputs; their band decisions at 0.25 on 433 of 440 (0.984); the
beta is about 3.5× cheaper per evaluation (2.05 vs 7.10 s, threaded). By the rule, the screen is
the beta at w = 0.25. Disclosed: the margin is narrow; w* is set by one near-threshold input at
d = 0.246 (0.312 with the grid). On the development inputs (information only) the two-stage
accuracy against the kernel's majority would be 0.992, as the block model alone.

## Samples (fresh)

- **Exclusions.** spectrum_ld.jl's `all_exclusions` (E6′'s pools and samples, the spectrum
  diagnosis's inputs, validate/published_sweep_results.csv, the table rows rerun in vignettes 09
  and 11), united with every input id in any file name or any .csv/.log/.txt/.md/.jl file under
  `validate/` of the main package and of every `wt-*` worktree except wt-paper, **including**
  wt-blocks-rep/validate/blocks_raw/ (blocks.jl's H, F and pool) and this task's development cache
  (validate/two_stage_dev_raw/); only this validation's own caches (validate/two_stage_raw/) and log
  are skipped, and files naming more than 50,000 ids (the deterministic `:rare_alleles` emulator on
  the whole table, which no kernel comparison used). The script prints the counts by directory,
  checks the samples' set intersection with the exclusions (must be 0), and separately with every
  id in blocks.jl's caches (must be 0).
- **P (primary).** E6′'s procedure (as blocks.jl's H) on a new pool: per stratum a
  MersenneTwister(20261141) permutation of the non-excluded inputs, prefiltered by
  0.3 ≤ MoveFitnessSD/kc ≤ 3 (kc = `drift_critical_rate`), 100 per stratum, their `:rare_alleles`
  k50 (`emulated_half_survival_rate`, burn-in and horizon 2000); picks with 0.3 ≤ MoveFitnessSD/k50
  ≤ 1.5, 3 per log-bin (4 bins), filled to **15 per stratum** from the rest, MersenneTwister(20261142):
  **120 inputs**. No outcome and no screen result is used in the selection.
  - Kernel: 40 runs per input, E6′'s protocol (init = :equilibrium, burn-in 1000, to generation
    12000, the move at 10000), seeds 51261116000 + 100j + r; outcome alive at +2000; f = fraction
    alive; the kernel's call is f > 0.5.
  - Emulators, per input in one job and in this order, each timed (wall time of the call in its
    thread, 8 threads): the two-stage emulator (seed `job_seed(20261144, id)`); the block model
    alone with the same seed (`emulate(p; loci = :blocks, blocks = 160, within = :gaussian,
    replicates = 40)`, so that inside the band the two-stage result equals it: common random
    numbers, checked); `:spectrum` with `ld = :dynamic` alone (grid, deterministic); the default
    `:spectrum`.
- **S (secondary).** 250 per stratum (2000), MersenneTwister(20261143) permutation, outside the
  exclusions and P's pool. The two-stage emulator (seed `job_seed(20261145, id)`), timed; on a
  stratified subsample of 32 (every 62nd of S in input order: 5/4/4/4/4/4/4/3 per stratum) also the block model alone
  and `ld = :dynamic` alone, timed in the same job. Truth: the table's outcomes (alive at 2000 =
  "Generations to extinction" > 2000).

## Tests

**Primary (P)** — calls P2000 > 1/2 against the kernel's majority f > 0.5.

- **(T1)** two-stage vs the block model alone: PASS iff (a) accuracy: (inputs the block model gets
  right and the two-stage wrong) − (the reverse) ≤ 1, **or** exact two-sided McNemar p > 0.05 on
  those discordant inputs (the direction is reported), **and** (b) AUC(P2000; f > 0.5) of the
  two-stage ≥ that of the block model − max(1/n₁, 1/n₀) (the most one input can move the AUC; n₁,
  n₀ the kernel-alive and kernel-dead counts).
- **(T2)** x50 (E6′d's `x50_rule`: d/√V_P, logistic regression with the kernel's f as outcome
  and P2000 as a fractional response, the separation rule otherwise): pooled ratio
  two-stage/kernel in [0.8, 1.25], **and** every stratum's ratio in [0.8, 1.25] (an undefined
  stratum fails).
- **(T3)** median over inputs of the two-stage's cost per input ≤ 1/3 of the median of the block
  model alone's.
- Reported without verdict: the stage-2 fraction (overall and per stratum), E6′a–d for the
  two-stage and the three comparators, their accuracy/AUC and pooled x50, the paired discordances
  with ld = :dynamic alone and the default `:spectrum`, the mean costs.

**Secondary (S)** — table outcomes as truth; E1–E4 of validate/emulator_criteria.md with these
rules, fixed now because blocks.jl found undefined statistics:

- **(E1)** AUC(P2000; alive at 2000) ≥ 0.85 over S, and ≥ 0.80 in every stratum where it is
  defined; a stratum with only one outcome has no AUC: it is listed and excluded. E1 also needs at
  least 6 of the 8 strata defined (else FAIL).
- **(E2)** deciles of P2000 (ties by input order): |observed − mean P2000| ≤ 0.10 in every decile;
  and the calibration slope (logistic regression of the outcome on logit P2000, P2000 clipped to
  [1e-6, 1 − 1e-6]) in [0.8, 1.25] **if the fit converges**. If it does not converge (separation:
  the maximum-likelihood slope does not exist), the slope is reported as undefined and E2 is scored
  on the deciles alone, and the log says so.
- **(E3)** unchanged: among inputs extinct within 2000 in the table with a predicted T (stage 2:
  the median of the block runs' extinction times, missing if half or more survive), Spearman ≥ 0.7
  and median T_pred/T_obs in [0.8, 1.25].
- **(E4)** x50 per stratum by E6′d's rule (`x50_rule`) for the table (outcomes) and for the
  emulator (P2000 as a fractional response), instead of the plain logistic `x50`, whose value is
  meaningless when the fit does not converge (the plain values are printed beside it). A stratum
  whose table x50 is undefined (all outcomes on one side, or no convergence without separation)
  is listed and excluded; one whose table x50 is defined and the emulator's not fails. PASS iff
  every defined stratum's ratio is in [0.8, 1.25] and at least 6 strata are defined.
- Reported: the stage-2 fraction with a 95% CI and per stratum; the screen alone's E1–E4 (its
  P2000 and T at k; no verdict); total cost — the two-stage's mean time per input on S, and the
  block model alone's and ld = :dynamic alone's on the subsample of 32, each × 119,984 (CPU-hours
  at threaded speed), and the ratios.

**Optional** (only if T1–T3 all pass): the wall time of a two-stage emulation of all 119,984 inputs
on 8 threads, estimated from S; run only if it fits in the remaining budget (cached with
`emulate_sweep`, git-ignored).

## Budget

About 15 CPU-hours for the task; development used about 4 CPU-hours (screen sweeps 3.1 grid + 0.9 beta). Estimated: pool k50s 0.5 CPU-h;
P kernel 4800 runs ≈ 3.3; P emulators ≈ 1.5 (block model alone ≈ 25 s per input, the rest ≈ 15 s);
S two-stage ≈ 2000 × 4 s (two beta evaluations) plus the block model on the stage-2 inputs (≈ 25 s each) ≈ 2.5; the
cost subsample ≈ 0.3. Total ≈ 8 CPU-h.

## Amendments

(none)

## Record before the fresh run

The `select` mode was run once before this file was committed (2026-11-16 11:13; no outcome
computed): 20,226 exclusions (spectrum_ld's 6,097 ∪ 19,132 scanned ids); P = 120 (pool 800; per
stratum and bin 3/3/3/3 except 2/3/3/3 in QT0100-Possion5.0 and 1/3/3/3 in QT0500-Possion5.0,
filled to 15 per stratum), S = 2000, all intersections with the exclusions and with the ids of
blocks.jl's caches (1,120) are 0.
