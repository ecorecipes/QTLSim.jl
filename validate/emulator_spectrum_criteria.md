# The emulator of the published sweep with the `:spectrum` closure: pre-specified criteria

Written and committed **before** the validation script was written or run (2026-11-13). The
script `validate/emulator_spectrum.jl` applies these criteria and writes
`validate/emulator_spectrum.log`. Any later change to a criterion, endpoint, closure, sample,
seed or setting is recorded below under "Amendments", dated, with whether it was made before
or after any result was seen, and the original verdict is still reported.

This is a validation, not a model change: `emulate` and the `:spectrum` closure are used as
they are on branch `emulator-spectrum` (from main 0c74137); nothing is tuned.

## Disclosure

- **`:spectrum` was chosen after seeing it outperform `:rare_alleles` on earlier samples**:
  on the 1000-input sample of `validate/emulator.jl` (seed 20261113; x50 within 0.96–1.12 of
  the table in 7 of 8 strata, no verdict) and on the 192 E6 inputs of
  `validate/emulator_variance.jl` (mean |f − P2000| 0.14 against 0.33, no verdict). Every
  input whose `:spectrum` result has been looked at is **excluded** from both samples below.
- The `:spectrum` closure's design (CLAUDE.md, `validate/vg_move.jl`) was partly made after
  studying vignette 10's scenarios; it was not designed on the published table.
- The published table motivated this test (E2 and E4 failed with `:rare_alleles`, twice), so
  E1–E4 are evaluated on a table that is not independent of the choice of closure, even on a
  fresh sample. **E6′ is the independent test**: fresh kernel simulations on inputs never run
  before.
- E6′'s window of rates (below) is placed lower, relative to the `:rare_alleles` k50, than
  E6's, because E4 and E6 showed that `:rare_alleles`' k50 is too high (seen results of the
  baseline closure, not of `:spectrum` on these inputs).
- Sizes were set from a timing pilot run before this file was written, on inputs that are
  excluded anyway (the E6 inputs): `:spectrum` took 4.4 thread-seconds per input
  (deterministic) and 7 with R = 8, a kernel run 1.3 / 2.1 / 2.5 / 4.6 thread-seconds at
  100 / 500 / 1000 / 5000 loci (8 threads). Only timings were printed.

## Model

- **Primary:** `emulate(p; loci = :spectrum, burnin = 2000, horizon = 2000, tail = 240_000)`,
  **deterministic** (`replicates = 0`). Justification: E6 found the kernel itself nearly
  all-or-nothing per input (138 / 37 / 17 of 192 inputs with none / all / some of 40 runs
  alive), so the emulator's near-binary P2000 is not by itself the problem; what E6 left open
  is whether the threshold (k50) is in the right place. Between-run noise did not move it
  (E6 slope 0.14 vs 0.11), and for `:spectrum` the noise model has the mean noise only (no
  haplotype moments, so no V_G noise).
- **Secondary (no verdicts):** `:spectrum` with the mean noise, `replicates = 8`, seed 20261116
  (per input by `QTLSim.job_seed`, as `emulate_sweep` does), on the E6′ inputs only.
- **Baseline (no verdicts):** the deterministic `:rare_alleles` emulator (the first run's model)
  on the same inputs, with the same statistics; its numbers are printed beside the primary's
  and marked "(as a baseline)".

## Exclusions

Excluded from every sample below: every input id that appears in
`validate/emulator_spectrum.csv` (the 1000-input sample, seed 20261113),
`validate/emulator_variance_spectrum.csv` (`:spectrum` on the E6 inputs),
`validate/emulator_variance_e6.csv` and `validate/emulator_variance_kernel.csv` (the 192 E6
inputs). The script checks that this union contains the 1000 inputs of
`sort(randperm(MersenneTwister(20261113), 119984)[1:1000])` and 192 E6 inputs.

## Sample for E1–E4: 2400 fresh inputs

Inputs in the order of `validate/emulator.jl` (stratum, then file name). One
`MersenneTwister(20261116)` is used in turn for the 8 strata in that order: in each stratum the
eligible inputs (not excluded), in input order, are permuted with `randperm(rng, n_eligible)`
and the first **300** taken; the sample is sorted. The strata are of equal size (14,998), so
this is an equal-probability sample of the non-excluded inputs. Matching and table conventions
exactly as in `validate/emulator_criteria.md` with its amendments 1, 2 and 4 (alive at 2000 =
"Generations to extinction" > 2000; T counted from the first generation of the move as 1;
deciles equal-count with ties broken by input order).

## Criteria on the published table (thresholds of `validate/emulator_criteria.md`)

- **E1 discrimination.** AUC of P2000 for alive-at-2000 ≥ 0.85 overall and ≥ 0.80 in each of
  the 8 strata.
- **E2 calibration.** In 10 equal-count deciles of P2000 the observed fraction alive is within
  0.10 of the mean P2000 in every decile; and the calibration slope (logistic regression of the
  outcome on logit(P2000), P2000 clipped to [1e-6, 1 − 1e-6]) lies in [0.8, 1.25].
- **E3 timing.** Among inputs extinct within 2000 generations in the table with a predicted T:
  Spearman ≥ 0.7 and median T_pred/T_obs in [0.8, 1.25].
- **E4 half-survival rate.** Per stratum, x50 (logistic regression on log(MoveFitnessSD/√V_P);
  fractional response for the emulator): emulator/table ratio in [0.8, 1.25] in every stratum.
- **E5 (exploratory, no verdict).** AUC of P240k for survived-240,000.
- Reported with E4, without verdicts: a bootstrap 95% interval of each stratum's ratio (1000
  resamples of inputs within the stratum, `MersenneTwister(1)`), and the table's x50 on all
  14,998 inputs of the stratum (a less noisy reference).

## E6′: a new independent test against fresh kernel simulations

**Inputs (160).** (Candidate source changed by amendment 1, before any result.) Chosen without looking at any table outcome, from the 2000 inputs of the
first run's parameter map (`validate/emulator_map.csv`, uniform sample of seed 20261114, which
holds each input's deterministic `:rare_alleles` k50), excluding the excluded inputs above and
the 2400 inputs of the E1–E4 sample. Let r = MoveFitnessSD/k50. Candidates: finite positive
k50 and **0.3 ≤ r ≤ 1.5**. In each stratum the candidates (sorted by input id) are split into 4
bins of equal width in log r over [0.3, 1.5]; **5** are drawn at random from each bin; a bin
with fewer than 5 contributes all of them and the stratum's remaining slots are filled at random
from its other undrawn candidates. One `MersenneTwister(20261117)` is used in turn for the
strata (in input order) and their bins (increasing r). **20 per stratum, 160 inputs.**

**Kernel runs.** Each input **40** times, exactly as E6: `simulate(Params(p; NoGenerations =
12000, Sample = 1, Seed = nothing); seed, init = :equilibrium, burnin = 1000)`; fresh seeds
`30261116000 + 100·j + r` for input j (1-based, in selection order) and replicate r = 1…40.
Alive at 2000 = not extinct by the end of the run (generation 11999). An extinction before the
move counts as not alive.

**Statistics.** Per input: the kernel's fraction alive f_i (of 40); the emulator's P2000_i;
rate_i = MoveFitnessSD/√V_P.

- **E6′a binned calibration** (as E6a): bins of P2000 [0, 0.1), [0.1, 0.3), [0.3, 0.5),
  [0.5, 0.7), [0.7, 0.9), [0.9, 1]; in every bin with ≥ 10 inputs |mean f − mean P2000| ≤ 0.10.
- **E6′b calibration slope** (as E6b): logistic regression of the 40 × 160 binary outcomes on
  logit(P2000) (clipped to [1e-6, 1 − 1e-6]); slope in [0.8, 1.25].
- **E6′c threshold test.** Kernel majority alive M_i = (f_i > 0.5); emulator call
  C_i = (P2000_i > 0.5). PASS iff the accuracy, the fraction of inputs with C_i = M_i, is
  ≥ 0.85 **and** the AUC of P2000_i for M_i (midranks) is ≥ 0.85. Reported: the 2 × 2 table,
  balanced accuracy, and the same by stratum (no verdict).
- **E6′d half-survival rate by stratum.** x50 of the kernel: logistic regression of the
  40 × 20 binary runs of the stratum on log rate; x50 of the emulator: the same regression
  with P2000 as a fractional response over the stratum's 20 inputs; x50 = exp(−β₀/β₁). If a
  regression does not converge (complete separation) and the responses ≥ 0.5 all have lower
  rates than the responses < 0.5 (for the kernel, f_i ≥ 0.5 vs < 0.5), x50 is the geometric
  mean of the largest rate with response ≥ 0.5 and the smallest with response < 0.5. x50 is
  **undefined** (the stratum fails) if the fitted slope β₁ ≥ 0, or if a regression does not
  converge and no such separation exists, or if all responses are on one side of 0.5. PASS iff
  emulator/kernel lies in [0.8, 1.25] in every stratum. Reported without verdict: the pooled
  ratio over all 160 inputs, and bootstrap 95% intervals (1000 resamples of inputs within the
  stratum, `MersenneTwister(2)`).
- **E6′ PASS** iff E6′a, E6′b, E6′c and E6′d all pass; each is also reported separately.
- Reported without verdicts: mean |f_i − P2000_i|, the Brier score, the binomial
  log-likelihood Σ log Bin(s_i; 40, P2000_i) (clipped), the intercept of E6′b, the AUC over
  runs, and the counts of inputs with none / all / some runs alive.

The baseline (`:rare_alleles`, deterministic) and the secondary (`:spectrum`, R = 8) get the
same statistics on the E6′ inputs, without verdicts.

## Budget

About 10 CPU-hours, on 8 threads: 2400 × 4.4 s (`:spectrum`) ≈ 2.9 h, 6400 kernel runs ×
~2.6 s ≈ 4.6 h, the secondaries and baseline ≈ 0.6 h.

## Amendments

1. **2026-11-13, after writing the script and running only its `select` mode (sample sizes and
   counts per bin, no emulation, no kernel run, no table outcome), BEFORE any result: the
   source of E6′'s candidates.** After the exclusions the parameter map has only 101 inputs
   with 0.3 ≤ r ≤ 1.5 (8, 8, 9, 9, 14, 13, 20, 20 per stratum; most near-k50 inputs of the map
   were taken by E6), so 160 cannot be drawn from it. Instead, the candidates come from a new
   pool, with the same k50 (deterministic `:rare_alleles`, `emulated_half_survival_rate` with
   its defaults, as the map) and the same rule afterwards:
   - In each stratum, the eligible inputs (not excluded, not in the E1–E4 sample), in input
     order, are permuted with one `MersenneTwister(20261117)` used in turn for the strata;
     walking the permutation, an input enters the pool if its MoveFitnessSD/kc lies in
     [0.3, 3.0], kc = `drift_critical_rate(p).k` (a prefilter; on the map 47% of such inputs
     have r in [0.3, 1.5], against 14% of all); the first **150** per stratum enter the pool.
   - k50 is computed for the 1200 pool inputs (about 2 s each; cache
     `validate/emulator_spectrum_pool.csv`); candidates are pool inputs with finite positive
     k50 and 0.3 ≤ r ≤ 1.5; then the bins, 5 per bin, the filling rule, 20 per stratum and
     160 in all, as above, with a second generator `MersenneTwister(20261118)` used in turn for
     the strata and bins (candidates sorted by input id).
   The map is no longer used (about 0.7 CPU-hours added to the budget). The pool's k50 is a
   selection quantity of the baseline closure, computed in the script's `select` mode before
   the run; it involves no table outcome and no `:spectrum` emulation. Kernel runs, seeds,
   statistics and thresholds are unchanged.
