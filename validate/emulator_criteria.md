# The emulator of the published sweep: pre-specified criteria (roadmap F4)

Written and committed **before** the emulator was implemented or run (2026-11-13). The
script `validate/emulator.jl` applies these criteria and writes `validate/emulator.log`;
any later change to a criterion, endpoint, closure, horizon or sample is recorded below
under "Amendments", dated, with the reason, and the original verdict is still reported.

## What is emulated

For each of the 119,984 published inputs (`Repository/Simulations/`, matched to rows of
`Repository/Results/Master Output file.xlsx` by parameters, as `validate/published_sweep.jl`
does), the moments lowering run along the published protocol's continuous move at the
input's own rate, with drift (`erosion = true`), started at the equilibrium genetic
variance with drift (`house_of_gauss_vg`, then a burn-in at a fixed optimum as in
`half_survival_rate`), and `extinction_distribution` along that trajectory (fixed size, as
QTSim) give for each input:

- P2000: the probability of being alive 2000 generations after the move starts;
- T: the median generation of extinction after the move (when P(extinct) reaches 0.5
  within the horizon);
- P240k: the probability of surviving 240,000 generations, extrapolated from P2000 with
  the per-generation hazard at the end of the horizon held constant (exploratory).

Primary closure: `loci = :rare_alleles` (the drift-regime closure; table/prediction of
the pre-move genetic variance 0.991 on 320 published inputs). Secondary: `loci = :spectrum`
on a uniform random sample of 1000 inputs (seed 20261113; it is about 4 s per input), and
`:multilocus` on the same 1000 for comparison. Burn-in and horizon are fixed at 2000 and
2000 generations; they are not tuned.

Observed outcomes, from the table: alive at 2000 = "Generations to extinction" > 2000;
extinction time = "Generations to extinction" (generation of death − 9999); survived =
"Extinct (0,1)" = 1.

## Criteria (primary closure, all inputs)

- **E1 discrimination.** The AUC of P2000 for the observed alive-at-2000 is ≥ 0.85 over
  all inputs and ≥ 0.80 in each of the 8 strata (loci × mean crossover number).
- **E2 calibration.** In 10 deciles of P2000 (all inputs), the observed fraction alive is
  within 0.10 of the mean P2000 in every decile; and the calibration slope (logistic
  regression of the observed outcome on logit(P2000), P2000 clipped to [1e-6, 1 − 1e-6]) lies
  in [0.8, 1.25].
- **E3 timing.** Among inputs extinct within 2000 generations in the table and with a
  predicted T, the Spearman correlation of T with the observed extinction time is ≥ 0.7,
  and the median of predicted/observed lies in [0.8, 1.25].
- **E4 the half-survival rate.** In each stratum, x50 — the rate (MoveFitnessSD/√V_P) at
  which half are alive at 2000, from a logistic regression of the outcome on log rate — is
  computed from the table's outcomes and from the emulator's P2000 (the same regression,
  outcome replaced by P2000 as a fractional response); their ratio lies in [0.8, 1.25] in
  every stratum.
- **E5 (exploratory, no verdict).** The AUC of P240k for survived-240,000.

Each criterion is reported PASS or FAIL with its numbers; the secondary closures are
reported with the same numbers on their sample, without verdicts.

## After the criteria: the parameter map (exploratory)

On a uniform random sample of 2000 inputs (seed 20261114), the emulated half-survival rate
k50 (bisection on MoveFitnessSD for P2000 = 0.5) and its elasticities with respect to the
continuous parameters (PhenotypeVar, Heritability, MutationScalar, MutationRatePerGamete,
the fitness width InitFitnessSD) by central differences with Richardson extrapolation
(`QTLSim.Para`, `sensitivities`); NoQTL and RecombPoissonMean are compared across strata.
Reported as medians and interquartile ranges, by stratum. Exploratory: no verdicts.

## Amendments

All recorded on 2026-11-13, after the code was written and timed on a few inputs and BEFORE
any emulation was compared with the table (no criterion, endpoint, closure, horizon, burn-in
or sample is changed; these settle points the criteria leave open):

1. **Counting T.** The emulator's median generation of extinction is counted as the table
   counts it: the first generation of the move is 1 (`T = median − burnin + 1`), so "alive
   at 2000" means the same event (death after the 2000th generation of the move) in both.
2. **Deciles (E2).** The 10 deciles are equal-count groups of the inputs sorted by P2000;
   ties (many inputs have P2000 exactly 0 or 1) are broken by input order (stratum, then
   file name).
3. **Finding k50 (parameter map).** "Bisection on MoveFitnessSD" is implemented as a
   bracketing search in log k (Illinois regula falsi on log(log 2) − log(−log P2000),
   safeguarded by bisection) to 1e-10 in log k (1e-9 for the perturbed parameters of the
   finite differences, from a bracket of ±0.2% around the base k50), because 12 bisection
   steps resolve log k only to about 1e-3, too coarse for central differences with
   rel = 1e-3. The starting bracket is `half_survival_rate`'s, kc/8 to 2kc
   (`drift_critical_rate`), extended by doubling up to 64kc and halving down to kc/512
   when it does not contain the root: timing runs on a few inputs found P2000 ≥ 0.5 at
   2kc for some of them. The root found is the same one bisection would converge to.
4. **Matching.** Rows are matched with `match_row`'s conditions through an index on
   Heritability (`row_index`, validate/published_table.jl), checked against `match_row`
   on every 600th input; all 119,984 inputs match a unique row.
