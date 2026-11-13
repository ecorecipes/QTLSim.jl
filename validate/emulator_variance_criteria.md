# Between-run variance in the emulator of the published sweep: pre-specified criteria

Written and committed **before** the change was implemented or run (2026-11-13). The script
`validate/emulator_variance.jl` applies these criteria and writes
`validate/emulator_variance.log`. Any later change to a criterion, endpoint, closure, number
of replicates, horizon or sample is recorded below under "Amendments", dated, with whether it
was made before or after any result was seen, and the original verdict is still reported.

## Disclosure

This change is motivated by the failure of E2 (calibration) and E4 (half-survival rate) of
the first validation of the emulator (`validate/emulator_criteria.md`, `validate/emulator.log`)
**on the same published table**. The diagnosis (post hoc): in the emulator only demographic
chance is random; the moments trajectory (mean, lag, V_G) is deterministic, so P2000 is
nearly a step function of the rate. E1–E4 below are therefore evaluated on a table that
motivated the design: they are not independent of it, even on a fresh sample of inputs.
**E6 is the independent test**: it compares the emulator with fresh simulations of the Julia
kernel, which were not used in the design.

## The noise model

`emulate(p; replicates = R)` (R > 0; `R = 0`, the default, is the deterministic emulator,
bitwise unchanged) runs R replicate trajectories of the moments lowering, each with its own
realised mean genotype:

- **Drift of the mean.** The drift-mode moments lowering (`erosion = true`) already carries,
  for its haplotype moments, `C`: the expected square of the deviation of the population's
  realised mean haplotype from the deterministic mean, so that `Ploidy²·C` is the variance
  among replicate runs of the mean genotype around the deterministic trajectory. In each
  replicate, every generation, the offspring's mean genotype receives a Gaussian innovation
  of variance `Ploidy²·(C_new − C_sel)`, the moments' own increase of `C` over reproduction
  and mutation (`C_sel` after selection, `C_new` in the offspring after mutation); to leading
  order this is `V_G/n + V_m/N` (n = N·W̄ the parents, the same N_e that the drift of the
  moments uses). Selection then acts on each replicate's realised lag, shrinking a deviation
  by `1 − V_G/s` per generation, the factor by which the moments shrink `C`, so that the
  variance of the replicates' mean about the deterministic one follows `Ploidy²·C` (exactly in
  the linear Gaussian recursion; to first order once the deviation changes W̄ and n). Each
  replicate also carries its own V_G, mean fitness and number of parents, which differ
  through the realised lag (n = N·W̄ depends on it). No constant is fitted; nothing is taken
  from the published table.
- **Burn-in.** The burn-in at a fixed optimum (2000 generations, as before) is run once,
  deterministically, and shared; at the start of the move each replicate starts from a mean
  genotype drawn from `N(z̄_det, Ploidy²·C)`, the moments' stationary spread of the mean. In
  the linear Gaussian recursion this is the law that running the noise through the burn-in
  would give, at half the cost; the second-order effect of the lag on W̄ and on the drift of
  V_G during the burn-in (of order 1/(4n) in W̄) is ignored.
- **For closures without haplotype moments** (`:spectrum`, `:aggregate`), the innovation is
  the leading-order `V_Gs/n + V_m/N` (V_Gs the survivors' V_G) and the starting spread the
  stationary variance of the AR(1) recursion, `σ²/(β(2 − β))`, β = V_G/s.
- **Demographic chance** is added as before: `extinction_distribution` along each
  replicate's trajectory of W̄ from the start of the move. P2000 is the mean over replicates
  of their probabilities of being alive 2000 generations after the move starts; T is the
  median of the mixture of their extinction-time distributions; P240k is the mean over
  replicates of each one's constant-end-hazard extrapolation.
- **Not included:** random fluctuation of V_G itself (V_G in each replicate is the moments'
  expected within-population variance, given that replicate's history of parents; the
  sampling of allele frequencies that makes V_G of a real run wander is a fourth-moment
  quantity the lowering does not carry), non-Gaussian phenotype distributions, and the
  feedback of the realised (rather than expected) number of survivors on drift.

The replicates draw from a `Xoshiro` generator seeded per input (`emulate_sweep` derives each
input's seed from a base seed and the input's id by a fixed hash), so results do not depend
on threads or job order and resume exactly from the cache.

## Settings

- Primary closure `loci = :rare_alleles` (the closure of the first run, so that the comparison
  with it isolates the variance), burn-in 2000, horizon 2000, tail 240,000, as before.
- **R = 32** replicates per input.
- Base seed of the emulator's replicates: 20261115.

## Sample for E1–E5

A fresh uniform random sample of **20,000** of the 119,984 published inputs, drawn with
`randperm(MersenneTwister(20261115), 119984)[1:20000]` over the inputs in the order of
`validate/emulator.jl` (stratum, then file name). Observed outcomes, matching and table
conventions exactly as in `validate/emulator_criteria.md` (with its amendments 1, 2 and 4).
The deterministic emulator (the first run's cache, `validate/emulator_runs.csv`) is reported on
the same sample as the baseline.

## Criteria on the published table (thresholds of `validate/emulator_criteria.md`)

- **E1 discrimination.** AUC of P2000 for alive-at-2000 ≥ 0.85 overall and ≥ 0.80 in each of
  the 8 strata.
- **E2 calibration.** In 10 equal-count deciles of P2000 (ties broken by input order) the
  observed fraction alive is within 0.10 of the mean P2000 in every decile; and the
  calibration slope (logistic regression of the outcome on logit(P2000), P2000 clipped to
  [1e-6, 1 − 1e-6]) lies in [0.8, 1.25].
- **E3 timing.** Among inputs extinct within 2000 generations in the table with a predicted T:
  Spearman ≥ 0.7 and median T_pred/T_obs in [0.8, 1.25].
- **E4 half-survival rate.** Per stratum, x50 (logistic regression on log(MoveFitnessSD/√V_P);
  fractional response for the emulator) emulator/table ratio in [0.8, 1.25] in every stratum.
- **E5 (exploratory, no verdict).** AUC of P240k for survived-240,000.

## E6: calibration against fresh kernel simulations (the independent test)

**Inputs.** Chosen without looking at any table outcome, from the 2000 inputs of the first
run's parameter map (`validate/emulator_map.csv`, uniform sample of seed 20261114), which
holds each one's deterministic emulated half-survival rate k50. Let r = MoveFitnessSD/k50.
Candidates are those with finite positive k50 and 0.4 ≤ r ≤ 2. In each of the 8 strata,
the candidates are split into 4 bins of equal width in log r over [0.4, 2]; 6 are drawn at
random from each bin (`MersenneTwister(20261115)`, candidates in map order, `randperm`); a bin
with fewer than 6 contributes all of them and the stratum's remaining slots are filled at
random from its other undrawn candidates. 24 per stratum, **192 inputs**. (Range chosen because
the first run's E4 put the table's x50 at 1.2–1.7 times the emulator's: the inputs extend
below the deterministic k50 more than above it.)

**Kernel runs.** Each input is run **R_k = 40** times in the Julia kernel under the published
protocol (continuous move from generation 10000), started from the stationary spectrum
(`simulate(p; init = :equilibrium, burnin = 1000)`, roadmap F1: generations 9000 to 11999
are run), `NoGenerations = 12000`; fresh seeds `20261115000 + 100·j + r` for input j
(1-based, in selection order) and replicate r = 1…40. Alive at 2000 = not extinct by the end
of generation 11999 (death in generation 12000 or later, the table's "Generations to
extinction" > 2000). An extinction before the move counts as not alive.

**Statistics.** Per input, the kernel's fraction alive f_i (of 40) and the emulator's P2000
(R = 32, primary closure); the deterministic emulator (R = 0) is reported on the same inputs as
the baseline, with the same statistics and the same thresholds.

- **E6a binned calibration.** Inputs are grouped by the emulator's P2000 into the fixed bins
  [0, 0.1), [0.1, 0.3), [0.3, 0.5), [0.5, 0.7), [0.7, 0.9), [0.9, 1]. In every bin with at
  least 10 inputs, |mean f_i − mean P2000| ≤ 0.10.
- **E6b calibration slope.** Logistic regression of the 40 × 192 binary kernel outcomes on
  logit(P2000) of their input (P2000 clipped to [1e-6, 1 − 1e-6]): slope in [0.8, 1.25].
- **E6 PASS** iff E6a and E6b both pass.
- Reported without verdicts: mean |f_i − P2000|, the Brier score of P2000 against the runs'
  outcomes, the binomial log-likelihood Σ log Bin(s_i; 40, P2000_i) (P clipped as above), the
  intercept of E6b, and the AUC of P2000 for the runs' outcomes.
- Diagnostics of the noise model (no verdicts): across each input's kernel replicates alive at
  generation 9999, the variance of the mean phenotype and the coefficient of variation of V_G,
  against the emulator's spread of the mean at the move (Ploidy²·C, plus V_E/N for a sample of
  phenotypes) and its V_G (no spread, by construction), summarised as medians of ratios.

## Secondary (no verdicts)

`loci = :spectrum` with the noise (R = 8; it is about 4 s per deterministic input) on the E6
inputs, with E6's statistics, and on the 1000-input sample of the first run (seed 20261113) with
E1–E4's statistics if the CPU budget allows; how much of it is run is decided from the pilot
timing before any result is seen and recorded as an amendment.

## Budget

About 10 CPU-hours in all, on 8 threads. Sizes are checked against a pilot timing before the
run; if they do not fit, the sample sizes are reduced **before** the run and the change
recorded under Amendments.

## Amendments

1. **2026-11-13, after implementing the mean noise and BEFORE any comparison with the table
   or with kernel simulations: noise in V_G added to the primary model.** A smoke test of the
   mean noise alone on the synthetic 20-locus input of `test/emulator.jl` (not a published
   input) gave P2000 = 1.00, 0.55 and 0.00 at 0.9, 1.0 and 1.1 times the deterministic k50:
   still nearly a step, since the lag's noise is linear and Gaussian and hardly changes the
   trajectory of W̄. (The same smoke test timed four published inputs, Input-00103 of four
   strata; their P2000 was 0 with and without noise. No table outcome was looked at.) The
   primary model therefore also has random fluctuation of V_G, from theory with no fitted
   constant: the infinite-sites Poisson random field of the stationary spectrum
   (`stationary_spectrum`, src/init.jl, the law `init = :equilibrium` samples), in which a
   site of effect a at frequency x adds 2a²x(1 − x) to V_G (`QTLSim.vg_fluctuation`):
   - at the start of the move each replicate's V_G is shifted by δ ~ N(0, Var V_G),
     Var V_G = 8nLμ E_a[a⁴ K(γ)], K(γ) = ∫₀¹ x(1−x) e^{−γx(1−x)} dx (the stationary variance
     of the field; neutral: CV² = 1/(4nLμ));
   - every generation of the move its offspring's V_G is shifted by δ ~ N(0, 4Lμ E_a[a⁴ M(γ)]),
     M(γ) = ∫₀¹ (1 − 2x)² e^{−γx(1−x)} dx, the variance of the one-generation change of
     Σ 2a²x(1 − x) by drift of the frequencies (Var Δx = x(1 − x)/(2n));
   - each δ is added to the genic class of the replicate's haplotype moments (so V_G changes by
     δ up to O(1/N)), clamped so that the genic within-haplotype variance keeps at least 10% of
     itself, and the lowering's own recursion relaxes the deviation (drift at 1/(2n),
     selection, the rare-alleles loss). Both variances are computed once per input, at the
     stationary state (lag 0), and held through the move.
   - Still not included: linkage disequilibrium and associations among linked sites in the
     fluctuation (they widen the law of V_G), the change of the field's variances during the
     move, non-Gaussian phenotypes. `:spectrum` (no haplotype moments) gets the mean noise only.

   The model as first specified (mean noise only, `vgnoise = false`) is reported as a
   **secondary**, without verdicts: on the same 20,000 inputs with E1–E5's statistics and on
   the E6 inputs with E6's. The noise-model diagnostics of E6 also compare the kernel's
   coefficient of variation of V_G at generation 9999 with `vg_fluctuation`'s `cv`.
   Clarification of E6's selection: "candidates in map order" means sorted by input id
   (stratum, then file name), the order of the map's sample; one `MersenneTwister(20261115)`
   is used in turn for the strata in the order of `validate/emulator.jl` and their bins in
   increasing r.

2. **2026-11-13, after the pilot timing and BEFORE the run (no result seen): sample sizes
   reduced to fit the budget.** The pilot (`validate/emulator_variance.jl pilot`, 8 threads)
   took 0.133 s wall per random input for the primary (R = 32), 0.258 s per input near k50,
   0.43 s wall per kernel run and 1.7 s wall per input for `:spectrum` (R = 8): the sizes
   above would take about 2.4 h wall, ~19 thread-hours, against a budget of ~10 CPU-hours.
   Changes: the sample for E1–E5 is the first **10,000** (not 20,000) inputs of the same
   permutation (`randperm(MersenneTwister(20261115), 119984)[1:10000]`, sorted), still a
   uniform random sample; the mean-only secondary is run on the first 5,000 of that
   permutation (and on all E6 inputs); `:spectrum` is run on the E6 inputs only, not on the
   1000-input sample. E6 (192 inputs × 40 kernel runs) is unchanged. Thresholds unchanged.
