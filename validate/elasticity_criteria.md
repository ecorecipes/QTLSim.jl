# Pre-specified criteria: are the emulator's gradients worth having?

Written 2026-11-18, **after** the table's own elasticities were fitted (validate/elasticity.jl
`table`, validate/elasticity_table.log; no emulator is called in that mode) and after a timing
run of single emulator calls (`timing`; seconds only, no outcome and no elasticity looked at), and
**before** any emulator elasticity or any half-survival rate of the sampled inputs was computed.
Committed before validate/elasticity.jl `emulator` is run. Output: validate/elasticity.log.
Phase A of the automatic-differentiation plan: no AD, no change to src/, no default changes
whatever the outcome.

## Question

Do the emulator's elasticities of the half-survival rate k50 (the rate of the continuous move at
which half the populations are alive 2000 generations into it) with respect to the inputs agree
with the elasticities implied by the published table?

## What the table varies (established in `table` mode)

All 119,984 input files match a unique table row. Per file, five inputs are drawn, independently
(largest |Spearman r| between two of them 0.019): the rate k = `MoveFitnessSD` (continuous move,
per generation; roughly uniform on 0–0.3), `PhenotypeVar` V_P (0.01–0.09), `Heritability` h²
(0.01–0.99), `MutationScalar` m (1e-4–0.05) and `InitFitnessSD` ω (0.003–10). **ω is given per
file; it is not the analytical width** h²V_P/√(m·V_E), so it does not move with the others.
`MutationRatePerGamete` (1e-4), `PopSize` (1000) and everything else are constant; `NoQTL` and
`RecombinationPoissonMean` define the eight strata. So the four inputs are θ = (V_P, h², m, ω),
each one `Params` field, and an emulator derivative perturbs that one field alone
(`Params(p; field = value·e^{±δ})`), all other fields held, as the table's inputs co-vary.

**Units.** The rate is raw, k = `MoveFitnessSD` per generation. (For the standardised rate
k/√V_P, the elasticity in V_P is lower by exactly 1/2 and the others are unchanged.)

**Outcome.** Alive at 2000 = the table's "Generations to extinction" > 2000 (as E1–E4 of
validate/emulator_criteria.md).

## What the table fits showed, and why the primary cells are terciles

Per stratum, a logistic regression of alive-at-2000 on log k and the logs of the four inputs
gives e_θ = −β_θ/β_k (bootstrap CIs, 1000 resamples; all 32 exclude 0). **These elasticities are
not constant**: a quadratic model improves the fit hugely (likelihood ratio 370–560 on 15 df in
every stratum), and fits within terciles of each input differ by more than 3 SE in 23 of 32
(stratum, input) cells. e_h² runs from about 0 at low h² to about −0.7 at high h², e_m from 0.2
to 0.6–0.7, and e_ω changes sign (about −0.2 at small ω, +0.3 at large ω); only e_V_P is roughly
constant (0.4–0.8, differences within noise). A single per-stratum number is then a
leverage-weighted average over the input's whole range (dominated by where log θ spreads most,
e.g. low h²), which no average of local gradients over a sample of inputs estimates: the table's
own quadratic-model local elasticities, averaged over near-boundary inputs, differ from its
linear-fit number by 0.1–0.17 for h² and ω. The quadratic model itself is not adequate as a local
truth either (adequacy 3 in the log: it compresses the dependence on h² and m and misses the sign
change in ω). Therefore, as the task allows ("compare the emulator's elasticities averaged over
the same sub-range"), **the primary comparison is by sub-range**:

**Primary cells (48).** For each locus number L ∈ {100, 500, 1000, 5000} (its strata u = 3 and
u = 5 pooled, with an indicator of u = 5; the per-stratum elasticities of the two differ little)
and each input θ: three terciles of θ, the terciles being those of the table's information weights
μ(1 − μ) (μ the per-stratum linear fit's probability), so each holds a third of the table's
information about the boundary. Table side: the linear logistic model fitted to the rows of the
tercile; e_θ of θ itself, with a bootstrap percentile 95% CI (600 resamples of the cell's rows).
All 48 cells are defined; 41 CIs exclude 0; 40 have |e| ≥ 0.1 (validate/elasticity_table.log,
cached in elasticity_raw/table_cells.tsv). Limitations, accepted in advance: within a tercile the
fit is still a leverage-weighted average over that tercile (and over the other inputs), while the
emulator's side is a median of local values; the table side has its own uncertainty (CIs of
half-width 0.03–0.15).

## The emulator's sample

Per stratum, 16 inputs drawn without replacement with probability ∝ μ(1 − μ) of the stratum's
linear table fit (Efraimidis–Spirakis keys, `MersenneTwister(20261119 + stratum)`), in a fixed
draw order: the inputs that carry the table's information about the boundary (a uniform sample
would be 92% inputs far from it). Ids: validate/elasticity_table.log, elasticity_raw/sample.tsv.
Nothing is fitted to the emulator on them, so freshness is not at issue; the table's outcomes of
these rows were used only through μ. Each input's elasticities are evaluated at the emulator's
**own** k50 for that input's (V_P, h², m, ω), not at the file's rate.

**Budget and the size actually run (decided on the timing run, before any result).** One screen
evaluation costs 1–4 s, the same closure on the grid 6–8.5 s, default `:spectrum` 4–6 s; an input
needs about 12 evaluations for k50 and 11 for the derivatives. With about 3 CPU-hours in all
(about 35 CPU-minutes already spent on the table fits and bootstraps), the brief's 25 per stratum
does not fit. So:

- **screen**: the first 12 of each stratum's 16, in draw order, processed in rounds of one input
  per stratum; after round 8, no new round starts once the summed job seconds of the screen exceed
  6000 (a rule on cost alone). So 64–96 inputs.
- **`:rare_alleles`** (cheap): all 128; its statistics are reported on the screen's inputs.
- **default `:spectrum`**: round 1 only (8 inputs, one per stratum): too few for cells; reported
  input by input beside the screen.

## Method (emulator)

Closures: (i) the two-stage screen `emulate(p; loci = :spectrum, spectrum = (; ld = :dynamic,
solver = :beta))`; (ii) default `emulate(p; loci = :spectrum)`; (iii) `emulate(p)` (`:rare_alleles`).
Burn-in and horizon 2000.

1. k50 by `emulated_half_survival_rate` (same `loci`/`spectrum`): `:rare_alleles` with the default
   bracket to `rtol = 1e-6`; the others from the bracket (k50_rare/4, k50_rare) (the function
   widens it if it misses) to `rtol = 1e-2` in log k. That is enough: the formula below is the
   elasticity of the level set through the point, and P2000 is nearly a step in k.
2. h(k, θ) = log(log 2) − log(−log P2000), the search's own function (the log cumulative hazard).
   At k50: e_θ = −(∂h/∂log θ)/(∂h/∂log k), each derivative a central difference with step
   δ = 1e-2 in the log (10 evaluations, plus one at k50 itself for the record).
3. **Undefined** (counted, never imputed): status 1, no k50 (0 or ∞); status 2, P2000 exactly 0 or
   1 at one of the ten points (h infinite) — retried once with δ = 1e-3, undefined if still so;
   status 3, ∂h/∂log k ≥ 0.
4. **Fallback of the beta solver** (a discontinuity). Recorded at every evaluation. A difference
   pair *straddles a switch* if one of its two points fell back to the grid and the other did not
   (`sstat`), and *shifts* if both fell back at different generations (`sgen`). An input with any
   straddling pair has all its ten difference evaluations recomputed with the grid solver
   (`spectrum = (; ld = :dynamic)`) at the screen's k50 (if the grid saturates there, at the grid's
   own k50, bracket k50·e^{±0.2}); up to 10 such inputs in sample order; any beyond 10 are
   **undefined** for the primary statistics (counted). Inputs with only shifts keep the beta
   solver's values; the first 3 of them (sample order) are also recomputed with the grid as a
   check of the size of the effect (reported). The primary ("consistent") screen value is the
   grid's for the straddling inputs and the beta's otherwise; the beta's own values for all
   inputs are reported as "screen as run" (no verdict).
5. **Numerical checks of the screen**, on round 1 (the first sample input of each stratum):
   - N1 step: the elasticities with δ = 1e-3 at the same k50. Stable if the median
     |e(1e-2) − e(1e-3)| ≤ 0.02 and its 90th percentile ≤ 0.05 over the 8 × 4 values.
   - N2 brute force: for the first sample input of strata 1, 3, 5, 7 and 8 (five inputs), k50
     re-solved at θ·e^{±0.05} for each input (bracket k50·e^{±0.15}, `rtol = 1e-3`), e = Δlog k50/0.1.
     Verified if the median |brute − implicit| ≤ 0.03 and the maximum ≤ 0.10 over the 5 × 4 values.
   If N1 or N2 fails, the verdicts below are still computed and reported, and flagged
   "numerically unreliable"; no step is changed after the fact.

## Statistics and verdicts

An input belongs to the cell (L, θ, tercile) of its own value of θ. The emulator's cell value is
the **median** over the cell's inputs with a defined elasticity (the mean is reported too). A cell
with fewer than 4 defined inputs is excluded from G1–G3 and listed.

- **G1 (sign).** Among cells whose table CI excludes 0: the screen's median has the table's sign.
  PASS if ≥ 90% of those cells.
- **G2 (size).** Among cells with |table e| ≥ 0.1: the screen's median lies in
  [lo − tol, hi + tol], with [lo, hi] the table's bootstrap CI and tol = max(0.15, 0.30·|table e|).
  PASS if ≥ 75% of those cells. Every failing cell is listed.
- **G3 (ranking).** Spearman correlation between the screen's medians and the table's e across
  all cells (not excluded) ≥ 0.8.
- **Overall gate: PASS iff G1 and G2 pass.** G3 is reported with its own PASS/FAIL.
- A table cell that is undefined (non-converged fit, βk ≥ 0, a coefficient beyond 1e3, or fewer
  than 90% usable bootstrap fits) would be excluded; there is none.

**Comparators, same statistics, no verdict:** `:rare_alleles` on the screen's inputs; the screen
"as run" (beta values for every input).

**Secondary, no verdict:**

- S1, the literal per-stratum reading (32 cells): the per-stratum linear-fit e_θ with its CI
  against the screen's median and mean over the stratum's inputs, with G1–G3's numbers for
  information. Expected in advance to be distorted by the non-constancy described above.
- S2, input by input: the screen's local elasticity against the table's quadratic-model local
  elasticity at the same input (elasticity_raw/sample.tsv): Spearman and median |difference| per
  θ. The quadratic model is known to be poor for ω.
- S3, default `:spectrum` beside the screen on the 8 round-1 inputs (differences per θ).
- S4, levels: the emulators' k50 against the table's boundary rate for the same inputs (quadratic
  fit), median ratio by L; and a least-squares regression of the screen's log k50 on the logs of
  the four inputs (per L, with a u indicator) beside the per-stratum linear table fits (a global
  slope from levels, not a gradient).
- Counts: undefined by status, fallbacks, straddles and shifts, grid recomputations, CPU seconds.

## Amendments

Any change after the commit of this file is recorded here with its reason and whether any
emulator result had been seen.

(none at the time of the first commit)

No criterion was amended after the commit. Notes on the run (2026-11-18, written after the
results were seen; they change no rule and no verdict):

- The screen stopped after round 8 (64 inputs) under the cost rule: its summed job seconds were
  7603 > 6000. Job seconds are wall time inside 8 concurrent threads and overstated the CPU time
  (the whole emulator mode took 94 CPU-minutes by `/usr/bin/time`), so the rule stopped earlier
  than the budget required. Not extended after the fact.
- N2 is logged "NOT VERIFIED" because the rule asks for all 5 × 4 values and one of the five
  brute-force inputs has an undefined (saturated) implicit elasticity; the 16 comparable values
  agree to 0.0065. The log's flag "numerically unreliable" is the rule applied as written.
- 12 of the 64 screen inputs are undefined (P2000 exactly 0 or 1 at a difference point even at
  δ = 1e-3). They all have ω < 2.2 (post hoc line in the log), so the undefined inputs are not a
  random subset: 17 of the 48 cells have fewer than 4 defined inputs and are excluded, among them
  all four low-ω cells. The verdicts cover 31 cells.
- The only code added after the results is the "post hoc" block at the end of `analyse`.
