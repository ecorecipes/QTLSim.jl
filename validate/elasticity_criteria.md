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

## Phase B addendum (2026-11-18): a target that does not saturate

Written after phase A's results (above), and **before any elasticity of a sampled input was
computed on the new target**. Committed, with the code it describes (src/k50_gradient.jl,
validate/elasticity_b.jl, test/k50_gradient.jl), before validate/elasticity_b.jl `run`. Output:
validate/elasticity_b.log. Phase A's log, caches and script are not touched. G1–G3, their
thresholds and tolerances, the cells, the terciles, the table side, the sample and its draw order,
and the minimum of 4 defined inputs per cell are **unchanged**.

**What had been looked at when this was written.** (i) A diagnosis of the saturation on five of
phase A's 12 undefined inputs (QT0500-Possion5.0/Input-09763, QT5000-Possion5.0/Input-11353,
QT0100-Possion5.0/Input-14634, QT0500-Possion3.0/Input-03628, QT1000-Possion5.0/Input-10781): the
screen's trajectory at rates within e^{±0.003} of phase A's k50 (mean fitness, per-generation
hazards, the generation of deterministic extinction, P2000); no derivative in any input and no
elasticity. (ii) Development of `k50_gradient` on synthetic parameters (test/k50_gradient.jl) and on
three table inputs **outside the sample** (QT0100-Possion3.0/Input-14211, QT0100-Possion5.0/
Input-01948 and Input-09336, ω = 0.65, 0.47, 1.10), where its elasticities in V_P and ω agreed with
a brute-force re-solve to 0.003. Nothing else.

**Where the saturation came from** (diagnosis). At a fixed size the probability of being alive is a
product over generations of 1 − h_t, h_t = P(Binomial(N, W̄_t) < 2). With narrow ω the trajectory
sits at a steady lag with h_t far below 1e-16 until it runs away, and then dies within a few
generations. So (a) P2000 is exactly 1 in Float64 while the cumulative hazard H is below 1.1e-16
(H from 1e-13 down to 1e-102 was seen at rates 3e-4 to 3e-3 below k50 in log k): the product does
not move; (b) P2000 is exactly 0 once the **deterministic trajectory itself dies inside the
horizon** (expected survivors below 2), which `extinction_distribution` takes as certain
extinction, although the hazard accumulated until then is finite and can be far below log 2 (0.026
to 5.6 in the five inputs): a jump of the target, not rounding; and (c) between the two, log H
moves through its whole range over 1e-3 to 1e-5 in log k (|∂h/∂log k| up to 1.6e5 on a development
input, against phase A's median 193), so a step of 1e-2 or 1e-3 lands on (a) on one side and (b) on
the other. Not a cause: the truncations of the exact size chain (1e-17, 1e-18), which is not used at
a fixed size; the breeders are not rounded in the hazard.

**The new target.** h(k, θ) = log(log 2) − log H, with H = Σ_t −log(1 − h_t) over the horizon summed
in logs (log h_t in closed form; `QTLSim.log_cumulative_hazard`), along the trajectory **continued**
through its deterministic extinction: log W̄ extrapolated quadratically from its last three
generations, first and second differences capped at zero (`QTLSim.emulated_log_hazard`). It equals
phase A's h wherever that was finite and the trajectory survived the horizon. Where the trajectory
dies inside the horizon the two differ by construction: phase A's P2000 = 1/2 was then the jump (b),
the new k50 is the root of the continued target, within a fraction of a generation of the runaway
of it (a relative difference in k of the order of 1/(|∂h/∂log k|), i.e. 1e-3 to 1e-6). Limits,
accepted in advance: the continuation is an extrapolation (it matters for one or two generations);
the report counts the inputs that rely on it.

**Method** (`k50_gradient(p; names = (V_P, h², m, ω))`, defaults, in place of method steps 1–3 of
phase A): k50 is the root of h in log k to |h| ≤ 0.02 (bracket: phase A's k50·e^{±0.006} for its 64
inputs — a search result, not an elasticity; (k50_rare/4, k50_rare) for the others, as in phase A;
moved toward the root if it misses it). e_θ = −(∂h/∂log θ)/(∂h/∂log k) by central differences whose
step is chosen for the difference in k to span about ±0.5 in h (at most 1e-2, at least 1e-10), the
same step for the four inputs (reduced for one whose difference spans more than ±2). **Undefined**
(counted, never imputed): no k50; no root (the bracket closes on a discontinuity); a difference
point not finite at the smallest step; ∂h/∂log k ≥ 0; a fallback switch as below; more than 200
evaluations.

**Fallback of the beta solver** (replaces step 4 of phase A; needed because the steps are now up to
10⁴ times smaller, so a jump inside a difference is no longer small beside it). A difference is used
only if its two points and the base point fell back to the grid in the same generation, or none
did; otherwise its step is divided by four, up to three times. An input for which a pair still
differs ("fallback switch"; this now includes phase A's "shifts") is recomputed with the grid solver
(`spectrum = (; ld = :dynamic)`) at the grid's own k50 (bracket k50·e^{±0.02}), up to 20 such inputs
in sample order (phase A's cap of 10 for 64 inputs, scaled to 128); beyond 20, or if the grid's
result is undefined, the input is **undefined**. The primary ("consistent") value is the grid's for
those inputs and the beta's otherwise. As a check of the size of the beta–grid difference (no
verdict), the first 3 defined inputs in sample order whose run at k50 fell back are also recomputed
with the grid.

**Sample.** The 64 inputs of phase A (rounds 1–8) are re-evaluated on the new target, same ids. The
sample is extended toward the pre-specified 16 per stratum = 128, with the ids already drawn
(validate/elasticity_table.log, elasticity_raw/sample.tsv), in draw order: rounds 9–12 are always
run (96 inputs); each of rounds 13–16 starts only if the process CPU time of the run (clock(), user
+ system over all threads) is then below 3.0 hours (a rule on cost alone; the whole of phase B has
about 5 CPU-hours). No more than 128. Six threads.

**Verdicts.** G1, G2, G3 and the overall gate are recomputed as written above over all 48 cells, on
the consistent values of **all inputs run** (the verdict), with the same exclusion of cells with
fewer than 4 defined inputs (listed). Reported beside it without verdict: the beta values alone;
phase A's 64 inputs alone on the new target; `:rare_alleles` (phase A's values, old target, not
recomputed). **The four low-ω tercile cells (ω, t1, one per locus number), which phase A could not
test, are listed one by one** with every input in them. These verdicts supersede nothing: phase A's
verdicts on 31 cells stand as logged; phase B's are the verdicts on the cells it could not see and
on the larger sample.

**Numerical checks** (same thresholds as phase A; if N1 or N2 fails the verdicts are reported and
flagged "numerically unreliable"):

- N1 (step), round 1 (8 inputs): the elasticities with the span ±0.5 against a span of ±0.125 (steps
  four times smaller) at the same k50. Stable if the median |difference| ≤ 0.02 and the 90th
  percentile ≤ 0.05 over the values defined at both.
- N2 (brute force), redone: the brute-force elasticities are **phase A's** (elasticity_raw/
  screen_brute.csv: k50 by the P2000 search re-solved at θ·e^{±0.05}; five inputs × four inputs, all
  20 finite; not recomputed). The rule, fixed now: over **every value that brute force defines**
  (finite), the implicit value on the new target (consistent value) must be defined, and the
  |brute − implicit| must have median ≤ 0.03 and maximum ≤ 0.10. One undefined implicit value
  where brute force is defined makes N2 NOT VERIFIED.

**Added descriptive checks (no verdict).**

- C1: on phase A's 52 defined inputs, the new elasticities against phase A's. Expected: agreement
  within 0.02, except possibly where phase A's pairs fell back in different generations (its 10
  "shift" inputs, whose phase A values carried a jump; reported separately).
- C2: e_V_P against e_m. Phase A's values nearly coincide. Analytically (stated now, checked in the
  report): the emulator reads V_P, h² and m only through V_E = (1 − h²)V_P (the environmental
  variance in V_s = ω² + V_E and in V_P = V_G + V_E) and V_m = m·V_E (through the variance of a
  mutation, m·V_E/(2Lμ)); the target h²V_P is not used, since V_G is set by mutation, selection and
  drift from the house-of-Gauss start. Hence two exact identities: e_h² = −h²/(1 − h²)·e_V_P, and
  e_V_P − e_m = the elasticity in V_E at fixed V_m and ω. The latter is small when ω² ≫ V_E (V_E ≤
  0.09 here, so for ω above about 1), and **should differ from zero at low ω**; if ω and V_E acted
  only through V_s it would equal e_ω·V_E/(2ω²) (they do not exactly: ω also enters the height
  ω/√(ω² + V_P) of mean fitness). The report gives the identity in h², e_V_P − e_m by ω class, and
  its relation to V_E/(ω² + V_E) and to that prediction.

**Secondary.** S1, S2 and the level comparison of S4 are recomputed on the new values (no verdict);
S3 (default `:spectrum`) and S4's regression are not rerun.

### Amendments to the addendum

(none at the time of its commit)

No rule was amended after the commit. Notes on the run (2026-11-18, written after the results were
seen; they change no rule and no verdict):

- All 16 rounds ran (128 inputs): the run took 1.49 CPU-hours (clock()), below the 3.0 that would
  have stopped it. No input kept a fallback switch after the step reductions, so no input was
  recomputed with the grid for the consistent values (they are the beta's); the grid check ran on 3.
- 126 of 128 inputs are defined, among them all 12 that phase A left undefined. The two undefined
  (status no_root, ω = 0.05 and 0.16) have a jump of the trajectory itself at k50: within 2e-9 in
  log k the target goes from h = 84–111 to a trajectory that dies more than 1000 generations before
  the end of the horizon (post hoc lines in the log). The continuation does not help there.
- 10 inputs (all with ω < 2.2) relied on the continuation.
- N2 is NOT VERIFIED under the rule fixed above: 19 of the 20 values agree within 0.0065, one does
  not (e_ω of QT5000-Possion5.0/Input-11353, ω = 0.37: implicit −0.326, phase A's brute force
  −1.087). The log's flag "numerically unreliable" is the rule applied as written. Post hoc (mode
  `posthoc`, added after the results): brute force on the new target over ω·e^{±0.05} reproduces
  −1.088, over ±0.02 and ±0.005 gives −0.328, so k50(ω) is not smooth over ±5% there and the
  brute-force comparator, not the implicit value, is off; this is an explanation, not a verdict.
- C2's identity in h², called exact above, holds to 4.5e-4 in 90% of the inputs but is off by
  0.007, 0.007 and 0.058 in three inputs with ω ≤ 0.37 (listed in the log); not diagnosed
  (numerical error at very steep targets, or a small dependence on h²V_P that the argument missed).
- Code added after the results: the `posthoc` mode of validate/elasticity_b.jl and the "post hoc"
  section of its report.
