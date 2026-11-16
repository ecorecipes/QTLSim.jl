# Linked directional terms: pre-specified criteria

Committed before any run or comparison on the confirmation set J (2026-11-15). Script
`validate/linked_dir.jl` (+ `validate/linked_dir_analysis.jl`), theory `validate/linked_dir_theory.jl`,
algebra `validate/symbolic/check_linked_dir.jl` (50 checks, PASS before this file was committed).
Running log `validate/linked_dir.log` (hypotheses before each test; development in sections 1–5).

**Disclosure.** Everything below was developed on D (linked_decay's development set) and H′
(linked_ld's confirmation set), both seen in earlier studies, re-run with the directional
instrumentation (bitwise the same runs). The choice of level B‡ (which terms to close) was made after
seeing D and H′ (log, section 5); nothing was fitted. On development B‡ overpredicts Σ E[yY] late in
the move and does not close the late overprediction of ΣC, so (a) and (d) are expected to fail on J.

## Level B‡ (closed directional terms)

linked_ld's level B′* (`FIXES[:Bps]`) plus, from `validate/linked_dir_theory.jl` (`levelBdag`):
(i) K3_s = ⟨w⟩_l m_s in the T recursion's source (biallelic K3 = w_l y_s, w = κ4/v, variance-weighted
group mean ⟨w⟩ = Σκ4/Σv) and the group's measured Σκ3 in its coupling −c_s κ3/S_c; (ii) the own-skew
part of the y directional term, Q^T_st = Σ_l E[T_s y_t], carried with selection
+β[(⟨s²⟩_l + ⟨w⟩_l)P_st + h K_s m_t] − 3v̄Q^T/S_c, sampling (1 − 5ε)Q^T − 2εQ^Tᵀ + εκ3_aΣ_B, recombination
(1 − ρ_s)(1 − ρ_t), driven by a zeroth-order B′* run (no feedback), adding β(Q^T + Q^Tᵀ) to P and 2× to
H, βT̄ to m. Dropped: the partner part J of the y term, the K term βκ₅, (β² − 1/S)κ₆/2, three-locus
cumulants. Drivers: each run's measured β, S_c, n, N, Σv, Σv², Σκ3, Σκ4, Σκ3²/v, c, Cbb (single-locus
spectrum moments and background covariances; no association statistic). Start: B′*'s stationary state
for the drivers averaged over t = −60..−41, T = Q^T = 0.

Comparators (reported, no verdict): B′* (linked_ld), BdT (B‡'s (i) only).

## Inputs, runs, seeds

- **J (confirmation):** `j_rows()`: inputs of the published sweep not in
  `validate/spectrum_diagnosis_inputs.csv` (so in none of D, F, G, C, Dctl, Cctl, E6′) and not in
  linked_ld's H′, per stratum (NoQTL × u) in the order of a MersenneTwister(20261118) permutation of its
  files, screened on parameters only (k²/V_G ≥ 0.01 with V_G = `house_of_gauss_vg`, k/k_c ≤ 1.4 with
  k_c = `drift_critical_rate`), the first 3 that pass per stratum: 24 inputs (all 8 strata filled;
  listed by `validate/linked_dir.jl jrows`; computed before this file from parameters only).
  Exclusion checked by set intersection with D (`ld_rows("D")`), F (`fresh_rows()`), G (`g_rows()`), H′
  (`hrows()`), the whole spectrum-diagnosis input file (160 ids) and every input with raw records in the
  linked/genic/diagnosis studies' raw directories (104 ids): all empty (log, section 6). Some J inputs
  may have had kernel runs in other studies of the published sweep with other seeds; their association
  statistics were never looked at.
- R = 3 runs per input, fresh seeds 60261116000 + 100j + r (j = 2001..2024); `init = :equilibrium`,
  1000 generations of burn-in, instrumented from MOVE − 60 to MOVE + 199 with `dir_run`.
- The task asks for at least 15 eligible inputs: read here as at least 15 inputs passing the
  parameter screen (24 do); the outcome eligibility below cannot be known before the runs.

## Statistics

As linked_source/linked_ld: window means over runs alive, then over generations; medians over inputs;
eligible inputs: measured cumulative 2βΣC over 0..99 ≤ −0.05 × genic variance before the move;
lr = log(predicted/measured), undefined (a failure of that input in the |lr| statistics) when of
opposite sign. Ratios are predicted/measured of window means.

## Predictions and verdict rules (on J; level B‡)

- **(a) Allele–background LD:** Σ E[yY] ratio, median over all J inputs, in 50..99 and in 100..199,
  within [0.75, 1.33]. PASS if both, PARTIAL if one, else FAIL.
- **(b) Variance heterogeneity:** Σ H ratio at 50..99, median over all inputs, within [0.75, 1.33]:
  PASS, else FAIL.
- **(c) Source:** Σ K ratio, median over all inputs, in −60..−1 and in 50..99, within [0.75, 1.33].
  PASS if both, PARTIAL if one, else FAIL.
- **(d) Closed association, eligible inputs:** ΣC in 25..49, 50..99 and 100..199 (100..199 is a verdict
  window here): median lr within ±0.3 and median |lr| ≤ 0.4 in each; cumulative loss 0..99, 2β ×
  predicted ΣC against the measured sign-symmetric loss 2βΣC_sym: median lr within ±0.3 and median
  |lr| ≤ 0.4. PASS if all four, PARTIAL if two or three, else FAIL.
- **(e) Decay rate (P1, level A, unchanged):** `:derived` median lr within ±0.2 and median |lr| ≤ 0.3 in
  each of 25..49, 50..99, 100..199, and its median |lr| below both comparators' (`:nosel`, `:spectrum`)
  in 50..99 and 100..199. PASS if all; PARTIAL if the accuracy conditions hold in 2 of 3 windows
  (comparator condition still required); else FAIL.
- Fewer than 6 eligible J inputs: (d) and (e) INCONCLUSIVE.

## Secondary pre-specified claim: the diagnosis (not a closure)

B′* driven additionally by the measured directional terms (β δP bin matrix, β δK rows, β δy) and the
measured K3 and Σκ3 in the T step (`levelBdir`, variant Bdt; log, sections 2–3): the rules of (d) on
J's eligible inputs. Expected (H′: 25..49 +0.11/0.27, 50..99 −0.12/0.37, 100..199 +0.03/0.49, loss
+0.17/0.55): PASS in 25..49 and 50..99, median within ±0.3 at 100..199. Reported as a verdict on the
diagnosis, with the same PASS/PARTIAL/FAIL rule; B′* driven by the measured directional terms only
(Bd) reported alongside.

Reported, not verdicts: the comparators; the measured size of the directional terms on J (as log
section 1); the exploratory driven-spectrum check only if (d) passes.

## Amendments

None. Any amendment will be recorded here, with its reason, before the J runs are analysed.
