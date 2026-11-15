# Linked LD: pre-specified criteria

Committed before any run or comparison on the confirmation set H′ (2026-11-15). Script
`validate/linked_ld.jl`, theory `validate/linked_ld_theory.jl`, algebra
`validate/symbolic/check_linked_ld.jl` (39 checks, PASS before this file was committed). Running
log `validate/linked_ld.log` (hypotheses stated before each test; development results in its
sections 1–4).

**Disclosure.** The diagnosis and the corrected recursion were developed on D, F and G, all seen
in earlier studies (linked_decay, linked_source): the stage accounting on D, F, G (log, section 1),
the class-space accounting and the full measured class moments on D (sections 2–3), the corrected
recursion driven by the runs of D, F, G (section 4). Nothing was fitted: every term is derived, but
which terms to keep was decided after seeing D, F, G. Only H′ gives verdicts.

## The corrected recursion (level B*; with linked_source's H coupling, level B′*)

`yfix_step!` with `between`, `pairing`, `mutfull`, `sind`, `vhet` all on (`FIX_ALL`), on top of
linked_decay's level B (`y_step!`) and linked_source's level B′ (`yH_step!`, `hsel = :bulmer`):

- `vhet` (the main term): the cross term of the mean LD and Bulmer's offset in the summed moments,
  Σ_l A E[y_l] b_lᵀ with b_l = −c v_l/S, uses Σ_l E[y_l] v_l = (Σv²/Σv) m (E[y_l] ∝ v_l under the
  recursion itself) instead of v̄ m;
- `between`: selection on individuals creates the same association between the two chromosomes
  of an individual as within each; recombination passes it to the gamete, so selection's change
  enters the gamete's LD in full and recombination acts only on the inherited part;
- `pairing`: the recombinant's cross-chromosome product under random pairing of the survivors'
  chromosomes, mean −ρ∘u/(2n − 1) and covariance [v_l Cov(R_s, R_t) + (ρ∘u)(ρ∘u)ᵀ]/(2n);
- `mutfull`: mutation in the background, v_l Σ_{k∈s} m_k/(2N) on the diagonal;
- `sind`: Bulmer's exact Gaussian curvature 1/(ω² + V_P) for 1/S_c.

Dropped: the directional third-cumulant term b κ(B_s, x_l, G), the fourth-cumulant term,
O(1/n²). Start, drivers and T recursion as linked_decay's level B and linked_source's level B′.

## Inputs, runs, seeds

- **H′ (confirmation):** `h_rows()`: inputs of the published sweep not in
  `validate/spectrum_diagnosis_inputs.csv` (so in none of D, F, G, C, Dctl, Cctl, E6′), per stratum
  (NoQTL × u) in the order of a MersenneTwister(20261117) permutation of its files, screened on
  parameters only as linked_source's G (k²/V_G ≥ 0.01 with V_G = `house_of_gauss_vg`, k/k_c ≤ 1.4
  with k_c = `drift_critical_rate`), the first 3 that pass per stratum: 24 inputs (listed by
  `julia --project=vignettes validate/linked_ld.jl hrows`; computed before this file, from parameters
  only). Disclosure: some of these inputs may have had kernel runs in other studies of the published
  sweep (published_sweep, vg_move, emulator checks) with other seeds; their association statistics
  were never looked at.
- R = 3 runs per input, fresh seeds 70261115000 + 100j + r (j = 1001..1024); `init = :equilibrium`,
  1000 generations of burn-in, instrumented from MOVE − 60 to MOVE + 199 with linked_source's
  `ls_run` (no decomposition pass, as G).

## Statistics

As linked_source: window means over runs alive, then over generations; eligible inputs: measured
cumulative 2βΣC over 0..99 ≤ −0.05 × genic variance before the move; lr = log(predicted/measured),
undefined (a failure for that input in the |lr| statistics) when of opposite sign. Medians over
inputs.

## Predictions and verdict rules (on H′)

- **(a) Allele–background LD, level B*:** predicted Σ E[yY] / measured Σ yY, over all H′ inputs, in
  windows −60..−1 and 50..99: median ratio in [0.80, 1.25] and median |lr| ≤ 0.30 in each window.
  PASS if both windows, PARTIAL if one, else FAIL. (Development: D 1.00/0.86, |lr| 0.07/0.16;
  G 1.02/0.91, 0.07/0.09; F 0.99/0.90, 0.11/0.16; linked_decay's level B on the same sets
  0.76/0.66, 0.90/0.71, 0.92/0.87.)
- **(b) Variance heterogeneity, level B′*:** predicted Σ H / measured Σ H over all H′ inputs, windows
  −60..−1 and 50..99: median ratio in [0.67, 1.5] and median |lr| ≤ 0.4 in each window. PASS if both,
  PARTIAL if one, else FAIL. (Development: D 0.87/0.67, G 0.62/0.54; linked_source's B′ 0.40/0.24,
  0.23/0.21.)
- **(c) Closed association, level B′*, eligible inputs (linked_source's (b) rules):** ΣC in 25..49
  and 50..99: median lr within ±0.3 and median |lr| ≤ 0.4; cumulative loss 0..99, 2β × predicted ΣC
  against the measured sign-symmetric loss 2βΣC_sym: median lr within ±0.3, median |lr| ≤ 0.4. PASS
  if all three, PARTIAL if two, else FAIL. (Development: D passes all three; G fails 25..49 and
  50..99 as linked_source's B′ did; the corrected y does not change C much.)
- **(d) Re-checks.** Source: predicted ΣK (level B′*) / measured ΣK, median over all inputs in
  −60..−1 and 50..99 within [0.75, 1.33]: PASS if both (linked_source's (a)). Decay rate: P1 of
  linked_decay/linked_source (level A, independent of y), unchanged: `:derived` median lr within ±0.2
  and median |lr| ≤ 0.3 in each of 25..49, 50..99, 100..199, and its median |lr| below both
  comparators' (`:nosel`, `:spectrum`) in 50..99 and 100..199; PASS if all, PARTIAL if the accuracy
  conditions hold in 2 of 3 windows (comparator condition still required), else FAIL.
- Fewer than 6 eligible H′ inputs: (c) and P1 INCONCLUSIVE.

Reported, not verdicts: window 100..199 of (a)–(d) (affected by the dropped fifth-order directional
term), level B* for ΣC, linked_decay's level B and linked_source's level B′ as comparators, each
correction alone on development data.

## Amendments

None. Any amendment will be recorded here, with its reason, before the H′ runs are analysed.
