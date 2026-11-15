# Linked source: pre-specified criteria

Committed before any comparison on the confirmation set G (2026-11-15). Script
`validate/linked_source.jl`, theory `validate/linked_source_theory.jl`, algebra
`validate/symbolic/check_linked_source.jl` (PASS before this file was committed). Running log
`validate/linked_source.log`.

**Disclosure.** The closure below was developed on set D (linked_decay's development inputs and
seeds, re-run with the extended instrumentation) after the decomposition of K on D and after
comparing two candidate closures with the kernel on D (log, sections 1–3). F (linked_decay's
confirmation set, whose results have been seen) may be analysed as further development and
is reported, not a verdict. Only G gives verdicts.

## The derived prediction (level B′)

Exact: K_s = κ(B_s, B, x_l, x_l) = H_s − 2y_sY_l, H_s = Cov(B_sB, d_l²) (variance heterogeneity).
Per group a of loci and class pair (s, t) (16 position bins × 7 split classes), state P_st = Σ_{l∈a}
E[y_s y_t], H_st = Σ_{l∈a} E[Cov(B_sB_t, d_l²)], K_st = H_st − 2P_st; source K_s = Σ_t K_st feeds
linked_decay's T recursion (decay unchanged). One generation:

1. selection (first order in 1/S_c): P by linked_decay's Bulmer step; the fourth cumulant
   K_st ← K_st − (c_s K_t· + c_t K_s·)/S_c (leading Leonov–Shiryaev terms of the stabilising
   selection change; c_s = Cov(B_s, G_chr), K_s· = Σ_t K_st); then H = K + 2P. Dropped: the
   directional β κ₅ term, κ₆, κ(B_s, l, l, l) (K3), products of third-order cumulants, the
   2Cov(x_l, G)κ(B_s, B, l, G) term (v_l ≪ Var(B)).
2. survival sampling, ε = 1/(2n) − 1/(2N): P ← (1 − 3ε)P + ε(H + Σ_B V), H ← (1 − 5ε)H + 2ε(P + Pᵀ)
   (exact to O(ε²) for iid sampling; Ohta & Kimura's E[D(1−2p)(1−2q)] equation in allele values).
3. gametes (single crossover): P ← RPR, H_st ← (1 − ρ_st)H_st, ρ_st the three-locus split
   probability (max on one side, sum across).
4. offspring sampling, ε = 1/(2N), as 2; mutation adds to P only (E[ΔH] = 0 exactly).

Start: joint stationary state of (P, H) for the drivers averaged over t = −60..−41, T = 0 at
t = −60; driven by each run's measured β, S_c, n, N and variances (as linked_decay's level B).

Comparators (reported, no verdict): level B (linked_decay, K = −2E[yY]); B′ without the
selection term on K (drift seeding only).

## Inputs, runs, seeds

- **G (confirmation):** from the 102 inputs of `validate/spectrum_diagnosis_inputs.csv` in no
  set and not in linked_decay's F, those passing a screen computed from parameters only (no run,
  no outcome after the move): k²/V_G ≥ 0.01 (k the rate of the move, V_G = `house_of_gauss_vg`)
  and k/k_c ≤ 1.4 (k_c = `drift_critical_rate`); at most 3 per stratum (NoQTL × u) in the order
  of MersenneTwister(20261115) (`g_rows()`): 22 inputs. The screen was chosen on D and F (it
  separated their eligible from their ineligible inputs, 16/16, with no ineligible input
  passing); disclosed. Disclosure as for F: these inputs' kernel survival (E6′) was measured in
  earlier studies with other seeds; their association statistics were never looked at.
- R = 3 runs per input, fresh seeds 90261115000 + 100j + r; `init = :equilibrium`, 1000
  generations of burn-in, instrumented from MOVE − 60 to MOVE + 199 (no π pass, no Gram).

## Statistics

As linked_decay: window means over runs alive, then over generations; eligible inputs: measured
cumulative 2βΣC over 0..99 ≤ −0.05 × genic variance before the move; lr = log(predicted/measured),
undefined (a failure for that input in the |lr| statistics) when of opposite sign.

## Predictions and verdict rules (on G)

- **(a) Source identity:** median over all G inputs of predicted ΣK (level B′) / measured ΣK in
  windows −60..−1 and 50..99 within [0.75, 1.33]. PASS if both, else FAIL.
- **(b) Closed level B′, eligible inputs:** ΣC in 25..49 and 50..99: median lr within ±0.3 and
  median |lr| ≤ 0.4; cumulative loss 0..99, 2β × predicted ΣC against the measured
  sign-symmetric loss 2βΣC_sym: median lr within ±0.3 and median |lr| ≤ 0.4. PASS if all three,
  PARTIAL if two, else FAIL.
- **P1 (decay rate, level A, unchanged from linked_decay):** for `:derived`, median lr within
  ±0.2 and median |lr| ≤ 0.3 in each of 25..49, 50..99, 100..199; and `:derived`'s median |lr|
  below both comparators' (`:nosel`, `:spectrum`) in 50..99 and in 100..199. PASS if all;
  PARTIAL if the accuracy conditions hold in 2 of 3 windows (comparator condition still
  required); else FAIL.
- Fewer than 6 eligible G inputs: (b) and P1 INCONCLUSIVE.

Reported, not verdicts: development results on D and F; window 100..199 of (a) and (b); the
comparators; any exploratory driven-spectrum check (only if (b) passes on G).

## Amendments

One amendment is allowed after the development results and before the G runs are analysed; it
will be recorded here with its reason. None so far.
