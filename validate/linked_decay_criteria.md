# Linked decay: pre-specified criteria

Committed before any comparison of prediction and kernel (2026-11-14). Script
`validate/linked_decay.jl`, theory `validate/linked_decay_theory.jl`, algebra
`validate/symbolic/check_linked_decay.jl` (76 checks, PASS before this file was committed).

## The derived prediction

Per locus l and class s of partner loci k (16 position bins × 7 split-probability classes), at the
adult stage, chromosome level:

- T_s = κ(B_s, x_l, x_l), C_l = Σ_s T_s = Cov(g_rest, d_l²) (cis); the genic loss through the
  association is 2β Σ_l C_l (genic_decline's other3; its sign-symmetric part is sym3);
- y_s = Cov(B_s, x_l) (allele–background LD), Y_l = Σ_s y_s; K_s = κ(B_s, B, x_l, x_l).

One generation, to first order in β and 1/S_c (S_c = ω² + V_E + Var(G_chr)):

    T_s ← f₃(1 − ρ_s)[T_s + β(K_s + K3_s) − (c_s(κ3_l + C_l) + 2 v_l T_s)/S_c]
    y_s ← f₁(1 − ρ_s)[y_s − (c_s(v_l + Y_l) + v_l y_s)/S_c] + drift noise (variance v_l Var(B_s)/(2n)) + mutation noise
    E[K_s] = −2 E[y_s Y_l]   (leading order in 1/n)

with f₁ ≈ 1 − 1/(2n), f₃ ≈ 1 − 3/(2n) (n survivors), ρ_s QTSim's split probability (single
crossover), c_s = Cov(B_s, G_chr). So the association of class s decays at
ρ_s + 3/(2n) + (its share of) Cov(B, G)/S_c per generation and is fed by β times the drift-seeded
squared allele–background LD, which itself decays at ρ_s + 1/(2n) + Cov(B, G)/S_c. With
δ_y = δ_T this reaches Barton's (1995) steady state; it starts from zero at the move.

## Inputs, runs, seeds

- **Development (D):** genic_decline's 13 discovery inputs (set D of
  `validate/spectrum_diagnosis_inputs.csv`), R = 3 runs each, seeds KSEED + 100j + r (r = 1..3),
  i.e. the same runs as genic_decline's runs 1–3 (instrumentation draws nothing).
- **Confirmation (F):** 16 fresh inputs, 2 per stratum (NoQTL × u), drawn by `fresh_rows()`
  (MersenneTwister(20261121)) from the 118 inputs of that file with no set (the E6′ inputs never
  used by genic_decline or the diagnosis's D/C sets); R = 3 runs each, fresh seeds
  80261114000 + 100j + r. Disclosure: these inputs' kernel survival (E6′) was measured in earlier
  studies with other seeds; their association statistics were never looked at.
- Every run: `init = :equilibrium`, 1000 generations of burn-in, instrumented from MOVE − 60 to
  MOVE + 199.
- Smoke test of the analysis code (crash check only, output discarded): one short run of a Dctl
  input, which is in neither D nor F.

## Statistics (per input: mean over the runs alive at each generation, then over the window)

- ΣC = Σ_l C_l, measured and predicted, in windows 25..49, 50..99, 100..199 after the move.
- Level A: the T recursion from the measured T at t = 0, driven by the run's measured K + K3, κ3,
  β, S_c, n, N and variances (bin source split over ρ classes in proportion to n_s/(ρ_s + 1/(2n))).
  Comparators: `:nosel` (no 1/S_c terms), `:spectrum` (the spectrum's δ: factor
  (1 − ρ)(1 − 1/(2n) − V_G/S), no coupling).
- Level B: K from the y recursion, started at its stationary state for the drivers averaged over
  t = −60..−41, T = 0 at t = −60; driven by the run's measured β, S_c, n, N and variances.
- Eligible inputs: measured cumulative 2βΣC over 0..99 ≤ −0.05 × genic variance before the move
  (mean over runs alive at +100).
- Log ratio lr = log(predicted / measured) (undefined if of opposite sign; counted as a failure of
  that input for the |lr| statistics, reported).

## Predictions and verdict rules (on F; D is reported the same way, as development)

- **H1 (source identity):** median over all inputs of ΣK / Σ(−2 y Y) in windows −60..−1 and
  50..99 within [0.75, 1.33]. PASS if both, else FAIL.
- **H2 (drift-seeded allele–background LD):** level B's Σ E[yY] / measured Σ yY: median ratio in
  [0.67, 1.5] and median |lr| ≤ 0.4, in windows −60..−1 and 50..99. PASS if both windows, else FAIL.
- **P1 (decay rate, level A, eligible inputs):** for `:derived`, median lr within ±0.2 and median
  |lr| ≤ 0.3 in each of 25..49, 50..99, 100..199; and `:derived`'s median |lr| below both
  comparators' in 50..99 and in 100..199. PASS if all; PARTIAL if the accuracy conditions hold in
  2 of 3 windows (comparator condition still required); else FAIL.
- **P2 (closed, level B, eligible inputs):** ΣC in 25..49 and 50..99: median lr within ±0.3 and
  median |lr| ≤ 0.4; cumulative loss 0..99, 2β × predicted ΣC against the measured sign-symmetric
  loss 2βΣC_sym: median lr within ±0.3, median |lr| ≤ 0.4. PASS if all three, PARTIAL if two, else
  FAIL.
- Fewer than 6 eligible F inputs: P1 and P2 INCONCLUSIVE.

Reported, not verdicts: development results; the cross-check of the measured 2β₂ΣC_sym with
genic_decline's sym3 on the same D runs; any `:spectrum` prototype (exploratory).

## Amendments

One amendment is allowed after the development results and before the confirmation runs are
analysed; it will be recorded here with its reason. None so far.
