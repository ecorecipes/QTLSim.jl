# Xu and Shen (2026): comparison targets and parameter mapping, written before the comparisons

This file is committed before `validate/xu_shen.jl` computes anything it is compared with;
the script prints it at the top of `validate/xu_shen.log`.

Source: Xu, K. and Shen, H. (2026), "Eco-evolutionary feedbacks generate bistability in
population persistence under gradual environmental change", bioRxiv preprint,
doi:10.64898/2026.08.05.743160, version 1 (posted 9 August 2026; not peer reviewed). The
preprint publishes equations and figures and, as far as we know, no code or data. `src/xu_shen.jl`
is a re-implementation from its equations 1–8 and S1–S21.

**Disclosure.** While the module was being written, and before this file was committed, the following
numbers had been computed with it and seen: the critical rate at ω² = 20, 50, 100 (0.1354,
0.1302, 0.1232) and at V_M = 0.05, 0.15, 0.2 (0.0871, 0.1631, 0.1904), and the equilibria at
v = 0, 0.08, 0.12 for the standard parameters ((V, N) = (2.10, 179), (2.07, 142), (1.97, 87)),
and a handful of basin-boundary values at (v, Δ₀) = (0.08, 3), (0.08, 5), (0.12, 3), (0.12, 5).
The read-offs below are the ones taken from the figures (rendered at 300 dpi) before those
computations and are left as they were read, including where they later looked imprecise
(A1 at V_M = 0.1). Everything else was computed after this file was committed.

## A. Read-offs from the preprint's figures (targets)

Precision is that of reading a rendered figure: about ±0.005 in a rate, ±0.05 in V, ±3 in N on
Figures 1–3 unless a range is given; "≈" marks values read between ticks.

**Figure 1** (standard parameters, V_M on the x axis, v on the y axis; boundary of the coloured
area = v_c).
- A1. boundary v_c: at V_M = 0.05 ≈ 0.08–0.09; at V_M = 0.10 ≈ 0.115–0.12; at V_M = 0.15 ≈ 0.16;
  at V_M = 0.20 ≈ 0.19.
- A2. V* is nearly independent of v (near-vertical contours) and runs from < 0.6 at the left
  edge (V_M ≈ 0.01) to > 3.0 at V_M = 0.2 and low v; the contour V* = 1.8 is near V_M ≈ 0.075–0.08.
- A3. N* at low v is in the top class (> 180) only for V_M below ≈ 0.085 at v → 0.01, and
  160–180 across the rest of the bottom edge; N* falls with v to < 60 (darkest class 40–60)
  just under the boundary at V_M ≥ 0.1.

**Figure 2** (x axis log₁₀(1/ω²) from −3 to 0).
- A4. panel (a), blue (standard r = 0.2, V_M = 0.1, K = 10³): v_c ≈ 0.087 at −3, peak ≈ 0.136
  near −1.2 (between −1.4 and −1.0), ≈ 0.093 at 0.
- A5. yellow (V_M = 0.5): ≈ 0.14 at −3, peak ≈ 0.19 near −1.5, ≈ 0.095 at 0.
- A6. green (K = 10⁴): ≈ 0.157 at −3, peak ≈ 0.17 near −1.9, ≈ 0.097 at 0.
- A7. red (r = 0.5): ≈ 0.19 at −3, peak ≈ 0.26 near −0.9, ≈ 0.24 at 0.
- A8. panels (b), (c) (caption: r = 0.2, V_M = 0.1, K = 1000, n_e = 1; v axis clipped at 0.12):
  the boundary is at v ≈ 0.045–0.05 at −3, reaches the top of the frame (0.12) between ≈ −1.2
  and ≈ −0.8, and is ≈ 0.075–0.08 at 0; V* > 6 at the far left and < 0.5–1 right of ≈ −0.5;
  N* is at most in the class "> 90" (bottom left), i.e. below about 100 everywhere.
  *Noted before computing:* these do not look like the same parameters as the blue line of
  panel (a) (0.087 at −3, peak 0.136 > 0.12 over a wider range) or as Figure 1b (N* up to
  180 at ω² = 50); to be examined.

**Figure 3** (standard parameters; V on x, N on y).
- A9. viable equilibrium (red dot): v = 0.08: (V*, N*) ≈ (2.06, 142); v = 0.12: ≈ (1.96, 88).
- A10. panels (a), (d) (Δ₀ = 0.5): every plotted start persists (V₀ from 0.2 to 3.0, N₀ from 10
  to 150); the boundary is drawn along the axes.
- A11. panel (b) (v = 0.08, Δ₀ = 3): every plotted start persists; boundary roughly horizontal
  at N₀ ≈ 2–5 for V₀ from ≈ 0.1 to 3.
- A12. panel (c) (v = 0.08, Δ₀ = 5): boundary through about (V₀, N₀) = (0.3, 180), (0.5, 75),
  (1.0, 38), (1.5, 25), (2.0, 17), (3.0, 11) (N₀ ± 5, more at the steep left end).
- A13. panel (e) (v = 0.12, Δ₀ = 3): boundary near-vertical at V₀ ≈ 0.02–0.05 above N₀ ≈ 50,
  then about (0.5, 15), (1.0, 10), (2.0, 7), (3.0, 6) (N₀ ± 3).
- A14. panel (f) (v = 0.12, Δ₀ = 5): boundary through about (1.7, 180), (2.0, 100), (2.5, 55),
  (3.0, 35) (N₀ ± 8).
- A15. qualitative: with Δ₀ = 0.5 trajectories overshoot in N (rise above the equilibrium, then
  fall); in (c) and (f) some persisting trajectories first fall in N and V; in (e) some
  trajectories that go extinct first rise in both N and V.

**Figure 4.**
- A16. panel (a) (critical N₀ with V₀ at mutation–selection–drift balance, floor of the axis at
  1): black solid (Δ₀ = 1, v = 0.05): 1 up to ≈ −0.65, ≈ 3 at 0. Red solid (Δ₀ = 1, v = 0.1):
  drawn from ≈ −2.6 to ≈ −0.1, leaves 1 at ≈ −1.35, ≈ 2 at −1.1, ≈ 5 at −0.6, ≈ 18 at −0.1.
  Black dashed (Δ₀ = 3, v = 0.05): leaves 1 at ≈ −1.65, ≈ 13 at −1.0, ≈ 40 at −0.7, ≈ 10³ at
  −0.15, 10⁴ at 0. Red dashed: leaves 1 at ≈ −2.1, ≈ 8 at −1.6, ≈ 30 at −1.1, ≈ 10² at −0.8,
  ≈ 10³ at −0.45, off the top (10⁴) at ≈ −0.35. (Factor ≈ 1.5 on the log axis.) The fourth
  legend entry reads "Δ₀ = 1, v = 0.1", the same as the second; from the order of the entries
  and the text we take the red dashed curve to be Δ₀ = 3, v = 0.1.
- A17. panel (b) (max v₁ against the abrupt shift δ): plateaus ≈ 0.1355 (ω² = 20) and ≈ 0.123
  (ω² = 100). Curves leave the plateau at δ ≈ 2.6 (ω² = 20, v₀ = 0.04), ≈ 1.9 (20, 0.08), ≈ 4.8
  (100, 0.04), ≈ 3.3 (100, 0.08), and reach 0 at δ ≈ 5.65, ≈ 4.95, ≈ 9.2, ≈ 7.7 (±0.15).
- A18. panel (c): a surface of boundary V₀ over (Δ, N) with the curve of equilibria for v from
  0 to v_c. The rate of the surface is not stated; the axes are only partly legible. Target is
  qualitative: boundary V₀ ≈ 0 for small lags, rising steeply beyond a lag of a few units, more
  steeply at small N; the curve of equilibria stays on the persisting side and meets the
  surface as v → v_c.
- A19. panel (d) (max duration T against v₁; solid v₀ = 0.02, dashed v₀ = 0.09): the curves
  diverge at v₁ ≈ 0.095–0.10 (ω² = 1), ≈ 0.123 (ω² = 100), ≈ 0.135–0.14 (ω² = 10). Solid
  curves: ω² = 100: ≈ 230 at v₁ = 0.15, ≈ 100 at 0.20, ≈ 60 at 0.25; ω² = 10: crosses the ω² = 100
  curve near v₁ ≈ 0.148 (T ≈ 220), ≈ 60 at 0.20, ≈ 30 at 0.25; ω² = 1: ≈ 190 at 0.11, ≈ 60 at
  0.15, ≈ 20 at 0.20. Dashed curves lie below the solid ones (ω² = 100: ≈ 180 at 0.15). (±20,
  more near the divergence.)

## B. Mapping to QTLSim's discrete-generation counterpart (`drift_critical_rate`)

Xu and Shen's eq. 6 has priority; the feedback itself is Bürger and Lynch (1995).
`drift_critical_rate` is a discrete-generation implementation matched to QTSim's life cycle.
The comparison is descriptive: no pass/fail.

QTLSim, per generation: N zygotes; viability W(z) = c₀ exp(−(z − θ)²/(2ω²)), c₀ = 1/(0.4√(2π))
≈ 0.9974; n = N W̄ breeders with W̄ = c₀ ω/√(ω² + V_P) · exp(−λ²/(2(ω² + V_P))); under
`Demography = :beverton_holt` Poisson(R n/(1 + (R − 1)n/K)) offspring, equilibrium breeders
n = K(RW̄ − 1)/(R − 1) and adults n/W̄; response V_G λ/(V_G + V_s), V_s = ω² + V_E;
N_e = `effective_size(N, n)`; V_G the house of Gauss: per gene c²/V_s + c(α²/(2V_s) + 1/(Ploidy·N_e)) = μα².

| Xu–Shen | QTLSim | status |
|---|---|---|
| time t, rate v per unit time | one generation, k = `MoveFitnessSD` (continuous move) | discrete vs continuous: agree to first order in the per-generation changes |
| ω² (trait fully heritable) | V_s = ω² + V_E | exact for selection on genotypes; the response V_Gλ/(V_G + V_s) equals VΔ/ω² only for V ≪ ω² |
| V (genetic = phenotypic) | V_G | – |
| V_M | m·V_E (`muteVar`), per individual per generation | exact |
| n_e | Ploidy·L under free recombination (the Gaussian equilibrium √(2L V_m V_s) = ω√(n_e V_M)); 1 for a haploid clone | only for the Gaussian part of the variance equation, see below |
| −V/N (haploid, census N) | −V/(Ploidy·N_e), N_e from the breeders | N ↔ Ploidy·N_e gene copies; N_e depends on the breeders n = N W̄, not on the census alone: not matchable exactly |
| r | ln R + ln c₀ − ½ ln(1 + V_E/ω²) | the growth rate at low density of a perfectly adapted population without genetic variance; the remaining load ½ ln(1 + V_G/V_s) + λ²/(2(V_s + V_G)) equals (V + Δ²)/(2ω²) only for V ≪ ω² |
| −N/K (logistic) | −ln(1 + (R − 1)n/K) (Beverton–Holt on breeders) | not matchable in general. Main choice (K1): K_XS = Ploidy·K/ln R, the same number of gene copies at no load. Alternative (K2): K_XS = Ploidy·K/(R − 1), the same slope at low density. They agree as R → 1 |

Where the two cannot be matched:
1. **The variance equation.** Summing the house of Gauss over Ploidy·L genes with V = Ploidy·L·c
   and V_M = Ploidy·L·μα² gives V_M = V²/(Ploidy·L·V_s) + V/(Ploidy·N_e) + V·α²/(2V_s): eq. 5 at
   equilibrium with n_e = Ploidy·L and N = Ploidy·N_e, plus a third loss term α²/(2V_s)·V (rare
   alleles of large effect, the house-of-cards limit). It has the form of extra drift, not of a
   different n_e, so no n_e reproduces it for all N. It vanishes as μ → ∞ at fixed V_M (α² → 0).
2. **Drift from breeders.** At a fixed size (`Demography = :fixed`) N has no dynamics at all; the
   feedback acts only through N_e(N, N W̄). Eq. 6 has no counterpart of that; fixed-size values
   are listed separately without a mapping.
3. **Weak selection.** Eq. 4 and the loads of eq. 1 are first order in V/ω² and Δ²/ω²; QTLSim keeps
   V_G/(V_G + V_s), ω/√(ω² + V_P) and exp(−λ²/(2(ω² + V_P))).
4. **Density regulation** (the last table row), and overlapping against discrete generations.
5. **Linkage.** Eq. 5 neglects LD; `drift_critical_rate` (house of Gauss) does too. QTSim's
   single crossover builds negative LD that neither has (`validate/spectrum_diagnosis.md`).

Planned evaluation (descriptive):
- (i) Xu–Shen-type inputs built in QTLSim's parameters: haploid clone (n_e = 1) and diploid with
  L loci (n_e = 2L), V_s = ω², V_M, r, K as in the preprint's standard set and one-at-a-time
  variations (ω² ∈ {20, 50, 100}, V_M ∈ {0.05, 0.1, 0.2}, K ∈ {500, 1000, 2000}, r ∈ {0.2, 0.5});
  mutation rate high enough that α²/(2V_s) is small. Report v_c(eq. 6), k_c(`drift_critical_rate`)
  and the ratio.
- (ii) Published-type QTSim inputs (diploid, N = 1000, L = 100 and 500) under Beverton–Holt
  R ∈ {1.5, 2, 5}: the same three numbers after mapping to (r, ω², V_M, K, n_e), plus the
  fixed-size k_c without a mapping.
- (iii) A ladder from `drift_critical_rate` to eq. 6, changing one ingredient at a time in a local
  copy of the steady-state feedback (no change to `src/critical.jl`): (1) drop the α²/(2V_s)
  term (their variance equation); (2) drift from the census of adults instead of N_e from the
  breeders; (3) first-order (continuous-time) response and loads; (4) logistic density
  regulation, which is eq. 6 at the mapped parameters. The ratios between successive rungs
  attribute the difference.
- (iv) The n_e that eq. 5 would need to reproduce a given equilibrium variance at v = 0 and
  N = Ploidy·N_e: n_e = V²/(V_s(V_M − V/(Ploidy·N_e))), for V from `house_of_gauss_vg`, from
  `multilocus_equilibrium(p; drift = true)` (rare alleles and LD under QTSim's single crossover)
  and its `loci = :multilocus` form, and from the published table's pre-move V_G of a sample of
  inputs (one QTSim realisation each; no new kernel runs). A rough estimate: where V is close to
  the neutral 2N_e V_M the denominator is a small difference and n_e is ill-conditioned.
