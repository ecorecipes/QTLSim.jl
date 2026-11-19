# Xu and Shen (2026) in QTLSim: what reproduces, and how it relates to `drift_critical_rate`

Xu, K. and Shen, H. (2026), "Eco-evolutionary feedbacks generate bistability in population
persistence under gradual environmental change", bioRxiv preprint,
doi:10.64898/2026.08.05.743160, version 1 (9 August 2026; not peer reviewed). The preprint
gives equations and figures and, as far as we know, no code or data. `src/xu_shen.jl` is a
re-implementation from its equations 1–8 and S1–S21. Its eq. 6 has priority over QTLSim's
`drift_critical_rate`; the feedback itself is Bürger and Lynch (1995).

Files: `src/xu_shen.jl`; `test/xu_shen.jl` (422 tests); `validate/symbolic/check_xu_shen.jl`
(77 checks); `validate/xu_shen_targets.md` (figure read-offs and the mapping, committed
before the comparisons); `validate/xu_shen.jl` → `validate/xu_shen.log` and
`validate/xu_shen_figs/fig{1,2,3,4}.png` (139 s of CPU time on one thread).

## 1. What reproduces

There is no published data, so every target is a value read off a figure rendered at
300 dpi; nothing is a pass/fail verdict. Numbers are "read-off → ours".

| figure | result |
|---|---|
| 1a boundary v_c(V_M) | 0.08–0.09 → 0.087 (V_M = 0.05); ≈ 0.16 → 0.163; ≈ 0.19 → 0.190; at V_M = 0.1 the first read-off (0.115–0.12) was too low: read again ≈ 0.13 → 0.130 |
| 1a, 1b V*, N* | V* < 0.6 … > 3.0 → 0.59 … 3.02; V* = 1.8 at V_M ≈ 0.075–0.08 → 0.076; N* = 180 contour at V_M ≈ 0.085 → 0.0855; N* < 60 under the boundary → 43–52 at v_c |
| 2a blue, green, red | ends and peaks agree to ≤ 0.002 in v_c and ≤ 0.13 in the peak's position (e.g. blue 0.087 / 0.136 near −1.2 / 0.093 → 0.0862 / 0.1359 at −1.15 / 0.0926) |
| 2a yellow ("V_M = 0.5") | does **not** match V_M = 0.5 (0.14 / 0.19 / 0.095 → 0.255 / 0.304 / 0.097); post hoc it matches V_M = 0.2 (0.140 / 0.192 at −1.45 / 0.096) |
| 2b, 2c | do **not** match the caption's parameters (boundary = panel a's blue line, N* up to 196, against a drawn boundary of 0.045–0.05 at −3 and N* < 100); post hoc K = 500 is the closest of 84 combinations (0.059 at −3, above 0.12 from −1.22 to −0.71, N* ≤ 98), not exact |
| 3 equilibria | (2.06, 142) → (2.067, 141.9); (1.96, 88) → (1.968, 87.4) |
| 3 basin boundaries | all six panels: (a), (b), (d) every plotted start persists (63/63); (c) (0.5, 75), (1.0, 38), (2.0, 17), (3.0, 11) → 82.6, 35.8, 16.5, 10.9; (e) 15, 10, 7, 6 → 16.5, 10.7, 6.9, 5.4; (f) (2.0, 100), (2.5, 55), (3.0, 35) → 96.9, 53.8, 36.9 |
| 3 transients | overshoot of N at Δ₀ = 0.5 (55 and 42 starts); persisting starts that first fall in N and V in (c), (f) (23, 13); extinct starts that first rise in both in (e) (4) |
| 4a | all four curves within the reading precision (e.g. Δ₀ = 3, v = 0.05: 13 / 40 / 10³ / 10⁴ → 13.9 / 39.1 / 1194 / 1.06·10⁴), taking the fourth legend entry as Δ₀ = 3 |
| 4b | plateaus 0.1355, 0.123 → 0.1354, 0.1232; zeros 5.65, 4.95, 9.2, 7.7 → 5.65, 4.96, 9.29, 7.72; the curves leave the plateau tangentially (exactly at 2.39, 1.72, 4.41, 2.87; 1% below it at 2.67, 1.99, 4.83, 3.29, against read-offs 2.6, 1.9, 4.8, 3.3) |
| 4c | qualitative only: the rate of the drawn surface is not stated. Ours (v = 0.08): boundary V₀ ≈ 0 up to a lag of 2.5–4, then steep, steeper at small N; all 79 equilibria for v′ ∈ [0, v_c) lie on the persisting side |
| 4d | the first read-offs misplaced the v₁ axis (the crop had no tick labels) and disagree by 25–45%; read again with the ticks, **not blind**: ω² = 100: 170, 89, 55 → 176, 84, 57; ω² = 10: 50, 30 → 48, 31, crossing at 0.143 → 0.141; ω² = 1: 185, 105, 40, 22 → 191, 107, 40, 22 |

Text claims checked: with no abrupt shift a population at the equilibrium of any v₀
persists at any v₁ < v_c (196 of 196 pairs); with V₀ at balance persistence can be
impossible from every N₀ although v < v_c (Δ₀ = 3, v = 0.1, ω² ≈ 1.4–1.6: the boundary V₀
exceeds V_c^M).

**The extinction rule.** In the equations N never reaches 0 and a small N can grow back, so
a rule is needed and the preprint does not state one. Ours is N < 10⁻³. The boundaries of
Figures 3 and 4a are identical to 4–6 digits for thresholds from 1 to 10⁻⁸, because the
boundary is the stable manifold of the saddle equilibrium (N = 11.5 and 28.0 at v = 0.08 and
0.12), not the threshold. Where the equations' own critical size is below one individual
(e.g. 0.014 at v = 0.08, Δ₀ = 0.5) Figure 4a's axis floor of 1 hides it.

## 2. Notes on the preprint's equations

Precise and neutral; none affects a result. Details and checks in
`validate/symbolic/check_xu_shen.log`.

1. **S17.** Printed: P₁ = J₂₃(J₁₁J₃₂ − J₃₁J₁₂) − J₃₃(J₁₁ + J₂₂)². Expanding a₁a₂ − a₃ with S16
   (and S18 as printed) gives −J₂₃(J₁₁J₃₂ − J₃₁J₁₂) − J₃₃(J₁₁ + J₂₂)². S19 follows the printed
   S17, so its terms V*²/(2N*ω⁴) and −v²/(N*V*) carry the opposite signs to the derived ones.
   S20–S21 are correct algebra given S19, but S20 bounds v² from above, which bounds the
   derived P₁ from above, so they do not establish P₁ > 0 as written (the derived P₁ is
   negative at some stable equilibria of random parameter sets). **The conclusion stands**:
   a₁a₂ − a₃ = J₂₃J₃₁J₁₂ + J₂₂J₂₃J₃₂ + J₃₃J₂₃J₃₂ − (J₁₁+J₂₂)(J₁₁J₂₂ + J₃₃(J₁₁+J₂₂+J₃₃)), every term
   positive under the signs of S8, so stability is decided by det J alone (S15).
2. Line 410: "(J₁₃, J₂₁, J₃₃) are zero"; J₃₃ = −N*/K is given two lines above and in S8.
3. Text after S5: the decreasing function is written with ½·V/(2ω²); S5 has V/(2ω²).
4. S1 calls extinction (N → 0, V/N → V_M, Δ → ∞) an equilibrium; it is an attracting state.
5. Figure 4a: the fourth legend entry repeats the second ("Δ₀ = 1, v = 0.1"); the curve
   matches Δ₀ = 3, v = 0.1.
6. Figure 2a, yellow curve, and Figure 2b, c: the drawn values do not match the stated
   parameters (section 1); V_M = 0.2 and K = 500 are our post hoc guesses, not established.

Everything else in eqs. 1–8 and S1–S16 was confirmed symbolically.

## 3. The mapping to `drift_critical_rate`, and its limits

Time unit = one generation, v = k; ω² ↔ V_s = ω² + V_E; V_M ↔ m·V_E; n_e ↔ Ploidy·L (1 for a
haploid clone); N ↔ Ploidy·N_e gene copies; r ↔ ln R + ln c₀ − ½ ln(1 + V_E/ω²); K ↔
Ploidy·K/ln R (K1, the same gene copies at no load) or Ploidy·K/(R − 1) (K2, the same
low-density slope). Summing the house of Gauss over genes gives eq. 5 at equilibrium plus a
third loss term V·α²/(2V_s).

| inputs | v_c(eq. 6)/k_c (K1) | with K2 |
|---|---|---|
| Xu–Shen-type, n_e = 1 (haploid clone), 8 parameter sets | 1.011–1.026 | 0.980–1.009 |
| Xu–Shen-type, n_e = 100 (diploid, 50 loci), 8 sets | 1.017–1.096 | 0.948–0.998 |
| published-type, Beverton–Holt R = 1.5 / 2 / 5 (6 inputs) | 1.006–1.169 / 1.018–1.184 / 1.090–1.267 | 0.82–0.96 / 0.71–0.84 / 0.44–0.53 |

Attribution (a ladder changing one ingredient at a time; the split depends on the order):

- the variance equation (dropping α²/(2V_s)): ×1.001–1.004 where α²/(2V_s) ≪ 1/(2N), up to
  ×1.16 for the test input with 100 loci (α²/(2V_s) = 0.28/(2N));
- drift from the census of adults instead of N_e from the breeders: ×1.07 (R = 1.5), ×1.14
  (R = 2), ×1.39 (R = 5);
- first-order (continuous-time) response and loads: ×1.000–1.004 for published-type inputs
  (V ≪ ω²), up to ×1.09 at r = 0.5;
- logistic instead of Beverton–Holt regulation (K1): ×0.94, 0.89, 0.78.

The last rung equals eq. 6 at the mapped parameters to 4 digits, which checks the mapping.

Limits. (1) Weak selection (ω² ≫ V): harmless for published inputs. (2) Overlapping against
discrete generations: only first-order agreement. (3) Density regulation cannot be matched
for all loads; K1 and K2 differ by up to a factor 2.5 at R = 5. (4) **A fixed size has no
counterpart**: N has no dynamics and the feedback acts through N_e(N, N W̄) alone; there
k_c is 1.8 times the Beverton–Holt (R = 2) value for the same input. (5) Linkage: neither
eq. 5 nor the house of Gauss has the LD of QTSim's single crossover.

Two observations about QTLSim's own functions made on the way (nothing changed):

- `drift_critical_rate` follows one branch of steady states (its inner iteration starts from
  the variance at N = K). Where the variance of a large population carries a load near or
  above the growth rate it returns less than the maximum over all steady states, or 0
  (3 of 16 Xu–Shen-type inputs with n_e = 100: 0 against 0.297, 0.296 against 0.434, 0.411
  against 0.421). Published-type inputs are far from that regime (0 of 18 differ).
- `effective_size(N, n)` is close to the harmonic mean 2Nn/(N + n). A plain neutral
  simulation of the same sampling scheme (N zygotes, n breeders, two distinct parents; not
  the kernel) gave N_e = n: 198, 100, 51 for n = 200, 100, 50 of N = 200, against 200, 134,
  81 from the formula. With N_e = n, k_c is 0.90–0.91 times its value for published-type
  inputs at R = 2. Not verified in the kernel; worth a separate check.

**The n_e that eq. 5 would need** for the equilibrium variance at v = 0 (a rough estimate:
close to the neutral 2N_eV_M the denominator is a small difference): from the house of
Gauss 45–2850, i.e. 0.23–0.29 of 2L; from `multilocus_equilibrium` (rare alleles and LD
under the single crossover) 39–159; from the published table's pre-move V_G (60 inputs per
stratum, one QTSim run each, 13–27% undefined) medians 14–21 (100 loci), 61–85 (500), 60–62
(1000), 102–151 (5000), i.e. 0.10 down to 0.01 of 2L. So n_e of order 20–150, far below the
number of loci and hardly growing with it; no single n_e represents the house-of-cards term.

## 4. Proposal (not executed): an individual-based test of bistability in the kernel

**Question.** Do finite populations in the kernel show the two basins: at one rate below
the critical rate, persistence from some initial states and extinction from others, with
the switch where the deterministic boundary is?

**Model cells.** `Demography = :beverton_holt` (a fixed size has no N dynamics), one
diploid trait, viability selection, `MoveType = :continuous` from generation 0,
`CrossoverModel = :all` and 50 loci at μ = 0.01–0.1 in a first, "Gaussian" arm (closest to
eq. 5), then QTSim's own regime (100 loci, μ = 10⁻⁴, single crossover) in a second arm.
Sizes must be larger than the preprint's standard set, whose equilibrium is a few tens of
diploid individuals after mapping (demographic noise would blur everything): K = 2000 and
5000, R = 1.5 and 2. Two rates per cell, 0.6 and 0.85 of the cell's critical rate; three
initial lags (0.5, 1 and 1.3 times the saddle's lag).

**Initial states.** N₀ and V₀ on a 6 × 6 grid spanning the predicted boundary (N₀ from 0.1
to 1 of the equilibrium size, V₀ from 0.2 to 1.2 of the equilibrium variance). To start the
kernel there: `World(p)` with capacity K; genomes from `equilibrium_genomes` of the
matching fixed-size `Params` at N₀ (the sampler rejects Beverton–Holt, so it is called on a
fixed-size copy), allele values scaled by √(V₀/V) and written with `setgenomes!`; the
current size set to N₀ (`w.N`, an internal field: `setgenomes!` leaves N unchanged, so a
small helper is needed and must be tested); `develop!`; the optimum `w.θ` set to the mean
genotype plus Δ₀; `init_env!`, then `run!(w, rng; resume = true)`. `equilibrium_world`
itself cannot be used (fixed size only, lag 0). Scaling allele values keeps the spectrum's
shape but not its relation to α², and the sampler has no LD: 20–50 generations at a fixed
optimum and fixed size before the move would restore LD at the cost of letting V drift from
V₀; record the realised (N, V, Δ) at the move and analyse against those.

**Predictions, fixed before the runs.** Three, each for every cell: (P1) the Xu–Shen
boundary at the mapped parameters (K1), (P2) the same from a discrete-generation version of
the three recursions with QTLSim's own ingredients (to be written and committed first),
(P3) `run_moments(…; erosion = true)` started from the same state where it can be. Outcome:
alive at 500 generations (the deterministic fates are settled within about 300), 40
replicates per state. Pre-specified statistics: the share of states called correctly
(alive in > 50% of replicates against the predicted fate), a logistic fit of survival
against the signed log distance to the predicted boundary (location within a factor 1.5 in
N₀, and a slope; a slope near 0 would mean no sharp basins), and, as the direct sign of
bistability, at least one pair of states in each cell with ≥ 90% and ≤ 10% survival at the
same rate. A null arm at 1.2 times the critical rate (all should die) and one at Δ₀ = 0 and
low rate (all should persist). Known in advance: stochasticity removes true bistability
(every finite population dies eventually), so the horizon is finite and the boundary is a
band; the preprint itself expects a lower critical rate and a larger extinction basin with
stochasticity.

**Cost.** 8 cells (2 K × 2 R × 2 arms) × 2 rates × 3 lags × 36 states × 40 replicates =
69,120 runs of ≤ 500 generations at 0.25–0.5 ms per generation (N ≤ 5000 scales it up to
about 1.5 ms): 3–6 CPU-hours, plus 0.5 h for the null arms and the helper's tests. A pilot
on one cell (4,320 runs, about 20 CPU-minutes) should come first, to check the helper and
the horizon, with its seeds excluded from the main run.
