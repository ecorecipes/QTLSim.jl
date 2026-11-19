# Bistability below the critical rate in an individual-based model: findings

2026-11-19. Criteria `validate/bistability_criteria.md` (committed after one disclosed
pilot and before every other run; no amendment afterwards), script
`validate/bistability.jl` + `bistability_analysis.jl`, log `validate/bistability.log`
(verdicts written by the script), figures `validate/bistability_figs/*.png` (one per cell:
survival over the (N₀, V₀) grid with the four predicted boundaries), post hoc diagnostic
`validate/bistability_posthoc.jl` / `.log`. Raw caches in `validate/bistability_raw/`
(git-ignored). `src/` untouched.

Xu and Shen (2026, bioRxiv preprint, doi:10.64898/2026.08.05.743160; not peer reviewed)
predict that below the critical rate persistence depends on the initial state. The test
asks whether QTLSim's kernel (explicit loci, drift, demographic stochasticity,
Beverton–Holt regulation) shows that, and where the boundary is. Their eq. 6 has priority
over `drift_critical_rate`; the feedback is Bürger and Lynch's (1995).

## Design as run

Two arms × three cells (K, R) = (1000, 2), (1000, 1.5), (3000, 2); the two (3000, 1.5)
cells were dropped for cost before the run. **G** ("Gaussian": 50 loci, 0.05 mutations per
locus per gamete, all of Poisson(5) crossovers, ω² = 25 ≫ V ≈ 1) and **Q** (QTSim-type:
100 loci, 10⁻⁴, single crossover, the input of `test/critical.jl`). Per cell: rates 0.6 and
0.85 of k_c (`drift_critical_rate` with N_e = breeders, bisection solver), initial lags 0,
0.5, 1 × the critical lag λ*, a 5 × 5 grid N₀ = 0.02–0.9 K by V₀ = 0.1–1.0 V_eq, 30
replicates per state, and two null arms (no move; 2 k_c). 6 × 200 states × 30 = 36,000
runs; horizon T = 500 generations (G), 2000 (Q, K = 1000), 1000 (Q, K = 3000), set per cell
from the deterministic prediction's latest death. Initial states: 30 pool populations per
cell (stationary-spectrum genomes, then 500 kernel generations at a fixed optimum), N₀
individuals subsampled, allele values rescaled to V₀, optimum set Δ₀ above the mean;
realised N₀ exact, |V/V₀ − 1| ≤ 2·10⁻¹⁵, lag exact. Pool V_G/V_eq 0.91–1.03.

Predictions fixed in advance, per state: P-XS (Xu and Shen's eqs. 1–5 at the K1 and K2
mappings) and P-D (a discrete-generation recursion with QTLSim's ingredients, with N_e =
breeders and with `effective_size`). All four put the boundary at V₀ ≈ 0.1–0.55 V_eq,
almost independent of N₀.

CPU: 3.2 CPU-hours for the experiment (11,483 s user), 10 CPU-minutes for the pilot, under
1 for analysis and diagnostic. One reporting crash after the runs (axis ticks of the
figures), fixed without touching runs or rules.

## Verdicts (as logged)

| | arm Q | arm G |
|---|---|---|
| null, no move (≥ 90% alive in every state with N₀ ≥ 0.05K) | PASS (and 100% at N₀ = 0.02K) | PASS (same) |
| null, 2 k_c (≤ 10% alive in every state) | PASS | PASS |
| **B1** bistability in ≥ 75% of units (cell × rate) | **FAIL (4 of 6)** | **FAIL (2 of 6)** |
| **B2** monotonicity | PASS (0 of 720 pairs) | PASS (0 of 720) |
| **B3** P-D, N_e = breeders | PASS (374/445, 84%) | FAIL (252/448, 56%) |
| **B3** P-D, N_e = `effective_size` | PASS (388/445, 87%) | FAIL (247/448, 55%) |
| **B3** P-XS, K1 | PASS (404/445, 91%) | FAIL (245/448, 55%) |
| **B3** P-XS, K2 | PASS (374/445, 84%) | FAIL (258/448, 58%) |

## What the kernel shows

**1. Dependence on the initial state exists, in the predicted form, in the QTSim-type
arm.** At 0.85 k_c all three Q cells have states that persist reliably and states that die
reliably at the same rate (30 / 8, 31 / 8 and 25 / 30 of 75 states), and Q-3000-2 also at
0.6 k_c (65 / 4). B1 fails in arm Q because at 0.6 k_c the two K = 1000 cells persist from
(almost) every state on the grid, 74 and 71 of 75 reliably, none dying reliably. That is
what the predictions say for most of those states (their boundaries there are below the
grid's lowest V₀ at lags 0 and 0.5 λ*, and at 0.07–0.13 V_eq at lag λ*), so the
B1 failure in arm Q is a grid that does not reach the extinction basin at the low rate,
not an absence of basins. Survival is monotone in N₀ and V₀ everywhere (no violation in
1440 pairs).

**2. The boundary is in V₀ (and Δ₀), not in N₀.** As predicted for these cells (fast
density regulation), the initial size hardly matters: from N₀ = 0.02K to 0.9K at fixed V₀
survival changes by +0.02 (Q) and +0.03 (G) on average, against +0.42 and +0.28 from V₀ =
0.1 to 1.0 V_eq at fixed N₀. Twenty founders with enough variance do as well as 900.

**3. In arm Q the deterministic boundaries are right to within a factor of two, and the
kernel's extinction basin is SMALLER than predicted.** Majority outcomes are matched in
84–91% of states. The logistic regressions have clear positive slopes (2.5 per ln V₀, 95%
intervals 2.4–2.6), and the kernel's half-survival variance is 0.56 (P-D breeders, P-XS
K2), 0.66 (P-D `effective_size`) and 0.77 (P-XS K1) of the predicted boundary variance:
populations persist from less initial variance than any of the deterministic models
allows. The direction is the opposite of the one Xu and Shen expect from stochasticity.
Likely reasons, not tested: a rescaled start has the equilibrium's number of segregating
alleles at smaller effects and regains variance faster than the house-of-Gauss recursion
(new mutations have the full effect α while standing effects are scaled down, and a
directional move raises the frequency of rare favourable alleles — the Barton–Turelli
rise of V_G that the package's spectrum closure has and these recursions do not); and the
boundary is sharp only for K = 3000 (slopes 4–8 per ln V₀ within a rate and lag) and
shallow for K = 1000 (1.1–3.6), where "alive at 2000" mixes fate with slow stochastic
loss. At the highest lag and rate the kernel's boundary meets the predictions (V50/V_eq
0.38–0.53 against 0.26–0.47).

**4. The two N_e variants and the two K mappings are not separated by this test.** In arm
Q P-XS K1 scores best (91%) and P-D with the correct N_e (breeders) worst together with
K2 (84%), because every prediction puts the boundary too high and the ones that put it
lowest win; that is not evidence for `effective_size` or for K1. Reading P-D at T
instead of at 10,000 generations changes no state.

**5. The "Gaussian" arm dies below the deterministic critical rate, from every initial
state.** At 0.85 k_c all 225 states of the three G cells die reliably, including the
equilibrium state itself; at 0.6 k_c all 75 states of G-3000-2 die, and in the two K =
1000 cells the kernel needs V₀ ≥ 0.73–1.04 V_eq (R = 2) or 0.44–0.66 V_eq (R = 1.5)
where the predictions say 0.16–0.26: x50 = 2.5–4.5 times the predicted boundary variance.
So B1 and B3 fail in arm G because the cells' critical rate is wrong, for all four
predictions alike (k_c from eq. 6, K1 and K2, is within −11% to +13% of the k_c used).
Survivors at 0.6 k_c hold 0.22–0.32 of V_eq and 26–44% of K, a state on its way down.

Post hoc (`bistability_posthoc.log`; exploratory, six populations per line, not a
verdict): under a move at 0.6 k_c the genic variance of the G arm falls to 0.94, 0.71,
0.45, 0.32 and 0.22 of its starting value after 100–500 generations, while without a move
it stays at 1.0, and drift alone among 900 breeders would leave 0.76 without any mutation.
Four times as many crossovers change nothing (0.90 … 0.21), so this is not linkage between
loci. With the same mutational variance and effect size spread over 500 loci instead of 50
the fall is to 0.80 (and the population keeps N ≈ 934). The loss is therefore within
loci: an allele at a locus is a non-recombining sum of many small mutations (2Nμ = 100 new
ones per locus per generation), and sustained directional selection sorts among whole
alleles faster than mutation rebuilds their spread. Neither eq. 5 of Xu and Shen nor the
house-of-Gauss recursion has a term by which the lag erodes the variance (in both, the lag
acts on the variance only through the population size), and a "Gaussian" regime built from
few loci with high mutation rates is exactly where that omission matters. The arm chosen
to be closest to the preprint's assumptions is the one farthest from its predictions; its
n_e ("effective number of recombining factors") would have to be well below 2L here. The
mechanism is a hypothesis from one diagnostic, not established.

## Answer to the question

- Persistence below the critical rate does depend on the initial state in an
  individual-based model, through the initial variance and lag, hardly through the initial
  size (in cells where density regulation is fast). Where the deterministic critical rate
  is about right (the QTSim-type arm), finite populations show both basins at one rate and
  the deterministic boundary predicts 84–91% of outcomes, with the kernel's boundary at
  0.56–0.77 of the predicted variance: the extinction basin is smaller than predicted, not
  larger.
- Where the variance erodes under the move (the 50-locus high-mutation arm), every
  prediction fails together, because the critical rate itself is overestimated by a factor
  of more than 1.2–1.7 (all states die at 0.85 k_c, and at 0.6 k_c for K = 3000): there the
  question of basins is moot at the rates tested.
- By the pre-set rules: B1 FAIL in both arms (Q 4/6, G 2/6), B2 PASS in both, B3 PASS for
  all four predictions in arm Q and FAIL for all four in arm G; all four null arms PASS.

## Limits

- Six cells, two genetic architectures, K ≤ 3000, R ≥ 1.5; the two K = 3000, R = 1.5 cells
  were not run. Slow density regulation (the preprint's r = 0.2), where N₀ should matter,
  was not tested.
- The rates are fractions of one particular k_c (N_e = breeders); in arm G that k_c is too
  high, so no rate below the kernel's own critical rate was tested there at K = 3000.
- Rescaled starts are not states a population reaches by itself (criteria file): low-V₀
  states have many alleles of small effect. A bottlenecked start (few alleles of full
  effect) could behave differently, and is the more natural test of the N₀ axis.
- Horizons are finite and differ between cells; in the Q arm the pilot's deaths had not
  levelled off by 1500 generations, so survival at T understates eventual extinction near
  the boundary (which would move the kernel's boundary up, toward the predictions).
- P-D is a recursion written for this test (house of Gauss per gene, no LD, no change of
  variance with the lag); the third prediction of the proposal (`run_moments` from the
  same state) was not made because `run_moments` cannot start at an arbitrary size.
- 30 replicates per state: a survival fraction has a standard error up to 0.09.
- The pilot (one cell, own seeds, 10 runs per state, times of death only) set the
  horizons; it is disclosed in the criteria file.
