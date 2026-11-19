# An individual-based test of bistability below the critical rate: criteria

Committed before any non-pilot run of `validate/bistability.jl`. The design refines the
proposal in `validate/xu_shen.md`, section 4.

Xu and Shen (2026, bioRxiv preprint, doi:10.64898/2026.08.05.743160; not peer reviewed)
predict that below the critical rate of a moving optimum a viable state and extinction are
both stable, so that persistence depends on the initial state (N₀, V₀, Δ₀): populations
that start on the wrong side of a basin boundary go extinct. Their model is deterministic,
in continuous time, with a Gaussian trait; they state that stochasticity should lower the
critical rate and enlarge the extinction basin. The feedback is Bürger and Lynch's (1995);
the preprint's eq. 6 has priority over QTLSim's `drift_critical_rate`.

**Question.** Does that dependence on the initial state hold in an individual-based model
with explicit loci, drift and demographic stochasticity (QTLSim's kernel under
Beverton–Holt regulation), and where is the boundary relative to deterministic
predictions?

## Model and cells

Kernel, one diploid trait, viability selection, `Demography = :beverton_holt`,
`MoveType = :continuous` from the first generation of a run.

- **Gaussian arm (G)**, close to the preprint's assumptions: 50 loci, 0.05 mutations per
  locus per gamete, `CrossoverModel = :all` with Poisson(5) crossovers, V_E = 0.5,
  m = 0.002 (V_m = 0.001), ω = 5 (ω² = 25 ≫ V ≈ 1). α²/(2V_s) = 4·10⁻⁶ ≪ 1/(2N).
- **QTSim-type arm (Q)**: published-type genetics, 100 loci, 10⁻⁴ mutations per locus per
  gamete, QTSim's single crossover (Poisson(3)), V_P = 0.03, h² = 0.4, m = 0.005, ω = 4
  (the input of `test/critical.jl`).
- (K, R): (1000, 2), (1000, 1.5), (3000, 2) in each arm: 6 cells ((3000, 1.5) was planned
  and dropped for cost before the run, see "Budget and trimming"). K = 1000 is the published
  size; the proposal's 2000 and 5000 are replaced by 1000 and 3000 for cost.

Per cell, with k_c and λ* the critical rate and critical lag of
`drift_critical_rate(p; ne = :breeders)` (N_e = the breeders, which stage 1 found to be
the effective size of this life cycle; the default formula gives a k_c 4–10% higher in
these cells) and (V_eq, N_eq) its steady state without a move:

- rates k = 0.6 k_c and 0.85 k_c;
- initial lags Δ₀ = 0, 0.5 λ*, λ* (the optimum is set Δ₀ above the mean genotype; the
  first generation's move adds k before the first selection, in the kernel and in P-D);
- initial sizes N₀ = 0.02, 0.05, 0.15, 0.4, 0.9 × K; initial variances V₀ = 0.1, 0.2, 0.35,
  0.6, 1.0 × V_eq (the grid was placed with the deterministic predictions below, which put
  their boundaries at V₀ ≈ 0.1–0.6 V_eq, nearly independent of N₀; no kernel run was used);
- 30 replicates per state; outcome: alive T generations after the start (T per cell: see
  "Pilot and horizon");
- null arms, on the same 25 (N₀, V₀) states at Δ₀ = 0: no move (k = 0), and k = 2 k_c.

2 × 3 × 25 + 2 × 25 = 200 states and 6000 runs per cell.

**Initial states** (a helper in the script, `start_world`; nothing in `src/`). For each
cell 30 pool populations are made, one per replicate index: genomes drawn from the
stationary spectrum (`equilibrium_genomes` of the fixed-size copy of the cell's
parameters), then 500 generations of the kernel at a fixed optimum under the cell's own
Beverton–Holt demography, which builds the linkage disequilibrium and genealogy of the
equilibrium. Replicate r of every state draws N₀ individuals without replacement from pool
population r, multiplies every allele value by √(V₀/V) (V the genetic variance of those N₀
individuals), and sets the optimum to their mean genotype plus Δ₀. Environmental
deviations are drawn afresh. So N₀, V₀ and Δ₀ are met exactly (checked in every run and
reported), genic variance and LD are scaled by the same factor (their ratio is that of the
equilibrium sample), and the 30 replicates of a state are 30 different populations, each
shared across states. What rescaling implies, disclosed: segregating effects are √(V₀/V)
times their equilibrium size while new mutations keep the scale α, so a start with V₀ <
V_eq has the equilibrium's number of segregating alleles with smaller effects, not fewer
alleles (as a bottleneck would leave), and regains variance from mutation at the normal
rate; in the Q arm (house-of-cards-like, variance in few alleles of large effect) such a
state is more Gaussian than any state the population would reach by itself. Subsampling
to a small N₀ loses rare alleles before the rescaling.

## Predictions (fixed here; computed by the script before it looks at survival)

For every state, persistence under four deterministic models:

- **P-XS (K1), P-XS (K2)**: Xu and Shen's eqs. 1–5 (`xu_shen_trajectory`) at the mapped
  parameters of `validate/xu_shen_targets.md`, part B (v = k, ω² ↔ ω² + V_E, V_M ↔ m V_E,
  n_e ↔ 2L, r ↔ ln R + ln c₀ − ½ln(1 + V_E/ω²), K ↔ 2K/ln R (K1) or 2K/(R − 1) (K2)), from
  (Δ₀, V₀, N = 2N₀ gene copies); extinct if k is at or above their own critical rate.
- **P-D (N_e = breeders), P-D (N_e = `effective_size`)**: a discrete-generation recursion
  with QTLSim's ingredients, in the kernel's order: λ ← λ + k; W̄ = c₀ω/√s·exp(−λ²/(2s)),
  s = ω² + V + V_E; breeders n = N W̄ (extinct if n < 2); λ ← λ(1 − V/s); V ← V + V_m −
  V²/(2L V_s) − V(α²/(2V_s) + 1/(2N_e)) (the house-of-Gauss recursion per gene, summed);
  N ← min(K, R n/(1 + (R − 1)n/K)). Its steady states are those of `drift_critical_rate`.
  Persistence = alive after 10,000 generations; read at T as a secondary line.
  (`run_moments` cannot be started at an arbitrary size, so the third prediction of the
  proposal is not made.)

For each prediction, state's N₀ and Δ₀, the **boundary** V_b is the smallest V₀ in
[V_eq/1000, 3V_eq] from which it persists (bisection on ln V₀; "none" if it persists from
V_eq/1000 or not from 3V_eq).

## Statistics and rules

A state "persists reliably" with ≥ 27 of 30 alive (≥ 90%), "dies reliably" with ≤ 3
(≤ 10%). Its majority outcome is alive if > 15 survive, dead if < 15; ties are excluded.

- **Null arms** (per arm, over its cells): no move — every state with N₀ ≥ 0.05K persists
  reliably (states with N₀ = 0.02K, 20 or 60 individuals, are reported, not judged);
  k = 2k_c — every state dies reliably. A failed null arm is reported with the verdicts and
  qualifies them; it does not change them.
- **B1, bistability exists.** Unit: cell × rate (2 per cell). A unit is bistable if among
  its 75 states (N₀ × V₀ × Δ₀) at least one persists reliably and at least one dies
  reliably. PASS for an arm if ≥ 75% of its units are bistable; INCONCLUSIVE if fewer than
  4 units were run.
- **B2, monotonicity.** Over all pairs of neighbouring states that differ in one step of
  N₀ or of V₀ (same cell, rate, lag): a violation is a fall of the survival fraction by
  more than 0.25 toward the larger value (the standard error of a difference of two
  fractions of 30 is at most 0.13). PASS for an arm if violations are ≤ 5% of its pairs.
- **B3, the deterministic boundary.** For each arm and each of the four predictions: the
  share of non-tied main-grid states whose majority outcome equals the prediction. PASS if
  ≥ 80%; INCONCLUSIVE if fewer than 50 states. Reported per cell as well.
- **Location and sharpness.** For each arm and prediction, a logistic regression of the
  survival counts on d = ln(V₀/V_b) over the states where V_b exists: the slope with a
  Wald 95% interval (a slope whose interval excludes 0 from above = a boundary in the
  predicted direction) and x50 = exp(−intercept/slope), the kernel's half-survival
  variance as a multiple of the predicted boundary. If survival is perfectly separated in
  d, the slope is reported as "separated" and x50 as the middle of the gap; with fewer than
  10 states or no variation the regression is undefined. x50 > 1 means the extinction basin
  is larger in the kernel than predicted, which is what Xu and Shen expect of
  stochasticity: it is reported as a direction, not as a failure of a prediction. Also per
  cell × rate × lag: the V₀/V_eq at which half survive, beside the four V_b/V_eq.
- Descriptive: how much survival changes with N₀ at fixed V₀ (the predictions' boundaries
  are nearly flat in N₀ in these cells); survivors' N, V_G and lag at T.

## Pilot and horizon (disclosed)

One pilot, `validate/bistability.jl pilot`, with its own seeds (80·10⁹ + …; the experiment
uses 81·10⁹ + …): the cell Q-K1000-R2, all 200 states × 10 runs to 1500 generations, and
timing runs in three other cells. It printed only the helper's errors, the pool's V_G, the
cost, and the distribution of the times of death pooled over the main grid and per null
arm, never survival by state. What it showed:

- helper: N₀ exact, |V/V₀ − 1| ≤ 1.2·10⁻¹⁵, lag exact; pool V_G 0.1235 (sd 0.022) against
  V_eq 0.1307, N 979;
- main grid: 167 deaths in 1500 runs, at times spread over the whole 1500 generations
  (median 688, 90% by 1267; 118 after generation 500, 43 after 1000); no-move null: no
  death in 250 runs; 2k_c null: 250 of 250 dead, median 268, latest 555;
- cost: 534 thread-seconds for the 2000 pilot runs; 0.30 s (G, K = 1000), 0.93 s (G, K =
  3000), 0.37 s (Q, K = 3000) per 500 generations of a surviving run.

So a horizon of 500 generations, as proposed, is too short for the Q arm. The horizon is
set per cell by a rule that uses only the deterministic prediction: **T = the smallest of
500, 1000, 2000 that exceeds the latest death of P-D (N_e = breeders) among the cell's 150
main-grid states** (computed before the run: G-1000-2 300, G-1000-1.5 441, G-3000-2 217,
Q-1000-2 1756, Q-1000-1.5 1472, Q-3000-2 726 generations) — T = 500 for the three G
cells, 2000 for the two Q cells with K = 1000, 1000 for Q-3000-2. The pilot's timing of
deaths was used to see that 500 is too short and for nothing else; the predictions and
the grid were fixed before it. Disclosed limit: the pilot's deaths did not level off by
1500 generations, so in the Q arm "alive at T" mixes the basin's fate with slow stochastic
extinction near the critical rate; survival at T is what is analysed, as fixed here.

## Budget and trimming

About 6 CPU-hours. With these horizons the pilot's costs put the eight cells at about 8.5
CPU-hours (Q-1000 36 min each, G-1000 24, Q-3000-2 60, G-3000-2 75, Q-3000-1.5 120,
G-3000-1.5 150). **Dropped before the run, by cost alone: Q-3000-1.5 and G-3000-1.5**
(the last two of the order below). The other six are run in this order and all are
finished: Q-1000-2, G-1000-2, Q-1000-1.5, G-1000-1.5, Q-3000-2, G-3000-2 (about 4.3
CPU-hours expected). Each arm keeps three cells, 6 units for B1.

Seeds: pool population r of cell i: 81·10⁹ + 10⁸ i + r; run r of state s: 81·10⁹ + 10⁸ i +
10⁴ s + 100 + r.
