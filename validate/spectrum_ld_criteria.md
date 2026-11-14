# Pre-specified criteria: the `:spectrum` closure with dynamic linkage disequilibrium

Written 2026-11-14, after development (validate/spectrum_ld_dev.log, validate/spectrum_ld_dev.jl)
and before any held-out result. Committed alone. The script is validate/spectrum_ld.jl; its
output goes to validate/spectrum_ld.log. No default changes whatever the outcome.

## The model under test

`emulate(p; loci = :spectrum, spectrum = (; ld = :dynamic))` (deterministic), i.e. the moments
lowering with `loci = :spectrum` and `spectrum = (; ld = :dynamic)` (defaults `ldbins = 20`,
`turnover = true`, `linkage = :barton`), src/spectrum.jl:

- V_G = g + D, g = Σ 2a²p(1 − p) the spectrum's genic variance, D = Σ_c d_c the linkage
  disequilibrium of the classes c of pairs of distinct loci, grouped by the bins (a, b) of 20
  bins of consecutive positions that hold their two loci and by their split probability ρ_c
  under `CrossoverModel` (`QTLSim.ld_classes`), ω_c their share of the L² ordered pairs.
- Each bin b has β_b = V_b/g_b (its genic share plus the LD of its pairs; half of a pair shared
  with another bin).
- Selection on a site of effect a at frequency p, lag λ, S = ω_e² + V_P:
  E[Δp] = pq[β̄ aλ/S − β²‾ (1 − λ²/S)(a² + 2c)(1 − 2p)/(2S)], with β̄ = Σ_b m_b β_b = V_G/g and
  β²‾ = Σ_b m_b β_b² (m_b the bin's share of the loci; one spectrum for all bins). The response
  is V_G λ/S.
- One generation of LD: d_c ← (1 − 1/(2n))[(1 − ρ_c)(1 + 2γ) d_c − ω_c β_a β_b g²/(2S)], n the
  parents, γ the relative change of g by selection in the generation (0 when neutral).
- Drift, mutation and Barton's linked-selection term as in the default `:spectrum`; mutations
  enter without LD. `run_moments(…; VG)` / `emulate`'s start sets the genic variance to VG,
  the LD at the equilibrium's ratio.

Nothing is fitted to kernel or table data: every quantity is a model quantity; `ldbins = 20`
is a numerical resolution (V_G/genic converged to < 3e-4 between 20 and 80 bins on five
inputs). The default `:spectrum` output is bit for bit unchanged.

## Disclosure: what development used

- The diagnosis's discovery inputs D and Dctl (21 inputs) and their kernel trajectories
  (10 runs each): static V_G, genic variance and LD at the move; V_G during the move; the
  genic variance and LD driven by each run's own lag and survivors; k50 on the diagnosis's six
  k50 inputs (its kernel cache).
- In step 2 of the dev log, the LD fixed point given the kernel's genic variance was also
  computed and printed for the C and Cctl inputs (static LD at the move only; their kernel runs
  had already been seen by the diagnosis). Nothing was tuned on them.
- The choices made in development: position bins (derived in step 2, not fitted), the turnover
  term (derived), Barton's term kept (an LD-consistent replacement was worse, step 6),
  genic-variance start semantics (step 5).

## Samples

- **H (held out).** E6′'s procedure (validate/emulator_spectrum.jl, amendment 1) on a new pool:
  - exclusions: every input excluded by validate/emulator_spectrum.jl (`excluded_ids`), its
    E1–E4 sample, its E6′ pool (1200, containing the 160 E6′ inputs), the diagnosis's 42
    diagnostic inputs and 64 fresh inputs, the inputs of validate/published_sweep_results.csv,
    and the inputs whose table row was rerun in vignettes 09 and 11;
  - pool: 100 per stratum, MersenneTwister(20261121) permutation, prefiltered by
    0.3 ≤ MoveFitnessSD/kc ≤ 3 (kc = `drift_critical_rate`), their `:rare_alleles` k50
    (`emulated_half_survival_rate`);
  - picks: 0.3 ≤ r = MoveFitnessSD/k50 ≤ 1.5, 5 per log-bin of r (4 bins), filled to 20 per
    stratum from the rest, MersenneTwister(20261122) — 160 inputs (the `select` mode was run
    before this file was committed: 6097 inputs excluded, 20 per stratum, 5 per bin but 4 in
    one bin of QT0500-Possion5.0). No table outcome is
    used in the selection.
  - kernel: 40 runs per input, exactly E6′'s protocol (init = :equilibrium, burn-in 1000,
    generations to 12000, the move at 10000), seeds 70261114000 + 100j + r (fresh), recording
    alive at +2000, V_G (mean over −200…−1), the genic variance (mean at −200, −100, −1), and
    V_G at +100 and +500 (mean over ±5 generations, runs alive throughout the window).
- **C (confirmation set of the diagnosis, 13 inputs; reported separately, no verdict).** Their
  E6′ kernel runs (40 each) and the diagnosis's 10 trajectories per input, already seen.
- **F (fresh table sample).** 150 per stratum (1200), MersenneTwister(20261123), outside the
  exclusions above and H's pool.

## Tests (model under test; each PASS/FAIL)

Model quantities: V_G before the move is the emulator trajectory's mean over −200…−1, the
genic variance and LD at the move are `spectrum_equilibrium(p; ld = :dynamic)`'s (the
trajectory approaches it during the 2000-generation burn-in), V_G at +100/+500 the
trajectory's mean over ±5 generations.

- **(a) the move, on H.**
  - (a1) median over inputs of model/kernel V_G before the move in [0.9, 1.1];
  - (a2) median of model/kernel genic variance in [0.9, 1.1];
  - (a3) on the inputs with kernel LD/genic ≤ −0.1, median of (model LD/genic)/(kernel LD/genic)
    in [0.8, 1.25].
  - (a) passes if a1, a2 and a3 pass.
- **(b) during the move, on H.** Median of model/kernel V_G at +100, and at +500, over the
  inputs with at least 10 of 40 kernel runs alive through the window and the model's trajectory
  alive, each in [0.8, 1.25]. (b) passes if both pass. (Development showed the genic decline
  under the move is still too late; (b), especially +500, is expected to be hard.)
- **(c) the threshold, on H (E6′'s statistics).**
  - E6′c: accuracy of P2000 > 0.5 for the kernel's majority alive ≥ 0.85 and AUC ≥ 0.85;
  - x50 per stratum: ratio emulator/kernel in [0.8, 1.25] in every stratum (E6′d's rule);
  - pooled x50 ratio (all H) in [0.8, 1.25].
  - E6′a/b are reported (no verdict here).
- **(d) the table, on F.** E1–E4 of validate/emulator_criteria.md, as validate/emulator_spectrum.jl
  applies them (E3 with the predicted median T of the deterministic trajectory; E5 not computed).

## Baselines (same inputs, no verdict)

The default `:spectrum`, the `ldratio` prototype (r from `multilocus_equilibrium(p; drift = true,
loci = :multilocus)`; H and C only) and the deterministic `:rare_alleles` emulator (H and F).

## Budget and resizing

About 12 CPU-hours for the whole task, of which development used about 2.5. Estimated: H kernel
≈ 4 CPU-h, emulations on H and C ≈ 1.5, F ≈ 3.5. A timing pilot (16 inputs of F) runs first; if
the total would exceed the budget, F is shrunk (per stratum, same seed order) by an amendment
recorded here before any held-out result.

## Amendments

1. **2026-11-14, before any held-out result** (after the timing pilot only: 16 inputs of F,
   emulator run time, no outcome looked at). The pilot gave ≈ 15 CPU-s per `ld = :dynamic`
   emulation (≈ 5 s of it the equilibrium) and ≈ 5 s for the default `:spectrum`; with F at 150
   per stratum the task would exceed its budget. F is shrunk to the first **100 per stratum** of
   the same MersenneTwister(20261123) permutation (800 inputs), and F's emulations skip
   `spectrum_equilibrium` (its genic variance and LD are used by (a) only, on H). Nothing else
   changes.
