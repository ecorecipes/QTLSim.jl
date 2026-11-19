# News

All notable changes to QTLSim.jl are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the package follows
[semantic versioning](https://semver.org/) (before 1.0, a minor release may break the API).

## [Unreleased]

### Changed (breaking for numbers)

- **The effective size is the number of breeders, and `drift_critical_rate` finds every steady
  state** (2026-11-19; the two options of validate/effective_size.md, added below as opt-in, are now
  the defaults).
  - *What changed.* `effective_size(N, n)` returns `n` (`scheme = :breeders`); it returned the
    Falconer-type form 4N/(V_k + 2), about the harmonic mean 2Nn/(N + n) of adults and breeders.
    The `ne` keyword of `house_of_gauss_vg`, `house_of_gauss_G`, `house_of_cards_vg`,
    `drift_critical_rate` and `half_survival_rate` defaults to `:breeders`.
    `drift_critical_rate` (and the bracket of `half_survival_rate`) uses `solver = :bisection`:
    where its iteration of the steady state fails it now finds the root by bisection. The earlier
    behaviours are `scheme`/`ne = :falconer` and `solver = :iteration`, with the old values exactly.
  - *Why.* For QTSim's life cycle (n of N adults survive, N offspring with two distinct parents
    among the survivors) the inbreeding, variance and mutation–drift effective sizes are n, whatever
    N (the eigenvalue size n + ½): derivation, 27 of 27 neutral kernel cells, Lean
    `EffectiveSize.lean` (`mutation_drift_balance`, `inbreeding_step`, `leading_bounds`,
    `falconerSize_gt_breeders`). The old form overstates N_e whenever fewer breed than the census
    (1.1 to 2 times with three quarters to a tenth surviving). The iteration alone returned less
    than the maximum over steady states, or 0, when the variance of a large population carries a load
    near ln R. When the breeders vary between generations the size is their harmonic mean; the
    functions are called with the expected number, which for binomial survival is less than one
    individual higher (not corrected; `effective_size`'s docstring).
  - *How much.* `drift_critical_rate`'s k_c is 0.72–0.89 times its earlier value at a fixed size
    (median 0.77 on 600 published inputs) and 0.83, 0.90 and 0.94 times (medians) under
    Beverton–Holt with R = 5, 2 and 1.5; its `M` falls too (1.75–1.79 → 1.21–1.25 at a fixed size
    in the test's sets). The solver changes 12 of 3360 published-type values (none at a fixed size or
    R = 5; all upward, five from 0) and 3 of 16 Xu–Shen-type ones. `house_of_gauss_vg(p)` and
    `house_of_cards_vg(p)` with everyone breeding use N_e = N instead of N + ½ (a relative change
    of the order of 1/(2N)); with `nsurv < PopSize` they are lower than before.
  - *What does not change.* The kernel; the moments lowering (`run_moments`,
    `multilocus_equilibrium`, `rare_alleles_vg`, `gaussian_vg_loci`, the spectrum) and `emulate`,
    bit for bit (validate/ne_default_unchanged.log: 486 outputs on 24 published inputs; the
    emulators' starting variance is kept as it was, `QTLSim.burnin_start_vg`); `half_survival_rate`'s
    and `emulated_half_survival_rate`'s k50 to the resolution of their searches (1e-11 relative for
    the latter), though the `kc` they return is the new one; `sustainable_rate`, `critical_rate`
    (the paper's recipe) and vignette 09's reproduction of the paper's table, which has its own code.
  - *The test of `drift_critical_rate`, rerun.* validate/critical_rate.jl on fresh seeds (addendum
    written before the run; tolerances unchanged; grid re-centred): the cost ratio k50(c)/k50(0) is
    within 0.06 of the prediction in **13 of 24** cells (run 1 with the old N_e: 24 of 24, largest
    error 0.054; now 0.135, the kernel's ratio above the prediction in 24 of 24), the fitted scale
    M_obs/M is 1.19–1.75 and within 0.8–1.25 in **2 of 8** groups (was 1.09–1.22, 8 of 8), the
    linear form still fits better than the square root in 8 of 8, `sustainable_rate`'s ratio is
    within tolerance in 0 of 24, and `half_survival_rate` is within 0.75–1.33 in 26 of 32 (as
    before; 100 loci 0.90–1.10, 500 loci 0.69–1.08). Verdicts by the pre-set rules: relative effect
    of a cost MIXED (was SUPPORTED), linear-margin form MIXED (was SUPPORTED),
    `half_survival_rate` SUPPORTED. The kernel's k50 is 0.96–1.37 × k_c with 100 loci and 0.56–0.88
    with 500 (was 0.77–1.03 and 0.44–0.75). So the steady-state theory gets the direction, the linear
    shape and the ordering of demographies, and overstates the size of a cost's effect (its M is
    19–75% low); the earlier agreement depended on the overstated N_e. Run 1 is kept as
    validate/critical_rate_falconer.log (and `critical_rate.jl falconer` reruns it). On the paper's
    710 populations k50 is 0.38–1.12 × k_c (was 0.29–1.06; validate/critical_rate_posthoc.log).
  - *Also updated.* Tests (the recursion's fixed point and eigenvalue, the closed forms at N_e = n,
    `:falconer` against recorded values, `drift_critical_rate` against an independent maximum over all
    steady states); symbolic checks (`check_critical.jl` 62 checks, both effective sizes;
    `check_effective_size.jl`); validate/xu_shen.log section 5 and xu_shen.md section 3 (eq. 6 over
    k_c is now 1.03–1.17 for Xu–Shen-type inputs and 1.07–1.51 for published-type ones; the first
    run is xu_shen_falconer.log); docs/src/plasticity_design.md; the prose of vignettes 02 and 14
    (their computed outputs do not change; the `.md`/`.html` were edited by hand, not re-rendered).
    Earlier validation scripts that used `drift_critical_rate` or the house of Gauss to choose inputs
    or starting values now name the settings they were run with (`ne = :falconer`,
    `solver = :iteration`, `QTLSim.burnin_start_vg`), so they still reproduce their logs; none was
    rerun. Entries below quote numbers as they were when written.

### Added

- **QTSimFast, a faster copy of QTSim with the same output bit for bit** (reference/QTSimFast; the
  package and its defaults are unchanged, and `find_qtsim()` never returns it: the original in
  reference/QTSim stays the reference). A modified version of QTSim's C source (GPL-3.0-or-later,
  © 2023 Andrew Meade; modified 2026-11, each file says how): for the same options file its standard
  output and exit code are the original's. Genomes are shared blocks of loci copied only on a
  crossover or a mutation; genotypes are summed in the original order, eight individuals side by
  side and optionally on threads (`QTSIM_THREADS`, default 1; the output does not depend on it);
  MT19937 is inlined and checked against `gsl_rng_mt19937` at start, every other draw is GSL's own
  in the original order; crossover maps are neither allocated nor sorted; statistics the original
  computes twice or never uses are computed once or not at all. Cached block sums, as in the Julia
  kernel, are ruled out by bit parity. `reference/QTSimFast/test_parity.sh` compares the two
  programs (`make -C reference/QTSimFast parity`, a step of CI's `QTSim parity` job; `parity-full`
  is the whole matrix: published inputs of every stratum, the single-locus protocol, a synthetic
  matrix of every option, thread counts 1, 2, 3, 4 and 8, block sizes, sanitizer builds).
  Timings are not reported yet: `make -C reference/QTSimFast bench` writes
  reference/QTSimFast/bench_results.md (original, QTSimFast at 1, 2, 4 and 8 threads, the Julia
  kernel) and is to be run on a quiet machine. See reference/QTSimFast/README.md.
- **An individual-based test of the bistability predicted by Xu and Shen (2026, bioRxiv preprint)**
  (validation only, nothing in src/; validate/bistability.md, criteria committed after one disclosed
  pilot and before every other run). 36,000 kernel runs under Beverton–Holt from prescribed initial
  states (N₀, V₀, Δ₀) at 0.6 and 0.85 of the critical rate, two genetic architectures × three (K, R),
  against four deterministic predictions (their eqs. 1–5 at two mappings; a discrete recursion with
  QTLSim's ingredients and either N_e). QTSim-type genetics: persistence depends on the initial
  variance and lag, hardly on the initial size; both basins occur at one rate in 4 of 6 cell × rate
  units (B1 FAIL by the pre-set 75%; the other two persist from the whole grid, as predicted); the
  deterministic boundaries match 84–91% of majority outcomes (B3 PASS for all four), with the kernel
  persisting from 0.56–0.77 of the predicted boundary variance (a smaller extinction basin, the
  opposite of the direction expected from stochasticity). A 50-locus high-mutation "Gaussian" arm
  dies from every state at 0.85 k_c and, for K = 3000, at 0.6 k_c (B1 and B3 FAIL for all four
  predictions): post hoc, its genic variance erodes under the move (to 0.22 in 500 generations; not
  with 500 loci, not changed by more crossovers), which neither variance equation contains.
  Monotonicity and all null arms pass. Figures in validate/bistability_figs/.
- **`effective_size(N, n; scheme)` and `drift_critical_rate(p; ne, solver)`: two options, both off by
  default** when added; **both are the defaults since 2026-11-19, see "Changed" above** (every default is bitwise unchanged; validate/effective_size.md, criteria committed first).
  (1) `effective_size(N, n)`, Falconer's 4N/(V_k + 2) with the breeders' V_k, is about the harmonic mean
  2Nn/(N + n) of adults and breeders. For QTSim's life cycle (n of N adults survive, N offspring with
  two distinct parents among the survivors) the inbreeding, variance, mutation–drift and eigenvalue
  sizes are all n (derivation in the criteria file and validate/symbolic/check_effective_size.jl): the
  formula puts a V_k whose mean is 2N/n into the expression for a mean of 2. Neutral kernel runs
  (N = 200 and 1000, a tenth to all surviving; decay of variance, mutation–drift balance, pair
  coalescence times from recorded pedigrees) reject the formula in 23 of 27 cells (it is 1.1 to 2 times
  the measured size when fewer than all survive) and fit the breeders' harmonic mean in 27 of 27.
  `scheme = :breeders` returns n; `house_of_gauss_vg`, `house_of_gauss_G`, `house_of_cards_vg`,
  `drift_critical_rate` and `half_survival_rate` take `ne = :breeders`. The moments lowering and the
  emulators never used the function (their drift is that of the n parents), so only closed forms called
  with `nsurv < PopSize` and `drift_critical_rate` are affected: its k_c falls to 0.72–0.98 of the
  default on published inputs (median 0.77 at a fixed size, 0.83, 0.90, 0.94 under Beverton–Holt with
  R = 5, 2, 1.5). With the breeders, validate/critical_rate.jl's agreement on the effect of a cost is
  worse (cost ratios within tolerance in 12 of 24 cells instead of 24, M_obs/M 1.18–1.77 instead of
  1.09–1.22): the earlier agreement owed something to the overstated N_e. Which default to use is left
  open.
  (2) `drift_critical_rate` iterates V = f(V) from the variance at N = K; where that variance carries a
  load the growth rate cannot pay, or |f′| > 1, the iteration reports no population or does not
  converge, and the rate is below the maximum over steady states or 0. `solver = :bisection` finds
  the (unique) root there and is otherwise the same computation: 3 of 16 Xu–Shen-type inputs change
  (0 → 0.297), no published input at a fixed size or with R = 5, 5 of 600 with R = 1.5 without a cost
  (one from 0 to 0.077, for which the kernel holds a population of 530 adults that follows the optimum
  at half that rate). Not the default, because published-type values under Beverton–Holt change.
  Tests: topics `theory`, `critical`.

- **The eco-evolutionary model of Xu and Shen (2026), re-implemented from its equations** (additive, a
  theory module beside the simulation; no default and no existing output changes; src/xu_shen.jl).
  Xu, K. and Shen, H. (2026), "Eco-evolutionary feedbacks generate bistability in population
  persistence under gradual environmental change", bioRxiv preprint, doi:10.64898/2026.08.05.743160
  (version 1; not peer reviewed; we know of no published code). The lag, the genetic variance and the
  population size as three coupled equations in continuous time: `XuShen` (parameters; the
  preprint's standard set by default), `xu_shen_field`, `xu_shen_critical_rate` (their eq. 6, by
  bisection on the monotone factor of h′), `xu_shen_equilibria`/`xu_shen_equilibrium` (eqs. 7–8: a
  saddle and a stable node below the critical rate), `xu_shen_jacobian`, `xu_shen_stability` (S8,
  Routh–Hurwitz), `xu_shen_trajectory` (an adaptive Rosenbrock integrator written here, no new
  dependency; extinction = N below a threshold, default 1e-3), `xu_shen_persists`,
  `xu_shen_boundary` (the basin boundary of their Figure 3), and their three scenarios
  (`xu_shen_balance_variance`, `xu_shen_critical_size`, `xu_shen_max_rate`, `xu_shen_max_duration`).
  Their eq. 6 has priority over `drift_critical_rate`, which is a discrete-generation counterpart
  matched to QTSim's life cycle; the feedback is Bürger and Lynch (1995).
  Tests: topic `xu_shen` (422). Symbolic checks: validate/symbolic/check_xu_shen.jl (77), with notes
  where a printed equation and the derivation differ: S17's first term has the opposite sign to the
  one that follows from S16 (so S19–S21 do not bound P₁ as written); the conclusion of S2 stands,
  because a₁a₂ − a₃ is a sum of positive terms at every equilibrium; no result is affected.
  Figures (validate/xu_shen.jl, validate/xu_shen.log, validate/xu_shen_figs/, summary
  validate/xu_shen.md; read-offs from the preprint's figures committed first in
  validate/xu_shen_targets.md): Figures 1, 2a (three of four curves), 3 and 4a, b, d agree with the
  read-offs within reading precision; the "V_M = 0.5" curve of Figure 2a matches V_M = 0.2 and
  panels 2b, c match K = 500 better than the caption's K = 1000 (both post hoc); Figure 4c's rate
  is not stated. The read-offs of Figure 4d were first taken against a misplaced axis and were read
  again after the comparison (disclosed in the log). Comparison with `drift_critical_rate`
  (descriptive): at mapped parameters eq. 6 is 1.01–1.10 times it for Xu–Shen-type inputs and
  1.01–1.27 times it for published-type inputs under Beverton–Holt (R = 1.5–5); the differences are
  drift from the census against N_e from the breeders (×1.07–1.39), logistic against Beverton–Holt
  regulation (×0.78–0.94) and the house-of-Gauss term α²/(2V_s) (×1.00–1.16), with the weak-selection
  forms below 0.5%. A fixed size has no counterpart in the preprint's model.

- **`k50_gradient`: the emulator's half-survival rate and its elasticities on a target that does not
  saturate** (additive; `emulate`, `extinction_distribution` and every default are unchanged;
  src/k50_gradient.jl, with `QTLSim.emulated_log_hazard` and `QTLSim.log_cumulative_hazard`). Phase A
  of the elasticity study left 12 of 64 near-boundary inputs undefined, all with ω < 2.2. The cause:
  under strong selection `P2000` is exactly 1 in `Float64` while the cumulative hazard is below
  1.1e-16 (the product of survival factors does not move), and exactly 0 once the deterministic
  moments trajectory dies inside the horizon (taken as certain extinction, although the hazard
  accumulated until then can be far below log 2), with the whole range of log H crossed over 1e-3
  to 1e-9 in log k. The target is now log H summed in logs along the trajectory continued through
  its deterministic extinction (quadratic extrapolation of log W̄, used over one or two
  generations); k50 is its root, and the elasticities are the implicit-function formula with a
  central difference spanning about ±0.5 in the target, refused across a fallback switch of the beta
  solver. Finite differences, not AD. Validation under phase A's unchanged criteria with a
  pre-committed addendum (validate/elasticity_criteria.md, validate/elasticity_b.jl,
  validate/elasticity_b.log; 128 inputs, 1.5 CPU-h): 126 of 128 defined (phase A's 12 undefined
  inputs all defined; the two undefined, ω = 0.05 and 0.16, have a true jump of the trajectory at
  k50). **G1 PASS** (sign 41/41), **G2 PASS** (size 34/40), **G3 PASS** (Spearman 0.921) over all 48
  cells, none excluded; but **all four low-ω cells fail G2**: the emulator's e_ω there is −0.43 to
  −0.63 against the table's −0.14 to −0.19 (right sign, about three times too strong). N1 stable
  (max 0.003); **N2 NOT VERIFIED** as written (19 of 20 values within 0.007; one, e_ω at ω = 0.37,
  −0.33 against phase A's brute force −1.09; post hoc, brute force over ±2% and ±0.5% gives −0.33,
  so the ±5% brute-force value is not local), hence the verdicts carry the flag "numerically
  unreliable" under the rule. On phase A's 52 defined inputs the new values agree within 0.02 on
  all 42 without a fallback shift and differ by up to 0.16 on the 10 with one (phase A's values
  carried the jump). The emulator reads V_P, h² and m only through V_E and m·V_E, so
  e_h² = −h²/(1 − h²)·e_V_P (to 5e-4 in 90% of inputs; three low-ω inputs off by 0.007–0.058,
  not diagnosed) and e_V_P − e_m is the elasticity in V_E: below 0.005 for ω ≥ 2.2, down to −0.69
  for ω < 1. Tests: topic `k50_gradient` (78).

- **A two-stage emulator** (off by default; `emulate(p; loci = :two_stage, band, screen, blocks,
  replicates, seed)`, `emulated_half_survival_rate(p; loci = :two_stage)`,
  `emulate_sweep(…; loci = :two_stage, seed)`; src/two_stage.jl). Stage 1 is the deterministic
  `:spectrum` emulator with `screen = (; ld = :dynamic, solver = :beta)` (`QTLSim.TWO_STAGE_SCREEN`)
  on every input. Stage 2 is the block model (160 blocks, Gaussian law, 40 replicates), run only
  where the input's rate is within `exp(±band)` of the screen's own half-survival rate
  (`band = QTLSim.TWO_STAGE_BAND = 0.25`). That condition is decided by one more screen evaluation,
  toward the flip (`QTLSim.two_stage_screen`). Results record `stage` and the screen's values. The
  half-survival rate brackets with the screen and bisects the block model from ±band around it,
  with common random numbers; `blocks_half_survival_rate` gained a `bracket`.
  Design, on development inputs only (validate/two_stage_dev.log: blocks.jl's 120 near-threshold
  and 320 table inputs with cached block results): w = 0.25 is the smallest width at which the screen
  disagrees with the block model on < 1% of inputs outside the band. The margin is narrow (one input
  at distance 0.246). The beta and grid screens make the same calls at k on all 440 inputs, their
  band decisions agree on 98%, and the grid would need w = 0.4; the beta is 3.5× cheaper per evaluation.
  Pre-specified validation (validate/two_stage_criteria.md, validate/two_stage.log; fresh inputs
  and seeds): on 120 near-threshold inputs, 40 kernel runs each, **(T1) PASS**: accuracy 0.992 and
  AUC 0.985, against 1.000/1.000 for the block model alone, one discordant input. ld = :dynamic alone
  scores 0.883 and the default `:spectrum` 0.817. **(T2) PASS**: pooled x50 ratio 0.967, strata
  0.88–1.07. **(T3) PASS**: median cost 4.7 s against 65 s for the block model alone (0.07), with 25%
  of these inputs in stage 2. On 2000 fresh table inputs, 4.05% go to stage 2. E1 PASS (AUC 0.999),
  E3 PASS and E4 PASS (0.95–1.06), but **E2 FAIL**: the deciles are within 0.003, yet the calibration
  slope is 0.50. Cost per input is 0.17× the block model alone and 0.93× ld = :dynamic (grid) alone. A
  two-stage run of the whole table is estimated at ≈ 27 h on 8 threads (not run). Not a default.

- **Haplotype-block lineages** (off by default; `run_blocks`, `emulate(p; loci = :blocks,
  blocks, within, partition, replicates, seed)`, `emulated_half_survival_rate(p; loci = :blocks)`,
  `emulate_sweep(…; loci = :blocks, seed)`; src/block_lineages.jl). The kernel's life cycle on
  `PopSize` individuals whose chromosomes are reduced to the sums of `blocks` blocks of
  consecutive loci (`partition = :map`, equal shares of QTSim's genetic map, or `:position`),
  so that associations between alleles and their linked backgrounds arise by construction
  rather than by a moment closure. QTSim's switch law is exact (`QTLSim.qtsim_switch_law`,
  P(c) = Σ_K Pois(K; u) C(M−c−1, K−1)/C(M, K)); a switch inside a block is either moved to the
  nearer boundary (`within = :snap`: exactly the kernel with a coarsened map, and the kernel
  itself with one block per locus) or kept, with Gaussian segregation of that block
  (`within = :gaussian`: ȳ + f d_A + (1 − f) d_B + √(2f(1 − f)V_b) ξ). Replicates are seeded by
  `QTLSim.job_seed`. Supported: one diploid trait, fixed size, viability selection, continuous
  move. `QTLSim.block_map`, `QTLSim.block_stats`, `QTLSim.block_equilibrium!`. Symbolics:
  validate/symbolic/check_blocks.jl (79 checks). Design (validate/blocks_dev.log, development
  inputs only): snapping does not converge in `blocks` up to 320 (V_G keeps rising: removing
  recombination inside a block matters even at map shares ≪ 1/B, Hill–Robertson within blocks);
  the Gaussian law meets the pre-stated convergence rule at 160 blocks (V_G 0.92–0.96 of the
  kernel's through +200 on the development inputs); its cost is only 1.3–2.9× below the kernel's
  and 4–9× above one `:spectrum` call.
  Pre-specified validation (validate/blocks_criteria.md, validate/blocks.jl, validate/blocks.log;
  160 blocks, Gaussian law; fresh inputs and seeds): on 120 inputs near the threshold, 40 kernel
  runs each, **(P1) threshold accuracy 0.992, AUC 1.000 PASS; (P2) pooled x50 ratio 0.980 PASS;
  (P3) x50 per stratum 0.81–1.02 PASS; (P4) not worse than `ld = :dynamic` PASS** (accuracy 0.975;
  default `:spectrum` 0.925, pooled x50 1.18; `:rare_alleles` 0.708). V_G at the move and at +50,
  +100, +200 at 0.97–0.98 of the kernel's (S1 PASS), the fall of block-genic variance under the
  move within 1% (S2 PASS). On 320 fresh table inputs (10 replicates) E3 PASS; E1, E2 and E4
  FAIL formally (a stratum with no survivor, separation, one 100-locus stratum at 1.32). Cost
  (S4) FAIL: 7.5 s per call against 5.8 s for `ld = :dynamic`; a run is only 1.1–1.8× cheaper
  than a kernel run. Not a default.

- **A beta closure for the `:spectrum` model** (off by default; `spectrum = (; solver = :beta)` in
  `run_moments`, `spectrum_equilibrium`, `emulate`, `emulated_half_survival_rate`; default
  `:grid` bit for bit unchanged). Per effect class the moments `P_k = ∫ p^k ψ`, `ψ = pqφ`,
  `k = 0, 1, 2`, closed by `ψ/P₀ = Beta(α, β)` (exact on the neutral spectrum `ψ ∝ 1 − p`), one
  explicit Euler step per generation of the exact ψ-moment equations (`QTLSim.beta_dP`,
  `QTLSim.beta_closure`) with Barton's term and the same mutation input, projected onto the
  Hausdorff conditions. Works with `ld = :dynamic`. A deterministic fallback hands the state to
  the grid (V_G kept exactly) for the rest of the run when a class carrying ≥ 0.1% of V_G has
  scaled directional selection `4n|aλ|/S > 100`, a U-shaped closed law, or needs the projection;
  results then carry `spectrum_fallback` (`run_moments`), `fallback` (`emulate`) or `fallbacks`
  (`emulated_half_survival_rate`); `QTLSim.spectrum_fallback`. The threshold was set on four
  inputs and the in-place handover chosen on a 40-input pilot (disclosed in the criteria).
  Pre-specified numerical-equivalence check against the grid: `validate/spectrum_beta_criteria.md`,
  `validate/spectrum_beta.jl`, `validate/spectrum_beta.log` (400 fresh published inputs):
  **11 of 13 PASS**. Equilibrium V_G median |r−1| 0.0005, 95th pct 0.009 PASS; V_G at +100
  0.0017 / 0.089 PASS; **V_G at +500 FAIL** (median 0.011, 95th pct 0.111 > 0.10; mostly the
  beta's slow downward drift, 0.80–0.94 of the grid, in runs that never fall back); P2000 calls
  400/400 PASS; T ratio median 1.000 PASS; half-survival rate (48 inputs) median ratio 0.998, 95th
  pct 0.027 PASS; with `ld = :dynamic` (200) the same, except **+500 FAIL** (95th pct 0.160) and
  P2000 calls 199/200; fallback rate 0.245 PASS; speed-up of `emulate` median 267× (min 0.90×,
  fallbacks 116×) PASS, but only 3.4× with `ld = :dynamic` (the LD-class step dominates). Use it
  for equilibria, emulator calls and half-survival rates; not for V_G hundreds of generations
  into a move. Symbolic checks in `check_spectrum.jl` (h), 79 checks pass.
- **Dynamic linkage disequilibrium in the `:spectrum` closure** (off by default;
  `spectrum = (; ld = :dynamic)` in `run_moments`, `emulate`, `emulated_half_survival_rate`,
  `spectrum_equilibrium`). `V_G` = the spectrum's genic variance + `D_t`, the LD of classes of
  pairs of distinct loci by position bins (`ldbins = 20`) and split probability under
  `CrossoverModel` (`QTLSim.ld_classes`, `QTLSim.LDClasses`), one generation
  `d_c ← (1 − 1/(2n))[(1 − ρ_c)(1 + 2γ)d_c − ω_c β_a β_b g²/(2S)]` (Bulmer's generation with each
  bin's regression `β_b` of the genotype on its loci, recombination, drift, and `γ` the sites'
  turnover; `QTLSim.spectrum_ld_step!`); each site's directional selection is scaled by
  `V_G/genic` (response `V_G λ/S`) and its stabilising selection by the bins' mean `β²`. Nothing
  fitted. Without drift, one bin and in the infinitesimal limit the fixed point is the
  multilocus model's (tested). Development (`validate/spectrum_ld_dev.jl`, `.log`): the
  multilocus model's LD is too strong mainly because it gives every pair the mean locus
  covariance with the genotype, while under QTSim's lowest-point crossover the tight pairs sit
  toward the chromosome's end, where that covariance is reduced by their own LD; resolving LD by
  position puts V_G/genic at the move at 0.997 of the kernel's (median, discovery inputs; mean
  field 0.914). The option does not explain the kernel's early fall of the genic variance under a
  move; Barton's linked-selection term is kept unchanged (an LD-consistent replacement was
  worse). Default `:spectrum` output is bit for bit unchanged; `run_moments(…; VG)` with the
  option sets the genic variance. Pre-specified validation: `validate/spectrum_ld_criteria.md`,
  `validate/spectrum_ld.jl`. Results (`validate/spectrum_ld.log`; 160 held-out inputs near the
  threshold, 40 fresh kernel runs each): **(a) PASS** (V_G at the move model/kernel median 1.000,
  genic 1.016, LD/genic 1.009), **(b) PASS** (V_G at +100 1.027, +500 1.103, IQR to 1.375),
  **(c)** E6′c PASS (accuracy 0.925, AUC 0.926; default `:spectrum` 0.863/0.870), pooled x50 ratio
  1.068 PASS (default 1.366), per-stratum x50 FAIL (QT0100-5 0.72, QT5000-3 1.77); **(d)** on 800
  fresh table inputs E1 PASS (0.983), E2 FAIL (slope 0.33; P2000 nearly binary), E3 PASS, E4 FAIL
  (QT1000-5 0.79, QT5000-3 1.66; the default fails E4 in QT5000-3 too). Not a default.
- **Parameters as an object, with sensitivities** (roadmap C3; `src/para.jl`).
  `QTLSim.Para(f, base, names)` sees a prediction `f(::Params)` as a function of the named
  `Float64` fields of `base` (integer, `Bool`, `Symbol` and array fields are an
  `ArgumentError`); `sensitivities(P; rel = 1e-3, log = true)` gives elasticities
  `∂log y/∂log θ` (or plain derivatives) by central differences with one level of
  Richardson extrapolation. They are finite differences: `Params` holds `Float64`s, so
  automatic differentiation would need it to be generic over the number type.
- **An emulator of the published sweep** (roadmap F4; `src/emulator.jl`). `emulate(p)` runs
  the moments lowering with drift (`erosion = true`, default `loci = :rare_alleles`) 2000
  generations at a fixed optimum from the house-of-Gauss `V_G`, then 2000 along the
  continuous move, and adds demographic stochasticity with `extinction_distribution`:
  `P2000`, the median extinction time `T` (counted as the table counts), `P240k` by
  constant end-hazard extrapolation, `V_G` at the move. `emulated_half_survival_rate(p)`
  solves `P2000 = 1/2` (bracketed regula falsi on the log cumulative hazard, to 1e-10 in
  log k); `emulate_sweep` runs a batch through `sweep` with a resumable cache.
  `validate/emulator.jl` applies the criteria fixed beforehand in
  `validate/emulator_criteria.md` (four clarifying amendments recorded before any
  comparison) to all 119,984 published inputs (`validate/emulator.log`). Primary closure
  `:rare_alleles`: **E1 PASS** (AUC of P2000 for alive at 2000 = 0.958; per stratum
  0.924–0.974); **E2 FAIL** (the emulator's P2000 is nearly binary, 0 or 1 for all but 213
  inputs; in the top decile, mean P2000 1.000 against 58.8% observed alive, largest
  decile gap 0.412; calibration slope 0.227); **E3 PASS** (Spearman 0.941 between
  predicted and observed extinction times, median ratio 1.048); **E4 FAIL** (the emulated
  half-survival rate is 1.19–1.74 times the table's, within [0.8, 1.25] only at 100 loci,
  and rising with the number of loci); E5 (exploratory) AUC of P240k for survival to
  240,000 = 0.941. On the pre-specified sample of 1000 (no verdicts): `:spectrum` has
  AUC 0.981, x50 ratios 0.96–1.12 in 7 of 8 strata (1.34 at 5000 loci, u = 5), median
  T ratio 1.000, but is as badly calibrated (top-decile gap 0.19, slope 0.27);
  `:multilocus` is no better than `:rare_alleles` (x50 ratios 1.05–3.46). Exploratory
  parameter map (2000 inputs): emulated x50 median 0.17 (0.16 at 100 loci, 0.16–0.19 at
  500–5000), elasticities of k50 about +0.6 in PhenotypeVar and in MutationScalar, −0.23
  in Heritability, −0.21 in InitFitnessSD (IQR spans 0) and near 0 in the mutation rate.
- `pair_classes` and `pair_recombination` are memoised by (loci, mean, crossover model):
  the same values, about 0.4 s saved per moments run at 5000 loci.
- `validate/published_table.jl`: the published table's loading and row matching, shared
  by `published_sweep.jl` and `emulator.jl`, with an index (`row_index`) for matching all
  119,984 inputs.

- **Between-run variance in the emulator** (roadmap F4, off by default; `src/emulator.jl`,
  `ext/QTLSimCatlabExt/moments_noise.jl`). `emulate(p; replicates = R, seed)` (or `rng`)
  shares the deterministic burn-in and runs `R` replicate moments trajectories through the
  move (`QTLSim.moments_replicates`), each with its own realised mean genotype and `V_G`:
  the mean drifts by the drift-mode lowering's own increase of `C` (Ploidy²·ΔC, about
  `V_G/n + V_m/N` per generation), starting from its stationary spread `Ploidy²·C`, and
  `V_G` fluctuates as the stationary Poisson random field of `stationary_spectrum` predicts
  (`QTLSim.vg_fluctuation`: Var V_G = 8nLμE[a⁴K(γ)], per-generation innovation
  4LμE[a⁴M(γ)]); `vgnoise = false` keeps the mean noise only. P2000 averages the
  replicates' `extinction_distribution`s. `emulate_sweep(...; replicates, seed)` seeds
  each input by `QTLSim.job_seed(seed, id)` (FNV-1a), so results do not depend on threads
  or the cache. `replicates = 0` is bitwise the deterministic emulator. Validation under
  `validate/emulator_variance_criteria.md` (written before the code; amendment 1, V_G noise,
  and amendment 2, smaller samples after the pilot, both before any result; disclosed: the
  change was motivated by E2 on the same table), `validate/emulator_variance.log` (about 6.5
  CPU-h): on 10,000 fresh inputs (seed 20261115), primary `:rare_alleles` with R = 32:
  **E1 PASS** (AUC 0.965; deterministic 0.958), **E2 FAIL** (slope 0.243, top-decile gap
  0.362; deterministic 0.228 and 0.425), **E3 PASS** (Spearman 0.940, median ratio 1.046),
  **E4 FAIL** (x50 ratio 1.20–1.91; deterministic 1.22–1.92). **E6 FAIL**, the independent
  test against 40 fresh kernel runs on each of 192 inputs near the deterministic k50
  (`init = :equilibrium`, burn-in 1000): slope 0.14 (deterministic 0.11), top bin mean
  P2000 0.997 vs 0.444 alive (deterministic 1.000 vs 0.418), log-likelihood −28,544 vs
  −34,499, AUC over runs 0.81 vs 0.76. The kernel is itself nearly binary per input (0 of
  40 alive for 138 inputs, 40 of 40 for 37, in between for 17), so the emulator's
  miscalibration is mainly bias (its k50 is too high), not missing variance. Diagnostics: the
  spread of the mean before the move matches the kernel (variance ratio 1.04) and so does
  V_G (1.03); the kernel's CV of V_G is 1.8 times `vg_fluctuation`'s (LD and linked sites are
  left out). Mean-only noise changes almost nothing. `:spectrum` with mean noise (R = 8) on
  the E6 inputs (no verdict): mean |f − P| 0.136 (primary 0.327), slope 0.16.
- **A fresh validation of the emulator with `loci = :spectrum`** (roadmap F4; no code
  change). `validate/emulator_spectrum.jl` under `validate/emulator_spectrum_criteria.md`
  (committed before the script; amendment 1, before any result: E6′'s candidates come from a
  new pool, because the parameter map had too few near-k50 inputs left), log
  `validate/emulator_spectrum.log`, about 7.5 CPU-h. Disclosed: `:spectrum` was chosen after it
  did better on earlier samples, and every input whose `:spectrum` result had been seen is
  excluded (1191 inputs). Primary: deterministic `emulate(p; loci = :spectrum)`. On 2400 fresh
  inputs (300 per stratum, seed 20261116): **E1 PASS** (AUC 0.977; `:rare_alleles` 0.953),
  **E2 FAIL** (slope 0.26, top-decile gap 0.24; `:rare_alleles` 0.22 and 0.47), **E3 PASS**
  (Spearman 0.952, median T ratio 1.000), **E4 FAIL** (x50 ratio 0.94–1.27, with 1000 loci and
  u = 3 at 1.251 and 5000 loci and u = 3 at 1.268 outside [0.8, 1.25]; `:rare_alleles` 1.10–2.27).
  **E6′ FAIL**, the independent test: 160 new inputs near the `:rare_alleles` k50, 40 fresh
  kernel runs each (`init = :equilibrium`, burn-in 1000): the kernel is again near-binary (91 /
  58 / 11 inputs with none / all / some alive). E6′a FAIL (top bin P2000 0.999 vs 0.711 alive),
  E6′b FAIL (slope 0.16), E6′c threshold test FAIL (accuracy 0.838, AUC 0.859, both against
  0.85; 25 inputs the emulator calls alive and the kernel dead, 1 the other way), E6′d FAIL
  (x50 ratio 0.88–3.12 by stratum, pooled 1.66; 20 inputs per stratum give wide bootstrap
  intervals). `:rare_alleles` on the same inputs: accuracy 0.638, pooled x50 ratio 3.38, mean
  |f − P| 0.36 (`:spectrum` 0.17). Post hoc: on the E6′ inputs the published table agrees with
  the kernel on 158 of 160, so near the threshold `:spectrum`'s k50 is still too high, despite
  its E4 ratios near 1. The default stays `:rare_alleles` (the user's decision).
- **Fixed** `stationary_mean_var` (F1 sampler): when the genic target exceeds twice the
  spectrum's S (narrow fitness width, e.g. published input QT1000-Possion3.0/Input-03628) it
  returned a negative variance and `equilibrium_genomes` threw a `DomainError`; S is now
  taken at the target's V_G in that case only (found during the run above).
- **The targets of the lowerings are categories** (roadmap C1; `ext/QTLSimCatlabExt/typed.jl`).
  The lowering functor now lands in `(Stage{B}, Morph{B})`: a `Morph{B}` is an image with
  its domain and codomain stages. Each generator's image is checked when the life cycle is
  lowered (a morphism of the target, acting on as many demes as its step has stages,
  callable as the backend calls it), so a wrong image is an `ArgumentError` at `lower`
  instead of a failure mid-run; composites check that their stages meet.
  `lower(f, target; typed = true)` returns the `Morph` (by default `lower` returns the image,
  as before). `LoweringCategory{B}()` declares each target (`:kernel`, `:eidos`,
  `:moments`) as a GATlab model of `ThSymmetricMonoidalCategory`. `test/laws.jl` checks the
  laws on random composites: associativity, units and functoriality exactly in all three
  targets (random bracketings of the steps' images with identities inserted, against the
  image of the whole composite: the same population bit for bit, the same Eidos text, the
  same moments); the braiding an involution; interchange and the braiding's naturality
  exactly in the moments lowering and, since they reorder the demes' random draws, in law
  in the kernel (Welch z over 150 replicates per side, |z| < 4, with a check that the test
  tells apart two composites that differ in law). Catlab normalises associativity in the
  syntax but not interchange or the braiding's naturality, so those tests are not
  syntactic. Also fixed on the way: an Eidos lowering given a non-Eidos image, or a moments
  lowering given a `DemeOp`, now reports why.

- **The naturality harness covers every validation family** (roadmap C2). A `:qtsim`
  backend runs the C reference (`ENV["QTSIM"]`, the bundled build or the paper's binary;
  `options.qtsim.qtsim`) on an options file written by `write_qtsim_options`, one seed per
  replicate, and parses its table into the kernel's columns (`gen`, `pheno_mean`,
  `geno_var`, `nsurv`, `theta`, `omega`, `popsize`, `extinct` on `Pop all died.`); it runs
  QTSim's generation of one deme only and rejects other life cycles, several demes and the
  Julia-only settings (drift, several traits, plasticity, fluctuating moves, `:all`
  crossovers, density regulation, fecundity selection) with an `ArgumentError` before
  anything runs. Two criteria: `:welch` (two stochastic backends, |t| < tol, the rule of
  `compare_to_qtsim.jl`) and `:diff` (the point difference alone, the rule of
  `reaction_norm.jl`). A statistic may return `missing` for a replicate
  (`QTLSim.surviving(s)`: extinct runs left out, as the scripts do); the result reports how
  many replicates were used. New statistics `QTLSim.lag_mean(window; squared)` and
  `QTLSim.persistence(start)`, and `QTLSim.welch_t`; a case's `note` is printed. Nine new
  cases in `validate/naturality_cases.jl`: `compare_to_qtsim.jl` A and B (kernel vs
  QTSim), `slim/compare_fluctuating.jl` periodic and red noise, the periodic moments
  against `expected_lag2` (test/fluctuating.jl), `slim/compare_demography.jl`'s moving
  optimum, `slim/compare_plasticity.jl` scenario 1, `reaction_norm.jl` C1 and `fst.jl`'s
  island k = 4, Nm = 1 against `neutral_deme_variances` fed `migrate!`'s own law of pairs'
  origins; every recorded verdict reproduced on the harness's own seeds
  (`validate/naturality.log`).

- **Naturality as the validation harness** (roadmap C2; `src/naturality.jl`,
  `ext/QTLSimCatlabExt/naturality.jl`, `docs/src/naturality_design.md`). Agreement between
  two backends is a naturality square of the lowering functors, σ ∘ K(f) = M(f) ∘ σ; a
  square is now data, a `NaturalityCase` (life cycle, parameters, the two backends among
  `:kernel`, `:slim`, `:moments`, `:theory`, statistics with their tolerances, the criterion
  — `:band`, `:z`, `:exact` or `:ks`, the scripts' rules — replicates, horizon, and the
  verdict the ported script recorded). `check_naturality(case)` (or
  `check_naturality(cycle, p; from, to, statistic, tolerance, …)`) runs it and says whether
  the verdict was reproduced; `naturality_steps(p; loci)` checks the square of every step of
  one generation from matched states (a linkage-equilibrium kernel population and its
  moments), so a disagreement is localised to the step whose square fails: with
  `loci = :aggregate` the square of `reproduce` fails, with `:multilocus` none does.
  `validate/naturality_cases.jl` holds 11 cases ported from `moments_loci.jl`,
  `moments_regimes.jl`, `moments_demes.jl`, `rare_alleles_traits.jl`,
  `local_adaptation.jl`, the neutral exact test and `slim/compare_demography.jl` (three of
  them recorded failures); `validate/naturality.jl` runs them (`naturality.log`). The
  scripts and their logs are unchanged.
- **Closures and model morphisms declared as data** (roadmap C7; `src/closures.jl`,
  `ext/QTLSimCatlabExt/morphisms.jl`). `QTLSim.closures()` declares each closure of the
  moments lowering with its rank in the order of refinements, its domain and its automatic
  domain; `loci = :auto` is now `QTLSim.auto_closure`, the finest closure whose automatic
  domain contains the context (model, demes, ploidy, one law of mutation effects,
  `between`). It gives exactly the earlier choices (a test enumerates 72 contexts against
  the earlier rule) except for one case, since fixed (see Fixed): `QTLSim.moments_loci(p)` of a heritable
  reaction norm said `:rare_local` while `run_moments` used `:multilocus`, as the
  declarations say. `model_morphisms()` and `QTLSim.closure_order()` declare 8 maps (one
  trait into two independent traits, the projection to trait 1, k demes to one, plastic to
  plain at rate (1 − φ)k, and four arrows between closures) with 21 reductions, from which
  a testset is generated; another checks the declared domains against `run_moments`'
  guards.
- **Support matrix and `explain`** (roadmap F11; `ext/QTLSimCatlabExt/support.jl`,
  `docs/src/support.md`, `scripts/support_matrix.jl`). `support_matrix()` probes 20
  features × 3 models × 7 columns (kernel, SLiM, five closures) on tiny populations, so
  "unsupported" is exactly where a guard throws, and reads the result against declared
  validation evidence; `support_markdown` writes the docs page, which a test keeps current.
  `explain(p; cycle)` says which model, closure and backends apply, and why.
- **The critical rate under drift** (`src/critical.jl`; roadmap F8 and the open problem of
  `validate/plasticity_erosion.jl`). `drift_critical_rate(p; cost)`: the lag lowers mean
  fitness, so fewer adults breed (`N W̄` at a fixed size, the Beverton–Holt equilibrium
  `K(RW̄ − 1)/(R − 1)`), drift lowers `V_G` (the house-of-Gauss variance at that `N_e`),
  and the lag grows; the critical rate is the largest steady response `V_G λ/(V_G + V_s)`
  over lags, above which the feedback runs away. The feedback is known (Bürger and Lynch
  1995, from simulations), and a critical rate as the maximum over steady states with
  drift-reduced variance is in Xu and Shen (2026, bioRxiv preprint, continuous time,
  logistic growth; spatial analogues in Polechová and Barton 2015 and Polechová 2026);
  `drift_critical_rate` is a discrete-generation implementation matched to QTSim's life
  cycle, tested against the kernel. A log-fitness cost `c` lowers the
  breeders as the lag load does, so by the envelope theorem
  `k_c(c)/k_c(0) = 1 − c/M + o(c)` (O(c²) when k_c is twice differentiable; proved in Lean as `ratio_linear`), `M = λ*²/s`, twice the lag load at the critical lag.
  `sustainable_rate`'s fixed-V_G threshold has the same first-order form,
  `1 − c/(2 ln(W̄(0)/W_min))`; drift makes the critical lag, and so `M`, 2.5–7 times
  smaller in the tested cells, so a cost bites several times harder. `half_survival_rate(p;
  cost, horizon)`: the rate at which half the populations survive `horizon` generations,
  from the moments lowering with drift (the allele-frequency spectrum by default) plus
  `extinction_distribution`, by bisection (needs Catlab). `sustainable_rate`'s docstring now
  says it is a fixed-V_G threshold, not the half-survival rate. No default or existing
  number changed (`house_of_gauss_vg` unchanged).
  - `validate/symbolic/check_critical.jl` (35 checks): the envelope theorem, the power-law
    case, the elasticities, and the code's `M` against `1/(κν)` and finite differences.
  - `validate/critical_rate_posthoc.log` (existing data only): on the erosion test's cells
    the predicted cost ratio is within 0.05 of the observed in 7 of 8; on the paper's 710
    populations in four demographies the half-survival rate is 0.29–1.06 ×
    `drift_critical_rate` (falling with the number of loci) against 0.07–0.56 ×
    `sustainable_rate`.
  - `validate/critical_rate.jl` (pre-specified, fresh seeds; 960 burn-ins × 13 rates, 100
    and 500 loci, fixed size and Beverton–Holt R = 5, 2, 1.5, log costs up to 1): the cost
    ratio `k50(c)/k50(0)` is within 0.054 of `drift_critical_rate`'s in 24 of 24 cells;
    over the tested range (observed ratios down to 0.25) a ratio linear in the cost fits
    better than a square root in 8 of 8 groups, with `M_obs/M` 1.09–1.22
    (the theory's `M` is about 18% low; its ratio errors lean the same way in 21 of 24);
    `sustainable_rate`'s ratio is within 0.06 in 0 of 24. `half_survival_rate` is within
    0.75–1.33 of the kernel in 26 of 32 cells (100 loci 0.89–1.11; 500 loci 0.68–1.09,
    failing under R = 5 without a cost and R = 2, 1.5 without or with a small cost). Verdicts SUPPORTED by the
    pre-set counts; the tolerances were set after a post hoc look at existing data (the
    script header says so), on fresh seeds.

- **Populations initialised from theory** (`src/init.jl`; roadmap F1). `stationary_spectrum(p)`
  is the stationary allele-frequency spectrum of mutation, stabilising selection at a fixed
  optimum and drift (infinite sites, linkage equilibrium, a closed form solved for V_G);
  `equilibrium_genomes` samples genomes from it (genic variance scaled by default to
  `multilocus_equilibrium(p; drift = true).genic`, the mean genotype drawn from its
  stationary law); `equilibrium_world(p, rng; burnin, at)` starts a `World` `burnin`
  generations before the first change of selection, keeping the protocol's generation
  numbers; `simulate`/`replicate(…; init = :equilibrium, burnin)` and
  `run!(w, rng; resume = true)`. Defaults unchanged (`init = :zero`). One diploid trait,
  fixed size, viability selection, no `StandardisePhenotype`, a reaction norm only without
  a cost; anything else is an `ArgumentError`.
  - `validate/init_equilibrium.jl` (criteria written before the first run; 200 populations
    per arm against 10,000 generations from zero, then a fast and a slow continuous move):
    the burn-in needed is B = 100 in the paper's 100-locus regime and in a house-of-cards
    cell, B = 1000 with 1000 loci (the missing LD and the narrow spread of V_G take
    longer), and B = 3000 in a Gaussian-regime cell, where the default genic target is
    5% high and the excess decays slowly. On 16 published inputs B = 300 and 1000 pass.
    Wall time for a run capped 2000 generations after the move: ×5.7, ×3.7, ×2.4 and ×5.6
    less. `validate/symbolic/check_init.jl`: 40 checks.

- **Ancestry: identity by descent from the block pool, pedigrees, coalescence times,
  tree-sequence tables** (`src/ancestry.jl`, `src/pedigree.jl`; roadmap F3).
  - *Block identity, free.* Two chromosomes that hold the same block id are identical by
    descent over that block, with no mutation and no crossover inside the block since
    their common ancestor; different ids imply nothing (`shared_blocks` states this
    exactly). `ibd_fraction`, `ibd_sharing` (all pairs within individuals, within and
    between demes, per block and genome-wide, in time linear in the chromosomes) and
    `ibd_fst`. `tag_founders!` gives each founder chromosome its own blocks, so that ids
    count descent from individual founders and not from the shared block of zeros.
  - *Pedigrees, off by default.* `record_pedigree!(w)` (or of the demes of a
    metapopulation) records every offspring chromosome's parental chromosomes and
    crossovers, the mutations with their effects, and the migrations, as tree-sequence
    tables (`Pedigree`). It draws no random numbers: the run is the same bit for bit, and
    without a pedigree a generation still allocates nothing. `simplify!`
    (Kelleher et al. 2018; `simplify_every = n`) bounds memory by keeping the genealogy
    of the living chromosomes. Every allele of the living population can be rebuilt from
    the pedigree alone, bit for bit (`allele`, `alleles`).
  - *Realised coalescence times.* `pair_coalescence` counts, for every pair of living gene
    copies at every locus, the time to their most recent common ancestor, within
    individuals, within demes and between demes (`coalescence_cdf`, `coalescence_mean`,
    `tmrca`, `marginal_tree`, `pedigree_ibd`): a direct test of `deme_coalescence`.
  - *Tree-sequence export.* `tskit_tables`, `write_tskit`, `read_tskit`, `check_tables`:
    nodes, edges, individuals, populations, sites, mutations and migrations in tskit's text
    format. **Not loaded into tskit** (no Python here): the tables are checked in Julia
    against tskit's documented requirements and round-tripped through the text files.
  - Cost (`validate/coalescence.log`, N = 1000, one thread, a busy machine): recording
    without simplification takes 1.2–1.9× the time of a plain generation and 138 kB per
    generation (one mutation and one crossover per gamete); simplifying every 100
    generations 2.5–3.2×, with tables of 1.3–1.8 MB after 500 generations.
  - `validate/coalescence.jl` (criteria in its header; amended once after a 6-run smoke
    test, before the first full run, and one script bug fixed after it, both described
    there): 300 runs each of one deme, islands of 4 and 10 demes at Nm = 0.1, 1, 10 and a
    ring of 8 demes, censuses at 50, 200 and 800 generations. 424 z tests of
    E[min(T, t)] and P(T ≤ s) against `deme_coalescence` fed `migrate!`'s own law of
    pairs' origins: max |z| 2.87 (PASS). Block identity never contradicts the pedigree
    (2,700 runs); with blocks of one locus it is equivalent to identity by descent without
    mutation in all of 900 runs, and its probability matches Σ P(T) e^{−2μT} (max |z| 1.41).
- **Phenotypic plasticity: a linear reaction norm** (`Plasticity`, `OptimumSlope`, `CueLag`,
  `PlasticityCostSD`; off by default, and then the kernel's output, its allocations and
  `slim_script`'s output are unchanged). `z = g + b·ε + e` with the environment measured by
  the optimum, `θ = B·ε`, and a cue `τ` generations old, so that the plastic part of the
  phenotype is `(b/B)·θ_{t−τ}` (after Chevin, Lande and Mace 2010); a cost of plasticity
  multiplies fitness by `exp(−b²/(2ω_b²))`. The norm is expressed in `select`, where the
  optimum of the generation is set, so its fields are parameters of `select`
  (`withparams(L.select; Plasticity = 0.5)`), and demes can have their own. All three
  lowerings: the kernel, `run_moments` and `slim_script`
  (`validate/slim/compare_plasticity.jl`, no discrepancy, smallest KS p 0.26). QTSim's `:phenotype` and `:fitness` moves
  re-anchor the optimum to the mean phenotype and are rejected with plasticity.
  Theory, derived from the moments recursion: `steady_lag`, `critical_rate`,
  `red_noise_lag_var`, `expected_lag2` and `lag_load` take the slope, the optimum's slope
  and the cue's delay and reduce exactly to their old values at `b = 0`;
  `red_noise_lag_autocov`; `sustainable_rate`, the rate at which the steady lag brings the
  mean fitness to the demography's threshold. Without a cost the lag scales by `1 − b/B`
  and the sustainable rate by `1/|1 − b/B|`; with a perfect cue and no cost this is exact
  (a plastic population at rate k is, in distribution, one without plasticity at rate
  (1 − b/B)k). `validate/plasticity_rates.jl` (pre-specified, 12,320 kernel runs in a
  published-like regime): without a cost the half-survival rate rose 1.33, 2.0 and 4.0-fold
  for b/B = 0.25, 0.5, 0.75 (12/12 cells within 2% of 1/(1 − b/B)); with a cost
  (ω_b = 0.5) the gain fell to 1.2–2.1-fold at a fixed size and vanished or reversed under
  Beverton–Holt (R = 2). The deterministic prediction of the cost's effect failed 4 of 12
  cells. Its loss of V_G before the move is not the reason (`validate/plasticity_erosion.jl`,
  pre-specified, fresh seeds: 3–9% of the effect); the cost acts after the move, 2–7 times
  more strongly than `sustainable_rate`'s fitness-margin term predicts.
- **A heritable reaction norm** (`ReactionNorm = true`, `NoTraits = 2`; a third model,
  `ReactionNormModel()`): elevation and slope are the two heritable values of the model of
  several traits, the phenotype is `z = a + b·ε + e`, selection acts on `z` and the cost on
  the slope. Kernel and moments lowerings; `reaction_norm_response`, `reaction_norm_lag`.
  The design note is `docs/src/plasticity_design.md`; the algebra is checked in
  `validate/symbolic/check_plasticity.jl`.
- **Vignette 14, plasticity and sustainable rates of change**
  (`vignettes/14_plasticity_and_sustainable_rates`): where the reaction norm sits in the
  life cycle, the exact correspondence in the kernel (same seed, same survivors and V_G),
  survival against the rate from `validate/plasticity_rates.csv` (read from the cache, not
  rerun) with the P2 failures, the cost trade-off from `sustainable_rate`, a heritable
  slope (kernel against moments, 8 runs per scenario), and a delayed cue under red noise
  (illustrative, 54 runs). About 45 s of code on three threads; run in CI's `vignettes` job.
- **Locally favoured alleles in metapopulations: `run_moments(...; loci = :rare_local)`**
  (`src/local_alleles.jl`). The rare-alleles closure takes every mutant to be deleterious
  wherever it is; under local adaptation with migration, demes lag behind different
  optima, and a rare mutant is favoured where its deme lags. To second order, its lineage
  grows at `π·σ + σᵀWσ` with `W = Π(I − P)^#P` (`P` the backward migration matrix), so
  selection against it weakens by `2Λ`, `Λ = Σ W_dc (ℓ_d − ℓ̄)(ℓ_c − ℓ̄)ᵀ`, `ℓ_d = s_d⁻¹(θ_d −
  z̄_d)` (`QTLSim.local_gain`; island model `(1 − m′)/m′·Var(ℓ)`). The closure's rare loss
  uses `[s⁻¹ − 2Λ]₊` in place of `s⁻¹`. No constant is fitted, and it equals
  `:rare_alleles` exactly for one deme, for a common optimum and without migration.
  `validate/local_adaptation.jl` (run 1, criteria written before any kernel run): the
  current default misses the kernel's within-deme variance by −45% (one trait, house of
  cards, migration 0.01) and by −148%/−85% (two traits, D1); `:rare_local` is within 0.8%
  and within 1–9%, passes the 15% accuracy criterion in all 10 cells, and is never worse
  than the multilocus Gaussian by more than 0.001. Its one failure, on the pre-specified
  rule "better than `:multilocus` in every house-of-cards cell", is a tie: where the gain
  switches the rare loss off, `:rare_local` is the multilocus model (0.008 against 0.007).
- `validate/symbolic/` (own environment with Symbolics.jl; CI job `Symbolic checks`, about
  30 s): checks the algebra of the closures and evaluates each identity against the code.
  Lean: `CategoricalQuantGen/LocalAdaptation.lean`.
- **Genetic variance under a moving optimum: `run_moments(p; loci = :spectrum)`**
  (`src/spectrum.jl`; one deme, one diploid trait, drift always on). The other moments
  models predict the genetic variance before the optimum moves but hold it roughly fixed
  afterwards, so they could not predict extinction at slow rates. This model follows the
  expected allele-frequency spectrum of the segregating mutations (a Poisson random field):
  selection moves each allele by Hayward and Sella's (2022) eq. 7 with the rare-alleles
  closure's Gaussian term, drift acts among the expected parents, and Barton's (1995)
  linked-selection factor lowers the effective size under directional selection. No
  constant is fitted. Directional selection sweeps or purges the alleles that carry the
  variance within a few hundred generations; the variance sweeps sustain scales with the
  survivors, and QTSim's single crossover per gamete adds linked selection, so a lagging
  population loses variance, lags further and dies. `spectrum_equilibrium` and `spectrum_vg`
  give its stationary state.
- `validate/vg_move.jl` (run 1, criteria written before any kernel run: 148 PASS, 17 FAIL,
  the failures kept):
  - It predicts extinction times in density-regulated populations at slow rates within
    0–20% of the kernel, where the rare-alleles model predicted survival or times up to
    1.7× too long. Median |log error| is 0.038, against 0.106.
  - It matches the published extinction times on 240 of the paper's inputs (median
    table/prediction 1.000).
  - It reproduces the rise of the variance when rare alleles of large effect sweep.
  - Caveat: the model was designed after studying vignette 10's runs and exploratory
    kernel runs of the same scenarios. The validation used fresh seeds, but its vignette 10
    cells are not independent of the design.
  - Where it fails:
    - It over-predicts the rise of the variance at intermediate rates in the
      house-of-cards regime.
    - In the Gaussian regime its variance is 1.5× too high, and it misses the kernel's
      collapse of the variance under a move.
    - The rate at which linked-selection associations decay is the open problem.
  - It is not the default. It is recommended for moving optima when 4Nμ ≲ 1 per locus,
    as in the paper's inputs and vignettes 10 and 11.

### Changed

- **Tests split by topic** (roadmap F11). `test/runtests.jl` is now a header and a runner; the
  75 testsets moved unchanged into 17 files `test/<topic>.jl`. Run a subset with
  `Pkg.test(test_args = ["traits", "demes"])`, as arguments to `test/runtests.jl`, or in
  `ENV["QTLSIM_TEST"]`; Aqua runs only with the full suite. Every topic passes alone.

- **Several traits: `loci = :auto` resolves to `:rare_local`** (was `:multilocus`), by the
  rule fixed before `validate/local_adaptation.jl` ran: the closure passes the accuracy and
  "no worse, better in the house of cards" criteria in D1 and D2, and for one population it
  is `:rare_alleles`, which already passed `validate/rare_alleles_traits.jl`'s T cells
  (up to 1000× closer to the kernel than the multilocus Gaussian). Pass `loci = :multilocus`
  for the old default.
- `loci = :rare_alleles` works for several traits (`rare_alleles_G`, `house_of_gauss_G`,
  `house_of_cards_G`) and for haploid clones, and `:auto` uses it for one-trait haploids.
  Several traits stay on `:multilocus` by default. The closure fails a pre-specified test
  of four demes with local optima in the house-of-cards regime: there it predicts under
  half the kernel's within-deme variance (`validate/rare_alleles_traits.jl`).

### Fixed

- Three gaps the support matrix found: `QTLSim.moments_loci` (so `multilocus_equilibrium`'s default) of a heritable reaction norm is `:multilocus`, as `run_moments` uses (was `:rare_local`); `loci = :spectrum` with several traits in several demes (or `run_trait_moments` directly) throws an `ArgumentError` instead of running the rare-alleles closure; `slim_script` checks the parameters (`checkparams`) before generating code, so two traits with `Plasticity` (or `StandardisePhenotype`) are refused cleanly instead of failing in the Eidos splice.

## [0.1.0] - unreleased

First release. QTLSim.jl is a Julia port of QTSim, the C simulation (Andrew Meade,
University of Reading) behind Pagel, Gardner and Meade, *Scientific Reports* (2025),
with a categorical description of its life cycle that lowers to the Julia kernel, to
SLiM and to deterministic moment equations, and several models that QTSim does not have.

### Kernel: the QTSim port

- `Params` (QTSim's options and defaults, a few renamed; `checkparams`), `World`,
  `simulate`, `run!`, `generation!` and the step operators (`develop!`, `select!`,
  `moveoptimum!`, `reproduce!`, `mutate!`, `gamete!`, `swapbuffers!`) reproduce QTSim's
  generation in QTSim's order: viability selection with QTSim's `0.4/ω` normalisation,
  random mating of survivors, Poisson recombination, Poisson–Gaussian mutation, and the
  `:phenotype`, `:fitness` and `:continuous` moves of the optimum.
- `CrossoverModel = :qtsim` (the default) reproduces QTSim's recombination, in which a
  gamete has at most one crossover; `:all` applies every crossover point.
- Genomes are shared, reference-counted blocks of about √`NoQTL` loci with cached sums,
  copied on write. One thread, 1000 diploids, before the move (`validate/benchmark.log`):

  | `NoQTL` | kernel (ms/gen) | dense array | QTSim | QTSim / kernel |
  |--:|--:|--:|--:|--:|
  | 1 | 0.10 | 0.12 | 0.21 | 2.0 |
  | 100 | 0.25 | 0.25 | 0.47 | 1.9 |
  | 500 | 0.33 | 0.47 | 1.41 | 4.3 |
  | 1000 | 0.39 | 1.10 | 2.86 | 7.3 |
  | 5000 | 0.78 | 4.95 | 13.2 | 16.9 |

  Allele values are bitwise those of the dense kernel for the same seed. `genomes(w)` and
  `setgenomes!(w, G)` read and write a dense `(NoQTL, Ploidy, N)` array.
- `read_qtsim_options` and `write_qtsim_options` read and write QTSim's options files.
- QTSim's C source is bundled unmodified in `reference/QTSim` (checked by `SHA256SUMS`),
  with a Makefile to build it.

### Validation against QTSim

- **Statistical**: `validate/compare_to_qtsim.jl`, Welch t-tests on equilibrium and
  directional scenarios; vignette 08 repeats it with KS tests.
- **Exact**: `validate/gsl_locked.jl`, a replica drawing QTSim's GSL random numbers in
  QTSim's order, is bit-identical to QTSim every generation, for one locus and for
  multi-locus inputs through extinction, with a negative control (`:all` diverges).
- **Published results**: `validate/published_sweep.jl` against the paper's per-run
  table. On two runs of about 590 inputs, the Julia/published ratio of time to extinction
  is 1.008 [0.990, 1.026] and 0.999 [0.993, 1.006], and of pre-move genetic variance
  1.007 [0.986, 1.029] and 0.999 [0.979, 1.020]; `CrossoverModel = :all` is rejected
  (extinction time ×1.12, V_G ×1.06) (`validate/published_sweep.log`).
- The SLiM model generated from the life cycle matches a hand-written SLiM port run for
  run and QTSim in distribution (`validate/slim/valid_generated.log`). A pre-specified
  test of 40,000 runs per program found mean V_G differences of −0.56% (SLiM − QTSim,
  p = 0.57) and +0.26% (kernel − QTSim, p = 0.79) (`validate/slim/three_way_confirm.log`).

### The categorical layer (package extension `QTLSimCatlabExt`, `import Catlab`)

- `lifecycle()`: the life cycle as a morphism of a free symmetric monoidal category
  presented in Catlab, with stages (Zygotes, Adults, Survivors, Offspring), steps
  (`develop`, `select`, `reproduce`, `record`, `rescale`, `mutate`, `swap`), `⋅` for
  composition, `⊗` for demes side by side and `migrate(k)` coupling them.
- Three lowerings by functors, `lower(f, target)`:
  - `:kernel`: `run_composed` runs any composite; QTSim's generation is bit-identical to
    `run!`. `operators()` and `operator(name, act)` give the kernel's operators directly.
  - `:eidos`: `slim_script` writes a nonWF SLiM model (any number of loci, one
    subpopulation per deme, per-deme parameters, migration matrices, density regulation,
    fecundity selection, several traits), validated against the kernel in
    `validate/slim/valid_*.log` (smallest KS p 0.029 to 0.263 across the checks, no
    discrepancy).
  - `:moments`: `run_moments` iterates the deterministic moments of a Gaussian
    population for any life cycle, with k demes, migration, drift (`erosion`,
    `between = true`) and a realised optimum sequence.
- `withparams(f; …)` puts parameters on steps; `L.migrate(k; rate)` carries a rate or a
  matrix; all three lowerings honour them. `params ⊗ demes(k)` and
  `[p1, …, pk] ⊗ demes(k)` build metapopulations.
- The population as an ACSet: `population_schema`, `lift`, `lower!`, `genotypes`,
  `phenotypes`, `genome_from_acset`.

### Models

- `OneTraitModel()` (QTSim's) and `TraitModel(n)` (n traits with pleiotropy) are two
  families of lowerings of the same life cycle; `model_of(p)` picks one from `Params`.
- Step parameters are typed by the model (`step_schema`): `θ` is a number for one trait
  and a vector for n traits, and `Ω` exists only for n traits.

### Biology beyond QTSim (Julia only, off by default)

- **Demography**: `Demography = :beverton_holt` with `GrowthRate` R (offspring
  Poisson(R n/(1+(R−1)n/K))); records carry `mean_fitness`, `sel_diff` and `popsize`.
  `validate/demography_sweep.jl` reruns 1200 published inputs under density regulation
  (`validate/demography_sweep_summary.md`).
- **Fecundity selection**: `Selection = :fecundity` (everyone survives, parents drawn with
  probability ∝ W), in the kernel, Eidos and moments (`validate/fecundity_summary.md`).
- **Fluctuating optima**: `MoveType = :periodic` and `:red_noise`.
- **Metapopulations and migration**: `island_migration`, `stepping_stone`, `migrate!`,
  `migration`; demes are visited in a random order, so migration is exchangeable.
- **Several traits**: `NoTraits = n` with mutational and selection covariances
  (`MutCorr`, `SelCorr`, `MutCov`, `SelCov`), directional and fluctuating moves along a
  vector (`MoveDirection`, `FluctDirection`, `MoveNoiseCov`); see the
  [design note](docs/src/traits_design.md).

### Theory

- Gaussian (Lande, Bürger) equilibria: `gaussian_vg`, `gaussian_vg_exact`,
  `selection_variance`, `effective_size`.
- Multilocus Gaussian moments with linkage disequilibrium and recombination:
  `gaussian_vg_loci`, `multilocus_equilibrium` (`validate/moments_loci.log`).
- Rare-alleles closure (`rare_alleles_vg`, `house_of_gauss_vg`), the moments default for
  one diploid trait; `house_of_cards_vg`. Across 60 cells of the Gaussian,
  house-of-cards and drift regimes, 122 pre-specified criteria pass and 4 fail
  (`validate/moments_regimes.log`).
- Lag theory: `steady_lag`, `critical_rate`, `mutation_adjusted`, `periodic_response`,
  `red_noise_lag_var`, `expected_lag2`, `lag_load`.
- Extinction and rescue: `extinction_distribution` (exact size chain, branching or
  binomial hazard), `rescue_probability` (Gomulkiewicz–Holt).
- Neutral differentiation among demes: `qst_island`, `island_variances`,
  `backward_migration`, `deme_coalescence`, `exchange_origins`, `neutral_deme_variances`,
  `fst_stepping_stone`, `pair_recombination`; 41 comparisons with simulation, largest
  |z| 2.57 (`validate/fst.log`).
- Several traits: `gaussian_G`, `gaussian_G_loci`, `selection_gradient`,
  `trait_response`, `steady_lag_vector`, `gmax`, `response_deflection` and the
  fluctuating-optimum covariances (`validate/traits_kernel_vs_moments.log`).

### Sweeps

- `sweep(f, jobs; cache, cost, ntasks)`: greedy threaded jobs with a resumable CSV cache,
  deterministic by seed; `replicate(p, n; seed, summary)`.

### Lean library

- The mathematical facts behind the decomposition (heritability algebra, stabilising
  selection, mutational input, unbiased recombination, composition and iteration,
  lowering from step images, demes and migration, monoidal naturality, the Gaussian
  moment recursion and more) are proved in Lean 4 with Mathlib in the companion library
  `CategoricalQuantGen`, outside this repository. Docstrings cite it as
  `lean/CategoricalQuantGen/…`.

### Vignettes

Thirteen Quarto vignettes, rendered into the documentation: getting started;
mutation–selection balance; a moving optimum; loci and recombination; the categorical
layer; recomposing the model; one life cycle, two backends; validation against QTSim;
replicating Pagel, Gardner and Meade (2025); demography and rescue; density regulation in
the paper's populations; metapopulations; several traits and pleiotropy.

### Continuous integration

- Tests on Julia 1.10 and 1 (Ubuntu) and 1 (macOS), with Aqua.
- QTSim parity: checks the bundled source's sums, builds it, runs `gsl_locked.jl` and
  `compare_to_qtsim.jl`.
- SLiM backend: builds SLiM 5.2, runs `check_generated.jl`, `compare_to_qtsim.jl`,
  `compare_demography.jl` and `compare_selection.jl` (more runs weekly).
- Vignette code (01–03, 05, 06, 10, 12, 13 on push; 04 weekly), the vignette sync check,
  and the Documenter site deployed to https://ecorecipes.github.io/QTLSim.jl/.

### Behaviour to know about

- `run_moments` and `multilocus_equilibrium` default to `loci = :auto`: the rare-alleles
  closure for one diploid trait, the multilocus model for several traits and for
  haploids. Pass `loci = :multilocus` or `loci = :aggregate` for the others.
- `migration_rate` is the forward rate at which an individual starts an exchange; since a
  migrant swaps places with an individual of another deme, the backward fraction of a
  deme's adults born elsewhere is about twice the rate (`backward_migration`).
- `Params` compares and hashes by value (`==`, `hash`), so equal parameter sets are equal
  dictionary keys.
- QTSim's literal `0.4` in the fitness normalisation (`W(θ) ≈ 0.9974`) and its
  single-crossover recombination are reproduced on purpose.
