# Pre-specified criteria: the haplotype-block lineage emulator

Written 2026-11-16, after development (validate/blocks_dev.jl, validate/blocks_dev.log) and before
any comparison on fresh data. Committed before validate/blocks.jl is run on fresh inputs (its
`select` mode, which only draws the samples and computes the pool's `:rare_alleles` k50s, may run
first; no outcome is looked at). The script is validate/blocks.jl; its output goes to
validate/blocks.log. No default changes whatever the outcome.

## The model under test

`emulate(p; loci = :blocks, blocks = 160, within = :gaussian)` and `run_blocks` (src/block_lineages.jl):

- The kernel's life cycle on `PopSize` diploid individuals (QTSim's fitness with its 0.4,
  viability selection, two distinct random survivors as parents, Poisson–Gaussian mutation with
  the kernel's variance, the continuous move of the optimum).
- Each chromosome is reduced to the sums y_b of B = min(160, L) blocks of consecutive loci,
  boundaries at equal steps of QTSim's genetic map (`partition = :map`; the law of the single
  switch point, `QTLSim.qtsim_switch_law`, P(c) = Σ_K Pois(K; u) C(M−c−1, K−1)/C(M, K), checked in
  validate/symbolic/check_blocks.jl).
- A gamete takes the kernel's switch point from its exact law. If the point is on a block
  boundary, the blocks are copied from the two strands. If it falls inside block b, leaving a
  fraction f of the block's loci on the starting strand, block b is drawn by Gaussian segregation:
  y = ȳ_b + f(y_A − ȳ_b) + (1 − f)(y_B − ȳ_b) + √(2f(1 − f)V_b) ξ, with ȳ_b and V_b the mean and
  variance of y_b over the parents' generation (all 2N chromosomes).
- Mutations: Poisson(2N·L·μ) per generation, each on a uniform chromosome and a uniform locus,
  added to that locus's block.
- Start: the kernel's equilibrium sampler (`equilibrium_genomes`'s law) summed over the blocks.

With L ≤ 160 every locus is its own block, every switch is on a boundary, and the model is the
kernel in law (the 100-locus strata test the implementation more than the reduction). Nothing is
fitted to kernel outcomes. B, the partition and the within-block law were chosen in development
by the convergence rule written in validate/blocks_dev.log before step 1 (B → 2B: median
|log V_G ratio| < 0.03 before the move and < 0.05 at +100 and +200, mean |ΔP| ≤ 0.10, ≤ 1 call of
21 changed); snapping crossovers to boundaries met it at no B ≤ 320 and was dropped.

## Disclosure: what development used

- Inputs: the spectrum diagnosis's D (13) and Dctl (8) sets, their 10 raw kernel trajectories and
  40 E6′ kernel runs per input, all seen by the earlier studies; 400 + 80 fresh kernel runs on two
  100-locus D inputs for the implementation check (step 0b).
- The comparisons with the kernel in development (step 2: V_G at B = 160 with the Gaussian law at
  0.95/0.96/0.92/0.93 of the kernel's before/+50/+100/+200; all 21 threshold calls agree) were
  made after the convergence rule had been written but before the choice of the within-block law
  was final: the choice of the Gaussian law over snapping was made on convergence (step 1b),
  and the kernel comparison agrees with it; B was set by the rule alone.
- Known before validation: V_G at B = 160 is 5–8% below the kernel's on the development inputs
  (more blocks → more variance, monotonically), and the cost is not much below the kernel's
  (validate/blocks_dev.log step 4).

## Samples (fresh)

- **Exclusions.** spectrum_ld.jl's `all_exclusions` (every input of validate/emulator_spectrum.jl's
  exclusions, E1–E4 sample and E6′ pool; the spectrum diagnosis's diagnostic and fresh inputs;
  validate/published_sweep_results.csv; the table rows rerun in vignettes 09 and 11), united with
  every input id that appears in any file name or in any .csv/.log/.txt/.md/.jl file under
  `validate/` of the main package and of every `wt-*` worktree except wt-paper (this covers the
  D, F, G, H′ and J sets of the linked-selection chain through their raw-run file names and logs,
  spectrum_ld's H pool, H and F, spectrum_beta's inputs, emulator_variance's samples). Files that
  name more than 50,000 ids are not used: they are the deterministic `:rare_alleles` emulator run
  on the whole table (validate/emulator.jl), which no kernel comparison of this design used. The
  script prints the counts and checks the samples' disjointness from the exclusions by set
  intersection (must be 0).
- **H (primary).** E6′'s procedure on a new pool: per stratum, a MersenneTwister(20261131)
  permutation of the non-excluded inputs, prefiltered by 0.3 ≤ MoveFitnessSD/kc ≤ 3
  (kc = `drift_critical_rate`), 100 per stratum, their `:rare_alleles` k50
  (`emulated_half_survival_rate`, burn-in and horizon 2000); picks with 0.3 ≤ r = MoveFitnessSD/k50
  ≤ 1.5, 3 per log-bin of r (4 bins) and filled to **15 per stratum** from the rest,
  MersenneTwister(20261132): **120 inputs**, stratified by locus count and Poisson mean. No table
  outcome is used in the selection.
  - Kernel: 40 runs per input, E6′'s protocol (init = :equilibrium, burn-in 1000, to generation
    12000, the move at 10000), fresh seeds 41261116000 + 100j + r. Recorded: alive at +2000; V_G
    (mean over −200…−1, and over ±5 generations around +50, +100, +200, runs alive throughout);
    the per-locus genic variance and the block-genic variance at B = 160 (genomes summed over the
    model's blocks; 2Σ_b Var y_b) at −200, −100, −1 (averaged), +50, +100, +200.
  - Model: 40 `run_blocks` runs per input with the same protocol (`within = :gaussian`,
    `blocks = 160`, burn-in 1000), seeds `Xoshiro(job_seed(20261135, "id#r"))`, the same
    statistics. P2000 = the fraction of these runs alive at +2000. (`emulate(p; loci = :blocks,
    burnin = 1000)` is the same model and protocol with generation numbers rebased; the 40 runs
    are used once for P2000 and the trajectories, to stay in budget.)
- **F (secondary).** 40 per stratum (320), MersenneTwister(20261133) permutation, outside the
  exclusions and H's pool. Model: `emulate_sweep(…; loci = :blocks, blocks = 160,
  within = :gaussian, replicates = 10, seed = 20261134, burnin = 2000, horizon = 2000)` (replicate r
  of input id from `Xoshiro(job_seed(job_seed(20261134, "F|id"), r))`).

## Tests

**Primary (H)**, E6′'s statistics as validate/spectrum_ld_criteria.md (c):

- **(P1)** E6′c: accuracy of P2000 > 0.5 for the kernel's majority alive (f > 0.5) ≥ 0.85 and
  AUC(P2000; f > 0.5) ≥ 0.85.
- **(P2)** pooled x50 ratio (emulator/kernel, `x50_rule` on all H) in [0.8, 1.25].
- **(P3)** x50 per stratum: ratio in [0.8, 1.25] in every stratum (E6′d's rule).
- **(P4)** paired comparison with `:spectrum, ld = :dynamic` on the same inputs: the number of
  inputs whose call the block model gets right and ld = :dynamic wrong is at least the number the
  other way (exact McNemar p reported). PASS = not worse.
- E6′a/b are printed (no verdict). Comparators on the same inputs, printed with the same
  statistics and no verdict: the default `:spectrum`, `ld = :dynamic` (deterministic, burn-in and
  horizon 2000, as in spectrum_ld) and `:rare_alleles`.

**Secondary.**

- **(S1)** V_G, model/kernel per input, median over inputs: before the move in [0.95, 1.05]; at
  +50, +100, +200 (inputs with ≥ 10 of 40 runs alive through the window in both) each in
  [0.9, 1.1].
- **(S2)** block-genic variance relative to its pre-move value, model/kernel, median at +50,
  +100, +200, each in [0.9, 1.1] (the fall of within-block variance under the move at the model's
  resolution; the kernel's per-locus genic ratio is printed beside it, no verdict, since the model
  has no per-locus state).
- **(S3)** on F: E1–E4 of validate/emulator_criteria.md as validate/emulator_spectrum.jl applies
  them (E3 with the median T of the replicates; E5 not computed). The default `:spectrum` and
  `:rare_alleles` on the same inputs as comparators (no verdict).
- **(S4)** cost (mode `cost`, one thread, 16 inputs of F, 2 per stratum): median wall time of one
  emulator call (10 replicates, burn-in 2000, horizon 2000) ≤ that of `:spectrum` with
  ld = :dynamic. Also printed: the default `:spectrum`'s time, and one run of the kernel's protocol
  in the block model and in the kernel (ratio by locus count).

Expected difficulties, stated now: the 5–8% V_G deficit at B = 160 may fail (S1) before the move;
(S3)'s E2 is coarse with 10 replicates; (S4) is expected to fail against the default
`:spectrum` and is uncertain against ld = :dynamic.

## Budget

About 15 CPU-hours for the task; development used about 5. Estimated: H pool k50s 0.5 CPU-h,
H kernel ≈ 3.5, H model ≈ 2, comparators on H ≈ 0.7, F model ≈ 2, F comparators ≈ 0.5. If a
timing pilot (mode `pilot`, 16 inputs of F, no outcome) shows that F would exceed 3 CPU-h, F is
shrunk per stratum (same permutation) by an amendment recorded here before any result.

## Amendments

(none)
