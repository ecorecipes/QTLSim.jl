# News

All notable changes to QTLSim.jl are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the package follows
[semantic versioning](https://semver.org/) (before 1.0, a minor release may break the API).

## [Unreleased]

### Added

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
