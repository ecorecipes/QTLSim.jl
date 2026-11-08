# Several traits: a second model through the same life cycle

QTLSim's categorical layer writes QTSim's life cycle once, as a morphism of a free
symmetric monoidal category (`lifecycle()`), and lowers it by functors to the kernel, to
SLiM and to the moments of a Gaussian population. The north star is a framework of which
QTSim is one instance. This note records what happened when a second model went through
the same syntax and lowerings: **n quantitative traits with pleiotropy** (Lande 1979,
1980; Jones, Arnold and Bürger 2003, 2007). It is a test of whether the decomposition
generalises.

## The model

- `NoTraits = n` traits (any `n ≥ 1`; `n = 1` is QTSim). Each locus carries a vector of
  `n` allelic effects, and the genotype `g ∈ ℝⁿ` is their sum.
- The phenotype is `z = g + e`, with `e ~ N(0, V_E I)`: every trait has QTSim's
  `V_E = (1 − h²) V_P`, independently.
- A mutation adds a draw from `N(0, M/(Ploidy·NoQTL·rate))` to a uniformly chosen locus,
  so that the expected input per individual per generation is the mutational covariance
  `M`: `m·V_E` on the diagonal and `r_μ·m·V_E` off it (`MutCorr`), or `MutCov`.
- Fitness is QTSim's Gaussian in `n` dimensions,
  `W(z) = c·exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ))`, with QTSim's `c = 1/(0.4√(2π)) ≈ 0.9974`, so that
  `n = 1` is QTSim's fitness exactly. `Ω = ω² S`: the scalar width `ω` (initial width,
  changed by `rescale` as in QTSim) times a shape `S` with unit diagonal and correlation
  `r_ω` (`SelCorr`), or `SelCov/ω₀²`. Selection is viability (or fecundity) selection, as
  for one trait.
- The optimum `θ ∈ ℝⁿ` moves along a unit direction `u` (`MoveDirection`, trait 1 by
  default): `θ += d·u` (`:continuous`), `θ = z̄ + √V_P·d·u` (`:phenotype`) or
  `θ = z̄ + ω·d·u` (`:fitness`). Selection on one trait then acts on the others through
  `G`'s and `Ω`'s correlations (Lande and Arnold 1983).

## Findings: the syntax did not change

**`lifecycle()` is unchanged.** The stages (Zygotes, Adults, Survivors, Offspring) and
the steps (develop, select, reproduce, record, rescale, mutate, swap) are the same
generators, QTSim's generation is the same composite, `migrate(k)` and `⊗` are the same.
No new generator was needed. Every lowering of the trait model is a functor out of the
*same* free category; what differs is where it sends the generators:

| | one trait (QTSim) | n traits |
|---|---|---|
| object image (an individual at a stage) | `g, e, z ∈ ℝ` | `g, e, z ∈ ℝⁿ` |
| `develop` | `z = Σ effects + e` | the same, per trait |
| `select` | `W = c·exp(−(z−θ)²/2ω²)`, θ moves by a rule | `W = c·exp(−½ (z−θ)ᵀΩ⁻¹(z−θ))`, θ moves along `u` |
| `reproduce` | gametes copy blocks of loci | the same code: a block holds `n` values per locus |
| `record` | QTSim's columns | the same, plus trait means, `G`, `P`, lag and selection differential vectors |
| `rescale` | `ω ← ChangeFitnessSD` | the same: `Ω = ω² S` scales with it |
| `mutate` | `+N(0, σ²)` at a random locus | `+N(0, M/(2L·rate))` at a random locus |
| `swap` | swap buffers | the same, with the trait environments |

`reproduce`, `rescale` and `swap` have literally the same images: recombination and
inheritance do not care what an allele is, only that it can be copied. `migrate`'s
image is the same too (a migrant's genome is block ids; only its environment vector had
to move with it). So the functor table changes only in `develop`, `select`, `mutate` and
`record`, and only because the object image changed from ℝ to ℝⁿ.

**The ACSet schema is unchanged too.** `SchPopulation` (Individual, Allele, `owner`,
attributes `effect`, `env`, `genotype`, `phenotype`, `fitness`) is the same presentation;
a population of `n` traits is an ACSet on it with vector-valued attribute types,
`PopulationACSet{Vector{Float64},Vector{Float64}}` (`TraitPopulation`). The genotype is
still the sum over the preimage `owner⁻¹(i)`, now in ℝⁿ. In categorical terms the model
changed the attribute functor, not the schema.

In the language of the north star: the n-trait model is a family of functors
`F_n : LifeCycle → Kern`, indexed by the object image `ℝⁿ`, with `F_1` QTSim's. The
lowerings are fixed by their images on the generators (`Lowering.lean`), so `F_n`'s
generation is the composite of its steps' images by the same theorem, with no new proof
obligation about the wiring. The new obligations are about the images: that `mutate`
injects `M` (the multivariate `mutationalVarianceInput_eq`) and the Gaussian recursion's
fixed point (see "For Lean" below).

## Where it did not work that way

These are the places where the generalisation was not free, as first found. Items 1–3
have since been addressed by making the model explicit (next section); the notes record
the original finding.

1. **Step parameters belong to the model, not to the syntax.** A generator carries
   parameters (`withparams(select; ω = 0.3)`, the table `STEP_PARAMS`). The bare
   presentation is model-independent, but the parameter vocabulary of each step is not:
   the trait model's `select` reads `SelCorr`/`SelCov`/`MoveDirection`, its `mutate`
   reads `MutCorr`/`MutCov`, and a fixed optimum `θ` is a vector. The parametrised
   generators `Step(name, params)` range over a model-specific parameter space. At first
   steps with parameters took the scalar parameters only. *Now addressed*: each model
   types its steps' parameters (`step_schema`), a vector `θ` and a matrix `Ω` included.
2. **The optimum's move is not a function of the syntax.** QTSim's moves are rules about a
   scalar; with several traits they need a direction, a new parameter of `select`'s image.
   The fluctuating moves (`:red_noise`, `:periodic`) are not generalised (a direction for
   the trend and a covariance for the noise would be needed) and `checkparams` rejects
   them with `NoTraits > 1`, as it does `StandardisePhenotype`.
3. **The kernel's images are one table with a branch, not two tables.** In the kernel the
   object image is a `World`, and a `World` with `traits` set is the ℝⁿ image; the operators
   `develop!`, `select!` and `mutate!` branch on it (and call `develop_traits!` etc.). So
   `lower(f, :kernel)` is one table whose images dispatch on the object image, rather than
   a second table. That is a Julia convenience: semantically they are two functors, and
   the trait moments lowering is written as a genuinely second table
   (`trait_moment_images`, passed as `images`) to show that the syntax is shared.
   *Now addressed*: `kernel_images(model)` is two tables, chosen by `lower(f, :kernel;
   model)`.
4. **QTSim's scalar outputs are a projection.** `record`'s QTSim columns (`pheno_mean`,
   `geno_var`, `theta`) describe trait 1; the trait columns are additional. The SLiM
   model's output format, the golden file and the parsers are keyed to the scalar columns.

## Models: step parameters and images belong to the model

The three findings above have one cause: the model (one trait or `n`) was implicit, read
off the `World` at run time, while step parameters were checked against one global table.
`src/models.jl` makes it explicit and keeps it minimal and concrete:

- `QTLSim.Model`, with two instances, `OneTraitModel()` (QTSim's) and `TraitModel(n)`;
  `model_of(p)` is the model of a `Params`.
- **Each model declares its own parameter schema per step** (`step_schema(model)`, from
  `step_fields`, the fields of `Params` a step's image reads, and `step_own_params`, the
  parameters that are not fields, each with a checker). For one trait, `select` takes a
  number `θ` and `ω` (and `StandardisePhenotype`); for `n` traits, `select` takes a
  vector `θ` of `n`, a fitness covariance matrix `Ω` (in place of `ω² S`), `ω`, `SelCorr`
  and `MoveDirection`, and `mutate` also takes `MutCorr` and `MutCov`. `STEP_PARAMS`
  survives as `step_schema(OneTraitModel())`.
- **`withparams` validates against a model**: `withparams(f; kw...)` against QTSim's,
  `withparams(f, model; kw...)` or `withparams(f, p; kw...)` against `model` or
  `model_of(p)`. **`lower` validates again against the model it targets** (`lower(f,
  target; model)`), so `withparams(g; θ = 0.1)` lowered for two traits, or `withparams(g,
  p2; θ = [0.1, 0.0])` lowered for one, is an `ArgumentError` at lowering, not a wrong run.
  The parametrised generators `Step(name, params)` are unchanged: the presentation stays
  model-independent, and the parameters are typed where they are interpreted.
- **Each model has its own table of images in each backend**: `kernel_images(model, rec)`,
  `moment_images(model, rec)` and `eidos_images(model)`, chosen by `lower(f, target;
  model)`. The two kernel tables share the images of `reproduce`, `record`, `rescale` and
  `swap` (`shared_kernel_images`) and differ in `develop`, `select` and `mutate`, which
  call `develop_one!`/`develop_traits!` etc. directly; each image refuses a population of
  the other model. `develop!`, `select!` and `mutate!` remain as the monolithic `run!`'s
  dispatchers. So `lower(f, :kernel; model = TraitModel(2))` is visibly a second functor
  from the same syntax. `run_composed`, `run_moments` and `slim_script` pass `model_of(p)`;
  a metapopulation must be of one model.
- **How a step's parameters act is per model** (`with_step_params(model, s, op)` for the
  kernel and the moments, `eidos_with_params(model, s, code)` for Eidos). For `n` traits a
  step's `θ` replaces the optimum vector, `Ω` the width (1) and the fitness shape, and the
  `Params` fields that define the state's derived matrices replace them for the step: the
  inverse shape `S⁻¹` (`SelCorr`), one mutation's covariance (`MutationScalar`,
  `MutationRatePerGamete`, `MutCorr`, `MutCov`) and the move direction. In Eidos they are
  literals (`THETA = c(…)`, `SINV`, `UDIR`, `MUTCOV`), so **per-deme trait parameters in a
  SLiM model go through steps**: `withparams(g, p; θ = [0.1, 0.0]) ⊗ withparams(g, p; Ω =
  …) ⋅ L.migrate(2)` gives each deme its own optimum vector or fitness covariance, while
  `slim_script(f, [p1, p2])` still requires the demes' trait fields of `Params` to agree.

What this shows about generalisation:

- The syntax needed no change; the change was to make the *codomain side* of each functor
  a first-class value. A model is exactly the data of a family of functors from the same
  free category: an object image, a table of step images per backend, and the parameter
  space of each parametrised generator. That is the shape a general framework would have
  (a model as an instance), and it is what a third model would implement.
- Parameters are not uniform across models even for the same step: `select`'s `θ` lives in
  the object image (ℝ or ℝⁿ), `Ω` exists only where the object image has more than one
  dimension, and `StandardisePhenotype` has no n-trait meaning. A step's parameter space
  is a function of the model, which is why one global table could not work.
- Typing at lowering, not at construction, keeps the presentation free: the same
  `withparams(L.select; ω = 0.3)` is valid in both models, since `ω` is in both schemas.
- It found a silent bug: before, a step's `MutationScalar` (or `MutationRatePerGamete`) on
  a population of traits was ignored in the kernel and in the moments, because the trait
  state's mutation matrix was derived from the deme's `Params` once, not from the step's.
  With the derived matrices part of what a step's parameters replace, a step's mutation
  fields now give the run of `Params` with those values, bit for bit (tested).
- One trait is bitwise unchanged: 36 hashes before and after (`simulate` in three regimes,
  `run_composed` with and without step parameters, two demes with and without, Beverton–Holt,
  fecundity, red noise, periodic, `run_moments` in six configurations including two demes
  with `between = true`, trait runs, and 13 `slim_script` outputs), the golden file and
  `check_generated.jl`. A generation still allocates nothing, now also with a step that
  overrides fields of `Params` (an `IdDict` keyed by the immutable `Params` boxed the key,
  272 bytes per generation; a typed vector searched by `===` does not).

What remains: `record` is shared (its trait columns are added when the `World` carries
traits) rather than per model; `slim_script`'s per-deme `Params` cannot vary the trait
fields (use steps); `checkparams` still knows both models (`checktraits`); and a model is
chosen from `NoTraits`, not declared, so `Params` remains the union of both models'
parameters. A general framework would let a model own its `Params`.

## What is implemented

- **Params** (`src/params.jl`): `NoTraits`, `MutCorr`, `SelCorr`, `MutCov`, `SelCov`,
  `MoveDirection`, checked by `checkparams`; `mutational_covariance`, `selection_shape`,
  `selection_covariance`, `move_direction`. Params now compares by value (it holds
  arrays).
- **Kernel** (`src/blocks.jl`, `src/model.jl`, `src/traits.jl`): one block pool per
  population as before, a block holding `n` values per locus (locus by locus) and caching
  `n` sums, so that a gamete still copies contiguous runs and shares blocks; with `n = 1`
  the layout and every operation are the old ones. `World.traits` (a `TraitState`) holds
  the `n × N` genotypes, phenotypes and environments, the optimum vector and the
  precomputed `S⁻¹` and `M^{1/2}`. QTSim's one-trait path is untouched: runs are bitwise
  identical to before (12 runs in four regimes compared by hash), the golden SLiM file
  is unchanged, and a generation still allocates nothing (nor does one of `n` traits).
- **Records**: `trait_mean`, `gmatrix`, `pmatrix`, `trait_lag` (θ − z̄),
  `trait_sel_diff`, one vector or matrix per recorded generation.
- **Theory** (`src/traits_theory.jl`): see the next section.
- **Moments lowering** (`ext/QTLSimCatlabExt/trait_moments.jl`): `trait_moment_images`,
  a second table for the seven steps on a `TraitMoments` state (mean vector, `G`, `P`,
  optimum vector), lowered from the unchanged `lifecycle()` by `lower(f, :moments; images
  = trait_moment_images(rec))`. `run_moments(p)` uses it when `NoTraits > 1`; with
  `NoTraits = 1` the table runs 1×1 matrices and reproduces moments.jl's one-trait
  recursion to 1e-10, with and without drift's erosion, with a move and a width change.
  One deme only: the moments of several demes (moments.jl) are not generalised.
- **Eidos lowering** (`ext/QTLSimCatlabExt/eidos_traits.jl`): `eidos_trait_images`, a
  second table of Eidos for the seven steps. An individual's `a0`/`a1` hold `L·n` values
  (the kernel's layout), `THETA` is a vector on the deme's subpopulation, genotypes and
  phenotypes are `n × N` matrices, fitness uses `matrixMult` with the constant `SINV`,
  mutation draws with `rmvnorm` from `MUTCOV = M/(2L·rate)`, and `record` appends the trait
  means, `G`'s upper triangle and the optimum to QTSim's columns. `slim_script` picks the
  table when `NoTraits > 1`; the one-trait model is byte-identical (golden file) and
  `check_generated.jl` still passes. Migration's image is unchanged (it swaps whole allele
  vectors), and the fecundity splice applies to the trait table unchanged, because the
  spliced segments are verbatim the one-trait ones. Limits: the demes of a SLiM model share
  the trait parameters, `M` must be positive definite (`rmvnorm`), and a step's own `θ` is
  rejected (as in the kernel and moments lowerings, where it would otherwise be a scalar
  for a vector); a step's own vector `θ`, `Ω`, `SelCorr`, `MoveDirection` and mutation
  fields lower to literals.

## The theory

With `P = G + V_E I`, `A = Ω + V_E I` (the matrix `V_s`) and `λ = θ − z̄`, Gaussian selection
on a Gaussian population gives (Lande 1979)

    Δz̄ = G (A + G)⁻¹ λ = G β,        G_survivors = G − G (A + G)⁻¹ G,

and with mutation `G' = G − G (A + G)⁻¹ G + M`. One trait is the scalar recursion of
`theory.jl` (`A + G = V_s + V_G`). This is the aggregate (one haploid locus) recursion,
`run_moments(p; loci = :aggregate)`; the default multilocus model applies the same
selection to the haplotype moments of src/multilocus.jl, as matrices.

- **Equilibrium G** (`gaussian_G`, of the aggregate recursion; `gaussian_G_loci` for
  `L` recombining diploid loci is `gaussian_G(2L·M, A)`): the unique positive semidefinite solution of
  `G (A + G)⁻¹ G = M` is `G = A^{1/2} X A^{1/2}` with `X = ½ (M̃ + (M̃² + 4M̃)^{1/2})`,
  `M̃ = A^{-1/2} M A^{-1/2}`, for any `M` and `A`, not only when they share eigenvectors.
  When they do, `G` has their eigenvectors and the one-trait formula on each axis. For
  weak mutation `G ≈ A # M`, the matrix geometric mean (Lande 1980). `gaussian_G_iterate`
  iterates the recursion as a numerical check, optionally with drift's erosion.
- **Steady lag** (`steady_lag_vector`): an optimum moving by `k = d·u` per generation is
  tracked with `G (A + G)⁻¹ λ = k`, so `λ = (A + G) G⁻¹ k`; for one trait
  `k (V_G + V_s)/V_G`. The mean moves exactly along `k` in the steady state, but the lag
  is deflected away from `g_max`: it is longest where genetic variance is least.
- **Line of least resistance** (`response_deflection`, `gmax`): the response `G β` turns
  from the gradient `β` toward `g_max` (Schluter 1996): in `G`'s eigenbasis each
  component of `β` is scaled by its eigenvalue, so the angle to the `g_max` axis can only
  shrink.

## Validation

`validate/traits_kernel_vs_moments.jl` (log `traits_kernel_vs_moments.log`) runs the kernel
against the moments lowering where the Gaussian recursion is exact, one haploid locus
with a continuum of alleles (2000 individuals, μ = 0.5, μV_s/α² ≈ 125), with `M` and `Ω`
that do not share eigenvectors, 12 replicates per scenario. Criteria fixed before the
first run, all passed:

| | kernel | moments (with drift) | criterion |
|---|---|---|---|
| G₁₁, G₁₂, G₂₂ (×10⁻³) | 4.09, 0.64, 2.80 | 4.18, 0.65, 2.85 | each within 0.10·√(G_ii G_jj): ratios 0.973–0.982 |
| genetic correlation | 0.188 | 0.189 | within 0.05 |
| g_max | | | 0.01° apart (≤ 10°) |
| steady lag, move 0.001 along (1, 2) | (−0.0018, 0.0306) | (−0.0017, 0.0293); `steady_lag_vector` with the kernel's G (−0.0016, 0.0298) | within 10% and 5° |

The kernel's G is about 2% below the moments', and the one-trait kernel in the same regime is
1.2% below the one-trait moments: the trait model misses the Gaussian approximation by
about as much as the scalar one. The lag shows Lande's deflection: the optimum moves at
63° from trait 1, the lag points at 93°, away from `g_max` (22°).

`validate/slim/compare_traits.jl` (log `valid_traits.log`, 40 runs per backend, KS tests,
fail if any p < 0.001, fixed beforehand) compares the kernel with the generated SLiM model
on two traits (`M = m·V_E·[1 0.4; 0.4 0.5]`, fitness correlation −0.3). A, one diploid
locus at equilibrium: G₁₁, G₁₂, G₂₂, genetic correlation and survivors, p = 0.26, 0.91,
0.26, 0.57, 0.76. B, ten loci with the optimum moving 0.001 per generation along (1, 2):
both lag components, the three entries of G and the final mean of trait 2, p = 0.57,
0.99, 0.91, 0.76, 0.26, 0.99. The smallest p is 0.26, so no discrepancy was found.

Many loci: with 50 diploid loci (N = 1000) the aggregate recursion `G' = G − G (A + G)⁻¹ G + M`
does not apply, for one trait as for two: the kernel's V_G is 7.0× `gaussian_vg_exact`
and its two-trait variances 6.5× `gaussian_G`. It treats all of selection's reduction of
G as inherited, which is right for one haploid locus only. The moments lowering's
default is now the multilocus model (src/multilocus.jl): each locus keeps only its share
of the reduction, the rest is linkage disequilibrium that recombination breaks up (by
class of pairs of loci), and drift enters with `erosion = true`. Its equilibrium is
`gaussian_G_loci(p)`, `gaussian_G(2L·M, A)` without drift, and
`gaussian_G_loci(p; drift = true)` with it. In the informational case of
`traits_kernel_vs_moments.log` the kernel is 0.985× the one-trait prediction with drift
and 0.93× the two-trait one, with genetic correlation 0.250 against 0.240. The
pre-specified test is criterion T of `validate/moments_loci.jl`: G₁₁ 0.940, G₁₂ 1.003,
G₂₂ 0.935, correlation 0.258 against 0.241, all within tolerance.

## For Lean

New identities the Julia code relies on, for the Lean layer:

1. *Multivariate mutation input*: if each of `K ~ Poisson(λ)` mutations adds an
   independent `N(0, C)` vector, the expected covariance added is `λ C`; with
   `λ = Ploidy·L·rate` and `C = M/(Ploidy·L·rate)` it is `M` (the matrix version of
   `mutationalVarianceInput_eq`; equivalently the scalar lemma for every `vᵀ M v`).
2. *Gaussian selection on a Gaussian vector*: for `z ~ N(z̄, P)` and
   `W ∝ exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ))`, the selected distribution is
   `N(z̄ + P (Ω + P)⁻¹ λ, P − P (Ω + P)⁻¹ P)`, and the regression of `g` on `z` gives the
   genotype's moments above (matrix algebra: `(P⁻¹ + Ω⁻¹)⁻¹ = P − P (P + Ω)⁻¹ P`).
3. *Equilibrium*: `X (I + X)⁻¹ X = M̃` with `X ⪰ 0` has the unique solution
   `X = ½ (M̃ + (M̃² + 4M̃)^{1/2})`, hence `gaussian_G`; with `n = 1` it is
   `gaussian_fixed_point`.
4. *Steady lag*: `G (A + G)⁻¹ λ = k` iff `λ = (A + G) G⁻¹ k` (`G`, `A + G` invertible);
   with `n = 1`, `steady_lag_response`.
5. *Deflection*: for `G` symmetric positive definite with eigenvalues `g_1 ≥ … ≥ g_n` and
   leading eigenvector `e_1`, `|cos∠(Gβ, e_1)| ≥ |cos∠(β, e_1)|`.
