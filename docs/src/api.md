# API reference

```@meta
CurrentModule = QTLSim
```

```@docs
QTLSim
```

## Simulation

```@docs
Params
checkparams
World
genomes
setgenomes!
simulate
run!
generation!
init_env!
records
record!
rescale!
develop!
select!
moveoptimum!
fluctuating_optimum
reproduce!
noffspring
mutate!
gamete!
swapbuffers!
```

## Populations initialised from theory

A run from zero genomes spends its first thousands of generations building its standing
variation. `simulate(p; init = :equilibrium, burnin = B)` (and `replicate` with the same
keywords) instead starts from a population drawn from the stationary allele-frequency
spectrum of mutation, stabilising selection and drift (`stationary_spectrum`, src/init.jl),
`B` generations before the first change of selection, and keeps the protocol's generation
numbers. The draw is in linkage equilibrium; the burn-in restores the rest. How long it
must be, and what it saves, is measured in `validate/init_equilibrium.jl`. One diploid
trait, fixed size, viability selection, kernel only.

```@docs
equilibrium_world
equilibrium_genomes
stationary_spectrum
QTLSim.stationary_sites
QTLSim.stationary_density
QTLSim.stationary_mean_var
QTLSim.sojourn_mass
QTLSim.sojourn_ratio
QTLSim.protocol_start
QTLSim.first_event
```

## Sweeps and replicates

```@docs
sweep
replicate
```

## QTSim options files

```@docs
read_qtsim_options
write_qtsim_options
```

## Derived quantities

```@docs
geneticVar
envVar
analyticalFitSD
initFitSD
muteVar
```

Several traits (`NoTraits > 1`; see [Several traits](traits_design.md)):

```@docs
mutational_covariance
selection_shape
selection_covariance
move_direction
fluct_direction
move_noise_covariance
traitvalues
```

Models: the families of lowerings of the one life cycle, each with its own step
parameters and tables of images:

```@docs
QTLSim.Model
OneTraitModel
TraitModel
ReactionNormModel
model_of
step_schema
QTLSim.step_fields
QTLSim.step_own_params
QTLSim.checkstep
```

## Theory

The predictions the simulations are compared with; `run_moments` iterates them for any
life cycle. The multilocus Gaussian model (`loci = :multilocus`: genic variance, linkage
disequilibrium broken up by recombination, optional drift) has the equilibrium
`gaussian_vg_loci`; `gaussian_vg_exact` is the equilibrium of the aggregate (one haploid
locus) recursion, `run_moments(p; loci = :aggregate)`. The rare-alleles closure
(`loci = :rare_alleles`, src/rare_alleles.jl), within 15% of the kernel in 46 of 47
weak-selection cells across the Gaussian, house-of-cards and drift regimes
(`validate/moments_regimes.jl`), has the equilibrium `rare_alleles_vg`;
`house_of_gauss_vg` is its closed form for loci in linkage equilibrium. For several traits
(src/rare_alleles_traits.jl) the equilibrium is `rare_alleles_G`, the linkage-equilibrium
form `house_of_gauss_G` and the multivariate house of cards `house_of_cards_G`. For
metapopulations whose demes lag behind different optima, `loci = :rare_local`
(src/local_alleles.jl, `validate/local_adaptation.jl`) weakens the closure by the gain
`QTLSim.local_gain` of locally favoured alleles; it equals `:rare_alleles` for one deme, a
common optimum and no migration. The default,
`loci = :auto`, is the rare-alleles closure for one trait (diploid, or haploid clones) and
`:rare_local` for several traits and for one trait in several demes, read from the declared
closures (`QTLSim.moments_loci`, `QTLSim.auto_closure`); `multilocus_equilibrium`
takes the same default.
Under a moving optimum, the allele-frequency spectrum (`loci = :spectrum`,
src/spectrum.jl, `validate/vg_move.jl`) follows how sweeps, drift and linked selection
change the variance; `spectrum_equilibrium` and `spectrum_vg` are its stationary state.

```@docs
selection_variance
gaussian_vg_loci
multilocus_equilibrium
gaussian_vg
gaussian_vg_exact
effective_size
house_of_cards_vg
house_of_gauss_vg
rare_alleles_vg
QTLSim.mutation_effect_variance
house_of_cards_G
house_of_gauss_G
rare_alleles_G
QTLSim.mutation_effect_covariance
QTLSim.rare_alleles_rate
QTLSim.rare_hc_matrix
QTLSim.hap_clonal
QTLSim.local_gain
QTLSim.local_weights
QTLSim.local_sinv
QTLSim.migration_classes
spectrum_equilibrium
spectrum_vg
QTLSim.AlleleSpectrum
QTLSim.spectrum_select
QTLSim.spectrum_step!
QTLSim.spectrum_linked
QTLSim.spectrum_mutate!
QTLSim.moments_loci
steady_lag
critical_rate
mutation_adjusted
extinction_distribution
rescue_probability
```

Fluctuating optima (`MoveType = :periodic` or `:red_noise`):

```@docs
periodic_response
red_noise_lag_var
expected_lag2
lag_load
```

Phenotypic plasticity, a linear reaction norm (Chevin, Lande and Mace 2010; Lande 2009; see
[Phenotypic plasticity](plasticity_design.md)). `steady_lag`, `critical_rate`,
`red_noise_lag_var`, `expected_lag2` and `lag_load` above take the slope, the optimum's
slope and the cue's delay:

```@docs
plasticity_ratio
plasticity_cost
QTLSim.isplastic
sustainable_rate
drift_critical_rate
half_survival_rate
red_noise_lag_autocov
reaction_norm_response
reaction_norm_lag
QTLSim.plastic_cue!
QTLSim.express_plasticity!
QTLSim.select_norm!
```

Several traits with pleiotropy (Lande 1979, 1980; Schluter 1996):

```@docs
selection_variance_matrix
gaussian_G_loci
gaussian_G
gaussian_G_iterate
selection_gradient
trait_response
steady_lag_vector
gmax
vector_angle
response_deflection
red_noise_lag_cov
periodic_lag_cov
expected_lag_cov
trait_lag_load
```

Neutral differentiation among demes (F_ST and Q_ST, `validate/fst.jl`):

```@docs
qst_island
island_variances
backward_migration
deme_coalescence
exchange_origins
neutral_deme_variances
fst_stepping_stone
pair_recombination
QTLSim.pair_classes
QTLSim.migration_mixing
```

## Ancestry

Identity by descent read off the block pool, which costs nothing, and recorded pedigrees,
which are off by default (`src/ancestry.jl`, `validate/coalescence.jl`). Two chromosomes
that hold the same block are identical by descent over it, without an intervening mutation
or crossover inside the block; chromosomes that hold different blocks may still be
identical by descent ([`shared_blocks`](@ref) says exactly what block identity means).

```@docs
tag_founders!
shared_blocks
ibd_fraction
ibd_sharing
ibd_fst
```

A pedigree records, for every offspring chromosome, the parental chromosomes it copies and
the crossovers between them, with the mutations and the migrations, as tree-sequence
tables. From it come the complete ancestry of every locus and the realised coalescence
times that [`deme_coalescence`](@ref) predicts.

```@docs
record_pedigree!
Pedigree
QTLSim.Ancestry
stop_recording!
simplify!
samples
pedigree_time
pedigree_bytes
marginal_tree
mrca
tmrca
allele
alleles
pedigree_ibd
pair_coalescence
PairCoalescence
coalescence_cdf
coalescence_mean
coalescence_transform
```

Tree-sequence tables in tskit's data model and text format. tskit is not a dependency and
the files have not been loaded into it; they are checked in Julia.

```@docs
tskit_tables
write_tskit
read_tskit
check_tables
table_tmrca
```

## The categorical layer

These functions are implemented by the package extension `QTLSimCatlabExt`, which Julia
loads when Catlab is loaded:

```julia
using QTLSim
import Catlab
using Catlab: ⋅, ⊗
L = lifecycle()                                                        # the life cycle as syntax
(; develop, select, reproduce, mutate, swap, generation) = operators()  # its image in the kernel
```

### The life cycle and its lowerings

```@docs
lifecycle
withparams
lower
slim_script
run_composed
run_moments
```

### Operators

```@docs
operators
operator
demes
migrate!
migration
island_migration
stepping_stone
QTLSim.check_migration
totalpop
ndemes
```

### The population as an ACSet

```@docs
population_schema
lift
lower!
genotypes
phenotypes
genome_from_acset
```

### Naturality: validation as data

A statistic σ maps one backend's states to another's; agreement between backends is a
naturality square σ ∘ K(f) = M(f) ∘ σ, exact for some life cycles and within a closure's
error for others. A square is a [`NaturalityCase`](@ref), a row of the table of cases in
`validate/naturality_cases.jl`; see [Naturality and model morphisms](naturality_design.md).

```@docs
NaturalityCase
NaturalityResult
check_naturality
naturality_steps
QTLSim.first_failure
QTLSim.Statistic
QTLSim.window_mean
QTLSim.metapop_mean
QTLSim.gmatrix_mean
QTLSim.qst_stat
QTLSim.ks_pvalue
```

### Closures and model morphisms, declared as data

The closures of the moments lowering, their domains and the order of refinements, and the
maps between models, each with the reductions it must satisfy (the tests are generated from
them). `loci = :auto` is the finest closure whose automatic domain contains the model and
the life cycle. The [support matrix](support.md) is generated from the same declarations
and the code's own guards.

```@docs
QTLSim.ClosureContext
QTLSim.closure_context
QTLSim.ClosureDecl
QTLSim.closures
QTLSim.auto_closure
QTLSim.closure_reasons
QTLSim.ModelMorphism
QTLSim.Reduction
QTLSim.check_reduction
model_morphisms
QTLSim.closure_order
explain
support_matrix
support_markdown
```

### Types of the extension

```@meta
CurrentModule = Main
```

```@docs
QTLSimCatlabExt
QTLSimCatlabExt.Stage
QTLSimCatlabExt.Step
QTLSimCatlabExt.STEP_PARAMS
QTLSimCatlabExt.kernel_images
QTLSimCatlabExt.moment_images
QTLSimCatlabExt.with_step_params
QTLSimCatlabExt.EidosCode
QTLSimCatlabExt.EidosMeta
QTLSimCatlabExt.DemeOp
QTLSimCatlabExt.Operator
QTLSimCatlabExt.Demes
QTLSimCatlabExt.Stratified
QTLSimCatlabExt.Population
QTLSimCatlabExt.Moments
QTLSimCatlabExt.DemeMoments
QTLSimCatlabExt.TraitMoments
QTLSimCatlabExt.TraitCov
QTLSimCatlabExt.TraitDemeMoments
QTLSimCatlabExt.trait_moment_images
QTLSimCatlabExt.run_trait_moments
```
