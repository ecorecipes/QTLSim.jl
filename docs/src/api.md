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
traitvalues
```

Models: the families of lowerings of the one life cycle, each with its own step
parameters and tables of images:

```@docs
QTLSim.Model
OneTraitModel
TraitModel
model_of
step_schema
QTLSim.step_fields
QTLSim.step_own_params
QTLSim.checkstep
```

## Theory

The predictions the simulations are compared with; `run_moments` iterates them for any
life cycle. Its default, the multilocus Gaussian model (genic variance, linkage
disequilibrium broken up by recombination, optional drift), has the equilibrium
`gaussian_vg_loci`; `gaussian_vg_exact` is the equilibrium of the aggregate (one haploid
locus) recursion, `run_moments(p; loci = :aggregate)`.

```@docs
selection_variance
gaussian_vg_loci
multilocus_equilibrium
gaussian_vg
gaussian_vg_exact
effective_size
house_of_cards_vg
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
QTLSim.migration_mixing
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
```
