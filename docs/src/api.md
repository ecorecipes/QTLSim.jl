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
reproduce!
noffspring
mutate!
gamete!
swapbuffers!
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

## Theory

The predictions the simulations are compared with; `run_moments` iterates them for any
life cycle.

```@docs
selection_variance
gaussian_vg
gaussian_vg_exact
effective_size
house_of_cards_vg
steady_lag
critical_rate
mutation_adjusted
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
QTLSimCatlabExt.EidosCode
QTLSimCatlabExt.EidosMeta
QTLSimCatlabExt.DemeOp
QTLSimCatlabExt.Operator
QTLSimCatlabExt.Demes
QTLSimCatlabExt.Stratified
QTLSimCatlabExt.Population
QTLSimCatlabExt.Moments
```
