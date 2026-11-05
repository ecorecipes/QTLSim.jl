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
World
simulate
run!
generation!
develop!
select!
moveoptimum!
reproduce!
mutate!
gamete!
swapbuffers!
```

## Derived quantities

```@docs
geneticVar
envVar
analyticalFitSD
initFitSD
muteVar
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
```

### Operators

```@docs
operators
operator
demes
migrate!
migration
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
QTLSimCatlabExt.Operator
QTLSimCatlabExt.Demes
QTLSimCatlabExt.Stratified
QTLSimCatlabExt.Population
```
