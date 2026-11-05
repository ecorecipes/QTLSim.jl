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
(; develop, select, reproduce, mutate, swap, generation) = operators()
```

```@docs
operators
operator
demes
run_composed
migrate!
migration
totalpop
ndemes
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
QTLSimCatlabExt.Operator
QTLSimCatlabExt.Demes
QTLSimCatlabExt.Stratified
QTLSimCatlabExt.Population
```
