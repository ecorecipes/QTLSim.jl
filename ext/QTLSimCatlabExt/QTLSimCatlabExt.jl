"""
    QTLSimCatlabExt

The categorical layer of QTLSim, loaded with Catlab. The population is an attributed
C-set (ACSet) on the schema of `schema.jl`; each step of the life cycle is an
`Operator`, composed with Catlab's diagrammatic composition `⋅` and stratified into
demes with its monoidal product `⊗` (`compose.jl`). The facts that justify the
decomposition, that composition is associative and that a run is the iterate of one
composed generation, are proved in lean/CategoricalQuantGen/Composition.lean.
"""
module QTLSimCatlabExt

using QTLSim
using QTLSim: World, Params, envsd
import QTLSim: operators, operator, demes, run_composed, migrate!, migration, totalpop,
               ndemes, population_schema, lift, lower!, genotypes, phenotypes,
               genome_from_acset
using Catlab: @present, @acset_type, FreeSchema, add_parts!, add_part!, set_subpart!,
              nparts, incident
import Catlab: ⋅, ⊗
using Statistics: mean, var

include("schema.jl")
include("compose.jl")

end # module
