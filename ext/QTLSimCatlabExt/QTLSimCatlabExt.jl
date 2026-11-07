"""
    QTLSimCatlabExt

The categorical layer of QTLSim, loaded with Catlab. The life cycle is a morphism of a
free category whose objects are the stages of a generation (`lifecycle.jl`), lowered by
functors to the simulation kernel, where each step is an `Operator` composed with
Catlab's `⋅` and stratified into demes with its monoidal product `⊗` (`compose.jl`), and
to a SLiM model in Eidos (`eidos.jl`), and to the deterministic recursions of
quantitative-genetic theory on a Gaussian population's moments (`moments.jl`). The population is an attributed C-set (ACSet) on
the schema of `schema.jl`. The facts that justify the decomposition, that composition is
associative, that a run is the iterate of one composed generation and that a lowering is
fixed by its action on the steps, are proved in lean/CategoricalQuantGen/Composition.lean
and Lowering.lean.
"""
module QTLSimCatlabExt

using QTLSim
using QTLSim: World, Params, envsd, records, record!, rescale!, init_env!
import QTLSim: operators, operator, demes, run_composed, migrate!, migration, totalpop,
               ndemes, population_schema, lift, lower!, genotypes, phenotypes,
               genome_from_acset, lifecycle, lower, slim_script, run_moments
using Catlab: @present, @acset_type, FreeSchema, FreeSymmetricMonoidalCategory, add_parts!,
              add_part!, set_subpart!, nparts, incident, functor, head, args, dom, codom, Hom,
              GATExpr, ObExpr
import Catlab: ⋅, ⊗, compose, id, otimes, braid, munit
using Statistics: mean, var
using Random: randperm

include("schema.jl")
include("compose.jl")
include("eidos.jl")
include("lifecycle.jl")
include("moments.jl")

end # module
