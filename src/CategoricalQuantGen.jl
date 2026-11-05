"""
    CategoricalQuantGen

Compositional reformulation of the Pagel et al. (Scientific Reports, 2025)
quantitative-genetics model — the executable counterpart of the Lean theory in
`lean/` and a compositional re-expression of the reference C simulation QTSim
(`Repository/Code/QTSim/`).

Two levels, in the AlgebraicJulia spirit:

* an **abstract** layer — the population as a Catlab/ACSets attributed C-set
  (`schema.jl`), and evolutionary processes as first-class `Operator`s composed
  with the unicode operators `⋅` (diagrammatic composition) and `⊗`
  (stratification) from Catlab (`compose.jl`);
* a **lowered** layer — the fast array kernel (`params.jl`, `model.jl`) whose
  numerics replicate QTSim so output is statistically consistent with the C
  reference (see `validate/`). `lift`/`lower!` translate between the two.

One generation is the composite `develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap`
and a run is its iterate (cf. Lean `Composition.lean`).
"""
module CategoricalQuantGen

using Random, Statistics, Distributions
using Catlab, Catlab.CategoricalAlgebra, ACSets

include("params.jl")
include("model.jl")
include("schema.jl")
include("compose.jl")

# lowered array kernel
export Params, World, simulate, run!, generation!,
       develop!, select!, moveoptimum!, reproduce!, mutate!, gamete!, swapbuffers!,
       geneticVar, envVar, analyticalFitSD, initFitSD, muteVar
# abstract ACSet data model (Catlab / ACSets)
export SchPopulation, Population, lift, lower!, genotypes, phenotypes, genome_from_acset
# categorical operator sugar (⋅ composition, ⊗ stratification)
export Operator, ⋅, ⊗, develop, select, reproduce, mutate, swap, generation,
       run_composed, Demes, Stratified, migrate!, migration, totalpop, ndemes

end # module
