"""
    QTLSim

A Julia port of QTSim, the quantitative-genetics simulation of Pagel, Gardner and Meade
(Scientific Reports, 2025; the C source calls itself QTLSim 1.0): an additive polygenic
trait under Gaussian stabilising viability selection toward an optimum that can move,
with random mating, Poisson recombination and Poisson–Gaussian mutation. The kernel
reproduces QTSim's numerics; `validate/` compares it with QTSim and with the paper's
published results.

The categorical layer, in which the population is an ACSet and the life cycle is a
composite of operators (`⋅`) that can be stratified into demes (`⊗`), is the package
extension `QTLSimCatlabExt`. Julia loads it when Catlab is loaded (`import Catlab`). Its
functions are declared, documented and exported here.
"""
module QTLSim

using Random, Statistics, Distributions

include("params.jl")
include("qtsim_io.jl")
include("model.jl")
include("theory.jl")
include("categorical.jl")

# the simulation kernel
export Params, World, simulate, run!, generation!, init_env!, checkparams,
       develop!, select!, moveoptimum!, reproduce!, mutate!, gamete!, swapbuffers!,
       geneticVar, envVar, analyticalFitSD, initFitSD, muteVar,
       read_qtsim_options, write_qtsim_options
# analytical predictions
export selection_variance, gaussian_vg, gaussian_vg_exact, effective_size, house_of_cards_vg,
       steady_lag, critical_rate, mutation_adjusted
# the categorical layer (implemented by QTLSimCatlabExt; load Catlab to use it)
export lifecycle, lower, slim_script, operators, operator, demes, run_composed, migrate!,
       migration, totalpop, ndemes, population_schema, lift, lower!, genotypes, phenotypes,
       genome_from_acset, run_moments

end # module
