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
using LinearAlgebra: I, kron, eigvals, eigen, Symmetric, cholesky, det, tr, dot, norm, schur, ordschur!, isposdef, Diagonal

include("params.jl")
include("models.jl")
include("qtsim_io.jl")
include("blocks.jl")
include("model.jl")
include("traits.jl")
include("theory.jl")
include("traits_theory.jl")
include("multilocus.jl")
include("rare_alleles.jl")
include("rare_alleles_traits.jl")
include("local_alleles.jl")
include("spectrum.jl")
include("init.jl")
include("migration.jl")
include("categorical.jl")
include("sweep.jl")

# the simulation kernel
export Params, World, genomes, setgenomes!, simulate, run!, generation!, init_env!, checkparams,
       develop!, select!, moveoptimum!, fluctuating_optimum, reproduce!, noffspring, mutate!, gamete!, swapbuffers!,
       geneticVar, envVar, analyticalFitSD, initFitSD, muteVar,
       mutational_covariance, selection_covariance, selection_shape, move_direction, traitvalues,
       read_qtsim_options, write_qtsim_options
# models: the families of lowerings of the one life cycle, and their step parameters
export OneTraitModel, TraitModel, model_of, step_schema
# analytical predictions
export selection_variance, gaussian_vg, gaussian_vg_exact, gaussian_vg_loci, multilocus_equilibrium, effective_size, house_of_cards_vg,
       house_of_gauss_vg, rare_alleles_vg, house_of_cards_G, house_of_gauss_G, rare_alleles_G,
       spectrum_equilibrium, spectrum_vg,
       stationary_spectrum, equilibrium_genomes, equilibrium_world,
       steady_lag, critical_rate, mutation_adjusted, extinction_distribution, rescue_probability
# fluctuating optima
export periodic_response, red_noise_lag_var, expected_lag2, lag_load
# several traits (Lande 1979, 1980)
export selection_variance_matrix, gaussian_G, gaussian_G_iterate, selection_gradient, trait_response,
       steady_lag_vector, gmax, vector_angle, response_deflection, gaussian_G_loci,
       red_noise_lag_cov, periodic_lag_cov, expected_lag_cov, trait_lag_load, fluct_direction, move_noise_covariance
# migration graphs among demes
export island_migration, stepping_stone
# neutral differentiation among demes (F_ST, Q_ST)
export backward_migration, deme_coalescence, neutral_deme_variances, exchange_origins, qst_island, island_variances, fst_stepping_stone,
       pair_recombination
# the categorical layer (implemented by QTLSimCatlabExt; load Catlab to use it)
export lifecycle, withparams, lower, slim_script, operators, operator, demes, run_composed, migrate!,
       migration, totalpop, ndemes, population_schema, lift, lower!, genotypes, phenotypes,
       genome_from_acset, run_moments
# sweeps and replicates on all threads, with a resumable CSV cache
export sweep, replicate

end # module
