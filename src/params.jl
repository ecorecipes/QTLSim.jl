# Simulation parameters and derived quantities.
#
# Field names and defaults mirror QTSim's `OPTIONS` (Repository/Code/QTSim/src/QTSim.h)
# and `SetDefOptions`. The derived quantities mirror `CaclOptParm` /
# `GetAnalyticalFitSD` and are the Julia counterparts of the Lean definitions in
# lean/CategoricalQuantGen/Heritability.lean.

"""
    Params

Configuration for one simulation. `Heritability` is the narrow-sense heritability
h² (a single scalar); `MutationScalar` is the mutation scalar m. See QTSim's
options file format for the meaning of each field.

`CrossoverModel` is the one field with no QTSim counterpart. QTSim draws
`Poisson(RecombPoissonMean)` crossover points per gamete but only the lowest one
switches strands (its point index never advances), so a gamete has at most one
crossover. `:qtsim` (default) reproduces that, as the published simulations did;
`:all` applies every drawn point (the evidently intended Poisson-crossover model).
With `NoQTL == 1` the two are identical.
"""
Base.@kwdef struct Params
    NoGenerations::Int          = 20000
    Sample::Int                 = 1
    PopSize::Int                = 1000
    NoQTL::Int                  = 1
    Ploidy::Int                 = 1        # 1 = haploid, 2 = diploid
    PhenotypeVar::Float64       = 1.0      # V_P
    Heritability::Float64       = 0.5      # h²
    MutationScalar::Float64     = 1/1000   # m
    MutationRatePerGamete::Float64 = 0.0002
    RecombPoissonMean::Float64  = 1.0
    CrossoverModel::Symbol      = :qtsim   # :qtsim (lowest crossover point only, as QTSim) | :all  — Julia-only
    InitFitnessSD::Float64      = -1.0     # <0 ⇒ use analytical width
    ChangeFitnessGen::Int       = -1       # generation to change fitness width
    ChangeFitnessSD::Float64    = 0.0      # new width at ChangeFitnessGen
    MoveFitnessSDGen::Int       = -1       # generation the optimum starts moving
    MoveFitnessSD::Float64      = 1.0      # move size (in the MoveType's units)
    MoveType::Symbol            = :phenotype  # :phenotype | :fitness | :continuous
    Drift::Bool                 = false
    StandardisePhenotype::Bool  = false
    Seed::Union{Nothing,Int}    = nothing
end

"Additive genetic variance target V_G = h²·V_P (QTSim `OptVarTarget`)."
geneticVar(p::Params) = p.Heritability * p.PhenotypeVar

"Environmental variance V_E = V_P − V_G (QTSim `EnvVar`)."
envVar(p::Params) = p.PhenotypeVar - geneticVar(p)

"""
    analyticalFitSD(p)

Analytical initial fitness-function width ω = V_G / √(m·V_E)
(QTSim `GetAnalyticalFitSD`; see Lean `analyticalFitSD_eq`).
"""
analyticalFitSD(p::Params) = geneticVar(p) / sqrt(p.MutationScalar * envVar(p))

"Initial fitness width: user `InitFitnessSD` if ≥0, else the analytical value."
initFitSD(p::Params) = p.InitFitnessSD < 0 ? analyticalFitSD(p) : p.InitFitnessSD

"""
    muteVar(p)

Expected per-generation mutational variance input, m·V_E (QTSim `MuteVar`).
The Poisson mutation operator is calibrated to this value; see Lean
`mutationalVarianceInput_eq`.
"""
muteVar(p::Params) = p.MutationScalar * envVar(p)

"""
    Params(p::Params; kw...)

A copy of `p` with the fields named in `kw` replaced, for example
`Params(p; NoGenerations = 5000, MoveFitnessSD = 0.01)`. An unknown field name is an
`ArgumentError`.
"""
function Params(p::Params; kw...)
    fields = Dict{Symbol,Any}(f => getfield(p, f) for f in fieldnames(Params))
    for (k, v) in kw
        haskey(fields, k) || throw(ArgumentError("Params has no field $k"))
        fields[k] = v
    end
    Params(; fields...)
end
