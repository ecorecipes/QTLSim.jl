# Simulation parameters and derived quantities.
#
# Field names and defaults mirror QTSim's `OPTIONS` (Repository/Code/QTSim/src/QTSim.h)
# and `SetDefOptions`. The derived quantities mirror `CaclOptParm` /
# `GetAnalyticalFitSD` and are the Julia counterparts of the Lean definitions in
# lean/CategoricalQuantGen/Heritability.lean.

"""
    Params

Configuration for one simulation, with QTSim's defaults. `Heritability` is the
narrow-sense heritability h² (a single scalar); `MutationScalar` is the mutation scalar
m. The documentation's Parameters page lists every field with its meaning and its
QTSim option. Values are checked when a `World` is made (see [`checkparams`](@ref)).

Two fields have no QTSim counterpart, or behave differently:

- `CrossoverModel`: QTSim draws `Poisson(RecombPoissonMean)` crossover points per
  gamete but only the lowest one switches strands (its point index never advances), so
  a gamete has at most one crossover. `:qtsim` (default) reproduces that, as the
  published simulations did; `:all` applies every drawn point (the evidently intended
  Poisson-crossover model). With `NoQTL == 1` the two are identical.
- `Drift`: QTSim has the option but switches it off at the start of every run, so it
  never has an effect there. Here `Drift = true` makes every fitness 1: a neutral run,
  for comparison.
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

"Initial fitness width: `InitFitnessSD`, or the analytical value when it is −1 (QTSim `CaclOptParm`)."
initFitSD(p::Params) = p.InitFitnessSD == -1 ? analyticalFitSD(p) : p.InitFitnessSD

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

"""
    checkparams(p)

Throw an `ArgumentError` if `p` is not a model QTSim would run, or would run
differently: `Ploidy` 1 or 2, at least one individual (two if diploid), locus, sampled
generation and mutation per gamete; `0 ≤ Heritability ≤ 1` and `PhenotypeVar > 0`;
`MoveType` one of `:phenotype`, `:fitness`, `:continuous` (QTSim falls back to the
phenotype move on an unknown name; a misspelt symbol here would silently change the
model); `CrossoverModel` `:qtsim` or `:all`; the schedules `MoveFitnessSDGen` and
`ChangeFitnessGen` either −1 (off) or a generation, with `ChangeFitnessSD > 0` when the
width changes; and `InitFitnessSD` either −1 (the analytical width, which needs
`V_E > 0` and `m > 0`) or positive. Called by `World(p)`.
"""
function checkparams(p::Params)
    bad(msg) = throw(ArgumentError("Params: " * msg))
    p.Ploidy in (1, 2) || bad("Ploidy must be 1 or 2, got $(p.Ploidy)")
    p.PopSize >= p.Ploidy || bad("PopSize must be at least $(p.Ploidy) for Ploidy $(p.Ploidy), got $(p.PopSize)")
    p.NoQTL >= 1 || bad("NoQTL must be at least 1, got $(p.NoQTL)")
    p.NoGenerations >= 0 || bad("NoGenerations must be non-negative, got $(p.NoGenerations)")
    p.Sample >= 1 || bad("Sample must be at least 1, got $(p.Sample)")
    p.PhenotypeVar > 0 || bad("PhenotypeVar must be positive, got $(p.PhenotypeVar)")
    0 <= p.Heritability <= 1 || bad("Heritability must be in [0, 1], got $(p.Heritability)")
    p.MutationScalar >= 0 || bad("MutationScalar must be non-negative, got $(p.MutationScalar)")
    p.MutationRatePerGamete > 0 || bad("MutationRatePerGamete must be positive, got $(p.MutationRatePerGamete)")
    p.RecombPoissonMean >= 0 || bad("RecombPoissonMean must be non-negative, got $(p.RecombPoissonMean)")
    p.CrossoverModel in (:qtsim, :all) || bad("CrossoverModel must be :qtsim or :all, got :$(p.CrossoverModel)")
    p.MoveType in (:phenotype, :fitness, :continuous) ||
        bad("MoveType must be :phenotype, :fitness or :continuous, got :$(p.MoveType)")
    (p.MoveFitnessSDGen == -1 || p.MoveFitnessSDGen >= 0) ||
        bad("MoveFitnessSDGen must be -1 (no move) or a generation, got $(p.MoveFitnessSDGen)")
    (p.ChangeFitnessGen == -1 || p.ChangeFitnessGen >= 0) ||
        bad("ChangeFitnessGen must be -1 (no change) or a generation, got $(p.ChangeFitnessGen)")
    (p.ChangeFitnessGen == -1 || p.ChangeFitnessSD > 0) ||
        bad("ChangeFitnessSD must be positive when ChangeFitnessGen is set, got $(p.ChangeFitnessSD)")
    if p.InitFitnessSD == -1
        envVar(p) > 0 && p.MutationScalar > 0 ||
            bad("the analytical fitness width (InitFitnessSD = -1) needs Heritability < 1 and MutationScalar > 0; set InitFitnessSD")
    else
        p.InitFitnessSD > 0 || bad("InitFitnessSD must be -1 (analytical) or positive, got $(p.InitFitnessSD)")
    end
    p
end
