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

Some fields have no QTSim counterpart, or behave differently:

- `CrossoverModel`: QTSim draws `Poisson(RecombPoissonMean)` crossover points per
  gamete but only the lowest one switches strands (its point index never advances), so
  a gamete has at most one crossover. `:qtsim` (default) reproduces that, as the
  published simulations did; `:all` applies every drawn point (the evidently intended
  Poisson-crossover model). With `NoQTL == 1` the two are identical.
- `Drift`: QTSim has the option but switches it off at the start of every run, so it
  never has an effect there. Here `Drift = true` makes every fitness 1: a neutral run,
  for comparison.
- `Demography` and `GrowthRate`: QTSim always refills the population to `PopSize`, so it
  dies only when fewer than two survive one generation (`:fixed`, the default). With
  `:beverton_holt` the population is density-regulated, with carrying capacity
  `PopSize` and `GrowthRate` offspring per survivor at low density (see
  [`noffspring`](@ref)), so that it shrinks when mean fitness is low and can be rescued
  by adaptation.
- `Selection`: QTSim's selection is on viability (`:viability`, the default): an
  individual survives iff `U < W(z)`, and survivors mate at random. With `:fecundity`
  everyone survives and each parent is drawn with probability `W/ΣW` (two distinct
  parents when diploid), and under `:beverton_holt` the offspring number uses the sum of
  fitnesses `ΣW = N·W̄` in place of the number of survivors (see [`noffspring`](@ref)).
  `ΣW` is the expected number of survivors under viability, so the two have the same
  expected selection and, at a fixed size, differ only in its variance. At a fixed size
  fecundity selection is soft: the population dies only if fewer than `Ploidy`
  individuals have positive fitness.
- `MoveType = :red_noise` or `:periodic`, with `MoveNoiseSD`, `MoveAutocorr`,
  `MoveAmplitude` and `MovePeriod`: QTSim's optimum moves by a rule; these two make it
  fluctuate, around a linear trend of `MoveFitnessSD` per generation, with stationary
  AR(1) noise or a sine wave (see [`moveoptimum!`](@ref)).
- `NoTraits`, `MutCorr`, `SelCorr`, `MutCov`, `SelCov` and `MoveDirection`: QTSim has one
  trait. With `NoTraits = n > 1` each locus has a vector of `n` effects (pleiotropy), a
  mutation adds a draw from `N(0, M/(Ploidy·NoQTL·MutationRatePerGamete))`, so that `M`
  ([`mutational_covariance`](@ref)) is injected per individual per generation, each trait
  has environmental variance `V_E`, independently, and fitness is
  `W(z) = c·exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ))` with QTSim's `c = 1/(0.4√(2π))`, so that one trait
  is QTSim's fitness. `Ω = ω² S`, where `ω` is the width (`InitFitnessSD` or the
  analytical width, changed by `ChangeFitnessSD`) and the shape `S` has unit diagonal
  and correlation `SelCorr`, or is `SelCov/ω²` initially. The optimum moves along the
  unit vector `u` ∝ `MoveDirection` ([`move_direction`](@ref)): by `d·u` per generation
  (`:continuous`), to `z̄ + √V_P·d·u` (`:phenotype`) or to `z̄ + ω d u` (`:fitness`). The
  fluctuating moves have the trend `d·u` per generation and fluctuate along the unit vector
  `f` ∝ `FluctDirection` (`u` by default, [`fluct_direction`](@ref)): `A sin(2π t/P)·f`
  (`:periodic`), or AR(1) noise with stationary covariance `MoveNoiseCov`, by default
  `MoveNoiseSD²·f fᵀ` ([`move_noise_covariance`](@ref)). `FluctDirection` and
  `MoveNoiseCov` are Julia-only, for several traits. See the design note
  `traits_design.md`.
- `Plasticity`, `OptimumSlope`, `CueLag` and `PlasticityCostSD` (Julia-only; off with
  `Plasticity = 0`, the default): phenotypic plasticity, a linear reaction norm
  `z = g + b·ε + e` with slope `b = Plasticity`, after Chevin, Lande and Mace (2010). The
  environment `ε` is measured by the optimum it sets, `θ = B·ε` with `B = OptimumSlope`, so
  the plastic part of the phenotype is `(b/B)·θ`: a fraction `b/B` of the optimum's
  displacement from the origin of the trait's scale is tracked without any genetic
  change. The cue is the environment `CueLag = τ` generations before selection
  (`τ = 0`: a perfect cue; the environment before the run is taken to be the first one
  seen), so `z_t = g + (b/B)·θ_{t−τ} + e`. A cost of plasticity multiplies every fitness
  by `exp(−b²/(2ω_b²))`, `ω_b = PlasticityCostSD` (`Inf`: no cost): stabilising selection
  on the slope around 0. The reaction norm is expressed in `select`, once the optimum of
  the generation is set ([`select!`](@ref)). It needs an optimum with its own trajectory:
  `:continuous`, `:red_noise`, `:periodic`, or none; QTSim's `:phenotype` and `:fitness`
  moves re-anchor the optimum to the mean phenotype every generation, which a plastic
  shift would only carry along, and are rejected.
- `ReactionNorm = true` with `NoTraits = 2` (Julia-only): a heritable reaction norm. Trait
  1 is the elevation `a` and trait 2 the slope: an individual's slope is `Plasticity + g₂`,
  its phenotype `z = g₁ + (Plasticity + g₂)·ε + e` with `ε = θ_{t−τ}/B`, selection acts on
  `z` (one-dimensional, width `ω`) and the cost on the individual's slope. `MutCov` gives
  the mutational covariance of elevation and slope. See `plasticity_design.md`.
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
    MoveType::Symbol            = :phenotype  # :phenotype | :fitness | :continuous | :red_noise | :periodic (last two Julia-only)
    MoveNoiseSD::Float64        = 0.0      # σ, stationary SD of the optimum's AR(1) noise (:red_noise)
    MoveAutocorr::Float64       = 0.0      # ρ, the noise's autocorrelation between generations (:red_noise)
    MoveAmplitude::Float64      = 0.0      # A, amplitude of the optimum's cycle (:periodic)
    MovePeriod::Float64         = 100.0    # P, period of the cycle in generations (:periodic)
    Drift::Bool                 = false
    StandardisePhenotype::Bool  = false
    Demography::Symbol          = :fixed   # :fixed (QTSim: refill to PopSize) | :beverton_holt — Julia-only
    GrowthRate::Float64         = 2.0      # R, offspring per survivor at low density (:beverton_holt)
    Selection::Symbol           = :viability  # :viability (QTSim: survive iff U < W) | :fecundity (mate ∝ W) — Julia-only
    # phenotypic plasticity, a linear reaction norm z = g + b·ε + e — Julia-only; off with Plasticity = 0
    Plasticity::Float64         = 0.0      # b, the slope of the reaction norm on the environment ε
    OptimumSlope::Float64       = 1.0      # B, the optimum's slope on the environment: θ = B·ε
    CueLag::Int                 = 0        # τ, generations from the cue that sets the plastic response to selection
    PlasticityCostSD::Float64   = Inf      # ω_b, width of stabilising selection on the slope (the cost of plasticity); Inf ⇒ no cost
    ReactionNorm::Bool          = false    # a heritable reaction norm: NoTraits = 2, trait 1 the elevation, trait 2 the slope
    # several traits with pleiotropy — Julia-only; with NoTraits = 1 the rest are unused
    NoTraits::Int               = 1        # n, traits each locus affects
    MutCorr::Float64            = 0.0      # r_μ, mutational correlation between every two traits
    SelCorr::Float64            = 0.0      # r_ω, correlation of the fitness function between every two traits
    MutCov::Matrix{Float64}     = zeros(0, 0)   # M, mutational covariance per generation (n×n); empty ⇒ m·V_E with MutCorr
    SelCov::Matrix{Float64}     = zeros(0, 0)   # Ω, initial fitness-function covariance (n×n); empty ⇒ ω² with SelCorr
    MoveDirection::Vector{Float64} = Float64[]  # u, the direction the optimum moves in; empty ⇒ trait 1
    FluctDirection::Vector{Float64} = Float64[] # f, the direction it fluctuates in (:periodic, default :red_noise); empty ⇒ u
    MoveNoiseCov::Matrix{Float64} = zeros(0, 0) # Σ, stationary covariance of the optimum's AR(1) noise (n×n); empty ⇒ MoveNoiseSD² f fᵀ
    Seed::Union{Nothing,Int}    = nothing
end

# Params holds arrays, so equality is by field values rather than by identity.
Base.:(==)(a::Params, b::Params) = all(f -> getfield(a, f) == getfield(b, f), fieldnames(Params))
Base.hash(p::Params, h::UInt) = foldl((h, f) -> hash(getfield(p, f), h), fieldnames(Params); init = hash(:Params, h))

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
`MoveType` one of `:phenotype`, `:fitness`, `:continuous`, `:red_noise`, `:periodic`
(QTSim falls back to the phenotype move on an unknown name; a misspelt symbol here would
silently change the model); `MoveNoiseSD ≥ 0`, `0 ≤ MoveAutocorr < 1`, a finite
`MoveAmplitude` and `MovePeriod > 0`; `CrossoverModel` `:qtsim` or `:all`; `Selection`
`:viability` or `:fecundity`; the schedules `MoveFitnessSDGen` and
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
    p.Demography in (:fixed, :beverton_holt) || bad("Demography must be :fixed or :beverton_holt, got :$(p.Demography)")
    p.GrowthRate > 0 || bad("GrowthRate must be positive, got $(p.GrowthRate)")
    p.Selection in (:viability, :fecundity) || bad("Selection must be :viability or :fecundity, got :$(p.Selection)")
    p.MoveType in (:phenotype, :fitness, :continuous, :red_noise, :periodic) ||
        bad("MoveType must be :phenotype, :fitness, :continuous, :red_noise or :periodic, got :$(p.MoveType)")
    p.MoveNoiseSD >= 0 || bad("MoveNoiseSD must be non-negative, got $(p.MoveNoiseSD)")
    0 <= p.MoveAutocorr < 1 || bad("MoveAutocorr must be in [0, 1), got $(p.MoveAutocorr)")
    isfinite(p.MoveAmplitude) || bad("MoveAmplitude must be finite, got $(p.MoveAmplitude)")
    p.MovePeriod > 0 && isfinite(p.MovePeriod) || bad("MovePeriod must be positive, got $(p.MovePeriod)")
    (p.MoveFitnessSDGen == -1 || p.MoveFitnessSDGen >= 0) ||
        bad("MoveFitnessSDGen must be -1 (no move) or a generation, got $(p.MoveFitnessSDGen)")
    (p.ChangeFitnessGen == -1 || p.ChangeFitnessGen >= 0) ||
        bad("ChangeFitnessGen must be -1 (no change) or a generation, got $(p.ChangeFitnessGen)")
    (p.ChangeFitnessGen == -1 || p.ChangeFitnessSD > 0) ||
        bad("ChangeFitnessSD must be positive when ChangeFitnessGen is set, got $(p.ChangeFitnessSD)")
    checktraits(p, bad)
    checkplasticity(p, bad)
    if p.InitFitnessSD == -1
        envVar(p) > 0 && p.MutationScalar > 0 ||
            bad("the analytical fitness width (InitFitnessSD = -1) needs Heritability < 1 and MutationScalar > 0; set InitFitnessSD")
    else
        p.InitFitnessSD > 0 || bad("InitFitnessSD must be -1 (analytical) or positive, got $(p.InitFitnessSD)")
    end
    p
end

# --- phenotypic plasticity ---------------------------------------------------

"""
    isplastic(p)

Whether `p` has a reaction norm: a fixed slope `Plasticity ≠ 0`, or a heritable one
(`ReactionNorm`).
"""
isplastic(p::Params) = p.Plasticity != 0 || p.ReactionNorm

"""
    plasticity_ratio(p) -> b/B

The fraction `Plasticity/OptimumSlope` of a displacement of the optimum that the plastic
response tracks (Chevin, Lande and Mace 2010: reaction norm slope `b`, environmental
sensitivity of selection `B`). 0 without plasticity, 1 for perfect plasticity.
"""
plasticity_ratio(p::Params) = p.Plasticity / p.OptimumSlope

"""
    plasticity_cost(p; b = p.Plasticity)

The factor `exp(−b²/(2ω_b²))` by which the cost of plasticity multiplies fitness, with
`ω_b = PlasticityCostSD`; 1 without a cost (`ω_b = Inf`) or without plasticity.
"""
plasticity_cost(p::Params; b::Real = p.Plasticity) = exp(-b^2 / (2 * p.PlasticityCostSD^2))

# checkparams for the plasticity fields
function checkplasticity(p::Params, bad)
    isfinite(p.Plasticity) || bad("Plasticity must be finite, got $(p.Plasticity)")
    isfinite(p.OptimumSlope) && p.OptimumSlope != 0 || bad("OptimumSlope must be finite and nonzero, got $(p.OptimumSlope)")
    p.CueLag >= 0 || bad("CueLag must be non-negative, got $(p.CueLag)")
    p.PlasticityCostSD > 0 || bad("PlasticityCostSD must be positive (Inf for no cost), got $(p.PlasticityCostSD)")
    if p.ReactionNorm
        p.NoTraits == 2 || bad("ReactionNorm needs NoTraits = 2 (trait 1 the elevation, trait 2 the slope), got NoTraits = $(p.NoTraits)")
        isempty(p.MoveDirection) && isempty(p.FluctDirection) && isempty(p.MoveNoiseCov) ||
            bad("with ReactionNorm the optimum is that of the phenotype, a number: leave MoveDirection, FluctDirection and MoveNoiseCov unset")
        p.SelCorr == 0 && isempty(p.SelCov) ||
            bad("with ReactionNorm selection acts on the one phenotype: leave SelCorr and SelCov unset (the cost is PlasticityCostSD)")
    elseif p.NoTraits > 1
        p.Plasticity == 0 || bad("Plasticity with NoTraits > 1 needs ReactionNorm = true (and NoTraits = 2)")
    end
    isplastic(p) && p.MoveFitnessSDGen != -1 && p.MoveType in (:phenotype, :fitness) &&
        bad("plasticity needs an optimum with its own trajectory (MoveType :continuous, :red_noise or :periodic); " *
            ":$(p.MoveType) re-anchors the optimum to the mean phenotype every generation")
    p
end

# --- several traits ---------------------------------------------------------

# the n×n matrix with unit diagonal and `r` off it
equicorrelation(n::Integer, r::Real) = [i == j ? 1.0 : Float64(r) for i in 1:n, j in 1:n]

issymmetric_approx(A) = size(A, 1) == size(A, 2) && all(isfinite, A) &&
                        all(abs(A[i, j] - A[j, i]) <= 1e-12 * max(1.0, abs(A[i, j])) for i in axes(A, 1), j in axes(A, 2))

mineigen(A) = minimum(eigvals(Symmetric(Matrix{Float64}(A))))

# checkparams for the trait fields: what a multi-trait model needs, and that a one-trait
# model does not silently ignore them
function checktraits(p::Params, bad)
    n = p.NoTraits
    n >= 1 || bad("NoTraits must be at least 1, got $n")
    if n == 1
        p.MutCorr == 0 && p.SelCorr == 0 && isempty(p.MutCov) && isempty(p.SelCov) && isempty(p.MoveDirection) &&
            isempty(p.FluctDirection) && isempty(p.MoveNoiseCov) ||
            bad("MutCorr, SelCorr, MutCov, SelCov, MoveDirection, FluctDirection and MoveNoiseCov are for NoTraits > 1; with one trait leave them unset")
        return p
    end
    -1 <= p.MutCorr <= 1 || bad("MutCorr must be in [-1, 1], got $(p.MutCorr)")
    -1 < p.SelCorr < 1 || bad("SelCorr must be in (-1, 1), got $(p.SelCorr)")
    if isempty(p.MutCov)
        mineigen(equicorrelation(n, p.MutCorr)) >= -1e-12 ||
            bad("MutCorr = $(p.MutCorr) is not a correlation of $n traits (it must be at least -1/(n-1))")
    else
        size(p.MutCov) == (n, n) && issymmetric_approx(p.MutCov) && mineigen(p.MutCov) >= -1e-12 ||
            bad("MutCov must be a symmetric positive semidefinite $n×$n matrix")
        p.MutCorr == 0 || bad("give MutCov or MutCorr, not both")
    end
    if isempty(p.SelCov)
        mineigen(equicorrelation(n, p.SelCorr)) > 0 ||
            bad("SelCorr = $(p.SelCorr) is not a correlation of $n traits (it must be above -1/(n-1))")
    else
        size(p.SelCov) == (n, n) && issymmetric_approx(p.SelCov) && mineigen(p.SelCov) > 0 ||
            bad("SelCov must be a symmetric positive definite $n×$n matrix")
        p.SelCorr == 0 || bad("give SelCov or SelCorr, not both")
    end
    isempty(p.MoveDirection) || (length(p.MoveDirection) == n && all(isfinite, p.MoveDirection) &&
                                 any(!iszero, p.MoveDirection)) ||
        bad("MoveDirection must be a nonzero vector of $n numbers, got $(p.MoveDirection)")
    isempty(p.FluctDirection) || (length(p.FluctDirection) == n && all(isfinite, p.FluctDirection) &&
                                  any(!iszero, p.FluctDirection)) ||
        bad("FluctDirection must be a nonzero vector of $n numbers, got $(p.FluctDirection)")
    if !isempty(p.MoveNoiseCov)
        size(p.MoveNoiseCov) == (n, n) && issymmetric_approx(p.MoveNoiseCov) && mineigen(p.MoveNoiseCov) >= -1e-12 ||
            bad("MoveNoiseCov must be a symmetric positive semidefinite $n×$n matrix")
        p.MoveNoiseSD == 0 || bad("give MoveNoiseCov or MoveNoiseSD, not both")
    end
    p.StandardisePhenotype && bad("StandardisePhenotype is for one trait")
    p
end

"""
    mutational_covariance(p)

The mutational covariance matrix `M` (n×n, n = `NoTraits`): the expected covariance of
the effects mutation adds to an individual's genotype per generation. `MutCov` when
given, otherwise `m·V_E` on the diagonal (as [`muteVar`](@ref) for one trait) and
`MutCorr·m·V_E` off it. Each mutation draws its effect vector from
`N(0, M/(Ploidy·NoQTL·MutationRatePerGamete))`, so the expected input is `M` (Lean
`mutationalVarianceInput_eq`, trait by trait and for every linear combination).
"""
mutational_covariance(p::Params) =
    isempty(p.MutCov) ? muteVar(p) .* equicorrelation(p.NoTraits, p.MutCorr) : Matrix{Float64}(p.MutCov)

"""
    selection_shape(p)

The shape `S` of the fitness function, `Ω = ω² S` (see [`selection_covariance`](@ref)):
unit diagonal and `SelCorr` off it, or `SelCov/ω₀²` with `ω₀` the initial width.
"""
selection_shape(p::Params) =
    isempty(p.SelCov) ? equicorrelation(p.NoTraits, p.SelCorr) : Matrix{Float64}(p.SelCov) ./ initFitSD(p)^2

"""
    selection_covariance(p; ω = initFitSD(p))

The covariance matrix `Ω` of the Gaussian fitness function `W(z) ∝ exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ))`
at width `ω`: `ω² S` with `S` the [`selection_shape`](@ref); initially `SelCov` when it is
given. For one trait it is `ω²`.
"""
selection_covariance(p::Params; ω::Real = initFitSD(p)) = ω^2 .* selection_shape(p)

"""
    move_direction(p)

The unit vector along which the optimum moves: `MoveDirection` normalised, or the first
trait's axis by default (so that one trait moves as in QTSim).
"""
function move_direction(p::Params)
    n = p.NoTraits
    isempty(p.MoveDirection) && return [i == 1 ? 1.0 : 0.0 for i in 1:n]
    u = Float64.(p.MoveDirection)
    u ./ sqrt(sum(abs2, u))
end

"""
    fluct_direction(p)

The unit vector along which a fluctuating optimum of several traits cycles (`:periodic`)
and, by default, wanders (`:red_noise`): `FluctDirection` normalised, or the move's
direction [`move_direction`](@ref) by default. For one trait, `[1.0]`.
"""
function fluct_direction(p::Params)
    isempty(p.FluctDirection) && return move_direction(p)
    f = Float64.(p.FluctDirection)
    f ./ sqrt(sum(abs2, f))
end

"""
    move_noise_covariance(p)

The stationary covariance `Σ` (n×n) of the AR(1) noise of a `:red_noise` optimum:
`MoveNoiseCov` when given, otherwise `MoveNoiseSD²·f fᵀ` with `f` the
[`fluct_direction`](@ref), noise along one direction. For one trait, `MoveNoiseSD²` as a
1×1 matrix.
"""
function move_noise_covariance(p::Params)
    isempty(p.MoveNoiseCov) || return Matrix{Float64}(p.MoveNoiseCov)
    f = fluct_direction(p)
    p.MoveNoiseSD^2 .* (f * f')
end
