# Models: the families of lowerings of the one life cycle. The life cycle (lifecycle(), in
# the categorical layer) is a morphism of a free symmetric monoidal category, and what its
# steps mean is fixed by a functor out of that category. QTLSim has two models, two
# families of functors from the same syntax: QTSim's one trait (`OneTraitModel`, the
# object image ℝ) and n traits with pleiotropy (`TraitModel(n)`, the object image ℝⁿ).
# A model fixes three things:
#
# - the parameters each step can carry (`step_fields`, `step_own_params`): a step with
#   parameters, `withparams(select, model; θ = …)`, is a generator whose parameters range
#   over the model's parameter space, so `θ` is a number for one trait and a vector of n
#   for n traits, and `Ω` (a fitness covariance) exists only for n traits;
# - the table of images of the steps in each backend (the kernel, Eidos and the moments,
#   in ext/QTLSimCatlabExt), which `lower(f, target; model)` chooses;
# - how a step's parameters act in the image (`with_step_params` in the extension).
#
# This is deliberately concrete: two models, written out. It is not the general framework
# (an abstract category of evolutionary systems of which these are instances), but it keeps
# the parts that differ between models in one place, so that a third model is a third
# set of methods. See docs/src/traits_design.md.

"""
    QTLSim.Model

A model of QTLSim: a family of lowerings of the one life cycle, which fixes the parameters
each step can carry ([`step_schema`](@ref)) and the images of the steps in each backend.
The models are [`OneTraitModel`](@ref) (QTSim's) and [`TraitModel`](@ref); `model_of(p)`
is the model of a `Params`.
"""
abstract type Model end

"""
    OneTraitModel()

QTSim's model: one trait, whose genotype, phenotype and optimum are numbers. The default
model of [`withparams`](@ref) and [`lower`](@ref).
"""
struct OneTraitModel <: Model end

"""
    TraitModel(n)

`n ≥ 2` traits with pleiotropy (`Params` `NoTraits = n`): genotypes, phenotypes and the
optimum are vectors of `n`, the fitness function has a covariance matrix `Ω`, and mutation
a covariance matrix `M`. See `docs/src/traits_design.md`.
"""
struct TraitModel <: Model
    n::Int
    function TraitModel(n::Integer)
        n >= 2 || throw(ArgumentError("a TraitModel has at least two traits; one trait is OneTraitModel()"))
        new(n)
    end
end

"""
    ReactionNormModel()

A heritable linear reaction norm (`Params` `ReactionNorm = true`, `NoTraits = 2`): an
individual carries two heritable values, as in [`TraitModel`](@ref)`(2)`, its elevation
and (the deviation from `Plasticity` of) its slope, but has one phenotype,
`z = a + b·ε + e`, on which selection acts, with the optimum a number as in
[`OneTraitModel`](@ref). Its images are `TraitModel(2)`'s for every step but `select`
(and, in the moments, `record`), which read the pair through the linear map
`(a, b) ↦ a + b·ε` of the generation's environment. See `docs/src/plasticity_design.md`.
"""
struct ReactionNormModel <: Model end

"""
    model_of(p::Params) -> Model

The model a `Params` belongs to: [`OneTraitModel`](@ref) when `NoTraits = 1`,
[`ReactionNormModel`](@ref) when `ReactionNorm` is set, otherwise
[`TraitModel`](@ref)`(NoTraits)`.
"""
model_of(p::Params) = p.NoTraits == 1 ? OneTraitModel() : p.ReactionNorm ? ReactionNormModel() : TraitModel(p.NoTraits)
model_of(m::Model) = m

"The number of traits of a model."
ntraits(::OneTraitModel) = 1
ntraits(m::TraitModel) = m.n
ntraits(::ReactionNormModel) = 2

# the fields of Params that place a default Params in a model (for checking a step's fields)
model_fields(m::Model) = (NoTraits = ntraits(m),)
model_fields(::ReactionNormModel) = (NoTraits = 2, ReactionNorm = true)

Base.show(io::IO, ::OneTraitModel) = print(io, "OneTraitModel()")
Base.show(io::IO, m::TraitModel) = print(io, "TraitModel(", m.n, ")")
Base.show(io::IO, ::ReactionNormModel) = print(io, "ReactionNormModel()")

# the fields of Params of a reaction norm, which `select` reads: the optimum of the generation,
# the measure of the environment, is set there
const PLASTICITY_FIELDS = (:Plasticity, :OptimumSlope, :CueLag, :PlasticityCostSD)

# the fields of Params a step's image reads, in every model
const STEP_FIELDS_SHARED = Dict{Symbol,Tuple}(
    :develop => (),
    :select => (:MoveType, :MoveFitnessSD, :MoveFitnessSDGen, :MoveNoiseSD, :MoveAutocorr,
                :MoveAmplitude, :MovePeriod, :Drift, :Selection),
    :reproduce => (:RecombPoissonMean, :CrossoverModel, :Demography, :GrowthRate),
    :record => (:Sample,),
    :rescale => (:ChangeFitnessGen, :ChangeFitnessSD),
    :mutate => (:MutationScalar, :MutationRatePerGamete),
    :swap => (),
    :migrate => ())

"""
    step_fields(model) -> Dict{Symbol,Tuple}

The fields of `Params` each step of the life cycle can override in `model` (see
[`step_schema`](@ref)).
"""
step_fields(::OneTraitModel) =
    merge(STEP_FIELDS_SHARED, Dict{Symbol,Tuple}(:select => (STEP_FIELDS_SHARED[:select]..., :StandardisePhenotype,
                                                              PLASTICITY_FIELDS...)))
step_fields(::ReactionNormModel) =
    merge(STEP_FIELDS_SHARED, Dict{Symbol,Tuple}(
        :select => (STEP_FIELDS_SHARED[:select]..., PLASTICITY_FIELDS...),
        :mutate => (STEP_FIELDS_SHARED[:mutate]..., :MutCorr, :MutCov)))
step_fields(::TraitModel) =
    merge(STEP_FIELDS_SHARED, Dict{Symbol,Tuple}(
        :select => (STEP_FIELDS_SHARED[:select]..., :SelCorr, :MoveDirection, :MoveNoiseCov, :FluctDirection),
        :mutate => (STEP_FIELDS_SHARED[:mutate]..., :MutCorr, :MutCov)))

# checkers of a step's own parameters (not fields of Params): each normalises a value or
# throws an ArgumentError naming the parameter
function check_real(k, v; positive = false)
    v isa Real && isfinite(v) || throw(ArgumentError("$k must be a finite number; got $v"))
    positive && v <= 0 && throw(ArgumentError("$k must be positive; got $v"))
    Float64(v)
end
function check_vector(k, v, n)
    v isa AbstractVector{<:Real} && length(v) == n && all(isfinite, v) ||
        throw(ArgumentError("$k must be a vector of $n finite numbers, one per trait; got $v"))
    Vector{Float64}(v)
end
function check_pd(k, v, n)
    v isa AbstractMatrix{<:Real} && size(v) == (n, n) && issymmetric_approx(v) && mineigen(v) > 0 ||
        throw(ArgumentError("$k must be a symmetric positive definite $n×$n matrix; got $v"))
    Matrix{Float64}((v + v') / 2)
end
check_rate(k, v) = v isa AbstractMatrix ? Matrix{Float64}(v) : check_real(k, v)

"""
    step_own_params(model) -> Dict{Symbol,NamedTuple}

The parameters each step can carry in `model` that are not fields of `Params`, each with
the function that checks and normalises its value: for one trait, `select`'s fixed width
`ω` and fixed optimum `θ`, numbers; for `n` traits, `ω`, a vector `θ` of `n` and a fitness
covariance matrix `Ω` (in place of `ω² S`); in both, `migrate`'s `rate`.
"""
step_own_params(::OneTraitModel) = Dict{Symbol,NamedTuple}(
    :select => (ω = v -> check_real(:ω, v; positive = true), θ = v -> check_real(:θ, v)),
    :migrate => (rate = v -> check_rate(:rate, v),))
step_own_params(::ReactionNormModel) = step_own_params(OneTraitModel())
step_own_params(m::TraitModel) = Dict{Symbol,NamedTuple}(
    :select => (ω = v -> check_real(:ω, v; positive = true), θ = v -> check_vector(:θ, v, m.n),
                Ω = v -> check_pd(:Ω, v, m.n)),
    :migrate => (rate = v -> check_rate(:rate, v),))

"""
    step_schema(model) -> Dict{Symbol,Tuple}

The parameters each step of the life cycle can carry in `model`, by step: the fields of
`Params` its image reads ([`step_fields`](@ref)) and its own parameters
([`step_own_params`](@ref)). [`withparams`](@ref) checks a step's parameters against the
schema of a model, and [`lower`](@ref) against the schema of the model it lowers for, so a
step written for one model is rejected by the other when its parameters do not fit (a
number `θ` for several traits, a matrix `Ω` for one).
"""
function step_schema(m::Model)
    f, o = step_fields(m), step_own_params(m)
    Dict{Symbol,Tuple}(s => (keys(get(o, s, (;)))..., f[s]...) for s in keys(f))
end

# the fields of Params that make up the optimum's move (a fixed θ replaces them)
const MOVE_FIELDS = (:MoveType, :MoveFitnessSD, :MoveFitnessSDGen, :MoveNoiseSD, :MoveAutocorr,
                     :MoveAmplitude, :MovePeriod, :MoveNoiseCov, :FluctDirection)

"""
    QTLSim.checkstep(model, name, params) -> NamedTuple

Check the parameters `params` of the step `name` against `model`'s schema and normalise
their values: own parameters by their checkers, fields of `Params` through `Params`' own
constructor and [`checkparams`](@ref) (so `ω = 1` becomes `1.0` and `MoveType = :bogus` is
an error when it is written). A fixed `θ` excludes the move's fields and `Ω` excludes `ω`
and the fitness shape `SelCorr`.
"""
function checkstep(m::Model, name::Symbol, params::NamedTuple)
    fields = get(step_fields(m), name, nothing)
    fields === nothing && throw(ArgumentError("the step $name has no parameters; only the steps of lifecycle() do"))
    own = get(step_own_params(m), name, (;))
    allowed = (keys(own)..., fields...)
    for k in keys(params)
        k in allowed ||
            throw(ArgumentError("in $(m), the step $name takes the parameters $(join(allowed, ", ")); got $k"))
    end
    haskey(params, :θ) && any(k -> k in MOVE_FIELDS, keys(params)) &&
        throw(ArgumentError("a fixed optimum θ replaces the move; give one or the other to $name"))
    haskey(params, :Ω) && any(k -> k in (:ω, :SelCorr), keys(params)) &&
        throw(ArgumentError("Ω is the whole fitness covariance; give it or ω and SelCorr to $name"))
    pf = NamedTuple(k => v for (k, v) in pairs(params) if k in fields)
    q = isempty(pf) ? nothing :
        checkparams(Params(Params(; InitFitnessSD = 1.0, ChangeFitnessSD = 1.0, model_fields(m)...); pf...))
    NamedTuple(k => (haskey(own, k) ? own[k](v) : getfield(q, k)) for (k, v) in pairs(params))
end
