# The life cycle, written once as a morphism of a free symmetric monoidal category, and
# its lowering to the simulation kernel (compose.jl) and to SLiM's Eidos (eidos.jl).
#
# The objects of the category are the stages of a generation and its generators are the
# steps between them:
#
#     develop:   Zygotes → Adults           record:  Offspring → Offspring
#     select:    Adults → Survivors         rescale: Offspring → Offspring
#     reproduce: Survivors → Offspring      mutate:  Offspring → Offspring
#                                           swap:    Offspring → Zygotes
#
# QTSim's generation is the composite develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅
# mutate ⋅ swap: Zygotes → Zygotes, in the order of its generation loop. The monoidal
# product ⊗ puts demes side by side: Zygotes⊗Zygotes is a metapopulation of two demes at
# the start of a generation, generation ⊗ generation lets each evolve on its own, and
# `migrate(k)`: Zygotes^⊗k → Zygotes^⊗k couples k of them. A functor out of a free
# category is fixed by where it sends the generators, so a backend is a table of images,
# one per step, and the image of any composite or product is the composite or product of
# the images (lean/CategoricalQuantGen/Lowering.lean, Demes.lean). Composites whose
# stages do not match are not morphisms at all: `select ⋅ develop` is a domain error.
#
# A step can carry parameters of its own, `withparams(select; ω = 0.3)`: a further
# generator, named by the step and its parameters (`Step`), whose image in each backend is
# the plain step's with those parameters in place of the deme's. Being a generator, it is
# covered by the same fact: a lowering is fixed by its images on the generators.

@present LifeCycle(FreeSymmetricMonoidalCategory) begin
    (Zygotes, Adults, Survivors, Offspring)::Ob
    develop::Hom(Zygotes, Adults)
    select::Hom(Adults, Survivors)
    reproduce::Hom(Survivors, Offspring)
    record::Hom(Offspring, Offspring)
    rescale::Hom(Offspring, Offspring)
    mutate::Hom(Offspring, Offspring)
    swap::Hom(Offspring, Zygotes)
end

const STEPS = (:develop, :select, :reproduce, :record, :rescale, :mutate, :swap)

"The number of stages, one per deme, in a product of stages."
nfactors(x) = head(x) === :otimes ? length(args(x)) : head(x) === :munit ? 0 : 1

"Zygotes⊗…⊗Zygotes with `k` factors: `k` demes at the start of a generation."
zygotes_word(k::Integer) = k == 1 ? LifeCycle[:Zygotes] : foldl(otimes, fill(LifeCycle[:Zygotes], k))

"""
    Step(name, params)

The name of a step that carries its own parameters: a generator `Hom(Step(name, params),
dom, codom)` of the free category, beside the presentation's generator `name`.
`params` is a named tuple of the fields of `Params` the step overrides, or of the
step's own parameters (`ω` and `θ` for `select`, `rate` for `migrate`). Two such
generators are equal when their names and parameters are, so the parameters are part
of the syntax: `withparams(L.select; ω = 0.3)` is a different morphism from `L.select`
and from `withparams(L.select; ω = 0.4)`.
"""
struct Step
    name::Symbol
    params::NamedTuple
end

Base.:(==)(a::Step, b::Step) = a.name == b.name && a.params == b.params
Base.hash(s::Step, h::UInt) = hash(s.params, hash(s.name, hash(:Step, h)))
showvalue(x) = x isa AbstractMatrix ? repr(x; context = :compact => true) : repr(x)
Base.show(io::IO, s::Step) =
    print(io, s.name, "(", join(("$k = $(showvalue(v))" for (k, v) in pairs(s.params)), ", "), ")")

"""
The parameters each step of the life cycle can carry: fields of `Params` that its
kernel image reads, and for `select` a fixed width `ω` and a fixed optimum `θ`.
"""
const STEP_PARAMS = Dict{Symbol,Tuple}(
    :develop => (),
    :select => (:ω, :θ, :MoveType, :MoveFitnessSD, :MoveFitnessSDGen, :MoveNoiseSD, :MoveAutocorr,
                :MoveAmplitude, :MovePeriod, :StandardisePhenotype, :Drift),
    :reproduce => (:RecombPoissonMean, :CrossoverModel, :Demography, :GrowthRate),
    :record => (:Sample,),
    :rescale => (:ChangeFitnessGen, :ChangeFitnessSD),
    :mutate => (:MutationScalar, :MutationRatePerGamete),
    :swap => (),
    :migrate => (:rate,))

const MOVE_FIELDS = (:MoveType, :MoveFitnessSD, :MoveFitnessSDGen, :MoveNoiseSD, :MoveAutocorr,
                     :MoveAmplitude, :MovePeriod)

"The name of a generator: a step's name, whether or not it carries parameters."
stepname(g) = (v = first(args(g)); v isa Step ? v.name : nameof(g))

"The parameters a generator carries (an empty named tuple for a plain step)."
stepparams(g) = (v = first(args(g)); v isa Step ? v.params : (;))

"The fields of `Params` among a step's parameters."
paramfields(s::Step) = NamedTuple(k => v for (k, v) in pairs(s.params) if k in fieldnames(Params))

# Check a step's parameters and normalise their types, through Params' own constructor and
# checks: `ω = 1` becomes 1.0, and `MoveType = :bogus` is an error when it is written.
function checkstep(name::Symbol, params::NamedTuple)
    allowed = get(STEP_PARAMS, name, nothing)
    allowed === nothing && throw(ArgumentError("the step $name has no parameters; only the steps of lifecycle() do"))
    for k in keys(params)
        k in allowed || throw(ArgumentError("the step $name takes the parameters $(join(allowed, ", ")); got $k"))
    end
    haskey(params, :θ) && any(k -> k in MOVE_FIELDS, keys(params)) &&
        throw(ArgumentError("a fixed optimum θ replaces the move; give one or the other to $name"))
    for k in (:ω, :θ)
        haskey(params, k) && !(params[k] isa Real && isfinite(params[k])) &&
            throw(ArgumentError("$k must be a finite number; got $(params[k])"))
    end
    haskey(params, :ω) && params.ω <= 0 && throw(ArgumentError("ω must be positive; got $(params.ω)"))
    fields = NamedTuple(k => v for (k, v) in pairs(params) if k in fieldnames(Params))
    q = isempty(fields) ? Params() : QTLSim.checkparams(Params(Params(; InitFitnessSD = 1.0, ChangeFitnessSD = 1.0); fields...))
    NamedTuple(k => (k in fieldnames(Params) ? getfield(q, k) :
                     k === :rate ? (v isa AbstractMatrix ? Matrix{Float64}(v) : Float64(v)) : Float64(v))
               for (k, v) in pairs(params))
end

function withparams(f::GATExpr; kw...)
    isempty(kw) && return f
    used = Set{Symbol}()
    gens = Dict{Any,Any}()
    for g in stepsof(f)
        g isa ObExpr && (gens[g] = g; continue)
        name = stepname(g)
        mine = NamedTuple(k => v for (k, v) in pairs(kw) if k in get(STEP_PARAMS, name, ()))
        if isempty(mine)
            gens[g] = g
        else
            union!(used, keys(mine))
            merged = merge(stepparams(g), mine)
            if name === :migrate && haskey(merged, :rate)
                rate = merged.rate
                rate isa AbstractMatrix ? QTLSim.check_migration(rate, nfactors(dom(g))) :
                    (0 <= rate <= 1 || throw(ArgumentError("a migration rate is in [0, 1]; got $rate")))
            end
            gens[g] = Hom(Step(name, checkstep(name, merged)), dom(g), codom(g))
        end
    end
    unused = setdiff(keys(kw), used)
    isempty(unused) ||
        throw(ArgumentError("no step of $f takes the parameter$(length(unused) == 1 ? "" : "s") $(join(unused, ", "))"))
    head(f) === :generator ? gens[f] :
        functor((FreeSymmetricMonoidalCategory.Ob, FreeSymmetricMonoidalCategory.Hom), f; generators = gens)
end

"The step that exchanges offspring among `k` demes, at its own `rate` when one is given."
function migrate_step(k::Integer; rate = nothing)
    k >= 2 || throw(ArgumentError("migration needs at least two demes; got $k"))
    g = Hom(:migrate, zygotes_word(k), zygotes_word(k))
    rate === nothing ? g : withparams(g; rate)
end

function lifecycle()
    steps = map(s -> LifeCycle[s], STEPS)
    NamedTuple{(:presentation, STEPS..., :generation, :migrate)}(
        (LifeCycle, steps..., foldl(⋅, steps), migrate_step))
end

"""
A product of stages, one per deme, as the image of an object under the lowering to
backend `B` (`:kernel`, `:eidos` or `:moments`). The kernel carries each deme's stages in
one `World` and the moments lowering in one `Moments`; in Eidos, a stage is the set of
variables a step can rely on (`EIDOS_SCOPE`).
"""
struct Stage{B}
    names::Vector{Symbol}
end

Stage{B}(name::Symbol) where {B} = Stage{B}([name])
nfactors(s::Stage) = length(s.names)

otimes(a::Stage{B}, b::Stage{B}) where {B} = Stage{B}(vcat(a.names, b.names))
otimes(a::Stage{B}, b::Stage{B}, c::Stage{B}, rest::Stage{B}...) where {B} = foldl(otimes, (a, b, c, rest...))
munit(::Type{Stage{B}}) where {B} = Stage{B}(Symbol[])

id(s::Stage{:kernel}) = nfactors(s) == 1 ? Operator(:id, (w, rng) -> w) :
                        DemeOp(nfactors(s), collect(1:nfactors(s)), (st, d, rng) -> st)
id(s::Stage{:eidos}) = nfactors(s) == 1 ? EidosCode(Symbol(:id_, s.names[1]), "") :
                       EidosMeta(nfactors(s), collect(1:nfactors(s)), pops -> "", :id)

"Where a braiding X⊗Y → Y⊗X sends the demes of X (m of them) and Y (n of them)."
braid_perm(m, n) = vcat(collect(m+1:m+n), collect(1:m))
braid(a::Stage{:kernel}, b::Stage{:kernel}) =
    DemeOp(nfactors(a) + nfactors(b), braid_perm(nfactors(a), nfactors(b)), (st, d, rng) -> st)
braid(a::Stage{:eidos}, b::Stage{:eidos}) =
    EidosMeta(nfactors(a) + nfactors(b), braid_perm(nfactors(a), nfactors(b)), pops -> "", :braid)

"The generators occurring in a life-cycle expression."
stepsof(f) = head(f) === :generator ? Any[f] :
             reduce(vcat, (stepsof(a) for a in args(f) if a isa GATExpr); init = Any[])

asimage(::Val{:kernel}, name, x::Union{Operator,DemeOp}) = x
asimage(::Val{:kernel}, name, x) = Operator(name, x)               # a function (state, rng) -> state
asimage(::Val{:eidos}, name, x::Union{EidosCode,EidosMeta}) = x
asimage(::Val{:eidos}, name, x::AbstractString) = EidosCode(name, x)
asimage(::Val{:moments}, name, x::Operator) = x
asimage(::Val{:moments}, name, x) = Operator(name, x)               # a function (moments, rng) -> moments

# The image of migration among `k` demes at `rate`, a number or a matrix.
migrate_image(::Val{:kernel}, k, rate) = migrate_op(k, rate)
migrate_image(::Val{:eidos}, k, rate) = eidos_migrate(k, rate isa AbstractMatrix)
migrate_image(::Val{:moments}, k, rate) = moments_migrate_op(k, rate)

"""
The image of a step that carries parameters (`s`, a `Step`): the one `table` gives for
`s` itself, if any; for `migrate`, migration at the step's rate; otherwise the image of
the plain step `s.name`, with the step's parameters in place of those of `Params` (see
`with_step_params` and `eidos_with_params`).
"""
function step_image(t::Val{T}, s::Step, table, k, migration_rate) where {T}
    haskey(table, s) && return asimage(t, Symbol(string(s)), table[s])
    if s.name === :migrate
        haskey(s.params, :rate) || return migrate_image(t, k, migration_rate)
        return T === :eidos ? eidos_migrate_at(k, s.params.rate) : migrate_image(t, k, s.params.rate)
    end
    haskey(table, s.name) ||
        throw(ArgumentError("no $T image for the step $s; give one with `images = Dict(:$(s.name) => …)`"))
    base = asimage(t, s.name, table[s.name])
    T === :eidos ? eidos_with_params(s, base) : with_step_params(s, base)
end

function lower(f, target::Symbol; images = Dict{Symbol,Any}(), record = nothing,
               migration_rate::Union{Real,AbstractMatrix} = 0.0, deme_constants = ())
    table = target === :kernel ? kernel_images(record) :
            target === :eidos ? eidos_images() :
            target === :moments ? moment_images(record) :
            throw(ArgumentError("unknown target :$target; the targets are :kernel, :eidos and :moments"))
    table = merge!(Dict{Any,Any}(table), Dict{Any,Any}(images))
    gens = Dict{Any,Any}()
    for g in stepsof(f)
        v = first(args(g))
        name = v isa Step ? v : nameof(g)
        if g isa ObExpr
            gens[g] = Stage{target}(name)
        elseif v isa Step
            gens[g] = step_image(Val(target), v, table, nfactors(dom(g)), migration_rate)
        elseif haskey(table, name)
            gens[g] = asimage(Val(target), name, table[name])
        elseif name === :migrate
            gens[g] = migrate_image(Val(target), nfactors(dom(g)), migration_rate)
        else
            throw(ArgumentError("no $target image for the step $name; give one with `images = Dict(:$name => …)`"))
        end
        # slim_script with demes of their own parameters: those constants are deme state
        target === :eidos && !isempty(deme_constants) && gens[g] isa EidosCode &&
            (gens[g] = EidosCode(gens[g].name, deme_state(gens[g].code, deme_constants)))
    end
    functor((Stage{target}, Any), f; generators = gens)
end

function operators()
    k = kernel_images()
    NamedTuple{(STEPS..., :generation)}((map(s -> k[s], STEPS)..., lower(lifecycle().generation, :kernel)))
end
