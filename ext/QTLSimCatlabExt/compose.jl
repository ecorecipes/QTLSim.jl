# Evolutionary processes as first-class operators: the image of the life cycle
# (lifecycle.jl) in the simulation kernel. An operator acts on a `World` and composes
# left to right with Catlab's `⋅`, `(f ⋅ g)(w) = g(f(w))`, which reads in the biological
# order and is the reverse of the `∘` in Composition.lean. Catlab's monoidal product `⊗`
# builds a metapopulation of demes over which the same operators act pointwise, with a
# migration operator coupling them; a life cycle on a tensor product of stages, one per
# deme, lowers to a `DemeOp` on such a metapopulation.

"""
    Operator(name, act, scope)

An evolutionary operator: `op(state, rng)` applies `act(state, rng)`. `scope` is
`:deme` for an operator on one population, applied to each deme of a metapopulation,
or `:meta` for one on the whole metapopulation (migration).
"""
struct Operator{F}
    name::Symbol
    act::F
    scope::Symbol
end

Operator(name::Symbol, act) = Operator(name, act, :deme)
operator(name::Symbol, act; scope::Symbol = :deme) = Operator(name, act, scope)

(op::Operator)(state, rng) = op.act(state, rng)

Base.show(io::IO, op::Operator) = print(io, "Operator(:", op.name, ")")

"""
    f ⋅ g

Apply `f`, then `g`; the rest of a composite is skipped once the population is extinct
(as in `run!`). A composite acts on the whole metapopulation if either factor does.
"""
function ⋅(f::Operator, g::Operator)
    sc = (f.scope === :meta || g.scope === :meta) ? :meta : :deme
    Operator(Symbol(f.name, :⋅, g.name), function (state, rng)
        f(state, rng)
        isextinct(state) && return state
        g(state, rng)
    end, sc)
end

# Catlab's `compose` is `⋅` on operators, so that a functor from the life-cycle category
# can land in them: GATlab's `functor` composes the images of the steps with `compose`.
compose(f::Operator, g::Operator) = f ⋅ g
compose(f::Operator, g::Operator, h::Operator, rest::Operator...) = foldl(⋅, (f, g, h, rest...))

isextinct(w::World) = w.extinct

# --- the kernel image of each step of the life cycle ---------------------
#
# Each model (QTLSim.Model) is a functor out of the same free category, so each has its
# own table of images. The tables share the images of the steps that do not care what an
# allele is (inheritance copies blocks of loci, the width changes on schedule, the
# buffers swap) and differ in develop, select and mutate, where the object image, ℝ or ℝⁿ,
# matters. An image checks that it got a population of its model.

# the steps whose kernel image is the same in both models
shared_kernel_images(rec) = Dict{Symbol,Any}(
    :reproduce => Operator(:reproduce, (w, rng) -> isextinct(w) ? w : QTLSim.reproduce!(w, rng)),
    :record => Operator(:record, (w, rng) -> record!(rec, w)),
    :rescale => Operator(:rescale, (w, rng) -> rescale!(w)),
    :swap => Operator(:swap, (w, rng) -> (QTLSim.swapbuffers!(w); w.gen += 1; w)))

"A `World` of the one-trait model, or an error that says which model it belongs to."
@inline function onetrait(w::World)
    w.traits === nothing ||
        throw(ArgumentError("the one-trait model's image of a step got a population of $(w.p.NoTraits) traits; " *
                            "lower with model = model_of(p), TraitModel($(w.p.NoTraits))"))
    w
end

"A `World` of `n` traits, or an error that says which model it belongs to."
@inline function ntraits_world(w::World, n::Int)
    t = w.traits
    w.p.ReactionNorm && throw(ArgumentError("the image of a step for $n traits got a population with a heritable reaction norm; lower with model = model_of(p), ReactionNormModel()"))
    t !== nothing && t.n == n ||
        throw(ArgumentError("the image of a step for $n traits got a population of $(w.p.NoTraits) trait" *
                            "$(w.p.NoTraits == 1 ? "" : "s"); lower with model = model_of(p)"))
    w
end

"""
    kernel_images(model, rec = nothing)

The kernel image of each step of the life cycle in `model`, by name: the kernel's
operators, with `record` appending to `rec` (see `QTLSim.records`), or discarding its
observations when `rec` is `nothing`. `OneTraitModel()` (QTSim's) and `TraitModel(n)` share
the images of reproduce, record, rescale and swap, and differ in develop, select and mutate.
"""
kernel_images(::OneTraitModel, rec = nothing) = merge!(shared_kernel_images(rec), Dict{Symbol,Any}(
    :develop => Operator(:develop, (w, rng) -> QTLSim.develop_one!(onetrait(w))),
    :select => Operator(:select, (w, rng) -> QTLSim.select_one!(onetrait(w), rng)),
    :mutate => Operator(:mutate, (w, rng) -> QTLSim.mutate_one!(onetrait(w), rng))))
kernel_images(m::TraitModel, rec = nothing) = merge!(shared_kernel_images(rec), Dict{Symbol,Any}(
    :develop => Operator(:develop, (w, rng) -> QTLSim.develop_traits!(ntraits_world(w, m.n))),
    :select => Operator(:select, (w, rng) -> QTLSim.select_traits!(ntraits_world(w, m.n), rng)),
    :mutate => Operator(:mutate, (w, rng) -> QTLSim.mutate_traits!(ntraits_world(w, m.n), rng))))
# a heritable reaction norm: TraitModel(2)'s images, but for select, which reads the pair
# (elevation, slope) through the environment, z = a + b·ε + e (QTLSim.select_norm!)
kernel_images(::ReactionNormModel, rec = nothing) = merge!(shared_kernel_images(rec), Dict{Symbol,Any}(
    :develop => Operator(:develop, (w, rng) -> QTLSim.develop_traits!(normworld(w))),
    :select => Operator(:select, (w, rng) -> QTLSim.select_norm!(normworld(w), rng)),
    :mutate => Operator(:mutate, (w, rng) -> QTLSim.mutate_traits!(normworld(w), rng))))
kernel_images(rec = nothing) = kernel_images(OneTraitModel(), rec)

"A `World` with a heritable reaction norm, or an error that says which model it belongs to."
@inline function normworld(w::World)
    w.p.ReactionNorm && w.traits !== nothing ||
        throw(ArgumentError("the reaction-norm model's image of a step got a population of $(model_of(w.p)); lower with model = model_of(p)"))
    w
end

"""
    with_step_params(model, s::Step, op::Operator)

The image of a step with parameters `s.params`, from the image `op` of the plain step,
in the kernel or the moments lowering of `model`: `op` runs with the deme's `Params`
overridden by the step's fields of `Params`, and with the step's own parameters in place of
the deme's state, all restored afterwards. A step's own `θ` replaces the move
(`MoveFitnessSDGen = -1` for that step). So the parameters act in that step only: the
deme's state keeps its own optimum and width, which is what `record` records, while a move
made by a step's own move parameters moves the deme's optimum.

For one trait the step's `ω` and `θ` (numbers) replace the deme's width and optimum. For
`n` traits `θ` is a vector, `Ω` a fitness covariance (in place of `ω² S`, so the width is
1 for that step), and the fields that define the state's derived matrices replace them:
`SelCorr` the fitness shape, `MutationScalar`, `MutationRatePerGamete`,
`MutCorr`/`MutCov` the mutational covariance, `MoveDirection` the move's direction,
`MoveNoiseSD`, `MoveNoiseCov` and `FluctDirection` the fluctuation's noise and direction.
"""
function with_step_params(::OneTraitModel, s, op::Operator)
    fields = paramfields(s)
    ω = get(s.params, :ω, nothing)
    θ = get(s.params, :θ, nothing)
    θ === nothing || (fields = merge(fields, (MoveFitnessSDGen = -1,)))
    cache = Tuple{Params,Params}[]           # the overridden Params, per deme's Params
    Operator(Symbol(string(s)), function (w, rng)
        p0, ω0, θ0 = w.p, w.ω, w.θ
        if !isempty(fields)
            i = cacheindex(cache, p0)
            i == 0 && (push!(cache, (p0, QTLSim.checkparams(Params(p0; fields...)))); i = length(cache))
            w.p = cache[i][2]
        end
        ω === nothing || (w.ω = ω)
        θ === nothing || (w.θ = θ)
        try
            op(w, rng)
        finally
            w.p = p0
            ω === nothing || (w.ω = ω0)
            θ === nothing || (w.θ = θ0)
        end
        w
    end, op.scope)
end

# A heritable reaction norm: its optimum and width are numbers, as for one trait, and its
# mutation is that of two traits. So a step's `ω` and `θ` (numbers) and its fields of Params
# act as for one trait, on the optimum of the phenotype (`norm_optimum`), and a `mutate`
# step's fields as for `TraitModel(2)`.
function with_step_params(::ReactionNormModel, s, op::Operator)
    s.name === :mutate && return with_step_params(TraitModel(2), s, op)
    fields = paramfields(s)
    ω = get(s.params, :ω, nothing)
    θ = get(s.params, :θ, nothing)
    θ === nothing || (fields = merge(fields, (MoveFitnessSDGen = -1,)))
    cache = Tuple{Params,Params}[]
    Operator(Symbol(string(s)), function (w, rng)
        p0, ω0, θ0 = w.p, w.ω, norm_optimum(w)
        if !isempty(fields)
            i = cacheindex(cache, p0)
            i == 0 && (push!(cache, (p0, QTLSim.checkparams(Params(p0; fields...)))); i = length(cache))
            w.p = cache[i][2]
        end
        ω === nothing || (w.ω = ω)
        θ === nothing || set_norm_optimum!(w, θ)
        try
            op(w, rng)
        finally
            w.p = p0
            ω === nothing || (w.ω = ω0)
            θ === nothing || set_norm_optimum!(w, θ0)
        end
        w
    end, op.scope)
end
norm_optimum(w::World) = w.θ
set_norm_optimum!(w::World, θ) = (w.θ = θ; w)

# the index in `cache`, a vector of tuples whose first element is a deme's Params, of the
# entry for `p` (the same object, ===), or 0; a linear search allocates nothing, where an
# IdDict boxes its immutable key
@inline function cacheindex(cache::Vector, p::Params)
    for i in eachindex(cache)
        cache[i][1] === p && return i
    end
    0
end

function with_step_params(::TraitModel, s, op::Operator)
    θ = get(s.params, :θ, nothing)
    Ω = get(s.params, :Ω, nothing)
    fields = θ === nothing ? paramfields(s) : merge(paramfields(s), (MoveFitnessSDGen = -1,))
    ω = Ω === nothing ? get(s.params, :ω, nothing) : 1.0      # Ω is the whole covariance: Ω = 1² Ω
    # per deme's Params: the overridden Params and the state's matrices for the step, all
    # derived from them (where the step overrides none of their fields, the deme's values)
    cache = Tuple{Params,Params,Tuple{Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Vector{Float64}}}[]
    Operator(Symbol(string(s)), function (w, rng)
        p0, ω0 = w.p, w.ω
        i = cacheindex(cache, p0)
        if i == 0
            q = isempty(fields) ? p0 : QTLSim.checkparams(Params(p0; fields...))
            push!(cache, (p0, q, derived_slots(w, q, Ω)))
            i = length(cache)
        end
        _, q, d = cache[i]
        saved = trait_slots(w)               # (θ, selection, mutation, direction, noise, cycle direction)
        w.p = q
        ω === nothing || (w.ω = ω)
        set_trait_slots!(w, (θ === nothing ? saved[1] : θ, d...))
        try
            op(w, rng)
        finally
            w.p = p0
            ω === nothing || (w.ω = ω0)
            # the optimum stays where the step's own move put it, unless the step fixed it
            set_trait_slots!(w, (θ === nothing ? trait_slots(w)[1] : saved[1], Base.tail(saved)...))
        end
        w
    end, op.scope)
end

# The state of n traits that a step's parameters replace, in the kernel: the optimum
# vector, the inverse of the fitness shape, a square root of one mutation's covariance, the
# move's direction, a square root of the fluctuation's noise covariance and the cycle's
# direction (fields of the World's TraitState). w.θ holds the optimum's trait 1.
trait_slots(w::World) = (t = w.traits; (t.θ, t.Sinv, t.Lmut, t.u, t.Lnoise, t.f))
function set_trait_slots!(w::World, (θ, Sinv, Lmut, u, Lnoise, f))
    t = w.traits
    t.θ, t.Sinv, t.Lmut, t.u, t.Lnoise, t.f = θ, Sinv, Lmut, u, Lnoise, f
    w.θ = θ[1]
    w
end
# the matrices of the step's Params (with Ω, width 1 and shape Ω) as the kernel holds them,
# computed as TraitState computes them
derived_slots(::World, q::Params, Ω) =
    (Matrix(inv(QTLSim.Symmetric(Ω === nothing ? QTLSim.selection_shape(q) : Ω))),
     QTLSim.psd_sqrt(QTLSim.mutational_covariance(q) ./ (q.Ploidy * q.NoQTL * q.MutationRatePerGamete)),
     QTLSim.move_direction(q), QTLSim.psd_sqrt(QTLSim.move_noise_covariance(q)), QTLSim.fluct_direction(q))

# --- a run driven by one lowered life cycle -------------------------------

function run_composed(w::World, rng; record::Bool=true, cycle = lifecycle().generation,
                      images = Dict{Symbol,Any}())
    zygotes = LifeCycle[:Zygotes]
    dom(cycle) == zygotes && codom(cycle) == zygotes ||
        throw(ArgumentError("a run iterates a life cycle Zygotes → Zygotes; got $(dom(cycle)) → $(codom(cycle))"))
    rec = records(; traits = w.p.NoTraits)
    step = lower(cycle, :kernel; images, record = record ? rec : nothing, model = model_of(w.p))
    p = w.p
    init_env!(w, rng)
    # whole cycles until NoGenerations generations are done: a cycle of several generations
    # advances w.gen at each swap, and one that does not advance it counts as one
    w.gen = 0
    while w.gen < p.NoGenerations
        it = w.gen
        step(w, rng)
        w.extinct && return (extinct=true, extinct_gen=w.gen, rec...)
        w.gen = max(w.gen, it + 1)
    end
    (extinct=false, extinct_gen=-1, rec...)
end

# --- stratification: ⊗ builds a metapopulation of demes ------------------

"Stratification into `k` demes, used as `p ⊗ demes(k)`."
struct Demes
    k::Int
end

demes(k::Integer) = Demes(k)

"A metapopulation: demes (`World`s) that operators act on one at a time."
mutable struct Stratified
    demes::Vector{World}
end

# The demes share one pool of blocks (see `World`), so that a migrant moves by swapping
# block ids; demes with different numbers of loci cannot share one.
⊗(p::Params, d::Demes) = (pool = QTLSim.BlockPool(p.NoQTL; traits = p.NoTraits); Stratified([World(p; pool) for _ in 1:d.k]))
function ⊗(w::World, d::Demes)
    w1 = deepcopy(w)
    w1.anc = nothing                 # the copies record no pedigree: record_pedigree!(w ⊗ demes(k))
    Stratified([w1; [QTLSim.shared_copy(w1) for _ in 2:d.k]])
end
function ⊗(ps::AbstractVector{Params}, d::Demes)
    length(ps) == d.k || throw(ArgumentError("$(length(ps)) parameter sets for $(d.k) demes"))
    all(q -> q.NoGenerations == ps[1].NoGenerations, ps) ||
        throw(ArgumentError("the demes of a metapopulation run for the same number of generations"))
    pools = Dict{Tuple{Int,Int},QTLSim.BlockPool}()
    Stratified([World(q; pool = get!(() -> QTLSim.BlockPool(q.NoQTL; traits = q.NoTraits), pools, (q.NoQTL, q.NoTraits)))
                for q in ps])
end

isextinct(s::Stratified) = all(w -> w.extinct, s.demes)

# the pedigree of a metapopulation, and identity by descent among its demes (src/ancestry.jl)
QTLSim.record_pedigree!(s::Stratified; kw...) = QTLSim.record_pedigree!(s.demes; kw...)
QTLSim.tag_founders!(s::Stratified) = QTLSim.tag_founders!(s.demes)
QTLSim.ibd_sharing(s::Stratified) = QTLSim.ibd_sharing(s.demes)
QTLSim.ibd_fst(s::Stratified) = QTLSim.ibd_fst(s.demes)
totalpop(s::Stratified) = sum(w -> w.N, s.demes; init=0)
ndemes(s::Stratified) = length(s.demes)

# a :meta operator acts on the whole metapopulation, a :deme operator on each deme
function (op::Operator)(s::Stratified, rng)
    if op.scope === :meta
        op.act(s, rng)
    else
        for w in s.demes
            w.extinct || op.act(w, rng)
        end
    end
    s
end

function migrate!(s::Stratified, rng; rate=0.01, demes=eachindex(s.demes))
    k = length(demes)
    k < 2 && return (count=0,)
    shape = QTLSim.genome_shape(s.demes[first(demes)])
    all(d -> QTLSim.genome_shape(s.demes[d]) == shape, demes) ||
        throw(ArgumentError("migration needs demes with the same loci and ploidy; got $([QTLSim.genome_shape(s.demes[d]) for d in demes])"))
    rate isa AbstractMatrix && return migrate_graph!(s, rng, QTLSim.check_migration(rate, k), demes)
    count = 0
    for di in randperm(rng, k)      # the demes in a uniformly random order, see the docstring
        w = s.demes[demes[di]]
        w.extinct && continue
        for i in 1:w.N
            rand(rng) < rate || continue
            dj = rand(rng, 1:k-1); dj >= di && (dj += 1)   # a different deme
            v = s.demes[demes[dj]]; v.extinct && continue
            count += swap_individuals!(w, i, v, rand(rng, 1:v.N))
        end
    end
    (count=count,)
end

# Exchange individual i of deme w and individual j of deme v: their genomes, by swapping
# block ids when the demes share a pool, as demes made by ⊗ do (QTLSim.exchange_genomes!),
# and their environments (every trait's, when there are several).
function swap_individuals!(w::World, i, v::World, j)
    QTLSim.exchange_genomes!(w, i, v, j)
    QTLSim.exchange_env!(w, i, v, j)
    1
end

# Migration along a graph: an individual of the deme in position `a` starts a swap with
# probability `M[a, b]` with a random individual of the deme in position `b`. One uniform
# draw decides both whether and where: `u < Σ_b M[a, b]`, and `b` is where the cumulative
# row first exceeds `u`. As in `migrate!`, the demes take their turns in a uniformly random
# order and the individuals of a deme in order.
function migrate_graph!(s::Stratified, rng, M::AbstractMatrix, demes)
    count = 0
    for a in randperm(rng, length(demes))
        w = s.demes[demes[a]]
        w.extinct && continue
        cum = cumsum(M[a, :])
        for i in 1:w.N
            u = rand(rng)
            u < cum[end] || continue
            b = findfirst(>(u), cum)
            v = s.demes[demes[b]]; v.extinct && continue
            count += swap_individuals!(w, i, v, rand(rng, 1:v.N))
        end
    end
    (count=count,)
end

migration(rate) = Operator(:migrate, (s, rng) -> (migrate!(s, rng; rate=rate); s), :meta)

# --- the image of a life cycle on several demes ---------------------------

"""
    DemeOp(n, perm, act)

The kernel image of a life-cycle morphism between tensor products of `n` stages, one
stage per deme: `act(s, demes, rng)` acts on the demes `s.demes[demes]` of a
metapopulation `s`, and output position `i` holds the deme at input position `perm[i]`
(a braiding exchanges demes; the others leave them in place). `⊗` places operators side
by side on consecutive demes, and `⋅` composes them.
"""
struct DemeOp{F}
    n::Int
    perm::Vector{Int}
    act::F
end

Base.show(io::IO, op::DemeOp) = print(io, "DemeOp(", op.n, " deme", op.n == 1 ? "" : "s", ")")

"An operator on one population as an operator on one deme, which it skips once extinct."
demeop(op::Operator) =
    DemeOp(1, [1], (s, d, rng) -> (w = s.demes[d[1]]; w.extinct || op(w, rng); s))
demeop(op::DemeOp) = op

const AnyOp = Union{Operator,DemeOp}

function compose(a::DemeOp, b::DemeOp)
    a.n == b.n || throw(ArgumentError("composing operators on $(a.n) and $(b.n) demes"))
    DemeOp(a.n, a.perm[b.perm], (s, d, rng) -> (a.act(s, d, rng); b.act(s, d[a.perm], rng); s))
end
compose(a::AnyOp, b::AnyOp) = compose(demeop(a), demeop(b))
compose(a::AnyOp, b::AnyOp, c::AnyOp, rest::AnyOp...) = foldl(compose, (a, b, c, rest...))
⋅(a::DemeOp, b::AnyOp) = compose(a, b)
⋅(a::Operator, b::DemeOp) = compose(a, b)

function otimes(a::DemeOp, b::DemeOp)
    na = a.n
    DemeOp(na + b.n, vcat(a.perm, na .+ b.perm),
           (s, d, rng) -> (a.act(s, d[1:na], rng); b.act(s, d[na+1:end], rng); s))
end
otimes(a::AnyOp, b::AnyOp) = otimes(demeop(a), demeop(b))
otimes(a::AnyOp, b::AnyOp, c::AnyOp, rest::AnyOp...) = foldl(otimes, (a, b, c, rest...))

"Migration among `k` demes at `rate`, as the kernel image of the step `migrate`."
migrate_op(k, rate) = DemeOp(k, collect(1:k), (s, d, rng) -> (migrate!(s, rng; rate, demes = d); s))

function run_composed(s::Stratified, rng; cycle, record::Bool=true, images = Dict{Symbol,Any}(),
                      migration_rate::Union{Real,AbstractMatrix} = 0.0)
    k = ndemes(s)
    zk = zygotes_word(k)
    dom(cycle) == zk && codom(cycle) == zk ||
        throw(ArgumentError("a run of $k demes iterates a life cycle on $k zygote stages, Zygotes⊗…⊗Zygotes; " *
                            "got $(dom(cycle)) → $(codom(cycle))"))
    m = model_of(first(s.demes).p)
    all(w -> model_of(w.p) == m, s.demes) ||
        throw(ArgumentError("the demes of a metapopulation belong to one model; got $(unique(model_of(w.p) for w in s.demes))"))
    rec = records(demes = s.demes, traits = first(s.demes).p.NoTraits)
    step = demeop(lower(cycle, :kernel; images, record = record ? rec : nothing, migration_rate, model = m))
    foreach(w -> init_env!(w, rng), s.demes)
    p = first(s.demes).p
    d = collect(1:k)
    extinct_gen = fill(-1, k)
    gen = 0                      # whole cycles, as for one World, with the demes in step
    while gen < p.NoGenerations
        foreach(w -> w.extinct || (w.gen = gen), s.demes)
        step.act(s, d, rng)
        for (i, w) in enumerate(s.demes)
            w.extinct && extinct_gen[i] < 0 && (extinct_gen[i] = w.gen)
        end
        isextinct(s) && break
        d = d[step.perm]
        gen = max(gen + 1, maximum(w.gen for w in s.demes if !w.extinct))
    end
    out = (extinct = isextinct(s), extinct_gen, gen = rec.gen, deme = rec.deme, pheno_mean = rec.pheno_mean,
           pheno_var = rec.pheno_var, geno_var = rec.geno_var, nsurv = rec.nsurv, theta = rec.theta, omega = rec.omega,
           mean_fitness = rec.mean_fitness, sel_diff = rec.sel_diff, popsize = rec.popsize)
    haskey(rec, :gmatrix) ? merge(out, NamedTuple(k => rec[k] for k in keys(QTLSim.trait_records()))) : out
end
