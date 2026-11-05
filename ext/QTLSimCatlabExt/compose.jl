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

"""
Empty columns for the statistics `record` collects, as `run!` returns them; with
`demes`, also the deme of each row, looked up in `demeof`.
"""
records(; demes = nothing) = demes === nothing ?
    (gen = Int[], pheno_mean = Float64[], pheno_var = Float64[], geno_var = Float64[],
     nsurv = Int[], theta = Float64[], omega = Float64[]) :
    (gen = Int[], deme = Int[], pheno_mean = Float64[], pheno_var = Float64[],
     geno_var = Float64[], nsurv = Int[], theta = Float64[], omega = Float64[],
     demeof = IdDict{World,Int}(w => i for (i, w) in enumerate(demes)))

"Append the parent generation's statistics to `rec` every `Sample` generations (QTSim `Output`)."
function record!(rec, w::World)
    (rec === nothing || w.gen % w.p.Sample != 0) && return w
    push!(rec.gen, w.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[w])
    push!(rec.pheno_mean, mean(w.phenotype))
    push!(rec.pheno_var, var(w.phenotype)); push!(rec.geno_var, var(w.genotype))
    push!(rec.nsurv, w.nsurv); push!(rec.theta, w.θ); push!(rec.omega, w.ω)
    w
end

"Change the fitness width at `ChangeFitnessGen`, effective from the next generation."
function rescale!(w::World)
    w.gen == w.p.ChangeFitnessGen && (w.ω = w.p.ChangeFitnessSD)
    w
end

"""
The kernel image of each step of the life cycle, by name. `record` appends to `rec`, or
discards its observations when `rec` is `nothing`.
"""
kernel_images(rec = nothing) = Dict{Symbol,Any}(
    :develop => Operator(:develop, (w, rng) -> QTLSim.develop!(w)),
    :select => Operator(:select, (w, rng) -> QTLSim.select!(w, rng)),
    :reproduce => Operator(:reproduce, (w, rng) -> isextinct(w) ? w : QTLSim.reproduce!(w, rng)),
    :record => Operator(:record, (w, rng) -> record!(rec, w)),
    :rescale => Operator(:rescale, (w, rng) -> rescale!(w)),
    :mutate => Operator(:mutate, (w, rng) -> QTLSim.mutate!(w, rng)),
    :swap => Operator(:swap, (w, rng) -> (QTLSim.swapbuffers!(w); w.gen += 1; w)))

# --- a run driven by one lowered life cycle -------------------------------

function run_composed(w::World, rng; record::Bool=true, cycle = lifecycle().generation,
                      images = Dict{Symbol,Any}())
    zygotes = LifeCycle[:Zygotes]
    dom(cycle) == zygotes && codom(cycle) == zygotes ||
        throw(ArgumentError("a run iterates a life cycle Zygotes → Zygotes; got $(dom(cycle)) → $(codom(cycle))"))
    rec = records()
    step = lower(cycle, :kernel; images, record = record ? rec : nothing)
    p = w.p
    esd = envsd(p)
    @inbounds for i in 1:p.PopSize
        w.env[i] = esd * randn(rng)
    end
    for it in 0:p.NoGenerations-1
        w.gen = it
        step(w, rng)
        w.extinct && return (extinct=true, extinct_gen=it, rec...)
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

⊗(p::Params, d::Demes) = Stratified([World(p) for _ in 1:d.k])
⊗(w::World, d::Demes) = Stratified([deepcopy(w) for _ in 1:d.k])

isextinct(s::Stratified) = all(w -> w.extinct, s.demes)
totalpop(s::Stratified) = sum(w -> w.p.PopSize, s.demes; init=0)
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
    count = 0
    for (di, d) in enumerate(demes)
        w = s.demes[d]
        w.extinct && continue
        for i in 1:w.p.PopSize
            rand(rng) < rate || continue
            dj = rand(rng, 1:k-1); dj >= di && (dj += 1)   # a different deme
            v = s.demes[demes[dj]]; v.extinct && continue
            j = rand(rng, 1:v.p.PopSize)
            for c in axes(w.genome, 2), l in axes(w.genome, 1)
                w.genome[l, c, i], v.genome[l, c, j] = v.genome[l, c, j], w.genome[l, c, i]
            end
            w.env[i], v.env[j] = v.env[j], w.env[i]
            count += 1
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
                      migration_rate::Real = 0.0)
    k = ndemes(s)
    zk = zygotes_word(k)
    dom(cycle) == zk && codom(cycle) == zk ||
        throw(ArgumentError("a run of $k demes iterates a life cycle on $k zygote stages, Zygotes⊗…⊗Zygotes; " *
                            "got $(dom(cycle)) → $(codom(cycle))"))
    rec = records(demes = s.demes)
    step = demeop(lower(cycle, :kernel; images, record = record ? rec : nothing, migration_rate))
    for w in s.demes
        esd = envsd(w.p)
        @inbounds for i in 1:w.p.PopSize
            w.env[i] = esd * randn(rng)
        end
    end
    p = first(s.demes).p
    d = collect(1:k)
    extinct_gen = fill(-1, k)
    for it in 0:p.NoGenerations-1
        foreach(w -> w.extinct || (w.gen = it), s.demes)
        step.act(s, d, rng)
        for (i, w) in enumerate(s.demes)
            w.extinct && extinct_gen[i] < 0 && (extinct_gen[i] = it)
        end
        isextinct(s) && break
        d = d[step.perm]
    end
    (extinct = isextinct(s), extinct_gen, gen = rec.gen, deme = rec.deme, pheno_mean = rec.pheno_mean,
     pheno_var = rec.pheno_var, geno_var = rec.geno_var, nsurv = rec.nsurv, theta = rec.theta, omega = rec.omega)
end
