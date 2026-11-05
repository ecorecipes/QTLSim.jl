# Syntactic sugar for the compositional model: unicode operators for
# diagrammatic composition (⋅) and stratification (⊗).
#
# The evolutionary operators of `model.jl` are morphisms on world state. Here we
# wrap them as first-class `Operator` values so a generation can be written as a
# left-to-right composite
#
#     generation = develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap
#
# `⋅` is DIAGRAMMATIC (left-to-right) composition: `(f ⋅ g)(w) = g(f(w))`. This
# reads in the biological order and is the reverse of the mathematical `∘` used
# in lean/CategoricalQuantGen/Composition.lean (`generation = mutate ∘ reproduce ∘
# select`). Stratification `⊗` builds a metapopulation of demes over which the
# same operators act pointwise, with a `migrate` operator coupling them.

# Extend Catlab's own categorical operators: `⋅` (\cdot, diagrammatic
# composition) and `⊗` (\otimes, monoidal product). Using Catlab's symbols makes
# our operators literally the AlgebraicJulia composition/stratification syntax.
import Catlab: ⋅, ⊗

# --- Operator: a named RNG-threaded morphism on world state --------------

"""
    Operator(name, act)

A first-class evolutionary operator wrapping an in-place kernel step. `act` has
signature `act(state, rng) -> state`. Operators are called like functions and
composed with [`⋅`](@ref).
"""
struct Operator{F}
    name::Symbol
    act::F
    scope::Symbol   # :deme — acts on a `World` (applied pointwise over `⊗`);
                    # :meta — acts on the whole `Stratified` metapopulation (e.g. migration).
end

# Default scope is :deme; most operators are per-individual/per-deme kernel steps.
Operator(name::Symbol, act) = Operator(name, act, :deme)

(op::Operator)(state, rng) = op.act(state, rng)

Base.show(io::IO, op::Operator) = print(io, "Operator(:", op.name, ")")

"""
    f ⋅ g

Diagrammatic (left-to-right) composition: apply `f` then `g`, i.e.
`(f ⋅ g)(w, rng) = g(f(w, rng), rng)`. Short-circuits once a `World` is extinct
so downstream operators are skipped (mirrors the `w.extinct` guards in `run!`).
"""
function ⋅(f::Operator, g::Operator)
    # A composite is :meta if any factor is — so migration composed into a
    # generation makes the whole generation act on the metapopulation, with each
    # factor re-dispatching (deme-scoped factors act pointwise, :meta ones whole).
    sc = (f.scope === :meta || g.scope === :meta) ? :meta : :deme
    Operator(Symbol(f.name, :⋅, g.name), function (state, rng)
        f(state, rng)
        isextinct(state) && return state
        g(state, rng)
    end, sc)
end

# extinction predicate — extended for stratified state below
isextinct(w::World) = w.extinct

# --- base operators wrapping the validated kernel ------------------------

"Genotype/phenotype development (`develop!`)."
const develop = Operator(:develop, (w, rng) -> develop!(w))

"Stabilizing selection + survivor sampling (`select!`)."
const select = Operator(:select, (w, rng) -> select!(w, rng))

"Random-mating reproduction into the offspring buffer (`reproduce!`)."
const reproduce = Operator(:reproduce, (w, rng) -> isextinct(w) ? w : reproduce!(w, rng))

"Poisson mutation of the offspring buffer (`mutate!`)."
const mutate = Operator(:mutate, (w, rng) -> mutate!(w, rng))

"Buffer swap advancing offspring to parents, incrementing the generation (`swapbuffers!`)."
const swap = Operator(:swap, (w, rng) -> (swapbuffers!(w); w.gen += 1; w))

"""
    generation

One generation as the composite `develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap`.
Equivalent (up to bookkeeping) to `generation!`; the categorical composite proved
associative and iterable in lean/CategoricalQuantGen/Composition.lean.
"""
const generation = develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap

# --- run driven by the composed operator ---------------------------------

"""
    run_composed(w::World, rng; record=true)

Drive a full run through the composed `generation` operator, reproducing `run!`
exactly (same RNG draw order, recording, and fitness-width schedule). Returns the
same statistics NamedTuple as `run!`. For identical seeds `run_composed` and
`run!` are bit-for-bit equal.
"""
function run_composed(w::World, rng; record::Bool=true)
    p = w.p
    esd = envsd(p)
    @inbounds for i in 1:p.PopSize
        w.env[i] = esd * randn(rng)
    end
    gen = Int[]; pmean = Float64[]; pvar = Float64[]
    gvar = Float64[]; nsurv = Int[]; θs = Float64[]; ωs = Float64[]
    for it in 0:p.NoGenerations-1
        w.gen = it
        develop(w, rng)
        select(w, rng)
        if w.extinct
            return (extinct=true, extinct_gen=it, gen=gen, pheno_mean=pmean,
                    pheno_var=pvar, geno_var=gvar, nsurv=nsurv, theta=θs, omega=ωs)
        end
        reproduce(w, rng)
        if record && (it % p.Sample == 0)
            push!(gen, it); push!(pmean, mean(w.phenotype)); push!(pvar, var(w.phenotype))
            push!(gvar, var(w.genotype)); push!(nsurv, w.nsurv)
            push!(θs, w.θ); push!(ωs, w.ω)
        end
        it == p.ChangeFitnessGen && (w.ω = p.ChangeFitnessSD)
        mutate(w, rng)
        swap(w, rng)
    end
    (extinct=false, extinct_gen=-1, gen=gen, pheno_mean=pmean,
     pheno_var=pvar, geno_var=gvar, nsurv=nsurv, theta=θs, omega=ωs)
end

# --- stratification: ⊗ builds a metapopulation of demes ------------------

"""
    Demes(k)

Stratification specifier: replicate the base population into `k` demes. Combine
with a base `World`/`Params` via [`⊗`](@ref).
"""
struct Demes
    k::Int
end

"""
    Stratified

A metapopulation: a vector of sub-`World`s (demes) sharing the same `Params`.
Evolutionary operators act on it pointwise (one deme at a time); a `migrate`
operator couples the demes. This is a discrete stratification — the categorical
`⊗` of `k` copies of the base dynamical system plus a migration morphism.
"""
mutable struct Stratified
    demes::Vector{World}
end

"""
    base ⊗ Demes(k)

Stratify a base population into `k` demes (`⊗` = `\\otimes`). `base` may be a
`Params` (each deme an independent fresh `World`) or a `World` (deep-copied into
each deme). Returns a [`Stratified`](@ref) metapopulation.
"""
⊗(p::Params, d::Demes) = Stratified([World(p) for _ in 1:d.k])
⊗(w::World, d::Demes) = Stratified([deepcopy(w) for _ in 1:d.k])

isextinct(s::Stratified) = all(w -> w.extinct, s.demes)

"Total number of individuals across all demes."
totalpop(s::Stratified) = sum(w -> w.p.PopSize, s.demes; init=0)
ndemes(s::Stratified) = length(s.demes)

# An Operator applied to a Stratified state: a :meta operator (e.g. migration, or
# any composite containing it) acts on the whole metapopulation; a :deme operator
# acts on each deme in turn. This is what lets a coupling operator be composed
# into a generation with `⋅` (e.g. `generation ⋅ migration(0.01)`).
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

"""
    migrate(s::Stratified, rng; rate=0.01)

Migration operator: each individual moves to a uniformly-chosen other deme with
probability `rate`, swapping genome columns and environments between demes.
Population sizes per deme are preserved by pairing each emigrant with a random
immigrant slot (a symmetric exchange), so `totalpop` is conserved. Returns a
`(count=...,)` NamedTuple of the number of exchanges performed.
"""
function migrate!(s::Stratified, rng; rate=0.01)
    k = length(s.demes)
    k < 2 && return (count=0,)
    count = 0
    for (di, w) in enumerate(s.demes)
        w.extinct && continue
        for i in 1:w.p.PopSize
            rand(rng) < rate || continue
            dj = rand(rng, 1:k-1); dj >= di && (dj += 1)   # a different deme
            v = s.demes[dj]; v.extinct && continue
            j = rand(rng, 1:v.p.PopSize)
            # symmetric exchange of individual i (deme di) ↔ individual j (deme dj)
            for c in axes(w.genome, 2), l in axes(w.genome, 1)
                w.genome[l, c, i], v.genome[l, c, j] = v.genome[l, c, j], w.genome[l, c, i]
            end
            w.env[i], v.env[j] = v.env[j], w.env[i]
            count += 1
        end
    end
    (count=count,)
end

"""
    migration(rate)

An [`Operator`](@ref) form of [`migrate!`](@ref) with a fixed `rate`, so it can
be composed with `⋅` into a stratified generation, e.g.
`generation ⋅ migration(0.01)`.
"""
migration(rate) = Operator(:migrate, (s, rng) -> (migrate!(s, rng; rate=rate); s), :meta)
