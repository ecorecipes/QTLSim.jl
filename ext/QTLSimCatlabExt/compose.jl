# Evolutionary processes as first-class operators. A generation is the left-to-right
# composite
#
#     generation = develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap
#
# `⋅` is Catlab's diagrammatic composition, `(f ⋅ g)(w) = g(f(w))`: it reads in the
# biological order and is the reverse of the `∘` in Composition.lean
# (`generation = mutate ∘ reproduce ∘ select`). Catlab's monoidal product `⊗` builds a
# metapopulation of demes over which the same operators act pointwise, with a migration
# operator coupling them.

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

isextinct(w::World) = w.extinct

# --- the steps of QTSim's generation loop, wrapping the kernel ----------------

const develop = Operator(:develop, (w, rng) -> QTLSim.develop!(w))
const select = Operator(:select, (w, rng) -> QTLSim.select!(w, rng))
const reproduce = Operator(:reproduce, (w, rng) -> isextinct(w) ? w : QTLSim.reproduce!(w, rng))
const mutate = Operator(:mutate, (w, rng) -> QTLSim.mutate!(w, rng))
const swap = Operator(:swap, (w, rng) -> (QTLSim.swapbuffers!(w); w.gen += 1; w))
const generation = develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap

operators() = (; develop, select, reproduce, mutate, swap, generation)

# --- a run driven by the composed operator -------------------------------

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
