# The compositional evolutionary model.
#
# Each biological process is an operator (a Kleisli morphism of the RNG monad):
# it maps the world state to a new world state, using the RNG. One generation is
# their composite
#
#     generation = mutate ∘ reproduce ∘ select ∘ develop
#
# and a run is the iterate of one generation — exactly the structure proved in
# lean/CategoricalQuantGen/Composition.lean (`generation`, `loop_eq_iterate`). The
# numerics replicate QTSim (Repository/Code/QTSim/src/QTSim.c) so that results
# are statistically consistent with the C reference.

"""
    World

Mutable population state carried through the operators. Holds the current
parents' genomes/environments, per-individual scratch (genotype/phenotype/
fitness), the survivor index set, the current selection optimum `θ` and width
`ω`, and a double buffer (`sgenome`/`senv`) for the offspring generation.

Genomes are stored as `(NoQTL, Ploidy, PopSize)`; the genotype of individual `i`
is the sum over all its loci and chromosomes (QTSim `CaclGenotype`).
"""
mutable struct World
    p::Params
    gen::Int
    genome::Array{Float64,3}
    env::Vector{Float64}
    genotype::Vector{Float64}
    phenotype::Vector{Float64}
    fitness::Vector{Float64}
    survivors::Vector{Int}
    nsurv::Int
    θ::Float64
    ω::Float64
    sgenome::Array{Float64,3}
    senv::Vector{Float64}
    extinct::Bool
end

function World(p::Params)
    p.CrossoverModel in (:qtsim, :all) ||
        throw(ArgumentError("CrossoverModel must be :qtsim or :all, got :$(p.CrossoverModel)"))
    N, L, P = p.PopSize, p.NoQTL, p.Ploidy
    World(p, 0,
          zeros(L, P, N), zeros(N),
          zeros(N), zeros(N), zeros(N),
          Int[], 0, 0.0, initFitSD(p),
          zeros(L, P, N), zeros(N), false)
end

# --- environment ---------------------------------------------------------

@inline envsd(p::Params) = sqrt(max(envVar(p), 0.0))

# --- operators (morphisms K → K) ----------------------------------------

"""
    develop!(w)

Genotype = Σ allele effects; phenotype = genotype + environment (QTSim
`CaclPopPhenotype`). Deterministic given the environments drawn at birth.
"""
function develop!(w::World)
    @inbounds for i in 1:w.p.PopSize
        w.genotype[i] = sum(@view w.genome[:, :, i])
        w.phenotype[i] = w.genotype[i] + w.env[i]
    end
    w
end

"Move the optimum `θ` if the schedule says so (QTSim `SetFitnessMean`)."
function moveoptimum!(w::World)
    p = w.p
    (p.MoveFitnessSDGen != -1 && w.gen >= p.MoveFitnessSDGen) || return w
    m = mean(w.phenotype)
    if p.MoveType == :phenotype
        w.θ = m + sqrt(p.PhenotypeVar) * p.MoveFitnessSD
    elseif p.MoveType == :fitness
        w.θ = m + w.ω * p.MoveFitnessSD
    else # :continuous
        w.θ += p.MoveFitnessSD
    end
    w
end

"""
    select!(w, rng)

Move the optimum, assign each individual a stabilizing-selection fitness, and
draw the survivor set (survive iff `rand() < fitness`). Mirrors QTSim `Fitness`
+ `SetSurviving`. Sets `w.extinct` if too few survive to reproduce.

Fitness is the normalized Gaussian `pdf(Normal(0,ω), (z−θ)/sd) / (0.4/ω)`,
replicating QTSim exactly (the `0.4 ≈ 1/√(2π)` constant included); the idealized
form `exp(−(z−θ)²/2ω²)` is what Lean `StabilizingSelection.lean` proves about.
"""
function select!(w::World, rng)
    p = w.p
    moveoptimum!(w)
    sd = p.StandardisePhenotype ? std(w.phenotype) : 1.0
    @inbounds for i in 1:p.PopSize
        if p.Drift
            w.fitness[i] = 1.0
        else
            x = (w.phenotype[i] - w.θ) / sd
            w.fitness[i] = pdf(Normal(0.0, w.ω), x) / (0.4 / w.ω)
        end
    end
    empty!(w.survivors)
    @inbounds for i in 1:p.PopSize
        rand(rng) < w.fitness[i] && push!(w.survivors, i)
    end
    w.nsurv = length(w.survivors)
    if (p.Ploidy == 1 && w.nsurv == 0) || (p.Ploidy == 2 && w.nsurv < 2)
        w.extinct = true
    end
    w
end

"""
    gamete!(dest, chrA, chrB, rng, recombMean, L, allpoints=false)

Build one gamete from a parent's two chromosomes (QTSim `CreateCombMap` +
`RecombinationPoissonChr`). Draw `Poisson(recombMean)` distinct crossover points
in `0:L-2` (none if the count is 0 or ≥ `L`), start on a random parental
chromosome, and switch strands after crossover points. QTSim never advances its
point index, so only the LOWEST point switches: a gamete has at most one
crossover. That is the default; `allpoints=true` switches at every point
(Params `CrossoverModel = :all`). With `L == 1` both reduce to picking one
parental chromosome at random (Mendelian segregation).
"""
function gamete!(dest, chrA, chrB, rng, recombMean, L, allpoints::Bool=false)
    npts = rand(rng, Poisson(recombMean))
    if npts == 0 || npts >= L
        src = rand(rng, Bool) ? chrA : chrB
        dest .= src
        return dest
    end
    pts = Int[]
    while length(pts) < npts
        c = rand(rng, 0:L-2)
        c in pts || push!(pts, c)
    end
    sort!(pts)
    allpoints || resize!(pts, 1)       # QTSim: only the lowest point switches
    cp = rand(rng, Bool) ? chrA : chrB
    other = cp === chrA ? chrB : chrA
    pi = 1
    @inbounds for idx in 0:L-1
        dest[idx+1] = cp[idx+1]
        if pi <= length(pts) && idx == pts[pi]
            cp, other = other, cp
            pi += 1
        end
    end
    dest
end

"""
    reproduce!(w, rng)

Fill the offspring buffer `sgenome`/`senv` by fitness-independent random mating
among survivors (QTSim `NewIn`). Haploid: clone a random survivor. Diploid: two
distinct survivors each contribute a (recombined) gamete. Offspring get a fresh
environment.
"""
function reproduce!(w::World, rng)
    p = w.p
    L = p.NoQTL
    esd = envsd(p)
    allpts = p.CrossoverModel === :all
    if p.Ploidy == 1
        @inbounds for j in 1:p.PopSize
            par = w.survivors[rand(rng, 1:w.nsurv)]
            @views w.sgenome[:, 1, j] .= w.genome[:, 1, par]
            w.senv[j] = esd * randn(rng)
        end
    else
        @inbounds for j in 1:p.PopSize
            p1 = w.survivors[rand(rng, 1:w.nsurv)]
            p2 = w.survivors[rand(rng, 1:w.nsurv)]
            while p2 == p1
                p2 = w.survivors[rand(rng, 1:w.nsurv)]
            end
            gamete!(@view(w.sgenome[:, 1, j]), @view(w.genome[:, 1, p1]),
                    @view(w.genome[:, 2, p1]), rng, p.RecombPoissonMean, L, allpts)
            gamete!(@view(w.sgenome[:, 2, j]), @view(w.genome[:, 1, p2]),
                    @view(w.genome[:, 2, p2]), rng, p.RecombPoissonMean, L, allpts)
            w.senv[j] = esd * randn(rng)
        end
    end
    w
end

"""
    mutate!(w, rng)

Add Gaussian mutations to the offspring buffer (QTSim `MutateIndividualChrNo`).
Each chromosome receives `Poisson(L·rate)` mutations, each of variance
`m/(ploidy·L·rate)·V_E`; in expectation this injects `m·V_E` of variance per
individual per generation (Lean `mutationalVarianceInput_eq`).
"""
function mutate!(w::World, rng)
    p = w.p
    L = p.NoQTL
    rate = p.MutationRatePerGamete
    mutsd = sqrt(max((p.MutationScalar / (p.Ploidy * L * rate)) * envVar(p), 0.0))
    λ = L * rate
    @inbounds for j in 1:p.PopSize
        for c in 1:p.Ploidy
            no = rand(rng, Poisson(λ))
            for _ in 1:no
                pos = rand(rng, 1:L)
                w.sgenome[pos, c, j] += mutsd * randn(rng)
            end
        end
    end
    w
end

"Swap the parent and offspring buffers (QTSim double-buffer `swap`)."
function swapbuffers!(w::World)
    w.genome, w.sgenome = w.sgenome, w.genome
    w.env, w.senv = w.senv, w.env
    w
end

"""
    generation!(w, rng)

One composed generation `mutate ∘ reproduce ∘ select ∘ develop`, advancing the
offspring into the current buffer. Provided as the standalone categorical
composite (used in tests); `run!` inlines the same sequence plus recording and
the fitness-width schedule. No-op past extinction.
"""
function generation!(w::World, rng)
    w.extinct && return w
    develop!(w)
    select!(w, rng)
    w.extinct && return w
    reproduce!(w, rng)
    mutate!(w, rng)
    swapbuffers!(w)
    w.gen += 1
    w
end

# --- driver --------------------------------------------------------------

"""
    run!(w, rng; record=true)

Iterate the generation operator for `NoGenerations`, replicating QTSim's
`RunSim` order: develop → select → reproduce → (record) → fitness-width change →
mutate → swap. Returns a NamedTuple of per-sampled-generation statistics
(`gen`, `pheno_mean`, `pheno_var`, `geno_var`, `nsurv`, `theta`, `omega`) plus
`extinct`/`extinct_gen`.
"""
function run!(w::World, rng; record::Bool=true)
    p = w.p
    esd = envsd(p)
    @inbounds for i in 1:p.PopSize
        w.env[i] = esd * randn(rng)
    end
    gen = Int[]; pmean = Float64[]; pvar = Float64[]
    gvar = Float64[]; nsurv = Int[]; θs = Float64[]; ωs = Float64[]
    for it in 0:p.NoGenerations-1
        w.gen = it
        develop!(w)
        select!(w, rng)
        if w.extinct
            return (extinct=true, extinct_gen=it, gen=gen, pheno_mean=pmean,
                    pheno_var=pvar, geno_var=gvar, nsurv=nsurv, theta=θs, omega=ωs)
        end
        reproduce!(w, rng)
        if record && (it % p.Sample == 0)
            push!(gen, it)
            push!(pmean, mean(w.phenotype))
            push!(pvar, var(w.phenotype))
            push!(gvar, var(w.genotype))
            push!(nsurv, w.nsurv)
            push!(θs, w.θ)
            push!(ωs, w.ω)
        end
        it == p.ChangeFitnessGen && (w.ω = p.ChangeFitnessSD)
        mutate!(w, rng)
        swapbuffers!(w)
    end
    (extinct=false, extinct_gen=-1, gen=gen, pheno_mean=pmean,
     pheno_var=pvar, geno_var=gvar, nsurv=nsurv, theta=θs, omega=ωs)
end

"""
    simulate(p::Params; seed=p.Seed)

Run one simulation from `Params`, returning the recorded statistics (see `run!`).
"""
function simulate(p::Params; seed=p.Seed)
    rng = seed === nothing ? MersenneTwister() : MersenneTwister(seed)
    run!(World(p), rng)
end
