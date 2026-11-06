# The compositional evolutionary model.
#
# Each biological process is an operator (a Kleisli morphism of the RNG monad):
# it maps the world state to a new world state, using the RNG. One generation is
# their composite, in QTSim's `RunSim` order (diagrammatic, first step first)
#
#     develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅ mutate ⋅ swap
#
# and a run is the iterate of one generation — exactly the structure proved in
# lean/CategoricalQuantGen/Composition.lean (`generation`, `loop_eq_iterate`). The
# categorical layer (ext/QTLSimCatlabExt) writes the same composite as a morphism of a
# free category and lowers it to these operators. The numerics replicate QTSim
# (reference/QTSim/src/QTSim.c) so that results are statistically consistent with the
# C reference.

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
    checkparams(p)
    N, L, P = p.PopSize, p.NoQTL, p.Ploidy
    World(p, 0,
          zeros(L, P, N), zeros(N),
          zeros(N), zeros(N), zeros(N),
          Int[], 0, 0.0, initFitSD(p),
          zeros(L, P, N), zeros(N), false)
end

# --- environment ---------------------------------------------------------

@inline envsd(p::Params) = sqrt(max(envVar(p), 0.0))

"""
    init_env!(w, rng)

Draw the founders' environmental deviations, `N(0, V_E)` each (QTSim `CreatePop`, via
`SetIndividualEnv`). `run!` and `run_composed` do this first; a `World` made with
`World(p)` starts with all deviations zero, so call it before driving a fresh `World`
with [`generation!`](@ref) or the operators directly.
"""
function init_env!(w::World, rng)
    esd = envsd(w.p)
    @inbounds for i in 1:w.p.PopSize
        w.env[i] = esd * randn(rng)
    end
    w
end

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
    gamete!(dest, chrA, chrB, rng, recombMean, L, allpoints=false, pts=Int[])

Build one gamete from a parent's two chromosomes (QTSim `CreateCombMap` +
`RecombinationPoissonChr`). Draw `Poisson(recombMean)` distinct crossover points
in `0:L-2` (none if the count is 0 or ≥ `L`), start on a random parental
chromosome, and switch strands after crossover points. QTSim never advances its
point index, so only the LOWEST point switches: a gamete has at most one
crossover. That is the default; `allpoints=true` switches at every point
(Params `CrossoverModel = :all`). With `L == 1` both reduce to picking one
parental chromosome at random (Mendelian segregation).
"""
function gamete!(dest, chrA, chrB, rng, recombMean, L, allpoints::Bool=false, pts::Vector{Int}=Int[])
    npts = rand(rng, Poisson(recombMean))
    if npts == 0 || npts >= L
        copyto!(dest, rand(rng, Bool) ? chrA : chrB)
        return dest
    end
    empty!(pts)                        # scratch for the crossover points, reused by reproduce!
    while length(pts) < npts
        c = rand(rng, 0:L-2)
        c in pts || push!(pts, c)
    end
    cp, other = rand(rng, Bool) ? (chrA, chrB) : (chrB, chrA)
    # Copy whole blocks between crossovers: loci 0:c (1-based 1:c+1) come from the current
    # strand, which switches after locus c. The draws are those of a locus-by-locus copy.
    if !allpoints                      # QTSim: only the lowest point switches
        c = minimum(pts)
        copyto!(dest, 1, cp, 1, c + 1)
        copyto!(dest, c + 2, other, c + 2, L - c - 1)
        return dest
    end
    sort!(pts)
    lo = 1
    for c in pts
        copyto!(dest, lo, cp, lo, c + 2 - lo)
        lo = c + 2
        cp, other = other, cp
    end
    copyto!(dest, lo, cp, lo, L - lo + 1)
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
    pts = Int[]                      # crossover-point scratch, shared by every gamete
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
                    @view(w.genome[:, 2, p1]), rng, p.RecombPoissonMean, L, allpts, pts)
            gamete!(@view(w.sgenome[:, 2, j]), @view(w.genome[:, 1, p2]),
                    @view(w.genome[:, 2, p2]), rng, p.RecombPoissonMean, L, allpts, pts)
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
    rescale!(w)

Change the fitness width to `ChangeFitnessSD` at the end of generation
`ChangeFitnessGen`, so that it takes effect from the next generation (QTSim `RunSim`).
"""
function rescale!(w::World)
    w.gen == w.p.ChangeFitnessGen && (w.ω = w.p.ChangeFitnessSD)
    w
end

"""
    records(; demes = nothing)

Empty columns for the statistics [`record!`](@ref) collects, as `run!` returns them;
with `demes`, a vector of `World`s, also the deme of each row.
"""
records(; demes = nothing) = demes === nothing ?
    (gen = Int[], pheno_mean = Float64[], pheno_var = Float64[], geno_var = Float64[],
     nsurv = Int[], theta = Float64[], omega = Float64[], mean_fitness = Float64[],
     sel_diff = Float64[]) :
    (gen = Int[], deme = Int[], pheno_mean = Float64[], pheno_var = Float64[],
     geno_var = Float64[], nsurv = Int[], theta = Float64[], omega = Float64[],
     mean_fitness = Float64[], sel_diff = Float64[],
     demeof = IdDict{World,Int}(w => i for (i, w) in enumerate(demes)))

"""
    record!(rec, w)

Append the parent generation's statistics to the columns `rec` (see [`records`](@ref))
every `Sample` generations (QTSim `Output`); a no-op when `rec` is `nothing`. Variances
are sample variances (n − 1), as GSL's. Besides QTSim's columns it records the mean
fitness `W̄` (one minus the load) and the selection differential, the survivors' mean
phenotype minus the mean phenotype.
"""
function record!(rec, w::World)
    (rec === nothing || w.gen % w.p.Sample != 0) && return w
    push!(rec.gen, w.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[w])
    push!(rec.pheno_mean, mean(w.phenotype))
    push!(rec.pheno_var, var(w.phenotype)); push!(rec.geno_var, var(w.genotype))
    push!(rec.nsurv, w.nsurv); push!(rec.theta, w.θ); push!(rec.omega, w.ω)
    push!(rec.mean_fitness, mean(w.fitness))
    push!(rec.sel_diff, mean(w.phenotype[i] for i in w.survivors) - mean(w.phenotype))
    w
end

"""
    generation!(w, rng; rec = nothing)

One generation of QTSim's life cycle, `develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅
mutate ⋅ swap`, advancing the offspring into the current buffer and `w.gen` by one;
`record` appends to `rec` (see [`records`](@ref)) when one is given. A no-op past
extinction. On a fresh `World`, call [`init_env!`](@ref) first: `run!` is `init_env!`
followed by `NoGenerations` of these.
"""
function generation!(w::World, rng; rec = nothing)
    w.extinct && return w
    develop!(w)
    select!(w, rng)
    w.extinct && return w
    reproduce!(w, rng)
    record!(rec, w)
    rescale!(w)
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
    init_env!(w, rng)
    rec = records()
    for it in 0:p.NoGenerations-1
        w.gen = it
        generation!(w, rng; rec = record ? rec : nothing)
        w.extinct && return (extinct=true, extinct_gen=it, rec...)
    end
    (extinct=false, extinct_gen=-1, rec...)
end

"""
    simulate(p::Params; seed=p.Seed)

Run one simulation from `Params` with a `MersenneTwister(seed)`, returning the recorded
statistics (see `run!`) and the `seed`. Without a seed one is drawn at random, and the
returned `seed` reproduces the run.

The same seed gives the same run with the same versions of Julia and Distributions.jl:
the streams of `rand(rng, 1:n)` and of the Poisson sampler are not guaranteed across
their releases.
"""
function simulate(p::Params; seed=p.Seed)
    seed === nothing && (seed = Int(rand(RandomDevice(), UInt32)))
    merge(run!(World(p), MersenneTwister(seed)), (; seed))
end
