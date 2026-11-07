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
is the sum over all its loci and chromosomes (QTSim `CaclGenotype`). The current
population is individuals `1:N`, with `N == PopSize` always under QTSim's `Demography =
:fixed`; under density regulation `N` varies, `PopSize` is the carrying capacity, and
`Noff` is the number of offspring in the buffer. `θ0` and `ε` hold the state of a
fluctuating optimum (`MoveType = :red_noise` or `:periodic`, see [`moveoptimum!`](@ref)).
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
    N::Int
    Noff::Int
    θ0::Float64          # the optimum when a fluctuating move starts
    ε::Float64           # the optimum's AR(1) noise, MoveType = :red_noise
end

function World(p::Params)
    checkparams(p)
    N, L, P = p.PopSize, p.NoQTL, p.Ploidy
    World(p, 0,
          zeros(L, P, N), zeros(N),
          zeros(N), zeros(N), zeros(N),
          Int[], 0, 0.0, initFitSD(p),
          zeros(L, P, N), zeros(N), false, N, N, 0.0, 0.0)
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
    @inbounds for i in 1:w.N
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
    @inbounds for i in 1:w.N
        w.genotype[i] = sum(@view w.genome[:, :, i])
        w.phenotype[i] = w.genotype[i] + w.env[i]
    end
    w
end

"""
    moveoptimum!(w, rng = nothing)

Move the optimum `θ` if the schedule says so (QTSim `SetFitnessMean`). From generation
`t0 = MoveFitnessSDGen` on, with `d = MoveFitnessSD`, `z̄` the mean phenotype and `t` the
generation:

- `:phenotype`: `θ = z̄ + √V_P·d`;
- `:fitness`: `θ = z̄ + ω·d`;
- `:continuous`: `θ += d`;
- `:red_noise` (Julia-only): `θ_t = θ0 + d (t − t0) + ε_t`, where `θ0` is the optimum
  before the move and `ε` is a stationary AR(1) process with standard deviation
  `σ = MoveNoiseSD` and autocorrelation `ρ = MoveAutocorr`: `ε_{t0} ~ N(0, σ²)` and
  `ε_t = ρ ε_{t−1} + σ √(1 − ρ²) ξ_t` with `ξ_t ~ N(0, 1)`. One normal draw per generation,
  from `rng`;
- `:periodic` (Julia-only): `θ_t = θ0 + d (t − t0) + A sin(2π (t − t0)/P)`, with
  `A = MoveAmplitude` and `P = MovePeriod`. Draws nothing.

The two fluctuating moves are exogenous: θ does not depend on the population. Their
trend is `d` per generation, starting from `θ0` at `t0` (so a step behind `:continuous`,
which already adds `d` at `t0`). With `ρ = 0` the noise is white; as `ρ → 1` it
wanders more slowly with the same variance, so `σ` and `ρ` set the size and the time
scale (`−1/ln ρ` generations) of the fluctuations separately.
"""
function moveoptimum!(w::World, rng = nothing)
    p = w.p
    (p.MoveFitnessSDGen != -1 && w.gen >= p.MoveFitnessSDGen) || return w
    if p.MoveType === :red_noise || p.MoveType === :periodic
        t = w.gen - p.MoveFitnessSDGen
        t == 0 && (w.θ0 = w.θ)
        if p.MoveType === :red_noise
            rng === nothing && throw(ArgumentError("MoveType = :red_noise draws from a random number generator"))
            σ, ρ = p.MoveNoiseSD, p.MoveAutocorr
            w.ε = t == 0 ? σ * randn(rng) : ρ * w.ε + σ * sqrt(1 - ρ^2) * randn(rng)
        end
        w.θ = fluctuating_optimum(p, w.θ0, t, w.ε)
        return w
    end
    m = mean(view(w.phenotype, 1:w.N))
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
    fluctuating_optimum(p, θ0, t, ε = 0.0)

The optimum of a fluctuating move `t` generations after it starts, from `θ0`:
`θ0 + d t + ε` for `MoveType = :red_noise`, with the current noise `ε`, and
`θ0 + d t + A sin(2π t/P)` for `:periodic` (see [`moveoptimum!`](@ref)). With `ε = 0` it
is the expected optimum under red noise.
"""
function fluctuating_optimum(p::Params, θ0::Real, t::Real, ε::Real = 0.0)
    trend = θ0 + p.MoveFitnessSD * t
    p.MoveType === :red_noise && return trend + ε
    p.MoveType === :periodic && return trend + p.MoveAmplitude * sin(2π * t / p.MovePeriod)
    throw(ArgumentError("fluctuating_optimum is for MoveType :red_noise or :periodic, got :$(p.MoveType)"))
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
    moveoptimum!(w, rng)
    sd = p.StandardisePhenotype ? std(view(w.phenotype, 1:w.N)) : 1.0
    @inbounds for i in 1:w.N
        if p.Drift
            w.fitness[i] = 1.0
        else
            x = (w.phenotype[i] - w.θ) / sd
            w.fitness[i] = pdf(Normal(0.0, w.ω), x) / (0.4 / w.ω)
        end
    end
    empty!(w.survivors)
    @inbounds for i in 1:w.N
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
    noffspring(w, rng)

The number of offspring the survivors produce. Under QTSim's `Demography = :fixed` it is
`PopSize`, whatever the number of survivors, and draws nothing. Under `:beverton_holt`,
with `n` survivors, carrying capacity `K = PopSize` and `R = GrowthRate` offspring per
survivor at low density, it is `Poisson(R n / (1 + (R − 1) n / K))`, at most `K`: the
population stays near `K` while survivors are plentiful, and declines when `R·W̄ < 1`.
"""
function noffspring(w::World, rng)
    p = w.p
    p.Demography === :fixed && return p.PopSize
    R, K, n = p.GrowthRate, p.PopSize, w.nsurv
    min(K, rand(rng, Poisson(R * n / (1 + (R - 1) * n / K))))
end

"""
    reproduce!(w, rng)

Fill the offspring buffer `sgenome`/`senv` by fitness-independent random mating
among survivors (QTSim `NewIn`). Haploid: clone a random survivor. Diploid: two
distinct survivors each contribute a (recombined) gamete. Offspring get a fresh
environment. How many offspring is [`noffspring`](@ref)'s: `PopSize`, as in QTSim, or
under density regulation a number that depends on the survivors.
"""
function reproduce!(w::World, rng)
    p = w.p
    L = p.NoQTL
    esd = envsd(p)
    allpts = p.CrossoverModel === :all
    n = w.Noff = noffspring(w, rng)
    n == 0 && (w.extinct = true; return w)
    pts = Int[]                      # crossover-point scratch, shared by every gamete
    if p.Ploidy == 1
        @inbounds for j in 1:n
            par = w.survivors[rand(rng, 1:w.nsurv)]
            @views w.sgenome[:, 1, j] .= w.genome[:, 1, par]
            w.senv[j] = esd * randn(rng)
        end
    else
        @inbounds for j in 1:n
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
    @inbounds for j in 1:w.Noff
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
    w.N = w.Noff
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
     sel_diff = Float64[], popsize = Int[]) :
    (gen = Int[], deme = Int[], pheno_mean = Float64[], pheno_var = Float64[],
     geno_var = Float64[], nsurv = Int[], theta = Float64[], omega = Float64[],
     mean_fitness = Float64[], sel_diff = Float64[], popsize = Int[],
     demeof = IdDict{World,Int}(w => i for (i, w) in enumerate(demes)))

"""
    record!(rec, w)

Append the parent generation's statistics to the columns `rec` (see [`records`](@ref))
every `Sample` generations (QTSim `Output`); a no-op when `rec` is `nothing`. Variances
are sample variances (n − 1), as GSL's. Besides QTSim's columns it records the mean
fitness `W̄` (one minus the load) and the selection differential, the survivors' mean
phenotype minus the mean phenotype, and the population size `N` (`PopSize`, unless the
population is density-regulated).
"""
function record!(rec, w::World)
    (rec === nothing || w.gen % w.p.Sample != 0) && return w
    push!(rec.gen, w.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[w])
    z = view(w.phenotype, 1:w.N)
    push!(rec.pheno_mean, mean(z))
    push!(rec.pheno_var, var(z)); push!(rec.geno_var, var(view(w.genotype, 1:w.N)))
    push!(rec.nsurv, w.nsurv); push!(rec.theta, w.θ); push!(rec.omega, w.ω)
    push!(rec.mean_fitness, mean(view(w.fitness, 1:w.N)))
    push!(rec.sel_diff, mean(w.phenotype[i] for i in w.survivors) - mean(z))
    push!(rec.popsize, w.N)
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
    w.extinct && return w            # no offspring (density regulation): no record, as at selection
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
