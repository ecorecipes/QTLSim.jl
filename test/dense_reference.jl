# A test-only copy of the dense kernel, src/model.jl as of commit 76c6058, before genomes
# were stored as shared blocks. Genomes are a dense (NoQTL, Ploidy, PopSize) array and an
# offspring double buffer. test/runtests.jl runs it beside the block kernel on the same
# seeds: allele values must agree bit for bit, genotypes to rounding. Do not change it,
# except to follow a change to the random draws of the kernel that is meant.
module DenseReference

using Random, Statistics, Distributions
using QTLSim: Params, checkparams, initFitSD, envVar

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

function init_env!(w::World, rng)
    esd = envsd(w.p)
    @inbounds for i in 1:w.N
        w.env[i] = esd * randn(rng)
    end
    w
end

# --- operators (morphisms K → K) ----------------------------------------

function develop!(w::World)
    @inbounds for i in 1:w.N
        w.genotype[i] = sum(@view w.genome[:, :, i])
        w.phenotype[i] = w.genotype[i] + w.env[i]
    end
    w
end

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

function fluctuating_optimum(p::Params, θ0::Real, t::Real, ε::Real = 0.0)
    trend = θ0 + p.MoveFitnessSD * t
    p.MoveType === :red_noise && return trend + ε
    p.MoveType === :periodic && return trend + p.MoveAmplitude * sin(2π * t / p.MovePeriod)
    throw(ArgumentError("fluctuating_optimum is for MoveType :red_noise or :periodic, got :$(p.MoveType)"))
end

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

function noffspring(w::World, rng)
    p = w.p
    p.Demography === :fixed && return p.PopSize
    R, K, n = p.GrowthRate, p.PopSize, w.nsurv
    min(K, rand(rng, Poisson(R * n / (1 + (R - 1) * n / K))))
end

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

function swapbuffers!(w::World)
    w.genome, w.sgenome = w.sgenome, w.genome
    w.env, w.senv = w.senv, w.env
    w.N = w.Noff
    w
end

function rescale!(w::World)
    w.gen == w.p.ChangeFitnessGen && (w.ω = w.p.ChangeFitnessSD)
    w
end

records(; demes = nothing) = demes === nothing ?
    (gen = Int[], pheno_mean = Float64[], pheno_var = Float64[], geno_var = Float64[],
     nsurv = Int[], theta = Float64[], omega = Float64[], mean_fitness = Float64[],
     sel_diff = Float64[], popsize = Int[]) :
    (gen = Int[], deme = Int[], pheno_mean = Float64[], pheno_var = Float64[],
     geno_var = Float64[], nsurv = Int[], theta = Float64[], omega = Float64[],
     mean_fitness = Float64[], sel_diff = Float64[], popsize = Int[],
     demeof = IdDict{World,Int}(w => i for (i, w) in enumerate(demes)))

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

end # module
