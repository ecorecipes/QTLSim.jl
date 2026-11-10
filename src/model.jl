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
    TraitState

The state of a population with `n > 1` traits (`Params` `NoTraits`), which a `World`
carries in `traits` (`nothing` for one trait): genotypes `G`, phenotypes `Z` and the
parents' and offspring's environmental deviations `E` and `SE`, `n × PopSize` each
(column `i` is individual `i`), the optimum vector `θ`, the inverse `Sinv` of the fitness
function's shape (`Ω = ω² S`), a square root `Lmut` of the covariance of one mutation's
effects (`Lmut·Lmutᵀ = M/(Ploidy·NoQTL·rate)`), the unit move direction `u`, scratch
vectors `ξ` and `d`, and for a fluctuating move the optimum `θ0` when it starts, the AR(1)
noise `ε`, a square root `Lnoise` of its stationary covariance and the unit direction `f`
of the cycle.
"""
mutable struct TraitState
    n::Int
    G::Matrix{Float64}
    Z::Matrix{Float64}
    E::Matrix{Float64}
    SE::Matrix{Float64}
    θ::Vector{Float64}
    Sinv::Matrix{Float64}
    Lmut::Matrix{Float64}
    u::Vector{Float64}
    ξ::Vector{Float64}
    d::Vector{Float64}
    θ0::Vector{Float64}
    ε::Vector{Float64}
    Lnoise::Matrix{Float64}
    f::Vector{Float64}
end

# a symmetric square root of a positive semidefinite matrix, which may be singular
function psd_sqrt(A::AbstractMatrix)
    F = eigen(Symmetric(Matrix{Float64}(A)))
    F.vectors * (sqrt.(max.(F.values, 0.0)) .* F.vectors')
end

function TraitState(p::Params, N::Integer)
    n = p.NoTraits
    Mper = mutational_covariance(p) ./ (p.Ploidy * p.NoQTL * p.MutationRatePerGamete)
    TraitState(n, zeros(n, N), zeros(n, N), zeros(n, N), zeros(n, N), zeros(n),
               inv(Symmetric(selection_shape(p))) |> Matrix, psd_sqrt(Mper), move_direction(p),
               zeros(n), zeros(n), zeros(n), zeros(n), psd_sqrt(move_noise_covariance(p)), fluct_direction(p))
end

"""
    World(p; blocksize = default_blocksize(p.NoQTL), pool = nothing)

Mutable population state carried through the operators. Holds the current parents'
genomes and environments, per-individual scratch (genotype, phenotype, fitness), the
survivor index set, the current selection optimum `θ` and width `ω`, and a double buffer
(`schr`/`senv`) for the offspring generation. Under `Selection = :fecundity`, [`select!`](@ref)
sets `weighted` and leaves the survivors' mating weights, their fitnesses, as cumulative
sums in `cumw` with total `wsum`, for [`reproduce!`](@ref) and [`noffspring`](@ref).

The current population is individuals `1:N`, with `N == PopSize` always under QTSim's
`Demography = :fixed`; under density regulation `N` varies, `PopSize` is the carrying
capacity, and `Noff` is the number of offspring in the buffer. `θ0` and `ε` hold the
state of a fluctuating optimum (`MoveType = :red_noise` or `:periodic`, see
[`moveoptimum!`](@ref)).

# Genomes

Each chromosome of `NoQTL` loci is split into `nb` blocks of `blocksize` loci (the last
may be shorter), `⌈√NoQTL⌉` by default, and stored as `nb` ids of blocks in a
pool (`QTLSim.BlockPool`): `w.chr[k, c, i]` is block `k` of chromosome `c` of
individual `i`, an `(nb, Ploidy, PopSize)` array, and `w.schr` the same for the
offspring. Blocks are shared and reference-counted. A gamete shares every block of its
parent's chromosomes except the one that holds a crossover, which it builds; a mutation
changes a block in place if nothing else holds it and copies it first otherwise. Each
block caches the sum of its values, and the genotype of individual `i` (QTSim
`CaclGenotype`) is the sum of its blocks' sums, chromosome by chromosome, block by
block. So a generation costs about `N·√L` operations rather than `N·L`.

Read the genomes with [`genomes`](@ref) and write them with [`setgenomes!`](@ref), which
convert to and from a dense `(NoQTL, Ploidy, n)` array.

Every slot of both buffers, all `PopSize` columns, holds one reference to a live block
at all times. A slot drops its reference only when it is overwritten, after it has
taken the new one: by [`reproduce!`](@ref) in the offspring buffer, which after a swap
holds the previous parents, by [`mutate!`](@ref) when it copies a block, by
`setgenomes!`, and by migration between demes that do not share a pool. So the previous
parents' blocks are released one offspring slot at a time during the next
`reproduce!`, never while a slot names them, and `swapbuffers!` moves no references.
Demes of a metapopulation made by `⊗` share one pool, so that a migrant is moved by
swapping its block ids.

The random draws, and so the allele values, are those of a dense array of alleles
(bit for bit, for the same seed). Genotypes are summed in another order, so they can
differ from a sum over the dense array in the last bits (about 1e-16 relative).
"""
mutable struct World
    p::Params
    gen::Int
    pool::BlockPool
    chr::Array{Int32,3}     # parents: block ids, (nb, Ploidy, PopSize)
    env::Vector{Float64}
    genotype::Vector{Float64}
    phenotype::Vector{Float64}
    fitness::Vector{Float64}
    survivors::Vector{Int}
    nsurv::Int
    θ::Float64
    ω::Float64
    schr::Array{Int32,3}    # offspring: block ids
    senv::Vector{Float64}
    extinct::Bool
    N::Int
    Noff::Int
    θ0::Float64          # the optimum when a fluctuating move starts
    ε::Float64           # the optimum's AR(1) noise, MoveType = :red_noise
    pts::Vector{Int}     # scratch for the crossover points of a gamete
    weighted::Bool       # the survivors mate in proportion to their fitness (Selection = :fecundity)
    wsum::Float64        # then: the sum of their fitnesses, ΣW
    cumw::Vector{Float64}   # and the cumulative sums of the fitnesses of 1:N
    traits::Union{Nothing,TraitState}   # NoTraits > 1: the traits (genotype etc. then hold trait 1)
    plastic::Float64     # the plastic part (b/B)·θ_{t−τ} of every phenotype, once select has expressed it (Plasticity ≠ 0)
    cues::Vector{Float64}   # the optima of the last CueLag + 1 generations, the latest first (a reaction norm's cue)
    cuegen::Int          # the generation whose optimum was last entered in `cues`
end

function World(p::Params; blocksize::Integer = default_blocksize(p.NoQTL),
               pool::Union{BlockPool,Nothing} = nothing)
    checkparams(p)
    N, L, P = p.PopSize, p.NoQTL, p.Ploidy
    if pool === nothing
        pool = BlockPool(L, blocksize; traits = p.NoTraits)
    else
        pool.L == L || throw(ArgumentError("a pool of chromosomes of $(pool.L) loci for NoQTL = $L"))
        pool.n == p.NoTraits || throw(ArgumentError("a pool of $(pool.n) traits for NoTraits = $(p.NoTraits)"))
    end
    nb, n = pool.nb, pool.n
    # every slot of both buffers starts on one block of zeros per block position
    chr = Array{Int32,3}(undef, nb, P, N)
    for k in 1:nb
        id = newblock!(pool)
        fill!(view(pool.vals, blockoffset(pool, id) .+ (1:n*blocklen(pool, k))), 0.0)
        fill!(view(pool.sums, (Int(id) - 1) * n .+ (1:n)), 0.0)
        pool.refs[id] = Int32(2 * P * N)
        chr[k, :, :] .= id
    end
    World(p, 0, pool, chr, zeros(N),
          zeros(N), zeros(N), zeros(N),
          sizehint!(Int[], N), 0, 0.0, initFitSD(p),
          copy(chr), zeros(N), false, N, N, 0.0, 0.0, sizehint!(Int[], 16), false, 0.0, zeros(N),
          p.NoTraits == 1 ? nothing : TraitState(p, N), 0.0, sizehint!(Float64[], p.CueLag + 1), -1)
end

"""
    shared_copy(w)

A copy of `w` whose genomes are blocks of the same pool: the copy's slots name the same
blocks, one more reference each. `World ⊗ demes(k)` makes its demes this way.
"""
function shared_copy(w::World)
    v = World(w.p; pool = w.pool)
    for f in fieldnames(World)
        f in (:p, :pool, :chr, :schr) && continue
        setfield!(v, f, deepcopy(getfield(w, f)))
    end
    for (dst, src) in ((v.chr, w.chr), (v.schr, w.schr)), i in eachindex(src)
        old = dst[i]
        dst[i] = retain!(w.pool, src[i])
        release!(w.pool, old)
    end
    v
end

"""
    genomes(w; offspring = false)

The genomes of the current population as a dense `(NoQTL, Ploidy, N)` array of allele
effects, a copy: `genomes(w)[l, c, i]` is the effect at locus `l` of chromosome `c` of
individual `i`. With `offspring = true`, those of the `Noff` offspring in the buffer.
With `n = NoTraits > 1` traits, an `(n, NoQTL, Ploidy, N)` array: `genomes(w)[r, l, c, i]`
is the effect on trait `r`.
"""
function genomes(w::World; offspring::Bool = false)
    pool, chr = w.pool, offspring ? w.schr : w.chr
    n = offspring ? w.Noff : w.N
    nt = pool.n
    G = nt == 1 ? Array{Float64,3}(undef, pool.L, size(chr, 2), n) :
                  Array{Float64,4}(undef, nt, pool.L, size(chr, 2), n)
    @inbounds for i in 1:n, c in axes(chr, 2), k in 1:pool.nb
        copyto!(G, nt * ((k - 1) * pool.B + pool.L * ((c - 1) + size(chr, 2) * (i - 1))) + 1,
                pool.vals, blockoffset(pool, chr[k, c, i]) + 1, nt * blocklen(pool, k))
    end
    G
end

"""
    setgenomes!(w, G; offspring = false)

Write the dense `(NoQTL, Ploidy, n)` array `G` into the genomes of individuals `1:n` of
the current population (with `offspring = true`, of the offspring buffer); `n` may be at
most `PopSize`. It changes neither `N` nor the genotypes: call [`develop!`](@ref) to
bring the genotypes and phenotypes up to date.
"""
function setgenomes!(w::World, G::AbstractArray{<:Real,3}; offspring::Bool = false)
    pool, chr = w.pool, offspring ? w.schr : w.chr
    pool.n == 1 || throw(ArgumentError("genomes of $(pool.n) traits are an (n, NoQTL, Ploidy, N) array"))
    size(G, 1) == pool.L && size(G, 2) == size(chr, 2) && size(G, 3) <= size(chr, 3) ||
        throw(ArgumentError("genomes of size $(size(G)) for $(pool.L) loci, ploidy $(size(chr, 2)) and at most $(size(chr, 3)) individuals"))
    for i in axes(G, 3), c in axes(G, 2)
        setchromosome!(w.pool, chr, c, i, view(G, :, c, i))
    end
    w
end

function setgenomes!(w::World, G::AbstractArray{<:Real,4}; offspring::Bool = false)
    pool, chr = w.pool, offspring ? w.schr : w.chr
    size(G, 1) == pool.n && size(G, 2) == pool.L && size(G, 3) == size(chr, 2) && size(G, 4) <= size(chr, 3) ||
        throw(ArgumentError("genomes of size $(size(G)) for $(pool.n) traits, $(pool.L) loci, ploidy $(size(chr, 2)) and at most $(size(chr, 3)) individuals"))
    for i in axes(G, 4), c in axes(G, 3)
        setchromosome!(w.pool, chr, c, i, view(G, :, :, c, i))
    end
    w
end

# Write the values `x` (length L, or n×L for n traits) into chromosome `c` of individual
# `i` of `chr`, as new blocks.
function setchromosome!(pool::BlockPool, chr::Array{Int32,3}, c, i, x::AbstractArray)
    n = pool.n
    for k in 1:pool.nb
        id = newblock!(pool)
        o, lo, len = blockoffset(pool, id), (k - 1) * pool.B * n, n * blocklen(pool, k)
        for t in 1:len
            pool.vals[o + t] = x[lo + t]
        end
        blocksum!(pool, id, blocklen(pool, k))
        setslot!(pool, chr, k, c, i, id)
    end
end

"The loci and ploidy of a `World`'s genomes (and the traits, when there are several)."
genome_shape(w::World) = w.p.NoTraits == 1 ? (w.p.NoQTL, w.p.Ploidy) : (w.p.NoQTL, w.p.Ploidy, w.p.NoTraits)

"""
    exchange_genomes!(w, i, v, j)

Exchange the genomes of individual `i` of `w` and individual `j` of `v`, which have the
same loci and ploidy. Swaps block ids, `nb·Ploidy` of them, when the two share a pool;
otherwise copies the blocks into the other pool, or, if the block sizes differ, the
alleles.
"""
function exchange_genomes!(w::World, i::Integer, v::World, j::Integer)
    genome_shape(w) == genome_shape(v) ||
        throw(ArgumentError("exchanging genomes of shapes $(genome_shape(w)) and $(genome_shape(v))"))
    pw, pv = w.pool, v.pool
    if pw === pv
        @inbounds for c in axes(w.chr, 2), k in axes(w.chr, 1)
            w.chr[k, c, i], v.chr[k, c, j] = v.chr[k, c, j], w.chr[k, c, i]
        end
    elseif pw.B == pv.B
        for c in axes(w.chr, 2), k in axes(w.chr, 1)
            a, b, len = w.chr[k, c, i], v.chr[k, c, j], blocklen(pw, k)
            setslot!(pv, v.chr, k, c, j, copyblock!(pv, pw, a, len))
            setslot!(pw, w.chr, k, c, i, copyblock!(pw, pv, b, len))
        end
    else
        pw.n == 1 || throw(ArgumentError("exchanging genomes of several traits between pools of different block sizes"))
        gi = genomes_of(w, i); gj = genomes_of(v, j)
        for c in axes(w.chr, 2)
            setchromosome!(pw, w.chr, c, i, view(gj, :, c))
            setchromosome!(pv, v.chr, c, j, view(gi, :, c))
        end
    end
    nothing
end

# the (NoQTL, Ploidy) genome of one parent
function genomes_of(w::World, i::Integer)
    pool = w.pool
    g = Matrix{Float64}(undef, pool.L, size(w.chr, 2))
    for c in axes(w.chr, 2), k in 1:pool.nb
        copyto!(g, (c - 1) * pool.L + (k - 1) * pool.B + 1,
                pool.vals, blockoffset(pool, w.chr[k, c, i]) + 1, blocklen(pool, k))
    end
    g
end

"""
    check_blocks(ws::World...)

Check the reference counts of the pool the `World`s `ws` share, which must be all the
`World`s that use it, against their buffers (see [`check_blocks`](@ref)). Returns the
number of live blocks. For tests.
"""
function check_blocks(ws::World...)
    pool = first(ws).pool
    all(w -> w.pool === pool, ws) || throw(ArgumentError("the Worlds do not share one pool"))
    check_blocks(pool, Iterators.flatten((w.chr, w.schr) for w in ws)...)
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
    t = w.traits
    if t !== nothing         # several traits: N(0, V_E) for each, independently
        @inbounds for i in 1:w.N, r in 1:t.n
            t.E[r, i] = esd * randn(rng)
        end
        return w
    end
    @inbounds for i in 1:w.N
        w.env[i] = esd * randn(rng)
    end
    w
end

# --- operators (morphisms K → K) ----------------------------------------

"""
    develop!(w)

Genotype = Σ allele effects; phenotype = genotype + environment (QTSim
`CaclPopPhenotype`). Deterministic given the environments drawn at birth. The sum of the
effects is the sum of the cached sums of the individual's blocks, chromosome by
chromosome, block by block (see [`World`](@ref)).

`develop!`, [`select!`](@ref) and [`mutate!`](@ref) serve both models: they call the
one-trait operator (`develop_one!` etc.) or the operator of several traits
(`develop_traits!` etc., src/traits.jl) by the `World`. The two models' tables of kernel
images (ext/QTLSimCatlabExt/compose.jl) call those directly.
"""
develop!(w::World) = w.traits === nothing ? develop_one!(w) : develop_traits!(w)

# develop! of QTSim's one-trait model
function develop_one!(w::World)
    sums, chr = w.pool.sums, w.chr
    nb, P = size(chr, 1), size(chr, 2)
    m = nb * P
    @inbounds for i in 1:w.N
        g = 0.0
        for t in m*(i-1)+1:m*i       # chromosome by chromosome, block by block
            g += sums[chr[t]]
        end
        w.genotype[i] = g
        w.phenotype[i] = g + w.env[i]
    end
    w.plastic = 0.0          # the reaction norm, if any, is expressed in select
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
    w.traits === nothing ? moveoptimum_one!(w, rng) : moveoptimum_traits!(w, rng)
end

# moveoptimum! of the one-trait model, from the move's first generation on
function moveoptimum_one!(w::World, rng)
    p = w.p
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

Under `Selection = :fecundity` (Julia-only) everyone survives and draws nothing here:
the survivors are all `N` individuals, and their fitnesses become mating weights
(`w.weighted`, `w.cumw`, `w.wsum`), which [`reproduce!`](@ref) samples parents by. The
population dies if fewer than `Ploidy` individuals have positive fitness (no two distinct
parents can be drawn), which at a fixed size happens only when fitness underflows to 0.

Fitness is the normalized Gaussian `pdf(Normal(0,ω), (z−θ)/sd) / (0.4/ω)`,
replicating QTSim exactly (the `0.4 ≈ 1/√(2π)` constant included); the idealized
form `exp(−(z−θ)²/2ω²)` is what Lean `StabilizingSelection.lean` proves about.

With phenotypic plasticity (`Plasticity = b ≠ 0`, Julia-only) the reaction norm is
expressed here, once the optimum of the generation is set, since the optimum is the
measure of the environment (`θ = B·ε`): every phenotype gains the plastic part
`(b/B)·θ_{t−τ}`, `τ = CueLag` (`QTLSim.express_plasticity!`), so that the recorded
phenotypes include it, and the cost of plasticity multiplies every fitness by
`exp(−b²/(2ω_b²))`. With a heritable reaction norm (`ReactionNorm`) the image is
`QTLSim.select_norm!`.
"""
select!(w::World, rng) = w.traits === nothing ? select_one!(w, rng) :
                         w.p.ReactionNorm ? select_norm!(w, rng) : select_traits!(w, rng)

# select! of the one-trait model: move the optimum, QTSim's fitness, the survivors
function select_one!(w::World, rng)
    p = w.p
    (p.MoveFitnessSDGen != -1 && w.gen >= p.MoveFitnessSDGen) && moveoptimum_one!(w, rng)
    p.Plasticity == 0 || express_plasticity!(w)
    sd = p.StandardisePhenotype ? std(view(w.phenotype, 1:w.N)) : 1.0
    @inbounds for i in 1:w.N
        if p.Drift
            w.fitness[i] = 1.0
        else
            x = (w.phenotype[i] - w.θ) / sd
            w.fitness[i] = pdf(Normal(0.0, w.ω), x) / (0.4 / w.ω)
        end
    end
    p.Plasticity == 0 || p.Drift || plasticity_cost!(w)
    survive!(w, rng)
end

# select! of several traits: move the optimum vector, the Gaussian fitness in n dimensions
# (src/traits.jl), the survivors
function select_traits!(w::World, rng)
    p = w.p
    (p.MoveFitnessSDGen != -1 && w.gen >= p.MoveFitnessSDGen) && moveoptimum_traits!(w, rng)
    fitness_traits!(w)
    survive!(w, rng)
end

# the survivors of viability selection (QTSim SetSurviving), or under fecundity selection
# everyone, with their fitnesses as mating weights; the same in both models
function survive!(w::World, rng)
    p = w.p
    empty!(w.survivors)
    if p.Selection === :fecundity
        return weigh!(w)
    end
    w.weighted = false
    @inbounds for i in 1:w.N
        rand(rng) < w.fitness[i] && push!(w.survivors, i)
    end
    w.nsurv = length(w.survivors)
    if (p.Ploidy == 1 && w.nsurv == 0) || (p.Ploidy == 2 && w.nsurv < 2)
        w.extinct = true
    end
    w
end

# Fecundity selection: all of 1:N survive, with their fitnesses as mating weights.
function weigh!(w::World)
    N = w.N
    append!(w.survivors, 1:N)
    w.nsurv = N
    length(w.cumw) < N && resize!(w.cumw, N)
    s, npos = 0.0, 0
    @inbounds for i in 1:N
        s += w.fitness[i]
        w.cumw[i] = s
        npos += w.fitness[i] > 0
    end
    w.wsum, w.weighted = s, true
    npos < w.p.Ploidy && (w.extinct = true)
    w
end

# A parent drawn with probability W_i/ΣW: the first i with cumw[i] > u, u ~ U(0, ΣW), so
# that an individual of fitness 0 is never drawn.
@inline function weighted_parent(w::World, rng)
    u = rand(rng) * w.wsum
    min(searchsortedlast(view(w.cumw, 1:w.N), u) + 1, w.N)
end

# A parent other than `p1`, drawn with probability W_i/Σ_{j≠p1} W_j, by a scan of the
# fitnesses (not of their cumulative sums, which lose the others' weight to rounding
# when p1's dominates). Redrawing until the parent differs from p1 has this distribution
# too, also after any number of failed redraws, so switching to it keeps the law.
function weighted_parent_except(w::World, rng, p1::Int)
    rest = 0.0
    @inbounds for i in 1:w.N
        i == p1 || (rest += w.fitness[i])
    end
    u = rand(rng) * rest
    acc, last = 0.0, 0
    @inbounds for i in 1:w.N
        (i == p1 || w.fitness[i] <= 0) && continue
        acc += w.fitness[i]; last = i
        acc > u && return i
    end
    last
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

This is the gamete on dense arrays of alleles. The kernel ([`reproduce!`](@ref)) builds
gametes from shared blocks, with the same random draws in the same order, so it copies the
same alleles.
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

# The kernel's gamete: `gamete!` on block ids. Writes chromosome `c` of offspring `j` (a
# slot column of `w.schr`) from the two chromosomes of parent `par` in `w.chr`, with the
# random draws of `gamete!`, in its order, so that the alleles are those `gamete!` would
# copy. Blocks wholly on one side of every switch are shared; a block with a switch inside
# it is built from the two strands, unless both strands hold the same block there.
function gamete_blocks!(w::World, c::Int, j::Int, par::Int, rng, recombMean::Float64,
                        allpoints::Bool)
    pool, chr, schr, pts = w.pool, w.chr, w.schr, w.pts
    L, B, nb, n = pool.L, pool.B, pool.nb, pool.n
    P = size(chr, 2)
    # linear offsets: block k of strand s of the parent is chr[pb[s] + k], and block k
    # of the gamete goes to schr[db + k]
    pb1 = nb * P * (par - 1)
    pb2 = pb1 + nb
    db = nb * ((c - 1) + P * (j - 1))
    npts = rand(rng, Poisson(recombMean))
    if npts == 0 || npts >= L
        pb = rand(rng, Bool) ? pb1 : pb2
        @inbounds for k in 1:nb
            share!(pool, schr, db + k, chr[pb + k])
        end
        return nothing
    end
    empty!(pts)
    while length(pts) < npts
        x = rand(rng, 0:L-2)
        x in pts || push!(pts, x)
    end
    pb, po = rand(rng, Bool) ? (pb1, pb2) : (pb2, pb1)   # the current strand and the other
    if allpoints
        sort!(pts)
    else                                 # QTSim: only the lowest point switches
        pts[1] = minimum(pts)
        npts = 1
    end
    # A switch after 0-based locus x: loci 1:x+1 (1-based) come from the current strand.
    # Block k holds loci lo:hi; a switch with x+1 < hi falls inside it, x+1 == hi on its
    # end, which changes the strand for the blocks after it only.
    ip = 1
    @inbounds for k in 1:nb
        hi = k == nb ? L : k * B
        a = chr[pb + k]
        if ip <= npts && pts[ip] + 1 < hi
            if a == chr[po + k]          # the same block on both strands: share it
                while ip <= npts && pts[ip] + 1 < hi
                    pb, po = po, pb; ip += 1
                end
                share!(pool, schr, db + k, a)
            else
                lo = (k - 1) * B + 1
                id = newblock!(pool)
                # locus l (its n values) goes to vals[o + (l - lo)n + 1 : o + (l - lo + 1)n]
                o = blockoffset(pool, id)
                l = lo
                while ip <= npts && pts[ip] + 1 < hi
                    e = pts[ip] + 1
                    copyto!(pool.vals, o + (l - lo) * n + 1, pool.vals,
                            blockoffset(pool, chr[pb + k]) + (l - lo) * n + 1, (e - l + 1) * n)
                    l = e + 1; pb, po = po, pb; ip += 1
                end
                copyto!(pool.vals, o + (l - lo) * n + 1, pool.vals,
                        blockoffset(pool, chr[pb + k]) + (l - lo) * n + 1, (hi - l + 1) * n)
                blocksum!(pool, id, hi - lo + 1)
                setslot!(pool, schr, db + k, id)
            end
        else
            share!(pool, schr, db + k, a)
        end
        if ip <= npts && pts[ip] + 1 == hi
            pb, po = po, pb; ip += 1
        end
    end
    nothing
end

"""
    noffspring(w, rng)

The number of offspring the survivors produce. Under QTSim's `Demography = :fixed` it is
`PopSize`, whatever the number of survivors, and draws nothing. Under `:beverton_holt`,
with `n` survivors, carrying capacity `K = PopSize` and `R = GrowthRate` offspring per
survivor at low density, it is `Poisson(R n / (1 + (R − 1) n / K))`, at most `K`: the
population stays near `K` while survivors are plentiful, and declines when `R·W̄ < 1`.

Under `Selection = :fecundity` nobody dies at selection, and `n` is the sum of the
fitnesses, `ΣW = N·W̄`, a real number: the expected number of survivors under viability
selection with the same phenotypes. So the expected offspring number is the viability
one's with `n` replaced by its expectation. As `n ↦ R n/(1 + (R − 1) n/K)` is concave,
viability's expectation is lower, by at most `R (R − 1) Var(n)/K ≤ R (R − 1) N/(4K)`
offspring (Var(n) = Σ W(1 − W) ≤ N/4), about 0.5 for `R = 2` and `N = K`: negligible
beside the Poisson noise.
"""
function noffspring(w::World, rng)
    p = w.p
    p.Demography === :fixed && return p.PopSize
    R, K = p.GrowthRate, p.PopSize
    n = w.weighted ? w.wsum : w.nsurv
    min(K, rand(rng, Poisson(R * n / (1 + (R - 1) * n / K))))
end

"""
    reproduce!(w, rng)

Fill the offspring buffer by fitness-independent random mating among survivors (QTSim
`NewIn`). Haploid: clone a random survivor. Diploid: two distinct survivors each
contribute a (recombined) gamete, drawn as by [`gamete!`](@ref). After fecundity
selection (`w.weighted`, see [`select!`](@ref)) each parent is drawn with probability
proportional to its fitness instead, and a second parent equal to the first is redrawn
in the same way. Offspring get a fresh environment. How many offspring is
[`noffspring`](@ref)'s: `PopSize`, as in QTSim, or under density regulation a number that
depends on the survivors.

Offspring share their parents' blocks (see [`World`](@ref)); each slot it fills releases
the block it named before, from the generation before the parents.
"""
function reproduce!(w::World, rng)
    p = w.p
    esd = envsd(p)
    allpts = p.CrossoverModel === :all
    n = w.Noff = noffspring(w, rng)
    n == 0 && (w.extinct = true; return w)
    pool, chr, schr = w.pool, w.chr, w.schr
    w.weighted && return reproduce_weighted!(w, rng, n, esd, allpts)
    if p.Ploidy == 1
        @inbounds for j in 1:n
            par = w.survivors[rand(rng, 1:w.nsurv)]
            nb = size(chr, 1)
            for k in 1:nb
                share!(pool, schr, nb * (j - 1) + k, chr[nb * (par - 1) + k])
            end
            newenv!(w, j, esd, rng)
        end
    else
        @inbounds for j in 1:n
            p1 = w.survivors[rand(rng, 1:w.nsurv)]
            p2 = w.survivors[rand(rng, 1:w.nsurv)]
            while p2 == p1
                p2 = w.survivors[rand(rng, 1:w.nsurv)]
            end
            gamete_blocks!(w, 1, j, p1, rng, p.RecombPoissonMean, allpts)
            gamete_blocks!(w, 2, j, p2, rng, p.RecombPoissonMean, allpts)
            newenv!(w, j, esd, rng)
        end
    end
    w
end

# the environmental deviation(s) of offspring j: N(0, V_E), for each trait when there are several
@inline function newenv!(w::World, j::Int, esd::Float64, rng)
    t = w.traits
    if t === nothing
        @inbounds w.senv[j] = esd * randn(rng)
    else
        @inbounds for r in 1:t.n
            t.SE[r, j] = esd * randn(rng)
        end
    end
    nothing
end

# reproduce! after fecundity selection: parents drawn in proportion to their fitness
function reproduce_weighted!(w::World, rng, n, esd, allpts)
    p, pool, chr, schr = w.p, w.pool, w.chr, w.schr
    nb = size(chr, 1)
    @inbounds for j in 1:n
        if p.Ploidy == 1
            par = weighted_parent(w, rng)
            for k in 1:nb
                share!(pool, schr, nb * (j - 1) + k, chr[nb * (par - 1) + k])
            end
        else
            p1 = weighted_parent(w, rng)
            p2 = weighted_parent(w, rng)
            tries = 1
            while p2 == p1
                # nearly all the weight on p1 (the rest can be 1e-300 when fitness nearly
                # underflows): after 100 redraws, draw from the others exactly
                p2 = tries < 100 ? weighted_parent(w, rng) : weighted_parent_except(w, rng, p1)
                tries += 1
            end
            gamete_blocks!(w, 1, j, p1, rng, p.RecombPoissonMean, allpts)
            gamete_blocks!(w, 2, j, p2, rng, p.RecombPoissonMean, allpts)
        end
        newenv!(w, j, esd, rng)
    end
    w
end

"""
    mutate!(w, rng)

Add Gaussian mutations to the offspring buffer (QTSim `MutateIndividualChrNo`).
Each chromosome receives `Poisson(L·rate)` mutations, each of variance
`m/(ploidy·L·rate)·V_E`; in expectation this injects `m·V_E` of variance per
individual per generation (Lean `mutationalVarianceInput_eq`).

A mutation changes its block in place if the offspring is the only holder, and a copy of
it otherwise; the block's sum is then recomputed from its values.
"""
mutate!(w::World, rng) = w.traits === nothing ? mutate_one!(w, rng) : mutate_traits!(w, rng)

# mutate! of the one-trait model
function mutate_one!(w::World, rng)
    p = w.p
    L = p.NoQTL
    rate = p.MutationRatePerGamete
    mutsd = sqrt(max((p.MutationScalar / (p.Ploidy * L * rate)) * envVar(p), 0.0))
    λ = L * rate
    pool, schr = w.pool, w.schr
    B = pool.B
    @inbounds for j in 1:w.Noff
        for c in 1:p.Ploidy
            no = rand(rng, Poisson(λ))
            for _ in 1:no
                pos = rand(rng, 1:L)
                x = mutsd * randn(rng)
                k = (pos - 1) ÷ B + 1
                len = blocklen(pool, k)
                id = schr[k, c, j]
                if pool.refs[id] > 1     # shared: copy on write
                    nid = copyblock!(pool, pool, id, len)
                    pool.refs[id] -= Int32(1)
                    schr[k, c, j] = id = nid
                end
                pool.vals[blockoffset(pool, id) + pos - (k - 1) * B] += x
                blocksum!(pool, id, len)
            end
        end
    end
    w
end

"""
    swapbuffers!(w)

Swap the parent and offspring buffers (QTSim double-buffer `swap`). No block references
change: the old parents' blocks stay referenced by the offspring buffer until the next
[`reproduce!`](@ref) overwrites it.
"""
function swapbuffers!(w::World)
    w.chr, w.schr = w.schr, w.chr
    w.env, w.senv = w.senv, w.env
    t = w.traits
    t === nothing || ((t.E, t.SE) = (t.SE, t.E))
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
function records(; demes = nothing, traits::Integer = 1)
    rec = demes === nothing ?
        (gen = Int[], pheno_mean = Float64[], pheno_var = Float64[], geno_var = Float64[],
         nsurv = Int[], theta = Float64[], omega = Float64[], mean_fitness = Float64[],
         sel_diff = Float64[], popsize = Int[]) :
        (gen = Int[], deme = Int[], pheno_mean = Float64[], pheno_var = Float64[],
         geno_var = Float64[], nsurv = Int[], theta = Float64[], omega = Float64[],
         mean_fitness = Float64[], sel_diff = Float64[], popsize = Int[],
         demeof = IdDict{World,Int}(w => i for (i, w) in enumerate(demes)))
    traits == 1 ? rec : merge(rec, trait_records())
end

# the columns of a model of several traits: per-trait means, the G and P matrices, the
# lag θ − z̄ and the selection differential vector, one vector or matrix per record
trait_records() = (trait_mean = Vector{Float64}[], gmatrix = Matrix{Float64}[], pmatrix = Matrix{Float64}[],
                   trait_lag = Vector{Float64}[], trait_sel_diff = Vector{Float64}[])

"""
    record!(rec, w)

Append the parent generation's statistics to the columns `rec` (see [`records`](@ref))
every `Sample` generations (QTSim `Output`); a no-op when `rec` is `nothing`. Variances
are sample variances (n − 1), as GSL's. Besides QTSim's columns it records the mean
fitness `W̄` (one minus the load) and the selection differential, the survivors' mean
phenotype minus the mean phenotype, and the population size `N` (`PopSize`, unless the
population is density-regulated). After fecundity selection, when all survive, the
selection differential is that of the mating weights, `Σ W z/ΣW − z̄`.
"""
function record!(rec, w::World)
    (rec === nothing || w.gen % w.p.Sample != 0) && return w
    push!(rec.gen, w.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[w])
    z = view(w.phenotype, 1:w.N)
    push!(rec.pheno_mean, mean(z))
    push!(rec.pheno_var, var(z)); push!(rec.geno_var, var(view(w.genotype, 1:w.N)))
    push!(rec.nsurv, w.nsurv); push!(rec.theta, w.θ); push!(rec.omega, w.ω)
    push!(rec.mean_fitness, mean(view(w.fitness, 1:w.N)))
    push!(rec.sel_diff, (w.weighted ? sum(w.fitness[i] * z[i] for i in 1:w.N) / w.wsum :
                         mean(w.phenotype[i] for i in w.survivors)) - mean(z))
    push!(rec.popsize, w.N)
    w.traits === nothing || record_traits!(rec, w)
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
    rec = records(; traits = p.NoTraits)
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
