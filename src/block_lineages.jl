# Haplotype-block lineages: a reduced representation of the kernel for the emulator.
#
# Each chromosome is coarsened into nb blocks of consecutive loci, and an individual carries
# only the sum of the allelic effects of each block of each of its chromosomes (no per-locus
# state). The life cycle is the kernel's, on N individuals: phenotype, viability selection
# toward the moving optimum, random mating of two distinct survivors, gametes with QTSim's
# crossover, Poisson–Gaussian mutation. Crossovers happen only at block boundaries: the
# kernel's crossover point is drawn from its exact law and moved to the nearer boundary of
# the block it falls in (nearer on the genetic map). Blocks are therefore non-recombining
# units, and for them the sum is a sufficient state: a mutation at any locus of a block adds
# its effect to the block's sum, whatever the alleles. So the block model is *exactly* the
# kernel with that coarsened crossover map, and with one block per locus it is the kernel
# itself (in law; the random draws differ). The only approximation is the map, which the
# number of blocks controls; associations between alleles and their linked backgrounds
# (Hill–Robertson, Barton 1995) arise by construction, not by a closure.
#
# Validation: validate/blocks_dev.jl (design, convergence in the number of blocks),
# validate/blocks_criteria.md and validate/blocks.jl (pre-specified, fresh inputs).

"""
    QTLSim.BlockMap

The block partition of a chromosome and the law of the switch point of a gamete on it
(see [`QTLSim.block_map`](@ref)). Fields:

- `L`, `nb`: loci and blocks;
- `hi`: the last locus (1-based) of each block, `hi[nb] = L`;
- `blockof`: the block of each locus;
- `nloci`: the loci of each block;
- `pnone`: the probability that a gamete is a copy of one parental chromosome (no
  crossover, or one moved to an end of the chromosome);
- `cum`: for `:qtsim`, the cumulative probabilities of a switch after block `j`,
  `j = 1:nb-1`, beyond `pnone` (`cum[end] = 1`): a uniform `u ≥ pnone` switches after block
  `searchsortedfirst(cum, u)`;
- `snap`: the boundary (`0:nb`) that each kernel crossover point `c ∈ 0:L-2` (a switch after
  0-based locus `c`) is moved to; `0` and `nb` are the ends (no switch);
- `model`: the `CrossoverModel`, `:qtsim` (QTSim's lowest point only) or `:all`;
- `recomb`: `RecombPoissonMean`;
- `kernel`: the kernel's law of the switch point under `:qtsim`, `P(c)` for `c ∈ 0:L-2`
  (`kernel[c+1]`), and the block map's `bmass`, the probability of a switch after each block.
"""
struct BlockMap
    L::Int
    nb::Int
    hi::Vector{Int}
    blockof::Vector{Int32}
    nloci::Vector{Int}
    pnone::Float64
    cum::Vector{Float64}
    snap::Vector{Int}
    model::Symbol
    recomb::Float64
    kernel::Vector{Float64}
    bmass::Vector{Float64}
    within::Symbol
    pnk::Float64
    gcum::Vector{Float64}
    gblock::Vector{Int32}
    gf::Vector{Float64}
end

"""
    QTLSim.qtsim_switch_law(L, u) -> (; P, pnone)

The law of the switch point of a gamete under QTSim's crossover (`CrossoverModel = :qtsim`,
`RecombPoissonMean = u`, `L` loci): `K ~ Poisson(u)` distinct points drawn uniformly from
the `M = L − 1` gaps `0:L-2`, none if `K = 0` or `K ≥ L`, and only the lowest point
switches. For `1 ≤ K ≤ M`, `P(min = c | K) = C(M−c−1, K−1)/C(M, K)`, so

    P(c) = Σ_{K=1}^{M} Pois(K; u) · C(M−c−1, K−1)/C(M, K),   c = 0, …, M−1,

and `pnone = P(K = 0) + P(K ≥ L)`. `P[c+1]` is `P(c)`. Checked symbolically and against
draws of the kernel's `gamete!` in validate/symbolic/check_blocks.jl.
"""
function qtsim_switch_law(L::Integer, u::Real)
    L >= 1 || throw(ArgumentError("need at least one locus, got $L"))
    M = L - 1
    P = zeros(M)
    M == 0 && return (; P, pnone = 1.0)
    pk = exp(-u)                          # Pois(0)
    pnone = pk
    Kmax = M
    tailsum = pk
    for K in 1:Kmax
        pk *= u / K                       # Pois(K)
        tailsum += pk
        pk < 1e-300 && break
        # r_c = C(M−c−1, K−1)/C(M, K): r_0 = K/M, r_{c+1} = r_c (M−c−K)/(M−c−1)
        r = K / M
        for c in 0:M-K
            P[c+1] += pk * r
            c + 1 <= M - K && (r *= (M - c - K) / (M - c - 1))
        end
    end
    pnone += max(0.0, 1 - tailsum)        # P(K ≥ L), up to rounding
    (; P, pnone)
end

"""
    QTLSim.block_map(p::Params, nb; partition = :map) -> BlockMap

Split the `NoQTL` loci of a chromosome into `min(nb, NoQTL)` blocks of consecutive loci and
move each crossover point of the kernel to a block boundary.

- `partition = :map` (default): boundaries at equal steps of the kernel's genetic map, the
  cumulative probability `F(c)` that a gamete switches at or before gap `c` (under
  `:qtsim`, the law [`QTLSim.qtsim_switch_law`](@ref); under `:all`, `c + 1` gaps of
  `RecombPoissonMean/(L−1)` expected switches each). Block `j` ends at the first gap with
  `F ≥ j·F(L−2)/nb`, so every block carries about `1/nb` of the switches: the error in the
  probability that two loci are split is at most about one block's share of the map,
  whatever their positions. Under QTSim's lowest-point rule the map is dense at the start
  of the chromosome and sparse at its end, so the blocks at the end hold many loci.
- `partition = :position`: blocks of equal numbers of loci.

A crossover point inside a block is moved to the block's boundary that is nearer on the
map (the lower one on a tie). Boundaries `0` and `nb` are the chromosome's ends, where a
switch is no switch (the starting strand is a fair coin). With `nb ≥ NoQTL` every locus is
a block and the map is the kernel's.

`within` is what a gamete does with a switch point inside a block: `:snap` (default) moves it
as above, so blocks never recombine inside; `:gaussian` keeps the kernel's switch point and
draws the recombinant block by Gaussian segregation, `ȳ_b + f(y_A − ȳ_b) + (1 − f)(y_B − ȳ_b) +
√(2f(1 − f)V_b) ξ`, `f` the fraction of the block's loci at or before the point, `ȳ_b` and `V_b`
the block's mean and variance over the parents' generation: exact in law if the two parts
of the block are Gaussian, independent and carry variance in proportion to their loci
(`:qtsim` only). Validation of the choice: validate/blocks_dev.log.
"""
function block_map(p::Params, nb::Integer; partition::Symbol = :map, within::Symbol = :snap)
    within in (:snap, :gaussian) || throw(ArgumentError("within is :snap or :gaussian, got :$within"))
    within === :gaussian && p.CrossoverModel !== :qtsim &&
        throw(ArgumentError("within = :gaussian is for CrossoverModel = :qtsim (one switch per gamete)"))
    L = p.NoQTL
    nb >= 1 || throw(ArgumentError("need at least one block, got $nb"))
    partition in (:map, :position) || throw(ArgumentError("partition is :map or :position, got :$partition"))
    model = p.CrossoverModel
    u = p.RecombPoissonMean
    M = L - 1
    law = model === :qtsim ? qtsim_switch_law(L, u) : (; P = fill(M > 0 ? u / M : 0.0, M), pnone = exp(-u))
    P = law.P
    F = cumsum(P)                          # F[c+1]: at or before gap c
    nb = min(nb, L)
    # block ends (0-based gap after the block's last locus); the last block ends at L-1
    ends = Int[]
    if nb > 1
        if partition === :position || M == 0 || F[end] <= 0 || nb >= L
            for j in 1:nb-1
                push!(ends, round(Int, j * L / nb) - 1)
            end
        else
            # the first gap with F ≥ j·Ftot/nb, but at least one locus per block: where one
            # gap carries more than a block's share of the map (the dense start of QTSim's
            # map at few loci) the blocks are single loci, and the rest follow the map
            Ftot = F[end]
            prev = -1
            for j in 1:nb-1
                c = searchsortedfirst(F, j * Ftot / nb) - 1
                c = clamp(c, prev + 1, M - (nb - j))
                push!(ends, c)
                prev = c
            end
        end
    end
    nbb = length(ends) + 1
    hi = [ends .+ 1; L]                    # last locus (1-based) of each block
    blockof = Vector{Int32}(undef, L)
    lo = 1
    for (b, h) in enumerate(hi)
        blockof[lo:h] .= b
        lo = h + 1
    end
    nloci = diff([0; hi])
    # snap each gap c ∈ 0:M-1 to a boundary 0:nbb; boundary j (1 ≤ j < nbb) is gap ends[j]
    snap = zeros(Int, M)
    Fat(c) = c < 0 ? 0.0 : F[c+1]          # F at gap c, 0 before the first
    for c in 0:M-1
        b = blockof[c+1]                   # the block of locus c+1 (1-based), whose last gap…
        if b < nbb && c == ends[b]
            snap[c+1] = b                  # on a boundary
            continue
        end
        # inside block b: its lower boundary is b−1 (gap ends[b−1], or the start), its upper b
        flo = b == 1 ? 0.0 : Fat(ends[b-1])
        fhi = b == nbb ? (M > 0 ? F[end] : 0.0) : Fat(ends[b])
        fc = Fat(c) - P[c+1] / 2           # the point's own mass at its middle
        snap[c+1] = (fc - flo <= fhi - fc) ? b - 1 : b
    end
    bmass = zeros(nbb - 1)
    pnone = law.pnone
    for c in 0:M-1
        j = snap[c+1]
        (j == 0 || j == nbb) ? (pnone += P[c+1]) : (bmass[j] += P[c+1])
    end
    rest = 1 - pnone
    cum = rest > 0 ? cumsum(bmass) ./ rest : Float64[]
    isempty(cum) || (cum[end] = 1.0)
    # within = :gaussian: the kernel's switch point itself, its block and the fraction of the
    # block's loci at or before it (1 on a boundary)
    gk = 1 - law.pnone
    gcum = gk > 0 ? cumsum(P) ./ gk : Float64[]
    isempty(gcum) || (gcum[end] = 1.0)
    lo = [1; hi[1:end-1] .+ 1]
    gblock = Int32[blockof[c+1] for c in 0:M-1]
    gf = [(c + 1 - lo[blockof[c+1]] + 1) / nloci[blockof[c+1]] for c in 0:M-1]
    BlockMap(L, nbb, hi, blockof, nloci, pnone, cum, snap, model, u, P, bmass, within, law.pnone, gcum, gblock, gf)
end

"""
    QTLSim.BlockPop

The state of one population of the block model (see [`run_blocks`](@ref)): block sums
`Y[b, 2i−1]` and `Y[b, 2i]` of individual `i`'s two chromosomes, the offspring buffer `Yo`,
and the per-generation scratch.
"""
mutable struct BlockPop
    Y::Matrix{Float64}
    Yo::Matrix{Float64}
    g::Vector{Float64}
    z::Vector{Float64}
    W::Vector{Float64}
    surv::Vector{Int}
    N::Int
    θ::Float64
    ω::Float64
    gen::Int
    extinct::Bool
    ybar::Vector{Float64}
    vb::Vector{Float64}
end

BlockPop(nb::Integer, N::Integer, ω::Real) =
    BlockPop(zeros(nb, 2N), zeros(nb, 2N), zeros(N), zeros(N), zeros(N), sizehint!(Int[], N), N, 0.0, float(ω), 0, false, zeros(nb), zeros(nb))

Base.copy(s::BlockPop) = BlockPop(copy(s.Y), copy(s.Yo), copy(s.g), copy(s.z), copy(s.W), copy(s.surv), s.N, s.θ, s.ω, s.gen, s.extinct, copy(s.ybar), copy(s.vb))

function check_blocks_params(p::Params)
    bad(msg) = throw(ArgumentError("the block model " * msg))
    p.Ploidy == 2 || bad("is diploid; Ploidy = $(p.Ploidy)")
    p.NoTraits == 1 || bad("has one trait; NoTraits = $(p.NoTraits)")
    (p.MoveFitnessSDGen == -1 || p.MoveType === :continuous) || bad("moves the optimum continuously (MoveType = :continuous); got :$(p.MoveType)")
    p.Selection === :viability || bad("has viability selection; Selection = :$(p.Selection)")
    p.Demography === :fixed || bad("has a fixed population size (its start is the kernel's equilibrium sampler); Demography = :$(p.Demography)")
    p.Plasticity == 0 || bad("has no plasticity")
    p.StandardisePhenotype && bad("does not standardise phenotypes")
    nothing
end

"""
    QTLSim.block_equilibrium!(s, p, bm, rng; spectrum = stationary_spectrum(p), target)

Fill the block sums of `s` from the kernel's equilibrium sampler ([`equilibrium_genomes`](@ref)),
summed over each block: the same law of genomes as `equilibrium_world`, without building the
dense `(L, 2, N)` array. `target` is the genic variance the effects are scaled to (as
`equilibrium_genomes`'s numeric `target`).
"""
function block_equilibrium!(s::BlockPop, p::Params, bm::BlockMap, rng::AbstractRNG;
                            spectrum = stationary_spectrum(p), target::Real = init_target(p, spectrum, :genic))
    st = spectrum
    N, L, nb = p.PopSize, p.NoQTL, bm.nb
    Y = s.Y
    fill!(Y, 0.0)
    if st.VG > 0
        κ = sqrt(target / st.VG)
        m = 2N
        cum = cumsum(1 / k for k in 1:m-1)
        sd = sqrt(st.α2)
        idx = collect(1:m)
        for _ in 1:rand(rng, Poisson(st.theta * cum[end]))
            k = min(searchsortedfirst(cum, rand(rng) * cum[end]), m - 1)
            x = rand(rng, Beta(k, m - k + 1))
            a = sd * randn(rng)
            rand(rng) < sojourn_ratio(spectrum_gamma(st, a^2), x) || continue
            b = bm.blockof[rand(rng, 1:L)]
            e = κ * a
            if 2k > m
                @inbounds for h in 1:m
                    Y[b, h] += e
                end
                k, e = m - k, -e
            end
            @inbounds for j in 1:k
                r = rand(rng, j:m)
                idx[j], idx[r] = idx[r], idx[j]
                Y[b, idx[j]] += e
            end
        end
        z = sqrt(stationary_mean_var(p, st; VG = target)) * randn(rng)
        shift = (z - sum(Y) / N) / (2L)
        @inbounds for h in 1:m, b in 1:nb
            Y[b, h] += shift * bm.nloci[b]
        end
    end
    s.N = N
    s.extinct = false
    s
end

# one generation of the block model (develop, select, reproduce, rescale, mutate, swap);
# returns the statistics of the parent generation
function block_generation!(s::BlockPop, p::Params, bm::BlockMap, rng::AbstractRNG, esd::Float64, mutsd::Float64)
    N, nb = s.N, bm.nb
    Y = s.Y
    # develop: genotype and phenotype
    gs = 0.0; gss = 0.0; zs = 0.0; zss = 0.0
    @inbounds for i in 1:N
        g = 0.0
        c1 = 2i - 1; c2 = 2i
        @simd for b in 1:nb
            g += Y[b, c1] + Y[b, c2]
        end
        s.g[i] = g
        z = g + esd * randn(rng)
        s.z[i] = z
    end
    # select: move the optimum, fitness, survivors
    if p.MoveFitnessSDGen != -1 && s.gen >= p.MoveFitnessSDGen
        s.θ += p.MoveFitnessSD
    end
    ω = s.ω
    c0 = 1 / (0.4 * sqrt(2π))
    empty!(s.surv)
    Wsum = 0.0
    @inbounds for i in 1:N
        x = s.z[i] - s.θ
        w = p.Drift ? 1.0 : c0 * exp(-x * x / (2ω * ω))
        s.W[i] = w
        Wsum += w
        rand(rng) < w && push!(s.surv, i)
    end
    ns = length(s.surv)
    gbar = sum(view(s.g, 1:N)) / N
    zbar = sum(view(s.z, 1:N)) / N
    VG = N > 1 ? sum(abs2, view(s.g, 1:N) .- gbar) / (N - 1) : 0.0
    VP = N > 1 ? sum(abs2, view(s.z, 1:N) .- zbar) / (N - 1) : 0.0
    stats = (; gen = s.gen, N, gbar, zbar, VG, VP, θ = s.θ, Wbar = Wsum / N, nsurv = ns)
    if ns < 2
        s.extinct = true
        return stats
    end
    # reproduce
    Noff = p.Demography === :fixed ? p.PopSize :
           min(p.PopSize, rand(rng, Poisson(p.GrowthRate * ns / (1 + (p.GrowthRate - 1) * ns / p.PopSize))))
    if Noff == 0
        s.extinct = true
        return stats
    end
    Yo = s.Yo
    allpts = bm.model === :all
    if bm.within === :gaussian            # the parents' block means and variances (all adults), column by column
        M2 = 2N
        ym, yv = s.ybar, s.vb
        fill!(ym, 0.0); fill!(yv, 0.0)
        @inbounds for h in 1:M2
            @simd for b in 1:nb
                ym[b] += Y[b, h]
            end
        end
        ym ./= M2
        @inbounds for h in 1:M2
            @simd for b in 1:nb
                d = Y[b, h] - ym[b]; yv[b] += d * d
            end
        end
        yv ./= M2
    end
    @inbounds for j in 1:Noff
        p1 = s.surv[rand(rng, 1:ns)]
        p2 = s.surv[rand(rng, 1:ns)]
        while p2 == p1
            p2 = s.surv[rand(rng, 1:ns)]
        end
        if bm.within === :gaussian
            _block_gamete_gauss!(Yo, 2j - 1, Y, p1, bm, rng, s)
            _block_gamete_gauss!(Yo, 2j, Y, p2, bm, rng, s)
        else
            _block_gamete!(Yo, 2j - 1, Y, p1, bm, rng, allpts)
            _block_gamete!(Yo, 2j, Y, p2, bm, rng, allpts)
        end
    end
    # rescale (takes effect next generation)
    s.gen == p.ChangeFitnessGen && (s.ω = p.ChangeFitnessSD)
    # mutate: Poisson(2·Noff·L·μ) mutations, each on a uniform chromosome and locus (the
    # same law as Poisson(L·μ) per chromosome)
    λ = 2Noff * p.NoQTL * p.MutationRatePerGamete
    if λ > 0 && mutsd > 0
        for _ in 1:rand(rng, Poisson(λ))
            h = rand(rng, 1:2Noff)
            b = bm.blockof[rand(rng, 1:p.NoQTL)]
            Yo[b, h] += mutsd * randn(rng)
        end
    end
    # swap
    s.Y, s.Yo = s.Yo, s.Y                 # Noff ≤ PopSize: the buffers hold 2·PopSize chromosomes
    s.N = Noff
    s.gen += 1
    stats
end

# within = :gaussian: the kernel's switch point; the blocks before its block from the starting
# strand, after it from the other, and its own block (unless the point is on its boundary) by
# Gaussian segregation: ȳ + f d_A + (1 − f) d_B + √(2f(1 − f)V_b) ξ
@inline function _block_gamete_gauss!(Yo, dst, Y, par, bm::BlockMap, rng, s::BlockPop)
    nb = bm.nb
    a, o = rand(rng, Bool) ? (2par - 1, 2par) : (2par, 2par - 1)
    u = rand(rng)
    if u < bm.pnk || isempty(bm.gcum)
        @inbounds @simd for b in 1:nb
            Yo[b, dst] = Y[b, a]
        end
        return nothing
    end
    c = min(searchsortedfirst(bm.gcum, (u - bm.pnk) / (1 - bm.pnk)), length(bm.gcum))
    k = Int(bm.gblock[c]); f = bm.gf[c]
    @inbounds for b in 1:k-1
        Yo[b, dst] = Y[b, a]
    end
    @inbounds if f >= 1
        Yo[k, dst] = Y[k, a]
    else
        m = s.ybar[k]
        Yo[k, dst] = m + f * (Y[k, a] - m) + (1 - f) * (Y[k, o] - m) + sqrt(2f * (1 - f) * s.vb[k]) * randn(rng)
    end
    @inbounds for b in k+1:nb
        Yo[b, dst] = Y[b, o]
    end
    nothing
end

# a gamete of parent `par` into column `dst` of Yo
@inline function _block_gamete!(Yo, dst, Y, par, bm::BlockMap, rng, allpts::Bool)
    nb = bm.nb
    a, o = rand(rng, Bool) ? (2par - 1, 2par) : (2par, 2par - 1)
    if !allpts
        u = rand(rng)
        if u < bm.pnone || nb == 1
            @inbounds @simd for b in 1:nb
                Yo[b, dst] = Y[b, a]
            end
            return nothing
        end
        j = searchsortedfirst(bm.cum, (u - bm.pnone) / (1 - bm.pnone))
        j = min(j, nb - 1)
        @inbounds for b in 1:j
            Yo[b, dst] = Y[b, a]
        end
        @inbounds for b in j+1:nb
            Yo[b, dst] = Y[b, o]
        end
        return nothing
    end
    # :all: the kernel's points, each moved to its boundary, every one switching
    L = bm.L
    K = rand(rng, Poisson(bm.recomb))
    if K == 0 || K >= L
        @inbounds for b in 1:nb
            Yo[b, dst] = Y[b, a]
        end
        return nothing
    end
    # distinct points as the kernel draws them, each moved to its boundary; a boundary
    # switches the strand if it received an odd number of points (boundaries 0 and nb are
    # the ends: a switch before the first block flips the whole gamete, which the fair
    # starting strand absorbs in law, and one after the last changes nothing)
    pts = Int[]
    while length(pts) < K
        c = rand(rng, 0:L-2)
        c in pts || push!(pts, c)
    end
    par = falses(nb)
    for c in pts
        j = bm.snap[c+1]
        1 <= j < nb && (par[j] = !par[j])
    end
    cur, oth = a, o
    @inbounds for b in 1:nb
        Yo[b, dst] = Y[b, cur]
        par[b] && ((cur, oth) = (oth, cur))
    end
    nothing
end

"""
    run_blocks(p::Params, rng; blocks = 40, partition = :map, within = :snap, burnin, at = QTLSim.first_event(p),
               record_from = nothing, spectrum = stationary_spectrum(p), map = nothing) -> NamedTuple

One run of the block model (src/block_lineages.jl): the kernel's life cycle on `PopSize`
diploid individuals, each chromosome reduced to the sums of its `blocks` blocks
([`QTLSim.block_map`](@ref)). It starts `burnin`
generations before generation `at` of the protocol (as [`equilibrium_world`](@ref)), from the
kernel's equilibrium sampler summed over the blocks ([`QTLSim.block_equilibrium!`](@ref)),
and runs to `NoGenerations − 1`, or until fewer than two individuals survive selection.

Returns `(; extinct, extinct_gen, start, gen, VG, VP, gbar, zbar, theta, Wbar, nsurv)`, the
statistics of the parent generation for every generation from `record_from` (default: the
start) on. Supported: one diploid trait, viability selection, `MoveType = :continuous` (or
no move), a fixed population size, `CrossoverModel` `:qtsim` or `:all`;
anything else is an `ArgumentError`.

With `blocks ≥ NoQTL` the map is the kernel's, so the run has the kernel's law. With fewer
blocks and `within = :snap` (default) only the map is coarsened (a block's sum is a sufficient
state while nothing recombines inside it), but V_G then converges slowly in `blocks`;
`within = :gaussian` keeps every switch point and closes the block that holds it by Gaussian
segregation ([`QTLSim.block_map`](@ref)). On the development inputs `within = :gaussian` with
160 blocks met the convergence rule: validate/blocks_dev.log, validate/blocks_criteria.md.
"""
function run_blocks(p::Params, rng::AbstractRNG; blocks::Integer = 40, partition::Symbol = :map, within::Symbol = :snap, burnin::Integer,
                    at::Integer = first_event(p), record_from = nothing, spectrum = nothing, map = nothing,
                    target = nothing, observe = nothing)
    checkparams(p); check_blocks_params(p)
    bm = map === nothing ? block_map(p, blocks; partition, within) : map
    st = spectrum === nothing ? stationary_spectrum(p) : spectrum
    tg = target === nothing ? init_target(p, st, :genic) : float(target)
    s = BlockPop(bm.nb, p.PopSize, initFitSD(p))
    block_equilibrium!(s, p, bm, rng; spectrum = st, target = tg)
    s.gen = protocol_start(p; burnin, at)
    _run_blocks!(s, p, bm, rng; record_from = record_from === nothing ? s.gen : record_from, observe)
end

function _run_blocks!(s::BlockPop, p::Params, bm::BlockMap, rng::AbstractRNG; record_from::Integer = s.gen,
                      observe = nothing, stop::Integer = p.NoGenerations)
    esd = envsd(p)
    L, μ = p.NoQTL, p.MutationRatePerGamete
    mutsd = μ > 0 ? sqrt(max((p.MutationScalar / (p.Ploidy * L * μ)) * envVar(p), 0.0)) : 0.0
    start = s.gen
    gen = Int[]; VG = Float64[]; VP = Float64[]; gbar = Float64[]; zbar = Float64[]
    theta = Float64[]; Wbar = Float64[]; nsurv = Int[]
    extinct_gen = -1
    while s.gen < stop && !s.extinct
        observe === nothing || observe(s, bm)
        st = block_generation!(s, p, bm, rng, esd, mutsd)
        if st.gen >= record_from
            push!(gen, st.gen); push!(VG, st.VG); push!(VP, st.VP); push!(gbar, st.gbar); push!(zbar, st.zbar)
            push!(theta, st.θ); push!(Wbar, st.Wbar); push!(nsurv, st.nsurv)
        end
        s.extinct && (extinct_gen = st.gen)
    end
    (; extinct = s.extinct, extinct_gen, start, gen, VG, VP, gbar, zbar, theta, Wbar, nsurv)
end

"""
    QTLSim.block_stats(s::BlockPop, bm) -> (; blockgenic, assoc, VG)

Block-level statistics of the current adults (population denominators over the `2N`
chromosomes): `blockgenic = 2 Σ_b Var(y_b)`, the variance of the block sums, which is the
genic variance when every block is one locus and otherwise also holds the linkage
disequilibrium within blocks; `assoc = Σ_b Cov(G − y_b, (y_b − ȳ_b)²)`, the association of
a block's squared deviation with the rest of its chromosome (`G` the chromosome's total),
the block-level counterpart of `Σ_l Cov(g_rest, d_l²)` of validate/genic_decline.md.
"""
function block_stats(Y::AbstractMatrix, nb::Integer, M::Integer)
    ybar = vec(sum(view(Y, :, 1:M); dims = 2)) ./ M
    G = vec(sum(view(Y, :, 1:M); dims = 1))
    Gbar = sum(G) / M
    bg = 0.0; assoc = 0.0
    @inbounds for b in 1:nb
        v = 0.0; c = 0.0; d2m = 0.0
        for h in 1:M
            d = Y[b, h] - ybar[b]
            v += d * d
        end
        d2m = v / M
        for h in 1:M
            d = Y[b, h] - ybar[b]
            c += (G[h] - Gbar - d) * (d * d - d2m)
        end
        bg += 2v / M
        assoc += c / M
    end
    (; blockgenic = bg, assoc)
end
block_stats(s::BlockPop, bm::BlockMap) = block_stats(s.Y, bm.nb, 2s.N)

"""
    QTLSim.blocks_emulate(p; blocks, replicates, seed, rng, burnin, horizon, tail, partition)

`emulate(p; loci = :blocks, …)`: see [`emulate`](@ref).
"""
function blocks_emulate(p::Params; blocks::Integer = 40, replicates::Integer = 40, seed = nothing,
                        rng::Union{Nothing,AbstractRNG} = nothing, burnin::Integer = 2000, horizon::Integer = 2000,
                        tail::Integer = 240_000, partition::Symbol = :map, within::Symbol = :snap)
    replicates >= 1 || throw(ArgumentError("the block model needs replicates ≥ 1, got $replicates"))
    seed === nothing && rng === nothing && throw(ArgumentError("the block model needs a `seed` or an `rng`"))
    q = Params(p; Sample = 1, NoGenerations = burnin + horizon, MoveFitnessSDGen = burnin, ChangeFitnessGen = -1)
    checkparams(q); check_blocks_params(q)
    bm = block_map(q, blocks; partition, within)
    st = stationary_spectrum(q)
    tg = init_target(q, st, :genic)
    seeds = rng === nothing ? [job_seed(seed, r) for r in 1:replicates] : [rand(rng, UInt64) for _ in 1:replicates]
    runs = [run_blocks(q, Xoshiro(sd); map = bm, burnin, at = burnin, spectrum = st, target = tg, record_from = burnin)
            for sd in seeds]
    _blocks_summary(runs, burnin, horizon, tail)
end

# the emulator's result from the replicate runs (each recorded from the move on)
function _blocks_summary(runs, burnin, horizon, tail)
    R = length(runs)
    last = burnin + horizon - 1
    # T counted as the table: the first generation of the move is 1; a run that dies in
    # generation g (its parents' generation g) has T = g − burnin + 1
    T = [r.extinct ? r.extinct_gen - burnin + 1 : typemax(Int) for r in runs]
    alive = count(>(horizon), T)
    P2000 = alive / R
    sT = sort(T)
    med = sT[cld(R, 2)] <= horizon && (isodd(R) || sT[R÷2+1] <= horizon) ? (isodd(R) ? sT[cld(R, 2)] : (sT[R÷2] + sT[R÷2+1]) / 2) : missing
    # end hazard: deaths in the last quarter of the horizon per run-generation at risk
    w0 = horizon - horizon ÷ 4
    atrisk = 0; deaths = 0
    for t in T
        t > w0 || continue
        deaths += t <= horizon
        atrisk += min(t, horizon) - w0
    end
    hazard = atrisk > 0 ? deaths / atrisk : 1.0
    P240k = P2000 == 0 ? 0.0 : P2000 * (1 - hazard)^(tail - horizon)
    at(r, f, g) = (k = searchsortedfirst(r.gen, g); k <= length(r.gen) && r.gen[k] == g ? getfield(r, f)[k] : NaN)
    mn(v) = (u = filter(isfinite, v); isempty(u) ? NaN : sum(u) / length(u))
    VG_move = mn([at(r, :VG, burnin) for r in runs])
    Wbar_move = mn([at(r, :Wbar, burnin) for r in runs])
    VG_end = mn([at(r, :VG, last) for r in runs])
    lag_end = mn([at(r, :theta, last) - at(r, :zbar, last) for r in runs])
    (; P2000, T = med, P240k, hazard, VG_move, Wbar_move, VG_end, lag_end, extinct = alive == 0 && all(r -> r.extinct, runs),
     alive, replicates = R)
end

# emulated_half_survival_rate(p; loci = :blocks): the replicates are burned in once (the
# burn-in does not depend on the rate), each rate continues copies of them with copies of
# their generators (common random numbers), and P2000(k) ≥ 1/2 is bisected in log k.
function blocks_half_survival_rate(p::Params; blocks::Integer = 40, replicates::Integer = 40, seed = nothing,
                                   burnin::Integer = 2000, horizon::Integer = 2000, partition::Symbol = :map, within::Symbol = :snap,
                                   rtol::Real = 1e-3, maxiter::Integer = 40, bracket = nothing)
    replicates >= 1 || throw(ArgumentError("the block model needs replicates ≥ 1, got $replicates"))
    seed === nothing && throw(ArgumentError("the block model needs a `seed`"))
    q = Params(p; Sample = 1, NoGenerations = burnin + horizon, MoveFitnessSDGen = burnin, ChangeFitnessGen = -1)
    checkparams(q); check_blocks_params(q)
    bm = block_map(q, blocks; partition, within)
    st = stationary_spectrum(q)
    tg = init_target(q, st, :genic)
    base = map(1:replicates) do r
        rng = Xoshiro(job_seed(seed, r))
        s = BlockPop(bm.nb, q.PopSize, initFitSD(q))
        block_equilibrium!(s, q, bm, rng; spectrum = st, target = tg)
        s.gen = protocol_start(q; burnin, at = burnin)
        _run_blocks!(s, q, bm, rng; stop = burnin)
        (s, rng)
    end
    evals = 0
    function P(k)
        evals += 1
        qk = Params(q; MoveFitnessSD = k)
        alive = 0
        for (s0, r0) in base
            s0.extinct && continue
            s = copy(s0); rng = copy(r0)
            _run_blocks!(s, qk, bm, rng; record_from = typemax(Int))
            s.extinct || (alive += 1)
        end
        alive / replicates
    end
    bisect(lo, hi) = begin
        for _ in 1:maxiter
            hi - lo <= rtol && break
            mid = (lo + hi) / 2
            P(exp(mid)) >= 0.5 ? (lo = mid) : (hi = mid)
        end
        exp((lo + hi) / 2)
    end
    if bracket !== nothing
        # a narrow bracket first (the two-stage emulator's band around the screen's k50): if it
        # misses the root it is moved toward it, each time four times wider in log k, up to four
        # times, before falling back to the wide search below (`bracketed = false`)
        a, b = log(bracket[1]), log(bracket[2])
        a < b || throw(ArgumentError("bracket must be (lo, hi) with 0 < lo < hi"))
        pa, pb = P(exp(a)) >= 0.5, P(exp(b)) >= 0.5
        for _ in 1:4
            (pa && !pb) && break
            w = 4(b - a)
            pb ? ((a, pa) = (b, pb); b = a + w; pb = P(exp(b)) >= 0.5) : ((b, pb) = (a, pa); a = b - w; pa = P(exp(a)) >= 0.5)
        end
        pa && !pb && return (; k50 = bisect(a, b), kc = NaN, evals, bracketed = true)
    end
    kc = drift_critical_rate(p).k
    bracket === nothing || return merge(_blocks_wide(P, kc, bisect), (; evals, bracketed = false))
    merge(_blocks_wide(P, kc, bisect), (; evals))
end

# the wide search of blocks_half_survival_rate: kc/8 to 2kc, the top doubled up to 64kc, the
# bottom halved down to kc/512
function _blocks_wide(P, kc, bisect)
    kc > 0 || return (; k50 = 0.0, kc)
    lo, hi = log(kc / 8), log(2kc)
    for _ in 1:5
        P(exp(hi)) >= 0.5 || break
        lo = hi; hi += log(2.0)
        hi > log(64kc) && return (; k50 = Inf, kc)
    end
    P(exp(hi)) >= 0.5 && return (; k50 = Inf, kc)
    for _ in 1:6
        P(exp(lo)) >= 0.5 && break
        hi = lo; lo -= log(2.0)
    end
    P(exp(lo)) >= 0.5 || return (; k50 = 0.0, kc)
    (; k50 = bisect(lo, hi), kc)
end
