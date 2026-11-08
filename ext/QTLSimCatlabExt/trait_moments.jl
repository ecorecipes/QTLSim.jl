# The moments lowering for several traits: a second table of images for the same steps
# of the same life cycle, on the moments of a multivariate Gaussian population (Lande
# 1979, 1980). It is kept apart from moments.jl, the one-trait table, to show that only
# the images differ: `run_trait_moments` lowers the unchanged `lifecycle()` with the table
# of `TraitModel(n)`, `trait_moment_images(rec)`. With one trait its recursion is
# moments.jl's; the fixed points are QTLSim's `gaussian_G_loci` (src/multilocus.jl;
# `gaussian_G` with `loci = :aggregate`) and `steady_lag_vector` (src/traits_theory.jl).
# By default it carries the haplotype moments of the multilocus model as matrices
# (`TraitCov`, the equations of src/multilocus.jl, shared with the one-trait table through
# `QTLSim.hap_select` and its siblings); with `loci = :aggregate` selection's reduction
# `G (Ω + P)⁻¹ G` is inherited whole.
#
# Several demes (k of them, `run_trait_moments(ps; cycle)`) are the one-trait k-deme
# moments (moments.jl, `DemeCov`) with every moment a matrix: the haplotype moments of
# each deme and, with drift (`between = true` or `erosion = true`), the cross-moments
# Q_ab of the demes, `E[(h − z̄_a/2)(h′ − z̄_b/2)ᵀ]`. Selection in deme a shrinks the
# deviation of its mean by `I − B_a`, `B_a = G_a s_a⁻¹`, so `Q_ab ← (I − B_a) Q_ab`; migration
# mixes means, haplotype moments and cross-moments with the outer products of the
# differences of the means in place of their squares.

"""
The state of the trait moments lowering: as `Moments`, with vectors for the means and the
optimum and matrices for the variances. `M` is the mutational covariance, `shape` the
fitness function's shape (`Ω = ω² S`), `u` the direction the optimum moves in, `θ0` the
optimum when a fluctuating move starts, `noise` the covariance of its AR(1) noise and `f`
the direction it cycles in; `cov` the haplotype moments of all the demes (`TraitCov`) and
`idx` this deme's index in them.
"""
mutable struct TraitMoments
    p::Params
    gen::Int
    zbar::Vector{Float64}
    G::Matrix{Float64}
    P::Matrix{Float64}
    θ::Vector{Float64}
    ω::Float64
    nsurv::Float64
    Wbar::Float64
    S::Vector{Float64}      # selection differential on the phenotype
    zs::Vector{Float64}     # survivors' mean genotype
    Gs::Matrix{Float64}     # survivors' G
    zoff::Vector{Float64}
    Goff::Matrix{Float64}
    extinct::Bool
    N::Float64
    Noff::Float64
    erosion::Bool
    VE::Float64             # environmental variance of each trait (a mixture of the demes' with migrants)
    weighted::Bool
    W2bar::Float64
    M::Matrix{Float64}
    shape::Matrix{Float64}
    u::Vector{Float64}
    cov::Any                # `nothing` (loci = :aggregate) or the haplotype moments, a `TraitCov`
    sinv::Matrix{Float64}   # (Ω + P)⁻¹, zero when neutral
    θ0::Vector{Float64}
    noise::Matrix{Float64}
    f::Vector{Float64}
    idx::Int
end

"""
The haplotype moments of the multilocus model by pair class, as matrices (src/multilocus.jl),
for k demes: `T[a, c]`, `R[a, c]` and `Q[a, b, c]` for deme a (and b) and pair class c, with
`Q[a, a, ⋅]` the moments of two distinct adults of one deme and `Q[a, b, ⋅] = Q[b, a, ⋅]ᵀ`
those of two demes. As moments.jl's `DemeCov`, with matrices.
"""
mutable struct TraitCov
    T::Matrix{Matrix{Float64}}
    R::Matrix{Matrix{Float64}}
    Q::Array{Matrix{Float64},3}
    ploidy::Int
    drift::Bool
    ω::Vector{Float64}      # QTLSim.pair_classes(p)
    ρ::Vector{Float64}
    origins::Any            # :exchange, :pool or a k² × k² matrix
    cache::Dict{Any,Any}    # (M, N, alive) => (P, J)
    α2::Vector{Any}         # per deme: loci = :rare_alleles, the per-mutation effect covariance C; 0.0: Gaussian
    clonal::Bool            # QTLSim.hap_clonal: the genome is inherited whole (haploid clones)
    gain::Any               # loci = :rare_local: per deme, the gain matrix Λ (src/local_alleles.jl); nothing otherwise
end

function TraitCov(p::Params, N::AbstractVector, G::AbstractVector; drift::Bool, origins = :exchange, α2 = 0.0,
                  localal::Bool = false)
    k, n = length(N), p.NoTraits
    ω, ρ = QTLSim.pair_classes(p)
    nc = length(ω)
    T = [zeros(n, n) for _ in 1:k, _ in 1:nc]
    R = [zeros(n, n) for _ in 1:k, _ in 1:nc]
    Q = [zeros(n, n) for _ in 1:k, _ in 1:k, _ in 1:nc]
    for a in 1:k
        st = QTLSim.hap_init(Matrix{Float64}(G[a]), drift ? float(N[a]) : Inf, p.Ploidy, nc)
        T[a, :], R[a, :], Q[a, a, :] = st
    end
    α2s = α2 isa AbstractVector ? Any[a for a in α2] : Any[α2 for _ in 1:k]
    length(α2s) == k || throw(ArgumentError("α2 gives $(length(α2s)) values for $k demes"))
    TraitCov(T, R, Q, p.Ploidy, drift, ω, ρ, origins, Dict{Any,Any}(), α2s, QTLSim.hap_clonal(p.Ploidy, ρ),
             localal ? [zeros(n, n) for _ in 1:k] : nothing)
end

_nd(c::TraitCov, N) = c.drift ? max(N, 2.0) : Inf
_state(c::TraitCov, a) = (c.T[a, :], c.R[a, :], c.Q[a, a, :])
_setstate!(c::TraitCov, a, st) = (c.T[a, :] = st[1]; c.R[a, :] = st[2]; c.Q[a, a, :] = st[3]; c)
_vpop(c::TraitCov, a, N) = sym(QTLSim.hap_vpop(_state(c, a), _nd(c, N), c.ploidy))
_cmean(c::TraitCov, a, N) = QTLSim.hap_cmean(sum(c.T[a, :]), sum(c.R[a, :]), sum(c.Q[a, a, :]), _nd(c, N), c.ploidy)

function TraitMoments(p::Params; zbar = zeros(p.NoTraits), G = zeros(p.NoTraits, p.NoTraits), erosion = false,
                      loci::Symbol = :multilocus, cov = :new, idx = 1)
    n, VE = p.NoTraits, QTLSim.envVar(p)
    Id = Matrix{Float64}(I, n, n)
    loci === :auto && (loci = QTLSim.moments_loci(p))
    cov === :new && (cov = loci in (:multilocus, :rare_alleles, :rare_local) ?
                     TraitCov(p, [float(p.PopSize)], [G]; drift = erosion, localal = loci === :rare_local,
                              α2 = loci in (:rare_alleles, :rare_local) ? QTLSim.mutation_effect_covariance(p) : 0.0) : nothing)
    TraitMoments(p, 0, collect(Float64, zbar), Matrix{Float64}(G), G + VE * Id, zeros(n), QTLSim.initFitSD(p),
                 p.PopSize, 1.0, zeros(n), collect(Float64, zbar), Matrix{Float64}(G), collect(Float64, zbar),
                 Matrix{Float64}(G), false, p.PopSize, p.PopSize, erosion, VE, false, 1.0,
                 QTLSim.mutational_covariance(p), QTLSim.selection_shape(p), QTLSim.move_direction(p),
                 cov, zeros(n, n), zeros(n), QTLSim.move_noise_covariance(p), QTLSim.fluct_direction(p), idx)
end

isextinct(m::TraitMoments) = m.extinct

sym(A) = (A + A') / 2

# the optimum: a given one (`θ`, from a realised sequence), the fluctuating moves at their
# expectation (the periodic move exactly, red noise as its trend; QTLSim.moveoptimum_fluct!),
# or QTSim's moves along u (moveoptimum_traits!)
function trait_moments_optimum!(m::TraitMoments, θ)
    p = m.p
    if θ !== nothing
        m.θ = collect(Float64, θ)
    elseif p.MoveFitnessSDGen != -1 && m.gen >= p.MoveFitnessSDGen
        d = p.MoveFitnessSD
        if p.MoveType === :red_noise || p.MoveType === :periodic
            t = m.gen - p.MoveFitnessSDGen
            t == 0 && (m.θ0 = copy(m.θ))
            c = p.MoveType === :periodic ? p.MoveAmplitude * sin(2π * t / p.MovePeriod) : 0.0
            m.θ = m.θ0 .+ (d * t) .* m.u .+ c .* m.f
        else
            m.θ = p.MoveType === :continuous ? m.θ .+ d .* m.u :
                  m.zbar .+ (p.MoveType === :phenotype ? sqrt(p.PhenotypeVar) : m.ω) * d .* m.u
        end
    end
    m
end

function trait_moments_select!(m::TraitMoments; θ = nothing)
    p = m.p
    trait_moments_optimum!(m, θ)
    n = length(m.zbar)
    if p.Drift
        m.Wbar, m.W2bar, m.S, m.zs, m.Gs = 1.0, 1.0, zeros(n), copy(m.zbar), copy(m.G)
        m.sinv = zeros(n, n)
    else
        # z ~ N(z̄, P), W ∝ exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ)): the parents are N(z̄ + P s⁻¹ λ, P − P s⁻¹ P)
        # with s = Ω + P and λ = θ − z̄, and g regresses on z with coefficient G P⁻¹
        Ω = m.ω^2 .* m.shape
        s = Ω + m.P
        λ = m.θ - m.zbar
        sλ = s \ λ
        m.Wbar = FITNESS_PEAK * sqrt(det(Ω) / det(s)) * exp(-dot(λ, sλ) / 2)
        m.S = m.P * sλ
        m.zs = m.zbar + m.G * sλ
        m.Gs = sym(m.G - m.G * (s \ m.G))
        s2 = Ω + 2m.P                                        # E[W²]: Ω/2 in place of Ω
        m.W2bar = FITNESS_PEAK^2 * sqrt(det(Ω) / det(s2)) * exp(-dot(λ, s2 \ λ))
        m.sinv = sym(inv(s))
    end
    m.weighted = p.Selection === :fecundity
    m.nsurv = m.weighted ? m.N : m.N * m.Wbar
    m.extinct = m.weighted ? m.N < p.Ploidy : m.nsurv < (p.Ploidy == 2 ? 2 : 1)
    c = m.cov
    if c !== nothing && !p.Drift
        a = m.idx
        st, B = QTLSim.hap_select(_state(c, a), m.sinv, _nd(c, m.N), c.ω, c.ploidy; α2 = c.α2[a], clonal = c.clonal,
                                  gain = c.gain === nothing ? nothing : c.gain[a])
        # the deviation of this deme's mean shrinks by I − B, and with it its cross-moments
        for b in axes(c.Q, 1)
            b == a && continue
            BQ = B * sum(c.Q[a, b, :])
            for k in eachindex(c.ω)
                c.Q[a, b, k] = c.Q[a, b, k] - c.ω[k] .* BQ
                c.Q[b, a, k] = Matrix(c.Q[a, b, k]')
            end
        end
        _setstate!(c, a, st)
    end
    m
end

# as moments_reproduce!: offspring inherit the parents' mean; their G comes from the
# haplotype moments, or with loci = :aggregate it is the parents', with optional drift
function trait_moments_reproduce!(m::TraitMoments)
    p = m.p
    m.zoff, m.Goff = copy(m.zs), copy(m.Gs)
    npar = m.weighted ? m.N * m.Wbar^2 / m.W2bar : m.nsurv
    c = m.cov
    if c === nothing
        m.erosion && (m.Goff .*= 1 - 1 / (2 * max(npar, 1.0)))
    else
        _setstate!(c, m.idx, QTLSim.hap_reproduce(_state(c, m.idx), c.drift ? max(npar, 1.0) : Inf, c.ρ, c.ploidy))
    end
    R, K = p.GrowthRate, p.PopSize
    nw = m.weighted ? m.N * m.Wbar : m.nsurv
    m.Noff = p.Demography === :fixed ? K : min(K, R * nw / (1 + (R - 1) * nw / K))
    m.Noff < 1 && (m.extinct = true)
    m
end

function trait_moments_record!(rec, m::TraitMoments)
    (rec === nothing || m.gen % m.p.Sample != 0) && return m
    push!(rec.gen, m.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[m])
    Gr = m.cov !== nothing && m.cov.drift ? m.G .* (m.N / (m.N - 1)) : m.G      # the expected sample covariance with drift
    push!(rec.pheno_mean, m.zbar[1]); push!(rec.pheno_var, m.P[1, 1]); push!(rec.geno_var, Gr[1, 1])
    push!(rec.nsurv, round(Int, m.nsurv)); push!(rec.theta, m.θ[1]); push!(rec.omega, m.ω)
    push!(rec.mean_fitness, m.Wbar); push!(rec.sel_diff, m.S[1]); push!(rec.popsize, round(Int, m.N))
    push!(rec.trait_mean, copy(m.zbar)); push!(rec.gmatrix, copy(Gr)); push!(rec.pmatrix, copy(m.P))
    push!(rec.trait_lag, m.θ - m.zbar); push!(rec.trait_sel_diff, copy(m.S))
    m
end

"""
    trait_moment_images(rec = nothing)

The image of each step of the life cycle in the moments lowering of several traits, by
name: the multivariate Gaussian recursion, `Δz̄ = G (Ω + P)⁻¹ (θ − z̄)`, with `G` from the
haplotype moments of the multilocus model (src/multilocus.jl), or with `loci = :aggregate`
`G' = G − G (Ω + P)⁻¹ G + M`. A second table for the same generators as `moment_images`.
"""
trait_moment_images(rec = nothing) = Dict{Symbol,Any}(
    :develop => Operator(:develop, (m, rng) -> (m.P = m.G + m.VE * Matrix{Float64}(I, size(m.G)...); m)),
    :select => Operator(:select, (m, rng) -> trait_moments_select!(m)),
    :reproduce => Operator(:reproduce, (m, rng) -> m.extinct ? m : trait_moments_reproduce!(m)),
    :record => Operator(:record, (m, rng) -> trait_moments_record!(rec, m)),
    :rescale => Operator(:rescale, (m, rng) -> (m.gen == m.p.ChangeFitnessGen && (m.ω = m.p.ChangeFitnessSD); m)),
    :mutate => Operator(:mutate, (m, rng) -> trait_moments_mutate!(m)),
    :swap => Operator(:swap, (m, rng) -> trait_moments_swap!(m)))

function trait_moments_mutate!(m::TraitMoments)
    m.Goff = m.Goff + m.M
    c = m.cov
    c === nothing || _setstate!(c, m.idx, QTLSim.hap_mutate(_state(c, m.idx), m.M, c.ploidy))
    m
end

function trait_moments_swap!(m::TraitMoments)
    m.zbar = m.zoff; m.G = m.Goff; m.N = m.Noff; m.VE = QTLSim.envVar(m.p); m.gen += 1
    c = m.cov
    c === nothing || m.extinct || (m.G = _vpop(c, m.idx, m.N))
    m
end

moment_images(::TraitModel, rec = nothing) = trait_moment_images(rec)

# the state of n traits that a step's parameters replace (compose.jl, with_step_params), in
# the moments: the optimum vector, the fitness shape, the mutational covariance, the move's
# direction, the fluctuation's noise covariance and the cycle's direction
trait_slots(m::TraitMoments) = (m.θ, m.shape, m.M, m.u, m.noise, m.f)
set_trait_slots!(m::TraitMoments, (θ, shape, M, u, noise, f)) =
    ((m.θ, m.shape, m.M, m.u, m.noise, m.f) = (θ, shape, M, u, noise, f); m)
derived_slots(::TraitMoments, q::Params, Ω) =
    (Ω === nothing ? QTLSim.selection_shape(q) : Matrix{Float64}(Ω), QTLSim.mutational_covariance(q),
     QTLSim.move_direction(q), QTLSim.move_noise_covariance(q), QTLSim.fluct_direction(q))

# --- several demes ------------------------------------------------------------

"The moments of the demes of a metapopulation of several traits, which a `DemeOp` updates."
struct TraitDemeMoments
    demes::Vector{TraitMoments}
end

isextinct(s::TraitDemeMoments) = all(isextinct, s.demes)

# The image of `migrate` on the moments of several traits: as moments_migrate! (moments.jl)
# with vectors for the means and outer products of their differences for their squares.
function moments_migrate!(s::TraitDemeMoments, demes, rate)
    ms = s.demes[demes]
    k = length(ms)
    M = rate isa AbstractMatrix ? QTLSim.check_migration(rate, k) :
        [a == b ? 0.0 : Float64(rate) / (k - 1) for a in 1:k, b in 1:k]
    alive = [!m.extinct for m in ms]
    c = ms[1].cov
    Nv = [m.N for m in ms]
    P, J = trait_mixing(c, M, Nv, alive)
    z, VE = [m.zbar for m in ms], [m.VE for m in ms]
    z2 = [alive[a] ? sum(P[a, cc] .* z[cc] for cc in 1:k) : z[a] for a in 1:k]
    if c === nothing
        G = [m.G for m in ms]
        for a in 1:k
            alive[a] || continue
            ms[a].G = sym(sum(P[a, cc] .* (G[cc] + (z[cc] - z2[a]) * (z[cc] - z2[a])') for cc in 1:k))
        end
    else
        ix = [m.idx for m in ms]
        h, n = c.ploidy, length(z[1])
        Δ = [(z[cc] - z2[a]) ./ h for a in 1:k, cc in 1:k]
        # the outer products of the differences of the means, which every pair class carries
        # in its share ω: Σ_c P_ac Δ_ac Δ_acᵀ, and Σ_cd J[ab, cd] Δ_ac Δ_bdᵀ = D_a J_ab D_bᵀ with
        # D_a = [Δ_a1 … Δ_ak] and J_ab the k × k block of J for the pair (a, b)
        D = [reduce(hcat, (Δ[a, cc] for cc in 1:k)) for a in 1:k]
        ET = [sum(P[a, cc] .* (Δ[a, cc] * Δ[a, cc]') for cc in 1:k) for a in 1:k]
        EQ = [D[a] * reshape(J[a + (b - 1) * k, :], k, k) * D[b]' for a in 1:k, b in 1:k]
        for j in eachindex(c.ω)
            κ = c.ω[j]
            T, R = c.T[ix, j], c.R[ix, j]
            Qm = reshape(reduce(hcat, (vec(c.Q[ix[cc], ix[d], j]) for cc in 1:k, d in 1:k)), n * n, k * k)
            QJ = Qm * J'                                 # column a + (b−1)k: Σ_cd J[ab, cd] vec(Q_cd)
            for a in 1:k
                alive[a] || continue
                c.T[ix[a], j] = sum(P[a, cc] .* T[cc] for cc in 1:k) + κ .* ET[a]
                c.R[ix[a], j] = sum(P[a, cc] .* R[cc] for cc in 1:k) + κ .* ET[a]
            end
            Q2 = [alive[a] && alive[b] ? reshape(QJ[:, a + (b - 1) * k], n, n) + κ .* EQ[a, b] : c.Q[ix[a], ix[b], j]
                  for a in 1:k, b in 1:k]
            for a in 1:k, b in 1:k
                c.Q[ix[a], ix[b], j] = (Q2[a, b] + Q2[b, a]') ./ 2
            end
        end
    end
    for a in 1:k
        alive[a] || continue
        ms[a].zbar = z2[a]
        ms[a].VE = sum(P[a, :] .* VE)
        c === nothing || (ms[a].G = _vpop(c, ms[a].idx, max(ms[a].N, 2.0)))
    end
    c === nothing || c.gain === nothing || set_gain!(c, [m.idx for m in ms], P, alive, [m.sinv * (m.θ - m.zbar) for m in ms])
    s
end

# the backward matrix and the law of the origins of pairs, as cov_mixing (moments.jl)
function trait_mixing(c, M, Nv, alive)
    (c === nothing || !c.drift) && (P = migration_mixing(M, Nv, alive); return (P, kron(P, P)))
    get!(c.cache, (M, Nv, alive)) do
        P = all(alive) ? QTLSim.backward_migration(M, Nv) : migration_mixing(M, Nv, alive)
        J = c.origins === :exchange ? QTLSim.exchange_origins(P, Nv) : c.origins === :pool ? kron(P, P) : c.origins
        size(J) == (length(Nv)^2, length(Nv)^2) ||
            throw(ArgumentError("origins must be $(length(Nv)^2) × $(length(Nv)^2); got $(size(J))"))
        (P, J)
    end
end

# the census statistics of `between = true` for several traits, over the living demes: as
# cov_summary (moments.jl), with matrices; `qst` is per trait
function trait_cov_summary(ms::AbstractVector{TraitMoments})
    ms = [m for m in ms if !m.extinct]
    k = length(ms)
    k == 0 && return nothing
    c = ms[1].cov
    Nv = [max(m.N, 2.0) for m in ms]
    within = [_vpop(c, m.idx, N) .* (N / (N - 1)) for (m, N) in zip(ms, Nv)]
    C = [_cmean(c, m.idx, N) for (m, N) in zip(ms, Nv)]
    h = c.ploidy
    pairs = [a == b ? zero(within[1]) :
             sym((ms[a].zbar - ms[b].zbar) * (ms[a].zbar - ms[b].zbar)' ./ 2 +
                 h^2 / 2 .* (C[a] + C[b] - sum(c.Q[ms[a].idx, ms[b].idx, :]) - sum(c.Q[ms[a].idx, ms[b].idx, :])'))
             for a in 1:k, b in 1:k]
    VW = sum(within) ./ k
    n = size(VW, 1)
    k == 1 && return (; within, pairs, VW, among = zero(VW), between = zero(VW), qst = zeros(n))
    among = sum(pairs) ./ (k * (k - 1))
    between = among - sum(within[a] ./ Nv[a] for a in 1:k) ./ k
    σ2 = (k - 1) / k .* between
    (; within, pairs, VW, among, between, qst = [σ2[i, i] / (σ2[i, i] + 2VW[i, i]) for i in 1:n])
end

trait_metapop_records() = (gen = Int[], within = Matrix{Float64}[], var_among = Matrix{Float64}[],
                           between = Matrix{Float64}[], qst = Vector{Float64}[], pairs = Matrix{Matrix{Float64}}[])

function record_trait_metapop!(mp, gen, ms)
    x = trait_cov_summary(ms)
    x === nothing && return mp
    push!(mp.gen, gen); push!(mp.within, x.VW); push!(mp.var_among, x.among); push!(mp.between, x.between)
    push!(mp.qst, x.qst); push!(mp.pairs, x.pairs)
    mp
end

# `loci` resolved for the demes `ps` of several traits and the life cycle `cycle`, with the
# per-mutation effect covariance C of each deme (0.0 unless :rare_alleles), as moments.jl's
# `resolve_loci`: the mutate steps of the life cycle must carry the same mutation fields
# (MutationScalar, MutationRatePerGamete, MutCorr, MutCov), which then set C; otherwise
# `:auto` is :multilocus and an explicit :rare_alleles an error.
function resolve_trait_loci(ps::AbstractVector{Params}, cycle, loci::Symbol)
    loci in (:multilocus, :aggregate) && return loci, Any[0.0 for _ in ps]
    fields = (:MutationScalar, :MutationRatePerGamete, :MutCorr, :MutCov)
    fs = unique([NamedTuple(k => v for (k, v) in pairs(stepparams(g)) if k in fields)
                 for g in stepsof(cycle) if !(g isa ObExpr) && stepname(g) === :mutate])
    if length(fs) > 1
        loci === :auto && return :multilocus, Any[0.0 for _ in ps]
        throw(ArgumentError("loci = :rare_alleles needs one law of mutation effects per deme; the mutate steps of the " *
                            "life cycle carry different mutation fields"))
    end
    qs = isempty(fs) ? ps : [Params(q; only(fs)...) for q in ps]
    if loci === :auto
        all(q -> QTLSim.moments_loci(q) === :rare_alleles, qs) || return :multilocus, Any[0.0 for _ in ps]
        loci = :rare_alleles
    end
    loci, Any[QTLSim.mutation_effect_covariance(q) for q in qs]
end

"""
    run_trait_moments(p; cycle = lifecycle().generation, images, zbar, G, erosion = false, optimum = nothing,
                      loci = :auto)
    run_trait_moments(ps::Vector{Params}; cycle, images, zbar, G, erosion = false, migration_rate = 0.0,
                      between = false, origins = :exchange, loci = :auto)

The moments lowering of a model of several traits: `cycle` lowered with the table of
`TraitModel(n)` (`trait_moment_images`) and iterated `NoGenerations` times from mean `zbar`
and genetic covariance `G` (zero by default; with several demes, one or one per deme).
`run_moments(p)` calls it when `p.NoTraits > 1`; called directly with one trait it runs the
multivariate recursion with 1×1 matrices, which is the one-trait recursion with `loci =
:multilocus`. `loci = :auto`, the default, resolves as `QTLSim.moments_loci` (the multilocus
Gaussian model while that is its resolution for several traits; `loci = :rare_alleles` is the
rare-alleles closure of src/rare_alleles_traits.jl, converging to `rare_alleles_G(p)`). The
multilocus model (converging to `gaussian_G_loci(p)`, or with `erosion = true` to
`gaussian_G_loci(p; drift = true)`); `loci = :aggregate` is the old recursion
(`gaussian_G(p)`). `optimum`, a vector of optimum vectors, one per generation, replaces the
move (for one deme; for example a realised red-noise sequence). A fluctuating move lowers to
its expected optimum: the periodic move exactly, red noise as its trend.

With `ps`, one `Params` per deme, a life cycle on `k = length(ps)` demes, coupled by
`migrate` steps at `migration_rate` (a number or a matrix) or their own rates: the k-deme
moments of moments.jl with matrices. With `between = true` (or `erosion = true`) the demes
are finite and drift; `between = true` also returns `metapop`, per recorded generation the
mean within-deme G (`within`, expected sample covariance), the expected sample covariance of
the deme mean genotypes (`var_among`), its part beyond sampling (`between`), `Q_ST` per trait
(`qst`) and the matrix of half the expected outer products of the differences of deme means
(`pairs`).
"""
function run_trait_moments(p::Params; cycle = lifecycle().generation, kw...)
    k = nfactors(dom(cycle))
    run_trait_moments(fill(p, max(k, 1)); cycle, kw...)
end

function run_trait_moments(ps::AbstractVector{Params}; cycle = lifecycle().generation, images = Dict{Symbol,Any}(),
                           zbar = zeros(ps[1].NoTraits), G = zeros(ps[1].NoTraits, ps[1].NoTraits), erosion = false,
                           optimum = nothing, loci::Symbol = :auto, migration_rate::Union{Real,AbstractMatrix} = 0.0,
                           between = false, origins::Union{Symbol,AbstractMatrix} = :exchange)
    k = length(ps)
    foreach(QTLSim.checkparams, ps)
    check_loci(loci, between)
    loci, α2 = resolve_trait_loci(collect(Params, ps), cycle, loci)
    between && check_between(ps, origins)
    p = ps[1]
    n = p.NoTraits
    all(q -> q.NoTraits == n, ps) || throw(ArgumentError("the demes of a metapopulation have the same number of traits"))
    all(q -> q.NoGenerations == p.NoGenerations, ps) ||
        throw(ArgumentError("the demes of a metapopulation run for the same number of generations"))
    loci !== :aggregate && any(q -> (q.Ploidy, q.NoQTL, q.RecombPoissonMean, q.CrossoverModel) !=
                                     (p.Ploidy, p.NoQTL, p.RecombPoissonMean, p.CrossoverModel), ps) &&
        throw(ArgumentError("the demes of a metapopulation share Ploidy, NoQTL, RecombPoissonMean and CrossoverModel"))
    zk = zygotes_word(k)
    dom(cycle) == zk && codom(cycle) == zk ||
        throw(ArgumentError("the moments of $k deme$(k == 1 ? "" : "s") of several traits iterate a life cycle on $k zygote " *
                            "stage$(k == 1 ? "" : "s"); got $(dom(cycle)) → $(codom(cycle))"))
    if optimum !== nothing          # a given optimum per generation, optimum[t + 1], replaces the move
        k == 1 || throw(ArgumentError("a given optimum sequence is for one deme"))
        length(optimum) >= p.NoGenerations ||
            throw(ArgumentError("optimum gives $(length(optimum)) generations; NoGenerations is $(p.NoGenerations)"))
        all(θ -> length(θ) == n, optimum) || throw(ArgumentError("each optimum is a vector of $n traits"))
        images = merge(Dict{Any,Any}(:select => Operator(:select, (m, rng) -> trait_moments_select!(m; θ = optimum[m.gen + 1]))),
                       Dict{Any,Any}(images))
    end
    rec = merge(records(), QTLSim.trait_records())
    k > 1 && (rec = merge((gen = rec.gen, deme = Int[]), Base.structdiff(rec, NamedTuple{(:gen,)}), (demeof = IdDict{TraitMoments,Int}(),)))
    # the model of n traits lowers with its own table; called with one trait, the same table runs 1×1 matrices
    model = model_of(p)
    step = demeop(lower(cycle, :moments; images = model isa TraitModel ? images : merge(trait_moment_images(rec), Dict{Any,Any}(images)),
                        record = rec, model, migration_rate))
    perdeme(x, name) = x isa AbstractVector && !(x isa AbstractVector{<:Real}) ?
        (length(x) == k ? x : throw(ArgumentError("$name gives $(length(x)) values for $k demes"))) : fill(x, k)
    zs, Gs = perdeme(zbar, "zbar"), perdeme(G, "G")
    drift = between || erosion
    cov = loci === :aggregate ? nothing : TraitCov(p, [float(q.PopSize) for q in ps], Gs; drift, origins, α2,
                                                            localal = loci === :rare_local)
    s = TraitDemeMoments([TraitMoments(q; zbar = z, G = g, erosion, loci, cov, idx = i)
                          for (i, (q, z, g)) in enumerate(zip(ps, zs, Gs))])
    k > 1 && for (i, m) in enumerate(s.demes); rec.demeof[m] = i; end
    mp = trait_metapop_records()
    d = collect(1:k)
    extinct_gen = fill(-1, k)
    gen = 0                      # whole cycles, as in run_composed, with the demes in step
    while gen < p.NoGenerations
        foreach(m -> m.extinct || (m.gen = gen), s.demes)
        between && gen % p.Sample == 0 && record_trait_metapop!(mp, gen, s.demes)
        step.act(s, d, nothing)
        for (i, m) in enumerate(s.demes)
            m.extinct && extinct_gen[i] < 0 && (extinct_gen[i] = m.gen)
        end
        isextinct(s) && break
        d = d[step.perm]
        gen = max(gen + 1, maximum(m.gen for m in s.demes if !m.extinct))
    end
    out = k == 1 ? (extinct = isextinct(s), extinct_gen = extinct_gen[1], rec...) :
                   (extinct = isextinct(s), extinct_gen, Base.structdiff(rec, NamedTuple{(:demeof,)})...)
    between ? (out..., metapop = mp) : out
end
