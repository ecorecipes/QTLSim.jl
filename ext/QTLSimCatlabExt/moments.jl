# The third lowering of the life cycle: to deterministic recursions on the moments of a
# Gaussian population, the infinite-population approximation of quantitative genetics
# (Lande 1975, 1976; Bürger and Lynch 1995). The same morphism of the life-cycle category
# that the kernel and SLiM simulate is lowered to its predicted means and variances, so
# a simulation can be checked against theory step for step, and any life cycle built
# from these steps has a prediction. The population is the multilocus Gaussian model
# (src/multilocus.jl, `loci = :multilocus`), whose fixed point is `gaussian_vg_loci`, by
# default with the rare-alleles closure (src/rare_alleles.jl, `loci = :rare_alleles`, fixed
# point `rare_alleles_vg`) for one diploid trait (`loci = :auto`, resolved by
# `QTLSim.moments_loci` and `resolve_loci` below); `loci = :aggregate` is the one-locus
# recursion of `gaussian_vg_exact`. The algebra is in
# lean/CategoricalQuantGen/Multilocus.lean and Theory.lean.

"""
The state of the moments lowering: the parents' mean genotype (the mean phenotype, as
environmental deviations have mean zero), genetic and phenotypic variances, the optimum
and width, the expected survivors and mean fitness, the survivors' moments and the
offspring's.
"""
mutable struct Moments
    p::Params
    gen::Int
    zbar::Float64
    VG::Float64
    VP::Float64
    θ::Float64
    ω::Float64
    nsurv::Float64
    Wbar::Float64
    S::Float64           # selection differential on the phenotype
    zs::Float64          # survivors' mean genotype
    VGs::Float64         # survivors' genetic variance
    zoff::Float64
    VGoff::Float64
    extinct::Bool
    N::Float64           # expected population size
    Noff::Float64
    θ0::Float64          # the optimum when a fluctuating move starts
    erosion::Bool        # drift erodes V_G by V_G/(2N_e) per generation, N_e = expected survivors
    VE::Float64          # environmental variance: V_E, or with migrants a mixture of the demes' V_E
    weighted::Bool       # Selection = :fecundity: all survive and mate in proportion to W
    W2bar::Float64       # the mean of W², for the number of parents under fecundity selection
    cov::Any             # `nothing`, or the demes' expected second moments (`DemeCov`, `between = true`)
    idx::Int             # this deme's index in `cov`
    sinv::Float64        # 1/s, s = ω² + V_P the width of selection on the phenotype (0 when neutral)
end

Moments(p::Params; zbar = 0.0, VG = 0.0, erosion = false, cov = nothing, idx = 1) =
    Moments(p, 0, zbar, VG, VG + QTLSim.envVar(p), 0.0, QTLSim.initFitSD(p), p.PopSize, 1.0, 0.0,
            zbar, VG, zbar, VG, false, p.PopSize, p.PopSize, 0.0, erosion, QTLSim.envVar(p), false, 1.0,
            cov, idx, 0.0)

# --- the multilocus model: expected second moments of the haplotypes ----------------
#
# Unless `loci = :aggregate`, the lowering carries, for each of its k demes, the
# expected second moments of the genes about the demes' deterministic means z̄_a (the
# `zbar` of each `Moments`), in units of a haplotype: one gamete's alleles summed over the
# loci, so that a diploid adult's genotype is h₁ + h₂ and its expected mean half the
# deme's z̄ (a haploid's genotype is its one haplotype).
#
#   T[a]    = E[(h − z̄_a/2)²]                  one haplotype of an adult of deme a
#   R[a]    = E[(h₁ − z̄_a/2)(h₂ − z̄_a/2)]      the two haplotypes of one adult (diploid)
#   Q[a, b] = E[(h − z̄_a/2)(h′ − z̄_b/2)]       haplotypes of two distinct adults of a and b
#
# each resolved by class of pairs of loci (`QTLSim.pair_classes`): class 1, a locus with
# itself, is the genic part (Σ_l over loci of the per-locus moments); the other classes are
# the covariances of alleles at different loci (linkage disequilibrium, LD), grouped by
# the probability ρ that a gamete splits the pair, so that each class recombines at its
# own rate (adjacent loci rarely, distant ones often). The full moment is the sum. The genetic
# variance within deme a is `V_G = 2T + 2R − 4C` (haploid `T − C`), with C the expected
# square of the deviation of the deme's mean haplotype from z̄_a/2.
#
# The images of select, reproduce and mutate on one deme's moments are QTLSim's
# `hap_select`, `hap_reproduce` and `hap_mutate`, whose equations are in
# src/multilocus.jl: Gaussian selection of width s lowers a haplotype's within-deme
# variance by D = H²/s (H = V_G/2, diploid), the genic parts by D/L (each locus carries
# 1/L of H; Lande's multilocus Gaussian), the rest is negative LD; recombination at the
# rates of `QTLSim.pair_classes(p)` turns LD back into variance, class by class; mutation adds V_m/Ploidy to
# each haplotype. With `loci = :rare_alleles` the genic part also loses β²α²/(2s) of
# itself (`hap_select`'s `α2`, one α² per deme). Without drift (the default) the demes are
# infinite, C = Q = 0, and one deme converges to `gaussian_vg_loci(p)` (`rare_alleles_vg(p)`
# with the rare-alleles closure). With drift (`erosion = true`, or `between =
# true`, which also records the statistics among demes) each deme has its expected N
# adults and n parents, and selection also shrinks the deviation of the deme's mean,
# C ← (1 − β)² C and Q_ab ← (1 − β) Q_ab with β = V_G/s. Each step's image (`cov_*!`) acts
# on these moments; for a neutral additive trait under the continuum-of-alleles mutation
# of `mutate!`, full = genic, and the expected squared difference of two genes is V_m
# times their expected coalescence time truncated at the start of the run: the steps
# composed are `deme_coalescence`'s recursion read forward, exactly.
#
# Migration (backward matrix P, law J of the origins of two distinct adults): with
# Δ_ac = (z̄_c − z̄′_a)/Ploidy, T′_a = Σ_c P_ac (T_c + Δ_ac²), R′_a = Σ_c P_ac (R_c + Δ_ac²),
# Q′_ab = Σ_cd J[ab, cd] (Q_cd + Δ_ac Δ_bd); in the genic parts Δ_ac Δ_bd is weighted by
# 1/L, the deterministic difference of the means spread equally over the loci. Without
# drift P is `migration_mixing` and J = P ⊗ P, so Q stays 0.
#
# With `loci = :aggregate` there are no haplotype moments (`cov = nothing`): selection's
# reduction V_G²/s is inherited whole, V ← V − V²/(V_s + V) + V_m, exact only for one
# haploid locus (or haploid clones), with fixed point `gaussian_vg_exact`; `erosion`
# multiplies V_G by 1 − 1/(2n). It is the recursion of QTLSim before the multilocus model.

"The expected second moments of the haplotypes of k demes (`loci = :multilocus` or `:rare_alleles`), by pair class; see the comment above."
mutable struct DemeCov
    T::Matrix{Float64}                   # k × classes
    R::Matrix{Float64}
    Q::Array{Float64,3}                  # k × k × classes
    ω::Vector{Float64}                   # the pair classes (QTLSim.pair_classes): weights,
    ρ::Vector{Float64}                   # and split probabilities; class 1 is genic
    origins::Any                         # :exchange, :pool or a k² × k² matrix
    cache::Dict{Any,Any}                 # (M, N, alive) => (P, J)
    ploidy::Int
    drift::Bool                          # finite demes (N adults, n parents); false: infinite
    α2::Vector{Float64}                  # per deme; loci = :rare_alleles: the mutation effect variance α²; 0: Gaussian
end

# demes start as Hardy–Weinberg populations in linkage equilibrium, with exactly their
# mean z̄ and population variance V
function DemeCov(N::AbstractVector, V::AbstractVector; origins = :exchange, ploidy::Integer = 2, drift::Bool = true,
                 classes = ([1.0, 0.0], [0.0, 0.0]), α2 = 0.0)
    k = length(N)
    ω, ρ = classes
    nc = length(ω)
    st = [QTLSim.hap_init(float(V[a]), drift ? float(N[a]) : Inf, ploidy, nc) for a in 1:k]
    T, R, Q = zeros(k, nc), zeros(k, nc), zeros(k, k, nc)
    for a in 1:k; T[a, :] = st[a][1]; Q[a, a, :] = st[a][3]; end
    α2s = α2 isa Real ? fill(float(α2), k) : collect(Float64, α2)
    length(α2s) == k || throw(ArgumentError("α2 gives $(length(α2s)) values for $k demes"))
    DemeCov(T, R, Q, collect(float.(ω)), collect(float.(ρ)), origins, Dict{Any,Any}(), ploidy, drift, α2s)
end
# α2: one value per deme, from `resolve_loci` (0 for the multilocus Gaussian model)
DemeCov(p::Params, N::AbstractVector, V::AbstractVector; α2 = 0.0, kw...) =
    DemeCov(N, V; ploidy = p.Ploidy, classes = QTLSim.pair_classes(p), α2, kw...)

# the number of adults (parents) of a deme in the moments: infinite without drift
_nd(c::DemeCov, N) = c.drift ? max(N, 2.0) : Inf
_np(c::DemeCov, n) = c.drift ? max(n, 1.0) : Inf
_state(c::DemeCov, a) = (c.T[a, :], c.R[a, :], c.Q[a, a, :])
_setstate!(c::DemeCov, a, st) = (c.T[a, :] = st[1]; c.R[a, :] = st[2]; c.Q[a, a, :] = st[3]; c)
# the expected square of the deviation of deme a's mean haplotype, C_aa
_cmean(c::DemeCov, a, N) = QTLSim.hap_cmean(sum(c.T[a, :]), sum(c.R[a, :]), sum(c.Q[a, a, :]), _nd(c, N), c.ploidy)
# the expected within-deme population variance (divisor N) of genotypes
_vpop(c::DemeCov, a, N) = QTLSim.hap_vpop(_state(c, a), _nd(c, N), c.ploidy)

function cov_select!(m::Moments)
    c, a = m.cov, m.idx
    m.sinv == 0 && return m                 # neutral: the identity
    st, β = QTLSim.hap_select(_state(c, a), m.sinv, _nd(c, m.N), c.ω, c.ploidy; α2 = c.α2[a],
                              clonal = QTLSim.hap_clonal(c.ploidy, c.ρ))
    for b in axes(c.Q, 1)
        b == a && continue
        Qab = sum(c.Q[a, b, :])
        for k in eachindex(c.ω)
            c.Q[a, b, k] -= c.ω[k] * β * Qab; c.Q[b, a, k] = c.Q[a, b, k]
        end
    end
    _setstate!(c, a, st)
    m
end

cov_reproduce!(m::Moments, npar) =
    (_setstate!(m.cov, m.idx, QTLSim.hap_reproduce(_state(m.cov, m.idx), _np(m.cov, npar), m.cov.ρ, m.cov.ploidy)); m)

function moments_mutate!(m::Moments)
    Vm = QTLSim.muteVar(m.p)
    m.VGoff += Vm
    c = m.cov
    c === nothing || _setstate!(c, m.idx, QTLSim.hap_mutate(_state(c, m.idx), Vm, c.ploidy))
    m
end

function moments_swap!(m::Moments)
    m.zbar = m.zoff; m.VG = m.VGoff; m.VE = QTLSim.envVar(m.p); m.N = m.Noff; m.gen += 1
    m.cov === nothing || m.extinct || (m.VG = _vpop(m.cov, m.idx, m.N))
    m
end

isextinct(m::Moments) = m.extinct

# QTSim's fitness: pdf(Normal(0, ω), x)/(0.4/ω) = c·exp(−x²/2ω²), with c = 1/(0.4√(2π))
const FITNESS_PEAK = 1 / (0.4 * sqrt(2π))

# A fluctuating move lowers to its expected optimum, `QTLSim.fluctuating_optimum` with no
# noise: the periodic move exactly, red noise as its trend. `θ`, when given, is the optimum
# instead (a realised sequence, `run_moments(p; optimum)`).
function moments_select!(m::Moments; θ = nothing)
    p = m.p
    if θ !== nothing
        m.θ = θ
    elseif p.MoveFitnessSDGen != -1 && m.gen >= p.MoveFitnessSDGen && p.MoveType in (:red_noise, :periodic)
        t = m.gen - p.MoveFitnessSDGen
        t == 0 && (m.θ0 = m.θ)
        m.θ = QTLSim.fluctuating_optimum(p, m.θ0, t)
    elseif p.MoveFitnessSDGen != -1 && m.gen >= p.MoveFitnessSDGen     # QTSim SetFitnessMean
        m.θ = p.MoveType === :phenotype ? m.zbar + sqrt(p.PhenotypeVar) * p.MoveFitnessSD :
              p.MoveType === :fitness ? m.zbar + m.ω * p.MoveFitnessSD : m.θ + p.MoveFitnessSD
    end
    if p.Drift
        m.Wbar, m.S, m.zs, m.VGs = 1.0, 0.0, m.zbar, m.VG
        m.W2bar, m.sinv = 1.0, 0.0
    else
        # Gaussian selection of width ωe on z ~ N(z̄, V_P): survivors are N(z̄ + V_P λ/s, V_P ωe²/s)
        # with λ = θ − z̄ and s = ωe² + V_P, and the genotype regresses on the phenotype
        ωe = p.StandardisePhenotype ? m.ω * sqrt(m.VP) : m.ω
        s = ωe^2 + m.VP
        λ = m.θ - m.zbar
        m.Wbar = FITNESS_PEAK * ωe / sqrt(s) * exp(-λ^2 / (2s))
        m.S = m.VP * λ / s
        m.zs = m.zbar + m.VG * λ / s
        m.VGs = m.VG - m.VG^2 / s
        m.sinv = 1 / s
        # E[W²]: W² is c² times a Gaussian of width ωe/√2
        s2 = ωe^2 + 2m.VP
        m.W2bar = FITNESS_PEAK^2 * ωe / sqrt(s2) * exp(-λ^2 / s2)
    end
    # Fecundity selection (Selection = :fecundity) weights parents by W instead of culling
    # by it: the same selection differential and parents' moments, as the parents'
    # distribution is W(z) p(z)/W̄ either way, but everyone survives
    m.weighted = p.Selection === :fecundity
    m.nsurv = m.weighted ? m.N : m.N * m.Wbar
    m.extinct = m.weighted ? m.N < p.Ploidy : m.nsurv < (p.Ploidy == 2 ? 2 : 1)
    m.cov === nothing || cov_select!(m)
    m
end

# Offspring inherit the survivors' mean. Their genetic variance comes from the haplotype
# moments (`cov_reproduce!`; the multilocus model), or with `loci = :aggregate` it is the
# survivors', all of selection's reduction inherited (exact for one haploid locus).
function moments_reproduce!(m::Moments)
    p = m.p
    m.zoff, m.VGoff = m.zs, m.VGs
    # optional drift: offspring are 2N gametes from the expected nsurv survivors, mated at
    # random, so two gametes share a parent with probability 1/nsurv and a locus is
    # identical by descent with probability 1/(2 nsurv). Under fecundity selection two
    # gametes share a parent with probability Σ(W_i/ΣW)² = E[W²]/(N W̄²), so the number of
    # parents is N W̄²/E[W²], at least N W̄ when W ≤ 1: fecundity selection drifts less.
    npar = m.weighted ? m.N * m.Wbar^2 / m.W2bar : m.nsurv
    m.cov === nothing ? (m.erosion && (m.VGoff *= 1 - 1 / (2 * max(npar, 1.0)))) : cov_reproduce!(m, npar)
    # density regulation counts the expected survivors N W̄ or, under fecundity selection,
    # the sum of the fitnesses, whose expectation is the same (QTLSim.noffspring)
    R, K = p.GrowthRate, p.PopSize
    n = m.weighted ? m.N * m.Wbar : m.nsurv
    # the expected number of offspring (the mean of QTLSim.noffspring), at most K
    m.Noff = p.Demography === :fixed ? K : min(K, R * n / (1 + (R - 1) * n / K))
    m.Noff < 1 && (m.extinct = true)
    m
end

function moments_record!(rec, m::Moments)
    (rec === nothing || m.gen % m.p.Sample != 0) && return m
    push!(rec.gen, m.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[m])
    push!(rec.pheno_mean, m.zbar); push!(rec.pheno_var, m.VP)
    push!(rec.geno_var, m.cov === nothing || !m.cov.drift ? m.VG : m.VG * m.N / (m.N - 1)); push!(rec.nsurv, round(Int, m.nsurv)); push!(rec.theta, m.θ)
    push!(rec.omega, m.ω); push!(rec.mean_fitness, m.Wbar); push!(rec.sel_diff, m.S)
    push!(rec.popsize, round(Int, m.N))
    m
end

"""
    moment_images(model, rec = nothing)

The image of each step of the life cycle in the moments lowering of `model`, by name: for
QTSim's one trait this table, for `TraitModel(n)` the table of `n` traits
(`trait_moment_images`, trait_moments.jl).
"""
moment_images(rec = nothing) = moment_images(OneTraitModel(), rec)
moment_images(::OneTraitModel, rec = nothing) = Dict{Symbol,Any}(
    :develop => Operator(:develop, (m, rng) -> (m.VP = m.VG + m.VE; m)),
    :select => Operator(:select, (m, rng) -> moments_select!(m)),
    :reproduce => Operator(:reproduce, (m, rng) -> m.extinct ? m : moments_reproduce!(m)),
    :record => Operator(:record, (m, rng) -> moments_record!(rec, m)),
    :rescale => Operator(:rescale, (m, rng) -> (m.gen == m.p.ChangeFitnessGen && (m.ω = m.p.ChangeFitnessSD); m)),
    :mutate => Operator(:mutate, (m, rng) -> moments_mutate!(m)),
    :swap => Operator(:swap, (m, rng) -> moments_swap!(m)))

function run_moments1(p::Params; cycle = lifecycle().generation, images = Dict{Symbol,Any}(),
                      zbar = 0.0, VG = 0.0, optimum = nothing, erosion = false, between = false, origins = :exchange,
                      loci::Symbol = :auto)
    QTLSim.checkparams(p)
    check_loci(loci, between)
    loci, α2 = resolve_loci([p], cycle, loci)
    between && check_between([p], origins)
    if optimum !== nothing          # a given optimum per generation, optimum[t + 1], replaces the move
        length(optimum) >= p.NoGenerations ||
            throw(ArgumentError("optimum gives $(length(optimum)) generations; NoGenerations is $(p.NoGenerations)"))
        images = merge(Dict{Symbol,Any}(:select => Operator(:select, (m, rng) -> moments_select!(m; θ = optimum[m.gen + 1]))),
                       Dict{Symbol,Any}(images))
    end
    zygotes = LifeCycle[:Zygotes]
    dom(cycle) == zygotes && codom(cycle) == zygotes ||
        throw(ArgumentError("the moments lowering iterates a life cycle Zygotes → Zygotes; got $(dom(cycle)) → $(codom(cycle))"))
    rec = records()
    step = lower(cycle, :moments; images, record = rec)
    cov = loci === :aggregate ? nothing : DemeCov(p, [float(p.PopSize)], [float(VG)]; α2, origins, drift = between || erosion)
    m = Moments(p; zbar, VG, erosion, cov)
    mp = metapop_records()
    m.gen = 0                    # whole cycles, as in run_composed
    while m.gen < p.NoGenerations
        it = m.gen
        between && it % p.Sample == 0 && record_metapop!(mp, it, [m])
        step(m, nothing)
        m.extinct && return between ? (extinct = true, extinct_gen = m.gen, rec..., metapop = mp) :
                                      (extinct = true, extinct_gen = m.gen, rec...)
        m.gen = max(m.gen, it + 1)
    end
    between ? (extinct = false, extinct_gen = -1, rec..., metapop = mp) : (extinct = false, extinct_gen = -1, rec...)
end

# --- several demes ----------------------------------------------------------

"The moments of the demes of a metapopulation, which the image of a morphism on several demes (a `DemeOp`) updates."
struct DemeMoments
    demes::Vector{Moments}
end

isextinct(s::DemeMoments) = all(isextinct, s.demes)

id(s::Stage{:moments}) = nfactors(s) == 1 ? Operator(:id, (m, rng) -> m) :
                         DemeOp(nfactors(s), collect(1:nfactors(s)), (st, d, rng) -> st)
braid(a::Stage{:moments}, b::Stage{:moments}) =
    DemeOp(nfactors(a) + nfactors(b), braid_perm(nfactors(a), nfactors(b)), (st, d, rng) -> st)

"""
    moments_migrate!(s, demes, rate)

The image of `migrate` in the moments lowering: the demes `s.demes[demes]` mix as the
kernel's migration mixes them, at `rate` (a number: `rate/(k − 1)` toward each other
deme; or a matrix), with expected composition `P` from `migration_mixing`. A deme becomes
a mixture of the demes it received from: its mean genotype is `Σ_c P[a, c] z̄_c`, its
genetic variance `Σ_c P[a, c] V_c + Σ_c P[a, c] (z̄_c − z̄_a)²`, the variance within the
demes plus that between their means, and its environmental variance `Σ_c P[a, c] V_E,c`,
since migrants keep the environmental deviations they were born with. These are the
mixture's first and second moments for the expected composition; the sampling variance
of the numbers of migrants, of order `rate/N`, is left out, as drift is in the rest of
the lowering. The mixture is not Gaussian, which the next selection assumes it is.
"""
function moments_migrate!(s, demes, rate)
    ms = s.demes[demes]
    k = length(ms)
    M = rate isa AbstractMatrix ? QTLSim.check_migration(rate, k) :
        [a == b ? 0.0 : Float64(rate) / (k - 1) for a in 1:k, b in 1:k]
    alive = [!m.extinct for m in ms]
    ms[1].cov === nothing || return cov_migrate!(s, ms, M, alive)
    P = migration_mixing(M, [m.N for m in ms], alive)
    z, V, VE = [m.zbar for m in ms], [m.VG for m in ms], [m.VE for m in ms]
    for a in 1:k
        alive[a] || continue
        w = P[a, :]
        za = sum(w .* z)
        ms[a].zbar = za
        ms[a].VG = sum(w .* V) + sum(w .* (z .- za) .^ 2)
        ms[a].VE = sum(w .* VE)
    end
    s
end

# The backward matrix P and the law J of the origins of pairs. Without drift P is
# `migration_mixing`, as for `loci = :aggregate`, and J = P ⊗ P. With drift P is
# `backward_migration` (`migration_mixing` averaged over the visiting orders; exactly
# symmetric for the island model) and J is `exchange_origins(P, N)` (the kernel's
# `migrate!`), `kron(P, P)` (`:pool`) or a given matrix. Cached per (M, N, alive).
function cov_mixing(c::DemeCov, M, Nv, alive)
    get!(c.cache, (M, Nv, alive)) do
        if !c.drift
            P = migration_mixing(M, Nv, alive)
            return (P, kron(P, P))
        end
        P = all(alive) ? backward_migration(M, Nv) : migration_mixing(M, Nv, alive)
        J = c.origins === :exchange ? exchange_origins(P, Nv) : c.origins === :pool ? kron(P, P) : c.origins
        size(J) == (length(Nv)^2, length(Nv)^2) ||
            throw(ArgumentError("origins must be $(length(Nv)^2) × $(length(Nv)^2); got $(size(J))"))
        (P, J)
    end
end

# migration of the means and of the second moments (see the comment on DemeCov)
function cov_migrate!(s, ms, M, alive)
    k = length(ms)
    c = ms[1].cov
    Nv = [m.N for m in ms]
    P, J = cov_mixing(c, M, Nv, alive)
    ix = [m.idx for m in ms]
    z, VE = [m.zbar for m in ms], [m.VE for m in ms]
    z2 = [alive[a] ? sum(P[a, :] .* z) : z[a] for a in 1:k]
    Δ = [(z[cc] - z2[a]) / c.ploidy for a in 1:k, cc in 1:k]
    function mix(T, R, Q, κ)
        T2, R2, Q2 = copy(T), copy(R), copy(Q)
        for a in 1:k
            alive[a] || continue
            T2[a] = sum(P[a, cc] * (T[cc] + κ * Δ[a, cc]^2) for cc in 1:k)
            R2[a] = sum(P[a, cc] * (R[cc] + κ * Δ[a, cc]^2) for cc in 1:k)
            for b in 1:k
                alive[b] || continue
                row = a + (b - 1) * k
                Q2[a, b] = sum(J[row, cc + (d - 1) * k] * (Q[cc, d] + κ * Δ[a, cc] * Δ[b, d]) for cc in 1:k, d in 1:k)
            end
        end
        T2, R2, (Q2 .+ Q2') ./ 2
    end
    for j in eachindex(c.ω)        # each pair class carries its share ω of the difference of the means
        c.T[ix, j], c.R[ix, j], c.Q[ix, ix, j] = mix(c.T[ix, j], c.R[ix, j], c.Q[ix, ix, j], c.ω[j])
    end
    for a in 1:k
        alive[a] || continue
        ms[a].zbar = z2[a]
        ms[a].VE = sum(P[a, :] .* VE)
        ms[a].VG = _vpop(c, ms[a].idx, max(ms[a].N, 2.0))
    end
    s
end

# the census statistics of `between = true`, over the living demes: `within` the expected
# within-deme sample variances, `pairs[a, b]` half the expected squared difference of the
# deme mean genotypes (the deterministic difference of the means included), and from them
# the statistics of `neutral_deme_variances`
function cov_summary(ms::AbstractVector{Moments})
    ms = [m for m in ms if !m.extinct]
    k = length(ms)
    k == 0 && return nothing
    c = ms[1].cov
    Nv = [max(m.N, 2.0) for m in ms]
    within = [_vpop(c, m.idx, N) * N / (N - 1) for (m, N) in zip(ms, Nv)]
    C = [_cmean(c, m.idx, N) for (m, N) in zip(ms, Nv)]
    pairs = [a == b ? 0.0 : (ms[a].zbar - ms[b].zbar)^2 / 2 + c.ploidy^2 / 2 * (C[a] + C[b] - 2sum(c.Q[ms[a].idx, ms[b].idx, :]))
             for a in 1:k, b in 1:k]
    VW = sum(within) / k
    k == 1 && return (; within, pairs, VW, among = 0.0, between = 0.0, qst = 0.0, qst_demes = 0.0)
    among = sum(pairs) / (k * (k - 1))
    between = among - sum(within ./ Nv) / k
    σ2 = (k - 1) / k * between
    (; within, pairs, VW, among, between, qst = σ2 / (σ2 + 2VW), qst_demes = between / (between + 2VW))
end

metapop_records() = (gen = Int[], within = Float64[], var_among = Float64[], between = Float64[], qst = Float64[],
                     qst_demes = Float64[], pairs = Matrix{Float64}[])

function record_metapop!(mp, gen, ms)
    x = cov_summary(ms)
    x === nothing && return mp
    push!(mp.gen, gen); push!(mp.within, x.VW); push!(mp.var_among, x.among); push!(mp.between, x.between)
    push!(mp.qst, x.qst); push!(mp.qst_demes, x.qst_demes); push!(mp.pairs, x.pairs)
    mp
end

function check_between(ps, origins)
    all(q -> q.Ploidy == 2, ps) || throw(ArgumentError("between = true needs diploid demes (Ploidy = 2)"))
    origins isa AbstractMatrix || origins in (:exchange, :pool) ||
        throw(ArgumentError("origins is :exchange, :pool or a k² × k² matrix; got $origins"))
end

function check_loci(loci, between)
    loci in (:auto, :multilocus, :rare_alleles, :aggregate) ||
        throw(ArgumentError("loci is :auto, :multilocus, :rare_alleles or :aggregate; got :$loci"))
    loci === :aggregate && between &&
        throw(ArgumentError("between = true tracks the haplotype moments of the multilocus model; it has no loci = :aggregate"))
end

"Migration among `k` demes at `rate`, as the image of the step `migrate` in the moments lowering."
moments_migrate_op(k, rate) = DemeOp(k, collect(1:k), (s, d, rng) -> moments_migrate!(s, d, rate))

# The mutation fields (`MutationScalar`, `MutationRatePerGamete`) that the mutate steps of
# `cycle` carry: a named tuple when every mutate step carries the same (empty when none
# does), `nothing` when they differ.
function mutate_fields(cycle)
    fs = unique([NamedTuple(k => v for (k, v) in pairs(stepparams(g)) if k in (:MutationScalar, :MutationRatePerGamete))
                 for g in stepsof(cycle) if !(g isa ObExpr) && stepname(g) === :mutate])
    isempty(fs) ? (;) : length(fs) == 1 ? only(fs) : nothing
end

# `loci` resolved for the demes `ps` and the life cycle `cycle`, with the rare-alleles
# effect variance α² of each deme (zeros unless :rare_alleles). `:auto` is
# `QTLSim.moments_loci` of every deme (with the mutate steps' fields), :rare_alleles if all
# of them are, and :multilocus if the mutate steps carry different mutation fields (then no
# one α² describes the mutations; an explicit :rare_alleles is an error there).
function resolve_loci(ps::AbstractVector{Params}, cycle, loci::Symbol)
    loci in (:multilocus, :aggregate) && return loci, zeros(length(ps))
    mf = mutate_fields(cycle)
    if mf === nothing
        loci === :auto && return :multilocus, zeros(length(ps))
        throw(ArgumentError("loci = :rare_alleles needs one law of mutation effects per deme; the mutate steps of the " *
                            "life cycle carry different MutationScalar or MutationRatePerGamete"))
    end
    qs = isempty(mf) ? ps : [Params(q; mf...) for q in ps]
    if loci === :auto
        all(q -> QTLSim.moments_loci(q) === :rare_alleles, qs) || return :multilocus, zeros(length(ps))
        loci = :rare_alleles
    end
    loci, [QTLSim.rare_alleles_effect(q, loci) for q in qs]
end

function run_moments(p::Params; cycle = lifecycle().generation, kw...)
    p.NoTraits > 1 && return run_trait_moments(p; cycle, kw...)    # trait_moments.jl
    k = nfactors(dom(cycle))
    k > 1 ? run_moments(fill(p, k); cycle, kw...) : run_moments1(p; cycle, kw...)
end

function run_moments(ps::AbstractVector{Params}; cycle, images = Dict{Symbol,Any}(), zbar = 0.0, VG = 0.0,
                     erosion = false, migration_rate::Union{Real,AbstractMatrix} = 0.0, between = false,
                     origins::Union{Symbol,AbstractMatrix} = :exchange, loci::Symbol = :auto)
    ps[1].NoTraits > 1 && return run_trait_moments(collect(Params, ps); cycle, images, erosion, migration_rate, between, origins,
                                                   loci,(zbar isa Real ? (;) : (; zbar))..., (VG isa Real ? (;) : (; G = VG))...)  # trait_moments.jl
    k = length(ps)
    k == 1 && return between ? run_moments1(ps[1]; cycle, images, zbar = only(vcat(zbar)), VG = only(vcat(VG)), erosion, between, origins, loci) :
                               run_moments1(ps[1]; cycle, images, zbar, VG, erosion, loci)
    foreach(QTLSim.checkparams, ps)
    check_loci(loci, between)
    between && check_between(ps, origins)
    loci, α2 = resolve_loci(collect(Params, ps), cycle, loci)
    loci !== :aggregate && any(q -> (q.Ploidy, q.NoQTL, q.RecombPoissonMean, q.CrossoverModel) !=
                                     (ps[1].Ploidy, ps[1].NoQTL, ps[1].RecombPoissonMean, ps[1].CrossoverModel), ps) &&
        throw(ArgumentError("the demes of a metapopulation share Ploidy, NoQTL, RecombPoissonMean and CrossoverModel"))
    all(q -> q.NoGenerations == ps[1].NoGenerations, ps) ||
        throw(ArgumentError("the demes of a metapopulation run for the same number of generations"))
    zk = zygotes_word(k)
    dom(cycle) == zk && codom(cycle) == zk ||
        throw(ArgumentError("the moments of $k demes iterate a life cycle on $k zygote stages, Zygotes⊗…⊗Zygotes; " *
                            "got $(dom(cycle)) → $(codom(cycle))"))
    rec = (gen = Int[], deme = Int[], pheno_mean = Float64[], pheno_var = Float64[], geno_var = Float64[],
           nsurv = Int[], theta = Float64[], omega = Float64[], mean_fitness = Float64[], sel_diff = Float64[],
           popsize = Int[], demeof = IdDict{Moments,Int}())
    perdeme(x, name) = x isa AbstractVector ? (length(x) == k ? x : throw(ArgumentError("$name gives $(length(x)) values for $k demes"))) : fill(x, k)
    zs, vs = perdeme(zbar, "zbar"), perdeme(VG, "VG")
    cov = loci === :aggregate ? nothing :
          DemeCov(ps[1], [float(q.PopSize) for q in ps], float.(vs); α2, origins, drift = between || erosion)
    s = DemeMoments([Moments(q; zbar = z, VG = v, erosion, cov, idx = i) for (i, (q, z, v)) in enumerate(zip(ps, zs, vs))])
    mp = metapop_records()
    for (i, m) in enumerate(s.demes); rec.demeof[m] = i; end
    step = demeop(lower(cycle, :moments; images, record = rec, migration_rate))
    d = collect(1:k)
    extinct_gen = fill(-1, k)
    gen = 0
    while gen < ps[1].NoGenerations
        foreach(m -> m.extinct || (m.gen = gen), s.demes)
        between && gen % ps[1].Sample == 0 && record_metapop!(mp, gen, s.demes)
        step.act(s, d, nothing)
        for (i, m) in enumerate(s.demes)
            m.extinct && extinct_gen[i] < 0 && (extinct_gen[i] = m.gen)
        end
        isextinct(s) && break
        d = d[step.perm]
        gen = max(gen + 1, maximum(m.gen for m in s.demes if !m.extinct))
    end
    out = (extinct = isextinct(s), extinct_gen, gen = rec.gen, deme = rec.deme, pheno_mean = rec.pheno_mean,
           pheno_var = rec.pheno_var, geno_var = rec.geno_var, nsurv = rec.nsurv, theta = rec.theta, omega = rec.omega,
           mean_fitness = rec.mean_fitness, sel_diff = rec.sel_diff, popsize = rec.popsize)
    between ? (out..., metapop = mp) : out
end
