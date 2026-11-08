# The third lowering of the life cycle: to deterministic recursions on the moments of a
# Gaussian population, the infinite-population approximation of quantitative genetics
# (Lande 1975, 1976; Bürger and Lynch 1995). The same morphism of the life-cycle category
# that the kernel and SLiM simulate is lowered to its predicted means and variances, so
# a simulation can be checked against theory step for step, and any life cycle built
# from these steps has a prediction. The fixed points are the formulas of
# src/theory.jl; the algebra is in lean/CategoricalQuantGen/Theory.lean.

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
end

Moments(p::Params; zbar = 0.0, VG = 0.0, erosion = false) =
    Moments(p, 0, zbar, VG, VG + QTLSim.envVar(p), 0.0, QTLSim.initFitSD(p), p.PopSize, 1.0, 0.0,
            zbar, VG, zbar, VG, false, p.PopSize, p.PopSize, 0.0, erosion, QTLSim.envVar(p), false, 1.0)

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
        m.W2bar = 1.0
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
    m
end

# Offspring inherit the survivors' mean and genetic variance: the Gaussian (continuum of
# alleles) approximation, in which selection's reduction of the variance is genic and
# lasting, and drift is ignored.
function moments_reproduce!(m::Moments)
    p = m.p
    m.zoff, m.VGoff = m.zs, m.VGs
    # optional drift: offspring are 2N gametes from the expected nsurv survivors, mated at
    # random, so two gametes share a parent with probability 1/nsurv and a locus is
    # identical by descent with probability 1/(2 nsurv). Under fecundity selection two
    # gametes share a parent with probability Σ(W_i/ΣW)² = E[W²]/(N W̄²), so the number of
    # parents is N W̄²/E[W²], at least N W̄ when W ≤ 1: fecundity selection drifts less.
    npar = m.weighted ? m.N * m.Wbar^2 / m.W2bar : m.nsurv
    m.erosion && (m.VGoff *= 1 - 1 / (2 * max(npar, 1.0)))
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
    push!(rec.geno_var, m.VG); push!(rec.nsurv, round(Int, m.nsurv)); push!(rec.theta, m.θ)
    push!(rec.omega, m.ω); push!(rec.mean_fitness, m.Wbar); push!(rec.sel_diff, m.S)
    push!(rec.popsize, round(Int, m.N))
    m
end

"The image of each step of the life cycle in the moments lowering, by name."
moment_images(rec = nothing) = Dict{Symbol,Any}(
    :develop => Operator(:develop, (m, rng) -> (m.VP = m.VG + m.VE; m)),
    :select => Operator(:select, (m, rng) -> moments_select!(m)),
    :reproduce => Operator(:reproduce, (m, rng) -> m.extinct ? m : moments_reproduce!(m)),
    :record => Operator(:record, (m, rng) -> moments_record!(rec, m)),
    :rescale => Operator(:rescale, (m, rng) -> (m.gen == m.p.ChangeFitnessGen && (m.ω = m.p.ChangeFitnessSD); m)),
    :mutate => Operator(:mutate, (m, rng) -> (m.VGoff += QTLSim.muteVar(m.p); m)),
    :swap => Operator(:swap, (m, rng) -> (m.zbar = m.zoff; m.VG = m.VGoff; m.VE = QTLSim.envVar(m.p); m.N = m.Noff; m.gen += 1; m)))

function run_moments1(p::Params; cycle = lifecycle().generation, images = Dict{Symbol,Any}(),
                      zbar = 0.0, VG = 0.0, optimum = nothing, erosion = false)
    QTLSim.checkparams(p)
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
    m = Moments(p; zbar, VG, erosion)
    m.gen = 0                    # whole cycles, as in run_composed
    while m.gen < p.NoGenerations
        it = m.gen
        step(m, nothing)
        m.extinct && return (extinct = true, extinct_gen = m.gen, rec...)
        m.gen = max(m.gen, it + 1)
    end
    (extinct = false, extinct_gen = -1, rec...)
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

"Migration among `k` demes at `rate`, as the image of the step `migrate` in the moments lowering."
moments_migrate_op(k, rate) = DemeOp(k, collect(1:k), (s, d, rng) -> moments_migrate!(s, d, rate))

function run_moments(p::Params; cycle = lifecycle().generation, kw...)
    p.NoTraits > 1 && return run_trait_moments(p; cycle, kw...)    # trait_moments.jl
    k = nfactors(dom(cycle))
    k > 1 ? run_moments(fill(p, k); cycle, kw...) : run_moments1(p; cycle, kw...)
end

function run_moments(ps::AbstractVector{Params}; cycle, images = Dict{Symbol,Any}(), zbar = 0.0, VG = 0.0,
                     erosion = false, migration_rate::Union{Real,AbstractMatrix} = 0.0)
    k = length(ps)
    k == 1 && return run_moments1(ps[1]; cycle, images, zbar, VG, erosion)
    foreach(QTLSim.checkparams, ps)
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
    s = DemeMoments([Moments(q; zbar = z, VG = v, erosion) for (q, z, v) in zip(ps, perdeme(zbar, "zbar"), perdeme(VG, "VG"))])
    for (i, m) in enumerate(s.demes); rec.demeof[m] = i; end
    step = demeop(lower(cycle, :moments; images, record = rec, migration_rate))
    d = collect(1:k)
    extinct_gen = fill(-1, k)
    gen = 0
    while gen < ps[1].NoGenerations
        foreach(m -> m.extinct || (m.gen = gen), s.demes)
        step.act(s, d, nothing)
        for (i, m) in enumerate(s.demes)
            m.extinct && extinct_gen[i] < 0 && (extinct_gen[i] = m.gen)
        end
        isextinct(s) && break
        d = d[step.perm]
        gen = max(gen + 1, maximum(m.gen for m in s.demes if !m.extinct))
    end
    (extinct = isextinct(s), extinct_gen, gen = rec.gen, deme = rec.deme, pheno_mean = rec.pheno_mean,
     pheno_var = rec.pheno_var, geno_var = rec.geno_var, nsurv = rec.nsurv, theta = rec.theta, omega = rec.omega,
     mean_fitness = rec.mean_fitness, sel_diff = rec.sel_diff, popsize = rec.popsize)
end
