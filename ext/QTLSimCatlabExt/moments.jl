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
    erosion::Bool        # drift erodes V_G by V_G/(2N_e) per generation, N_e = expected survivors
end

Moments(p::Params; zbar = 0.0, VG = 0.0, erosion = false) =
    Moments(p, 0, zbar, VG, VG + QTLSim.envVar(p), 0.0, QTLSim.initFitSD(p), p.PopSize, 1.0, 0.0,
            zbar, VG, zbar, VG, false, p.PopSize, p.PopSize, erosion)

isextinct(m::Moments) = m.extinct

# QTSim's fitness: pdf(Normal(0, ω), x)/(0.4/ω) = c·exp(−x²/2ω²), with c = 1/(0.4√(2π))
const FITNESS_PEAK = 1 / (0.4 * sqrt(2π))

function moments_select!(m::Moments)
    p = m.p
    if p.MoveFitnessSDGen != -1 && m.gen >= p.MoveFitnessSDGen     # QTSim SetFitnessMean
        m.θ = p.MoveType === :phenotype ? m.zbar + sqrt(p.PhenotypeVar) * p.MoveFitnessSD :
              p.MoveType === :fitness ? m.zbar + m.ω * p.MoveFitnessSD : m.θ + p.MoveFitnessSD
    end
    if p.Drift
        m.Wbar, m.S, m.zs, m.VGs = 1.0, 0.0, m.zbar, m.VG
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
    end
    m.nsurv = m.N * m.Wbar
    m.extinct = m.nsurv < (p.Ploidy == 2 ? 2 : 1)
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
    # identical by descent with probability 1/(2 nsurv)
    m.erosion && (m.VGoff *= 1 - 1 / (2 * max(m.nsurv, 1.0)))
    R, K, n = p.GrowthRate, p.PopSize, m.nsurv
    # the expected number of offspring (the mean of QTLSim.noffspring), at most K
    m.Noff = p.Demography === :fixed ? K : min(K, R * n / (1 + (R - 1) * n / K))
    m.Noff < 1 && (m.extinct = true)
    m
end

function moments_record!(rec, m::Moments)
    (rec === nothing || m.gen % m.p.Sample != 0) && return m
    push!(rec.gen, m.gen); push!(rec.pheno_mean, m.zbar); push!(rec.pheno_var, m.VP)
    push!(rec.geno_var, m.VG); push!(rec.nsurv, round(Int, m.nsurv)); push!(rec.theta, m.θ)
    push!(rec.omega, m.ω); push!(rec.mean_fitness, m.Wbar); push!(rec.sel_diff, m.S)
    push!(rec.popsize, round(Int, m.N))
    m
end

"The image of each step of the life cycle in the moments lowering, by name."
moment_images(rec = nothing) = Dict{Symbol,Any}(
    :develop => Operator(:develop, (m, rng) -> (m.VP = m.VG + QTLSim.envVar(m.p); m)),
    :select => Operator(:select, (m, rng) -> moments_select!(m)),
    :reproduce => Operator(:reproduce, (m, rng) -> m.extinct ? m : moments_reproduce!(m)),
    :record => Operator(:record, (m, rng) -> moments_record!(rec, m)),
    :rescale => Operator(:rescale, (m, rng) -> (m.gen == m.p.ChangeFitnessGen && (m.ω = m.p.ChangeFitnessSD); m)),
    :mutate => Operator(:mutate, (m, rng) -> (m.VGoff += QTLSim.muteVar(m.p); m)),
    :swap => Operator(:swap, (m, rng) -> (m.zbar = m.zoff; m.VG = m.VGoff; m.N = m.Noff; m.gen += 1; m)))

id(s::Stage{:moments}) = nfactors(s) == 1 ? Operator(:id, (m, rng) -> m) :
    throw(ArgumentError("the moments lowering covers one deme"))

function run_moments(p::Params; cycle = lifecycle().generation, images = Dict{Symbol,Any}(),
                     zbar = 0.0, VG = 0.0, erosion = false)
    QTLSim.checkparams(p)
    zygotes = LifeCycle[:Zygotes]
    dom(cycle) == zygotes && codom(cycle) == zygotes ||
        throw(ArgumentError("the moments lowering iterates a life cycle Zygotes → Zygotes; got $(dom(cycle)) → $(codom(cycle))"))
    rec = records()
    step = lower(cycle, :moments; images, record = rec)
    m = Moments(p; zbar, VG, erosion)
    for it in 0:p.NoGenerations-1
        m.gen = it
        step(m, nothing)
        m.extinct && return (extinct = true, extinct_gen = it, rec...)
    end
    (extinct = false, extinct_gen = -1, rec...)
end
