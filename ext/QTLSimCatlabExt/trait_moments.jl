# The moments lowering for several traits: a second table of images for the same steps
# of the same life cycle, on the moments of a multivariate Gaussian population (Lande
# 1979, 1980). It is kept apart from moments.jl, the one-trait table, to show that only
# the images differ: `run_trait_moments` lowers the unchanged `lifecycle()` with
# `images = trait_moment_images(rec)`. With one trait its recursion is moments.jl's; the
# fixed points are QTLSim's `gaussian_G` and `steady_lag_vector` (src/traits_theory.jl).

"""
The state of the trait moments lowering: as `Moments`, with vectors for the means and the
optimum and matrices for the variances. `M` is the mutational covariance, `S` the fitness
function's shape (`Ω = ω² S`) and `u` the direction the optimum moves in.
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
    VE::Float64
    weighted::Bool
    W2bar::Float64
    M::Matrix{Float64}
    shape::Matrix{Float64}
    u::Vector{Float64}
end

function TraitMoments(p::Params; zbar = zeros(p.NoTraits), G = zeros(p.NoTraits, p.NoTraits), erosion = false)
    n, VE = p.NoTraits, QTLSim.envVar(p)
    Id = Matrix{Float64}(I, n, n)
    TraitMoments(p, 0, collect(Float64, zbar), Matrix{Float64}(G), G + VE * Id, zeros(n), QTLSim.initFitSD(p),
                 p.PopSize, 1.0, zeros(n), collect(Float64, zbar), Matrix{Float64}(G), collect(Float64, zbar),
                 Matrix{Float64}(G), false, p.PopSize, p.PopSize, erosion, VE, false, 1.0,
                 QTLSim.mutational_covariance(p), QTLSim.selection_shape(p), QTLSim.move_direction(p))
end

isextinct(m::TraitMoments) = m.extinct

sym(A) = (A + A') / 2

function trait_moments_select!(m::TraitMoments)
    p = m.p
    if p.MoveFitnessSDGen != -1 && m.gen >= p.MoveFitnessSDGen       # moveoptimum_traits!
        d = p.MoveFitnessSD
        m.θ = p.MoveType === :continuous ? m.θ .+ d .* m.u :
              m.zbar .+ (p.MoveType === :phenotype ? sqrt(p.PhenotypeVar) : m.ω) * d .* m.u
    end
    n = length(m.zbar)
    if p.Drift
        m.Wbar, m.W2bar, m.S, m.zs, m.Gs = 1.0, 1.0, zeros(n), copy(m.zbar), copy(m.G)
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
    end
    m.weighted = p.Selection === :fecundity
    m.nsurv = m.weighted ? m.N : m.N * m.Wbar
    m.extinct = m.weighted ? m.N < p.Ploidy : m.nsurv < (p.Ploidy == 2 ? 2 : 1)
    m
end

# as moments_reproduce!: offspring inherit the parents' moments, with optional drift
function trait_moments_reproduce!(m::TraitMoments)
    p = m.p
    m.zoff, m.Goff = copy(m.zs), copy(m.Gs)
    npar = m.weighted ? m.N * m.Wbar^2 / m.W2bar : m.nsurv
    m.erosion && (m.Goff .*= 1 - 1 / (2 * max(npar, 1.0)))
    R, K = p.GrowthRate, p.PopSize
    nw = m.weighted ? m.N * m.Wbar : m.nsurv
    m.Noff = p.Demography === :fixed ? K : min(K, R * nw / (1 + (R - 1) * nw / K))
    m.Noff < 1 && (m.extinct = true)
    m
end

function trait_moments_record!(rec, m::TraitMoments)
    (rec === nothing || m.gen % m.p.Sample != 0) && return m
    push!(rec.gen, m.gen)
    push!(rec.pheno_mean, m.zbar[1]); push!(rec.pheno_var, m.P[1, 1]); push!(rec.geno_var, m.G[1, 1])
    push!(rec.nsurv, round(Int, m.nsurv)); push!(rec.theta, m.θ[1]); push!(rec.omega, m.ω)
    push!(rec.mean_fitness, m.Wbar); push!(rec.sel_diff, m.S[1]); push!(rec.popsize, round(Int, m.N))
    push!(rec.trait_mean, copy(m.zbar)); push!(rec.gmatrix, copy(m.G)); push!(rec.pmatrix, copy(m.P))
    push!(rec.trait_lag, m.θ - m.zbar); push!(rec.trait_sel_diff, copy(m.S))
    m
end

"""
    trait_moment_images(rec = nothing)

The image of each step of the life cycle in the moments lowering of several traits, by
name: the multivariate Gaussian recursion, `Δz̄ = G (Ω + P)⁻¹ (θ − z̄)`,
`G' = G − G (Ω + P)⁻¹ G + M`. A second table for the same generators as `moment_images`.
"""
trait_moment_images(rec = nothing) = Dict{Symbol,Any}(
    :develop => Operator(:develop, (m, rng) -> (m.P = m.G + m.VE * Matrix{Float64}(I, size(m.G)...); m)),
    :select => Operator(:select, (m, rng) -> trait_moments_select!(m)),
    :reproduce => Operator(:reproduce, (m, rng) -> m.extinct ? m : trait_moments_reproduce!(m)),
    :record => Operator(:record, (m, rng) -> trait_moments_record!(rec, m)),
    :rescale => Operator(:rescale, (m, rng) -> (m.gen == m.p.ChangeFitnessGen && (m.ω = m.p.ChangeFitnessSD); m)),
    :mutate => Operator(:mutate, (m, rng) -> (m.Goff = m.Goff + m.M; m)),
    :swap => Operator(:swap, (m, rng) -> (m.zbar = m.zoff; m.G = m.Goff; m.N = m.Noff; m.gen += 1; m)))

"""
    run_trait_moments(p; cycle = lifecycle().generation, images, zbar, G, erosion = false)

The moments lowering of a model of several traits: `cycle`, one deme's life cycle, lowered
with the images `trait_moment_images` and iterated `NoGenerations` times from mean `zbar`
and genetic covariance `G` (zero by default). `run_moments(p)` calls it when
`p.NoTraits > 1`; called directly with one trait it runs the multivariate recursion with
1×1 matrices, which is the one-trait recursion.
"""
function run_trait_moments(p::Params; cycle = lifecycle().generation, images = Dict{Symbol,Any}(),
                           zbar = zeros(p.NoTraits), G = zeros(p.NoTraits, p.NoTraits), erosion = false,
                           optimum = nothing)
    QTLSim.checkparams(p)
    optimum === nothing || throw(ArgumentError("a given optimum sequence is for one trait"))
    zygotes = LifeCycle[:Zygotes]
    dom(cycle) == zygotes && codom(cycle) == zygotes ||
        throw(ArgumentError("the moments of several traits are for one deme, a life cycle Zygotes → Zygotes; " *
                            "got $(dom(cycle)) → $(codom(cycle))"))
    rec = merge(records(), QTLSim.trait_records())
    step = lower(cycle, :moments; images = merge(trait_moment_images(rec), Dict{Symbol,Any}(images)), record = rec)
    m = TraitMoments(p; zbar, G, erosion)
    while m.gen < p.NoGenerations
        it = m.gen
        step(m, nothing)
        m.extinct && return (extinct = true, extinct_gen = m.gen, rec...)
        m.gen = max(m.gen, it + 1)
    end
    (extinct = false, extinct_gen = -1, rec...)
end
hastraits(::TraitMoments) = true
