# The moments lowering of a heritable reaction norm (`ReactionNormModel`, src/plasticity.jl):
# the table of `TraitModel(2)` (trait_moments.jl) for the pair (elevation, slope), with the
# images of `select` and `record` replaced. The pair is read through the linear map
# (a, b) ↦ a + b·ε of the generation's environment ε = θ_{t−τ}/B: selection is Gaussian on
# the one phenotype z = a + b ε + e, with the cost of plasticity on the slope
# (`QTLSim.reaction_norm_response`; Lande 2009), and acts on the haplotype moments of the
# multilocus model through the matrix K with G_s = G − G K G, as (Ω + P)⁻¹ does for two
# traits. The fixed slope of the one-trait model is in moments.jl (`moments_select!`).

norm_optimum(m::TraitMoments) = m.θ[1]
set_norm_optimum!(m::TraitMoments, θ) = (m.θ = [Float64(θ), 0.0]; m)

function norm_moments_select!(m::TraitMoments; θ = nothing)
    p = m.p
    trait_moments_optimum!(m, θ === nothing ? nothing : [Float64(θ), 0.0])
    ε = moments_cue!(m, m.θ[1]) / p.OptimumSlope
    m.env = ε
    v = [1.0, ε]
    Gv = m.G * v
    Vz = dot(v, Gv) + m.VE
    # the covariance matrix of the phenotype and the slope
    m.P = [Vz Gv[2]; Gv[2] m.G[2, 2]]
    if p.Drift
        m.Wbar, m.W2bar, m.S, m.zs, m.Gs = 1.0, 1.0, zeros(2), copy(m.zbar), copy(m.G)
        m.sinv = zeros(2, 2)
    else
        r = QTLSim.reaction_norm_response(m.G, ε; ω = m.ω, VE = m.VE, θ = m.θ[1], abar = m.zbar[1],
                                          bbar = p.Plasticity + m.zbar[2], ωb = p.PlasticityCostSD)
        m.Wbar, m.W2bar = FITNESS_PEAK * r.Wbar, FITNESS_PEAK^2 * r.W2bar
        m.zs, m.Gs, m.sinv = m.zbar + r.Δ, r.Gs, r.K
        # the differentials of the phenotype (its environmental part responds too) and of the slope
        m.S = [dot(v, r.Δ) + m.VE * r.lag / (m.ω^2 + Vz), r.Δ[2]]
    end
    trait_moments_survive!(m)
end

function norm_moments_record!(rec, m::TraitMoments)
    (rec === nothing || m.gen % m.p.Sample != 0) && return m
    push!(rec.gen, m.gen); haskey(rec, :deme) && push!(rec.deme, rec.demeof[m])
    Gr = m.cov !== nothing && m.cov.drift ? m.G .* (m.N / (m.N - 1)) : m.G
    bbar = m.p.Plasticity + m.zbar[2]
    zbar = m.zbar[1] + bbar * m.env
    push!(rec.pheno_mean, zbar); push!(rec.pheno_var, m.P[1, 1]); push!(rec.geno_var, Gr[1, 1])
    push!(rec.nsurv, round(Int, m.nsurv)); push!(rec.theta, m.θ[1]); push!(rec.omega, m.ω)
    push!(rec.mean_fitness, m.Wbar); push!(rec.sel_diff, m.S[1]); push!(rec.popsize, round(Int, m.N))
    push!(rec.trait_mean, [zbar, bbar]); push!(rec.gmatrix, copy(Gr)); push!(rec.pmatrix, copy(m.P))
    push!(rec.trait_lag, [m.θ[1] - zbar, -bbar]); push!(rec.trait_sel_diff, copy(m.S))
    m
end

"""
    moment_images(ReactionNormModel(), rec = nothing)

The moments lowering of a heritable reaction norm: `trait_moment_images` for the pair
(elevation, slope), with `select` the Gaussian selection on the phenotype `z = a + b·ε + e`
and the cost on the slope, and `record` reporting the mean phenotype and the mean slope
(`trait_mean = (z̄, b̄)`, as the kernel's records).
"""
moment_images(::ReactionNormModel, rec = nothing) = merge!(trait_moment_images(rec), Dict{Symbol,Any}(
    :select => Operator(:select, (m, rng) -> norm_moments_select!(m)),
    :record => Operator(:record, (m, rng) -> norm_moments_record!(rec, m))))

eidos_images(::ReactionNormModel) =
    throw(ArgumentError("a heritable reaction norm (ReactionNorm) has no lowering to Eidos yet; the kernel and the moments have"))
