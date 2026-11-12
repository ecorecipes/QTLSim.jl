# Phenotypic plasticity: a linear reaction norm (Chevin, Lande and Mace 2010; Lande 2009).
#
# The environment ε is measured by the optimum it sets, θ = B·ε (B = `OptimumSlope`, the
# environmental sensitivity of selection), and an individual's phenotype is
#
#     z = a + b·ε_cue + e,     ε_cue = θ_{t−τ}/B,
#
# with elevation a (the genotype at the reference environment ε = 0, where the optimum is
# the origin of the trait's scale), slope b, and a cue that is the environment τ = `CueLag`
# generations before selection. Two models:
#
# - a fixed slope (`Plasticity = b`, one trait): every individual has the same slope, so the
#   plastic part (b/B)·θ_{t−τ} is common to all, a = g. A parameter of `select`, where the
#   optimum of the generation is set (QTSim's order is develop, then the optimum moves, then
#   fitness), expressed by `express_plasticity!`.
# - a heritable reaction norm (`ReactionNorm = true`, `NoTraits = 2`): the two heritable
#   values of the model of several traits are the elevation, a = g₁, and the deviation of the
#   slope from `Plasticity`, b = Plasticity + g₂. Selection acts on the one phenotype z, and
#   the cost on the individual's slope. The kernel image of `select` is `select_norm!`; the
#   other steps are those of `TraitModel(2)`. See docs/src/plasticity_design.md.
#
# The cost of plasticity is stabilising selection on the slope around 0: fitness is
# multiplied by exp(−b²/(2ω_b²)), ω_b = `PlasticityCostSD`.

"""
    QTLSim.plastic_cue!(w) -> θ_{t−τ}

Enter the optimum of the current generation in the deme's record of its last `CueLag + 1`
optima and return the one `τ = CueLag` generations ago, the cue of the reaction norm (the
environment is `θ/B`). Before the first generation the environment is taken to have been
the first one seen. Called once per generation by `select`; a second call in the same
generation replaces the entry instead of shifting the record.
"""
function plastic_cue!(w::World)
    c, n = w.cues, w.p.CueLag + 1
    if isempty(c)
        for _ in 1:n
            push!(c, w.θ)
        end
    else
        while length(c) < n               # a step with a longer lag than the deme's
            push!(c, c[end])
        end
        if w.cuegen != w.gen
            @inbounds for i in length(c):-1:2
                c[i] = c[i-1]
            end
        end
        @inbounds c[1] = w.θ
    end
    w.cuegen = w.gen
    @inbounds c[n]
end

"""
    QTLSim.express_plasticity!(w)

Add the plastic part `(b/B)·θ_{t−τ}` of the reaction norm to every phenotype (`b =
Plasticity`, `B = OptimumSlope`, `τ = CueLag`), and keep it in `w.plastic`.
[`develop!`](@ref) gives the phenotypes at the reference environment, `g + e`, and
[`select!`](@ref) calls this once the optimum of the generation is set.
"""
function express_plasticity!(w::World)
    shift = plasticity_ratio(w.p) * plastic_cue!(w)
    Δ = shift - w.plastic
    if Δ != 0
        @inbounds for i in 1:w.N
            w.phenotype[i] += Δ
        end
    end
    w.plastic = shift
    w
end

# the cost of a fixed slope: every fitness times exp(−b²/(2ω_b²))
function plasticity_cost!(w::World)
    c = plasticity_cost(w.p)
    c == 1 && return w
    @inbounds for i in 1:w.N
        w.fitness[i] *= c
    end
    w
end

"""
    QTLSim.select_norm!(w, rng)

`select` for a heritable reaction norm (`ReactionNorm = true`, `NoTraits = 2`): move the
optimum of the phenotype (a number: `w.θ`, as for one trait), express each individual's
norm, `z = g₁ + (b₀ + g₂)·ε + e₁` with `b₀ = Plasticity` and `ε = θ_{t−τ}/B`, and select on
`z` with QTSim's Gaussian fitness times the cost of the individual's slope `b = b₀ + g₂`,

    W = c·exp(−(z − θ)²/(2ω²))·exp(−b²/(2ω_b²)),   c = 1/(0.4√(2π)).

The trait columns then hold the phenotype and the slope: `Z[1, i] = z`, `Z[2, i] = b`, so
that the records' `trait_mean` is `(z̄, b̄)`, `gmatrix` the genetic covariance matrix of
elevation and slope, and `pmatrix` the covariance of phenotype and slope. The second
trait's environmental deviation is not used.
"""
function select_norm!(w::World, rng)
    p, t = w.p, w.traits::TraitState
    (p.MoveFitnessSDGen != -1 && w.gen >= p.MoveFitnessSDGen) && moveoptimum_one!(w, rng)
    t.θ[1] = w.θ; t.θ[2] = 0.0
    ε = plastic_cue!(w) / p.OptimumSlope
    b0, θ = p.Plasticity, w.θ
    iω2 = 1 / (2 * w.ω^2)
    ic = 1 / (2 * p.PlasticityCostSD^2)          # 0 without a cost
    @inbounds for i in 1:w.N
        b = b0 + t.G[2, i]
        z = t.G[1, i] + b * ε + t.E[1, i]
        t.Z[1, i] = z; t.Z[2, i] = b
        w.phenotype[i] = z
        w.fitness[i] = p.Drift ? 1.0 : TRAIT_FITNESS_PEAK * exp(-(z - θ)^2 * iω2 - b^2 * ic)
    end
    w.plastic = b0 * ε
    survive!(w, rng)
end

# --- theory -------------------------------------------------------------------------------
#
# The lag and the critical rate with a fixed slope are in theory.jl (`steady_lag`,
# `critical_rate`, `red_noise_lag_var`, `expected_lag2`, `lag_load`, each with `b`, `B` and
# the cue's delay). Below: the rate a population can sustain, from the moments recursion,
# and the response of a heritable reaction norm.

"""
    sustainable_rate(p; VG, ω = initFitSD(p), b = p.Plasticity)

The largest rate of change of the optimum (per generation, `MoveType = :continuous`) at
which the deterministic recursion of the moments lowering persists, for a population with
genetic variance `VG`: at the steady lag `λ = (1 − b/B) k (V_G + V_s)/V_G`
([`steady_lag`](@ref)) the mean fitness is

    W̄(λ) = c·exp(−b²/(2ω_b²))·ω/√s·exp(−λ²/(2s)),   s = ω² + V_G + V_E,  c = 1/(0.4√(2π)),

and the population persists while `W̄ > W_min`: `1/R` under `Demography = :beverton_holt`
(growth when rare, `R W̄ > 1`, Lean `fitness_extinct`), and `Ploidy/PopSize` at a fixed size
(QTSim dies when fewer than `Ploidy` survive; here, in expectation). So

    k_max = V_G/s · √(2 s ln(W̄(0)/W_min)) / |1 − b/B|,

0 when `W̄(0) ≤ W_min` and `Inf` at `b = B`. It is Bürger and Lynch's (1995) critical rate
with the population's own growth factor in place of `B₀` and without their allowance for
the variance of the lag ([`critical_rate`](@ref) is the form of the paper's table), and
the discrete-generation counterpart of Chevin, Lande and Mace's (2010) critical rate with
plasticity: without a cost, plasticity multiplies it by `1/|1 − b/B|`; a cost lowers
`ln(W̄(0)/W_min)` by `b²/(2ω_b²)`, so that an intermediate slope can be best.

What it is: a deterministic threshold for a population whose genetic variance stays at
`VG` whatever the lag. What it is not: the rate at which half of finite populations
survive. Under a moving optimum the lag removes breeders, drift then lowers `V_G`, and the
lag grows, a feedback this threshold leaves out (described by Bürger and Lynch 1995, who
note that their fixed-variance critical rate overestimates persistence for this reason);
in QTSim's regime the kernel's half-survival rate within 2000 generations is 0.07–0.56
times this rate (0.14–0.39 in `validate/plasticity_rates.log`; the paper's 710 populations
in vignettes 09 and 11, in `validate/critical_rate_posthoc.log`), and a cost lowers it
several times more than the fitness-margin term here says
(`validate/plasticity_erosion.log`). To first order in a log cost `c` this threshold falls
by the fraction `c/(2 ln(W̄(0)/W_min))`, the same linear form as the drift critical rate's
`c/M`, but with a margin several times larger, so it predicts too small an effect of the
cost. The feedback's critical rate is [`drift_critical_rate`](@ref), the half-survival rate
[`half_survival_rate`](@ref).
"""
function sustainable_rate(p::Params; VG::Real, ω::Real = initFitSD(p), b::Real = p.Plasticity)
    s = ω^2 + VG + envVar(p)
    W0 = TRAIT_FITNESS_PEAK * plasticity_cost(p; b) * ω / sqrt(s)
    Wmin = p.Demography === :beverton_holt ? 1 / p.GrowthRate : p.Ploidy / p.PopSize
    W0 > Wmin || return 0.0
    VG / s * sqrt(2s * log(W0 / Wmin)) / abs(1 - b / p.OptimumSlope)
end

"""
    reaction_norm_response(G, ε; ω, VE, θ, abar, bbar, ωb = Inf)

One generation of Gaussian selection on a heritable linear reaction norm (Lande 2009;
Gavrilets and Scheiner 1993): individuals have elevation `a` and slope `b`, with means
`abar`, `bbar` and additive genetic covariance matrix `G` (2×2, elevation first), phenotype
`z = a + b ε + e`, `e ~ N(0, VE)`, in the environment `ε`, fitness
`W ∝ exp(−(z − θ)²/(2ω²))·exp(−b²/(2ω_b²))`. With `v = (1, ε)`, the genotypes `(a, b)` of
the selected parents are Gaussian with

    G_s = (I + G Q)⁻¹ G,     Q = v vᵀ/(ω² + V_E) + e₂ e₂ᵀ/ω_b²,
    Δ(ā, b̄) = G_s (v (θ − z̄)/(ω² + V_E) − e₂ b̄/ω_b²),     z̄ = ā + b̄ ε.

Without a cost this is `Δ = G v (θ − z̄)/s`, `s = ω² + V_E + vᵀ G v`: the slope responds in
proportion to `ε`, so plasticity evolves only away from the reference environment, and the
genetic variance of the phenotype, `vᵀ G v = G_aa + 2 G_ab ε + G_bb ε²`, grows with `ε²`.
Returns `(; Δ, Gs, K, zbar, lag, Vz, Wbar, W2bar)`: the response, the parents' covariance,
the matrix `K = Q (I + G Q)⁻¹` with `G_s = G − G K G`, the mean phenotype, the lag `θ − z̄`,
the phenotypic variance `vᵀ G v + V_E`, the mean fitness as a fraction of the peak and the
mean squared fitness likewise.
"""
function reaction_norm_response(G::AbstractMatrix, ε::Real; ω::Real, VE::Real, θ::Real, abar::Real, bbar::Real,
                                ωb::Real = Inf)
    size(G) == (2, 2) || throw(ArgumentError("G is the 2×2 covariance matrix of elevation and slope; got size $(size(G))"))
    v = [1.0, Float64(ε)]
    c = 1 / ωb^2
    Q = (v * v') ./ (ω^2 + VE) .+ [0.0 0.0; 0.0 c]
    K = Q / (I + G * Q)
    K = (K + K') / 2
    Gs = G - G * K * G
    Gs = (Gs + Gs') / 2
    zbar = abar + bbar * ε
    lag = θ - zbar
    Δ = Gs * (v .* (lag / (ω^2 + VE)) .- [0.0, c * bbar])
    Vz = dot(v, G * v) + VE
    (; Δ, Gs, K, zbar, lag, Vz, Wbar = norm_mean_fitness(G, v, VE, lag, bbar, ω^2, ωb^2),
     W2bar = norm_mean_fitness(G, v, VE, lag, bbar, ω^2 / 2, ωb^2 / 2))
end

# E[exp(−(z−θ)²/(2ω²))·exp(−b²/(2ω_b²))] for a Gaussian population: selection on the phenotype,
# then the cost on the slopes of those selected
function norm_mean_fitness(G, v, VE, lag, bbar, ω2, ωb2)
    s = ω2 + dot(v, G * v) + VE
    W = sqrt(ω2 / s) * exp(-lag^2 / (2s))
    isfinite(ωb2) || return W
    Gv2 = G[2, 1] * v[1] + G[2, 2] * v[2]
    b1 = bbar + Gv2 * lag / s
    V1 = G[2, 2] - Gv2^2 / s
    W * sqrt(ωb2 / (ωb2 + V1)) * exp(-b1^2 / (2 * (ωb2 + V1)))
end

"""
    reaction_norm_lag(G, ε; k, B = 1, bbar, ω, VE, ωb = Inf)

The quasi-steady lag `θ − z̄` of a population with a heritable reaction norm whose optimum
moves `k` per generation (`θ = B ε`, so the environment moves `k/B`), when the environment
is `ε` and the mean slope `bbar`: the mean phenotype changes by `vᵀ Δ(ā, b̄) + b̄ k/B` per
generation (to first order in `k`), which equals `k` at

    λ = ((1 − b̄/B) k + (G_ab + G_bb ε) b̄/ω_b²) · s/(vᵀ G v),   s = ω² + V_E + vᵀ G v,

with `v = (1, ε)` (the cost term to first order in `1/ω_b²`). It generalises
[`steady_lag`](@ref), which it is for `G_bb = G_ab = 0`; it is quasi-steady because `ε`,
`b̄` and `vᵀ G v` change as the environment moves: far from the reference environment the
slope's variance dominates, `λ → (1 − b̄/B) k s/(G_bb ε²)`, and the slope itself approaches
`B` as `d b̄/dt ≈ (B − b̄)/t`.
"""
function reaction_norm_lag(G::AbstractMatrix, ε::Real; k::Real, B::Real = 1, bbar::Real, ω::Real, VE::Real,
                           ωb::Real = Inf)
    Vg = G[1, 1] + 2 * G[1, 2] * ε + G[2, 2] * ε^2
    ((1 - bbar / B) * k + (G[1, 2] + G[2, 2] * ε) * bbar / ωb^2) * (ω^2 + VE + Vg) / Vg
end
