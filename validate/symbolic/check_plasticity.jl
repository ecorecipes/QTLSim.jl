# Symbolic checks of the theory of phenotypic plasticity (src/plasticity.jl, the plasticity
# keywords of src/theory.jl, docs/src/plasticity_design.md), and numerical spot checks
# against the code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_plasticity.jl
# Writes validate/symbolic/check_plasticity.log; exits 1 if any check fails.
#
# A reaction norm z = g + b ε + e with the environment measured by the optimum, θ = B ε, a
# cue τ generations old and φ = b/B: the mean genotype follows ḡ′ = ḡ + β (θ_t − ḡ − φ θ_{t−τ}).
#
# (a) the steady lag behind θ_t = k t is (1 − φ) k/β whatever τ, and steady_lag, expected_lag2,
#     lag_load, critical_rate and red_noise_lag_var reduce exactly (===) to their forms
#     without plasticity at b = 0;
# (b) the sustainable rate: at k_max the steady lag is the lag at which W̄ = W_min, and
#     sustainable_rate, critical_rate scale by 1/|1 − φ| without a cost;
# (c) red noise: the lag L obeys L′ = a L + θ′ − θ (a = 1 − β); c₀ = Cov(θ, L), Γ(0), the
#     recursion and closed form of Γ(h), and Var(L_t − φ L_{t−τ}) = (1 + φ²)Γ(0) − 2φΓ(τ),
#     (1 − φ)²Γ(0) for τ = 0; a sine wave: |1 − φ e^{−iντ}|² = 1 + φ² − 2φ cos ντ;
# (d) a heritable reaction norm: G_s = (I + G Q)⁻¹ G = G − G K G with K = Q (I + G Q)⁻¹;
#     without a cost G_s = G − G v vᵀ G/s and Δ = G v λ/s; the fixed points in a constant
#     environment (with a cost: λ = 0 and b̄ = 0, the slope is lost; without: the line λ = 0);
#     the quasi-steady lag, and ε·Δb̄ → (1 − b̄/B) k far from the reference environment;
#     the mean fitness with the cost, by quadrature.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_plasticity.log")
const rng = Xoshiro(20261110)
runif(a, b) = runif(rng, a, b)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_plasticity.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the steady lag with a fixed slope, and the reductions at b = 0")
    @variables k β φ τ t c
    g(t) = (1 - φ) * k * t + c                 # a mean genotype moving at (1 − φ) k
    lag(t) = k * t - g(t) - φ * k * (t - τ)    # θ_t − ḡ_t − φ θ_{t−τ}
    check("the lag θ − z̄ of a genotype moving at (1 − φ)k is constant in t", iszero_sym(lag(t + 1) - lag(t)))
    csol = φ * k * τ - (1 - φ) * k / β         # the offset at which the response keeps pace
    check("ḡ′ − ḡ = β·lag at lag = (1 − φ)k/β, for every delay τ of the cue",
          iszero_sym(Symbolics.substitute(g(t + 1) - g(t) - β * lag(t), Dict(c => csol))) &&
          iszero_sym(Symbolics.substitute(lag(t), Dict(c => csol)) - (1 - φ) * k / β))
    for trial in 1:3
        k0, VG, Vs, b, B, τ0 = runif(0.001, 0.1), runif(0.05, 1.0), runif(0.5, 20.0), runif(-0.5, 1.5), runif(0.5, 3.0), rand(rng, 0:6)
        β0, φ0 = VG / (VG + Vs), b / B
        θ(s) = s < 0 ? 0.0 : k0 * s
        gb = 0.0; λ = 0.0
        for s in 0:200_000
            λ = θ(s) - gb - φ0 * θ(s - τ0)
            gb += β0 * λ
        end
        check(@sprintf("recursion iterated, b/B=%.2f τ=%d: lag → steady_lag(k; VG, Vs, b, B)", φ0, τ0),
              isapprox(λ, steady_lag(k0; VG, Vs, b, B); rtol = 1e-8))
        p = Params(MoveType = :continuous, MoveFitnessSD = k0, PhenotypeVar = 1.0, Heritability = 0.5, InitFitnessSD = sqrt(Vs),
                   Plasticity = b, OptimumSlope = B, CueLag = τ0)
        check("    steady_lag(p; VG) is that with p's Plasticity and OptimumSlope",
              steady_lag(p; VG) == steady_lag(k0; VG, Vs = selection_variance(p), b, B))
    end
    VG, Vs = 0.13, 2.7
    p0 = Params(MoveType = :continuous, MoveFitnessSD = 0.02, PhenotypeVar = 1.0, Heritability = 0.4, InitFitnessSD = 1.5)
    kw = (n = 420, VG = 0.13, ω = 1.5, VE = 0.6, sd50 = 0.9, N = 1000)
    red = Params(p0; MoveType = :red_noise, MoveNoiseSD = 0.3, MoveAutocorr = 0.7)
    per = Params(p0; MoveType = :periodic, MoveAmplitude = 0.4, MovePeriod = 25.0)
    check("b = 0: steady_lag === k (V_G + V_s)/V_G, also with B ≠ 1",
          steady_lag(0.02; VG, Vs, b = 0.0, B = 2.5) === 0.02 * (VG + Vs) / VG && steady_lag(0.02; VG, Vs) === 0.02 * (VG + Vs) / VG)
    check("b = 0: critical_rate === the table's form, whatever B and ω_b",
          critical_rate(; kw..., b = 0.0, B = 2.5, ωb = 0.3) === critical_rate(; kw...))
    check("b = 0: red_noise_lag_var === 2σ²(1 − ρ)/((2 − β)(1 − ρ(1 − β))), whatever τ",
          red_noise_lag_var(0.3, 0.7; VG, Vs, b = 0.0, τ = 5) === red_noise_lag_var(0.3, 0.7; VG, Vs) &&
          red_noise_lag_var(0.3, 0.7; VG, Vs) === (β0 = VG / (VG + Vs); 2 * 0.3^2 * (1 - 0.7) / ((2 - β0) * (1 - 0.7 * (1 - β0)))))
    check("b = 0: expected_lag2 and lag_load ignore OptimumSlope, CueLag and PlasticityCostSD (===), all three moves",
          all(q -> expected_lag2(Params(q; OptimumSlope = 2.0, CueLag = 4, PlasticityCostSD = 0.7); VG) === expected_lag2(q; VG) &&
                   lag_load(Params(q; OptimumSlope = 2.0, CueLag = 4, PlasticityCostSD = 0.7); VG) === lag_load(q; VG), (p0, red, per)))

    # ===================================================================================
    say("\n(b) the sustainable rate and the critical rate")
    @variables VGs s r ℓ               # r² = 2 s ℓ, ℓ = ln(W̄(0)/W_min)
    kmax = VGs / s * r / (1 - φ)
    λc = (1 - φ) * kmax * s / VGs       # the steady lag at k_max, with V_G + V_s = s
    check("the steady lag at k_max = V_G/s·√(2 s ℓ)/(1 − φ) is √(2 s ℓ): W̄ = W_min there", iszero_sym(λc - r))
    for trial in 1:4
        b, B, ωb = runif(0.0, 0.9), runif(0.8, 2.0), trial <= 2 ? Inf : runif(0.8, 3.0)
        bh = trial % 2 == 0
        p = Params(PopSize = 1000, Ploidy = 2, PhenotypeVar = runif(0.5, 2.0), Heritability = runif(0.2, 0.6), InitFitnessSD = runif(1.0, 4.0),
                   Plasticity = b, OptimumSlope = B, PlasticityCostSD = ωb, Demography = bh ? :beverton_holt : :fixed, GrowthRate = runif(1.5, 5.0))
        V = runif(0.05, 0.5)
        km = sustainable_rate(p; VG = V)
        sv = initFitSD(p)^2 + V + envVar(p)
        λ = steady_lag(km; VG = V, Vs = selection_variance(p), b, B)
        Wb = QTLSim.TRAIT_FITNESS_PEAK * plasticity_cost(p) * initFitSD(p) / sqrt(sv) * exp(-λ^2 / (2sv))
        Wmin = bh ? 1 / p.GrowthRate : 2 / 1000
        check(@sprintf("sustainable_rate (%s, b/B=%.2f, ω_b=%.2g): W̄ at its steady lag is W_min", bh ? "Beverton–Holt" : "fixed size", b / B, ωb),
              isapprox(Wb, Wmin; rtol = 1e-10))
        ωb == Inf && check("    without a cost it is 1/|1 − b/B| times the rate without plasticity",
                           isapprox(km, sustainable_rate(p; VG = V, b = 0.0) / abs(1 - b / B); rtol = 1e-12))
        ωb == Inf && check("    and so is critical_rate",
                           isapprox(critical_rate(; kw..., b, B), critical_rate(; kw...) / abs(1 - b / B); rtol = 1e-12))
    end
    # critical_rate with a cost: ln B₀ less b²/(2ω_b²)
    No = 1000 * (2 * 420 - 3) / (2 * 420 * 418); Vk = 2No * (1 - 2No / 1000); Ne = 4000 / (Vk + 2)
    Vss = 1.5^2 + 0.6; Vλ = Vss + 0.13 + Vss / (2Ne); B0 = No * 1.5 / sqrt(Vλ)
    check("critical_rate with a cost: V_G/(V_G + V_s)·√(2V_λ(ln B₀ − b²/(2ω_b²)))/|1 − b/B|/sd50",
          isapprox(critical_rate(; kw..., b = 0.5, B = 2.0, ωb = 0.9),
                   0.13 / (0.13 + Vss) * sqrt(2Vλ * (log(B0) - 0.25 / (2 * 0.81))) / 0.75 / 0.9; rtol = 1e-12))

    # ===================================================================================
    say("\n(c) fluctuating optima: red noise and a cycle, with a delayed cue")
    @variables σ ρ a Γ0
    βa = 1 - a
    c0 = σ^2 * (1 - ρ) / (1 - a * ρ)
    check("c₀ = Cov(θ_t, L_t) = σ²(1 − βρ/(1 − aρ)) = σ²(1 − ρ)/(1 − aρ)", iszero_sym(σ^2 * (1 - βa * ρ / (1 - a * ρ)) - c0))
    G0 = 2σ^2 * (1 - ρ) / ((1 + a) * (1 - a * ρ))
    check("Γ(0): (1 − a²)Γ(0) = 2σ²(1 − ρ) − 2a(1 − ρ)c₀, the variance of L′ = aL + θ′ − θ",
          iszero_sym((1 - a^2) * G0 - 2σ^2 * (1 - ρ) + 2a * (1 - ρ) * c0))
    Γ(h) = a^h * G0 - (1 - ρ) * c0 * sum(a^(h - 1 - j) * ρ^j for j in 0:h-1; init = 0)
    check("Γ(h + 1) = aΓ(h) − (1 − ρ)ρ^h c₀ for h = 0…5, with Σ_j a^{h−1−j}ρ^j = (a^h − ρ^h)/(a − ρ)",
          all(h -> iszero_sym(Γ(h + 1) - a * Γ(h) + (1 - ρ) * ρ^h * c0) &&
                   iszero_sym((a - ρ) * sum(a^(h - j) * ρ^j for j in 0:h) - (a^(h + 1) - ρ^(h + 1))), 0:5))
    check("τ = 0: (1 + φ²)Γ(0) − 2φΓ(0) = (1 − φ)²Γ(0)", iszero_sym((1 + φ^2) * G0 - 2φ * Γ(0) - (1 - φ)^2 * G0))
    @variables cs
    check("a cycle: |1 − φ e^{−iντ}|² = (1 − φ cos x)² + φ² sin² x = 1 + φ² − 2φ cos x (sin² = 1 − cos²)",
          iszero_sym((1 - φ * cs)^2 + φ^2 * (1 - cs^2) - (1 + φ^2 - 2φ * cs)))
    for trial in 1:3
        σ0, ρ0, V, Vs0, b, B, τ0 = runif(0.1, 1.0), runif(0.0, 0.95), runif(0.05, 1.0), runif(0.5, 10.0), runif(0.1, 1.2), runif(0.8, 2.0), rand(rng, 1:6)
        β0, φ0 = V / (V + Vs0), b / B
        # the lag as a moving average of the optimum, L_t = θ_t − β Σ_{j≥1} (1 − β)^{j−1} θ_{t−j}
        cf = [1.0; [-β0 * (1 - β0)^(j - 1) for j in 1:4000]]
        Γn(h) = σ0^2 * sum(cf[i] * sum(cf[j] * ρ0^abs(i - j + h) for j in eachindex(cf)) for i in eachindex(cf))
        check(@sprintf("red_noise_lag_autocov(σ, ρ=%.2f, h=%d) against the moving-average sum", ρ0, τ0),
              isapprox(red_noise_lag_autocov(σ0, ρ0, τ0; VG = V, Vs = Vs0), Γn(τ0); rtol = 1e-7))
        check(@sprintf("red_noise_lag_var with b/B=%.2f, τ=%d: (1 + φ²)Γ(0) − 2φΓ(τ)", φ0, τ0),
              isapprox(red_noise_lag_var(σ0, ρ0; VG = V, Vs = Vs0, b, B, τ = τ0), (1 + φ0^2) * Γn(0) - 2φ0 * Γn(τ0); rtol = 1e-7))
        # a cycle, iterated: the lag's amplitude
        P0, A0 = runif(15.0, 200.0), runif(0.1, 1.0)
        θ(s) = A0 * sin(2π * s / P0)
        gb = 0.0; lags = Float64[]
        for s in 0:60_000
            λ = θ(s) - gb - φ0 * θ(s - τ0)
            s > 40_000 && push!(lags, λ)
            gb += β0 * λ
        end
        q = Params(MoveType = :periodic, MoveAmplitude = A0, MovePeriod = P0, MoveFitnessSD = 0.0, PhenotypeVar = 1.0, Heritability = 0.5,
                   InitFitnessSD = sqrt(Vs0 - 0.5 + 1e-12 > 0 ? Vs0 - 0.5 : 0.1), Plasticity = b, OptimumSlope = B, CueLag = τ0)
        check(@sprintf("a cycle of period %.0f iterated: mean squared lag = expected_lag2(p; VG)", P0),
              isapprox(sum(abs2, lags) / length(lags), expected_lag2(q; VG = V); rtol = 2e-3))
    end

    # ===================================================================================
    say("\n(d) a heritable reaction norm")
    @variables Gaa Gab Gbb ε ω2 VE cb λs abar bbar θs kk Bs
    G = [Gaa Gab; Gab Gbb]
    v = [1, ε]
    sv = ω2 + VE + (v' * G * v)
    Q0 = (v * v') ./ (ω2 + VE)
    Gs0 = G - (G * v) * (G * v)' ./ sv
    Id = [1 0; 0 1]
    check("no cost: (I + G Q)(G − G v vᵀ G/s) = G, Q = v vᵀ/(ω² + V_E), s = ω² + V_E + vᵀGv", iszero_sym((Id + G * Q0) * Gs0 - G))
    check("no cost: Δ = G_s v λ/(ω² + V_E) = G v λ/s", iszero_sym(Gs0 * v .* (λs / (ω2 + VE)) - G * v .* (λs / sv)))
    lagq = (1 - bbar / Bs) * kk * sv / (v' * G * v)
    Δq = G * v .* (lagq / sv)
    check("the quasi-steady lag: vᵀΔ + b̄ k/B = k at λ = (1 − b̄/B) k s/(vᵀGv)", iszero_sym(v' * Δq + bbar * kk / Bs - kk))
    @variables x
    far = Symbolics.substitute(Δq[2] / x, Dict(ε => 1 / x))      # ε·Δb̄ at ε = 1/x
    far = Symbolics.simplify_fractions(far)
    check("far from the reference environment ε·Δb̄ → (1 − b̄/B) k (so d b̄/dt ≈ (B − b̄)/t)",
          iszero_sym(Symbolics.substitute(far, Dict(x => 0)) - (1 - bbar / Bs) * kk))
    # with a cost, exactly, in rational arithmetic
    okK, okfix, okline = true, true, true
    for trial in 1:20
        rr() = big(rand(rng, 1:40)) // big(rand(rng, 1:40))
        A = [rr() rr(); rr() rr()]
        Gr = A * A' + (1 // 10) * Id
        e, w2, ve, c1 = rr() - 1 // 2, rr(), rr(), rr()
        vr = [1 // 1, e]
        Qr = (vr * vr') .// (w2 + ve) + [0 0; 0 c1]
        Kr = Qr * inv(Id + Gr * Qr)
        Gsr = inv(Id + Gr * Qr) * Gr
        okK &= Gsr == Gr - Gr * Kr * Gr && Kr == Kr' && Gsr == Gsr'
        # the response at ā = θ − b̄ε … with a cost the only fixed point has λ = 0 and b̄ = 0
        th = rr()
        Δ(ab, bb) = Gsr * (vr .* ((th - ab - bb * e) // (w2 + ve)) - [0, c1 * bb])
        okfix &= Δ(th, 0 // 1) == [0, 0] && Δ(th - e // 2, 1 // 2) != [0, 0]
        # without a cost every (ā, b̄) with ā + b̄ε = θ is a fixed point
        Q1 = (vr * vr') .// (w2 + ve)
        G1 = inv(Id + Gr * Q1) * Gr
        bb = rr()
        okline &= G1 * (vr .* ((th - (th - bb * e) - bb * e) // (w2 + ve))) == [0, 0]
    end
    check("with a cost (20 rational cases): G_s = (I + GQ)⁻¹G = G − G K G, K = Q(I + GQ)⁻¹, both symmetric", okK)
    check("with a cost, a constant environment: Δ = 0 at λ = 0, b̄ = 0 (the slope is lost), not at λ = 0, b̄ ≠ 0", okfix)
    check("without a cost: Δ = 0 on the whole line ā + b̄ε = θ", okline)
    for trial in 1:3
        A = randn(rng, 2, 2); Gn = A * A' ./ 4 + 0.02I
        e, ω0, ve, θ0, ab, bb, ωb = runif(-3.0, 3.0), runif(0.7, 3.0), runif(0.2, 1.0), runif(-1.0, 1.0), runif(-1.0, 1.0), runif(-0.5, 1.5), runif(0.5, 3.0)
        vn = [1.0, e]
        r0 = reaction_norm_response(Gn, e; ω = ω0, VE = ve, θ = θ0, abar = ab, bbar = bb)
        sn = ω0^2 + ve + vn' * Gn * vn
        check("reaction_norm_response without a cost: Δ = G v λ/s, G_s = G − G v vᵀ G/s, W̄ = ω/√s·exp(−λ²/2s)",
              isapprox(r0.Δ, Gn * vn .* ((θ0 - ab - bb * e) / sn); rtol = 1e-10) && isapprox(r0.Gs, Gn - (Gn * vn) * (Gn * vn)' ./ sn; rtol = 1e-10) &&
              isapprox(r0.Wbar, ω0 / sqrt(sn) * exp(-(θ0 - ab - bb * e)^2 / (2sn)); rtol = 1e-12))
        rc = reaction_norm_response(Gn, e; ω = ω0, VE = ve, θ = θ0, abar = ab, bbar = bb, ωb)
        Qn = (vn * vn') ./ (ω0^2 + ve) + [0 0; 0 1 / ωb^2]
        check("    with a cost: G_s = (I + GQ)⁻¹G and Δ = G_s (v λ/(ω² + V_E) − e₂ b̄/ω_b²)",
              isapprox(rc.Gs, (I + Gn * Qn) \ Gn; rtol = 1e-10) &&
              isapprox(rc.Δ, rc.Gs * (vn .* ((θ0 - ab - bb * e) / (ω0^2 + ve)) .- [0.0, bb / ωb^2]); rtol = 1e-10))
        # the mean fitness with the cost by quadrature over (z, b): z | b is Gaussian
        Vz, Vb, Czb = vn' * Gn * vn + ve, Gn[2, 2], (Gn * vn)[2]
        zb = ab + bb * e
        bs = range(bb - 9sqrt(Vb), bb + 9sqrt(Vb); length = 4001)
        f(b) = (μ = zb + Czb / Vb * (b - bb); V = Vz - Czb^2 / Vb;
                exp(-(b - bb)^2 / (2Vb)) / sqrt(2π * Vb) * exp(-b^2 / (2ωb^2)) * sqrt(ω0^2 / (ω0^2 + V)) * exp(-(μ - θ0)^2 / (2 * (ω0^2 + V))))
        quad = sum(f, bs) * step(bs)
        check("    its mean fitness E[exp(−(z−θ)²/2ω²)·exp(−b²/2ω_b²)], by quadrature", isapprox(rc.Wbar, quad; rtol = 1e-8))
        k0, B0 = runif(0.005, 0.05), runif(0.8, 2.0)
        λq = reaction_norm_lag(Gn, e; k = k0, B = B0, bbar = bb, ω = ω0, VE = ve)
        rq = reaction_norm_response(Gn, e; ω = ω0, VE = ve, θ = ab + bb * e + λq, abar = ab, bbar = bb)
        check("    reaction_norm_lag: the mean phenotype then moves k per generation (vᵀΔ + b̄ k/B = k)",
              isapprox(vn' * rq.Δ + bb * k0 / B0, k0; rtol = 1e-10))
        check("    reaction_norm_lag without genetic variance in the slope is steady_lag",
              isapprox(reaction_norm_lag([Gn[1, 1] 0.0; 0.0 0.0], e; k = k0, B = B0, bbar = bb, ω = ω0, VE = ve),
                       steady_lag(k0; VG = Gn[1, 1], Vs = ω0^2 + ve, b = bb, B = B0); rtol = 1e-12))
    end

    finish(LOG, t0)
end

runmain(main, @__FILE__)
