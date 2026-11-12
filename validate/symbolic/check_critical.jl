# Symbolic checks of the critical rate under drift (src/critical.jl, drift_critical_rate),
# and numerical spot checks against the code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_critical.jl
# Writes validate/symbolic/check_critical.log; exits 1 if any check fails.
#
# (a) the envelope theorem: if ln k(λ, c) = ln λ + H(u), u = λ²/(2s) + c, then at a stationary
#     point in λ, ∂ ln k/∂c = −s/λ²; so k_c(c)/k_c(0) = 1 − c/M + O(c²), M = λ*²/s; for the
#     power law H = −εu (V_G ∝ n^ε, V_G ≪ V_s) exactly: λ* = √(s/ε), k_c(c) = k_c(0)e^{−εc};
# (b) the elasticities: of the house-of-Gauss V_G in the number of genes ng = Ploidy·N_e,
#     1/(ng(2c/V_s + b)); of effective_size(N, n) in n, 2N(n − 2)/(n²(V_k + 2)); of the
#     Beverton–Holt breeders K(RW − 1)/(R − 1) in W, RW/(RW − 1); so M = 1/(κν) with
#     κ = V_s/(V_G + V_s)·ε_V·ε_N and ν = 1 (fixed size) or RW/(RW − 1);
# (c) numerically: drift_critical_rate's k against a brute-force maximum, its M against the
#     finite-difference slope in the cost, the envelope (total = partial derivative in c) and
#     M = 1/(κν) at the critical lag.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_critical.log")
const rng = Xoshiro(20261111)
runif(a, b) = runif(rng, a, b)

function randparams()
    dem = rand(rng, (:fixed, :beverton_holt))
    Params(PopSize = rand(rng, (500, 1000, 2000)), Ploidy = 2, NoQTL = rand(rng, (50, 100, 500, 1000)),
           PhenotypeVar = runif(0.02, 0.1), Heritability = runif(0.2, 0.6), MutationScalar = runif(0.002, 0.012),
           MutationRatePerGamete = 1e-4, InitFitnessSD = runif(1.5, 5.0), MoveType = :continuous,
           Demography = dem, GrowthRate = dem === :fixed ? 1.5 : runif(1.5, 5.0))
end

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_critical.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the envelope theorem and the linear-margin form")
    @variables λ s c h1 ε A
    u = λ^2 / (2s) + c
    dλ = 1 / λ + h1 * λ / s            # ∂ ln k/∂λ with h1 = H′(u)
    h1star = -s / λ^2                  # the root of dλ = 0 in h1
    check("∂ ln k/∂λ = 0 ⇔ H′(u) = −s/λ²", iszero_sym(Symbolics.substitute(dλ, Dict(h1 => h1star))))
    check("then ∂ ln k/∂c = H′(u) = −s/λ²: M = λ*²/s", iszero_sym(h1star + 1 / (λ^2 / s)))
    lnk = log(A) + log(λ) - ε * u      # the power law H = −εu
    d = Symbolics.expand_derivatives(Differential(λ)(lnk))
    check("power law: ∂ ln k/∂λ = 0 at λ* = √(s/ε)", iszero_sym(Symbolics.substitute(d, Dict(λ => sqrt(s / ε)))))
    kmax(cc) = Symbolics.substitute(lnk, Dict(λ => sqrt(s / ε), c => cc))
    check("power law: ln k_c(c) − ln k_c(0) = −εc exactly (M = 1/ε = λ*²/s)", iszero_sym(kmax(c) - kmax(0) + ε * c))
    d2 = Symbolics.expand_derivatives(Differential(λ)(d))
    check("power law: the stationary point is a maximum (∂² ln k/∂λ² = −2ε/s at λ*)",
          iszero_sym(Symbolics.substitute(d2, Dict(λ => sqrt(s / ε))) + 2ε / s))

    # ===================================================================================
    say("\n(b) the elasticities")
    @variables cg Vs b μ α2 ng N n K R W
    F = cg^2 / Vs + cg * b - μ * α2      # the house-of-Gauss equation per gene
    dcdb = -Symbolics.expand_derivatives(Differential(b)(F)) / Symbolics.expand_derivatives(Differential(cg)(F))
    check("implicit derivative: dc/db = −c/(2c/V_s + b)", iszero_sym(dcdb + cg / (2cg / Vs + b)))
    # b = α²/(2V_s) + 1/ng, so db/d ln ng = −1/ng and d ln c/d ln ng = (−dc/db)/(c ng)
    check("ε_V = d ln V_G/d ln ng = 1/(ng (2c/V_s + b))", iszero_sym(-dcdb / (cg * ng) - 1 / (ng * (2cg / Vs + b))))
    Vk = 2N * (n - 1) / n^2
    Ne = 4N / (Vk + 2)
    εN = Symbolics.expand_derivatives(Differential(n)(log(Ne))) * n
    check("ε_N = d ln N_e/d ln n = 2N(n − 2)/(n²(V_k + 2))", iszero_sym(εN - 2N * (n - 2) / (n^2 * (Vk + 2))))
    nbh = K * (R * W - 1) / (R - 1)
    ν = Symbolics.expand_derivatives(Differential(W)(log(nbh))) * W
    check("Beverton–Holt: d ln n/d ln W̄ = RW̄/(RW̄ − 1)", iszero_sym(ν - R * W / (R * W - 1)))

    # ===================================================================================
    say("\n(c) numerical checks against drift_critical_rate")
    for trial in 1:6
        p = randparams()
        r = drift_critical_rate(p)
        r.k > 0 || continue
        Vs = selection_variance(p)
        rate(λ, cc) = (st = QTLSim._feedback_state(p, λ, cc); st.VG * λ / (st.VG + Vs))
        brute = maximum(rate(l, 0.0) for l in range(0, 4r.λ, length = 20_001))
        tag = @sprintf("%s L=%d N=%d", p.Demography, p.NoQTL, p.PopSize)
        check("$tag: k_c is the maximum over λ (brute force, 2e-4)", isapprox(r.k, brute; rtol = 2e-4) && r.k >= brute * (1 - 1e-9))
        h = 1e-3
        total = (log(drift_critical_rate(p; cost = h).k) - log(r.k)) / h
        partial = (log(rate(r.λ, h)) - log(rate(r.λ, 0.0))) / h
        check(@sprintf("%s: envelope, d ln k_c/dc %.4f = ∂ ln k/∂c at λ* %.4f (1%%)", tag, total, partial),
              isapprox(total, partial; rtol = 0.01))
        check(@sprintf("%s: d ln k_c/dc %.4f ≈ −1/M %.4f (s varies with V_G: 5%%)", tag, total, -1 / r.M),
              isapprox(total, -1 / r.M; rtol = 0.05))
        # M = 1/(κν) from the elasticities at the critical point
        α2 = QTLSim.mutation_effect_variance(p)
        ngv = p.Ploidy * effective_size(r.N, r.n)
        bb = α2 / (2Vs) + 1 / ngv
        cgene = r.VG / (2p.NoQTL)
        εV = 1 / (ngv * (2cgene / Vs + bb))
        Vk0 = 2r.N * (r.n - 1) / r.n^2
        εN0 = 2r.N * (r.n - 2) / (r.n^2 * (Vk0 + 2))
        R = p.GrowthRate
        ν0 = p.Demography === :fixed ? 1.0 : R * r.Wbar / (R * r.Wbar - 1)
        # under Beverton–Holt N = n/W̄ also moves with W̄: d ln N_e/d ln W̄ = ε_N ν + (∂ ln N_e/∂ ln N)(ν − 1)
        εNN = p.Demography === :fixed ? 0.0 : 2 / (Vk0 + 2)          # ∂ ln N_e/∂ ln N
        κν = Vs / (r.VG + Vs) * εV * (εN0 * ν0 + εNN * (ν0 - 1))
        check(@sprintf("%s: M = λ*²/s %.4f ≈ 1/(κν) %.4f (5%%)", tag, r.M, 1 / κν), isapprox(r.M, 1 / κν; rtol = 0.05))
    end
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0, PhenotypeVar = 0.03, Heritability = 0.4, MutationScalar = 0.005,
               MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous)
    check("_house_of_gauss(p, 2N) = house_of_gauss_vg(p) (nsurv = N)",
          isapprox(QTLSim._house_of_gauss(p, 2effective_size(1000, 1000)), house_of_gauss_vg(p); rtol = 1e-14))
    check("a cost lowers the critical rate; plasticity without a cost divides it by 1 − b/B",
          drift_critical_rate(p; cost = 0.3).k < drift_critical_rate(p).k &&
          isapprox(drift_critical_rate(Params(p; Plasticity = 0.5)).k, 2drift_critical_rate(p).k; rtol = 1e-12))

    finish(LOG, t0)
end

runmain(main, @__FILE__)
