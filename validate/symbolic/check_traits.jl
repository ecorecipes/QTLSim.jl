# Symbolic checks of the theory for several traits (src/traits_theory.jl,
# src/rare_alleles_traits.jl), and numerical spot checks against its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_traits.jl
# Writes validate/symbolic/check_traits.log; exits 1 if any check fails.
#
# (a) gaussian_G: G = A^½ X A^½ turns G(A + G)⁻¹G = M into X(I + X)⁻¹X = M̃ (symbolic,
#     general symmetric 2 × 2), X = (M̃ + (M̃² + 4M̃)^½)/2 solves it on each eigenvalue, and
#     so does (m + √(m² + 4ma))/2 on each shared axis; numerically for general M, A
#     (n = 2, 3), against gaussian_G_iterate, and gaussian_G_loci = gaussian_G(2L·M, A);
# (b) steady_lag_vector: λ = (A + G)G⁻¹k solves G(A + G)⁻¹λ = k (symbolic 2 × 2; numerical);
# (c) red_noise_lag_cov: the lag–optimum recursion, the n = 1 Lyapunov solution, which is
#     red_noise_lag_var, and numerically X = F X Fᵀ + Q for n = 2;
# (d) periodic_lag_cov and periodic_response: the transfer function's modulus D, the lag
#     amplitude 2 sin(ν/2)/D, the n = 1 reduction, and simulations of the recursion;
# (e) the deflection inequality (the response is never further from the g_max axis than
#     the gradient), symbolically for 2 × 2 and on random examples;
# (f) the rare-alleles Riccati equation c A⁻¹ c + J c + c Jᵀ = μ_u C at house_of_gauss_G,
#     its stabilising root, its one-trait reduction to house_of_gauss_vg and its limits.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_traits.log")
const rng = Xoshiro(20261110)
runif(a, b) = runif(rng, a, b)

# the Lyapunov solution by its series Σ F^k Q F^kᵀ (F stable)
function lyap_series(F, Q; tol = 1e-15)
    X, T = copy(Q), copy(Q)
    for _ in 1:10^6
        T = F * T * F'
        X += T
        maximum(abs, T) <= tol * maximum(abs, X) && return X
    end
    error("lyap_series did not converge")
end

twotrait(; L = 10, N = 500, μ = 1e-3, kw...) =
    Params(; Ploidy = 2, NoQTL = L, PopSize = N, MutationRatePerGamete = μ, RecombPoissonMean = 2.0, NoTraits = 2,
           MutCorr = runif(-0.6, 0.6), SelCorr = runif(-0.6, 0.6), PhenotypeVar = runif(0.5, 2.0),
           Heritability = runif(0.3, 0.7), MutationScalar = 10^runif(-3, -2), kw...)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_traits.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) gaussian_G, the multivariate Gaussian balance G(A + G)⁻¹G = M")
    @variables h11 h12 h22 x11 x12 x22 μ r m a
    H = [h11 h12; h12 h22]
    X = [x11 x12; x12 x22]
    Gs = H * X * H
    As = H * H
    Id2 = [1 0; 0 1]
    # A + G = H(I + X)H, so (A + G)⁻¹ = H⁻¹(I + X)⁻¹H⁻¹ and G(A + G)⁻¹G = H X(I + X)⁻¹X H;
    # X commutes with I + X, so X(I + X)⁻¹X = X²(I + X)⁻¹ is a function of X
    check("G = A^½XA^½: A + G = A^½(I + X)A^½ (general symmetric 2 × 2), so G(A + G)⁻¹G = A^½X(I + X)⁻¹XA^½",
          iszero_sym(As + Gs - H * (Id2 + X) * H))
    check("X commutes with I + X, so X(I + X)⁻¹X = X²(I + X)⁻¹ = M̃ is solved eigenvalue by eigenvalue",
          iszero_sym(X * (Id2 + X) - (Id2 + X) * X))
    x = (μ + r) / 2                                  # r = √(μ² + 4μ)
    check("on an eigenvalue μ of M̃: x = (μ + √(μ² + 4μ))/2 gives x²/(1 + x) = μ",
          iszero_mod_sqrt(x^2 - μ * (1 + x), r, μ^2 + 4μ))
    g = (m + r) / 2                                  # r = √(m² + 4ma)
    check("shared axes: g = (m + √(m² + 4ma))/2 gives g²/(a + g) = m",
          iszero_mod_sqrt(g^2 - m * (a + g), r, m^2 + 4m * a))

    for n in (2, 3), trial in 1:2
        M = randspd(rng, n; scale = 10^runif(-3, -1)); A = randspd(rng, n; scale = runif(1.0, 20.0))
        G = gaussian_G(M, A)
        check(@sprintf("gaussian_G, n = %d: G(A + G)⁻¹G = M, G ⪰ 0, = gaussian_G_iterate", n),
              isapprox(G * ((A + G) \ G), M; rtol = 1e-9) && minimum(eigvals(Symmetric(G))) >= -1e-14 &&
              isapprox(G, QTLSim.gaussian_G_iterate(M, A); rtol = 1e-8))
    end
    U = Matrix(qr(randn(rng, 3, 3)).Q)
    mv, av = 10 .^ (-3 .+ 2rand(rng, 3)), 1 .+ 10rand(rng, 3)
    Gc = gaussian_G(U * Diagonal(mv) * U', U * Diagonal(av) * U')
    check("gaussian_G, shared eigenvectors: eigenvalues (m + √(m² + 4ma))/2 on each axis",
          isapprox(Gc, U * Diagonal((mv .+ sqrt.(mv .^ 2 .+ 4mv .* av)) ./ 2) * U'; rtol = 1e-10))
    pr = twotrait()
    check("gaussian_G_loci = gaussian_G(2L·M, A) (diploid, recombining)",
          isapprox(gaussian_G_loci(pr), gaussian_G(2pr.NoQTL .* mutational_covariance(pr), selection_variance_matrix(pr)); rtol = 1e-12))
    pr1 = twotrait(L = 1)
    check("gaussian_G_loci = gaussian_G(2·M, A) without recombination (one locus)",
          isapprox(gaussian_G_loci(pr1), gaussian_G(2 .* mutational_covariance(pr1), selection_variance_matrix(pr1)); rtol = 1e-12))

    # ===================================================================================
    say("\n(b) steady_lag_vector")
    @variables g11 g12 g22 a11 a12 a22 k1 k2
    Gm, Am, kv = [g11 g12; g12 g22], [a11 a12; a12 a22], [k1, k2]
    @variables u1 u2
    uv = [u1, u2]                                    # u = G⁻¹k, so k = Gu and λ = (A + G)u
    check("λ = (A + G)G⁻¹k: with u = G⁻¹k, (A + G)⁻¹λ = u and G(A + G)⁻¹λ = Gu = k (symbolic 2 × 2)",
          iszero_sym((Am + Gm) * uv - (Am * uv + Gm * uv)))
    for n in (2, 3)
        G = randspd(rng, n; scale = 0.1); A = randspd(rng, n; scale = 10.0); k = randn(rng, n) .* 0.01
        λ = steady_lag_vector(k; G, A)
        check("steady_lag_vector, n = $n: trait_response(λ) = k", isapprox(trait_response(λ; G, A), k; rtol = 1e-10))
    end
    check("steady_lag_vector, n = 1: steady_lag k(V_G + V_s)/V_G",
          isapprox(steady_lag_vector([0.02]; G = fill(0.3, 1, 1), A = fill(9.0, 1, 1))[1], steady_lag(0.02; VG = 0.3, Vs = 9.0); rtol = 1e-12))

    # ===================================================================================
    say("\n(c) red_noise_lag_cov, red_noise_lag_var")
    @variables β ρ σ2 λt θt η
    # one step: z̄' = z̄ + βλ, θ' = ρθ + η, so λ' = θ' − z̄' = (1 − β)λ + (ρ − 1)θ + η
    zt = θt - λt
    check("λ_{t+1} = (1 − β)λ_t + (ρ − 1)θ_t + √(1 − ρ²)η_{t+1}",
          iszero_sym((ρ * θt + η) - (zt + β * λt) - ((1 - β) * λt + (ρ - 1) * θt + η)))
    F = [1-β ρ-1; 0 ρ]
    Q = (1 - ρ^2) * σ2 .* [1 1; 1 1]
    X12 = σ2 * (1 - ρ) / (1 - ρ * (1 - β))
    X11 = 2σ2 * (1 - ρ) / ((2 - β) * (1 - ρ * (1 - β)))
    Xs = [X11 X12; X12 σ2]
    check("n = 1: X = [2σ²(1−ρ)/((2−β)(1−ρ(1−β))), σ²(1−ρ)/(1−ρ(1−β)); ·, σ²] solves X = FXFᵀ + Q",
          iszero_sym(F * Xs * F' + Q - Xs))
    check("white noise (ρ = 0): 2σ²/(2 − β)", iszero_sym(Symbolics.substitute(X11, Dict(ρ => 0)) - 2σ2 / (2 - β)))
    for trial in 1:3
        σ0, ρ0, VG0, Vs0 = runif(0.1, 2.0), runif(-0.5, 0.95), runif(0.05, 1.0), runif(1.0, 30.0)
        check(@sprintf("red_noise_lag_cov (n = 1) = red_noise_lag_var = symbolic X₁₁, ρ = %.2f", ρ0),
              isapprox(red_noise_lag_cov(fill(σ0^2, 1, 1), ρ0; G = fill(VG0, 1, 1), A = fill(Vs0, 1, 1))[1],
                       red_noise_lag_var(σ0, ρ0; VG = VG0, Vs = Vs0); rtol = 1e-10) &&
              isapprox(red_noise_lag_var(σ0, ρ0; VG = VG0, Vs = Vs0),
                       num(X11, Dict(σ2 => σ0^2, ρ => ρ0, β => VG0 / (VG0 + Vs0))); rtol = 1e-12))
    end
    for trial in 1:2
        n = 2
        G = randspd(rng, n; scale = 0.3); A = randspd(rng, n; scale = 5.0); Σ = randspd(rng, n; scale = 0.5); ρ0 = runif(0.0, 0.9)
        B = G / (A + G)
        Fn = [I-B (ρ0 - 1) .* I(n); zeros(n, n) ρ0 .* I(n)]
        Qn = (1 - ρ0^2) .* [Σ Σ; Σ Σ]
        Xn = lyap_series(Fn, Qn)
        check(@sprintf("red_noise_lag_cov, n = 2, ρ = %.2f: the lag block of X = FXFᵀ + Q (series)", ρ0),
              isapprox(red_noise_lag_cov(Σ, ρ0; G, A), Xn[1:n, 1:n]; rtol = 1e-9))
    end
    gd, ad, sd = [0.2, 0.05], [4.0, 9.0], [0.3, 1.1]
    check("red_noise_lag_cov, shared axes: each axis the one-trait value with its own β",
          isapprox(red_noise_lag_cov(Diagonal(sd .^ 2) |> Matrix, 0.6; G = Diagonal(gd) |> Matrix, A = Diagonal(ad) |> Matrix),
                   Diagonal([red_noise_lag_var(sd[i], 0.6; VG = gd[i], Vs = ad[i]) for i in 1:2]); rtol = 1e-10, atol = 1e-14))

    # ===================================================================================
    say("\n(d) periodic_response, periodic_lag_cov")
    @variables C Sn sh ch
    D2 = 1 - 2(1 - β) * C + (1 - β)^2
    check("|e^{iν} − (1 − β)|² = 1 − 2(1 − β)cos ν + (1 − β)² = D²",
          iszero_mod_sqrt((C - (1 - β))^2 + Sn^2 - D2, Sn, 1 - C^2))
    # half angles: cos ν = ch² − sh², sin ν = 2 sh ch, ch² = 1 − sh²
    check("|e^{iν} − 1|² = (2 sin(ν/2))², so lag amplitude 2 sin(ν/2)/D",
          iszero_mod_sqrt((ch^2 - sh^2 - 1)^2 + (2sh * ch)^2 - 4sh^2, ch, 1 - sh^2))
    # the transfer function: Z e^{iν} = (1 − β)Z + βΘ ⇒ Z/Θ = β/(e^{iν} − (1 − β)), λ/Θ = 1 − Z/Θ
    check("lag/optimum = (e^{iν} − 1)/(e^{iν} − (1 − β)) = 1 − β/(e^{iν} − (1 − β))",
          iszero_sym((C - 1) / (C - (1 - β)) - (1 - β / (C - (1 - β)))))   # with e^{iν} as one symbol
    for P0 in (7.0, 30.0, 200.0)
        VG0, Vs0 = runif(0.1, 1.0), runif(2.0, 20.0)
        r = periodic_response(P0; VG = VG0, Vs = Vs0)
        β0, ν = VG0 / (VG0 + Vs0), 2π / P0
        z = 0.0; T = round(Int, 40 / β0 + 20P0)
        zs = Float64[]
        for t in 0:T
            θ = sin(ν * t)
            t > T - 5P0 && push!(zs, z - r.gain * sin(ν * t - r.phase))
            z += β0 * (θ - z)
        end
        check(@sprintf("periodic_response, P = %.0f: z̄_t = gain·sin(νt − phase) in the recursion", P0), maximum(abs, zs) < 1e-9)
        lagd = [sin(ν * t) - (r.gain * sin(ν * t - r.phase)) for t in 0:round(Int, 4P0)]
        check(@sprintf("… and the lag's amplitude is lag_amplitude (%.3f)", r.lag_amplitude),
              isapprox(sqrt(2 * sum(abs2, lagd[1:end-1]) / (length(lagd) - 1)), r.lag_amplitude; rtol = 0.02) &&
              isapprox(abs(1 - β0 / (exp(im * ν) - (1 - β0))), r.lag_amplitude; rtol = 1e-12))
        check("periodic_lag_cov, n = 1 = lag_amplitude²/2",
              isapprox(periodic_lag_cov(P0, [1.0]; G = fill(VG0, 1, 1), A = fill(Vs0, 1, 1))[1], r.lag_amplitude^2 / 2; rtol = 1e-12))
    end
    r = periodic_response(1e6; VG = 0.2, Vs = 9.8)
    check("slow cycle: gain → 1, delay → 1/β", isapprox(r.gain, 1; rtol = 1e-6) && isapprox(r.delay, 1 / 0.02; rtol = 1e-4))
    # n = 2: the mean outer product over whole cycles of the simulated recursion
    G = randspd(rng, 2; scale = 0.3); A = randspd(rng, 2; scale = 5.0); f = randn(rng, 2); P0 = 25
    B = G / (A + G)
    zv = zeros(2); acc = zeros(2, 2); T = 4000; ncyc = 40
    for t in 0:T+ncyc*P0-1
        θ = sin(2π * t / P0) .* f
        t >= T && (acc .+= (θ - zv) * (θ - zv)')
        zv = zv + B * (θ - zv)
    end
    check("periodic_lag_cov, n = 2: = the cycle average of λλᵀ in the recursion",
          isapprox(periodic_lag_cov(P0, f; G, A), acc ./ (ncyc * P0); rtol = 1e-8))

    # ===================================================================================
    say("\n(e) the deflection of the response toward g_max")
    @variables gg1 gg2 b1 b2
    # cos² of the angle to e₁: response g₁²β₁²/(g₁²β₁² + g₂²β₂²), gradient β₁²/(β₁² + β₂²)
    check("cos²(response, e₁) − cos²(gradient, e₁) has the sign of (g₁² − g₂²)β₁²β₂²",
          iszero_sym((gg1^2 * b1^2 * (b1^2 + b2^2) - b1^2 * (gg1^2 * b1^2 + gg2^2 * b2^2)) - (gg1^2 - gg2^2) * b1^2 * b2^2))
    ok, worst = true, -Inf
    for n in (2, 3, 4), _ in 1:100
        G = randspd(rng, n; scale = runif(0.01, 1.0), ridge = 0.01); A = randspd(rng, n; scale = 10.0); λ = randn(rng, n)
        d = response_deflection(λ; G, A)
        worst = max(worst, d.angle_response_gmax - d.angle_gradient_gmax)
        ok &= d.angle_response_gmax <= d.angle_gradient_gmax + 1e-9 && isapprox(d.response, G * d.β; rtol = 1e-12)
    end
    check(@sprintf("response_deflection, 300 random examples (n = 2–4): angle to g_max never grows (max %.1e°)", worst), ok)

    # ===================================================================================
    say("\n(f) the rare-alleles Riccati equation, house_of_gauss_G")
    @variables cc Aa α2 ng μu
    J1 = (α2 / (2Aa) + 1 / ng) / 2               # one trait: 𝓡 = α²/(2A)
    check("one trait: c²/A + 2Jc = μα² is c²/V_s + c(α²/(2V_s) + 1/ng) = μα²",
          iszero_sym((cc^2 / Aa + 2J1 * cc - μu * α2) - (cc^2 / Aa + cc * (α2 / (2Aa) + 1 / ng) - μu * α2)))
    for trial in 1:3
        pr = twotrait(L = rand(rng, 5:30), N = rand(rng, 100:2000), μ = 10^runif(-4, -2))
        nu, μ0 = QTLSim.rare_units(pr)
        Cm = QTLSim.mutation_effect_covariance(pr)
        A = selection_variance_matrix(pr)
        Ai = inv(A)
        ng0 = pr.Ploidy * effective_size(pr.PopSize, pr.PopSize)
        J = (QTLSim.rare_alleles_rate(Cm, Ai) + (1 / ng0) * I) ./ 2
        G = house_of_gauss_G(pr)
        c = G ./ nu
        res = c * Ai * c + J * c + c * J' - μ0 .* Cm
        check(@sprintf("house_of_gauss_G solves c A⁻¹ c + Jc + cJᵀ = μ_u C (L=%d μ=%.1e)", pr.NoQTL, μ0),
              maximum(abs, res) <= 1e-9 * maximum(abs, μ0 .* Cm))
        check("… the stabilising root: J + c A⁻¹ has eigenvalues of positive real part, c ⪰ 0",
              minimum(real, eigvals(J + c * Ai)) > 0 && minimum(eigvals(Symmetric(c))) >= -1e-15)
    end
    for trial in 1:2
        p1 = Params(Ploidy = 2, NoQTL = rand(rng, 5:40), PopSize = rand(rng, 100:2000), MutationRatePerGamete = 10^runif(-4, -2),
                    RecombPoissonMean = 2.0, MutationScalar = 10^runif(-3, -2), PhenotypeVar = runif(0.5, 2.0))
        check("house_of_gauss_G, one trait = house_of_gauss_vg (Hermisson–Wagner A4)",
              isapprox(house_of_gauss_G(p1)[1, 1], house_of_gauss_vg(p1); rtol = 1e-9))
    end
    # limits: rare mutations without drift → house_of_cards_G (error O(μ)); many small
    # mutations without drift → G A⁻¹ G = 2L M (s = A: the weak-mutation form of gaussian_G_loci)
    pr = twotrait(L = 20, μ = 1e-7, MutationScalar = 1e-8, InitFitnessSD = 3.0)   # α² ≪ A: no rate capped at 1
    hc = house_of_cards_G(pr)
    check(@sprintf("μ → 0, no drift: house_of_gauss_G → house_of_cards_G (rel. diff. %.1e at μ = 1e-7, α² ≪ A)",
                   maximum(abs, house_of_gauss_G(pr; drift = false) - hc) / maximum(abs, hc)),
          isapprox(house_of_gauss_G(pr; drift = false), hc; rtol = 1e-2))
    pr = twotrait(L = 20, μ = 0.4, MutationScalar = 1e-3)
    A = selection_variance_matrix(pr); Ah = QTLSim.symsqrt(A); Aih = inv(Ah)
    Nm = 2pr.NoQTL .* mutational_covariance(pr)
    Gw = Ah * QTLSim.symsqrt(Aih * Nm * Aih) * Ah            # G A⁻¹ G = 2L M
    Gh = house_of_gauss_G(pr; drift = false)
    check(@sprintf("many small mutations, no drift: house_of_gauss_G → G with G A⁻¹ G = 2L·M (rel. diff. %.1e)",
                   maximum(abs, Gh - Gw) / maximum(abs, Gw)), isapprox(Gh, Gw; rtol = 0.02))
    Gl = gaussian_G_loci(pr)
    check(@sprintf("… which is gaussian_G_loci (G(A + G)⁻¹G = 2L·M) only to O(G/A) (rel. diff. %.1e)",
                   maximum(abs, Gw - Gl) / maximum(abs, Gl)), isapprox(Gw, Gl; rtol = 0.2))

    finish(LOG, t0)
end

runmain(main, @__FILE__)
