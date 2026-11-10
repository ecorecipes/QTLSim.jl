# Symbolic checks of the allele-frequency spectrum model (`loci = :spectrum`, src/spectrum.jl
# and ext/QTLSimCatlabExt/spectrum.jl), and numerical spot checks against its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_spectrum.jl
# Writes validate/symbolic/check_spectrum.log; exits 1 if any check fails.
#
# (a) the per-allele selection E[Δp] = pq[aλ/S − (1 − λ²/S)(a² + 2c)(1 − 2p)/(2S)]:
#     Wright's formula with ln W̄ = −λ²/(2S) − ln(S)/2 gives it exactly (at c = 0), and it is
#     the exact change under viability selection on a diploid locus in a Gaussian
#     background to second order in the effect a (Hayward and Sella 2022, eq. 7), at any
#     lag; near p = 0 it is Kopp and Hermisson's s(α) = α(2λ − α)/(2S) up to a²λ²/(2S²);
#     the 2c term gives each gene the Gaussian loss c²/S of the rare-alleles closure, and
#     the a² term the fourth-moment loss Σa⁴p/(2S), which is α²c/(2S) at the house of cards;
# (b) the neutral Poisson random field: E[pq] decays by 1 − 1/(2N) a generation, the
#     diffusion's Green's function gives ∫pq dt = 2N p₀q₀, so V_G = (2N − 1)V_m ≈ 2N V_m;
#     the lowering's implicit Euler step gives (2N − 1/(2N))V_m;
# (c) the symmetric-spectrum identity Σ 2a E[Δp] = V_G λ/S when φ(a,·) = φ(−a,·), and the
#     directional part of Σ 2a E[Δp] being V_G λ/S in general;
# (d) the sweep variance: fixation probability 2aλ/S times ∫2a²pq dt = 2a S/λ is 4a²,
#     whatever λ, so the sweeps hold 2N_e V_m;
# (e) Barton's linked-selection factor N/N_e = 1 + b²V_G Σ_c ω_c/(2(ρ_c + δ)²);
# (f) numerical spot checks of spectrum_seff, spectrum_select, spectrum_vg, spectrum_linked,
#     spectrum_mutate!, the Gauss–Hermite effects and spectrum_equilibrium.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_spectrum.log")
const rng = Xoshiro(20261109)
runif(a, b) = runif(rng, a, b)

# the documented E[Δp]
dp_formula(a, p, λ, S, c) = p * (1 - p) * (a * λ / S - (1 - λ^2 / S) * (a^2 + 2c) * (1 - 2p) / (2S))

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_spectrum.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the selection of one allele, E[Δp]")
    @variables a p θ zb0 Vr ω2 ε λ Sb S c
    q = 1 - p
    # Wright's formula Δp = (pq/2) ∂ln W̄/∂p, with z̄ = z̄₀ + 2ap and V_P = V_r + 2a²pq (Hardy–Weinberg)
    zbar = zb0 + 2a * p
    Sx = ω2 + Vr + 2a^2 * p * q
    lnW = -(θ - zbar)^2 / (2Sx) - log(Sx) / 2
    wright = p * q / 2 * Symbolics.expand_derivatives(Differential(p)(lnW))
    check("Wright's formula on ln W̄ = −λ²/(2S) − ln S/2 is E[Δp] (c = 0), exactly",
          iszero_sym(wright - dp_formula(a, p, θ - zbar, Sx, 0)))
    # the exact change under viability selection: genotypes with g = 0, 1, 2 copies of an
    # allele of effect εa in a Gaussian background of lag L_b and width S_b, fitness
    # w_g ∝ exp(−(L_b − gεa)²/(2S_b)); with λ = L_b − 2εap the population's lag and
    # S = S_b + 2ε²a²pq its S. w_g/w_0 = exp(x_g), x_g = O(ε), expanded to x³.
    xg(g) = (2g * ε * a * (λ + 2ε * a * p) - g^2 * ε^2 * a^2) / (2Sb)
    rg(g) = (x = xg(g); 1 + x + x^2 / 2 + x^3 / 6)
    exact = (p * (rg(2) - rg(1)) + q * (rg(1) - rg(0))) / (p^2 * rg(2) + 2p * q * rg(1) + q^2 * rg(0))   # Δp/(pq)
    form = dp_formula(ε * a, p, λ, Sb + 2ε^2 * a^2 * p * q, 0) / (p * q)
    dif = exact - form
    check("diploid viability selection, exact: Δp/(pq) = E[Δp]/(pq) + O(ε³) (order ε⁰, ε¹)",
          iszero_sym(coeffj(dif, ε, 0)) && iszero_sym(coeffj(dif, ε, 1)))
    check("… order ε²: −(1 − λ²/S)a²(1 − 2p)/(2S), at any lag (Hayward–Sella eq. 7)", iszero_sym(coeffj(dif, ε, 2)))
    check("near p = 0 (c = 0): Kopp–Hermisson a(2λ − a)/(2S) + a²λ²/(2S²)",
          iszero_sym(Symbolics.substitute(dp_formula(a, p, λ, S, 0) / p, Dict(p => 0)) -
                     (a * (2λ - a) / (2S) + a^2 * λ^2 / (2S^2))))
    # the 2c term: two sites of one gene (a₁, p₁), (a₂, p₂), c = Σ a²pq; selection at the optimum
    # changes c by Σ a²(1 − 2p)Δp; its part linear in each p_i is −Σa⁴p/(2S) − c²/S
    @variables a1 a2 p1 p2
    cg = a1^2 * p1 * (1 - p1) + a2^2 * p2 * (1 - p2)
    Δc = a1^2 * (1 - 2p1) * dp_formula(a1, p1, 0, S, cg) + a2^2 * (1 - 2p2) * dp_formula(a2, p2, 0, S, cg)
    target = -(a1^4 * p1 + a2^4 * p2) / (2S) - (a1^2 * p1 + a2^2 * p2)^2 / S
    D = Δc - target
    mixed = Symbolics.expand_derivatives(Differential(p1)(Differential(p2)(D)))
    check("rare alleles, 2c term: Δc = −Σa⁴p/(2S) − c²/S to first order in each p",
          iszero_sym(Symbolics.substitute(D, Dict(p1 => 0, p2 => 0))) &&
          iszero_sym(Symbolics.substitute(Symbolics.expand_derivatives(Differential(p1)(D)), Dict(p1 => 0, p2 => 0))) &&
          iszero_sym(Symbolics.substitute(Symbolics.expand_derivatives(Differential(p2)(D)), Dict(p1 => 0, p2 => 0))) &&
          iszero_sym(Symbolics.substitute(mixed, Dict(p1 => 0, p2 => 0))))
    # at the house of cards p(a) ∝ g(a)/a², so Σ a⁴p / Σ a²p = E[a²] = α²: the a² term is α²c/(2S)
    @variables K g1 g2
    pa1, pa2 = K * g1 / a1^2, K * g2 / a2^2
    check("house of cards p ∝ g/a²: Σa⁴p = (Σga²/Σg)·Σa²p, i.e. α²c",
          iszero_sym((a1^4 * pa1 + a2^4 * pa2) * (g1 + g2) - (g1 * a1^2 + g2 * a2^2) * (a1^2 * pa1 + a2^2 * pa2)))

    # ===================================================================================
    say("\n(b) the neutral Poisson random field")
    @variables N L μ α2 p0
    Ep2 = p^2 + p * q / (2N)                   # E[p′²] under binomial sampling of 2N genes
    check("Wright–Fisher: E[p′q′] = (1 − 1/(2N)) pq", iszero_sym((p - Ep2) - (1 - 1 / (2N)) * p * q))
    check("Σ_t (1 − 1/(2N))^t pq = 2N pq", iszero_sym(p * q / (1 - (1 - 1 / (2N))) - 2N * p * q))
    u = 2N * p * q
    check("diffusion: u = 2N pq solves (pq/(4N))u″ = −pq, u(0) = u(1) = 0",
          iszero_sym(p * q / (4N) * Symbolics.expand_derivatives(Differential(p)(Differential(p)(u))) + p * q) &&
          iszero_sym(Symbolics.substitute(u, Dict(p => 0))) && iszero_sym(Symbolics.substitute(u, Dict(p => 1))))
    Vm = 2L * μ * α2                           # V_m = Ploidy·L·μ·α², diploid
    nmut = 2N * L * μ                          # new mutations per generation, each at p₀ = 1/(2N)
    VG = nmut * 2α2 * 2N * Symbolics.substitute(p * q, Dict(p => 1 / (2N)))
    check("V_G = 2NLμ · 2α² · 2N p₀q₀ = (2N − 1) V_m", iszero_sym(VG - (2N - 1) * Vm))
    # implicit Euler: E[pq] is divided by 1 + 1/(2N) a step, so Σ_t = (2N + 1) p₀q₀
    VGie = nmut * 2α2 * Symbolics.substitute(p * q, Dict(p => 1 / (2N))) / (1 - 1 / (1 + 1 / (2N)))
    check("implicit Euler step: V_G = (2N − 1/(2N)) V_m", iszero_sym(VGie - (2N - 1 / (2N)) * Vm))

    # ===================================================================================
    say("\n(c) the response of the mean")
    @variables pp
    pair = 2a * dp_formula(a, pp, λ, S, c) + 2(-a) * dp_formula(-a, pp, λ, S, c)
    check("symmetric pair (±a at one p): Σ 2a E[Δp] = (Σ 2a²pq) λ/S",
          iszero_sym(pair - 2 * 2a^2 * pp * (1 - pp) * λ / S))
    dirpart = 2a1 * p1 * (1 - p1) * a1 * λ / S + 2a2 * p2 * (1 - p2) * a2 * λ / S
    check("directional part of Σ 2a E[Δp] is V_G λ/S, V_G = Σ 2a²pq",
          iszero_sym(dirpart - (2a1^2 * p1 * (1 - p1) + 2a2^2 * p2 * (1 - p2)) * λ / S))

    # ===================================================================================
    say("\n(d) the sweep variance")
    @variables s t Ne
    pt = p0 * exp(s * t) / (1 - p0 + p0 * exp(s * t))
    dpt = Symbolics.expand_derivatives(Differential(t)(pt))
    check("logistic sweep: dp/dt = s pq", iszero_sym(simp(dpt - s * pt * (1 - pt))))
    check("so d(p/s)/dt = pq: ∫ pq dt over the sweep = (1 − p₀)/s ≈ 1/s", iszero_sym(simp(dpt / s - pt * (1 - pt))))
    sa = a * λ / S                             # the heterozygote's selection coefficient
    check("∫ 2a²pq dt = 2a²/s = 2a S/λ", iszero_sym(2a^2 / sa - 2a * S / λ))
    check("Haldane 2s × 2aS/λ = 4a², whatever λ and S", iszero_sym(2sa * 2a * S / λ - 4a^2))
    # summed over effects: 2NLμ mutants, the half with aλ > 0 fixing at 2s·N_e/N
    check("2NLμ · (N_e/N) · E[4a²; a > 0] = 2N_e V_m, with E[a²; a > 0] = α²/2",
          iszero_sym(2N * L * μ * (Ne / N) * 4 * (α2 / 2) - 2Ne * Vm))

    # ===================================================================================
    say("\n(e) linked selection")
    @variables n b V σω cdrift r
    check("Var(Δp) = pq/(2n)·N/N_e with N/N_e = 1 + b²VΣ ⇔ pq(1/(2n) + c), c = b²VΣ/(2n)",
          iszero_sym((1 + b^2 * V * σω) / (2n) - (1 / (2n) + b^2 * V * σω / (2n))))
    # Barton: L − 1 other loci, each with var(W) = b²V/L, at split probability r: Σ var(W)/(2r²)
    check("Barton per locus: (L − 1)·(b²V/L)/(2r²) = b²V·(1 − 1/L)/(2r²) (weights Σω = 1 − 1/L)",
          iszero_sym((L - 1) * (b^2 * V / L) / (2r^2) - b^2 * V * (1 - 1 / L) / (2r^2)))
    # the selection gradient of Gaussian selection, ∂ln W̄/∂z̄ = λ/S: relative fitness w ≈ 1 + b(z − z̄)
    @variables zz
    lnWz = -(θ - zz)^2 / (2S) - log(S) / 2
    check("∂ln W̄/∂z̄ = λ/S = b (additive variance in relative fitness b²V)",
          iszero_sym(Symbolics.expand_derivatives(Differential(zz)(lnWz)) - (θ - zz) / S))

    # ===================================================================================
    say("\n(f) numerical spot checks against src/spectrum.jl")
    for trial in 1:2
        L0 = rand(rng, 10:40); N0 = rand(rng, 100:400); μ0 = 10^runif(-3.5, -2.5)
        pr = Params(Ploidy = 2, NoQTL = L0, PopSize = N0, MutationRatePerGamete = μ0, MutationScalar = 10^runif(-3, -2),
                    RecombPoissonMean = runif(0.5, 3.0), PhenotypeVar = runif(0.5, 2.0), Heritability = runif(0.3, 0.7))
        eq = spectrum_equilibrium(pr)
        sp = eq.spectrum
        α20 = QTLSim.mutation_effect_variance(pr)
        Vg = sum(2sp.a[j]^2 * sp.p[i] * (1 - sp.p[i]) * sp.n[j, i] for j in eachindex(sp.a), i in eachindex(sp.p))
        check(@sprintf("spectrum_vg = Σ 2a²pq n (L=%d N=%d μ=%.1e)", L0, N0, μ0), isapprox(QTLSim.spectrum_vg(sp), Vg; rtol = 1e-12))
        check("Gauss–Hermite effects: Σg = 1, Σga² = α², Σga⁴ = 3α⁴, symmetric",
              isapprox(sum(sp.g), 1; rtol = 1e-12) && isapprox(sum(sp.g .* sp.a .^ 2), α20; rtol = 1e-10) &&
              isapprox(sum(sp.g .* sp.a .^ 4), 3α20^2; rtol = 1e-9) && isapprox(sort(sp.a), -reverse(sort(sp.a)); atol = 1e-12 * sqrt(α20)))
        check("Σ_{a>0} g a² = α²/2 (the half that fixes under a lag)", isapprox(sum(sp.g[j] * sp.a[j]^2 for j in eachindex(sp.a) if sp.a[j] > 0), α20 / 2; rtol = 1e-10))
        λ0, S0 = runif(-1.0, 1.0), runif(5.0, 50.0)
        c0 = Vg / (2L0)
        a0, p00 = sp.a[end], 0.3
        check("spectrum_seff = E[Δp]/(pq)", isapprox(QTLSim.spectrum_seff(a0, p00, λ0, S0, c0) * p00 * (1 - p00),
                                                     dp_formula(a0, p00, λ0, S0, c0); rtol = 1e-12))
        # the equilibrium spectrum at λ = 0 is symmetric in a (up to rounding): Σ 2aE[Δp] = V_G λ/S
        sym = maximum(abs, sp.n .- reverse(sp.n; dims = 1)) / maximum(sp.n)
        dz = QTLSim.spectrum_select(sp, λ0, S0)
        check(@sprintf("symmetric spectrum (asymmetry %.0e): spectrum_select = V_G λ/S", sym), isapprox(dz, Vg * λ0 / S0; rtol = 1e-8))
        # an asymmetric spectrum: Σ 2a E[Δp] = V_G λ/S minus the skew term
        sk = deepcopy(sp); sk.n = sp.n .* (1 .+ 0.5 .* sign.(sp.a))
        Vk = QTLSim.spectrum_vg(sk); ck = Vk / (2L0)
        skew = sum(2sk.a[j] * sk.p[i] * (1 - sk.p[i]) * (1 - λ0^2 / S0) * (sk.a[j]^2 + 2ck) * (1 - 2sk.p[i]) / (2S0) * sk.n[j, i]
                   for j in eachindex(sk.a), i in 2:length(sk.p)-1)
        check("asymmetric spectrum: spectrum_select = V_G λ/S − Σ 2a pq (1−λ²/S)(a²+2c)(1−2p)/(2S) n",
              isapprox(QTLSim.spectrum_select(sk, λ0, S0), Vk * λ0 / S0 - skew; rtol = 1e-10))
        # linked selection: the code against the formula, and against the sum over the pairs themselves
        npar = runif(50.0, N0)
        QTLSim.spectrum_select(sp, λ0, S0)
        δ = 1 / (2npar) + Vg / S0
        b2V = (λ0 / S0)^2 * Vg
        form = b2V * sum(sp.ω[k] / (2(sp.ρ[k] + δ)^2) for k in 2:length(sp.ω)) / (2npar)
        check("spectrum_linked = b²V Σ_{c≥2} ω_c/(2(ρ_c + δ)²)/(2n)", isapprox(QTLSim.spectrum_linked(sp, npar), form; rtol = 1e-12))
        check("pair classes of distinct loci weigh 1 − 1/L", isapprox(sum(sp.ω[2:end]), 1 - 1 / L0; rtol = 1e-12))
        ρp, wp = QTLSim.pair_splits(L0, pr.RecombPoissonMean, pr.CrossoverModel)
        perpair = b2V * (1 - 1 / L0) * sum(wp .* (1 ./ (2 .* (ρp .+ δ) .^ 2))) / (2npar)
        check(@sprintf("… = Barton's sum over the pairs themselves (binned classes: %.1e)", perpair / form - 1),
              isapprox(form, perpair; rtol = 0.03))
        # mutation: number and mean frequency of the new sites exact
        before = copy(sp.n)
        QTLSim.spectrum_mutate!(sp, N0)
        add = sp.n .- before
        check("spectrum_mutate!: 2NLμ sites, mean frequency 1/(2N)",
              isapprox(sum(add), 2N0 * L0 * μ0; rtol = 1e-12) &&
              isapprox(sum(vec(sum(add; dims = 1)) .* sp.p) / sum(add), 1 / (2N0); rtol = 1e-10))
        # equilibrium close to the stochastic house of Gauss (loci in linkage equilibrium): an
        # approximation (infinite sites against the per-gene closure), so a loose tolerance
        hg = house_of_gauss_vg(pr)
        check(@sprintf("spectrum_vg ≈ house_of_gauss_vg (approximate; ratio %.4f)", eq.VG / hg), isapprox(eq.VG, hg; rtol = 0.08))
    end
    # neutral: the implicit Euler equilibrium (2N − 1/(2N))V_m, as an expected sample variance ×N/(N − 1)
    for N0 in (100, 250)
        pr = Params(Ploidy = 2, NoQTL = 20, PopSize = N0, MutationRatePerGamete = 1e-3, RecombPoissonMean = 2.0, Drift = true)
        Vm0 = muteVar(pr)
        v = spectrum_vg(pr)
        check(@sprintf("neutral, N = %d: spectrum_vg·(N−1)/N = (2N − 1/(2N))V_m (ratio %.6f)", N0, v * (N0 - 1) / N0 / ((2N0 - 1 / (2N0)) * Vm0)),
              isapprox(v * (N0 - 1) / N0, (2N0 - 1 / (2N0)) * Vm0; rtol = 1e-4))
    end

    finish(LOG, t0)
end

runmain(main, @__FILE__)
