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
#     spectrum_mutate!, the Gauss–Hermite effects and spectrum_equilibrium;
# (g) the dynamic LD (`ld = :dynamic`): the class recursion's fixed point (Bulmer's with no drift),
#     its derivation from the haplotype moments, D' w̄² = D(1 + s)(1 + t) and Δln D = Δln p_Aq_A +
#     Δln p_Bq_B, the bins' Σ m_b β_b = V_G/g, and spot checks of ld_classes, _ld_betas!,
#     spectrum_ld_step!, spectrum_select and spectrum_equilibrium;
# (h) the beta closure (`solver = :beta`): the ψ-moment equations against the generator, its
#     exactness on the neutral ψ = θ(1 − p), the beta moment map and its inverse, the Hausdorff
#     conditions for k ≤ 4, and spot checks of the beta solver.

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

    # ===================================================================================
    say("\n(g) dynamic linkage disequilibrium (spectrum = (; ld = :dynamic))")
    @variables dc ρc ωc Vv Sv nn γv gg Ff x1 x2 x3 x4 sA tB
    # one class: d' = dr[(1 − ρ)(1 + 2γ) d − ω V²/(2S)]; fixed point
    dr = 1 - 1 / (2nn)
    stepd(d) = dr * ((1 - ρc) * (1 + 2γv) * d - ωc * Vv^2 / (2Sv))
    dstar = -dr * ωc * Vv^2 / (2Sv) / (1 - dr * (1 - ρc) * (1 + 2γv))
    check("LD fixed point d* = −dr ω V²/(2S)/(1 − dr(1 − ρ)(1 + 2γ)), exactly", iszero_sym(stepd(dstar) - dstar))
    stepb(d) = (1 - ρc) * d - ωc * Vv^2 / (2Sv)
    check("no drift, γ = 0: Bulmer's d* = −ω V²/(2Sρ); ρ = ½, ω = 1: d* = −V²/S",
          iszero_sym(stepb(-ωc * Vv^2 / (2Sv * ρc)) + ωc * Vv^2 / (2Sv * ρc)) &&
          iszero_sym(Symbolics.substitute(-ωc * Vv^2 / (2Sv * ρc), Dict(ωc => 1, ρc => 1 // 2)) + Vv^2 / Sv))
    # derivation of the generation term: selection removes ω V²/(4S) from T and from R (ω D, D = (V/2)²/S);
    # random mating resets R to 0, recombination T ← (1 − ρ)T + ρR: d' = 2T' = (1 − ρ)d − ω V²/(2S)
    @variables Tc
    Ts = Tc - ωc * Vv^2 / (4Sv); Rs = -ωc * Vv^2 / (4Sv)      # D = H s⁻¹ H, H = V/2
    check("haplotype moments: 2[(1 − ρ)T_s + ρR_s] = (1 − ρ)(2T) − ω V²/(2S)",
          iszero_sym(2 * ((1 - ρc) * Ts + ρc * Rs) - ((1 - ρc) * 2Tc - ωc * Vv^2 / (2Sv))))
    # the genic variance given: V = g + D with D = −V² F/(2S) has the positive root
    g_, F_, S_ = 0.37, 41.0, 3.3
    vr = (-1 + sqrt(1 + 4g_ * F_ / (2S_))) / (2F_ / (2S_))
    check("V = g − V² F/(2S): the root (−1 + √(1 + 4gF/(2S)))/(F/S) (numeric)",
          isapprox(vr, g_ - vr^2 * F_ / (2S_); rtol = 1e-12))
    # the infinitesimal limit: the genic class loses β²g²/(2LS) = V²/(2LS) (the multilocus genic loss)
    @variables βv Lv
    check("genic loss with β = V/g: β² g²/(2LS) = V²/(2LS)", iszero_sym(Symbolics.substitute(βv^2 * gg^2 / (2Lv * Sv), Dict(βv => Vv / gg)) - Vv^2 / (2Lv * Sv)))
    # two-locus haploid selection with multiplicative fitnesses (haplotypes AB, Ab, aB, ab):
    # x_i' = x_i w_i/w̄, so D' w̄² = x₁w₁x₄w₄ − x₂w₂x₃w₃ = D (1 + s)(1 + t)
    @variables ε
    w = [(1 + sA) * (1 + tB), 1 + sA, 1 + tB, 1]
    check("selection without epistasis: D' w̄² = D(1 + s)(1 + t), exactly",
          iszero_sym(x1 * w[1] * x4 * w[4] - x2 * w[2] * x3 * w[3] - (x1 * x4 - x2 * x3) * (1 + sA) * (1 + tB)))
    # to first order (s, t → εs, εt; x₄ = 1 − x₁ − x₂ − x₃): D'/D = (1 + s)(1 + t)/w̄², so
    # Δln D = s(1 − 2p_A) + t(1 − 2p_B), and p_A' = p_A(1 + s)(…)/w̄ gives Δln(p_Aq_A) = s(1 − 2p_A)
    X4 = 1 - x1 - x2 - x3
    we = [(1 + ε * sA) * (1 + ε * tB), 1 + ε * sA, 1 + ε * tB, 1]
    wbar = x1 * we[1] + x2 * we[2] + x3 * we[3] + X4 * we[4]
    pA = x1 + x2; pB = x1 + x3
    c1 = coeffj((1 + ε * sA) * (1 + ε * tB) / wbar^2, ε, 1)
    check("first order: Δln D = s(1 − 2p_A) + t(1 − 2p_B)", iszero_sym(c1 - (sA * (1 - 2pA) + tB * (1 - 2pB))))
    pAn = (x1 * we[1] + x2 * we[2]) / wbar
    c1pq = coeffj(Symbolics.substitute(pAn * (1 - pAn) / (pA * (1 - pA)), Dict(tB => 0)), ε, 1)
    check("first order: Δln(p_A q_A) = s(1 − 2p_A), so Δln D = Δln p_Aq_A + Δln p_Bq_B", iszero_sym(c1pq - sA * (1 - 2pA)))
    # bins: Σ_b m_b β_b = V_G/g with β_b = (g m_b + LD within b + ½ LD with other bins)/(g m_b) (two bins)
    @variables m1 d11 d12 d22
    β1 = (gg * m1 + d11 + d12 / 2) / (gg * m1); β2 = (gg * (1 - m1) + d22 + d12 / 2) / (gg * (1 - m1))
    check("two bins: m₁β₁ + m₂β₂ = (g + D)/g", iszero_sym(m1 * β1 + (1 - m1) * β2 - (gg + d11 + d12 + d22) / gg))

    say("  numerical spot checks against src/spectrum.jl")
    pl = Params(PopSize = 1000, Ploidy = 2, NoQTL = 300, PhenotypeVar = 0.05, Heritability = 0.4, MutationScalar = 0.005,
                MutationRatePerGamete = 1e-4, RecombPoissonMean = 5.0, InitFitnessSD = 3.0)
    c = QTLSim.ld_classes(pl; bins = 7)
    d = -0.2 .* rand(rng, length(c.ρ)) ./ length(c.ρ); d0 = copy(d)        # LD below each bin's genic share (no clamping)
    β = zeros(7); g0, S0, n0, γ0 = 0.21, 4.1, 870.0, -3e-4
    QTLSim._ld_betas!(β, c, d, g0)
    hand = zeros(7)
    for k in eachindex(d)
        if c.a[k] == c.b[k]; hand[c.a[k]] += d[k]; else hand[c.a[k]] += d[k] / 2; hand[c.b[k]] += d[k] / 2; end
    end
    check("_ld_betas!: β_b = (g m_b + LD within + ½ LD across)/(g m_b)", isapprox(β, (g0 .* c.m .+ hand) ./ (g0 .* c.m); rtol = 1e-12))
    check("Σ_b m_b β_b = (g + D)/g", isapprox(sum(c.m .* β), (g0 + sum(d)) / g0; rtol = 1e-12))
    QTLSim.spectrum_ld_step!(d, c, β, g0, S0, n0, γ0)
    want = [(1 - 1 / (2n0)) * ((1 - c.ρ[k]) * (1 + 2γ0) * d0[k] - c.ω[k] * β[c.a[k]] * β[c.b[k]] * g0^2 / (2S0)) for k in eachindex(d0)]
    check("spectrum_ld_step! = the recursion, class by class", isapprox(d, want; rtol = 1e-12))
    ω1, ρ1 = QTLSim.pair_classes(pl); c1c = QTLSim.ld_classes(pl; bins = 1)
    check("ld_classes with one bin = pair_classes of distinct loci", isapprox(c1c.ρ, ρ1[2:end]; rtol = 1e-12) && isapprox(c1c.ω, ω1[2:end]; rtol = 1e-12))
    check("ld_classes weights: Σω = 1 − 1/L, Σm = 1 (7 bins)", isapprox(sum(c.ω), 1 - 1 / 300; rtol = 1e-12) && isapprox(sum(c.m), 1; rtol = 1e-12))
    el = spectrum_equilibrium(pl; ld = :dynamic)
    check(@sprintf("spectrum_equilibrium(ld = :dynamic): VG = genic + ld (%.4f = %.4f %+.4f), ld < 0", el.VG, el.genic, el.ld),
          isapprox(el.VG, el.genic + el.ld; rtol = 1e-12) && el.ld < 0)
    sl = el.spectrum
    Vt = QTLSim.spectrum_vtotal(sl); Sl = initFitSD(pl)^2 + envVar(pl) + Vt
    check("the response at the start of a move is V_G λ/S, V_G = genic + LD (symmetric spectrum)",
          isapprox(QTLSim.spectrum_select(sl, 0.4, Sl), Vt * 0.4 / Sl; rtol = 1e-8))
    # at the equilibrium every class is at its fixed point given β, g, S, n and γ of the last generation
    # (to the equilibrium's tolerance)
    QTLSim.spectrum_select(sl, 0.0, Sl)
    gl = QTLSim.spectrum_vg(sl); QTLSim._ld_betas!(sl.β, sl.ldc, sl.d, gl)
    nl = el.nparents; drl = 1 - 1 / (2nl)
    fp = [-drl * sl.ldc.ω[k] * sl.β[sl.ldc.a[k]] * sl.β[sl.ldc.b[k]] * gl^2 / (2Sl) / (1 - drl * (1 - sl.ldc.ρ[k]) * (1 + 2sl.γ))
          for k in eachindex(sl.d)]
    check(@sprintf("equilibrium LD = Σ class fixed points (ratio %.4f)", sum(sl.d) / sum(fp)), isapprox(sum(sl.d), sum(fp); rtol = 0.01))

    # ===================================================================================
    say("\n(h) the beta closure (solver = :beta): ψ-moment equations, exactness, moment map, Hausdorff bounds")
    @variables pp s0 s1 d4 θn
    qq = 1 - pp
    Dp(e) = Symbolics.expand_derivatives(Differential(pp)(e))
    for m in 0:4
        f = pp^(m + 1) * qq                                    # f_m, ∫ f_m φ = P_m; f_m(0) = f_m(1) = 0
        gen = pp * qq * (s0 + s1 * pp) * Dp(f) + pp * qq * d4 * Dp(Dp(f))   # M f′ + D f″, M = pq s, D = pq d₄
        rhs = pp * qq * QTLSim.beta_dP(m, s0, s1, d4, m == 0 ? 0 : pp^(m - 1), pp^m, pp^(m + 1), pp^(m + 2))
        check("generator on f_m = p^(m+1)(1−p) is pq·beta_dP(m, …) with P_k → p^k, m = $m", iszero_sym(gen - rhs))
    end
    # the neutral PRF ψ = θ(1 − p): P_k = θ/((k+1)(k+2)), drift alone leaves P_m (m ≥ 1) unchanged
    # and takes θ d₄ from P₀ (what mutation restores), and the beta closure is exact on it
    Pn(k) = θn / ((k + 1) * (k + 2))
    check("neutral ψ = θ(1−p): drift dP_m = 0 for m = 1…4, dP_0 = −θ d₄",
          all(iszero_sym(QTLSim.beta_dP(m, 0, 0, d4, Pn(m - 1), Pn(m), Pn(m + 1), Pn(m + 2))) for m in 1:4) &&
          iszero_sym(QTLSim.beta_dP(0, 0, 0, d4, 0, Pn(0), Pn(1), Pn(2)) + θn * d4))
    μn3, μn4, αn, βn = QTLSim.beta_closure(1 // 3, 1 // 6)
    check("beta_closure(1/3, 1/6) = (1/10, 1/15, α = 1, β = 2) exactly (neutral μ₃ = 1/10, μ₄ = 1/15)",
          (μn3, μn4, αn, βn) == (1 // 10, 1 // 15, 1 // 1, 2 // 1))
    # the moment map (α, β) ↦ (μ₁, μ₂) and its inverse; μ₃, μ₄ the beta law's
    @variables α_ β_
    ν_ = α_ + β_
    bm(k) = prod((α_ + i) / (ν_ + i) for i in 0:k-1)
    b3, b4, αi, βi = QTLSim.beta_closure(bm(1), bm(2))
    check("beta_closure inverts the moment map: α(μ₁(α,β), μ₂(α,β)) = α, β likewise", iszero_sym(αi - α_) && iszero_sym(βi - β_))
    check("beta_closure's μ₃, μ₄ = Π_{i<k}(α+i)/(ν+i), k = 3, 4", iszero_sym(b3 - bm(3)) && iszero_sym(b4 - bm(4)))
    check("for a beta law μ₂ − μ₁² = αβ/(ν²(ν+1)) and μ₁ − μ₂ = αβ/(ν(ν+1)) (inside the k = 2 bounds)",
          iszero_sym(bm(2) - bm(1)^2 - α_ * β_ / (ν_^2 * (ν_ + 1))) && iszero_sym(bm(1) - bm(2) - α_ * β_ / (ν_ * (ν_ + 1))))
    # Hausdorff conditions for k ≤ 4 (Hankel determinants of μ and of the shifted (μ_k − μ_{k+1})):
    # k = 2: μ₂ − μ₁² ≥ 0, μ₁ − μ₂ ≥ 0; k = 3: μ₁μ₃ − μ₂² ≥ 0, (1 − μ₁)(μ₂ − μ₃) − (μ₁ − μ₂)² ≥ 0;
    # k = 4: det[1 μ₁ μ₂; μ₁ μ₂ μ₃; μ₂ μ₃ μ₄] ≥ 0, det[μ₁−μ₂ μ₂−μ₃; μ₂−μ₃ μ₃−μ₄] ≥ 0
    haus(μ) = (μ[2] - μ[1]^2, μ[1] - μ[2], μ[1] * μ[3] - μ[2]^2, (1 - μ[1]) * (μ[2] - μ[3]) - (μ[1] - μ[2])^2,
               det([1 μ[1] μ[2]; μ[1] μ[2] μ[3]; μ[2] μ[3] μ[4]]), (μ[1] - μ[2]) * (μ[3] - μ[4]) - (μ[2] - μ[3])^2)
    okm, okb = true, true
    for _ in 1:2000
        k = rand(rng, 1:4); x = rand(rng, k); w = rand(rng, k); w ./= sum(w)          # a measure of k atoms on [0, 1]
        okm &= all(h -> h >= -1e-12, haus([sum(w .* x .^ j) for j in 1:4]))
        al, be = exp(runif(-3, 3)), exp(runif(-3, 3))
        m1, m2 = al / (al + be), al * (al + 1) / ((al + be) * (al + be + 1))
        c3, c4, _, _ = QTLSim.beta_closure(m1, m2)
        okb &= all(h -> h >= -1e-12, haus([m1, m2, c3, c4]))
    end
    check("the Hausdorff conditions (k ≤ 4) hold for 2000 random atomic measures on [0, 1]", okm)
    check("the closed moments (μ₁, μ₂, μ₃, μ₄) of 2000 random beta laws satisfy them", okb)
    # numerical spot checks of the beta solver against its definition
    pb = Params(PopSize = 1000, Ploidy = 2, NoQTL = 200, PhenotypeVar = 0.05, Heritability = 0.4, MutationScalar = 0.005,
                MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 3.0)
    eb = spectrum_equilibrium(pb; solver = :beta); sb = eb.spectrum
    gb = QTLSim.spectrum_vg(sb)
    check("spectrum_vg (beta) = Σ 2a² P₀", isapprox(gb, sum(2 .* sb.a .^ 2 .* sb.P[:, 1]); rtol = 1e-14))
    Sb_, λb, cb = initFitSD(pb)^2 + envVar(pb) + gb, 0.37, gb / (2pb.NoQTL)
    dzb = QTLSim.spectrum_select(sb, λb, Sb_)
    hand = sum(begin
                   e0 = QTLSim.spectrum_seff(sb.a[j], 0.0, λb, Sb_, cb); e1 = QTLSim.spectrum_seff(sb.a[j], 1.0, λb, Sb_, cb) - e0
                   2sb.a[j] * (e0 * sb.P[j, 1] + e1 * sb.P[j, 2])
               end for j in eachindex(sb.a))
    check("spectrum_select (beta) = Σ 2a ∫ s ψ with spectrum_seff's s(p) = s₀ + s₁p", isapprox(dzb, hand; rtol = 1e-12))
    eg = spectrum_equilibrium(pb)
    check(@sprintf("beta equilibrium V_G / grid's = %.5f (within 1%%)", eb.VG / eg.VG), abs(eb.VG / eg.VG - 1) < 0.01)
    sc = deepcopy(sb); QTLSim._beta_to_grid!(sc)
    check("handing the state to the grid keeps V_G exactly (relative 1e-12)", sc.solver === :grid && isapprox(QTLSim.spectrum_vg(sc), gb; rtol = 1e-12))
    s2 = deepcopy(sb); P0 = copy(s2.P); QTLSim.spectrum_mutate!(s2, 1000)
    x0 = 1 / 2000; inc = 2000 * pb.NoQTL * pb.MutationRatePerGamete .* s2.g .* x0 * (1 - x0)
    check("spectrum_mutate! (beta) adds 2N L μ g f_m(1/(2N)) to P_m", isapprox(s2.P, P0 .+ inc .* [1 x0 x0^2]; rtol = 1e-12))

    finish(LOG, t0)
end

runmain(main, @__FILE__)
