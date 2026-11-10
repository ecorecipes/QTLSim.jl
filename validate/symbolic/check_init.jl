# Symbolic checks of the initialisation from theory (src/init.jl: `stationary_spectrum`,
# `equilibrium_genomes`, `equilibrium_world`), and numerical spot checks against its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_init.jl
# Writes validate/symbolic/check_init.log; exits 1 if any check fails.
#
# (a) the density sampled, φ(x) ∝ e^{−γx(1−x)} I(x)/(x(1−x)) with I(x) = ∫ₓ¹e^{γy(1−y)}dy, is a
#     stationary solution of the diffusion with the code's drift and selection terms,
#     M = −σx(1−x)(1−2x), D = x(1−x)/(4n), γ = 4nσ: its probability flux Mφ − (Dφ)′ is the
#     constant 1/(4n) (not zero: sites flow to fixation), the flux of fixations equals the
#     input of mutations times their fixation probability to first order in the entry
#     frequency 1/(2N), the fixation probability solves the backward equation, and at γ = 0
#     the density is the neutral 1/x;
# (b) its variance: F_γ(x) + F_γ(1−x) = 1, so ∫2a²x(1−x)φ dx = 4nLμ g(a) a² J(γ),
#     J = ∫₀¹e^{−γx(1−x)}dx; J = 1 − γ/6 + γ²/60 − … and J ≈ (2/γ)(1 + 2/γ + 12/γ² + 120/γ³);
# (c) the limits of V = 4nLμ E[a²J]: neutral 2nV_m exactly; with J → 2/γ (leading order in
#     1/γ) the house of cards 4LμS when c ≪ a² and the Gaussian V² = 2L V_m S when a² ≪ c;
#     replacing J by 1/(1 + γ/2), which has the same two limits, and a² by α² gives
#     Hermisson and Wagner's equation A4 (`house_of_gauss_vg`) exactly;
# (d) the sampler: sites at k copies among m chromosomes have the neutral mean θ/k and
#     x | k ~ Beta(k, m − k + 1); with k copies placed at random the expected sample variance
#     of the genotypes is 2a²·p̂q̂·m/(m−1), whose expectation over k ~ Bin(m, x) is 2a²x(1−x),
#     exactly;
# (e) the mean genotype's AR(1) around the optimum has stationary variance S/(n(2 − b));
# (f) numerical: `sojourn_mass` and `sojourn_ratio` against brute-force quadrature, the
#     acceptance ratio ≤ 1, `stationary_spectrum` against its own equation, its limits,
#     `spectrum_vg` (the same model on a grid, by implicit Euler steps) and
#     `house_of_gauss_vg`, and the sampler's Monte Carlo means (sites, genic variance, mean).

include("common.jl")

const LOG = joinpath(@__DIR__, "check_init.log")
const rng = Xoshiro(20261110)
runif(a, b) = runif(rng, a, b)

# composite Simpson
function simpson(f, a, b; n = 20_000)
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:n-1
        s += (isodd(i) ? 4 : 2) * f(a + i * h)
    end
    s * h / 3
end
# J(γ) by brute force: 2∫₀^½ e^{−γx(1−x)}dx, the boundary layer of width 1/γ apart
function J_brute(γ)
    f(x) = exp(-γ * x * (1 - x))
    b = min(0.5, 60 / max(γ, 1e-300))
    2 * (simpson(f, 0.0, b) + (b < 0.5 ? simpson(f, b, 0.5) : 0.0))
end
# F_γ(x) = ∫ₓ¹e^{γy(1−y)}dy/∫₀¹, with e^{−γ/4} taken out of both
function F_brute(γ, x; n = 200_000)
    f(y) = exp(γ * (y * (1 - y) - 0.25))
    simpson(f, x, 1.0; n) / simpson(f, 0.0, 1.0; n)
end
R_brute(γ, x) = exp(-γ * x * (1 - x)) * F_brute(γ, x) / (1 - x)

# ∫₀¹ of a polynomial in x (exact, term by term)
polyint01(e, x, deg) = sum(coeffj(e, x, j) / (j + 1) for j in 0:deg)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_init.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the sampled density is stationary for the code's drift and selection")
    # E stands for e^{γx(1−x)} and Iv for I(x) = ∫ₓ¹E dy: dE/dx = γ(1−2x)E, dI/dx = −E
    @variables x γ n σ E Iv I0 N L μ g
    Dx = Differential(x); DE = Differential(E); DI = Differential(Iv)
    dtot(e) = Symbolics.expand_derivatives(Dx(e)) + Symbolics.expand_derivatives(DE(e)) * γ * (1 - 2x) * E -
              Symbolics.expand_derivatives(DI(e)) * E
    φ = Iv / (E * x * (1 - x))
    M = -(γ / (4n)) * x * (1 - x) * (1 - 2x)           # E[Δx], σ = γ/(4n): spectrum_seff at λ = 0
    Dc = x * (1 - x) / (4n)                              # Var(Δx)/2 = x(1−x)/(4n)
    flux = M * φ - dtot(Dc * φ)
    check("probability flux Mφ − (Dφ)′ of φ = e^{−γx(1−x)}I(x)/(x(1−x)) is the constant 1/(4n)", iszero_sym(flux - 1 / (4n)))
    check("so φ is stationary: −(Mφ)′ + (Dφ)″ = 0", iszero_sym(-dtot(M * φ) + dtot(dtot(Dc * φ))))
    check("the code's selection: spectrum_seff(a, x, 0, S, c)·x(1−x) = −σx(1−x)(1−2x), σ = (a²+2c)/(2S)",
          let a = runif(-2, 2), xx = runif(0, 1), S = runif(1, 50), c = runif(0, 1)
              isapprox(QTLSim.spectrum_seff(a, xx, 0.0, S, c) * xx * (1 - xx), -(a^2 + 2c) / (2S) * xx * (1 - xx) * (1 - 2xx); rtol = 1e-13)
          end)
    # the fixation probability u(x) = ∫₀ˣE dy/I0: u′ = E/I0, u″ = γ(1−2x)E/I0
    u1, u2 = E / I0, γ * (1 - 2x) * E / I0
    check("fixation probability u = ∫₀ˣe^{γy(1−y)}dy/I(0) solves the backward equation M u′ + D u″ = 0", iszero_sym(M * u1 + Dc * u2))
    # input 2NLμ g sites per generation at 1/(2N), each fixing with probability u(1/(2N)) = (1/(2N))·u′(0) + O(N⁻²), u′(0) = 1/I0;
    # the density 4nLμ g φ/I0 carries the flux 4nLμ g/I0 · 1/(4n)
    check("flux of fixations Lμg/I(0) = input 2NLμg × fixation probability, to first order in 1/(2N)",
          iszero_sym(4n * L * μ * g / I0 * (1 / (4n)) - 2N * L * μ * g * (1 / (2N)) * Symbolics.substitute(u1, Dict(E => 1))))
    check("neutral (γ = 0, E = 1, I = 1 − x): φ = 1/x, the density 4nLμg/x", iszero_sym(Symbolics.substitute(φ, Dict(E => 1, Iv => 1 - x)) - 1 / x))

    # ===================================================================================
    say("\n(b) the variance of the density")
    @variables y a
    check("y(1−y) is symmetric about ½, so I(x) + I(1−x) = I(0) and F_γ(x) + F_γ(1−x) = 1",
          iszero_sym(Symbolics.substitute(y * (1 - y), Dict(y => 1 - y)) - y * (1 - y)))
    ok = true
    for γv in (0.3, 4.0, 37.0, 150.0)
        xs = [runif(0.02, 0.98) for _ in 1:3]
        ok &= all(abs(F_brute(γv, xv) + F_brute(γv, 1 - xv) - 1) < 1e-9 for xv in xs)
        v = simpson(xv -> 2 * exp(-γv * xv * (1 - xv)) * F_brute(γv, xv; n = 10_000), 0.0, 1.0; n = 4000)
        ok &= isapprox(v, J_brute(γv); rtol = 1e-6)
    end
    check("numerically: ∫₀¹ 2x(1−x)·e^{−γx(1−x)}F_γ/(x(1−x)) dx = J(γ) = ∫₀¹e^{−γx(1−x)}dx", ok)
    uu = x * (1 - x)
    series = polyint01(Symbolics.expand(1 - γ * uu + γ^2 * uu^2 / 2), x, 4)
    check("J(γ) = 1 − γ/6 + γ²/60 + O(γ³) (the integral of the expanded exponential)", iszero_sym(series - (1 - γ / 6 + γ^2 / 60)))
    check("QTLSim.sojourn_mass(γ) = 1 − γ/6 + γ²/60 − γ³/840 for small γ (to 1e-10 at γ = 0.01)",
          abs(QTLSim.sojourn_mass(0.01) - (1 - 0.01 / 6 + 0.01^2 / 60 - 0.01^3 / 840)) < 1e-10)
    check("QTLSim.sojourn_mass(γ) ≈ (2/γ)(1 + 2/γ + 12/γ² + 120/γ³) for large γ (to 1e-6 at γ = 400, 1e4)",
          all(γv -> isapprox(QTLSim.sojourn_mass(γv), 2 / γv * (1 + 2 / γv + 12 / γv^2 + 120 / γv^3); rtol = 1e-6), (400.0, 1e4)))

    # ===================================================================================
    say("\n(c) the limits of V = 4nLμ E_a[a² J(γ(a))], γ = 2n(a² + 2c)/S, c = V/(2L)")
    @variables α2 S V Vm c a2
    check("neutral (J = 1): V = 4nLμα² = 2nV_m with V_m = 2Lμα²", iszero_sym(4n * L * μ * α2 - 2n * (2L * μ * α2)))
    check("house of cards (J → 2/γ, c = 0): a²·2/γ = S/n for every a, so V = 4LμS (leading order in 1/γ)",
          iszero_sym(4n * L * μ * (a2 * 2 / (2n * a2 / S)) - 4L * μ * S))
    γg = 2n * (V / L) / S
    check("Gaussian (J → 2/γ, a² ≪ 2c): V·(4nLμα²·2/γ) = 2L V_m S, i.e. V = √(2L V_m S) (leading order)",
          iszero_sym(V * (4n * L * μ * α2 * 2 / γg) - 2L * (2L * μ * α2) * S))
    γc = 2n * (α2 + 2c) / S
    check("J → 1/(1 + γ/2), a² → α²: c = 2nμα²/(1 + γ/2) is Hermisson–Wagner A4, c²/S + c(α²/(2S) + 1/(2n)) = μα²",
          iszero_sym((c * (1 + γc / 2) - 2n * μ * α2) / (2n) - (c^2 / S + c * (α2 / (2S) + 1 / (2n)) - μ * α2)))

    # ===================================================================================
    say("\n(d) the sampler")
    @variables m k
    # k copies among m = 2N chromosomes, an individual holds two of them: X ∈ {0, 1, 2}
    EX = 2k / m
    EX2 = EX + 2 * k * (k - 1) / (m * (m - 1))
    varX = EX2 - EX^2
    Nn = m / 2
    check("k copies at random: Var(copies per individual) = 2p̂q̂(m−2)/(m−1), p̂ = k/m", iszero_sym(varX - 2 * (k / m) * (1 - k / m) * (m - 2) / (m - 1)))
    check("expected sample variance (divisor N − 1) of the genotypes: 2a²·p̂q̂·m/(m−1)",
          iszero_sym(a^2 * varX * Nn / (Nn - 1) - 2a^2 * (k / m) * (1 - k / m) * m / (m - 1)))
    # k ~ Bin(m, x): E[k] = mx, E[k²] = mx(1−x) + m²x²
    Ek, Ek2 = m * x, m * x * (1 - x) + m^2 * x^2
    check("its expectation over k ~ Bin(m, x) is 2a²x(1−x), exactly", iszero_sym(2a^2 * (Ek / m - Ek2 / m^2) * m / (m - 1) - 2a^2 * x * (1 - x)))
    ok = true
    for kk in 1:5
        integrand = Symbolics.expand(binomial(6, kk) * x^(kk - 1) * (1 - x)^(6 - kk))
        ok &= iszero_sym(polyint01(integrand, x, 6) - 1 // kk)
    end
    check("neutral sites at k copies of m = 6: ∫₀¹(1/x)·Bin(k; m, x)dx = 1/k, k = 1..5", ok)
    mm = big(2000)
    check("… and for m = 2000, every k: C(m,k)·B(k, m−k+1) = 1/k (exact rationals); x | k ~ Beta(k, m−k+1)",
          all(binomial(mm, kk) * (factorial(big(kk - 1)) * factorial(mm - kk) // factorial(mm)) == 1 // kk for kk in (1, 2, 17, 500, 1999)))

    # ===================================================================================
    say("\n(e) the mean genotype around the optimum")
    @variables v b
    vstat = S / (n * (2 - V / S))
    check("AR(1) z̄′ = (1−b)z̄ + ε, b = V/S, Var ε = V/n: the variance S/(n(2−b)) is stationary",
          iszero_sym(Symbolics.substitute((1 - b)^2 * v + V / n - v, Dict(v => vstat, b => V / S))))

    # ===================================================================================
    say("\n(f) numerical checks of src/init.jl")
    γs = (1e-6, 0.01, 0.5, 3.0, 20.0, 79.0, 81.0, 300.0, 5000.0, 2e5)
    check("sojourn_mass(γ) = ∫₀¹e^{−γx(1−x)}dx by brute force, γ from 1e-6 to 2e5 (rel. 1e-9)",
          all(γv -> isapprox(QTLSim.sojourn_mass(γv), J_brute(γv); rtol = 1e-9), γs))
    worst, rmax = 0.0, 0.0
    for γv in (1e-6, 0.01, 0.5, 3.0, 20.0, 300.0, 5000.0), xv in (1e-4, 0.003, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999)
        rb = R_brute(γv, xv)
        rmax = max(rmax, rb)
        rb > 1e-200 && (worst = max(worst, abs(QTLSim.sojourn_ratio(γv, xv) / rb - 1)))
    end
    say(@sprintf("    worst relative difference %.2g; largest ratio %.12f", worst, rmax))
    check("sojourn_ratio(γ, x) = e^{−γx(1−x)}F_γ(x)/(1−x) by brute force (rel. 1e-7)", worst < 1e-7)
    check("the acceptance ratio φ/φ_neutral is at most 1 (brute force, 56 points)", rmax <= 1 + 1e-9)
    fine = maximum(QTLSim.sojourn_ratio(γv, xv) for γv in exp.(range(log(1e-8), log(1e6); length = 141)), xv in range(1e-6, 1 - 1e-6; length = 301))
    check("… and on a 141 × 301 grid of the code's own ratio before clamping would be ≤ 1 (max is 1 at γ → 0)", fine <= 1.0)
    check("sojourn_ratio → 1 as γ → 0, and its x → 1 limit is continuous",
          QTLSim.sojourn_ratio(0.0, 0.3) == 1.0 && isapprox(QTLSim.sojourn_ratio(7.0, 1 - 1e-13), QTLSim.sojourn_ratio(7.0, 1 - 1e-9); rtol = 1e-6))

    base = Params(Ploidy = 2, PopSize = 1000, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01, RecombPoissonMean = 2.0)
    cells = [("Gaussian (μ = 0.01, L = 100)", Params(base; NoQTL = 100, MutationRatePerGamete = 0.01)),
             ("house of cards (μ = 2.5e-5, L = 100)", Params(base; NoQTL = 100, MutationRatePerGamete = 2.5e-5)),
             ("paper-like (4Nμ = 0.4, L = 1000, ω = 5.13)", Params(base; NoQTL = 1000, PhenotypeVar = 0.0348, Heritability = 0.2325,
                                                                     MutationScalar = 0.00268, MutationRatePerGamete = 1e-4, InitFitnessSD = 5.13)),
             ("small (N = 200, L = 20)", Params(base; PopSize = 200, NoQTL = 20, MutationRatePerGamete = 1e-3))]
    # the equation of stationary_spectrum, by brute force in a and x
    function V_rhs(p, st)
        α = sqrt(st.α2)
        f(t) = t^2 * J_brute(2 * st.nparents * (st.α2 * t^2 + st.VG / p.NoQTL) / st.S) * 2 * exp(-t^2 / 2) / sqrt(2π)
        4 * st.nparents * p.NoQTL * p.MutationRatePerGamete * st.α2 * simpson(f, 0.0, 12.0; n = 4000)
    end
    okeq = okn = okgrid = true
    gridratio = Float64[]
    for (name, p) in cells
        st = stationary_spectrum(p)
        rhs = V_rhs(p, st)
        sv = spectrum_vg(p) * (p.PopSize - 1) / p.PopSize      # spectrum_equilibrium reports Σ2a²pq·N/(N−1)
        hog = house_of_gauss_vg(p; nsurv = st.nparents)
        say(@sprintf("    %-44s V %.6g  rhs/V − 1 %+.1e  spectrum_vg/V %.4f  house_of_gauss/V %.3f  sites %.1f", name, st.VG, rhs / st.VG - 1, sv / st.VG, hog / st.VG, QTLSim.stationary_sites(st)))
        okeq &= abs(rhs / st.VG - 1) < 1e-6
        okn &= isapprox(st.nparents, p.PopSize * QTLSim.QTSIM_FITNESS_PEAK * sqrt(initFitSD(p)^2 / st.S); rtol = 1e-12) &&
               isapprox(st.S, initFitSD(p)^2 + st.VG + envVar(p); rtol = 1e-12)
        push!(gridratio, sv / st.VG)
    end
    check("stationary_spectrum(p).VG solves V = 4nLμE[a²J(γ)] (brute-force quadrature, rel. 1e-6), 4 cells", okeq)
    check("its n = N·c₀√(ω²/S) and S = ω² + V + V_E are those of the code's selection", okn)
    check("spectrum_vg(p) (the same model on a frequency grid, implicit Euler) within 1% of it: Gaussian, paper-like, small",
          all(i -> abs(gridratio[i] - 1) < 0.01, (1, 3, 4)))
    # in the house-of-cards cell the grid model's 20 effect classes (Gauss–Hermite nodes) do not
    # resolve a²/(a² + 2c) near a = 0: it is 16% above the closed form by default and converges
    # towards it with more classes (40: 4%; 80 and 160: 1.2%)
    sv80 = spectrum_vg(cells[2][2]; classes = 80) * 999 / 1000
    say(@sprintf("    house of cards: spectrum_vg/V = %.4f with the default 20 effect classes, %.4f with 80", gridratio[2], sv80 / stationary_spectrum(cells[2][2]).VG))
    check("house-of-cards cell: spectrum_vg(p; classes = 80) within 2% of it (the default 20 classes: 16% above)",
          abs(sv80 / stationary_spectrum(cells[2][2]).VG - 1) < 0.02 && 1.1 < gridratio[2] < 1.2)
    pn = Params(base; NoQTL = 50, MutationRatePerGamete = 1e-3, Drift = true)
    stn = stationary_spectrum(pn)
    check("neutral (Drift = true): V = 2N V_m and sites = 4NLμ Σ_{k<2N} 1/k",
          isapprox(stn.VG, 2 * 1000 * muteVar(pn); rtol = 1e-10) && isapprox(QTLSim.stationary_sites(stn), 4 * 1000 * 50 * 1e-3 * sum(1 / k for k in 1:1999); rtol = 1e-12))
    pw = Params(base; NoQTL = 50, MutationRatePerGamete = 1e-3, InitFitnessSD = 3000.0)
    stw = stationary_spectrum(pw)
    check("very weak selection (ω = 3000): V within 0.1% of 2nV_m", isapprox(stw.VG, 2 * stw.nparents * muteVar(pw); rtol = 1e-3))
    ph = Params(base; NoQTL = 10_000, MutationRatePerGamete = 1e-7, MutationScalar = 0.01, InitFitnessSD = 3.0, PopSize = 100_000)
    sth = stationary_spectrum(ph)
    say(@sprintf("    house-of-cards limit: V/(4LμS) = %.4f (γ at a = α: %.3g; c/α² = %.2g)", sth.VG / (4 * 10_000 * 1e-7 * sth.S),
                 2 * sth.nparents * sth.α2 / sth.S, sth.gene / sth.α2))
    check("rare, strongly selected mutations (Nα²/S large, c ≪ α²): V within 2% of the house of cards 4LμS",
          isapprox(sth.VG, 4 * 10_000 * 1e-7 * sth.S; rtol = 0.02))
    # the density integrates to the variance
    let (name, p) = cells[4], st = stationary_spectrum(p), α = sqrt(st.α2)
        inner(av) = simpson(u -> (xv = exp(u); 2 * av^2 * xv * (1 - xv) * QTLSim.stationary_density(st, av, xv) * xv), log(1e-9), 0.0; n = 2000)
        v = 2 * simpson(inner, 0.0, 10α; n = 400)
        check("∫∫ 2a²x(1−x)·stationary_density da dx = VG (rel. 1e-4)", isapprox(v, st.VG; rtol = 1e-4))
    end
    # the sampler's Monte Carlo means, fixed seed: the small cell
    let (name, p) = cells[4], st = stationary_spectrum(p), R = 3000, nsites = QTLSim.stationary_sites(st)
        r = Xoshiro(7)
        genic(G) = sum(QTLSim.var(view(G, l, 1, :) .+ view(G, l, 2, :)) for l in axes(G, 1))
        ns, vs, gs = Float64[], Float64[], Float64[]
        for _ in 1:R
            s = QTLSim.sample_equilibrium(p, r; target = :spectrum, mean = 0.0, spectrum = st)
            push!(ns, s.sites); push!(gs, genic(s.G)); push!(vs, QTLSim.var(vec(sum(s.G, dims = (1, 2)))))
        end
        sem(v) = QTLSim.std(v) / sqrt(length(v))
        say(@sprintf("    sampler, %d draws of the small cell: sites %.3f ± %.3f (expected %.3f); genic %.5f ± %.5f, V_G %.5f ± %.5f (expected %.5f)",
                     R, QTLSim.mean(ns), sem(ns), nsites, QTLSim.mean(gs), sem(gs), QTLSim.mean(vs), sem(vs), st.VG))
        check("sampler: mean number of segregating sites = QTLSim.stationary_sites (within 4 SE and 2%)",
              abs(QTLSim.mean(ns) - nsites) < 4 * sem(ns) && abs(QTLSim.mean(ns) / nsites - 1) < 0.02)
        check("sampler: number of sites is Poisson (variance/mean within 10% of 1)", abs(QTLSim.var(ns) / QTLSim.mean(ns) - 1) < 0.1)
        check("sampler: mean genic variance = VG (within 4 SE and 3%)", abs(QTLSim.mean(gs) - st.VG) < 4 * sem(gs) && abs(QTLSim.mean(gs) / st.VG - 1) < 0.03)
        check("sampler: mean variance of the genotypes = VG (within 4 SE and 3%)", abs(QTLSim.mean(vs) - st.VG) < 4 * sem(vs) && abs(QTLSim.mean(vs) / st.VG - 1) < 0.03)
        ms = [(G = equilibrium_genomes(p, r; target = :spectrum, spectrum = st); sum(G) / p.PopSize) for _ in 1:R]
        check("sampler: the mean genotype is N(0, S/(n(2−b))) (mean within 4 SE of 0, variance within 10%)",
              abs(QTLSim.mean(ms)) < 4 * sem(ms) && abs(QTLSim.var(ms) / QTLSim.stationary_mean_var(p, st) - 1) < 0.1)
        g2 = [genic(equilibrium_genomes(p, r; target = 2 * st.VG, mean = 0.3, spectrum = st)) for _ in 1:R]
        G = equilibrium_genomes(p, r; target = 2 * st.VG, mean = 0.3, spectrum = st)
        check("sampler: target = 2V doubles the genic variance (within 4 SE), mean = 0.3 sets the mean genotype",
              abs(QTLSim.mean(g2) - 2 * st.VG) < 4 * sem(g2) && isapprox(sum(G) / p.PopSize, 0.3; atol = 1e-10))
    end

    finish(LOG, t0)
end

runmain(main, @__FILE__)
