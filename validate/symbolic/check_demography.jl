# Symbolic checks of the demography theory (src/theory.jl, docs/src/parameters.md), and
# numerical spot checks against the code. (The fluctuating optimum's periodic_response
# and red_noise_lag_var are in check_traits.jl, with their matrix forms.)
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_demography.jl
# Writes validate/symbolic/check_demography.log; exits 1 if any check fails.
#
# (a) Beverton–Holt: with n = N W̄ survivors and R n/(1 + (R − 1)n/K) offspring, the adults
#     settle at K(RW̄ − 1)/((R − 1)W̄), the survivors at K(RW̄ − 1)/(R − 1), K when W̄ = 1,
#     and the fixed point is stable (slope 1/(RW̄)) when RW̄ > 1;
# (b) the Galton–Watson process of extinction_distribution(:branching) and
#     rescue_probability: the thinned Poisson pgf φ(s) = 1 − W + W e^{r(s − 1)} has φ(1) = 1,
#     mean rW and variance rW(1 + r − rW), and its extinction probability is the least
#     root of φ(q) = q, below 1 iff rW > 1; numerically against rescue_probability and the
#     branching extinction probabilities.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_demography.log")
const rng = Xoshiro(20261112)
runif(a, b) = runif(rng, a, b)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_demography.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) Beverton–Holt")
    @variables R W K N
    bh(n) = R * n / (1 + (R - 1) * n / K)
    Nstar = K * (R * W - 1) / ((R - 1) * W)
    check("adults: N* = K(RW̄ − 1)/((R − 1)W̄) solves N = R N W̄/(1 + (R − 1)N W̄/K)", iszero_sym(bh(Nstar * W) - Nstar))
    check("survivors: n* = N*W̄ = K(RW̄ − 1)/(R − 1)", iszero_sym(Nstar * W - K * (R * W - 1) / (R - 1)))
    check("W̄ = 1: N* = K", iszero_sym(Symbolics.substitute(Nstar, Dict(W => 1)) - K))
    slope = Symbolics.expand_derivatives(Differential(N)(bh(N * W)))
    check("slope of the map at N*: 1/(RW̄) (stable iff RW̄ > 1)", iszero_sym(Symbolics.substitute(slope, Dict(N => Nstar)) - 1 / (R * W)))
    for trial in 1:3
        R0, W0, K0 = runif(1.2, 5.0), runif(0.5, 1.0), rand(rng, 100:5000)
        R0 * W0 > 1.05 || continue
        p = Params(PopSize = K0, Demography = :beverton_holt, GrowthRate = R0)
        Nn = float(K0)
        for _ in 1:20_000
            Nn = QTLSim.expected_offspring(p, Nn * W0)
        end
        check(@sprintf("expected_offspring iterated, R=%.2f W̄=%.2f K=%d: → K(RW̄ − 1)/((R − 1)W̄)", R0, W0, K0),
              isapprox(Nn, K0 * (R0 * W0 - 1) / ((R0 - 1) * W0); rtol = 1e-9))
    end

    # ===================================================================================
    say("\n(b) the Galton–Watson process")
    @variables s r w
    φ = 1 - w + w * exp(r * (s - 1))
    d1 = Symbolics.expand_derivatives(Differential(s)(φ))
    d2 = Symbolics.expand_derivatives(Differential(s)(d1))
    at1(e) = Symbolics.substitute(e, Dict(s => 1))
    # exp(0) is not folded by the simplifier: the derivatives are symbolic, evaluated at random (r, W)
    okφ = all(1:5) do _
        r0, w0 = runif(0.5, 5.0), runif(0.1, 1.0)
        ev(e) = Symbolics.build_function(e, s, r, w; expression = Val(false))(1.0, r0, w0)
        isapprox(ev(φ), 1; rtol = 1e-14) && isapprox(ev(d1), r0 * w0; rtol = 1e-14) && isapprox(ev(d2), r0^2 * w0; rtol = 1e-14)
    end
    check("φ(1) = 1, φ′(1) = rW (the mean), φ″(1) = r²W (variance rW(1 + r − rW))", okφ)
    for trial in 1:3
        R0, ω0, VG0, VE0 = runif(1.5, 4.0), runif(1.0, 5.0), runif(0.1, 1.0), runif(0.5, 2.0)
        N0 = rand(rng, 1:20)
        res = rescue_probability(; N0, R = R0, d0 = 0.0, VG = VG0, VE = VE0, ω = ω0)
        Winf = ω0 / sqrt(ω0^2 + VG0 + VE0)
        f(q) = 1 - Winf + Winf * exp(R0 * (q - 1)) - q
        if R0 * Winf <= 1
            check("rescue_probability = 0 when RW̄ ≤ 1 (subcritical)", res.rescue == 0)
            continue
        end
        lo, hi = 0.0, 1 - 1e-12                         # f(0) > 0, f(1⁻) < 0 for RW̄ > 1
        for _ in 1:200
            mid = (lo + hi) / 2
            f(mid) > 0 ? (lo = mid) : (hi = mid)
        end
        check(@sprintf("rescue_probability, no shift (RW̄=%.2f): 1 − q^N₀, q the least root of φ(q) = q", R0 * Winf),
              isapprox(res.rescue, 1 - lo^N0; rtol = 1e-9, atol = 1e-12))
        # extinction_distribution(:branching) with a constant W̄ and K so large that r_t = R: q^N₀ in the limit
        p = Params(PopSize = 10^12, Demography = :beverton_holt, GrowthRate = R0)
        d = extinction_distribution(p, fill(Winf, 60); N0 = N0, method = :branching)
        check("extinction_distribution(:branching), constant W̄, K → ∞: P(extinct) → q^N₀",
              isapprox(d.pext[end], lo^N0; rtol = 1e-6, atol = 1e-10))
    end

    finish(LOG, t0)
end

runmain(main, @__FILE__)
