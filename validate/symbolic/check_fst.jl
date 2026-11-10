# Symbolic checks of neutral differentiation among demes (src/theory.jl, src/migration.jl):
# Q_ST in the island model and the within-deme variance, and numerical spot checks
# against the code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_fst.jl
# Writes validate/symbolic/check_fst.log; exits 1 if any check fails.
#
# (a) the island model reduces deme_coalescence's recursion to three variables (W, the
#     two genes of one adult; S₀, two adults of one deme; S₁, of two demes); its fixed
#     point gives V_W = 2NkV_m whatever m, σ²_B = 2kV_m λ/(1 − λ) and
#     Q_ST = λ/(λ + 2N(1 − λ)) (:demes; ×(k − 1)/k on the ratio for :total) with
#     λ = (1 − mk/(k − 1))², the difference of the same-origin probabilities of two adults
#     of one deme and of two demes; Wright's k = ∞ form and Takahata's small-m forms;
# (b) the kernel's exchanges (exchange_origins): the same recursion with the two-deme
#     same-origin probability lowered by 2m/((k − 1)N) gives V_W = 2NV_m(kN − 1)/(N − 1),
#     and two distinct adults share a deme with probability (N − 1)/(kN − 1);
# (c) Slatkin's stepping-stone F_ST from t₀ and t_i;
# (d) numerical spot checks of qst_island, island_variances, deme_coalescence,
#     neutral_deme_variances, exchange_origins and backward_migration.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_fst.log")
const rng = Xoshiro(20261111)
runif(a, b) = runif(rng, a, b)

# the reduced recursion: residuals of S₀ and S₁ (W = 1 + S₀) for same-origin probabilities π₀, π₁
function reduced_residuals(π0, π1, N, S0, S1)
    W = 1 + S0
    O_same = 1 + W / (2N) + (1 - 1 / N) * S0
    [S0 - (π0 * O_same + (1 - π0) * (1 + S1)), S1 - (π1 * O_same + (1 - π1) * (1 + S1))]
end

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_fst.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the island model: the reduced coalescence recursion")
    # first in the abstract: any same-origin probabilities π₀ (one deme) and π₁ (two demes)
    @variables P0 P1 N m k
    s0 = 1 + 2N * (1 - P0 + P1) / P1
    s1 = s0 + (P0 - P1) / P1
    w = 1 + s0
    check("fixed point: S₀ = 1 + 2N(1 − π₀ + π₁)/π₁, S₁ = S₀ + (π₀ − π₁)/π₁, W = 1 + S₀",
          iszero_sym(reduced_residuals(P0, P1, N, s0, s1)))
    within = 2s0 - w                                            # /V_m
    half = (1 - 1 / N) * s0 + w / (2N)
    pairs = 2s1 - 2half
    between = pairs - within / N
    check("so V_W = V_m(S₀ − 1) = 2NV_m(1 − π₀ + π₁)/π₁ and σ²_B = 2V_m(π₀ − π₁)/π₁",
          iszero_sym(within - 2N * (1 - P0 + P1) / P1) && iszero_sym(between - 2(P0 - P1) / P1))
    # the island model
    π0 = (1 - m)^2 + m^2 / (k - 1)
    π1 = 2(1 - m) * m / (k - 1) + (k - 2) * m^2 / (k - 1)^2
    λ = (1 - m * k / (k - 1))^2
    check("island: π₀ − π₁ = (1 − mk/(k − 1))² = λ", iszero_sym(π0 - π1 - λ))
    check("island: π₁ = (1 − λ)/k", iszero_sym(π1 - (1 - λ) / k))
    # hence V_W = 2NkV_m and σ²_B = 2kV_m λ/(1 − λ); written in λ and k
    @variables lam
    Wk = 2N * (1 - lam)/((1 - lam) / k)
    Bk = 2lam / ((1 - lam) / k)
    check("V_W = 2NkV_m, whatever m", iszero_sym(Wk - 2N * k))
    check("σ²_B = 2kV_m λ/(1 − λ)", iszero_sym(Bk - 2k * lam / (1 - lam)))
    R = Bk / (2Wk)
    check("σ²_B/(2V_W) = λ/(2N(1 − λ)), so Q_ST (:demes) = λ/(λ + 2N(1 − λ))",
          iszero_sym(R - lam / (2N * (1 - lam))) && iszero_sym(R / (R + 1) - lam / (lam + 2N * (1 - lam))))
    Rt = (k - 1) / k * R
    check("Q_ST (:total) = (k−1)/k·R/((k−1)/k·R + 1) = (k−1)λ/((k−1)λ + 2Nk(1 − λ))",
          iszero_sym(Rt / (Rt + 1) - (k - 1) * lam / ((k - 1) * lam + 2N * k * (1 - lam))))
    λ∞ = (1 - m)^2
    check("k = ∞: λ/(λ + 2N(1 − λ)) = (1 − m)²/(2N − (2N − 1)(1 − m)²) (Wright)",
          iszero_sym(λ∞ / (λ∞ + 2N * (1 - λ∞)) - λ∞ / (2N - (2N - 1) * λ∞)))
    # small m: 1/Q_ST − 1 = 2N(1 − λ)/λ; 1 − λ = 2mk/(k − 1) − m²k²/(k − 1)², λ = 1 + O(m)
    oneminus = 1 - λ
    check("small m: 1 − λ = 2mk/(k−1) + O(m²), so 1/Q_ST − 1 = 4Nmk/(k−1) + O(m²) (:demes, Takahata)",
          iszero_sym(coeffj(oneminus, m, 0)) && iszero_sym(coeffj(oneminus, m, 1) - 2k / (k - 1)))
    check("… and ×k/(k − 1): 4Nm(k/(k − 1))² + O(m²) (:total)",
          iszero_sym(2N * coeffj(oneminus, m, 1) * k / (k - 1) - 4N * (k / (k - 1))^2))
    check("panmixia m = (k − 1)/k: λ = 0, Q_ST = 0", iszero_sym(Symbolics.substitute(λ, Dict(m => (k - 1) / k))))

    # ===================================================================================
    say("\n(b) the kernel's exchanges")
    # exchange_origins lowers the two-deme same-origin probability by 2δ, δ = P_ab/N = m/((k − 1)N).
    # With u = m/(k − 1): π₁ = u(2 − ku), 1 − π₀ + π₁ = kπ₁, 2δ = 2u/N, so by the abstract fixed point
    # V_W = 2NV_m (kπ₁ − 2δ)/(π₁ − 2δ) = 2NV_m (k(2 − ku) − 2/N)/(2 − ku − 2/N)
    @variables u
    check("island: π₁ = u(2 − ku), u = m/(k − 1)", iszero_sym(Symbolics.substitute(π1, Dict(m => u * (k - 1))) - u * (2 - k * u)))
    Vx = 2N * (k * (2 - k * u) - 2 / N) / (2 - k * u - 2 / N)
    check("exchanges: V_W → 2NV_m(kN − 1)/(N − 1) as m → 0", iszero_sym(coeffj(Vx, u, 0) - 2N * (k * N - 1) / (N - 1)))
    check("… only to leading order: ∂V_W/∂u = 2NV_m·2k(k − 1)/(N(2 − 2/N)²) at m = 0 (relative O(m/N))",
          iszero_sym(coeffj(Vx, u, 1) - 2N * 2k * (k - 1) / (N * (2 - 2 / N)^2)))
    check("two distinct adults of the metapopulation share a deme w.p. kN(N−1)/(kN(kN−1)) = (N−1)/(kN−1)",
          iszero_sym(k * N * (N - 1) / (k * N * (k * N - 1)) - (N - 1) / (k * N - 1)))

    # ===================================================================================
    say("\n(c) Slatkin's stepping stone")
    @variables i
    t0s = 2N * k
    ti = 2N * k + (k - i) * i / (2m)
    check("F_ST = (t_i − t₀)/(t_i + t₀) = 1/(1 + 8Nmk/((k − i)i)) (Slatkin 1991, eqs 16–17)",
          iszero_sym((ti - t0s) / (ti + t0s) - 1 / (1 + 8N * m * k / ((k - i) * i))))

    # ===================================================================================
    say("\n(d) numerical spot checks against the code")
    for (k0, N0, rate) in ((2, 20, 0.01), (3, 50, 0.002), (5, 30, 0.03), (8, 12, 0.1))
        P = backward_migration(k0, rate, N0)
        m0 = 1 - P[1, 1]
        π0n = (1 - m0)^2 + m0^2 / (k0 - 1); π1n = 2(1 - m0) * m0 / (k0 - 1) + (k0 - 2) * m0^2 / (k0 - 1)^2
        dd = Dict(P0 => π0n, P1 => π1n, N => N0)
        c = deme_coalescence(P, N0)
        check(@sprintf("deme_coalescence, k=%d N=%d m=%.4f: W, S₀, S₁ = the symbolic fixed point", k0, N0, m0),
              isapprox(c.W[1], num(w, dd); rtol = 1e-9) && isapprox(c.S[1, 1], num(s0, dd); rtol = 1e-9) &&
              isapprox(c.S[1, 2], num(s1, dd); rtol = 1e-9) && all(x -> isapprox(x, c.S[1, 2]; rtol = 1e-9), (c.S[a, b] for a in 1:k0, b in 1:k0 if a != b)))
        Vm0 = 10^runif(-4, -2)
        v = neutral_deme_variances(P, N0, Vm0)
        iv = island_variances(Vm0, N0, m0, k0)
        check("neutral_deme_variances = island_variances: within 2NkV_m, between 2kV_m λ/(1 − λ)",
              all(x -> isapprox(x, iv.within; rtol = 1e-9), v.within) && isapprox(v.between, iv.between; rtol = 1e-8) &&
              isapprox(iv.within, 2N0 * k0 * Vm0; rtol = 1e-12))
        check("qst_island (:total, :demes) = neutral_deme_variances' qst, qst_demes",
              isapprox(qst_island(N0, m0; k = k0), v.qst; rtol = 1e-8) && isapprox(qst_island(N0, m0; k = k0, among = :demes), v.qst_demes; rtol = 1e-8))
        vx = neutral_deme_variances(P, N0, Vm0; origins = :exchange)
        u0 = m0 / (k0 - 1)
        vxe = 2N0 * Vm0 * (k0 * (2 - k0 * u0) - 2 / N0) / (2 - k0 * u0 - 2 / N0)
        check(@sprintf("origins = :exchange: within = 2NV_m(k(2 − ku) − 2/N)/(2 − ku − 2/N) (÷ (kN−1)/(N−1) form: 1 + %.1e)",
                       vxe / (2N0 * Vm0 * (k0 * N0 - 1) / (N0 - 1)) - 1),
              all(x -> isapprox(x, vxe; rtol = 1e-9), vx.within) && isapprox(vxe, 2N0 * Vm0 * (k0 * N0 - 1) / (N0 - 1); rtol = 2m0 / N0))
        J = exchange_origins(P, N0)
        ix(a, b) = a + (b - 1) * k0
        okm = all(isapprox(sum(J[ix(a, b), ix(c1, d1)] for d1 in 1:k0), P[a, c1]; atol = 1e-14) &&
                  isapprox(sum(J[ix(a, b), ix(d1, c1)] for d1 in 1:k0), P[b, c1]; atol = 1e-14) for a in 1:k0, b in 1:k0, c1 in 1:k0 if a != b)
        δ = P[1, 2] / N0
        check("exchange_origins: rows sum to 1, marginals P, (a,b)→(a,a) is P_aa P_ba − P_ab/N",
              all(r -> isapprox(sum(r), 1; atol = 1e-13), eachrow(J)) && okm &&
              isapprox(J[ix(1, 2), ix(1, 1)], P[1, 1] * P[2, 1] - δ; rtol = 1e-12) &&
              isapprox(J[ix(1, 2), ix(2, 1)], P[1, 2] * P[2, 1] + δ; rtol = 1e-12))
        if k0 <= 5
            Pall = QTLSim.migration_mixing(island_migration(k0, rate), fill(float(N0), k0), trues(k0))
            check("backward_migration (island, one order) = migration_mixing over all k! orders",
                  isapprox(P, Pall; rtol = 1e-12) && all(r -> isapprox(sum(r), 1; atol = 1e-13), eachrow(P)))
        end
        check(@sprintf("backward m ≈ 2·rate for small rates (m/(2·rate) = %.4f)", m0 / (2rate)), isapprox(m0, 2rate; rtol = 2rate + 2 / N0))
        cT = deme_coalescence(P, N0; t = [1, 20000])
        check("deme_coalescence(t): t = 1 gives W = S = 1; large t approaches the equilibrium",
              all(x -> isapprox(x, 1; atol = 1e-12), [cT[1].W; vec(cT[1].S)]) && isapprox(cT[2].S[1, 2], c.S[1, 2]; rtol = 1e-3))
    end
    check("qst_island, k = ∞: (1 − m)²/(2N − (2N − 1)(1 − m)²)",
          isapprox(qst_island(40, 0.02), 0.98^2 / (80 - 79 * 0.98^2); rtol = 1e-12))

    finish(LOG, t0)
end

runmain(main, @__FILE__)
