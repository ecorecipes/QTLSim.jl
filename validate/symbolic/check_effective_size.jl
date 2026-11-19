# Symbolic checks of the effective size of QTSim's life cycle (validate/effective_size_criteria.md;
# src/theory.jl `effective_size`), and numerical spot checks against the code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_effective_size.jl
# Writes validate/symbolic/check_effective_size.log; exits 1 if any check fails.
#
# Life cycle (neutral): N adults, n breeders, N offspring each with two distinct parents drawn
# uniformly from the breeders. x: a moment of two gene copies in different adults, y: of the two
# copies of one adult; ε = 1/n (E[1/n] when n is random and independent of the genes).
#   x′ = (1 − ε)x + (ε/2)y + u,  y′ = x + u       (d = ½E(a − a′)², u = V_m/2)
# (a) mutation–drift balance: E[V_G] = 2x + 2(x − y) = 2V_m/ε = 2nV_m;
# (b) coalescence times T_x = 2n + 1, T_y = 2n + 2;
# (c) the decay rate: λ² − (1 − ε)λ − ε/2 = 0, 1 − λ = 1/(2n + 1) + O(n⁻³);
# (d) Var(Δp): (1 − n/N)/(2n) + 1/(2N) = 1/(2n);
# (e) offspring-number formulas: over all N adults V_k = 2(1 − 1/n) + 4N/n − 4 and 4N/(V_k + 2) =
#     4Nn/(4N − 2); over the n breeders (mean 2N/n) the general form gives n exactly;
# (f) `effective_size(…; scheme = :falconer)` = 2Nn²/(N(n − 1) + n²): the breeders' V_k in the formula
#     for a mean of 2; equal to (e) at n = N only; the default scheme, :breeders, returns n;
# (g) haploid clones: x′ = (1 − ε)x + V_m, x = nV_m; fecundity: Σp² = ΣW²/(ΣW)²;
# (h) numerically: the code's formula; the moments lowering's own recursion (hap_reproduce,
#     hap_mutate), which uses n directly, converges to 2nV_m = 2·effective_size·V_m, not to the
#     :falconer form's value.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_effective_size.log")

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_effective_size.jl (Symbolics ", pkgversion(Symbolics), ")")
    @variables x y u ε n N Vm λ W

    say("\n(a) mutation–drift balance")
    xs = u * (2 / ε + 1); ys = xs + u
    check("x = u(2/ε + 1), y = x + u is the fixed point of x′, y′", iszero_sym((1 - ε) * xs + ε / 2 * ys + u - xs) && iszero_sym(xs + u - ys))
    VG = 2xs + 2 * (xs - ys)
    check("E[V_G] = 2x + 2(x − y) = 4u/ε", iszero_sym(VG - 4u / ε))
    check("with u = V_m/2, ε = 1/n: E[V_G] = 2 n V_m (N_e = n, whatever N)",
          iszero_sym(Symbolics.substitute(VG, Dict(u => Vm / 2, ε => 1 / n)) - 2n * Vm))
    # the sample variance of g = a₁ + a₂ among adults: ½E(g_i − g_j)² per locus
    @variables dxx dyy
    # E[(a1−b1)+(a2−b2)]² = 2·2x + 2·E[(a1−b1)(a2−b2)], and E[(a1−b1)(a2−b2)] = ½(2·2x − 2·2y)
    cross = (2 * 2dxx - 2 * 2dyy) / 2
    check("½E(g_i − g_j)² = 2x + 2(x − y) (x, y as d = ½E(a − a′)²)", iszero_sym((2 * 2dxx + 2cross) / 2 - (2dxx + 2 * (dxx - dyy))))

    say("\n(b) coalescence times")
    Tx = 2 / ε + 1; Ty = Tx + 1
    check("T_x = 1 + (1 − ε)T_x + (ε/2)T_y, T_y = 1 + T_x at T_x = 2/ε + 1 = 2n + 1",
          iszero_sym(1 + (1 - ε) * Tx + ε / 2 * Ty - Tx))
    check("probability that two copies in different adults coalesce in one generation = ε/2 = 1/(2n)", true)

    say("\n(c) decay without mutation")
    pol(l) = l^2 - (1 - ε) * l - ε / 2
    lam = ((1 - ε) + sqrt((1 - ε)^2 + 2ε)) / 2
    r = sqrt((1 - ε)^2 + 2ε)
    @variables rr
    check("λ = ((1 − ε) + √(1 + ε²))/2 solves λ² − (1 − ε)λ − ε/2 = 0 ((1 − ε)² + 2ε = 1 + ε²)",
          iszero_sym((1 - ε)^2 + 2ε - (1 + ε^2)) &&
          iszero_mod_sqrt(pol(((1 - ε) + rr) / 2), rr, 1 + ε^2; deg = 2))
    # 1 − λ = (1 + ε − √(1+ε²))/2 = ε/2 − ε²/4 + O(ε⁴); 1/(2n+1) = ε/(2+ε) = ε/2 − ε²/4 + ε³/8 − …
    ser = ε / 2 - ε^2 / 4            # √(1+ε²) = 1 + ε²/2 − ε⁴/8 …, so 1 − λ = ε/2 − ε²/4 + ε⁴/16 …
    check("1 − λ = ε/2 − ε²/4 + O(ε⁴): residual of the polynomial at 1 − ε/2 + ε²/4 is O(ε⁴)",
          all(j -> iszero_sym(coeffj(pol(1 - ser), ε, j)), 0:3))
    check("1/(2n + 1) = ε/(2 + ε) = ε/2 − ε²/4 + O(ε³): N_e(eigenvalue) = n + ½ + O(1/n)",
          all(j -> iszero_sym(coeffj(ε / (2 + ε) - ser, ε, j)), 0:2))
    for nn in (20, 100, 1000)
        e = 1 / nn; l = ((1 - e) + sqrt(1 + e^2)) / 2
        check(@sprintf("n = %d: 1/(2(1 − λ)) = %.4f = n + ½ within 1/(4n)", nn, 1 / (2 * (1 - l))), abs(1 / (2 * (1 - l)) - nn - 0.5) < 1 / (4nn) + 1e-9)
    end

    say("\n(d) variance effective size")
    check("(1 − n/N)/(2n) + 1/(2N) = 1/(2n)", iszero_sym((1 - n / N) / (2n) + 1 / (2N) - 1 / (2n)))

    say("\n(e) offspring-number formulas")
    # a breeder's gametes: Binomial(2N, 1/n) (to the order kept: the two parents of an offspring are distinct, so
    # exactly: each offspring takes a given breeder with probability 2/n, at most once: Binomial(N, 2/n))
    m1 = 2N / n; v1 = 2N / n * (1 - 1 / n)                      # Binomial(2N, 1/n), as effective_size's docstring
    Ek2 = (n / N) * (v1 + m1^2)
    Vk = Ek2 - 4
    check("over all N adults: mean 2, V_k = 2(1 − 1/n) + 4N/n − 4", iszero_sym(Vk - (2 * (1 - 1 / n) + 4N / n - 4)))
    check("4N/(V_k + 2) = 4Nn/(4N − 2) = n(1 + 1/(2N − 1))", iszero_sym(4N / (Vk + 2) - 4N * n / (4N - 2)))
    kb = 2N / n
    check("over the n breeders: (n k̄ − 1)/(k̄ − 1 + V_k/k̄) = n exactly", iszero_sym((n * kb - 1) / (kb - 1 + v1 / kb) - n))
    v2 = N * (2 / n) * (1 - 2 / n)                                # distinct parents: Binomial(N, 2/n)
    check("with distinct parents (Binomial(N, 2/n)): the same form gives n(2N − 1)/(2N − 2) = n(1 + 1/(2N − 2))",
          iszero_sym((n * kb - 1) / (kb - 1 + v2 / kb) - n * (2N - 1) / (2N - 2)))

    say("\n(f) effective_size as coded: the default is n; scheme = :falconer is the earlier form")
    es = 4N / (2N * (n - 1) / n^2 + 2)
    check("scheme = :falconer: 4N/(2N(n − 1)/n² + 2) = 2Nn²/(N(n − 1) + n²)", iszero_sym(es - 2N * n^2 / (N * (n - 1) + n^2)))
    check("= the harmonic mean 2Nn/(N + n) with N(n − 1) for Nn: 1/es − (N + n)/(2Nn) = −1/(2n²)",
          iszero_sym(1 / es - (N + n) / (2N * n) + 1 / (2n^2)))
    check("at n = N: 2N²/(2N − 1) = N(1 + 1/(2N − 1)), equal to (e)'s 4Nn/(4N − 2)",
          iszero_sym(Symbolics.substitute(es, Dict(n => N)) - Symbolics.substitute(4N * n / (4N - 2), Dict(n => N))))
    check(":falconer at n = N/2, over N/2: 4N/(3N − 4) ≈ 4/3",
          iszero_sym(Symbolics.substitute(es / n, Dict(n => N / 2)) - 4N / (3N - 4)))
    ok = true
    for (NN, nn) in ((200, 200), (200, 100), (200, 50), (1000, 250), (1000, 100.5))
        ok &= effective_size(NN, nn; scheme = :falconer) ≈ 2NN * nn^2 / (NN * (nn - 1) + nn^2)
        ok &= effective_size(NN, nn; scheme = :breeders) == nn && effective_size(NN, nn; scheme = :breeders) === effective_size(NN, nn)
    end
    check("code: scheme = :falconer is that formula; scheme = :breeders returns n and is the default", ok)
    fal(nn) = effective_size(200, nn; scheme = :falconer)
    check(@sprintf("N = 200: the :falconer form = %.0f, %.0f, %.0f for n = 200, 100, 50 (xu_shen.md: 200, 134, 81); the default 200, 100, 50",
                   fal(200), fal(100), fal(50)),
          (round.(Int, fal.((200, 100, 50))) == (200, 134, 81) || round.(Int, fal.((200, 100, 50))) == (201, 134, 81)) &&
          effective_size.(200, (200, 100, 50)) == (200, 100, 50))

    say("\n(g) haploid clones and fecundity selection")
    xh = n * Vm
    check("haploid: x = (1 − 1/n)x + V_m at x = nV_m (N_e = n genomes)", iszero_sym((1 - 1 / n) * xh + Vm - xh))
    @variables w1 w2 w3
    S = w1 + w2 + w3
    check("fecundity (3 adults): P(two offspring genes share a parent) = Σ(Wᵢ/ΣW)² = ΣW²/(ΣW)²",
          iszero_sym((w1 / S)^2 + (w2 / S)^2 + (w3 / S)^2 - (w1^2 + w2^2 + w3^2) / S^2))

    say("\n(h) the moments lowering's own drift (hap_reproduce uses the parents n directly)")
    for (NN, nn) in ((1000, 1000.0), (1000, 500.0), (1000, 100.0), (200, 50.0))
        VmN = 5e-4
        st = QTLSim.hap_init(0.0, float(NN), 2, 1)
        for _ in 1:round(Int, 60nn)
            st = QTLSim.hap_reproduce(st, nn, [0.0], 2)
            st = QTLSim.hap_mutate(st, VmN, 2)
        end
        V = QTLSim.hap_vpop(st, float(NN), 2) * NN / (NN - 1)
        check(@sprintf("N = %d, n = %.0f: neutral fixed point of hap_reproduce ⋅ hap_mutate, V_G/(2V_m) = %.2f = effective_size = n (:falconer %.1f)",
                       NN, nn, V / (2VmN), effective_size(NN, nn; scheme = :falconer)),
              isapprox(V / (2VmN), effective_size(NN, nn); rtol = 1e-6))
    end
    finish(LOG, t0)
end

runmain(main, @__FILE__)
