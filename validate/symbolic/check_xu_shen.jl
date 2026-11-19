# Symbolic checks of the algebra of Xu and Shen (2026, bioRxiv preprint,
# doi:10.64898/2026.08.05.743160, version 1), Supplement S1–S2, as used by src/xu_shen.jl,
# and numerical spot checks against the code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_xu_shen.jl
# Writes validate/symbolic/check_xu_shen.log; exits 1 if any check fails.
#
# w stands for ω², M for V_M, n for n_e, c = n w.
# (a) the equilibrium: S2 (eq. 7) and S4 (eq. 8, the argument of eq. 6) from S1 = 0;
# (b) S5, S6a, S6b;
# (c) the Jacobian S8 from the derivatives of S1 at the equilibrium;
# (d) the determinant: S10–S15;
# (e) the third Routh–Hurwitz condition: S16–S21, and a decomposition of a₁a₂ − a₃ into
#     positive terms;
# (f) the invariances used by the tests;
# (g) numerically: the Julia functions against these expressions.
#
# Where a printed equation differs from what the derivation gives, the check is of the
# derived form and the difference is recorded as a NOTE (not a failure), with whether it
# affects a result. The notes are collected at the end of the log.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_xu_shen.log")
const rng = Xoshiro(20261118)
const notes = String[]
note(s) = (push!(notes, s); say("  NOTE: ", s))

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_xu_shen.jl (Symbolics ", pkgversion(Symbolics), ")")
    @variables Δ V N v r w M K n
    c = n * w
    D(x) = e -> Symbolics.expand_derivatives(Differential(x)(e))
    FΔ = v - V * Δ / w
    FV = M - V^2 / c - V / N
    FN = N * (r - (V + Δ^2) / (2w) - N / K)
    N1 = V * c / (M * c - V^2)
    f = r - V / (2w) - N1 / K
    g = 2V^2 / w
    h = g * f
    Δs = v * w / V
    # r from S4 (v² = h(V)): the equilibrium condition solved for r, so that identities at an
    # equilibrium are rational in (v, V, w, M, K, n)
    rS4 = v^2 * w / (2V^2) + V / (2w) + N1 / K
    ateq(e) = Symbolics.substitute(e, Dict(Δ => Δs, N => N1, r => rS4))

    # ===================================================================================
    say("\n(a) the equilibrium (S2, S4; eqs. 7, 8)")
    check("dΔ/dt = 0 at Δ* = vω²/V (S2, eq. 7)", iszero_sym(Symbolics.substitute(FΔ, Dict(Δ => Δs))))
    check("dV/dt = 0 at N* = V n_e ω²/(V_M n_e ω² − V²) (S2, eq. 7)", iszero_sym(Symbolics.substitute(FV, Dict(N => N1))))
    N2 = K * (r - (V + Δs^2) / (2w))
    check("dN/dt = 0 (N > 0) ⇔ N = N₂(V) = K(r − (V + Δ*²)/(2ω²)) (S3, eq. 8)",
          iszero_sym(Symbolics.substitute(FN / N, Dict(Δ => Δs, N => N2))))
    check("N₁ = N₂ ⇔ v² = (2V²/ω²)[r − V/(2ω²) − N₁/K] (S4; the argument of eq. 6)",
          iszero_sym((N2 - N1) * 2V^2 / (K * w) - (h - v^2)))
    check("all three right-hand sides vanish at (Δ*, V, N₁) with r from S4", iszero_sym([ateq(FΔ), ateq(FV), ateq(FN)]))

    # ===================================================================================
    say("\n(b) S5, S6a, S6b")
    dN1 = D(V)(N1)
    ddN1 = D(V)(dN1)
    check("S6a: N₁′ = n_e ω²(V_M n_e ω² + V²)/(V_M n_e ω² − V²)²", iszero_sym(dN1 - c * (M * c + V^2) / (M * c - V^2)^2))
    check("S6b: N₁″ = 2V n_e ω²(3V_M n_e ω² + V²)/(V_M n_e ω² − V²)³", iszero_sym(ddN1 - 2V * c * (3M * c + V^2) / (M * c - V^2)^3))
    dh = D(V)(h)
    φ = 2f - V * dN1 / K - V / (2w)
    check("S5: h′ = (2V/ω²)[2f − V N₁′/K − V/(2ω²)]", iszero_sym(dh - 2V / w * φ))
    check("g′ = 4V/ω² and f′ = −1/(2ω²) − N₁′/K", iszero_sym(D(V)(g) - 4V / w) && iszero_sym(D(V)(f) + 1 / (2w) + dN1 / K))
    dφ = D(V)(φ)
    check("φ′ = −3/(2ω²) − 3N₁′/K − V N₁″/K (each term negative on (0, V_c^M): φ decreases)",
          iszero_sym(dφ + 3 / (2w) + 3dN1 / K + V * ddN1 / K))
    check("φ(0) = 2r", iszero_sym(Symbolics.substitute(φ, Dict(V => 0)) - 2r))
    note("text after S5 (line 394): the decreasing function is written with the term ½·V/(2ω²); S5 itself has V/(2ω²) " *
         "(checked above). Either is decreasing, so the argument and eq. 6 are unaffected.")
    note("S1 (line 377) describes extinction as N → 0, V/N → V_M, Δ → +∞; this is an attracting state rather than an " *
         "equilibrium of S1 (Δ grows without bound), which is how src/xu_shen.jl treats it (a threshold on N).")

    # ===================================================================================
    say("\n(c) the Jacobian (S7, S8)")
    Fs, ys = [FΔ, FV, FN], [Δ, V, N]
    Jgen = [D(y)(F) for F in Fs, y in ys]
    Jcode = [-V/w -Δ/w 0; 0 (-2V / c-1 / N) V/N^2; -N*Δ/w -N/(2w) (r - (V + Δ^2) / (2w)-2N / K)]
    check("the general Jacobian of S1 is xu_shen_jacobian's matrix", iszero_sym(Jgen .- Jcode))
    S8 = [-V/w -v/V 0; 0 (-2V / c-1 / N1) V/N1^2; -N1*v/V -N1/(2w) -N1/K]
    check("at the equilibrium it is S8 (J₁₂ = −v/V, J₃₁ = −N*v/V, J₃₃ = −N*/K)", iszero_sym(ateq.(Jgen) .- S8))
    note("line 410: \"All other terms (J13, J21, J33) are zero\"; J₃₃ = −N*/K is given two lines above and in S8, " *
         "so the list should read (J₁₃, J₂₁). No effect.")

    # ===================================================================================
    say("\n(d) the determinant (S10–S15)")
    @variables J11 J12 J22 J23 J31 J32 J33
    Jsym = [J11 J12 0; 0 J22 J23; J31 J32 J33]
    det3(A) = A[1, 1] * (A[2, 2] * A[3, 3] - A[2, 3] * A[3, 2]) - A[1, 2] * (A[2, 1] * A[3, 3] - A[2, 3] * A[3, 1]) +
              A[1, 3] * (A[2, 1] * A[3, 2] - A[2, 2] * A[3, 1])
    X = J11 * J32 - J31 * J12
    check("S10: det J = J₁₁J₂₂J₃₃ − J₂₃(J₁₁J₃₂ − J₃₁J₁₂)", iszero_sym(det3(Jsym) - (J11 * J22 * J33 - J23 * X)))
    sub8 = Dict(J11 => S8[1, 1], J12 => S8[1, 2], J22 => S8[2, 2], J23 => S8[2, 3], J31 => S8[3, 1], J32 => S8[3, 2], J33 => S8[3, 3])
    at8(e) = Symbolics.substitute(e, sub8)
    check("S11: J₁₁J₃₂ − J₃₁J₁₂ = V*N*/(2ω⁴) − v²N*/V*²", iszero_sym(at8(X) - (V * N1 / (2w^2) - v^2 * N1 / V^2)))
    check("S12: −J₂₃(J₁₁J₃₂ − J₃₁J₁₂) = [v² − V*³/(2ω⁴)]/(N*V*)", iszero_sym(at8(-J23 * X) - (v^2 - V^3 / (2w^2)) / (N1 * V)))
    check("S13: J₁₁J₂₂J₃₃ = −(N*/(Kω²))[V_M + V*²/(n_e ω²)]", iszero_sym(at8(J11 * J22 * J33) + N1 / (K * w) * (M + V^2 / c)))
    dheq = Symbolics.substitute(dh, Dict(r => rS4))
    check("S14, first line: h′(V*) = 2v²/V* − V*²/ω⁴ − 2V*²N₁′/(Kω²)", iszero_sym(dheq - (2v^2 / V - V^2 / w^2 - 2V^2 * dN1 / (K * w))))
    check("S14, second line: = (2/V*)[v² − V*³/(2ω⁴)] − (2N*²/(Kω²))[V_M + V*²/(n_e ω²)]",
          iszero_sym(dheq - (2 / V * (v^2 - V^3 / (2w^2)) - 2N1^2 / (K * w) * (M + V^2 / c))))
    check("S15: det J = h′(V*)/(2N*)", iszero_sym(det3(S8) - dheq / (2N1)))

    # ===================================================================================
    say("\n(e) the third Routh–Hurwitz condition (S16–S21)")
    a1 = -(J11 + J22 + J33)
    a2 = J11 * J22 + J22 * J33 - J23 * J32 + J11 * J33
    a3 = -J11 * J22 * J33 + J23 * X
    @variables λ
    charpoly = det3(Jsym .- λ .* [1 0 0; 0 1 0; 0 0 1])
    check("S16: det(J − λI) = −(λ³ + a₁λ² + a₂λ + a₃) with the printed a₁, a₂, a₃", iszero_sym(charpoly + λ^3 + a1 * λ^2 + a2 * λ + a3))
    P1print = J23 * X - J33 * (J11 + J22)^2
    P2 = -(J11 + J22) * (J11 * J22 - J23 * J32) + J33 * J23 * J32 - J33^2 * (J11 + J22)
    P1der = -J23 * X - J33 * (J11 + J22)^2
    check("a₁a₂ − a₃ = P₁ + P₂ with P₁ = −J₂₃(J₁₁J₃₂ − J₃₁J₁₂) − J₃₃(J₁₁ + J₂₂)² and P₂ as in S18", iszero_sym(a1 * a2 - a3 - (P1der + P2)))
    check("with P₁ as printed in S17 the sum differs from a₁a₂ − a₃ by 2J₂₃(J₁₁J₃₂ − J₃₁J₁₂)", iszero_sym(P1print + P2 - (a1 * a2 - a3) - 2J23 * X))
    @variables aa bb
    S19 = V^2 / (2N1 * w^2) + N1 / K * (aa * V / w + 1 / N1)^2 - v^2 / (N1 * V)
    asub = Dict(aa => 1 + 2 / n, bb => 1 + 2 / n + 4 / n^2)
    check("S19 is S17 as printed, evaluated with S8 (a = 1 + 2/n_e)", iszero_sym(at8(P1print) - Symbolics.substitute(S19, asub)))
    P1eq = -V^2 / (2N1 * w^2) + N1 / K * (aa * V / w + 1 / N1)^2 + v^2 / (N1 * V)
    check("the derived P₁ at the equilibrium: −V*²/(2N*ω⁴) + (N*/K)(aV*/ω² + 1/N*)² + v²/(N*V*)",
          iszero_sym(at8(P1der) - Symbolics.substitute(P1eq, asub)))
    note("S17: P₁ is printed as J₂₃(J₁₁J₃₂ − J₃₁J₁₂) − J₃₃(J₁₁ + J₂₂)². Expanding a₁a₂ − a₃ with S16 gives " *
         "−J₂₃(J₁₁J₃₂ − J₃₁J₁₂) − J₃₃(J₁₁ + J₂₂)² (P₂ of S18 as printed): the first term has the opposite sign. S19 " *
         "follows S17 as printed, so its first and third terms, V*²/(2N*ω⁴) and −v²/(N*V*), have the opposite signs " *
         "to the derived ones. S20 and S21 are correct algebra given S19 as printed (checked below), but S20 bounds v² " *
         "from above, which bounds the derived P₁ from above, not below, so S19–S21 do not establish P₁ > 0 as written.")
    # S20, S21 as printed
    check("S20: h′(V*) < 0 ⇔ v² < V*³/(2ω⁴) + V*³N₁′/(Kω²) (h′ is linear in v²: S14)",
          iszero_sym(dheq * V / 2 - (v^2 - V^3 / (2w^2) - V^3 * dN1 / (K * w))))
    v2max = V^3 / (2w^2) + V^3 * dN1 / (K * w)
    S21_1 = N1 / K * (aa * V / w + 1 / N1)^2 - V^2 / (K * N1 * w) * dN1
    check("S21 line 1: S19 at v² = the bound of S20", iszero_sym(V^2 / (2N1 * w^2) + N1 / K * (aa * V / w + 1 / N1)^2 - v2max / (N1 * V) - S21_1))
    S21_2 = N1 / K * (aa * V / w + 1 / N1)^2 - N1 / (K * w) * (M + V^2 / c)
    check("S21 line 2: V*²N₁′ = N*²[V_M + V*²/(n_e ω²)]", iszero_sym(S21_1 - S21_2))
    S21_3 = N1 / K * (aa * V / w + 1 / N1)^2 - N1 / (K * w) * (2V^2 / c + V / N1)
    check("S21 line 3: V_M = V*²/(n_e ω²) + V*/N*", iszero_sym(S21_2 - S21_3))
    S21_4 = N1 / K * (bb * V^2 / w^2 + (2aa - 1) * V / (w * N1) + 1 / N1^2)
    check("S21 line 4: b = 1 + 2/n_e + 4/n_e² (a² − 2/n_e)", iszero_sym(Symbolics.substitute(S21_3 - S21_4, asub)))
    # the conclusion holds: a decomposition of a₁a₂ − a₃ into terms of one sign
    s = J11 + J22
    pos = J23 * J31 * J12 + J22 * J23 * J32 + J33 * J23 * J32 - s * (J11 * J22 + J33 * (s + J33))
    check("a₁a₂ − a₃ = J₂₃J₃₁J₁₂ + J₂₂J₂₃J₃₂ + J₃₃J₂₃J₃₂ − (J₁₁+J₂₂)(J₁₁J₂₂ + J₃₃(J₁₁+J₂₂+J₃₃))", iszero_sym(a1 * a2 - a3 - pos))
    note("The conclusion of S2 (the equilibrium V₂ is locally stable whenever it exists) stands: in the identity " *
         "above, with J₁₁, J₂₂, J₃₃, J₁₂, J₃₁, J₃₂ < 0 < J₂₃ (S8, v > 0), every term is positive (the J₁₁J₂₃J₃₂ parts of " *
         "P₁ and P₂ cancel), so a₁a₂ − a₃ > 0 at every equilibrium with N* > 0, and stability is decided by the sign of " *
         "det J alone (S15). No result of the preprint depends on the sign in S17.")

    # ===================================================================================
    say("\n(f) invariances")
    @variables a b
    fld(d) = [Symbolics.substitute(F, d) for F in Fs]
    tr = fld(Dict(Δ => a * Δ, V => a^2 * V, v => a * v, w => a^2 * w, M => a^2 * M))
    check("trait units: (Δ, v) × a, (V, ω², V_M) × a² scales (dΔ, dV, dN) by (a, a², 1)", iszero_sym(tr .- [a, a^2, 1] .* Fs))
    tm = fld(Dict(r => b * r, w => w / b, M => b * M, K => K / b^2, v => b * v, N => N / b))
    check("time units: (r, V_M, v) × b, ω²/b, K/b², N/b scales (dΔ, dV, dN) by (b, b, 1)", iszero_sym(tm .- [b, b, 1] .* Fs))
    check("K → ∞ at fixed V: h → (2V²/ω²)(r − V/(2ω²)), maximal at V = 4rω²/3",
          iszero_sym(Symbolics.substitute(D(V)(2V^2 / w * (r - V / (2w))), Dict(V => 4r * w / 3))))
    check("xu_shen_balance_variance solves eq. 5 at fixed N (dV/dt = 0 numerically, 20 points)",
          all(1:20) do _
              mm = XuShen(VM = 10.0^runif(rng, -3, 0), ne = rand(rng, (1.0, 5.0, 100.0)), ω2 = 10.0^runif(rng, 0, 3))
              Nn = 10.0^runif(rng, -1, 6)
              abs(num(FV, Dict(N => Nn, M => mm.VM, n => mm.ne, w => mm.ω2, V => xu_shen_balance_variance(mm, Nn)))) < 1e-12 * mm.VM
          end)

    # ===================================================================================
    say("\n(g) numerical checks against src/xu_shen.jl")
    nstab = npos = nXpos = neq = 0
    P1min = Inf
    for trial in 1:60
        m = XuShen(r = 0.05 + rand(rng), ω2 = 10.0^runif(rng, 0, 3), VM = 10.0^runif(rng, -3, 0), K = 10.0^runif(rng, 2, 5),
                   ne = rand(rng, (1, 2, 10, 100)))
        d0 = Dict(r => m.r, w => m.ω2, M => m.VM, K => m.K, n => m.ne)
        cr = xu_shen_critical_rate(m)
        if trial <= 8
            tag = @sprintf("r=%.2f ω²=%.3g V_M=%.3g K=%.3g n_e=%g", m.r, m.ω2, m.VM, m.K, m.ne)
            Vt = runif(rng, 0.05, 0.95) * xu_shen_vmax(m)
            dV = merge(d0, Dict(V => Vt))
            check("$tag: _xs_h, _xs_phi, _xs_dN1 equal h, φ, N₁′",
                  isapprox(QTLSim._xs_h(m, Vt), num(h, dV); rtol = 1e-10) && isapprox(QTLSim._xs_phi(m, Vt), num(φ, dV); rtol = 1e-9, atol = 1e-12) &&
                  isapprox(QTLSim._xs_dN1(m, Vt), num(dN1, dV); rtol = 1e-10))
            check(@sprintf("%s: h′(V_c) = 0 at xu_shen_critical_rate's V_c (φ = %.1e) and v_c² = h(V_c)", tag, num(φ, merge(d0, Dict(V => cr.V)))),
                  abs(num(φ, merge(d0, Dict(V => cr.V)))) < 1e-8 * m.r && isapprox(cr.v^2, num(h, merge(d0, Dict(V => cr.V))); rtol = 1e-10))
            brute = maximum(x -> num(h, merge(d0, Dict(V => x))), range(0, xu_shen_vmax(m), length = 4002)[2:end-1])
            check("$tag: v_c² is the maximum of h (grid of 4000, 1e-5)", isapprox(cr.v^2, brute; rtol = 1e-5) && cr.v^2 >= brute * (1 - 1e-12))
            y = (runif(rng, 0, 5), Vt, 10.0^runif(rng, 0, 3))
            dy = merge(d0, Dict(Δ => y[1], V => y[2], N => y[3], v => 0.1))
            check("$tag: xu_shen_jacobian equals the symbolic Jacobian at a random state",
                  isapprox(xu_shen_jacobian(m, y...), [num(e, dy) for e in Jgen]; rtol = 1e-10))
            check("$tag: xu_shen_field equals S1", isapprox(collect(xu_shen_field(m, y...; v = 0.1)), [num(e, dy) for e in Fs]; rtol = 1e-12))
        end
        vv = runif(rng, 0.05, 0.98) * cr.v
        for e in xu_shen_equilibria(m, vv)
            neq += 1
            de = merge(d0, Dict(V => e.V, v => vv))
            st = xu_shen_stability(m, e.Δ, e.V, e.N)
            Jn = xu_shen_jacobian(m, e.Δ, e.V, e.N)
            Xn = Jn[1, 1] * Jn[3, 2] - Jn[3, 1] * Jn[1, 2]
            P1n = -Jn[2, 3] * Xn - Jn[3, 3] * (Jn[1, 1] + Jn[2, 2])^2
            nXpos += Xn > 0
            e.stable && (P1min = min(P1min, P1n / st.hurwitz))
            npos += st.hurwitz > 0 && st.a1 > 0
            nstab += (e.stable == (e.V > cr.V)) && (st.a3 > 0) == (num(dh, merge(de, Dict(r => m.r))) < 0) &&
                     isapprox(-st.a3, num(dh, merge(de, Dict(r => m.r))) / (2e.N); rtol = 1e-6, atol = 1e-14) &&
                     abs(num(h, merge(de, Dict(r => m.r))) - vv^2) < 1e-9 * vv^2
        end
    end
    check("60 random parameter sets, $neq equilibria: a₁ > 0 and a₁a₂ − a₃ > 0 at all", npos == neq && neq == 120)
    check("… stable ⇔ V* > V_c ⇔ h′(V*) < 0, det J = h′/(2N*) and v² = h(V*) (S4, S15) at all", nstab == neq)
    say(@sprintf("  J₁₁J₃₂ − J₃₁J₁₂ > 0 (where the sign in S17 matters) at %d of %d equilibria; smallest derived P₁/(a₁a₂ − a₃) at a stable one: %.3g",
                 nXpos, neq, P1min))
    P1min < 0 && note("the derived P₁ is negative at some stable equilibria of the random sets (smallest P₁/(a₁a₂ − a₃) = " * @sprintf("%.3g", P1min) *
                      "); there P₂ makes up the difference, as the decomposition in (e) shows.")

    say("\nNotes on the preprint's equations (", length(notes), "):")
    for (i, s) in enumerate(notes)
        say("  ", i, ". ", s)
    end
    finish(LOG, t0)
end

runmain(main, @__FILE__)
