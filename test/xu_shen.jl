# Tests: xu_shen. Included by runtests.jl; run alone with
# `julia --project=vignettes test/runtests.jl xu_shen`.
#
# The continuous-time model of Xu and Shen (2026, bioRxiv preprint,
# doi:10.64898/2026.08.05.743160), re-implemented from its equations (src/xu_shen.jl).

@testset "Xu and Shen (2026): equilibria, critical rate, stability" begin
    m = XuShen()
    @test (m.r, m.ω2, m.VM, m.K, m.ne) == (0.2, 50.0, 0.1, 1000.0, 1.0)             # the preprint's standard parameters
    @test XuShen(m; K = 500).K == 500 && XuShen(m; K = 500).r == m.r
    @test_throws ArgumentError XuShen(r = -0.1)
    @test_throws ArgumentError XuShen(ω2 = Inf)
    @test xu_shen_vmax(m) ≈ sqrt(50 * 0.1)
    # the field (eqs. 2, 4, 5)
    @test collect(xu_shen_field(m, 2.0, 1.5, 80.0; v = 0.1)) ≈
          [0.1 - 1.5 * 2 / 50, 0.1 - 1.5^2 / 50 - 1.5 / 80, 80 * (0.2 - (1.5 + 4) / 100 - 80 / 1000)]
    c = xu_shen_critical_rate(m)
    h(V) = QTLSim._xs_h(m, V)
    # eq. 6: v_c² is the maximum of h over (0, V_c^M), by brute force
    Vs = range(0, xu_shen_vmax(m), length = 200_002)[2:end-1]
    @test c.v^2 ≈ maximum(h, Vs) rtol = 1e-8
    @test all(h(V) <= c.v^2 * (1 + 1e-12) for V in Vs)
    @test c.Δ ≈ c.v * m.ω2 / c.V && c.N ≈ QTLSim._xs_N1(m, c.V)
    # equilibria satisfy the ODE to 1e-10; two below v_c, none above, one at v = 0
    for mm in (m, XuShen(ω2 = 5.0), XuShen(VM = 0.5, ne = 10), XuShen(r = 0.5, K = 1e4), XuShen(ω2 = 1000.0, K = 200.0))
        vc = xu_shen_critical_rate(mm).v
        for v in (0.0, 0.2vc, 0.7vc, 0.99vc, (1 - 1e-9) * vc)
            es = xu_shen_equilibria(mm, v)
            @test length(es) == (v == 0 ? 1 : 2)
            for e in es
                @test maximum(abs, xu_shen_field(mm, e.Δ, e.V, e.N; v)) < 1e-10
                @test 0 < e.V < xu_shen_vmax(mm) && e.N > 0
                @test e.Δ ≈ v * mm.ω2 / e.V                                                    # eq. 7
                @test e.N ≈ mm.K * (mm.r - (e.V + e.Δ^2) / (2mm.ω2)) rtol = 1e-8               # eq. 8
            end
            @test issorted(es; by = e -> e.V)
            @test es[end].stable && (v == 0 || !es[1].stable)
            @test xu_shen_equilibrium(mm, v).V == es[end].V
        end
        # the equilibrium disappears exactly at the rate from the maximisation
        @test length(xu_shen_equilibria(mm, vc * (1 + 1e-9))) == 0
        @test xu_shen_equilibrium(mm, 1.01vc) === nothing
        # at v_c the two roots merge at the maximiser
        e2 = xu_shen_equilibria(mm, (1 - 1e-9) * vc)
        cc = xu_shen_critical_rate(mm)
        @test e2[1].V ≈ cc.V rtol = 1e-3
        @test e2[2].V ≈ cc.V rtol = 1e-3
    end
    @test_throws ArgumentError xu_shen_equilibria(m, -0.1)
    @test_throws ArgumentError xu_shen_equilibria(XuShen(K = Inf), 0.05)

    # the Jacobian (S8) against finite differences of the field, at equilibria and elsewhere
    function fdjac(mm, y, v)
        J = zeros(3, 3)
        for j in 1:3
            d = 1e-6 * max(abs(y[j]), 1e-3)
            yp = collect(y); ym = collect(y); yp[j] += d; ym[j] -= d
            J[:, j] = (collect(xu_shen_field(mm, yp...; v)) .- collect(xu_shen_field(mm, ym...; v))) ./ 2d
        end
        J
    end
    for v in (0.05, 0.08, 0.12), e in xu_shen_equilibria(m, v)
        J = xu_shen_jacobian(m, e.Δ, e.V, e.N)
        @test J ≈ fdjac(m, (e.Δ, e.V, e.N), v) rtol = 1e-6 atol = 1e-9
        @test J[1, 2] ≈ -v / e.V && J[3, 1] ≈ -e.N * v / e.V && J[3, 3] ≈ -e.N / m.K                # the forms printed in S8
        @test J[1, 3] == 0 && J[2, 1] == 0
        s = xu_shen_stability(m, e.Δ, e.V, e.N)
        # Routh–Hurwitz agrees with the eigenvalues
        @test s.stable == all(real(λ) < 0 for λ in s.eigenvalues)
        @test s.a1 ≈ -sum(s.eigenvalues) && s.a3 ≈ -real(prod(s.eigenvalues))
        @test s.a1 > 0 && s.hurwitz > 0                                # these two always hold; a₃ decides
        # S15: det J = h′(V*)/(2N*)
        dh = (h(e.V * (1 + 1e-6)) - h(e.V * (1 - 1e-6))) / (2e-6 * e.V)
        @test -s.a3 ≈ dh / (2e.N) rtol = 1e-5
        @test e.stable == (e.V > c.V)                                  # the falling side of h
    end
    @test xu_shen_jacobian(m, 1.0, 1.5, 80.0) ≈ fdjac(m, (1.0, 1.5, 80.0), 0.07) rtol = 1e-6       # away from equilibrium
    # a scan: Routh–Hurwitz and the eigenvalues agree at every equilibrium of random parameter sets
    rng = Xoshiro(20261118)
    nagree = 0
    for _ in 1:200
        mm = XuShen(r = 0.05 + rand(rng), ω2 = 10.0^(3rand(rng)), VM = 10.0^(-3 + 3rand(rng)), K = 10.0^(2 + 3rand(rng)),
                    ne = rand(rng, (1, 2, 10, 100)))
        v = rand(rng) * xu_shen_critical_rate(mm).v
        for e in xu_shen_equilibria(mm, v)
            s = xu_shen_stability(mm, e.Δ, e.V, e.N)
            nagree += s.stable == all(real(λ) < 0 for λ in s.eigenvalues) && s.hurwitz > 0 && s.a1 > 0
        end
    end
    @test nagree == 400
end

@testset "Xu and Shen (2026): limits and invariances" begin
    m = XuShen()
    # K → ∞: no drift, V → ω√(n_e V_M) and the classical lag vω²/V; v_c → the classical critical rate
    big = XuShen(m; K = 1e13)
    e = xu_shen_equilibrium(big, 0.05)
    @test e.V ≈ xu_shen_vmax(m) rtol = 1e-6
    @test e.Δ ≈ 0.05 * m.ω2 / xu_shen_vmax(m) rtol = 1e-6
    @test e.N > 1e11
    ci = xu_shen_critical_rate(XuShen(m; K = Inf))
    Vm = xu_shen_vmax(m)
    @test ci.V == Vm && ci.v ≈ Vm / sqrt(m.ω2) * sqrt(2m.r - Vm / m.ω2) && ci.N == Inf
    @test xu_shen_critical_rate(big).v ≈ ci.v rtol = 1e-4
    @test xu_shen_critical_rate(m).v < ci.v                              # drift lowers the critical rate
    # K = ∞ with a wide fitness function: the maximum of (2V²/ω²)(r − V/(2ω²)) is interior, at 4rω²/3
    wide = XuShen(r = 0.01, ω2 = 2.0, VM = 0.5, K = Inf)
    @test xu_shen_critical_rate(wide).V ≈ 4 * 0.01 * 2 / 3
    # balance variance: the root of eq. 5 at fixed N; limits
    for N in (0.5, 10.0, 1e3, 1e9)
        V = xu_shen_balance_variance(m, N)
        @test abs(m.VM - V^2 / (m.ne * m.ω2) - V / N) < 1e-12
    end
    @test xu_shen_balance_variance(m, 1e15) ≈ Vm
    @test xu_shen_balance_variance(m, 1e-6) ≈ m.VM * 1e-6 rtol = 1e-6
    # invariance 1, units of the trait: (Δ, v) × a, (V, ω², V_M) × a²
    a = 3.7
    ma = XuShen(m; ω2 = a^2 * m.ω2, VM = a^2 * m.VM)
    c, ca = xu_shen_critical_rate(m), xu_shen_critical_rate(ma)
    @test ca.v ≈ a * c.v && ca.V ≈ a^2 * c.V && ca.N ≈ c.N
    e, ea = xu_shen_equilibrium(m, 0.08), xu_shen_equilibrium(ma, 0.08a)
    @test ea.Δ ≈ a * e.Δ && ea.V ≈ a^2 * e.V && ea.N ≈ e.N
    # invariance 2, units of time: (r, V_M, v) × b, ω²/b, K/b², N/b
    b = 2.5
    mb = XuShen(r = b * m.r, ω2 = m.ω2 / b, VM = b * m.VM, K = m.K / b^2, ne = m.ne)
    cb = xu_shen_critical_rate(mb)
    @test cb.v ≈ b * c.v && cb.V ≈ c.V && cb.N ≈ c.N / b
    eb = xu_shen_equilibrium(mb, 0.08b)
    @test eb.Δ ≈ e.Δ && eb.V ≈ e.V && eb.N ≈ e.N / b
    f, fb = xu_shen_field(m, 3.0, 1.2, 60.0; v = 0.08), xu_shen_field(mb, 3.0, 1.2, 60.0 / b; v = 0.08b)
    @test collect(fb) ≈ [b * f[1], b * f[2], f[3]]                       # d/dt × b, and N/b
    # n_e enters only through n_e ω² in eq. 5: V_c^M ∝ √n_e
    @test xu_shen_vmax(XuShen(m; ne = 4)) ≈ 2Vm
end

@testset "Xu and Shen (2026): trajectories, basins, scenarios" begin
    m = XuShen()
    v = 0.08
    es = xu_shen_equilibria(m, v)
    sad, e = es[1], es[2]
    # from near the stable equilibrium the trajectory converges to it
    tr = xu_shen_trajectory(m, 1.1e.Δ, 0.9e.V, 1.2e.N; v)
    @test tr.fate === :persists && tr.equilibrium.V == e.V
    @test tr.Δ[end] ≈ e.Δ rtol = 2e-3
    @test tr.V[end] ≈ e.V rtol = 2e-3
    @test tr.N[end] ≈ e.N rtol = 2e-3
    @test issorted(tr.t) && tr.t[1] == 0 && tr.N[1] ≈ 1.2e.N
    # at the equilibrium it stays; with save = false only the end point
    t0 = xu_shen_trajectory(m, e.Δ, e.V, e.N; v, save = false)
    @test t0.fate === :persists && length(t0.t) == 1
    # the integrator against the exact solution of a decoupled case: N at V, Δ held near equilibrium is
    # checked through a tolerance study instead: tighter tolerance, same trajectory
    ta = xu_shen_trajectory(m, 5.0, 2.0, 100.0; v, rtol = 1e-6, tmax = 50.0)
    tb = xu_shen_trajectory(m, 5.0, 2.0, 100.0; v, rtol = 1e-10, tmax = 50.0)
    @test ta.t[end] == tb.t[end] == 50.0 && ta.fate === :undecided
    @test ta.N[end] ≈ tb.N[end] rtol = 1e-4
    @test ta.V[end] ≈ tb.V[end] rtol = 1e-4
    @test ta.Δ[end] ≈ tb.Δ[end] rtol = 1e-4
    # logistic growth alone has a closed form: Δ = 0 = v and V at its balance is not constant, so
    # test the N equation with a huge ω² (no load): N(t) = rK/(1 + (rK/N0 − 1)e^{−rt})
    flat = XuShen(r = 0.3, ω2 = 1e12, VM = 1e-6, K = 500.0)
    tl = xu_shen_trajectory(flat, 0.0, 1e-3, 5.0; v = 0.0, tmax = 20.0, rtol = 1e-9)
    @test tl.N[end] ≈ 150 / (1 + (150 / 5 - 1) * exp(-0.3 * 20)) rtol = 1e-5
    # bistability: inside and outside the basin (cf. Figure 3c: v = 0.08, Δ₀ = 5)
    @test xu_shen_persists(m, 5.0, 2.0, 100.0; v)
    @test !xu_shen_persists(m, 5.0, 0.5, 30.0; v)
    @test xu_shen_trajectory(m, 5.0, 0.5, 30.0; v).fate === :extinct
    Vb = xu_shen_boundary(m, 50.0; Δ0 = 5.0, v)
    @test 0.5 < Vb < 2.0
    @test xu_shen_persists(m, 5.0, 1.02Vb, 50.0; v) && !xu_shen_persists(m, 5.0, 0.98Vb, 50.0; v)
    # the saddle is on the boundary: at its lag and size the critical V₀ is its variance
    @test xu_shen_boundary(m, sad.N; Δ0 = sad.Δ, v) ≈ sad.V rtol = 1e-3
    @test xu_shen_persists(m, sad.Δ, 1.01sad.V, sad.N; v) && !xu_shen_persists(m, sad.Δ, 0.99sad.V, sad.N; v)
    # the extinction threshold does not move the boundary here
    @test xu_shen_boundary(m, 50.0; Δ0 = 5.0, v, Next = 1.0) ≈ Vb rtol = 1e-6
    @test xu_shen_boundary(m, 50.0; Δ0 = 5.0, v, Next = 1e-9) ≈ Vb rtol = 1e-6
    # above the critical rate every start dies
    vc = xu_shen_critical_rate(m).v
    @test !xu_shen_persists(m, 0.0, 3.0, 500.0; v = 1.05vc)
    @test xu_shen_boundary(m, 50.0; Δ0 = 1.0, v = 1.05vc) == Inf
    # a small lag: the whole plotted range persists (Figure 3a), the boundary is at the floor
    @test xu_shen_boundary(m, 10.0; Δ0 = 0.5, v) == 1e-6 * xu_shen_vmax(m)
    @test_throws ArgumentError xu_shen_trajectory(m, 1.0, -1.0, 10.0; v)
    # scenario 1: V₀ at balance, the critical size (Figure 4a)
    strong = XuShen(m; ω2 = 2.0)
    Nc = xu_shen_critical_size(strong; Δ0 = 3.0, v = 0.05)
    @test 5 < Nc < 1e4
    @test xu_shen_persists(strong, 3.0, xu_shen_balance_variance(strong, 1.05Nc), 1.05Nc; v = 0.05)
    @test !xu_shen_persists(strong, 3.0, xu_shen_balance_variance(strong, 0.95Nc), 0.95Nc; v = 0.05)
    @test xu_shen_critical_size(m; Δ0 = 1.0, v = 0.05) == 1.0                      # the floor
    @test xu_shen_critical_size(m; Δ0 = 1.0, v = 1.05vc) == Inf
    @test xu_shen_critical_size(m; Δ0 = 5.0, v, V0 = 2.0, Nmin = 5.0) > 5
    # scenario 2: a change of rate with an abrupt shift (Figure 4b)
    @test xu_shen_max_rate(m; v0 = 0.04) == vc                                      # no shift: any v₁ < v_c
    v1 = xu_shen_max_rate(m; v0 = 0.04, δ = 5.0, steps = 25)
    @test 0 < v1 < vc
    e0 = xu_shen_equilibrium(m, 0.04)
    @test xu_shen_persists(m, e0.Δ + 5, e0.V, e0.N; v = 0.97v1) && !xu_shen_persists(m, e0.Δ + 5, e0.V, e0.N; v = 1.03v1)
    @test xu_shen_max_rate(m; v0 = 0.04, δ = 30.0) == 0
    @test xu_shen_max_rate(m; v0 = 0.08, δ = 5.0, steps = 25) < v1                  # a faster original rate, less room
    @test_throws ArgumentError xu_shen_max_rate(m; v0 = 2vc)
    # scenario 3: a temporary rise of the rate (Figure 4d)
    @test xu_shen_max_duration(m; v0 = 0.02, v1 = 0.9vc, Tmax = 5000.0) == Inf
    T = xu_shen_max_duration(m; v0 = 0.02, v1 = 0.2, steps = 25)
    @test 10 < T < 1000
    @test xu_shen_trajectory(m, e.Δ, e.V, e.N; v = 0.02, phases = [(0.2, 0.0)], save = false).fate === :persists
    e2 = xu_shen_equilibrium(m, 0.02)
    @test xu_shen_persists(m, e2.Δ, e2.V, e2.N; v = 0.02, phases = [(0.2, 0.95T)])
    @test !xu_shen_persists(m, e2.Δ, e2.V, e2.N; v = 0.02, phases = [(0.2, 1.05T)])
    @test xu_shen_max_duration(m; v0 = 0.09, v1 = 0.2, steps = 25) < T              # closer to the edge, less time
end
