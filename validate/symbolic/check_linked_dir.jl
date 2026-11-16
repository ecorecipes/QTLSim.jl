# Symbolic and exact checks of the directional terms of the allele–background closure
# (validate/linked_dir.jl, validate/linked_dir_theory.jl), and numerical spot checks of its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_linked_dir.jl
# Writes validate/symbolic/check_linked_dir.log; exits 1 if any check fails.
#
# Notation (chromosome level, centred moments): A = B_s (a class of partners), R the rest of the
# background, B = A + R, X = x_l, d = x_l − x̄_l, G = B + X (chromosome genotype), y_A = Cov(A,X),
# Y = Cov(B,X), v = Var(X), T_A = κ(A,X,X), C = κ(B,X,X), J_A = κ(A,B,X), H_A = Cov(AB, d²),
# K_A = κ(A,B,X,X) = H_A − 2y_AY, K3_A = κ(A,X,X,X).
#
# (a) the measurement's derivatives (dir_measure!): under weights 1 + εh (E[h] = 0) the derivatives at
#     ε = 0 of y, T, J, H, K, yY equal the central-moment formulas used in the code (exact rationals);
# (b) one biallelic locus (values 0, a; frequency p): d² = s d + v, d³ = (s² + v)d + sv (s = κ3/v),
#     κ4 = v(s² − 2v), κ5 = v s(s² − 8v) (symbolic);
# (c) biallelic X, any background (exact rational populations): T_A = s y_A, H_A = s J_A, K3_A = w y_A
#     (w = s² − 2v), and s κ(A,B,B,X) = Cov(AB², d²) − 2Cov(A,B) C − Var(B) T_A;
# (d) the directional change (†): under weights 1 + ε(G − Ḡ), biallelic X: dH_A/dε = w J_A − 2Y u_A +
#     s κ(A,B,B,X), ds/dε = −2(v + Y), dT_A/dε = K_A + K3_A, dy_A/dε = T_A + J_A, dJ_A/dε = K_A + κ(A,B,B,X),
#     dK_A/dε = κ₅(A,B,X,X,G) (exact rationals);
# (e) Gaussian background: (A, B) jointly Gaussian and E[d² | B] = v + γ(B − B̄) give
#     Cov(AB², d²) − 2Cov(A,B) C − Var(B) T_A = 0 (Isserlis, symbolic), so with a Gaussian background
#     (†) reduces to w J_A − 2Y u_A;
# (f) two biallelic loci, exact joint law P_kP_l(1 + D d_kd_l/(v_kv_l)) (all cumulants exact): κ(k,k,l) = s_k D;
#     κ(k,k,l,l) = s_ks_lD − 2D², κ(k,l,l,l) = w_l D, κ(k,k,k,l,l) = s_l w_k D − 6 s_k D²
#     (exact polynomials in D at random rational p_k, p_l, a_k, a_l; three trials);
# (g) iid sampling of M: E[T̂_A ŷ_B] = (1 − 5/M) T_A y_B − (2/M) y_A T_B + (μ(A,B,X,X,X) − v μ(A,B,X))/M
#     + O(1/M²), and the biallelic seed μ(A,B,X³) − vμ(A,B,X) = s v Cov(A,B) + s² μ(A,B,X);
# (h) spot checks of validate/linked_dir_theory.jl: `dir_drivers`, `k3_source`, `qt_increments`,
#     `qt_step!` against an explicit evaluation of the recursion.

include("common.jl")
include(joinpath(@__DIR__, "..", "linked_decay_theory.jl"))
include(joinpath(@__DIR__, "..", "linked_source_theory.jl"))
include(joinpath(@__DIR__, "..", "linked_dir_theory.jl"))

const LOG = joinpath(@__DIR__, "check_linked_dir.log")
const rng = Xoshiro(20261118)

# ---------------------------------------------------------------------------------------
# expectations of products of sample means under M iid draws (as check_linked_source.jl)

const NV = 3
const Mon = NTuple{NV,Int}
const Poly = Dict{Vector{Mon},Any}
Amon(e::Mon) = Poly([e] => 1)
pone() = Poly(Mon[] => 1)
function padd(p::Poly, q::Poly, c = 1)
    r = copy(p)
    for (k, v) in q
        r[k] = get(r, k, 0) + c * v
    end
    r
end
function pmul(p::Poly, q::Poly)
    r = Poly()
    for (k1, v1) in p, (k2, v2) in q
        k = sort(vcat(k1, k2))
        r[k] = get(r, k, 0) + v1 * v2
    end
    r
end
ppow(p, n) = n == 0 ? pone() : pmul(p, ppow(p, n - 1))
pscale(p, c) = Poly(k => c * v for (k, v) in p)
unitm(i) = ntuple(j -> j == i ? 1 : 0, NV)
function central(e::Mon)
    acc = Poly()
    for idx in Iterators.product((0:e[i] for i in 1:NV)...)
        c = prod(binomial(e[i], idx[i]) * (-1)^(e[i] - idx[i]) for i in 1:NV)
        term = pscale(Amon(Mon(idx)), c)
        for i in 1:NV
            term = pmul(term, ppow(Amon(unitm(i)), e[i] - idx[i]))
        end
        acc = padd(acc, term)
    end
    acc
end
function setpartitions(n)
    n == 0 && return [Vector{Vector{Int}}()]
    out = Vector{Vector{Vector{Int}}}()
    for p in setpartitions(n - 1)
        for i in eachindex(p)
            q = deepcopy(p); push!(q[i], n); push!(out, q)
        end
        q = deepcopy(p); push!(q, [n]); push!(out, q)
    end
    out
end
falling(M, r) = prod((M - i for i in 0:r-1); init = 1)
function expect(p::Poly, M, μ)
    tot = 0
    for (fac, c) in p
        k = length(fac)
        s = 0
        for π in setpartitions(k)
            term = falling(M, length(π)) / M^k
            for B in π
                term *= μ(reduce((x, y) -> x .+ y, fac[B]))
            end
            s += term
        end
        tot += c * s
    end
    tot
end

# ---------------------------------------------------------------------------------------
# exact rational populations and joint cumulants (any order, by the partition formula)

randq() = Rational{BigInt}(rand(rng, -9:9), rand(rng, 1:5))
wmean(w, f) = sum(w .* f) / sum(w)
cm(w, vs...) = (c = [v .- wmean(w, v) for v in vs]; wmean(w, reduce(.*, c)))
function kap(w, vs...)
    n = length(vs)
    c = [v .- wmean(w, v) for v in vs]
    tot = zero(eltype(w))
    for π in setpartitions(n)
        length(π) >= 1 || continue
        any(B -> length(B) == 1, π) && continue            # centred: blocks of size 1 vanish
        k = length(π)
        term = (-1)^(k - 1) * factorial(k - 1) * one(eltype(w))
        for B in π
            term *= wmean(w, reduce(.*, (c[i] for i in B)))
        end
        tot += term
    end
    tot
end
Hof(w, A, B, X) = cm(w, A, B, X, X) - cm(w, A, B) * cm(w, X, X)

# a random population of n chromosomes: A, R arbitrary rationals, X biallelic (values 0, a) if `bial`
function randpop(n; bial = true)
    a = randq(); a == 0 && (a = 1 // 1)
    A = [randq() for _ in 1:n]; R = [randq() for _ in 1:n]
    X = bial ? [rand(rng, Bool) ? a : zero(a) for _ in 1:n] : [randq() for _ in 1:n]
    allequal(X) && (X[1] = X[1] == 0 ? a : zero(a))
    w = [Rational{BigInt}(rand(rng, 1:6)) for _ in 1:n]
    (; w, A, R, B = A .+ R, X)
end

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_linked_dir.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the measurement's derivatives under weights 1 + εh")
    @variables ε
    for trial in 1:2
        n = 7
        P = randpop(n; bial = false)
        h0 = [randq() for _ in 1:n]
        h = h0 .- wmean(P.w, h0)
        wε = P.w .* (1 .+ ε .* h)
        A, B, X = P.A, P.B, P.X
        stat(w) = (y = cm(w, A, X), T = cm(w, A, X, X), J = cm(w, A, B, X), H = Hof(w, A, B, X),
                   K = Hof(w, A, B, X) - 2cm(w, A, X) * cm(w, B, X), yY = cm(w, A, X) * cm(w, B, X))
        Sε = stat(wε)
        der(e) = simp(Symbolics.substitute(Symbolics.expand_derivatives(Differential(ε)(e)), Dict(ε => 0)))
        w = P.w
        Eh(f) = wmean(w, h .* f)
        c(v) = v .- wmean(w, v)
        a_, b_, d_ = c(A), c(B), c(X)
        v = cm(w, X, X); y = cm(w, A, X); Y = cm(w, B, X); T = cm(w, A, X, X); J = cm(w, A, B, X)
        Cl = cm(w, B, X, X); CBB = cm(w, A, B)
        δx = Eh(d_); δv = Eh(d_ .^ 2); δB = Eh(b_); δY = Eh(b_ .* d_); δc = Eh(a_)
        δy = Eh(a_ .* d_)
        δT = Eh(a_ .* d_ .^ 2) - δc * v - 2δx * y
        δJ = Eh(a_ .* b_ .* d_) - δc * Y - δB * y - δx * CBB
        δCBB = Eh(a_ .* b_)
        δH = Eh(a_ .* b_ .* d_ .^ 2) - δc * Cl - δB * T - 2δx * J - δCBB * v - CBB * δv
        δK = δH - 2δy * Y - 2y * δY
        δyY = δy * Y + y * δY
        ok = all(iszero_sym(der(getfield(Sε, s)) - f) for (s, f) in ((:y, δy), (:T, δT), (:J, δJ), (:H, δH), (:K, δK), (:yY, δyY)))
        check("trial $trial: d/dε of y, T, J, H, K, yY = dir_measure!'s formulas (exact)", ok)
    end

    # ===================================================================================
    say("\n(b) one biallelic locus: d² = s d + v, d³ = (s² + v)d + sv, κ4 = v(s² − 2v), κ5 = v s(s² − 8v)")
    @variables a p
    vals = [a, 0]; pr = [p, 1 - p]
    mu1 = a * p
    m(k) = sum(pr[i] * (vals[i] - mu1)^k for i in 1:2)
    v = m(2); κ3 = m(3); s = κ3 / v
    ok2 = all(iszero_sym((vals[i] - mu1)^2 - (s * (vals[i] - mu1) + v)) for i in 1:2)
    ok3 = all(iszero_sym((vals[i] - mu1)^3 - ((s^2 + v) * (vals[i] - mu1) + s * v)) for i in 1:2)
    check("d² = s d + v at both allele values", ok2)
    check("d³ = (s² + v) d + s v at both allele values", ok3)
    κ4 = m(4) - 3m(2)^2
    κ5 = m(5) - 10m(3) * m(2)
    check("κ4 = v (s² − 2v)", iszero_sym(κ4 - v * (s^2 - 2v)))
    check("κ5 = v s (s² − 8v)", iszero_sym(κ5 - v * s * (s^2 - 8v)))
    check("s = a(1 − 2p), v = a²p(1 − p)", iszero_sym(s - a * (1 - 2p)) && iszero_sym(v - a^2 * p * (1 - p)))

    # ===================================================================================
    say("\n(c) biallelic X, any background (exact rational populations)")
    for trial in 1:3
        P = randpop(9)
        w, A, B, X = P.w, P.A, P.B, P.X
        v = cm(w, X, X); s = cm(w, X, X, X) / v; wl = s^2 - 2v
        y = cm(w, A, X); T = kap(w, A, X, X); J = kap(w, A, B, X); H = Hof(w, A, B, X)
        K3 = kap(w, A, X, X, X); C = kap(w, B, X, X)
        check("trial $trial: T_A = s y_A, H_A = s J_A, K3_A = w y_A", T == s * y && H == s * J && K3 == wl * y)
        d = X .- wmean(w, X)
        a_ = A .- wmean(w, A); b_ = B .- wmean(w, B)
        CovAB2 = wmean(w, a_ .* b_ .^ 2 .* d .^ 2) - wmean(w, a_ .* b_ .^ 2) * v
        lhs = s * kap(w, A, B, B, X)
        check("trial $trial: s κ(A,B,B,X) = Cov(AB², d²) − 2Cov(A,B)C − Var(B)T_A", lhs == CovAB2 - 2cm(w, A, B) * C - cm(w, B, B) * T)
    end

    # ===================================================================================
    say("\n(d) the directional change under weights 1 + ε(G − Ḡ), biallelic X (exact)")
    for trial in 1:3
        P = randpop(8)
        w, A, B, X = P.w, P.A, P.B, P.X
        G = B .+ X
        h = G .- wmean(w, G)
        wε = w .* (1 .+ ε .* h)
        der(e) = simp(Symbolics.substitute(Symbolics.expand_derivatives(Differential(ε)(e)), Dict(ε => 0)))
        v = cm(w, X, X); s = cm(w, X, X, X) / v; wl = s^2 - 2v; Y = cm(w, B, X); yA = cm(w, A, X)
        J = kap(w, A, B, X); T = kap(w, A, X, X); K = kap(w, A, B, X, X); K3 = kap(w, A, X, X, X)
        κABBX = kap(w, A, B, B, X)
        Hε = Hof(wε, A, B, X)
        vε = cm(wε, X, X); sε = cm(wε, X, X, X) / vε
        check("trial $trial: dH_A/dε = w J_A − 2Y(T_A + J_A) + s κ(A,B,B,X)  (†)", iszero_sym(der(Hε) - (wl * J - 2Y * (T + J) + s * κABBX)))
        check("trial $trial: (the form without −2Y T_A is wrong)", !iszero_sym(der(Hε) - (wl * J - 2Y * J + s * κABBX)))
        check("trial $trial: ds/dε = −2(v + Y)", iszero_sym(der(sε) - (-2(v + Y))))
        check("trial $trial: dT_A/dε = K_A + K3_A; dy_A/dε = T_A + J_A", iszero_sym(der(cm(wε, A, X, X)) - (K + K3)) &&
                                                                       iszero_sym(der(cm(wε, A, X)) - (T + J)))
        check("trial $trial: dJ_A/dε = K_A + κ(A,B,B,X)", iszero_sym(der(cm(wε, A, B, X)) - (K + κABBX)))
        Kε = Hof(wε, A, B, X) - 2cm(wε, A, X) * cm(wε, B, X)
        check("trial $trial: dK_A/dε = κ₅(A,B,X,X,G)", iszero_sym(der(Kε) - kap(w, A, B, X, X, G)))
    end

    # ===================================================================================
    say("\n(e) Gaussian background: the cubic heterogeneity cancels −2Cov(A,B)C − Var(B)T_A")
    @variables σA σB ρ γ vv
    # (A, B) = (σA(ρ Z + √(1−ρ²) W), σB Z), E[d² | B] = vv + γ B; Isserlis: E[Z⁴] = 3, E[Z²] = 1, odd = 0
    covAB = ρ * σA * σB; varB = σB^2
    C = γ * varB                      # Cov(B, d²)
    TA = γ * covAB                    # Cov(A, d²)
    CovAB2 = γ * 3 * covAB * varB     # E[A B² (d² − vv)] = γ E[A B³] = 3γ Cov(A,B) Var(B), E[AB²] = 0
    check("Cov(AB², d²) − 2Cov(A,B)C − Var(B)T_A = 0", iszero_sym(CovAB2 - 2covAB * C - varB * TA))

    # ===================================================================================
    say("\n(f) two biallelic loci, joint law P_kP_l(1 + D d_k d_l/(v_k v_l)): cumulants in powers of D")
    @variables D
    for trial in 1:3
        pk = Rational{BigInt}(rand(rng, 1:9), 10); pl = Rational{BigInt}(rand(rng, 1:9), 10)
        ak = randq(); ak == 0 && (ak = 1 // 1); al = randq(); al == 0 && (al = -1 // 2)
        xk = [ak, 0]; xl = [al, 0]; prk = [pk, 1 - pk]; prl = [pl, 1 - pl]
        mk = ak * pk; ml = al * pl
        vk = ak^2 * pk * (1 - pk); vl = al^2 * pl * (1 - pl)
        sk = ak * (1 - 2pk); sl = al * (1 - 2pl); wk = sk^2 - 2vk; wl = sl^2 - 2vl
        pj(i, j) = prk[i] * prl[j] * (1 + D * (xk[i] - mk) * (xl[j] - ml) / (vk * vl))
        Ejk(e1, e2) = sum(pj(i, j) * (xk[i] - mk)^e1 * (xl[j] - ml)^e2 for i in 1:2, j in 1:2)   # the marginals stay P_k, P_l
        check("trial $trial: κ(k,k,l) = s_k D exactly", iszero_sym(Ejk(2, 1) - sk * D))
        κkkll = Ejk(2, 2) - Ejk(2, 0) * Ejk(0, 2) - 2 * Ejk(1, 1)^2
        check("trial $trial: κ(k,k,l,l) = s_k s_l D − 2D² exactly", iszero_sym(κkkll - (sk * sl * D - 2D^2)))
        κklll = Ejk(1, 3) - 3 * Ejk(1, 1) * Ejk(0, 2)
        check("trial $trial: κ(k,l,l,l) = w_l D exactly", iszero_sym(κklll - wl * D))
        # centred five-variable cumulant k,k,k,l,l: μ32 − 3μ20μ12 − 6μ11μ21 − μ02μ30
        κkkkll = Ejk(3, 2) - 3 * Ejk(2, 0) * Ejk(1, 2) - 6 * Ejk(1, 1) * Ejk(2, 1) - Ejk(0, 2) * Ejk(3, 0)
        check("trial $trial: κ(k,k,k,l,l) = s_l w_k D − 6 s_k D² exactly", iszero_sym(κkkkll - (sl * wk * D - 6sk * D^2)))
    end

    # ===================================================================================
    say("\n(g) iid sampling of M: E[T̂_A ŷ_B]")
    @variables M ε2
    μs = Dict{Mon,Any}()
    μ(e) = sum(e) == 0 ? 1 : sum(e) == 1 ? 0 : get!(μs, e) do
        Symbolics.variable(Symbol("μ", join(e)))
    end
    TA_ = central((1, 0, 2)); yB_ = central((0, 1, 1))
    ETy = expect(pmul(TA_, yB_), M, μ)
    lead(e) = simp(Symbolics.substitute(e, Dict(M => 1 / ε2)))
    Tm = μ((1, 0, 2)); yBm = μ((0, 1, 1)); yAm = μ((1, 0, 1)); TBm = μ((0, 1, 2)); vm = μ((0, 0, 2))
    want = -5Tm * yBm - 2yAm * TBm + μ((1, 1, 3)) - vm * μ((1, 1, 1))
    check("E[T̂_A ŷ_B] = T_A y_B + [−5T_Ay_B − 2y_AT_B + μ(A,B,X³) − vμ(A,B,X)]/M + O(1/M²)",
          iszero_sym(coeffj(lead(ETy), ε2, 0) - Tm * yBm) && iszero_sym(coeffj(lead(ETy), ε2, 1) - want))
    # biallelic seed (exact rational populations)
    okseed = true
    for trial in 1:3
        P = randpop(9)
        w, A, B, X = P.w, P.A, P.B, P.X
        v = cm(w, X, X); s = cm(w, X, X, X) / v
        okseed &= cm(w, A, B, X, X, X) - v * cm(w, A, B, X) == s * v * cm(w, A, B) + s^2 * cm(w, A, B, X)
    end
    check("biallelic seed μ(A,B,X³) − vμ(A,B,X) = s v Cov(A,B) + s² μ(A,B,X) (3 populations)", okseed)

    # ===================================================================================
    say("\n(h) spot checks of validate/linked_dir_theory.jl")
    dd = dir_drivers([2.0, 0.0], [1.0, 0.0], [0.6, 0.0], [0.8, 0.0])
    check("dir_drivers: ⟨w⟩ = Σκ4/Σv, ⟨s²⟩ = Σ(κ3²/v)/Σv, h = Σv²/(Σv)² (empty group 0)",
          dd.w ≈ [0.3, 0.0] && dd.s2 ≈ [0.4, 0.0] && dd.h ≈ [0.25, 0.0])
    check("k3_source = ⟨w⟩_l m", k3_source([1.0, -2.0], dd, 1) ≈ [0.3, -0.6])
    G2 = LDGroups(20, [3.0 2.0; 3.0 2.0], [0.01 0.2; 0.01 0.2], [10, 10], [1, 2])
    d = LDDrivers(β = 0.7, Sc = 5.0, dS = 0.0, n = 400.0, N = 1000.0, Vm = 0.0, V = [2.0, 1.0], V2 = [1.0, 0.3],
                  k3 = [0.05, -0.02], g = [2.0, 1.0], c = [0.3, 0.2], Cbb = [2.5 0.1; 0.1 1.2])
    QT0 = [0.01 0.002; -0.003 0.004]; QT = copy(QT0)
    P = [0.2 0.05; 0.05 0.1]; m = [-0.1, -0.05]; Ks = [-0.3, -0.1]
    qt_step!(QT, P, m, Ks, G2, d, dd, 1)
    Σ = class_cov(G2, d, 1)
    q = QT0 .+ 0.7 .* ((0.4 + 0.3) .* P .+ 0.25 .* (Ks * m')) .- (3 * (2.0 / 10) / 5.0) .* QT0
    R = [0.99, 0.8]
    εs = 1 / 800 - 1 / 2000; εo = 1 / 2000
    q = (1 - 5εs) .* q .- 2εs .* q' .+ (εs * 0.05) .* Σ
    q = R .* q .* R'
    q = (1 - 5εo) .* q .- 2εo .* q' .+ (εo * 0.05) .* Σ
    check("qt_step!: selection, two sampling stages (with −2εQᵀ and the κ3 seed), recombination", QT ≈ q)
    inc = qt_increments(QT, [0.1, 0.2], d)
    check("qt_increments: δP = β(Q + Qᵀ), δH = 2δP, δm = β T̄", inc.δP ≈ 0.7 .* (QT .+ QT') && inc.δH ≈ 2 .* inc.δP && inc.δm ≈ [0.07, 0.14])

    finish(LOG, t0)
end

runmain(main, @__FILE__)
