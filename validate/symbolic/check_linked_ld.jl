# Symbolic and exact checks of the corrected allele–background LD recursion (validate/linked_ld.jl,
# validate/linked_ld_theory.jl), and numerical spot checks of its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_linked_ld.jl
# Writes validate/symbolic/check_linked_ld.log; exits 1 if any check fails.
#
# Notation (chromosome level): locus l with allele value x, background class B_s = Σ_{k∈s} x_k,
# y_s = Cov(B_s, x), Y = Σ_s y_s, v = Var(x), ρ_kl the probability that a gamete's crossover separates
# k and l (QTSim: at most one switch, after the lowest of Poisson(u) points).
#
# (a) gamete pool of n individuals under a single switch point (any law of the switch point):
#     Cov_pool(x_k, x_l) = (1 − ρ_kl) D_kl + ρ_kl Z_kl exactly, D the survivors' within-chromosome
#     covariance, Z_kl = (1/2n) Σ_i (d_k^{i2} d_l^{i1} + d_k^{i1} d_l^{i2}) the cross-chromosome one;
#     summed over a class: y_gam = y − Cov(R, x) + Z, R = Σ_k ρ_kl x_k (as `gamete_pool`);
# (b) random pairing of 2n chromosomes into n individuals: E[Z] = −Cov(R, x)/(2n − 1) exactly
#     (enumeration of all pairings), and 2n E[(Z − E Z)²] → v Var(R) + Cov(R, x)² (leading order);
# (c) selection on individuals, weights 1 + h(G₁ + G₂), chromosomes independent: the between-
#     chromosome association created is Cov_w(B^{(1)}, x^{(2)}) = −Cov(B, G)Cov(x, G)/S at first
#     order (quadratic h = −(G − Ḡ)²/(2S) + b(G − Ḡ)); the within-chromosome change under the tilt
#     1 + b u + q(u² − V) is b κ(B,x,G) + q(κ(B,x,G,G) + 2c_B c_x) − b² c_B c_x + O(bq, q², b³),
#     i.e. with q = (b² − 1/S)/2: −c_B c_x/S + b κ(B,x,G) + (b² − 1/S)κ(B,x,G,G)/2;
# (d) mutation adding independent δ (Var m_k per locus) to each of M chromosomes:
#     E[(ΔŶ)²] = (Var(B) m_l + v Σ_{k≠l} m_k)/M + (M − 1) m_l Σ_{k≠l} m_k / M² exactly;
# (e) the composite gamete map: diag(1 − ρ) A + diag(ρ)(A − I) = diag(1 − ρ) + (A − I): selection's
#     change passes to the gamete in full, recombination acts on the inherited LD only;
# (g) heterogeneity among loci: under y_l ↦ A y_l + b v_l the stationary mean is ∝ v_l, so the summed
#     moments' cross term needs Σ_l E[y_l] v_l = (Σv²/Σv) m, not v̄ m (exact identities);
# (f) spot checks of the code: `rsum!`, `class_seg`/`class_sums!`, `gamete_pool` against an explicit
#     enumeration of the pool, `crossmom` (rank-one) against the dense form, `yfix_step!` with every
#     correction off equal to `y_step!` (level B) and to `yH_step!` (level B′), `class_covR`.

include("common.jl")
include(joinpath(@__DIR__, "..", "linked_decay_theory.jl"))
include(joinpath(@__DIR__, "..", "linked_source_theory.jl"))
include(joinpath(@__DIR__, "..", "linked_ld_theory.jl"))

const LOG = joinpath(@__DIR__, "check_linked_ld.log")
const rng = Xoshiro(20261116)

randq() = Rational{BigInt}(rand(rng, -9:9), rand(rng, 1:5))
mean_(v) = sum(v) / length(v)
cov_(a, b) = mean_((a .- mean_(a)) .* (b .- mean_(b)))

# all perfect matchings of 1:2n
function matchings(v::Vector{Int})
    isempty(v) && return [Tuple{Int,Int}[]]
    a = v[1]; out = Vector{Vector{Tuple{Int,Int}}}()
    for j in 2:length(v)
        rest = [v[k] for k in 2:length(v) if k != j]
        for m in matchings(rest)
            push!(out, [(a, v[j]); m])
        end
    end
    out
end

# Minimal stand-ins so that the measurement code of linked_ld.jl can be checked without a World:
# the pieces it needs from validate/linked_ld.jl are re-included below from a header-free copy.
function load_measurement_code()
    src = read(joinpath(@__DIR__, "..", "linked_ld.jl"), String)
    # take the definitions of STBuf, moments, rsum!, gamete_pool, class_seg, class_sums!
    pick(name_start, name_end) = (i = findfirst(name_start, src); j = findnext(name_end, src, last(i)); src[first(i):first(j)-1])
    code = pick("mutable struct STBuf", "function load_pool!") *
           pick("\"\"\"\nPer-locus moments", "\"\"\"\nThe recursion's class-averaged") *
           pick("struct CSeg", "mutable struct CSBuf")
    include_string(@__MODULE__, code)
end

load_measurement_code()

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_linked_ld.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the gamete pool under a single switch point")
    for trial in 1:3
        L = 6; n = 4
        chrom = [[[randq() for _ in 1:L] for _ in 1:2] for _ in 1:n]      # chrom[i][c][k]
        # any law of the switch point: P(no switch) p0, P(switch after locus c) pc[c], c = 1..L−1
        wts = [Rational{BigInt}(rand(rng, 1:5)) for _ in 1:L]
        wts ./= sum(wts)
        p0 = wts[1]; pc = wts[2:end]
        # ρ_kl = P(switch point between k and l)
        ρ(k, l) = k == l ? 0 // 1 : sum(pc[min(k, l):max(k, l)-1])
        # enumerate the pool: each individual, start strand, switch point
        pool = Tuple{Rational{BigInt},Vector{Rational{BigInt}}}[]
        for i in 1:n, st in 1:2
            o = 3 - st
            push!(pool, (p0 / (2n), chrom[i][st]))
            for c in 1:L-1
                push!(pool, (pc[c] / (2n), [k <= c ? chrom[i][st][k] : chrom[i][o][k] for k in 1:L]))
            end
        end
        Wp = sum(first, pool)
        pm(k) = sum(w * x[k] for (w, x) in pool) / Wp
        pcov(k, l) = sum(w * (x[k] - pm(k)) * (x[l] - pm(l)) for (w, x) in pool) / Wp
        flat = [chrom[i][c] for i in 1:n for c in 1:2]
        sm(k) = mean_([x[k] for x in flat])
        D(k, l) = mean_([(x[k] - sm(k)) * (x[l] - sm(l)) for x in flat])
        Z(k, l) = sum((chrom[i][2][k] - sm(k)) * (chrom[i][1][l] - sm(l)) + (chrom[i][1][k] - sm(k)) * (chrom[i][2][l] - sm(l)) for i in 1:n) / (2n)
        ok = all(pcov(k, l) == (1 - ρ(k, l)) * D(k, l) + ρ(k, l) * Z(k, l) for k in 1:L, l in 1:L if k != l)
        check("trial $trial: Cov_pool(x_k, x_l) = (1 − ρ_kl) D_kl + ρ_kl Z_kl, all pairs (exact)", ok)
        okm = all(pm(k) == sm(k) for k in 1:L)
        check("trial $trial: the pool's allele means are the survivors'", okm)
        l = 3
        Ygam = sum(pcov(k, l) for k in 1:L if k != l)
        Y = sum(D(k, l) for k in 1:L if k != l)
        CRx = sum(ρ(k, l) * D(k, l) for k in 1:L if k != l)
        Zl = sum(ρ(k, l) * Z(k, l) for k in 1:L if k != l)
        check("trial $trial: summed over partners, Y_gam = Y − Cov(R, x) + Σρ Z", Ygam == Y - CRx + Zl)
    end

    # ===================================================================================
    say("\n(b) random pairing of the survivors' chromosomes")
    for (trial, n) in enumerate((2, 3, 4))
        d = [randq() for _ in 1:2n]; e = [randq() for _ in 1:2n]
        d .-= mean_(d); e .-= mean_(e)
        C = mean_(d .* e)
        Ms = matchings(collect(1:2n))
        Zs = [sum(d[a] * e[b] + d[b] * e[a] for (a, b) in m) / (2n) for m in Ms]
        EZ = mean_(Zs)
        check("n = $n: E[Z] = −Cov(R, x)/(2n − 1) over all $(length(Ms)) pairings (exact)", EZ == -C / (2n - 1))
    end
    # leading-order variance: 2n Var(Z) → v Var(R) + Cov(R, x)² (random populations of growing size)
    errs = Float64[]
    for n in (50, 200, 800)
        r = Xoshiro(7 + n)
        x = randn(r, 2n); R = 0.4 .* x .+ randn(r, 2n)
        x .-= mean_(x); R .-= mean_(R)
        v = mean_(x .^ 2); VR = mean_(R .^ 2); CRx = mean_(x .* R)
        # exact Var over random pairings, by the pair-moment formula: Var(Σ_p z_p) with z_p = x_a R_b + x_b R_a
        # computed by Monte Carlo over many pairings
        K = 20000; acc = 0.0; acc2 = 0.0
        for _ in 1:K
            pm = randperm(r, 2n)
            z = 0.0
            for i in 1:n
                a, b = pm[2i-1], pm[2i]
                z += x[a] * R[b] + x[b] * R[a]
            end
            z /= 2n
            acc += z; acc2 += z^2
        end
        varZ = acc2 / K - (acc / K)^2
        push!(errs, abs(2n * varZ / (v * VR + CRx^2) - 1))
    end
    check(@sprintf("2n Var(Z)/(v Var R + Cov²) → 1 (|error| %.3f, %.3f, %.3f at n = 50, 200, 800; Monte Carlo)", errs...),
          errs[end] < 0.05 && errs[end] < errs[1] + 0.02)

    # ===================================================================================
    say("\n(c) selection on individuals: between- and within-chromosome changes")
    @variables b q S ε
    for trial in 1:3
        K = 5
        # chromosome 1 values (B, x, G1 = B + x + rest), chromosome 2 independent: product measure
        B1 = [randq() for _ in 1:K]; X1 = [randq() for _ in 1:K]; O1 = [randq() for _ in 1:K]
        G1 = B1 .+ X1 .+ O1
        B2 = [randq() for _ in 1:K]; X2 = [randq() for _ in 1:K]; O2 = [randq() for _ in 1:K]
        G2 = B2 .+ X2 .+ O2
        pr = [(i, j) for i in 1:K for j in 1:K]           # individuals: (chromosome 1 = i, chromosome 2 = j), uniform
        Gm = mean_([G1[i] + G2[j] for (i, j) in pr]); Vg = mean_([(G1[i] + G2[j] - Gm)^2 for (i, j) in pr])
        h(i, j) = b * (G1[i] + G2[j] - Gm) - ((G1[i] + G2[j] - Gm)^2 - Vg) / (2S)
        wt(i, j) = 1 + ε * h(i, j)
        Wt = sum(wt(i, j) for (i, j) in pr)
        Ew(f) = sum(wt(i, j) * f(i, j) for (i, j) in pr) / Wt
        mB1 = Ew((i, j) -> B1[i]); mX2 = Ew((i, j) -> X2[j])
        covb = Ew((i, j) -> (B1[i] - mB1) * (X2[j] - mX2))
        first_ = coeffj(covb, ε, 1)
        cB = cov_(B1, G1); cx = cov_(X2, G2)
        check("trial $trial: d/dε Cov_w(B^(1), x^(2)) = −Cov(B,G)Cov(x,G)/S (b drops out; exact)", iszero_sym(first_ + cB * cx / S))
        check("trial $trial: Cov(B^(1), x^(2)) = 0 before selection (independent chromosomes)", iszero_sym(coeffj(covb, ε, 0)))
    end
    for trial in 1:3
        K = 7
        Bv = [randq() for _ in 1:K]; Xv = [randq() for _ in 1:K]; Ov = [randq() for _ in 1:K]
        Gv = Bv .+ Xv .+ Ov
        Gm = mean_(Gv); u = Gv .- Gm; Vv = mean_(u .^ 2)
        w = [1 + b * u[i] + q * (u[i]^2 - Vv) for i in 1:K]
        Wt = sum(w)
        Ew(f) = sum(w[i] * f(i) for i in 1:K) / Wt
        mB = Ew(i -> Bv[i]); mX = Ew(i -> Xv[i])
        covw = Ew(i -> (Bv[i] - mB) * (Xv[i] - mX))
        Y0 = cov_(Bv, Xv)
        cB = cov_(Bv, Gv); cx = cov_(Xv, Gv)
        dB = Bv .- mean_(Bv); dX = Xv .- mean_(Xv)
        k3 = mean_(dB .* dX .* u)
        k4 = mean_(dB .* dX .* u .^ 2) - Y0 * Vv - 2cB * cx
        # Taylor coefficients in (b, q)
        e = covw - Y0
        c_b = Symbolics.substitute(coeffj(e, b, 1), Dict(q => 0))
        c_q = Symbolics.substitute(coeffj(e, q, 1), Dict(b => 0))
        c_bb = Symbolics.substitute(coeffj(e, b, 2), Dict(q => 0))
        check("trial $trial: tilt 1 + bu + q(u² − V): ∂_b = κ(B,x,G), ∂_q = κ(B,x,G,G) + 2c_B c_x, ∂²_b/2 = −c_B c_x (exact)",
              iszero_sym(c_b - k3) && iszero_sym(c_q - (k4 + 2cB * cx)) && iszero_sym(c_bb + cB * cx))
    end
    # with q = (b² − 1/S)/2 (the Gaussian tilt e^{bu − u²/2S} to second order) the b² terms cancel
    @variables κ3 κ4 cc
    ΔY = b * κ3 + q * (κ4 + 2cc) - b^2 * cc
    check("q = (b² − 1/S)/2: ΔY = −c_B c_x/S + b κ3 + (b² − 1/S) κ4/2",
          iszero_sym(Symbolics.substitute(ΔY, Dict(q => (b^2 - 1 / S) / 2)) - (-cc / S + b * κ3 + (b^2 - 1 / S) * κ4 / 2)))

    # ===================================================================================
    say("\n(d) mutation")
    for trial in 1:2
        M = 4
        eB = [randq() for _ in 1:M]; dx = [randq() for _ in 1:M]
        eB .-= mean_(eB); dx .-= mean_(dx)
        sB = Rational{BigInt}(rand(rng, 1:3)); sl = Rational{BigInt}(rand(rng, 1:3))   # δ_B = ±sB, δ_l = ±sl
        Y0 = mean_(eB .* dx)
        tot = 0 // 1; cnt = 0
        for sgB in Iterators.product(fill((-1, 1), M)...), sgl in Iterators.product(fill((-1, 1), M)...)
            Bp = eB .+ sB .* collect(sgB); xp = dx .+ sl .* collect(sgl)
            tot += (cov_(Bp, xp) - Y0)^2; cnt += 1
        end
        EΔ2 = tot / cnt
        mB = sB^2; ml = sl^2; VB = mean_(eB .^ 2); v = mean_(dx .^ 2)
        check("trial $trial: E[(ΔŶ)²] = (Var(B) m_l + v m_B)/M + (M − 1) m_B m_l/M² (exact enumeration)",
              EΔ2 == (VB * ml + v * mB) / M + (M - 1) * mB * ml / M^2)
    end

    # ===================================================================================
    say("\n(e) the composite gamete map")
    @variables r1 r2 r3 vb c1 c2 c3 Sx
    ρv = [r1, r2, r3]; cv = [c1, c2, c3]
    Id = Matrix{Num}(I, 3, 3)
    A = Id .- (cv * ones(Num, 3)' .+ vb .* Id) ./ Sx
    T = Diagonal(1 .- ρv) * A .+ Diagonal(ρv) * (A .- Id)
    check("diag(1 − ρ) A + diag(ρ)(A − I) = diag(1 − ρ) + (A − I)", iszero_sym(T .- (Diagonal(1 .- ρv) .+ (A .- Id))))

    # ===================================================================================
    say("\n(g) heterogeneity among loci in Bulmer's step (`vhet`)")
    @variables v1 v2 v3 a11 a12 a21 a22 b1 b2 μ1 μ2 p11 p12 p22 q11 q12 q22 w11 w12 w22
    A2 = [a11 a12; a21 a22]; bb = [b1, b2]
    # the stationary mean of y_l ← A y_l + b v_l is (I − A)⁻¹ b v_l: proportional to v_l (exact rational A, b)
    Aq = Rational{BigInt}[9/10 -1/20; 1/30 4/5]; bq = Rational{BigInt}[-1/7, -2/9]
    ml(v) = inv(Matrix{Rational{BigInt}}(I, 2, 2) .- Aq) * (bq .* v)
    ms = [ml(v) for v in (v1, v2, v3)]
    check("stationary mean E[y_l] = (I − A)⁻¹ b v_l ∝ v_l, so Σ_l E[y_l] v_l = (Σv²/Σv) Σ_l E[y_l]",
          iszero_sym(sum(ms[i] .* [v1, v2, v3][i] for i in 1:3) .- ((v1^2 + v2^2 + v3^2) / (v1 + v2 + v3)) .* sum(ms)))
    # one step of the summed moments with E[y_l] = v_l μ and arbitrary per-locus second moments P_l
    μv = [μ1, μ2]
    Pl = [[p11 p12; p12 p22], [q11 q12; q12 q22], [w11 w12; w12 w22]]
    vs = [v1, v2, v3]
    exact = sum(A2 * Pl[i] * A2' .+ (A2 * (vs[i] .* μv)) * (vs[i] .* bb)' .+ (vs[i] .* bb) * (A2 * (vs[i] .* μv))' .+ (vs[i]^2) .* (bb * bb') for i in 1:3)
    P = sum(Pl); m = sum(vs[i] .* μv for i in 1:3)
    V = v1 + v2 + v3; V2s = v1^2 + v2^2 + v3^2
    withw(w) = A2 * P * A2' .+ (A2 * m) * (w .* bb)' .+ (w .* bb) * (A2 * m)' .+ (bb * bb') .* V2s
    check("Σ_l E[(A y_l + b v_l)(A y_l + b v_l)ᵀ] = A P Aᵀ + (Am)(w b)ᵀ + (w b)(Am)ᵀ + b bᵀ Σv² with w = Σv²/Σv (exact)",
          iszero_sym(exact .- withw(V2s / V)))
    check("  … and not with the old w = v̄ = Σv/n unless the v_l are equal (difference ∝ Σv² − (Σv)²/n)",
          !iszero_sym(exact .- withw(V / 3)) && iszero_sym(Symbolics.substitute.(exact .- withw(V / 3), Ref(Dict(v2 => v1, v3 => v1)))))
    check("  n Σv²/(Σv)² = 1 + CV² ≥ 1 (Cauchy–Schwarz form: n Σv² − (Σv)² = Σ_{i<j}(v_i − v_j)²)",
          iszero_sym(3V2s - V^2 - ((v1 - v2)^2 + (v1 - v3)^2 + (v2 - v3)^2)))

    # ===================================================================================
    say("\n(f) spot checks of the code")
    let L = 60, u = 3.0
        S = qtsim_survival(L, u)
        x = randn(rng, L); R = zeros(L)
        rsum!(R, x, S)
        check("rsum!: R_l = Σ_k |S_k − S_l| x_k (prefix sums)", maximum(abs.(R .- [sum(abs(S[k] - S[l]) * x[k] for k in 1:L) for l in 1:L])) < 1e-12)
        cs = class_seg(L, u)
        Bv = zeros(length(cs.cls)); Rv = zeros(length(cs.cls)); P0 = zeros(L + 1); P1 = zeros(L + 1)
        class_sums!(Bv, Rv, x, cs, S, P0, P1)
        G = ld_groups(L, u)
        ok = true
        for l in 1:L, (q, s) in enumerate(cs.cls[cs.coff[l]:cs.coff[l+1]-1])
            ks = [k for k in 1:L if k != l && (ld_rclass(abs(S[k] - S[l])) - 1) * LD_NB + ld_bin(k - 1, L) == s]
            ok &= abs(Bv[cs.coff[l]+q-1] - sum(x[ks])) < 1e-12 && abs(Rv[cs.coff[l]+q-1] - sum(abs(S[k] - S[l]) * x[k] for k in ks)) < 1e-12
        end
        check("class_seg/class_sums!: class sums and ρ-weighted class sums of every locus", ok)
        # gamete_pool against an explicit enumeration of the pool (QTSim's switch law)
        n = 5
        Bf = STBuf(L, n, u)
        Bf.X[:, 1:2n] .= randn(rng, L, 2n)
        for m in 1:2n
            Bf.G[m] = sum(view(Bf.X, :, m))
        end
        ms = moments(Bf, 1:n, ones(n))
        gp = gamete_pool(Bf, 1:n, ms)
        pc = [S[c] - S[c+1] for c in 1:L-1]; p0 = 1 - S[1]
        pool = Tuple{Float64,Vector{Float64}}[]
        for i in 1:n, st in 1:2
            a = Bf.X[:, 2(i-1)+st]; o = Bf.X[:, 2(i-1)+3-st]
            push!(pool, (p0 / (2n), a))
            for c in 1:L-1
                push!(pool, (pc[c] / (2n), [k <= c ? a[k] : o[k] for k in 1:L]))
            end
        end
        Wp = sum(first, pool)
        xm = [sum(w * x[k] for (w, x) in pool) / Wp for k in 1:L]
        Gm = sum(w * sum(x) for (w, x) in pool) / Wp
        Yp = [sum(w * (sum(x) - x[l] - (Gm - xm[l])) * (x[l] - xm[l]) for (w, x) in pool) / Wp for l in 1:L]
        check("gamete_pool: Y^Gp = Y^S − Cov(R,x) + Z equals the enumerated pool's Cov(B_l, x_l)", maximum(abs.(gp.Ygp .- Yp)) < 1e-10)
    end
    let ns = 9
        P = randn(rng, ns, ns); P = P * P'; m = randn(rng, ns); c = abs.(randn(rng, ns)); vbar = 0.3; V2 = 0.7; Sx = 4.0
        αA = 1 - vbar / Sx; Ad = Matrix(1.0I, ns, ns) .- (c * ones(ns)' .+ vbar .* Matrix(1.0I, ns, ns)) ./ Sx
        Dd = Ad .- Matrix(1.0I, ns, ns)
        β = -c ./ Sx
        dense = Ad * P * Dd' .+ (Ad * m) * (vbar .* β') .+ (vbar .* β) * (Dd * m)' .+ (β * β') .* V2
        r1 = crossmom((αA, -c ./ Sx), (-vbar / Sx, -c ./ Sx), β, β, m, P, vbar, V2)
        check("crossmom (rank-one maps) equals the dense Σ_l (A y + β v)(D y + β v)ᵀ", maximum(abs.(r1 .- dense)) < 1e-12)
    end
    let L = 120, u = 3.0
        G = ld_groups(L, u); ρ2 = class_rho2(L, u)
        ns = size(G.n, 2)
        r = Xoshiro(11)
        Cbb = randn(r, LD_NB, LD_NB); Cbb = 0.01 .* (Cbb * Cbb') .+ Diagonal(fill(0.05, LD_NB))
        d = LDDrivers(β = 0.01, Sc = 3.0, dS = 0.05, n = 800.0, N = 1000.0, Vm = 1e-3, V = fill(0.02, LD_NB), V2 = fill(1e-4, LD_NB),
                      k3 = zeros(LD_NB), g = fill(0.02, LD_NB), c = fill(0.03, LD_NB), Cbb = Cbb)
        st1 = LDState(G); y_stationary!(st1, G, d; maxit = 300)
        st2 = LSState(G)
        for a in 1:LD_NB
            st2.y.m[a] .= st1.m[a]; st2.y.P[a] .= st1.P[a]
        end
        st3 = deepcopy(st2); st3.H .= [0.3 .* p for p in st3.y.P]; st4 = deepcopy(st3)
        for a in 1:LD_NB
            G.nb[a] == 0 && continue
            y_step!(st1, G, d, a); yfix_step!(st2, G, d, a; ρ2)
            yH_step!(st3, G, d, a; hsel = :bulmer); yfix_step!(st4, G, d, a; ρ2, couple = true, hsel = :bulmer)
        end
        e1 = maximum(maximum(abs.(st1.P[a] .- st2.y.P[a])) / max(maximum(abs.(st1.P[a])), 1e-300) for a in 1:LD_NB if G.nb[a] > 0)
        e2 = maximum(maximum(abs.(st3.y.P[a] .- st4.y.P[a])) / max(maximum(abs.(st3.y.P[a])), 1e-300) for a in 1:LD_NB if G.nb[a] > 0)
        e3 = maximum(maximum(abs.(st3.H[a] .- st4.H[a])) / max(maximum(abs.(st3.H[a])), 1e-300) for a in 1:LD_NB if G.nb[a] > 0)
        check(@sprintf("yfix_step! with no correction = y_step! (rel. %.1e)", e1), e1 < 1e-12)
        stA = LDState(G); y_stationary!(stA, G, d; tol = 1e-12)
        stB = LSState(G); yfix_stationary!(stB, G, d; ρ2, tol = 1e-12)
        e4 = maximum(maximum(abs.(stA.P[a] .- stB.y.P[a])) / max(maximum(abs.(stA.P[a])), 1e-300) for a in 1:LD_NB if G.nb[a] > 0)
        check(@sprintf("yfix_stationary! (active classes only) = y_stationary! (full) for level B (rel. %.1e)", e4), e4 < 1e-6)
        check(@sprintf("yfix_step!(couple, hsel = :bulmer) with no correction = yH_step! (P %.1e, H %.1e)", e2, e3), e2 < 1e-12 && e3 < 1e-12)
        a = 8
        Σ = class_cov(G, d, a); CR = class_covR(G, d, a, ρ2, Σ)
        ρ = G.ρ[a, :]
        check("class_covR: off-diagonal ρ_s ρ_t Σ_B[s,t]", all(isapprox(CR[s, t], ρ[s] * ρ[t] * Σ[s, t]; atol = 1e-15) for s in 1:ns, t in 1:ns if s != t))
        okd = all(begin
            bb = G.bofs[s]; fs = G.n[a, s] / max(G.nb[bb], 1)
            isapprox(CR[s, s], ρ[s]^2 * Σ[s, s] + (ρ2[a, s] - ρ[s]^2) * fs * d.g[bb]; atol = 1e-15)
        end for s in 1:ns)
        check("class_covR: diagonal adds the within-class spread (ρ2 − ρ̄²) f_s g_b ≥ 0", okd && all(ρ2[a, s] >= ρ[s]^2 - 1e-15 for s in 1:ns))
    end

    finish(LOG, t0)
end

runmain(main, @__FILE__)
