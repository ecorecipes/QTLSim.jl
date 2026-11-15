# Symbolic checks of the source of the association between a locus's squared deviation and its
# linked background (validate/linked_source.jl, validate/linked_source_theory.jl), and numerical
# spot checks of its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_linked_source.jl
# Writes validate/symbolic/check_linked_source.log; exits 1 if any check fails.
#
# Notation (chromosome level, centred moments): A = B_s (a class of partner loci), B the whole
# background of locus l, X = x_l, d_l = x_l − x̄_l, y_A = Cov(A, X), Y = Cov(B, X), v = Var(X),
# H = Cov(AB, d_l²) = μ(A,B,X,X) − Cov(A,B)v (variance heterogeneity), K = κ(A,B,X,X) (the source).
#
# (a) the exact decomposition: K = H − 2y_A Y; with A, B sums of loci, H = H2 + H3 and K = (two-locus
#     cumulants) + (three-locus cumulants), where H2 = Σ_k Cov(d_k², d_l²) and y_A Y = Σ_k D_kl² + cross;
#     a biallelic X has d² − v = s d (s = κ3/v), so H = s κ(A,B,X) and T = κ(A,X,X) = s y_A;
# (b) iid sampling of M (offspring): E[Ĥ] = (1 − 5/M)H + 4 y_A Y/M + O(1/M²) (exact rational form
#     checked), E[ŷ_A Ŷ] = (1 − 3/M) y_A Y + (H + Cov(A,B) v)/M + O(1/M²), E[K̂] = (1 − 7/M)K −
#     (2/M)(Cov(A,B)v + 2y_A Y) + O(1/M²); two loci (A = x_k): E[Ĥ] = (1 − 5/M)H + 4y²/M — the allele-value
#     form of Ohta & Kimura's E[D(1−2p)(1−2q)] equation; biallelic loci: exactly theirs (check on
#     haplotype frequencies);
# (c) survival as a subset of m of M: to leading order iid with 1/M → 1/m − 1/M (exact enumeration,
#     the error falls as the sizes grow);
# (d) single-crossover recombination of three loci: h = μ(k,j,l,l) − D_kj v_l ← p₀ h with p₀ the
#     probability of no split among them (= 1 − max(ρ_kl, ρ_jl) on one side, 1 − ρ_kl − ρ_jl on both);
#     two loci H ← (1 − ρ)H;
# (e) mutation on a random one of M chromosomes: E[ΔH] = 0 exactly;
# (f) selection: under weights 1 + εh, dH/dε = κ₅(A,B,X,X,h) + 2[κ(A,X,h)Y + y_A κ(B,X,h)] exactly;
#     biallelic: Δs = −2β(v + Y) at first order (every locus's skew falls under directional selection);
# (g) the neutral fixed point of one pair: H/(2y²) → 2ε/(5ε + ρ) (one stage), ≤ 2/5;
#     `pair_fixed_point` solves the two-stage map exactly;
# (h) spot checks of validate/linked_source_theory.jl: sample_stage! against (b), span_matrix,
#     yH_step! without selection against pair_fixed_point (one class).

include("common.jl")
include(joinpath(@__DIR__, "..", "linked_decay_theory.jl"))
include(joinpath(@__DIR__, "..", "linked_source_theory.jl"))

const LOG = joinpath(@__DIR__, "check_linked_source.log")
const rng = Xoshiro(20261115)

# ---------------------------------------------------------------------------------------
# expectations of products of sample means under M iid draws (as check_linked_decay.jl, for any
# number of variables NV)

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
# exact rational populations

randq() = Rational{BigInt}(rand(rng, -9:9), rand(rng, 1:5))
wmean(w, f) = sum(w .* f) / sum(w)
cm(w, vs...) = (c = [v .- wmean(w, v) for v in vs]; wmean(w, reduce(.*, c)))      # central mixed moment
function kap(w, vs...)
    n = length(vs)
    c = [v .- wmean(w, v) for v in vs]
    E(ix...) = wmean(w, reduce(.*, (c[i] for i in ix)))
    n == 2 && return E(1, 2)
    n == 3 && return E(1, 2, 3)
    n == 4 && return E(1, 2, 3, 4) - E(1, 2) * E(3, 4) - E(1, 3) * E(2, 4) - E(1, 4) * E(2, 3)
    if n == 5
        s = E(1, 2, 3, 4, 5)
        for pr in ((1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5))
            s -= E(pr...) * E(Tuple(setdiff(1:5, pr))...)
        end
        return s
    end
    error("order")
end
Hof(w, A, B, X) = cm(w, A, B, X, X) - cm(w, A, B) * cm(w, X, X)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_linked_source.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) exact decomposition of K")
    for trial in 1:3
        K = 9
        w = [Rational{BigInt}(rand(rng, 1:4)) for _ in 1:K]
        xs = [[randq() for _ in 1:K] for _ in 1:4]     # loci k1, k2 (class s), j (rest), l
        A = xs[1] .+ xs[2]; Bv = A .+ xs[3]; X = xs[4]
        yA = cm(w, A, X); Y = cm(w, Bv, X)
        H = Hof(w, A, Bv, X)
        check("trial $trial: K = κ(A,B,X,X) = H − 2 y_A Y", kap(w, A, Bv, X, X) == H - 2yA * Y)
        H2 = sum(cm(w, xs[k], xs[k], X, X) - cm(w, xs[k], xs[k]) * cm(w, X, X) for k in 1:2)
        H3 = H - H2
        H3direct = sum(cm(w, xs[k], xs[j], X, X) - cm(w, xs[k], xs[j]) * cm(w, X, X) for k in 1:2 for j in 1:3 if j != k)
        check("trial $trial: H = H2 + H3, H3 = Σ_{k∈s} Σ_{j≠k} Cov(d_k d_j, d_l²)", H3 == H3direct)
        Dsq = sum(cm(w, xs[k], X)^2 for k in 1:2)
        cross = sum(cm(w, xs[k], X) * cm(w, xs[j], X) for k in 1:2 for j in 1:3 if j != k)
        check("trial $trial: y_A Y = Σ_k D_kl² + Σ_{k≠j} D_kl D_jl", yA * Y == Dsq + cross)
        K2 = sum(kap(w, xs[k], xs[k], X, X) for k in 1:2)
        K3 = sum(kap(w, xs[k], xs[j], X, X) for k in 1:2 for j in 1:3 if j != k)
        check("trial $trial: K = two-locus + three-locus cumulants; two-locus = H2 − 2ΣD²", K2 + K3 == kap(w, A, Bv, X, X) && K2 == H2 - 2Dsq)
    end
    for trial in 1:3
        K = 8
        w = [Rational{BigInt}(rand(rng, 1:4)) for _ in 1:K]
        a = randq(); a == 0 && (a = 1 // 1)
        X = [rand(rng, Bool) ? a : 0 // 1 for _ in 1:K]
        X[1] = a; X[2] = 0
        A = [randq() for _ in 1:K]; Bv = A .+ [randq() for _ in 1:K]
        v = cm(w, X, X); s = cm(w, X, X, X) / v
        d = X .- wmean(w, X)
        check("trial $trial: biallelic X: d² − v = s d exactly", all(d .^ 2 .- v .== s .* d))
        check("trial $trial: biallelic X: H = s κ(A,B,X)", Hof(w, A, Bv, X) == s * cm(w, A, Bv, X))
        check("trial $trial: biallelic X: T = κ(A,X,X) = s y_A", cm(w, A, X, X) == s * cm(w, A, X))
    end

    # ===================================================================================
    say("\n(b) iid sampling of M (offspring)")
    @variables M ε
    μs = Dict{Mon,Any}()
    μ(e) = sum(e) == 0 ? 1 : sum(e) == 1 ? 0 : get!(μs, e) do
        Symbolics.variable(Symbol("μ", join(e)))
    end
    Hp = padd(central((1, 1, 2)), pmul(central((1, 1, 0)), central((0, 0, 2))), -1)
    yy = pmul(central((1, 0, 1)), central((0, 1, 1)))
    Hm = μ((1, 1, 2)) - μ((1, 1, 0)) * μ((0, 0, 2)); yYm = μ((1, 0, 1)) * μ((0, 1, 1)); Cv = μ((1, 1, 0)) * μ((0, 0, 2))
    EH = expect(Hp, M, μ); EyY = expect(yy, M, μ)
    exactH = ((M^3 - 5M^2 + 8M - 4) * μ((1, 1, 2)) - (M^3 - 5M^2 + 8M - 4) * Cv + (4M^2 - 12M + 8) * yYm) / M^3
    check("E[Ĥ] = (M−1)(M−2)²/M³ · H + 4(M−1)(M−2) y_A Y/M³ (exact)", iszero_sym(EH - exactH))
    check("  that is (M−1)(M−2)[(M−2)H + 4 y_A Y]/M³", iszero_sym(exactH - (M - 1) * (M - 2) * ((M - 2) * Hm + 4yYm) / M^3))
    lead(e) = simp(Symbolics.substitute(e, Dict(M => 1 / ε)))
    check("E[Ĥ] = (1 − 5/M) H + 4 y_A Y/M + O(1/M²)", iszero_sym(coeffj(lead(EH), ε, 1) - (-5Hm + 4yYm)) &&
                                                          iszero_sym(coeffj(lead(EH), ε, 0) - Hm))
    check("E[ŷ_A Ŷ] = (1 − 3/M) y_A Y + (H + Cov(A,B) v)/M + O(1/M²)",
          iszero_sym(coeffj(lead(EyY), ε, 1) - (-3yYm + Hm + Cv)) && iszero_sym(coeffj(lead(EyY), ε, 0) - yYm))
    EK = simp(EH - 2EyY)
    Km = Hm - 2yYm
    check("E[K̂] = (1 − 7/M) K − (2/M)(Cov(A,B) v + 2 y_A Y) + O(1/M²)",
          iszero_sym(coeffj(lead(EK), ε, 1) - (-7Km - 2(Cv + 2yYm))))
    # two loci: A = B = x_k: variables (k, l) via (1, 0, ·) with B ≡ A: set μ(e) for e = (i, j, m) to ν(i+j, m)
    νs = Dict{Tuple{Int,Int},Any}()
    ν(i, j) = i + j == 0 ? 1 : i + j == 1 ? 0 : get!(νs, (i, j)) do
        Symbolics.variable(Symbol("ν", i, j))
    end
    μ2(e) = ν(e[1] + e[2], e[3])
    EH2 = expect(Hp, M, μ2); Ey2 = expect(yy, M, μ2)
    H2m = ν(2, 2) - ν(2, 0) * ν(0, 2); y2m = ν(1, 1)^2
    check("two loci: E[Ĥ] = (1 − 5/M)H + 4y²/M + O(1/M²)", iszero_sym(coeffj(lead(EH2), ε, 1) - (-5H2m + 4y2m)))
    check("two loci: E[ŷ²] = (1 − 3/M)y² + (H + v_k v_l)/M + O(1/M²)",
          iszero_sym(coeffj(lead(Ey2), ε, 1) - (-3y2m + H2m + ν(2, 0) * ν(0, 2))))
    # biallelic loci (indicators I, J; frequencies p, q, LD D): central moments as polynomials
    @variables p q D
    function bial(i, j)
        # E[(I − p)^i (J − q)^j] over haplotypes AB (p q + D), Ab (p(1−q) − D), aB ((1−p)q − D), ab
        f = [(1 - p, 1 - q, p * q + D), (1 - p, -q, p * (1 - q) - D), (-p, 1 - q, (1 - p) * q - D), (-p, -q, (1 - p) * (1 - q) + D)]
        sum(x^i * y^j * w for (x, y, w) in f)
    end
    μb(e) = bial(e[1] + e[2], e[3])
    check("biallelic: H = (1 − 2p)(1 − 2q) D (Ohta & Kimura's moment)", iszero_sym(bial(2, 2) - bial(2, 0) * bial(0, 2) - (1 - 2p) * (1 - 2q) * D))
    EHb = expect(Hp, M, μb); Eyb = expect(yy, M, μb)
    check("biallelic: E[D̂(1−2p̂)(1−2q̂)] = (1 − 5/M)D(1−2p)(1−2q) + 4D²/M + O(1/M²)",
          iszero_sym(coeffj(lead(EHb), ε, 1) - (-5(1 - 2p) * (1 - 2q) * D + 4D^2)))
    check("biallelic: E[D̂²] = (1 − 3/M)D² + (D(1−2p)(1−2q) + p(1−p)q(1−q))/M + O(1/M²)",
          iszero_sym(coeffj(lead(Eyb), ε, 1) - (-3D^2 + (1 - 2p) * (1 - 2q) * D + p * (1 - p) * q * (1 - q))))

    # ===================================================================================
    say("\n(c) survival: a subset of m of M, to leading order iid with 1/m − 1/M")
    errs = Float64[]
    for (Mp, mp) in ((8, 4), (12, 6), (16, 8))
        # a population with LD and heterogeneity: X depends on A
        A = [randq() for _ in 1:Mp]; Bv = A .+ [randq() for _ in 1:Mp]
        X = A .* [randq() for _ in 1:Mp] .+ [randq() for _ in 1:Mp]
        w1 = ones(Rational{BigInt}, Mp)
        H0 = Hof(w1, A, Bv, X); yY0 = cm(w1, A, X) * cm(w1, Bv, X)
        tot = 0 // 1; nsub = 0
        for s in 0:2^Mp-1
            count_ones(s) == mp || continue
            ix = findall(digits(s; base = 2, pad = Mp) .== 1)
            tot += Hof(ones(Rational{BigInt}, mp), A[ix], Bv[ix], X[ix]); nsub += 1
        end
        e = 1 // mp - 1 // Mp
        pred = H0 * (1 - 5e) + 4e * yY0
        push!(errs, Float64(abs(tot / nsub - pred) / (abs(H0) + abs(yY0))))
    end
    check("subset error falls with size (relative error $(round.(errs; sigdigits = 2)))", errs[3] < errs[2] < errs[1] && errs[3] < 0.05)

    # ===================================================================================
    say("\n(d) single-crossover recombination of three loci")
    @variables p0 p1 p2 p3 μkjll Dkj vl
    # gamete from two independent homologs; partitions: none (p0), l | kj (p1), k | jl (p2), j | kl (p3)
    μn = p0 * μkjll + p1 * Dkj * vl + p2 * 0 + p3 * 0
    Dn = (p0 + p1) * Dkj
    check("h′ = μ′ − D′_kj v_l = p₀ (μ − D_kj v_l)", iszero_sym(μn - Dn * vl - p0 * (μkjll - Dkj * vl)))
    for (L, u) in ((9, 3.0), (20, 5.0))
        S = qtsim_survival(L, u)
        ρ(i, j) = abs(S[i+1] - S[j+1])
        ok = true
        for (k, j, l) in ((0, 2, 5), (6, 1, 3), (2, 3, 8))
            lo, hi = minimum((k, j, l)), maximum((k, j, l))
            p0v = 1 - (S[lo+1] - S[hi+1])           # no crossover among them
            same = (k < l) == (j < l)
            ok &= isapprox(1 - p0v, same ? max(ρ(k, l), ρ(j, l)) : ρ(k, l) + ρ(j, l); atol = 1e-14)
        end
        check("L = $L, u = $u: P(split among k, j, l) = max(ρ_kl, ρ_jl) on one side, ρ_kl + ρ_jl across", ok)
    end

    # ===================================================================================
    say("\n(e) mutation on a random one of M chromosomes")
    for trial in 1:2
        Mp = 6
        A = [randq() for _ in 1:Mp]; Bv = A .+ [randq() for _ in 1:Mp]; X = [randq() for _ in 1:Mp]
        w1 = ones(Rational{BigInt}, Mp)
        H0 = Hof(w1, A, Bv, X)
        @variables ee α2
        EHm = 0
        for i in 1:Mp
            x2 = [j == i ? X[j] + ee : X[j] + 0 * ee for j in 1:Mp]
            cx = x2 .- sum(x2) / Mp; cA = A .- sum(A) / Mp; cB = Bv .- sum(Bv) / Mp
            Hm_ = sum(cA .* cB .* cx .^ 2) / Mp - (sum(cA .* cB) / Mp) * (sum(cx .^ 2) / Mp)
            Hm_ = Symbolics.expand(Hm_)
            EHm += (Symbolics.substitute(Hm_, Dict(ee => 0)) + coeffj(Hm_, ee, 2) * α2) / Mp
        end
        # the ε² parts of μ(A,B,X,X) and of Cov(A,B)·v both grow by α²Cov(A,B)(M − 1)/M²
        check("trial $trial: E[ΔH] = 0 exactly", iszero_sym(EHm - H0))
    end

    # ===================================================================================
    say("\n(f) selection")
    for trial in 1:4
        K = 7
        w = [Rational{BigInt}(rand(rng, 1:6)) for _ in 1:K]
        A = [randq() for _ in 1:K]; Bv = A .+ [randq() for _ in 1:K]; X = [randq() for _ in 1:K]
        h = [randq() for _ in 1:K]
        @variables e
        we = w .* (1 .+ e .* h)
        cen(v) = v .- sum(we .* v) / sum(we)
        Ew(f) = sum(we .* f) / sum(we)
        He = Ew(cen(A) .* cen(Bv) .* cen(X) .^ 2) - Ew(cen(A) .* cen(Bv)) * Ew(cen(X) .^ 2)
        dH = Symbolics.substitute(Symbolics.expand_derivatives(Differential(e)(He)), Dict(e => 0))
        val(x) = Rational{BigInt}(Symbolics.value(simp(x)))
        rhs = kap(w, A, Bv, X, X, h) + 2 * (kap(w, A, X, h) * kap(w, Bv, X) + kap(w, A, X) * kap(w, Bv, X, h))
        check("trial $trial: dH/dε = κ₅(A,B,X,X,h) + 2[κ(A,X,h)Y + y_A κ(B,X,h)] exactly", val(dH) == rhs)
    end
    @variables a β Dl qf
    # biallelic: Δq = β Cov(1_A, u), u = a(1_A − q) + B: Cov = a q(1−q) + Cov(1_A, B) =: a q(1−q) + Dl
    dq = β * (a * qf * (1 - qf) + Dl)
    s0 = a * (1 - 2qf); s1 = a * (1 - 2(qf + dq))
    v0 = a^2 * qf * (1 - qf); Yl = a * Dl
    check("biallelic: Δs = −2β(v + Y) at first order", iszero_sym(s1 - s0 + 2β * (v0 + Yl)))

    # ===================================================================================
    say("\n(g) neutral fixed point of one pair")
    @variables ρ e1
    # one stage, then recombination: P ← (1−ρ)²[(1−3ε)P + ε(H + Cv)], H ← (1−ρ)[(1−5ε)H + 4εP]
    @variables Pp Hh Cvv
    # leading order (ρ = r e1, ε = e e1): fixed point H/P = 4ε/(5ε + ρ)
    @variables r ee2
    Hfix = 4ee2 * Pp / (5ee2 + r)
    resid = (1 - r * e1) * ((1 - 5ee2 * e1) * Hfix + 4ee2 * e1 * Pp) - Hfix
    check("H fixed point: H/y² = 4ε/(5ε + ρ) to leading order, i.e. H/(2y²) = 2ε/(5ε + ρ)", iszero_sym(coeffj(resid, e1, 1)))
    check("  at ρ = 0 the ratio H/(2y²) is 2/5", iszero_sym(Symbolics.substitute(Hfix / (2Pp), Dict(r => 0)) - 2 // 5))
    fp = pair_fixed_point(0.0, 1e-4, 1e-4, 1.0)
    check("pair_fixed_point at ρ = 0, small ε: H/(2P) ≈ 0.4", isapprox(fp.H / (2fp.P), 0.4; rtol = 1e-3))
    fp2 = pair_fixed_point(0.05, 2.5e-4, 2.5e-4, 1.0)
    check("pair_fixed_point at ρ = 0.05: H/(2P) ≈ 2ε/(5ε+ρ) with ε = 5e-4 per generation",
          isapprox(fp2.H / (2fp2.P), 2 * 5e-4 / (5 * 5e-4 + 0.05); rtol = 0.05))

    # ===================================================================================
    say("\n(h) spot checks of validate/linked_source_theory.jl")
    G1 = LDGroups(10, fill(9.0, 1, 1), fill(0.0, 1, 1), [10], [1])
    st = LSState(G1)
    st.y.P[1] .= 2.0; st.H[1] .= 0.7; st.y.m[1] .= 0.3
    Σ = fill(1.5, 1, 1)
    sample_stage!(st, Σ, 0.2, 0.01, 1)
    check("sample_stage!: P ← (1 − 3ε)P + ε(H + ΣV)", isapprox(st.y.P[1][1], 0.97 * 2 + 0.01 * (0.7 + 0.3)))
    check("sample_stage!: H ← (1 − 5ε)H + 4εP", isapprox(st.H[1][1], 0.95 * 0.7 + 0.04 * 2))
    # span matrix of 3 bins
    Gs = LDGroups(30, fill(1.0, 3, 3), [0.1 0.2 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3], [10, 10, 10], [1, 2, 3])
    sp = span_matrix(Gs, 2)
    check("span_matrix: opposite sides add, same side max", sp[1, 3] ≈ 0.4 && sp[1, 1] ≈ 0.1 && sp[1, 2] ≈ 0.2 && sp[3, 2] ≈ 0.3)
    # yH_step! without selection, one class, Σ = V·g: stationary = pair_fixed_point (Vm = 0)
    Gp = LDGroups(2, fill(1.0, 1, 1), fill(0.02, 1, 1), [1], [1])
    d = LDDrivers(β = 0.0, Sc = 50.0, dS = 0.0, n = 500.0, N = 1000.0, Vm = 0.0, V = [0.5], V2 = [0.25], k3 = [0.0],
                  g = [0.5], c = [0.5], Cbb = fill(0.5, 1, 1))
    st2 = LSState(Gp)
    for t in 1:200000
        yH_step!(st2, Gp, d, 1; select = false)
    end
    Σp = class_cov(Gp, d, 1)
    fp3 = pair_fixed_point(0.02, 1 / 1000 - 1 / 2000, 1 / 2000, Σp[1, 1] * 0.5)
    check("yH_step!(select = false) stationary P, H = pair_fixed_point", isapprox(st2.y.P[1][1], fp3.P; rtol = 1e-6) &&
                                                                         isapprox(st2.H[1][1], fp3.H; rtol = 1e-6))
    check("source_K = Σ_t (H − 2P)", source_K(st2, 1)[1] ≈ st2.H[1][1] - 2st2.y.P[1][1])

    finish(LOG, t0)
end

runmain(main, @__FILE__)
