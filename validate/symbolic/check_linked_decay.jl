# Symbolic checks of the recursion for the association between a locus's squared deviation and
# its linked background (validate/linked_decay.jl, validate/linked_decay_theory.jl), and
# numerical spot checks of its code.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_linked_decay.jl
# Writes validate/symbolic/check_linked_decay.log; exits 1 if any check fails.
#
# Notation (chromosome level): x_k allele value at locus k, d_l = x_l − x̄_l, B_s a background sum,
# T = κ(B_s, l, l) (third order), y = Cov(B_s, x_l) (second), K = κ(B_s, B, l, l) (fourth).
#
# (a) sampling (offspring = M iid draws from the parents' pool): E[m₁₁] = (1 − 1/M)μ₁₁, every
#     third-order central moment E[m] = (1 − 1/M)(1 − 2/M)μ (so drift multiplies T by ≈ 1 − 3/M and
#     y by 1 − 1/M), Var(m₁₁) = (μ₂₂ − μ₁₁²)/M + O(1/M²); with x_l independent of (B_s, B):
#     E[m(B_s,B,l,l) − m(B_s,B)m(l,l)] = 0 exactly and E[m(B_s,l)m(B,l)] = (M − 1)Cov(B_s,B)v_l/M²,
#     so the sample's K is −2·(the sample's y_s Y) in expectation: drift seeds K only through the
#     squared allele–background LD (variance heterogeneity and three-locus parts are not seeded);
# (b) survival as a subset of m of M (no replacement): E[m₃] factor ((m−1)(m−2)/m²)(M²/((M−1)(M−2)))
#     = 1 − 3/m + 3/M + …, E[m₂] factor (m−1)M/(m(M−1)); exact enumeration of all subsets;
# (c) selection: under weights w = 1 + εh the first-order change of any joint cumulant is
#     κ(A, h) (exact derivative); Leonov–Shiryaev for h containing u² (u centred):
#     κ(k,l,l,u²) = κ₅(k,l,l,u,u) + 2[κ(k,u)κ(l,l,u) + 2κ(l,u)κ(k,l,u)],
#     κ(k,l,u²) = κ₄(k,l,u,u) + 2κ(k,u)κ(l,u); hence for h = a u − u²/(2S_c), to first order,
#     ΔT = a κ(k,l,l,u) − [κ(k,u)κ(l,l,u) + 2κ(l,u)κ(k,l,u)]/S_c − κ₅/(2S_c),
#     Δy = a κ(k,l,u) − κ(k,u)κ(l,u)/S_c − κ₄(k,l,u,u)/(2S_c);
#     chromosome-level weight E[π | G₁] for G₂ ~ N(Ḡ₂, V₂): width S_c = ω² + V_E + V₂ (exact Gaussian integral);
# (d) recombination (gamete from two independent homologs, split with probability ρ):
#     T′ = (1 − ρ)T, y′ = (1 − ρ)y, κ₂₂′ = (1 − ρ)κ₂₂ + 2ρ(1 − ρ)y² (why y² is carried separately);
# (e) mutation (effect ε, mean 0, variance α², on a random one of M chromosomes): E[ΔT] = 0,
#     E[Δy] = 0, E[(Δy)²] = α² Var(B)/M²;
# (f) one biallelic locus: κ(B, x_l, x_l) = a²(1 − 2p)Cov(B, 1_A), and Δ(a²pq) for Δp = βD is
#     β κ(B,l,l) − β² Cov(B, x_l)²: the first-order loss is βC, its autocovariance part the
#     linked change of p;
# (g) the recursion's closed forms: T_t = βKλ(1 − λᵗ)/(1 − λ) from T₀ = 0; stationary
#     E[y²] = ν/(1 − λ_y²); the steady-state loss per locus, −β² v V_chr Σω/(n(ρ+δ_y)(ρ+δ_T)) to
#     leading order, is Barton's (1995) excess drift loss v b² V_G Σω/(2n(ρ+δ)²) when δ_y = δ_T = δ
#     and V_chr = V_G/2;
# (h) spot checks of validate/linked_decay_theory.jl: qtsim_survival against QTLSim.pair_splits,
#     ld_groups, the rank-one form of A P Aᵀ in y_step!, T_step! against the closed form.

include("common.jl")
include(joinpath(@__DIR__, "..", "linked_decay_theory.jl"))

const LOG = joinpath(@__DIR__, "check_linked_decay.log")
const rng = Xoshiro(20261114)

# ---------------------------------------------------------------------------------------
# (a) expectations of products of sample means under M iid draws: a polynomial in sample
# means A(x^e) (e an exponent tuple over nv variables) is a Dict(sorted factor list => coefficient)

const Poly = Dict{Vector{NTuple{3,Int}},Any}
Amon(e::NTuple{3,Int}) = Poly([e] => 1)
pone() = Poly(NTuple{3,Int}[] => 1)
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

# central sample moment m_e = A(Π (x_i − X̄_i)^{e_i})
function central(e::NTuple{3,Int})
    acc = Poly()
    for i in 0:e[1], j in 0:e[2], k in 0:e[3]
        c = binomial(e[1], i) * binomial(e[2], j) * binomial(e[3], k) * (-1)^(e[1] - i + e[2] - j + e[3] - k)
        term = pscale(Amon((i, j, k)), c)
        term = pmul(term, ppow(Amon((1, 0, 0)), e[1] - i))
        term = pmul(term, ppow(Amon((0, 1, 0)), e[2] - j))
        term = pmul(term, ppow(Amon((0, 0, 1)), e[3] - k))
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

# E[poly] for M iid draws from a population with central moments μ(e) (means zero)
function expect(p::Poly, M, μ)
    tot = 0
    for (fac, c) in p
        k = length(fac)
        s = 0
        for π in setpartitions(k)
            term = falling(M, length(π)) / M^k
            for B in π
                e = reduce((x, y) -> x .+ y, fac[B])
                term *= μ(e)
            end
            s += term
        end
        tot += c * s
    end
    tot
end

# ---------------------------------------------------------------------------------------
# (b), (c), (e) exact rational populations

randq() = Rational{BigInt}(rand(rng, -9:9), rand(rng, 1:5))
wmean(w, f) = sum(w .* f) / sum(w)
# joint cumulants of centred variables (weights w): orders 2–5
function kap(w, vs...)
    c = [v .- wmean(w, v) for v in vs]
    E(ix...) = wmean(w, reduce(.*, (c[i] for i in ix)))
    n = length(vs)
    n == 2 && return E(1, 2)
    n == 3 && return E(1, 2, 3)
    n == 4 && return E(1, 2, 3, 4) - E(1, 2) * E(3, 4) - E(1, 3) * E(2, 4) - E(1, 4) * E(2, 3)
    if n == 5
        s = E(1, 2, 3, 4, 5)
        for pr in ((1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5))
            rest = Tuple(setdiff(1:5, pr))
            s -= E(pr...) * E(rest...)
        end
        return s
    end
    error("order")
end

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_linked_decay.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) offspring sampling: M iid draws from the parents' pool")
    @variables M
    μs = Dict{NTuple{3,Int},Any}()
    μ(e) = e == (0, 0, 0) ? 1 : sum(e) == 1 ? 0 : get!(μs, e) do
        Symbolics.variable(Symbol("μ", join(e)))
    end
    check("E[m₁₁] = (1 − 1/M) μ₁₁", iszero_sym(expect(central((1, 1, 0)), M, μ) - (1 - 1 / M) * μ((1, 1, 0))))
    f3 = (1 - 1 / M) * (1 - 2 / M)
    for e in ((1, 2, 0), (2, 1, 0), (3, 0, 0), (1, 1, 1))
        check("E[m_$(join(e))] = (1 − 1/M)(1 − 2/M) μ_$(join(e)) (third order)",
              iszero_sym(expect(central(e), M, μ) - f3 * μ(e)))
    end
    m11 = central((1, 1, 0))
    varm11 = expect(pmul(m11, m11), M, μ) - expect(m11, M, μ)^2
    @variables ε
    lead = simp(Symbolics.substitute(varm11 / ε, Dict(M => 1 / ε)))
    check("Var(m₁₁) = (μ₂₂ − μ₁₁²)/M + O(1/M²)",
          iszero_sym(Symbolics.substitute(lead, Dict(ε => 0)) - (μ((2, 2, 0)) - μ((1, 1, 0))^2)))
    # x_l (variable 3) independent of (B_s, B) (variables 1, 2): its mixed central moments factor
    μind(e) = e == (0, 0, 0) ? 1 : (e[3] == 1 || (e[1] + e[2] == 1)) ? 0 :
              (e[1] + e[2] == 0 ? μ((0, 0, e[3])) : e[3] == 0 ? μ((e[1], e[2], 0)) : μ((e[1], e[2], 0)) * μ((0, 0, e[3])))
    mBBll = central((1, 1, 2)); mBB = central((1, 1, 0)); mll = central((0, 0, 2))
    mBsl = central((1, 0, 1)); mBl = central((0, 1, 1))
    check("independent x_l: E[m(B_s,B,l,l) − m(B_s,B)m(l,l)] = 0 exactly (no variance heterogeneity seeded)",
          iszero_sym(expect(padd(mBBll, pmul(mBB, mll), -1), M, μind)))
    check("independent x_l: E[m(B_s,l) m(B,l)] = (M − 1)Cov(B_s,B) v_l/M²",
          iszero_sym(expect(pmul(mBsl, mBl), M, μind) - (M - 1) / M^2 * μ((1, 1, 0)) * μ((0, 0, 2))))
    Khat = padd(padd(mBBll, pmul(mBB, mll), -1), pmul(mBsl, mBl), -2)
    check("independent x_l: E[K̂] = −2 E[ŷ_s Ŷ] = −2(M − 1)Cov(B_s,B)v_l/M²",
          iszero_sym(expect(Khat, M, μind) + 2 * (M - 1) / M^2 * μ((1, 1, 0)) * μ((0, 0, 2))))

    # ===================================================================================
    say("\n(b) survival as a subset of m of M without replacement (exact enumeration)")
    for (Mp, mp) in ((7, 4), (8, 3))
        X = [randq() for _ in 1:Mp]; Y = [randq() for _ in 1:Mp]
        w1 = ones(Rational{BigInt}, Mp)
        μ12 = kap(w1, X, Y, Y); μ11 = kap(w1, X, Y)
        subs = [findall(digits(s; base = 2, pad = Mp) .== 1) for s in 0:2^Mp-1]
        subs = filter(s -> length(s) == mp, subs)
        e3 = sum(kap(ones(Rational{BigInt}, mp), X[s], Y[s], Y[s]) for s in subs) / length(subs)
        e2 = sum(kap(ones(Rational{BigInt}, mp), X[s], Y[s]) for s in subs) / length(subs)
        check("M = $Mp, m = $mp: E[m₁₂] = ((m−1)(m−2)/m²)(M²/((M−1)(M−2))) μ₁₂ exactly",
              e3 == (mp - 1) * (mp - 2) // mp^2 * (Mp^2 // ((Mp - 1) * (Mp - 2))) * μ12)
        check("M = $Mp, m = $mp: E[m₁₁] = ((m−1)/m)(M/(M−1)) μ₁₁ exactly", e2 == (mp - 1) // mp * (Mp // (Mp - 1)) * μ11)
    end
    # m = 1/ε, M = 1/ε²: ((m−1)(m−2)/m²)(M²/((M−1)(M−2))) written in ε
    fsub = (1 - ε) * (1 - 2ε) / ((1 - ε^2) * (1 - 2ε^2))
    check("subset factor = 1 − 3/m + O(1/m², 1/M) (m = 2n survivors' genes)",
          iszero_sym(coeffj(fsub, ε, 1) + 3))
    @variables n N
    check("drift_f3 = 1 − 3/(2n) + O(1/n², 1/N)",
          iszero_sym(coeffj(simp(drift_f3(1 / ε, 1 / ε^2)), ε, 1) + 3 // 2))
    check("drift_f1 = 1 − 1/(2n) + O(1/n², 1/N)",
          iszero_sym(coeffj(simp(drift_f1(1 / ε, 1 / ε^2)), ε, 1) + 1 // 2))
    check("drift_f3 at n = N is (1 − 1/(2N))(1 − 2/(2N)) (iid offspring stage only)",
          iszero_sym(drift_f3(N, N) - (1 - 1 / (2N)) * (1 - 2 / (2N))))

    # ===================================================================================
    say("\n(c) selection: first-order cumulant changes and Leonov–Shiryaev")
    for trial in 1:4
        K = 7
        w = [Rational{BigInt}(rand(rng, 1:6)) for _ in 1:K]
        xk = [randq() for _ in 1:K]; xl = [randq() for _ in 1:K]; G = [randq() for _ in 1:K]
        u = G .- wmean(w, G)
        h = [randq() for _ in 1:K]
        # exact derivative at ε = 0 of the cumulants under weights w(1 + εh)
        @variables e
        we = w .* (1 .+ e .* h)
        cen(v) = v .- sum(we .* v) / sum(we)
        Ew(f) = sum(we .* f) / sum(we)
        T3 = Ew(cen(xk) .* cen(xl) .^ 2); C2 = Ew(cen(xk) .* cen(xl))
        dT = Symbolics.substitute(Symbolics.expand_derivatives(Differential(e)(T3)), Dict(e => 0))
        dC = Symbolics.substitute(Symbolics.expand_derivatives(Differential(e)(C2)), Dict(e => 0))
        val(x) = Rational{BigInt}(Symbolics.value(simp(x)))
        check("trial $trial: d/dε κ(k,l,l)|₀ = κ(k,l,l,h) exactly", val(dT) == kap(w, xk, xl, xl, h))
        check("trial $trial: d/dε κ(k,l)|₀ = κ(k,l,h) exactly", val(dC) == kap(w, xk, xl, h))
        u2 = u .^ 2
        check("trial $trial: κ(k,l,l,u²) = κ₅(k,l,l,u,u) + 2[κ(k,u)κ(l,l,u) + 2κ(l,u)κ(k,l,u)]",
              kap(w, xk, xl, xl, u2) == kap(w, xk, xl, xl, u, u) +
                                        2 * (kap(w, xk, u) * kap(w, xl, xl, u) + 2 * kap(w, xl, u) * kap(w, xk, xl, u)))
        check("trial $trial: κ(k,l,u²) = κ₄(k,l,u,u) + 2κ(k,u)κ(l,u)",
              kap(w, xk, xl, u2) == kap(w, xk, xl, u, u) + 2 * kap(w, xk, u) * kap(w, xl, u))
        # the combined first-order change for h = a u − u²/(2S_c)
        a = randq(); Sc = abs(randq()) + 1
        hh = a .* u .- u2 ./ (2Sc)
        lhs = kap(w, xk, xl, xl, hh)
        rhs = a * kap(w, xk, xl, xl, u) - (kap(w, xk, u) * kap(w, xl, xl, u) + 2 * kap(w, xl, u) * kap(w, xk, xl, u)) / Sc -
              kap(w, xk, xl, xl, u, u) / (2Sc)
        check("trial $trial: ΔT = aκ(k,l,l,u) − [κ(k,u)κ(l,l,u) + 2κ(l,u)κ(k,l,u)]/S_c − κ₅/(2S_c)", lhs == rhs)
        lhs2 = kap(w, xk, xl, hh)
        rhs2 = a * kap(w, xk, xl, u) - kap(w, xk, u) * kap(w, xl, u) / Sc - kap(w, xk, xl, u, u) / (2Sc)
        check("trial $trial: Δy = aκ(k,l,u) − κ(k,u)κ(l,u)/S_c − κ₄(k,l,u,u)/(2S_c)", lhs2 == rhs2)
    end
    # chromosome-level width: ∫ exp(−(g₁+g₂−θ)²/(2S′)) N(g₂; m, V₂) dg₂ ∝ exp(−(g₁+m−θ)²/(2(S′+V₂)))
    @variables g1 m2 V2 Sp θs
    # completing the square: the exponent's second derivative in g₁ after integrating g₂
    # (precision 1/S′ − (1/S′)²/(1/S′ + 1/V₂)) equals 1/(S′ + V₂)
    prec = 1 / Sp - (1 / Sp)^2 / (1 / Sp + 1 / V2)
    check("E[π | G₁] for Gaussian G₂: curvature 1/S_c, S_c = ω² + V_E + V₂", iszero_sym(prec - 1 / (Sp + V2)))

    # ===================================================================================
    say("\n(d) recombination")
    @variables ρ μ12 μ11 μ22 vk vl
    # homologs A, B independent copies; central moments; split w.p. ρ takes l from B
    Tn = (1 - ρ) * μ12 + ρ * 0 * vl           # E[x_k^A (x_l^B)²] = E[x_k]E[x_l²] = 0
    check("T′ = (1 − ρ) T", iszero_sym(Tn - (1 - ρ) * μ12))
    check("y′ = (1 − ρ) y (E[x_k^A x_l^B] = 0)", iszero_sym(((1 - ρ) * μ11 + ρ * 0) - (1 - ρ) * μ11))
    yn = (1 - ρ) * μ11
    m22n = (1 - ρ) * μ22 + ρ * vk * vl
    κ22 = μ22 - vk * vl - 2μ11^2
    κ22n = m22n - vk * vl - 2yn^2
    check("κ₂₂′ = (1 − ρ)κ₂₂ + 2ρ(1 − ρ)y² (cumulant not simply decaying)", iszero_sym(κ22n - ((1 - ρ) * κ22 + 2ρ * (1 - ρ) * μ11^2)))

    # ===================================================================================
    say("\n(e) mutation on a random chromosome")
    for trial in 1:2
        Mp = 6
        Bv = [randq() for _ in 1:Mp]; xv = [randq() for _ in 1:Mp]
        w1 = ones(Rational{BigInt}, Mp)
        T0 = kap(w1, Bv, xv, xv); y0 = kap(w1, Bv, xv); VB = kap(w1, Bv, Bv)
        @variables ee α2
        ET = 0; Ey = 0; Ey2 = 0
        for i in 1:Mp
            x2 = [j == i ? xv[j] + ee : xv[j] + 0 * ee for j in 1:Mp]
            cx = x2 .- sum(x2) / Mp; cB = Bv .- sum(Bv) / Mp
            Tm = sum(cB .* cx .^ 2) / Mp; ym = sum(cB .* cx) / Mp
            Eε(f) = (f0 = Symbolics.substitute(f, Dict(ee => 0)); f0 + coeffj(f, ee, 2) * α2)   # E ε = 0, E ε² = α²
            ET += Eε(Tm) / Mp; Ey += Eε(ym) / Mp; Ey2 += Eε(Symbolics.expand((ym - y0)^2)) / Mp
        end
        check("trial $trial: E[ΔT] = 0", iszero_sym(ET - T0))
        check("trial $trial: E[Δy] = 0", iszero_sym(Ey - y0))
        check("trial $trial: E[(Δy)²] = α² Var(B)/M²", iszero_sym(Ey2 - α2 * VB / Mp^2))
    end

    # ===================================================================================
    say("\n(f) one biallelic locus")
    @variables a p D β EB
    # E[(B − B̄)(a1_A − ap)²] = a² E[(B − B̄)(1_A(1 − 2p) + p²)] with E[(B − B̄)1_A] = D
    κBll = a^2 * ((1 - 2p) * D + p^2 * 0)
    check("κ(B, x_l, x_l) = a²(1 − 2p) Cov(B, 1_A)", iszero_sym(κBll - a^2 * (1 - 2p) * D))
    dv = a^2 * ((p + β * D) * (1 - p - β * D) - p * (1 - p))
    check("Δ(a²pq) for Δp = βD is βκ(B,l,l) − β²Cov(B,x_l)²", iszero_sym(dv - (β * κBll - β^2 * (a * D)^2)))

    # ===================================================================================
    say("\n(g) closed forms of the recursion")
    @variables λ Kc ν λy δy δT VG v Vc ω
    T = 0
    for t in 1:6
        T = λ * (T + β * Kc)
        check("T_$t = βKλ(1 − λ^$t)/(1 − λ) from T₀ = 0", iszero_sym(T - β * Kc * λ * (1 - λ^t) / (1 - λ)))
    end
    Q = ν / (1 - λy^2)
    check("stationary E[y²] = ν/(1 − λ_y²)", iszero_sym(λy^2 * Q + ν - Q))
    # leading order: 1 − λ_y² ≈ 2(ρ + δ_y), λ_T/(1 − λ_T) ≈ 1/(ρ + δ_T); ν = v V_s/(2n), K = −2Q
    loss_mine = 2β * (β * (-2 * (v * Vc / (2n)) / (2(ρ + δy)))) / (ρ + δT)
    check("steady-state loss per locus = −β² v V_s/(n(ρ + δ_y)(ρ + δ_T))",
          iszero_sym(loss_mine + β^2 * v * Vc / (n * (ρ + δy) * (ρ + δT))))
    barton = -v * β^2 * VG * ω / (2n * (ρ + δy)^2)       # Barton (1995): excess drift loss of 2v at N/N_e − 1 = b²V_G ω/(2(ρ+δ)²)
    mine_eq = Symbolics.substitute(loss_mine, Dict(δT => δy, Vc => ω * VG / 2))
    check("with δ_y = δ_T = δ and V_chr = V_G/2 it is Barton's excess drift loss", iszero_sym(mine_eq - barton))
    # exact factors vs their leading forms
    lam_y = (1 - ρ) * (1 - δy); lam_T = (1 - ρ) * (1 - δT)
    check("1 − λ_y² = 2(ρ + δ_y) + O(ρ², ρδ, δ²)",
          iszero_sym(coeffj(Symbolics.substitute(1 - lam_y^2, Dict(ρ => ε * ρ, δy => ε * δy)), ε, 1) - 2(ρ + δy)))
    check("1 − λ_T = ρ + δ_T + O(ρδ)",
          iszero_sym(coeffj(Symbolics.substitute(1 - lam_T, Dict(ρ => ε * ρ, δT => ε * δT)), ε, 1) - (ρ + δT)))

    # ===================================================================================
    say("\n(h) spot checks of validate/linked_decay_theory.jl")
    for (L, u) in ((7, 3.0), (12, 5.0), (30, 3.0))
        ρp, _ = QTLSim.pair_splits(L, u, :qtsim)
        S = qtsim_survival(L, u)
        mine = [S[i+1] - S[j+1] for i in 0:L-2 for j in i+1:L-1]
        check("qtsim_survival = pair_splits, L = $L, u = $u", maximum(abs.(mine .- ρp)) < 1e-14)
        G = ld_groups(L, u; nb = 4)
        check("ld_groups: partners per locus sum to L − 1 (L = $L)", all(a -> G.nb[a] == 0 || isapprox(sum(G.n[a, :]), L - 1), 1:4))
        tot = sum(G.nb[a] * sum(G.n[a, :] .* G.ρ[a, :]) for a in 1:4) / 2
        check("ld_groups: Σ ρ over pairs = pair_splits total (L = $L)", isapprox(tot, sum(ρp); rtol = 1e-12))
    end
    # rank-one form of A P Aᵀ
    ns = 5; c = rand(rng, ns); vb = 0.3; Sc = 7.0; P = randspd(rng, ns)
    A = I - (c * ones(ns)' + vb * I) / Sc
    AP = P .- (c * sum(P; dims = 1) .+ vb .* P) ./ Sc
    APA = AP .- (sum(AP; dims = 2) * c' .+ vb .* AP) ./ Sc
    check("rank-one form of A P Aᵀ (y_step!)", maximum(abs.(APA .- A * P * A')) < 1e-12)
    # T_step! :nosel with constant K and one class = closed form λ(1 − λᵗ)/(1 − λ)
    G1 = LDGroups(10, fill(9.0, 1, 1), fill(0.05, 1, 1), [10], [1])
    d = LDDrivers(β = 0.2, Sc = 50.0, dS = 0.01, n = 500.0, N = 1000.0, Vm = 0.0, V = [1.0], V2 = [0.1], k3 = [0.0],
                  g = [1.0], c = [1.0], Cbb = fill(1.0, 1, 1))
    Tv = [0.0]; Kv = [-0.3]
    lam = drift_f3(500.0, 1000.0) * 0.95
    for t in 1:20
        T_step!(Tv, G1, d, 1, Kv; mode = :nosel)
    end
    check("T_step!(:nosel) = βKλ(1 − λᵗ)/(1 − λ)", isapprox(Tv[1], 0.2 * -0.3 * lam * (1 - lam^20) / (1 - lam); rtol = 1e-12))
    # y_step! without selection, one class: stationary P = ν/(1 − λ²) with the two-stage noise
    st = LDState(G1)
    for t in 1:20000
        y_step!(st, G1, d, 1; select = false)
    end
    fS = 1 - 1 / 1000 + 1 / 2000; fO = 1 - 1 / 2000
    νtot = 0.9 * ((1 / 1000 - 1 / 2000) * fO^2 * 0.95^2 + 1 / 2000)        # Σ_B = (9/10)·g: the class holds 9 of the bin's 10 loci
    check("y_step!(select = false) stationary P = ν/(1 − λ²)", isapprox(st.P[1][1, 1], νtot / (1 - (fS * fO * 0.95)^2); rtol = 1e-8))

    finish(LOG, t0)
end

runmain(main, @__FILE__)
