# The source of the association T_s = κ(B_s, x_l, x_l) (validate/linked_source.jl; algebra checked in
# validate/symbolic/check_linked_source.jl). Extends validate/linked_decay_theory.jl.
#
# Exact decomposition (centred moments of the adults' chromosomes, any population):
#   K_s = κ(B_s, B, x_l, x_l) = H_s − 2 y_s Y_l,   H_s = Cov(B_s B, d_l²) = μ(B_s, B, l, l) − Cov(B_s, B) v_l.
# linked_decay closed K_s ≈ −2 y_s Y_l, i.e. it dropped H_s, the variance heterogeneity of locus l
# across background states (two-locus H2 = Σ_{k∈s} Cov(d_k², d_l²) and three-locus H3).
#
# The moments of one sampling stage (M chromosomes drawn iid from a pool; exact to O(1/M²)):
#   E[Ĥ_st]     = (1 − 5/M) H_st + (2/M)(y_s y_t + y_t y_s)          (H_st = Cov(B_s B_t, d_l²))
#   E[ŷ_s ŷ_t]  = (1 − 3/M) y_s y_t + (1/M)(H_st + Cov(B_s, B_t) v_l)
# i.e. drift seeds the variance heterogeneity from the squared allele–background LD (the
# allele-value form of Hill & Robertson's (1968) / Ohta & Kimura's (1969) E[D(1−2p)(1−2q)] moment)
# and H feeds back into the LD's sampling variance. A subset of m of M (survival) acts to leading
# order as iid sampling with 1/M replaced by 1/m − 1/M. Single-crossover recombination:
#   y_s ← (1 − ρ_s) y_s,  H_st ← (1 − ρ_{st}) H_st,  ρ_{st} the probability that the crossover falls
# among the three loci (k ∈ s, j ∈ t, l): max(ρ_s, ρ_t) when k and j lie on the same side of l,
# ρ_s + ρ_t on opposite sides. Mutation leaves E[H] unchanged (exactly).
# Neutral fixed point of one pair at split probability ρ (one stage, ε = 1/M, ρ, ε small):
#   H/(2y²) = 2ε/(5ε + ρ) ≤ 0.4, so for ρ ≫ 1/M the closure K = −2yY holds and for the tightly
#   linked pairs of QTSim's map (ρ ≲ 1/M near the chromosome's end) it overstates |K| by up to
#   a factor 1/(1 − 0.4) = 1.67.
#
# Selection: the first-order change of H_st under weights 1 + h is (exact)
#   ΔH_st = κ₅(B_s, B_t, l, l, h) + 2[κ(B_s, l, h) y_t + y_s κ(B_t, l, h)];
# what the closure keeps of it is set by the `hsel` option below (see linked_source.log).

"""
Per group a: the three-locus split probability ρ_{st} of classes s, t (side by position bin of
the partner relative to the group; the group's own bin counts as the same side).
"""
function span_matrix(G::LDGroups, a)
    ns = size(G.n, 2)
    S = zeros(ns, ns)
    for s in 1:ns, t in 1:ns
        bs, bt = G.bofs[s], G.bofs[t]
        opp = (bs < a && bt > a) || (bs > a && bt < a)
        ρs, ρt = G.ρ[a, s], G.ρ[a, t]
        S[s, t] = opp ? min(1.0, ρs + ρt) : max(ρs, ρt)
    end
    S
end

"""
State of level B′: per group a the y moments (`m`, `P`, as LDState), the variance heterogeneity
H (class × class, summed over l ∈ a) and the association T.
"""
mutable struct LSState
    y::LDState
    H::Vector{Matrix{Float64}}
    span::Vector{Matrix{Float64}}
end
LSState(G::LDGroups) = (nb = size(G.n, 1); ns = size(G.n, 2);
    LSState(LDState(G), [zeros(ns, ns) for _ in 1:nb], [span_matrix(G, a) for a in 1:nb]))

"""
One sampling stage of rate ε for group a: P ← (1 − 3ε)P + ε(H + Σ_B V_a), H ← (1 − 5ε)H + 2ε(P + Pᵀ),
m ← (1 − ε)m (all simultaneous; P the raw second moment E[y yᵀ] summed over l ∈ a).
`couple = false` drops H's feedback and seeding (linked_decay's y recursion, apart from P's factor).
"""
function sample_stage!(st::LSState, Σ, Va, ε, a; couple = true)
    m, P, H = st.y.m[a], st.y.P[a], st.H[a]
    Pold = copy(P)
    if couple
        P .= (1 - 3ε) .* P .+ ε .* (H .+ Σ .* Va)
        H .= (1 - 5ε) .* H .+ (2ε) .* (Pold .+ Pold')
    else
        P .= (1 - ε)^2 .* P .+ ε .* Σ .* Va
    end
    m .*= (1 - ε)
    st
end

"""
One generation of the y and H moments for group a (level B′): selection on y (Bulmer, as
`y_step!`), `hsel` on H, survival stage (ε = 1/(2n) − 1/(2N)), gametes (single crossover:
P ← RPR, H ← H ∘ (1 − ρ_{st})), offspring stage (ε = 1/(2N)), mutation (P only).
"""
function yH_step!(st::LSState, G::LDGroups, d::LDDrivers, a; select = true, couple = true, hsel = :none, Tvec = nothing,
                  c = class_c(G, d, a), Σ = class_cov(G, d, a))
    m, P = st.y.m[a], st.y.P[a]
    H = st.H[a]
    na = max(G.nb[a], 1)
    vbar = d.V[a] / na
    ns = length(m)
    if select
        isc = 1 / d.Sc
        Pold_sel = copy(P)
        Am = m .- (c .* sum(m) .+ vbar .* m) .* isc
        bm = -c .* (d.V[a] * isc)
        AP = P .- (c * sum(P; dims = 1) .+ vbar .* P) .* isc
        APA = AP .- (sum(AP; dims = 2) * c' .+ vbar .* AP) .* isc
        Amb = Am * (-c' .* (vbar * isc))
        P .= APA .+ Amb .+ Amb' .+ (c * c') .* (d.V2[a] * isc^2)
        m .= Am .+ bm
        if hsel === :bulmer
            # stabilising selection on the fourth cumulant K = H − 2P (leading Leonov–Shiryaev terms):
            # ΔK_st = −[c_s K_t· + c_t K_s·]/S_c  (K_s· = Σ_t K_st; K3, κ₆, products of third-order
            # cumulants and the directional βκ₅ dropped); H then follows from the selected P
            Kst = H .- 2 .* Pold_sel
            rs = vec(sum(Kst; dims = 2))
            Kst .-= (c * rs' .+ rs * c') .* isc
            H .= Kst .+ 2 .* P
        end
    end
    εs = 1 / (2d.n) - 1 / (2d.N)
    sample_stage!(st, Σ, d.V[a], εs, a; couple)
    R = [1 - G.ρ[a, s] for s in 1:ns]
    m .*= R
    P .= R .* P .* R'
    H .*= (1 .- st.span[a])
    sample_stage!(st, Σ, d.V[a], 1 / (2d.N), a; couple)
    P .+= Σ .* ((d.Vm / 2) * na / G.L / (2d.N))
    st
end

"""
Stationary y and H moments for fixed drivers (iterate `yH_step!`).
"""
function yH_stationary!(st::LSState, G::LDGroups, d::LDDrivers; maxit = 40000, tol = 1e-8, kw...)
    for a in axes(G.n, 1)
        G.nb[a] == 0 && continue
        c = class_c(G, d, a); Σ = class_cov(G, d, a)
        for it in 1:maxit
            chk = it % 50 == 0
            chk && (P0 = copy(st.y.P[a]); H0 = copy(st.H[a]))
            yH_step!(st, G, d, a; c, Σ, kw...)
            chk || continue
            sc = max(maximum(abs, st.y.P[a]), 1e-300)
            max(maximum(abs, st.y.P[a] .- P0), maximum(abs, st.H[a] .- H0)) <= tol * sc && break
        end
    end
    st
end

"""
Predicted source vector of group a: K_s = Σ_t (H_st − 2P_st).
"""
source_K(st::LSState, a) = vec(sum(st.H[a]; dims = 2)) .- 2 .* vec(sum(st.y.P[a]; dims = 2))

"""
Neutral fixed point of one pair (one class, both stages, recombination ρ, background variance
product Cv): returns (P, H). Solves the 2×2 linear map exactly.
"""
function pair_fixed_point(ρ, εs, εo, Cv)
    # one stage: [P; H] ← A(ε)[P; H] + ε[Cv; 0]
    A(ε) = [1-3ε ε; 4ε 1-5ε]
    Rm = [(1 - ρ)^2 0; 0 1-ρ]
    M = A(εo) * Rm * A(εs)
    b = A(εo) * Rm * [εs * Cv; 0] .+ [εo * Cv; 0]
    x = (I - M) \ b
    (P = x[1], H = x[2])
end
