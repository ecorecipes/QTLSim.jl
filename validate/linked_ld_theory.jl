# The corrected recursion for the allele–background LD y_s = Cov(B_s, x_l) (validate/linked_ld.jl;
# algebra checked in validate/symbolic/check_linked_ld.jl). Extends validate/linked_decay_theory.jl
# (level B's y recursion, `y_step!`) and validate/linked_source_theory.jl (level B′, `yH_step!`).
#
# The stage accounting of linked_ld.jl (sections 1–3 of its log) names what the old recursion misses.
# One generation, chromosome level, first order in the selection coefficients and in 1/n:
#
# 0. (the main term) Heterogeneity among loci in Bulmer's step. Selection maps y_l ↦ A y_l + b_l with
#    the offset b_l = −c v_l/S proportional to the locus's variance v_l. The moments are carried summed
#    over the loci l of a group, m = Σ_l E[y_l], P = Σ_l E[y_l y_lᵀ], so the cross term
#    Σ_l A E[y_l] b_lᵀ needs Σ_l E[y_l] v_l. The old recursion used v̄ m (v̄ = Σv/n_a: averaging over loci
#    before the product). Under the recursion itself E[y_l] ∝ v_l exactly (every source of the mean is
#    ∝ v_l and the map does not depend on l), so Σ_l E[y_l] v_l = (Σv²/Σv) m (`vhet`). The cross term is
#    positive (mean and offset are both the negative Bulmer LD), and Σv²/Σv ≫ v̄ when the loci's
#    variances are unequal (n_a Σv²/(Σv)² = 1 + CV² of v_l): the old recursion under-weighted it.
# 1. Selection on individuals (weights h(G₁ + G₂)) changes the LD within each chromosome by
#    Δ = D y + b_l (D = −(c 1ᵀ + v̄ I)/S, b_l = −c v_l/S: Bulmer's step) AND creates the same
#    association between the two chromosomes of an individual, Cov(B_s^{(1)}, x_l^{(2)}) = Δ
#    (exact at first order: the mixed cumulants of independent chromosomes vanish, so only the
#    product term −Cov(B_s, G)Cov(x_l, G)/S survives, and Cov(B_s^{(1)}, G) = Cov(B_s, G_chr)).
# 2. A gamete takes x_l and B_s from the same chromosome with probability 1 − ρ_s, from the two
#    chromosomes with probability ρ_s. So the gamete pool's LD is
#        y_gam = (1 − ρ)∘u + ρ∘(Δ + Z̃),
#    u the survivors' within-chromosome LD, Z̃ the survivors' between-chromosome product beyond Δ:
#    under random pairing E[Z̃ | survivors] = −ρ∘u/(2n − 1) and
#        E[Z̃_s Z̃_t] = [v_l Cov(R_s, R_t) + (ρ∘u)_s (ρ∘u)_t]/(2n),   R_s = Σ_{k∈s} ρ_kl x_k
#    (n survivors). The old recursion had y_gam = (1 − ρ)∘(u_with Δ): it attenuated selection's change by
#    1 − ρ and had no recombinant noise. Together: y_gam = (1 − ρ)∘(y + noise) + Δ + ρ∘Z̃ — selection
#    acts on the gamete's LD in full, recombination only on the inherited part.
# 3. Mutation adds Var = (v_l Σ_{k∈s} m_k + δ_st m_l Var(B_s))/(2N) to E[y_s y_t]; the old recursion
#    kept only the second (Var(B_s) m_l, all s, t) term: the first is added on the diagonal.
# 4. The curvature 1/S of exact Gaussian conditioning (Bulmer's ω² + V_P), S = ω² + V_E + V_G
#    (individual level) instead of the first-order S_c = ω² + V_E + V_G,chr (`sind`).
#
# Dropped (stated in the log): the directional third-cumulant term b κ(B_s, x_l, G) of selection
# (≈ b(s_l y_s + Σ_{k∈s} s_k D_kl) for biallelic loci; needs per-locus skew-weighted moments), the
# fourth-cumulant term −κ(B_s, x_l, G, G)/(2S) (small before the move), O(1/n²), products of
# selection and drift terms, the cross-pair LD among background loci beyond Σ_B (class_cov).

"""
Class mean of ρ² per group a and class s (`ld_groups` layout), for the recombinant noise's diagonal
(the spread of ρ within a class).
"""
function class_rho2(L::Integer, u::Real; nb = LD_NB)
    S = qtsim_survival(L, u)
    ns = nb * LD_NR
    r2 = zeros(nb, ns); cnt = zeros(nb, ns)
    for i in 0:L-1, j in 0:L-1
        i == j && continue
        r = abs(S[i+1] - S[j+1])
        a = ld_bin(i, L; nb); s = (ld_rclass(r) - 1) * nb + ld_bin(j, L; nb)
        r2[a, s] += r^2; cnt[a, s] += 1
    end
    r2 ./ max.(cnt, 1)
end

"""
Covariance of the recombinant sums R_s = Σ_{k∈s} ρ_kl x_k per locus of group a (class form):
ρ_s ρ_t Σ_B[s, t] off the within-class diagonal, plus the within-class spread of ρ on the diagonal:
Σ_{k∈s} ρ_kl² v_k = ρ2_s f_s g_b (f_s = n[a,s]/nb[b]).
"""
function class_covR(G::LDGroups, d::LDDrivers, a, ρ2, Σ)
    ns = size(G.n, 2)
    ρ = G.ρ[a, :]
    C = (ρ * ρ') .* Σ
    for s in 1:ns
        b = G.bofs[s]
        fs = G.n[a, s] / max(G.nb[b], 1)
        C[s, s] += (ρ2[a, s] - ρ[s]^2) * fs * d.g[b]
    end
    C
end

# The selection maps are rank-one updates of a multiple of I: y ↦ α y + u (1ᵀy) (A: α = 1 − v̄/S,
# D = A − I: α = −v̄/S; u = −c/S in both), applied to moments in O(ns²).
lmul1(α, u, P) = α .* P .+ u * sum(P; dims = 1)
rmul1(P, α, u) = α .* P .+ sum(P; dims = 2) * u'
apply1(α, u, m) = α .* m .+ u .* sum(m)
# Σ_l (X y_l + x v_l)(W y_l + w v_l)ᵀ with Σ_l y_l v_l ≈ v̄ m (as y_step!), X = (αx, ux), W = (αw, uw)
crossmom(X, W, x, w, m, P, vbar, V2) = rmul1(lmul1(X..., P), W...) .+ apply1(X..., m) * (vbar .* w') .+
                                        (vbar .* x) * apply1(W..., m)' .+ (x * w') .* V2

"""
    yfix_step!(st, G, d, a; between, pairing, mutfull, sind, ρ2, couple, hsel)

One generation of the y moments (mean `m`, raw second moment `P`, summed over l ∈ a) for group a,
with the corrections above switched by keywords; all `false` is `y_step!` exactly (level B), and with
`couple = true` (H in `st.H`, an `LSState`) and `hsel` it is `yH_step!` exactly (level B′).
`st` is an `LSState` (its `.y` an `LDState`).
"""
function yfix_step!(st::LSState, G::LDGroups, d::LDDrivers, a; between = false, pairing = false, mutfull = false,
                    sind = false, vhet = false, ρ2 = nothing, couple = false, hsel = :none,
                    c = class_c(G, d, a), Σ = class_cov(G, d, a))
    m, P, H = st.y.m[a], st.y.P[a], st.H[a]
    na = max(G.nb[a], 1)
    vbar = d.V[a] / na
    ns = length(m)
    Sx = sind ? d.Sc + sum(d.Cbb) : d.Sc
    isc = 1 / Sx
    A = (1 - vbar * isc, -c .* isc)
    Dm = (-vbar * isc, -c .* isc)
    β = -c .* isc
    # weight of the mean × Bulmer-offset cross term Σ_l E[y_l] v_l: v̄ m (old) or, as E[y_l] ∝ v_l
    # under the recursion itself, (Σv²/Σv) m (`vhet`)
    vw = vhet ? d.V2[a] / max(d.V[a], 1e-300) : vbar
    # selection: within-chromosome (Bulmer) and, for `between`, the between-chromosome association Δ
    Pw = crossmom(A, A, β, β, m, P, vw, d.V2[a]); mw = apply1(A..., m) .+ β .* d.V[a]
    if couple && hsel === :bulmer
        Kst = H .- 2 .* P
        rs = vec(sum(Kst; dims = 2))
        Kst .-= (c * rs' .+ rs * c') .* isc
        H .= Kst .+ 2 .* Pw
    end
    if between
        Pwd = crossmom(A, Dm, β, β, m, P, vw, d.V2[a]); Pdd = crossmom(Dm, Dm, β, β, m, P, vw, d.V2[a])
        md = apply1(Dm..., m) .+ β .* d.V[a]
    end
    # survival sampling, ε = 1/(2n) − 1/(2N)
    εs = 1 / (2d.n) - 1 / (2d.N)
    fS = 1 - εs
    if couple
        Pu = (1 - 3εs) .* Pw .+ εs .* (H .+ Σ .* d.V[a])
        H .= (1 - 5εs) .* H .+ (2εs) .* (Pw .+ Pw')
    else
        Pu = fS^2 .* Pw .+ Σ .* (d.V[a] * εs)
    end
    mu = fS .* mw
    # gametes (single crossover)
    ρ = G.ρ[a, :]
    R = 1 .- ρ
    Rp = pairing ? R .- ρ ./ (2d.n - 1) : R
    mg = Rp .* mu
    Pg = Rp .* Pu .* Rp'
    if between
        mg .+= ρ .* (fS .* md)
        X = (Rp .* (fS^2 .* Pwd)) .* ρ'
        Pg .+= X .+ X' .+ (ρ .* (fS^2 .* Pdd)) .* ρ'
    end
    if pairing
        CR = class_covR(G, d, a, ρ2, Σ)
        Pg .+= (d.V[a] .* CR .+ (ρ .* Pu) .* ρ') ./ (2d.n)
    end
    couple && (H .*= (1 .- st.span[a]))
    # offspring sampling, ε = 1/(2N)
    εo = 1 / (2d.N)
    if couple
        Pn = (1 - 3εo) .* Pg .+ εo .* (H .+ Σ .* d.V[a])
        H .= (1 - 5εo) .* H .+ (2εo) .* (Pg .+ Pg')
        mg .*= (1 - εo)
    else
        Pn = (1 - εo)^2 .* Pg .+ Σ .* (d.V[a] * εo)
        mg .*= (1 - εo)
    end
    # mutation (m per locus per chromosome = (V_m/2)/L)
    mloc = (d.Vm / 2) / G.L
    Pn .+= Σ .* (mloc * na / (2d.N))
    if mutfull
        for s in 1:ns
            Pn[s, s] += d.V[a] * G.n[a, s] * mloc / (2d.N)
        end
    end
    m .= mg; P .= Pn
    st
end

"""
Stationary y (and H) moments for fixed drivers, iterating `yfix_step!` per group on its active
classes only (n[a, s] > 0; the others stay zero exactly), then writing them back into `st`.
"""
function yfix_stationary!(st::LSState, G::LDGroups, d::LDDrivers; maxit = 40000, tol = 1e-9, ρ2 = nothing, kw...)
    for a in axes(G.n, 1)
        G.nb[a] == 0 && continue
        idx = findall(>(0), G.n[a, :])
        Gr = LDGroups(G.L, G.n[:, idx], G.ρ[:, idx], G.nb, G.bofs[idx])
        sr = LSState(Gr)
        sr.y.m[a] .= st.y.m[a][idx]; sr.y.P[a] .= st.y.P[a][idx, idx]; sr.H[a] .= st.H[a][idx, idx]
        ρ2r = ρ2 === nothing ? nothing : ρ2[:, idx]
        c = class_c(Gr, d, a); Σ = class_cov(Gr, d, a)
        for it in 1:maxit
            chk = it % 25 == 0
            chk && (P0 = copy(sr.y.P[a]); H0 = copy(sr.H[a]))
            yfix_step!(sr, Gr, d, a; c, Σ, ρ2 = ρ2r, kw...)
            chk || continue
            sc = max(maximum(abs, sr.y.P[a]), 1e-300)
            max(maximum(abs, sr.y.P[a] .- P0), maximum(abs, sr.H[a] .- H0)) <= tol * sc && break
        end
        fill!(st.y.m[a], 0.0); fill!(st.y.P[a], 0.0); fill!(st.H[a], 0.0)
        st.y.m[a][idx] .= sr.y.m[a]; st.y.P[a][idx, idx] .= sr.y.P[a]; st.H[a][idx, idx] .= sr.H[a]
    end
    st
end
