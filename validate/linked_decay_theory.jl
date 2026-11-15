# The derived recursion for the association between a locus's squared deviation and its
# linked genomic background (validate/linked_decay.jl; algebra checked in
# validate/symbolic/check_linked_decay.jl). Plain Julia, no dependencies beyond QTLSim, so the
# symbolic check can include it.
#
# Chromosome-level quantities (the adults, i.e. zygotes from random union):
#   x_k        allele value at locus k on a chromosome; v_k its variance (gene level)
#   B_s(l)     Σ_{k ∈ s, k ≠ l} x_k: the background of locus l in class s (partner bin × split-probability class)
#   y_{l,s}    Cov(B_s, x_l)                          (second order: allele–background LD)
#   T_{l,s}    κ(B_s, x_l, x_l) = E[ΔB_s d_l²]        (third order: the association this study is about)
#   C_l        Σ_s T_{l,s} = Cov(g_rest, d_l²) (cis)  (the genic accounting's other3 is 2β Σ_l C_l)
#   K_{l,s}    κ(B_s, B, x_l, x_l)                    (fourth order: the source)
#
# One generation (selection on adults, survival, gametes with QTSim's single crossover, offspring
# sampling, mutation), to first order in the selection gradient β and the chromosome-level
# curvature 1/S_c (S_c = ω² + V_E + Var(G_chr), the width of E[π | one chromosome]):
#
#   T_s ← f₃ (1 − ρ_s) [ T_s + β (K_s + K3_s) − (c_s (κ3_l + C_l) + 2 v_l T_s)/S_c ]
#   y_s ← f₁ (1 − ρ_s) [ y_s − (c_s (v_l + Y_l) + v_l y_s)/S_c ] + drift noise (Var ≈ v_l Var(B_s)/(2n))
#
# with c_s = Cov(B_s, G_chr), Y_l = Σ_s y_s, K3_s = κ(B_s, l, l, l), f₁ ≈ 1 − 1/(2n), f₃ ≈ 1 − 3/(2n)
# (n survivors), and, to leading order in 1/n, E[K_s] = −2 E[y_s Y_l] (variance-heterogeneity and
# three-locus parts of K are not seeded by sampling from a population in linkage equilibrium).
# Hence the decay rate of the association of class s:  ρ_s + 3/(2n) + (its share of) Cov(B,G)/S_c,
# and of the allele–background LD that feeds it:      ρ_s + 1/(2n) + Cov(B,G)/S_c.

using LinearAlgebra

const LD_NB = 16                                     # position bins (partner bins b and locus groups a)
const LD_REDGES = [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3] # split-probability classes r (upper edges; last class open)
const LD_NR = length(LD_REDGES) + 1

ld_bin(i0, L; nb = LD_NB) = min(nb, fld(i0 * nb, L) + 1)       # 0-based locus i0
ld_rclass(r) = (k = findfirst(e -> r < e, LD_REDGES); k === nothing ? LD_NR : k)

"""
    qtsim_survival(L, u)

S[c+1] = Σ_{K ≥ 1, K < L} P(K) P(lowest of K distinct points in 0:L−2 is ≥ c), c = 0..L−1: QTSim's
single crossover (`gamete!`): pair (i, j), 0-based i < j, is split iff the lowest point lies in
i..j−1, so ρ_ij = S[i+1] − S[j+1]. As `QTLSim.pair_splits` (checked).
"""
function qtsim_survival(L::Integer, u::Real)
    pois = QTLSim.Poisson(u)
    Ks = [K for K in 1:L-1 if QTLSim.pdf(pois, K) > 1e-17 || K <= u]
    S = zeros(L)
    for K in Ks
        w = QTLSim.pdf(pois, K)
        for c in 0:L-1
            s = 1.0
            for i in 0:K-1
                s *= max(L - 1 - c - i, 0) / (L - 1 - i)
            end
            S[c+1] += w * s
        end
    end
    S
end

"""
Classes of partners for loci of each group a: for s = (b, r) (partner bin, ρ class), `n[a, s]` the
mean number of partners k ≠ l per locus l ∈ a, `ρ[a, s]` their mean split probability, `nb[b]` the
loci in bin b, `na[a]` the loci in group a.
"""
struct LDGroups
    L::Int
    n::Matrix{Float64}      # (NB, NB·NR), s = (r − 1)·NB + b
    ρ::Matrix{Float64}
    nb::Vector{Int}
    bofs::Vector{Int}       # partner bin of class s
end

function ld_groups(L::Integer, u::Real; nb = LD_NB)
    S = qtsim_survival(L, u)
    ns = nb * LD_NR
    cnt = zeros(nb, ns); rsum = zeros(nb, ns)
    na = zeros(Int, nb)
    for i in 0:L-1
        a = ld_bin(i, L; nb)
        na[a] += 1
        for j in 0:L-1
            j == i && continue
            r = i < j ? S[i+1] - S[j+1] : S[j+1] - S[i+1]
            s = (ld_rclass(r) - 1) * nb + ld_bin(j, L; nb)
            cnt[a, s] += 1; rsum[a, s] += r
        end
    end
    ρ = [cnt[a, s] > 0 ? rsum[a, s] / cnt[a, s] : 0.0 for a in 1:nb, s in 1:ns]
    n = cnt ./ max.(na, 1)
    bofs = [(s - 1) % nb + 1 for s in 1:ns]
    LDGroups(L, n, ρ, na, bofs)
end

"""
Drivers of the recursion for one generation (measured in the kernel run, or theory):
β (regression of π/π̄ on g), Sc (chromosome-level width), VG/S (for the spectrum's δ comparator),
n (survivors), N (offspring), Vm (genic mutational input per generation, genic units),
per group a: V[a] = Σ_{l∈a} v_l, V2[a] = Σ v_l², k3[a] = Σ κ3_l; per bin b: g[b] = Σ_{k∈b} v_k,
c[b] = Cov(B_b, G_chr), Cbb = Cov(B_b, B_b′).
"""
Base.@kwdef struct LDDrivers
    β::Float64; Sc::Float64; dS::Float64; n::Float64; N::Float64; Vm::Float64
    V::Vector{Float64}; V2::Vector{Float64}; k3::Vector{Float64}
    g::Vector{Float64}; c::Vector{Float64}; Cbb::Matrix{Float64}
end

# per-locus-in-a class vectors: c_s and the background covariance Σ_B (s × s)
function class_c(G::LDGroups, d::LDDrivers, a)
    [d.c[G.bofs[s]] * G.n[a, s] / max(G.nb[G.bofs[s]], 1) for s in axes(G.n, 2)]
end
function class_cov(G::LDGroups, d::LDDrivers, a)
    ns = size(G.n, 2)
    Σ = zeros(ns, ns)
    for s in 1:ns, t in 1:ns
        b, e = G.bofs[s], G.bofs[t]
        fs = G.n[a, s] / max(G.nb[b], 1); ft = G.n[a, t] / max(G.nb[e], 1)
        Σ[s, t] = fs * ft * (d.Cbb[b, e] - (b == e ? d.g[b] : 0.0)) + (s == t ? fs * d.g[b] : 0.0)
    end
    Σ
end

drift_f1(n, N) = (1 - 1 / (2n) + 1 / (2N)) * (1 - 1 / (2N))
drift_f3(n, N) = (1 - 3 / (2n) + 3 / (2N)) * (1 - 1 / (2N)) * (1 - 2 / (2N))

"""
State of the LD part (level B): per group a the mean `m[a]` and raw second moment `P[a]` of the
class vector y_l summed over l ∈ a, and the association `T[a]` (class vector, summed over l ∈ a).
"""
mutable struct LDState
    m::Vector{Vector{Float64}}
    P::Vector{Matrix{Float64}}
    T::Vector{Vector{Float64}}
end
LDState(G::LDGroups) = (nb = size(G.n, 1); ns = size(G.n, 2);
    LDState([zeros(ns) for _ in 1:nb], [zeros(ns, ns) for _ in 1:nb], [zeros(ns) for _ in 1:nb]))

"""
One generation of the y (allele–background LD) moments for group a, in place. Selection
(Bulmer; first order in 1/S_c): y ← A y + b_l, A = I − (c 1ᵀ + v̄ I)/S_c, b_l = −c v_l/S_c; survival
sampling (factor and noise of a subset of n of N); gametes (1 − ρ_s); offspring sampling (2N draws);
mutation (noise M_a Σ_B/(2N)).
"""
function y_step!(st::LDState, G::LDGroups, d::LDDrivers, a; select = true)
    m, P = st.m[a], st.P[a]
    na = max(G.nb[a], 1)
    vbar = d.V[a] / na
    c = class_c(G, d, a); Σ = class_cov(G, d, a)
    ns = length(m)
    if select
        isc = 1 / d.Sc
        # A x = x − (c (1ᵀx) + v̄ x)/S_c
        Am = m .- (c .* sum(m) .+ vbar .* m) .* isc
        bm = -c .* (d.V[a] * isc)                    # Σ_l b_l
        AP = P .- (c * sum(P; dims = 1) .+ vbar .* P) .* isc      # A P
        APA = AP .- (sum(AP; dims = 2) * c' .+ vbar .* AP) .* isc  # A P Aᵀ
        Amb = Am * (-c' .* (vbar * isc))             # Σ_l A E[y_l] b_lᵀ ≈ v̄-weighted
        P .= APA .+ Amb .+ Amb' .+ (c * c') .* (d.V2[a] * isc^2)
        m .= Am .+ bm
    end
    fS = 1 - 1 / (2d.n) + 1 / (2d.N)
    fO = 1 - 1 / (2d.N)
    m .*= fS; P .*= fS^2
    P .+= Σ .* (d.V[a] * (1 / (2d.n) - 1 / (2d.N)))
    R = [1 - G.ρ[a, s] for s in 1:ns]
    m .*= fO .* R
    P .= (fO^2) .* (R .* P .* R')
    P .+= Σ .* (d.V[a] / (2d.N))
    P .+= Σ .* ((d.Vm / 2) * na / G.L / (2d.N))     # mutation: Σ_{l∈a} m_l = (V_m/2)(n_a/L)
    st
end

"""
One generation of the association T for group a. `K` (class vector, Σ_{l∈a}) is the source:
level B passes −2 P 1 (theory), level A the measured K + K3 split over ρ classes. `mode`:
:derived (the recursion above), :nosel (no selection term), :spectrum (the spectrum's δ:
T ← (1 − ρ)(1 − 1/(2n) − V_G/S)(T + βK)). `k3` adds the κ3 coupling −c_s κ3/S_c (level A).
"""
function T_step!(T::Vector{Float64}, G::LDGroups, d::LDDrivers, a, K; mode = :derived, k3 = 0.0)
    ns = length(T)
    R = [1 - G.ρ[a, s] for s in 1:ns]
    if mode === :derived
        na = max(G.nb[a], 1)
        vbar = d.V[a] / na
        c = class_c(G, d, a)
        C = sum(T)
        T .= T .+ d.β .* K .- (c .* (k3 + C) .+ 2vbar .* T) ./ d.Sc
        T .*= drift_f3(d.n, d.N) .* R
    elseif mode === :nosel
        T .= (T .+ d.β .* K) .* drift_f3(d.n, d.N) .* R
    elseif mode === :spectrum
        T .= (T .+ d.β .* K) .* (1 - 1 / (2d.n) - d.dS) .* R
    else
        error("mode $mode")
    end
    T
end

"""
Stationary y moments for fixed drivers (β irrelevant: the y recursion has no β term at this
order): iterate `y_step!` until the change is below `tol` relative.
"""
function y_stationary!(st::LDState, G::LDGroups, d::LDDrivers; maxit = 20000, tol = 1e-10)
    for a in axes(G.n, 1)
        for it in 1:maxit
            P0 = copy(st.P[a])
            y_step!(st, G, d, a)
            maximum(abs, st.P[a] .- P0) <= tol * max(maximum(abs, st.P[a]), 1e-300) && break
        end
    end
    st
end

"""
Split weights of a bin-level source over the ρ classes of that bin, for group a: proportional to
the uncoupled stationary E[y_s²] shape, n_s/(ρ_s + δ_y), δ_y = 1/(2n).
"""
function split_weights(G::LDGroups, a, n)
    nb = size(G.n, 1); ns = size(G.n, 2)
    w = [G.n[a, s] / (G.ρ[a, s] + 1 / (2n)) for s in 1:ns]
    tot = zeros(nb)
    for s in 1:ns
        tot[G.bofs[s]] += w[s]
    end
    [tot[G.bofs[s]] > 0 ? w[s] / tot[G.bofs[s]] : 0.0 for s in 1:ns]
end
