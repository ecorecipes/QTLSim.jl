# The directional (first order in the selection gradient β) terms of the allele–background closure
# (validate/linked_dir.jl; algebra checked in validate/symbolic/check_linked_dir.jl). Extends
# validate/linked_decay_theory.jl (T recursion), validate/linked_source_theory.jl (H coupling) and
# validate/linked_ld_theory.jl (level B′*, `yfix_step!`). Plain Julia.
#
# Notation (chromosome level, adults; locus l, partner class s, B_s the partner sum, B = G_chr − x_l):
#   y_s = Cov(B_s, x_l), T_s = κ(B_s,l,l), J_s = κ(B_s,B,l), u_s = T_s + J_s = κ(B_s,l,G_chr),
#   H_s = Cov(B_sB, d_l²), K_s = κ(B_s,B,l,l) = H_s − 2y_sY_l, K3_s = κ(B_s,l,l,l).
# Exact first-order changes under the weights 1 + β(G − Ḡ) (the linear part of selection):
#   Δy_s = β u_s,  ΔT_s = β(K_s + K3_s),  ΔJ_s = β(K_s + κ(B_s,B,B,l)),  ΔK_s = β κ₅(B_s,B,l,l,G),
#   ΔH_s = ΔK_s + 2β(u_sY_l + y_sU_l),  Δκ3_l = β(κ4_l + κ(l,l,l,B)).
# The odd quantities (T, J, u, κ₅, κ3_l) vanish in expectation before the move (stabilising selection
# around the optimum is symmetric under G → −G); under the move their means are O(β).
#
# Biallelic projection (exact for a biallelic locus l, any background): d² = s d + v and
# d³ = (s² + v)d + sv (s = κ3/v), hence T_s = s_l y_s, H_s = s_l J_s, K3_s = w_l y_s (w = s² − 2v = κ4/v)
# and, with Δs_l = −2β(v_l + Y_l),
#   ΔH_s = β[w_l J_s − 2Y_l u_s + s_l κ(B_s,B,B,l)]     (u_s = T_s + J_s)                      (†)
#
# Level B‡ (lowest controlled order; second order in β for the even states, first order for the odd
# ones, which are driven by the zeroth-order (level B′*) even states — no feedback):
#  (i)  T recursion: the source gets K3_s = ⟨w⟩_l m_s (biallelic K3 = w_l y_s with the locus's w
#       decoupled from its LD, variance weights: ⟨w⟩_l = Σκ4/Σv of the group), the coupling
#       −c_s κ3/S_c gets the group's Σκ3 (a single-locus spectrum moment, a measured driver like Σv²,
#       as linked_decay's level A);
#  (ii) the own-skew part of the y term, Q^T_st = Σ_l E[T_s y_t] = Σ_l E[s_l y_s y_t] (biallelic), carried:
#       selection  Q^T += β[(⟨s²⟩_l + ⟨w⟩_l) P_st + h K_s m_t] − 3v̄ Q^T/S_c
#                  (sources E[T_s T_t] = ⟨s²⟩_l P_st, E[K3_s y_t] = ⟨w⟩_l P_st, E[K_s y_t] = h K_s m_t,
#                   h = Σv²/(Σv)², ⟨s²⟩_l = Σ(κ3²/v)/Σv; decay: T's 2v̄/S_c and y's v̄/S_c),
#       each sampling stage of rate ε  Q^T ← (1 − 5ε)Q^T − 2ε Q^Tᵀ + ε κ3_a Σ_B  (E[T̂_sŷ_t] to O(ε): the
#                  seed μ(B_s,B_t,l,l,l) − v_l μ(B_s,B_t,l) is, biallelic, s_l v_l Cov(B_s,B_t) + s_l² μ(B_s,B_t,l);
#                  the three-locus moment μ(B_s,B_t,l) is dropped),
#       recombination  Q^T ← (1 − ρ_s)(1 − ρ_t) Q^T (both pool statistics);
#       increments of the even states at selection: P += β(Q^T + Q^Tᵀ), H += 2β(Q^T + Q^Tᵀ) (K unchanged),
#       m += β T̄ (the mean of Δy = β u, own-skew part);
#  (iii) not closed at this order (dropped, quantified in validate/linked_dir.log): the partner part
#       J of u (dominated by three-locus third cumulants Σ_{k≠j} κ(k,j,l): the pairwise part Σ_k s_k D_kl
#       is −0.4 of J), hence E[T_sJ_t], E[J_sJ_t], Σ_t E[J_s y_t]; the K term β κ₅ (by (†) a ~5% residual
#       of terms −β Var(B)T_s − 2β Cov(B_s,B)C_l that the cubic heterogeneity Cov(B_sB², d_l²) cancels);
#       the quadratic-tilt term (β² − 1/S)κ₆/2 on K; drift seeding of Q^T by other odd fifth moments; the
#       stabilising coupling terms c_s(·) of Q^T.

"""
Single-locus spectrum drivers of one generation (per group): `w` = ⟨w⟩ = Σκ4/Σv, `s2` = ⟨s²⟩ =
Σ(κ3²/v)/Σv, `h` = Σv²/(Σv)².
"""
struct DirDrivers
    w::Vector{Float64}; s2::Vector{Float64}; h::Vector{Float64}
end
function dir_drivers(V, V2, sk4, ss2v)
    NB = length(V)
    DirDrivers([V[a] > 0 ? sk4[a] / V[a] : 0.0 for a in 1:NB], [V[a] > 0 ? ss2v[a] / V[a] : 0.0 for a in 1:NB],
               [V[a] > 0 ? V2[a] / V[a]^2 : 0.0 for a in 1:NB])
end

"""
One generation of Q^T (group a, active classes), driven by the zeroth-order even states `P`, `m` and
source `Ks` = Σ_t(H − 2P). In place.
"""
function qt_step!(QT::Matrix{Float64}, P, m, Ks, g::LDGroups, d::LDDrivers, dd::DirDrivers, a; Σ = class_cov(g, d, a))
    ns = size(QT, 1)
    na = max(g.nb[a], 1)
    vbar = d.V[a] / na
    QT .+= d.β .* ((dd.s2[a] + dd.w[a]) .* P .+ dd.h[a] .* (Ks * m')) .- (3vbar / d.Sc) .* QT
    R = [1 - g.ρ[a, s] for s in 1:ns]
    εs = 1 / (2d.n) - 1 / (2d.N); εo = 1 / (2d.N)
    QT .= (1 - 5εs) .* QT .- (2εs) .* QT' .+ (εs * d.k3[a]) .* Σ
    QT .= R .* QT .* R'
    QT .= (1 - 5εo) .* QT .- (2εo) .* QT' .+ (εo * d.k3[a]) .* Σ
    QT
end

"""
The increments of level B‡'s even states at the selection stage: (δP, δH, δm) from Q^T and the
zeroth-order T̄.
"""
function qt_increments(QT, T0, d::LDDrivers)
    δP = d.β .* (QT .+ QT')
    (; δP, δH = 2 .* δP, δm = d.β .* T0)
end

"""
The K3 source of the T recursion (biallelic closure): K3_s = ⟨w⟩_l m_s.
"""
k3_source(m, dd::DirDrivers, a) = dd.w[a] .* m
