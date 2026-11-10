# The rare-alleles closure for several traits (`run_moments(p; loci = :rare_alleles)` with
# `NoTraits = n > 1`), its equilibrium `rare_alleles_G`, the multivariate house of cards
# `house_of_cards_G` and the linkage-equilibrium closed form `house_of_gauss_G`; and the
# closure for genomes that do not recombine (haploid clones). The one-trait closure is
# src/rare_alleles.jl; the comparison with the kernel is validate/rare_alleles_traits.jl.
#
# The per-locus view with n traits. At one allele position, allelic values are vectors
# x ∈ ℝⁿ with within-population covariance c (n × n). Under weak Gaussian selection of
# covariance s = Ω + P on the genotype, an allele's marginal fitness is
# 1 − ½ yᵀ s⁻¹ y with y = x − x̄, so selection changes c by
#
#   Δc_ij = −½ Σ_kl s⁻¹_kl (E[y_i y_j y_k y_l] − c_ij c_kl)
#         = −(c s⁻¹ c)_ij − ½ Σ_kl s⁻¹_kl κ_ijkl,
#
# the multivariate weak-selection moment recursion (Turelli and Barton 1990, one trait;
# the n-trait form follows the same way from the covariance of y yᵀ with fitness), with κ
# the fourth cumulant tensor of y. The Gaussian closure κ = 0 is Lande's per-locus
# recursion c s⁻¹ c, the multilocus model of src/multilocus.jl (`hap_select`'s D = H s⁻¹ H).
#
# Rare alleles. A mutant of effect vector a, a ~ N(0, C) with C = M/(Ploidy·L·μ) the
# per-mutation effect covariance (`mutation_effect_covariance`; the kernel's `mutate!`), is
# removed at rate ½ aᵀ s⁻¹ a, so at the house-of-cards balance it has frequency
# q(a) = μ g(a)·2/(aᵀ s⁻¹ a) and the locus holds
#
#   c_HC = Σ q a aᵀ = 2μ X,      X = E_g[a aᵀ / (aᵀ s⁻¹ a)],
#
# while selection's fourth-moment term removes ½ Σ q a aᵀ (aᵀ s⁻¹ a) = μ C, the input.
# X is explicit: with C^½ s⁻¹ C^½ = U Λ Uᵀ (Λ = diag λ), X = C^½ U diag(e) Uᵀ C^½ with
#
#   e_i = E[w_i² / Σ_j λ_j w_j²] (w ~ N(0, I)) = ∫₀^∞ (1 + 2tλ_i)⁻¹ Π_j (1 + 2tλ_j)^(−½) dt,
#
# so Σ λ_i e_i = 1 (tr(s⁻¹X) = 1), e = 1/λ for one trait (X = s), and for two traits
# e_i = 1/(√λ_i (√λ_1 + √λ_2)). This is the multivariate house of cards: for spherical
# mutation and selection (C ∝ s) every λ_i = λ, X = s/n, and each trait holds 1/n of
# the one-trait house-of-cards variance 4Lμ V_s, the dilution by pleiotropy that Turelli
# (1985) and Wagner (1989) describe. (We could not read either paper's equations; the
# closed forms here are derived, and `house_of_cards_G` is checked against Monte Carlo
# averages over a in test/runtests.jl.)
#
# The closure. A linear closure for the fourth-moment term must give μ C at c = 2μX for
# every μ and must fix the shape of c as well as its size. (Two tempting forms fail:
# ½ tr(s⁻¹c) C removes the right amount but leaves the shape of c undetermined at the
# house-of-cards balance, and ¼(C s⁻¹ c + c s⁻¹ C) has c_HC = 2μ s, n times too large for
# spherical pleiotropy.) Taking each principal direction w_i of the pair (C, s⁻¹) to lose
# its variance at its house-of-cards rate ρ_i = 1/(2e_i), and a covariance between
# directions at the mean of their rates, gives
#
#   ½ Σ_kl s⁻¹_kl κ_ijkl  ≈  ½ (𝓡 c + c 𝓡ᵀ),   𝓡 = ½ C X⁻¹ = C^½ U diag(ρ) Uᵀ C^(−½),
#
# exact at the house-of-cards balance (½(𝓡 X + X 𝓡ᵀ)·2μ = μ C), and for one trait
# 𝓡 = α²/(2s), the one-trait closure κ₄ = α² c. Selection then changes the genic
# covariance by
#
#   Δc = −c s⁻¹ c − ½ (𝓡 c + c 𝓡ᵀ).
#
# As for one trait, LD changes a mutant's average effect on the genotype from a to β a,
# β = H H_genic⁻¹ (the regression of the genotype on a locus's alleles), so 𝓡 is taken
# with s̃⁻¹ = βᵀ s⁻¹ β in place of s⁻¹ (one trait: β²/s), and ρ_i is capped at 1.
#
# Linkage equilibrium, s = A = Ω + V_E I fixed, drift among ng genes: the per-gene fixed
# point is the stabilising solution of the continuous algebraic Riccati equation
#
#   c A⁻¹ c + J c + c Jᵀ = μ C,     J = ½ (𝓡 + I/ng),
#
# `house_of_gauss_G` (solved by the stable invariant subspace of its Hamiltonian; G = ng_u·c
# over the allele positions). For one trait it is Hermisson and Wagner's (2004) equation
# A4, `house_of_gauss_vg`; with μ → 0 it is the stochastic house of cards (a Lyapunov
# equation), with ng = ∞ and μ → 0 the house of cards above (while α² ≪ A, so that no rate
# ρ_i reaches its cap of 1); for many small mutations without drift it is the Gaussian
# G A⁻¹ G = 2L M (s = A), the weak-mutation form of `gaussian_G_loci`'s G (A + G)⁻¹ G = 2L M
# (validate/symbolic/check_traits.jl). No closed form in elementary
# functions exists unless C and A share eigenvectors.
#
# Genomes that do not recombine (haploid clones, as the kernel's haploid `reproduce!`; or a
# diploid with RecombPoissonMean = 0). A clone is one unit of inheritance: a mutant clone
# differs from its parent by the mutation it carries whatever locus it hit, and it is
# selected on its whole genotype. The genic/LD split of the multilocus model is then
# bookkeeping (without recombination the genic class grows and the LD class falls without
# bound while their sum converges, so β = H/H_genic → 0 would switch the closure off).
# The closure acts on the whole within-population variance of a haplotype, class by class,
# with s⁻¹ itself (β = 1): per unit, ΔV = U α² − V/ng − α² V/(2s) − V²/s with U = L μ the
# mutation rate per haplotype; at equilibrium V²/s + V (α²/(2s) + 1/ng) = L μ α² = V_m
# for haploids (`house_of_gauss_vg` with `gaussian_units`).

"""
    QTLSim.mutation_effect_covariance(p) -> C

The covariance of one mutation's effect vector, `C = M/(Ploidy·L·μ)` with `M` =
[`mutational_covariance`](@ref)`(p)` (an `n × n` matrix, `1 × 1` for one trait, whose entry is
[`QTLSim.mutation_effect_variance`](@ref)`(p)`).
"""
function mutation_effect_covariance(p::Params)
    p.MutationRatePerGamete > 0 || throw(ArgumentError("no mutations (MutationRatePerGamete = 0): C is undefined"))
    mutational_covariance(p) ./ (p.Ploidy * p.NoQTL * p.MutationRatePerGamete)
end

# The house-of-cards weights e_i = E[w_i²/Σ λ_j w_j²], w ~ N(0, I), for λ ≥ 0, as rates
# ρ_i = 1/(2e_i) (0 for a direction without selection): closed forms for n ≤ 2, otherwise
# the integral ½∫₀^∞ τ^(n/2−1) (τ+λ_i)⁻¹ Π_j (τ+λ_j)^(−½) dτ by the trapezoid rule in log τ,
# which converges exponentially for this analytic, exponentially decaying integrand.
function rare_hc_rates(λ::AbstractVector{<:Real})
    n = length(λ)
    n == 1 && return [λ[1] / 2]
    r = sqrt.(max.(λ, 0.0))
    n == 2 && return [r[1] * (r[1] + r[2]) / 2, r[2] * (r[1] + r[2]) / 2]
    pos = filter(>(0), λ)
    isempty(pos) && return zeros(n)
    lo, hi = log(minimum(pos)) - 80 / n, log(maximum(pos)) + 40
    h = 0.02
    e = zeros(n)
    for y in lo:h:hi
        τ = exp(y)
        base = 0.5 * τ^(n / 2) * prod(1 / sqrt(τ + l) for l in λ)
        for i in 1:n
            e[i] += h * base / (τ + λ[i])
        end
    end
    [λ[i] > 0 ? 1 / (2e[i]) : 0.0 for i in 1:n]
end

# C^½ U and its pseudo-inverse Uᵀ C^(−½) on the range of C, with the eigenvalues λ of
# C^½ sinv C^½ = U Λ Uᵀ
function rare_basis(C::AbstractMatrix, sinv::AbstractMatrix)
    F = eigen(Symmetric(Matrix{Float64}(C)))
    keep = F.values .> 1e-12 * max(maximum(F.values), 0.0)
    V, g = F.vectors[:, keep], sqrt.(F.values[keep])
    K = Symmetric((g .* (V' * sinv * V)) .* g')
    E = eigen(K)
    T = V * (g .* E.vectors)
    Tp = (E.vectors' ./ g') * V'
    T, Tp, max.(E.values, 0.0)
end

"""
    QTLSim.rare_alleles_rate(C, sinv) -> 𝓡

The rate operator of the rare-alleles closure for several traits (src/rare_alleles_traits.jl):
selection's fourth-cumulant term removes `½(𝓡 c + c 𝓡ᵀ)` of a locus's covariance `c`, with
`𝓡 = ½ C X⁻¹ = C^½ U diag(ρ) Uᵀ C^(−½)`, `C` the per-mutation effect covariance, `sinv` the
inverse selection covariance, `C^½ sinv C^½ = U Λ Uᵀ` and `ρ_i = 1/(2e_i)` the house-of-cards
rates (capped at 1). For one trait `α²·sinv/2`.
"""
function rare_alleles_rate(C::AbstractMatrix, sinv::AbstractMatrix)
    T, Tp, λ = rare_basis(C, sinv)
    ρ = min.(rare_hc_rates(λ), 1.0)
    (T .* ρ') * Tp
end

"""
    QTLSim.rare_hc_matrix(C, sinv) -> X

`X = E[a aᵀ/(aᵀ sinv a)]` for `a ~ N(0, C)`: a locus at the house-of-cards balance holds
`2μX` (src/rare_alleles_traits.jl). For one trait `s`; `tr(sinv X) = 1`.
"""
function rare_hc_matrix(C::AbstractMatrix, sinv::AbstractMatrix)
    T, _, λ = rare_basis(C, sinv)
    ρ = rare_hc_rates(λ)
    all(>(0), ρ) || throw(ArgumentError("a direction of mutation is not selected: the house of cards has no equilibrium"))
    _sym((T .* (1 ./ (2 .* ρ))') * T')
end

# selection's rare-alleles loss from a locus class of within-population covariance c:
# one trait (α2 a number, as the scalar code always did) or several (α2 = C a matrix),
# with the inverse selection covariance sinv (already including β when there is LD)
rare_loss(α2::Real, sinv::Real, c::Real) = min(α2 * sinv / 2, 1.0) * c
rare_loss(C::AbstractMatrix, sinv::AbstractMatrix, c::AbstractMatrix) = _sym(rare_alleles_rate(C, sinv) * c)

# the inverse selection covariance felt by a rare mutant through its average effect β a,
# β = H H_genic⁻¹: βᵀ s⁻¹ β (one trait β² s⁻¹; β = 1 when the genic part is not positive)
function rare_sinv(H::Real, Hg::Real, sinv::Real)
    β2 = Hg > 0 ? (H / Hg)^2 : 1.0
    β2 * sinv
end
function rare_sinv(H::AbstractMatrix, Hg::AbstractMatrix, sinv::AbstractMatrix)
    isposdef(Symmetric(Matrix{Float64}(Hg))) || return sinv
    β = H / Hg
    _sym(β' * sinv * β)
end

"""
    QTLSim.hap_clonal(ploidy, ρ) -> Bool

Whether the genome is inherited whole (haploid clones, or no pair of loci is ever split):
the rare-alleles closure then acts on a haplotype's whole variance (src/rare_alleles_traits.jl).
"""
hap_clonal(ploidy::Integer, ρ::AbstractVector) = ploidy == 1 || all(iszero, ρ)

"""
    house_of_cards_G(p)

The multivariate house-of-cards equilibrium G matrix of `L` loci (no drift, rare
mutations, `s = A` = [`selection_variance_matrix`](@ref)`(p)`): `n_u · 2μ_u X` with
`X = E[a aᵀ/(aᵀ A⁻¹ a)]`, `a ~ N(0, C)`, `C` = [`QTLSim.mutation_effect_covariance`](@ref)`(p)`,
over `n_u` units of inheritance with `μ_u` mutations each (`2L` genes with `μ` for a
recombining diploid, one clone with `Lμ` for a haploid). For one trait `2·n_u μ_u V_s`,
Turelli's `4Lμ V_s`; for spherical mutation and selection in `n` traits, `1/n` of that per
trait. See src/rare_alleles_traits.jl.
"""
function house_of_cards_G(p::Params)
    nu, μu = rare_units(p)
    nu * 2μu .* rare_hc_matrix(mutation_effect_covariance(p), inv(selection_variance_matrix(p)))
end

# the units of inheritance whose variance balances mutation separately, and the mutation
# rate of each: 2L genes at μ (recombining diploid), 2 haplotypes at Lμ (no recombination),
# 1 clone at Lμ (haploid)
function rare_units(p::Params)
    nu = gaussian_units(p)
    nu, nu == p.Ploidy * p.NoQTL ? p.MutationRatePerGamete : p.Ploidy * p.NoQTL * p.MutationRatePerGamete / nu
end

"""
    house_of_gauss_G(p; nsurv = p.PopSize, drift = true)

The multivariate "stochastic house of Gauss": the equilibrium G matrix of the rare-alleles
closure for loci in linkage equilibrium with `s = A` = [`selection_variance_matrix`](@ref)`(p)`,
`G = n_u c` with `c` the stabilising solution of the algebraic Riccati equation

    c A⁻¹ c + J c + c Jᵀ = μ_u C,    J = ½ (𝓡 + I/ng),

`C` = [`QTLSim.mutation_effect_covariance`](@ref)`(p)`, `𝓡` = [`QTLSim.rare_alleles_rate`](@ref)`(C, A⁻¹)`,
`ng = Ploidy·N_e` genes ([`effective_size`](@ref) with `nsurv` survivors; `∞` with
`drift = false`) and units as in [`house_of_cards_G`](@ref). For one trait it is
[`house_of_gauss_vg`](@ref) (Hermisson and Wagner 2004, eq. A4); for rare mutations without
drift [`house_of_cards_G`](@ref) (when `α² ≪ A`, as the rates `ρ_i` are capped at 1); for many
small mutations without drift the solution of `G A⁻¹ G = 2L·M`, the weak-mutation form of
[`gaussian_G_loci`](@ref) (`s = A` rather than `A + G`).
[`rare_alleles_G`](@ref) is the same closure with LD and `s = A + G`.
"""
function house_of_gauss_G(p::Params; nsurv::Real = p.PopSize, drift::Bool = true)
    n = p.NoTraits
    nu, μu = rare_units(p)
    C = mutation_effect_covariance(p)
    A = selection_variance_matrix(p)
    Ainv = inv(A)
    ng = drift ? p.Ploidy * effective_size(p.PopSize, nsurv) : Inf
    J = (rare_alleles_rate(C, Ainv) + (1 / ng) * I) ./ 2
    nu .* care_stabilising(J, Ainv, μu .* C)
end

# the stabilising solution X ⪰ 0 of X B X + J X + X Jᵀ = Q, from the stable invariant
# subspace of the Hamiltonian [−Jᵀ −B; −Q J] (Laub's Schur method)
function care_stabilising(J::AbstractMatrix, B::AbstractMatrix, Q::AbstractMatrix)
    n = size(J, 1)
    Ham = [-Matrix(J') -Matrix(B); -Matrix(Q) Matrix(J)]
    S = schur(Ham)
    ordschur!(S, real.(S.values) .< 0)
    Z = S.Z[:, 1:n]
    _sym(real.(Z[n+1:end, :] / Z[1:n, :]))
end

"""
    rare_alleles_G(p; drift = false)

The equilibrium G matrix of the trait moments lowering with `loci = :rare_alleles`: the
multilocus model of [`gaussian_G_loci`](@ref) whose genic covariance also loses
`½(𝓡 c + c 𝓡ᵀ)` to selection per generation (src/rare_alleles_traits.jl), with LD and
`s = Ω + P`; [`multilocus_equilibrium`](@ref)`(p; loci = :rare_alleles).VG`. For one trait
[`rare_alleles_vg`](@ref) (as a 1 × 1 matrix); for loci in linkage equilibrium
[`house_of_gauss_G`](@ref).
"""
function rare_alleles_G(p::Params; drift::Bool = false)
    G = multilocus_equilibrium(p; drift, loci = :rare_alleles).VG
    G isa Number ? fill(G, 1, 1) : G
end
