# The rare-alleles closure with locally favoured alleles (`loci = :rare_local`): the
# rare-alleles closure of src/rare_alleles.jl and src/rare_alleles_traits.jl for a
# metapopulation whose demes lag behind different optima. The algebra is checked with
# Symbolics.jl in validate/symbolic/check_closures.jl, the theorems are in
# lean/CategoricalQuantGen/LocalAdaptation.lean, and the comparison with the kernel is
# validate/local_adaptation.jl.
#
# The problem. The rare-alleles closure removes a locus's genic variance at the
# house-of-cards rate: a rare mutant of effect a is selected against at s(a) = ½aᵀs⁻¹a,
# the curvature of the fitness surface at the deme's mean (src/rare_alleles.jl), so that
# at the balance it holds q(a) = μ g(a)/s(a), and the closure's extra loss ½(𝓡c + c𝓡ᵀ),
# 𝓡 = ½CX⁻¹ (one trait α²/(2s)), is exact there. In a deme whose mean lags behind its
# optimum by λ_d = θ_d − z̄_d, a mutant's marginal fitness is
#
#   w_d(a) ≈ 1 + aᵀℓ_d − ½ aᵀ s_d⁻¹ a,       ℓ_d = s_d⁻¹ λ_d,
#
# (expand exp(−½(y + a − λ)ᵀs⁻¹(y + a − λ)) averaged over the background y ~ N(0, P) to
# second order in a), so a mutant of the right sign is favoured where the deme lags.
# One population, or demes that lag alike, have a common lag; its directional part
# moves the mean and, like the current closure for one deme, we leave it out (it averages
# out over the signs of a, and the closure is unchanged there). Under local adaptation
# with migration the lags differ among demes (migration pulls each deme's mean toward
# the others'), and an allele's frequency is shared through migration.
#
# A rare mutant lineage's frequencies x_d in the k demes. Selection then migration
# (backward matrix P, P_dc the fraction of deme d's adults born in c; `migration_mixing`):
#
#   x′ = P (I + Σ(a)) x,      Σ = diag(σ_d),   σ_d(a) = aᵀℓ_d − ½ aᵀ s_d⁻¹ a,
#
# whose growth rate is the leading eigenvalue of P(I + Σ). For a migration class
# (irreducible P, stationary π, πP = π), perturbation in σ of the simple eigenvalue 1 of P
# (right eigenvector 𝟙, left π) gives
#
#   r(a) = π·σ + σᵀ Π (I − P)^# P σ + O(σ³),     Π = diag(π),
#
# with (I − P)^# = (I − P + 𝟙π)⁻¹ − 𝟙π the group inverse. As (I − P)^#𝟙 = 0 and
# π(I − P)^# = 0, the second-order term depends only on the deviations of the lags from
# their π-mean, ℓ_d − ℓ̄; the first-order term π·σ is the common lag (left out, as above)
# and the stabilising part −½aᵀ(Σπ_d s_d⁻¹)a. To second order in a,
#
#   r(a) = aᵀℓ̄ − ½ aᵀ (s̄⁻¹ − 2Λ) a,     Λ = Σ_dc W_dc (ℓ_d − ℓ̄)(ℓ_c − ℓ̄)ᵀ (symmetrised),
#   W = Π (I − P)^# P,
#
# so the variation of the lags among demes coupled by migration weakens selection against
# rare mutants by 2Λ (the covariance between a mutant's local frequency and its local
# advantage: it is commoner where it is favoured). For the island model (P = (1 − m′)I +
# m′𝟙πᵀ on the deviations, π uniform) W acts as (1 − m′)/m′ times the π-variance:
# Λ = (1 − m′)/m′ · Var_π(ℓ), the classical gain of a rare allele from spatially varying
# selection under migration (the second-order term of the protected-polymorphism
# analysis; Nagylaki 1992, ch. 6; Bürger 2014, "A survey of migration–selection models in
# population genetics", Discrete Contin. Dyn. Syst. B 19:883, section 3; we checked the
# expansion numerically against the eigenvalue, not against those texts' equations).
# When it exceeds the stabilising curvature, rare mutants of small effect invade: they
# are locally favoured alleles held by migration–selection balance, not rare alleles
# held by mutation, and their variance behaves as the Gaussian model's (Barton 1999,
# "Clines in polygenic traits", Genet. Res. 74:223: gene flow inflates the genetic variance
# as mutations of small effect would; we read the abstract only).
#
# The closure. A rare mutant is selected against at ½aᵀ(s_d⁻¹ − 2Λ)a in deme d, so the
# house-of-cards balance, and the closure exact there, are the current ones with
#
#   s_d⁻¹  →  [s_d⁻¹ − 2Λ]₊     (the positive part: directions in which mutants are
#                                 favoured are not held by rare alleles; 0 for one trait),
#
# in the rare-alleles loss only (with the LD regression β as before, βᵀ(s⁻¹ − 2Λ)β); the
# Gaussian part of selection, the means and everything else are unchanged. One trait:
# the loss is β²α² max(0, s⁻¹ − 2Λ)/2 of the genic variance. No constant is fitted: Λ
# comes from the deme means, optima and widths of the moments and from P, the migration
# matrix the moments already use. Exact reductions (validate/symbolic, Lean
# LocalAdaptation.lean, test/runtests.jl): Λ = 0, so the closure is the current one,
#
# - for one deme (k = 1: no deviations),
# - for demes that lag alike (a common optimum and identical demes: ℓ_d all equal), and
# - without migration (P = I: every deme is its own migration class, Λ = 0),
#
# and Λ → 0 as migration becomes complete (m′ → 1: the island W vanishes).
# Simplifications: second order in σ (accurate when σ ≪ m′, i.e. migration is strong
# against selection on one mutant; for weaker migration the expansion overstates the
# gain, but there the lags, which migration causes, are small); the lags are taken from
# the deterministic means after migration and the widths of the last selection; Λ is one
# matrix per migration class; the skew of the alleles is ignored, as before.

"""
    QTLSim.migration_classes(P) -> Vector{Vector{Int}}

The migration classes of the demes: the connected components of the graph with an edge
between `d` and `c` when `P[d, c] > 0` or `P[c, d] > 0`.
"""
function migration_classes(P::AbstractMatrix)
    k = size(P, 1)
    seen = falses(k)
    comps = Vector{Int}[]
    for d0 in 1:k
        seen[d0] && continue
        comp, stack = Int[], [d0]
        seen[d0] = true
        while !isempty(stack)
            d = pop!(stack)
            push!(comp, d)
            for c in 1:k
                if !seen[c] && c != d && (P[d, c] > 0 || P[c, d] > 0)
                    seen[c] = true; push!(stack, c)
                end
            end
        end
        push!(comps, sort(comp))
    end
    comps
end

"""
    QTLSim.local_weights(P) -> W

`W = Π (I − P)^# P` for an irreducible backward migration matrix `P` with stationary
distribution `π` (`Π = diag(π)`, `(I − P)^# = (I − P + 𝟙π)⁻¹ − 𝟙π` the group inverse): the
second-order coefficient of a rare lineage's growth rate, `r = π·σ + σᵀWσ + O(σ³)`, the
leading eigenvalue of `P (I + diag σ)`. See src/local_alleles.jl.
"""
function local_weights(P::AbstractMatrix)
    k = size(P, 1)
    A = Matrix{Float64}(I, k, k) - P
    π = [A'; ones(1, k)] \ [zeros(k); 1.0]
    one_π = ones(k) * π'
    Diagonal(π) * (inv(A + one_π) - one_π) * P
end

"""
    QTLSim.local_gain(P, ℓ) -> Vector

The gain `Λ_d` of a rare mutant from the variation of the lags among the demes coupled by
migration, for each deme: `Λ = Σ_dc W_dc (ℓ_d − ℓ̄)(ℓ_c − ℓ̄)ᵀ` (symmetrised) over the
deme's migration class ([`QTLSim.migration_classes`](@ref)), `W` =
[`QTLSim.local_weights`](@ref)`(P)` on the class, `ℓ_d = s_d⁻¹(θ_d − z̄_d)` (a number per deme
for one trait, a vector for several) and `ℓ̄` their `π`-mean. Exactly zero (`0.0` or a
zero matrix) for a class of one deme or whose `ℓ_d` are all equal. A rare mutant of effect
`a` is then selected against at `½aᵀ(s⁻¹ − 2Λ)a`. See src/local_alleles.jl.
"""
function local_gain(P::AbstractMatrix, ℓ::AbstractVector)
    k = length(ℓ)
    size(P) == (k, k) || throw(ArgumentError("P is $(size(P)); expected $k × $k"))
    z = _sym(zero(ℓ[1] * ℓ[1]'))
    out = Any[z for _ in 1:k]
    for comp in migration_classes(P)
        length(comp) > 1 || continue
        ls = ℓ[comp]
        all(l -> l == ls[1], ls) && continue
        Pc = Matrix{Float64}(P[comp, comp])
        Pc ./= sum(Pc; dims = 2)
        W = local_weights(Pc)
        A = Matrix{Float64}(I, length(comp), length(comp)) - Pc
        π = [A'; ones(1, length(comp))] \ [zeros(length(comp)); 1.0]
        lbar = sum(π[i] * ls[i] for i in eachindex(ls))
        dev = [l - lbar for l in ls]
        Λ = _sym(sum(W[i, j] * (dev[i] * dev[j]') for i in eachindex(dev), j in eachindex(dev)))
        for d in comp
            out[d] = Λ
        end
    end
    ℓ[1] isa Number ? Float64[x for x in out] : [Matrix{Float64}(x) for x in out]
end

# the positive part of a symmetric matrix (a number: max(x, 0))
psd_part(x::Real) = max(x, 0.0)
function psd_part(x::AbstractMatrix)
    F = eigen(Symmetric(Matrix{Float64}(x)))
    _sym(F.vectors * Diagonal(max.(F.values, 0.0)) * F.vectors')
end

"""
    QTLSim.local_sinv(sinv, Λ) -> s̃⁻¹

The inverse selection covariance a rare mutant feels in a deme of a metapopulation with
locally favoured alleles: the positive part of `sinv − 2Λ` ([`QTLSim.local_gain`](@ref)). For
`Λ = 0` it is `sinv` itself. See src/local_alleles.jl.
"""
local_sinv(sinv, Λ) = iszero(Λ) ? sinv : psd_part(_sym(sinv - 2 * Λ))
