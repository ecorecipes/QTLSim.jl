# Mutation–selection–drift balance across regimes: the rare-alleles closure of the moments
# lowering (`run_moments(p; loci = :rare_alleles)`, the default for one diploid trait,
# `moments_loci`), its equilibrium, and the closed form
# it reduces to for loci in linkage equilibrium, Hermisson and Wagner's (2004) "stochastic
# house of Gauss". The comparison with the kernel is validate/moments_regimes.jl.
#
# The per-locus view. Take one allele position (a gene: L of them per haplotype) with
# allelic values x of within-population variance c and fourth central moment m₄. Under
# weak Gaussian selection of width s on the genotype (s = ω² + V_P, so V_s + V_G), an
# allele's marginal fitness is quadratic in its deviation, 1 − (x − x̄)²/(2s), and
# selection changes c by −(m₄ − c²)/(2s) (Turelli and Barton 1990, the weak-selection
# moment recursion). Mutation adds μα² per gene per generation, α² = V_m/(Ploidy·L·μ) the
# variance of one mutation's effect (`mutation_effect_variance`; QTSim's MutVar), and drift
# among 2N genes removes c/(2N_e). The Gaussian closure m₄ = 3c² gives Δc = −c²/s, Lande's
# per-locus recursion and the multilocus model of src/multilocus.jl, which is right when a
# locus carries many mutations' worth of variance (c ≫ α²). When mutations are rare and
# large, a locus is mostly one allele with a few rare mutants, and m₄ is dominated by them:
# a mutant of effect a at frequency q contributes q a² to c and q a⁴ to m₄, and at the
# house-of-cards balance q(a) = μ g(a)·2s/a² (g the law of effects), so m₄/c = E_g[a²] = α²
# exactly. The rare-alleles closure keeps both parts,
#
#   m₄ = 3c² + α² c         (fourth cumulant κ₄ = α² c),
#
# so that selection changes the genic variance by
#
#   Δc = −c²/s − α² c/(2s),
#
# the Gaussian loss plus a loss of α²/(2s) of the variance itself: each rare mutant of
# effect a is removed at its own rate a²/(2s), whatever the variance around it. The
# per-gene recursion is then
#
#   Δc = μα² − c/(2N_e) − α² c/(2s) − c²/s,
#
# whose fixed point, with u = μ, v = α², Hermisson and Wagner's s = 1/(2V_s) and their
# N_e the number of genes (2N_e diploid), is exactly their equation A4,
#
#   c/v = (√(s_e² + 2(1 + Θ)s_e + 1) − s_e − 1)/(4 s_e),   s_e = 2N_e α²/(2V_s),  Θ = 4·2N_e μ,
#
# the "stochastic house of Gauss" (their equation 14 is an algebraic variant): with
# V_G = 2L c it is Lynch and Hill's (1986) neutral 2N_e V_m when selection is weak
# against drift (s_e → 0), Bürger's stochastic house of cards
# 4Lμ V_s/(1 + V_s/(N_e α²)) (Bürger, Wagner and Stettinger 1989; Bürger and Lande 1994)
# when Θ s_e ≪ (1 + s_e)², hence Turelli's (1984) house of cards 4Lμ V_s for large N_e,
# and the Gaussian √(2L V_m V_s) (Lande 1975; Kimura 1965) when Θ s_e ≫ (1 + s_e)².
# `house_of_gauss_vg` is that closed form (loci in linkage equilibrium, s = V_s).
#
# The stochastic house of cards is more than an interpolation here: for independent
# mutations (infinite sites) with Gaussian effects, the diffusion of one mutant of effect
# a under drift and the selection a locus feels at the optimum, Δp = −S pq(1 − 2p) with
# S = a²/(2s) (selection against the rarer allele), gives an expected lifetime Σ pq of
# h(γ) = 2∫₀¹ e^{−γpq} ∫ₚ¹ e^{γy(1−y)} dy dp / ∫₀¹ e^{γy(1−y)} dy, γ = 4N_e S, and the
# average over a ~ N(0, α²) of a² h(γ_a) is α²/(1 + 2N_e α²/(2s)) to quadrature precision
# (checked numerically at γ_α = 0.5–50 while building this model; not in Lean), which is
# Bürger's form. The rare-alleles recursion has the same fixed point in that regime.
#
# In the moments lowering, `loci = :rare_alleles` adds the loss to the genic class of the
# multilocus model (`hap_select`'s `α2`), so it keeps linkage disequilibrium,
# recombination by pair class, drift through the haplotype moments, s = V_s + V_G, moving
# optima, demography and demes. As in the Gaussian part, where a locus loses
# (H/L) s⁻¹ (H/L) with H/L its alleles' covariance with the genotype rather than their
# genic variance, a rare mutant is selected through its average effect on the genotype,
# β a with β = H/H_genic the regression of the genotype on a locus's alleles (LD included,
# exchangeable loci); the genic class loses β² α²/(2s) of its within-population variance,
# with β = 1 in linkage equilibrium. No constant is fitted: α² is the kernel's mutation
# effect variance.
#
# Simplifications: selection is weak per mutation (a²/(2s) ≪ 1; the loss is capped at the
# whole variance); κ₄ = α² c is exact at the deterministic house-of-cards balance, and
# negligible against 3c² in the Gaussian regime, but only approximate in between and under
# drift, where the alleles at a locus are a few clades of shared mutations; the mean's
# response to directional selection ignores the skew of the alleles (Barton and Turelli 1987).

"""
    QTLSim.mutation_effect_variance(p) -> α²

The variance of one mutation's effect on the trait, `α² = V_m/(Ploidy·L·μ)`, with `V_m` =
[`muteVar`](@ref)`(p)`, `L = NoQTL` and `μ = MutationRatePerGamete` per locus per gamete
(QTSim's `MutVar`, the kernel's `mutate!`).
"""
function mutation_effect_variance(p::Params)
    p.MutationRatePerGamete > 0 || throw(ArgumentError("no mutations (MutationRatePerGamete = 0): α² is undefined"))
    muteVar(p) / (p.Ploidy * p.NoQTL * p.MutationRatePerGamete)
end

"""
    QTLSim.moments_loci(p, loci = :auto) -> Symbol

The model of the loci that the moments lowering ([`run_moments`](@ref)) and
[`multilocus_equilibrium`](@ref) use when asked for `loci`: any value but `:auto` is
returned as it is, and `:auto`, the default, resolves to
- `:rare_alleles`, the rare-alleles closure (src/rare_alleles.jl), for one trait, diploid or
  haploid (clones; src/rare_alleles_traits.jl, within 4% of the kernel in all 8 haploid
  cells of `validate/rare_alleles_traits.jl` on two seed sets, up to 3000× closer than the
  Gaussian). Diploid:
  within 15% of the kernel in 46 of 47 weak-selection cells of
  `validate/moments_regimes.jl`, never farther from it than the multilocus Gaussian, and
  up to 10× closer in the house-of-cards and drift regimes, to which it reduces to the
  Gaussian when mutations are many and small;
- `:rare_local`, the closure for several traits with locally favoured alleles
  (src/rare_alleles_traits.jl, src/local_alleles.jl), when `NoTraits > 1`: for one
  population it is `:rare_alleles`, far better than the multilocus Gaussian in the
  house-of-cards regime (`validate/rare_alleles_traits.jl`, T cells); for demes with local
  optima and migration it is within 15% of the kernel's within-deme G where `:rare_alleles`
  gives under half of it (`validate/local_adaptation.jl`, D1; D2 too).
`run_moments` resolves `:auto` to `:multilocus` when the `mutate` steps of the life
cycle carry different `MutationScalar` or `MutationRatePerGamete`, as no one `α²` then
describes the mutations; when they all carry the same, `α²` is theirs. The demes of a
metapopulation each have their own `α²`.
"""
function moments_loci(p::Params, loci::Symbol = :auto)
    loci === :auto || return loci
    p.NoTraits > 1 && return :rare_local
    :rare_alleles
end

# α² for the rare-alleles closure, 0 for the multilocus Gaussian model
function rare_alleles_effect(p::Params, loci::Symbol)
    loci = moments_loci(p, loci)
    loci in (:multilocus, :rare_alleles, :rare_local) ||
        throw(ArgumentError("loci is :auto, :multilocus, :rare_alleles, :rare_local or :aggregate; got :$loci"))
    loci === :multilocus && return 0.0
    p.NoTraits == 1 || return mutation_effect_covariance(p)      # several traits: src/rare_alleles_traits.jl
    mutation_effect_variance(p)
end

"""
    house_of_gauss_vg(p; nsurv = p.PopSize, drift = true)

The equilibrium genetic variance of `L` loci in linkage equilibrium under mutation,
stabilising selection and drift, interpolating the neutral, house-of-cards and Gaussian
regimes without a fitted constant: Hermisson and Wagner's (2004, equation A4) "stochastic
house of Gauss", `V_G = Ploidy·L·c` with each gene's variance `c` the positive root of

    c²/V_s + c (α²/(2V_s) + 1/(Ploidy·N_e)) = μα²,

i.e. `c = 2μα²/(b + √(b² + 4μα²/V_s))`, `b = α²/(2V_s) + 1/(Ploidy·N_e)`, with `μ =
MutationRatePerGamete`, `α²` = [`QTLSim.mutation_effect_variance`](@ref)`(p)`, `V_s` =
[`selection_variance`](@ref)`(p)` and `N_e` from [`effective_size`](@ref) with `nsurv`
survivors (`N_e = ∞` with `drift = false`). A haploid genome, inherited whole by its clones
(or a diploid's haplotype when no pair of loci recombines), is one unit with `Lμ` in place
of `μ` and `V_G = Ploidy·c` (src/rare_alleles_traits.jl). Its limits are Lynch and Hill's neutral
`2N_e V_m` (`V_s → ∞`), Bürger's stochastic house of cards [`house_of_cards_vg`](@ref)
(rare mutations, `μ → 0`), Turelli's house of cards `4Lμ V_s` (rare mutations, no drift)
and the Gaussian `√(2L V_m V_s)` (many mutations, [`gaussian_vg_loci`](@ref)'s weak-mutation
form). It is the fixed point of the per-gene recursion behind the rare-alleles closure
(src/rare_alleles.jl); [`rare_alleles_vg`](@ref) is the same closure in the multilocus
model, with linkage disequilibrium and `s = V_s + V_G`.
"""
function house_of_gauss_vg(p::Params; nsurv::Real = p.PopSize, drift::Bool = true)
    α2, Vs = mutation_effect_variance(p), selection_variance(p)
    nu, μ = rare_units(p)          # 2L genes at μ; a haploid clone (or a haplotype that never recombines) at Lμ
    ng = drift ? p.Ploidy * effective_size(p.PopSize, nsurv) : Inf
    b = α2 / (2Vs) + 1 / ng
    c = 2μ * α2 / (b + sqrt(b^2 + 4μ * α2 / Vs))
    nu * c
end

"""
    rare_alleles_vg(p; drift = false)

The equilibrium genetic variance of the moments lowering with `loci = :rare_alleles`,
the default for one diploid trait ([`QTLSim.moments_loci`](@ref)):
`run_moments(p)` converges to it (with `drift = true`, `run_moments(p; erosion = true)`). The
multilocus model of [`gaussian_vg_loci`](@ref) with the rare-alleles closure, in which the
genic variance also loses `β²α²/(2s)` of itself to selection per generation, `α²` =
[`QTLSim.mutation_effect_variance`](@ref)`(p)`, `s = V_s + V_G` and `β = V_G/V_genic` the
regression of the genotype on a locus's alleles (`β = 1` without LD; src/rare_alleles.jl). It tends to
[`gaussian_vg_loci`](@ref) when mutations are many and small and to the (stochastic) house
of cards when they are rare and large; for loci in linkage equilibrium it is
[`house_of_gauss_vg`](@ref). See [`multilocus_equilibrium`](@ref) for the details.
"""
rare_alleles_vg(p::Params; drift::Bool = false) = multilocus_equilibrium(p; drift, loci = :rare_alleles).VG
