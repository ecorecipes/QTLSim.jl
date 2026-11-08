# The multilocus Gaussian model: the recursion the moments lowering iterates by default
# (ext/QTLSimCatlabExt/moments.jl and trait_moments.jl) and its equilibria. The algebra is
# in lean/CategoricalQuantGen/Multilocus.lean; the comparison with the kernel in
# validate/moments_loci.jl.
#
# A population of `NoQTL` additive loci, diploid or haploid, whose allelic values are
# jointly Gaussian (Lande 1975; Bulmer 1971, 1980; Turelli and Barton 1990, 1994). Its
# state is the expected second moments of the haplotypes about the deterministic mean,
# in units of one haplotype (one gamete's alleles summed over the loci):
#
#   T = E[(h − z̄/2)²]                 one haplotype (diploid; h − z̄ when haploid)
#   R = E[(h₁ − z̄/2)(h₂ − z̄/2)]       the two haplotypes of one adult (diploid only)
#   Q = E[(h − z̄/2)(h′ − z̄/2)]        haplotypes of two distinct adults
#
# each a sum over ordered pairs of loci (l, m) of the covariance of the alleles at l and
# m, and kept resolved by class of pairs (`pair_classes`): class 1, l = m, of weight
# ω₁ = 1/L, is the genic part; the other classes group the pairs of distinct loci by the
# probability ρ_c that a gamete takes their alleles from different haplotypes of its
# parent, with ω_c their share of the L² pairs. A class of distinct loci is linkage
# disequilibrium (LD). Each moment is a number for one trait and a matrix
# (`(h − z̄/2)(h − z̄/2)ᵀ`) for several. With N adults and T, R, Q summed over the classes,
#
#   C  = Q + (T + R − 2Q)/(2N)        (haploid: Q + (T − Q)/N), the expected square of the mean's deviation,
#   H  = T + R − 2C                   (haploid: T − C), the within-population covariance of a haplotype and the genotype,
#   V_G = 2H                          (haploid: H), the within-population genetic variance.
#
# Without drift N = ∞, so C = Q = 0, and V_G = 2T + 2R. One generation:
#
# - select, Gaussian of covariance s = Ω + P (one trait: s = ω² + V_P = V_s + V_G): the
#   survivors' allelic values are again Gaussian, each allele's regressed on the
#   phenotype through its covariance with the genotype. A haplotype's within-population
#   variance falls by D = H s⁻¹ H, and so does the covariance of an adult's two
#   haplotypes (an adult's h₁ − h₂ does not change). By exchangeability each locus
#   carries 1/L of H, so the pair (l, m) loses (H/L) s⁻¹ (H/L), and class c loses ω_c D:
#   the genic part D/L, the rest, D(1 − 1/L), negative LD. With B = V_G s⁻¹ the mean's
#   deviation shrinks to (I − B)δ, C ← C − ΔC with ΔC = B C + C Bᵀ − B C Bᵀ, and class c's
#   part of C by ω_c ΔC.
# - reproduce, n parents (diploid: an offspring's two gametes come from distinct parents;
#   haploid: an offspring is a copy of one parent): in every class R_c ← Q_c,
#   Q_c ← (T_c + R_c)/(2n) + (1 − 1/n) Q_c (haploid Q_c ← T_c/n + (1 − 1/n) Q_c), and
#   recombination T_c ← (1 − ρ_c) T_c + ρ_c R_c (ρ₁ = 0: a locus is never split from
#   itself). A haploid offspring inherits its parent's haplotype whole, T ← T. Without
#   drift n = ∞.
# - mutate: the genic class of T gains M/Ploidy, M = m·V_E (`muteVar`) per individual.
#
# Without drift the genic part of a diploid population's variance changes by
# −V_G²/(2Ls) + V_m per generation, so at equilibrium V_G²/(V_s + V_G) = 2L·V_m whatever
# the recombination (as long as every ρ_c > 0), the fixed point `gaussian_vg_loci`: LD
# lowers V_G below the genic variance, by L V_m Σ_c ω_c/ρ_c at equilibrium, and the
# genic variance rises until V_G itself balances mutation (Lande 1975). For one haploid
# locus, or haploid clones of any number of loci, all of selection's reduction is
# inherited and the recursion is `gaussian_vg_exact`'s `V ← V − V²/(V_s + V) + V_m`.
# Simplifications: the loci are exchangeable (equal shares of H; Lean
# `genic_reduction_ge`: unequal shares, as when a few large mutations hold the variance,
# lose more), the LD of a class is shared equally by its pairs, and a nonlinear function
# of the random within-population variance is taken at its expected value
# (`D = E[H] s⁻¹ E[H]`, as is Gaussian selection with drift on).

# QTSim's fitness peak: pdf(Normal(0, ω), x)/(0.4/ω) = c·exp(−x²/2ω²)
const QTSIM_FITNESS_PEAK = 1 / (0.4 * sqrt(2π))

# a number's transpose is itself, so the formulas below serve one trait and several
_sym(x::Number) = x
_sym(x::AbstractMatrix) = (x + x') / 2
_eye(x::Number) = one(x)
_eye(x::AbstractMatrix) = Matrix{Float64}(I, size(x)...)

"""
    QTLSim.pair_classes(p) -> (ω, ρ)

The classes of pairs of loci over which the multilocus model's haplotype moments are
resolved: class 1 is a locus with itself (the genic part, weight `ω₁ = 1/L`, never split,
`ρ₁ = 0`); the other classes group the ordered pairs of distinct loci by the probability
`ρ` that the kernel's crossover (`gamete!`, `CrossoverModel`) puts their alleles on
different haplotypes of a gamete, with `ω` the fraction of the `L²` ordered pairs in each.
Under `:qtsim` (the lowest crossover point only) the pair `(i, j)`, `i < j`, is split when
the lowest of `K ~ Poisson(mean)` distinct points in `0:L−2` lies in `i:j−1`; under `:all`
when an odd number of the points does, which depends only on `j − i`. Pairs whose `ρ`
differ by less than 2% share a class (their weighted mean `ρ`), so a class's LD,
which at equilibrium is proportional to `1/ρ`, is resolved to about 0.1%. The weighted
mean of `ρ` over the pairs of distinct loci is [`pair_recombination`](@ref).
"""
function pair_classes(p::Params)
    L = p.NoQTL
    L <= 1 && return ([1.0], [0.0])
    ρpairs, wpairs = pair_splits(L, p.RecombPoissonMean, p.CrossoverModel)
    # bins of 2% in ρ; pairs never split form a class of their own
    key(r) = r <= 0 ? typemin(Int) : floor(Int, log(r) / log(1.02))
    acc = Dict{Int,Tuple{Float64,Float64}}()
    for (r, w) in zip(ρpairs, wpairs)
        a = get(acc, key(r), (0.0, 0.0))
        acc[key(r)] = (a[1] + w, a[2] + w * r)
    end
    ks = sort(collect(keys(acc)))
    ω = [1 / L; [acc[k][1] for k in ks]]
    ρ = [0.0; [acc[k][2] / acc[k][1] for k in ks]]
    ω[2:end] .*= (1 - 1 / L) / sum(ω[2:end])          # exactly 1 − 1/L in all
    (ω, ρ)
end

# the split probability of every pair class of distinct loci and its weight among them:
# under :all by distance d (L − d pairs), under :qtsim every pair i < j
function pair_splits(L::Integer, mean::Real, model::Symbol)
    model in (:qtsim, :all) || throw(ArgumentError("model is :qtsim or :all; got :$model"))
    pois = Poisson(mean)
    Ks = [K for K in 1:L-1 if pdf(pois, K) > 1e-17 || K <= mean]
    if model === :all
        ρ = zeros(L - 1)
        for K in Ks, d in 1:L-1
            h = Hypergeometric(d, L - 1 - d, K)
            ρ[d] += pdf(pois, K) * sum(pdf(h, x) for x in 1:2:min(d, K); init = 0.0)
        end
        return ρ, [(L - d) / (L * (L - 1) / 2) for d in 1:L-1]
    end
    # P(min ≥ c) for K points, then pair (i, j) (0-based i < j) is split iff min ∈ i:j−1
    S = zeros(L)                                   # Σ_K pK P(min ≥ c), c = 0:L−1
    for K in Ks
        w = pdf(pois, K)
        for c in 0:L-1
            s = 1.0
            for i in 0:K-1; s *= max(L - 1 - c - i, 0) / (L - 1 - i); end
            S[c+1] += w * s
        end
    end
    ρ = Float64[]; wt = Float64[]
    np = L * (L - 1) / 2
    for i in 0:L-2, j in i+1:L-1
        push!(ρ, S[i+1] - S[j+1]); push!(wt, 1 / np)
    end
    ρ, wt
end

"""
    QTLSim.hap_init(V, N, ploidy, nclass)

The haplotype moments `(T, R, Q)` (each a vector over the `nclass` pair classes of
[`pair_classes`](@ref)) of a Hardy–Weinberg population in linkage equilibrium with
within-population genetic variance `V` (a number or a matrix) among `N` adults (`Inf`
without drift): all of it genic, class 1. See src/multilocus.jl.
"""
function hap_init(V, N::Real, ploidy::Integer, nclass::Integer = 2)
    T0 = ploidy == 2 ? V ./ 2 : V .* 1
    Q0 = ploidy == 2 ? -T0 ./ (2N - 2) : -T0 ./ (N - 1)
    z = zero(T0)
    ([k == 1 ? T0 : z for k in 1:nclass], [z for _ in 1:nclass], [k == 1 ? Q0 : z for k in 1:nclass])
end

hap_cmean(T, R, Q, N::Real, ploidy::Integer) = ploidy == 2 ? Q + (T + R - 2Q) / (2N) : Q + (T - Q) / N

"The within-population (divisor N) genetic variance of the haplotype moments `st`."
function hap_vpop(st, N::Real, ploidy::Integer)
    T, R, Q = sum(st[1]), sum(st[2]), sum(st[3])
    C = hap_cmean(T, R, Q, N, ploidy)
    ploidy == 2 ? 2T + 2R - 4C : T - C
end

"The genic part (class 1) of the within-population genetic variance."
hap_genic(st, N::Real, ploidy::Integer) = hap_vpop(([st[1][1]], [st[2][1]], [st[3][1]]), N, ploidy)

"""
    QTLSim.hap_select(st, sinv, N, ω, ploidy; α2 = 0) -> (st′, B)

Gaussian selection of covariance `s` (`sinv` = s⁻¹) on the haplotype moments `st` of a
population of `N` adults, with pair-class weights `ω` ([`pair_classes`](@ref)); `B = V_G s⁻¹`
is the factor by which the deviation of the mean shrinks, for the moments with other
populations. See src/multilocus.jl. With `α2 = α²` > 0 (one trait), the rare-alleles
closure of src/rare_alleles.jl: the genic part also loses `α²/(2s)` of itself.
"""
function hap_select(st, sinv, N::Real, ω::AbstractVector, ploidy::Integer; α2::Real = 0.0)
    T, R, Q = st
    Tf, Rf, Qf = sum(T), sum(R), sum(Q)
    C = hap_cmean(Tf, Rf, Qf, N, ploidy)
    H = ploidy == 2 ? Tf + Rf - 2C : Tf - C
    B = ploidy * H * sinv
    D = _sym(H * sinv * H)
    ΔC = _sym(B * C + C * B' - B * C * B')          # (I − B) C (I − B)ᵀ = C − ΔC
    Tn, Rn, Qn = similar(T), similar(R), similar(Q)
    for k in eachindex(T)
        Ck = hap_cmean(T[k], R[k], Q[k], N, ploidy)
        Cn = Ck - ω[k] * ΔC
        σ2 = T[k] - Ck - ω[k] * D
        # rare alleles (α2 = α², `loci = :rare_alleles`): a locus's alleles have fourth
        # cumulant α² times their variance, so the genic part also loses β²α²/(2s) of the
        # within-population genic variance, β = H/H_genic the regression of the genotype
        # on a locus's alleles (LD shrinks a mutant's effective effect; src/rare_alleles.jl)
        if k == 1 && α2 > 0
            Hg = ploidy == 2 ? T[k] + R[k] - 2Ck : T[k] - Ck      # the genic part of H
            β2 = Hg > 0 ? (H / Hg)^2 : 1.0
            σ2 -= min(β2 * α2 * sinv / 2, 1.0) * (T[k] - Ck)
        end
        if ploidy == 2
            r = R[k] - Ck - ω[k] * D
            Tn[k], Rn[k], Qn[k] = Cn + σ2, Cn + r, Cn - (σ2 + r) / (2N - 2)
        else
            Tn[k], Rn[k], Qn[k] = Cn + σ2, R[k], Cn - σ2 / (N - 1)
        end
    end
    (Tn, Rn, Qn), B
end

"""
    QTLSim.hap_reproduce(st, n, ρ, ploidy)

Random mating of `n` parents (`Inf` without drift), with `ρ` the split probabilities of
the pair classes: the haplotype moments of the offspring. See src/multilocus.jl.
"""
function hap_reproduce(st, n::Real, ρ::AbstractVector, ploidy::Integer)
    T, R, Q = st
    ploidy == 1 && return (T, R, [T[k] / n + (1 - 1 / n) * Q[k] for k in eachindex(T)])
    ([(1 - ρ[k]) * T[k] + ρ[k] * R[k] for k in eachindex(T)], copy(Q),
     [(T[k] + R[k]) / (2n) + (1 - 1 / n) * Q[k] for k in eachindex(T)])
end

"Mutation adds `M/ploidy` to each haplotype, locus by locus (the genic class)."
function hap_mutate(st, M, ploidy::Integer)
    T = copy(st[1]); T[1] = T[1] + M / ploidy
    (T, st[2], st[3])
end


# the number of allele positions whose variance balances mutation separately: 2L for a
# diploid that recombines, 2 for a diploid without recombination (one locus, or none),
# 1 for a haploid (clones)
function gaussian_units(p::Params)
    p.Ploidy == 1 && return 1
    pair_recombination(p) > 0 ? 2p.NoQTL : 2
end

"""
    gaussian_vg_loci(p; drift = false)

The equilibrium genetic variance of the multilocus Gaussian model, the fixed point that
the moments lowering (`run_moments`, by default) converges to. Without drift, the
positive root of `V²/(V_s + V) = n·V_m`,

    V = (n V_m + √(n² V_m² + 4 n V_m V_s))/2,

with `V_s` = [`selection_variance`](@ref)`(p)`, `V_m` = [`muteVar`](@ref)`(p)` and `n = 2L`
for a diploid of `L` loci that recombine (`n = 2` without recombination, `n = 1` for a
haploid, whose offspring are clones). Each allele position balances mutation against the
part of selection's reduction that its share of the variance inherits, `V²/(n s)`, while
the rest is linkage disequilibrium, which recombination breaks up; at equilibrium the
genic variance exceeds `V` by the LD, and `V` itself does not depend on the rate of
recombination (Lande 1975; Bulmer 1980; Turelli and Barton 1990). It is
[`gaussian_vg_exact`](@ref) with `n·V_m` in place of `V_m`, and `√(n V_m V_s)` for weak
mutation (Lande's per-locus sum).

With `drift = true`, the stationary expected within-population sample variance of the
same recursion with drift among `N = PopSize` adults and their parents (survivors, or
`N W̄²/E[W²]` under fecundity selection; [`multilocus_equilibrium`](@ref)), which is what
`run_moments(p; erosion = true)` converges to and the kernel's `geno_var` estimates.

The Gaussian assumption fails when few mutations of large effect hold the variance at a
locus (`α² = V_m/(2Lμ)` large against the locus's share `V/(2L)`, the house-of-cards
regime, [`house_of_cards_vg`](@ref)); there the kernel holds less variance than either.
"""
gaussian_vg_loci(p::Params; drift::Bool = false) =
    drift ? multilocus_equilibrium(p; drift).VG : gaussian_vg_closed(gaussian_units(p) * muteVar(p), selection_variance(p))

gaussian_vg_closed(Vm, Vs) = (Vm + sqrt(Vm^2 + 4Vm * Vs)) / 2

"""
    multilocus_equilibrium(p; drift = false, loci = :multilocus, tol = 1e-14, maxiter = 10^7)

The equilibrium of the multilocus Gaussian recursion of the moments lowering, found by
iterating one generation (selection, reproduction, mutation; see src/multilocus.jl) from
[`gaussian_vg_loci`](@ref)'s deterministic value, at a fixed optimum, until `V_G` changes
by less than `tol` relative to itself. With `drift = true` the population has `N =
PopSize` adults (the Beverton–Holt equilibrium size under `Demography = :beverton_holt`)
and `n` parents, `N W̄` survivors (`N W̄²/E[W²]` under fecundity selection). For several
traits the variances are matrices. Returns `(; VG, genic, ld, N, nparents, Wbar)`: `VG`
the within-population genetic variance (the expected sample variance, divisor `N − 1`,
with drift), `genic` its genic part and `ld = VG − genic` the contribution of linkage
disequilibrium. With `loci = :rare_alleles` (one trait), the recursion of
`run_moments(p; loci = :rare_alleles)`, with the rare-alleles closure of
src/rare_alleles.jl ([`rare_alleles_vg`](@ref)).
"""
function multilocus_equilibrium(p::Params; drift::Bool = false, loci::Symbol = :multilocus, tol::Real = 1e-14, maxiter::Integer = 10^7)
    p.Drift && throw(ArgumentError("no equilibrium without selection (Drift = true)"))
    nt, L, h = p.NoTraits, p.NoQTL, p.Ploidy
    α2 = rare_alleles_effect(p, loci)
    ω, ρ = pair_classes(p)
    VE = envVar(p)
    if nt == 1
        M, Ω, Id = muteVar(p), initFitSD(p)^2, 1.0
        V0 = gaussian_vg_loci(p)
    else
        M, Ω = mutational_covariance(p), selection_covariance(p)
        Id = Matrix{Float64}(I, nt, nt)
        V0 = gaussian_G_loci(p)
    end
    N = float(p.PopSize)
    st = hap_init(V0, drift ? N : Inf, h, length(ω))
    VG, Wbar, n = V0, 1.0, Inf
    c = QTSIM_FITNESS_PEAK
    for _ in 1:maxiter
        Nd = drift ? max(N, 2.0) : Inf
        V = hap_vpop(st, Nd, h)
        P = V + VE * Id
        if nt == 1
            ωe2 = p.StandardisePhenotype ? Ω * P : Ω
            s = ωe2 + P
            sinv = 1 / s
            Wbar = c * sqrt(ωe2 / s)
            W2 = c^2 * sqrt(ωe2 / (ωe2 + 2P))
        else
            s = Ω + P
            sinv = inv(s)
            Wbar = c * sqrt(det(Ω) / det(s))
            W2 = c^2 * sqrt(det(Ω) / det(Ω + 2P))
        end
        st, _ = hap_select(st, sinv, Nd, ω, h; α2)
        npar = p.Selection === :fecundity ? N * Wbar^2 / W2 : N * Wbar
        n = drift ? max(npar, 1.0) : Inf
        st = hap_reproduce(st, n, ρ, h)
        st = hap_mutate(st, M, h)
        nw = N * Wbar
        N = p.Demography === :fixed ? float(p.PopSize) : min(p.PopSize, p.GrowthRate * nw / (1 + (p.GrowthRate - 1) * nw / p.PopSize))
        N < 1 && throw(ArgumentError("the population cannot persist at the optimum"))
        Vn = hap_vpop(st, drift ? max(N, 2.0) : Inf, h)
        done = maximum(abs, Vn - VG) <= tol * maximum(abs, Vn)
        VG = Vn
        done && break
    end
    Nd = drift ? max(N, 2.0) : Inf
    g = hap_genic(st, Nd, h)
    k = drift ? Nd / (Nd - 1) : 1.0
    (VG = VG * k, genic = g * k, ld = (VG - g) * k, N, nparents = n, Wbar)
end

"""
    gaussian_G_loci(p; drift = false)

The equilibrium G matrix of the multilocus Gaussian model for several traits, the fixed
point of the trait moments lowering (`run_moments`, by default): without drift the
positive semidefinite solution of `G (A + G)⁻¹ G = n M`, [`gaussian_G`](@ref)`(n M, A)`,
with `n` allele positions as in [`gaussian_vg_loci`](@ref) (`2L` for a diploid that
recombines). Each position's share of `G`, `G/n`, loses `(G/n)(A + G)⁻¹(G/n)` to selection
for good, and the rest of selection's reduction is linkage disequilibrium that
recombination breaks up. For one trait it is `gaussian_vg_loci`; for weak mutation `G ≈
A # (n M)`, the matrix geometric mean. With `drift = true`, the stationary expected sample
covariance with drift ([`multilocus_equilibrium`](@ref)).
"""
gaussian_G_loci(p::Params; drift::Bool = false) =
    drift ? multilocus_equilibrium(p; drift).VG : gaussian_G(gaussian_units(p) .* mutational_covariance(p), selection_variance_matrix(p))
