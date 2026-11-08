# The multilocus Gaussian model: the recursion the moments lowering iterates by default
# (ext/QTLSimCatlabExt/moments.jl and trait_moments.jl) and its equilibria. The algebra is
# in lean/CategoricalQuantGen/Multilocus.lean.
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
# and their genic parts Tg, Rg, Qg, the same sums over pairs of alleles at one locus. The
# full moment minus the genic part is linkage disequilibrium (LD). Each is a number for
# one trait and a matrix (`(h − z̄/2)(h − z̄/2)ᵀ`) for several. With N adults,
#
#   C  = Q + (T + R − 2Q)/(2N)        (haploid: Q + (T − Q)/N), the expected square of the mean's deviation,
#   H  = T + R − 2C                   (haploid: T − C), the within-population covariance of a haplotype and the genotype,
#   V_G = 2H                          (haploid: H), the within-population genetic variance.
#
# Without drift N = ∞, so C = Q = 0, and V_G = 2T + 2R. One generation:
#
# - select, Gaussian of covariance s = Ω + P (one trait: s = ω² + V_P = V_s + V_G): the
#   survivors' allelic values are again Gaussian, each allele's regressed on the
#   phenotype. A haplotype's within-population variance falls by D = H s⁻¹ H, and so does
#   the covariance of an adult's two haplotypes (an adult's h₁ − h₂ does not change). By
#   exchangeability each locus carries 1/L of H, so its alleles' variance falls by
#   (H/L) s⁻¹ (H/L), and summed over the loci the genic parts fall by D/L. The rest of
#   the reduction, D(1 − 1/L), is negative LD, which recombination breaks up. With
#   B = V_G s⁻¹ the mean's deviation shrinks to (I − B)δ, so C ← (I − B) C (I − B)ᵀ, and
#   each locus's share by (B/L)δ, so Cg ← Cg − (B C + C Bᵀ)/L + B C Bᵀ/L.
# - reproduce, n parents (diploid: an offspring's two gametes come from distinct parents;
#   haploid: an offspring is a copy of one parent): R ← Q, Q ← (T + R)/(2n) + (1 − 1/n) Q
#   (haploid Q ← T/n + (1 − 1/n) Q), the genic parts alike, and recombination
#   T ← Tg + (1 − ρ)(T − Tg) + ρ (R − Rg), with ρ = `pair_recombination(p)` the probability
#   that a gamete's alleles at two random distinct loci come from different haplotypes of
#   its parent. A haploid offspring inherits its parent's haplotype whole, T ← T. Without
#   drift n = ∞.
# - mutate: T and Tg gain M/Ploidy, M = m·V_E (`muteVar`) per individual.
#
# Without drift the genic part of a diploid population's variance changes by
# −V_G²/(2Ls) + V_m per generation, so at equilibrium V_G²/(V_s + V_G) = 2L·V_m whatever
# the recombination (as long as ρ > 0), the fixed point `gaussian_vg_loci`: LD lowers V_G
# below the genic variance, and the genic variance rises until V_G itself balances
# mutation (Lande 1975). For one haploid locus, or haploid clones of any number of loci,
# all of selection's reduction is inherited and the recursion is `gaussian_vg_exact`'s
# `V ← V − V²/(V_s + V) + V_m`. Two simplifications: every pair of loci recombines at the
# mean rate ρ, and a nonlinear function of the random within-population variance is taken
# at its expected value (`D = E[H] s⁻¹ E[H]`, as is Gaussian selection with drift on).

# QTSim's fitness peak: pdf(Normal(0, ω), x)/(0.4/ω) = c·exp(−x²/2ω²)
const QTSIM_FITNESS_PEAK = 1 / (0.4 * sqrt(2π))

# a number's transpose is itself, so the formulas below serve one trait and several
_sym(x::Number) = x
_sym(x::AbstractMatrix) = (x + x') / 2
_eye(x::Number) = one(x)
_eye(x::AbstractMatrix) = Matrix{Float64}(I, size(x)...)

"""
    QTLSim.hap_init(V, N, ploidy)

The haplotype moments `(T, R, Q, Tg, Rg, Qg)` of a Hardy–Weinberg population in linkage
equilibrium with within-population genetic variance `V` (a number or a matrix) among `N`
adults (`Inf` without drift). See src/multilocus.jl.
"""
function hap_init(V, N::Real, ploidy::Integer)
    T = ploidy == 2 ? V ./ 2 : V .* 1
    Q = ploidy == 2 ? -T ./ (2N - 2) : -T ./ (N - 1)
    z = zero(T)
    (T, z, Q, copy(T), z, copy(Q))
end

hap_cmean(T, R, Q, N::Real, ploidy::Integer) = ploidy == 2 ? Q + (T + R - 2Q) / (2N) : Q + (T - Q) / N

"The within-population (divisor N) genetic variance of the haplotype moments `st`."
function hap_vpop(st, N::Real, ploidy::Integer)
    T, R, Q = st
    C = hap_cmean(T, R, Q, N, ploidy)
    ploidy == 2 ? 2T + 2R - 4C : T - C
end

"""
    QTLSim.hap_select(st, sinv, N, L, ploidy) -> (st′, B)

Gaussian selection of covariance `s` (`sinv` = s⁻¹) on the haplotype moments `st` of a
population of `N` adults at `L` loci; `B = V_G s⁻¹` is the factor by which the deviation of
the mean shrinks, for the moments with other populations. See src/multilocus.jl.
"""
function hap_select(st, sinv, N::Real, L::Integer, ploidy::Integer)
    T, R, Q, Tg, Rg, Qg = st
    C, Cg = hap_cmean(T, R, Q, N, ploidy), hap_cmean(Tg, Rg, Qg, N, ploidy)
    σ2, σ2g = T - C, Tg - Cg
    r, rg = ploidy == 2 ? (R - C, Rg - Cg) : (zero(R), zero(Rg))
    H = σ2 + r
    B = ploidy * H * sinv
    D = _sym(H * sinv * H)
    σ2, σ2g = σ2 - D, σ2g - D / L
    A = _eye(B) - B
    Cn = _sym(A * C * A')
    Cgn = _sym(Cg - (B * C + C * B') / L + B * C * B' / L)
    if ploidy == 2
        r, rg = r - D, rg - D / L
        return (Cn + σ2, Cn + r, Cn - (σ2 + r) / (2N - 2), Cgn + σ2g, Cgn + rg, Cgn - (σ2g + rg) / (2N - 2)), B
    end
    (Cn + σ2, R, Cn - σ2 / (N - 1), Cgn + σ2g, Rg, Cgn - σ2g / (N - 1)), B
end

"""
    QTLSim.hap_reproduce(st, n, ρ, ploidy)

Random mating of `n` parents (`Inf` without drift) with pair recombination probability
`ρ`: the haplotype moments of the offspring. See src/multilocus.jl.
"""
function hap_reproduce(st, n::Real, ρ::Real, ploidy::Integer)
    T, R, Q, Tg, Rg, Qg = st
    ploidy == 1 && return (T, R, T / n + (1 - 1 / n) * Q, Tg, Rg, Tg / n + (1 - 1 / n) * Qg)
    (Tg + (1 - ρ) * (T - Tg) + ρ * (R - Rg), Q, (T + R) / (2n) + (1 - 1 / n) * Q,
     Tg, Qg, (Tg + Rg) / (2n) + (1 - 1 / n) * Qg)
end

"Mutation adds `M/ploidy` to each haplotype, locus by locus."
hap_mutate(st, M, ploidy::Integer) = (st[1] + M / ploidy, st[2], st[3], st[4] + M / ploidy, st[5], st[6])

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
    multilocus_equilibrium(p; drift = false, tol = 1e-14, maxiter = 10^7)

The equilibrium of the multilocus Gaussian recursion of the moments lowering, found by
iterating one generation (selection, reproduction, mutation; see src/multilocus.jl) from
[`gaussian_vg_loci`](@ref)'s deterministic value, at a fixed optimum, until `V_G` changes
by less than `tol` relative to itself. With `drift = true` the population has `N =
PopSize` adults (the Beverton–Holt equilibrium size under `Demography = :beverton_holt`)
and `n` parents, `N W̄` survivors (`N W̄²/E[W²]` under fecundity selection). For several
traits the variances are matrices. Returns `(; VG, genic, ld, N, nparents, Wbar)`: `VG`
the within-population genetic variance (the expected sample variance, divisor `N − 1`,
with drift), `genic` its genic part and `ld = VG − genic` the contribution of linkage
disequilibrium.
"""
function multilocus_equilibrium(p::Params; drift::Bool = false, tol::Real = 1e-14, maxiter::Integer = 10^7)
    p.Drift && throw(ArgumentError("no equilibrium without selection (Drift = true)"))
    nt, L, h = p.NoTraits, p.NoQTL, p.Ploidy
    ρ = pair_recombination(p)
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
    st = hap_init(V0, drift ? N : Inf, h)
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
        st, _ = hap_select(st, sinv, Nd, L, h)
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
    g = hap_vpop(st[4:6], Nd, h)
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
