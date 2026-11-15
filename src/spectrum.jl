# The allele-frequency spectrum of the segregating mutations: the moments lowering with
# `loci = :spectrum` (ext/QTLSimCatlabExt/spectrum.jl), a model of the genetic variance
# that follows how directional selection, stabilising selection and drift move each
# mutation's frequency, so that the variance can rise when a moving optimum drives rare
# alleles of large effect to intermediate frequency (Barton and Turelli 1987; Bürger 1999)
# and fall when selection fixes or purges the alleles that carried it faster than drift
# would, while the population, and with it the number of parents, shrinks. The comparison
# with the kernel is validate/vg_move.jl.
#
# The state is the expected number of segregating mutations of effect a at frequency p,
# φ(a, p), the Poisson random field of Sawyer and Hartl (1992): every mutation is a new
# site (infinite sites, which the continuum-of-alleles kernel approaches when a locus
# rarely carries more than a few segregating mutations, 4Nμ ≲ 1 per locus) in linkage
# equilibrium with the others. A diploid's genotype is Σ 2a p over the sites, so
#
#   z̄ = z̄₀ + Σ 2a p,   V_G = Σ 2a² p(1 − p)    (Hardy–Weinberg, linkage equilibrium).
#
# One generation, with the population's lag λ = θ − z̄ and S = ω_e² + V_P (s⁻¹ of the
# lowering, `Moments.sinv`):
#
# - select: a mutant's frequency changes by Wright's formula Δp = (pq/2) ∂ln W̄/∂p, with
#   ln W̄ = −λ²/(2S) − ln(S)/2 for Gaussian selection (W̄ ∝ √(ω_e²/S) exp(−λ²/2S)):
#
#     E[Δp] = pq [a λ/S − (1 − λ²/S)(a² + 2c)(1 − 2p)/(2S)],
#
#   Hayward and Sella (2022, eLife 11:e66697, eq. 7, with V_S = S) without the 2c: a
#   directional (additive) term of selection coefficient aλ/S and a stabilising
#   (underdominant) one, which near p = 0 is Kopp and Hermisson's (2009, Genetics
#   183:1453, eq. 7) s(α) = α(2δ − α)/(2V_s), the selection coefficient of a new mutation
#   behind a moving optimum. The extra 2c, with c = V_G/(2L) the variance of one gene
#   (locus × haplotype), is the Gaussian part of the rare-alleles closure
#   (src/rare_alleles.jl): there a gene of variance c loses c²/s + α²c/(2s) per
#   generation, and giving each of its mutations the underdominance (a² + 2c)/(2S)
#   reproduces both terms when the mutations are rare, so that at a fixed optimum the
#   spectrum's equilibrium is Hermisson and Wagner's (2004) stochastic house of Gauss,
#   `house_of_gauss_vg`, which reduces to the neutral, house-of-cards and Gaussian
#   balances. Selection moves the mean by E[Δz̄] = Σ 2a E[Δp] = V_G λ/S − (skew terms):
#   the alleles' third moment, which the Gaussian response V_G λ/S ignores, enters
#   through (1 − 2p).
# - reproduce: drift among n parents (the expected survivors; N W̄²/E[W²] under fecundity
#   selection), Var(Δp) = pq/(2n), as in the haplotype moments of src/multilocus.jl,
#   times N/N_e = 1 + b² V_G Σ_c ω_c/(2(ρ_c + δ)²) under directional selection of
#   gradient b = λ/S: linked selection (Robertson 1961; Barton 1995, Genetics 140:821,
#   N/N_e = 1 + var(W)/(2r²) for a site at recombination r from a locus of additive
#   variance in relative fitness var(W), here b² V_G/L per locus), summed over the pair
#   classes of distinct loci of the kernel's crossover (`pair_classes`), with drift
#   and selection dissolving the association at 1/(2n) + V_G/S (Santiago and Caballero
#   1998: recombination, drift and the reduction of the selected variance by selection). `spectrum_linked`.
# - mutate: 2N L μ new mutations (N offspring, L loci, μ per locus per gamete), each at
#   frequency 1/(2N), of effect a ~ N(0, α²), α² = V_m/(2Lμ).
#
# Selection and drift are the Wright–Fisher diffusion ∂φ/∂t = −∂(Mφ)/∂p + ∂²(Dφ)/∂p²,
# M = E[Δp], D = Var(Δp)/2 (Kimura 1964), which a generation advances by one implicit Euler
# step on a Chebyshev grid of frequencies (fine near 0, where mutations enter, and near 1,
# where sweeps end), with Scharfetter–Gummel (Chang–Cooper) fluxes, exact for a locally
# constant M/D = 4n s(p) and positivity-preserving; sites reaching 0 or 1 leave the
# variance (lost or fixed). Effects a are Gauss–Hermite nodes of N(0, α²). Nothing is
# fitted to the kernel: α², L, μ, N, n, S, λ and the crossover's pair classes are the model's.
#
# Under a moving optimum this is the origin–fixation picture of Kopp and Hermisson in a
# population: a mutation of effect a fixes with probability ≈ 2aλ/S (Haldane) and, during
# its sweep, contributes ∫2a²pq dt = 2a/(aλ/S) to V_G, so that the sweeps' stationary
# variance is ≈ 2N_e V_m whatever the lag (Hill 1982, PNAS 79:142: asymptotic response
# 2N_e V_m times the selection gradient), but it turns over on the time scale of a sweep,
# ~ln(4N s)/s with s = aλ/S, not of drift, 2N_e: when the lag grows and the survivors,
# hence N_e, fall, V_G follows within a few hundred generations, which lowers the response,
# raises the lag and lowers N_e further. With little recombination (QTSim applies at most
# one crossover per gamete) the linked-selection term lowers N_e further as b² V_G grows,
# and the variance collapses faster still; with free recombination it vanishes. The haplotype moments of src/multilocus.jl and
# src/rare_alleles.jl erode V_G only at 1/(2N_e) per generation, and their directional
# selection moves the mean without using up variance.
#
# Dynamic linkage disequilibrium (`ld = :dynamic`, off by default; validate/spectrum_ld_dev.log,
# validate/spectrum_ld_criteria.md). Stabilising selection with tight linkage builds negative
# LD (Bulmer 1971), which the spectrum above leaves out (V_G = genic). With the option,
# V_G = g + D, g the spectrum's genic variance and D = Σ_c d_c the LD of classes c of pairs of
# distinct loci, grouped by the bins (a, b) of `ldbins` bins of consecutive positions that
# hold their loci and by their split probability ρ_c (`ld_classes`; ω_c their share of the L²
# ordered pairs). One generation, from the multilocus model's class recursion (src/multilocus.jl)
# with the locus's covariance with the genotype resolved by position:
#
#   d_c ← (1 − 1/(2n)) [(1 − ρ_c)(1 + 2γ) d_c − ω_c β_a β_b g²/(2S)],
#
# Bulmer's generation (a pair of loci gains −V_iV_j/S in the gametes, V_i = β_b g/L its
# covariance with the genotype), recombination, drift among n parents (Hill and Robertson
# 1968), and γ, the relative change of g by selection: under selection without epistasis
# Δln D_ij = Δln p_iq_i + Δln p_jq_j, so the LD of a pair follows its sites' variances
# (mean field over effects). β_b = V_b/g_b, V_b the bin's genic share plus the LD of its
# pairs. Each site's selection uses the regression of the genotype on it: directional
# β̄ aλ/S with β̄ = Σ_b m_b β_b = V_G/g (the response is V_G λ/S), stabilising
# β²‾ (a² + 2c)(1 − 2p)/(2S) with β²‾ = Σ_b m_b β_b² (the multilocus genic loss (H/L)²/s and
# the rare-alleles closure's β²). Mutations enter without LD. Without drift and with γ = 0,
# one bin and the infinitesimal limit, the fixed point is the multilocus model's: V_G =
# gaussian_vg_loci, D = −L V_m Σ_c ω_c/ρ_c. Resolving the LD by position matters under
# QTSim's crossover (only the lowest point switches): the tightly linked pairs lie toward the
# chromosome's end, where the loci's covariance with the genotype is reduced by their own LD,
# so LD is generated most weakly where it decays most slowly; the mean field over loci (the
# multilocus model, `ldbins = 1`) makes the LD too strong (dev log, step 2).
# Not modelled: the frequency distribution of sites by position (one spectrum for all bins),
# LD between sites of one locus, LD heterogeneity by effect size (one γ), the random part of
# the LD (left to Barton's term, which is kept as in the default: the expected LD D and
# Barton's extra variance of Δp are different moments), and the early fall of the genic
# variance under a move that the kernel shows (dev log, step 5: the option does not explain it).
#
# Simplifications: linkage equilibrium apart from the linked-selection term (no Bulmer
# LD unless `ld = :dynamic`; Barton's loose-linkage N_e for every pair class, and none for the mutations of the
# same locus); infinite sites (two
# mutations of one locus are separate sites, not one allele); weak selection per
# mutation (Wright's formula); the expected spectrum along the deterministic trajectory
# of λ and n (the nonlinear feedback through λ and n is taken at the expectation, as in
# the rest of the moments lowering).

"""
    QTLSim.LDClasses

The classes of unordered pairs of distinct loci over which the linkage disequilibrium of the
spectrum (`ld = :dynamic`, src/spectrum.jl) is resolved: the loci are cut into `length(m)`
bins of consecutive positions (`m[b]` the fraction of the loci in bin `b`), and the pairs
with one locus in bin `a[c]` and the other in bin `b[c]` (`a[c] ≤ b[c]`) are grouped by the
probability `ρ[c]` that the kernel's crossover splits them (as [`QTLSim.pair_classes`](@ref):
2% bins of `ρ`, their weighted mean), `ω[c]` the fraction of the `L²` ordered pairs in the
class. With one bin the classes are `pair_classes(p)`'s of distinct loci. See [`ld_classes`](@ref).
"""
struct LDClasses
    a::Vector{Int}
    b::Vector{Int}
    ρ::Vector{Float64}
    ω::Vector{Float64}
    m::Vector{Float64}
end

"""
    QTLSim.AlleleSpectrum

The expected allele-frequency spectrum of the segregating mutations of a diploid
population (src/spectrum.jl): `n[j, i]` the expected number of mutations of effect `a[j]`
at frequency `p[i]` (nodes `p[1] = 0` and `p[end] = 1` are the absorbing boundaries and
stay empty). `b`, `κ` and `S` are the selection the next `spectrum_step!` applies.

`ldratio` (keyword of the constructor and of [`spectrum_equilibrium`](@ref), through
`run_moments(…; spectrum = (; ldratio))`; default 1, the model as described above) is an
exploratory prototype of validate/spectrum_diagnosis.jl, not a validated option: a factor
`r = V_G/genic` for the negative linkage disequilibrium the spectrum leaves out. The
lowering's `V_G` is `r` times the spectrum's genic variance, and the directional term of
selection on every mutation is `r·aλ/S`, so that the response is `r` times the genic one.

`ld = :dynamic` (keyword of the constructor and of [`spectrum_equilibrium`](@ref), through
`run_moments(…; spectrum = (; ld = :dynamic))` or `emulate(…; spectrum = (; ld = :dynamic))`;
default `:none`, the model above) adds the dynamic linkage disequilibrium described in
src/spectrum.jl: `d` the LD of the classes `ldc` ([`QTLSim.ld_classes`](@ref), `ldbins = 20`
bins of positions), the lowering's `V_G` is [`QTLSim.spectrum_vtotal`](@ref) = genic + LD,
and each site's selection uses the bins' regressions `β`. `turnover = false` drops the
`(1 + 2γ)` factor. Exclusive with `ldratio`.

`solver = :beta` (keyword of the constructor and of [`spectrum_equilibrium`](@ref), through
`run_moments(…; spectrum = (; solver = :beta))`, `emulate`, `emulated_half_survival_rate`;
default `:grid`, the frequency grid, bit for bit unchanged) replaces the grid by the beta closure
of src/spectrum.jl: per effect class the moments `P_k = ∫ p^k pq φ`, `k = 0, 1, 2` (`P`), closed by a
beta law, one explicit Euler step of the exact moment equations per generation. **Regime of
validity:** drift-shaped spectra, i.e. a fixed optimum and weak sweeps (population-scaled
directional selection `4n|aλ|/S ≲ 100` on the classes that carry the variance); not strong sweeps,
where a class's spectrum turns into a drift part plus a sweep part. **Fallback:** after every
step the classes carrying at least 0.1% of `V_G` are checked, and if one has `4n|r₁aλ|/S > 100`, a
U-shaped closed law (`α < 1`, `β < 1`) or needed the Hausdorff projection beyond 1e-6, the step is
discarded, the state is handed to the grid (`V_G` kept exactly) and the run continues on the grid;
`fallback` records the step ([`QTLSim.spectrum_fallback`](@ref)). Deterministic. Works with
`ld = :dynamic`, not with `ldratio`. Equivalence with the grid: validate/spectrum_beta_criteria.md,
validate/spectrum_beta.log.
"""
mutable struct AlleleSpectrum
    p::Vector{Float64}          # frequency grid, 0 and 1 included
    h::Vector{Float64}          # the frequency width of each node (0 at the boundaries)
    a::Vector{Float64}          # mutation effects (Gauss–Hermite nodes of N(0, α²))
    g::Vector{Float64}          # their probabilities
    n::Matrix{Float64}          # length(a) × length(p): expected numbers of sites
    L::Int
    μ::Float64
    α2::Float64
    λ::Float64                  # the lag θ − z̄ at the last selection
    S::Float64                  # ω_e² + V_P at the last selection (Inf: neutral)
    ω::Vector{Float64}          # pair classes of loci (pair_classes): weights,
    ρ::Vector{Float64}          # and split probabilities (class 1: the same locus)
    linkage::Symbol             # :barton (linked selection lowers N_e) or :none
    cd::Float64                 # the extra drift of linked selection in the last step (Var(Δp) = pq(1/(2n) + cd))
    ldratio::Float64            # prototype (validate/spectrum_diagnosis.jl): V_G/genic, the linkage disequilibrium factor (1: none)
    ld::Symbol                  # :none (linkage equilibrium, the default) or :dynamic (Bulmer LD, `d`)
    ldc::LDClasses              # ld = :dynamic: the classes of pairs of distinct loci by position bins and split probability
    d::Vector{Float64}          # ld = :dynamic: the LD of each class of `ldc`, genotype units
    β::Vector{Float64}          # ld = :dynamic: each bin's regression of the genotype on its loci at the last selection
    γ::Float64                  # ld = :dynamic: the relative change of the genic variance by the last selection
    gsel::Float64               # ld = :dynamic: the genic variance at the last selection
    turnover::Bool              # ld = :dynamic: LD also scales with its sites' variances under selection (γ)
    solver::Symbol              # the current solver: :grid (the frequency grid) or :beta (the beta closure of the moments)
    requested::Symbol           # the solver asked for (stays :beta after a fallback to the grid)
    P::Matrix{Float64}          # solver = :beta: length(a) × 3, the moments P_k = ∫ p^k ψ of ψ = pq φ, k = 0, 1, 2
    Pw::Matrix{Float64}         # solver = :beta: work array of P
    nsteps::Int                 # spectrum_step! calls since construction or the last reset (generations)
    fallback::Int               # solver = :beta: the step at which the closure fell back to the grid (−1: never)
    # work arrays
    lo::Vector{Float64}; di::Vector{Float64}; up::Vector{Float64}; rhs::Vector{Float64}; cp::Vector{Float64}
end

"Gauss–Hermite nodes and weights for the standard normal (Golub–Welsch)."
function gauss_hermite_normal(K::Integer)
    A = zeros(K, K)
    for i in 1:K-1
        A[i, i+1] = A[i+1, i] = sqrt(i)          # probabilists' Hermite recurrence
    end
    E = eigen(Symmetric(A))
    E.values, E.vectors[1, :] .^ 2
end

function AlleleSpectrum(p::Params; grid::Integer = 240, classes::Integer = 20, linkage::Symbol = :barton, ldratio::Real = 1.0,
                        ld::Symbol = :none, turnover::Bool = true, ldbins::Integer = 20, solver::Symbol = :grid)
    ld in (:none, :dynamic) || throw(ArgumentError("ld is :none or :dynamic; got :$ld"))
    solver in (:grid, :beta) || throw(ArgumentError("solver is :grid or :beta; got :$solver"))
    solver === :beta && ldratio != 1 && throw(ArgumentError("solver = :beta does not implement the prototype ldratio"))
    ld === :dynamic && ldratio != 1 && throw(ArgumentError("ld = :dynamic and the prototype ldratio are exclusive"))
    0 < ldratio <= 1 || throw(ArgumentError("ldratio must be in (0, 1], got $ldratio"))
    linkage in (:barton, :none) || throw(ArgumentError("linkage is :barton or :none; got :$linkage"))
    p.Ploidy == 2 || throw(ArgumentError("loci = :spectrum is for diploids; Ploidy = $(p.Ploidy)"))
    p.NoTraits == 1 || throw(ArgumentError("loci = :spectrum is for one trait; NoTraits = $(p.NoTraits)"))
    iseven(classes) || throw(ArgumentError("classes must be even (no neutral node)"))
    α2 = mutation_effect_variance(p)
    M = grid
    x = [(1 - cos(π * i / (M + 1))) / 2 for i in 0:M+1]
    x[1], x[end] = 0.0, 1.0
    h = zeros(M + 2)
    for i in 2:M+1
        h[i] = (x[i+1] - x[i-1]) / 2
    end
    z, w = gauss_hermite_normal(classes)
    ω, ρ = pair_classes(p)
    ldbins >= 1 || throw(ArgumentError("ldbins must be at least 1, got $ldbins"))
    ldc = ld === :dynamic ? ld_classes(p; bins = ldbins) : LDClasses(Int[], Int[], Float64[], Float64[], Float64[])
    AlleleSpectrum(x, h, sqrt(α2) .* z, w, zeros(classes, M + 2), p.NoQTL, p.MutationRatePerGamete, α2,
                   0.0, Inf, ω, ρ, linkage, 0.0, float(ldratio), ld, ldc, zeros(length(ldc.ρ)), ones(length(ldc.m)), 0.0, 0.0, turnover,
                   solver, solver, zeros(classes, 3), zeros(classes, 3), 0, -1, zeros(M + 2), zeros(M + 2), zeros(M + 2), zeros(M + 2), zeros(M + 2))
end

"The within-population genetic variance `Σ 2a² p(1 − p)` of the spectrum (`Σ 2a² P₀` with `solver = :beta`)."
function spectrum_vg(sp::AlleleSpectrum)
    sp.solver === :beta && return _beta_vg(sp)
    v = 0.0
    for i in 2:length(sp.p)-1
        pq = sp.p[i] * (1 - sp.p[i])
        for j in eachindex(sp.a)
            v += 2sp.a[j]^2 * pq * sp.n[j, i]
        end
    end
    v
end

# the expected change of a frequency per generation divided by pq, for effect a: the
# directional term aλ/S and the stabilising one, (1 − λ²/S)(a² + 2c)(1 − 2p)/(2S)
@inline spectrum_seff(a, p, λ, S, c, r = 1.0, r2 = 1.0) =
    isinf(S) ? 0.0 : r * a * λ / S - r2 * (1 - λ^2 / S) * (a^2 + 2c) * (1 - 2p) / (2S)

"The linkage disequilibrium `Σ_c d_c` of the spectrum (`ld = :dynamic`; 0 otherwise)."
spectrum_ld(sp::AlleleSpectrum) = sp.ld === :dynamic ? sum(sp.d) : 0.0

"""
    QTLSim.ld_classes(p; bins = 20) -> LDClasses

The classes of pairs of distinct loci of [`QTLSim.LDClasses`](@ref) for `p`'s loci and crossover
(`NoQTL`, `RecombPoissonMean`, `CrossoverModel`), with the loci cut into `min(bins, NoQTL)`
bins of consecutive positions. Memoised.
"""
function ld_classes(p::Params; bins::Integer = 20)
    L = p.NoQTL
    B = clamp(bins, 1, max(L, 1))
    key = (L, Float64(p.RecombPoissonMean), p.CrossoverModel, B)
    hit = lock(() -> get(_LD_CLASSES, key, nothing), _PAIR_CLASSES_LOCK)
    hit === nothing || return hit
    c = _ld_classes(L, p.RecombPoissonMean, p.CrossoverModel, B)
    lock(() -> (_LD_CLASSES[key] = c), _PAIR_CLASSES_LOCK)
    c
end

const _LD_CLASSES = Dict{Tuple{Int,Float64,Symbol,Int},LDClasses}()

function _ld_classes(L::Integer, mean::Real, model::Symbol, B::Integer)
    binof(i) = min(B, fld(i * B, L) + 1)                  # 0-based locus i
    m = zeros(B)
    for i in 0:L-1
        m[binof(i)] += 1 / L
    end
    L <= 1 && return LDClasses(Int[], Int[], Float64[], Float64[], m)
    ρp, _ = pair_splits(L, mean, model)                   # :qtsim every pair i < j in order; :all by distance
    rkey(r) = r <= 0 ? typemin(Int) : floor(Int, log(r) / log(1.02))
    acc = Dict{Tuple{Int,Int,Int},Tuple{Float64,Float64}}()
    k = 0
    for i in 0:L-2, j in i+1:L-1
        k += 1
        r = model === :qtsim ? ρp[k] : ρp[j-i]
        kk = (binof(i), binof(j), rkey(r))
        a = get(acc, kk, (0.0, 0.0))
        acc[kk] = (a[1] + 1, a[2] + r)
    end
    ks = sort(collect(keys(acc)))
    LDClasses([x[1] for x in ks], [x[2] for x in ks], [acc[x][2] / acc[x][1] for x in ks],
              [2acc[x][1] / L^2 for x in ks], m)
end

# each bin's regression of the genotype on its loci, β_b = V_b/g_b: V_b the bin's share of
# V_G, its genic share g m_b plus the LD of its pairs (all of a pair within the bin, half of
# a pair with another bin); written into β
function _ld_betas!(β::AbstractVector, ldc::LDClasses, d::AbstractVector, g::Real)
    B = length(ldc.m)
    @inbounds for b in 1:B
        β[b] = g * ldc.m[b]
    end
    @inbounds for c in eachindex(d)
        a, b = ldc.a[c], ldc.b[c]
        if a == b
            β[a] += d[c]
        else
            β[a] += d[c] / 2; β[b] += d[c] / 2
        end
    end
    @inbounds for b in 1:B
        β[b] = g > 0 && ldc.m[b] > 0 ? max(β[b], 0.0) / (g * ldc.m[b]) : 1.0
    end
    β
end

"""
    QTLSim.spectrum_vtotal(sp)

The genetic variance of the lowering: with `ld = :dynamic` the genic variance plus the
linkage disequilibrium, `V_G = Σ 2a² p(1 − p) + Σ_c d_c` (at least 0); otherwise `ldratio`
times the genic variance (the genic variance itself by default).
"""
function spectrum_vtotal(sp::AlleleSpectrum)
    g = spectrum_vg(sp)
    sp.ld === :dynamic ? max(g + sum(sp.d), 0.0) : sp.ldratio * g
end

# the factors of the directional and the stabilising terms of selection on a site: with
# ld = :dynamic, the mean over the bins of β_b and of β_b² (a site's covariance with the
# genotype is β_b times its variance; the sites of every bin share one spectrum), where
# Σ_b m_b β_b = V_G/genic; else the prototype's (ldratio, 1). Updates sp.β.
function _ld_factors(sp::AlleleSpectrum, g)
    sp.ld === :dynamic || return (sp.ldratio, 1.0)
    _ld_betas!(sp.β, sp.ldc, sp.d, g)
    r1 = 0.0; r2 = 0.0
    for b in eachindex(sp.β)
        r1 += sp.ldc.m[b] * sp.β[b]; r2 += sp.ldc.m[b] * sp.β[b]^2
    end
    (r1, r2)
end

"""
    QTLSim.spectrum_select(sp, λ, S) -> Δz̄

Record the selection of one generation (lag `λ`, `S = ω_e² + V_P`; `S = Inf` when
neutral) for the next [`QTLSim.spectrum_step!`](@ref), and return the expected change of the mean
genotype among the survivors, `Σ 2a E[Δp]`.
"""
function spectrum_select(sp::AlleleSpectrum, λ::Real, S::Real)
    sp.λ, sp.S = λ, S
    isinf(S) && return 0.0
    g = spectrum_vg(sp)
    c = g / (2sp.L)
    r1, r2 = _ld_factors(sp, g)
    dz = 0.0
    dg = 0.0
    sp.solver === :beta && ((dz, dg) = _beta_select(sp, λ, S, c, r1, r2))
    for i in (sp.solver === :beta ? (2:1) : (2:length(sp.p)-1))
        pi = sp.p[i]; pq = pi * (1 - pi)
        for j in eachindex(sp.a)
            e = spectrum_seff(sp.a[j], pi, λ, S, c, r1, r2)
            dz += 2sp.a[j] * pq * e * sp.n[j, i]          # (this association keeps the default bit for bit)
            dg += 2sp.a[j]^2 * pq * (1 - 2pi) * e * sp.n[j, i]
        end
    end
    if sp.ld === :dynamic
        sp.γ = g > 0 ? dg / g : 0.0
        sp.gsel = g
    end
    dz
end

@inline _bern(x) = abs(x) < 1e-8 ? 1 - x / 2 : x / expm1(x)       # x/(eˣ − 1)

"""
    QTLSim.spectrum_linked(sp, n) -> c

The extra coalescence rate per generation of a site whose loosely linked loci carry
additive variance in relative fitness (Robertson 1961; Barton 1995, Genetics 140:821:
`N/N_e = 1 + var(W)/(2r²)` for a locus at recombination fraction `r` from one with
variance `var(W)`): under directional selection of gradient `b = λ/S` each other locus
carries `var(W) = b² V_G/L`, so a random site has `N/N_e = 1 + b² V_G Σ_c ω_c/(2(ρ_c + δ)²)`
over the pair classes of distinct loci ([`QTLSim.pair_classes`](@ref)), with `δ = 1/(2n) + V_G/S`
the rate at which drift and selection (the Gaussian reduction of the selected variance,
Bulmer 1971) dissolve the association of a site with its background, the two factors
besides recombination of Santiago and Caballero (1998); it keeps the sum finite for
tightly linked loci and small when selection is strong. The
locus itself is not counted (Barton's approximation is for loose linkage). Returns
`c = (N/N_e − 1)/(2n)`, so that `Var(Δp) = pq(1/(2n) + c)`. See src/spectrum.jl.
"""
function spectrum_linked(sp::AlleleSpectrum, n::Real)
    (isinf(sp.S) || sp.λ == 0) && return 0.0
    V = spectrum_vg(sp)
    b2V = (sp.λ / sp.S)^2 * V
    δ = 1 / (2n) + V / sp.S
    acc = 0.0
    for k in 2:length(sp.ω)
        acc += sp.ω[k] / (2 * (sp.ρ[k] + δ)^2)
    end
    b2V * acc / (2n)
end

"""
    QTLSim.spectrum_step!(sp, npar)

One generation of selection (as last recorded by [`QTLSim.spectrum_select`](@ref)) and drift
among `npar` parents: one implicit Euler step of the Wright–Fisher diffusion of every
effect class. See src/spectrum.jl.
"""
function spectrum_step!(sp::AlleleSpectrum, npar::Real)
    sp.nsteps += 1
    if sp.solver === :beta
        _beta_step!(sp, npar) && return sp
        _beta_to_grid!(sp)                       # out of the closure's domain: this step and the rest on the grid
    end
    x, h = sp.p, sp.h
    M = length(x) - 2
    nn = max(npar, 1.0)
    gnc = spectrum_vg(sp)
    c = isinf(sp.S) ? 0.0 : gnc / (2sp.L)
    r1, r2 = _ld_factors(sp, gnc)
    sp.cd = sp.linkage === :barton ? spectrum_linked(sp, nn) : 0.0
    d4 = 1 / (4nn) + sp.cd / 2              # D / pq: Var(Δp) = pq (1/(2n) + c_d)
    lo, di, up, rhs, cp = sp.lo, sp.di, sp.up, sp.rhs, sp.cp
    for j in eachindex(sp.a)
        a = sp.a[j]
        # unknowns n_i (i = 2:M+1); u_i = D_i n_i/h_i, D = pq/(4n); flux F_{i+½} = (B(−cΔ)u_i − B(cΔ)u_{i+1})/Δ
        fill!(lo, 0.0); fill!(di, 0.0); fill!(up, 0.0)
        for k in 1:M+1                          # interface between nodes k and k+1
            Δ = x[k+1] - x[k]
            pm = (x[k] + x[k+1]) / 2
            q = spectrum_seff(a, pm, sp.λ, sp.S, c, r1, r2) * Δ / d4
            Bm, Bp = _bern(-q) / Δ, _bern(q) / Δ
            # contribution of u_k (left) and u_{k+1} (right) to the flux F = Bm u_k − Bp u_{k+1}
            # node k loses F, node k+1 gains F (dn/dt)
            rk = k >= 2 ? x[k] * (1 - x[k]) * d4 / h[k] : 0.0
            rk1 = k + 1 <= M + 1 ? x[k+1] * (1 - x[k+1]) * d4 / h[k+1] : 0.0
            if k >= 2                           # node k: dn_k/dt −= Bm r_k n_k − Bp r_{k+1} n_{k+1}
                di[k] += Bm * rk
                k + 1 <= M + 1 && (up[k] -= Bp * rk1)
            end
            if k + 1 <= M + 1                   # node k+1: dn/dt += Bm r_k n_k − Bp r_{k+1} n_{k+1}
                di[k+1] += Bp * rk1
                k >= 2 && (lo[k+1] -= Bm * rk)
            end
        end
        # (I + A) n_new = n_old, A = [lo di up] (dt = 1 generation)
        for i in 2:M+1
            rhs[i] = sp.n[j, i]
            di[i] += 1.0
        end
        # Thomas algorithm on i = 2:M+1
        cp[2] = up[2] / di[2]; rhs[2] /= di[2]
        for i in 3:M+1
            m = di[i] - lo[i] * cp[i-1]
            cp[i] = up[i] / m
            rhs[i] = (rhs[i] - lo[i] * rhs[i-1]) / m
        end
        for i in M:-1:2
            rhs[i] -= cp[i] * rhs[i+1]
        end
        for i in 2:M+1
            sp.n[j, i] = max(rhs[i], 0.0)
        end
    end
    sp.ld === :dynamic && spectrum_ld_step!(sp.d, sp.ldc, sp.β, isinf(sp.S) ? 0.0 : sp.gsel, sp.S, nn,
                                            sp.turnover && !isinf(sp.S) ? sp.γ : 0.0)
    sp
end

"""
    QTLSim.spectrum_ld_step!(d, ldc, β, g, S, n, γ) -> d

One generation of the linkage disequilibrium of the classes `ldc` of pairs of distinct loci
([`QTLSim.LDClasses`](@ref); `ld = :dynamic`, src/spectrum.jl): for every class `c` of pairs
between bins `a` and `b` (weight `ω_c`, split probability `ρ_c`),

    d_c ← (1 − 1/(2n)) [(1 − ρ_c)(1 + 2γ) d_c − ω_c β_a β_b g²/(2S)],

Bulmer's generation of negative LD by Gaussian selection of `S = ω_e² + V_P` (a pair of
loci gains −V_i V_j/S, `V_i = β_b g/L` a locus's covariance with the genotype, the half of
selection's change that a gamete inherits; `g` the genic variance, `g = 0` when neutral),
recombination at the class's rate, drift among `n` parents (`n = Inf`: none) and, through
`γ`, the change of the LD of a pair of sites with the variances of its sites under
selection (`Δ ln D = Δ ln pq_i + Δ ln pq_j`; `γ` the relative change of the genic variance
by selection in the generation). With every `β_b = β` the generation is `ω_c V²/(2S)`,
`V = βg`, the multilocus model's (src/multilocus.jl).
"""
function spectrum_ld_step!(d::AbstractVector, ldc::LDClasses, β::AbstractVector, g::Real, S::Real, n::Real, γ::Real)
    dr = isinf(n) ? 1.0 : 1 - 1 / (2n)
    G = max(1 + 2γ, 0.0)
    gen = isinf(S) || g == 0 ? 0.0 : g^2 / (2S)
    @inbounds for c in eachindex(d)
        d[c] = dr * ((1 - ldc.ρ[c]) * G * d[c] - ldc.ω[c] * β[ldc.a[c]] * β[ldc.b[c]] * gen)
    end
    d
end

"""
    QTLSim.spectrum_mutate!(sp, N)

`2N L μ` new mutations among the gametes of `N` offspring, each at frequency `1/(2N)`,
spread over the effect classes by their probabilities and over the two grid nodes around
`1/(2N)` so that their number and mean frequency are exact.
"""
function spectrum_mutate!(sp::AlleleSpectrum, N::Real)
    N < 1 && return sp
    sp.solver === :beta && return _beta_mutate!(sp, N)
    p0 = 1 / (2N)
    x = sp.p
    k = searchsortedlast(x, p0)
    k = clamp(k, 2, length(x) - 2)              # p0 between x[k] and x[k+1], interior nodes when possible
    total = 2N * sp.L * sp.μ
    if x[k] > p0                                # below the first interior node: put it there
        wl, wr = 1.0, 0.0
    else
        wr = (p0 - x[k]) / (x[k+1] - x[k]); wl = 1 - wr
    end
    for j in eachindex(sp.a)
        sp.n[j, k] += total * sp.g[j] * wl
        sp.n[j, k+1] += total * sp.g[j] * wr
    end
    sp
end

# --- the beta closure (solver = :beta; off by default) ------------------------------------
#
# A faster stand-in for the frequency grid when sweeps are weak (validate/spectrum_beta_criteria.md,
# validate/spectrum_beta.jl). Per effect class a_j, the pq-weighted spectrum ψ = pq φ is a finite
# measure on [0, 1] without atoms at 0 and 1 (the neutral one is ψ ∝ 1 − p), represented by its
# moments P_k = ∫ p^k ψ, k = 0, 1, 2; V_G = Σ 2a² P₀ and Δz̄ = Σ 2a ∫ s ψ. For any f with
# f(0) = f(1) = 0, the diffusion of φ (drift M = pq s(p), D = pq d₄, s(p) = s₀ + s₁p the
# selection of spectrum_seff, linear in p, d₄ = 1/(4n) + c_d/2 with Barton's c_d) gives
#
#   d/dt ∫ f φ = ∫ [s f′ + d₄ f″] ψ,
#
# and with f_m = p^{m+1}(1 − p), ∫ f_m φ = P_m:
#
#   dP_m = s₀[(m+1)P_m − (m+2)P_{m+1}] + s₁[(m+1)P_{m+1} − (m+2)P_{m+2}]
#          + d₄[m(m+1)P_{m−1} − (m+1)(m+2)P_m],
#
# exact (absorption at 0 and 1 included: f vanishes there), but dP₂ needs P₃ and P₄. The closure
# takes ψ/P₀ = Beta(α, β) with the same mean μ₁ and second moment μ₂: ν = μ₁(1 − μ₁)/(μ₂ − μ₁²) − 1,
# α = μ₁ν, β = (1 − μ₁)ν, μ₃ = μ₂(α + 2)/(ν + 2), μ₄ = μ₃(α + 3)/(ν + 3). It is exact on the neutral
# spectrum, Beta(1, 2). One generation is one explicit Euler step of these equations (the grid
# takes one implicit step: both are first order in the generation; the grid's implicit step adds
# an O(Δt²) term to the λ² response of V_G, 2(A_{a⁴} − 2Q₄)/S², which is not "fixed" here or
# there: this solver is checked against the grid as it is), followed by the projection of
# (μ₁, μ₂) onto the Hausdorff moment conditions of a measure on [0, 1], 0 < μ₁ < 1,
# μ₁² < μ₂ < μ₁. Mutations add 2N L μ g_j f_m(1/(2N)) to P_m.
#
# Domain and fallback. The closure is faithful while drift shapes each class's ψ, as at a fixed
# optimum (stabilising selection and drift) and under weak sweeps; when directional selection on
# a class is strong against drift, its sites sweep deterministically, ψ becomes a mixture of a
# drift part near 0 and a sweep part (U-shaped in the limit), and two moments no longer carry its
# dynamics. After every step, each class carrying at least BETA_SHARE of V_G (and not an atom of
# new mutations, α + β < BETA_NUMAX) is checked; the step is outside the domain if for such a
# class (i) its population-scaled directional selection Γ = 4n|r₁ a λ|/S exceeds BETA_GMAX = 100,
# (ii) the closed law is U-shaped, α < 1 and β < 1 (this includes μ₂ → μ₁, ν → 0), or (iii) the
# Hausdorff projection (or P₀ < 0) moved (μ₁, μ₂) by more than BETA_HTOL of μ₁(1 − μ₁).
# Criterion (i) and its threshold were set on four published inputs before the validation
# (validate/spectrum_beta_criteria.md): the shape criteria (ii) and (iii) alone fire too late
# (04943, 4Nα²/S ≈ 12: V_G 38% above the grid's 100 generations into the move, U-shape at 126)
# or never (06159: 18% below at +600), and a tracked third moment compared with the beta law's
# saturated in accurate and inaccurate runs alike. Then that step is discarded, the last state inside the domain is converted to the grid (each
# class's beta law integrated over the grid's cells, so V_G is kept exactly) and the step and the
# rest of the run use the grid; the run never switches back. `fallback` records the step. The
# rule depends only on the state, so it is deterministic. Per run rather than per class, because
# the classes are coupled through V_G, λ, n and Barton's term: a mixed state would need two
# solvers in one generation for no gain (once a class leaves, the run is in the strong-sweep
# regime the grid is for).

"Classes carrying less than this share of `V_G` are not checked by the beta closure's domain test."
const BETA_SHARE = 1e-3
"The relative move of `(μ₁, μ₂)` by the Hausdorff projection beyond which the beta closure falls back to the grid."
const BETA_HTOL = 1e-6
"The population-scaled directional selection `4n|a λ|/S` of a class beyond which the beta closure falls back to the grid."
const BETA_GMAX = 100.0
"Classes whose beta law is nearly an atom (`α + β` above this; new mutations) are not checked."
const BETA_NUMAX = 1e6

"""
    QTLSim.beta_closure(μ1, μ2) -> (μ3, μ4, α, β)

The third and fourth raw moments of the beta law with mean `μ1` and raw second moment `μ2`
(`0 < μ1 < 1`, `μ1² < μ2 < μ1`), and its parameters: `ν = μ1(1 − μ1)/(μ2 − μ1²) − 1`,
`α = μ1ν`, `β = (1 − μ1)ν`, `μ3 = μ2(α + 2)/(ν + 2)`, `μ4 = μ3(α + 3)/(ν + 3)`. The closure of
`solver = :beta` (src/spectrum.jl); exact for the neutral spectrum `ψ ∝ 1 − p` (`α = 1`, `β = 2`).
"""
@inline function beta_closure(μ1::Real, μ2::Real)
    ν = μ1 * (1 - μ1) / (μ2 - μ1^2) - 1
    α = μ1 * ν
    μ3 = μ2 * (α + 2) / (ν + 2)
    (μ3, μ3 * (α + 3) / (ν + 3), α, (1 - μ1) * ν)
end

# project (μ₁, μ₂) of row j of Q onto 0 < μ₁ < 1, μ₁² < μ₂ < μ₁; returns the relative move
# (Inf when P₀ < 0) and whether the closed beta law is U-shaped (α < 1 and β < 1)
function _beta_project!(Q::AbstractMatrix, j::Integer)
    P0 = Q[j, 1]
    if P0 <= 0
        Q[j, 1] = Q[j, 2] = Q[j, 3] = 0.0
        return (P0 < 0 ? Inf : 0.0), false
    end
    m1, m2 = Q[j, 2] / P0, Q[j, 3] / P0
    μ1 = clamp(m1, 1e-12, 1 - 1e-12)
    μ2 = clamp(m2, μ1^2 * (1 + 1e-9), μ1 * (1 - 1e-9))
    Q[j, 2], Q[j, 3] = μ1 * P0, μ2 * P0
    _, _, α, β = beta_closure(μ1, μ2)
    max(abs(μ1 - m1), abs(μ2 - m2)) / (μ1 * (1 - μ1)), α < 1 && β < 1
end

_beta_vg(sp::AlleleSpectrum) = (v = 0.0; for j in eachindex(sp.a); v += 2sp.a[j]^2 * sp.P[j, 1]; end; v)

# Δz̄ = Σ 2a ∫ s ψ and Σ 2a² ∫ (1 − 2p) s ψ with s = e₀ + e₁p (spectrum_seff)
function _beta_select(sp::AlleleSpectrum, λ, S, c, r1, r2)
    dz = 0.0; dg = 0.0
    for j in eachindex(sp.a)
        a = sp.a[j]
        P0, P1, P2 = sp.P[j, 1], sp.P[j, 2], sp.P[j, 3]
        κ = r2 * (1 - λ^2 / S) * (a^2 + 2c) / (2S)
        e0, e1 = r1 * a * λ / S - κ, 2κ
        dz += 2a * (e0 * P0 + e1 * P1)
        dg += 2a^2 * (e0 * (P0 - 2P1) + e1 * (P1 - 2P2))
    end
    dz, dg
end

"""
    QTLSim.beta_dP(m, s0, s1, d4, Pm1, Pm, Pp1, Pp2)

The change in one generation of `P_m = ∫ p^m ψ`, `ψ = pq φ`, under selection `s(p) = s0 + s1 p`
and drift `d4` (`Var(Δp) = 2 d4 pq`), from `P_{m−1}, P_m, P_{m+1}, P_{m+2}`:
`s0[(m+1)P_m − (m+2)P_{m+1}] + s1[(m+1)P_{m+1} − (m+2)P_{m+2}] + d4[m(m+1)P_{m−1} − (m+1)(m+2)P_m]`
(src/spectrum.jl; checked against the generator in validate/symbolic/check_spectrum.jl).
"""
@inline beta_dP(m, s0, s1, d4, Pm1, Pm, Pp1, Pp2) =
    s0 * ((m + 1) * Pm - (m + 2) * Pp1) + s1 * ((m + 1) * Pp1 - (m + 2) * Pp2) + d4 * (m * (m + 1) * Pm1 - (m + 1) * (m + 2) * Pm)

# one generation of selection and drift on the moments (explicit Euler); false, leaving the
# state as it was, if the step leaves the closure's domain
function _beta_step!(sp::AlleleSpectrum, npar::Real)
    nn = max(npar, 1.0)
    g = _beta_vg(sp)
    neutral = isinf(sp.S)
    c = neutral ? 0.0 : g / (2sp.L)
    r1, r2 = _ld_factors(sp, g)
    cd = sp.linkage === :barton ? spectrum_linked(sp, nn) : 0.0
    d4 = 1 / (4nn) + cd / 2
    P, Q = sp.P, sp.Pw
    ok = true
    for j in eachindex(sp.a)
        a = sp.a[j]
        P0 = P[j, 1]
        if P0 <= 0
            Q[j, 1] = Q[j, 2] = Q[j, 3] = 0.0
            continue
        end
        μ1 = clamp(P[j, 2] / P0, 1e-12, 1 - 1e-12)              # (an atom, as the first mutations, is a limit of beta laws)
        μ2 = clamp(P[j, 3] / P0, μ1^2 * (1 + 1e-9), μ1 * (1 - 1e-9))
        μ3, μ4, α, β = beta_closure(μ1, μ2)
        P1, P2, P3, P4 = μ1 * P0, μ2 * P0, μ3 * P0, μ4 * P0
        if neutral
            s0 = s1 = 0.0
        else
            κ = r2 * (1 - sp.λ^2 / sp.S) * (a^2 + 2c) / (2sp.S)
            s0, s1 = r1 * a * sp.λ / sp.S - κ, 2κ
        end
        Q[j, 1] = P0 + beta_dP(0, s0, s1, d4, 0.0, P0, P1, P2)
        Q[j, 2] = P1 + beta_dP(1, s0, s1, d4, P0, P1, P2, P3)
        Q[j, 3] = P2 + beta_dP(2, s0, s1, d4, P1, P2, P3, P4)
        bind, ushape = _beta_project!(Q, j)
        if g > 0 && 2a^2 * P0 >= BETA_SHARE * g && α + β < BETA_NUMAX
            Γ = neutral ? 0.0 : 4nn * abs(r1 * a * sp.λ) / sp.S
            (Γ > BETA_GMAX || ushape || bind > BETA_HTOL) && (ok = false)
        end
    end
    ok || return false
    sp.cd = cd
    P .= Q
    sp.ld === :dynamic && spectrum_ld_step!(sp.d, sp.ldc, sp.β, neutral ? 0.0 : sp.gsel, sp.S, nn,
                                            sp.turnover && !neutral ? sp.γ : 0.0)
    true
end

function _beta_mutate!(sp::AlleleSpectrum, N::Real)
    p0 = 1 / (2N)
    total = 2N * sp.L * sp.μ
    for j in eachindex(sp.a)
        f = total * sp.g[j] * p0 * (1 - p0)
        sp.P[j, 1] += f; sp.P[j, 2] += f * p0; sp.P[j, 3] += f * p0^2
    end
    sp
end

# hand the state to the grid: each class's beta law ψ/P₀ integrated over the cells of the grid's
# interior nodes (edges halfway between nodes, the first at 0 and the last at 1), n_i = ψ(cell)/(p_i q_i),
# so that the grid's V_G is the closure's exactly
function _beta_to_grid!(sp::AlleleSpectrum)
    x = sp.p
    M = length(x) - 2
    fill!(sp.n, 0.0)
    for j in eachindex(sp.a)
        P0 = sp.P[j, 1]
        P0 > 0 || continue
        μ1 = clamp(sp.P[j, 2] / P0, 1e-12, 1 - 1e-12)
        _, _, α, β = beta_closure(μ1, clamp(sp.P[j, 3] / P0, μ1^2 * (1 + 1e-9), μ1 * (1 - 1e-9)))
        D = Beta(α, β)
        Fl = 0.0
        for i in 2:M+1
            Fu = i == M + 1 ? 1.0 : cdf(D, (x[i] + x[i+1]) / 2)
            sp.n[j, i] = P0 * (Fu - Fl) / (x[i] * (1 - x[i]))
            Fl = Fu
        end
    end
    sp.solver = :grid
    sp.fallback = sp.nsteps
    sp
end

"""
    QTLSim.spectrum_fallback(sp) -> Int

For a spectrum made with `solver = :beta`, the step (generation, counted from the start of the
run) at which the beta closure left its domain and the run fell back to the frequency grid:
`-1` if it never did, `0` if it fell back before the run (in the equilibrium that started it).
`-1` for the grid solver. See src/spectrum.jl.
"""
spectrum_fallback(sp::AlleleSpectrum) = sp.fallback

"Scale the spectrum so that its genetic variance is `V` (a spectrum without sites stays empty); with `ld = :dynamic` its genic variance, the linkage disequilibrium scaled in proportion."
function spectrum_scale!(sp::AlleleSpectrum, V::Real)
    if sp.ld === :dynamic                       # the genic variance to V, the LD in proportion
        v = spectrum_vg(sp)
        v > 0 && (sp.n .*= V / v; sp.d .*= V / v; sp.P .*= V / v)
        return sp
    end
    v = spectrum_vg(sp)
    v > 0 && (sp.n .*= V / v; sp.P .*= V / v)
    sp
end

"""
    spectrum_equilibrium(p; tol = 1e-9, maxgen = 40 * p.PopSize + 20_000, grid = 240, classes = 20)

The stationary allele-frequency spectrum of the moments lowering with `loci = :spectrum`
at a fixed optimum (src/spectrum.jl): the recursion of one generation (Gaussian
selection with the lag 0, drift among the expected survivors, mutation; the Beverton–Holt
equilibrium size under density regulation) iterated from no variation until `V_G`
changes by less than `tol` relative to itself over 100 generations. Returns `(; VG, N,
nparents, spectrum)`. For mutations in linkage equilibrium and weak selection its `VG` is
close to [`house_of_gauss_vg`](@ref). With `ld = :dynamic` (and `ldbins`, `turnover`; see
[`QTLSim.AlleleSpectrum`](@ref)) `VG` is genic + LD, and the result also has `genic` and `ld`
(expected sample variances, ×N/(N − 1)). With `solver = :beta` the equilibrium is the beta
closure's (milliseconds instead of seconds), the grid's if the closure falls back (the
spectrum's `fallback ≥ 0`); see [`QTLSim.AlleleSpectrum`](@ref).
"""
function spectrum_equilibrium(p::Params; tol::Real = 1e-9, maxgen::Integer = 40 * p.PopSize + 20_000,
                              grid::Integer = 240, classes::Integer = 20, linkage::Symbol = :barton, ldratio::Real = 1.0,
                              ld::Symbol = :none, turnover::Bool = true, ldbins::Integer = 20, solver::Symbol = :grid)
    sp = AlleleSpectrum(p; grid, classes, linkage, ldratio, ld, turnover, ldbins, solver)
    dyn = ld === :dynamic
    VE = envVar(p)
    Ω = initFitSD(p)^2
    c0 = QTSIM_FITNESS_PEAK
    N = float(p.PopSize)
    VG, npar = 0.0, N
    last = -1.0
    for t in 1:maxgen
        P = (dyn ? spectrum_vtotal(sp) : sp.ldratio * VG) + VE
        ωe2 = p.StandardisePhenotype ? Ω * P : Ω
        S = p.Drift ? Inf : ωe2 + P
        Wbar = p.Drift ? 1.0 : c0 * sqrt(ωe2 / S)
        W2 = p.Drift ? 1.0 : c0^2 * sqrt(ωe2 / (ωe2 + 2P))
        spectrum_select(sp, 0.0, S)
        npar = p.Selection === :fecundity ? N * Wbar^2 / W2 : N * Wbar
        spectrum_step!(sp, npar)
        nw = N * Wbar
        N = p.Demography === :fixed ? float(p.PopSize) : min(p.PopSize, p.GrowthRate * nw / (1 + (p.GrowthRate - 1) * nw / p.PopSize))
        N < 1 && throw(ArgumentError("the population cannot persist at the optimum"))
        spectrum_mutate!(sp, N)
        VG = dyn ? spectrum_vtotal(sp) : spectrum_vg(sp)
        if t % 100 == 0
            last > 0 && abs(VG - last) <= tol * VG * 100 && break
            last = VG
        end
    end
    dyn && return (; VG = VG * N / (N - 1), genic = spectrum_vg(sp) * N / (N - 1), ld = spectrum_ld(sp) * N / (N - 1),
                   N, nparents = npar, spectrum = sp)
    (; VG = VG * N / (N - 1), N, nparents = npar, spectrum = sp)
end

"""
    spectrum_vg(p; kw...)

The equilibrium genetic variance (expected sample variance) of the allele-frequency
spectrum model, `spectrum_equilibrium(p; kw...).VG`; see [`spectrum_equilibrium`](@ref).
"""
spectrum_vg(p::Params; kw...) = spectrum_equilibrium(p; kw...).VG
