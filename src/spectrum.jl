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
# Simplifications: linkage equilibrium apart from the linked-selection term (no Bulmer
# LD; Barton's loose-linkage N_e for every pair class, and none for the mutations of the
# same locus); infinite sites (two
# mutations of one locus are separate sites, not one allele); weak selection per
# mutation (Wright's formula); the expected spectrum along the deterministic trajectory
# of λ and n (the nonlinear feedback through λ and n is taken at the expectation, as in
# the rest of the moments lowering).

"""
    QTLSim.AlleleSpectrum

The expected allele-frequency spectrum of the segregating mutations of a diploid
population (src/spectrum.jl): `n[j, i]` the expected number of mutations of effect `a[j]`
at frequency `p[i]` (nodes `p[1] = 0` and `p[end] = 1` are the absorbing boundaries and
stay empty). `b`, `κ` and `S` are the selection the next `spectrum_step!` applies.
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

function AlleleSpectrum(p::Params; grid::Integer = 240, classes::Integer = 20, linkage::Symbol = :barton)
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
    AlleleSpectrum(x, h, sqrt(α2) .* z, w, zeros(classes, M + 2), p.NoQTL, p.MutationRatePerGamete, α2,
                   0.0, Inf, ω, ρ, linkage, 0.0, zeros(M + 2), zeros(M + 2), zeros(M + 2), zeros(M + 2), zeros(M + 2))
end

"The within-population genetic variance `Σ 2a² p(1 − p)` of the spectrum."
function spectrum_vg(sp::AlleleSpectrum)
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
@inline spectrum_seff(a, p, λ, S, c) = isinf(S) ? 0.0 : a * λ / S - (1 - λ^2 / S) * (a^2 + 2c) * (1 - 2p) / (2S)

"""
    QTLSim.spectrum_select(sp, λ, S) -> Δz̄

Record the selection of one generation (lag `λ`, `S = ω_e² + V_P`; `S = Inf` when
neutral) for the next [`QTLSim.spectrum_step!`](@ref), and return the expected change of the mean
genotype among the survivors, `Σ 2a E[Δp]`.
"""
function spectrum_select(sp::AlleleSpectrum, λ::Real, S::Real)
    sp.λ, sp.S = λ, S
    isinf(S) && return 0.0
    c = spectrum_vg(sp) / (2sp.L)
    dz = 0.0
    for i in 2:length(sp.p)-1
        pi = sp.p[i]; pq = pi * (1 - pi)
        for j in eachindex(sp.a)
            dz += 2sp.a[j] * pq * spectrum_seff(sp.a[j], pi, λ, S, c) * sp.n[j, i]
        end
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
    x, h = sp.p, sp.h
    M = length(x) - 2
    nn = max(npar, 1.0)
    c = isinf(sp.S) ? 0.0 : spectrum_vg(sp) / (2sp.L)
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
            q = spectrum_seff(a, pm, sp.λ, sp.S, c) * Δ / d4
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
    sp
end

"""
    QTLSim.spectrum_mutate!(sp, N)

`2N L μ` new mutations among the gametes of `N` offspring, each at frequency `1/(2N)`,
spread over the effect classes by their probabilities and over the two grid nodes around
`1/(2N)` so that their number and mean frequency are exact.
"""
function spectrum_mutate!(sp::AlleleSpectrum, N::Real)
    N < 1 && return sp
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

"Scale the spectrum so that its genetic variance is `V` (a spectrum without sites stays empty)."
function spectrum_scale!(sp::AlleleSpectrum, V::Real)
    v = spectrum_vg(sp)
    v > 0 && (sp.n .*= V / v)
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
close to [`house_of_gauss_vg`](@ref).
"""
function spectrum_equilibrium(p::Params; tol::Real = 1e-9, maxgen::Integer = 40 * p.PopSize + 20_000,
                              grid::Integer = 240, classes::Integer = 20, linkage::Symbol = :barton)
    sp = AlleleSpectrum(p; grid, classes, linkage)
    VE = envVar(p)
    Ω = initFitSD(p)^2
    c0 = QTSIM_FITNESS_PEAK
    N = float(p.PopSize)
    VG, npar = 0.0, N
    last = -1.0
    for t in 1:maxgen
        P = VG + VE
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
        VG = spectrum_vg(sp)
        if t % 100 == 0
            last > 0 && abs(VG - last) <= tol * VG * 100 && break
            last = VG
        end
    end
    (; VG = VG * N / (N - 1), N, nparents = npar, spectrum = sp)
end

"""
    spectrum_vg(p; kw...)

The equilibrium genetic variance (expected sample variance) of the allele-frequency
spectrum model, `spectrum_equilibrium(p; kw...).VG`; see [`spectrum_equilibrium`](@ref).
"""
spectrum_vg(p::Params; kw...) = spectrum_equilibrium(p; kw...).VG
