# Analytical predictions the simulations are compared with: the equilibrium genetic
# variance under mutation–selection balance (Gaussian and house-of-cards approximations),
# the lag behind an optimum moving at a constant rate, Bürger and Lynch's critical rate of
# environmental change, and the paper's adjustment for the mutation scalar. The forms are
# those the paper's results table uses (see vignettes 02, 03 and 09); the algebra is in
# lean/CategoricalQuantGen/Theory.lean.

"""
    selection_variance(p) -> V_s

The width of stabilising selection on the genotype, `V_s = ω² + V_E`, for the initial
fitness width ω of `p` ([`initFitSD`](@ref)). Selection on a phenotype `z = G + E` with
`E ~ N(0, V_E)` acts on `G` as Gaussian selection of variance `ω² + V_E`.
"""
selection_variance(p::Params) = initFitSD(p)^2 + envVar(p)

"""
    gaussian_vg(p)

The Gaussian approximation to the equilibrium genetic variance under mutation–selection
balance, `√(V_m V_s)`, with `V_m = m·V_E` ([`muteVar`](@ref)) and `V_s` from
[`selection_variance`](@ref) (Lande 1975). It ignores drift, and does not depend on the
number of loci. The moments lowering (`run_moments`) has the slightly larger fixed point
[`gaussian_vg_exact`](@ref).
"""
gaussian_vg(p::Params) = sqrt(muteVar(p) * selection_variance(p))

"""
    gaussian_vg_exact(p)

The fixed point of the Gaussian recursion `V ← V − V²/(V_s + V) + V_m`, which the
moments lowering iterates: `(V_m + √(V_m² + 4 V_m V_s))/2`. It is
`√(V_m V_s)·(1 + O(√(V_m/V_s)))`, close to [`gaussian_vg`](@ref) when mutation is weak
relative to selection (Lean `gaussian_fixed_point`).
"""
function gaussian_vg_exact(p::Params)
    Vm, Vs = muteVar(p), selection_variance(p)
    (Vm + sqrt(Vm^2 + 4Vm * Vs)) / 2
end

"""
    effective_size(N, nsurv)

Falconer's effective size `N_e = 4N/(V_k + 2)` when `nsurv` of the `N` individuals
survive and each of the next generation's `2N` gametes comes from a survivor chosen at
random: the number a survivor contributes is binomial, with variance
`V_k = 2N (nsurv − 1)/nsurv²`. With every individual surviving, `N_e ≈ N`.
"""
function effective_size(N::Real, nsurv::Real)
    Vk = 2N * (nsurv - 1) / nsurv^2
    4N / (Vk + 2)
end

"""
    house_of_cards_vg(p; nsurv = p.PopSize)

Bürger's stochastic house-of-cards approximation to the equilibrium genetic variance,
`4nμ V_s / (1 + V_s/(N_e α²))`, with `n` loci, mutation rate `μ` per locus per gamete,
`α² = V_m/(2nμ)` the variance of one mutation's effect, and `N_e` from
[`effective_size`](@ref) with `nsurv` survivors per generation (Bürger 1989; the paper's
Supplementary Figure S1). Unlike [`gaussian_vg`](@ref), it grows with the number of loci.
"""
function house_of_cards_vg(p::Params; nsurv::Real = p.PopSize)
    n, μ = p.NoQTL, p.MutationRatePerGamete
    Vs = selection_variance(p)
    α2 = muteVar(p) / (2n * μ)
    Ne = effective_size(p.PopSize, nsurv)
    4n * μ * Vs / (1 + Vs / (Ne * α2))
end

"""
    steady_lag(k; VG, Vs)
    steady_lag(p; VG)

The lag `θ − z̄` at which a population whose mean responds to Gaussian selection keeps
pace with an optimum moving `k` per generation: `k (V_G + V_s)/V_G` (Lande; Bürger and
Lynch 1995). The response to a lag `λ` is `V_G λ/(V_G + V_s)` per generation, which equals
`k` at that lag (Lean `steady_lag_response`). With `p`, `k = MoveFitnessSD` and
`V_s` = [`selection_variance`](@ref)`(p)`, for `MoveType = :continuous`.
"""
steady_lag(k::Real; VG::Real, Vs::Real) = k * (VG + Vs) / VG
steady_lag(p::Params; VG::Real) = steady_lag(p.MoveFitnessSD; VG, Vs = selection_variance(p))

"""
    critical_rate(; n, VG, ω, VE, sd50, N)

Bürger and Lynch's critical rate of environmental change, as the paper's results table
computes it, in units of the phenotypic standard deviation `sd50` at 50% mortality:

    k_c = V_G/(V_G + V_s) · √(2 V_λ ln B₀),   V_s = ω² + V_E,  V_λ = V_s + V_G + V_s/(2N_e),
    B₀ = N₀ ω/√V_λ,   N₀ = N(2n − 3)/(2n(n − 2)),   V_k = 2N₀(1 − 2N₀/N),  N_e = 4N/(V_k + 2),

where `n` is the number of survivors at the table's generation 10,047 and `V_G` the
genetic variance just before the move. `NaN` where the table has none: `n < 3`, or
`B₀ ≤ 1`. (The paper's equation 2 is printed without the factor 2 and with `V_s + V_G` for
`V_λ`; this form reproduces its table to machine precision, see vignette 09.)
"""
function critical_rate(; n::Real, VG::Real, ω::Real, VE::Real, sd50::Real, N::Real)
    n >= 3 || return NaN
    No = N * (2n - 3) / (2n * (n - 2))          # offspring per surviving parent
    Vk = 2No * (1 - 2No / N)                    # variance in offspring number
    Ne = 4N / (Vk + 2)
    Vs = ω^2 + VE
    Vλ = Vs + VG + Vs / (2Ne)
    B0 = No * ω / sqrt(Vλ)
    B0 > 1 ? VG / (VG + Vs) * sqrt(2Vλ * log(B0)) / sd50 : NaN
end

"""
    mutation_adjusted(rate, m, dist; w = ones(length(m)))

The paper's equation 1: each rate in `rate` weighted by the density of its population's
mutation scalar `m` under `dist` (the paper's gamma distribution: the maximum-likelihood
fit to all 119,984 mutation scalars, shape 1.201 and scale 0.004635), relative to the
(`w`-weighted) mean density.
"""
function mutation_adjusted(rate::AbstractVector, m::AbstractVector, dist; w = ones(length(m)))
    d = pdf.(Ref(dist), m)
    rate .* d ./ (sum(w .* d) / sum(w))
end
