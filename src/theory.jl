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

# --- fluctuating optima ----------------------------------------------------
#
# Under Gaussian selection with a constant genetic variance V_G, the mean responds to the
# lag λ_t = θ_t − z̄_t as z̄_{t+1} = z̄_t + β λ_t, with β = V_G/(V_G + V_s): a linear filter
# of the optimum (the moments lowering's recursion once V_G has settled). Its response
# to a linear trend is `steady_lag`; to a sine wave and to AR(1) noise it is below, after
# Charlesworth (1993) and Lande and Shannon (1996), in this discrete-generation form.

"""
    periodic_response(P; VG, Vs)

The steady response of the mean phenotype to an optimum that cycles with period `P`
generations, `θ_t = A sin(2π t/P)`, when the mean follows `z̄_{t+1} = z̄_t + β (θ_t − z̄_t)`
with `β = V_G/(V_G + V_s)`. With `ν = 2π/P`, the mean cycles as
`z̄_t = gain·A sin(ν t − phase)`, where

    D = √(1 − 2(1 − β) cos ν + (1 − β)²),   gain = β/D,
    phase = atan(sin ν, cos ν − (1 − β)),

and `delay = phase/ν` is the phase lag in generations. The lag `θ − z̄` cycles with
amplitude `lag_amplitude·A`, `lag_amplitude = 2 sin(ν/2)/D`. A slow cycle (`P → ∞`) is
tracked with `gain → 1` and `delay → 1/β`, the delay of the steady lag behind a trend; a
fast one is averaged out (`gain → 0`). Returns `(; gain, phase, delay, lag_amplitude)`.
"""
function periodic_response(P::Real; VG::Real, Vs::Real)
    β = VG / (VG + Vs)
    ν = 2π / P
    D = sqrt(1 - 2(1 - β) * cos(ν) + (1 - β)^2)
    phase = atan(sin(ν), cos(ν) - (1 - β))
    (gain = β / D, phase = phase, delay = phase / ν, lag_amplitude = 2sin(ν / 2) / D)
end

"""
    red_noise_lag_var(σ, ρ; VG, Vs)

The stationary variance of the lag `θ − z̄` when the optimum is AR(1) noise with
standard deviation `σ` and autocorrelation `ρ`, `θ_t = ρ θ_{t−1} + σ√(1 − ρ²) ξ_t`, and the
mean follows `z̄_{t+1} = z̄_t + β (θ_t − z̄_t)` with `β = V_G/(V_G + V_s)`:

    E[(θ − z̄)²] = 2σ² (1 − ρ) / ((2 − β)(1 − ρ (1 − β))).

White noise (`ρ = 0`) gives `2σ²/(2 − β)`, more than `σ²`, as the mean chases noise
that has already gone; slowly wandering noise (`ρ → 1`) is tracked and the lag vanishes.
"""
function red_noise_lag_var(σ::Real, ρ::Real; VG::Real, Vs::Real)
    β = VG / (VG + Vs)
    2σ^2 * (1 - ρ) / ((2 - β) * (1 - ρ * (1 - β)))
end

"""
    expected_lag2(p; VG)

The expected squared lag `E[(θ − z̄)²]` of a population with genetic variance `VG`
under the move of `p`, once the response is steady, with `V_s` =
[`selection_variance`](@ref)`(p)`: the squared [`steady_lag`](@ref) of the trend
`d = MoveFitnessSD`, plus, for `MoveType = :red_noise`, [`red_noise_lag_var`](@ref), or,
for `:periodic`, the cycle's mean squared lag `(lag_amplitude·A)²/2`
([`periodic_response`](@ref)). For `:continuous` it is the squared steady lag. The
trend and the fluctuation add because the response is linear and the fluctuation has
mean zero (over a cycle, for `:periodic`).
"""
function expected_lag2(p::Params; VG::Real)
    Vs = selection_variance(p)
    p.MoveType in (:continuous, :red_noise, :periodic) ||
        throw(ArgumentError("expected_lag2 is for an exogenous optimum (:continuous, :red_noise, :periodic), got :$(p.MoveType)"))
    trend = steady_lag(p.MoveFitnessSD; VG, Vs)^2
    p.MoveType === :red_noise && return trend + red_noise_lag_var(p.MoveNoiseSD, p.MoveAutocorr; VG, Vs)
    p.MoveType === :periodic && return trend + (periodic_response(p.MovePeriod; VG, Vs).lag_amplitude * p.MoveAmplitude)^2 / 2
    trend
end

"""
    lag_load(p; VG)

The expected lag load on the log scale: the mean of `ln W̄(0) − ln W̄(λ)` for a lag λ,
`E[λ²]/(2(ω² + V_P))` with `V_P = VG + V_E`, ω the initial width and `E[λ²]` from
[`expected_lag2`](@ref). Gaussian selection of width ω on `z ~ N(z̄, V_P)` gives
`ln W̄(λ) = ln W̄(0) − λ²/(2(ω² + V_P))` (Lande and Shannon 1996).
"""
lag_load(p::Params; VG::Real) = expected_lag2(p; VG) / (2 * (initFitSD(p)^2 + VG + envVar(p)))
