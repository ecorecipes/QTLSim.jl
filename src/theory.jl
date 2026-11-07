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

# --- demographic stochasticity along a deterministic trajectory --------------------------

# the expected number of offspring of n survivors (the mean of `noffspring`, at most K)
expected_offspring(p::Params, n::Real) = p.Demography === :fixed ? float(p.PopSize) :
    min(float(p.PopSize), p.GrowthRate * n / (1 + (p.GrowthRate - 1) * n / p.PopSize))

# P(fewer than `nmin` of N survive), each with probability w: the population cannot mate
_too_few(N::Integer, w::Real, nmin::Integer) =
    w >= 1 ? float(N < nmin) : cdf(Binomial(N, clamp(w, 0, 1)), nmin - 1)

# Add c·pmf(k) of Binomial(N, w) to acc[k+1] for every k where the pmf is not negligible,
# by the ratio recurrence outward from the mode.
function _add_binomial!(acc, c, N::Integer, w::Real)
    w >= 1 && (acc[N+1] += c; return acc)
    w <= 0 && (acc[1] += c; return acc)
    mode = min(N, floor(Int, (N + 1) * w))
    p0 = exp(logpdf(Binomial(N, w), mode))
    ρ = w / (1 - w)
    pk = p0
    for k in mode:N
        acc[k+1] += c * pk
        k < N || break
        pk *= (N - k) / (k + 1) * ρ
        pk < 1e-17 * p0 && break
    end
    pk = p0
    for k in mode:-1:1
        pk *= k / ((N - k + 1) * ρ)
        pk < 1e-17 * p0 && break
        acc[k] += c * pk
    end
    acc
end

# Add c·P(min(K, Poisson(μ)) = k) to acc[k+1], k ≥ 1; return the mass at k = 0.
function _add_capped_poisson!(acc, c, μ::Real, K::Integer)
    μ <= 0 && return c
    mode = floor(Int, μ)
    p0 = exp(logpdf(Poisson(μ), mode))
    below = 0.0                                   # the mass placed at 1:K-1
    pk = p0
    for k in mode:K-1
        k >= 1 && (acc[k+1] += c * pk; below += pk)
        pk *= μ / (k + 1)
        pk < 1e-17 * p0 && break
    end
    pk = p0
    for k in mode-1:-1:1
        pk *= (k + 1) / μ
        pk < 1e-17 * p0 && break
        k <= K - 1 && (acc[k+1] += c * pk; below += pk)
    end
    p0zero = exp(-μ)
    acc[K+1] += c * max(0.0, 1 - below - p0zero)  # the cap: at least K offspring become K
    c * p0zero
end

# QTSim's fixed size: an independent trial each generation
function _fixed_pext(Wbar, gens, K, nmin)
    alive, out = 1.0, Float64[]
    for t in gens
        alive *= 1 - _too_few(K, Wbar[t+1], nmin)
        push!(out, 1 - alive)
    end
    out
end

# the Galton–Watson process in a varying environment, by backward iteration of the
# per-adult generating functions φ_t(s) = 1 − W̄_t + W̄_t exp(r_t (s − 1))
function _branching_pext(Wbar, N, gens, K, R)
    from = first(gens)
    φ(t, s) = (w = Wbar[t+1]; r = R / (1 + (R - 1) * N[t+1] * w / K); 1 - w + w * exp(r * (s - 1)))
    map(gens) do T
        u = 0.0
        for t in T:-1:from
            u = φ(t, u)
        end
        u^N[from+1]
    end
end

# the exact chain on sizes 0:K: binomial survival, then min(K, Poisson) offspring
function _chain_pext(Wbar, N, gens, K, R, nmin)
    prob = zeros(K + 1); prob[clamp(round(Int, N[first(gens)+1]), 0, K)+1] = 1.0
    surv = zeros(K + 1); out = Float64[]; lost = 0.0
    for t in gens
        fill!(surv, 0.0)
        for n in 0:K
            c = prob[n+1]
            c > 1e-18 && _add_binomial!(surv, c, n, Wbar[t+1])
        end
        fill!(prob, 0.0)
        lost += sum(@view surv[1:nmin])
        for n in nmin:K
            c = surv[n+1]
            c > 1e-18 || continue
            lost += _add_capped_poisson!(prob, c, R * n / (1 + (R - 1) * n / K), K)
        end
        push!(out, min(lost, 1.0))
    end
    out
end

"""
    extinction_distribution(p, Wbar; N0 = p.PopSize, from = 0, method = :chain)
    extinction_distribution(p, m; from = 0, method = :chain)

The distribution of the time to extinction from demographic stochasticity alone, along a
deterministic trajectory of mean fitness: `Wbar[t+1]` is the mean fitness `W̄_t` in
generation `t`, or `m` is the result of `run_moments(p)` (with `Sample = 1`), whose
`mean_fitness` it takes. The population starts with `N0` adults in generation 0; with
`from`, the stochastic process starts in generation `from`, with the deterministic size of
that generation, and is certain to survive before it.

In generation `t` each of the `N_t` adults survives with probability `W̄_t` (viability
selection; the binomial ignores the variation of fitness among individuals), and the
population is lost if fewer than two survive (one, if haploid). Under QTSim's `Demography
= :fixed` the survivors are replaced by `PopSize` offspring, so each generation is an
independent trial, lost with probability `P(Binomial(PopSize, W̄_t) < 2)`. Under
`:beverton_holt` the `n` survivors have `min(K, Poisson(R n/(1 + (R − 1) n/K)))` offspring
([`noffspring`](@ref)), and `method` chooses how this is followed:

- `:branching`: the density dependence is linearised along the deterministic trajectory.
  Each survivor has `Poisson(r_t)` offspring independently, with `r_t = R/(1 + (R − 1)
  n_t/K)` the per-capita rate at the expected number of survivors `n_t = N_t W̄_t`, so the
  adults form a Galton–Watson process in a varying environment, with per-adult
  generating function `φ_t(s) = 1 − W̄_t + W̄_t exp(r_t (s − 1))` (survival thins the
  Poisson). Then `P(extinct by T) = (φ_from ∘ ⋯ ∘ φ_T)(0)^{N_from}`, by backward iteration.
  This ignores the cap and the need for two parents, and, near `K`, the regulation that
  pulls the size back to the trajectory, so it overstates the spread of sizes and
  brings extinction forward (vignette 10);
- `:chain` (the default): the exact distribution of the size given the trajectory, a
  Markov chain on `0:K` driven by `W̄_t`, with the binomial survival, the cap at `K` and
  the need for two parents. It takes about a second for `K = 1000` and 2000 generations.

Returns `(; gen, pext, median, psurvive)`: the generations `from:T−1`, `P(extinct by the
end of generation t)` for each, the median generation of extinction (`missing` if `pext`
stays below one half), and the probability of surviving the whole trajectory. A
trajectory `m` that `run_moments` ends with a deterministic extinction in generation
`m.extinct_gen` is certain to be extinct by then (the last `gen`, with `pext = 1`).
"""
function extinction_distribution(p::Params, Wbar::AbstractVector{<:Real}; N0::Real = p.PopSize,
                                 from::Integer = 0, method::Symbol = :chain)
    T = length(Wbar)
    0 <= from < T || throw(ArgumentError("from must be in 0:$(T - 1), got $from"))
    method in (:branching, :chain) ||
        throw(ArgumentError("method must be :branching or :chain, got :$method"))
    method === :chain && p.Demography !== :fixed && p.PopSize > 10^6 &&
        throw(ArgumentError("the chain on 0:PopSize is too large for PopSize = $(p.PopSize); use method = :branching"))
    K, nmin = p.PopSize, p.Ploidy == 2 ? 2 : 1
    N = Vector{Float64}(undef, T)                 # the deterministic sizes, as run_moments
    N[1] = N0
    for t in 1:T-1
        N[t+1] = expected_offspring(p, N[t] * Wbar[t])
    end
    gens = from:T-1
    pext = p.Demography === :fixed ? _fixed_pext(Wbar, gens, K, nmin) :
           method === :branching ? _branching_pext(Wbar, N, gens, K, p.GrowthRate) :
           _chain_pext(Wbar, N, gens, K, p.GrowthRate, nmin)
    med = findfirst(>=(0.5), pext)
    (gen = collect(gens), pext, median = med === nothing ? missing : gens[med],
     psurvive = 1 - pext[end])
end

function extinction_distribution(p::Params, m::NamedTuple; from::Integer = 0,
                                 method::Symbol = :chain)
    m.gen == 0:length(m.gen)-1 ||
        throw(ArgumentError("the trajectory must record every generation from 0 (Sample = 1)"))
    d = extinction_distribution(p, m.mean_fitness; from, method)
    m.extinct || return d
    # the trajectory stops at its deterministic extinction, which is taken as certain
    gen, pext = [d.gen; m.extinct_gen], [d.pext; 1.0]
    (gen, pext, median = gen[findfirst(>=(0.5), pext)], psurvive = 0.0)
end

"""
    rescue_probability(; N0, R, d0, VG, VE, ω, Wmax = 1, T = 10_000)

Gomulkiewicz and Holt's (1995) evolutionary rescue after an abrupt shift of the optimum
by `d0`, on QTLSim's life cycle without density regulation. Fitness is
`Wmax·exp(−(z − θ)²/(2ω²))`, so with phenotypic variance `V_P = V_G + V_E` and lag `d_t`
the mean fitness is `W̄_t = Wmax·ω/√s · exp(−d_t²/(2s))`, `s = ω² + V_P`, and the lag
shrinks as `d_{t+1} = (1 − V_G/s) d_t` (`V_G` constant). Each adult survives with
probability `W̄_t` and has `Poisson(R)` offspring, so the population is a Galton–Watson
process with generating function `1 − W̄_t + W̄_t exp(R(u − 1))`, and the probability
that `N0` adults leave descendants for ever is one minus its extinction probability,
computed by backward iteration from the fixed point of the limiting generation (after
`T` generations). (`Wmax = 1/(0.4√(2π)) ≈ 0.997` is the peak of QTSim's fitness.)

Returns `(; rescue, Nmin, tmin)`: the probability of rescue, and the minimum of the
deterministic size `N_{t+1} = R W̄_t N_t` and the generation it is reached (Gomulkiewicz
and Holt's criterion compares `Nmin` with a critical size). `rescue = 0` when the
limiting growth `R W̄_∞ ≤ 1`.
"""
function rescue_probability(; N0::Real, R::Real, d0::Real, VG::Real, VE::Real, ω::Real,
                            Wmax::Real = 1.0, T::Integer = 10_000)
    s = ω^2 + VG + VE
    Wt(t) = Wmax * ω / sqrt(s) * exp(-(d0 * (1 - VG / s)^t)^2 / (2s))
    φ(w, u) = 1 - w + w * exp(R * (u - 1))
    W∞ = Wmax * ω / sqrt(s)
    R * W∞ <= 1 && return (rescue = 0.0, Nmin = 0.0, tmin = -1)
    q = 0.0                                       # the limiting extinction probability
    for _ in 1:10^7
        q′ = φ(W∞, q)
        abs(q′ - q) < 1e-15 && (q = q′; break)
        q = q′
    end
    for t in T-1:-1:0
        q = φ(Wt(t), q)
    end
    N, Nmin, tmin = float(N0), float(N0), 0
    for t in 0:T-1
        N *= R * Wt(t)
        N < Nmin && ((Nmin, tmin) = (N, t + 1))
    end
    (rescue = 1 - q^N0, Nmin, tmin)
end

# --- neutral differentiation among demes -------------------------------------
#
# A neutral additive trait in a metapopulation: k demes of N diploids, each running
# QTSim's generation with every fitness 1 (`Drift = true`), then `migrate(k)`. Genes are
# continuum-of-alleles: every copy of a locus gains, per generation, an independent
# increment of variance s = V_m/(2L) (the Poisson–Gaussian mutation of `mutate!`, V_m =
# `muteVar`). Two copies of a locus whose lineages have been apart for T generations
# differ by 2sT in expected square, whatever the mutation rate, so every expected
# quadratic form of the genotypes is linear in pairwise coalescence times (Slatkin 1991;
# Whitlock 1999, eqs 2–6), and copies at different loci are uncorrelated. The census is
# QTSim's record: the adults of a generation, after the previous generation's migration.

"""
    backward_migration(M, N)
    backward_migration(k, rate, N)

The backward migration matrix of the kernel's migration (`migrate!`) at a census of
equal-or-not deme sizes `N` (a number or a vector): `P[a, c]` is the probability that a
random individual of deme `a` after migration was in deme `c` before it. `M` is a
forward matrix as `migrate!` takes it ([`island_migration`](@ref),
[`stepping_stone`](@ref)); `(k, rate)` is the island model at a scalar `rate`.

It is [`QTLSim.migration_mixing`](@ref), the expected composition after migration
averaged over the random order in which `migrate!` visits the demes, exactly. For the
island model with equal deme sizes one visiting order is enough: relabelling the demes
leaves `M` unchanged, so the mean of the diagonal is the same for every order, and the
average over orders is `1 − m` on the diagonal and `m/(k − 1)` off it, with `m` that mean
emigrant fraction. Since a migrant exchanges places with an individual of another deme,
`m ≈ 2·rate` for small rates.
"""
function backward_migration(M::AbstractMatrix, N)
    k = size(M, 1)
    check_migration(M, k)
    Nv = N isa Real ? fill(float(N), k) : float.(collect(N))
    length(Nv) == k || throw(ArgumentError("$(length(Nv)) deme sizes for $k demes"))
    alive = trues(k)
    offdiag = [M[a, b] for a in 1:k for b in 1:k if a != b]
    if k >= 2 && allequal(offdiag) && allequal(Nv)          # the island model: one order
        P1 = migration_mixing(M, Nv, alive, collect(1:k))
        m = 1 - sum(P1[a, a] for a in 1:k) / k
        return [a == c ? 1 - m : m / (k - 1) for a in 1:k, c in 1:k]
    end
    migration_mixing(M, Nv, alive)
end
backward_migration(k::Integer, rate::Real, N) = backward_migration(island_migration(k, rate), N)

# One generation of the pairwise coalescence recursion, census to census, as an affine
# map on x = [W; vec(S)]: W[a] for the two genes of one adult of deme a, S[a, b] for
# genes of two distinct adults of demes a and b. Backwards from the census: migration
# (one adult came from deme c with probability P[a, c]; two distinct adults of demes a
# and b came from distinct offspring of demes c and d with probability
# J[a + (b−1)k, c + (d−1)k]), then reproduction in deme c (the genes of one offspring
# come from two distinct parents; genes of two offspring from one parent with
# probability 1/N_c, and then from the same copy with probability 1/2). Each step back
# adds one generation.
function _coalescence_step(x::AbstractVector, P::AbstractMatrix, J::AbstractMatrix, Nv::AbstractVector)
    k = length(Nv)
    W = view(x, 1:k)
    S = reshape(view(x, k+1:k+k^2), k, k)
    Ow = [1 + S[c, c] for c in 1:k]
    O = [c == d ? 1 + W[c] / (2Nv[c]) + (1 - 1 / Nv[c]) * S[c, c] : 1 + S[c, d] for c in 1:k, d in 1:k]
    vcat(P * Ow, J * vec(O))
end

"""
    exchange_origins(P, N)

The law of the deme origins of two distinct adults after the kernel's migration
(`migrate!`), to first order in the migration rates, as a `k² × k²` matrix for
[`deme_coalescence`](@ref)'s `origins`. `migrate!` exchanges two individuals of demes
`a` and `b`; with `V = (N_a P[a,b] + N_b P[b,a])/2` exchanges between them, each puts an
adult of `b` into `a` and one of `a` into `b`, so for two adults of demes `a ≠ b` the
origins `(b, a)` and `(a, b)` are each `V/(N_a N_b)` more likely, and `(a, a)` and
`(b, b)` as much less likely, than under independent origins `kron(P, P)`. The marginals
stay `P`. Chains of exchanges (an immigrant exchanged again) are of second order.
"""
function exchange_origins(P::AbstractMatrix, N)
    k = size(P, 1)
    Nv = N isa Real ? fill(float(N), k) : float.(collect(N))
    J = kron(P, P)
    ix(a, b) = a + (b - 1) * k
    for a in 1:k, b in 1:k
        a == b && continue
        δ = (Nv[a] * P[a, b] + Nv[b] * P[b, a]) / 2 / (Nv[a] * Nv[b])
        J[ix(a, b), ix(a, a)] -= δ
        J[ix(a, b), ix(b, b)] -= δ
        J[ix(a, b), ix(a, b)] += δ
        J[ix(a, b), ix(b, a)] += δ
    end
    J
end

"""
    deme_coalescence(P, N; t = nothing, origins = :pool) -> (; W, S)

Expected pairwise coalescence times, in generations, of genes in a metapopulation of
diploid demes of sizes `N` (a number or a vector) running QTSim's neutral generation and
then migration with backward matrix `P` ([`backward_migration`](@ref)), at the census
(the adults, after the previous migration). `W[a]`: the two genes of one adult of deme
`a`; `S[a, b]`: genes of two distinct adults of demes `a` and `b` (also for `a == b`).

With `t = nothing` the equilibrium, the solution of the linear fixed-point equations;
with an integer `t`, `E[min(T, t)]` at generation `t` of a run whose genomes all start
equal (as `World`s do, at zero): the recursion `x_t = A x_{t−1} + b` from `x_0 = 0`,
iterated exactly; with a sorted vector of generations, a vector of results, from one
pass.

`origins` is the law of the deme origins of two distinct adults, before migration: a
`k² × k²` matrix, `origins[a + (b−1)k, c + (d−1)k]` the probability that two distinct
adults of demes `a` and `b` after migration were in demes `c` and `d` before it; or
`:pool` (the default), `kron(P, P)`, the two origins independent, as when migrants are
drawn from a pool (the textbook island model); or `:exchange`, the first-order law of
the kernel's `migrate!` ([`exchange_origins`](@ref)). `migrate!` *exchanges*
individuals, which correlates origins: when an adult of deme `b` came from deme `a`, its
partner went from `b` to `a`. Two adults of different demes then share a deme of origin
less often than under `:pool`, by a relative `1/N`, as the exchange keeps every deme's
size: two distinct adults of the whole metapopulation are in one deme with probability
exactly `(N − 1)/(kN − 1)` (equal demes), not `1/k`. That raises the within-deme
variance from `2NkV_m` to `2NV_m (kN − 1)/(N − 1)` and Q_ST by a fraction of a percent,
at N = 50. `validate/fst.jl` estimates the law from `migrate!` itself and passes it here.
"""
function deme_coalescence(P::AbstractMatrix, N; t::Union{Nothing,Integer,AbstractVector{<:Integer}} = nothing,
                          origins::Union{Symbol,AbstractMatrix} = :pool)
    k = size(P, 1)
    size(P) == (k, k) || throw(ArgumentError("P must be square; got $(size(P))"))
    all(r -> isapprox(sum(r), 1; atol = 1e-9), eachrow(P)) || throw(ArgumentError("rows of a backward migration matrix sum to one"))
    Nv = N isa Real ? fill(float(N), k) : float.(collect(N))
    length(Nv) == k || throw(ArgumentError("$(length(Nv)) deme sizes for $k demes"))
    all(>=(2), Nv) || throw(ArgumentError("diploid demes need at least two individuals"))
    origins = origins === :pool ? kron(P, P) : origins === :exchange ? exchange_origins(P, Nv) :
              origins isa Symbol ? throw(ArgumentError("origins is :pool, :exchange or a matrix; got :$origins")) : origins
    size(origins) == (k^2, k^2) || throw(ArgumentError("origins must be $(k^2) × $(k^2); got $(size(origins))"))
    all(r -> isapprox(sum(r), 1; atol = 1e-6), eachrow(origins)) || throw(ArgumentError("rows of `origins` sum to one"))
    n = k + k^2
    unpack(x) = (W = x[1:k], S = reshape(x[k+1:end], k, k))
    if t === nothing
        b = _coalescence_step(zeros(n), P, origins, Nv)
        A = Matrix{Float64}(undef, n, n)
        e = zeros(n)
        for j in 1:n
            e[j] = 1.0
            A[:, j] = _coalescence_step(e, P, origins, Nv) .- b
            e[j] = 0.0
        end
        return unpack((I(n) - A) \ b)
    end
    ts = t isa Integer ? [t] : collect(t)
    all(>=(0), ts) && issorted(ts) || throw(ArgumentError("t must be non-negative (and sorted, for several); got $t"))
    out = Vector{NamedTuple}(undef, length(ts))
    x, now = zeros(n), 0
    for (j, tj) in enumerate(ts)
        while now < tj
            x = _coalescence_step(x, P, origins, Nv)
            now += 1
        end
        out[j] = unpack(x)
    end
    t isa Integer ? out[1] : out
end

"""
    neutral_deme_variances(P, N, Vm; t = nothing, origins = :pool)

Expected genetic variances of a neutral additive trait in a metapopulation of diploid
demes (sizes `N`, backward migration `P`, mutational variance `Vm` per individual per
generation, [`muteVar`](@ref)), at the census, from [`deme_coalescence`](@ref) (with the
same `origins` and
same `t`; for a vector of generations, a vector of these):

- `within[a]`: the expected within-deme sample variance (n − 1) of genotypes,
  `Vm (2S[a,a] − W[a])`;
- `pairs[a, b]`: half the expected squared difference of the mean genotypes of demes `a`
  and `b`, `Vm (2S[a,b] − (1 − 1/N_a) S[a,a] − W[a]/(2N_a) − (same for b))`;
- `among`: the expected sample variance (k − 1) of the k deme means, the mean of
  `pairs` over `a ≠ b`;
- `between`: the parametric between-deme variance σ²_B, `among` less the sampling
  variance of a deme mean, mean(`within ./ N`); it is `2Vm (S_ab − S_aa)` averaged over
  pairs of demes;
- `qst`: `σ²/(σ² + 2 V_W)` with `V_W` = mean(`within`) and σ² the among-deme variance
  relative to the whole metapopulation, `(k − 1)/k · between` (Wright's F_ST, as in
  Takahata 1983, Slatkin 1991 eq 8, Whitlock 1999 eqs 7–10, where Q_ST = F_ST for a
  neutral additive trait);
- `qst_demes`: the same with σ² = `between`, the variance among demes as a random sample
  of demes (an ANOVA's variance component);
- `W`, `S`: the coalescence times.
"""
function neutral_deme_variances(P::AbstractMatrix, N, Vm::Real; t::Union{Nothing,Integer,AbstractVector{<:Integer}} = nothing,
                                origins::Union{Symbol,AbstractMatrix} = :pool)
    k = size(P, 1)
    Nv = N isa Real ? fill(float(N), k) : float.(collect(N))
    c = deme_coalescence(P, Nv; t, origins)
    c isa AbstractVector ? [_deme_variances(x.W, x.S, Nv, Vm) for x in c] : _deme_variances(c.W, c.S, Nv, Vm)
end

function _deme_variances(W, S, Nv, Vm)
    k = length(W)
    within = [Vm * (2S[a, a] - W[a]) for a in 1:k]
    half(a) = (1 - 1 / Nv[a]) * S[a, a] + W[a] / (2Nv[a])
    pairs = [a == b ? 0.0 : Vm * (2S[a, b] - half(a) - half(b)) for a in 1:k, b in 1:k]
    among = k >= 2 ? sum(pairs) / (k * (k - 1)) : 0.0
    between = among - sum(within ./ Nv) / k
    VW = sum(within) / k
    σ2 = (k - 1) / k * between
    (; within, pairs, among, between, qst = σ2 / (σ2 + 2VW), qst_demes = between / (between + 2VW), W, S)
end

"""
    qst_island(N, m; k = Inf, among = :total)

The equilibrium Q_ST (= F_ST, Whitlock 1999) of a neutral additive trait in the island
model of `k` demes of `N` diploids with backward migration rate `m` (the fraction of a
deme's adults that were born elsewhere; [`backward_migration`](@ref) gives it for the
kernel's `migration_rate`), for QTSim's generation followed by migration and censused
after migration:

    Q_ST = λ / (λ + 2N (1 − λ)) · c,   λ = (1 − m k/(k − 1))²,

exactly, from the coalescence recursion of [`deme_coalescence`](@ref), where `c`
accounts for the convention: with `among = :total` (Wright's F_ST relative to the whole
metapopulation; Takahata 1983, Slatkin 1991 eq 14) the ratio `R = λ/(2N(1 − λ))` of
between- to twice the within-deme variance is scaled by `(k − 1)/k`, so that
`Q_ST = (k−1)/k·R / ((k−1)/k·R + 1)`; with `among = :demes` (variance among demes as a
sample, an ANOVA component) `Q_ST = R/(R + 1) = λ/(λ + 2N(1 − λ))`. For small `m` both
are close to Takahata's `1/(1 + 4Nm (k/(k − 1))²)` (`:total`) and to `1/(1 + 4Nm k/(k−1))`
(`:demes`); with `k = Inf` both are Wright's discrete-generation
`(1 − m)² / (2N − (2N − 1)(1 − m)²) ≈ 1/(1 + 4Nm)`. `Q_ST = 0` at panmixia,
`m = (k − 1)/k`, and → 1 as `m → 0`. Selfing is excluded (two distinct parents), which
changes the within-individual coalescence time but not Q_ST.
"""
function qst_island(N::Real, m::Real; k::Real = Inf, among::Symbol = :total)
    among in (:total, :demes) || throw(ArgumentError("among is :total or :demes; got :$among"))
    k > 1 || throw(ArgumentError("the island model needs k > 1 demes; got $k"))
    0 <= m <= 1 || throw(ArgumentError("m must be in [0, 1]; got $m"))
    λ = isinf(k) ? (1 - m)^2 : (1 - m * k / (k - 1))^2
    λ == 1 && return 1.0
    R = λ / (2N * (1 - λ))
    among === :total && !isinf(k) && (R *= (k - 1) / k)
    R / (R + 1)
end

"""
    island_variances(Vm, N, m, k)

The equilibrium within- and between-deme genetic variances of a neutral additive trait in
the island model ([`qst_island`](@ref)): `within = 2NkVm`, the variance of an undivided
population of the same total size whatever the migration rate (Lande 1992; Whitlock
1999 eq 15; Strobeck 1987 for coalescence times), `between = 2kVm λ/(1 − λ)`, the
parametric variance among deme means as a sample of demes (σ²_B of
[`neutral_deme_variances`](@ref)), and `total_among = (k−1)/k·between`, which is
Whitlock's `4NkF_ST Vm/(1 − F_ST)` (eq 15ff.); with `qst = qst_island(N, m; k)`. The
expected sample variance among the k deme means is `between + within/N`.
"""
function island_variances(Vm::Real, N::Real, m::Real, k::Integer)
    k >= 2 || throw(ArgumentError("the island model needs at least two demes; got $k"))
    λ = (1 - m * k / (k - 1))^2
    within = 2N * k * Vm
    between = λ == 1 ? Inf : 2k * Vm * λ / (1 - λ)
    (; within, between, total_among = (k - 1) / k * between, qst = qst_island(N, m; k))
end

"""
    fst_stepping_stone(N, m, k, i)

Slatkin's (1991, eq 17) approximation to the pairwise F_ST between two demes `i` steps
apart on a ring of `k` demes of `N` diploids, emigration rate `m` (half to each
neighbour): `1/(1 + 8Nmk/((k − i) i))`, from the mean coalescence times
`t_0 = 2Nk` and `t_i = 2Nk + (k − i) i/(2m)` (eq 16), for small `m`. The exact value for
the kernel's ring is [`neutral_deme_variances`](@ref) with the two demes' `pairs` entry,
as `σ² = pairs[a,b] − (within[a]/N_a + within[b]/N_b)/2` halved (k = 2 total
convention) over `σ² + 2V_W`.
"""
function fst_stepping_stone(N::Real, m::Real, k::Integer, i::Integer)
    0 <= i <= k || throw(ArgumentError("distance i must be in 0:$k; got $i"))
    (i == 0 || i == k) && return 0.0
    1 / (1 + 8N * m * k / ((k - i) * i))
end
