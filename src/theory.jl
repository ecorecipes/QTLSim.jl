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
