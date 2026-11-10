# Populations initialised from theory: genomes drawn from the stationary allele-frequency
# spectrum of mutation, stabilising selection and drift at a fixed optimum, so that a run
# need not spend thousands of generations building its standing variation from zero
# genomes. The check of the algebra is validate/symbolic/check_init.jl, the comparison with
# fully burned-in populations validate/init_equilibrium.jl.
#
# THE DENSITY. The model is that of src/spectrum.jl at the lag 0: every mutation is a
# site of its own (infinite sites, Poisson random field), in linkage equilibrium with the
# others; a site of effect a at frequency x changes by
#
#   E[Δx] = −σ x(1 − x)(1 − 2x),  σ = (a² + 2c)/(2S),      Var(Δx) = x(1 − x)/(2n),
#
# with S = ω² + V_P, c = V_G/(2L) the variance of one gene and n = N·W̄ the expected
# number of survivors, the parents (`spectrum_seff` with λ = 0, and `spectrum_step!`'s
# drift; the linked-selection term of `spectrum_linked` vanishes at λ = 0). Mutations enter
# at 2NLμ sites per generation, each at frequency 1/(2N), of effect a ~ N(0, α²). The
# stationary density of sites of effect a at frequency x above 1/(2N) is the sojourn-time
# density of the diffusion (Ewens 2004, eq. 4.23; Kimura 1969) times the input,
#
#   φ(a, x) = 4nLμ g(a) · e^{−γ x(1−x)}/(x(1 − x)) · F_γ(x),   γ = 4nσ = 2n(a² + 2c)/S,
#   F_γ(x) = ∫ₓ¹ e^{γ y(1−y)} dy / ∫₀¹ e^{γ y(1−y)} dy,
#
# to first order in the entry frequency 1/(2N) (g the density of N(0, α²)). It carries a
# constant flux towards fixation, Lμ g(a)/∫₀¹e^{γy(1−y)}dy sites per generation, and
# reduces to the neutral 4nLμ g(a)/x at γ = 0. Its variance is, by the symmetry
# F_γ(x) + F_γ(1 − x) = 1,
#
#   V = ∫∫ 2a² x(1 − x) φ = 4nLμ E_a[a² J(γ(a))],   J(γ) = ∫₀¹ e^{−γ x(1−x)} dx,
#
# in which c = V/(2L), S = ω² + V + V_E and n = N W̄(S) depend on V: `stationary_spectrum`
# solves this for V. Its limits are 2nV_m (neutral, J = 1), Turelli's house of cards
# 4LμS (γ → ∞ with c ≪ a²: J → 2/γ) and the Gaussian √(2L V_m S)-type balance when c ≫ a².
#
# THE SAMPLER (`equilibrium_genomes`). The number of sites that are at k copies among the
# 2N chromosomes of the population is Poisson with mean ∫∫ φ(a, x) Bin(k; 2N, x) dx da.
# For the neutral density that is 4nLμ/k, whatever a; so: draw a Poisson number of neutral
# sites with mean 4nLμ Σ_{k=1}^{2N−1} 1/k, give each a count k with probability ∝ 1/k, a
# frequency x ~ Beta(k, 2N − k + 1) (the law of x given k) and an effect a ~ N(0, α²), and
# keep it with probability φ/φ_neutral = e^{−γx(1−x)} F_γ(x)/(1 − x) ≤ 1. A kept site goes
# to a uniformly random locus and to k chromosomes drawn at random without replacement, and
# the allele of a chromosome at a locus is the sum of the effects of the sites it carries
# there. Sites at k = 2N (fixed) are a constant and are left out. The expected sample
# variance of the genotypes is then exactly V (hypergeometric counts per individual).
#
# WHAT IT LEAVES OUT, for a burn-in to restore:
# - linkage disequilibrium: the sites are placed independently, while stabilising
#   selection keeps a negative LD (multilocus_equilibrium's `ld`). By default the effects
#   are scaled so that the expected genic variance is `multilocus_equilibrium(p; drift =
#   true).genic` rather than the spectrum's own V, so that V_G starts above its
#   equilibrium by about −ld and comes down as the LD builds up (a few generations for
#   loosely linked loci), instead of starting at it and falling below;
# - the genealogy within a locus (a continuum-of-alleles locus's mutations are nested in
#   its alleles' genealogy, not placed independently), and the associations among linked
#   neutral sites (which widen the distribution of V_G among populations);
# - anything the model itself gets wrong (see src/spectrum.jl and src/rare_alleles.jl).
#
# THE MEAN. The mean genotype is not put on the optimum: with weak selection it wanders
# around it, as an AR(1) process z̄′ − θ = (1 − b)(z̄ − θ) + ε, b = V_G/S, Var ε = V_G/n
# (Lande 1976), with stationary variance S/(n(2 − b)), and returns only at rate b per
# generation (thousands of generations for the paper's populations). The sampler shifts
# all alleles so that the mean genotype is a draw from N(θ, S/(n(2 − b))).

# Gauss–Legendre nodes and weights on [−1, 1] (Golub–Welsch)
function gauss_legendre(K::Integer)
    A = zeros(K, K)
    for i in 1:K-1
        A[i, i+1] = A[i+1, i] = i / sqrt(4i^2 - 1)
    end
    E = eigen(Symmetric(A))
    E.values, 2 .* E.vectors[1, :] .^ 2
end

const _GL32 = gauss_legendre(32)

function _glint(f, lo::Real, hi::Real)
    x, w = _GL32
    h, m = (hi - lo) / 2, (hi + lo) / 2
    s = 0.0
    for i in eachindex(x)
        s += w[i] * f(m + h * x[i])
    end
    s * h
end

"""
    QTLSim.sojourn_mass(γ) -> J

`J(γ) = ∫₀¹ e^{−γ x(1−x)} dx`, the factor by which stabilising (underdominant) selection of
scaled strength `γ = 2n(a² + 2c)/S` reduces the variance that the segregating sites of
effect `a` hold, against the neutral spectrum (src/init.jl): `J(0) = 1`, `J ≈ 1 − γ/6` for
small `γ` and `J ≈ 2/γ` for large. By Gauss–Legendre quadrature, to about 1e-13.
"""
function sojourn_mass(γ::Real)
    γ <= 0 && return 1.0
    f(x) = exp(-γ * x * (1 - x))
    b = min(0.5, 40 / γ)            # the integrand falls from 1 by e⁻⁴⁰ within 40/γ of 0
    2 * (_glint(f, 0.0, b) + (b < 0.5 ? _glint(f, b, 0.5) : 0.0))
end

"""
    QTLSim.sojourn_ratio(γ, x) -> φ/φ_neutral

The stationary density of segregating sites at frequency `x` under stabilising selection
of scaled strength `γ` relative to the neutral one (src/init.jl),

    e^{−γ x(1−x)} F_γ(x)/(1 − x),   F_γ(x) = ∫ₓ¹ e^{γ y(1−y)} dy / ∫₀¹ e^{γ y(1−y)} dy,

by the normal distribution function (`e^{γ y(1−y)} = e^{γ/4} e^{−γ(y−½)²}`). It is 1 at
`γ = 0` and at most 1, which makes it the acceptance probability of
[`equilibrium_genomes`](@ref)'s sampler.
"""
function sojourn_ratio(γ::Real, x::Real)
    γ < 1e-10 && return 1.0
    r = sqrt(2γ)
    z1 = r / 2
    den = 1 - 2 * ccdf(Normal(), z1)
    1 - x < 1e-12 && return r * pdf(Normal(), z1) / den          # the limit x → 1
    num = ccdf(Normal(), r * (x - 0.5)) - ccdf(Normal(), z1)
    min(1.0, exp(-γ * x * (1 - x)) * num / (den * (1 - x)))
end

# E[h(|t|)] for t ~ N(0, 1), for h smooth on a logarithmic scale of t: Gauss–Legendre panels
const _EFFECT_BREAKS = (0.0, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 8.0, 10.0)
function _effect_mean(h)
    s = 0.0
    f(t) = h(t) * 2 * pdf(Normal(), t)
    for i in 1:length(_EFFECT_BREAKS)-1
        s += _glint(f, _EFFECT_BREAKS[i], _EFFECT_BREAKS[i+1])
    end
    s
end

# what the initialisation from theory covers
function check_init(p::Params)
    bad(msg) = throw(ArgumentError("initialisation from the stationary spectrum " * msg))
    p.NoTraits == 1 || bad("is for one trait; NoTraits = $(p.NoTraits)")
    p.Ploidy == 2 || bad("is for diploids; Ploidy = $(p.Ploidy)")
    p.Demography === :fixed || bad("is for a fixed population size; Demography = :$(p.Demography)")
    p.Selection === :viability || bad("is for viability selection; Selection = :$(p.Selection)")
    p.StandardisePhenotype && bad("does not cover StandardisePhenotype")
    p
end

"""
    stationary_spectrum(p) -> (; VG, nparents, S, Wbar, gene, α2, theta, L, chromosomes)

The stationary state of the allele-frequency spectrum of segregating mutations under
mutation, stabilising selection at a fixed optimum and drift, in closed form (src/init.jl):
the model of [`spectrum_equilibrium`](@ref) in continuous frequencies and effects, which
[`equilibrium_genomes`](@ref) samples. The density of sites of effect `a` at frequency `x` is

    φ(a, x) = 4nLμ g(a) e^{−γ x(1−x)} F_γ(x)/(x(1 − x)),   γ = 2n(a² + 2c)/S,

with `g` the density of `N(0, α²)`, `F_γ(x) = ∫ₓ¹e^{γy(1−y)}dy/∫₀¹e^{γy(1−y)}dy`, `n = N·W̄`
parents, `S = ω² + V_G + V_E` and `c = V_G/(2L)`; `VG` is the root of

    V = 4nLμ E_a[a² J(γ(a))],   J(γ) = ∫₀¹ e^{−γ x(1−x)} dx

(found by bisection; `n`, `S` and `γ` depend on `V`). Returns that `VG` (the expected
sample variance of the genotypes of a population drawn from the spectrum, in linkage
equilibrium), `nparents = n`, `S`, `Wbar`, `gene = c`, `α2`, the neutral scale `theta =
4nLμ`, the loci `L` and the `chromosomes`, `2N`; [`QTLSim.stationary_sites`](@ref) gives the
expected number of segregating sites.
Its limits are the neutral `2nV_m`, the house of cards `4LμS` and a Gaussian balance; it
is close to [`house_of_gauss_vg`](@ref) and to `spectrum_vg(p)`. With `Drift = true`
(no selection) it is the neutral spectrum. One diploid trait, fixed size, viability
selection.
"""
function stationary_spectrum(p::Params)
    checkparams(p); check_init(p)
    N, L, μ = p.PopSize, p.NoQTL, p.MutationRatePerGamete
    α2, VE, Ω = mutation_effect_variance(p), envVar(p), p.Drift ? Inf : initFitSD(p)^2
    function state(V)
        S = p.Drift ? Inf : Ω + V + VE
        Wbar = p.Drift ? 1.0 : QTSIM_FITNESS_PEAK * sqrt(Ω / S)
        (; S, Wbar, n = max(N * Wbar, 1.0))
    end
    γ(a2, V, s) = isinf(s.S) ? 0.0 : 2 * s.n * (a2 + V / L) / s.S
    function f(V)
        s = state(V)
        4 * s.n * L * μ * α2 * _effect_mean(t -> t^2 * sojourn_mass(γ(α2 * t^2, V, s)))
    end
    lo, hi = 0.0, 2.0 * N * muteVar(p) * 1.001 + floatmin()     # f ≤ 4nLμα² = 2nV_m ≤ 2N V_m
    V = hi
    if α2 > 0
        for _ in 1:200
            V = (lo + hi) / 2
            f(V) > V ? (lo = V) : (hi = V)
            hi - lo <= 1e-14 * hi && break
        end
    else
        V = 0.0
    end
    s = state(V)
    (; VG = V, nparents = s.n, S = s.S, Wbar = s.Wbar, gene = V / (2L), α2, theta = 4 * s.n * L * μ, L, chromosomes = 2N)
end

# γ = 2n(a² + 2c)/S of a site of effect a in the stationary spectrum st
@inline spectrum_gamma(st, a2::Real) = isinf(st.S) ? 0.0 : 2 * st.nparents * (a2 + 2 * st.gene) / st.S

"""
    QTLSim.stationary_sites(st) -> expected number of segregating sites

The expected number of sites segregating among the `2N` chromosomes of a population drawn
from the stationary spectrum `st` = [`stationary_spectrum`](@ref)`(p)`: `Σ_k ∫∫ φ(a, x)
Bin(k; 2N, x) dx da` over `k = 1, …, 2N − 1`, by the sampler's decomposition (the neutral
`theta/k` sites at `k` copies, thinned by [`QTLSim.sojourn_ratio`](@ref)); `theta·Σ 1/k`
when neutral.
"""
function stationary_sites(st)
    m = st.chromosomes
    (st.α2 > 0 && spectrum_gamma(st, st.α2) > 0) || return st.theta * sum(1 / k for k in 1:m-1)
    # Σ_k Bin(k; m, x) = 1 − (1−x)^m − x^m, against the density theta·ratio/x: integrate in ln x
    acc(x) = _effect_mean(t -> sojourn_ratio(spectrum_gamma(st, st.α2 * t^2), x))
    g(u) = (x = exp(u); (1 - (1 - x)^m - x^m) * acc(x))
    edges = range(log(1e-6 / m), 0.0; length = 49)
    st.theta * sum(_glint(g, edges[i], edges[i+1]) for i in 1:48)
end

"""
    QTLSim.stationary_density(st, a, x) -> φ(a, x)

The density of segregating sites of effect `a` at frequency `x` of the stationary
spectrum `st` = [`stationary_spectrum`](@ref)`(p)`, per unit effect and frequency, summed
over the `L` loci: `st.theta · g(a) · sojourn_ratio(γ(a), x)/x`. For checks.
"""
function stationary_density(st, a::Real, x::Real)
    st.theta * pdf(Normal(0.0, sqrt(st.α2)), a) * sojourn_ratio(spectrum_gamma(st, a^2), x) / x
end

"""
    QTLSim.stationary_mean_var(p, st = stationary_spectrum(p); VG = st.VG) -> Var(z̄ − θ)

The stationary variance of the mean genotype around a fixed optimum: the mean follows
`z̄′ − θ = (1 − b)(z̄ − θ) + ε` with `b = V_G/S` (the Gaussian response to stabilising
selection) and `Var ε = V_G/n` (drift among `n` parents), so its variance settles at
`S/(n(2 − b))`, about `(ω² + V_P)/(2N_e)` (Lande 1976). 0 without selection (`Drift =
true`), where the mean has no stationary law.
"""
function stationary_mean_var(p::Params, st = stationary_spectrum(p); VG::Real = st.VG)
    (isinf(st.S) || VG <= 0) && return 0.0
    st.S / (st.nparents * (2 - VG / st.S))
end

# the genic variance the sampler aims at
function init_target(p::Params, st, target)
    target isa Real && return float(target)
    target === :spectrum && return st.VG
    target === :genic || throw(ArgumentError("target is :genic, :spectrum or a variance; got :$target"))
    (p.Drift || st.VG <= 0) && return st.VG              # neutral: no LD expected
    multilocus_equilibrium(p; drift = true).genic
end

"""
    equilibrium_genomes(p, rng; target = :genic, mean = :stationary, spectrum = stationary_spectrum(p)) -> G

Genomes of a population near its mutation–selection–drift equilibrium at the optimum 0: a
dense `(NoQTL, 2, PopSize)` array of allele values (as [`genomes`](@ref) gives and
[`setgenomes!`](@ref) takes), drawn from the stationary allele-frequency spectrum
[`stationary_spectrum`](@ref)`(p)` (src/init.jl): a Poisson number of segregating sites,
each with a number of copies `k`, a frequency and an effect from the stationary density
(neutral proposals thinned by [`QTLSim.sojourn_ratio`](@ref)), at a random locus, carried
by `k` chromosomes drawn at random; an allele is the sum of the effects of the sites its
chromosome carries at its locus.

The sites are independent: the population is in Hardy–Weinberg and linkage equilibrium,
without the negative linkage disequilibrium that stabilising selection maintains and
without the genealogy of the alleles within a locus. A burn-in restores those (see
[`equilibrium_world`](@ref)).

- `target`: the expected genic variance. `:genic` (default) scales the effects so that it
  is `multilocus_equilibrium(p; drift = true).genic`, the genic variance of the moments
  lowering's equilibrium, above its `V_G` by the linkage disequilibrium that the burn-in
  will build; `:spectrum` leaves the spectrum's own variance, `stationary_spectrum(p).VG`;
  a number is used as it is (for example the genic variance measured in a few populations
  burned in from zero).
- `mean`: `:stationary` (default) shifts all alleles so that the mean genotype is a draw
  from its stationary law around the optimum, `N(0, S/(n(2 − V_G/S)))`
  ([`QTLSim.stationary_mean_var`](@ref)); a number puts the mean genotype there.

`spectrum` takes a `stationary_spectrum(p)` computed before, and a numeric `target` avoids
solving the moments equilibrium, for many draws of one `p` (as [`replicate`](@ref) does).

One diploid trait, fixed size, viability selection, no `StandardisePhenotype`; anything
else is an `ArgumentError`.
"""
equilibrium_genomes(p::Params, rng::AbstractRNG; target = :genic, mean = :stationary, spectrum = stationary_spectrum(p)) =
    sample_equilibrium(p, rng; target, mean, spectrum).G

# the sampler: the genomes, the number of segregating sites placed and the scale κ of their effects
function sample_equilibrium(p::Params, rng::AbstractRNG; target = :genic, mean = :stationary, spectrum = stationary_spectrum(p))
    st = spectrum
    check_init(p)
    (st.L == p.NoQTL && st.chromosomes == 2 * p.PopSize) || throw(ArgumentError("the spectrum is not that of these Params"))
    mean isa Real || mean === :stationary || throw(ArgumentError("mean is :stationary or a number; got :$mean"))
    N, L = p.PopSize, p.NoQTL
    G = zeros(L, 2, N)
    st.VG > 0 || return (; G, sites = 0, scale = 1.0)
    Vt = init_target(p, st, target)
    Vt >= 0 || throw(ArgumentError("the target genic variance must be non-negative, got $Vt"))
    κ = sqrt(Vt / st.VG)
    m = 2N
    cum = cumsum(1 / k for k in 1:m-1)
    sd = sqrt(st.α2)
    idx = collect(1:m)
    sites = 0
    for _ in 1:rand(rng, Poisson(st.theta * cum[end]))
        k = min(searchsortedfirst(cum, rand(rng) * cum[end]), m - 1)
        x = rand(rng, Beta(k, m - k + 1))
        a = sd * randn(rng)
        rand(rng) < sojourn_ratio(spectrum_gamma(st, a^2), x) || continue
        l = rand(rng, 1:L)
        sites += 1
        e = κ * a
        if 2k > m                              # common: all chromosomes but a random (m − k)-subset
            @inbounds for h in 1:m
                G[l + L * (h - 1)] += e
            end
            k, e = m - k, -e
        end
        @inbounds for j in 1:k                 # a random k-subset of the chromosomes
            r = rand(rng, j:m)
            idx[j], idx[r] = idx[r], idx[j]
            G[l + L * (idx[j] - 1)] += e
        end
    end
    z = mean isa Real ? float(mean) : sqrt(stationary_mean_var(p, st; VG = Vt)) * randn(rng)
    G .+= (z - sum(G) / N) / (2L)
    (; G, sites, scale = κ)
end

"""
    QTLSim.protocol_start(p; burnin, at = QTLSim.first_event(p)) -> g0

The generation at which a run initialised from theory starts, `max(0, at − burnin)`: it
replaces generations `0:g0−1` of the protocol, and the generations from `g0` on keep
their numbers, so that the optimum moves, the width changes and the output is sampled at
the protocol's generations. `at` is the generation by which the population must be at
equilibrium, by default [`QTLSim.first_event`](@ref).
"""
function protocol_start(p::Params; burnin::Integer, at::Integer = first_event(p))
    burnin >= 0 || throw(ArgumentError("burnin must be non-negative, got $burnin"))
    at >= 0 || throw(ArgumentError("at must be a generation, got $at"))
    max(0, at - burnin)
end

"""
    QTLSim.first_event(p)

The first generation at which the protocol of `p` changes the selection: the earlier of
`MoveFitnessSDGen` and `ChangeFitnessGen` (those that are set), or `NoGenerations` if
neither is.
"""
function first_event(p::Params)
    t = p.NoGenerations
    p.MoveFitnessSDGen != -1 && (t = min(t, p.MoveFitnessSDGen))
    p.ChangeFitnessGen != -1 && (t = min(t, p.ChangeFitnessGen))
    t
end

"""
    equilibrium_world(p, rng; burnin, at = QTLSim.first_event(p), target = :genic, mean = :stationary,
                      spectrum = stationary_spectrum(p))

A `World` started near its mutation–selection–drift equilibrium instead of from zero
genomes: genomes from [`equilibrium_genomes`](@ref), environmental deviations drawn
([`init_env!`](@ref)), and the generation counter at `g0 = max(0, at − burnin)`
([`QTLSim.protocol_start`](@ref)). Continue it with `run!(w, rng; resume = true)`, which
runs generations `g0:NoGenerations−1` of the protocol: the first `burnin` of them, before
generation `at` (by default the generation of the first change of selection, the move of
the optimum in the published protocol), are the burn-in that restores what the sampler
leaves out, and all keep the protocol's generation numbers, so output rows, the move and
the width change are where they are in a run from zero; only generations before `g0` are
missing. `simulate(p; init = :equilibrium, burnin)` does both.

The optimum starts at 0 and the width at the initial one, as in `World(p)`. Supported:
one diploid trait, fixed size, viability selection (see `equilibrium_genomes`). How much
burn-in is needed depends on the regime: `validate/init_equilibrium.jl`.
"""
function equilibrium_world(p::Params, rng::AbstractRNG; burnin::Integer, at::Integer = first_event(p),
                           target = :genic, mean = :stationary, spectrum = stationary_spectrum(p))
    g0 = protocol_start(p; burnin, at)
    w = World(p)
    setgenomes!(w, equilibrium_genomes(p, rng; target, mean, spectrum))
    init_env!(w, rng)
    w.gen = g0
    w
end
