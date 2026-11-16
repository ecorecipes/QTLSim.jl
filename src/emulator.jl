# An emulator of the published sweep (ROADMAP F4): for a published-protocol input (a
# continuous move of the optimum), the moments lowering with drift along the move plus
# demographic stochasticity along its trajectory, in place of a QTSim run. The pipeline is
# `half_survival_rate`'s (critical.jl); validate/emulator.jl tests it against the
# published table under the pre-specified criteria of validate/emulator_criteria.md.

"""
    emulate(p::Params; loci = :rare_alleles, burnin = 2000, horizon = 2000, tail = 240_000,
            VG = house_of_gauss_vg(p))

Emulate a run of the published sweep's protocol (`MoveType = :continuous`, the optimum
moving `MoveFitnessSD` per generation) without simulating it. The moments lowering
([`run_moments`](@ref), `erosion = true` for drift, closure `loci`) is run `burnin`
generations at a fixed optimum from genetic variance `VG` (the house-of-Gauss variance
with drift by default), then `horizon` generations with the optimum moving; demographic
stochasticity along that deterministic trajectory is added by
[`extinction_distribution`](@ref) from the start of the move (at a fixed size, as QTSim,
each generation is lost with probability `P(Binomial(PopSize, W̄_t) < Ploidy)`). This is
the pipeline of [`half_survival_rate`](@ref) at one rate. Requires Catlab (the moments
lowering is in the extension).

Returns a `NamedTuple`:

- `P2000`: the probability of being alive `horizon` generations after the move starts
  (named for the default horizon), i.e. that the generation of death counted as the
  published table counts it is after `horizon`;
- `T`: the median generation of extinction after the move, counted as the table's
  "Generations to extinction" (the first generation of the move is 1), or `missing` if
  the probability of extinction stays below one half within the horizon (`T ≤ 0`, and
  `P2000 = 0`, if the deterministic trajectory dies before the move, as when selection is
  too strong for the population to persist at all);
- `P240k`: the probability of surviving `tail` generations after the move, extrapolated
  from `P2000` with the per-generation hazard of the last generation of the horizon held
  constant, `P2000·(1 − h)^(tail − horizon)` (0 when the trajectory itself goes extinct);
- `hazard`: that end hazard `h`;
- `VG_move`, `Wbar_move`: the genetic variance and mean fitness of the trajectory in the
  first generation of the move (`NaN` if it died before);
- `VG_end`, `lag_end`: genetic variance and lag `θ − z̄` in the last generation (`NaN` if it
  died before);
- `extinct`: whether the deterministic trajectory itself died (its mean number of
  survivors fell below `Ploidy`), in which case extinction is taken as certain then.

How well it predicts the published runs: `validate/emulator.jl`, `validate/emulator.log`.

**Haplotype-block lineages** (`loci = :blocks`; off by default). Instead of a moments
closure, `replicates` (40 when `replicates = 0`) stochastic runs of the block model
([`run_blocks`](@ref), src/block_lineages.jl): the kernel's life cycle on `PopSize`
individuals whose chromosomes are reduced to the sums of `blocks` blocks of consecutive
loci (`partition = :map`, equal shares of QTSim's genetic map, or `:position`), with
crossovers moved to block boundaries; with one block per locus it is the kernel in law.
Each run starts `burnin` generations before the move from the kernel's equilibrium sampler
summed over the blocks. Replicate `r` draws from `Xoshiro(QTLSim.job_seed(seed, r))` (or
from seeds drawn from `rng`), so a result is a function of `seed`. `P2000` is the fraction
of runs alive `horizon` generations into the move, `T` the median of the runs' extinction
generations (counted as the table counts; `missing` if half or more survive), `hazard` the
deaths per run-generation at risk in the last quarter of the horizon (a crude end hazard),
`P240k` its extrapolation as above, and `VG_move`, `Wbar_move`, `VG_end` and `lag_end` means
over the runs alive then; the result also has `alive` and `replicates`. `VG` and `spectrum`
are ignored. Validation: `validate/blocks_dev.log` (design, convergence in `blocks`),
`validate/blocks_criteria.md`, `validate/blocks.log`.

`spectrum` is passed to [`run_moments`](@ref) for `loci = :spectrum`: e.g.
`spectrum = (; ld = :dynamic)`, the dynamic linkage disequilibrium of src/spectrum.jl
(off by default; validation `validate/spectrum_ld_criteria.md`, `validate/spectrum_ld.log`), or
`spectrum = (; solver = :beta)`, the beta closure of the spectrum (off by default; much faster,
for weak sweeps, with a deterministic fallback to the grid; [`QTLSim.AlleleSpectrum`](@ref)),
when the result also has `fallback`, the generation (from the start of the burn-in) at which the
run fell back to the grid, −1 if it did not (with replicates the earliest among them).

**Between-run variance** (`replicates = R > 0`; off by default). With `replicates = 0` the
trajectory is deterministic and only demographic chance is random, so `P2000` is nearly a
step function of the rate. With `R > 0` the burn-in is run once and `R` replicate
trajectories continue it through the move, each with its own realised mean genotype
([`QTLSim.moments_replicates`](@ref), which states the noise included and not included):
the mean drifts by the moments' own drift variance (to leading order `V_G/n + V_m/N` per
generation, `n = N·W̄` the parents), starting from the stationary spread at the move, and
selection acts on each replicate's lag. `extinction_distribution` is applied along each
replicate's `W̄`. Then `P2000` is the mean of the replicates' probabilities of being alive,
`T` the median of the mixture of their extinction-time distributions, `hazard` the
mixture's end hazard, `P240k` the mean of the replicates' extrapolations, `VG_move` the
common `V_G` at the move, `Wbar_move`, `VG_end` and `lag_end` means over the replicates
(the last two over those whose trajectory is alive at the end), and `extinct` whether every
replicate's trajectory died. The replicates draw from `rng`, or from `Xoshiro(seed)`; one of
the two is required, and the result is a function of it. Validation:
`validate/emulator_variance_criteria.md`, `validate/emulator_variance.log`.
"""
function emulate(p::Params; loci::Symbol = :rare_alleles, burnin::Integer = 2000, horizon::Integer = 2000,
                 tail::Integer = 240_000, VG::Real = house_of_gauss_vg(p), replicates::Integer = 0, vgnoise::Bool = true,
                 rng::Union{Nothing,AbstractRNG} = nothing, seed::Union{Nothing,Integer} = nothing, spectrum = (;),
                 blocks::Integer = 40, partition::Symbol = :map, within::Symbol = :snap)
    loci === :blocks && return blocks_emulate(p; blocks, replicates = replicates == 0 ? 40 : replicates, seed, rng,
                                              burnin, horizon, tail, partition, within)
    p.MoveType === :continuous || throw(ArgumentError("emulate is for the published protocol, MoveType = :continuous; got :$(p.MoveType)"))
    p.NoTraits == 1 || throw(ArgumentError("emulate is for one trait; NoTraits = $(p.NoTraits)"))
    burnin >= 0 && horizon >= 2 || throw(ArgumentError("need burnin ≥ 0 and horizon ≥ 2, got $burnin and $horizon"))
    tail >= horizon || throw(ArgumentError("tail ($tail) must be at least the horizon ($horizon)"))
    replicates >= 0 || throw(ArgumentError("replicates must be non-negative, got $replicates"))
    q = Params(p; Sample = 1, NoGenerations = burnin + horizon, MoveFitnessSDGen = burnin)
    if replicates > 0
        rng === nothing && seed === nothing &&
            throw(ArgumentError("replicates > 0 needs an `rng` or a `seed`"))
        return _emulate_replicates(q, rng === nothing ? Xoshiro(seed) : rng; replicates, loci, burnin, horizon, tail, VG, vgnoise, spectrum)
    end
    m = run_moments(q; VG, erosion = true, loci, spectrum)
    n = length(m.gen)
    at(v, g) = g + 1 <= n ? Float64(v[g+1]) : NaN
    if n <= burnin                     # the trajectory died before the move (T ≤ 0)
        return _with_fallback((; P2000 = 0.0, T = m.extinct_gen - burnin + 1, P240k = 0.0, hazard = 1.0, VG_move = NaN,
                               Wbar_move = NaN, VG_end = NaN, lag_end = NaN, extinct = true), m)
    end
    d = extinction_distribution(q, m; from = burnin)
    P2000 = d.psurvive
    pe = d.pext
    hazard = m.extinct ? 1.0 :
             length(pe) >= 2 && pe[end-1] < 1 ? (pe[end] - pe[end-1]) / (1 - pe[end-1]) : 1.0
    P240k = m.extinct || P2000 == 0 ? 0.0 : P2000 * (1 - hazard)^(tail - horizon)
    T = d.median === missing ? missing : d.median - burnin + 1
    last = burnin + horizon - 1
    _with_fallback((; P2000, T, P240k, hazard, VG_move = at(m.geno_var, burnin), Wbar_move = at(m.mean_fitness, burnin),
                    VG_end = at(m.geno_var, last), lag_end = at(m.theta, last) - at(m.pheno_mean, last), extinct = m.extinct), m)
end

# spectrum = (; solver = :beta): the emulator's result also carries `fallback`, the generation
# (from the start of the burn-in) at which the beta closure fell back to the grid, −1 if never
# (with replicates the earliest among them); the default result is unchanged
_with_fallback(r, m) = haskey(m, :spectrum_fallback) ? merge(r, (; fallback = m.spectrum_fallback)) : r
function _with_fallback(r, reps::AbstractVector)
    (isempty(reps) || !haskey(first(reps), :spectrum_fallback)) && return r
    f = [x.spectrum_fallback for x in reps if x.spectrum_fallback >= 0]
    merge(r, (; fallback = isempty(f) ? -1 : minimum(f)))
end

# emulate with between-run variance: R replicate trajectories (moments_replicates), each
# with demographic chance along its W̄, mixed with equal weights. `q` is the protocol's Params.
function _emulate_replicates(q::Params, rng::AbstractRNG; replicates, loci, burnin, horizon, tail, VG, vgnoise, spectrum = (;))
    vg = vgnoise && loci in (:rare_alleles, :multilocus, :rare_local, :auto) ? vg_fluctuation(q) : nothing
    r = moments_replicates(q; replicates, rng, from = burnin, VG, loci, vg, spectrum)
    if r.det.extinct                      # the shared burn-in died before the move (T ≤ 0)
        return (; P2000 = 0.0, T = r.det.extinct_gen - burnin + 1, P240k = 0.0, hazard = 1.0, VG_move = NaN,
                Wbar_move = NaN, VG_end = NaN, lag_end = NaN, extinct = true)
    end
    grid = burnin:burnin+horizon-1
    mix = zeros(length(grid))
    P240 = 0.0
    VGm, Wm, VGe, lage = NaN, Float64[], Float64[], Float64[]
    for rep in r.replicates
        mf = [r.det.mean_fitness; rep.mean_fitness]
        if length(mf) <= burnin               # died in the first generation of the move
            mix .+= 1.0
            continue
        end
        isnan(VGm) && (VGm = rep.geno_var[1])
        push!(Wm, rep.mean_fitness[1])
        d = extinction_distribution(q, (gen = collect(0:length(mf)-1), mean_fitness = mf, extinct = rep.extinct,
                                        extinct_gen = rep.extinct_gen); from = burnin)
        for (j, t) in enumerate(grid)
            i = searchsortedlast(d.gen, t)
            mix[j] += i == 0 ? 0.0 : d.pext[i]
        end
        pe = d.pext
        h = rep.extinct ? 1.0 : length(pe) >= 2 && pe[end-1] < 1 ? (pe[end] - pe[end-1]) / (1 - pe[end-1]) : 1.0
        P = d.psurvive
        P240 += rep.extinct || P == 0 ? 0.0 : P * (1 - h)^(tail - horizon)
        if !rep.extinct
            push!(VGe, rep.geno_var[end]); push!(lage, rep.theta[end] - rep.pheno_mean[end])
        end
    end
    R = length(r.replicates)
    mix ./= R
    P2000 = 1 - mix[end]
    hazard = mix[end-1] < 1 ? (mix[end] - mix[end-1]) / (1 - mix[end-1]) : 1.0
    med = findfirst(>=(0.5), mix)
    m(v) = isempty(v) ? NaN : sum(v) / length(v)
    _with_fallback((; P2000, T = med === nothing ? missing : grid[med] - burnin + 1, P240k = P240 / R, hazard, VG_move = VGm,
                    Wbar_move = m(Wm), VG_end = m(VGe), lag_end = m(lage), extinct = isempty(VGe)), r.replicates)
end

"""
    emulated_half_survival_rate(p::Params; loci = :rare_alleles, burnin = 2000, horizon = 2000,
                                bracket = nothing, rtol = 1e-10, maxiter = 100)

The rate of the continuous move (`MoveFitnessSD`, per generation) at which the emulator
([`emulate`](@ref)) predicts that half the populations are alive `horizon` generations
after the move starts: the root of `P2000(k) = 1/2`. The root is bracketed as in
[`half_survival_rate`](@ref), between [`drift_critical_rate`](@ref)`(p).k/8` and twice it,
extended (the top doubled up to `64kc`, the bottom halved down to `kc/512`) when that does
not bracket it: over a finite horizon the emulator can keep half the populations alive
above the critical rate,
(or by `bracket = (lo, hi)`, a narrower bracket to try first when the root is known to be
near, e.g. for finite differences; if it misses the root it is moved toward it, each time
four times wider in `log k`, up to four times, before falling back to the wide bracket), and found in `log k` by a bracketing search on
`log(log 2) − log(−log P2000)` (the log cumulative hazard, nearly linear in `log k` where
`P2000` is nearly a step): regula falsi (the Illinois variant) safeguarded by bisection,
stopped when the bracket is narrower than `rtol` in `log k` or the function is within
`1e-13` of zero. Bisection alone
(12 halvings in `half_survival_rate`) resolves `log k` to about 1e-3, too coarse for the
finite differences of [`sensitivities`](@ref), which this is meant to feed.

Returns `(; k50, kc, evals)`: the rate, `drift_critical_rate(p).k` (`NaN` when the given
`bracket` sufficed and it was not needed), and the number of emulator evaluations. `k50 = Inf` if half are alive at `64kc`, 0 if not at `kc/512` (and when `kc = 0`).
With `loci = :spectrum, spectrum = (; solver = :beta)` (the beta closure; see [`emulate`](@ref)) the
result also has `fallbacks`, the number of evaluations in which the closure fell back to the grid.

With `loci = :blocks` (the block model; see [`emulate`](@ref)) `replicates` runs (default 40)
are burned in once from `Xoshiro(QTLSim.job_seed(seed, r))` (`seed` required), every rate
continues copies of them and their generators (common random numbers), and the largest
`k` with `P2000(k) ≥ 1/2` is bisected in `log k` to `1e-3` within the same bracket as
above. `P2000(k)` is then a fraction of runs, so `k50` locates a step of a noisy function.
"""
function emulated_half_survival_rate(p::Params; loci::Symbol = :rare_alleles, burnin::Integer = 2000,
                                     horizon::Integer = 2000, bracket = nothing, rtol::Real = 1e-10,
                                     maxiter::Integer = 100, spectrum = (;), blocks::Integer = 40, replicates::Integer = 40,
                                     seed::Union{Nothing,Integer} = nothing, partition::Symbol = :map, within::Symbol = :snap)
    p.MoveType === :continuous || throw(ArgumentError("emulated_half_survival_rate is for MoveType = :continuous"))
    loci === :blocks && return blocks_half_survival_rate(p; blocks, replicates, seed, burnin, horizon, partition, within)
    kc = NaN
    VG = house_of_gauss_vg(p)
    evals = 0
    nfall = 0
    beta = loci === :spectrum && get(spectrum, :solver, :grid) === :beta
    # the root of P2000 = 1/2 in log k, as log(log 2) − log(−log P2000): the log of the
    # cumulative hazard is close to linear in log k, where P2000 itself is nearly a step
    function f(lk)
        evals += 1
        e = emulate(Params(p; MoveFitnessSD = exp(lk)); loci, burnin, horizon, VG, spectrum)
        beta && e.fallback >= 0 && (nfall += 1)
        P = e.P2000
        log(log(2.0)) - log(-log(P))           # +Inf at P = 1, −Inf at P = 0
    end
    lo = hi = flo = fhi = NaN
    found = false
    if bracket !== nothing
        a, b = log(bracket[1]), log(bracket[2])
        a < b || throw(ArgumentError("bracket must be (lo, hi) with 0 < lo < hi"))
        fa, fb = f(a), f(b)
        for _ in 1:4                      # widen a bracket that misses the root, toward it
            (fa >= 0 && fb < 0) && break
            w = 4(b - a)
            fb >= 0 ? ((a, fa) = (b, fb); b = a + w; fb = f(b)) : ((b, fb) = (a, fa); a = b - w; fa = f(a))
        end
        fa >= 0 && fb < 0 && ((lo, hi, flo, fhi) = (a, b, fa, fb); found = true)
    end
    if !found
        kc = drift_critical_rate(p).k
        kc > 0 || return _hsr_result((k50 = 0.0, kc, evals), beta, nfall)
        a, b, fa = log(kc / 8), log(2kc), NaN
        fb = f(b)
        for _ in 1:5                      # half alive at 2kc: double the top, up to 64kc
            fb < 0 && break
            a, fa = b, fb
            b += log(2.0); fb = f(b)
        end
        fb >= 0 && return _hsr_result((k50 = Inf, kc, evals), beta, nfall)
        isnan(fa) && (fa = f(a))
        for _ in 1:6                      # not half alive at kc/8: halve the bottom, down to kc/512
            fa >= 0 && break
            b, fb = a, fa
            a -= log(2.0); fa = f(a)
        end
        fa < 0 && return _hsr_result((k50 = 0.0, kc, evals), beta, nfall)
        lo, hi, flo, fhi = a, b, fa, fb
    end
    # P2000 decreases in k: f(lo) ≥ 0 > f(hi)
    side = 0
    for _ in 1:maxiter
        hi - lo <= rtol && break
        x = lo + (hi - lo) * flo / (flo - fhi)      # regula falsi (flo ≥ 0 > fhi)
        (x <= lo || x >= hi || !isfinite(x)) && (x = (lo + hi) / 2)
        fx = f(x)
        abs(fx) < 1e-13 && return _hsr_result((k50 = exp(x), kc, evals), beta, nfall)
        if fx >= 0
            lo, flo = x, fx
            side == -1 && (fhi /= 2)          # Illinois: halve the weight of an end kept twice
            side = -1
        else
            hi, fhi = x, fx
            side == 1 && (flo /= 2)
            side = 1
        end
    end
    _hsr_result((k50 = exp((lo + hi) / 2), kc, evals), beta, nfall)
end

"""
    emulate_sweep(ps; ids = eachindex(ps), cache = nothing, progress = 0, ntasks = Threads.nthreads(),
                  seed = nothing, kw...)

[`emulate`](@ref) on every `Params` in `ps`, on all threads through [`sweep`](@ref), with
its resumable CSV `cache` keyed by `ids` (strings or numbers, one per input). `kw` are
passed to `emulate`. Returns one `NamedTuple` per input, in order, with the fields of
`emulate`'s result; `T` is `NaN` where `emulate` gives `missing` (so the rows
can be cached).

With `replicates > 0` among `kw`, `seed` (an integer) is required: each input's replicates
draw from `Xoshiro(QTLSim.job_seed(seed, id))`, a fixed (FNV-1a) hash of the base seed and
the input's id, so the results do not depend on the threads, the order of the jobs or
which of them come from the cache.
"""
function emulate_sweep(ps::AbstractVector{Params}; ids = eachindex(ps), cache = nothing, progress::Integer = 0,
                       ntasks::Integer = Threads.nthreads(), seed::Union{Nothing,Integer} = nothing, kw...)
    length(ids) == length(ps) || throw(ArgumentError("one id per input: $(length(ids)) ids, $(length(ps)) inputs"))
    noisy = get(kw, :replicates, 0) > 0 || get(kw, :loci, :rare_alleles) === :blocks
    noisy && seed === nothing && throw(ArgumentError("replicates > 0 (and loci = :blocks) need a base `seed`"))
    jobs = collect(zip(ids, ps))
    sweep(jobs; cache, id = first, progress, ntasks, cost = j -> j[2].NoQTL) do (i, p)
        r = noisy ? emulate(p; seed = job_seed(seed, i), kw...) : emulate(p; kw...)
        (; P2000 = r.P2000, T = r.T === missing ? NaN : Float64(r.T), r.P240k, r.hazard, r.VG_move,
         r.Wbar_move, r.VG_end, r.lag_end, r.extinct, (haskey(r, :fallback) ? (; fallback = r.fallback) : (;))...,
         (haskey(r, :alive) ? (; alive = r.alive) : (;))...)
    end
end

"""
    QTLSim.job_seed(seed, id) -> UInt64

A seed for job `id` of a sweep with base `seed`: the 64-bit FNV-1a hash of the bytes of
`string(seed, "/", id)`. Fixed across Julia versions and platforms (unlike `hash`).
"""
function job_seed(seed::Integer, id)
    h = 0xcbf29ce484222325
    for b in codeunits(string(seed, "/", id))
        h = (h ⊻ b) * 0x00000100000001b3
    end
    h
end

# ∫₀¹ w(x) e^{−γx(1−x)} dx for w symmetric about ½, as sojourn_mass integrates it
function _sojourn_moment(γ::Real, w)
    f(x) = w(x) * exp(-max(γ, 0.0) * x * (1 - x))
    b = γ <= 0 ? 0.5 : min(0.5, 40 / γ)
    2 * (_glint(f, 0.0, b) + (b < 0.5 ? _glint(f, b, 0.5) : 0.0))
end

"""
    QTLSim.vg_fluctuation(p, st = stationary_spectrum(p)) -> (; VG, stationary, innovation, cv)

The random fluctuation of the genetic variance of one population (genic, in linkage
equilibrium) under the stationary allele-frequency spectrum `st` of
[`stationary_spectrum`](@ref) (src/init.jl): the infinite-sites Poisson random field in
which a site of effect `a` at frequency `x` contributes `2a²x(1 − x)` to `V_G`, with density
`φ(a, x) = 4nLμ g(a) e^{−γx(1−x)} F_γ(x)/(x(1 − x))`, `γ = 2n(a² + 2c)/S`.

- `stationary`: `Var V_G = ∫∫ (2a²x(1−x))² φ = 8nLμ E_a[a⁴ K(γ)]`, `K(γ) = ∫₀¹ x(1−x)e^{−γx(1−x)}dx`
  (sites are Poisson, so the variance is the integral of the square; `F_γ(x) + F_γ(1−x) = 1`
  halves it);
- `innovation`: the variance of the change of `V_G` in one generation by drift of the
  frequencies, `Σ (2a²)²(1 − 2x)² Var(Δx)` with `Var(Δx) = x(1 − x)/(2n)`,
  `= 4Lμ E_a[a⁴ M(γ)]`, `M(γ) = ∫₀¹ (1 − 2x)² e^{−γx(1−x)}dx`;
- `cv = √stationary/VG`.

Neutral (`γ = 0`): `stationary = 4nLμα⁴`, `cv² = 1/(4nLμ)`, and `innovation = 4Lμα⁴`, so a
deviation decaying at `1/(2n)` per generation (heterozygosity under drift) keeps the
stationary variance. Leaves out linkage disequilibrium, associations among linked sites
(which widen the law of V_G) and the input of new sites (`O(1/N)` of the innovation). Used
by [`emulate`](@ref)`(p; replicates)` for the noise in `V_G`.
"""
function vg_fluctuation(p::Params, st = stationary_spectrum(p))
    α2, L, μ, n = st.α2, st.L, p.MutationRatePerGamete, st.nparents
    α2 > 0 || return (; VG = st.VG, stationary = 0.0, innovation = 0.0, cv = 0.0)
    γa(t) = spectrum_gamma(st, α2 * t^2)
    stat = 8 * n * L * μ * α2^2 * _effect_mean(t -> t^4 * _sojourn_moment(γa(t), x -> x * (1 - x)))
    innov = 4 * L * μ * α2^2 * _effect_mean(t -> t^4 * _sojourn_moment(γa(t), x -> (1 - 2x)^2))
    (; VG = st.VG, stationary = stat, innovation = innov, cv = sqrt(stat) / st.VG)
end

# spectrum = (; solver = :beta): also `fallbacks`, the number of evaluations in which the closure fell back to the grid
_hsr_result(r, beta, nfall) = beta ? merge(r, (; fallbacks = nfall)) : r
