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
"""
function emulate(p::Params; loci::Symbol = :rare_alleles, burnin::Integer = 2000, horizon::Integer = 2000,
                 tail::Integer = 240_000, VG::Real = house_of_gauss_vg(p))
    p.MoveType === :continuous || throw(ArgumentError("emulate is for the published protocol, MoveType = :continuous; got :$(p.MoveType)"))
    p.NoTraits == 1 || throw(ArgumentError("emulate is for one trait; NoTraits = $(p.NoTraits)"))
    burnin >= 0 && horizon >= 2 || throw(ArgumentError("need burnin ≥ 0 and horizon ≥ 2, got $burnin and $horizon"))
    tail >= horizon || throw(ArgumentError("tail ($tail) must be at least the horizon ($horizon)"))
    q = Params(p; Sample = 1, NoGenerations = burnin + horizon, MoveFitnessSDGen = burnin)
    m = run_moments(q; VG, erosion = true, loci)
    n = length(m.gen)
    at(v, g) = g + 1 <= n ? Float64(v[g+1]) : NaN
    if n <= burnin                     # the trajectory died before the move (T ≤ 0)
        return (; P2000 = 0.0, T = m.extinct_gen - burnin + 1, P240k = 0.0, hazard = 1.0, VG_move = NaN,
                Wbar_move = NaN, VG_end = NaN, lag_end = NaN, extinct = true)
    end
    d = extinction_distribution(q, m; from = burnin)
    P2000 = d.psurvive
    pe = d.pext
    hazard = m.extinct ? 1.0 :
             length(pe) >= 2 && pe[end-1] < 1 ? (pe[end] - pe[end-1]) / (1 - pe[end-1]) : 1.0
    P240k = m.extinct || P2000 == 0 ? 0.0 : P2000 * (1 - hazard)^(tail - horizon)
    T = d.median === missing ? missing : d.median - burnin + 1
    last = burnin + horizon - 1
    (; P2000, T, P240k, hazard, VG_move = at(m.geno_var, burnin), Wbar_move = at(m.mean_fitness, burnin),
     VG_end = at(m.geno_var, last), lag_end = at(m.theta, last) - at(m.pheno_mean, last), extinct = m.extinct)
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
"""
function emulated_half_survival_rate(p::Params; loci::Symbol = :rare_alleles, burnin::Integer = 2000,
                                     horizon::Integer = 2000, bracket = nothing, rtol::Real = 1e-10,
                                     maxiter::Integer = 100)
    p.MoveType === :continuous || throw(ArgumentError("emulated_half_survival_rate is for MoveType = :continuous"))
    kc = NaN
    VG = house_of_gauss_vg(p)
    evals = 0
    # the root of P2000 = 1/2 in log k, as log(log 2) − log(−log P2000): the log of the
    # cumulative hazard is close to linear in log k, where P2000 itself is nearly a step
    function f(lk)
        evals += 1
        P = emulate(Params(p; MoveFitnessSD = exp(lk)); loci, burnin, horizon, VG).P2000
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
        kc > 0 || return (k50 = 0.0, kc, evals)
        a, b, fa = log(kc / 8), log(2kc), NaN
        fb = f(b)
        for _ in 1:5                      # half alive at 2kc: double the top, up to 64kc
            fb < 0 && break
            a, fa = b, fb
            b += log(2.0); fb = f(b)
        end
        fb >= 0 && return (k50 = Inf, kc, evals)
        isnan(fa) && (fa = f(a))
        for _ in 1:6                      # not half alive at kc/8: halve the bottom, down to kc/512
            fa >= 0 && break
            b, fb = a, fa
            a -= log(2.0); fa = f(a)
        end
        fa < 0 && return (k50 = 0.0, kc, evals)
        lo, hi, flo, fhi = a, b, fa, fb
    end
    # P2000 decreases in k: f(lo) ≥ 0 > f(hi)
    side = 0
    for _ in 1:maxiter
        hi - lo <= rtol && break
        x = lo + (hi - lo) * flo / (flo - fhi)      # regula falsi (flo ≥ 0 > fhi)
        (x <= lo || x >= hi || !isfinite(x)) && (x = (lo + hi) / 2)
        fx = f(x)
        abs(fx) < 1e-13 && return (k50 = exp(x), kc, evals)
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
    (k50 = exp((lo + hi) / 2), kc, evals)
end

"""
    emulate_sweep(ps; ids = eachindex(ps), cache = nothing, progress = 0, ntasks = Threads.nthreads(), kw...)

[`emulate`](@ref) on every `Params` in `ps`, on all threads through [`sweep`](@ref), with
its resumable CSV `cache` keyed by `ids` (strings or numbers, one per input). `kw` are
passed to `emulate`. Returns one `NamedTuple` per input, in order, with the fields of
`emulate`'s result; `T` is `NaN` where `emulate` gives `missing` (so the rows
can be cached).
"""
function emulate_sweep(ps::AbstractVector{Params}; ids = eachindex(ps), cache = nothing, progress::Integer = 0,
                       ntasks::Integer = Threads.nthreads(), kw...)
    length(ids) == length(ps) || throw(ArgumentError("one id per input: $(length(ids)) ids, $(length(ps)) inputs"))
    jobs = collect(zip(ids, ps))
    sweep(jobs; cache, id = first, progress, ntasks, cost = j -> j[2].NoQTL) do (i, p)
        r = emulate(p; kw...)
        (; P2000 = r.P2000, T = r.T === missing ? NaN : Float64(r.T), r.P240k, r.hazard, r.VG_move,
         r.Wbar_move, r.VG_end, r.lag_end, r.extinct)
    end
end
