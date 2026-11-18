# The gradient of the emulator's half-survival rate (AD plan, phase B): a target that does not
# saturate, and the elasticities of k50 by the implicit-function formula on it.
#
# Why `emulate`'s P2000 cannot be differentiated near k50 when selection is strong (narrow ω):
# the moments trajectory then sits at a steady lag with a hazard far below 1e-16 per generation
# until it runs away, and it dies within a few generations of doing so. P2000 as a function of the
# rate is then (i) exactly 1 in Float64 while the cumulative hazard H is below 1.1e-16 (the product
# of the (1 − hazard) factors does not move; H down to 1e-100 was seen), and (ii) exactly 0 once the
# deterministic trajectory itself dies inside the horizon (its expected survivors fall below
# `Ploidy`, which `extinction_distribution` takes as certain extinction), although the hazard
# accumulated until then can be far below log 2 (0.03 was seen). Between the two, log H changes
# by its whole range over a relative change of the rate of 1e-3 to 1e-6.
# Not a cause: the truncation of the size chain (1e-17, 1e-18), which is not used at a fixed size.
#
# The target here: log H, the log cumulative hazard over the horizon, summed in logs (so it is
# finite from H = 1e-300 and below up to H of thousands), along the trajectory continued through
# its deterministic extinction (a quadratic extrapolation of log W̄). It is a smooth function of
# the rate and the parameters; its level set log H = log log 2 is the emulator's P2000 = 1/2
# exactly when the trajectory survives the horizon there, and within a fraction of a generation
# of the runaway otherwise. Validation: validate/elasticity_criteria.md (phase B addendum),
# validate/elasticity_b.jl, validate/elasticity_b.log.

# log P(fewer than `nmin` of K survive), each with probability w, in closed form (nmin 1 or 2):
# the log of `_too_few`, without its underflow
function _log_too_few(K::Integer, w::Real, nmin::Integer)
    w <= 0 && return 0.0
    w >= 1 && return K < nmin ? 0.0 : -Inf
    nmin == 1 ? K * log1p(-w) : (K - 1) * log1p(-w) + log1p((K - 1) * w)
end

# log of the integrated hazard η = −log(1 − h) of one generation, h = P(fewer than nmin survive)
function _log_generation_hazard(K::Integer, w::Real, nmin::Integer)
    lh = _log_too_few(K, w, nmin)
    lh == -Inf && return -Inf
    lh < -36.0 && return lh                         # η = h + h²/2 + …
    if lh > -1e-8                                   # h → 1: 1 − h = P(at least nmin survive), to leading order in K w
        w <= 0 && return Inf
        l1mh = nmin == 1 ? log(K) + log(w) : log(K) + log(K - 1) - log(2.0) + 2 * log(w)
        return log(-l1mh)
    end
    l1mh = lh > -log(2.0) ? log(-expm1(lh)) : log1p(-exp(lh))
    log(-l1mh)
end

"""
    QTLSim.log_cumulative_hazard(p::Params, Wbar; from = 0, generations = length(Wbar) - from)

The log of the cumulative hazard of extinction from demographic stochasticity along a trajectory
of mean fitness at a fixed size (QTSim's `Demography = :fixed`), for the first `generations`
generations from `from`: `log H`, `H = Σ_t −log(1 − h_t)`, `h_t = P(Binomial(PopSize, W̄_t) <`
`Ploidy)`, so that the probability of surviving them is `exp(−H)` (what
[`extinction_distribution`](@ref) returns as `psurvive`). Everything is carried in logs
(`log h_t` in closed form, `log1p`/`expm1`, a running log-sum-exp), so the result is finite and
smooth where `psurvive` is exactly 1 (`H < 1.1e-16`) or exactly 0 in `Float64`.

If `Wbar` ends before `from + generations` (a moments trajectory stops at its deterministic
extinction, when the expected survivors fall below `Ploidy`), it is **continued**: `log W̄` is
extrapolated quadratically from its last three values, with its first and second differences
capped at zero (so it keeps falling at least as fast). This replaces `extinction_distribution`'s
convention that the deterministic extinction is certain, which makes `psurvive` jump to 0; the
continuation matters only for the first generation or two after the trajectory ends, after which
`H` is far above any level of interest.

Returns `(; logH, continued)`: the log cumulative hazard and the number of generations taken
from the continuation. Only for `Demography = :fixed`.
"""
function log_cumulative_hazard(p::Params, Wbar::AbstractVector{<:Real}; from::Integer = 0,
                               generations::Integer = length(Wbar) - from)
    p.Demography === :fixed ||
        throw(ArgumentError("log_cumulative_hazard is for Demography = :fixed, got :$(p.Demography)"))
    n = length(Wbar)
    0 <= from < n || throw(ArgumentError("from must be in 0:$(n - 1), got $from"))
    generations >= 1 || throw(ArgumentError("generations must be positive, got $generations"))
    K, nmin = p.PopSize, p.Ploidy == 2 ? 2 : 1
    acc = -Inf                                      # log H so far
    add(l) = l == -Inf ? acc : acc == -Inf ? l : (a = max(acc, l); a == Inf ? Inf : a + log1p(exp(min(acc, l) - a)))
    last = min(n, from + generations)
    for t in from+1:last
        acc = add(_log_generation_hazard(K, Float64(Wbar[t]), nmin))
    end
    continued = from + generations - last
    if continued > 0
        lw = log(max(Float64(Wbar[n]), 0.0))
        d1 = n >= 2 ? min(lw - log(Float64(Wbar[n-1])), 0.0) : 0.0
        d2 = n >= 3 ? min(lw - 2 * log(Float64(Wbar[n-1])) + log(Float64(Wbar[n-2])), 0.0) : 0.0
        for _ in 1:continued
            (acc == Inf || !isfinite(lw)) && (acc = Inf; break)
            d1 += d2
            lw += d1
            acc = add(_log_generation_hazard(K, exp(lw), nmin))
        end
    end
    (; logH = acc, continued)
end

"""
    QTLSim.emulated_log_hazard(p::Params; loci = :spectrum, spectrum = QTLSim.TWO_STAGE_SCREEN,
                               burnin = 2000, horizon = 2000, VG = house_of_gauss_vg(p))

The deterministic emulator ([`emulate`](@ref) with `replicates = 0`, a moments closure `loci`)
with its probability of being alive `horizon` generations into the move returned in logs, so
that it does not saturate: the same moments trajectory, and [`QTLSim.log_cumulative_hazard`](@ref)
along it from the start of the move.

Returns `(; logH, h, logP, log1mP, P2000, extinct, continued, fallback)`:

- `logH`: the log of the cumulative hazard `H` over the horizon; `logP = −H` is the log of the
  probability of being alive and `log1mP = log(1 − e^{−H})` that of being extinct;
- `h = log(log 2) − logH`, the function whose root in the rate is the half-survival rate (as in
  [`emulated_half_survival_rate`](@ref): positive where more than half are alive, decreasing in
  the rate);
- `P2000 = exp(−H)`. It equals `emulate`'s `P2000` (to rounding) when the deterministic
  trajectory survives the horizon. When the trajectory dies inside it (`extinct`), `emulate`
  gives exactly 0 and this the value along the continued trajectory (`continued` generations of
  it), usually tiny but not always: the hazard accumulated before a sudden runaway can be small;
- `fallback`: the beta solver's fallback generation (−1 if none, or no beta solver).

A trajectory that dies before the move has `logH = Inf` (`P2000 = 0`). Requires Catlab; only for
`MoveType = :continuous`, one trait, `Demography = :fixed`.
"""
function emulated_log_hazard(p::Params; loci::Symbol = :spectrum, spectrum = TWO_STAGE_SCREEN, burnin::Integer = 2000,
                             horizon::Integer = 2000, VG::Real = house_of_gauss_vg(p))
    p.MoveType === :continuous || throw(ArgumentError("emulated_log_hazard is for MoveType = :continuous; got :$(p.MoveType)"))
    p.NoTraits == 1 || throw(ArgumentError("emulated_log_hazard is for one trait; NoTraits = $(p.NoTraits)"))
    p.Demography === :fixed || throw(ArgumentError("emulated_log_hazard is for Demography = :fixed, got :$(p.Demography)"))
    burnin >= 0 && horizon >= 2 || throw(ArgumentError("need burnin ≥ 0 and horizon ≥ 2, got $burnin and $horizon"))
    q = Params(p; Sample = 1, NoGenerations = burnin + horizon, MoveFitnessSDGen = burnin)
    m = run_moments(q; VG, erosion = true, loci, spectrum)
    fallback = haskey(m, :spectrum_fallback) ? Int(m.spectrum_fallback) : -1
    if length(m.gen) <= burnin                 # died before the move
        return (; logH = Inf, h = -Inf, logP = -Inf, log1mP = 0.0, P2000 = 0.0, extinct = true, continued = horizon, fallback)
    end
    c = log_cumulative_hazard(q, m.mean_fitness; from = burnin, generations = horizon)
    H = exp(c.logH)
    log1mP = H > log(2.0) ? log1p(-exp(-H)) : log(-expm1(-H))
    (; logH = c.logH, h = log(log(2.0)) - c.logH, logP = -H, log1mP, P2000 = exp(-H), extinct = m.extinct,
     continued = c.continued, fallback)
end

# one evaluation of the target for k50_gradient: h and the fallback generation at log k and
# the named fields multiplied by exp(shift)
# (the starting variance of the burn-in is each point's own, as in `emulate`)
_kg_point(p::Params, lk::Real, kw; field::Union{Nothing,Symbol} = nothing, shift::Real = 0.0) =
    emulated_log_hazard(field === nothing ? Params(p; MoveFitnessSD = exp(lk)) :
                        Params(p; MoveFitnessSD = exp(lk), (field => getfield(p, field) * exp(shift),)...); kw...)

"""
    k50_gradient(p::Params; loci = :spectrum, spectrum = QTLSim.TWO_STAGE_SCREEN,
                 names = (:PhenotypeVar, :Heritability, :MutationScalar, :InitFitnessSD),
                 rel = 1e-2, dh = 0.5, htol = 0.02, bracket = nothing, onswitch = :grid,
                 burnin = 2000, horizon = 2000, rtol = 1e-2, minstep = 1e-10, maxevals = 200)

The emulator's half-survival rate `k50` (the rate of the continuous move at which half the
populations are alive `horizon` generations into it; [`emulated_half_survival_rate`](@ref)) and
its elasticities `∂log k50/∂log θ` with respect to the `Params` fields `names` (positive
`Float64` fields, as for [`QTLSim.Para`](@ref)), each field varied alone.

**Target.** `h(k, θ) = log(log 2) − log H`, the log cumulative hazard over the horizon carried in
logs along the continued trajectory ([`QTLSim.emulated_log_hazard`](@ref)): finite and smooth
across the whole neighbourhood of `k50`, where `emulate`'s `P2000` is exactly 0 or 1 in `Float64`
on one side or both when selection is strong. `k50` is the root of `h` in `log k`; the
elasticities are the implicit-function formula `e_θ = −(∂h/∂log θ)/(∂h/∂log k)` there, each
derivative a central difference.

**Method.** (1) A bracket: `bracket = (lo, hi)` if given (moved toward the root, each time four
times wider, up to four times, if it misses it), else `emulated_half_survival_rate` to `rtol` in
`log k`. (2) The root of `h` in the bracket, by regula falsi (Illinois) where `h` is finite and
moderate and bisection elsewhere, until `|h| ≤ htol`. (3) The step: `h` can change by 1 over a
relative change of the rate of 1e-2 to 1e-6, so the step in `log k` is chosen for the central
difference to span about `±dh` in `h` (never more than `rel`; at least `minstep`), and the same
step is used for each parameter, reduced for one whose difference spans more than `±4dh` (or
which would leave its range, as `Heritability > 1`). (4) **The beta solver's fallback** to the
grid is a discontinuity of the trajectory in the generation at which it happens. A difference is
valid only if its two points and the base point fell back in the same generation (or none did);
otherwise its step is divided by four, up to three times. If a pair still differs, then with
`onswitch = :grid` (and `spectrum.solver = :beta`) the whole computation is repeated with the grid
solver, from a bracket around the beta's `k50`, and that is the result (`solver = :grid`, with the
beta's values in `beta`); with `onswitch = :flag` the result is the beta's, with `status =
:fallback_switch`.

Returns a `NamedTuple`:

- `k50`; `elasticities`, a `NamedTuple` keyed by `names` (`NaN` where undefined);
- `defined` (`status === :ok`) and `status`: `:ok`; `:no_k50` (no root: `k50` is 0 or `Inf`);
  `:no_root` (the bracket closed without `|h| ≤ htol`: a discontinuity at the root);
  `:not_finite` (a difference point not finite at the smallest step); `:not_decreasing`
  (`∂h/∂log k ≥ 0`); `:fallback_switch`; `:max_evals`;
- `h0` (`h` at `k50`), `dh_dlogk`, `step` (in `log k`), `steps` (per parameter), `asymmetry`
  (`|h₊ + h₋ − 2h0|/|h₊ − h₋|` of the rate's difference: the relative size of the curvature
  over the step);
- `fallback`: the fallback generation at `k50` and at the two points of each difference,
  `(; base, k = (plus, minus), <name> = (plus, minus), …)` (−1 = none); `straddle` (a pair, or a
  pair and the base, of which one fell back and another did not, at the final steps) and
  `shift` (all fell back, in different generations);
- `continued`: whether the trajectory of the base or a difference point died inside the horizon,
  so that the target relied on the continuation (then `emulate`'s `P2000` jumps to 0 next to
  `k50`, and its own half-survival rate is that jump, within a fraction of a generation of the
  runaway from this `k50`);
- `solver` (`:beta`, `:grid` or `:none` for other closures), `beta` (the beta solver's result
  when the grid replaced it, else `nothing`), `evals` (emulator evaluations, all included).

**Limits.** Finite differences, not automatic differentiation (see [`sensitivities`](@ref)).
`k50` is the root of the smooth target, not of `emulate`'s `P2000` where they differ (see
`continued`). The continuation is an extrapolation, used over one or two generations. Each field
is varied alone: fields that the emulator reads only through a combination have dependent
elasticities (it reads `PhenotypeVar`, `Heritability` and `MutationScalar` only through
`V_E = (1 − h²)V_P` and `m·V_E`, so `e_h² = −h²/(1 − h²)·e_V_P` exactly, and `e_V_P − e_m` is the
elasticity in `V_E` at a fixed mutational variance, small when `ω² ≫ V_E`). A deterministic
closure only (`replicates = 0`); `Demography = :fixed`. Requires Catlab.
"""
function k50_gradient(p::Params; loci::Symbol = :spectrum, spectrum = TWO_STAGE_SCREEN,
                      names = (:PhenotypeVar, :Heritability, :MutationScalar, :InitFitnessSD),
                      rel::Real = 1e-2, dh::Real = 0.5, htol::Real = 0.02, bracket = nothing, onswitch::Symbol = :grid,
                      burnin::Integer = 2000, horizon::Integer = 2000, rtol::Real = 1e-2, minstep::Real = 1e-10,
                      maxevals::Integer = 200)
    onswitch in (:grid, :flag) || throw(ArgumentError("onswitch must be :grid or :flag, got :$onswitch"))
    rel > 0 && dh > 0 && htol > 0 && 0 < minstep <= rel ||
        throw(ArgumentError("need rel, dh, htol > 0 and 0 < minstep ≤ rel"))
    loci in (:blocks, :two_stage) && throw(ArgumentError("k50_gradient needs a deterministic closure, not loci = :$loci"))
    P = Para(identity, p, names)                    # checks the names (Float64 fields of Params)
    for (n, θ) in zip(P.names, P.x0)
        θ > 0 || throw(ArgumentError("an elasticity needs a positive parameter; $n = $θ"))
    end
    isbeta = loci === :spectrum && get(spectrum, :solver, :grid) === :beta
    r = _k50_gradient(p, P.names, loci, spectrum, isbeta ? :beta : loci === :spectrum ? :grid : :none;
                      rel, dh, htol, bracket, burnin, horizon, rtol, minstep, maxevals)
    (r.status === :fallback_switch && onswitch === :grid && isbeta) || return r
    grid = NamedTuple(k => v for (k, v) in pairs(spectrum) if k !== :solver)
    g = _k50_gradient(p, P.names, loci, grid, :grid; rel, dh, htol, bracket = (r.k50 * exp(-0.02), r.k50 * exp(0.02)),
                      burnin, horizon, rtol, minstep, maxevals)
    merge(g, (; beta = r, evals = g.evals + r.evals))
end

function _k50_gradient(p::Params, names::Vector{Symbol}, loci, spectrum, solver; rel, dh, htol, bracket, burnin, horizon, rtol,
                       minstep, maxevals)
    kw = (; loci, spectrum, burnin, horizon)
    evals = 0
    nan = NamedTuple{Tuple(names)}(ntuple(_ -> NaN, length(names)))
    nofb = merge((; base = -1, k = (-1, -1)), NamedTuple{Tuple(names)}(ntuple(_ -> (-1, -1), length(names))))
    result(status; k50 = NaN, e = nan, h0 = NaN, hk = NaN, step = NaN, steps = nan, asym = NaN, fb = nofb, straddle = false,
           shift = false, continued = false) =
        (; k50, elasticities = e, defined = status === :ok, status, h0, dh_dlogk = hk, step, steps, asymmetry = asym,
         fallback = fb, straddle, shift, continued, solver, beta = nothing, evals)
    ev(lk; field = nothing, shift = 0.0) = (evals += 1; _kg_point(p, lk, kw; field, shift))
    # (1) a bracket a < b in log k with h(a) ≥ 0 > h(b)
    if bracket === nothing
        s = emulated_half_survival_rate(p; loci, spectrum, burnin, horizon, rtol)
        evals += s.evals
        0 < s.k50 < Inf || return result(:no_k50; k50 = s.k50)
        a, b = log(s.k50) - 0.51rtol, log(s.k50) + 0.51rtol
    else
        a, b = log(bracket[1]), log(bracket[2])
        a < b || throw(ArgumentError("bracket must be (lo, hi) with 0 < lo < hi"))
    end
    ea, eb = ev(a), ev(b)
    for _ in 1:4
        (ea.h >= 0 && eb.h < 0) && break
        w = 4(b - a)
        eb.h >= 0 ? ((a, ea) = (b, eb); b = a + w; eb = ev(b)) : ((b, eb) = (a, ea); a = b - w; ea = ev(a))
    end
    (ea.h >= 0 && eb.h < 0) || return result(:no_k50; k50 = eb.h >= 0 ? Inf : 0.0)
    # (2) the root of h: regula falsi (Illinois) where h is finite and moderate, else bisection
    x, e0 = a, ea
    abs(eb.h) < abs(ea.h) && ((x, e0) = (b, eb))
    fa, fb, side = ea.h, eb.h, 0
    slope = NaN
    while abs(e0.h) > htol
        (b - a < 1e-14 || evals >= maxevals) && return result(evals >= maxevals ? :max_evals : :no_root; k50 = exp((a + b) / 2), h0 = e0.h)
        usable = isfinite(fa) && isfinite(fb) && max(abs(fa), abs(fb)) < 30
        x = usable ? a + (b - a) * fa / (fa - fb) : (a + b) / 2
        (x <= a || x >= b || !isfinite(x)) && (x = (a + b) / 2)
        e0 = ev(x)
        if e0.h >= 0
            isfinite(e0.h) && isfinite(eb.h) && (slope = (eb.h - e0.h) / (b - x))
            a, ea, fa = x, e0, e0.h
            side == -1 && (fb /= 2); side = -1
        else
            isfinite(e0.h) && isfinite(ea.h) && (slope = (e0.h - ea.h) / (x - a))
            b, eb, fb = x, e0, e0.h
            side == 1 && (fa /= 2); side = 1
        end
    end
    isnan(slope) && isfinite(ea.h) && isfinite(eb.h) && b > a && (slope = (eb.h - ea.h) / (b - a))
    lk0, h0, fb0 = x, e0.h, e0.fallback
    k50 = exp(lk0)
    continued = e0.extinct
    same(u, v) = u.fallback == v.fallback
    # (3) the step in log k: a central difference spanning about ±dh in h
    δ = isfinite(slope) && slope < 0 ? clamp(dh / abs(slope), minstep, rel) : rel
    ep = em = e0
    ok = false
    cuts = 0
    for _ in 1:12
        evals >= maxevals && return result(:max_evals; k50, h0)
        ep, em = ev(lk0 + δ), ev(lk0 - δ)
        fin = isfinite(ep.h) && isfinite(em.h)
        span = fin ? abs(ep.h - em.h) / 2 : Inf
        switch = !(same(ep, e0) && same(em, e0))
        if fin && span <= 3dh && !switch
            ok = true; break
        end
        δ <= minstep && break
        if fin && span <= 3dh                       # only a fallback switch: a smaller step, up to three times
            cuts += 1
            cuts > 3 && break
            δ = max(δ / 4, minstep)
        else
            δ = max(fin ? δ * dh / span : δ / 8, minstep)
        end
    end
    fbs = Dict{Symbol,Tuple{Int,Int}}(:k => (ep.fallback, em.fallback))
    all3(u, v) = (fb0, u, v)
    function flags(fbs)
        st = sh = false
        for (u, v) in values(fbs)
            t = all3(u, v)
            any(>=(0), t) && any(<(0), t) && (st = true)
            all(>=(0), t) && !(t[1] == t[2] == t[3]) && (sh = true)
        end
        st, sh
    end
    fbnt(fbs) = merge((; base = fb0, k = fbs[:k]), NamedTuple{Tuple(names)}(Tuple(get(fbs, n, (-1, -1)) for n in names)))
    if !ok
        fin = isfinite(ep.h) && isfinite(em.h)
        st, sh = flags(fbs)
        return result(!fin ? :not_finite : (st || sh) ? :fallback_switch : :not_finite; k50, h0, step = δ, fb = fbnt(fbs),
                      straddle = st, shift = sh, continued)
    end
    continued |= ep.extinct | em.extinct
    hk = (ep.h - em.h) / 2δ
    asym = abs(ep.h + em.h - 2h0) / abs(ep.h - em.h)
    hk < 0 || return result(:not_decreasing; k50, h0, hk, step = δ, asym, fb = fbnt(fbs), continued)
    # (4) the parameters, with the same step unless a difference is too wide or straddles a switch
    es = Float64[]; ss = Float64[]
    status = :ok
    for n in names
        θ = getfield(p, n)
        d = n === :Heritability && θ < 1 ? min(δ, -log(θ) / 2) : δ
        val = NaN
        cuts = 0
        for _ in 1:12
            evals >= maxevals && (status = :max_evals; break)
            up, um = ev(lk0; field = n, shift = d), ev(lk0; field = n, shift = -d)
            fbs[n] = (up.fallback, um.fallback)
            fin = isfinite(up.h) && isfinite(um.h)
            span = fin ? abs(up.h - um.h) / 2 : Inf
            switch = !(same(up, e0) && same(um, e0))
            if fin && span <= 4dh && !switch
                val = -((up.h - um.h) / 2d) / hk
                continued |= up.extinct | um.extinct
                break
            end
            if d <= minstep
                status === :ok && (status = switch && fin ? :fallback_switch : :not_finite); break
            end
            if fin && span <= 4dh
                cuts += 1
                cuts > 3 && (status === :ok && (status = :fallback_switch); break)
                d = max(d / 4, minstep)
            else
                d = max(fin ? d * dh / span : d / 8, minstep)
            end
        end
        push!(es, val); push!(ss, d)
        status === :max_evals && break
    end
    while length(es) < length(names)
        push!(es, NaN); push!(ss, NaN)
    end
    st, sh = flags(fbs)
    status === :ok && any(isnan, es) && (status = :not_finite)
    T = NamedTuple{Tuple(names)}
    result(status; k50, e = T(Tuple(es)), h0, hk, step = δ, steps = T(Tuple(ss)), asym, fb = fbnt(fbs), straddle = st,
           shift = sh, continued)
end
