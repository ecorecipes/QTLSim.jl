# A two-stage emulator (off by default): a cheap deterministic screen (the :spectrum closure,
# by default with dynamic linkage disequilibrium and the beta solver) on every input, and the
# haplotype-block lineage model (src/block_lineages.jl) only where the screen is uncertain,
# i.e. where the input's rate lies within a band of the screen's own half-survival rate.
#
# Design: validate/two_stage_dev.jl, validate/two_stage_dev.log (development inputs only).
# Validation: validate/two_stage_criteria.md, validate/two_stage.jl, validate/two_stage.log.

"""
    QTLSim.TWO_STAGE_BAND

The default half-width of the two-stage emulator's band, in `log k`: an input goes to the
block model when its rate `k` is within a factor `exp(band)` of the screen's own half-survival
rate. Set on development inputs only (validate/two_stage_dev.log).
"""
const TWO_STAGE_BAND = 0.2

"The default screen of the two-stage emulator: `:spectrum` with dynamic LD and the beta solver."
const TWO_STAGE_SCREEN = (; ld = :dynamic, solver = :beta)

"""
    QTLSim.two_stage_screen(p; band, screen, burnin, horizon, tail) -> NamedTuple

Stage 1 of the two-stage emulator ([`emulate`](@ref)`(p; loci = :two_stage)`): the
deterministic emulator `emulate(p; loci = :spectrum, spectrum = screen)` at the input's rate
`k`, and whether `k` lies within `exp(±band)` of the screen's half-survival rate `k50`, decided
without finding `k50`: the screen is evaluated once more, at `k·e^band` if it calls the input
alive at `k` (`P2000 > 1/2`) and at `k·e^−band` if it calls it dead; the input is inside the band
iff the two calls differ (`P2000` decreases in `k`, so `k50` is then between the two rates).
With `band = 0` nothing is inside, with `band = Inf` everything; a screen trajectory that dies
before the move (`T ≤ 0`) does so at every rate (the burn-in does not depend on `k`) and is
outside without the second evaluation.

Returns `(; result, inband, probe, probe_P2000, evals, fallback)`: the screen's `emulate`
result at `k`, the decision, the probed rate (`NaN` if none) and its `P2000`, the number of
screen evaluations, and the beta solver's fallback generation (−1 if none, or if no beta).
"""
function two_stage_screen(p::Params; band::Real = TWO_STAGE_BAND, screen = TWO_STAGE_SCREEN, burnin::Integer = 2000,
                          horizon::Integer = 2000, tail::Integer = 240_000)
    band >= 0 || throw(ArgumentError("band must be non-negative, got $band"))
    e = emulate(p; loci = :spectrum, spectrum = screen, burnin, horizon, tail)
    fb = get(e, :fallback, -1)
    alive = e.P2000 > 0.5
    before = e.T !== missing && e.T <= 0
    (band == 0 || before) && return (; result = e, inband = false, probe = NaN, probe_P2000 = NaN, evals = 1, fallback = fb)
    isinf(band) && return (; result = e, inband = true, probe = NaN, probe_P2000 = NaN, evals = 1, fallback = fb)
    k = p.MoveFitnessSD * exp(alive ? band : -band)
    ep = emulate(Params(p; MoveFitnessSD = k); loci = :spectrum, spectrum = screen, burnin, horizon, tail)
    fp = get(ep, :fallback, -1)
    fb = fb < 0 ? fp : fp < 0 ? fb : min(fb, fp)
    (; result = e, inband = (ep.P2000 > 0.5) != alive, probe = k, probe_P2000 = ep.P2000, evals = 2, fallback = fb)
end

# emulate(p; loci = :two_stage, …)
function two_stage_emulate(p::Params; band::Real = TWO_STAGE_BAND, screen = TWO_STAGE_SCREEN, blocks::Integer = 160,
                           replicates::Integer = 40, seed = nothing, rng::Union{Nothing,AbstractRNG} = nothing,
                           burnin::Integer = 2000, horizon::Integer = 2000, tail::Integer = 240_000,
                           partition::Symbol = :map, within::Symbol = :gaussian)
    seed === nothing && rng === nothing && throw(ArgumentError("loci = :two_stage needs a `seed` or an `rng` (for the block model)"))
    s = two_stage_screen(p; band, screen, burnin, horizon, tail)
    e = s.result
    sv = (; inband = s.inband, screen_P2000 = e.P2000, screen_T = e.T === missing ? NaN : Float64(e.T),
          screen_probe = s.probe, screen_probe_P2000 = s.probe_P2000, screen_evals = s.evals, screen_fallback = s.fallback)
    if !s.inband
        return merge((; e.P2000, e.T, e.P240k, e.hazard, e.VG_move, e.Wbar_move, e.VG_end, e.lag_end, e.extinct,
                      alive = -1, replicates = 0, stage = 1), sv)
    end
    b = blocks_emulate(p; blocks, replicates, seed, rng, burnin, horizon, tail, partition, within)
    merge((; b.P2000, b.T, b.P240k, b.hazard, b.VG_move, b.Wbar_move, b.VG_end, b.lag_end, b.extinct,
           b.alive, b.replicates, stage = 2), sv)
end

# emulated_half_survival_rate(p; loci = :two_stage, …): the screen's k50 (bracketed as the
# deterministic emulator's), then the block model's, bisected from a bracket of ±band around it
function two_stage_half_survival_rate(p::Params; band::Real = TWO_STAGE_BAND, screen = TWO_STAGE_SCREEN, blocks::Integer = 160,
                                      replicates::Integer = 40, seed = nothing, burnin::Integer = 2000, horizon::Integer = 2000,
                                      partition::Symbol = :map, within::Symbol = :gaussian, rtol_screen::Real = 1e-3,
                                      rtol::Real = 1e-3)
    seed === nothing && throw(ArgumentError("loci = :two_stage needs a `seed` (for the block model)"))
    band > 0 && isfinite(band) || throw(ArgumentError("the half-survival rate needs a finite band > 0, got $band"))
    s = emulated_half_survival_rate(p; loci = :spectrum, spectrum = screen, burnin, horizon, rtol = rtol_screen)
    ks = s.k50
    if !(0 < ks < Inf)                     # the screen finds no root (dies before the move, or survives at 64kc)
        return (; k50 = ks, k50_screen = ks, kc = s.kc, stage = 1, evals_screen = s.evals, evals_blocks = 0, bracketed = false)
    end
    b = blocks_half_survival_rate(p; blocks, replicates, seed, burnin, horizon, partition, within, rtol,
                                  bracket = (ks * exp(-band), ks * exp(band)))
    (; k50 = b.k50, k50_screen = ks, kc = isnan(s.kc) ? b.kc : s.kc, stage = 2, evals_screen = s.evals, evals_blocks = b.evals,
     bracketed = b.bracketed)
end
