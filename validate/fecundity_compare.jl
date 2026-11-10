# Viability against fecundity selection in the published sweep, paired by input and seed.
#
# Reads the caches of validate/demography_sweep.jl, run with SELECTION = viability (QTSim's)
# and fecundity at a fixed size and at one or more growth rates, and compares the two
# selection modes on the inputs both ran: the share alive CAP generations after the move,
# the median time to extinction, and the paired ratio of times to extinction (fecundity /
# viability), and each mode's time to extinction against QTSim's fixed-size (viability)
# run of the same input. Inputs are compared in two groups, NoQTL ≤ 1000 and NoQTL = 5000,
# each with the inputs that every set compared ran. Strata (NoQTL × u) count equally, as
# in demography_sweep.jl. Confidence intervals are 95% percentile intervals from 2000
# bootstrap resamples of inputs within strata (fixed seed), so they cover the sampling of
# inputs and of the runs' randomness, which is independent between the two modes (same
# seed, but the two modes draw differently).
#
# Run from the package root, after the sweeps (any environment with QTLSim):
#   julia --project=. validate/fecundity_compare.jl [R=2.0] [CAP=2000]
# R may be a comma-separated list (e.g. 5,2,1.5); a growth rate or group without runs
# under both modes is skipped. Writes the tables to stdout.

using Printf, Random, Statistics

const GROWTH = length(ARGS) >= 1 ? parse.(Float64, split(ARGS[1], ',')) : [2.0]
const CAP = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
const DIR = joinpath(@__DIR__, "demography_sweep")
const GROUPS = (("NoQTL ≤ 1000", L -> L <= 1000), ("NoQTL = 5000", L -> L == 5000))

function readruns(file)
    path = joinpath(DIR, file)
    isfile(path) || return Dict{Int,Any}()
    lines = readlines(path)
    h = Symbol.(split(lines[1], ','))
    conv(v) = v == "true" ? true : v == "false" ? false : (x = tryparse(Int, v); x === nothing ? parse(Float64, v) : x)
    Dict(r.row => r for r in (NamedTuple{Tuple(h)}(Tuple(conv.(split(l, ',')))) for l in lines[2:end]))
end

alive(r) = !r.premove && !r.extinct
died(r) = !r.premove && r.extinct

wquantile(v, w, p) = (o = sortperm(v); c = cumsum(w[o]); v[o[searchsortedfirst(c, p * c[end])]])
wmedian(v, w) = isempty(v) ? NaN : wquantile(v, w, 0.5)

"The statistic `stat(rows, weights)` and its stratified bootstrap 95% interval."
function boot(stat, rows, stratum; B = 2000)
    groups = Dict{Any,Vector{Int}}()
    for r in rows; push!(get!(groups, stratum(r), Int[]), r); end
    w(rs) = (n = Dict(k => count(==(k), stratum.(rs)) for k in keys(groups)); [1 / n[stratum(r)] for r in rs])
    est = stat(rows, w(rows))
    rng = Xoshiro(20261107)
    bs = Float64[]
    for _ in 1:B
        rs = reduce(vcat, [g[rand(rng, 1:length(g), length(g))] for g in values(groups)])
        v = stat(rs, w(rs))
        isfinite(v) && push!(bs, v)
    end
    est, quantile(bs, 0.025), quantile(bs, 0.975)
end
fmt(t; d = 1, pct = false) = pct ? @sprintf("%.*f%% (%.*f–%.*f)", d, 100t[1], d, 100t[2], d, 100t[3]) :
                                   @sprintf("%.*f (%.*f–%.*f)", d, t[1], d, t[2], d, t[3])

"McNemar's exact test on the discordant pairs: a ~ Binomial(a + b, 1/2) if the modes do not differ."
mcnemar(a, b) = (n = a + b; n == 0 ? 1.0 : Float64(min(1.0, 2 * sum(binomial(big(n), k) for k in 0:min(a, b)) / big(2)^n)))

"One demography (`:fixed` or a growth rate) on one group of inputs; returns a summary row for the overview."
function compare(name, V, F, Vfixed, rows, stratum)
    println("### ", name, "\n")
    println("| | viability | fecundity |\n|:--|--:|--:|")
    share(S) = (rs, w) -> sum(w[j] for j in eachindex(rs) if alive(S[rs[j]]); init = 0.0) / sum(w)
    @printf("| died before the move | %d | %d |\n", count(i -> V[i].premove, rows), count(i -> F[i].premove, rows))
    sV, sF = boot(share(V), rows, stratum), boot(share(F), rows, stratum)
    @printf("| alive %d generations after the move | %s | %s |\n", CAP, fmt(sV; pct = true), fmt(sF; pct = true))
    med(S) = (rs, w) -> (k = [j for j in eachindex(rs) if died(S[rs[j]])];
                         wmedian(Float64[S[rs[j]].T for j in k], w[k]))
    @printf("| median generations to extinction (those that died) | %s | %s |\n",
            fmt(boot(med(V), rows, stratum); d = 0), fmt(boot(med(F), rows, stratum); d = 0))
    @printf("| died after the move (of %d) | %d | %d |\n", length(rows), count(i -> died(V[i]), rows), count(i -> died(F[i]), rows))
    med50(S) = (rs, w) -> (k = [j for j in eachindex(rs) if !S[rs[j]].premove]; wmedian(Float64[S[rs[j]].T50 for j in k], w[k]))
    @printf("| median generations to 50%% mortality (W̄ ≈ 0.5 for fecundity) | %s | %s |\n",
            fmt(boot(med50(V), rows, stratum); d = 0), fmt(boot(med50(F), rows, stratum); d = 0))
    medlag(S) = (rs, w) -> (k = [j for j in eachindex(rs) if !S[rs[j]].premove && isfinite(S[rs[j]].lag)];
                            wmedian(Float64[S[rs[j]].lag for j in k], w[k]))
    @printf("| median standardised lag | %s | %s |\n",
            fmt(boot(medlag(V), rows, stratum); d = 3), fmt(boot(medlag(F), rows, stratum); d = 3))
    fixedsize = V === Vfixed
    if fixedsize
        tw = (rs, w) -> sum(w[j] for j in eachindex(rs) if !F[rs[j]].premove && F[rs[j]].Tw < 0; init = 0.0) / sum(w)
        @printf("| N·W̄ never below 2 within %d generations | – | %s |\n", CAP, fmt(boot(tw, rows, stratum); pct = true))
        mtw = (rs, w) -> (k = [j for j in eachindex(rs) if !F[rs[j]].premove && F[rs[j]].Tw >= 0];
                          wmedian(Float64[F[rs[j]].Tw for j in k], w[k]))
        @printf("| median generations until N·W̄ < 2 (those that got there) | – | %s |\n", fmt(boot(mtw, rows, stratum); d = 0))
    end
    println()
    println("Paired by input:\n")
    both = [i for i in rows if died(V[i]) && died(F[i])]
    ratio = (rs, w) -> wmedian([F[i].T / V[i].T for i in rs], w)
    rFV = length(both) >= 10 ? boot(ratio, both, stratum) : (NaN, NaN, NaN)
    length(both) >= 10 &&
        @printf("- time to extinction, fecundity / viability (median over %d inputs where both died after the move): %s\n",
                length(both), fmt(rFV; d = 2))
    rV = rF = (NaN, NaN, NaN); nall = 0
    if fixedsize
        pw = [i for i in rows if died(V[i]) && !F[i].premove && F[i].Tw >= 0]
        twr = (rs, w) -> wmedian([F[i].Tw / V[i].T for i in rs], w)
        length(pw) >= 10 && @printf("- generations until N·W̄ < 2 under fecundity / time to extinction under viability (median over %d inputs): %s\n",
                                    length(pw), fmt(boot(twr, pw, stratum); d = 2))
    else
        all3 = [i for i in rows if died(V[i]) && died(F[i]) && died(Vfixed[i])]
        nall = length(all3)
        if nall >= 10
            rV = boot((rs, w) -> wmedian([V[i].T / Vfixed[i].T for i in rs], w), all3, stratum)
            rF = boot((rs, w) -> wmedian([F[i].T / Vfixed[i].T for i in rs], w), all3, stratum)
            @printf("- time to extinction against QTSim's fixed-size (viability) run of the same input, median over %d inputs where all three died after the move: viability %s, fecundity %s\n",
                    nall, fmt(rV; d = 2), fmt(rF; d = 2))
        end
    end
    a = count(i -> alive(V[i]) && !alive(F[i]), rows); b = count(i -> !alive(V[i]) && alive(F[i]), rows)
    pmc = mcnemar(a, b)
    @printf("- alive under viability but dead under fecundity: %d; the reverse: %d (McNemar exact p = %.3g)\n\n", a, b, pmc)
    (; n = length(rows), sV, sF, rFV, nboth = length(both), rV, rF, nall, a, b, pmc)
end

function main()
    vfixed, ffixed = readruns("runs_fixed.csv"), readruns("runs_fecundity_fixed.csv")
    stratum(i) = (vfixed[i].L, vfixed[i].u)
    println("# Viability against fecundity selection\n")
    println("$CAP generations after the move; strata weighted equally; 95% stratified bootstrap intervals ",
            "(2000 resamples of inputs within strata).\n")
    overview = []
    for (gname, ingroup) in GROUPS
        base = sort([i for i in keys(vfixed) if ingroup(vfixed[i].L)])
        fbase = filter(i -> haskey(ffixed, i), base)
        any(R -> any(i -> haskey(readruns("runs_fecundity_R$(R).csv"), i), base), GROWTH) || isempty(fbase) && continue
        println("## ", gname, "\n")
        if !isempty(fbase)
            strata = sort(unique(stratum.(fbase)))
            println("Fixed size: $(length(fbase)) inputs (",
                    join(["NoQTL $(L), u $(Int(u)): $(count(i -> stratum(i) == (L, u), fbase))" for (L, u) in strata], "; "), ").\n")
            compare("Fixed size (QTSim's demography)", vfixed, ffixed, vfixed, fbase, stratum)
        end
        for R in GROWTH
            V, F = readruns("runs_R$(R).csv"), readruns("runs_fecundity_R$(R).csv")
            rows = sort([i for i in base if haskey(V, i) && haskey(F, i)])
            isempty(rows) && continue
            s = compare("Beverton–Holt, R = $R ($(length(rows)) inputs)", V, F, vfixed, rows, stratum)
            push!(overview, (; gname, R, s...))
        end
    end
    isempty(overview) && return
    println("## Overview\n")
    println("Beverton–Holt runs; time ratios are medians over inputs where the runs compared died after the move.\n")
    println("| group | R | inputs | alive, viability | alive, fecundity | discordant (V only / F only) | T fecundity / viability | T viability / fixed | T fecundity / fixed |")
    println("|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for o in overview
        @printf("| %s | %g | %d | %s | %s | %d / %d | %s | %s | %s |\n", o.gname, o.R, o.n, fmt(o.sV; pct = true),
                fmt(o.sF; pct = true), o.a, o.b, fmt(o.rFV; d = 2), fmt(o.rV; d = 2), fmt(o.rF; d = 2))
    end
end

main()
