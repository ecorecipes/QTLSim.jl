# Viability against fecundity selection in the published sweep, paired by input and seed.
#
# Reads the caches of validate/demography_sweep.jl, run with SELECTION = viability (QTSim's)
# and fecundity at a fixed size and at one growth rate, and compares the two selection
# modes on the inputs both ran: the share alive CAP generations after the move, the median
# time to extinction, and the paired ratio of times to extinction (fecundity / viability).
# Strata (NoQTL × u) count equally, as in demography_sweep.jl. Confidence intervals are
# 95% percentile intervals from 2000 bootstrap resamples of inputs within strata (fixed
# seed), so they cover the sampling of inputs and of the runs' randomness, which is
# independent between the two modes (same seed, but the two modes draw differently).
#
# Run from the package root, after the sweeps (any environment with QTLSim):
#   julia --project=. validate/fecundity_compare.jl [R=2.0] [CAP=2000]
# Writes the tables to stdout.

using Printf, Random, Statistics

const R = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 2.0
const CAP = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
const DIR = joinpath(@__DIR__, "demography_sweep")

function readruns(file)
    lines = readlines(joinpath(DIR, file))
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

function main()
    sets = Dict((:viability, :fixed) => readruns("runs_fixed.csv"), (:viability, :R) => readruns("runs_R$(R).csv"),
                (:fecundity, :fixed) => readruns("runs_fecundity_fixed.csv"),
                (:fecundity, :R) => readruns("runs_fecundity_R$(R).csv"))
    rows = sort(collect(intersect((keys(s) for s in values(sets))...)))
    stratum(i) = (sets[(:viability, :fixed)][i].L, sets[(:viability, :fixed)][i].u)
    strata = sort(unique(stratum.(rows)))
    println("# Viability against fecundity selection\n")
    println("$(length(rows)) inputs run under both modes, in $(length(strata)) strata ",
            "(", join(["NoQTL $(L), u $(Int(u)): $(count(i -> stratum(i) == (L, u), rows))" for (L, u) in strata], "; "),
            "); $CAP generations after the move; strata weighted equally; 95% stratified bootstrap intervals.\n")
    for (dem, name) in ((:fixed, "fixed size (QTSim)"), (:R, "Beverton–Holt, R = $R"))
        V, F = sets[(:viability, dem)], sets[(:fecundity, dem)]
        println("## ", name, "\n")
        println("| | viability | fecundity |\n|:--|--:|--:|")
        share(S) = (rs, w) -> sum(w[j] for j in eachindex(rs) if alive(S[rs[j]]); init = 0.0) / sum(w)
        @printf("| died before the move | %d | %d |\n", count(i -> V[i].premove, rows), count(i -> F[i].premove, rows))
        @printf("| alive %d generations after the move | %s | %s |\n", CAP,
                fmt(boot(share(V), rows, stratum); pct = true), fmt(boot(share(F), rows, stratum); pct = true))
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
        if dem === :fixed
            tw = (rs, w) -> sum(w[j] for j in eachindex(rs) if !F[rs[j]].premove && F[rs[j]].Tw < 0; init = 0.0) / sum(w)
            @printf("| N·W̄ never below 2 within %d generations | – | %s |\n", CAP, fmt(boot(tw, rows, stratum); pct = true))
            mtw = (rs, w) -> (k = [j for j in eachindex(rs) if !F[rs[j]].premove && F[rs[j]].Tw >= 0];
                              wmedian(Float64[F[rs[j]].Tw for j in k], w[k]))
            @printf("| median generations until N·W̄ < 2 (those that got there) | – | %s |\n", fmt(boot(mtw, rows, stratum); d = 0))
        end
        println()
        # paired
        both = [i for i in rows if died(V[i]) && died(F[i])]
        ratio = (rs, w) -> wmedian([F[i].T / V[i].T for i in rs], w)
        println("Paired by input:\n")
        if length(both) >= 10
            @printf("- time to extinction, fecundity / viability (median over %d inputs where both died after the move): %s\n",
                    length(both), fmt(boot(ratio, both, stratum); d = 2))
        end
        if dem === :fixed
            pw = [i for i in rows if died(V[i]) && !F[i].premove && F[i].Tw >= 0]
            twr = (rs, w) -> wmedian([F[i].Tw / V[i].T for i in rs], w)
            length(pw) >= 10 && @printf("- generations until N·W̄ < 2 under fecundity / time to extinction under viability (median over %d inputs): %s\n",
                                        length(pw), fmt(boot(twr, pw, stratum); d = 2))
        end
        a = count(i -> alive(V[i]) && !alive(F[i]), rows); b = count(i -> !alive(V[i]) && alive(F[i]), rows)
        # McNemar's exact test on the discordant pairs: a ~ Binomial(a + b, 1/2) if the modes do not differ
        n = a + b
        pmc = n == 0 ? 1.0 : min(1.0, 2 * sum(binomial(big(n), k) for k in 0:min(a, b)) / big(2)^n)
        @printf("- alive under viability but dead under fecundity: %d; the reverse: %d (McNemar exact p = %.3g)\n\n", a, b, Float64(pmc))
    end
end

main()
