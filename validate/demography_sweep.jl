# Rerun the paper's published sweep under Beverton–Holt density regulation.
#
# Vignette 11 reran a weighted stratified sample of 710 of the paper's 119,984 inputs
# with `Demography = :beverton_holt` and growth rates R = 5, 2 and 1.5, for 2000
# generations after the optimum starts to move. This script runs a larger sample, up
# to the whole sweep, and runs each input at QTSim's fixed size too, so that every
# comparison is paired by input.
#
# The sample takes a fraction SCALE of each of the 8 strata (NoQTL × u), after shuffling
# each stratum once with a fixed generator, so a smaller sample is the start of a larger
# one and the caches of a small run are reused by a larger one. Strata count equally in
# every summary (weight 1/n per run, n the stratum's sample size).
#
# Seeds: run of table row i has seed SEED_BASE + i, for every demography (common random
# numbers within this replication). Vignettes 09 and 11 used seed = i (i ≤ 119,984), and
# published_sweep.jl seeds of at most 6.1e8, so these runs are independent of theirs.
#
# Parameters come from the published table (`Repository/Results/Master Output file.xlsx`),
# as in vignettes 09 and 11. Each run is cached as one row of
# validate/demography_sweep/runs_<demography>.csv as soon as it finishes, so an
# interrupted sweep resumes. Populations that die before the move are recorded
# (premove = true), not dropped; they count as dead in the share alive.
#
# Run from the package root, in the vignettes environment (XLSX):
#   julia -t 10 --project=vignettes validate/demography_sweep.jl [R=5,2,1.5] [SCALE=0.01] [CAP=2000] [MAXL]
# R is a comma-separated list of growth rates; SCALE the fraction of each stratum
# (1.0 = all 119,984 inputs); CAP the generations after the move; MAXL skips strata with
# more loci. The summary is printed and written to validate/demography_sweep/summary.md,
# with the mean seconds per run and the CPU-hours they imply at scales 0.01, 0.05 and 1.
#
# Cost (2026-11-07, 12-core Apple M3 Pro, Julia 1.12; validate/benchmark.jl): most runs die
# soon after the move, so a run costs about 10,000 generations, at 0.25, 0.33, 0.39 and
# 0.78 ms per generation on one thread for 100, 500, 1000 and 5000 loci. With the fixed
# size and three growth rates (four runs per input) that is about 6 CPU-hours at scale
# 0.01, 30 at 0.05 and 600 for the whole sweep, 45% of it in the 5000-locus strata
# (MAXL = 1000 cuts it to 330). Several threads take up to about twice the CPU time of
# one, so scale 0.01 takes roughly an hour on 10 threads (estimated).

using QTLSim, XLSX, Printf, Random, Statistics

const GROWTH = length(ARGS) >= 1 ? Tuple(parse.(Float64, split(ARGS[1], ','))) : (5.0, 2.0, 1.5)
const SCALE  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.01
const CAP    = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2000
const MAXL   = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : typemax(Int)
const MOVE_GEN = 10_000
const POP = 1000
const SEED_BASE = 2_026_110_700_000           # 2026-11-07, then five digits for the row
const SAMPLE_SEED = 20261107
const STRATA = [(L, u) for L in (100, 500, 1000, 5000) for u in (3.0, 5.0) if L <= MAXL]
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MASTER = joinpath(ROOT, "Repository", "Results", "Master Output file.xlsx")
const OUT = joinpath(@__DIR__, "demography_sweep")

function read_master(path)
    dt = XLSX.readtable(path, "Master output file")
    lab = strip.(string.(dt.column_labels))
    num(v) = v isa AbstractString ? parse(Float64, v) : v === missing ? NaN : Float64(v)
    f(name) = num.(dt.data[findfirst(==(name), lab)])
    (; L = round.(Int, f("No. Loci")), u = f("recombination u=3,5"), h2 = f("Heritability"),
       VP = f("PhenotypeVar"), m = f("MutationScalar Vm/Ve"), d = f("d_env"), ω = f("w, fitness function sd"),
       T = round.(Int, f("Generations to extinction")))
end

"A fraction `scale` of each stratum; the sample of a smaller scale is the start of a larger one."
function nested_sample(t, scale)
    rng = Xoshiro(SAMPLE_SEED)
    rows = Int[]
    for (L, u) in [(L, u) for L in (100, 500, 1000, 5000) for u in (3.0, 5.0)]   # all strata, so MAXL keeps the draws
        s = shuffle(rng, findall(i -> t.L[i] == L && t.u[i] == u, eachindex(t.L)))
        L <= MAXL && append!(rows, sort(s[1:clamp(round(Int, scale * length(s)), 1, length(s))]))
    end
    rows
end

# Vignette 11's measures: 50% mortality is the first generation after the move at which
# the fraction surviving selection is closest to one half. N_eq is the mean population
# size over the 100 generations before the move, N50 the size at 50% mortality.
function measures(o, d; movegen = MOVE_GEN)
    gen, zbar, sd, ns, N = o.gen, o.pheno_mean, sqrt.(o.pheno_var), o.nsurv, o.popsize
    i0 = findfirst(==(movegen - 1), gen)
    post = findall(>=(movegen), gen)
    k = isempty(post) ? lastindex(gen) : post[argmin(abs.(ns[post] ./ N[post] .- 0.5))]
    T = o.extinct ? gen[end] + 2 - movegen : length(post)
    T50 = gen[k] + 2 - movegen
    sd50 = sd[k]
    x = d / sd50
    rate_ext = (zbar[end] - zbar[i0]) / T / sd50
    (; T, extinct = o.extinct, T50, sd50, x, rate_ext, lag = x - rate_ext, VG = o.geno_var[i0],
       N_eq = mean(N[i0-99:i0]), N50 = Float64(N[k]))
end

function run_one(t, i, R)
    dem = R === nothing ? (; Demography = :fixed) : (; Demography = :beverton_holt, GrowthRate = R)
    p = Params(; NoGenerations = MOVE_GEN + CAP, Sample = 1, PopSize = POP, Ploidy = 2, NoQTL = t.L[i],
               PhenotypeVar = t.VP[i], Heritability = t.h2[i], MutationScalar = t.m[i],
               MutationRatePerGamete = 1e-4, RecombPoissonMean = t.u[i], InitFitnessSD = t.ω[i],
               MoveFitnessSDGen = MOVE_GEN, MoveFitnessSD = t.d[i], MoveType = :continuous, dem...)
    t0 = time()
    o = simulate(p; seed = SEED_BASE + i)
    s = time() - t0
    head = (; row = i, L = t.L[i], u = t.u[i])
    if findfirst(==(MOVE_GEN - 1), o.gen) === nothing             # died before the move
        return (; head..., premove = true, died_gen = o.extinct_gen, T = -1, extinct = true, T50 = -1,
                sd50 = NaN, x = NaN, rate_ext = NaN, lag = NaN, VG = NaN, N_eq = NaN, N50 = NaN, seconds = s)
    end
    (; head..., premove = false, died_gen = o.extinct_gen, measures(o, t.d[i])..., seconds = s)
end

# --- summaries --------------------------------------------------------------
function wquantile(v, w, p)
    o = sortperm(v); c = cumsum(w[o])
    v[o[searchsortedfirst(c, p * c[end])]]
end
wmed(v, w) = (ok = isfinite.(v); any(ok) ? wquantile(v[ok], w[ok], 0.5) : NaN)

"Weighted logistic regression of `y` on `[1 log x]`, by Newton's method."
function logistic(x, y, w)
    X = [ones(length(x)) log.(x)]; b = zeros(2)
    for _ in 1:50
        μ = 1 ./ (1 .+ exp.(-X * b))
        g = X' * (w .* (y .- μ)); H = X' * (X .* (w .* μ .* (1 .- μ)))
        b += H \ g
    end
    b
end

function summarise(io, t, rows, sets)
    nstr = Dict(s => count(i -> (t.L[i], t.u[i]) == s, rows) for s in STRATA)
    wt(r) = 1 / nstr[(r.L, r.u)]
    W = sum(wt, sets[1][2])
    alive(r) = !r.premove && !r.extinct
    tab(names) = (println(io, "| | ", join(names, " | "), " |"); println(io, "|:--|", repeat("--:|", length(names))))
    line(name, vals) = println(io, "| ", name, " | ", join(vals, " | "), " |")
    col(f) = [f(rs) for (_, rs) in sets]
    post(rs) = filter(r -> !r.premove, rs)
    wm(rs, f) = (p = post(rs); @sprintf("%.0f", wmed(Float64[f(r) for r in p], wt.(p))))

    println(io, "## Each demography\n")
    pub = [min(t.T[r.row], CAP) for r in sets[1][2]]
    tab([first.(sets); "published (QTSim)"])
    line("population size before the move (median)", [col(rs -> wm(rs, r -> r.N_eq)); "$POP"])
    line("died before the move (of $(length(rows)))", [col(rs -> string(count(r -> r.premove, rs))); "0"])
    line("alive $CAP generations after the move",
         [col(rs -> @sprintf("%.1f%%", 100 * sum(wt(r) for r in rs if alive(r); init = 0.0) / W));
          @sprintf("%.1f%%", 100 * sum(wt(r) for r in sets[1][2] if t.T[r.row] > CAP; init = 0.0) / W)])
    line("generations to extinction (median)",
         [col(rs -> wm(rs, r -> r.T)); @sprintf("%.0f", wmed(Float64.(pub), wt.(sets[1][2])))])
    line("generations to 50% mortality (median)", [col(rs -> wm(rs, r -> r.T50)); "–"])
    f3(rs, f) = (p = post(rs); @sprintf("%.3f", wmed(Float64[f(r) for r in p], wt.(p))))
    line("rate of change to extinction (median)", [col(rs -> f3(rs, r -> r.rate_ext)); "–"])
    line("standardised lag (median)", [col(rs -> f3(rs, r -> r.lag)); "–"])
    line("rate of change, populations alive at $CAP (median)",
         [col(rs -> (a = filter(alive, rs); isempty(a) ? "–" : f3(a, r -> r.rate_ext))); "–"])
    xs = col(rs -> begin
        p = filter(r -> !r.premove && isfinite(r.x) && r.x > 0, rs)
        length(p) < 10 || all(alive, p) || !any(alive, p) ? "–" :
            (b = logistic([r.x for r in p], Float64[alive(r) for r in p], wt.(p)); @sprintf("%.4f", exp(-b[1] / b[2])))
    end)
    line("standardised rate at which half are alive (logistic x50)", [xs; "–"])

    println(io, "\n## Paired by input, against the same input at a fixed size\n")
    fixed = Dict(r.row => r for r in sets[1][2])
    dsets = sets[2:end]
    tab(first.(dsets))
    pairs(rs) = [(r, fixed[r.row]) for r in rs if haskey(fixed, r.row)]
    line("time to extinction / fixed (median, both died after the move)", [begin
        b = [(r, f) for (r, f) in pairs(rs) if !r.premove && !f.premove && r.extinct && f.extinct]
        isempty(b) ? "–" : @sprintf("%.2f", wmed([r.T / f.T for (r, f) in b], [wt(r) for (r, _) in b]))
    end for (_, rs) in dsets])
    line("alive at a fixed size, dead here", [string(count(((r, f),) -> alive(f) && !alive(r), pairs(rs))) for (_, rs) in dsets])
    line("dead at a fixed size, alive here", [string(count(((r, f),) -> !alive(f) && alive(r), pairs(rs))) for (_, rs) in dsets])
    println(io, "\n(The fixed-size column `published (QTSim)` uses the table's own runs of the same inputs, censored at $CAP.)")
end

function cost_report(io, t, sets)
    println(io, "\n## Cost\n")
    println(io, "Mean seconds per run (on this machine, with the threads it ran on), and the CPU-hours")
    println(io, "they imply for every demography run here (fixed size and each R) at three scales:\n")
    tab = Dict{Tuple{Int,Float64},Float64}()
    for (L, u) in STRATA
        tot = 0.0
        for (_, rs) in sets
            s = [r.seconds for r in rs if (r.L, r.u) == (L, u)]
            tot += isempty(s) ? NaN : mean(s)
        end
        tab[(L, u)] = tot
        @printf(io, "- NoQTL %4d, u %.0f: %.1f s per input over %d demographies\n", L, u, tot, length(sets))
    end
    for sc in (0.01, 0.05, 1.0)
        h = sum(tab[(L, u)] * round(Int, sc * count(i -> t.L[i] == L && t.u[i] == u, eachindex(t.L))) for (L, u) in STRATA) / 3600
        @printf(io, "- scale %.2f: about %.0f CPU-hours\n", sc, h)
    end
end

function main()
    t0 = time()
    isfile(MASTER) || error("this script needs the paper's results table at $MASTER")
    t = read_master(MASTER)
    rows = nested_sample(t, SCALE)
    mkpath(OUT)
    @printf("%d inputs (scale %g of each stratum%s), %d generations after the move, R = %s, %d threads\n",
            length(rows), SCALE, MAXL < typemax(Int) ? ", NoQTL ≤ $MAXL" : "", CAP, join(GROWTH, ", "), Threads.nthreads())
    sets = Pair{String,Vector}[]
    for R in (nothing, GROWTH...)
        name = R === nothing ? "fixed (QTSim)" : "R = $R"
        cache = joinpath(OUT, R === nothing ? "runs_fixed.csv" : "runs_R$(R).csv")
        @printf("%s: cache %s (%.0f min so far)\n", name, relpath(cache, pwd()), (time() - t0) / 60)
        rs = sweep(i -> run_one(t, i, R), rows; cache, cost = i -> t.L[i],
                   progress = max(1, length(rows) ÷ 20))
        push!(sets, name => rs)
    end
    io = IOBuffer()
    println(io, "# Density regulation in the published sweep\n")
    println(io, "$(length(rows)) inputs, scale $SCALE of each stratum, $CAP generations after the move, ",
            "seeds $SEED_BASE + row. Strata weighted equally.\n")
    summarise(io, t, rows, [first(s) => last(s) for s in sets])
    cost_report(io, t, sets)
    s = String(take!(io))
    print(s)
    write(joinpath(OUT, "summary.md"), s)
    @printf("\ntotal %.1f min\n", (time() - t0) / 60)
end

main()
