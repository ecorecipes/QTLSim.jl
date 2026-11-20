# QTSimFast's fast modes (QTSIM_FAST=1, 2) against its exact mode: the comparison
# pre-specified in validate/qtsimfast_fast_criteria.md (committed before this script).
#
#     julia -t 6 --project=vignettes validate/qtsimfast_fast.jl [BINARY]
#
# BINARY defaults to reference/QTSimFast/build/QTSimFast (build it first: make -C
# reference/QTSimFast). Writes validate/qtsimfast_fast.log, the per-run statistics
# validate/qtsimfast_fast_runs.csv and eCDF plots validate/qtsimfast_fast_figs/*.png.
# Exit status 1 if level 1 or level 2 fails. About 3 CPU-hours.

using HypothesisTests, Statistics, Printf, Dates, CairoMakie

const HERE = @__DIR__
const BIN = abspath(get(ARGS, 1, joinpath(HERE, "..", "reference", "QTSimFast", "build", "QTSimFast")))
const EX = joinpath(HERE, "..", "reference", "QTSimFast", "examples")
isfile(BIN) || error("no binary at $BIN")

# ------------------------------------------------------------------ scenarios
base(; L, pois, gens, move = "100000000 0.0") = """
Sample 1
PopSize 1000
NoQTL $L
Generations $gens
Ploidy 2
MoveFitnessSD $move
PhenotypeVar 0.02
MutationScalar 0.001
Heritability 0.6
InitFitnessSD 1.0
OptimumMoveType continuous
MutationRatePerGamete 0.0001
RecombinationPoissonMean $pois
"""
example(name) = join(filter(l -> !startswith(l, "#") && !startswith(lowercase(l), "seed"),
                            readlines(joinpath(EX, name))), "\n") * "\n"

const SCEN = [
    (id = "A", R = 1000, opts = example("single_locus.txt")),
    (id = "B", R = 1000, opts = base(L = 100, pois = 3, gens = 4000)),
    (id = "C", R = 1000, opts = base(L = 1000, pois = 5, gens = 4000)),
    (id = "D", R = 500, opts = base(L = 5000, pois = 3, gens = 2500)),
    (id = "E", R = 1000, opts = example("haploid_0257_loci.txt")),
    (id = "F", R = 1000, opts = base(L = 500, pois = 5, gens = 3000, move = "2000 0.05")),
    (id = "G", R = 1000, opts = base(L = 100, pois = 3, gens = 3000, move = "2000 0.009")),
]
# QTF_SMOKE=1: a check of this script, 24 runs per scenario on other seeds, written to
# *_smoke files (its numbers are not part of the comparison)
const SMOKE = get(ENV, "QTF_SMOKE", "0") == "1"
const SEED_EXACT = SMOKE ? 990_500_001 : 710_000_001
const SEED_FAST2 = SMOKE ? 990_600_001 : 720_000_001
const SFX = SMOKE ? "_smoke" : ""

# ------------------------------------------------------------------ one run
function runqt(opts, seed, fast)
    path, io = mktemp()
    write(io, opts, "Seed $seed\n"); close(io)
    cmd = fast == 0 ? `$BIN $path` : addenv(`$BIN $path`, "QTSIM_FAST" => string(fast))
    out = read(pipeline(cmd; stderr = devnull), String)
    rm(path)
    out
end

# columns: 1 Generation, 2 Phenotype Mean, 3 Phenotype Ratio, 4 Genotype Var, 5 P Survived,
# 6 Fitness Function Mean, 7 Fitness Function SD, 8 Surviving Slope, 9 Std Slope, 10 w
function parse_out(out)
    rows = Dict{Int,Vector{Float64}}()
    for l in eachline(IOBuffer(out))
        f = split(l, '\t')
        length(f) >= 9 || continue
        g = tryparse(Int, f[1]); g === nothing && continue
        v = [something(tryparse(Float64, x), NaN) for x in f[2:9]]
        rows[g] = [Float64(g); v]
    end
    died = occursin("Pop all died", out)
    (; rows, died, last = isempty(rows) ? -1 : maximum(keys(rows)))
end

col(r, g, c) = haskey(r.rows, g) ? r.rows[g][c] : NaN
function wmean(r, lo, hi, c)
    all(g -> haskey(r.rows, g), lo:hi) || return NaN
    mean(r.rows[g][c] for g in lo:hi)
end

# The statistics of one run: name => value (NaN = no value); proportions as 0/1 with a "p:" prefix
function stats(id, r)
    VG, PM, PS, FM, SL = 4, 2, 5, 6, 8
    if id == "A"
        ["a1 VG 4000-4999" => wmean(r, 4000, 4999, VG), "a2 PSurv 4000-4999" => wmean(r, 4000, 4999, PS),
         "a3 PhenMean @5049" => col(r, 5049, PM), "p:a4 alive at end" => r.died ? 0.0 : 1.0,
         "a5 time of death" => r.died ? Float64(r.last) : NaN]
    elseif id in ("B", "C", "D")
        lo, hi = id == "D" ? (2000, 2499) : (3000, 3999)
        k = lowercase(id)
        ["$(k)1 VG window" => wmean(r, lo, hi, VG), "$(k)2 VG @last" => col(r, hi, VG),
         "$(k)3 PhenMean window" => wmean(r, lo, hi, PM), "$(k)4 PSurv window" => wmean(r, lo, hi, PS),
         "$(k)5 SurvSlope window" => wmean(r, lo, hi, SL)]
    elseif id == "E"
        ["e1 VG 200-299" => wmean(r, 200, 299, VG), "e2 PSurv 200-299" => wmean(r, 200, 299, PS),
         "p:e3 alive at end" => r.died ? 0.0 : 1.0, "e4 PhenMean @last" => r.died ? NaN : col(r, 599, PM)]
    elseif id == "F"
        ["f1 VG 1500-1999" => wmean(r, 1500, 1999, VG), "f2 time of death" => r.died ? Float64(r.last) : NaN,
         "p:f3 alive at end" => r.died ? 0.0 : 1.0]
    else
        ["g1 VG 1500-1999" => wmean(r, 1500, 1999, VG), "p:g2 alive at end" => r.died ? 0.0 : 1.0,
         "g3 time of death" => r.died ? Float64(r.last) : NaN,
         "g4 lag @2199" => col(r, 2199, FM) - col(r, 2199, PM)]
    end
end

function firstdiff(a, b)
    la, lb = split(a, '\n'), split(b, '\n')
    for i in 1:min(length(la), length(lb))
        la[i] == lb[i] || return (f = split(la[i], '\t'); tryparse(Int, f[1]) === nothing ? "line $i" : "generation $(f[1])")
    end
    "length"
end

# ------------------------------------------------------------------ runs
logio = open(joinpath(HERE, "qtsimfast_fast$(SFX).log"), "w")
say(x...) = (println(x...); println(logio, x...); flush(logio))
say("QTSimFast fast modes against the exact mode (criteria: validate/qtsimfast_fast_criteria.md)")
say("started ", now(), "; binary ", BIN, "; ", Threads.nthreads(), " threads")

results = Dict{String,Any}()
csv = open(joinpath(HERE, "qtsimfast_fast_runs$(SFX).csv"), "w")
println(csv, "scenario,program,seed,statistic,value")
csvlock = ReentrantLock()
t0 = time()
for s in SCEN
    R = SMOKE ? 24 : s.R
    ex = Vector{Any}(undef, R); f2 = Vector{Any}(undef, R)
    ident = falses(R); diffs = Vector{String}(undef, R)
    Threads.@threads :greedy for i in 1:R
        seed = SEED_EXACT + i - 1
        oe = runqt(s.opts, seed, 0)
        o1 = runqt(s.opts, seed, 1)
        ident[i] = oe == o1
        diffs[i] = ident[i] ? "" : firstdiff(oe, o1)
        ex[i] = stats(s.id, parse_out(oe))
        seed2 = SEED_FAST2 + i - 1
        f2[i] = stats(s.id, parse_out(runqt(s.opts, seed2, 2)))
        lock(csvlock) do
            for (k, v) in ex[i]; println(csv, s.id, ",exact,", seed, ",", k, ",", v); end
            for (k, v) in f2[i]; println(csv, s.id, ",fast2,", seed2, ",", k, ",", v); end
        end
    end
    results[s.id] = (; ex, f2, ident, diffs, R)
    say(@sprintf("scenario %s: %d runs per program done (%.0f s elapsed)", s.id, R, time() - t0))
end
close(csv)

# ------------------------------------------------------------------ level 1
say("\n== Level 1 (QTSIM_FAST=1, the exact runs' seeds): byte-identical outputs")
l1pass = true
for s in SCEN
    r = results[s.id]; n = count(r.ident)
    frac = n / r.R
    ok = frac >= 0.95
    global l1pass &= ok
    d = [r.diffs[i] for i in 1:r.R if !r.ident[i]]
    say(@sprintf("  %s: %d of %d identical (%.1f%%)%s", s.id, n, r.R, 100frac, ok ? "" : "  BELOW 95%"),
        isempty(d) ? "" : "; first differences at: " * join(first(d, 10), ", ") * (length(d) > 10 ? ", …" : ""))
end
say("Level 1: ", l1pass ? "PASS" : "FAIL")

# ------------------------------------------------------------------ level 2
say("\n== Level 2 (QTSIM_FAST=2, own seeds) against the exact runs")
tests = []   # (scenario, name, p, D or NaN, note)
figdir = joinpath(HERE, "qtsimfast_fast_figs$(SFX)"); mkpath(figdir)
for s in SCEN
    r = results[s.id]
    names = first.(r.ex[1])
    fig = Figure(size = (300 * length(names), 280), fontsize = 11)
    for (j, nm) in enumerate(names)
        a = filter(!isnan, [r.ex[i][j][2] for i in 1:r.R])
        b = filter(!isnan, [r.f2[i][j][2] for i in 1:r.R])
        if startswith(nm, "p:")
            ka, kb = count(==(1.0), a), count(==(1.0), b)
            p = pvalue(FisherExactTest(ka, length(a) - ka, kb, length(b) - kb))
            push!(tests, (s.id, nm[3:end], p, NaN, @sprintf("exact %d/%d, fast2 %d/%d", ka, length(a), kb, length(b))))
            continue
        end
        if length(a) < 20 || length(b) < 20
            push!(tests, (s.id, nm, NaN, NaN, @sprintf("NOT TESTED: %d and %d values", length(a), length(b))))
            continue
        end
        ks = ApproximateTwoSampleKSTest(a, b)
        D = ks.δ
        ad = try pvalue(KSampleADTest(a, b)) catch; NaN end
        ma, mb = mean(a), mean(b)
        se = sqrt(var(a) / length(a) + var(b) / length(b))
        rel = ma != 0 ? @sprintf("mean diff %+.4f%% [%+.4f, %+.4f]", 100(mb - ma) / abs(ma), 100(mb - ma - 1.96se) / abs(ma), 100(mb - ma + 1.96se) / abs(ma)) : ""
        push!(tests, (s.id, nm, pvalue(ks), D, @sprintf("n %d/%d, D %.4f, AD p %.3f, %s", length(a), length(b), D, ad, rel)))
        ax = Axis(fig[1, j], title = "$(s.id): $nm", titlesize = 10)
        sa, sb = sort(a), sort(b)
        stairs!(ax, sa, (1:length(sa)) ./ length(sa), color = :black, label = "exact")
        stairs!(ax, sb, (1:length(sb)) ./ length(sb), color = :orangered, linestyle = :dash, label = "fast 2")
        j == 1 && axislegend(ax, position = :rb, framevisible = false, labelsize = 9)
    end
    save(joinpath(figdir, "scenario_$(s.id).png"), fig)
end

# Holm over the tests performed
done = [t for t in tests if !isnan(t[3])]
m = length(done)
order = sortperm([t[3] for t in done])
rejected = falses(m)
for (rank, k) in enumerate(order)
    done[k][3] <= 0.05 / (m - rank + 1) || break
    rejected[k] = true
end
for t in tests
    k = findfirst(x -> x === t, done)
    verdict = isnan(t[3]) ? "" : (rejected[k] ? "  REJECTED (Holm)" : "")
    say(@sprintf("  %s %-24s p %s  %s%s", t[1], t[2], isnan(t[3]) ? "  –   " : @sprintf("%.4f", t[3]), t[5], verdict))
end
nt = count(t -> isnan(t[3]), tests)
l2pass = !any(rejected)
say(@sprintf("\n%d primary tests performed (%d not testable); smallest p %.4f; Holm rejections: %d", m, nt, minimum(t[3] for t in done), count(rejected)))
say("Level 2: ", l2pass ? "PASS" : "FAIL: " * join([done[k][1] * " " * done[k][2] for k in 1:m if rejected[k]], "; "))
say("\nfinished ", now(), @sprintf(" (%.0f s)", time() - t0))
close(logio)
exit(l1pass && l2pass ? 0 : 1)
