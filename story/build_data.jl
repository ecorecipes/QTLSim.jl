# The data of the story (story/index.html): every number and curve it shows, as JSON.
#
#     julia -t 4 --project=vignettes story/build_data.jl
#
# Reads the runs the vignettes store (vignettes/09_replicating_pagel_2025/*.csv,
# vignettes/11_demography_in_the_papers_populations/*.csv), the published-sweep
# comparison (validate/published_sweep_results.csv) and, for the paper's own numbers,
# its results tables in ../Repository/Results (not distributed with this package).
# From the paper's tables only aggregates are written: recomputed medians, binned
# curves, counts; never its rows. Two short QTLSim runs illustrate the model. The
# definitions are vignette 09's and vignette 11's, copied, so the numbers are theirs.
# Writes story/data/story.json.

using QTLSim
using DelimitedFiles, Distributions, Printf, Random, Statistics, XLSX, Dates

const PKG = normpath(joinpath(@__DIR__, ".."))
const V09 = joinpath(PKG, "vignettes", "09_replicating_pagel_2025")
const V11 = joinpath(PKG, "vignettes", "11_demography_in_the_papers_populations")
const RESULTS = normpath(joinpath(PKG, "..", "Repository", "Results"))
const MASTER = joinpath(RESULTS, "Master Output file.xlsx")
const FIGDATA = joinpath(RESULTS, "Figures 1-2 and SI Figures 1-3 data.xlsx")
const OUT = joinpath(@__DIR__, "data", "story.json")
isfile(MASTER) || error("needs the paper's results tables in $RESULTS")

const MOVE_GEN = 10_000
const CAP = 2_000
const POP = 1000
const SAMPLE_SIZE = Dict(100 => 125, 500 => 125, 1000 => 75, 5000 => 30)

# ------------------------------------------------------------------ vignette 09's definitions

function read_master(path)
    dt = XLSX.readtable(path, "Master output file")
    lab = strip.(string.(dt.column_labels))
    num(v) = v isa AbstractString ? parse(Float64, v) : v === missing ? NaN : Float64(v)
    f(name) = num.(dt.data[findfirst(==(name), lab)])
    (; L = round.(Int, f("No. Loci")), u = f("recombination u=3,5"), h2 = f("Heritability"),
       VP = f("PhenotypeVar"), m = f("MutationScalar Vm/Ve"), d = f("d_env"), ω = f("w, fitness function sd"),
       T = f("Generations to extinction"), T50 = f("Generations to 50% mortality"),
       sd50 = f("Phenotype SD at 50% mortality"), x = f("d_env/trait sd at 50% mortality"),
       rate_ext = f("Standardised phenotypic change to extinction"),
       rate50 = f("standardised phenotypic change to 50% mortality"), lag = f("standardised lag"),
       P47 = f("P Survived at generation 10047"), VG = f("Genotype Var at 10000"),
       shc = f("Burger Vg(SHC) eq 21a"), kc = f("standardised kc"),
       adjusted = f("mutation adjusted standardised change"), survived = f("Extinct (0,1)") .== 1)
end

function wquantile(v, w, p)
    o = sortperm(v); c = cumsum(w[o])
    v[o[searchsortedfirst(c, p * c[end])]]
end
"median and quartiles, weighted, of the finite values"
function q3(v, w = ones(length(v)))
    ok = isfinite.(v)
    [wquantile(v[ok], w[ok], p) for p in (0.5, 0.25, 0.75)]
end
med(v, w = ones(length(v))) = q3(v, w)[1]

function logfit(lag, T, w = ones(length(T)))
    ok = (lag .> 0) .& (T .> 0) .& isfinite.(lag)
    x = log.(lag[ok]); y = log.(T[ok]); s = sqrt.(w[ok])
    X = [ones(length(x)) x x .^ 2]
    b = (s .* X) \ (s .* y)
    ȳ = sum(w[ok] .* y) / sum(w[ok])
    (; intercept = b[1], slope = b[2], quadratic = b[3],
       R2 = 1 - sum(w[ok] .* (y - X * b) .^ 2) / sum(w[ok] .* (y .- ȳ) .^ 2), n = count(ok))
end

function critical_rate(; n, VG, ω, VE, sd50, N = POP)
    n >= 3 || return NaN
    No = N * (2n - 3) / (2n * (n - 2))
    Vk = 2No * (1 - 2No / N)
    Ne = 4N / (Vk + 2)
    Vs = ω^2 + VE
    Vλ = Vs + VG + Vs / (2Ne)
    B0 = No * ω / sqrt(Vλ)
    B0 > 1 ? VG / (VG + Vs) * sqrt(2Vλ * log(B0)) / sd50 : NaN
end

function stratified_sample(t; n = SAMPLE_SIZE, seed = 2025)
    rng = Xoshiro(seed)
    rows = Int[]
    for L in (100, 500, 1000, 5000), u in (3.0, 5.0)
        s = findall(i -> t.L[i] == L && t.u[i] == u, eachindex(t.L))
        append!(rows, sort(shuffle(rng, s)[1:n[L]]))
    end
    rows
end

function read_csv(path)
    a, h = readdlm(path, ',', header = true)
    h = vec(String.(h))
    Dict(h[j] => a[:, j] for j in eachindex(h))
end

function read_runs(path)
    c = read_csv(path)
    g(k, T) = T.(c[k])
    (; row = g("row", Int), L = g("L", Int), u = g("u", Float64), T = g("T", Int), extinct = Bool.(c["extinct"]),
       T50 = g("T50", Int), sd50 = g("sd50", Float64), x = g("x", Float64), rate_ext = g("rate_ext", Float64),
       rate50 = haskey(c, "rate50") ? g("rate50", Float64) : fill(NaN, length(c["row"])),
       lag = g("lag", Float64), P47 = haskey(c, "P47") ? g("P47", Float64) : fill(NaN, length(c["row"])),
       VG = g("VG", Float64))
end

function binned(x, y, w = ones(length(x)); nbins = 24)
    o = sortperm(x); c = cumsum(w[o]) ./ sum(w)
    bin = min.(max.(ceil.(Int, c .* nbins), 1), nbins)
    [(wquantile(x[o][bin .== b], w[o][bin .== b], 0.5), wquantile(y[o][bin .== b], w[o][bin .== b], 0.5))
     for b in 1:nbins if any(==(b), bin)]
end
function mean_in_bins(x, y, w; nbins)
    o = sortperm(x); c = cumsum(w[o]) ./ sum(w)
    bin = min.(max.(ceil.(Int, c .* nbins), 1), nbins)
    [(wquantile(x[o][bin .== b], w[o][bin .== b], 0.5),
      sum(w[o][bin .== b] .* y[o][bin .== b]) / sum(w[o][bin .== b])) for b in 1:nbins if any(==(b), bin)]
end
mortality(P47) = replace(1 .- P47 ./ POP, NaN => 1.0)

function logistic(x, y, w)
    X = [ones(length(x)) log.(x)]; b = zeros(2)
    for _ in 1:50
        μ = 1 ./ (1 .+ exp.(-X * b))
        g = X' * (w .* (y .- μ)); H = X' * (X .* (w .* μ .* (1 .- μ)))
        b += H \ g
    end
    b
end

"ratio of geometric means a/b (weights w), with a 95% interval"
function gratio(a, b, w)
    ok = (a .> 0) .& (b .> 0) .& isfinite.(a) .& isfinite.(b)
    δ = log.(a[ok]) .- log.(b[ok]); ω = w[ok] ./ sum(w[ok])
    m = sum(ω .* δ); se = sqrt(sum(ω .^ 2 .* (δ .- m) .^ 2))
    [exp(m), exp(m - 1.96se), exp(m + 1.96se)]
end

# ------------------------------------------------------------------ a small JSON writer

r6(x) = isfinite(x) ? round(x; sigdigits = 6) : nothing
json(io, x::Nothing) = print(io, "null")
json(io, x::Bool) = print(io, x ? "true" : "false")
json(io, x::Integer) = print(io, x)
json(io, x::AbstractFloat) = (v = r6(x); v === nothing ? print(io, "null") : print(io, v == round(v) && abs(v) < 1e15 ? Int(v) : v))
json(io, x::AbstractString) = (print(io, '"'); print(io, escape_string(x)); print(io, '"'))
json(io, x::Symbol) = json(io, String(x))
json(io, x::Tuple) = json(io, collect(x))
function json(io, x::AbstractVector)
    print(io, '[')
    for (i, v) in enumerate(x)
        i > 1 && print(io, ',')
        json(io, v)
    end
    print(io, ']')
end
function json(io, x::Union{AbstractDict,NamedTuple})
    print(io, '{')
    for (i, (k, v)) in enumerate(pairs(x))
        i > 1 && print(io, ',')
        json(io, String(Symbol(k))); print(io, ':'); json(io, v)
    end
    print(io, '}')
end
pts(v) = [[r6(a), r6(b)] for (a, b) in v]

# ------------------------------------------------------------------ the paper and the reruns

t = read_master(MASTER)
rows = stratified_sample(t)
runs = read_runs(joinpath(V09, "replication_runs.csv"))
w = [1 / SAMPLE_SIZE[L] for L in runs.L]
pub(f) = getfield(t, f)[runs.row]
capT(T) = min.(Float64.(T), CAP)
share(alive, w = ones(length(alive))) = 100 * sum(w .* alive) / sum(w)

MGAMMA = fit_mle(Gamma, t.m)
mutation_adjusted(rate, m, w = ones(length(m))) = rate .* pdf.(MGAMMA, m) ./ (sum(w .* pdf.(MGAMMA, m)) / sum(w))
adj = mutation_adjusted(t.rate_ext[t.survived], t.m[t.survived])
kc = filter(isfinite, t.kc[t.survived])

paper = [
    (; quantity = "Generations to extinction", paper = "141 (79–326)", table = q3(t.T)),
    (; quantity = "Generations to 50% mortality", paper = "47 (25–113)", table = q3(t.T50)),
    (; quantity = "Rate of phenotypic change to extinction", paper = "0.056 (0.030–0.103)", table = q3(t.rate_ext)),
    (; quantity = "Rate of phenotypic change to 50% mortality", paper = "0.033 (0.017–0.058)", table = q3(t.rate50)),
    (; quantity = "Rate of phenotypic change, survivors", paper = "0.031 (0.011–0.066)", table = q3(t.rate_ext[t.survived])),
    (; quantity = "The same, adjusted for the mutation scalar", paper = "0.019 (0.006–0.055)", table = q3(adj)),
    (; quantity = "Critical rate k_c, survivors", paper = "0.022 (0.008–0.056)", table = q3(t.kc[t.survived])),
]
fitp = logfit(t.lag, t.T)
kc_pct = [100 * mean(kc .<= r) for r in (0.044, 0.088)]

# design: strata and the spread of the drawn parameters
design = (;
    populations = length(t.L), survived = count(t.survived),
    strata = [(; L, u, n = count(i -> t.L[i] == L && t.u[i] == u, eachindex(t.L)), sample = SAMPLE_SIZE[L])
              for L in (100, 500, 1000, 5000) for u in (3.0, 5.0)],
    parameters = [(; name, q = [r6(quantile(filter(isfinite, v), p)) for p in (0.05, 0.5, 0.95)])
                  for (name, v) in (("phenotypic variance V_P", t.VP), ("heritability h²", t.h2),
                                    ("mutation scalar m", t.m), ("fitness width ω", t.ω),
                                    ("rate of change d", t.d))])

# Figure 1 and survival
fig1a = XLSX.readtable(FIGDATA, "Figure 1a,d")
paper1a = sort(unique(collect(zip(Float64.(fig1a.data[1]), Float64.(fig1a.data[3])))))
fig1bc = XLSX.readtable(FIGDATA, "Figures 1b,c")
paper1bc = [Float64.(c) for c in fig1bc.data]
okp = paper1bc[2] .> 0
pos_p = pub(:lag) .> 0
pos_r = runs.lag .> 0

fig1 = (;
    a = (; paper = pts([(x, min(T, CAP)) for (x, T) in paper1a]),
         published = pts(binned(pub(:x), capT(pub(:T)), w)), reruns = pts(binned(runs.x, capT(runs.T), w))),
    b = (; paper = pts(zip(paper1bc[1][okp], paper1bc[2][okp])),
         published = pts(binned(pub(:x)[pos_p], pub(:lag)[pos_p], w[pos_p])),
         reruns = pts(binned(runs.x[pos_r], runs.lag[pos_r], w[pos_r]))),
    c = (; paper = pts(zip(paper1bc[2][okp], paper1bc[3][okp])),
         published = pts(binned(pub(:lag)[pos_p], mortality(pub(:P47))[pos_p], w[pos_p])),
         reruns = pts(binned(runs.lag[pos_r], mortality(runs.P47)[pos_r], w[pos_r]))),
    fit = (; intercept = fitp.intercept, slope = fitp.slope, quadratic = fitp.quadratic, R2 = fitp.R2))

survival = (;
    all = pts(mean_in_bins(t.x, Float64.(t.T .>= CAP), ones(length(t.x)); nbins = 50)),
    published = pts(mean_in_bins(pub(:x), Float64.(pub(:T) .>= CAP), w; nbins = 24)),
    reruns = pts(mean_in_bins(runs.x, Float64.(.!runs.extinct), w; nbins = 24)))

# genetic variance against the stochastic house of cards (vignette 09)
function house_of_cards(i, ns)
    N, n, μ = POP, t.L[i], 1e-4
    VE = (1 - t.h2[i]) * t.VP[i]
    Vs = t.ω[i]^2 + VE
    α2 = t.m[i] * VE / (2n * μ)
    Ne = 4N / (2N * (ns - 1) / ns^2 + 2)
    4n * μ * Vs / (1 + Vs / (Ne * α2))
end
ns_eq = Float64.(read_csv(joinpath(V09, "replication_runs.csv"))["ns_eq"])
shc_jl = [house_of_cards(i, ns) for (i, ns) in zip(runs.row, ns_eq)]
vg = (; published = med(pub(:VG) ./ pub(:shc), w), reruns = med(runs.VG ./ shc_jl, w))

# the reruns, population by population (QTLSim's own runs)
reruns = [(; L = runs.L[i], u = runs.u[i], x = r6(runs.x[i]), T = runs.T[i], extinct = runs.extinct[i],
             lag = r6(runs.lag[i]), VG = r6(runs.VG[i]), shc = r6(shc_jl[i]), rate = r6(runs.rate_ext[i]))
          for i in eachindex(runs.row)]

agreement = (;
    vg = gratio(runs.VG, pub(:VG), w), time = gratio(capT(runs.T), capT(pub(:T)), w),
    alive_reruns = share(.!runs.extinct, w), alive_published = share(pub(:T) .>= CAP, w),
    fit_reruns = logfit(runs.lag, capT(runs.T), w), fit_published = logfit(pub(:lag), capT(pub(:T)), w))

# survivors, run on to 240,000 generations
long = read_runs(joinpath(V09, "survivor_runs.csv"))
k = indexin(long.row, runs.row)
survived_jl = falses(length(runs.row)); survived_jl[k] .= .!long.extinct
survived_pub = pub(:survived)
sv = findall(.!long.extinct)
rows_jl = long.row[sv]; w_jl = [1 / SAMPLE_SIZE[L] for L in long.L[sv]]
rate_jl = long.rate_ext[sv]
kc_jl = [critical_rate(n = long.P47[i], VG = long.VG[i], ω = t.ω[r], VE = (1 - t.h2[r]) * t.VP[r], sd50 = long.sd50[i])
         for (i, r) in zip(sv, rows_jl)]
rows_pub = runs.row[survived_pub]; w_pub = w[survived_pub]
adj_pub = mutation_adjusted(t.rate_ext[rows_pub], t.m[rows_pub], w_pub)
adj_jl = mutation_adjusted(rate_jl, t.m[rows_jl], w_jl)
survivors = (;
    run_on = length(long.row),
    both = count(survived_pub .& survived_jl), published_only = count(survived_pub .& .!survived_jl),
    reruns_only = count(.!survived_pub .& survived_jl), neither = count(.!survived_pub .& .!survived_jl),
    share_reruns = share(survived_jl, w), share_published = share(survived_pub, w), share_all = share(t.survived),
    stats = [(; quantity = "Rate of phenotypic change", paper = "0.031 (0.011–0.066)",
               all = q3(t.rate_ext[t.survived]), published = q3(t.rate_ext[rows_pub], w_pub), reruns = q3(rate_jl, w_jl)),
             (; quantity = "Adjusted for the mutation scalar", paper = "0.019 (0.006–0.055)",
               all = q3(adj), published = q3(adj_pub, w_pub), reruns = q3(adj_jl, w_jl)),
             (; quantity = "Critical rate k_c", paper = "0.022 (0.008–0.056)",
               all = q3(t.kc[t.survived]), published = q3(t.kc[rows_pub], w_pub), reruns = q3(kc_jl, w_jl))])

summary = [
    (; quantity = "Generations to extinction", paper = "141", table = med(t.T), published = med(capT(pub(:T)), w), reruns = med(capT(runs.T), w)),
    (; quantity = "Generations to 50% mortality", paper = "47", table = med(t.T50), published = med(pub(:T50), w), reruns = med(Float64.(runs.T50), w)),
    (; quantity = "Rate of phenotypic change to extinction", paper = "0.056", table = med(t.rate_ext), published = med(pub(:rate_ext), w), reruns = med(runs.rate_ext, w)),
    (; quantity = "Rate of phenotypic change to 50% mortality", paper = "0.033", table = med(t.rate50), published = med(pub(:rate50), w), reruns = med(runs.rate50, w)),
    (; quantity = "Alive 2000 generations after the move (%)", paper = "", table = share(t.T .>= CAP), published = share(pub(:T) .>= CAP, w), reruns = share(.!runs.extinct, w)),
    (; quantity = "Survived 240,000 generations (%)", paper = "5.1", table = share(t.survived), published = share(survived_pub, w), reruns = share(survived_jl, w)),
    (; quantity = "Survivors' rate of phenotypic change", paper = "0.031", table = med(t.rate_ext[t.survived]), published = med(t.rate_ext[rows_pub], w_pub), reruns = med(rate_jl, w_jl)),
    (; quantity = "Survivors' critical rate k_c", paper = "0.022", table = med(t.kc[t.survived]), published = med(t.kc[rows_pub], w_pub), reruns = med(kc_jl, w_jl)),
]

# ------------------------------------------------------------------ what moves the result

# QTSim's single crossover: the published sweep run with QTSim's crossover model and with
# every crossover applied (validate/published_sweep.jl), against the published runs
ps = read_csv(joinpath(PKG, "validate", "published_sweep_results.csv"))
num(k) = [v isa Number ? Float64(v) : NaN for v in ps[k]]   # empty cells: not run (every crossover at 5000 loci)
wps = ones(length(ps["row"]))
pubT = min.(num("pub_T"), CAP)
crossover = (;
    inputs = length(ps["row"]),
    qtsim = (; vg = gratio(num("jl_qtsim_gv10000"), num("pub_gv10000"), wps), time = gratio(num("jl_qtsim_Tcens"), pubT, wps)),
    all = (; vg = gratio(num("jl_all_gv10000"), num("pub_gv10000"), wps), time = gratio(num("jl_all_Tcens"), pubT, wps)),
    by_loci = [(; L, qtsim = gratio(num("jl_qtsim_Tcens")[ps["NoQTL"] .== L], pubT[ps["NoQTL"] .== L], wps[ps["NoQTL"] .== L])[1],
                all = gratio(num("jl_all_Tcens")[ps["NoQTL"] .== L], pubT[ps["NoQTL"] .== L], wps[ps["NoQTL"] .== L])[1])
               for L in (100, 500, 1000, 5000) if any(==(L), ps["NoQTL"])])

# density regulation (vignette 11): the same 710 populations, Beverton–Holt with R = 5, 2, 1.5
function read_demo(R)
    c = read_csv(joinpath(V11, "demography_runs_R$(R).csv"))
    (; row = Int.(c["row"]), L = Int.(c["L"]), T = Int.(c["T"]), extinct = Bool.(c["extinct"]),
       x = Float64.(c["x"]), rate_ext = Float64.(c["rate_ext"]))
end
weight(L) = 1 / SAMPLE_SIZE[L]
WALL = sum(weight.(t.L[rows]))
sets = [("Fixed size (QTSim)", (; row = runs.row, L = runs.L, T = runs.T, extinct = runs.extinct, x = runs.x, rate_ext = runs.rate_ext));
        [("Regulated, R = $R", read_demo(R)) for R in (5.0, 2.0, 1.5)]]
demography = map(sets) do (name, r)
    ok = isfinite.(r.x) .& (r.x .> 0)
    xs = r.x[ok]; alive = .!r.extinct[ok]; wv = weight.(r.L[ok])
    b = logistic(xs, Float64.(alive), wv)
    edges = exp.(range(log(minimum(xs)), log(maximum(xs)), length = 13))
    bins = [[r6(exp(mean(log.(xs[i])))), r6(sum(wv[i] .* alive[i]) / sum(wv[i]))]
            for i in (findall(v -> edges[b] <= v <= edges[b+1], xs) for b in 1:12) if length(i) >= 5]
    g = exp.(range(log(minimum(xs)), log(maximum(xs)), length = 60))
    jf = indexin(r.row, runs.row)
    both = [i for i in eachindex(r.row) if jf[i] !== nothing && r.extinct[i] && runs.extinct[jf[i]]]
    tratio = isempty(both) ? NaN : med(r.T[both] ./ runs.T[jf[both]], weight.(r.L[both]))
    (; name, alive = 100 * sum(weight.(r.L) .* .!r.extinct) / WALL, x50 = exp(-b[1] / b[2]),
       died_before = length(rows) - length(r.row), time_ratio = tratio, bins,
       curve = [[r6(v), r6(1 / (1 + exp(-(b[1] + b[2] * log(v)))))] for v in g])
end

# ------------------------------------------------------------------ the model, by example

function example(d; seed = 12)
    p = Params(NoGenerations = MOVE_GEN + CAP, Sample = 1, PopSize = POP, Ploidy = 2, NoQTL = 100,
               PhenotypeVar = 0.02, Heritability = 0.6, MutationScalar = 0.002, MutationRatePerGamete = 1e-4,
               RecombPoissonMean = 3.0, InitFitnessSD = 1.0, MoveFitnessSDGen = MOVE_GEN, MoveFitnessSD = d,
               MoveType = :continuous)
    o = simulate(p; seed)
    keep = [i for i in eachindex(o.gen) if (o.gen[i] < MOVE_GEN - 200 ? o.gen[i] % 20 == 0 : o.gen[i] % 2 == 0) || i == lastindex(o.gen)]
    (; d, seed, extinct = o.extinct, last = o.gen[end],
       params = (; N = POP, L = 100, u = 3, VP = 0.02, h2 = 0.6, m = 0.002, mu = 1e-4, omega = 1.0, move = MOVE_GEN),
       gen = o.gen[keep], z = r6.(o.pheno_mean[keep]), theta = r6.(o.theta[keep]),
       vg = r6.(o.geno_var[keep]), vp = r6.(o.pheno_var[keep]), nsurv = o.nsurv[keep])
end
examples = [example(0.008), example(0.015)]

# ------------------------------------------------------------------ write

data = (;
    meta = (; generated = string(Dates.today()), qtlsim = string(pkgversion(QTLSim)),
            sources = ["vignettes/09_replicating_pagel_2025", "vignettes/11_demography_in_the_papers_populations",
                       "validate/published_sweep_results.csv", "the paper's results tables (aggregates only)"],
            move = MOVE_GEN, cap = CAP, N = POP),
    paper, kc_percentiles = kc_pct, design, examples, survival, fig1, reruns, agreement, vg,
    survivors, summary, crossover, demography)

mkpath(dirname(OUT))
open(io -> json(io, data), OUT, "w")
@printf("wrote %s (%.0f KB)\n", OUT, filesize(OUT) / 1024)
