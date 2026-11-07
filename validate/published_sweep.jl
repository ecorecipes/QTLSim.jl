# Reproduce the paper's published simulation sweep with the compositional model.
#
# Ground truth is `Repository/Results/Master Output file.xlsx`: one row per shipped
# input file in `Repository/Simulations/QT<NoQTL>-Possion<u>/` (8 strata × 14998
# QTSim runs: NoQTL ∈ {100,500,1000,5000}, RecombinationPoissonMean u ∈ {3,5},
# PopSize 1000, diploid, optimum moving continuously from generation 10000,
# 250000 generations; other parameters drawn per file). Rows are matched to files
# by their parameters. In the table, "Generations to extinction" counts from
# the move start (240000 = survived) and "Extinct (0,1)" is 1 for SURVIVORS.
#
# Three checks:
#   1. ANCHOR: the local QTSim binary, run on shipped files, reproduces the
#      published per-run numbers exactly (pre-move state, extinction time).
#   2. EQUILIBRIUM: the Julia kernel (its own RNG) on a random sample of inputs per
#      stratum. Genotypic variance and survivor count at generation 10000 are
#      compared with the published values for the same inputs (paired tests).
#   3. EXTINCTION: time to extinction after the move starts, censored at CAP
#      generations, paired by input, plus binned medians against the standardised
#      rate of environmental change d_env/√V_P (the paper's headline relationship).
# Checks 2–3 run for CrossoverModel = :qtsim (QTSim's single-crossover gametes,
# the default) and :all (every Poisson crossover point applied).
#
# Run from the package root (needs the vignettes env for XLSX + HypothesisTests, and the QTSim binary):
#   julia -t 10 --project=vignettes validate/published_sweep.jl [CAP=2000] [SCALE=1.0] [MAXL]
# SCALE multiplies the per-stratum sample sizes; MAXL skips strata with more loci. Per-run results are written to
# validate/published_sweep_results.csv.
#
# Cost (2026-11-07, 12-core Apple M3 Pro, Julia 1.12): the default sample (589 inputs, 1154
# kernel runs of up to 12,000 generations) is about 1 CPU-hour at the single-thread speeds
# of validate/benchmark.jl; a 4-thread run of `2000 0.1` (58 inputs) took 4 minutes with
# other jobs on the machine, so the whole sample takes roughly half an hour on 4 threads
# (extrapolated). A 5000-locus World holds about 20 MB (dense genomes took 160 MB). The
# last full run, in published_sweep.log, took 146 minutes on 8 threads with the dense
# kernel.

include(joinpath(@__DIR__, "gsl_locked.jl"))     # QTSIM, SIMDIR, ROOT, read_qtsim_opts
using XLSX, HypothesisTests, Statistics, Random, Printf

const CAP   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
const SCALE = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.0
const MAXL  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : typemax(Int)   # skip strata with more loci
const MASTER = joinpath(ROOT, "Repository", "Results", "Master Output file.xlsx")
const MOVE_GEN = 10000
# The table numbers generations from 1: its "at 10000" columns are QTSim's output
# line for (0-based) generation 9999, the last generation before the optimum moves.
const REF_GEN = MOVE_GEN - 1
# (stratum dir, NoQTL, u, sample size, crossover models to run)
const STRATA = [("QT0100-Possion3.0", 100, 3, 150, (:qtsim, :all)), ("QT0100-Possion5.0", 100, 5, 150, (:qtsim, :all)),
                ("QT0500-Possion3.0", 500, 3,  90, (:qtsim, :all)), ("QT0500-Possion5.0", 500, 5,  90, (:qtsim, :all)),
                ("QT1000-Possion3.0", 1000, 3, 45, (:qtsim, :all)), ("QT1000-Possion5.0", 1000, 5, 45, (:qtsim, :all)),
                ("QT5000-Possion3.0", 5000, 3, 12, (:qtsim,)),      ("QT5000-Possion5.0", 5000, 5, 12, (:qtsim,))]

# --- published table -----------------------------------------------------
function load_master()
    dt = XLSX.readtable(MASTER, "Master output file")
    lab = string.(dt.column_labels)
    c(name) = dt.data[findfirst(==(name), lab)]
    num(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)     # some cells are stored as text
    I(name) = round.(Int, num.(c(name))); F(name) = num.(c(name))
    (L = I("No. Loci"), u = I("recombination u=3,5"), H = F("Heritability"),
     PV = F("PhenotypeVar"), MS = F("MutationScalar Vm/Ve"), d = F("d_env"),
     T = I("Generations to extinction"), surv = I("Extinct (0,1)"),
     pm = F("Phenotype Mean at 10000"), gv = F("Genotype Var at 10000"), ns = I("P Survived at 10000"))
end

function match_row(tab, L, u, q::Params)
    hits = findall(i -> tab.L[i] == L && tab.u[i] == u &&
                        abs(tab.H[i] - q.Heritability) < 1e-9 && abs(tab.PV[i] - q.PhenotypeVar) < 1e-9 &&
                        # some rows round MutationScalar or d_env to 3 significant figures
                        abs(tab.d[i] - q.MoveFitnessSD) <= 5e-3 * q.MoveFitnessSD &&
                        abs(tab.MS[i] - q.MutationScalar) <= 5e-3 * q.MutationScalar,
                   eachindex(tab.L))
    length(hits) == 1 ? hits[1] : 0
end

# --- QTSim run on an unmodified shipped file -------------------------------
function qtsim_shipped(path)
    g = Int[]; pm = Float64[]; gv = Float64[]; ns = Int[]; b = false; ext = false
    for ln in eachline(`$QTSIM $path`)
        startswith(ln, "Pop all died") && (ext = true; break)
        startswith(ln, "Generation\t") && (b = true; continue)
        b || continue
        f = split(ln, '\t'); length(f) < 5 && continue
        push!(g, parse(Int, f[1])); push!(pm, parse(Float64, f[2])); push!(gv, parse(Float64, f[4])); push!(ns, parse(Int, f[5]))
    end
    (gen=g, pm=pm, gv=gv, ns=ns, extinct=ext, tdie = ext ? last(g) + 1 : -1)   # shipped files use Sample 1
end

# --- one Julia run ---------------------------------------------------------
function julia_run(q::Params, model::Symbol, seed::Int)
    kw = Dict{Symbol,Any}(k => getfield(q, k) for k in fieldnames(Params))
    kw[:NoGenerations] = MOVE_GEN + CAP; kw[:Sample] = REF_GEN; kw[:CrossoverModel] = model
    o = simulate(Params(; kw...); seed=seed)
    i = findfirst(==(REF_GEN), o.gen)
    (gv = i === nothing ? NaN : o.geno_var[i], ns = i === nothing ? -1 : o.nsurv[i],
     pm = i === nothing ? NaN : o.pheno_mean[i], tdie = o.extinct ? o.extinct_gen : -1)
end

function sweep_main()
    isfile(QTSIM) || error(qtsim_missing())
    t0 = time()
    tab = load_master()
    @printf("published table: %d runs loaded (%.0f s); CAP=%d generations after the move; SCALE=%.2f; %d threads\n",
            length(tab.L), time() - t0, CAP, SCALE, Threads.nthreads())

    # sample files per stratum and match them to table rows
    rng = MersenneTwister(20261104)
    jobs = NamedTuple[]
    for (dir, L, u, n, models) in filter(s -> s[2] <= MAXL, STRATA)
        files = sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, dir))))
        pick = sort(shuffle(rng, files)[1:max(1, round(Int, n * SCALE))])
        for f in pick
            q, _ = read_qtsim_opts(joinpath(SIMDIR, dir, f))
            (q.NoQTL == L && q.RecombPoissonMean == u) || error("$dir/$f does not match its stratum")
            r = match_row(tab, L, u, q)
            r == 0 && (println("  no unique table row for $dir/$f; skipped"); continue)
            push!(jobs, (dir=dir, file=f, L=L, u=u, q=q, row=r, models=models))
        end
    end
    println("matched $(length(jobs)) sampled input files to table rows")

    # 1. ANCHOR: QTSim on shipped files that died quickly, vs the published row
    println("\n── 1. anchor: local QTSim on unmodified shipped files vs the published row ──")
    anchors = NamedTuple[]                    # first two quick-dying runs in three strata
    for (dir, L, u) in filter(s -> s[2] <= MAXL, [("QT0100-Possion3.0", 100, 3), ("QT0500-Possion5.0", 500, 5), ("QT1000-Possion3.0", 1000, 3)])
        nfound = 0
        for f in sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, dir))))
            q, _ = read_qtsim_opts(joinpath(SIMDIR, dir, f)); r = match_row(tab, L, u, q)
            (r > 0 && tab.surv[r] == 0 && tab.T[r] < 300) || continue
            push!(anchors, (dir=dir, file=f, L=L, u=u, q=q, row=r)); nfound += 1
            nfound == 2 && break
        end
    end
    ares = Vector{Any}(undef, length(anchors))
    Threads.@threads :dynamic for k in eachindex(anchors)
        ares[k] = qtsim_shipped(joinpath(SIMDIR, anchors[k].dir, anchors[k].file))
    end
    offsets = Int[]; anchor_ok = true
    for (a, o) in zip(anchors, ares)
        i = findfirst(==(REF_GEN), o.gen)
        ok = i !== nothing && f6(o.pm[i]) == f6(tab.pm[a.row]) && f6(o.gv[i]) == f6(tab.gv[a.row]) && o.ns[i] == tab.ns[a.row]
        off = tab.T[a.row] - (o.tdie - MOVE_GEN); push!(offsets, off)
        anchor_ok &= ok
        @printf("  %s/%s row %-6d pre-move (mean, V_G, survivors) %s; dies at gen %d, table says %d after the move (offset %+d)\n",
                a.dir, a.file, a.row, ok ? "IDENTICAL" : "DIFFERENT", o.tdie, tab.T[a.row], off)
    end
    off = length(unique(offsets)) == 1 ? offsets[1] : error("inconsistent extinction-time offsets $offsets")
    println(anchor_ok ? "  ⇒ the local QTSim reproduces the published per-run values exactly; table time = (dying generation − $MOVE_GEN) + $off" :
                        "  ⇒ MISMATCH: the local QTSim does not reproduce the published rows")

    # 2–3. Julia kernel on the sample, both crossover models
    tasks = [(k, m) for (k, j) in enumerate(jobs) for m in j.models]
    out = Dict{Tuple{Int,Symbol},Any}(); lk = ReentrantLock(); done = Threads.Atomic{Int}(0)
    println("\nrunning $(length(tasks)) Julia simulations …")
    Threads.@threads :dynamic for t in tasks
        k, m = t
        res = julia_run(jobs[k].q, m, 1_000_003 * k + (m === :all ? 7 : 0))
        lock(lk) do; out[(k, m)] = res; end
        d = Threads.atomic_add!(done, 1) + 1
        d % 50 == 0 && @printf("  %d/%d done (%.0f min)\n", d, length(tasks), (time() - t0) / 60)
    end

    Tc_pub(j) = min(tab.T[j.row], CAP)
    Tc_jl(r)  = r.tdie < 0 ? CAP : min(r.tdie - MOVE_GEN + off, CAP)
    premove(r) = 0 <= r.tdie < MOVE_GEN          # died before the optimum moved
    npre = count(premove, values(out))
    npre > 0 && println("  $npre Julia run(s) died before the move; excluded from the extinction comparison")

    # per-run results first, so a failure in the summaries below loses nothing
    csv = joinpath(@__DIR__, "published_sweep_results.csv")
    open(csv, "w") do io
        println(io, "dir,file,row,NoQTL,u,d_env,PhenotypeVar,Heritability,MutationScalar,InitFitnessSD,",
                    "pub_gv10000,pub_ns10000,pub_T,jl_qtsim_gv10000,jl_qtsim_ns10000,jl_qtsim_Tcens,jl_all_gv10000,jl_all_ns10000,jl_all_Tcens")
        for (k, j) in enumerate(jobs)
            q = j.q; a = get(out, (k, :qtsim), nothing); b = get(out, (k, :all), nothing)
            fmt(r) = r === nothing ? ",," : "$(r.gv),$(r.ns),$(Tc_jl(r))"
            println(io, join((j.dir, j.file, j.row, j.L, j.u, q.MoveFitnessSD, q.PhenotypeVar, q.Heritability,
                              q.MutationScalar, q.InitFitnessSD, tab.gv[j.row], tab.ns[j.row], tab.T[j.row]), ","),
                    ",", fmt(a), ",", fmt(b))
        end
    end
    println("per-run results written to ", relpath(csv, ROOT))

    println("\n── 2. equilibrium just before the move (table generation 10000; paired by input; Julia vs published) ──")
    @printf("  %-18s %-6s %4s  %-30s %-28s\n", "stratum", "model", "n", "log V_G ratio: mean ± SE (p)", "survivors diff: mean ± SE (p)")
    groups = [(s[1], [k for (k, j) in enumerate(jobs) if j.dir == s[1]], s[5]) for s in STRATA if s[2] <= MAXL]
    push!(groups, ("all strata", collect(eachindex(jobs)), (:qtsim, :all)))
    for (name, ks, models) in groups, m in models
        ks2 = [k for k in ks if haskey(out, (k, m)) && out[(k, m)].gv > 0 && tab.gv[jobs[k].row] > 0]
        length(ks2) < 3 && continue
        lr = [log(out[(k, m)].gv / tab.gv[jobs[k].row]) for k in ks2]
        dn = [out[(k, m)].ns - tab.ns[jobs[k].row] for k in ks2]
        @printf("  %-18s %-6s %4d  %+.3f ± %.3f (p=%.2f)%-8s %+6.1f ± %4.1f (p=%.2f)\n", name, m, length(ks2),
                mean(lr), std(lr) / sqrt(length(lr)), pvalue(OneSampleTTest(lr)), "",
                mean(dn), std(dn) / sqrt(length(dn)), pvalue(OneSampleTTest(Float64.(dn))))
    end

    println("\n── 3. extinction after the move, censored at $CAP generations (paired by input) ──")
    @printf("  %-18s %-6s %4s  %-22s %-22s %-24s\n", "stratum", "model", "n", "median time pub / jl", "survive≥CAP pub / jl", "paired log-time p (signed rank)")
    for (name, ks, models) in groups, m in models
        ks2 = [k for k in ks if haskey(out, (k, m)) && !premove(out[(k, m)])]
        length(ks2) < 3 && continue
        tp = [Tc_pub(jobs[k]) for k in ks2]; tj = [Tc_jl(out[(k, m)]) for k in ks2]
        d = log.(tj) .- log.(tp)
        p = (length(d) < 10 || all(iszero, d)) ? NaN : pvalue(ApproximateSignedRankTest(d))
        @printf("  %-18s %-6s %4d  %7.0f / %-12.0f %5.1f%% / %-14s p=%.2f\n", name, m, length(ks2),
                median(tp), median(tj), 100mean(tp .>= CAP), @sprintf("%.1f%%", 100mean(tj .>= CAP)), p)
    end

    println("\n── headline: median time to extinction (censored at $CAP) by standardised rate d_env/√V_P ──")
    xs = [jobs[k].q.MoveFitnessSD / sqrt(jobs[k].q.PhenotypeVar) for k in eachindex(jobs)]
    edges = quantile(xs, range(0, 1; length=9))
    @printf("  %-17s %4s  %10s  %12s  %10s\n", "d_env/√V_P bin", "n", "published", "Julia :qtsim", "Julia :all")
    for b in 1:8
        ks = [k for k in eachindex(jobs) if (b == 1 ? xs[k] >= edges[b] : xs[k] > edges[b]) && xs[k] <= edges[b+1]]
        isempty(ks) && continue
        med(m) = (v = [Tc_jl(out[(k, m)]) for k in ks if haskey(out, (k, m)) && !premove(out[(k, m)])]; isempty(v) ? "-" : @sprintf("%.0f", median(v)))
        @printf("  %.4f–%-9.4f %4d  %10.0f  %12s  %10s\n", edges[b], edges[b+1], length(ks),
                median([Tc_pub(jobs[k]) for k in ks]), med(:qtsim), med(:all))
    end

    @printf("\ntotal %.0f min\n", (time() - t0) / 60)
end

sweep_main()
