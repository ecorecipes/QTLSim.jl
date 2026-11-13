# The emulator of the published sweep (ROADMAP F4) against the published table, under the
# pre-specified criteria of validate/emulator_criteria.md (written and committed before the
# emulator was implemented or run). Writes validate/emulator.log.
#
# For each of the 119,984 published inputs (Repository/Simulations/QT<L>-Possion<u>/,
# matched to rows of Repository/Results/Master Output file.xlsx by parameters, as
# validate/published_sweep.jl does), `emulate` (src/emulator.jl) gives P2000 (alive 2000
# generations after the move), T (median extinction generation after the move, counted as
# the table counts) and P240k (constant end-hazard extrapolation). Then:
#   E1 AUC of P2000 for alive-at-2000 (table T > 2000): ≥ 0.85 overall, ≥ 0.80 per stratum;
#   E2 deciles of P2000: |observed − mean P2000| ≤ 0.10 in every decile, and calibration
#      slope (logistic regression of the outcome on logit P2000, P2000 clipped to
#      [1e-6, 1 − 1e-6]) in [0.8, 1.25];
#   E3 among table T ≤ 2000 with a predicted T: Spearman ≥ 0.7, median T_pred/T_obs in [0.8, 1.25];
#   E4 per stratum, x50 (rate d/√V_P at which half are alive at 2000, logistic regression
#      on log rate) from the table and from P2000 (fractional response): ratio in [0.8, 1.25];
#   E5 (exploratory) AUC of P240k for survived-240,000.
# Primary closure :rare_alleles on all inputs; :spectrum and :multilocus (and :rare_alleles
# for comparison) on a uniform sample of 1000 (seed 20261113), no verdicts. Then the
# exploratory parameter map on a uniform sample of 2000 (seed 20261114):
# `emulated_half_survival_rate` and its elasticities by `QTLSim.Para`/`sensitivities`.
#
# Run from the package root, in the vignettes environment (XLSX; Catlab for the moments):
#   julia -t 8 --project=vignettes validate/emulator.jl          # everything; resumes from the caches
#   julia -t 8 --project=vignettes validate/emulator.jl pilot    # timing only, no comparison with the table
# Caches (resumable `sweep` CSVs, keyed by "stratum/file"): validate/emulator_runs.csv
# (primary), emulator_spectrum.csv, emulator_multilocus.csv, emulator_map.csv.

using QTLSim, Statistics, Random, Printf, LinearAlgebra
import Catlab
include(joinpath(@__DIR__, "published_table.jl"))   # load_master, row_index

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SIMDIR = joinpath(ROOT, "Repository", "Simulations")
const MASTER = joinpath(ROOT, "Repository", "Results", "Master Output file.xlsx")
const STRATA = [("QT0100-Possion3.0", 100, 3), ("QT0100-Possion5.0", 100, 5), ("QT0500-Possion3.0", 500, 3),
                ("QT0500-Possion5.0", 500, 5), ("QT1000-Possion3.0", 1000, 3), ("QT1000-Possion5.0", 1000, 5),
                ("QT5000-Possion3.0", 5000, 3), ("QT5000-Possion5.0", 5000, 5)]
const MODE = isempty(ARGS) ? "full" : ARGS[1]
const BURNIN, HORIZON, TAIL = 2000, 2000, 240_000          # fixed by the criteria
const MAP_NAMES = (:PhenotypeVar, :Heritability, :MutationScalar, :MutationRatePerGamete, :InitFitnessSD)

# --- output to the terminal and the log ----------------------------------------
const LOG = Ref{IO}(stdout)
say(s...) = (println(stdout, s...); LOG[] === stdout || (println(LOG[], s...); flush(LOG[])))
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))

# --- the inputs ------------------------------------------------------------------
function read_inputs()
    ids, ps, strat = String[], Params[], Int[]
    for (s, (dir, L, u)) in enumerate(STRATA)
        for f in sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, dir))))
            q = Params(read_qtsim_options(joinpath(SIMDIR, dir, f)); Seed = nothing)
            (q.NoQTL == L && q.RecombPoissonMean == u && q.MoveType === :continuous) ||
                error("$dir/$f does not match its stratum or protocol")
            push!(ids, "$dir/$f"); push!(ps, q); push!(strat, s)
        end
    end
    ids, ps, strat
end

# --- statistics (no packages) ----------------------------------------------------
function midranks(x)
    o = sortperm(x); r = similar(x, Float64); i = 1; n = length(x)
    while i <= n
        j = i
        while j < n && x[o[j+1]] == x[o[i]]; j += 1; end
        for k in i:j; r[o[k]] = (i + j) / 2; end
        i = j + 1
    end
    r
end

"AUC = P(score of a positive > score of a negative) + ½P(tie), by midranks."
function auc(score, y)
    r = midranks(score); n1 = count(y); n0 = length(y) - n1
    (n1 == 0 || n0 == 0) && return NaN
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
end

spearman(x, y) = cor(midranks(x), midranks(y))

"Logistic regression of y ∈ [0, 1] (binary or fractional) on the columns of X, by Newton–Raphson with step halving."
function logistic(X, y; maxiter = 200)
    β = zeros(size(X, 2))
    ll(b) = (η = X * b; sum(@. y * η - log1p(exp(-abs(η))) - max(η, 0)))
    cur = ll(β); conv = false
    for _ in 1:maxiter
        η = X * β; μ = @. 1 / (1 + exp(-η)); w = @. μ * (1 - μ)
        H = X' * (w .* X); g = X' * (y .- μ)
        Δ = try
            (Symmetric(H) + 1e-12 * opnorm(H) * I) \ g
        catch
            break
        end
        t = 1.0; new = ll(β + Δ)
        while new < cur - 1e-12 && t > 1e-8
            t /= 2; new = ll(β + t * Δ)
        end
        β += t * Δ
        conv = abs(new - cur) < 1e-10 * (1 + abs(cur)) && maximum(abs, t * Δ) < 1e-8 * (1 + maximum(abs, β))
        cur = new
        conv && break
    end
    β, conv
end

logit(p) = log(p / (1 - p))
q25(v) = quantile(v, 0.25); q75(v) = quantile(v, 0.75)
pf(b) = b ? "PASS" : "FAIL"

# x50 from a logistic regression of y on log rate
function x50(rate, y)
    β, conv = logistic([ones(length(rate)) log.(rate)], Float64.(y))
    exp(-β[1] / β[2]), conv
end

# the criteria's numbers for one set of emulations `em` against outcomes of the inputs `ix`
function criteria(em, ix, tab, rows, ps, strat; verdicts::Bool)
    P = [e.P2000 for e in em]; P240 = [e.P240k for e in em]; Tp = [e.T for e in em]
    Tobs = [tab.T[rows[i]] for i in ix]
    alive = Tobs .> 2000
    survived = [tab.surv[rows[i]] == 1 for i in ix]
    st = strat[ix]
    res = Dict{Symbol,Bool}()
    v(b) = verdicts ? pf(b) : "-"

    say("  n = $(length(ix)); observed alive at 2000: $(count(alive)) ($(round(100mean(alive), digits = 2))%); " *
        "survived 240,000: $(count(survived)); mean P2000 = $(round(mean(P), sigdigits = 4)); " *
        "P2000 exactly 0: $(count(==(0), P)), exactly 1: $(count(==(1), P)), in (0.01, 0.99): $(count(x -> 0.01 < x < 0.99, P))")

    # E1
    a = auc(P, alive)
    as = [auc(P[st.==s], alive[st.==s]) for s in eachindex(STRATA)]
    e1 = a >= 0.85 && all(x -> x >= 0.80, filter(!isnan, as)) && !any(isnan, as)
    res[:E1] = e1
    sayf("  E1 discrimination: AUC(P2000; alive at 2000) = %.4f (≥ 0.85)  %s", a, v(e1))
    for s in eachindex(STRATA)
        sayf("       %-18s n = %6d  alive %6d  AUC = %.4f (≥ 0.80)", STRATA[s][1], count(st .== s), count(alive[st.==s]), as[s])
    end

    # E2
    o = sortperm(P; alg = Base.Sort.DEFAULT_STABLE)          # ties broken by input order
    n = length(P); worst = 0.0
    say("  E2 calibration: deciles of P2000 (ties broken by input order)")
    say("       decile     n   mean P2000   observed   |diff|")
    for d in 1:10
        k = o[(div((d - 1) * n, 10) + 1):div(d * n, 10)]
        mp, ob = mean(P[k]), mean(alive[k])
        worst = max(worst, abs(mp - ob))
        sayf("       %6d %6d   %10.4f   %8.4f   %6.4f", d, length(k), mp, ob, abs(mp - ob))
    end
    Pc = clamp.(P, 1e-6, 1 - 1e-6)
    β, conv = logistic([ones(n) logit.(Pc)], Float64.(alive))
    e2 = worst <= 0.10 && 0.8 <= β[2] <= 1.25
    res[:E2] = e2
    sayf("     largest |observed − mean P2000| = %.4f (≤ 0.10); calibration slope = %.4f, intercept %.4f%s (in [0.8, 1.25])  %s",
         worst, β[2], β[1], conv ? "" : " (not converged)", v(e2))

    # E3
    k = findall(i -> Tobs[i] <= 2000 && !isnan(Tp[i]), eachindex(Tobs))
    nobs = count(<=(2000), Tobs)
    if length(k) >= 3
        ρ = spearman(Tp[k], Float64.(Tobs[k])); mr = median(Tp[k] ./ Tobs[k])
        e3 = ρ >= 0.7 && 0.8 <= mr <= 1.25
        sayf("  E3 timing: %d of %d inputs extinct within 2000 in the table have a predicted T; Spearman = %.4f (≥ 0.7), " *
             "median T_pred/T_obs = %.4f (IQR %.3f–%.3f; in [0.8, 1.25])  %s",
             length(k), nobs, ρ, mr, q25(Tp[k] ./ Tobs[k]), q75(Tp[k] ./ Tobs[k]), v(e3))
    else
        e3 = false
        say("  E3 timing: fewer than 3 inputs with both an observed and a predicted T  $(v(e3))")
    end
    res[:E3] = e3

    # E4
    say("  E4 half-survival rate x50 = d/√V_P (logistic regression on log rate), per stratum")
    say("       stratum                n   x50 table   x50 emulator   ratio   (in [0.8, 1.25])")
    e4 = true
    for s in eachindex(STRATA)
        k = findall(==(s), st)
        isempty(k) && continue
        rate = [ps[ix[i]].MoveFitnessSD / sqrt(ps[ix[i]].PhenotypeVar) for i in k]
        xt, ct = x50(rate, alive[k]); xe, ce = x50(rate, P[k])
        ok = 0.8 <= xe / xt <= 1.25
        e4 &= ok
        sayf("       %-18s %6d   %9.5f   %12.5f   %5.3f%s%s", STRATA[s][1], length(k), xt, xe, xe / xt,
             ok ? "" : "  outside", (ct && ce) ? "" : "  (a regression did not converge)")
    end
    res[:E4] = e4
    say("     E4: $(v(e4))")

    # E5
    sayf("  E5 (exploratory) AUC(P240k; survived 240,000) = %.4f; mean P240k %.4g vs observed survival %.4f",
         auc(P240, survived), mean(P240), mean(survived))
    res
end

# --- the parameter map -----------------------------------------------------------
function map_job((id, p))
    t = @elapsed h = emulated_half_survival_rate(p; burnin = BURNIN, horizon = HORIZON)
    k = h.k50
    el = if 0 < k < Inf
        P = QTLSim.Para(q -> emulated_half_survival_rate(q; burnin = BURNIN, horizon = HORIZON,
                                                         bracket = (k * exp(-0.002), k * exp(0.002)), rtol = 1e-9).k50,
                        p, MAP_NAMES)
        sensitivities(P; rel = 1e-3)
    else
        NamedTuple{MAP_NAMES}(ntuple(_ -> NaN, length(MAP_NAMES)))
    end
    (; k50 = k, x50 = k / sqrt(p.PhenotypeVar), kc = h.kc, el..., seconds = t)
end

# --- main ------------------------------------------------------------------------
function main()
    t0 = time()
    ids, ps, strat = read_inputs()
    n = length(ps)
    kw = (; burnin = BURNIN, horizon = HORIZON, tail = TAIL)

    if MODE == "pilot"
        rng = MersenneTwister(1)
        pick = vcat([shuffle(rng, findall(==(s), strat))[1:40] for s in eachindex(STRATA)]...)
        say("pilot: $(length(pick)) inputs (40 per stratum), $(Threads.nthreads()) threads; timings only")
        for loci in (:rare_alleles, :multilocus)
            emulate_sweep(ps[pick[1:8]]; loci, kw...)                           # compile
            t = @elapsed emulate_sweep(ps[pick]; loci, kw...)
            sayf("  %-13s %.1f s wall for %d inputs: %.4f s wall per input", loci, t, length(pick), t / length(pick))
        end
        sp = pick[1:5:end]
        t = @elapsed emulate_sweep(ps[sp]; loci = :spectrum, kw...)
        sayf("  %-13s %.1f s wall for %d inputs: %.3f s wall per input", :spectrum, t, length(sp), t / length(sp))
        mp = pick[3:20:end]
        map_job((ids[mp[1]], ps[mp[1]]))
        t = @elapsed rs = sweep(map_job, collect(zip(ids[mp], ps[mp])))
        sayf("  parameter map %.1f s wall for %d inputs: %.2f s wall per input; per-input single-thread seconds for k50 alone %s",
             t, length(mp), t / length(mp), join([@sprintf("%.2f", r.seconds) for r in rs], " "))
        return
    end

    LOG[] = open(joinpath(@__DIR__, "emulator.log"), "w")
    say("validate/emulator.jl — the emulator of the published sweep (roadmap F4) under validate/emulator_criteria.md")
    say("burn-in $BURNIN, horizon $HORIZON, tail $TAIL; $(Threads.nthreads()) threads; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    tab = load_master(MASTER)
    find = row_index(tab)
    rows = [find(STRATA[strat[i]][2], STRATA[strat[i]][3], ps[i]) for i in 1:n]
    matched = findall(>(0), rows)
    dup = length(matched) - length(unique(rows[matched]))
    chk = 1:600:n                                     # the indexed matching is match_row's
    all(i -> rows[i] == match_row(tab, STRATA[strat[i]][2], STRATA[strat[i]][3], ps[i]), chk) ||
        error("row_index disagrees with match_row")
    say("inputs: $n files; matched to a unique table row: $(length(matched)); unmatched: $(n - length(matched)); " *
        "rows matched twice: $dup; table rows $(length(tab.L))")

    say("\n== primary closure :rare_alleles, all matched inputs ==")
    t = @elapsed prim = emulate_sweep(ps[matched]; ids = ids[matched], loci = :rare_alleles, progress = 10_000,
                                      cache = joinpath(@__DIR__, "emulator_runs.csv"), kw...)
    sayf("  (%.0f s wall this session; cached rows are not re-run)", t)
    res = criteria(prim, matched, tab, rows, ps, strat; verdicts = true)
    say("\n  VERDICTS (primary): ", join(["$k $(pf(res[k]))" for k in (:E1, :E2, :E3, :E4)], ", "))

    # secondary closures on the pre-specified sample of 1000 inputs (drawn from all n files)
    s1 = sort(randperm(MersenneTwister(20261113), n)[1:1000])
    s1m = filter(i -> rows[i] > 0, s1)
    say("\n== secondary closures on the uniform sample of 1000 inputs (seed 20261113; $(length(s1m)) matched); no verdicts ==")
    pos = Dict(i => k for (k, i) in enumerate(matched))
    for loci in (:rare_alleles, :multilocus, :spectrum)
        say("\n-- loci = :$loci --")
        em = loci === :rare_alleles ? prim[[pos[i] for i in s1m]] :
             emulate_sweep(ps[s1m]; ids = ids[s1m], loci, progress = 100,
                           cache = joinpath(@__DIR__, "emulator_$(loci).csv"), kw...)
        criteria(em, s1m, tab, rows, ps, strat; verdicts = false)
    end

    # the parameter map
    s2 = sort(randperm(MersenneTwister(20261114), n)[1:2000])
    say("\n== the parameter map (exploratory) on the uniform sample of 2000 inputs (seed 20261114) ==")
    say("  k50: emulated_half_survival_rate (rate per generation at which P2000 = 1/2; x50 = k50/√V_P);")
    say("  elasticities ∂log k50/∂log θ by sensitivities (central differences in log θ, rel = 1e-3, one Richardson level)")
    t = @elapsed mp = sweep(map_job, collect(zip(ids[s2], ps[s2])); id = first, progress = 100,
                            cache = joinpath(@__DIR__, "emulator_map.csv"))
    sayf("  (%.0f s wall this session)", t)
    fin = [0 < r.k50 < Inf for r in mp]
    say("  k50 finite and positive for $(count(fin)) of $(length(mp)); 0 for $(count(r -> r.k50 == 0, mp)), Inf for $(count(r -> r.k50 == Inf, mp))")
    # the table's own x50 per stratum for reference (from E4's regression on all inputs)
    say("\n  medians [IQR] by stratum (finite k50 only)")
    cols = (:x50, :k50, MAP_NAMES...)
    say("  ", rpad("stratum", 18), " ", lpad("n", 5), join([lpad(string(c), 28) for c in cols]))
    for s in [eachindex(STRATA); 0]
        k = [j for (j, i) in enumerate(s2) if (s == 0 || strat[i] == s) && fin[j]]
        isempty(k) && continue
        cells = map(cols) do c
            x = filter(isfinite, [getfield(mp[j], c) for j in k])
            isempty(x) ? lpad("-", 28) : lpad(@sprintf("%.4g [%.4g, %.4g]", median(x), q25(x), q75(x)), 28)
        end
        say("  ", rpad(s == 0 ? "all" : STRATA[s][1], 18), " ", lpad(string(length(k)), 5), join(cells))
    end
    say("  (the elasticity of x50 = k50/√V_P in PhenotypeVar is that of k50 minus 1/2; NoQTL and RecombPoissonMean")
    say("   enter only through the strata above.)")
    sayf("\ntotal %.0f min wall", (time() - t0) / 60)
    close(LOG[])
end

main()
