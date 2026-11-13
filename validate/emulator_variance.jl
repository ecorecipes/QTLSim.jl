# Between-run variance in the emulator of the published sweep (roadmap F4), under the
# pre-specified criteria of validate/emulator_variance_criteria.md (written and committed
# before the change was implemented; amendment 1, also before any result, added the noise
# in V_G). Writes validate/emulator_variance.log.
#
#   E1–E5: the criteria of validate/emulator_criteria.md, on a fresh uniform sample of 10,000 (amendment 2)
#          published inputs (seed 20261115), for `emulate(p; replicates = 32)` (primary:
#          :rare_alleles, noise in the mean and in V_G), with verdicts; the deterministic
#          emulator (the first run's cache validate/emulator_runs.csv) and the mean-only noise
#          (vgnoise = false) on the first 5,000 of them as baseline and secondary.
#   E6:    calibration against fresh kernel simulations: 192 inputs chosen from the first
#          run's parameter map by r = MoveFitnessSD/k50 (no table outcome used), each run 40
#          times in the kernel (init = :equilibrium, burnin = 1000, generations 9000–11999).
#   Secondary: :spectrum with the mean noise (R = 8) on the E6 inputs (amendment 2; NSPEC1000 > 0 also runs part of the
#          1000-input sample of the first run, seed 20261113).
#
# Run from the package root, in the vignettes environment (XLSX; Catlab for the moments):
#   julia -t 8 --project=vignettes validate/emulator_variance.jl pilot   # timings only
#   julia -t 8 --project=vignettes validate/emulator_variance.jl         # everything; resumes from the caches
# Caches (resumable `sweep` CSVs): validate/emulator_variance_runs.csv (primary),
# emulator_variance_meanonly.csv, emulator_variance_kernel.csv (E6 kernel runs),
# emulator_variance_e6.csv (emulator on the E6 inputs), emulator_variance_spectrum.csv.

using QTLSim, Statistics, Random, Printf, LinearAlgebra
using Distributions: Binomial, logpdf
import Catlab
include(joinpath(@__DIR__, "published_table.jl"))   # load_master, row_index, match_row

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SIMDIR = joinpath(ROOT, "Repository", "Simulations")
const MASTER = joinpath(ROOT, "Repository", "Results", "Master Output file.xlsx")
const STRATA = [("QT0100-Possion3.0", 100, 3), ("QT0100-Possion5.0", 100, 5), ("QT0500-Possion3.0", 500, 3),
                ("QT0500-Possion5.0", 500, 5), ("QT1000-Possion3.0", 1000, 3), ("QT1000-Possion5.0", 1000, 5),
                ("QT5000-Possion3.0", 5000, 3), ("QT5000-Possion5.0", 5000, 5)]
const MODE = isempty(ARGS) ? "full" : ARGS[1]
const BURNIN, HORIZON, TAIL = 2000, 2000, 240_000
const R, SEED, NSAMPLE, NMEANONLY = 32, 20261115, 10_000, 5_000   # amendment 2 (was 20,000)
const RK, KBURN, NE6PER = 40, 1000, 24
const RSPEC = 8
const NSPEC1000 = parse(Int, get(ENV, "NSPEC1000", "0"))   # how many of the 1000-sample get :spectrum (amendment 2)
const KW = (; burnin = BURNIN, horizon = HORIZON, tail = TAIL)
cache(name) = joinpath(@__DIR__, name)

const LOG = Ref{IO}(stdout)
say(s...) = (println(stdout, s...); LOG[] === stdout || (println(LOG[], s...); flush(LOG[])))
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))

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

# --- statistics, as validate/emulator.jl -------------------------------------------
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
function auc(score, y)
    r = midranks(score); n1 = count(y); n0 = length(y) - n1
    (n1 == 0 || n0 == 0) && return NaN
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
end
spearman(x, y) = cor(midranks(x), midranks(y))
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
function x50(rate, y)
    β, conv = logistic([ones(length(rate)) log.(rate)], Float64.(y))
    exp(-β[1] / β[2]), conv
end

# E1–E5 for emulations `em` of the inputs `ix` (as validate/emulator.jl's `criteria`)
function criteria(em, ix, tab, rows, ps, strat; label)
    P = [e.P2000 for e in em]; P240 = [e.P240k for e in em]; Tp = [e.T for e in em]
    Tobs = [tab.T[rows[i]] for i in ix]
    alive = Tobs .> 2000
    survived = [tab.surv[rows[i]] == 1 for i in ix]
    st = strat[ix]
    res = Dict{Symbol,Bool}()
    say("  n = $(length(ix)); observed alive at 2000: $(count(alive)) ($(round(100mean(alive), digits = 2))%); " *
        "survived 240,000: $(count(survived)); mean P2000 = $(round(mean(P), sigdigits = 4)); " *
        "P2000 exactly 0: $(count(==(0), P)), exactly 1: $(count(==(1), P)), in (0.01, 0.99): $(count(x -> 0.01 < x < 0.99, P))")
    a = auc(P, alive)
    as = [auc(P[st.==s], alive[st.==s]) for s in eachindex(STRATA)]
    res[:E1] = a >= 0.85 && !any(isnan, as) && all(>=(0.80), as)
    sayf("  E1 discrimination: AUC(P2000; alive at 2000) = %.4f (≥ 0.85)  %s", a, label(res[:E1]))
    for s in eachindex(STRATA)
        sayf("       %-18s n = %6d  alive %6d  AUC = %.4f (≥ 0.80)", STRATA[s][1], count(st .== s), count(alive[st.==s]), as[s])
    end
    o = sortperm(P; alg = Base.Sort.DEFAULT_STABLE)
    n = length(P); worst = 0.0
    say("  E2 calibration: deciles of P2000 (ties broken by input order)")
    say("       decile     n   mean P2000   observed   |diff|")
    for d in 1:10
        k = o[(div((d - 1) * n, 10) + 1):div(d * n, 10)]
        mp, ob = mean(P[k]), mean(alive[k])
        worst = max(worst, abs(mp - ob))
        sayf("       %6d %6d   %10.4f   %8.4f   %6.4f", d, length(k), mp, ob, abs(mp - ob))
    end
    β, conv = logistic([ones(n) logit.(clamp.(P, 1e-6, 1 - 1e-6))], Float64.(alive))
    res[:E2] = worst <= 0.10 && 0.8 <= β[2] <= 1.25
    sayf("     largest |observed − mean P2000| = %.4f (≤ 0.10); calibration slope = %.4f, intercept %.4f%s (in [0.8, 1.25])  %s",
         worst, β[2], β[1], conv ? "" : " (not converged)", label(res[:E2]))
    k = findall(i -> Tobs[i] <= 2000 && !isnan(Tp[i]), eachindex(Tobs))
    nobs = count(<=(2000), Tobs)
    if length(k) >= 3
        ρ = spearman(Tp[k], Float64.(Tobs[k])); mr = median(Tp[k] ./ Tobs[k])
        res[:E3] = ρ >= 0.7 && 0.8 <= mr <= 1.25
        sayf("  E3 timing: %d of %d inputs extinct within 2000 in the table have a predicted T; Spearman = %.4f (≥ 0.7), " *
             "median T_pred/T_obs = %.4f (IQR %.3f–%.3f; in [0.8, 1.25])  %s",
             length(k), nobs, ρ, mr, q25(Tp[k] ./ Tobs[k]), q75(Tp[k] ./ Tobs[k]), label(res[:E3]))
    else
        res[:E3] = false
        say("  E3 timing: fewer than 3 inputs with both an observed and a predicted T  $(label(false))")
    end
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
    say("     E4: $(label(e4))")
    sayf("  E5 (exploratory) AUC(P240k; survived 240,000) = %.4f; mean P240k %.4g vs observed survival %.4f",
         auc(P240, survived), mean(P240), mean(survived))
    res
end

# --- E6 ---------------------------------------------------------------------------------
# the inputs: from the first run's parameter map, by r = MoveFitnessSD/k50 (criteria, E6)
function e6_inputs(ids, ps, strat)
    pos = Dict(id => i for (i, id) in enumerate(ids))
    lines = readlines(cache("emulator_map.csv"))
    hdr = split(lines[1], ","); ck = findfirst(==("k50"), hdr)
    cand = Tuple{Int,Float64}[]
    for l in lines[2:end]
        f = split(l, ",")
        k50 = parse(Float64, f[ck]); i = pos[f[1]]
        0 < k50 < Inf || continue
        r = ps[i].MoveFitnessSD / k50
        0.4 <= r <= 2 && push!(cand, (i, r))
    end
    sort!(cand; by = c -> ids[c[1]])                     # input order (stratum, then file name)
    edges = exp.(range(log(0.4), log(2.0); length = 5))
    bin(r) = min(4, searchsortedlast(edges, r))
    rng = MersenneTwister(SEED)
    pick = Int[]
    for s in eachindex(STRATA)
        cs = [c for c in cand if strat[c[1]] == s]
        got = Int[]
        for b in 1:4
            cb = [c[1] for c in cs if bin(c[2]) == b]
            append!(got, cb[randperm(rng, length(cb))[1:min(6, length(cb))]])
        end
        rest = [c[1] for c in cs if !(c[1] in got)]
        need = NE6PER - length(got)
        need > 0 && append!(got, rest[randperm(rng, length(rest))[1:min(need, length(rest))]])
        append!(pick, got)
    end
    pick, Dict(c[1] => c[2] for c in cand)
end

kernel_params(p) = Params(p; NoGenerations = 12000, Sample = 1, Seed = nothing)

function kernel_job(j)
    (; key, p, seed) = j
    o = simulate(kernel_params(p); seed, init = :equilibrium, burnin = KBURN)
    k = findfirst(==(9999), o.gen)
    (; alive = !o.extinct, extinct_gen = o.extinct_gen,
     zbar9999 = k === nothing ? NaN : o.pheno_mean[k], VG9999 = k === nothing ? NaN : o.geno_var[k])
end

# E6's statistics for the emulator's P2000 `P` against kernel successes `s` of RK runs
function e6stats(P, s; label)
    f = s ./ RK; n = length(P)
    res = true
    say("       bin            inputs   mean P2000   kernel fraction   |diff|")
    edges = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]
    for b in 1:6
        lo = edges[b]; hi = b < 6 ? edges[b+1] : Inf
        k = findall(x -> lo <= x < hi || (b == 6 && x == 1.0), P)
        if isempty(k)
            sayf("       [%.1f, %s) %8d", lo, b < 6 ? @sprintf("%.1f", hi) : "1]", 0)
            continue
        end
        d = abs(mean(f[k]) - mean(P[k]))
        ok = length(k) < 10 || d <= 0.10
        res &= ok
        sayf("       [%.1f, %s %8d   %10.4f   %15.4f   %6.4f%s", lo, b < 6 ? @sprintf("%.1f)", hi) : "1]  ", length(k),
             mean(P[k]), mean(f[k]), d, length(k) < 10 ? "  (fewer than 10: no verdict)" : ok ? "" : "  outside")
    end
    Pc = clamp.(P, 1e-6, 1 - 1e-6)
    β, conv = logistic([ones(n) logit.(Pc)], f)
    e6b = 0.8 <= β[2] <= 1.25
    sayf("     E6a binned (|diff| ≤ 0.10 in every bin with ≥ 10 inputs): %s", label(res))
    sayf("     E6b calibration slope = %.4f (intercept %.4f%s; in [0.8, 1.25]): %s", β[2], β[1], conv ? "" : ", not converged", label(e6b))
    brier = mean(@. f * (1 - P)^2 + (1 - f) * P^2)
    ll = sum(logpdf(Binomial(RK, Pc[i]), s[i]) for i in 1:n)
    y = vcat([[trues(s[i]); falses(RK - s[i])] for i in 1:n]...)
    sc = vcat([fill(P[i], RK) for i in 1:n]...)
    sayf("     mean |f − P2000| = %.4f; Brier = %.4f; binomial log-likelihood = %.1f; AUC over runs = %.4f; mean P2000 %.4f vs kernel %.4f",
         mean(abs.(f .- P)), brier, ll, auc(sc, y), mean(P), mean(f))
    say("     E6: $(label(res && e6b))")
    res && e6b
end

# --- main ------------------------------------------------------------------------------
function main()
    t0 = time()
    ids, ps, strat = read_inputs()
    n = length(ps)
    perm = randperm(MersenneTwister(SEED), n)
    S = sort(perm[1:NSAMPLE])
    Smo = sort(perm[1:NMEANONLY])                     # mean-only secondary (amendment 2)

    if MODE == "pilot"
        say("pilot: $(Threads.nthreads()) threads; timings only, no comparison with the table or the kernel")
        rng = MersenneTwister(1)
        pick = vcat([shuffle(rng, findall(==(s), strat))[1:10] for s in eachindex(STRATA)]...)
        emulate_sweep(ps[pick[1:8]]; KW..., replicates = R, seed = 1)
        t = @elapsed emulate_sweep(ps[pick]; KW..., replicates = R, seed = 1)
        sayf("  primary R = %d: %.1f s wall for %d random inputs: %.3f s wall per input", R, t, length(pick), t / length(pick))
        e6, _ = e6_inputs(ids, ps, strat)
        sp = e6[1:6:end]
        t = @elapsed emulate_sweep(ps[sp]; KW..., replicates = R, seed = 1)
        sayf("  primary R = %d on %d E6 inputs (near k50): %.3f s wall per input", R, length(sp), t / length(sp))
        kp = e6[3:24:end]
        sweep(kernel_job, [(; key = 0, p = ps[i], seed = 1) for i in kp[1:2]])
        t = @elapsed rs = sweep(kernel_job, [(; key = k, p = ps[i], seed = 100 + k) for (k, i) in enumerate(repeat(kp, 2))])
        sayf("  kernel: %.1f s wall for %d runs (2 per stratum, near k50): %.2f s wall per run; alive %d", t, length(rs),
             t / length(rs), count(r -> r.alive, rs))
        t = @elapsed emulate_sweep(ps[e6[5:48:end]]; KW..., loci = :spectrum, replicates = RSPEC, seed = 1)
        sayf("  :spectrum R = %d: %.1f s wall per input (E6 inputs)", RSPEC, t / length(e6[5:48:end]))
        sayf("  total %.1f min", (time() - t0) / 60)
        return
    end

    LOG[] = open(joinpath(@__DIR__, "emulator_variance.log"), "w")
    say("validate/emulator_variance.jl — between-run variance in the emulator (roadmap F4), under validate/emulator_variance_criteria.md")
    say("burn-in $BURNIN, horizon $HORIZON, tail $TAIL; R = $R replicates, base seed $SEED; $(Threads.nthreads()) threads; ",
        Libc.strftime("%Y-%m-%d %H:%M", time()))
    tab = load_master(MASTER)
    find = row_index(tab)
    rows = [find(STRATA[strat[i]][2], STRATA[strat[i]][3], ps[i]) for i in 1:n]
    all(i -> rows[i] == match_row(tab, STRATA[strat[i]][2], STRATA[strat[i]][3], ps[i]), S[1:200:end]) ||
        error("row_index disagrees with match_row")
    all(>(0), rows[S]) || error("an input of the sample has no table row")
    say("sample: $(length(S)) of $n inputs (randperm(MersenneTwister($SEED)), sorted); all matched to a table row")
    vl(b) = pf(b); nov(b) = "-"; base(b) = "($(pf(b)) as a baseline)"

    # primary
    say("\n== primary: :rare_alleles, noise in the mean and in V_G, R = $R ==")
    t = @elapsed prim = emulate_sweep(ps[S]; ids = ids[S], loci = :rare_alleles, KW..., replicates = R, seed = SEED,
                                      progress = 2000, cache = cache("emulator_variance_runs.csv"))
    sayf("  (%.0f s wall this session; cached rows are not re-run)", t)
    res = criteria(prim, S, tab, rows, ps, strat; label = vl)

    say("\n== baseline: the deterministic emulator (replicates = 0; the first run's cache) on the same inputs ==")
    det = emulate_sweep(ps[S]; ids = ids[S], loci = :rare_alleles, KW..., cache = cache("emulator_runs.csv"), progress = 5000)
    resd = criteria(det, S, tab, rows, ps, strat; label = base)

    say("\n== secondary: the noise in the mean only (vgnoise = false; the model as first specified), R = $R, on the first $NMEANONLY inputs of the permutation; no verdicts ==")
    t = @elapsed mo = emulate_sweep(ps[Smo]; ids = ids[Smo], loci = :rare_alleles, KW..., replicates = R, seed = SEED,
                                    vgnoise = false, progress = 2000, cache = cache("emulator_variance_meanonly.csv"))
    sayf("  (%.0f s wall this session)", t)
    criteria(mo, Smo, tab, rows, ps, strat; label = nov)

    # E6
    say("\n== E6: calibration against fresh kernel simulations ==")
    e6, rof = e6_inputs(ids, ps, strat)
    say("  $(length(e6)) inputs from the first run's map (0.4 ≤ r = MoveFitnessSD/k50 ≤ 2): per stratum ",
        join([count(i -> strat[i] == s, e6) for s in eachindex(STRATA)], ", "),
        "; r quartiles ", join([@sprintf("%.3f", quantile([rof[i] for i in e6], q)) for q in (0.25, 0.5, 0.75)], ", "))
    jobs = [(; key = "$(ids[i])#$r", p = ps[i], seed = 20261115000 + 100j + r) for (j, i) in enumerate(e6) for r in 1:RK]
    t = @elapsed kr = sweep(kernel_job, jobs; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 400,
                            cache = cache("emulator_variance_kernel.csv"))
    sayf("  kernel: %d runs (%d per input, init = :equilibrium, burnin = %d), %.0f s wall this session", length(kr), RK, KBURN, t)
    s = [count(kr[(j-1)*RK+r].alive for r in 1:RK) for j in eachindex(e6)]
    sayf("  kernel fraction alive at 2000: mean %.4f; inputs with 0 alive %d, all alive %d, in between %d",
         mean(s ./ RK), count(==(0), s), count(==(RK), s), count(x -> 0 < x < RK, s))
    e6ids = ids[e6]
    runs(name; kw...) = emulate_sweep(ps[e6]; ids = ["$name|" * x for x in e6ids], loci = :rare_alleles, KW..., kw...,
                                      cache = cache("emulator_variance_e6.csv"))
    Pp = [x.P2000 for x in runs("primary"; replicates = R, seed = SEED)]
    Pd = [x.P2000 for x in emulate_sweep(ps[e6]; ids = e6ids, loci = :rare_alleles, KW..., cache = cache("emulator_runs.csv"))]
    Pm = [x.P2000 for x in runs("meanonly"; replicates = R, seed = SEED, vgnoise = false)]
    say("\n  -- primary (noise in the mean and in V_G, R = $R) --")
    e6res = e6stats(Pp, s; label = vl)
    say("\n  -- baseline: deterministic emulator --")
    e6stats(Pd, s; label = base)
    say("\n  -- secondary: mean-only noise (no verdict) --")
    e6stats(Pm, s; label = nov)

    # diagnostics of the noise model at generation 9999 (the generation before the move)
    say("\n  -- noise-model diagnostics (no verdicts): kernel replicates alive at generation 9999 vs the emulator at the move --")
    rz, rv, cvk, cve = Float64[], Float64[], Float64[], Float64[]
    for (j, i) in enumerate(e6)
        k = [kr[(j-1)*RK+r] for r in 1:RK if !isnan(kr[(j-1)*RK+r].VG9999)]
        length(k) >= 10 || continue
        p = ps[i]
        q = Params(p; Sample = 1, NoGenerations = BURNIN + HORIZON, MoveFitnessSDGen = BURNIN)
        mr = QTLSim.moments_replicates(q; replicates = 400, rng = Xoshiro(j), from = BURNIN, VG = house_of_gauss_vg(p),
                                       loci = :rare_alleles, vg = QTLSim.vg_fluctuation(q))
        live = [x for x in mr.replicates if !isempty(x.pheno_mean)]
        ze = var([x.pheno_mean[1] for x in live]) + QTLSim.envVar(p) / p.PopSize      # the mean of N phenotypes
        vk = [x.VG9999 for x in k]
        push!(rz, var([x.zbar9999 for x in k]) / ze)
        push!(rv, mean(vk) / mr.det.geno_var[end])
        push!(cvk, std(vk) / mean(vk)); push!(cve, QTLSim.vg_fluctuation(q).cv)
    end
    sayf("  %d inputs with ≥ 10 kernel replicates alive at 9999", length(rz))
    sayf("  variance of the mean phenotype among replicates, kernel/emulator: median %.3f [IQR %.3f, %.3f]", median(rz), q25(rz), q75(rz))
    sayf("  mean V_G, kernel/emulator (deterministic V_G at the move): median %.3f [IQR %.3f, %.3f]", median(rv), q25(rv), q75(rv))
    sayf("  CV of V_G among replicates: kernel median %.3f [IQR %.3f, %.3f]; vg_fluctuation cv median %.3f [IQR %.3f, %.3f]; " *
         "ratio kernel/theory median %.3f", median(cvk), q25(cvk), q75(cvk), median(cve), q25(cve), q75(cve), median(cvk ./ cve))
    say("  (the kernel's CV includes the sampling error of a variance from 40 replicates, about 1/√78 ≈ 0.11 relative)")

    # secondary :spectrum
    say("\n== secondary: :spectrum with the noise in the mean (R = $RSPEC); no verdicts ==")
    t = @elapsed Psp = [x.P2000 for x in emulate_sweep(ps[e6]; ids = ["spectrum|" * x for x in e6ids], loci = :spectrum, KW...,
                                                         replicates = RSPEC, seed = SEED, progress = 50,
                                                         cache = cache("emulator_variance_spectrum.csv"))]
    sayf("  on the E6 inputs (%.0f s wall this session):", t)
    e6stats(Psp, s; label = nov)
    if NSPEC1000 > 0
        s1 = sort(randperm(MersenneTwister(20261113), n)[1:1000])[1:NSPEC1000]
        t = @elapsed sp = emulate_sweep(ps[s1]; ids = ids[s1], loci = :spectrum, KW..., replicates = RSPEC, seed = SEED,
                                        progress = 100, cache = cache("emulator_variance_spectrum.csv"))
        say("\n  on the first $NSPEC1000 inputs (in input order) of the first run's 1000-input sample (seed 20261113):")
        sayf("  (%.0f s wall this session)", t)
        criteria(sp, s1, tab, rows, ps, strat; label = nov)
    end

    say("\nVERDICTS (primary): ", join(["$k $(pf(res[k]))" for k in (:E1, :E2, :E3, :E4)], ", "), ", E6 $(pf(e6res))")
    say("baseline (deterministic) on the same table sample: ", join(["$k $(pf(resd[k]))" for k in (:E1, :E2, :E3, :E4)], ", "))
    sayf("total %.0f min wall", (time() - t0) / 60)
    close(LOG[])
end

main()
