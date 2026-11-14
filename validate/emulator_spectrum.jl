# The emulator of the published sweep with the `:spectrum` closure (roadmap F4), under the
# pre-specified criteria of validate/emulator_spectrum_criteria.md (committed before this
# script was written or anything was run). Writes validate/emulator_spectrum.log.
#
#   E1–E5: the criteria of validate/emulator_criteria.md on 2400 fresh inputs (300 per
#          stratum, MersenneTwister(20261116)), excluding every input whose :spectrum result
#          was looked at before; primary `emulate(p; loci = :spectrum)` (deterministic), the
#          deterministic :rare_alleles emulator as a baseline.
#   E6′:   160 new inputs (0.3 ≤ MoveFitnessSD/k50(:rare_alleles) ≤ 1.5, amendment 1;
#          MersenneTwister(20261117/8)), disjoint from all of the above, each run 40 times in the
#          kernel exactly as E6 of validate/emulator_variance.jl (inputs from a new pool, amendment 1); calibration (E6′a, b), a
#          threshold test (E6′c) and x50 by stratum (E6′d). Baseline :rare_alleles and
#          secondary :spectrum with the mean noise (R = 8) on the same inputs, no verdicts.
#
# Run from the package root, in the vignettes environment (XLSX; Catlab for the moments):
#   julia -t 8 --project=vignettes validate/emulator_spectrum.jl select   # the samples only (incl. the pool k50s), no emulation of the samples
#   julia -t 8 --project=vignettes validate/emulator_spectrum.jl          # everything; resumes from the caches
# Caches (resumable `sweep` CSVs): validate/emulator_spectrum_runs.csv (emulations, ids
# prefixed by model), emulator_spectrum_pool.csv (E6′ pool k50s) and emulator_spectrum_kernel.csv (E6′ kernel runs).

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
const KW = (; burnin = BURNIN, horizon = HORIZON, tail = TAIL)
const NPER, SAMPLE_SEED = 300, 20261116                      # E1–E4 sample
const E6SEED, E6SEED2, NPOOL, NBIN_PER, NE6PER, RLO, RHI = 20261117, 20261118, 150, 5, 20, 0.3, 1.5   # amendment 1
const RK, KBURN, KSEED = 40, 1000, 30261116000
const RSPEC, NOISE_SEED = 8, 20261116
cache(name) = joinpath(@__DIR__, name)
const RUNS = cache("emulator_spectrum_runs.csv")

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

# --- statistics, as validate/emulator.jl and validate/emulator_variance.jl ----------------
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

"""
x50 under E6′d's rule: the logistic regression's exp(−β₀/β₁); on non-convergence with a
complete separation (responses ≥ 0.5 all at lower rates than responses < 0.5) the geometric
mean of the two rates at the gap; NaN (undefined) otherwise, or if β₁ ≥ 0, or if all
responses are on one side of 0.5. Returns (x50, how).
"""
function x50_rule(rate, y)
    hi = y .>= 0.5
    (all(hi) || !any(hi)) && return NaN, "all on one side of 0.5"
    β, conv = logistic([ones(length(rate)) log.(rate)], Float64.(y))
    if conv
        β[2] >= 0 && return NaN, "slope ≥ 0"
        return exp(-β[1] / β[2]), "regression"
    end
    a, b = maximum(rate[hi]), minimum(rate[.!hi])
    a < b ? (sqrt(a * b), "separation") : (NaN, "not converged, no separation")
end

function bootstrap_ratio(rate, ya, yb, rng; B = 1000, f = x50)
    n = length(rate); out = Float64[]
    for _ in 1:B
        k = rand(rng, 1:n, n)
        xa = f(rate[k], ya[k])[1]; xb = f(rate[k], yb[k])[1]
        isfinite(xa) && isfinite(xb) && xb > 0 && push!(out, xa / xb)
    end
    length(out) < B ÷ 2 ? (NaN, NaN, length(out)) : (quantile(out, 0.025), quantile(out, 0.975), length(out))
end

# E1–E5 for emulations `em` of the inputs `ix` (as validate/emulator.jl's `criteria`), plus
# E4's bootstrap intervals and the table's x50 on the whole stratum
function criteria(em, ix, tab, rows, ps, strat; label, alltab = nothing)
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
    say("       stratum                n   x50 table   x50 emulator   ratio   (in [0.8, 1.25])   bootstrap 95%     table x50, whole stratum")
    e4 = true
    brng = MersenneTwister(1)
    for s in eachindex(STRATA)
        k = findall(==(s), st)
        isempty(k) && continue
        rate = [ps[ix[i]].MoveFitnessSD / sqrt(ps[ix[i]].PhenotypeVar) for i in k]
        xt, ct = x50(rate, alive[k]); xe, ce = x50(rate, P[k])
        ok = 0.8 <= xe / xt <= 1.25
        e4 &= ok
        lo, hi, _ = bootstrap_ratio(rate, P[k], alive[k], brng)
        whole = alltab === nothing ? NaN : alltab[s]
        sayf("       %-18s %6d   %9.5f   %12.5f   %5.3f%s   [%.3f, %.3f]   %.5f%s", STRATA[s][1], length(k), xt, xe, xe / xt,
             ok ? "          " : "  outside ", lo, hi, whole, (ct && ce) ? "" : "  (a regression did not converge)")
    end
    res[:E4] = e4
    say("     E4: $(label(e4))")
    sayf("  E5 (exploratory) AUC(P240k; survived 240,000) = %.4f; mean P240k %.4g vs observed survival %.4f",
         auc(P240, survived), mean(P240), mean(survived))
    res
end

# --- exclusions and samples ----------------------------------------------------------------
function excluded_ids()
    ex = Set{String}()
    for f in ("emulator_spectrum.csv", "emulator_variance_spectrum.csv", "emulator_variance_e6.csv")
        for l in readlines(cache(f))[2:end]
            push!(ex, String(last(split(split(l, ",")[1], "|"))))
        end
    end
    for l in readlines(cache("emulator_variance_kernel.csv"))[2:end]
        push!(ex, String(first(split(split(l, ",")[1], "#"))))
    end
    ex
end

function main_sample(ids, strat, ex)
    rng = MersenneTwister(SAMPLE_SEED)
    S = Int[]
    for s in eachindex(STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        append!(S, el[randperm(rng, length(el))[1:NPER]])
    end
    sort(S)
end

# E6′'s inputs (criteria, amendment 1): a pool of NPOOL per stratum prefiltered by
# MoveFitnessSD/kc in [0.3, 3], their :rare_alleles k50, then NBIN_PER per bin of r = MoveFitnessSD/k50
function e6p_pool(ids, ps, strat, ex)
    rng = MersenneTwister(E6SEED)
    pool = Int[]
    for s in eachindex(STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        got = 0
        for i in el[randperm(rng, length(el))]
            kc = drift_critical_rate(ps[i]).k
            (kc > 0 && 0.3 <= ps[i].MoveFitnessSD / kc <= 3.0) || continue
            push!(pool, i); got += 1
            got == NPOOL && break
        end
    end
    pool
end

function e6p_inputs(ids, ps, strat, ex)
    pool = e6p_pool(ids, ps, strat, ex)
    t = @elapsed k50s = sweep(i -> (; k50 = emulated_half_survival_rate(ps[i]; burnin = BURNIN, horizon = HORIZON).k50), pool;
                              id = i -> ids[i], cost = i -> ps[i].NoQTL, progress = 300,
                              cache = cache("emulator_spectrum_pool.csv"))
    sayf("E6′ pool: k50 (:rare_alleles) of %d inputs, %.0f s wall this session", length(pool), t)
    cand = Tuple{Int,Float64}[]
    for (i, k) in zip(pool, k50s)
        k50 = k.k50
        0 < k50 < Inf || continue
        r = ps[i].MoveFitnessSD / k50
        RLO <= r <= RHI && push!(cand, (i, r))
    end
    sort!(cand; by = c -> ids[c[1]])
    edges = exp.(range(log(RLO), log(RHI); length = 5))
    bin(r) = min(4, searchsortedlast(edges, r))
    rng = MersenneTwister(E6SEED2)
    pick = Int[]; nb = zeros(Int, length(STRATA), 4); ncand = zeros(Int, length(STRATA), 4)
    for s in eachindex(STRATA)
        cs = [c for c in cand if strat[c[1]] == s]
        got = Int[]
        for b in 1:4
            cb = [c[1] for c in cs if bin(c[2]) == b]
            ncand[s, b] = length(cb)
            g = cb[randperm(rng, length(cb))[1:min(NBIN_PER, length(cb))]]
            nb[s, b] = length(g)
            append!(got, g)
        end
        rest = [c[1] for c in cs if !(c[1] in got)]
        need = NE6PER - length(got)
        need > 0 && append!(got, rest[randperm(rng, length(rest))[1:min(need, length(rest))]])
        append!(pick, got)
    end
    pick, Dict(c[1] => c[2] for c in cand), nb, ncand
end

kernel_params(p) = Params(p; NoGenerations = 12000, Sample = 1, Seed = nothing)
function kernel_job(j)
    (; key, p, seed) = j
    o = simulate(kernel_params(p); seed, init = :equilibrium, burnin = KBURN)
    k = findfirst(==(9999), o.gen)
    (; alive = !o.extinct, extinct_gen = o.extinct_gen,
     zbar9999 = k === nothing ? NaN : o.pheno_mean[k], VG9999 = k === nothing ? NaN : o.geno_var[k])
end

# E6′'s statistics for P2000 `P` against kernel successes `s` (of RK) on inputs with rates `rate` and strata `st`
function e6pstats(P, s, rate, st; label)
    f = s ./ RK; n = length(P)
    res = Dict{Symbol,Bool}()
    # E6′a
    a = true
    say("    E6′a binned calibration")
    say("       bin            inputs   mean P2000   kernel fraction   |diff|")
    edges = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]
    for b in 1:6
        lo = edges[b]; hi = b < 6 ? edges[b+1] : Inf
        k = findall(x -> lo <= x < hi || (b == 6 && x == 1.0), P)
        if isempty(k)
            sayf("       [%.1f, %s %8d", lo, b < 6 ? @sprintf("%.1f)", hi) : "1]  ", 0)
            continue
        end
        d = abs(mean(f[k]) - mean(P[k]))
        ok = length(k) < 10 || d <= 0.10
        a &= ok
        sayf("       [%.1f, %s %8d   %10.4f   %15.4f   %6.4f%s", lo, b < 6 ? @sprintf("%.1f)", hi) : "1]  ", length(k),
             mean(P[k]), mean(f[k]), d, length(k) < 10 ? "  (fewer than 10: no verdict)" : ok ? "" : "  outside")
    end
    res[:a] = a
    sayf("     E6′a (|diff| ≤ 0.10 in every bin with ≥ 10 inputs): %s", label(a))
    # E6′b
    Pc = clamp.(P, 1e-6, 1 - 1e-6)
    β, conv = logistic([ones(n) logit.(Pc)], f)
    res[:b] = 0.8 <= β[2] <= 1.25
    sayf("     E6′b calibration slope = %.4f (intercept %.4f%s; in [0.8, 1.25]): %s", β[2], β[1], conv ? "" : ", not converged", label(res[:b]))
    # E6′c
    M = f .> 0.5; C = P .> 0.5
    acc = mean(C .== M); ac = auc(P, M)
    tp, fn, fp, tn = count(C .& M), count(.!C .& M), count(C .& .!M), count(.!C .& .!M)
    bacc = (tp / max(tp + fn, 1) + tn / max(tn + fp, 1)) / 2
    res[:c] = acc >= 0.85 && !isnan(ac) && ac >= 0.85
    say("    E6′c threshold test: kernel majority alive (f > 0.5) vs emulator P2000 > 0.5")
    sayf("       emulator alive & kernel alive %d, emulator alive & kernel dead %d, emulator dead & kernel alive %d, both dead %d",
         tp, fp, fn, tn)
    sayf("       accuracy %.4f (≥ 0.85), AUC(P2000; f > 0.5) %.4f (≥ 0.85), balanced accuracy %.4f: %s", acc, ac, bacc, label(res[:c]))
    say("       by stratum (no verdict): ", join([@sprintf("%s %.2f", STRATA[s][1][3:6] * "/" * STRATA[s][1][end-2], mean(C[st.==s] .== M[st.==s]))
                                             for s in eachindex(STRATA) if any(st .== s)], ", "))
    # E6′d
    say("    E6′d x50 = d/√V_P by stratum (kernel: runs; emulator: P2000 as a fractional response)")
    say("       stratum              n   x50 kernel (how)            x50 emulator (how)          ratio   (in [0.8, 1.25])   bootstrap 95%")
    d = true
    brng = MersenneTwister(2)
    for s in eachindex(STRATA)
        k = findall(==(s), st)
        isempty(k) && continue
        xk, hk = x50_rule(rate[k], f[k]); xe, he = x50_rule(rate[k], P[k])
        ok = isfinite(xk) && isfinite(xe) && 0.8 <= xe / xk <= 1.25
        d &= ok
        lo, hi, nb = bootstrap_ratio(rate[k], P[k], f[k], brng; f = x50_rule)
        sayf("       %-18s %3d   %8.5f (%-15s   %8.5f (%-15s   %6.3f%s   [%.3f, %.3f] (%d of 1000 defined)", STRATA[s][1], length(k),
             xk, hk * ")", xe, he * ")", xe / xk, ok ? "          " : "  outside ", lo, hi, nb)
    end
    res[:d] = d
    xk, hk = x50_rule(rate, f); xe, he = x50_rule(rate, P)
    sayf("       pooled (no verdict) %3d   %8.5f (%s)   %8.5f (%s)   ratio %.3f", n, xk, hk, xe, he, xe / xk)
    sayf("     E6′d (ratio in [0.8, 1.25] in every stratum): %s", label(d))
    brier = mean(@. f * (1 - P)^2 + (1 - f) * P^2)
    ll = sum(logpdf(Binomial(RK, Pc[i]), s[i]) for i in 1:n)
    y = vcat([[trues(s[i]); falses(RK - s[i])] for i in 1:n]...)
    sc = vcat([fill(P[i], RK) for i in 1:n]...)
    sayf("     mean |f − P2000| = %.4f; Brier = %.4f; binomial log-likelihood = %.1f; AUC over runs = %.4f; mean P2000 %.4f vs kernel %.4f",
         mean(abs.(f .- P)), brier, ll, auc(sc, y), mean(P), mean(f))
    all_ = res[:a] && res[:b] && res[:c] && res[:d]
    say("     E6′: $(label(all_))  (a $(label(res[:a])), b $(label(res[:b])), c $(label(res[:c])), d $(label(res[:d])))")
    res[:all] = all_
    res
end

# --- main ------------------------------------------------------------------------------
function main()
    t0 = time()
    ids, ps, strat = read_inputs()
    n = length(ps)
    ex = excluded_ids()
    s1000 = ids[sort(randperm(MersenneTwister(20261113), n)[1:1000])]
    all(in(ex), s1000) || error("the exclusions do not contain the 1000-input sample of seed 20261113")
    ne6 = length(Set(String(first(split(split(l, ",")[1], "#"))) for l in readlines(cache("emulator_variance_kernel.csv"))[2:end]))
    ne6 == 192 || error("expected 192 E6 inputs, found $ne6")
    S = main_sample(ids, strat, ex)
    MODE == "select" || (LOG[] = open(joinpath(@__DIR__, "emulator_spectrum.log"), "w"))
    say("validate/emulator_spectrum.jl — the emulator with loci = :spectrum, under validate/emulator_spectrum_criteria.md")
    say("burn-in $BURNIN, horizon $HORIZON, tail $TAIL; $(Threads.nthreads()) threads; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    say("excluded inputs: $(length(ex)) (contains the 1000-input sample of seed 20261113 and the 192 E6 inputs)")
    say("E1–E4 sample: $(length(S)) inputs ($NPER per stratum, MersenneTwister($SAMPLE_SEED)); overlap with the exclusions $(count(i -> ids[i] in ex, S))")
    E, rof, nb, ncand = e6p_inputs(ids, ps, strat, union(ex, Set(ids[S])))
    say("E6′ inputs: $(length(E)) (amendment 1: a pool of $NPOOL per stratum with MoveFitnessSD/kc in [0.3, 3], MersenneTwister($E6SEED); " *
        "$RLO ≤ r = MoveFitnessSD/k50(:rare_alleles) ≤ $RHI, MersenneTwister($E6SEED2)); " *
        "overlap with the exclusions $(count(i -> ids[i] in ex, E)), with the E1–E4 sample $(length(intersect(E, S)))")
    say("  candidates per stratum and bin of r: " * join([join(ncand[s, :], "/") for s in eachindex(STRATA)], ", "))
    say("  drawn per stratum and bin: " * join([join(nb[s, :], "/") for s in eachindex(STRATA)], ", ") *
        "; per stratum " * join([count(i -> strat[i] == s, E) for s in eachindex(STRATA)], ", "))
    say("  r quartiles " * join([@sprintf("%.3f", quantile([rof[i] for i in E], q)) for q in (0.25, 0.5, 0.75)], ", "))
    MODE == "select" && return

    tab = load_master(MASTER)
    find = row_index(tab)
    rows = [find(STRATA[strat[i]][2], STRATA[strat[i]][3], ps[i]) for i in 1:n]
    all(i -> rows[i] == match_row(tab, STRATA[strat[i]][2], STRATA[strat[i]][3], ps[i]), S[1:50:end]) ||
        error("row_index disagrees with match_row")
    all(>(0), rows[S]) || error("an input of the sample has no table row")
    alltab = map(eachindex(STRATA)) do s
        k = findall(==(s), strat)
        x50([ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in k], [tab.T[rows[i]] > 2000 for i in k])[1]
    end
    vl(b) = pf(b); base(b) = "($(pf(b)) as a baseline)"; nov(b) = "-"

    # E1–E5
    say("\n== primary: loci = :spectrum, deterministic, on the $(length(S)) fresh inputs ==")
    t = @elapsed prim = emulate_sweep(ps[S]; ids = ["spectrum|" * x for x in ids[S]], loci = :spectrum, KW...,
                                      progress = 200, cache = RUNS)
    sayf("  (%.0f s wall this session; cached rows are not re-run)", t)
    res = criteria(prim, S, tab, rows, ps, strat; label = vl, alltab)

    say("\n== baseline: loci = :rare_alleles, deterministic, on the same inputs ==")
    t = @elapsed ra = emulate_sweep(ps[S]; ids = ["rare_alleles|" * x for x in ids[S]], loci = :rare_alleles, KW...,
                                    progress = 1000, cache = RUNS)
    sayf("  (%.0f s wall this session)", t)
    resb = criteria(ra, S, tab, rows, ps, strat; label = base, alltab)

    # E6′
    say("\n== E6′: against fresh kernel simulations on $(length(E)) new inputs ==")
    jobs = [(; key = "$(ids[i])#$r", p = ps[i], seed = KSEED + 100j + r) for (j, i) in enumerate(E) for r in 1:RK]
    t = @elapsed kr = sweep(kernel_job, jobs; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 400,
                            cache = cache("emulator_spectrum_kernel.csv"))
    sayf("  kernel: %d runs (%d per input, init = :equilibrium, burnin = %d), %.0f s wall this session", length(kr), RK, KBURN, t)
    s = [count(kr[(j-1)*RK+r].alive for r in 1:RK) for j in eachindex(E)]
    sayf("  kernel fraction alive at 2000: mean %.4f; inputs with 0 alive %d, all alive %d, in between %d",
         mean(s ./ RK), count(==(0), s), count(==(RK), s), count(x -> 0 < x < RK, s))
    rate = [ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in E]
    st = strat[E]
    em(model; kw...) = [x.P2000 for x in emulate_sweep(ps[E]; ids = ["e6p-$model|" * x for x in ids[E]], KW..., kw...,
                                                         progress = 50, cache = RUNS)]
    t = @elapsed begin
        Ps = em("spectrum"; loci = :spectrum)
        Pr = em("rare_alleles"; loci = :rare_alleles)
        Pn = em("spectrum-noise"; loci = :spectrum, replicates = RSPEC, seed = NOISE_SEED)
    end
    sayf("  emulations on the E6′ inputs: %.0f s wall this session", t)
    say("\n  -- primary: :spectrum, deterministic --")
    r6 = e6pstats(Ps, s, rate, st; label = vl)
    say("\n  -- baseline: :rare_alleles, deterministic --")
    r6b = e6pstats(Pr, s, rate, st; label = base)
    say("\n  -- secondary: :spectrum with the mean noise, R = $RSPEC (no verdict) --")
    e6pstats(Pn, s, rate, st; label = nov)

    say("\nVERDICTS (primary :spectrum): ", join(["$k $(pf(res[k]))" for k in (:E1, :E2, :E3, :E4)], ", "),
        ", E6′ $(pf(r6[:all])) (E6′a $(pf(r6[:a])), E6′b $(pf(r6[:b])), E6′c $(pf(r6[:c])), E6′d $(pf(r6[:d])))")
    say("baseline :rare_alleles (no verdicts): ", join(["$k $(pf(resb[k]))" for k in (:E1, :E2, :E3, :E4)], ", "),
        ", E6′ $(pf(r6b[:all])) (a $(pf(r6b[:a])), b $(pf(r6b[:b])), c $(pf(r6b[:c])), d $(pf(r6b[:d])))")
    sayf("total %.0f min wall this session", (time() - t0) / 60)
    close(LOG[])
end

main()
