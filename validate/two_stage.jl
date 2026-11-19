# The two-stage emulator (`emulate(p; loci = :two_stage)`, src/two_stage.jl) under the
# pre-specified criteria of validate/two_stage_criteria.md (committed before this was run on
# fresh data). Writes validate/two_stage.log.
#
#   P:  fresh published inputs near the threshold (E6′'s procedure on a new pool), 40 fresh kernel
#       runs each: (T1) threshold accuracy/AUC against the block model alone, (T2) x50, (T3) cost.
#       Comparators on the same inputs: the block model alone, ld = :dynamic alone, default :spectrum.
#   S:  a fresh random sample of the published table (250 per stratum), table outcomes as truth:
#       the stage-2 fraction, E1–E4 (well-defined rules, see the criteria), total cost.
#
# Run from the package root (vignettes environment: XLSX, Catlab):
#   julia -t 8 --project=vignettes validate/two_stage.jl select   # exclusions and samples only (incl. the pool's k50s)
#   julia -t 8 --project=vignettes validate/two_stage.jl          # everything; resumes from the caches
# Caches (resumable `sweep` CSVs, git-ignored): validate/two_stage_raw/.

let src = read(joinpath(@__DIR__, "blocks.jl"), String)
    # blocks.jl's functions (and through it spectrum_ld.jl's and E6′'s statistics), without its main
    include_string(Main, replace(src, r"\nbmain\(\)\s*$" => "\n"), joinpath(@__DIR__, "blocks.jl"))
end

# 2026-11-19: `ne = :falconer`, `solver = :iteration` and `QTLSim.burnin_start_vg` below keep the effective size, the solver and
# the starting variance this script was run with (the package's defaults are now `:breeders` and `:bisection`; NEWS.md), so that
# it still reproduces its log.
const TMODE = isempty(ARGS) ? "full" : ARGS[1]
# the configuration under test (validate/two_stage_dev.log; fixed in validate/two_stage_criteria.md)
const TBAND = QTLSim.TWO_STAGE_BAND
const TSCREEN = QTLSim.TWO_STAGE_SCREEN
const TPOOLSEED, TPICKSEED, TSSEED = 20261141, 20261142, 20261143
const TKSEED = 51261116000                       # kernel seeds TKSEED + 100j + r
const TBSEED, TSBSEED = 20261144, 20261145       # job_seed bases: the block model on P, on S
const TNPOOL, TNS, TRK, TRB = 100, 250, 40, 40   # pool per stratum, S per stratum, kernel runs per P input, block replicates
const TRAW = joinpath(@__DIR__, "two_stage_raw")
tcache(name) = (mkpath(TRAW); joinpath(TRAW, name))

# exclusions: every input id in any file or raw-run name under validate/ of the main package and of
# every wt-* worktree except wt-paper (incl. wt-blocks-rep's blocks_raw/: blocks.jl's H, F and pool,
# and this task's development inputs), except files naming > 50,000 ids, and this validation's own
# caches and log; united with spectrum_ld's exclusions (E6′'s pools and samples, the diagnosis,
# the published sweep comparison, vignettes 09/11)
function t_scanned_ids(; maxids = 50_000)
    ex = Set{String}(); src = Dict{String,Int}(); skipped = String[]
    dirs = [d for d in readdir(ROOTWT; join = true) if isdir(joinpath(d, "validate")) &&
            (startswith(basename(d), "wt-") || basename(d) == "julia") && basename(d) != "wt-paper"]
    for d in dirs, (root, _, files) in walkdir(joinpath(d, "validate"))
        for f in files
            path = joinpath(root, f)
            (occursin("two_stage_raw", path) || basename(path) == "two_stage.log") && continue
            ids = Set{String}()
            for m in eachmatch(IDRX, f)
                push!(ids, normid(m.match))
            end
            if any(endswith(f, e) for e in (".csv", ".log", ".txt", ".md", ".jl"))
                filesize(path) < 200_000_000 || continue
                for m in eachmatch(IDRX, read(path, String))
                    push!(ids, normid(m.match))
                end
            end
            isempty(ids) && continue
            if length(ids) > maxids
                push!(skipped, "$(relpath(path, ROOTWT)) ($(length(ids)) ids)")
                continue
            end
            union!(ex, ids)
            key = relpath(root, ROOTWT)
            src[key] = get(src, key, 0) + length(ids)
        end
    end
    ex, src, skipped
end

# P: E6′'s procedure on a new pool (as blocks.jl's H, new seeds)
function tp_inputs(ids, ps, strat, ex)
    rng = MersenneTwister(TPOOLSEED)
    pool = Int[]
    for s in eachindex(M6.STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        got = 0
        for i in el[randperm(rng, length(el))]
            kc = drift_critical_rate(ps[i]; ne = :falconer, solver = :iteration).k
            (kc > 0 && 0.3 <= ps[i].MoveFitnessSD / kc <= 3.0) || continue
            push!(pool, i); got += 1
            got == TNPOOL && break
        end
    end
    k50s = sweep(i -> (; k50 = emulated_half_survival_rate(ps[i]; burnin = BURNIN, horizon = HORIZON).k50), pool;
                 id = i -> ids[i], cost = i -> ps[i].NoQTL, progress = 200, cache = tcache("two_stage_pool.csv"))
    cand = Tuple{Int,Float64}[]
    for (i, k) in zip(pool, k50s)
        0 < k.k50 < Inf || continue
        r = ps[i].MoveFitnessSD / k.k50
        M6.RLO <= r <= M6.RHI && push!(cand, (i, r))
    end
    sort!(cand; by = c -> ids[c[1]])
    edges = exp.(range(log(M6.RLO), log(M6.RHI); length = 5))
    bin(r) = min(4, searchsortedlast(edges, r))
    rng = MersenneTwister(TPICKSEED)
    pick = Int[]; nb = zeros(Int, 8, 4)
    for s in eachindex(M6.STRATA)
        cs = [c for c in cand if strat[c[1]] == s]
        got = Int[]
        for b in 1:4
            cb = [c[1] for c in cs if bin(c[2]) == b]
            g = cb[randperm(rng, length(cb))[1:min(BNBIN, length(cb))]]
            nb[s, b] = length(g); append!(got, g)
        end
        rest = [c[1] for c in cs if !(c[1] in got)]
        need = BNPER - length(got)
        need > 0 && append!(got, rest[randperm(rng, length(rest))[1:min(need, length(rest))]])
        append!(pick, got)
    end
    pick, pool, nb
end

function ts_sample(ids, strat, ex)
    rng = MersenneTwister(TSSEED)
    S = Int[]
    for s in eachindex(M6.STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        append!(S, el[randperm(rng, length(el))[1:TNS]])
    end
    sort(S)
end

tk_job(j) = (o = simulate(kernel_params(j.p); seed = j.seed, init = :equilibrium, burnin = KBURN); (; alive = !o.extinct, extinct_gen = o.extinct_gen))

# one P input: the four emulators in turn, each timed (wall time of the call in its thread)
function tp_job(id, p)
    sd = job_seed(TBSEED, id)
    t2 = @elapsed e2 = emulate(p; loci = :two_stage, band = TBAND, screen = TSCREEN, seed = sd)
    tb = @elapsed eb = emulate(p; loci = :blocks, blocks = BSTAR, within = BWITHIN, replicates = TRB, seed = sd)
    tl = @elapsed el = emulate(p; loci = :spectrum, spectrum = SPEC)
    ts = @elapsed es = emulate(p; loci = :spectrum)
    nan(x) = x === missing ? NaN : Float64(x)
    (; P2 = e2.P2000, T2 = nan(e2.T), stage = e2.stage, inband = e2.inband, sP = e2.screen_P2000, sT = e2.screen_T,
     probe = e2.screen_probe, probeP = e2.screen_probe_P2000, sfb = e2.screen_fallback, t2,
     Pb = eb.P2000, Tb = nan(eb.T), tb, Pl = el.P2000, Tl = nan(el.T), tl, Ps = es.P2000, Ts = nan(es.T), ts)
end

# S: the two-stage emulator on every input (timed); on a stratified subsample of 32 (4 per
# stratum) also the block model alone and ld = :dynamic alone, timed in the same job
function ts_job(id, p, sub::Bool)
    sd = job_seed(TSBSEED, id)
    t2 = @elapsed e2 = emulate(p; loci = :two_stage, band = TBAND, screen = TSCREEN, seed = sd)
    nan(x) = x === missing ? NaN : Float64(x)
    tb = sub ? (@elapsed emulate(p; loci = :blocks, blocks = BSTAR, within = BWITHIN, replicates = TRB, seed = sd)) : NaN
    tl = sub ? (@elapsed emulate(p; loci = :spectrum, spectrum = SPEC)) : NaN
    (; P2000 = e2.P2000, T = nan(e2.T), stage = e2.stage, inband = e2.inband, sP = e2.screen_P2000, sT = e2.screen_T,
     sfb = e2.screen_fallback, t2, tb, tl)
end

# --- E1–E4 with the well-defined rules of validate/two_stage_criteria.md
function criteria2(P, Tp, ix, tab, rows, ps, strat; label)
    Tobs = [tab.T[rows[i]] for i in ix]
    alive = Tobs .> 2000
    st = strat[ix]
    res = Dict{Symbol,Bool}()
    lsay("  n = $(length(ix)); observed alive at 2000: $(count(alive)) ($(round(100mean(alive), digits = 2))%); mean P2000 = " *
         "$(round(mean(P), sigdigits = 4)); P2000 exactly 0: $(count(==(0), P)), exactly 1: $(count(==(1), P)), in (0.01, 0.99): " *
         "$(count(x -> 0.01 < x < 0.99, P))")
    a = M6.auc(P, alive)
    as = [M6.auc(P[st.==s], alive[st.==s]) for s in eachindex(M6.STRATA)]
    nd = count(!isnan, as)
    res[:E1] = a >= 0.85 && nd >= 6 && all(x -> isnan(x) || x >= 0.80, as)
    lsayf("  E1 discrimination: AUC(P2000; alive at 2000) = %.4f (≥ 0.85); per stratum ≥ 0.80 where defined, %d of 8 defined (≥ 6)  %s",
          a, nd, label(res[:E1]))
    for s in eachindex(M6.STRATA)
        lsayf("       %-18s n = %4d  alive %4d  AUC = %s", M6.STRATA[s][1], count(st .== s), count(alive[st.==s]),
              isnan(as[s]) ? "undefined (one outcome only; excluded)" : @sprintf("%.4f", as[s]))
    end
    o = sortperm(P; alg = Base.Sort.DEFAULT_STABLE)
    n = length(P); worst = 0.0
    lsay("  E2 calibration: deciles of P2000 (ties broken by input order)")
    for d in 1:10
        k = o[(div((d - 1) * n, 10) + 1):div(d * n, 10)]
        mp, ob = mean(P[k]), mean(alive[k])
        worst = max(worst, abs(mp - ob))
        lsayf("       %6d %6d   %10.4f   %8.4f   %6.4f", d, length(k), mp, ob, abs(mp - ob))
    end
    β, conv = M6.logistic([ones(n) M6.logit.(clamp.(P, 1e-6, 1 - 1e-6))], Float64.(alive))
    slope_ok = conv ? 0.8 <= β[2] <= 1.25 : true
    res[:E2] = worst <= 0.10 && slope_ok
    lsayf("     largest |observed − mean P2000| = %.4f (≤ 0.10); calibration slope = %.4f, intercept %.4f %s  %s", worst, β[2], β[1],
          conv ? "(converged; in [0.8, 1.25])" : "(NOT converged: separation, slope undefined; E2 scored on the deciles alone)",
          label(res[:E2]))
    k = findall(i -> Tobs[i] <= 2000 && !isnan(Tp[i]), eachindex(Tobs))
    if length(k) >= 3
        ρ = M6.spearman(Tp[k], Float64.(Tobs[k])); mr = median(Tp[k] ./ Tobs[k])
        res[:E3] = ρ >= 0.7 && 0.8 <= mr <= 1.25
        lsayf("  E3 timing: %d of %d inputs extinct within 2000 in the table have a predicted T; Spearman = %.4f (≥ 0.7), " *
              "median T_pred/T_obs = %.4f (IQR %.3f–%.3f; in [0.8, 1.25])  %s", length(k), count(<=(2000), Tobs), ρ, mr,
              quantile(Tp[k] ./ Tobs[k], 0.25), quantile(Tp[k] ./ Tobs[k], 0.75), label(res[:E3]))
    else
        res[:E3] = false
        lsay("  E3 timing: fewer than 3 inputs with both an observed and a predicted T  $(label(false))")
    end
    lsay("  E4 x50 = d/√V_P per stratum by E6′d's rule (x50_rule; plain logistic x50 in brackets, no verdict)")
    e4 = true; ndef = 0
    for s in eachindex(M6.STRATA)
        kk = findall(==(s), st)
        rate = [ps[ix[i]].MoveFitnessSD / sqrt(ps[ix[i]].PhenotypeVar) for i in kk]
        xt, ht = M6.x50_rule(rate, Float64.(alive[kk])); xe, he = M6.x50_rule(rate, P[kk])
        pt, pe = M6.x50(rate, alive[kk])[1], M6.x50(rate, P[kk])[1]
        if !isfinite(xt)
            lsayf("       %-18s table x50 undefined (%s): stratum excluded   [plain: table %.5f, emulator %.5f]", M6.STRATA[s][1], ht, pt, pe)
            continue
        end
        ndef += 1
        ok = isfinite(xe) && 0.8 <= xe / xt <= 1.25
        e4 &= ok
        lsayf("       %-18s table %.5f (%s)  emulator %.5f (%s)  ratio %.3f%s   [plain: %.5f, %.5f]", M6.STRATA[s][1], xt, ht, xe, he,
              xe / xt, ok ? "" : "  outside", pt, pe)
    end
    res[:E4] = e4 && ndef >= 6
    lsayf("     E4 (every defined stratum in [0.8, 1.25], %d of 8 defined, ≥ 6): %s", ndef, label(res[:E4]))
    res
end

function tmain()
    t0 = time()
    ids, ps, strat = inv_(() -> M6.read_inputs())
    tab = inv_(() -> M6.load_master(M6.MASTER))
    find = inv_(() -> M6.row_index(tab))
    rows = [find(M6.STRATA[strat[i]][2], M6.STRATA[strat[i]][3], ps[i]) for i in eachindex(ps)]
    ex, nv = all_exclusions(ids, ps, strat, tab, rows)
    n0 = length(ex)
    sc, src, skipped = t_scanned_ids()
    union!(ex, sc)
    TMODE == "select" || (LLOG[] = open(joinpath(HERE, "two_stage.log"), "w"); M6.LOG[] = LLOG[])
    lsay("validate/two_stage.jl — the two-stage emulator (loci = :two_stage, band = $TBAND, screen = $TSCREEN), under validate/two_stage_criteria.md")
    lsay("$(Threads.nthreads()) threads; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    lsay("excluded inputs: $(length(ex)) = spectrum_ld's exclusions $n0 (of which $nv from vignettes 09/11) ∪ $(length(sc)) ids found in files of the worktrees' validate/ directories")
    lsay("  id occurrences found (summed over files), by directory: " * join(["$k $v" for (k, v) in sort(collect(src))], "; "))
    lsay("  whole-table files not used for exclusion: " * join(skipped, "; "))
    P, pool, nb = tp_inputs(ids, ps, strat, ex)
    S = ts_sample(ids, strat, union(ex, Set(ids[pool])))
    sub = S[1:div(length(S), 32):end][1:32]
    lsay("P: $(length(P)) inputs (pool $(length(pool)), MersenneTwister($TPOOLSEED); picks MersenneTwister($TPICKSEED)); per stratum and bin " *
         join([join(nb[s, :], "/") for s in 1:8], ", ") * "; per stratum " * join([count(i -> strat[i] == s, P) for s in 1:8], ", "))
    lsay("S: $(length(S)) inputs ($TNS per stratum, MersenneTwister($TSSEED)); cost subsample $(length(sub)) (" *
         join([count(i -> strat[i] == s, sub) for s in 1:8], "/") * " per stratum)")
    lsay("disjointness (set intersection with the exclusions): P ∩ ex = $(length(intersect(Set(ids[P]), ex))), pool ∩ ex = " *
         "$(length(intersect(Set(ids[pool]), ex))), S ∩ ex = $(length(intersect(Set(ids[S]), ex))), S ∩ pool = $(length(intersect(S, pool)))")
    blkraw = normpath(joinpath(ROOTWT, "wt-blocks-rep", "validate", "blocks_raw"))
    if isdir(blkraw)
        prev = Set{String}()
        for f in readdir(blkraw; join = true), m in eachmatch(IDRX, read(f, String))
            push!(prev, normid(m.match))
        end
        lsay("  check: ids in blocks.jl's caches (H, F, pool): $(length(prev)); ∩ P = $(length(intersect(prev, Set(ids[P])))), " *
             "∩ pool = $(length(intersect(prev, Set(ids[pool])))), ∩ S = $(length(intersect(prev, Set(ids[S]))))")
    end
    TMODE == "select" && return

    # --- P
    lsay("\n== P: fresh inputs near the threshold, $TRK fresh kernel runs each ==")
    jobs = [(; key = "$(ids[i])#$r", p = ps[i], seed = TKSEED + 100j + r) for (j, i) in enumerate(P) for r in 1:TRK]
    t = @elapsed kr = sweep(tk_job, jobs; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 400, cache = tcache("two_stage_kernel.csv"))
    lsayf("  kernel: %d runs, %.0f s wall this session", length(kr), t)
    s = [count(r -> r.alive, kr[(j-1)*TRK+1:j*TRK]) for j in eachindex(P)]
    lsayf("  kernel fraction alive at 2000: mean %.4f; inputs with 0 alive %d, all alive %d, in between %d",
          mean(s ./ TRK), count(==(0), s), count(==(TRK), s), count(x -> 0 < x < TRK, s))
    t = @elapsed E = sweep(i -> tp_job(ids[i], ps[i]), P; id = i -> "P|$(ids[i])", cost = i -> ps[i].NoQTL, progress = 20,
                           cache = tcache("two_stage_P.csv"))
    lsayf("  emulations on P (two-stage, blocks alone, ld = :dynamic, default :spectrum, in turn per input): %.0f s wall this session", t)
    rate = [ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in P]
    st = strat[P]
    f = s ./ TRK; M = f .> 0.5
    vl(b) = pf(b); base(b) = "($(pf(b)) as a comparator)"
    xk, _ = inv_(() -> M6.x50_rule(rate, f))
    pooled(Pv) = inv_(() -> M6.x50_rule(rate, Pv))[1] / xk
    P2 = [e.P2 for e in E]; Pb = [e.Pb for e in E]; Pl = [e.Pl for e in E]; Ps = [e.Ps for e in E]
    inb = [e.stage == 2 for e in E]
    # implementation check: in the band the two-stage result is the block model's (same seed)
    nbad = count(j -> inb[j] && P2[j] != Pb[j], eachindex(E))
    lsayf("  stage 2 (in the band): %d of %d (%.4f); per stratum %s; implementation check, inputs in the band whose two-stage P2000 ≠ the block model's: %d (must be 0)",
          count(inb), length(inb), mean(inb), join([count(inb[st.==k]) for k in 1:8], "/"), nbad)
    lsayf("  screen calls (P2000 > 1/2 at k) vs kernel majority: %d of %d agree", count((e.sP > 0.5) == m for (e, m) in zip(E, M)), length(E))
    lsay("\n  -- model under test: two-stage --")
    r2 = inv_(() -> M6.e6pstats(P2, s, rate, st; label = vl))
    rp2 = pooled(P2)
    for (name, Pv) in (("block model alone", Pb), ("ld = :dynamic alone (grid)", Pl), ("default :spectrum", Ps))
        lsay("\n  -- comparator: $name (no verdict) --")
        inv_(() -> M6.e6pstats(Pv, s, rate, st; label = base))
        lsayf("    pooled x50 ratio %.3f", pooled(Pv))
    end
    # T1
    acc(Pv) = mean((Pv .> 0.5) .== M)
    auc(Pv) = M6.auc(Pv, M)
    c2 = (P2 .> 0.5) .== M; cb = (Pb .> 0.5) .== M
    b_only, two_only = count(cb .& .!c2), count(c2 .& .!cb)
    pmc = b_only + two_only == 0 ? 1.0 : min(1.0, 2 * Distributions.cdf(Binomial(b_only + two_only, 0.5), min(b_only, two_only)))
    n1, n0 = count(M), count(.!M)
    tol = max(1 / max(n1, 1), 1 / max(n0, 1))
    accok = (b_only - two_only <= 1) || pmc > 0.05
    aucok = auc(P2) >= auc(Pb) - tol
    T1 = accok && aucok
    lsay("\n  (T1) threshold calls against the kernel's majority")
    lsayf("    accuracy: two-stage %.4f, blocks alone %.4f, ld = :dynamic alone %.4f, default :spectrum %.4f", acc(P2), acc(Pb), acc(Pl), acc(Ps))
    lsayf("    AUC:      two-stage %.4f, blocks alone %.4f, ld = :dynamic alone %.4f, default :spectrum %.4f", auc(P2), auc(Pb), auc(Pl), auc(Ps))
    lsayf("    discordant with blocks alone: blocks right & two-stage wrong %d, two-stage right & blocks wrong %d (exact McNemar p %.3g; direction: %s)",
          b_only, two_only, pmc, b_only > two_only ? "blocks better" : b_only < two_only ? "two-stage better" : "tie")
    lsayf("    accuracy part (difference ≤ 1 input, or McNemar p > 0.05): %s; AUC part (two-stage ≥ blocks − max(1/n₁, 1/n₀) = %.4f − %.4f): %s",
          pf(accok), auc(Pb), tol, pf(aucok))
    lsay("    (T1) $(pf(T1))")
    for (name, Pv) in (("ld = :dynamic alone", Pl), ("default :spectrum", Ps))
        cc = (Pv .> 0.5) .== M
        a_, b_ = count(c2 .& .!cc), count(cc .& .!c2)
        lsayf("    vs %s (no verdict): two-stage right & it wrong %d, it right & two-stage wrong %d", name, a_, b_)
    end
    # T2
    lsayf("\n  (T2) pooled x50 ratio two-stage/kernel %.3f (in [0.8, 1.25]); per stratum (E6′d above): %s", rp2, pf(r2[:d]))
    T2 = inband(rp2, 0.8, 1.25) && r2[:d]
    lsay("    (T2) $(pf(T2))")
    # T3
    t2 = [e.t2 for e in E]; tb = [e.tb for e in E]; tl = [e.tl for e in E]; ts_ = [e.ts for e in E]
    T3 = median(t2) <= median(tb) / 3
    lsay("\n  (T3) cost per input (wall time of each call in its thread, 8 threads, the four calls in turn in one job)")
    lsayf("    median: two-stage %.2f s (stage 1 only %.2f s, stage 2 %.2f s), blocks alone %.2f s, ld = :dynamic alone %.2f s, default :spectrum %.2f s",
          median(t2), any(.!inb) ? median(t2[.!inb]) : NaN, any(inb) ? median(t2[inb]) : NaN, median(tb), median(tl), median(ts_))
    lsayf("    mean:   two-stage %.2f s, blocks alone %.2f s, ld = :dynamic alone %.2f s, default :spectrum %.2f s", mean(t2), mean(tb), mean(tl), mean(ts_))
    lsayf("    median two-stage / median blocks alone = %.3f (≤ 1/3): %s", median(t2) / median(tb), pf(T3))

    # --- S
    lsay("\n== S: $(length(S)) fresh table inputs, the table's outcomes as truth ==")
    subset = Set(sub)
    t = @elapsed ES = sweep(i -> ts_job(ids[i], ps[i], i in subset), S; id = i -> "S|$(ids[i])", cost = i -> ps[i].NoQTL, progress = 200,
                            cache = tcache("two_stage_S.csv"))
    lsayf("  emulations on S: %.0f s wall this session", t)
    stS = strat[S]
    inS = [e.stage == 2 for e in ES]
    fS = mean(inS); ci = 1.96sqrt(fS * (1 - fS) / length(S))
    lsayf("  stage-2 fraction: %d of %d = %.4f (95%% CI %.4f–%.4f); per stratum %s", count(inS), length(S), fS, max(fS - ci, 0), fS + ci,
          join([@sprintf("%.3f", mean(inS[stS.==k])) for k in 1:8], ", "))
    lsay("\n  -- model under test: two-stage --")
    rS = criteria2([e.P2000 for e in ES], [e.T for e in ES], S, tab, rows, ps, strat; label = vl)
    lsay("\n  -- comparator: the screen alone ($(TSCREEN), its P2000 and T at k; no verdict) --")
    criteria2([e.sP for e in ES], [e.sT for e in ES], S, tab, rows, ps, strat; label = base)
    # cost
    ks = [j for (j, i) in enumerate(S) if i in subset]
    t2S = [e.t2 for e in ES]
    lsay("\n  cost (wall time per call in its thread on 8 threads)")
    lsayf("    two-stage on all S: mean %.2f s, median %.2f s (stage 1 mean %.2f s, stage 2 mean %.2f s)", mean(t2S), median(t2S),
          mean(t2S[.!inS]), any(inS) ? mean(t2S[inS]) : NaN)
    sb = [ES[j].tb for j in ks]; sl = [ES[j].tl for j in ks]; s2 = t2S[ks]
    lsayf("    subsample of %d: mean two-stage %.2f s, blocks alone %.2f s, ld = :dynamic alone %.2f s", length(ks), mean(s2), mean(sb), mean(sl))
    N = 119_984
    lsayf("    extrapolated to the %d inputs (CPU-hours, threaded wall time per call summed): two-stage %.1f (from all S), blocks alone %.1f, ld = :dynamic alone %.1f (from the subsample)",
          N, mean(t2S) * N / 3600, mean(sb) * N / 3600, mean(sl) * N / 3600)
    lsayf("    ratios: two-stage/blocks alone %.3f, two-stage/ld = :dynamic alone %.3f; wall time on 8 threads for the table: two-stage ≈ %.1f h",
          mean(t2S) / mean(sb), mean(t2S) / mean(sl), mean(t2S) * N / 3600 / 8)

    lsay("\nVERDICTS (loci = :two_stage, band $TBAND, screen $TSCREEN): primary (T1) $(pf(T1)), (T2) $(pf(T2)), (T3) $(pf(T3)); " *
         "secondary on S: stage-2 fraction $(round(fS, digits = 4)); " * join(["$k $(pf(rS[k]))" for k in (:E1, :E2, :E3, :E4)], ", "))
    lsayf("total %.0f min wall this session", (time() - t0) / 60)
    LLOG[] === stdout || close(LLOG[])
end

tmain()
