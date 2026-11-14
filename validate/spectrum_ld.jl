# The :spectrum closure with dynamic linkage disequilibrium (`spectrum = (; ld = :dynamic)`,
# src/spectrum.jl) under the pre-specified criteria of validate/spectrum_ld_criteria.md
# (committed before this was run). Writes validate/spectrum_ld.log.
#
#   H:   ~160 held-out published inputs near the threshold (E6′'s procedure on a new pool),
#        40 fresh kernel runs each, with the kernel's V_G and genic variance before the move
#        and V_G at +100 and +500: tests (a), (b), (c).
#   C:   the diagnosis's confirmation set (13 inputs; their E6′ kernel runs and the diagnosis's
#        10 trajectories per input), reported separately.
#   F:   a fresh table sample (150 per stratum): E1–E4 of validate/emulator_criteria.md (d).
#   Baselines (no verdict): the default :spectrum, the `ldratio` prototype, :rare_alleles.
#
# Run from the package root (vignettes environment: XLSX, Catlab):
#   julia -t 8 --project=vignettes validate/spectrum_ld.jl select   # the samples only (incl. the pool's k50s)
#   julia -t 8 --project=vignettes validate/spectrum_ld.jl pilot    # timing of the emulations on 16 inputs of F
#   julia -t 8 --project=vignettes validate/spectrum_ld.jl          # everything; resumes from the caches
# Caches (resumable `sweep` CSVs): validate/spectrum_ld_pool.csv, spectrum_ld_kernel.csv,
# spectrum_ld_runs.csv. The diagnosis's raw trajectories (for C's (a) and (b)) are read from
# ENV RAWDIR (default ../wt-spectrum-diagnosis/validate/spectrum_diagnosis_raw).

ENV["RAWDIR"] = get(ENV, "RAWDIR", normpath(joinpath(@__DIR__, "..", "..", "wt-spectrum-diagnosis", "validate",
                                                     "spectrum_diagnosis_raw")))
include(joinpath(@__DIR__, "spectrum_diagnosis.jl"))      # e6p_module, readp, load_runs, genic_stats, spectrum_traj, ld_ratio, …
using Distributions: Binomial, logpdf

const M6 = e6p_module()                                   # validate/emulator_spectrum.jl's functions
const LMODE = isempty(ARGS) ? "full" : ARGS[1]
const SPEC = (; ld = :dynamic)                            # the model under test
const POOLSEED, PICKSEED, NPOOL2 = 20261121, 20261122, 100
const FSEED, NF = 20261123, 100          # amendment 1: 100 per stratum (the first 100 of the same permutation)
const HK, HKSEED = 40, 70261114000
const LLOG = Ref{IO}(stdout)
lsay(s...) = (println(stdout, s...); LLOG[] === stdout || (println(LLOG[], s...); flush(LLOG[])))
lsayf(fmt, args...) = lsay(Printf.format(Printf.Format(fmt), args...))
lcache(name) = joinpath(HERE, name)
const LRUNS = lcache("spectrum_ld_runs.csv")
inv_(f) = Base.invokelatest(f)
pf(b) = b ? "PASS" : "FAIL"

# --- exclusions ---------------------------------------------------------------------------
function all_exclusions(ids, ps, strat, tab, rows)
    ex = inv_(() -> M6.excluded_ids())
    S = inv_(() -> M6.main_sample(ids, strat, ex))
    union!(ex, ids[S])
    pool = inv_(() -> M6.e6p_pool(ids, ps, strat, union(ex, Set(ids[S]))))
    union!(ex, ids[pool])                                  # E6′'s pool, which contains E6′'s inputs
    union!(ex, readlines(joinpath(HERE, "spectrum_diagnosis_fresh_inputs.txt")))
    for r in read_inputs_csv()
        push!(ex, r.id)
    end
    for l in readlines(joinpath(HERE, "published_sweep_results.csv"))[2:end]
        f = split(l, ","); push!(ex, "$(f[1])/$(f[2])")
    end
    # vignettes 09 and 11 (table rows): the inputs whose table row was rerun
    vrows = Set{Int}()
    vig = joinpath(HERE, "..", "vignettes")
    for f in [joinpath(vig, "09_replicating_pagel_2025", "replication_runs.csv");
              [joinpath(vig, "11_demography_in_the_papers_populations", "demography_runs_R$R.csv") for R in ("1.5", "2.0", "5.0")]]
        isfile(f) || continue
        for l in readlines(f)[2:end]
            x = tryparse(Int, split(l, ",")[1]); x === nothing || push!(vrows, x)
        end
    end
    nv = 0
    for i in eachindex(ids)
        rows[i] in vrows && (push!(ex, ids[i]); nv += 1)
    end
    ex, nv
end

# H: E6′'s procedure on a new pool
function h_inputs(ids, ps, strat, ex)
    rng = MersenneTwister(POOLSEED)
    pool = Int[]
    for s in eachindex(M6.STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        got = 0
        for i in el[randperm(rng, length(el))]
            kc = drift_critical_rate(ps[i]).k
            (kc > 0 && 0.3 <= ps[i].MoveFitnessSD / kc <= 3.0) || continue
            push!(pool, i); got += 1
            got == NPOOL2 && break
        end
    end
    t = @elapsed k50s = sweep(i -> (; k50 = emulated_half_survival_rate(ps[i]; burnin = BURNIN, horizon = HORIZON).k50), pool;
                              id = i -> ids[i], cost = i -> ps[i].NoQTL, progress = 200, cache = lcache("spectrum_ld_pool.csv"))
    lsayf("H pool: k50 (:rare_alleles) of %d inputs, %.0f s wall this session", length(pool), t)
    cand = Tuple{Int,Float64}[]
    for (i, k) in zip(pool, k50s)
        0 < k.k50 < Inf || continue
        r = ps[i].MoveFitnessSD / k.k50
        M6.RLO <= r <= M6.RHI && push!(cand, (i, r))
    end
    sort!(cand; by = c -> ids[c[1]])
    edges = exp.(range(log(M6.RLO), log(M6.RHI); length = 5))
    bin(r) = min(4, searchsortedlast(edges, r))
    rng = MersenneTwister(PICKSEED)
    pick = Int[]; nb = zeros(Int, 8, 4)
    for s in eachindex(M6.STRATA)
        cs = [c for c in cand if strat[c[1]] == s]
        got = Int[]
        for b in 1:4
            cb = [c[1] for c in cs if bin(c[2]) == b]
            g = cb[randperm(rng, length(cb))[1:min(M6.NBIN_PER, length(cb))]]
            nb[s, b] = length(g); append!(got, g)
        end
        rest = [c[1] for c in cs if !(c[1] in got)]
        need = M6.NE6PER - length(got)
        need > 0 && append!(got, rest[randperm(rng, length(rest))[1:min(need, length(rest))]])
        append!(pick, got)
    end
    pick, pool, nb
end

# F: a fresh table sample, NF per stratum, outside the exclusions and H's pool
function f_sample(ids, strat, ex)
    rng = MersenneTwister(FSEED)
    S = Int[]
    for s in eachindex(M6.STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        append!(S, el[randperm(rng, length(el))[1:NF]])
    end
    sort(S)
end

# --- the kernel: E6′'s runs (simulate's draws) with V_G each generation and the genic variance
# at MOVE − 200, − 100, − 1 (genic_stats draws nothing)
function h_kernel_job(j)
    q = kernel_params(j.p)
    rng = MersenneTwister(j.seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    vpre = 0.0; npre = 0; gpre = 0.0; ngp = 0
    v100 = 0.0; n100 = 0; v500 = 0.0; n500 = 0
    extinct_gen = -1
    for it in w.gen:q.NoGenerations-1
        w.gen = it
        if it in (MOVE - 200, MOVE - 100, MOVE - 1)
            gpre += genic_stats(w; seg = false)[1]; ngp += 1
        end
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        t = it - MOVE
        if -200 <= t <= -1 || 95 <= t <= 105 || 495 <= t <= 505
            v = var(view(w.genotype, 1:w.N))
            -200 <= t <= -1 ? (vpre += v; npre += 1) : 95 <= t <= 105 ? (v100 += v; n100 += 1) : (v500 += v; n500 += 1)
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    (; alive = extinct_gen < 0, extinct_gen, VGpre = vpre / max(npre, 1), genicpre = gpre / max(ngp, 1),
     VG100 = n100 == 11 ? v100 / 11 : NaN, VG500 = n500 == 11 ? v500 / 11 : NaN)
end

# --- the emulator, with the trajectory: P2000 (emulate's pipeline), V_G at the move (mean over
# −200…−1, as the kernel's), at +100 and +500 (mean over ±5), and the equilibrium's genic and LD
function em_job(p, kw; eq = true)
    st = spectrum_traj(p; spectrum = kw)
    m = st.m
    if length(m.gen) <= BURNIN                             # died before the move
        P, T = 0.0, Float64(m.extinct_gen - BURNIN + 1)
    else
        d = extinction_distribution(st.q, m; from = BURNIN)
        P, T = d.psurvive, d.median === missing ? NaN : Float64(d.median - BURNIN + 1)
    end
    w(lo, hi) = (k = findall(t -> lo <= t <= hi, st.t); length(k) == hi - lo + 1 ? mean(st.VG[k]) : NaN)
    e = eq && haskey(kw, :ld) && kw.ld === :dynamic ? spectrum_equilibrium(p; kw...) : nothing
    (; P2000 = P, T, VGpre = w(-200, -1), VG100 = w(95, 105), VG500 = w(495, 505),
     genic_eq = e === nothing ? NaN : e.genic, ld_eq = e === nothing ? NaN : e.ld)
end

function em_sweep(model, kwf, ix, ids, ps; eq = true)
    sweep(collect(ix); id = i -> "$model|$(ids[i])", cost = i -> ps[i].NoQTL, progress = 100, cache = LRUNS) do i
        em_job(ps[i], kwf(ps[i]); eq)
    end
end
ra_sweep(tag, ix, ids, ps) = emulate_sweep(ps[ix]; ids = ["$tag|" * x for x in ids[ix]], loci = :rare_alleles,
                                           burnin = BURNIN, horizon = HORIZON, tail = M6.TAIL, progress = 200, cache = lcache("spectrum_ld_runs_ra.csv"))

# --- (a) and (b) ----------------------------------------------------------------------------
med(v) = (u = filter(isfinite, v); isempty(u) ? NaN : median(u))
qq(v, q) = (u = filter(isfinite, v); isempty(u) ? NaN : quantile(u, q))
inband(x, lo, hi) = isfinite(x) && lo <= x <= hi
function ab_tests(label, K, E; verdict = true)
    # K: per input (VGpre, genicpre, VG100, VG100n, VG500, VG500n); E: em_job results (model under test)
    rv = [E[i].VGpre / K[i].VGpre for i in eachindex(K)]
    rg = [E[i].genic_eq / K[i].genicpre for i in eachindex(K)]
    kr = [(K[i].VGpre - K[i].genicpre) / K[i].genicpre for i in eachindex(K)]
    mr = [E[i].ld_eq / E[i].genic_eq for i in eachindex(K)]
    sel = findall(x -> x <= -0.1, kr)
    rl = [mr[i] / kr[i] for i in sel]
    a1 = inband(med(rv), 0.9, 1.1); a2 = inband(med(rg), 0.9, 1.1); a3 = inband(med(rl), 0.8, 1.25)
    v(b) = verdict ? pf(b) : "-"
    lsayf("    (a1) V_G at the move, model/kernel: median %.3f (IQR %.3f–%.3f; n %d; in [0.9, 1.1]) %s",
          med(rv), qq(rv, 0.25), qq(rv, 0.75), count(isfinite, rv), v(a1))
    lsayf("    (a2) genic variance at the move, model (equilibrium)/kernel: median %.3f (IQR %.3f–%.3f; in [0.9, 1.1]) %s",
          med(rg), qq(rg, 0.25), qq(rg, 0.75), v(a2))
    lsayf("    (a3) LD/genic, model/kernel, on the %d inputs with kernel LD/genic ≤ −0.1: median %.3f (IQR %.3f–%.3f; in [0.8, 1.25]) %s",
          length(sel), med(rl), qq(rl, 0.25), qq(rl, 0.75), v(a3))
    b = Bool[]
    for (t, f, fn) in ((100, :VG100, :n100), (500, :VG500, :n500))
        ok = [i for i in eachindex(K) if getfield(K[i], fn) >= 10 && isfinite(getfield(E[i], f))]
        r = [getfield(E[i], f) / getfield(K[i], f) for i in ok]
        bb = inband(med(r), 0.8, 1.25)
        push!(b, bb)
        lsayf("    (b) V_G at +%d, model/kernel (kernel mean over ≥ 10 runs alive; model alive): median %.3f (IQR %s; n %d; in [0.8, 1.25]) %s",
              t, med(r), isempty(r) ? "-" : @sprintf("%.3f–%.3f", quantile(r, 0.25), quantile(r, 0.75)), length(r), v(bb))
    end
    (; a = a1 && a2 && a3, a1, a2, a3, b = all(b), b100 = b[1], b500 = b[2])
end

# kernel summaries per input from the H kernel runs
function k_summary(kr, nin)
    map(1:nin) do j
        rs = kr[(j-1)*HK+1:j*HK]
        f(x) = (u = filter(isfinite, x); isempty(u) ? NaN : mean(u))
        (; s = count(r -> r.alive, rs), VGpre = f([r.VGpre for r in rs]), genicpre = f([r.genicpre for r in rs]),
         VG100 = f([r.VG100 for r in rs]), n100 = count(r -> isfinite(r.VG100), rs),
         VG500 = f([r.VG500 for r in rs]), n500 = count(r -> isfinite(r.VG500), rs))
    end
end

# C: the diagnosis's 10 trajectories per input for (a), (b), the E6′ 40 runs for (c)
function c_summary(rows)
    kr, _ = load_e6p()
    map(rows) do r
        runs = load_runs(r.id, 10)
        win(x, f, xf, lo, hi) = (v = [getfield(x, f)[k] for k in eachindex(getfield(x, xf)) if lo <= getfield(x, xf)[k] - MOVE <= hi]; v)
        f(v) = isempty(v) ? NaN : mean(v)
        vpre = f([mean(win(x, :VG, :gen, -200, -1)) for x in runs])
        gpre = f([mean(win(x, :genic, :xg, -200, -10)) for x in runs])
        a100 = [mean(win(x, :VG, :gen, 95, 105)) for x in runs if length(win(x, :VG, :gen, 95, 105)) == 11]
        a500 = [mean(win(x, :VG, :gen, 495, 505)) for x in runs if length(win(x, :VG, :gen, 495, 505)) == 11]
        (; s = count(x -> x[2], kr[r.id]), VGpre = vpre, genicpre = gpre, VG100 = f(a100), n100 = 4length(a100),
         VG500 = f(a500), n500 = 4length(a500))       # n scaled to 40 runs: "≥ 10 of 40" is ≥ 3 of 10
    end
end

function main()
    t0 = time()
    ids, ps, strat = inv_(() -> M6.read_inputs())
    tab = inv_(() -> M6.load_master(M6.MASTER))
    find = inv_(() -> M6.row_index(tab))
    rows = [find(M6.STRATA[strat[i]][2], M6.STRATA[strat[i]][3], ps[i]) for i in eachindex(ps)]
    ex, nv = all_exclusions(ids, ps, strat, tab, rows)
    LMODE in ("select", "pilot") || (LLOG[] = open(joinpath(HERE, "spectrum_ld.log"), "w"); M6.LOG[] = LLOG[])
    lsay("validate/spectrum_ld.jl — the :spectrum closure with dynamic LD, under validate/spectrum_ld_criteria.md")
    lsay("burn-in $BURNIN, horizon $HORIZON; $(Threads.nthreads()) threads; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    lsay("excluded inputs: $(length(ex)) (of which $nv from vignettes 09/11's table rows)")
    H, pool, nb = h_inputs(ids, ps, strat, ex)
    lsay("H: $(length(H)) inputs (pool $(length(pool)), MersenneTwister($POOLSEED); picks MersenneTwister($PICKSEED)); per stratum and bin " *
         join([join(nb[s, :], "/") for s in 1:8], ", ") * "; per stratum " * join([count(i -> strat[i] == s, H) for s in 1:8], ", "))
    F = f_sample(ids, strat, union(ex, Set(ids[pool])))
    lsay("F: $(length(F)) inputs ($NF per stratum, MersenneTwister($FSEED)); overlap with H's pool $(length(intersect(F, pool)))")
    LMODE == "select" && return
    if LMODE == "pilot"
        P = F[1:div(length(F), 16):end][1:16]
        for (name, kw) in (("ld", SPEC), ("spectrum", (;)))
            t = @elapsed em_sweep("pilot-$name", _ -> kw, P, ids, ps)
            lsayf("pilot %s: 16 inputs in %.0f s wall on %d threads", name, t, Threads.nthreads())
        end
        return
    end

    # --- H: kernel and emulations
    lsay("\n== H: held-out inputs near the threshold, $(HK) fresh kernel runs each ==")
    jobs = [(; key = "$(ids[i])#$r", p = ps[i], seed = HKSEED + 100j + r) for (j, i) in enumerate(H) for r in 1:HK]
    t = @elapsed kr = sweep(h_kernel_job, jobs; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 400,
                            cache = lcache("spectrum_ld_kernel.csv"))
    lsayf("  kernel: %d runs, %.0f s wall this session", length(kr), t)
    K = k_summary(kr, length(H))
    s = [k.s for k in K]
    lsayf("  kernel fraction alive at 2000: mean %.4f; inputs with 0 alive %d, all alive %d, in between %d",
          mean(s ./ HK), count(==(0), s), count(==(HK), s), count(x -> 0 < x < HK, s))
    rate = [ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in H]
    st = strat[H]
    t = @elapsed begin
        Eld = em_sweep("ld", _ -> SPEC, H, ids, ps)
        Esp = em_sweep("spectrum", _ -> (;), H, ids, ps)
        Elr = em_sweep("ldratio", p -> (; ldratio = ld_ratio(p)), H, ids, ps)
        Era = ra_sweep("rare_alleles", H, ids, ps)
    end
    lsayf("  emulations on H: %.0f s wall this session", t)
    vl(b) = pf(b); base(b) = "($(pf(b)) as a baseline)"
    lsay("\n  -- model under test: :spectrum, ld = :dynamic --")
    rab = ab_tests("H", K, Eld)
    lsay("    (c) E6′'s statistics:")
    r6 = inv_(() -> M6.e6pstats([e.P2000 for e in Eld], s, rate, st; label = vl))
    xk, _ = inv_(() -> M6.x50_rule(rate, s ./ HK)); xe, _ = inv_(() -> M6.x50_rule(rate, [e.P2000 for e in Eld]))
    rc_pool = inband(xe / xk, 0.8, 1.25)
    lsayf("    (c) pooled x50 ratio %.3f (in [0.8, 1.25]) %s", xe / xk, pf(rc_pool))
    for (name, E) in (("default :spectrum", Esp), ("ldratio prototype", Elr))
        lsay("\n  -- baseline: $name (no verdict) --")
        ab_tests("H", K, E; verdict = false)
        inv_(() -> M6.e6pstats([e.P2000 for e in E], s, rate, st; label = base))
    end
    lsay("\n  -- baseline: :rare_alleles, deterministic (no verdict) --")
    inv_(() -> M6.e6pstats([e.P2000 for e in Era], s, rate, st; label = base))

    # --- C
    lsay("\n== C: the diagnosis's confirmation set (13 inputs; kernel runs already seen; reported separately) ==")
    crow = filter(r -> r.set == "C", read_inputs_csv())
    ci = [findfirst(==(r.id), ids) for r in crow]
    KC = c_summary(crow)
    sC = [k.s for k in KC]
    rateC = [ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in ci]; stC = strat[ci]
    for (name, kwf) in (("ld", _ -> SPEC), ("spectrum", _ -> (;)), ("ldratio", p -> (; ldratio = ld_ratio(p))))
        E = em_sweep("C-$name", kwf, ci, ids, ps)
        lsay("  -- $name --")
        ab_tests("C", KC, E; verdict = false)
        P = [e.P2000 for e in E]
        M = sC ./ RK .> 0.5; Cc = P .> 0.5
        lsayf("    threshold calls: accuracy %.3f (emulator alive & kernel dead %d, emulator dead & kernel alive %d)",
              mean(M .== Cc), count(Cc .& .!M), count(.!Cc .& M))
    end

    # --- F: E1–E4 on the table
    lsay("\n== (d) F: E1–E4 of validate/emulator_criteria.md on $(length(F)) fresh table inputs ==")
    alltab = map(1:8) do s
        k = findall(==(s), strat)
        inv_(() -> M6.x50([ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in k], [tab.T[rows[i]] > 2000 for i in k]))[1]
    end
    t = @elapsed begin
        Fld = em_sweep("F-ld", _ -> SPEC, F, ids, ps; eq = false)
        Fsp = em_sweep("F-spectrum", _ -> (;), F, ids, ps; eq = false)
        Fra = ra_sweep("F-rare_alleles", F, ids, ps)
    end
    lsayf("  emulations on F: %.0f s wall this session", t)
    asem(E) = [(; P2000 = e.P2000, P240k = NaN, T = e.T) for e in E]
    lsay("\n  -- model under test: ld = :dynamic --")
    rd = inv_(() -> M6.criteria(asem(Fld), F, tab, rows, ps, strat; label = vl, alltab))
    lsay("\n  -- baseline: default :spectrum --")
    inv_(() -> M6.criteria(asem(Fsp), F, tab, rows, ps, strat; label = base, alltab))
    lsay("\n  -- baseline: :rare_alleles --")
    inv_(() -> M6.criteria(Fra, F, tab, rows, ps, strat; label = base, alltab))

    lsay("\nVERDICTS (ld = :dynamic): (a) $(pf(rab.a)) (a1 $(pf(rab.a1)), a2 $(pf(rab.a2)), a3 $(pf(rab.a3))); " *
         "(b) $(pf(rab.b)) (+100 $(pf(rab.b100)), +500 $(pf(rab.b500))); (c) E6′c $(pf(r6[:c])), x50 per stratum $(pf(r6[:d])), " *
         "pooled x50 $(pf(rc_pool)); (d) " * join(["$k $(pf(rd[k]))" for k in (:E1, :E2, :E3, :E4)], ", "))
    lsayf("total %.0f min wall this session", (time() - t0) / 60)
    LLOG[] === stdout || close(LLOG[])
end

main()
