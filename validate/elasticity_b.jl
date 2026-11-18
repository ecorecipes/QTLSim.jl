# Phase B of the elasticity study: the emulator's elasticities of the half-survival rate on a
# target that does not saturate (`k50_gradient`, src/k50_gradient.jl), re-checked against the
# published table under the criteria of phase A (validate/elasticity_criteria.md, G1–G3 and
# tolerances unchanged) as extended by its "Phase B addendum" (committed before this script
# computed any elasticity of a sampled input). Phase A's script, caches and log are not touched;
# the table side (cells, terciles, sample) is phase A's.
#
# Run from the package root (vignettes environment):
#   julia -t 6 --project=vignettes validate/elasticity_b.jl run       # the emulator (resumable caches)
#   julia --project=vignettes validate/elasticity_b.jl report         # comparison and verdicts → validate/elasticity_b.log
# Caches (git-ignored): validate/elasticity_raw/b_*.csv.

const BMODE = isempty(ARGS) ? "report" : ARGS[1]
let                                             # phase A's helpers without running any of its modes
    empty!(ARGS); push!(ARGS, "none")
    include(joinpath(@__DIR__, "elasticity.jl"))
end

const BGRID = (; ld = :dynamic)                 # the screen's closure on the frequency grid
const BMINROUNDS, BCPULIMIT = 12, 3.0 * 3600    # rounds always run; CPU seconds after which no further round starts
const BGRIDCAP, BGRIDCHECK = 20, 3              # fallback-switch inputs recomputed with the grid; ok inputs with a fallback, as a check
const BSTATUS = (:ok, :no_k50, :no_root, :not_finite, :not_decreasing, :fallback_switch, :max_evals)
cpu_seconds() = ccall(:clock, Clong, ()) / 1e6  # process CPU time (user + system, all threads; macOS CLOCKS_PER_SEC)

"One input: `k50_gradient` as a flat row for the cache."
function b_job(p::Params; kw...)
    t = @elapsed g = k50_gradient(p; names = GIN, kw...)
    e = collect(values(g.elasticities)); s = collect(values(g.steps)); f = g.fallback
    (; k50 = g.k50, evals = g.evals, status = findfirst(==(g.status), BSTATUS) - 1, step = g.step, h0 = g.h0, hk = g.dh_dlogk,
     e1 = e[1], e2 = e[2], e3 = e[3], e4 = e[4], s1 = s[1], s2 = s[2], s3 = s[3], s4 = s[4], asym = g.asymmetry,
     fb0 = f.base, fbmin = minimum(minimum, values(f)[2:end]), fbmax = maximum(maximum, values(f)[2:end]),
     straddle = g.straddle, shift = g.shift, continued = g.continued, seconds = t)
end

function run_mode()
    t0 = time(); c0 = cpu_seconds()
    S = g_read_sample()
    ps = Dict(id => g_params(id) for id in S.id)
    say("phase B run: $(Threads.nthreads()) threads; QTLSim at ", pathof(QTLSim))
    occursin("wt-ad", pathof(QTLSim)) || error("not the worktree's QTLSim")
    let p = ps[S.id[1]]                              # compile
        for sp in (GSCREEN, BGRID)
            QTLSim.emulated_log_hazard(Params(p; MoveFitnessSD = 1e-3); spectrum = sp, burnin = 20, horizon = 20)
        end
    end
    A, rare = g_cache("screen.csv"), g_cache("rare.csv")
    # the bracket: phase A's k50 where it has one (a search result, not an elasticity), else from :rare_alleles as phase A did
    function bracket(id)
        haskey(A, id) && 0 < A[id][:k50] < Inf && return (A[id][:k50] * exp(-0.006), A[id][:k50] * exp(0.006))
        haskey(rare, id) && 0 < rare[id][:k50] < Inf && return (rare[id][:k50] / 4, rare[id][:k50])
        nothing
    end
    job(id) = b_job(ps[id]; bracket = bracket(id), onswitch = :flag)
    nround = 0
    for r in 1:GNPER
        (r > BMINROUNDS && cpu_seconds() - c0 > BCPULIMIT) && break
        ids = r == 1 ? S.id[S.order.<=BMINROUNDS] : r <= BMINROUNDS ? String[] : S.id[S.order.==r]
        isempty(ids) || sweep(job, ids; cache = gcache("b_screen.csv"), progress = 16)
        nround = r
        (r == 1 || r >= BMINROUNDS) && sayf("  screen through round %2d: CPU %.0f s, wall %.0f s", max(r, BMINROUNDS), cpu_seconds() - c0, time() - t0)
        r == 1 && (nround = BMINROUNDS)
    end
    nround = max(nround, BMINROUNDS)
    open(io -> println(io, nround), gcache("b_rounds.txt"), "w")
    run = [i for i in eachindex(S.id) if S.order[i] <= nround]
    scr = sweep(id -> error("not cached"), S.id[run]; cache = gcache("b_screen.csv"))
    k = Dict(S.id[i] => scr[a].k50 for (a, i) in enumerate(run))
    # fallback switches that smaller steps did not avoid: the grid solver at its own k50
    sw = [S.id[i] for (a, i) in enumerate(run) if scr[a].status == 5]
    ck = [S.id[i] for (a, i) in enumerate(run) if scr[a].status == 0 && scr[a].fb0 >= 0]
    todo = [sw[1:min(end, BGRIDCAP)]; ck[1:min(end, BGRIDCHECK)]]
    grid(id) = b_job(ps[id]; spectrum = BGRID, bracket = (k[id] * exp(-0.02), k[id] * exp(0.02)))
    rg = sweep(grid, todo; cache = gcache("b_grid.csv"))
    sayf("  fallback switches: %d inputs (recomputed with the grid: %d); grid check on %d inputs: CPU %.0f s", length(sw), min(length(sw), BGRIDCAP),
         min(length(ck), BGRIDCHECK), cpu_seconds() - c0)
    # N1: a four times smaller span on round 1
    r1 = [id for id in S.id[S.order.==1] if 0 < k[id] < Inf]
    sweep(id -> b_job(ps[id]; bracket = (k[id] * exp(-1e-4), k[id] * exp(1e-4)), dh = 0.125, onswitch = :flag), r1; cache = gcache("b_step.csv"))
    sayf("phase B run done: CPU %.0f s (%.2f h), wall %.0f s", cpu_seconds() - c0, (cpu_seconds() - c0) / 3600, time() - t0)
    open(io -> println(io, cpu_seconds() - c0, " ", time() - t0, " ", Threads.nthreads()), gcache("b_cpu.txt"), "a")
end

BMODE == "run" && run_mode()

b_e(r) = r[:status] == 0 ? [r[:e1], r[:e2], r[:e3], r[:e4]] : fill(NaN, 4)
b_fmt(v) = join([@sprintf("% .3f", x) for x in v], " ")

function report_mode()
    LOG[] = open(joinpath(@__DIR__, "elasticity_b.log"), "w")
    say("validate/elasticity_b.jl report — phase B: the emulator's elasticities of the half-survival rate on the non-saturating target")
    say("criteria: validate/elasticity_criteria.md (G1–G3 unchanged) with its Phase B addendum; table side: validate/elasticity_table.log; ",
        Libc.strftime("%Y-%m-%d %H:%M", time()))
    S = g_read_sample(); cells = g_read_cells()
    nround = parse(Int, strip(read(gcache("b_rounds.txt"), String)))
    run = [i for i in eachindex(S.id) if S.order[i] <= nround]
    old = [i for i in eachindex(S.id) if S.order[i] <= 8]
    scr, grd, stp = g_cache("b_screen.csv"), g_cache("b_grid.csv"), g_cache("b_step.csv")
    A, rare, brt = g_cache("screen.csv"), g_cache("rare.csv"), g_cache("screen_brute.csv")
    R = [scr[S.id[i]] for i in run]
    VE = S.VP .* (1 .- S.h2)

    say("\n== the runs ==")
    say("screen (two-stage screen, k50_gradient): $(nround) rounds of one input per stratum = $(length(run)) inputs; the first 8 rounds are phase A's 64 inputs")
    tsec(d) = sum(r[:seconds] for r in values(d); init = 0.0)
    sayf("summed job seconds (wall inside threads, overstating CPU): screen %.0f, grid %.0f, step check %.0f", tsec(scr), tsec(grd), tsec(stp))
    if isfile(gcache("b_cpu.txt"))
        for l in readlines(gcache("b_cpu.txt"))
            f = split(l)
            sayf("process CPU time of a `run` invocation (clock(): user + system over all threads): %.0f s = %.2f h; wall %.0f s; %s threads",
                 parse(Float64, f[1]), parse(Float64, f[1]) / 3600, parse(Float64, f[2]), f[3])
        end
    end
    st = [Int(r[:status]) for r in R]
    say("  status of the screen as run (beta solver): ", join(["$(BSTATUS[c+1]) $(count(==(c), st))" for c in 0:6], ", "))
    sayf("  evaluations per input: median %d, max %d; step in log k: median %.2g, range %.2g to %.2g", median([r[:evals] for r in R]),
         maximum(r[:evals] for r in R), median(fin([r[:step] for r in R])), minimum(fin([r[:step] for r in R])), maximum(fin([r[:step] for r in R])))
    ok = st .== 0
    sayf("  at k50 (ok inputs): |h| median %.4f, max %.4f; ∂h/∂log k median %.0f (range %.0f to %.0f); asymmetry median %.3f, 90th pct %.3f, max %.3f",
         median(abs(r[:h0]) for r in R[ok]), maximum(abs(r[:h0]) for r in R[ok]), median(r[:hk] for r in R[ok]), minimum(r[:hk] for r in R[ok]),
         maximum(r[:hk] for r in R[ok]), median(r[:asym] for r in R[ok]), quantile([r[:asym] for r in R[ok]], 0.9), maximum(r[:asym] for r in R[ok]))
    nc = count(r -> r[:continued] == 1, R)
    sayf("  inputs whose target relied on the continuation of a trajectory that died inside the horizon (base or a difference point): %d of %d", nc, length(R))
    lowω = [a for (a, i) in enumerate(run) if S.omega[i] < 2.2]
    sayf("    of them with ω < 2.2: %d of %d; with ω ≥ 2.2: %d of %d", count(a -> R[a][:continued] == 1, lowω), length(lowω),
         nc - count(a -> R[a][:continued] == 1, lowω), length(R) - length(lowω))

    say("\n== fallback of the beta solver ==")
    say("  inputs whose run at k50 fell back: $(count(r -> r[:fb0] >= 0, R)) of $(length(R)); status fallback_switch (a pair and the base in different ",
        "fallback generations at every step tried): $(count(==(5), st))")
    Eas = fill(NaN, length(run), 4); Eco = fill(NaN, length(run), 4); Era = fill(NaN, length(run), 4)
    nre = 0; nlost = 0
    for (a, i) in enumerate(run)
        id = S.id[i]; r = R[a]
        Eas[a, :] .= b_e(r)
        haskey(rare, id) && (Era[a, :] .= g_e(rare[id]))
        if r[:status] == 5
            haskey(grd, id) && grd[id][:status] == 0 ? (Eco[a, :] .= b_e(grd[id]); nre += 1) : (nlost += 1)
        else
            Eco[a, :] .= Eas[a, :]
        end
    end
    say("  fallback-switch inputs recomputed with the grid and used in the consistent values: $nre; left undefined (beyond the cap of $BGRIDCAP, or the grid undefined): $nlost")
    for (a, i) in enumerate(run)
        id = S.id[i]; haskey(grd, id) || continue
        g = grd[id]; r = R[a]
        if r[:status] == 5
            sayf("    %-34s switch   grid: status %s  e %s   k50 grid/beta %.4f", id, BSTATUS[Int(g[:status])+1], b_fmt(b_e(g)), g[:k50] / r[:k50])
        else
            sayf("    %-34s check    beta e %s   grid e %s   max |diff| %.3f   k50 grid/beta %.4f", id, b_fmt(b_e(r)), b_fmt(b_e(g)),
                 maximum(abs.(b_e(r) .- b_e(g))), g[:k50] / r[:k50])
        end
    end

    say("\n== undefined inputs (consistent values) ==")
    und = [a for a in eachindex(run) if !all(isfinite, Eco[a, :])]
    say("  defined $(length(run) - length(und)) of $(length(run)); undefined $(length(und))")
    for a in und
        i = run[a]
        sayf("    %-34s ω %.3g  status %s  k50 %.4g  step %.2g", S.id[i], S.omega[i], BSTATUS[Int(R[a][:status])+1], R[a][:k50], R[a][:step])
    end
    sayf("  inputs with ω < 2.2: %d, defined %d (phase A: 17 of its 64 inputs, 5 defined)", length(lowω), count(a -> all(isfinite, Eco[a, :]), lowω))
    oa = [a for (a, i) in enumerate(run) if S.order[i] <= 8]
    wasund = [a for a in oa if A[S.id[run[a]]][:status] != 0]
    sayf("  phase A's %d undefined inputs: now defined %d", length(wasund), count(a -> all(isfinite, Eco[a, :]), wasund))
    for a in wasund
        i = run[a]
        sayf("    %-34s ω %.3g  k50 %.5g (phase A %.5g)  ∂h/∂log k %.0f  step %.2g  continued %d  e %s", S.id[i], S.omega[i], R[a][:k50], A[S.id[i]][:k50],
             R[a][:hk], R[a][:step], R[a][:continued], b_fmt(Eco[a, :]))
    end

    say("\n== numerical checks ==")
    d1 = Float64[]
    for id in S.id[S.order.==1]
        haskey(stp, id) || continue
        a, b = b_e(scr[id]), b_e(stp[id])
        append!(d1, abs.(a .- b))
        sayf("  N1 %-34s span ±0.5 (step %.2g): %s   span ±0.125 (step %.2g): %s", id, scr[id][:step], b_fmt(a), stp[id][:step], b_fmt(b))
    end
    nd1 = length(d1); d1 = fin(d1)
    n1 = !isempty(d1) && median(d1) <= 0.02 && quantile(d1, 0.9) <= 0.05
    isempty(d1) || sayf("  N1 step: |e(±0.5) − e(±0.125)| over %d values (%d undefined at either span): median %.4f (≤ 0.02), 90th percentile %.4f (≤ 0.05), max %.4f: %s",
                        length(d1), nd1 - length(d1), median(d1), quantile(d1, 0.9), maximum(d1), n1 ? "stable" : "NOT STABLE")
    d2 = Float64[]; nb = 0; miss = 0
    for id in sort(collect(keys(brt)))
        a = Eco[findfirst(i -> S.id[i] == id, run), :]; b = [brt[id][:b1], brt[id][:b2], brt[id][:b3], brt[id][:b4]]
        for j in 1:4
            isfinite(b[j]) || continue
            nb += 1
            isfinite(a[j]) ? push!(d2, abs(a[j] - b[j])) : (miss += 1)
        end
        sayf("  N2 %-34s implicit (new target): %s   brute force (phase A): %s", id, b_fmt(a), b_fmt(b))
    end
    n2 = nb > 0 && miss == 0 && median(d2) <= 0.03 && maximum(d2) <= 0.10
    sayf("  N2 brute force: %d values defined by brute force, %d of them with an undefined implicit value; |brute − implicit| median %.4f (≤ 0.03), max %.4f (≤ 0.10): %s",
         nb, miss, isempty(d2) ? NaN : median(d2), isempty(d2) ? NaN : maximum(d2), n2 ? "verified" : "NOT VERIFIED")
    say("  numerical status of the verdicts below: ", n1 && n2 ? "reliable (N1 and N2 pass)" : "NUMERICALLY UNRELIABLE (N1 or N2 failed)")

    say("\n== primary comparison: 48 cells (locus number × input × information tercile of the input) ==")
    g_gstats("the screen on the new target, consistent values (grid where a fallback switch remained), all $(length(run)) inputs", Eco, S, run, cells; verdict = true)
    g_gstats("the screen as run (beta values; fallback-switch inputs undefined)", Eas, S, run, cells; verdict = false)
    g_gstats("phase A's 64 inputs alone on the new target (consistent values)", Eco[oa, :], S, run[oa], cells; verdict = false)
    g_gstats(":rare_alleles, phase A's values (old target) on the inputs run", Era, S, run, cells; verdict = false)

    say("\n== the low-ω tercile cells, one by one (the cells phase A could not test) ==")
    for c in cells
        (c.j == 4 && c.t == 1) || continue
        ins = [a for (a, i) in enumerate(run) if (S.strat[i] + 1) ÷ 2 == c.li && g_tercile(cells, c.li, 4, S.omega[i]) == 1]
        v = fin(Eco[ins, 4])
        tol = max(0.15, 0.30 * abs(c.e))
        sayf("  L %4d omega t1 (ω ≤ %.3g): table % .3f [% .3f, % .3f]; emulator n %d defined of %d: median % .3f, mean % .3f, quartiles % .3f % .3f; sign %s, size %s",
             c.L, c.hi_θ, c.e, c.lo, c.hi, length(v), length(ins), med(v), avg(v), isempty(v) ? NaN : quantile(v, 0.25), isempty(v) ? NaN : quantile(v, 0.75),
             length(v) < 4 ? "—" : sign(median(v)) == sign(c.e) ? "ok" : "WRONG", length(v) < 4 ? "—" : c.lo - tol <= median(v) <= c.hi + tol ? "ok" : "FAIL")
        for a in ins
            i = run[a]
            sayf("      %-34s ω %.3f  k50 %.4g  e_VP % .3f e_h2 % .3f e_m % .3f e_ω % .3f  %s", S.id[i], S.omega[i], R[a][:k50], Eco[a, 1], Eco[a, 2], Eco[a, 3],
                 Eco[a, 4], all(isfinite, Eco[a, :]) ? (R[a][:continued] == 1 ? "(continued)" : "") : "UNDEFINED: " * string(BSTATUS[Int(R[a][:status])+1]))
        end
    end

    say("\n== C1 (descriptive): the new target against phase A's values on its defined inputs ==")
    Eold = fill(NaN, length(oa), 4)
    for (b, a) in enumerate(oa)
        Eold[b, :] .= g_e(A[S.id[run[a]]])
    end
    both = [b for b in eachindex(oa) if all(isfinite, Eold[b, :]) && all(isfinite, Eco[oa[b], :])]
    sayf("  phase A defined %d inputs; of them defined now %d", count(b -> all(isfinite, Eold[b, :]), eachindex(oa)), length(both))
    D = abs.(Eco[oa[both], :] .- Eold[both, :])
    sayf("  |new − old| over %d values: median %.4f, 90th pct %.4f, max %.4f; values within 0.02: %d of %d; inputs with all four within 0.02: %d of %d",
         length(D), median(D), quantile(vec(D), 0.9), maximum(D), count(<=(0.02), D), length(D), count(r -> maximum(r) <= 0.02, eachrow(D)), length(both))
    shifted = [b for b in both if A[S.id[run[oa[b]]]][:sgen] > 0]
    plain = setdiff(both, shifted)
    for (nm, set) in (("phase A pairs in different fallback generations (\"shift\", $(length(shifted)) inputs)", shifted), ("the others ($(length(plain)) inputs)", plain))
        isempty(set) && continue
        Dx = abs.(Eco[oa[set], :] .- Eold[set, :])
        sayf("    %s: median %.4f, max %.4f, values within 0.02: %d of %d", nm, median(Dx), maximum(Dx), count(<=(0.02), Dx), length(Dx))
    end
    for b in both
        maximum(abs.(Eco[oa[b], :] .- Eold[b, :])) > 0.02 || continue
        id = S.id[run[oa[b]]]
        sayf("    > 0.02: %-34s ω %.3g  old %s (step %.0e, shift pairs %d)   new %s (step %.1e)   k50 new/old %.5f", id, S.omega[run[oa[b]]], b_fmt(Eold[b, :]),
             A[id][:step], A[id][:sgen], b_fmt(Eco[oa[b], :]), R[oa[b]][:step], R[oa[b]][:k50] / A[id][:k50])
    end

    say("\n== C2 (descriptive): e_V_P against e_m, and the identity in h² ==")
    def = [a for a in eachindex(run) if all(isfinite, Eas[a, :]) && R[a][:status] == 0]
    i_ = run[def]
    x = S.h2[i_] ./ (1 .- S.h2[i_])
    idd = abs.(Eas[def, 2] .+ x .* Eas[def, 1])
    sayf("  identity e_h² = −h²/(1 − h²)·e_V_P (beta values, %d inputs): |difference| median %.2g, 90th pct %.2g, max %.2g; relative to |e_h²|: median %.2g, max %.2g",
         length(def), median(idd), quantile(idd, 0.9), maximum(idd), median(idd ./ abs.(Eas[def, 2])), maximum(idd ./ abs.(Eas[def, 2])))
    d = Eas[def, 1] .- Eas[def, 3]                     # e_V_P − e_m = the elasticity in V_E at fixed m·V_E and ω
    u = VE[i_] ./ (S.omega[i_] .^ 2 .+ VE[i_])
    pred = Eas[def, 4] .* VE[i_] ./ (2 .* S.omega[i_] .^ 2)   # if ω and V_E acted only through V_s = ω² + V_E
    sayf("  e_V_P − e_m: median % .4f, range % .4f to % .4f; |e_V_P − e_m| ≤ 0.01 in %d of %d inputs", median(d), minimum(d), maximum(d), count(<=(0.01), abs.(d)), length(d))
    for (nm, sel) in (("ω ≥ 2.2", S.omega[i_] .>= 2.2), ("1 ≤ ω < 2.2", (S.omega[i_] .< 2.2) .& (S.omega[i_] .>= 1)), ("ω < 1", S.omega[i_] .< 1))
        any(sel) || continue
        sayf("    %-12s n %3d  V_E/(ω² + V_E) median %.2g   e_V_P − e_m: median % .4f, min % .4f, max % .4f   predicted e_ω·V_E/(2ω²): median % .4f; |observed − predicted| median %.4f, max %.4f",
             nm, count(sel), median(u[sel]), median(d[sel]), minimum(d[sel]), maximum(d[sel]), median(pred[sel]), median(abs.(d[sel] .- pred[sel])),
             maximum(abs.(d[sel] .- pred[sel])))
    end
    sayf("  Spearman of |e_V_P − e_m| with V_E/(ω² + V_E): %.3f; of (e_V_P − e_m) with the prediction e_ω·V_E/(2ω²): %.3f", spearman(abs.(d), u), spearman(d, pred))

    say("\n== S1: the literal per-stratum reading (32 cells; no verdict; distorted by the non-constancy, see the criteria) ==")
    T, _ = readdlm(gcache("table_elasticities.tsv"), '\t'; header = true)
    let te = Float64[], em = Float64[], g1 = Bool[], g2 = Bool[]
        for s in eachindex(STRATA)
            ins = [a for (a, i) in enumerate(run) if S.strat[i] == s]
            line = ""
            for j in 1:4
                r = findfirst(q -> Int(T[q, 1]) == s && T[q, 2] == GLAB[j], axes(T, 1))
                e, lo, hi = Float64(T[r, 3]), Float64(T[r, 5]), Float64(T[r, 6])
                v = fin(Eco[ins, j]); isempty(v) && continue
                m = median(v); push!(te, e); push!(em, m)
                (lo > 0 || hi < 0) && push!(g1, sign(m) == sign(e))
                tol = max(0.15, 0.3abs(e))
                abs(e) >= 0.1 && push!(g2, lo - tol <= m <= hi + tol)
                line *= @sprintf("  %s table % .3f emu median % .3f mean % .3f (n %d) |", GLAB[j], e, m, mean(v), length(v))
            end
            say("    ", rpad(STRATA[s][1], 18), line)
        end
        sayf("    for information: sign %d/%d, size %d/%d, Spearman %.3f over %d cells", count(g1), length(g1), count(g2), length(g2), spearman(te, em), length(te))
    end

    say("\n== S2: input by input against the table's quadratic-model local elasticity (no verdict; the quadratic is poor for ω) ==")
    for j in 1:4
        k = [a for a in eachindex(run) if isfinite(Eco[a, j])]
        q = S.q[run[k], j]
        sayf("  e_%-6s n %3d  Spearman %.3f  median |emu − table| %.3f  median (emu − table) % .3f   emu quartiles % .3f % .3f % .3f   table % .3f % .3f % .3f",
             GLAB[j], length(k), spearman(Eco[k, j], q), median(abs.(Eco[k, j] .- q)), median(Eco[k, j] .- q),
             quantile(Eco[k, j], 0.25), median(Eco[k, j]), quantile(Eco[k, j], 0.75), quantile(q, 0.25), median(q), quantile(q, 0.75))
    end

    say("\n== S4: levels (no verdict) ==")
    say("  k50 of the screen (new target) / the table's boundary rate at the same inputs (quadratic fit): median [quartiles] by locus number")
    say("    ", join(map(1:4) do li
        v = [R[a][:k50] / exp(S.lkq[i]) for (a, i) in enumerate(run) if (S.strat[i] + 1) ÷ 2 == li && 0 < R[a][:k50] < Inf]
        @sprintf("L %d: %.2f [%.2f, %.2f]", GLS[li], median(v), quantile(v, 0.25), quantile(v, 0.75))
    end, "   "))
    kr = [R[a][:k50] / A[S.id[run[a]]][:k50] for a in oa if 0 < R[a][:k50] < Inf && 0 < A[S.id[run[a]]][:k50] < Inf]
    sayf("  k50 on the new target / phase A's k50 (to 1e-2 in log k) on its %d inputs: median %.5f, range %.5f to %.5f", length(kr), median(kr), minimum(kr), maximum(kr))
    say("\n== post hoc (added after the results were seen; not part of the criteria; no verdict changed) ==")
    ph = g_cache("b_posthoc.csv")
    if !isempty(ph)
        say("  N2's one large difference, e_ω of QT5000-Possion5.0/Input-11353.txt (ω = 0.372; implicit ", @sprintf("%.3f", Eco[findfirst(i -> S.id[i] == "QT5000-Possion5.0/Input-11353.txt", run), 4]),
            ", phase A's brute force over ω·e^{±0.05}: ", @sprintf("%.3f", brt["QT5000-Possion5.0/Input-11353.txt"][:b4]), "). Brute force on the new target, k50 re-solved to |h| ≤ 1e-3:")
        for s in sort(collect(keys(ph)); by = x -> -parse(Float64, x))
            r = ph[s]
            sayf("    ω·e^{±%s}: Δlog k50/Δlog ω = % .3f   (local implicit e_ω at ω·e^{+%s}: % .3f, at ω·e^{−%s}: % .3f)", s, r[:brute], s, r[:eplus], s, r[:eminus])
        end
    end
    nr = g_cache("b_posthoc_noroot.csv")
    for id in sort(collect(keys(nr)))
        r = nr[id]
        sayf("  no_root %-34s at k50·e^{−1e-9}: h % .3f (fallback %d, trajectory died in the horizon %d, continued %d); at k50·e^{+1e-9}: h % .3f (fallback %d, died %d, continued %d)",
             id, r[:hlo], r[:fblo], r[:extlo], r[:contlo], r[:hhi], r[:fbhi], r[:exthi], r[:conthi])
    end
    bad = [a for a in def if abs(Eas[a, 2] + S.h2[run[a]] / (1 - S.h2[run[a]]) * Eas[a, 1]) > 0.005]
    say("  C2, inputs where the identity in h² is off by more than 0.005 (", length(bad), "):")
    for a in bad
        i = run[a]
        sayf("    %-34s h² %.3f ω %.3g  e_V_P % .4f  e_h² % .4f  −h²/(1−h²)·e_V_P % .4f  steps V_P %.2g h² %.2g (log k %.2g)", S.id[i], S.h2[i], S.omega[i], Eas[a, 1], Eas[a, 2],
             -S.h2[i] / (1 - S.h2[i]) * Eas[a, 1], R[a][:s1], R[a][:s2], R[a][:step])
    end
    close(LOG[])
end

# Post hoc (added after the results were seen; no verdict, changes none): N2's one large difference,
# e_ω of QT5000-Possion5.0/Input-11353 (ω = 0.37). Brute force on the new target with smaller steps.
function posthoc_mode()
    id = "QT5000-Possion5.0/Input-11353.txt"
    p = g_params(id)
    k0 = g_cache("b_screen.csv")[id][:k50]
    rows = sweep([0.05, 0.02, 0.005]; cache = gcache("b_posthoc.csv")) do s
        ks = map((s, -s)) do d
            k50_gradient(Params(p; InitFitnessSD = p.InitFitnessSD * exp(d)); names = (:InitFitnessSD,), htol = 1e-3,
                         bracket = (k0 * exp(-0.1), k0 * exp(0.1)), onswitch = :flag)
        end
        (; kplus = ks[1].k50, kminus = ks[2].k50, brute = (log(ks[1].k50) - log(ks[2].k50)) / 2s, eplus = ks[1].elasticities[1],
         eminus = ks[2].elasticities[1], evals = ks[1].evals + ks[2].evals)
    end
    foreach(println, rows)
    # the inputs with status no_root: the target on either side of the closed bracket
    S = g_read_sample(); scr = g_cache("b_screen.csv")
    nr = [x for x in S.id if haskey(scr, x) && scr[x][:status] == 2]
    r2 = sweep(nr; cache = gcache("b_posthoc_noroot.csv")) do x
        q = g_params(x); k = scr[x][:k50]
        lo = QTLSim.emulated_log_hazard(Params(q; MoveFitnessSD = k * exp(-1e-9)))
        hi = QTLSim.emulated_log_hazard(Params(q; MoveFitnessSD = k * exp(1e-9)))
        (; hlo = lo.h, hhi = hi.h, fblo = lo.fallback, fbhi = hi.fallback, extlo = lo.extinct, exthi = hi.extinct, contlo = lo.continued, conthi = hi.continued)
    end
    foreach(println, r2)
end
BMODE == "posthoc" && posthoc_mode()

BMODE == "report" && report_mode()
