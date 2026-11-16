# Design of the two-stage emulator (`emulate(p; loci = :two_stage)`, src/two_stage.jl) on
# DEVELOPMENT inputs only: the band rule's width and the choice of the screen's solver.
# Writes validate/two_stage_dev.log.
#
# Development data (disclosed): the inputs of validate/blocks.jl's validation, whose block-model
# results are cached (git-ignored) in ../wt-blocks-rep/validate/blocks_raw/ (ENV BLOCKSRAW):
#   Hd: blocks.jl's H, 120 near-threshold inputs; block model = fraction alive at +2000 among 40
#       run_blocks runs (blocks_traj.csv; 160 blocks, Gaussian law, the kernel's protocol);
#       kernel = 40 kernel runs each (blocks_kernel.csv; printed for information only, not used
#       to set anything);
#   Fd: blocks.jl's F, 320 table inputs (40 per stratum, a random sample of the table);
#       block model = `emulate(…; loci = :blocks)` with 10 replicates (blocks_emulate.csv).
# For every development input and each screen (`:spectrum` with ld = :dynamic, solver :grid or
# :beta; deterministic, burn-in and horizon 2000) the distance d = |log(k/k50_screen)| of the
# input's rate from the screen's own half-survival rate, to ±0.01 up to 0.7 (beyond: d > 0.7),
# by evaluating the screen's call P2000 > 1/2 at k e^{±x}.
#
#   julia -t 8 --project=vignettes validate/two_stage_dev.jl
# Cache (resumable sweep CSV, git-ignored): validate/two_stage_dev_raw/.

using QTLSim, Statistics, Random, Printf
import Catlab
using QTLSim: sweep

const HERE = @__DIR__
const SIMDIR = normpath(joinpath(HERE, "..", "..", "Repository", "Simulations"))
const BRAW = get(ENV, "BLOCKSRAW", normpath(joinpath(HERE, "..", "..", "wt-blocks-rep", "validate", "blocks_raw")))
const RAW = joinpath(HERE, "two_stage_dev_raw")
const WMAX, RES = 0.7, 0.01
const WIDTHS = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7]
const SCREENS = ("grid" => (; ld = :dynamic), "beta" => (; ld = :dynamic, solver = :beta))
const STRATA = ["QT0100-Possion3.0", "QT0100-Possion5.0", "QT0500-Possion3.0", "QT0500-Possion5.0",
                "QT1000-Possion3.0", "QT1000-Possion5.0", "QT5000-Possion3.0", "QT5000-Possion5.0"]
const LOG = Ref{IO}(stdout)
say(s...) = (println(stdout, s...); LOG[] === stdout || (println(LOG[], s...); flush(LOG[])))
sayf(fmt, a...) = say(Printf.format(Printf.Format(fmt), a...))

function readcsv(path)
    l = readlines(path)
    h = split(l[1], ",")
    [Dict(zip(h, split(x, ","))) for x in l[2:end]]
end

# the development inputs and the block model's (and, for Hd, the kernel's) calls
function dev_inputs()
    traj = readcsv(joinpath(BRAW, "blocks_traj.csv"))
    ker = readcsv(joinpath(BRAW, "blocks_kernel.csv"))
    em = readcsv(joinpath(BRAW, "blocks_emulate.csv"))
    frac(rows) = begin
        d = Dict{String,Vector{Bool}}()
        for r in rows
            id = split(r["id"], "#")[1]
            push!(get!(d, id, Bool[]), r["alive"] == "true")
        end
        Dict(k => (mean(v), length(v)) for (k, v) in d)
    end
    tb, tk = frac(traj), frac(ker)
    H = sort(collect(keys(tb)))
    @assert all(tb[i][2] == 40 && tk[i][2] == 40 for i in H)
    F = sort([String(split(r["id"], "|")[2]) for r in em])
    fb = Dict(String(split(r["id"], "|")[2]) => parse(Float64, r["P2000"]) for r in em)
    (; H, F, Pb_H = [tb[i][1] for i in H], Pk_H = [tk[i][1] for i in H], Pb_F = [fb[i] for i in F])
end

# the screen's call and the distance of the input's rate from the screen's k50
function dist_job(p, spec)
    ev = 0
    call(k) = (ev += 1; e = emulate(Params(p; MoveFitnessSD = k); loci = :spectrum, spectrum = spec); (e.P2000 > 0.5, e))
    k = p.MoveFitnessSD
    c0, e0 = call(k)
    fb = get(e0, :fallback, -1)
    if e0.T !== missing && e0.T <= 0           # dies before the move at every rate
        return (; d = Inf, alive = c0, P2000 = e0.P2000, T = Float64(e0.T), evals = ev, fallback = fb)
    end
    s = c0 ? 1.0 : -1.0                        # toward the flip: up if alive, down if dead
    cw, _ = call(k * exp(s * WMAX))
    cw == c0 && return (; d = Inf, alive = c0, P2000 = e0.P2000, T = e0.T === missing ? NaN : Float64(e0.T), evals = ev, fallback = fb)
    lo, hi = 0.0, WMAX                         # call(k e^{s lo}) == c0, call(k e^{s hi}) != c0
    while hi - lo > 2RES
        mid = (lo + hi) / 2
        call(k * exp(s * mid))[1] == c0 ? (lo = mid) : (hi = mid)
    end
    (; d = (lo + hi) / 2, alive = c0, P2000 = e0.P2000, T = e0.T === missing ? NaN : Float64(e0.T), evals = ev, fallback = fb)
end

function main()
    mkpath(RAW)
    LOG[] = open(joinpath(HERE, "two_stage_dev.log"), "w")
    say("validate/two_stage_dev.jl — design of the two-stage emulator on development inputs only")
    say("$(Threads.nthreads()) threads; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    D = dev_inputs()
    ids = [D.H; D.F]
    ps = Dict(i => read_qtsim_options(joinpath(SIMDIR, i)) for i in ids)
    say("development inputs: Hd $(length(D.H)) (blocks.jl's H: block model 40 run_blocks runs, kernel 40 runs), " *
        "Fd $(length(D.F)) (blocks.jl's F: block model 10 replicates); Hd ∩ Fd = $(length(intersect(D.H, D.F)))")
    say("""
RULES, written before any distance was computed (this header is printed before the sweep):
 - The screen's call is P2000 > 1/2 at the input's rate k; an input is inside the band of half-width w
   iff d = |log(k/k50_screen)| < w, i.e. iff the screen's call at k·e^{±w} (toward the flip) differs
   from its call at k (the implemented rule, QTLSim.two_stage_screen, decides this with one more
   evaluation). Here d is resolved to ±$RES up to $WMAX; candidate widths $(WIDTHS).
 - Target: the screen disagrees with the block model (its call P2000 > 1/2) on < 1% of the
   development inputs OUTSIDE the band, separately on Hd (near threshold) and on Fd (random table
   sample). w* = the smallest candidate width meeting both, for each screen.
 - Screen choice: the beta solver is adequate if it meets the target at a width whose stage-2
   fraction on Fd (the table estimate) is at most 5 percentage points above the grid's at the grid's
   own w*, and if the two screens' decisions (inside/outside the band at the chosen width) agree on
   ≥ 95% of Hd ∪ Fd. Then the beta (≈ 3× cheaper per evaluation) is the default screen; otherwise
   the grid. Cost per evaluation is measured below on the development inputs.
 - The chosen width is the beta's (or grid's) w*, unless the stage-2 fraction of Hd at w* exceeds
   1/2 (then T3, cost ≤ 1/3 of the block model on near-threshold inputs, cannot pass; reported, and
   the smallest width meeting the target is still chosen).
 - The kernel's outcomes on Hd are printed for information only; they set nothing.
""")
    res = Dict{String,Any}()
    tim = Dict{String,Float64}()
    for (name, spec) in SCREENS
        jobs = collect(ids)
        t = @elapsed r = sweep(i -> begin
                                   t0 = time(); x = dist_job(ps[i], spec); merge(x, (; secs = time() - t0))
                               end, jobs; id = i -> "$name|$i", cost = i -> ps[i].NoQTL, progress = 50,
                               cache = joinpath(RAW, "two_stage_dev_dist.csv"))
        res[name] = Dict(zip(jobs, r))
        tim[name] = t
        ev = sum(x.evals for x in r); se = sum(x.secs for x in r)
        sayf("screen %s: %d inputs, %d evaluations, %.0f s summed over jobs (%.2f s per evaluation, threaded), %.0f s wall this session",
             name, length(r), ev, se, se / ev, t)
    end
    nH, nF = length(D.H), length(D.F)
    blkH = D.Pb_H .> 0.5; blkF = D.Pb_F .> 0.5; kerH = D.Pk_H .> 0.5
    wstar = Dict{String,Float64}()
    for (name, _) in SCREENS
        R = res[name]
        dH = [R[i].d for i in D.H]; dF = [R[i].d for i in D.F]
        cH = [R[i].alive for i in D.H]; cF = [R[i].alive for i in D.F]
        say("\n== screen $name: calls at k ==")
        sayf("  Hd: screen ≠ blocks %d of %d; screen ≠ kernel majority %d (info); blocks ≠ kernel %d (info)",
             count(cH .!= blkH), nH, count(cH .!= kerH), count(blkH .!= kerH))
        sayf("  Fd: screen ≠ blocks %d of %d", count(cF .!= blkF), nF)
        sayf("  d of the disagreeing inputs: Hd %s; Fd %s", string(round.(sort(dH[cH.!=blkH]), digits = 3)),
             string(round.(sort(dF[cF.!=blkF]), digits = 3)))
        say("  width   Hd in band   Fd in band (table est., 95% CI)   Hd outside: ≠blocks (rate)   Fd outside: ≠blocks (rate)   target   Hd 2-stage acc vs kernel (info)")
        ws = NaN
        for w in WIDTHS
            iH = dH .< w; iF = dF .< w
            oH = .!iH; oF = .!iF
            eH = count(cH[oH] .!= blkH[oH]); eF = count(cF[oF] .!= blkF[oF])
            rH = eH / max(count(oH), 1); rF = eF / max(count(oF), 1)
            ok = rH < 0.01 && rF < 0.01
            isnan(ws) && ok && (ws = w)
            f = mean(iF); ci = 1.96sqrt(f * (1 - f) / nF)
            two = ifelse.(iH, blkH, cH)
            sayf("  %5.2f   %3d (%.3f)   %3d (%.3f, [%.3f, %.3f])            %2d of %3d (%.4f)            %2d of %3d (%.4f)        %s     %.4f",
                 w, count(iH), mean(iH), count(iF), f, max(f - ci, 0), f + ci, eH, count(oH), rH, eF, count(oF), rF, ok ? "met " : "not ", mean(two .== kerH))
        end
        wstar[name] = ws
        sayf("  w* (%s) = %s", name, isnan(ws) ? "none of the candidates" : string(ws))
        sayf("  blocks alone on Hd vs kernel (info): accuracy %.4f; screen alone %.4f", mean(blkH .== kerH), mean(cH .== kerH))
    end
    # the decision with both screens
    say("\n== beta vs grid: the decision (inside/outside the band) ==")
    for w in WIDTHS
        ag = [(res["grid"][i].d < w) == (res["beta"][i].d < w) for i in ids]
        sayf("  w = %.2f: decisions agree on %d of %d (%.4f); calls at k agree on %d", w, count(ag), length(ids), mean(ag),
             count(res["grid"][i].alive == res["beta"][i].alive for i in ids))
    end
    dd = [abs(res["grid"][i].d - res["beta"][i].d) for i in ids if isfinite(res["grid"][i].d) && isfinite(res["beta"][i].d)]
    sayf("  |d_grid − d_beta| where both ≤ %.1f: n %d, median %.3f, 90th pct %.3f, max %.3f", WMAX, length(dd), median(dd), quantile(dd, 0.9), maximum(dd))
    fbk = count(res["beta"][i].fallback >= 0 for i in ids)
    sayf("  beta fell back to the grid (at k) in %d of %d inputs", fbk, length(ids))
    say("\n== decision ==")
    wg, wb = wstar["grid"], wstar["beta"]
    fF(name, w) = mean(res[name][i].d < w for i in D.F)
    fH(name, w) = mean(res[name][i].d < w for i in D.H)
    if isnan(wb)
        say("  beta: no candidate width meets the target: screen = grid")
        choice, w = "grid", wg
    else
        agree = mean((res["grid"][i].d < wb) == (res["beta"][i].d < wb) for i in ids)
        okf = isnan(wg) || fF("beta", wb) <= fF("grid", wg) + 0.05
        sayf("  beta w* %.2f (Fd stage-2 fraction %.4f, Hd %.4f); grid w* %s (Fd %.4f, Hd %.4f); decisions agree at the beta's w* on %.4f",
             wb, fF("beta", wb), fH("beta", wb), string(wg), isnan(wg) ? NaN : fF("grid", wg), isnan(wg) ? NaN : fH("grid", wg), agree)
        choice, w = okf && agree >= 0.95 ? ("beta", wb) : ("grid", wg)
    end
    say("  CHOSEN: screen = $choice, band = $w")
    isnan(w) || fH(choice, w) > 0.5 && say("  note: more than half of Hd is in the band at this width, so T3 cannot pass")
    close(LOG[])
end

main()
