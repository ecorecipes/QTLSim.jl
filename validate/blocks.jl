# The haplotype-block lineage emulator (`emulate(p; loci = :blocks)`, src/block_lineages.jl)
# under the pre-specified criteria of validate/blocks_criteria.md (committed before this was
# run on fresh data). Writes validate/blocks.log.
#
#   H:  fresh published inputs near the deterministic (:rare_alleles) k50, E6′'s procedure on a
#       new pool, 40 fresh kernel runs each (alive at +2000; V_G, genic and block-genic variance
#       before and during the move): the primary threshold tests and the trajectory ratios.
#   F:  a fresh table sample: E1–E4 of validate/emulator_criteria.md.
#   Comparators (same inputs): the default :spectrum, :spectrum with ld = :dynamic, :rare_alleles.
#   Cost: wall time per emulator call on one thread.
#
# Run from the package root (vignettes environment: XLSX, Catlab):
#   julia -t 8 --project=vignettes validate/blocks.jl select     # exclusions and samples only (incl. the pool's k50s)
#   julia -t 8 --project=vignettes validate/blocks.jl pilot      # timing on 16 inputs of F (no outcome)
#   julia -t 8 --project=vignettes validate/blocks.jl            # everything; resumes from the caches
#   julia -t 1 --project=vignettes validate/blocks.jl cost       # cost per call, one thread
# Caches (resumable `sweep` CSVs, git-ignored): validate/blocks_raw/.

let src = read(joinpath(@__DIR__, "spectrum_ld.jl"), String)
    # spectrum_ld.jl's functions (exclusions, E6′'s statistics through M6, em_job), without its main()
    include_string(Main, replace(src, r"\nmain\(\)\s*$" => "\n"), joinpath(@__DIR__, "spectrum_ld.jl"))
end
using QTLSim: job_seed, block_map, block_stats
import Distributions

const BMODE = isempty(ARGS) ? "full" : ARGS[1]
const BSTAR, BWITHIN = 160, :gaussian               # blocks and within-block law (validate/blocks_dev.log: set by convergence)
const BBURN_H, BBURN_F = 1000, 2000                # block emulator burn-in on H (the kernel's protocol) and on F (emulate's default)
const BPOOLSEED, BPICKSEED, BFSEED = 20261131, 20261132, 20261133
const BHKSEED = 41261116000                        # kernel seeds BHKSEED + 100j + r
const BLSEED, BTSEED = 20261134, 20261135          # job_seed bases: emulator replicates, trajectory runs
const BNPOOL, BNF, BRK, BRF = 100, 40, 40, 10    # pool per stratum, F per stratum, runs per H input, replicates per F input
const BNBIN, BNPER = 3, 15                         # H: per bin of r, per stratum
const BRAW = joinpath(@__DIR__, "blocks_raw")
bcache(name) = (mkpath(BRAW); joinpath(BRAW, name))
const TIMES = (-200, -100, -1, 50, 100, 200)
const ROOTWT = normpath(joinpath(@__DIR__, "..", ".."))

# --- exclusions: spectrum_ld's (E6′'s pools and samples, the diagnosis's inputs, the published
# sweep comparison, vignettes 09/11) plus every input id found in any file or raw-run name of any
# worktree's validate/ (the linked-selection chain's D, F, G, H′, J; spectrum_ld's H pool, held-out
# set and F; spectrum_beta's inputs; …), except whole-table files (the deterministic emulator's
# map of all 119,984 inputs, which no kernel comparison of this design used)
const IDRX = r"QT\d{4}-Possion\d\.\d(?:/|__)Input-\d{5}"
normid(m) = replace(m, "__" => "/") * ".txt"
function scanned_ids(; maxids = 50_000)
    ex = Set{String}(); src = Dict{String,Int}(); skipped = String[]
    dirs = [d for d in readdir(ROOTWT; join = true) if isdir(joinpath(d, "validate")) &&
            (startswith(basename(d), "wt-") || basename(d) == "julia") && basename(d) != "wt-paper"]   # wt-paper is not read
    for d in dirs, (root, _, files) in walkdir(joinpath(d, "validate"))
        for f in files
            path = joinpath(root, f)
            ids = Set{String}()
            for m in eachmatch(IDRX, f)                    # raw-run file names
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

function b_exclusions(ids, ps, strat, tab, rows)
    ex, nv = all_exclusions(ids, ps, strat, tab, rows)
    sc, src, skipped = scanned_ids()
    n0 = length(ex)
    union!(ex, sc)
    ex, nv, n0, length(sc), src, skipped
end

# H: E6′'s procedure on a new pool (seeds BPOOLSEED, BPICKSEED)
function bh_inputs(ids, ps, strat, ex)
    rng = MersenneTwister(BPOOLSEED)
    pool = Int[]
    for s in eachindex(M6.STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        got = 0
        for i in el[randperm(rng, length(el))]
            kc = drift_critical_rate(ps[i]).k
            (kc > 0 && 0.3 <= ps[i].MoveFitnessSD / kc <= 3.0) || continue
            push!(pool, i); got += 1
            got == BNPOOL && break
        end
    end
    k50s = sweep(i -> (; k50 = emulated_half_survival_rate(ps[i]; burnin = BURNIN, horizon = HORIZON).k50), pool;
                 id = i -> ids[i], cost = i -> ps[i].NoQTL, progress = 200, cache = bcache("blocks_pool.csv"))
    cand = Tuple{Int,Float64}[]
    for (i, k) in zip(pool, k50s)
        0 < k.k50 < Inf || continue
        r = ps[i].MoveFitnessSD / k.k50
        M6.RLO <= r <= M6.RHI && push!(cand, (i, r))
    end
    sort!(cand; by = c -> ids[c[1]])
    edges = exp.(range(log(M6.RLO), log(M6.RHI); length = 5))
    bin(r) = min(4, searchsortedlast(edges, r))
    rng = MersenneTwister(BPICKSEED)
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

function bf_sample(ids, strat, ex)
    rng = MersenneTwister(BFSEED)
    S = Int[]
    for s in eachindex(M6.STRATA)
        el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in ex)]
        append!(S, el[randperm(rng, length(el))[1:BNF]])
    end
    sort(S)
end

# --- the kernel on H: E6′'s protocol, V_G each generation in the windows, and at the TIMES the
# per-locus genic variance and the block-genic variance at BSTAR blocks (genomes summed over the
# blocks); reading the World draws nothing
function genic_and_blocks(w, bm)
    G = genomes(w)
    L, P, N = size(G)
    genic = 0.0
    buf = Vector{Float64}(undef, P * N)
    for l in 1:L
        k = 0
        for i in 1:N, c in 1:P
            k += 1; buf[k] = G[l, c, i]
        end
        genic += P * var(buf; corrected = false)
    end
    Y = zeros(bm.nb, P * N)
    for i in 1:N, c in 1:P, l in 1:L
        Y[bm.blockof[l], P*(i-1)+c] += G[l, c, i]
    end
    genic, block_stats(Y, bm.nb, P * N).blockgenic
end

function bk_job(j)
    q = kernel_params(j.p)
    bm = block_map(q, BSTAR; within = BWITHIN)
    rng = MersenneTwister(j.seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    vw = Dict(t => (0.0, 0) for t in (:pre, 50, 100, 200))
    gen = Dict{Int,Float64}(); bgen = Dict{Int,Float64}()
    extinct_gen = -1
    for it in w.gen:q.NoGenerations-1
        w.gen = it
        t = it - MOVE
        if t in TIMES
            g, b = genic_and_blocks(w, bm); gen[t] = g; bgen[t] = b
        end
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        key = -200 <= t <= -1 ? :pre : abs(t - 50) <= 5 ? 50 : abs(t - 100) <= 5 ? 100 : abs(t - 200) <= 5 ? 200 : nothing
        if key !== nothing
            v = var(view(w.genotype, 1:w.N))
            a, n = vw[key]; vw[key] = (a + v, n + 1)
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    vm(k, n) = vw[k][2] == n ? vw[k][1] / n : NaN
    pre(d) = mean(d[t] for t in (-200, -100, -1))
    at(d, t) = get(d, t, NaN)
    (; alive = extinct_gen < 0, extinct_gen, VGpre = vm(:pre, 200), VG50 = vm(50, 11), VG100 = vm(100, 11), VG200 = vm(200, 11),
     genicpre = pre(gen), genic50 = at(gen, 50), genic100 = at(gen, 100), genic200 = at(gen, 200),
     bgpre = pre(bgen), bg50 = at(bgen, 50), bg100 = at(bgen, 100), bg200 = at(bgen, 200))
end

# --- the block model's trajectories on H (run_blocks with the kernel's protocol, the same
# statistics; separate from the emulator's replicates, which give P2000)
function bt_job(j)
    q = kernel_params(j.p)
    bg = Dict{Int,Float64}()
    obs(s, bm) = (s.gen - MOVE in TIMES && (bg[s.gen-MOVE] = block_stats(s, bm).blockgenic); nothing)
    r = run_blocks(q, Xoshiro(j.seed); blocks = BSTAR, within = BWITHIN, burnin = KBURN, record_from = MOVE - 200, observe = obs)
    win(lo, hi) = (k = findall(g -> lo <= g - MOVE <= hi, r.gen); length(k) == hi - lo + 1 ? mean(r.VG[k]) : NaN)
    pre = all(t -> haskey(bg, t), (-200, -100, -1)) ? mean(bg[t] for t in (-200, -100, -1)) : NaN
    (; alive = !r.extinct, VGpre = win(-200, -1), VG50 = win(45, 55), VG100 = win(95, 105), VG200 = win(195, 205),
     bgpre = pre, bg50 = get(bg, 50, NaN), bg100 = get(bg, 100, NaN), bg200 = get(bg, 200, NaN))
end

fmean(v) = (u = filter(isfinite, v); isempty(u) ? NaN : mean(u))
function per_input(res, nin, R, fields)
    map(1:nin) do j
        rs = res[(j-1)*R+1:j*R]
        merge((; s = count(r -> r.alive, rs)),
              NamedTuple{fields}(Tuple(fmean([getfield(r, f) for r in rs]) for f in fields)),
              (; n50 = count(r -> isfinite(r.VG50), rs), n100 = count(r -> isfinite(r.VG100), rs), n200 = count(r -> isfinite(r.VG200), rs)))
    end
end

# --- (S1)/(S2): trajectory ratios
function traj_tests(K, T; verdict = true)
    v(b) = verdict ? pf(b) : "-"
    out = Dict{Symbol,Bool}()
    rpre = [T[i].VGpre / K[i].VGpre for i in eachindex(K)]
    ok = inband(med(rpre), 0.95, 1.05); out[:pre] = ok
    lsayf("    (S1) V_G before the move, model/kernel: median %.3f (IQR %.3f–%.3f; n %d; in [0.95, 1.05]) %s",
          med(rpre), qq(rpre, 0.25), qq(rpre, 0.75), count(isfinite, rpre), v(ok))
    for t in (50, 100, 200)
        f = Symbol("VG$t"); nf = Symbol("n$t")
        ix = [i for i in eachindex(K) if getfield(K[i], nf) >= 10 && getfield(T[i], nf) >= 10]
        r = [getfield(T[i], f) / getfield(K[i], f) for i in ix]
        ok = inband(med(r), 0.9, 1.1); out[Symbol("vg$t")] = ok
        lsayf("    (S1) V_G at +%d, model/kernel (≥ 10 of 40 runs alive in each): median %.3f (IQR %.3f–%.3f; n %d; in [0.9, 1.1]) %s",
              t, med(r), qq(r, 0.25), qq(r, 0.75), length(r), v(ok))
    end
    for t in (50, 100, 200)
        f = Symbol("bg$t"); nf = Symbol("n$t")
        ix = [i for i in eachindex(K) if getfield(K[i], nf) >= 10 && getfield(T[i], nf) >= 10]
        r = [(getfield(T[i], f) / T[i].bgpre) / (getfield(K[i], f) / K[i].bgpre) for i in ix]
        ok = inband(med(r), 0.9, 1.1); out[Symbol("bg$t")] = ok
        kg = [getfield(K[i], Symbol("genic$t")) / K[i].genicpre for i in ix]
        kb = [getfield(K[i], f) / K[i].bgpre for i in ix]
        lsayf("    (S2) block-genic(+%d)/block-genic(pre), model/kernel: median %.3f (IQR %.3f–%.3f; n %d; in [0.9, 1.1]) %s;" *
              " kernel: block-genic ratio median %.3f, per-locus genic ratio median %.3f", t, med(r), qq(r, 0.25), qq(r, 0.75),
              length(r), v(ok), med(kb), med(kg))
    end
    out
end

function cost_mode(ids, ps, F)
    P = F[1:div(length(F), 16):end][1:16]
    lsay("cost per emulator call on one thread ($(Threads.nthreads()) thread), 16 inputs of F (2 per stratum)")
    p0 = ps[P[1]]
    emulate(p0; loci = :blocks, blocks = 10, replicates = 2, seed = 1, burnin = 50, horizon = 50)
    emulate(p0; loci = :spectrum, burnin = 50, horizon = 50); emulate(p0; loci = :spectrum, burnin = 50, horizon = 50, spectrum = SPEC)
    tb = Float64[]; ts = Float64[]; tl = Float64[]; t1 = Float64[]; tk = Float64[]
    simulate(Params(kernel_params(p0); NoGenerations = 20); seed = 1)
    for i in P
        p = ps[i]
        push!(tb, @elapsed emulate(p; loci = :blocks, blocks = BSTAR, within = BWITHIN, replicates = BRF, seed = 1, burnin = BBURN_F, horizon = HORIZON))
        push!(ts, @elapsed emulate(p; loci = :spectrum, burnin = BURNIN, horizon = HORIZON))
        push!(tl, @elapsed emulate(p; loci = :spectrum, burnin = BURNIN, horizon = HORIZON, spectrum = SPEC))
        q = kernel_params(p)
        push!(t1, @elapsed run_blocks(q, Xoshiro(1); blocks = BSTAR, within = BWITHIN, burnin = KBURN))
        push!(tk, @elapsed simulate(q; seed = 1, init = :equilibrium, burnin = KBURN))
        lsayf("  %-34s L %4d: blocks (%d replicates) %.2f s, :spectrum %.2f s, ld = :dynamic %.2f s; one run of the kernel's protocol: blocks %.2f s, kernel %.2f s",
              ids[i], p.NoQTL, BRF, tb[end], ts[end], tl[end], t1[end], tk[end])
    end
    lsayf("  median: blocks %.2f s, :spectrum %.2f s, ld = :dynamic %.2f s; median ratio blocks/:spectrum %.2f, blocks/ld %.2f",
          median(tb), median(ts), median(tl), median(tb ./ ts), median(tb ./ tl))
    for L in (100, 500, 1000, 5000)
        k = [j for (j, i) in enumerate(P) if ps[i].NoQTL == L]
        isempty(k) || lsayf("  L %4d: one run, kernel/blocks time ratio median %.2f", L, median(tk[k] ./ t1[k]))
    end
    lsay("  (S4) median cost per call of the block emulator ≤ that of :spectrum with ld = :dynamic: $(pf(median(tb) <= median(tl)))")
end

function bmain()
    t0 = time()
    ids, ps, strat = inv_(() -> M6.read_inputs())
    tab = inv_(() -> M6.load_master(M6.MASTER))
    find = inv_(() -> M6.row_index(tab))
    rows = [find(M6.STRATA[strat[i]][2], M6.STRATA[strat[i]][3], ps[i]) for i in eachindex(ps)]
    ex, nv, n0, nsc, src, skipped = b_exclusions(ids, ps, strat, tab, rows)
    BMODE in ("select", "pilot", "cost") || (LLOG[] = open(joinpath(HERE, "blocks.log"), "w"); M6.LOG[] = LLOG[])
    lsay("validate/blocks.jl — the haplotype-block lineage emulator (loci = :blocks, blocks = $BSTAR), under validate/blocks_criteria.md")
    lsay("$(Threads.nthreads()) threads; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    lsay("excluded inputs: $(length(ex)) = spectrum_ld's exclusions $n0 (of which $nv from vignettes 09/11) ∪ $nsc ids found in files of the worktrees' validate/ directories")
    lsay("  id occurrences found (summed over files), by directory: " * join(["$k $v" for (k, v) in sort(collect(src))], "; "))
    lsay("  whole-table files not used for exclusion: " * join(skipped, "; "))
    H, pool, nb = bh_inputs(ids, ps, strat, ex)
    F = bf_sample(ids, strat, union(ex, Set(ids[pool])))
    lsay("H: $(length(H)) inputs (pool $(length(pool)), MersenneTwister($BPOOLSEED); picks MersenneTwister($BPICKSEED)); per stratum and bin " *
         join([join(nb[s, :], "/") for s in 1:8], ", ") * "; per stratum " * join([count(i -> strat[i] == s, H) for s in 1:8], ", "))
    lsay("F: $(length(F)) inputs ($BNF per stratum, MersenneTwister($BFSEED))")
    lsay("disjointness (set intersection with the exclusions): H ∩ ex = $(length(intersect(Set(ids[H]), ex))), pool ∩ ex = " *
         "$(length(intersect(Set(ids[pool]), ex))), F ∩ ex = $(length(intersect(Set(ids[F]), ex))), F ∩ pool = $(length(intersect(F, pool)))")
    BMODE == "select" && return
    BMODE == "cost" && return cost_mode(ids, ps, F)
    if BMODE == "pilot"
        P = F[1:div(length(F), 16):end][1:16]
        t = @elapsed emulate_sweep(ps[P]; ids = ["pilot|" * x for x in ids[P]], loci = :blocks, blocks = BSTAR, within = BWITHIN, replicates = BRF,
                                   seed = BLSEED, burnin = BBURN_F, horizon = HORIZON, cache = bcache("blocks_pilot.csv"))
        lsayf("pilot: blocks emulator on 16 inputs of F in %.0f s wall on %d threads", t, Threads.nthreads())
        return
    end

    # --- H
    lsay("\n== H: fresh inputs near the threshold, $BRK fresh kernel runs each ==")
    jobs = [(; key = "$(ids[i])#$r", p = ps[i], seed = BHKSEED + 100j + r) for (j, i) in enumerate(H) for r in 1:BRK]
    t = @elapsed kr = sweep(bk_job, jobs; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 400, cache = bcache("blocks_kernel.csv"))
    lsayf("  kernel: %d runs, %.0f s wall this session", length(kr), t)
    fl = (:VGpre, :VG50, :VG100, :VG200, :genicpre, :genic50, :genic100, :genic200, :bgpre, :bg50, :bg100, :bg200)
    K = per_input(kr, length(H), BRK, fl)
    s = [k.s for k in K]
    lsayf("  kernel fraction alive at 2000: mean %.4f; inputs with 0 alive %d, all alive %d, in between %d",
          mean(s ./ BRK), count(==(0), s), count(==(BRK), s), count(x -> 0 < x < BRK, s))
    t = @elapsed begin
        tj = [(; key = "$(ids[i])#$r", p = ps[i], seed = job_seed(BTSEED, "$(ids[i])#$r")) for i in H for r in 1:BRK]
        tr = sweep(bt_job, tj; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 800, cache = bcache("blocks_traj.csv"))
        Esp = sweep(i -> em_job(ps[i], (;); eq = false), H; id = i -> "spectrum|$(ids[i])", cost = i -> ps[i].NoQTL, progress = 40,
                    cache = bcache("blocks_spectrum.csv"))
        Eld = sweep(i -> em_job(ps[i], SPEC; eq = false), H; id = i -> "ld|$(ids[i])", cost = i -> ps[i].NoQTL, progress = 40,
                    cache = bcache("blocks_spectrum.csv"))
        Era = emulate_sweep(ps[H]; ids = ["H-ra|" * x for x in ids[H]], loci = :rare_alleles, burnin = BURNIN, horizon = HORIZON,
                            tail = M6.TAIL, progress = 200, cache = bcache("blocks_ra.csv"))
    end
    lsayf("  emulations on H: %.0f s wall this session", t)
    T = per_input(tr, length(H), BRK, (:VGpre, :VG50, :VG100, :VG200, :bgpre, :bg50, :bg100, :bg200))
    rate = [ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in H]
    st = strat[H]
    vl(b) = pf(b); base(b) = "($(pf(b)) as a comparator)"
    xk, _ = inv_(() -> M6.x50_rule(rate, s ./ BRK))
    pooled(E) = (xe = inv_(() -> M6.x50_rule(rate, [e.P2000 for e in E]))[1]; xe / xk)
    lsay("\n  -- model under test: loci = :blocks, blocks = $BSTAR, within :$BWITHIN, $BRK runs, burn-in $BBURN_H --")
    Pb = [t_.s / BRK for t_ in T]
    r6 = inv_(() -> M6.e6pstats(Pb, s, rate, st; label = vl))
    rp = pooled([(; P2000 = x) for x in Pb]); P2 = inband(rp, 0.8, 1.25)
    lsayf("    (P2) pooled x50 ratio %.3f (in [0.8, 1.25]) %s", rp, pf(P2))
    acc(P) = mean((P .> 0.5) .== (s ./ BRK .> 0.5))
    comp = Dict{String,Any}()
    for (name, E) in (("default :spectrum", Esp), ("ld = :dynamic", Eld), (":rare_alleles", Era))
        lsay("\n  -- comparator: $name (deterministic; no verdict) --")
        P = [e.P2000 for e in E]
        inv_(() -> M6.e6pstats(P, s, rate, st; label = base))
        lsayf("    pooled x50 ratio %.3f", pooled(E))
        comp[name] = P
    end
    # (P4) paired comparison of threshold calls with ld = :dynamic (exact McNemar)
    M = s ./ BRK .> 0.5
    cb = (Pb .> 0.5) .== M; cl = (comp["ld = :dynamic"] .> 0.5) .== M
    b10, b01 = count(cb .& .!cl), count(.!cb .& cl)
    pmc = b10 + b01 == 0 ? 1.0 : min(1.0, 2 * Distributions.cdf(Binomial(b10 + b01, 0.5), min(b10, b01)))
    P4 = b10 >= b01
    lsayf("\n    (P4) threshold accuracy blocks %.3f vs ld = :dynamic %.3f vs default :spectrum %.3f; discordant pairs: blocks right %d, ld right %d (exact McNemar p %.3g): blocks not worse %s",
          acc(Pb), acc(comp["ld = :dynamic"]), acc(comp["default :spectrum"]), b10, b01, pmc, pf(P4))
    lsay("\n  -- trajectories (model: $BRK run_blocks runs per input; kernel: the $BRK runs above) --")
    tt = traj_tests(K, T)

    # --- F
    lsay("\n== F: E1–E4 of validate/emulator_criteria.md on $(length(F)) fresh table inputs ==")
    alltab = map(1:8) do s_
        k = findall(==(s_), strat)
        inv_(() -> M6.x50([ps[i].MoveFitnessSD / sqrt(ps[i].PhenotypeVar) for i in k], [tab.T[rows[i]] > 2000 for i in k]))[1]
    end
    t = @elapsed begin
        Fb = emulate_sweep(ps[F]; ids = ["F|" * x for x in ids[F]], loci = :blocks, blocks = BSTAR, within = BWITHIN, replicates = BRF, seed = BLSEED,
                           burnin = BBURN_F, horizon = HORIZON, tail = M6.TAIL, progress = 50, cache = bcache("blocks_emulate.csv"))
        Fsp = sweep(i -> em_job(ps[i], (;); eq = false), F; id = i -> "F-spectrum|$(ids[i])", cost = i -> ps[i].NoQTL, progress = 100,
                    cache = bcache("blocks_spectrum.csv"))
        Fra = emulate_sweep(ps[F]; ids = ["F-ra|" * x for x in ids[F]], loci = :rare_alleles, burnin = BURNIN, horizon = HORIZON,
                            tail = M6.TAIL, progress = 200, cache = bcache("blocks_ra.csv"))
    end
    lsayf("  emulations on F: %.0f s wall this session", t)
    asem(E) = [(; P2000 = e.P2000, P240k = NaN, T = e.T) for e in E]
    lsay("\n  -- model under test: loci = :blocks --")
    rd = inv_(() -> M6.criteria(asem(Fb), F, tab, rows, ps, strat; label = vl, alltab))
    lsay("\n  -- comparator: default :spectrum --")
    inv_(() -> M6.criteria(asem(Fsp), F, tab, rows, ps, strat; label = base, alltab))
    lsay("\n  -- comparator: :rare_alleles --")
    inv_(() -> M6.criteria(Fra, F, tab, rows, ps, strat; label = base, alltab))

    lsay("\nVERDICTS (loci = :blocks, blocks = $BSTAR, within :$BWITHIN): primary (P1) E6′c $(pf(r6[:c])), (P2) pooled x50 $(pf(P2)), " *
         "(P3) x50 per stratum $(pf(r6[:d])), (P4) not worse than ld = :dynamic $(pf(P4)); secondary (S1) V_G pre $(pf(tt[:pre])), " *
         "+50 $(pf(tt[:vg50])), +100 $(pf(tt[:vg100])), +200 $(pf(tt[:vg200])); (S2) block-genic +50 $(pf(tt[:bg50])), " *
         "+100 $(pf(tt[:bg100])), +200 $(pf(tt[:bg200])); (S3) " * join(["$k $(pf(rd[k]))" for k in (:E1, :E2, :E3, :E4)], ", ") *
         "; (S4) see the cost mode's output")
    lsayf("total %.0f min wall this session", (time() - t0) / 60)
    LLOG[] === stdout || close(LLOG[])
end

bmain()
