# Design phase of the haplotype-block lineage representation (src/block_lineages.jl,
# `run_blocks`, `emulate(p; loci = :blocks)`). Development inputs only: the spectrum
# diagnosis's D and Dctl sets (validate/spectrum_diagnosis_inputs.csv; seen by every study of
# the linked-selection chain), their 10 raw kernel trajectories per input
# (../wt-spectrum-diagnosis/validate/spectrum_diagnosis_raw, ENV RAWDIR) and their 40 E6′ kernel
# runs (validate/emulator_spectrum_kernel.csv). Running log with the hypotheses written before
# each step: validate/blocks_dev.log.
#
# Run from the package root (vignettes environment):
#   julia -t 8 --project=vignettes validate/blocks_dev.jl check     # step 0: one block per locus vs the kernel (100 loci)
#   julia -t 8 --project=vignettes validate/blocks_dev.jl check400  # step 0b: the same, pre-move V_G, 400 runs each
#   julia -t 8 --project=vignettes validate/blocks_dev.jl conv      # step 1: convergence in the number of blocks (non-recombining blocks)
#   julia -t 8 --project=vignettes validate/blocks_dev.jl position  # step 1, :position partition
#   julia -t 8 --project=vignettes validate/blocks_dev.jl gauss     # step 1b: the same with Gaussian segregation within blocks
#   julia -t 8 --project=vignettes validate/blocks_dev.jl kcompare  # step 2: the cached runs against the kernel on D + Dctl
#   julia -t 8 --project=vignettes validate/blocks_dev.jl assoc     # step 3: block-level association, kernel vs model
#   julia -t 8 --project=vignettes validate/blocks_dev.jl cost      # step 4: cost per emulator call
# Caches (resumable `sweep` CSVs, git-ignored): validate/blocks_dev_raw/.

ENV["RAWDIR"] = get(ENV, "RAWDIR", normpath(joinpath(@__DIR__, "..", "..", "wt-spectrum-diagnosis", "validate",
                                                     "spectrum_diagnosis_raw")))
include(joinpath(@__DIR__, "spectrum_diagnosis.jl"))   # readp, kernel_params, load_runs, load_e6p, read_inputs_csv, inset, genic_stats, MOVE, KBURN, KSEED
using QTLSim: job_seed, block_map, block_stats

const DRAW = joinpath(@__DIR__, "blocks_dev_raw")
const DEVSEED = 51116
const BS = (5, 10, 20, 40, 80, 160, 320)
const RDEV = 40
const BSG = (10, 20, 40, 80, 160)
devrows() = inset(read_inputs_csv(), "D+Dctl")

# one block-model run with the kernel's protocol (equilibrium start 1000 generations before
# the move, to generation 12000): alive at +2000, death generation, V_G before and during the move
function bjob(p, B, seed; partition = :map, within = :snap)
    q = kernel_params(p)
    r = run_blocks(q, Xoshiro(seed); blocks = B, burnin = KBURN, record_from = MOVE - 200, partition, within)
    win(lo, hi) = (k = findall(g -> lo <= g - MOVE <= hi, r.gen); length(k) == hi - lo + 1 ? mean(r.VG[k]) : NaN)
    (; alive = !r.extinct, extinct_gen = r.extinct_gen, VGpre = win(-200, -1), VG50 = win(45, 55), VG100 = win(95, 105),
     VG200 = win(195, 205))
end

function bsweep(tag, rows, Bs; R = RDEV, partition = :map, within = :snap)
    mkpath(DRAW)
    jobs = [(; key = "$tag|$(r.id)|$B|$k", p = readp(r.id), B, seed = job_seed(DEVSEED, "$tag|$(r.id)|$B|$k"))
            for r in rows for B in Bs if B <= readp(r.id).NoQTL for k in 1:R]
    res = sweep(j -> bjob(j.p, j.B, j.seed; partition, within), jobs; id = j -> j.key, cost = j -> j.B, progress = 500,
                cache = joinpath(DRAW, "blocks_dev_$tag.csv"))
    Dict(j.key => x for (j, x) in zip(jobs, res))
end

# per-input summaries of a set of block runs
function bsum(res, tag, id, B; R = RDEV)
    rs = [res["$tag|$id|$B|$k"] for k in 1:R if haskey(res, "$tag|$id|$B|$k")]
    f(x) = (u = filter(isfinite, x); isempty(u) ? NaN : mean(u))
    (; P = count(r -> r.alive, rs) / length(rs), VGpre = f([r.VGpre for r in rs]), VG50 = f([r.VG50 for r in rs]),
     VG100 = f([r.VG100 for r in rs]), VG200 = f([r.VG200 for r in rs]), n100 = count(r -> isfinite(r.VG100), rs),
     n200 = count(r -> isfinite(r.VG200), rs))
end

# the kernel's trajectories (10 raw runs) and E6′ survival (40 runs) of an input
function ksum(id, kr)
    runs = load_runs(id, 10)
    win(x, lo, hi) = (k = findall(g -> lo <= g - MOVE <= hi, x.gen); length(k) == hi - lo + 1 ? mean(x.VG[k]) : NaN)
    gat(x, t) = (k = findfirst(==(MOVE + t), x.xg); k === nothing ? NaN : x.genic[k])
    f(v) = (u = filter(isfinite, v); isempty(u) ? NaN : mean(u))
    gpre = f([mean([x.genic[k] for k in eachindex(x.xg) if -200 <= x.xg[k] - MOVE <= -10]) for x in runs])
    (; P = count(x -> x[2], kr[id]) / length(kr[id]), VGpre = f([win(x, -200, -1) for x in runs]),
     VG50 = f([win(x, 45, 55) for x in runs]), VG100 = f([win(x, 95, 105) for x in runs]), VG200 = f([win(x, 195, 205) for x in runs]),
     n100 = count(x -> isfinite(win(x, 95, 105)), runs), genicpre = gpre,
     genic50 = f([gat(x, 50) for x in runs]), genic100 = f([gat(x, 100) for x in runs]), genic200 = f([gat(x, 200) for x in runs]))
end

med(v) = (u = filter(isfinite, v); isempty(u) ? NaN : median(u))
lr(a, b) = a > 0 && b > 0 && isfinite(a) && isfinite(b) ? log(a / b) : NaN

# --- step 0 ----------------------------------------------------------------------------------
function step_check()
    rows = filter(r -> startswith(r.id, "QT0100"), devrows())
    res = bsweep("check", rows, (100,); R = RDEV)
    kr, _ = load_e6p()
    println("step 0: one block per locus (B = L = 100) vs the kernel, $(length(rows)) inputs, $RDEV block runs, kernel: 10 trajectories + 40 E6′ runs")
    for r in rows
        b = bsum(res, "check", r.id, 100); k = ksum(r.id, kr)
        bp = binomial_p(round(Int, b.P * RDEV), RDEV, round(Int, k.P * 40), 40)
        @printf("  %-32s alive@2000 blocks %.3f kernel %.3f (Fisher-type p %.3f); V_G pre blocks %.5f kernel %.5f (ratio %.3f); +100 %.5f / %.5f\n",
                r.id, b.P, k.P, bp, b.VGpre, k.VGpre, b.VGpre / k.VGpre, b.VG100, k.VG100)
    end
end

# two-sample test of proportions (Fisher's exact, two-sided, by enumeration)
function binomial_p(a, n, b, m)
    t = a + b
    lp(x) = logpdf(Hypergeometric(n, m, t), x)
    lo, hi = max(0, t - m), min(n, t)
    p0 = lp(a)
    sum(exp(lp(x)) for x in lo:hi if lp(x) <= p0 + 1e-9)
end

# step 0b: the pre-move V_G at B = L, many runs (R = 80: kernel seeds 99001.., block job_seed(77, r);
# R = 400: kernel seeds 99101.., block job_seed(78, r))
function step_check_powered(R = 400; kseed = R == 80 ? 99_000 : 99_100, bseed = R == 80 ? 77 : 78)
    for id in ("QT0100-Possion5.0/Input-07602.txt", "QT0100-Possion5.0/Input-06159.txt")
        q = Params(kernel_params(readp(id)); NoGenerations = MOVE)
        kv = zeros(R); bv = zeros(R)
        Threads.@threads :greedy for r in 1:R
            o = simulate(q; seed = kseed + r, init = :equilibrium, burnin = KBURN)
            kv[r] = mean(o.geno_var[end-199:end])
            b = run_blocks(q, Xoshiro(job_seed(bseed, r)); blocks = 100, burnin = KBURN)
            bv[r] = mean(b.VG[end-199:end])
        end
        t = (mean(bv) - mean(kv)) / sqrt(var(bv) / R + var(kv) / R)
        @printf("%s: V_G pre (mean of −200..−1), %d runs each: kernel %.5f (sd %.5f), blocks B=L %.5f (sd %.5f), ratio %.3f, Welch t %.2f\n",
                id, R, mean(kv), std(kv), mean(bv), std(bv), mean(bv) / mean(kv), t)
    end
end

# --- step 1 ----------------------------------------------------------------------------------
function step_conv(; partition = :map, tag = "conv", within = :snap, Bs = BS)
    rows = devrows()
    res = bsweep(tag, rows, Bs; partition, within)
    println("step 1 ($tag, partition :$partition, within :$within): per input and B, P(alive at +2000), V_G before the move, V_G at +100, over $RDEV runs")
    tab = Dict{Tuple{String,Int},Any}()
    for r in rows, B in BS
        B <= readp(r.id).NoQTL || continue
        tab[(r.id, B)] = bsum(res, tag, r.id, B)
    end
    for r in rows
        L = readp(r.id).NoQTL
        line = join([@sprintf("B%d %.2f/%.4f/%.4f", B, tab[(r.id, B)].P, tab[(r.id, B)].VGpre, tab[(r.id, B)].VG100) for B in Bs if B <= L], "  ")
        @printf("  %-32s %s\n", r.id, line)
    end
    println("\n  successive doublings B → 2B (over inputs with both): median |log V_G ratio| before / at +100 / at +200, mean |ΔP|, calls P > 1/2 that differ")
    for (B, B2) in zip(Bs[1:end-1], Bs[2:end])
        ids = [r.id for r in rows if B2 <= readp(r.id).NoQTL]
        isempty(ids) && continue
        a = [tab[(i, B)] for i in ids]; b = [tab[(i, B2)] for i in ids]
        @printf("  %4d → %4d (n %2d): %.3f / %.3f / %.3f, %.3f, %d\n", B, B2, length(ids),
                med(abs.(lr.([x.VGpre for x in a], [x.VGpre for x in b]))), med(abs.(lr.([x.VG100 for x in a], [x.VG100 for x in b]))),
                med(abs.(lr.([x.VG200 for x in a], [x.VG200 for x in b]))),
                mean(abs.([x.P for x in a] .- [x.P for x in b])), count(k -> (a[k].P > 0.5) != (b[k].P > 0.5), eachindex(a)))
    end
    println("\n  the same by locus count (median |log V_G ratio| at +100; mean |ΔP|)")
    for L in (100, 500, 1000, 5000)
        ids = [r.id for r in rows if readp(r.id).NoQTL == L]
        s = String[]
        for (B, B2) in zip(Bs[1:end-1], Bs[2:end])
            B2 <= L || continue
            a = [tab[(i, B)] for i in ids]; b = [tab[(i, B2)] for i in ids]
            push!(s, @sprintf("%d→%d %.3f %.2f", B, B2, med(abs.(lr.([x.VG100 for x in a], [x.VG100 for x in b]))),
                              mean(abs.([x.P for x in a] .- [x.P for x in b]))))
        end
        println("  L = $L ($(length(ids)) inputs): ", join(s, "; "))
    end
    tab
end

# --- step 2 ----------------------------------------------------------------------------------
function step_compare()
    rows = devrows()
    tab = step_conv()
    kr, _ = load_e6p()
    K = Dict(r.id => ksum(r.id, kr) for r in rows)
    println("\nstep 2: against the kernel (10 trajectories, 40 E6′ runs) on D + Dctl")
    for B in BS
        ids = [r.id for r in rows if B <= readp(r.id).NoQTL]
        length(ids) == length(rows) || continue
        Pb = [tab[(i, B)].P for i in ids]; Pk = [K[i].P for i in ids]
        @printf("  B %3d: threshold accuracy %.3f (%d of %d), mean |P − f| %.3f; median V_G ratio pre %.3f, +50 %.3f, +100 %.3f, +200 %.3f\n",
                B, mean((Pb .> 0.5) .== (Pk .> 0.5)), count((Pb .> 0.5) .== (Pk .> 0.5)), length(ids), mean(abs.(Pb .- Pk)),
                med([tab[(i, B)].VGpre / K[i].VGpre for i in ids]), med([tab[(i, B)].VG50 / K[i].VG50 for i in ids]),
                med([tab[(i, B)].VG100 / K[i].VG100 for i in ids]), med([tab[(i, B)].VG200 / K[i].VG200 for i in ids]))
    end
    println("\n  baselines on the same inputs: :spectrum (default) and ld = :dynamic (deterministic emulators; P2000, V_G ratios)")
    for (name, kw) in (("spectrum", (;)), ("ld = :dynamic", (; ld = :dynamic)))
        E = Dict(r.id => spectrum_traj(readp(r.id); spectrum = kw) for r in rows)
        P = Float64[]; vr = [Float64[] for _ in 1:4]
        for r in rows
            st = E[r.id]
            d = length(st.m.gen) <= BURNIN ? 0.0 : extinction_distribution(st.q, st.m; from = BURNIN).psurvive
            push!(P, d)
            w(lo, hi) = (k = findall(t -> lo <= t <= hi, st.t); length(k) == hi - lo + 1 ? mean(st.VG[k]) : NaN)
            for (j, (lo, hi), kf) in zip(1:4, ((-200, -1), (45, 55), (95, 105), (195, 205)), (:VGpre, :VG50, :VG100, :VG200))
                push!(vr[j], w(lo, hi) / getfield(K[r.id], kf))
            end
        end
        Pk = [K[r.id].P for r in rows]
        @printf("  %-14s threshold accuracy %.3f (%d of %d), mean |P − f| %.3f; median V_G ratio pre %.3f, +50 %.3f, +100 %.3f, +200 %.3f\n",
                name, mean((P .> 0.5) .== (Pk .> 0.5)), count((P .> 0.5) .== (Pk .> 0.5)), length(rows), mean(abs.(P .- Pk)),
                med(vr[1]), med(vr[2]), med(vr[3]), med(vr[4]))
    end
    println("\n  per input: kernel f, blocks P at B = 20/40/80/160, kernel V_G +100, blocks V_G +100 at 40/80/160; kernel genic pre/+100")
    for r in rows
        k = K[r.id]; L = readp(r.id).NoQTL
        g(B, f) = B <= L ? getfield(tab[(r.id, B)], f) : NaN
        @printf("  %-32s %-5s f %.3f  P %.2f/%.2f/%.2f/%.2f  VG100 %.5f | %.5f %.5f %.5f  genic %.5f → %.5f\n", r.id, r.set, k.P,
                g(20, :P), g(40, :P), g(80, :P), g(160, :P), k.VG100, g(40, :VG100), g(80, :VG100), g(160, :VG100), k.genicpre, k.genic100)
    end
end

# step 2′: the cached runs of steps 1 and 1b against the kernel (10 trajectories, 40 E6′ runs) on
# D + Dctl, and :spectrum's threshold calls on the same inputs (its P2000 from the diagnosis's
# inputs file)
function step_kcompare()
    rows = devrows(); kr, _ = load_e6p()
    K = Dict(r.id => ksum(r.id, kr) for r in rows)
    Pk = Dict(r.id => K[r.id].P for r in rows)
    @printf("  default :spectrum (deterministic, as in E6′): threshold calls agree with the kernel on %d of %d (D: %d of %d)\n",
            count(r -> (r.P > 0.5) == (Pk[r.id] > 0.5), rows), length(rows),
            count(r -> r.set == "D" && (r.P > 0.5) == (Pk[r.id] > 0.5), rows), count(r -> r.set == "D", rows))
    for (tag, Bs, w) in (("conv", BS, :snap), ("gmap", (10, 20, 40, 80, 160, 320), :gaussian))
        res = bsweep(tag, rows, Bs; within = w)
        for B in Bs
            ids = [r.id for r in rows if B <= readp(r.id).NoQTL]
            T = Dict(i => bsum(res, tag, i, B) for i in ids)
            @printf("  %-5s B %3d (n %d): blocks/kernel V_G median pre %.3f +50 %.3f +100 %.3f +200 %.3f; threshold calls agree %d/%d; mean |P − f| %.3f\n",
                    w, B, length(ids), med([T[i].VGpre / K[i].VGpre for i in ids]), med([T[i].VG50 / K[i].VG50 for i in ids]),
                    med([T[i].VG100 / K[i].VG100 for i in ids]), med([T[i].VG200 / K[i].VG200 for i in ids]),
                    count(i -> (T[i].P > 0.5) == (K[i].P > 0.5), ids), length(ids), mean(abs(T[i].P - K[i].P) for i in ids))
        end
    end
end

# --- step 3: the block-level association, kernel (genomes summed over the same blocks) vs model
const ATS = (-1, 25, 50, 100, 200)
function kassoc_job(p, B, seed)
    q = kernel_params(p)
    bm = block_map(q, B)
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    out = Dict{Int,Tuple{Float64,Float64,Float64}}()
    for it in w.gen:MOVE+200
        w.gen = it
        if it - MOVE in ATS
            G = genomes(w)
            L, P, N = size(G)
            Y = zeros(bm.nb, 2N)
            for i in 1:N, c in 1:2, l in 1:L
                Y[bm.blockof[l], 2(i-1)+c] += G[l, c, i]
            end
            s = block_stats(Y, bm.nb, 2N)
            out[it-MOVE] = (s.blockgenic, s.assoc, genic_stats(w; seg = false)[1])
        end
        generation!(w, rng)
        w.extinct && break
    end
    (; [Symbol("bg$(t < 0 ? "m1" : t)") => (haskey(out, t) ? out[t][1] : NaN) for t in ATS]...,
     [Symbol("as$(t < 0 ? "m1" : t)") => (haskey(out, t) ? out[t][2] : NaN) for t in ATS]...,
     [Symbol("ge$(t < 0 ? "m1" : t)") => (haskey(out, t) ? out[t][3] : NaN) for t in ATS]...)
end
function bassoc_job(p, B, seed)
    q = kernel_params(p)
    out = Dict{Int,Tuple{Float64,Float64}}()
    obs(s, bm) = (s.gen - MOVE in ATS && (st = block_stats(s, bm); out[s.gen-MOVE] = (st.blockgenic, st.assoc)); nothing)
    run_blocks(Params(q; NoGenerations = MOVE + 201), Xoshiro(seed); blocks = B, burnin = KBURN, record_from = typemax(Int), observe = obs,
               within = Symbol(get(ENV, "BWITHIN", "gaussian")))
    (; [Symbol("bg$(t < 0 ? "m1" : t)") => (haskey(out, t) ? out[t][1] : NaN) for t in ATS]...,
     [Symbol("as$(t < 0 ? "m1" : t)") => (haskey(out, t) ? out[t][2] : NaN) for t in ATS]...)
end
function step_assoc(; B = parse(Int, get(ENV, "BASSOC", "80")), RK_ = 4, RB = 16)
    rows = inset(read_inputs_csv(), "D")
    mkpath(DRAW)
    kj = [(; key = "k|$(r.id)|$B|$k", p = readp(r.id), seed = KSEED + 100r.j + k) for r in rows for k in 1:RK_]
    K = sweep(j -> kassoc_job(j.p, B, j.seed), kj; id = j -> j.key, cost = j -> j.p.NoQTL, progress = 10,
              cache = joinpath(DRAW, "blocks_dev_kassoc.csv"))
    bj = [(; key = "b|$(r.id)|$B|$k", p = readp(r.id), seed = job_seed(DEVSEED, "assoc|$(r.id)|$B|$k")) for r in rows for k in 1:RB]
    Bm = sweep(j -> bassoc_job(j.p, B, j.seed), bj; id = j -> j.key, cost = _ -> B, progress = 50,
               cache = joinpath(DRAW, "blocks_dev_bassoc.csv"))
    println("step 3: block-level statistics at B = $B on D ($(length(rows)) inputs; kernel $RK_ runs, model $RB runs, means over runs alive)")
    println("  per input: block-genic(t)/block-genic(−1) kernel | model at t = 50, 100, 200;  Σ assoc / block-genic² kernel | model at 50, 100, 200;  kernel per-locus genic(t)/genic(−1)")
    rat = Dict(t => Float64[] for t in (50, 100, 200)); ras = Dict(t => Float64[] for t in (50, 100, 200))
    for r in rows
        m(v, f) = (u = filter(isfinite, [getfield(x, f) for x in v]); isempty(u) ? NaN : mean(u))
        ks = [K[i] for i in eachindex(kj) if kj[i].key |> k -> occursin("|$(r.id)|", k)]
        bs = [Bm[i] for i in eachindex(bj) if occursin("|$(r.id)|", bj[i].key)]
        kb0, bb0 = m(ks, :bgm1), m(bs, :bgm1)
        s = String[]
        for t in (50, 100, 200)
            kb, bb = m(ks, Symbol("bg$t")) / kb0, m(bs, Symbol("bg$t")) / bb0
            ka, ba = m(ks, Symbol("as$t")) / kb0^2, m(bs, Symbol("as$t")) / bb0^2
            push!(s, @sprintf("%.3f|%.3f  %+.3f|%+.3f", kb, bb, ka, ba))
            push!(rat[t], lr(bb, kb)); isfinite(ka) && isfinite(ba) && push!(ras[t], ba - ka)
        end
        ge = [m(ks, Symbol("ge$t")) / m(ks, :gem1) for t in (50, 100, 200)]
        @printf("  %-32s %s   genic %.3f %.3f %.3f\n", r.id, join(s, "   "), ge...)
    end
    for t in (50, 100, 200)
        @printf("  t = %3d: median log(model/kernel) of block-genic(t)/(−1) %.3f (n %d); median model − kernel association/bg² %+.3f\n",
                t, med(rat[t]), count(isfinite, rat[t]), med(ras[t]))
    end
end

# --- step 4: cost ------------------------------------------------------------------------------
function step_cost()
    rows = devrows()
    println("step 4: wall time of one emulator call on one thread (R = 40, burn-in 1000, horizon 2000), by locus count; kernel: 40 runs of the same protocol")
    for L in (100, 500, 1000, 5000)
        r = first(r for r in rows if readp(r.id).NoQTL == L)
        p = readp(r.id)
        emulate(p; loci = :blocks, blocks = 10, replicates = 2, seed = 1, burnin = 100, horizon = 100)  # compile
        for (B, w) in ((160, :gaussian), (320, :gaussian), (320, :snap))
            t = @elapsed e = emulate(p; loci = :blocks, blocks = B, within = w, replicates = 40, seed = 2, burnin = 1000, horizon = 2000)
            @printf("  L %4d  %-30s B %3d %-8s: %.1f s for 40 replicates (P2000 %.2f)\n", L, r.id, B, w, t, e.P2000)
        end
        t = @elapsed emulate(p; loci = :spectrum, burnin = 2000, horizon = 2000)
        t2 = @elapsed emulate(p; loci = :spectrum, burnin = 2000, horizon = 2000, spectrum = (; ld = :dynamic))
        q = kernel_params(p)
        tk = @elapsed for k in 1:4
            simulate(q; seed = 1000 + k, init = :equilibrium, burnin = KBURN)
        end
        @printf("  L %4d  :spectrum %.1f s, ld = :dynamic %.1f s; kernel 4 runs %.1f s (×10 for 40)\n", L, t, t2, tk)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "conv" : ARGS[1]
    mode == "check" ? step_check() : mode == "check400" ? step_check_powered(400) : mode == "conv" ? step_conv() :
    mode == "position" ? step_conv(; partition = :position, tag = "pos") :
    mode == "gauss" ? (step_conv(; within = :gaussian, tag = "gmap", Bs = BSG); step_conv(; partition = :position, within = :gaussian, tag = "gpos", Bs = BSG)) : mode == "compare" ? step_compare() :
    mode == "assoc" ? step_assoc() : mode == "kcompare" ? step_kcompare() : mode == "cost" ? step_cost() : error("unknown mode $mode")
end
