# Realised coalescence times and identity by descent against the coalescence theory.
#
# `deme_coalescence` predicts the pairwise coalescence times of gene copies in the
# kernel's neutral metapopulations; `validate/fst.jl` tests it through genetic variances
# only. Here the times themselves are read off recorded pedigrees (`record_pedigree!`,
# `pair_coalescence`, src/ancestry.jl) and compared with it, together with the identity
# of blocks (`ibd_sharing`) and the tree-sequence tables.
#
# Usage (from the package root; the vignettes environment has Catlab):
#     julia -t 3 --project=vignettes validate/coalescence.jl [R]
# R replicate runs per scenario (default 300). Writes validate/coalescence.log.
#
# The model, as in fst.jl: k demes of N = 50 diploids with 20 loci, each running QTSim's
# generation with every fitness 1 (`Drift = true`), then `migrate!` with the island or the
# ring stepping-stone matrix, written out by hand (`generation!` on each deme, then
# `migrate!`; the test suite checks that this is the lowered life cycle, bit for bit).
# Scenarios: one deme; islands of k = 4 and 10 demes at Nm = 0.1, 1, 10 (m the backward
# migrant fraction); a ring of 8 demes at Nm = 1, 10. The census is after migration. The
# pedigree is recorded from generation 0 and simplified every 50 generations.
#
# Predictions: with x_t = `deme_coalescence(P, N; t)` (W for the two copies in an
# individual, S[a, b] for copies in distinct individuals of demes a and b), for a census t
# generations after recording began,
#     E[min(T, t)] = x_t,     P(T ≤ s) = 1 − (x_{s+1} − x_s)   for s ≤ t,
# the second because x_{s+1} − x_s = P(T > s) and the backward recursion does not depend
# on the date of the census. The law of the deme origins of two distinct adults over one
# `migrate!` (`origins`) is estimated from `migrate!` itself, 40,000 calls in 4 batches, as
# in fst.jl; the batches give the prediction's Monte Carlo error. One deme needs no
# such estimate. Nothing else is estimated or fitted.
#
# PRE-SPECIFIED CRITERIA (written before the first run of this script):
#
# A. Coalescence times. At censuses t = 50, 200 and 800, for each class of pairs (within
#    an individual; distinct individuals of one deme; different demes, on the ring by
#    distance) the statistics E[min(T, t)] and P(T ≤ s), s ∈ {1, 10, 100, t} with s ≤ t,
#    each the mean over runs of the run's mean over all pairs of the class and all loci.
#    Pairs within a run are not independent, so SE is the standard deviation across runs
#    over √R, combined with the prediction's Monte Carlo SE; z = (simulated − predicted)/SE.
#    A P(T ≤ s) with prediction p is tested only if the predicted number of coalesced
#    pairs, summed over runs (p · pairs · loci · R), is at least 100 and p ≤ 0.99, so that
#    z is roughly normal; otherwise it is listed as "sparse" and not tested. Where the
#    prediction is exactly 0 or 1, the simulation must equal it exactly. FAIL if any tested
#    |z| ≥ 4. With about 450 tested comparisons the chance of a false alarm is about 3% if
#    the z were independent and normal.
# B. Identity of blocks against the pedigree. (1) For every run of A, in a second run of
#    the same seed whose pedigree is never simplified, at t = 200, for every block position
#    and every locus in it, chromosomes that hold the same block must have the same
#    `pedigree_ibd` label (identical by descent with no mutation in the block since the
#    common ancestor): zero violations. (2) In separate runs with blocks
#    of one locus (`blocksize = 1`), mutation rate μ = 0.005 per locus, census t = 200:
#    holding the same block must be *equivalent* to the same label, at every locus: zero
#    mismatches. (3) In those runs the probabilities of block identity by class
#    (`ibd_sharing`) are tested against Σ_{T ≤ t} P(T = T) e^{−2μT}, the same z test.
# C. Exact checks. In every run, the pairwise coalescence times at t = 200 from the
#    pedigree that was never simplified (the second run of B(1)) equal those from the
#    simplified one. On the first run of every scenario: the unsimplified tables pass
#    `check_tables` with complete cover; every allele of the final population rebuilt from the pedigree equals
#    the kernel's, bit for bit; the tree-sequence tables pass `check_tables`, survive the
#    text round trip (`write_tskit`, `read_tskit`), and give the pedigree's coalescence
#    times (`table_tmrca`) for 200 random pairs and loci. Any failure is a FAIL.
#
# The script exits 1 on any FAIL and keeps every comparison in the log. tskit itself is
# not available here: the text files are checked in Julia only (see `write_tskit`).
# Also reported, not tested: the cost of recording, the `:pool` and `:exchange` laws'
# predictions, F_ST from block identity.
#
# AMENDED after a smoke test of 6 runs per scenario (2026-11-10), before the first full
# run; nothing was changed after a full run. The smoke test showed two faults of the
# criteria as first written, both in the test and not in what it tests:
# - B(1) was to be checked at t = 800 on the simplified pedigree, and found 5 "violations"
#   in 54 runs. They were artefacts: simplification keeps a mutation on the lineage of its
#   own locus only, so the label of a block of several loci cannot be computed from a
#   simplified pedigree (`pedigree_ibd` now refuses to). B(1) is therefore checked on a
#   pedigree that is never simplified, at t = 200, and C's comparison of simplified and
#   unsimplified coalescence times, which uses the same second run, covers every run
#   instead of the first.
# - A's rule for testing a P(T ≤ s) was "at least 100 predicted pairs on the rarer side".
#   For p near 1 the rare pairs are the uncoalesced ones, which sit in a few runs with a
#   deep split (one deme, t = 800: p = 0.9996, and all 6 runs gave exactly 1, SE 0), so
#   their mean over runs is far from normal. The rule is now p · pairs · loci · R ≥ 100 and
#   p ≤ 0.99. Recent coalescences (small p) are many nearly independent events and stay
#   tested.
# With 6 runs the z are t-distributed with 5 degrees of freedom, so that run's |z| ≥ 4
# (6 of 410) say nothing; it is not kept.
#
# RUN HISTORY. Run 1 (coalescence_run1.log, R = 300): PASS, 328 z tests, max |z| 2.87. But
# the script had a bug, made with the amendment above: the rule p ≤ 0.99 was also applied to
# the means E[min(T, t)], which are all above 0.99, so none of them was tested (they are the
# 81 of that log's 101 "sparse" rows that are means), against criterion A. Also, a
# prediction below 1e-9 was taken as exactly 0 and required a simulated 0, which is wrong for
# the tiny but positive probabilities between distant demes of the ring; only P(T ≤ 1) within
# an individual is exactly 0. Both were fixed in the script (the criteria above are
# unchanged) and it was run again with the same seeds, so with the same simulated values:
# coalescence.log.

using QTLSim, Random, Statistics, Printf
import Catlab
using Catlab: ⊗

const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
const LOG = joinpath(@__DIR__, "coalescence.log")
const ZMAX = 4.0
const CENSUS = (50, 200, 800)
const SS = (1, 10, 100)
const EVERY = 50
const MINEVENTS = 100
const NDEME = 50
const MU1 = 0.005                 # mutation rate per locus in the one-locus-block runs
const T1 = 200                    # and their census

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); flush(stdout); println(out, s))
const zs = Tuple{String,Float64}[]
const exact = Tuple{String,Bool}[]
const sparse = Ref(0)

ringdist(a, b, k) = min(mod(a - b, k), mod(b - a, k))

# the forward rate whose backward migrant fraction 1 − P[a, a] is m, by bisection (fst.jl)
function forward_rate(Mof, k, N, m)
    lo, hi = 0.0, 0.99
    f(r) = 1 - backward_migration(Mof(k, r), N)[1, 1] - m
    f(hi) < 0 && error("backward fraction $m is out of reach")
    for _ in 1:60
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    (lo + hi) / 2
end

# the law of the deme origins of two distinct adults over one migrate!, from `ncalls` calls
# in `nbatch` batches (fst.jl's migration_check): each individual carries its deme of origin
# in its environment slot, which migrate! exchanges with it
function origins_law(M, k, N, ncalls, seed; nbatch = 4)
    p = Params(PopSize = N, NoQTL = 1, Ploidy = 2, NoGenerations = 1, Drift = true)
    s = p ⊗ demes(k)
    rng = Xoshiro(seed)
    Q = [zeros(k^2, k^2) for _ in 1:nbatch]
    nsum = [zeros(k, k) for _ in 1:nbatch]
    n = zeros(k, k)
    for t in 1:ncalls
        for (a, w) in enumerate(s.demes); fill!(view(w.env, 1:w.N), a); end
        migrate!(s, rng; rate = M)
        fill!(n, 0)
        for (a, w) in enumerate(s.demes), i in 1:w.N
            n[a, Int(w.env[i])] += 1
        end
        bt = mod1(t, nbatch)
        v = vec(n)
        Q[bt] .+= v .* v'
        nsum[bt] .+= n
    end
    map(1:nbatch) do bt
        nc = length(bt:nbatch:ncalls)
        Jb = zeros(k^2, k^2)
        for a in 1:k, b in 1:k, c in 1:k, d in 1:k
            x = Q[bt][a + (c - 1) * k, b + (d - 1) * k]
            Jb[a + (b - 1) * k, c + (d - 1) * k] = a == b ? (x - (c == d) * nsum[bt][a, c]) / (nc * N * (N - 1)) :
                                                            x / (nc * N^2)
        end
        Jb
    end
end

# --- classes of pairs and statistics ---------------------------------------------------------

# the classes of pairs of a scenario: (label, function of (W, S) -> the class mean, pairs per locus per run)
function classes(k, N, ring)
    n = 2N
    cl = Any[("within individual", (W, S) -> mean(W), k * N),
             ("within deme", (W, S) -> mean(S[a, a] for a in 1:k), k * (n * (n - 1) ÷ 2 - N))]
    if ring
        for d in 1:k÷2
            prs = [(a, b) for a in 1:k for b in a+1:k if ringdist(a, b, k) == d]
            push!(cl, ("demes at distance $d", (W, S) -> mean(S[a, b] for (a, b) in prs), length(prs) * n^2))
        end
    elseif k > 1
        push!(cl, ("between demes", (W, S) -> mean(S[a, b] for a in 1:k for b in a+1:k), k * (k - 1) ÷ 2 * n^2))
    end
    cl
end

# the statistics of one census, in a fixed order: per class, E[min(T, t)] and P(T ≤ s)
statlist(t) = [:mean; [s for s in (SS..., t) if s <= t && (s != t || !(t in SS))]]

function census_stats(pc, cl, t)
    m = coalescence_mean(pc)
    v = Float64[]
    for (_, f, _) in cl, st in statlist(t)
        if st === :mean
            push!(v, f(m.W, m.S))
        else
            c = coalescence_cdf(pc, st)
            push!(v, f(c.W, c.S))
        end
    end
    v
end

# the predictions for the same statistics from the sequence x_0, x_1, … of deme_coalescence
function census_pred(xs, cl, t)
    v = Float64[]
    for (_, f, _) in cl, st in statlist(t)
        x(u) = f(xs[u+1].W, xs[u+1].S)
        push!(v, st === :mean ? x(t) : 1 - (x(st + 1) - x(st)))
    end
    v
end

# P(same block) for blocks of one locus: Σ_{T ≤ t} P(T = T) e^{−2μT}
function ibd_pred(xs, cl, t, μ)
    map(cl) do (_, f, _)
        x(u) = f(xs[u+1].W, xs[u+1].S)
        sum(exp(-2μ * T) * ((x(T) - x(T - 1)) - (x(T + 1) - x(T))) for T in 1:t)
    end
end

# --- one run ------------------------------------------------------------------------------------

chrs(ws) = [(d, i, c) for (d, w) in enumerate(ws) for i in 1:w.N for c in 1:2]   # in the order of samples(ped)

# block identity against the pedigree at the present: (violations of "same block ⇒ same label",
# mismatches of "same block ⇔ same label", pairs with the same block, pairs with the same label)
function block_check(ped, ws)
    pool = ws[1].pool
    ch = chrs(ws)
    viol = mism = nblock = nlabel = 0
    pairs(counts) = sum(c * (c - 1) ÷ 2 for c in values(counts); init = 0)
    for kb in 1:pool.nb
        lo, hi = (kb - 1) * pool.B + 1, min(kb * pool.B, pool.L)
        ids = [ws[d].chr[kb, c, i] for (d, i, c) in ch]
        idcount = Dict{Int32,Int}()
        for id in ids; idcount[id] = get(idcount, id, 0) + 1; end
        for l in lo:hi
            lab = pedigree_ibd(ped, l; block = lo:hi)
            labof = Dict{Int32,Int32}(); idof = Dict{Int32,Int32}()
            labcount = Dict{Int32,Int}()
            for (id, lb) in zip(ids, lab)
                get!(labof, id, lb) == lb || (viol += 1)        # one block, two labels
                get!(idof, lb, id) == id || (mism += 1)         # one label, two blocks
                labcount[lb] = get(labcount, lb, 0) + 1
            end
            nblock += pairs(idcount); nlabel += pairs(labcount)
        end
    end
    (; viol, mism, nblock, nlabel)
end

function newmeta(p, k; blocksize = nothing)
    blocksize === nothing && return p ⊗ demes(k)
    pool = QTLSim.BlockPool(p.NoQTL, blocksize)
    s = p ⊗ demes(k)
    empty!(s.demes); append!(s.demes, [World(p; pool) for _ in 1:k])
    s
end

function advance!(s, rng, M, k)
    foreach(w -> generation!(w, rng), s.demes)
    k > 1 && migrate!(s, rng; rate = M)
end

function runA(pl, rep, i)
    k = pl.k
    s = newmeta(pl.p, k)
    ped = record_pedigree!(s; simplify_every = EVERY)
    rng = Xoshiro(200_000 * i + rep)
    foreach(w -> init_env!(w, rng), s.demes)
    stats = Vector{Float64}[]
    maxbytes = bytes200 = 0
    pc200 = pair_coalescence(ped)
    for t in 1:last(CENSUS)
        advance!(s, rng, pl.M, k)
        t % EVERY == EVERY - 1 && (maxbytes = max(maxbytes, pedigree_bytes(ped)))
        if t in CENSUS
            pc = pair_coalescence(ped)
            push!(stats, census_stats(pc, pl.cl, t))
            t == 200 && (pc200 = pc; bytes200 = pedigree_bytes(ped))
        end
    end
    # the same run with a pedigree that is never simplified, to t = 200: B(1) and C
    s2 = newmeta(pl.p, k); full = record_pedigree!(s2)
    rng = Xoshiro(200_000 * i + rep)
    foreach(w -> init_env!(w, rng), s2.demes)
    for t in 1:200; advance!(s2, rng, pl.M, k); end
    pcf = pair_coalescence(full)
    simp = pcf.W == pc200.W && pcf.S == pc200.S
    bc = block_check(full, s2.demes)
    complete = rep == 1 ? (try check_tables(full; complete = true) catch e; false end) : true
    checks = rep == 1 ? exact_checks(pl, i, ped, s) : nothing
    (; stats, bc, simp, complete, maxbytes, bytes = pedigree_bytes(ped), nodes = length(ped.node_time),
     full_bytes = pedigree_bytes(full), bytes200, checks)
end

# criterion C, on the first run of a scenario
function exact_checks(pl, i, ped, s)
    # every allele from the pedigree
    smp = samples(ped)
    ch = chrs(s.demes)
    G = [genomes(w) for w in s.demes]
    alle = all(1:pl.p.NoQTL) do l
        x = alleles(ped, l)
        all(x[smp.node[j]] === G[d][l, c, ii] for (j, (d, ii, c)) in enumerate(ch))
    end
    # the tables
    tb = tskit_tables(ped)
    valid = try check_tables(tb) catch e; say("  check_tables: ", sprint(showerror, e)); false end
    dir = mktempdir()
    write_tskit(dir, tb)
    back = read_tskit(dir)
    round = back == tb
    rng = Xoshiro(i)
    trees = all(1:200) do _
        l = rand(rng, 1:pl.p.NoQTL); a, b = rand(rng, smp.node, 2)
        a == b || table_tmrca(back, l, a - 1, b - 1) == tmrca(ped, l, a, b)
    end
    (; alle, valid, round, trees,
     rows = (length(tb.nodes.time), length(tb.edges.left), length(tb.individuals.flags), length(tb.mutations.site),
             length(tb.migrations.node)))
end

# criterion B(2), B(3): blocks of one locus
function runB(pl, rep, i)
    k = pl.k
    p = Params(pl.p; MutationRatePerGamete = MU1)
    s = newmeta(p, k; blocksize = 1)
    ped = record_pedigree!(s; simplify_every = EVERY)
    rng = Xoshiro(900_000 + 200_000 * i + rep)
    foreach(w -> init_env!(w, rng), s.demes)
    for t in 1:T1; advance!(s, rng, pl.M, k); end
    sh = ibd_sharing(s.demes)
    (; bc = block_check(ped, s.demes), ibd = [f(sh.ind, sh.between) for (_, f, _) in pl.cl],
     fst = k > 1 ? ibd_fst(s.demes).fst : NaN)
end

# --- the cost of recording -----------------------------------------------------------------------

function cost()
    say("\n== 0. The cost of recording (reported, not tested) ==")
    say(@sprintf("1-minute load average %.1f; times are the fastest of 3 runs of 500 generations after 1500, on one thread", Sys.loadavg()[1]))
    say(@sprintf("  %-41s %7s %12s %14s %16s", "N = 1000 diploids", "ms/gen", "× plain", "bytes/gen", "table at end, MB"))
    for (L, u) in ((100, 3.0), (1000, 3.0))
        p = Params(PopSize = 1000, NoQTL = L, Ploidy = 2, NoGenerations = 10^6, RecombPoissonMean = u, PhenotypeVar = 0.03,
                   Heritability = 0.3, MutationScalar = 0.003, MutationRatePerGamete = 1 / L, InitFitnessSD = 1.0)
        plain = 0.0
        for (name, every) in (("plain", -1), ("pedigree, never simplified", 0), ("pedigree, simplified every 100", 100),
                              ("pedigree, simplified every 20", 20))
            best, bytes, total, founders = Inf, 0.0, 0, 0
            for seed in 1:3
                w = World(p); rng = Xoshiro(seed); init_env!(w, rng)
                for _ in 1:1500; generation!(w, rng); end
                ped = every < 0 ? nothing : record_pedigree!(w; simplify_every = every)
                founders = ped === nothing ? 0 : pedigree_bytes(ped)
                t = @elapsed for _ in 1:500; generation!(w, rng); end
                best = min(best, t)
                ped === nothing || (bytes = (pedigree_bytes(ped) - founders) / 500; total = pedigree_bytes(ped))
            end
            every < 0 && (plain = best)
            say(@sprintf("  L = %-5d %-31s %7.3f %12.2f %14s %16s", L, name, 2best, best / plain,
                         every == 0 ? @sprintf("%.0f", bytes) : "", every < 0 ? "" : @sprintf("%.2f", total / 2^20)))
            every == 0 && say(@sprintf("            (of which the founders, with their %d non-zero alleles as mutations at time 0: %.2f MB)",
                                       (founders - 36 * p.PopSize) ÷ 20, founders / 2^20))
        end
    end
    say("  (one mutation per gamete, crossovers Poisson(3) of which QTSim applies one; recording starts at generation 1500, from")
    say("  genomes that are not zero; the simplified tables' size is at generation 500 of the recording, just simplified)")
end

# --- main ----------------------------------------------------------------------------------------

function scenarios()
    sc = Any[(name = "one deme", k = 1, m = 0.0, ring = false, Mof = nothing)]
    for k in (4, 10), Nm in (0.1, 1.0, 10.0)
        push!(sc, (name = "island k=$k Nm=$Nm", k, m = Nm / NDEME, ring = false, Mof = island_migration))
    end
    for Nm in (1.0, 10.0)
        push!(sc, (name = "ring k=8 Nm=$Nm", k = 8, m = Nm / NDEME, ring = true, Mof = (k, r) -> stepping_stone(k, r)))
    end
    sc
end

# `events` is the predicted number of coalesced pairs over all runs for a P(T ≤ s), `nothing`
# for a mean, which is always tested; `isexact` marks the one prediction that is exactly 0,
# P(T ≤ 1) within an individual (its two copies come from two different parents)
function ztest!(name, sim, se, pred, pse, events; isexact = false)
    if isexact
        ok = abs(pred) < 1e-9 && all(==(0), sim)
        push!(exact, (name, ok))
        return ok ? "exact" : "FAIL (not exact)"
    end
    if events !== nothing && (events < MINEVENTS || pred > 0.99)
        sparse[] += 1
        return "sparse"
    end
    z = (mean(sim) - pred) / sqrt(se^2 + pse^2)
    push!(zs, (name, z))
    @sprintf("%6.2f", z)
end

function main()
    say("QTLSim validate/coalescence.jl  ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    say("Julia threads: ", Threads.nthreads(), "; replicate runs per scenario R = ", R)
    say("Pre-specified criteria: see the header of validate/coalescence.jl. FAIL if any tested |z| ≥ $ZMAX,")
    say("any exact prediction is missed, any block identity contradicts the pedigree, or any exact check fails.")
    cost()

    p0 = Params(Ploidy = 2, NoQTL = 20, RecombPoissonMean = 1.0, PhenotypeVar = 1.0, Heritability = 0.5,
                MutationScalar = 0.002, MutationRatePerGamete = 0.05, Drift = true, PopSize = NDEME, NoGenerations = last(CENSUS))
    N = NDEME
    Tmax = last(CENSUS)
    plans = map(enumerate(scenarios())) do (i, s)
        k = s.k
        if k == 1
            r, M, P = 0.0, zeros(1, 1), ones(1, 1)
            Js = [ones(1, 1)]
        else
            r = forward_rate(s.Mof, k, N, s.m)
            M = s.Mof(k, r)
            P = backward_migration(M, N)
            Js = origins_law(M, k, N, 40_000, 7000 + i)
        end
        J = sum(Js) / length(Js)
        ts = collect(0:Tmax+1)
        xs = deme_coalescence(P, N; t = ts, origins = J)
        xb = [deme_coalescence(P, N; t = ts, origins = Jb) for Jb in Js]
        xpool = deme_coalescence(P, N; t = ts, origins = :pool)
        xexch = deme_coalescence(P, N; t = ts, origins = :exchange)
        eq = deme_coalescence(P, N; origins = J)
        (; s..., r, M, P, xs, xb, xpool, xexch, eq, p = p0, cl = classes(k, N, s.ring))
    end

    say("\n== A. Realised coalescence times against deme_coalescence ==")
    say(@sprintf("N = %d diploids per deme, %d loci, recombination mean %g; censuses at t = %s; pedigrees simplified every %d generations",
                 N, p0.NoQTL, p0.RecombPoissonMean, join(CENSUS, ", "), EVERY))
    jobs = [(i, rep) for i in eachindex(plans) for rep in 1:R]
    sort!(jobs; by = j -> -plans[j[1]].k)
    res = Vector{Any}(undef, length(jobs))
    t0 = time()
    Threads.@threads :greedy for j in eachindex(jobs)
        res[j] = runA(plans[jobs[j][1]], jobs[j][2], jobs[j][1])
    end
    say(@sprintf("%d runs in %.0f s", length(jobs), time() - t0))

    viol = 0
    for (i, pl) in enumerate(plans)
        rs = [res[j] for j in eachindex(jobs) if jobs[j][1] == i]
        k = pl.k
        say("\n", pl.name, k == 1 ? "" : @sprintf(": forward rate %.5f, backward m = %.5f", pl.r, 1 - pl.P[1, 1]),
            @sprintf("; equilibrium mean coalescence times: within individual %.1f, within deme %.1f", mean(pl.eq.W),
                     mean(pl.eq.S[a, a] for a in 1:k)), k == 1 ? "" : @sprintf(", between demes %.1f", mean(pl.eq.S[a, b] for a in 1:k for b in 1:k if a != b)))
        say(@sprintf("  %-22s %5s %-12s %11s %9s %11s %9s %8s %11s %11s", "pairs", "t", "statistic", "predicted", "pred. SE",
                     "simulated", "SE", "z", "pool", "exchange"))
        for (ic, t) in enumerate(CENSUS)
            pred = census_pred(pl.xs, pl.cl, t)
            pb = reduce(hcat, census_pred(x, pl.cl, t) for x in pl.xb)
            pse = size(pb, 2) > 1 ? vec(std(pb; dims = 2)) ./ sqrt(size(pb, 2)) : zeros(length(pred))
            ppool, pexch = census_pred(pl.xpool, pl.cl, t), census_pred(pl.xexch, pl.cl, t)
            sim = reduce(hcat, r.stats[ic] for r in rs)
            j = 0
            for (label, _, npairs) in pl.cl, st in statlist(t)
                j += 1
                x = view(sim, j, :)
                se = std(x) / sqrt(length(x))
                stname = st === :mean ? "E[min(T,t)]" : "P(T ≤ $st)"
                events = st === :mean ? nothing : pred[j] * npairs * p0.NoQTL * length(x)
                z = ztest!("$(pl.name), $label, t=$t, $stname", x, se, pred[j], pse[j], events;
                           isexact = st == 1 && label == "within individual")
                say(@sprintf("  %-22s %5d %-12s %11.6g %9.2g %11.6g %9.2g %8s %11.6g %11.6g", label, t, stname, pred[j], pse[j],
                             mean(x), se, z, ppool[j], pexch[j]))
            end
        end
        v = sum(r.bc.viol for r in rs); viol += v
        say(@sprintf("  B(1) blocks of %d loci at t = 200, pedigree never simplified, %d runs: %d violations of 'same block ⇒ identical by descent without mutation'; blocks find %.1f%% of such pairs",
                     QTLSim.default_blocksize(p0.NoQTL), length(rs), v,
                     100 * sum(r.bc.nblock for r in rs) / sum(r.bc.nlabel for r in rs)))
        say(@sprintf("  pedigree tables: at most %.2f MB between simplifications, %.2f MB and %d nodes (%d living chromosomes) at t = %d (mean over runs)",
                     maximum(r.maxbytes for r in rs) / 2^20, mean(r.bytes for r in rs) / 2^20, round(Int, mean(r.nodes for r in rs)), 2N * k, Tmax))
        c = only(r.checks for r in rs if r.checks !== nothing)
        nsimp = count(r -> r.simp, rs)
        for (name, ok) in (("simplified = unsimplified coalescence times at t = 200, every run", nsimp == length(rs)),
                           ("unsimplified tables complete (run 1)", all(r -> r.complete, rs)), ("alleles rebuilt from the pedigree, bit for bit", c.alle),
                           ("check_tables", c.valid), ("text round trip", c.round), ("table_tmrca = tmrca, 200 pairs", c.trees))
            push!(exact, ("$(pl.name), C: $name", ok))
        end
        say(@sprintf("  C: simplified = unsimplified coalescence times at t = 200 in %d of %d runs (tables %.2f MB against %.2f MB unsimplified, means); run 1: unsimplified tables complete: %s; alleles from the pedigree: %s;",
                     nsimp, length(rs), mean(r.bytes200 for r in rs) / 2^20, mean(r.full_bytes for r in rs) / 2^20,
                     all(r -> r.complete, rs) ? "yes" : "NO", c.alle ? "identical" : "DIFFERENT"))
        say(@sprintf("    tables (%d nodes, %d edges, %d individuals, %d mutations, %d migration rows): check_tables %s, text round trip %s, table_tmrca = tmrca %s",
                     c.rows..., c.valid ? "passes" : "FAILS", c.round ? "identical" : "DIFFERS", c.trees ? "yes" : "NO"))
    end
    push!(exact, ("B(1): no block identity without identity by descent", viol == 0))

    say("\n== B. Blocks of one locus: block identity ⇔ identity by descent without mutation, and its probability ==")
    say(@sprintf("blocksize = 1, μ = %g per locus per generation, census t = %d, founders tagged; prediction Σ_T P(T = T) e^{−2μT}", MU1, T1))
    bi = [i for (i, pl) in enumerate(plans) if pl.name in ("one deme", "island k=4 Nm=1.0", "ring k=8 Nm=1.0")]
    jobsB = [(i, rep) for i in bi for rep in 1:R]
    resB = Vector{Any}(undef, length(jobsB))
    t0 = time()
    Threads.@threads :greedy for j in eachindex(jobsB)
        resB[j] = runB(plans[jobsB[j][1]], jobsB[j][2], jobsB[j][1])
    end
    say(@sprintf("%d runs in %.0f s", length(jobsB), time() - t0))
    for i in bi
        pl = plans[i]
        rs = [resB[j] for j in eachindex(jobsB) if jobsB[j][1] == i]
        v, mm = sum(r.bc.viol for r in rs), sum(r.bc.mism for r in rs)
        push!(exact, ("$(pl.name), B(2): same block ⇔ same pedigree label", v == 0 && mm == 0))
        say("\n", pl.name, @sprintf(": %d runs; %d pairs of chromosomes × loci with the same block; mismatches with the pedigree: %d + %d",
                                     length(rs), sum(r.bc.nblock for r in rs), v, mm))
        say(@sprintf("  %-22s %11s %9s %11s %9s %8s", "pairs", "predicted", "pred. SE", "block ids", "SE", "z"))
        pred = ibd_pred(pl.xs, pl.cl, T1, MU1)
        pb = reduce(hcat, ibd_pred(x, pl.cl, T1, MU1) for x in pl.xb)
        pse = size(pb, 2) > 1 ? vec(std(pb; dims = 2)) ./ sqrt(size(pb, 2)) : zeros(length(pred))
        for (j, (label, _, npairs)) in enumerate(pl.cl)
            x = [r.ibd[j] for r in rs]
            se = std(x) / sqrt(length(x))
            z = ztest!("$(pl.name), blocks of one locus, $label", x, se, pred[j], pse[j], pred[j] * npairs * p0.NoQTL * length(x))
            say(@sprintf("  %-22s %11.6g %9.2g %11.6g %9.2g %8s", label, pred[j], pse[j], mean(x), se, z))
        end
        if pl.k > 1
            say(@sprintf("  F_ST from block identity (ibd_fst): mean over runs %.4f (SE %.4f); reported, not tested", mean(r.fst for r in rs),
                         std([r.fst for r in rs]) / sqrt(length(rs))))
        end
    end

    bad = [(n, z) for (n, z) in zs if !(abs(z) < ZMAX)]
    badx = [n for (n, ok) in exact if !ok]
    say("\n", length(zs), " z tests, max |z| = ", @sprintf("%.2f", maximum(abs(z) for (_, z) in zs)), "; ",
        count(abs(z) >= 3 for (_, z) in zs), " with |z| ≥ 3, ", count(abs(z) >= 2 for (_, z) in zs), " with |z| ≥ 2 (expected about ",
        round(Int, 0.0455 * length(zs)), " if independent and normal); ", length(exact), " exact checks; ", sparse[], " sparse comparisons not tested")
    if isempty(bad) && isempty(badx)
        say("PASS: every tested |z| < $ZMAX and every exact check holds")
    else
        say("FAIL: ", length(bad), " z test(s) with |z| ≥ $ZMAX and ", length(badx), " exact check(s) failed:")
        foreach(((n, z),) -> say(@sprintf("  %-70s z = %.2f", n, z)), bad)
        foreach(n -> say("  ", n), badx)
    end
    write(LOG, String(take!(out)))
    (isempty(bad) && isempty(badx)) || exit(1)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
