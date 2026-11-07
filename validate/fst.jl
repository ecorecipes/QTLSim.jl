# Neutral differentiation among demes: the kernel's metapopulations against
# population-genetic theory, through Q_ST (= F_ST for a neutral additive trait,
# Whitlock 1999) and the within- and between-deme genetic variances.
#
# Usage (from the package root; it needs Catlab, which the vignettes environment has):
#     julia -t 4 --project=vignettes validate/fst.jl [R]
# R replicates per scenario (default 400). Writes validate/fst.log.
#
# The model: k demes of N diploids, each running QTSim's generation with every fitness 1
# (`Drift = true`, as vignette 07 makes demes neutral), then `migrate(k)` with the island
# or the ring stepping-stone matrix. The trait is a continuum-of-alleles additive trait
# with mutational variance V_m = `muteVar(p)` per individual per generation. Statistics
# are on genotypes (not phenotypes) of the adults at QTSim's census, after migration:
# the within-deme sample variance V_W (mean over demes), the sample variance among the k
# deme means, Q_ST = σ²/(σ² + 2V_W) with σ² = (k−1)/k·(among − V_W/N), and for the ring
# half the mean squared difference of deme means at each distance. Each replicate
# averages them over every 10th generation of the second half of its run; replicates
# are independent, so standard errors are over replicates (Q_ST's by a bootstrap).
#
# The predictions are exact expectations for the same generations, from the pairwise
# coalescence recursion of `deme_coalescence` (genomes start at zero, as in `World`),
# given the law of the deme origins of two distinct adults over one migration step.
# Section 1 checks `migrate!`'s marginals against `backward_migration` and shows that
# its exchanges correlate the origins of pairs (so the textbook independent-origins
# "pool" law is off by a relative 1/N); section 2 therefore estimates the pair law from
# `migrate!` itself, by Monte Carlo of the migration step alone, and feeds it to the
# recursion. Everything else in the predictions (drift, mutation, the life-cycle order)
# is exact.
#
# PRE-SPECIFIED CRITERION: every tested statistic (the marginals of section 1; V_W, the
# variance among deme means, Q_ST and, on the ring, the squared differences by distance
# in section 2) gets z = (simulated − predicted)/SE, SE combining the simulation's and
# the prediction's Monte Carlo error; the script exits 1 if any |z| ≥ 4. With about 40
# comparisons a false alarm has probability about 0.003. The closed forms at
# equilibrium (qst_island, Takahata's and Wright's approximations, Slatkin's ring), the
# pool law's predictions and the moments lowering are reported, not tested.

using QTLSim, Random, Statistics, Printf
import Catlab
using Catlab: ⋅, ⊗

const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 400
const LOG = joinpath(@__DIR__, "fst.log")
const ZMAX = 4.0
const L = lifecycle()
const G = L.generation

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s))
const zs = Tuple{String,Float64}[]
check!(name, z) = (push!(zs, (name, z)); z)

ringdist(a, b, k) = min(mod(a - b, k), mod(b - a, k))

# the forward rate whose backward migrant fraction 1 − P[a, a] is m, by bisection
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

# --- 1. migrate! against backward_migration, and the law of pairs' origins ---

"""
Empirical composition after `migrate!`: each individual carries its deme of origin in
its environment slot (which `migrate!` exchanges with it). Returns, per call, the
fraction of deme a from deme c and the probabilities that two distinct individuals of
one deme, or of two different demes, share an origin; and, in `nbatch` batches of
calls, the law of the origins of two distinct individuals as `deme_coalescence` takes
it (`origins[a + (b−1)k, c + (d−1)k]`), from the counts n[a, c] of deme a's individuals
that came from deme c: E[n_ac n_bd]/N² for a ≠ b, E[n_ac (n_ad − [c = d])]/(N(N − 1))
for a = b.
"""
function migration_check(M, k, N, ncalls, seed; nbatch = 4)
    p = Params(PopSize = N, NoQTL = 1, Ploidy = 2, NoGenerations = 1, Drift = true)
    s = p ⊗ demes(k)
    rng = Xoshiro(seed)
    Pf = zeros(ncalls, k, k)
    qs, qd = zeros(ncalls), zeros(ncalls)
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
        Pf[t, :, :] = n ./ N
        qs[t] = sum(n[a, c] * (n[a, c] - 1) for a in 1:k, c in 1:k) / (k * N * (N - 1))
        qd[t] = sum(n[a, c] * n[b, c] for a in 1:k, b in 1:k, c in 1:k if a != b) / (k * (k - 1) * N^2)
        bt = mod1(t, nbatch)
        v = vec(n)
        Q[bt] .+= v .* v'
        nsum[bt] .+= n
    end
    J = map(1:nbatch) do bt
        nc = length(bt:nbatch:ncalls)
        Jb = zeros(k^2, k^2)
        for a in 1:k, b in 1:k, c in 1:k, d in 1:k
            x = Q[bt][a + (c - 1) * k, b + (d - 1) * k]
            Jb[a + (b - 1) * k, c + (d - 1) * k] = a == b ? (x - (c == d) * nsum[bt][a, c]) / (nc * N * (N - 1)) :
                                                            x / (nc * N^2)
        end
        Jb
    end
    (; Pf, qs, qd, J)
end

qpair(J, k) = (mean(sum(J[a + (a - 1) * k, c + (c - 1) * k] for c in 1:k) for a in 1:k),
               mean(sum(J[a + (b - 1) * k, c + (c - 1) * k] for c in 1:k) for a in 1:k, b in 1:k if a != b))

function report_migration(name, M, k, N; ncalls = 20_000, seed = 1)
    P = backward_migration(M, N)
    o = migration_check(M, k, N, ncalls, seed)
    say("\n", name, ": k = $k, N = $N, $ncalls calls of migrate!")
    say(@sprintf("  %-34s %12s %12s %10s %7s", "quantity", "predicted", "migrate!", "SE", "z"))
    island = all(M[a, b] == M[1, 2] for a in 1:k, b in 1:k if a != b)
    classes = island ? [0, 1] : sort(unique(ringdist(a, c, k) for a in 1:k, c in 1:k))
    cls(a, c) = island ? Int(a != c) : ringdist(a, c, k)
    for d in classes                             # tested: the marginals P
        idx = [(a, c) for a in 1:k, c in 1:k if cls(a, c) == d]
        pred = mean(P[a, c] for (a, c) in idx)
        per = [mean(o.Pf[t, a, c] for (a, c) in idx) for t in 1:ncalls]
        se = std(per) / sqrt(ncalls)
        z = check!("$name P(distance $d)", (mean(per) - pred) / se)
        say(@sprintf("  %-34s %12.6f %12.6f %10.2e %7.2f", "P[a,c], distance $d", pred, mean(per), se, z))
    end
    # reported, not tested: pairs of distinct individuals sharing an origin, against the
    # pool law (independent origins) and the first-order exchange law
    say(@sprintf("  %-34s %12s %12s %12s %10s %8s", "pairs of distinct individuals", "pool", "exchange", "migrate!", "SE", "z(pool)"))
    qp, qe = qpair(kron(P, P), k), qpair(exchange_origins(P, N), k)
    for (j, lab, v) in ((1, "same deme, same origin", o.qs), (2, "different demes, same origin", o.qd))
        se = std(v) / sqrt(ncalls)
        say(@sprintf("  %-34s %12.6f %12.6f %12.6f %10.2e %8.2f", lab, qp[j], qe[j], mean(v), se, (mean(v) - qp[j]) / se))
    end
    # deme sizes are kept, so (N−1)q + (k−1)N q′ = N − 1 in every call (pairs counted exactly)
    ident = maximum(abs((N - 1) * o.qs[t] + (k - 1) * N * o.qd[t] - (N - 1)) for t in 1:ncalls)
    say(@sprintf("  size conservation (N−1)q + (k−1)Nq′ = N−1: max deviation over calls %.1e; the pool law gives %.4f",
                 ident, (N - 1) * qp[1] + (k - 1) * N * qp[2]))
    if island
        m = 1 - P[1, 1]
        λe = mean(o.qs) - mean(o.qd)
        r = λe / (2N * (1 - λe)) * (k - 1) / k
        say(@sprintf("  λ = q − q′: pool %.6f, migrate! %.6f; equilibrium Q_ST %.5f (pool) vs %.5f (migrate!'s λ)",
                     qp[1] - qp[2], λe, qst_island(N, m; k), r / (r + 1)))
    end
end

# --- 2. metapopulation runs ---------------------------------------------------

function runmeta(p, k, M, cycle, seed, census)
    s = p ⊗ demes(k)
    idx = IdDict(w => i for (i, w) in enumerate(s.demes))
    mg = zeros(k, length(census)); vg = zeros(k, length(census))
    rec = operator(:record, function (w, rng)
        j = searchsortedfirst(census, w.gen)
        if j <= length(census) && census[j] == w.gen
            g = view(w.genotype, 1:w.N); a = idx[w]
            mg[a, j] = mean(g); vg[a, j] = var(g)
        end
        w
    end)
    o = run_composed(s, Xoshiro(seed); cycle, record = false, images = Dict(:record => rec), migration_rate = M)
    any(o.extinct_gen .>= 0) && error("a neutral deme died")
    (; mg, vg)
end

# per-replicate summaries, averaged over the census generations
function summaries(r, k, N, ring)
    nc = size(r.mg, 2)
    within = mean(r.vg)
    among = mean(var(view(r.mg, :, j)) for j in 1:nc)
    dists = ring ? (1:k÷2) : (1:0)
    pairs = [mean(mean((r.mg[a, j] - r.mg[b, j])^2 / 2 for a in 1:k for b in 1:k if a < b && ringdist(a, b, k) == d)
                  for j in 1:nc) for d in dists]
    (; within, among, pairs)
end

runjob(plans, (i, rep)) = (pl = plans[i];
    summaries(runmeta(pl.p, pl.k, pl.M, pl.cycle, 100_000 * i + rep, pl.census), pl.k, pl.N, pl.ring))

qst_of(within, among, k, N) = (σ2 = (k - 1) / k * (among - within / N); σ2 / (σ2 + 2within))

function scenarios()
    N = 50
    sc = []
    for k in (4, 10), Nm in (0.1, 1.0, 10.0)
        push!(sc, (name = "island k=$k Nm=$Nm", k, N, m = Nm / N, ring = false, Mof = island_migration))
    end
    for Nm in (1.0, 10.0)
        push!(sc, (name = "ring k=8 Nm=$Nm", k = 8, N, m = Nm / N, ring = true, Mof = (k, r) -> stepping_stone(k, r)))
    end
    sc
end

function main()
    say("QTLSim validate/fst.jl  ", Dates_now())
    say("Julia threads: ", Threads.nthreads(), "; replicates per scenario R = ", R)
    say("Pre-specified criterion: exit 1 if any |z| ≥ $ZMAX, z = (simulated − predicted)/SE.")

    say("\n== 1. migrate! against backward_migration, and the origins of pairs ==")
    report_migration("island r=0.05", island_migration(4, 0.05), 4, 50)
    report_migration("island r=0.005", island_migration(10, 0.005), 10, 50)
    report_migration("ring r=0.1", stepping_stone(8, 0.1), 8, 50)

    say("\n== 2. neutral metapopulations against the exact recursion ==")
    p0 = Params(Ploidy = 2, NoQTL = 20, RecombPoissonMean = 1.0, PhenotypeVar = 1.0, Heritability = 0.5,
                MutationScalar = 0.002, MutationRatePerGamete = 0.05, Drift = true)
    Vm = muteVar(p0)
    say(@sprintf("N = 50 diploids per deme, %d loci, %g mutations per chromosome per generation, V_m = %g", p0.NoQTL,
                 p0.NoQTL * p0.MutationRatePerGamete, Vm))
    say("Predictions use migrate!'s own law of pairs' origins, estimated from 40,000 calls in 4 batches;")
    say("their SE (from the batches) is combined with the simulation's in z. 'pool' is the prediction with")
    say("independent origins (the textbook island model), shown for comparison.")
    plans = map(enumerate(scenarios())) do (i, s)
        r = forward_rate(s.Mof, s.k, s.N, s.m)
        M = s.Mof(s.k, r)
        P = backward_migration(M, s.N)
        Js = migration_check(M, s.k, s.N, 40_000, 7000 + i).J
        J = sum(Js) / length(Js)
        eq = neutral_deme_variances(P, s.N, Vm; origins = J)
        T = maximum(eq.S)                         # the slowest mean coalescence time
        gens = ceil(Int, 5T / 100) * 100
        census = collect(gens ÷ 2:10:gens-1)
        p = Params(p0; PopSize = s.N, NoGenerations = gens)
        cycle = foldl(⊗, fill(G, s.k)) ⋅ L.migrate(s.k)
        (; s..., r, M, P, Js, J, eq, gens, census, p, cycle)
    end
    jobs = [(i, rep) for i in eachindex(plans) for rep in 1:R]
    sort!(jobs; by = j -> -plans[j[1]].gens * plans[j[1]].k)
    res = Vector{Any}(undef, length(jobs))
    t0 = time()
    Threads.@threads :greedy for j in eachindex(jobs)
        res[j] = runjob(plans, jobs[j])    # a function: `i` and `pl` are also variables of main, which a
                                           # closure would share between threads
    end
    say(@sprintf("%d runs in %.0f s", length(jobs), time() - t0))

    brng = Xoshiro(20261107)
    for (i, pl) in enumerate(plans)
        rs = [res[j] for j in eachindex(jobs) if jobs[j][1] == i]
        k, N = pl.k, pl.N
        dists = pl.ring ? (1:k÷2) : (1:0)
        bydist(x, d) = mean(x.pairs[a, b] for a in 1:k for b in 1:k if a < b && ringdist(a, b, k) == d)
        # the predicted statistics, averaged over the census generations, for a law of origins
        function predict(origins)
            th = neutral_deme_variances(pl.P, N, Vm; t = pl.census, origins)
            tw, ta = mean(mean(x.within) for x in th), mean(x.among for x in th)
            [tw; ta; qst_of(tw, ta, k, N); [mean(bydist(x, d) for x in th) for d in dists]]
        end
        pred = predict(pl.J)
        pbatch = reduce(hcat, predict(J) for J in pl.Js)
        pse = vec(std(pbatch; dims = 2)) ./ sqrt(size(pbatch, 2))
        pool = predict(:pool)
        w, a = [x.within for x in rs], [x.among for x in rs]
        boot = [(ix = rand(brng, 1:length(rs), length(rs)); qst_of(mean(w[ix]), mean(a[ix]), k, N)) for _ in 1:2000]
        sims = [mean(w); mean(a); qst_of(mean(w), mean(a), k, N); [mean(x.pairs[d] for x in rs) for d in dists]]
        ses = [std(w); std(a); std(boot) * sqrt(length(rs)); [std([x.pairs[d] for x in rs]) for d in dists]] ./ sqrt(length(rs))
        eqs = [mean(pl.eq.within); pl.eq.among; pl.eq.qst; [bydist(pl.eq, d) for d in dists]]
        labels = ["V_W (within demes)"; "var of deme means"; "Q_ST"; ["½(Ḡa−Ḡb)², distance $d" for d in dists]]
        m = 1 - pl.P[1, 1]
        say("\n", pl.name, @sprintf(": forward rate %.5f, backward m = %.5f, Nm = %.3g; %d generations, census %d:10:%d, %d replicates",
                                     pl.r, m, N * m, pl.gens, first(pl.census), last(pl.census), length(rs)))
        say(@sprintf("  %-26s %10s %9s %10s %9s %21s %6s %10s %10s", "statistic", "predicted", "pred. SE", "simulated", "SE",
                     "95% CI", "z", "pool", "equilib."))
        for j in eachindex(labels)
            se = sqrt(ses[j]^2 + pse[j]^2)
            z = check!("$(pl.name) $(labels[j])", (sims[j] - pred[j]) / se)
            say(@sprintf("  %-26s %10.5g %9.2g %10.5g %9.2g [%9.4g, %9.4g] %6.2f %10.5g %10.5g", labels[j], pred[j], pse[j],
                         sims[j], ses[j], sims[j] - 1.96ses[j], sims[j] + 1.96ses[j], z, pool[j], eqs[j]))
        end
        if pl.ring
            say("  pairwise F_ST at equilibrium (two demes, total convention): exact vs Slatkin 1991 eq 17:")
            for d in dists
                σ2 = (bydist(pl.eq, d) - pl.eq.within[1] / N) / 2
                say(@sprintf("    distance %d: %.4f exact, %.4f Slatkin", d, σ2 / (σ2 + 2pl.eq.within[1]), fst_stepping_stone(N, m, k, d)))
            end
        else
            say(@sprintf("  equilibrium Q_ST: %.4f (migrate!'s law), %.4f (qst_island, pool law), Takahata %.4f, Wright 1/(1+4Nm) %.4f",
                         pl.eq.qst, qst_island(N, m; k), 1 / (1 + 4N * m * (k / (k - 1))^2), 1 / (1 + 4N * m)))
            say(@sprintf("  equilibrium V_W: %.5g (migrate!'s law), 2NkV_m = %.5g (Lande 1992), 2NV_m(kN−1)/(N−1) = %.5g (exchange)",
                         mean(pl.eq.within), 2N * k * Vm, 2N * Vm * (k * N - 1) / (N - 1)))
        end
        say(@sprintf("  at the census the expected V_W is %.1f%% and the variance among means %.1f%% of equilibrium",
                     100pred[1] / eqs[1], 100pred[2] / eqs[2]))
    end

    # --- 3. the moments lowering ---
    say("\n== 3. the moments lowering (run_moments, k demes, migration_mixing, erosion) ==")
    pl = plans[2]                                 # island k = 4, Nm = 1
    mo = run_moments(fill(pl.p, pl.k); cycle = pl.cycle, erosion = true, migration_rate = pl.M)
    last_gen = maximum(mo.gen)
    sel = mo.gen .== last_gen
    say(@sprintf("%s, generation %d: moments V_G within %.4g (exact expectation %.4g), variance among deme means %.4g (exact %.4g)",
                 pl.name, last_gen, mean(mo.geno_var[sel]), mean(neutral_deme_variances(pl.P, pl.N, Vm; t = last_gen, origins = pl.J).within),
                 var(mo.pheno_mean[sel]), neutral_deme_variances(pl.P, pl.N, Vm; t = last_gen, origins = pl.J).among))
    say("The moments lowering is deterministic: deme means that start equal stay equal, so it has no between-deme")
    say("variance, and erosion removes V_G/(2N) per deme without the variance migrants bring back from other demes;")
    say("its within-deme fixed point is an isolated deme's 2N·V_m, not the metapopulation's 2Nk·V_m.")

    bad = [(n, z) for (n, z) in zs if !(abs(z) < ZMAX)]
    say("\n", length(zs), " comparisons; max |z| = ", @sprintf("%.2f", maximum(abs(z) for (_, z) in zs)))
    if isempty(bad)
        say("PASS: every |z| < $ZMAX")
    else
        say("FAIL: ", length(bad), " comparison(s) with |z| ≥ $ZMAX:")
        foreach(((n, z),) -> say(@sprintf("  %-50s z = %.2f", n, z)), bad)
    end
    write(LOG, String(take!(out)))
    isempty(bad) || exit(1)
end

Dates_now() = (t = time(); Libc.strftime("%Y-%m-%d %H:%M", t))

abspath(PROGRAM_FILE) == (@__FILE__) && main()
