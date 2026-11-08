# The moments lowering of a metapopulation with drift among demes (`run_moments(...;
# between = true)`) against the kernel: within-deme genetic variance, variance among deme
# means and Q_ST, neutral and under weak stabilising selection, with one optimum for all
# demes and with an optimum per deme (local adaptation).
#
# Usage (from the package root; it needs Catlab, which the vignettes environment has):
#     julia -t 3 --project=vignettes validate/moments_demes.jl [R]
# R replicates per scenario (default 300). Writes validate/moments_demes.log.
#
# The model: k demes of N = 100 diploids, 50 loci, each running QTSim's generation (a
# step `withparams(g; θ = θ_a)` gives deme a its optimum), then `migrate(k)` at forward
# rate 0.01 (island; backward migrant fraction m ≈ 0.02, Nm ≈ 2). V_E = 1, V_m = 1e-3
# per individual per generation, selection width ω = 5 (V_s = ω² + V_E ≈ 26, β = V_G/s
# about 0.005 per generation: weak). Genomes start at zero. Statistics are on the adults'
# genotypes at QTSim's census (after migration): the within-deme sample variance V_W
# (mean over demes), the sample variance of the k deme means, and Q_ST = σ²/(σ² + 2V_W),
# σ² = (k−1)/k·(among − V_W/N). Each replicate averages them over generations
# 2000:10:3990 of a 4000-generation run; standard errors are over replicates (Q_ST's by
# a bootstrap). Predictions are the same averages of `run_moments(...; between = true)`
# (law of origins `:exchange`, the default), and, for comparison, of the old lowering
# (`erosion = true`, deterministic means; skipped for k = 10, whose `migration_mixing`
# over 20,000 visiting orders per generation is slow).
#
# PRE-SPECIFIED CRITERIA (written before the kernel runs):
# - Neutral scenarios (A, D0): the prediction is exact up to the first-order exchange
#   law, so z = (kernel − predicted)/SE must satisfy |z| < 4 for V_W, the variance among
#   means and Q_ST.
# - Selection scenarios (B, C, D, E): the prediction is the Gaussian approximation, so a
#   relative tolerance: V_W within 10%, the variance among means and Q_ST within 15% of
#   the kernel. A statistic FAILS only if the kernel's whole 95% CI lies outside the
#   prediction's tolerance band, i.e. |kernel − pred| − 1.96 SE > tol·pred.
# - Signatures: under local adaptation (C, E) Q_ST exceeds the neutral F_ST of the same
#   metapopulation (A, D0), and under a common optimum (B, D) it falls below it, in the
#   kernel (95% CIs separated) and in the moments.
# The script exits 1 if any criterion fails.
#
# HISTORY (for honesty about what was fixed before the run): the first version of
# `between = true` reduced the within-deme variance by selection as the one-deme lowering
# does (V_G²/s per generation, genic and lasting, i.e. all variance at one locus) and
# carried migrants' deviations as genic variance. A 6-replicate smoke test of this
# script (seeds 1_000_001…) missed the kernel by a factor 2–5 under selection. The model
# was then rederived for L loci (each locus responds with its share β/L of the
# selection; deterministic mean differences spread over the loci; LD tracked as the
# difference between the full and the genic moments and broken up at the rate
# `pair_recombination(p)` of the kernel's crossover). It has no fitted constant. The
# tolerances above were not changed, and the run below uses fresh seeds (7_000_000 + …).
#
# Rerun on branch gauss (same seeds and criteria) after two changes to the model, which
# is now also the one-deme default (src/multilocus.jl): the genic part of selection's
# reduction uses each locus's share of the full covariance of a haplotype with the
# genotype (LD included), where the first version used its genic part, and LD is resolved
# by class of pairs of loci, each recombining at its own rate (`QTLSim.pair_classes`),
# not at the mean rate. The "old" column is `loci = :aggregate` with `erosion = true`.

using QTLSim, Random, Statistics, Printf
import Catlab
using Catlab: ⋅, ⊗

const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
const LOG = joinpath(@__DIR__, "moments_demes.log")
const L = lifecycle()
const G = L.generation
const RATE = 0.01
const GENS = 4000
const CENSUS = collect(2000:10:GENS-10)

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s))
const failures = String[]

base(N) = Params(Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                 MutationRatePerGamete = 0.02, PopSize = N, NoGenerations = GENS, Sample = 10, InitFitnessSD = 5.0)

function scenarios()
    sc = []
    for (k, tag) in ((5, ""), (10, "10"))
        N = 100
        push!(sc, (name = k == 5 ? "A neutral k=5" : "D0 neutral k=10", k, N, neutral = true, θ = zeros(k), tag = "neutral$tag"))
        push!(sc, (name = k == 5 ? "B common optimum k=5" : "D common optimum k=10", k, N, neutral = false, θ = zeros(k), tag = "common$tag"))
        push!(sc, (name = k == 5 ? "C local optima k=5" : "E local optima k=10", k, N, neutral = false,
                   θ = collect(range(-2.0, 2.0, k)), tag = "local$tag"))
    end
    sc
end

params(s) = s.neutral ? Params(base(s.N); Drift = true) : base(s.N)
cycle(s) = foldl(⊗, [withparams(G; θ = s.θ[a]) for a in 1:s.k]) ⋅ L.migrate(s.k)

function runmeta(s, seed)
    p, k = params(s), s.k
    st = p ⊗ demes(k)
    idx = IdDict(w => i for (i, w) in enumerate(st.demes))
    mg = zeros(k, length(CENSUS)); vg = zeros(k, length(CENSUS))
    rec = operator(:record, function (w, rng)
        j = searchsortedfirst(CENSUS, w.gen)
        if j <= length(CENSUS) && CENSUS[j] == w.gen
            gt = view(w.genotype, 1:w.N); a = idx[w]
            mg[a, j] = mean(gt); vg[a, j] = var(gt)
        end
        w
    end)
    o = run_composed(st, Xoshiro(seed); cycle = cycle(s), record = false, images = Dict(:record => rec),
                     migration_rate = RATE)
    any(o.extinct_gen .>= 0) && error("a deme died in $(s.name)")
    (within = mean(vg), among = mean(var(view(mg, :, j)) for j in axes(mg, 2)))
end

qst_of(within, among, k, N) = (σ2 = (k - 1) / k * (among - within / N); σ2 / (σ2 + 2within))

function predict(s)
    p = params(s)
    mo = run_moments(p; cycle = cycle(s), migration_rate = RATE, between = true)
    sel = findall(in(CENSUS), mo.metapop.gen)
    w, a = mean(mo.metapop.within[sel]), mean(mo.metapop.var_among[sel])
    new = [w, a, qst_of(w, a, s.k, s.N)]
    old = fill(NaN, 3)
    if s.k <= 5
        mo2 = run_moments(p; cycle = cycle(s), migration_rate = RATE, erosion = true, loci = :aggregate)
        ws, as = Float64[], Float64[]
        for t in CENSUS
            rows = findall(==(t), mo2.gen)
            push!(ws, mean(mo2.geno_var[rows])); push!(as, var(mo2.pheno_mean[rows]))
        end
        old = [mean(ws), mean(as), qst_of(mean(ws), mean(as), s.k, s.N)]
    end
    (; new, old, mo)
end

function main()
    say("QTLSim validate/moments_demes.jl")
    say("Julia threads: ", Threads.nthreads(), "; replicates per scenario R = ", R)
    p0 = base(100)
    P = backward_migration(island_migration(5, RATE), 100)
    say(@sprintf("N = 100 per deme, %d loci, V_E = %g, V_m = %g, ω = %g, forward rate %g (backward m = %.4f at k = 5)",
                 p0.NoQTL, envVar(p0), muteVar(p0), p0.InitFitnessSD, RATE, 1 - P[1, 1]))
    say("census generations $(first(CENSUS)):10:$(last(CENSUS)) of $GENS")
    sc = scenarios()
    jobs = [(i, r) for i in eachindex(sc) for r in 1:R]
    sort!(jobs; by = j -> -sc[j[1]].k)
    res = Vector{Any}(undef, length(jobs))
    t0 = time()
    Threads.@threads :greedy for j in eachindex(jobs)
        i, r = jobs[j]
        res[j] = runmeta(sc[i], 7_000_000 + 10_000 * i + r)
    end
    say(@sprintf("%d kernel runs in %.0f s", length(jobs), time() - t0))
    brng = Xoshiro(20261107)
    qst = Dict{String,Tuple{Float64,Float64,Float64}}()      # tag => (kernel, SE, moments)
    labels = ["V_W (within demes)", "var of deme means", "Q_ST"]
    for (i, s) in enumerate(sc)
        rs = [res[j] for j in eachindex(jobs) if jobs[j][1] == i]
        w, a = [x.within for x in rs], [x.among for x in rs]
        n = length(rs)
        boot = [(ix = rand(brng, 1:n, n); qst_of(mean(w[ix]), mean(a[ix]), s.k, s.N)) for _ in 1:2000]
        sims = [mean(w), mean(a), qst_of(mean(w), mean(a), s.k, s.N)]
        ses = [std(w) / sqrt(n), std(a) / sqrt(n), std(boot)]
        pr = predict(s)
        say("\n", s.name, s.neutral ? " (Drift = true)" : "", "; optima ", s.θ == zeros(s.k) ? "all 0" : "range(-2, 2, $(s.k))",
            "; ", n, " replicates")
        say(@sprintf("  %-20s %10s %10s %9s %21s %8s %10s %8s", "statistic", "moments", "kernel", "SE", "95% CI",
                     s.neutral ? "z" : "rel.err", "old", "verdict"))
        for j in 1:3
            d = sims[j] - pr.new[j]
            if s.neutral
                z = d / ses[j]
                ok = abs(z) < 4
                stat = @sprintf("%8.2f", z)
            else
                tol = j == 1 ? 0.10 : 0.15
                ok = abs(d) - 1.96ses[j] <= tol * pr.new[j]
                stat = @sprintf("%+7.1f%%", 100d / pr.new[j])
            end
            ok || push!(failures, "$(s.name) $(labels[j])")
            say(@sprintf("  %-20s %10.5g %10.5g %9.2g [%9.4g, %9.4g] %8s %10.5g %8s", labels[j], pr.new[j], sims[j], ses[j],
                         sims[j] - 1.96ses[j], sims[j] + 1.96ses[j], stat, pr.old[j], ok ? "ok" : "FAIL"))
        end
        qst[s.tag] = (sims[3], ses[3], pr.new[3])
        if !s.neutral && s.θ != zeros(s.k)
            mp = pr.mo.metapop
            sel = findall(in(CENSUS), pr.mo.gen)
            zb = [mean(pr.mo.pheno_mean[r] for r in sel if pr.mo.deme[r] == d) for d in 1:s.k]
            say(@sprintf("  moments: deterministic deme means %s (optima %s)", string(round.(zb; digits = 3)), string(s.θ)))
        end
    end
    say("\n== signatures: Q_ST against the neutral F_ST of the same metapopulation ==")
    for (k10, neu, com, loc) in (("", "neutral", "common", "local"), ("10", "neutral10", "common10", "local10"))
        fk, fse, fm = qst[neu]
        for (tag, want) in ((com, :below), (loc, :above))
            qk, qse, qm = qst[tag]
            sep = want === :above ? (qk - 1.96qse > fk + 1.96fse) : (qk + 1.96qse < fk - 1.96fse)
            mok = want === :above ? qm > fm : qm < fm
            ok = sep && mok
            ok || push!(failures, "signature $tag $want F_ST")
            say(@sprintf("  %-10s Q_ST %s F_ST: kernel %.4f ± %.4f vs %.4f ± %.4f; moments %.4f vs %.4f  %s",
                         tag, want === :above ? ">" : "<", qk, 1.96qse, fk, 1.96fse, qm, fm, ok ? "ok" : "FAIL"))
        end
    end
    say("")
    if isempty(failures)
        say("PASS: every pre-specified criterion holds")
    else
        say("FAIL: ", join(failures, "; "))
    end
    write(LOG, String(take!(out)))
    isempty(failures) || exit(1)
end

main()
