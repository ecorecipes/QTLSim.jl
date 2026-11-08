# The moments lowering of a metapopulation of several traits with drift among demes
# (`run_moments(...; between = true)` for `NoTraits = 2`, ext/QTLSimCatlabExt/trait_moments.jl)
# against the kernel: the within-deme G matrix, the covariance matrix of the deme mean
# genotypes, Q_ST per trait, and the deme means under local optima along different
# directions of the two traits.
#
# Usage (from the package root; it needs Catlab, which the vignettes environment has):
#     julia -t 3 --project=vignettes validate/trait_moments_demes.jl [R]
# R replicates per scenario (default 200). Writes validate/trait_moments_demes.log.
#
# The model: k = 4 demes of N = 100 diploids, 50 loci, two traits with mutational
# correlation 0.3 (`MutCorr`) and fitness-function correlation 0.2 (`SelCorr`), each deme
# running QTSim's generation (a step `withparams(g, p; θ = θ_a)` gives deme a its optimum
# vector), then `migrate(4)` at forward rate 0.01 (island). V_E = 1 per trait, V_m = 1e-3,
# width ω = 5 (Ω = 25·S). Genomes start at zero. Scenarios:
#   T0  neutral (Drift = true): the moments' prediction is exact to first order;
#   T1  local optima θ_a = 2·(cos φ_a, sin φ_a), φ = 0°, 60°, 120°, 180°: four directions
#       in the plane of the traits, so that each deme's response is deflected by G and the
#       demes differ along both traits.
# Statistics, on the adults' genotypes at QTSim's census (after migration), per replicate
# averaged over generations 2000:10:3990 of a 4000-generation run: the within-deme sample
# covariance matrix G_W (mean over demes; entries G11, G12, G22), the sample covariance
# matrix of the 4 deme mean genotype vectors (A11, A12, A22), Q_ST of each trait,
# σ²_ii/(σ²_ii + 2 G_W,ii) with σ² = (k−1)/k·(A − G_W/N), and (T1) each deme's mean
# genotype vector. Standard errors over replicates (Q_ST's by a bootstrap). Predictions
# are the same averages of `run_moments(p; cycle, migration_rate, between = true)`.
#
# PRE-SPECIFIED CRITERIA (written before the kernel runs):
# - T0 (neutral): |z| < 4 with z = (kernel − predicted)/SE, for G11, G12, G22, A11, A12,
#   A22 and both Q_ST.
# - T1 (selection; the Gaussian approximation): tolerances relative to the scale of each
#   matrix, sqrt(X_ii X_jj) of the prediction for entry (i, j): G_W entries within 10%,
#   A entries within 15%, Q_ST within 15% of the prediction; each component of each
#   deme's mean within 0.1·|θ| = 0.2 of the moments' deterministic mean. As in
#   validate/moments_demes.jl, a statistic FAILS only if the kernel's whole 95% CI lies
#   outside the tolerance band: |kernel − pred| − 1.96 SE > tol.
# - Signature: in T1 the moments' and the kernel's deme means are deflected from the
#   direction of their optimum the same way: for each deme with an optimum off the
#   traits' axes (60°, 120°), the angle of the kernel's mean is within 10° of the
#   moments'.
# The script exits 1 if any criterion fails.

using QTLSim, Random, Statistics, Printf
import Catlab
using Catlab: ⋅, ⊗

const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200
const LOG = joinpath(@__DIR__, "trait_moments_demes.log")
const L = lifecycle()
const G = L.generation
const RATE = 0.01
const GENS = 4000
const CENSUS = collect(2000:10:GENS-10)
const K = 4
const N = 100
const PHI = [0.0, 60.0, 120.0, 180.0]

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s))
const failures = String[]

const BASE = Params(Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                    MutationRatePerGamete = 0.02, PopSize = N, NoGenerations = GENS, Sample = 10, InitFitnessSD = 5.0,
                    NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)

scenarios() = [(name = "T0 neutral", neutral = true, θ = [zeros(2) for _ in 1:K]),
               (name = "T1 local optima along four directions", neutral = false,
                θ = [2.0 .* [cosd(φ), sind(φ)] for φ in PHI])]

params(s) = s.neutral ? Params(BASE; Drift = true) : BASE
cycle(s) = foldl(⊗, [withparams(G, params(s); θ = s.θ[a]) for a in 1:K]) ⋅ L.migrate(K)

# entries (1,1), (1,2), (2,2) of a 2×2 matrix
ent(X) = [X[1, 1], X[1, 2], X[2, 2]]
qst_of(W, A) = (σ2 = (K - 1) / K .* (A .- W ./ N); [σ2[i, i] / (σ2[i, i] + 2W[i, i]) for i in 1:2])

function runmeta(s, seed)
    p = params(s)
    st = p ⊗ demes(K)
    idx = IdDict(w => i for (i, w) in enumerate(st.demes))
    nc = length(CENSUS)
    mg = [zeros(2) for _ in 1:K, _ in 1:nc]
    gw = [zeros(2, 2) for _ in 1:K, _ in 1:nc]
    rec = operator(:record, function (w, rng)
        j = searchsortedfirst(CENSUS, w.gen)
        if j <= nc && CENSUS[j] == w.gen
            gt = traitvalues(w; what = :genotype)
            a = idx[w]
            mg[a, j] = vec(mean(gt; dims = 2)); gw[a, j] = cov(gt; dims = 2)
        end
        w
    end)
    o = run_composed(st, Xoshiro(seed); cycle = cycle(s), record = false, images = Dict(:record => rec),
                     migration_rate = RATE)
    any(o.extinct_gen .>= 0) && error("a deme died in $(s.name)")
    W = mean(mean(gw[:, j]) for j in 1:nc)
    A = mean(cov(reduce(hcat, mg[:, j]); dims = 2) for j in 1:nc)
    means = [mean(mg[a, :]) for a in 1:K]
    (; W, A, means)
end

function predict(s)
    mo = run_moments(params(s); cycle = cycle(s), migration_rate = RATE, between = true, loci = :multilocus)
    sel = findall(in(CENSUS), mo.metapop.gen)
    W, A = mean(mo.metapop.within[sel]), mean(mo.metapop.var_among[sel])
    rows = findall(in(CENSUS), mo.gen)
    means = [mean(mo.trait_mean[r] for r in rows if mo.deme[r] == a) for a in 1:K]
    (; W, A, qst = qst_of(W, A), means)
end

function main()
    say("QTLSim validate/trait_moments_demes.jl")
    say("Julia threads: ", Threads.nthreads(), "; replicates per scenario R = ", R)
    say(@sprintf("k = %d demes of N = %d, %d loci, 2 traits (MutCorr %.1f, SelCorr %.1f), V_E = %g, V_m = %g, ω = %g, forward rate %g",
                 K, N, BASE.NoQTL, BASE.MutCorr, BASE.SelCorr, envVar(BASE), muteVar(BASE), BASE.InitFitnessSD, RATE))
    say("census generations $(first(CENSUS)):10:$(last(CENSUS)) of $GENS")
    sc = scenarios()
    jobs = [(i, r) for i in eachindex(sc) for r in 1:R]
    res = Vector{Any}(undef, length(jobs))
    t0 = time()
    Threads.@threads :greedy for j in eachindex(jobs)
        i, r = jobs[j]
        res[j] = runmeta(sc[i], 9_000_000 + 10_000 * i + r)
    end
    say(@sprintf("%d kernel runs in %.0f s", length(jobs), time() - t0))
    brng = Xoshiro(20261107)
    for (i, s) in enumerate(sc)
        rs = [res[j] for j in eachindex(jobs) if jobs[j][1] == i]
        n = length(rs)
        pr = predict(s)
        say("\n", s.name, s.neutral ? " (Drift = true)" : "", "; ", n, " replicates")
        say(@sprintf("  %-22s %10s %10s %9s %21s %9s %8s", "statistic", "moments", "kernel", "SE", "95% CI",
                     s.neutral ? "z" : "dev/scale", "verdict"))
        function check(label, pred, sims, scale, tol)
            m, se = mean(sims), std(sims) / sqrt(n)
            d = m - pred
            ok = s.neutral ? abs(d / se) < 4 : abs(d) - 1.96se <= tol * scale
            stat = s.neutral ? @sprintf("%8.2f", d / se) : @sprintf("%+8.1f%%", 100d / scale)
            ok || push!(failures, "$(s.name) $label")
            say(@sprintf("  %-22s %10.5g %10.5g %9.2g [%9.4g, %9.4g] %9s %8s", label, pred, m, se, m - 1.96se, m + 1.96se,
                         stat, ok ? "ok" : "FAIL"))
        end
        for (X, Xs, tol, nm) in ((pr.W, [x.W for x in rs], 0.10, "G_W"), (pr.A, [x.A for x in rs], 0.15, "A"))
            for (lab, (a, b)) in zip(("11", "12", "22"), ((1, 1), (1, 2), (2, 2)))
                check(nm * lab, X[a, b], [x[a, b] for x in Xs], sqrt(X[a, a] * X[b, b]), tol)
            end
        end
        # Q_ST per trait, with a bootstrap SE
        for t in 1:2
            boot = [(ix = rand(brng, 1:n, n); qst_of(mean(rs[ix[q]].W for q in 1:n), mean(rs[ix[q]].A for q in 1:n))[t]) for _ in 1:2000]
            qk = qst_of(mean(x.W for x in rs), mean(x.A for x in rs))[t]
            se = std(boot)
            d = qk - pr.qst[t]
            ok = s.neutral ? abs(d / se) < 4 : abs(d) - 1.96se <= 0.15 * pr.qst[t]
            ok || push!(failures, "$(s.name) Q_ST trait $t")
            say(@sprintf("  %-22s %10.5g %10.5g %9.2g [%9.4g, %9.4g] %9s %8s", "Q_ST trait $t", pr.qst[t], qk, se,
                         qk - 1.96se, qk + 1.96se, s.neutral ? @sprintf("%8.2f", d / se) : @sprintf("%+8.1f%%", 100d / pr.qst[t]),
                         ok ? "ok" : "FAIL"))
        end
        s.neutral && continue
        say("  deme means (moments vs kernel), tolerance 0.2 per component:")
        for a in 1:K
            for t in 1:2
                check(@sprintf("deme %d (φ=%3.0f°) z̄_%d", a, PHI[a], t), pr.means[a][t], [x.means[a][t] for x in rs], 2.0, 0.1)
            end
        end
        say("  signature: deflection of the deme means from their optima (angles in degrees):")
        for a in 1:K
            km = mean(x.means[a] for x in rs)
            angm, angk = atand(pr.means[a][2], pr.means[a][1]), atand(km[2], km[1])
            offaxis = PHI[a] in (60.0, 120.0)
            ok = !offaxis || abs(angk - angm) <= 10
            ok || push!(failures, "$(s.name) deflection deme $a")
            say(@sprintf("    deme %d: optimum %5.1f°, moments mean %6.1f°, kernel mean %6.1f°%s", a, PHI[a], angm, angk,
                         offaxis ? (ok ? "  ok" : "  FAIL") : ""))
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
