# Several traits in two backends: QTLSim's kernel against the SLiM model generated from
# the same life cycle with the Eidos images of several traits (QTSim has one trait, so it
# cannot be the reference). Two traits, M and Ω not sharing eigenvectors (mutational
# covariance m·V_E·[1 0.4; 0.4 0.5], fitness correlation −0.3). Two-sample KS tests on:
#   A, one diploid locus, stabilising selection: the G matrix's entries, the genetic
#      correlation and the number of survivors, averaged over generations 1000–1999;
#   B, ten loci, an optimum moving by 0.001 per generation along (1, 2) from generation
#      1000: the lag vector's components and the G matrix's entries averaged over
#      generations 1500–1999, and the mean of trait 2 at the end.
# Pre-specified: exits 1 if any p < 0.001. Writes valid_traits.log.
#
# Run from the package root (R runs per backend and scenario, default 40; at most 3 threads here):
#   julia -t 3 --project=validate/slim validate/slim/compare_traits.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 40 : parse(Int, ARGS[1])
const LOG = open(joinpath(@__DIR__, "valid_traits.log"), "w")
say(s) = (print(stdout, s); print(LOG, s); flush(LOG))

base = Params(PopSize = 500, NoQTL = 1, Ploidy = 2, NoGenerations = 2000, PhenotypeVar = 0.02, Heritability = 0.5,
              MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)
pA = Params(base; NoTraits = 2, MutCov = muteVar(base) .* [1.0 0.4; 0.4 0.5], SelCorr = -0.3)
pB = Params(pA; NoQTL = 10, RecombPoissonMean = 2.0, MoveFitnessSDGen = 1000, MoveType = :continuous,
            MoveFitnessSD = 0.001, MoveDirection = [1.0, 2.0])

"A run's columns, from the kernel or from SLiM: gen, nsurv, the trait means, G's upper triangle and the optimum."
function kernel_run(p, seed)
    o = simulate(p; seed)
    (; gen = o.gen, nsurv = o.nsurv, z = o.trait_mean, G = [[g[1, 1], g[1, 2], g[2, 2]] for g in o.gmatrix],
       theta = [o.trait_mean[i] + o.trait_lag[i] for i in eachindex(o.gen)], extinct = o.extinct)
end
function slim_run(path, seed)
    out = read(`$SLIM -s $seed $path`, String)
    rows = [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    f(r, i) = parse(Float64, r[i])
    (; gen = [parse(Int, r[1]) for r in rows], nsurv = [parse(Int, r[4]) for r in rows],
       z = [[f(r, 7), f(r, 8)] for r in rows], G = [[f(r, 9), f(r, 10), f(r, 11)] for r in rows],
       theta = [[f(r, 12), f(r, 13)] for r in rows], extinct = occursin("Pop all died.", out))
end

window(r, lo, hi) = findall(g -> lo <= g <= hi, r.gen)
wmean(v, r, lo, hi) = mean(v[window(r, lo, hi)])

function both(p, seed0)
    path = tempname() * ".slim"; write(path, slim_script(lifecycle().generation, p))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 1:2R
        j <= R ? (k[j] = kernel_run(p, seed0 + j)) : (s[j - R] = slim_run(path, seed0 + j))
    end
    rm(path; force = true)
    any(r -> r.extinct, k) || any(r -> r.extinct, s) ? @warn("extinctions in a scenario meant to survive") : nothing
    k, s
end

const PMIN = Ref(1.0)
function report(name, a, b)
    pv = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
    PMIN[] = min(PMIN[], pv)
    say(@sprintf("  %-38s kernel %10.4g ± %-9.2g SLiM %10.4g ± %-9.2g KS p = %.3f\n", name, mean(a), std(a), mean(b), std(b), pv))
end

say("QTLSim kernel vs generated SLiM model, two traits, R = $R runs each, $(Threads.nthreads()) threads\n")
say(@sprintf("M = m·V_E·[1 0.4; 0.4 0.5], Ω = ω²·[1 -0.3; -0.3 1], ω = %.2f, V_E = %.3g\n\n", initFitSD(pA), envVar(pA)))

say("── A: one diploid locus, stabilising selection, generations 1000-1999 ──\n")
t0 = time(); k, s = both(pA, 10_000); t = time() - t0
for (i, nm) in enumerate(("G11", "G12", "G22"))
    report(nm, [wmean(getindex.(r.G, i), r, 1000, 1999) for r in k], [wmean(getindex.(r.G, i), r, 1000, 1999) for r in s])
end
corr(r) = wmean([g[2] / sqrt(g[1] * g[3]) for g in r.G], r, 1000, 1999)
report("genetic correlation", corr.(k), corr.(s))
report("survivors", [wmean(r.nsurv, r, 1000, 1999) for r in k], [wmean(r.nsurv, r, 1000, 1999) for r in s])
say(@sprintf("  (%.0f s; theory: gaussian_G for one locus, per-locus sum gaussian_G_loci G11 = %.4g)\n\n",
             t, gaussian_G_loci(pA)[1, 1]))

say("── B: ten loci, the optimum moving 0.001 per generation along (1, 2) from 1000, generations 1500-1999 ──\n")
t0 = time(); k, s = both(pB, 20_000); t = time() - t0
lag(r, i) = wmean([r.theta[j][i] - r.z[j][i] for j in eachindex(r.gen)], r, 1500, 1999)
report("lag, trait 1", [lag(r, 1) for r in k], [lag(r, 1) for r in s])
report("lag, trait 2", [lag(r, 2) for r in k], [lag(r, 2) for r in s])
for (i, nm) in enumerate(("G11", "G12", "G22"))
    report(nm, [wmean(getindex.(r.G, i), r, 1500, 1999) for r in k], [wmean(getindex.(r.G, i), r, 1500, 1999) for r in s])
end
report("mean of trait 2, last generation", [r.z[end][2] for r in k], [r.z[end][2] for r in s])
Gk = mean([wmean(r.G, r, 1500, 1999) for r in k])
λpred = steady_lag_vector(pB; G = [Gk[1] Gk[2]; Gk[2] Gk[3]])
say(@sprintf("  (%.0f s; steady_lag_vector with the kernel's G: (%.4f, %.4f))\n", t, λpred...))

say(@sprintf("\nsmallest p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy"))
close(LOG)
PMIN[] < 0.001 && exit(1)
