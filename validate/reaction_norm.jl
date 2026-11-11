# A heritable reaction norm under a sustained move: the kernel against the moments lowering
# and the quasi-steady lag (ROADMAP F2, stage 2). `ReactionNorm = true`, `NoTraits = 2`:
# elevation a and slope b heritable, z = a + b ε + e, ε = θ/B with B = 1 and a perfect cue.
#
# Written before the first run of this script (the same parameter set had been run at 8
# replicates in test/runtests.jl, whose tolerances these are, tightened where stated).
#
# Regime: 1000 diploids, 50 loci, μ = 0.01 per locus (the Gaussian regime, where the
# multilocus moments are within a few per cent of the kernel for one trait), V_P = 1,
# h² = 0.5, ω = 1, u = 3, CrossoverModel = :qtsim; mutational covariance diag(V_m, V_m/4) for
# (elevation, slope), V_m = m·V_E; initial slope 0.2. 1000 generations of stabilising
# selection at ε = 0, then the optimum moves 0.02 per generation for 1500 generations.
# Scenario A: no cost (the slope is neutral before the move). Scenario B: cost ω_b = 1.
# 40 runs each; the moments lowering with drift (`erosion = true`), loci = :multilocus.
#
# Criteria, at 250, 500, 1000 and 1500 generations after the move (means over the runs;
# lags and G averaged over the 50 generations up to each time):
#   C1  mean slope: |kernel − moments| ≤ 0.1 (the slope goes from 0.2 toward B = 1).
#   C2  lag θ − z̄: kernel/quasi-steady within 0.75–1.25, the quasi-steady lag being
#       reaction_norm_lag at each run's own G, mean slope and environment (no moments
#       involved; scenario A only, the formula's cost term being first-order).
#   C3  lag: |kernel − moments| ≤ 0.1 phenotypic SD of the kernel at the move.
#   C4  the slope's genetic variance G_bb: kernel/moments within a factor 2.
#   C5  (once) in scenario A the slope ends above 0.8 in both; in scenario B it ends lower
#       than in scenario A in both, by at least 0.05.
# Failures stay in the log (validate/reaction_norm.log).
#
# Run from the package root:  julia -t 3 --project=vignettes validate/reaction_norm.jl [SEED0]

using QTLSim, Statistics, Printf, Random
import Catlab

const SEED0 = isempty(ARGS) ? 31_000 : parse(Int, ARGS[1])
const NRUN = 40
const T0 = 1000
const TIMES = (250, 500, 1000, 1500)
base = Params(NoGenerations = T0 + 1500, PopSize = 1000, Ploidy = 2, NoQTL = 50, MutationRatePerGamete = 0.01,
              MoveType = :continuous, MoveFitnessSDGen = T0, MoveFitnessSD = 0.02, RecombPoissonMean = 3.0,
              InitFitnessSD = 1.0, PhenotypeVar = 1.0, Heritability = 0.5)
Vm = muteVar(base)
norm(cost) = Params(base; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2, MutCov = [Vm 0.0; 0.0 Vm / 4], PlasticityCostSD = cost)

const NPASS = Ref(0); const NALL = Ref(0)
verdict(ok) = (NALL[] += 1; NPASS[] += ok; ok ? "PASS" : "FAIL")
win(t) = T0 + t - 49:T0 + t        # record rows are generation + 1

function scenario(name, p, seedoff)
    println("\n── $name ──")
    rs = Vector{Any}(undef, NRUN)
    Threads.@threads :greedy for i in 1:NRUN
        rs[i] = simulate(p; seed = SEED0 + seedoff + i)
    end
    m = run_moments(p; erosion = true, loci = :multilocus)
    sdP = sqrt(mean(r.pheno_var[T0] for r in rs))
    @printf("  at the move: kernel mean slope %.3f (SD among runs %.3f), G = [%.3f %.3f; %.3f %.3f]; moments slope %.3f, G = [%.3f %.3f; %.3f %.3f]; phenotypic SD %.3f\n",
            mean(r.trait_mean[T0][2] for r in rs), std([r.trait_mean[T0][2] for r in rs]),
            mean(r.gmatrix[T0][1, 1] for r in rs), mean(r.gmatrix[T0][1, 2] for r in rs), mean(r.gmatrix[T0][2, 1] for r in rs),
            mean(r.gmatrix[T0][2, 2] for r in rs), m.trait_mean[T0][2], m.gmatrix[T0]..., sdP)
    ends = Float64[]
    for t in TIMES
        w = win(t)
        bk, bm = mean(r.trait_mean[T0 + t][2] for r in rs), m.trait_mean[T0 + t][2]
        lagk = mean(mean(x[1] for x in r.trait_lag[w]) for r in rs)
        lagm = mean(x[1] for x in m.trait_lag[w])
        gk = mean(mean(G[2, 2] for G in r.gmatrix[w]) for r in rs)
        gm = mean(G[2, 2] for G in m.gmatrix[w])
        @printf("  +%4d  slope kernel %.3f moments %.3f  C1 %s | lag kernel %.4f moments %.4f  C3 %s | G_bb kernel %.4g moments %.4g  C4 %s\n",
                t, bk, bm, verdict(abs(bk - bm) <= 0.1), lagk, lagm, verdict(abs(lagk - lagm) <= 0.1sdP), gk, gm, verdict(0.5 <= gk / gm <= 2))
        if p.PlasticityCostSD == Inf
            q = mean(mean(reaction_norm_lag(r.gmatrix[g] .* ((p.PopSize - 1) / p.PopSize), r.theta[g]; k = p.MoveFitnessSD,
                                            bbar = r.trait_mean[g][2], ω = 1.0, VE = envVar(p)) for g in w) for r in rs)
            @printf("         quasi-steady lag at the kernel's own G, slope and environment %.4f; kernel/it %.3f  C2 %s\n",
                    q, lagk / q, verdict(0.75 <= lagk / q <= 1.25))
        end
        t == TIMES[end] && append!(ends, (bk, bm))
    end
    surv = count(r -> !r.extinct, rs)
    println("  runs alive at the end: $surv/$NRUN")
    ends
end

println("A heritable reaction norm, kernel vs moments and theory; $NRUN runs, seeds from $SEED0, $(Threads.nthreads()) threads")
t0 = time()
a = scenario("A: no cost", norm(Inf), 0)
b = scenario("B: cost ω_b = 1", norm(1.0), 1000)
println("\n── C5 ──")
@printf("  scenario A: final slope kernel %.3f, moments %.3f, both above 0.8  %s\n", a..., verdict(all(>(0.8), a)))
@printf("  scenario B: final slope kernel %.3f, moments %.3f, each at least 0.05 below scenario A's  %s\n", b...,
        verdict(b[1] <= a[1] - 0.05 && b[2] <= a[2] - 0.05))
@printf("\n%d of %d criteria pass (%.1f min)\n", NPASS[], NALL[], (time() - t0) / 60)
NPASS[] == NALL[] || exit(1)
