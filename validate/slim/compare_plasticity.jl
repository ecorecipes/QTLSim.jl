# Phenotypic plasticity in two backends: QTLSim's kernel against the SLiM model generated
# from the same life cycle, with a reaction norm of fixed slope (`Plasticity`; QTSim has no
# plasticity, so it cannot be the reference here). The parameter set is
# compare_selection.jl's (100 loci, stabilising selection), with an optimum that then moves.
# Two-sample KS tests on run-level statistics, in four scenarios chosen to exercise each
# part of the Eidos splice:
#   1. b/B = 0.5, a perfect cue, the optimum moving 0.01 per generation: the lag, the genetic
#      variance and the survivors while it moves (500 generations, the lag still growing);
#   2. b = 1 with B = 2 and a cue 5 generations old, a periodic optimum (amplitude 0.3,
#      period 40) on a trend: the lag's mean and its standard deviation over the cycle (the
#      cue register CUES together with the fluctuating move);
#   3. b = 0.6 with the cost ω_b = 1 (fitness × 0.835), Beverton–Holt R = 2, the optimum moving
#      0.03 per generation: the population size before the move and the time to extinction;
#   4. two demes without migration, one plastic (b = 0.5, a step's own slope) and one not
#      (per-deme and per-step parameters): each deme's lag.
# Exits 1 if any p < 0.001. Logged in valid_plasticity.log.
#
# Run from the package root (R runs per backend and scenario, default 40):
#   julia -t 3 --project=validate/slim validate/slim/compare_plasticity.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
using Catlab: ⊗
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 40 : parse(Int, ARGS[1])

base = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 2500, PhenotypeVar = 0.022,
              Heritability = 0.32, MutationScalar = 0.004, MutationRatePerGamete = 1e-4,
              RecombPoissonMean = 3.0, InitFitnessSD = 5.0, MoveType = :continuous, MoveFitnessSDGen = 2000)

"A run's columns, from the kernel or from SLiM: gen, mean phenotype, geno_var, survivors, θ, popsize, and whether it died."
function kernel_run(p, seed)
    o = simulate(p; seed)
    (; gen = o.gen, zbar = o.pheno_mean, geno_var = o.geno_var, nsurv = o.nsurv, theta = o.theta,
       popsize = o.popsize, extinct = o.extinct)
end
function slim_rows(path, seed)
    out = read(`$SLIM -s $seed $path`, String)
    [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)], occursin("Pop all died.", out)
end
function slim_run(path, p, seed)
    rows, dead = slim_rows(path, seed)
    col(j, T) = [parse(T, r[j]) for r in rows]
    (; gen = col(1, Int), zbar = col(2, Float64), geno_var = col(3, Float64), nsurv = col(4, Int),
       theta = col(5, Float64), popsize = p.Demography === :fixed ? fill(p.PopSize, length(rows)) : col(7, Int),
       extinct = dead)
end

after_move(r, p) = r.gen[end] + 1 - p.MoveFitnessSDGen
window(r, lo, hi) = findall(g -> lo <= g <= hi, r.gen)
wmean(r, f, lo, hi) = mean(f(r)[window(r, lo, hi)])
wstd(r, f, lo, hi) = std(f(r)[window(r, lo, hi)])
lag(r) = r.theta .- r.zbar

function both(p, seed0)
    path = tempname() * ".slim"; write(path, slim_script(lifecycle().generation, p))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 1:2R
        j <= R ? (k[j] = kernel_run(p, seed0 + j)) : (s[j - R] = slim_run(path, p, seed0 + j))
    end
    rm(path; force = true)
    k, s
end

const PMIN = Ref(1.0)
function report(name, a, b)
    pv = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
    PMIN[] = min(PMIN[], pv)
    @printf("  %-40s kernel %9.4g ± %-8.2g SLiM %9.4g ± %-8.2g KS p = %.3f\n", name, mean(a), std(a), mean(b), std(b), pv)
end
extinctions(name, k, s) = @printf("  %-40s kernel %d/%d  SLiM %d/%d\n", name, count(r -> r.extinct, k), R, count(r -> r.extinct, s), R)

println("QTLSim kernel vs generated SLiM model, phenotypic plasticity, R = $R runs each, $(Threads.nthreads()) threads\n")
t0 = time()

println("── 1. b/B = 0.5, perfect cue; the optimum moving 0.01 per generation from 2000 ──")
p1 = Params(base; MoveFitnessSD = 0.01, Plasticity = 0.5)
k, s = both(p1, 110_000)
extinctions("extinctions", k, s)
report("genetic variance, gens 1500-1999", [wmean(r, x -> x.geno_var, 1500, 1999) for r in k], [wmean(r, x -> x.geno_var, 1500, 1999) for r in s])
report("θ − z̄, gens 2250-2499", [wmean(r, lag, 2250, 2499) for r in k], [wmean(r, lag, 2250, 2499) for r in s])
report("genetic variance, gens 2250-2499", [wmean(r, x -> x.geno_var, 2250, 2499) for r in k], [wmean(r, x -> x.geno_var, 2250, 2499) for r in s])
report("survivors, gens 2250-2499", [wmean(r, x -> x.nsurv, 2250, 2499) for r in k], [wmean(r, x -> x.nsurv, 2250, 2499) for r in s])

println("\n── 2. b = 1, B = 2, cue 5 generations old; periodic optimum (A = 0.3, P = 40) on a trend of 0.005 ──")
p2 = Params(base; MoveType = :periodic, MoveAmplitude = 0.3, MovePeriod = 40.0, MoveFitnessSD = 0.005,
            Plasticity = 1.0, OptimumSlope = 2.0, CueLag = 5)
k, s = both(p2, 120_000)
extinctions("extinctions", k, s)
report("θ − z̄, mean over gens 2100-2499", [wmean(r, lag, 2100, 2499) for r in k], [wmean(r, lag, 2100, 2499) for r in s])
report("θ − z̄, SD over gens 2100-2499", [wstd(r, lag, 2100, 2499) for r in k], [wstd(r, lag, 2100, 2499) for r in s])
report("survivors, gens 2100-2499", [wmean(r, x -> x.nsurv, 2100, 2499) for r in k], [wmean(r, x -> x.nsurv, 2100, 2499) for r in s])

println("\n── 3. b = 0.6 with the cost ω_b = 1; Beverton-Holt R = 2; the optimum moving 0.03 per generation ──")
p3 = Params(base; Demography = :beverton_holt, GrowthRate = 2.0, MoveFitnessSD = 0.03, Plasticity = 0.6,
            PlasticityCostSD = 1.0, NoGenerations = 4000)
k, s = both(p3, 130_000)
extinctions("extinct within 2000 generations", k, s)
report("generations after the move", [after_move(r, p3) for r in k], [after_move(r, p3) for r in s])
report("population size, gens 1500-1999", [wmean(r, x -> x.popsize, 1500, 1999) for r in k], [wmean(r, x -> x.popsize, 1500, 1999) for r in s])
report("survivors, gens 1500-1999", [wmean(r, x -> x.nsurv, 1500, 1999) for r in k], [wmean(r, x -> x.nsurv, 1500, 1999) for r in s])

println("\n── 4. two demes, no migration: deme 1 with a step's own slope 0.5 (cue 2 generations old), deme 2 without ──")
L = lifecycle(); g = L.generation
p4 = Params(base; MoveFitnessSD = 0.01)
cyc = withparams(g; Plasticity = 0.5, CueLag = 2) ⊗ g
path = tempname() * ".slim"; write(path, slim_script(cyc, p4))
kk = Vector{Any}(undef, R); ss = Vector{Any}(undef, R)
Threads.@threads :greedy for j in 1:2R
    if j <= R
        o = run_composed(p4 ⊗ demes(2), MersenneTwister(140_000 + j); cycle = cyc)
        kk[j] = [(gen = o.gen[o.deme .== d], zbar = o.pheno_mean[o.deme .== d], theta = o.theta[o.deme .== d]) for d in 1:2]
    else
        rows, _ = slim_rows(path, 140_000 + j)
        ss[j - R] = [(rs = [r for r in rows if r[7] == string(d)];
                      (gen = [parse(Int, r[1]) for r in rs], zbar = [parse(Float64, r[2]) for r in rs], theta = [parse(Float64, r[5]) for r in rs]))
                     for d in 1:2]
    end
end
rm(path; force = true)
for d in 1:2
    report("deme $d: θ − z̄, gens 2250-2499", [wmean(r[d], lag, 2250, 2499) for r in kk], [wmean(r[d], lag, 2250, 2499) for r in ss])
end

@printf("\nsmallest p = %.3g: %s  (%.1f min)\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy", (time() - t0) / 60)
PMIN[] < 0.001 && exit(1)
