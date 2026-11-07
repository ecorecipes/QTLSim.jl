# Density regulation in two backends: QTLSim's kernel against the SLiM model generated
# from the same life cycle, with `Demography = :beverton_holt` (QTSim has a fixed size,
# so it cannot be the reference here). The parameter set is vignette 10's (100 loci,
# 5000 generations of stabilising selection, then a continuous move). Two-sample KS tests
# on: the equilibrium population size and genetic variance, and the time to extinction
# after the move at two rates and growth rates. Exits 1 if any p < 0.001.
#
# Run from the package root (R runs per backend and scenario, default 40; ~10 min on 8 threads):
#   julia -t auto --project=validate/slim validate/slim/compare_demography.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 40 : parse(Int, ARGS[1])

base = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 7000, PhenotypeVar = 0.022,
              Heritability = 0.32, MutationScalar = 0.004, MutationRatePerGamete = 1e-4,
              RecombPoissonMean = 3.0, InitFitnessSD = 5.0, MoveFitnessSDGen = 5000,
              MoveType = :continuous, Demography = :beverton_holt)

"A run's columns, from the kernel or from SLiM: gen, geno_var, popsize, and whether it died."
function kernel_run(p, seed)
    o = simulate(p; seed)
    (; gen = o.gen, geno_var = o.geno_var, popsize = o.popsize, extinct = o.extinct)
end
function slim_run(path, seed)
    out = read(`$SLIM -s $seed $path`, String)
    rows = [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    (; gen = [parse(Int, r[1]) for r in rows], geno_var = [parse(Float64, r[3]) for r in rows],
       popsize = [parse(Int, r[7]) for r in rows], extinct = occursin("Pop all died.", out))
end

# the same definition for both: generations after the move to the last recorded generation, plus one
after_move(r, p) = r.gen[end] + 1 - p.MoveFitnessSDGen
window(r, lo, hi) = findall(g -> lo <= g <= hi, r.gen)

function both(p, seed0)
    path = tempname() * ".slim"; write(path, slim_script(lifecycle().generation, p))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 1:2R
        j <= R ? (k[j] = kernel_run(p, seed0 + j)) : (s[j - R] = slim_run(path, seed0 + j))
    end
    rm(path; force = true)
    k, s
end

const PMIN = Ref(1.0)
function report(name, a, b)
    pv = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
    PMIN[] = min(PMIN[], pv)
    @printf("  %-34s kernel %9.4g ± %-8.2g SLiM %9.4g ± %-8.2g KS p = %.3f\n", name, mean(a), std(a), mean(b), std(b), pv)
end

println("QTLSim kernel vs generated SLiM model, Beverton-Holt demography, R = $R runs each\n")

println("── equilibrium (stabilising selection, 5000 generations), growth rate 1.5 ──")
p1 = Params(base; GrowthRate = 1.5, MoveFitnessSDGen = -1, NoGenerations = 5000)
k, s = both(p1, 10_000)
report("population size, gens 4000-4999", [mean(r.popsize[window(r, 4000, 4999)]) for r in k],
       [mean(r.popsize[window(r, 4000, 4999)]) for r in s])
report("genetic variance, gens 4000-4999", [mean(r.geno_var[window(r, 4000, 4999)]) for r in k],
       [mean(r.geno_var[window(r, 4000, 4999)]) for r in s])
@printf("  %-34s kernel %d/%d  SLiM %d/%d\n", "extinctions", count(r -> r.extinct, k), R, count(r -> r.extinct, s), R)

for (d, gr, seed0) in ((0.02, 2.0, 20_000), (0.01, 1.5, 30_000))
    println("\n── the optimum moving $d per generation from 5000, growth rate $gr ──")
    p = Params(base; MoveFitnessSD = d, GrowthRate = gr)
    kk, ss = both(p, seed0)
    @printf("  %-34s kernel %d/%d  SLiM %d/%d\n", "extinct within 2000 generations",
            count(r -> r.extinct, kk), R, count(r -> r.extinct, ss), R)
    report("generations after the move", [after_move(r, p) for r in kk], [after_move(r, p) for r in ss])
    report("population size before the move", [mean(r.popsize[window(r, 4000, 4999)]) for r in kk],
           [mean(r.popsize[window(r, 4000, 4999)]) for r in ss])
end

@printf("\nsmallest p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy")
PMIN[] < 0.001 && exit(1)
