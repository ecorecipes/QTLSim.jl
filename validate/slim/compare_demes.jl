# Demes with parameters of their own in two backends: QTLSim's kernel against the SLiM
# model generated from the same life cycle (QTSim has no demes, so it cannot be the
# reference here). Three checks:
#
# A. three demes with their own Params, [p1, p2, p3] ⊗ demes(3): widths ω = 0.15, 0.3 and
#    0.6 and sizes 400, 500 and 600, coupled by migration at rate 0.01;
# B. two demes whose steps carry the parameters, withparams(g; θ = 0.1, ω = 0.2) ⊗
#    withparams(g; θ = -0.1, ω = 0.4) ⋅ migrate(2; rate = 0.02): fixed optima of their
#    own, so the demes diverge, and migration holds them together;
# C. run for run in SLiM, a life cycle whose steps carry a continuous move and a mutation
#    scalar is the same model as the plain life cycle with those values in Params.
#
# Two-sample KS tests in A and B on each deme's mean genetic variance and mean phenotype
# over the last 500 generations, and on the difference of the deme means in B. Exits 1 if
# any p < 0.001 or any run of C differs.
#
# Run from the package root (R runs per backend and scenario, default 50; ~3 min on 3 threads):
#   julia -t 3 --project=validate/slim validate/slim/compare_demes.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
using Catlab: ⋅, ⊗
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 50 : parse(Int, ARGS[1])

L = lifecycle(); g = L.generation
base = Params(PopSize = 500, NoQTL = 1, Ploidy = 2, NoGenerations = 1000, PhenotypeVar = 0.02,
              Heritability = 0.5, MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)

"Per deme, the mean phenotype and mean genetic variance over generations ≥ `from`."
function summarise(gen, deme, zbar, vg, k; from = 500)
    map(1:k) do d
        i = findall(j -> deme[j] == d && gen[j] >= from, eachindex(gen))
        (zbar = mean(zbar[i]), vg = mean(vg[i]))
    end
end
function kernel_run(ps, cycle, seed)
    o = run_composed(ps ⊗ demes(length(ps)), MersenneTwister(seed); cycle)
    summarise(o.gen, o.deme, o.pheno_mean, o.geno_var, length(ps))
end
slim_out(path, seed) = read(`$SLIM -s $seed $path`, String)
function slim_run(path, seed, k)
    rows = [split(l, '\t') for l in split(slim_out(path, seed), '\n') if occursin(r"^\d+\t", l)]
    summarise([parse(Int, r[1]) for r in rows], [parse(Int, r[7]) for r in rows],
              [parse(Float64, r[2]) for r in rows], [parse(Float64, r[3]) for r in rows], k)
end

function both(ps, cycle, seed0)
    path = tempname() * ".slim"; write(path, slim_script(cycle, ps))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 1:2R
        j <= R ? (k[j] = kernel_run(ps, cycle, seed0 + j)) : (s[j - R] = slim_run(path, seed0 + j, length(ps)))
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

println("QTLSim kernel vs generated SLiM model, demes with parameters of their own, R = $R runs each\n")

println("── A: three demes with their own Params (ω 0.15, 0.3, 0.6; N 400, 500, 600), migration 0.01 ──")
psA = [Params(base; InitFitnessSD = 0.15, PopSize = 400), base, Params(base; InitFitnessSD = 0.6, PopSize = 600)]
cycleA = g ⊗ g ⊗ g ⋅ L.migrate(3; rate = 0.01)
k, s = both(psA, cycleA, 10_000)
for d in 1:3
    report("deme $d genetic variance", [r[d].vg for r in k], [r[d].vg for r in s])
    report("deme $d mean phenotype", [r[d].zbar for r in k], [r[d].zbar for r in s])
end

println("\n── B: two demes, steps with θ = ±0.1 and ω = 0.2, 0.4, migration 0.02 ──")
cycleB = withparams(g; θ = 0.1, ω = 0.2) ⊗ withparams(g; θ = -0.1, ω = 0.4) ⋅ L.migrate(2; rate = 0.02)
k, s = both([base, base], cycleB, 20_000)
for d in 1:2
    report("deme $d genetic variance", [r[d].vg for r in k], [r[d].vg for r in s])
    report("deme $d mean phenotype", [r[d].zbar for r in k], [r[d].zbar for r in s])
end
report("difference of the deme means", [r[1].zbar - r[2].zbar for r in k], [r[1].zbar - r[2].zbar for r in s])

println("\n── C: steps with parameters = Params with those values, run for run in SLiM ──")
pC = Params(base; NoGenerations = 400)
over = (MoveType = :continuous, MoveFitnessSDGen = 200, MoveFitnessSD = 0.002, MutationScalar = 0.03)
a = tempname() * ".slim"; write(a, slim_script(withparams(g; over...), pC))
b = tempname() * ".slim"; write(b, slim_script(g, Params(pC; over...)))
same = [slim_out(a, seed) == slim_out(b, seed) for seed in 1:10]
rm(a); rm(b)
@printf("  identical output in %d of %d runs\n", count(same), length(same))

@printf("\nsmallest p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy")
(PMIN[] < 0.001 || !all(same)) && exit(1)
