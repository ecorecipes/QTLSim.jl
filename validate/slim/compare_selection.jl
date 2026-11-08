# Fecundity selection in two backends: QTLSim's kernel against the SLiM model generated
# from the same life cycle, with `Selection = :fecundity` (QTSim has only viability
# selection, so it cannot be the reference here). The parameter set is vignette 10's and
# compare_demography.jl's (100 loci, stabilising selection, then a continuous move).
# Two-sample KS tests, at a fixed size and under Beverton–Holt density regulation, on: the
# equilibrium genetic variance and lag, the lag and genetic variance while the optimum
# moves, the equilibrium population size, and the time to extinction after the move.
# Exits 1 if any p < 0.001. Logged in valid_selection.log.
#
# Run from the package root (R runs per backend and scenario, default 40; 80 min on 1 thread;
# CI runs 8 on push, a gross-breakage check, and 40 weekly):
#   julia -t 4 --project=validate/slim validate/slim/compare_selection.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 40 : parse(Int, ARGS[1])

base = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 4000, PhenotypeVar = 0.022,
              Heritability = 0.32, MutationScalar = 0.004, MutationRatePerGamete = 1e-4,
              RecombPoissonMean = 3.0, InitFitnessSD = 5.0, MoveType = :continuous, Selection = :fecundity)

"A run's columns, from the kernel or from SLiM: gen, mean phenotype, geno_var, θ, popsize, and whether it died."
function kernel_run(p, seed)
    o = simulate(p; seed)
    (; gen = o.gen, zbar = o.pheno_mean, geno_var = o.geno_var, nsurv = o.nsurv, theta = o.theta,
       popsize = o.popsize, extinct = o.extinct)
end
function slim_run(path, p, seed)
    out = read(`$SLIM -s $seed $path`, String)
    rows = [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    col(j, T) = [parse(T, r[j]) for r in rows]
    (; gen = col(1, Int), zbar = col(2, Float64), geno_var = col(3, Float64), nsurv = col(4, Int),
       theta = col(5, Float64), popsize = p.Demography === :fixed ? fill(p.PopSize, length(rows)) : col(7, Int),
       extinct = occursin("Pop all died.", out))
end

after_move(r, p) = r.gen[end] + 1 - p.MoveFitnessSDGen
window(r, lo, hi) = findall(g -> lo <= g <= hi, r.gen)
wmean(r, f, lo, hi) = mean(f(r)[window(r, lo, hi)])
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
    @printf("  %-36s kernel %9.4g ± %-8.2g SLiM %9.4g ± %-8.2g KS p = %.3f\n", name, mean(a), std(a), mean(b), std(b), pv)
end
extinctions(name, k, s) = @printf("  %-36s kernel %d/%d  SLiM %d/%d\n", name, count(r -> r.extinct, k), R, count(r -> r.extinct, s), R)

println("QTLSim kernel vs generated SLiM model, Selection = :fecundity, R = $R runs each, $(Threads.nthreads()) threads\n")
t0 = time()

println("── fixed size: equilibrium, then the optimum moving 0.01 per generation from 3000 ──")
p1 = Params(base; MoveFitnessSDGen = 3000, MoveFitnessSD = 0.01)
k, s = both(p1, 40_000)
extinctions("extinctions", k, s)
allsurv = (count(r -> all(r.nsurv .== r.popsize), k), count(r -> all(r.nsurv .== r.popsize), s))
@printf("  %-36s kernel %d/%d  SLiM %d/%d\n", "runs in which everyone survives", allsurv[1], R, allsurv[2], R)
allsurv == (R, R) || (PMIN[] = 0.0)        # fecundity selection kills nobody: a discrepancy if it does
report("genetic variance, gens 2000-2999", [wmean(r, x -> x.geno_var, 2000, 2999) for r in k],
       [wmean(r, x -> x.geno_var, 2000, 2999) for r in s])
report("θ − z̄, gens 2000-2999", [wmean(r, lag, 2000, 2999) for r in k], [wmean(r, lag, 2000, 2999) for r in s])
report("θ − z̄, gens 3500-3999", [wmean(r, lag, 3500, 3999) for r in k], [wmean(r, lag, 3500, 3999) for r in s])
report("genetic variance, gens 3500-3999", [wmean(r, x -> x.geno_var, 3500, 3999) for r in k],
       [wmean(r, x -> x.geno_var, 3500, 3999) for r in s])

println("\n── Beverton-Holt, growth rate 1.5: equilibrium (stabilising selection, 3000 generations) ──")
p2 = Params(base; Demography = :beverton_holt, GrowthRate = 1.5, NoGenerations = 3000)
k, s = both(p2, 50_000)
extinctions("extinctions", k, s)
report("population size, gens 2000-2999", [wmean(r, x -> x.popsize, 2000, 2999) for r in k],
       [wmean(r, x -> x.popsize, 2000, 2999) for r in s])
report("genetic variance, gens 2000-2999", [wmean(r, x -> x.geno_var, 2000, 2999) for r in k],
       [wmean(r, x -> x.geno_var, 2000, 2999) for r in s])

for (d, gr, seed0) in ((0.02, 2.0, 60_000), (0.01, 1.5, 70_000))
    println("\n── Beverton-Holt: the optimum moving $d per generation from 3000, growth rate $gr ──")
    p = Params(base; Demography = :beverton_holt, GrowthRate = gr, MoveFitnessSDGen = 3000, MoveFitnessSD = d,
               NoGenerations = 5000)
    kk, ss = both(p, seed0)
    extinctions("extinct within 2000 generations", kk, ss)
    report("generations after the move", [after_move(r, p) for r in kk], [after_move(r, p) for r in ss])
    report("population size before the move", [wmean(r, x -> x.popsize, 2000, 2999) for r in kk],
           [wmean(r, x -> x.popsize, 2000, 2999) for r in ss])
end

@printf("\nsmallest p = %.3g: %s  (%.1f min)\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy", (time() - t0) / 60)
PMIN[] < 0.001 && exit(1)
