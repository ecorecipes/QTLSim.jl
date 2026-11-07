# Migration in two backends: QTLSim's kernel against the SLiM model generated from the same
# life cycle, `g ⊗ g ⊗ g ⋅ migrate(3)`, on three neutral demes (QTSim has no demes, so it
# cannot be the reference here). Both visit the demes in a random order at each migration
# step (Lean `randomOrderSweep`, `lean/CategoricalQuantGen/Migration.lean`). Three
# migration models: a scalar rate (the island model, `migRate`), the island matrix and a
# line of stepping stones (`MIG`). Two-sample KS tests on, at the last generation: the
# squared differences between deme means, end demes and neighbours, and the genetic
# variance within a deme. Exits 1 if any p < 0.001.
#
# Run from the package root (R runs per backend and model, default 300; about an hour on 3 threads):
#   julia -t auto --project=validate/slim validate/slim/compare_migration.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
using Catlab: ⋅, ⊗
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 300 : parse(Int, ARGS[1])

neutral = Params(PopSize = 200, NoQTL = 10, RecombPoissonMean = 1.0, Ploidy = 2, NoGenerations = 300,
                 PhenotypeVar = 0.02, Heritability = 0.5, MutationScalar = 0.02,
                 MutationRatePerGamete = 0.01, Drift = true)
L = lifecycle(); g = L.generation
const CYCLE = g ⊗ g ⊗ g ⋅ L.migrate(3)
const LAST = neutral.NoGenerations - 1

"Deme means and genetic variances at the last generation, from either backend."
function at_end(gen, pheno_mean, geno_var, deme)
    i = [findfirst((gen .== LAST) .& (deme .== d)) for d in 1:3]
    (; mean = pheno_mean[i], vg = geno_var[i])
end
function kernel_run(rate, seed)
    o = run_composed(neutral ⊗ demes(3), Xoshiro(seed); cycle = CYCLE, migration_rate = rate)
    at_end(o.gen, o.pheno_mean, o.geno_var, o.deme)
end
function slim_run(path, seed)
    out = read(`$SLIM -s $seed $path`, String)
    rows = [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    at_end([parse(Int, r[1]) for r in rows], [parse(Float64, r[2]) for r in rows],
           [parse(Float64, r[3]) for r in rows], [parse(Int, r[7]) for r in rows])
end

function both(rate, seed0)
    path = tempname() * ".slim"; write(path, slim_script(CYCLE, neutral; migration_rate = rate))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 1:2R
        j <= R ? (k[j] = kernel_run(rate, seed0 + j)) : (s[j - R] = slim_run(path, seed0 + j))
    end
    rm(path; force = true)
    k, s
end

const PMIN = Ref(1.0)
function report(name, a, b)
    pv = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
    PMIN[] = min(PMIN[], pv)
    @printf("  %-40s kernel %9.4g  SLiM %9.4g  KS p = %.3f\n", name, mean(a), mean(b), pv)
end

ends(r) = (r.mean[1] - r.mean[3])^2
neighbours(r) = ((r.mean[1] - r.mean[2])^2 + (r.mean[2] - r.mean[3])^2) / 2
within(r) = mean(r.vg)

println("QTLSim kernel vs generated SLiM model, three neutral demes of $(neutral.PopSize), ",
        "$(neutral.NoGenerations) generations, R = $R runs each\n")
for (name, rate, seed0) in (("island, scalar rate 0.05", 0.05, 10_000),
                            ("island, island_migration(3, 0.05)", island_migration(3, 0.05), 20_000),
                            ("line, stepping_stone(3, 0.05; circular = false)",
                             stepping_stone(3, 0.05; circular = false), 30_000))
    println("── ", name, " ──")
    k, s = both(rate, seed0)
    report("squared difference, end demes", ends.(k), ends.(s))
    report("squared difference, neighbours", neighbours.(k), neighbours.(s))
    report("genetic variance within a deme", within.(k), within.(s))
    @printf("  %-40s kernel %9.3f  SLiM %9.3f\n\n", "ratio of the mean squared differences",
            mean(ends.(k)) / mean(neighbours.(k)), mean(ends.(s)) / mean(neighbours.(s)))
end

@printf("smallest p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy")
PMIN[] < 0.001 && exit(1)
