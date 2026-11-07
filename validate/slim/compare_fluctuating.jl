# Fluctuating optima in two backends: QTLSim's kernel against the SLiM model generated
# from the same life cycle, with `MoveType = :periodic` and `:red_noise` (QTSim has
# neither, so it cannot be the reference here). One locus, scenario C's selection
# (ω = 0.3, V_E = 0.01, m·V_E = 2e-4), 1000 generations of stabilising selection, then
# 2000 of a fluctuating optimum with no trend. Two-sample KS tests, per run, over
# generations 1500-2999, on: the mean squared lag θ − z̄, the genetic variance, the
# number of survivors and, for red noise, the optimum's SD and lag-1 autocorrelation.
# The theory column is `expected_lag2` at the kernel's mean genetic variance (it ignores
# drift and the fluctuation of V_G). Exits 1 if any p < 0.001.
#
# Run from the package root (R runs per backend and scenario, default 40; ~3 min on 3 threads):
#   julia -t 3 --project=validate/slim validate/slim/compare_fluctuating.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 40 : parse(Int, ARGS[1])

base = Params(PopSize = 1000, Ploidy = 2, NoQTL = 1, NoGenerations = 3000, PhenotypeVar = 0.02,
              Heritability = 0.5, MutationScalar = 0.02, MutationRatePerGamete = 0.5,
              InitFitnessSD = 0.3, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.0)

"A run's columns, from the kernel or from SLiM: gen, pheno_mean, geno_var, nsurv, theta, and whether it died."
function kernel_run(p, seed)
    o = simulate(p; seed)
    (; o.gen, o.pheno_mean, o.geno_var, o.nsurv, o.theta, o.extinct)
end
function slim_run(path, seed)
    out = read(`$SLIM -s $seed $path`, String)
    rows = [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    col(i, T) = [parse(T, r[i]) for r in rows]
    (; gen = col(1, Int), pheno_mean = col(2, Float64), geno_var = col(3, Float64), nsurv = col(4, Int),
       theta = col(5, Float64), extinct = occursin("Pop all died.", out))
end

function both(p, seed0)
    path = tempname() * ".slim"; write(path, slim_script(lifecycle().generation, p))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 1:2R
        j <= R ? (k[j] = kernel_run(p, seed0 + j)) : (s[j - R] = slim_run(path, seed0 + j))
    end
    rm(path; force = true)
    k, s
end

window(r) = findall(g -> 1500 <= g <= 2999, r.gen)
lag2(r) = (w = window(r); mean((r.theta[w] .- r.pheno_mean[w]) .^ 2))
stat(f, runs) = [f(r) for r in runs if !r.extinct]

const PMIN = Ref(1.0)
function report(name, a, b)
    pv = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
    PMIN[] = min(PMIN[], pv)
    @printf("  %-30s kernel %9.4g ± %-8.2g SLiM %9.4g ± %-8.2g KS p = %.3f\n", name, mean(a), std(a), mean(b), std(b), pv)
end

println("QTLSim kernel vs generated SLiM model, fluctuating optima, R = $R runs each\n")

scenarios = (("periodic, A = 0.1, P = 200", Params(base; MoveType = :periodic, MoveAmplitude = 0.1, MovePeriod = 200.0), 40_000),
             ("periodic, A = 0.15, P = 50", Params(base; MoveType = :periodic, MoveAmplitude = 0.15, MovePeriod = 50.0), 50_000),
             ("red noise, σ = 0.1, ρ = 0.9", Params(base; MoveType = :red_noise, MoveNoiseSD = 0.1, MoveAutocorr = 0.9), 60_000),
             ("red noise, σ = 0.1, ρ = 0", Params(base; MoveType = :red_noise, MoveNoiseSD = 0.1, MoveAutocorr = 0.0), 70_000))

for (name, p, seed0) in scenarios
    println("── $name ──")
    k, s = both(p, seed0)
    @printf("  %-30s kernel %d/%d  SLiM %d/%d\n", "extinctions", count(r -> r.extinct, k), R, count(r -> r.extinct, s), R)
    report("mean squared lag", stat(lag2, k), stat(lag2, s))
    report("genetic variance", stat(r -> mean(r.geno_var[window(r)]), k), stat(r -> mean(r.geno_var[window(r)]), s))
    report("survivors", stat(r -> mean(r.nsurv[window(r)]), k), stat(r -> mean(r.nsurv[window(r)]), s))
    if p.MoveType === :red_noise
        report("SD of θ", stat(r -> std(r.theta[window(r)]), k), stat(r -> std(r.theta[window(r)]), s))
        report("lag-1 autocorrelation of θ", stat(r -> (θ = r.theta[window(r)]; cor(θ[1:end-1], θ[2:end])), k),
               stat(r -> (θ = r.theta[window(r)]; cor(θ[1:end-1], θ[2:end])), s))
    end
    VG = mean(stat(r -> mean(r.geno_var[window(r)]), k))
    @printf("  %-30s %.4g (Gaussian theory at the kernel's V_G)\n\n", "expected squared lag", expected_lag2(p; VG))
end

@printf("smallest p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy")
PMIN[] < 0.001 && exit(1)
