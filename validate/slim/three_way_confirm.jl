# Pre-specified confirmatory test of QTSim against the SLiM model generated from QTLSim's
# life cycle and against QTLSim's kernel, in the rare-mutation regime (N = 1000, one locus,
# m = 0.001, mutation rate 0.0002 per gamete, analytic ω, 500 generations), where earlier
# smaller batches had suggested QTSim's mean genetic variance ran ~2% high.
#
# Primary endpoint: mean V_G at generation 499; Welch t-tests QTSim vs SLiM and QTSim vs
# kernel, α = 0.025 each (Bonferroni). Secondary: KS, Anderson–Darling, upper quantiles.
# Random seeds per program (no shared seeds). Last result in `three_way_confirm.log`.
#
# Run from the package root; R runs per program (default 40000, ~2.5 h on 8 threads):
#   julia -t 8 --project=validate/slim validate/slim/three_way_confirm.jl [R]
using QTLSim, Printf, Random, Statistics, HypothesisTests
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
include(joinpath(@__DIR__, "..", "qtsim.jl"))
const SLIM = find_slim(); SLIM === nothing && error(slim_missing())
const QTSIM = find_qtsim(); QTSIM === nothing && error(qtsim_missing())

function qtsim_options(q::Params, seed)
    """
    Sample $(q.Sample)
    PopSize $(q.PopSize)
    NoQTL $(q.NoQTL)
    Generations $(q.NoGenerations)
    Ploidy $(q.Ploidy)
    PhenotypeVar $(q.PhenotypeVar)
    MutationScalar $(q.MutationScalar)
    Heritability $(q.Heritability)
    MutationRatePerGamete $(q.MutationRatePerGamete)
    Seed $seed
    """
end

"Genetic variance in QTSim's last output line."
function qtsim_end_vg(q::Params, seed)
    path = tempname() * ".txt"; write(path, qtsim_options(q, seed))
    text = try read(`$QTSIM $path`, String) finally rm(path; force = true) end
    occursin("Pop all died", text) && error("QTSim run died (seed $seed)")
    rows = [split(l, '\t') for l in split(text, '\n') if occursin(r"^\d+\t", l)]
    parse(Float64, rows[end][4])
end

"Genetic variance in the SLiM model's last output line."
function slim_end_vg(model::AbstractString, seed)
    path = tempname() * ".slim"; write(path, model)
    text = try read(`$SLIM -s $seed $path`, String) finally rm(path; force = true) end
    rows = [split(l, '\t') for l in split(text, '\n') if occursin(r"^\d+\t", l)]
    parse(Float64, rows[end][3])
end

kernel_end_vg(q::Params, seed) = simulate(q; seed).geno_var[end]

R = isempty(ARGS) ? 40_000 : parse(Int, ARGS[1])
base = Params(PopSize = 1000, NoQTL = 1, Ploidy = 2, PhenotypeVar = 0.02, Heritability = 0.5, Sample = 1,
              MutationScalar = 0.001, MutationRatePerGamete = 0.0002, InitFitnessSD = -1.0, NoGenerations = 500)
model = slim_script(lifecycle().generation, base)
programs = (:QTSim, :SLiM, :kernel)
seedlist = Dict(p => unique(rand(Xoshiro(9000 + i), 2^30:2^31, 2R))[1:R] for (i, p) in enumerate(programs))
jobs = [(p, k) for p in programs for k in 1:R]
out = Vector{Float64}(undef, length(jobs))
t = @elapsed Threads.@threads :greedy for j in eachindex(jobs)
    p, k = jobs[j]; seed = seedlist[p][k]
    out[j] = p === :QTSim ? qtsim_end_vg(base, seed) : p === :SLiM ? slim_end_vg(model, seed) : kernel_end_vg(base, seed)
end
d = Dict(p => [out[j] for j in eachindex(jobs) if jobs[j][1] === p] for p in programs)

@printf("%d runs per program (%.0f min); genetic variance at generation %d\n", R, t / 60, base.NoGenerations - 1)
for p in programs
    x = d[p]
    @printf("  %-7s mean %.7f ± %.7f   median %.7f   q90 %.6f   q99 %.6f\n", p, mean(x), std(x) / sqrt(length(x)),
            median(x), quantile(x, 0.9), quantile(x, 0.99))
end
for (a, b) in ((:SLiM, :QTSim), (:kernel, :QTSim), (:SLiM, :kernel))
    u, v = d[a], d[b]
    diff = mean(u) - mean(v); se = sqrt(var(u) / length(u) + var(v) / length(v))
    @printf("  %-6s − %-6s %+.2f%% ± %.2f%% (95%%)  Welch p = %.3g  KS p = %.3g  AD p = %.3g\n", a, b,
            100 * diff / mean(v), 196 * se / mean(v), pvalue(UnequalVarianceTTest(u, v)),
            pvalue(ApproximateTwoSampleKSTest(u, v)), pvalue(KSampleADTest(u, v)))
end
