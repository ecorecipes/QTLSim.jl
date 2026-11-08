# Fluctuating optima with several traits (`MoveType = :red_noise` and `:periodic` with
# `NoTraits = 2`), in three lowerings of the same life cycle: QTLSim's kernel against the
# theory of src/traits_theory.jl (`red_noise_lag_cov`, `periodic_lag_cov`, given the
# kernel's G) and against the SLiM model generated with the Eidos images of several traits.
#
# Two traits, 10 diploid loci, N = 500, V_E = 0.5 per trait, M = m·V_E·[1 0.4; 0.4 1]
# (m = 0.01), Ω = ω²·[1 −0.3; −0.3 1] with ω = 1.5. From generation 1000 the optimum
# fluctuates with no trend:
#   F1, red noise with stationary covariance Σ = [0.2 0.05; 0.05 0.1] and autocorrelation
#       ρ = 0.7 (`MoveNoiseCov`, `MoveAutocorr`);
#   F2, a cycle of amplitude 0.6 and period 25 generations along (1, −1)/√2
#       (`MoveAmplitude`, `MovePeriod`, `FluctDirection`).
# Per run, over generations 1500–2499: the mean of λλᵀ (λ = θ − z̄, the recorded lag; its
# entries L11, L12, L22) and the mean G matrix (G11, G12, G22).
#
# PRE-SPECIFIED CRITERIA (written before any run):
# - Theory: the kernel's mean of λλᵀ over runs against `red_noise_lag_cov(Σ, ρ; G, A)` (F1)
#   or `0.6²·periodic_lag_cov(25, f; G, A)` (F2), with G the kernel's mean G over runs and
#   A = `selection_variance_matrix(p)`. Each entry (i, j) must be within 10% of
#   sqrt(Λ_ii Λ_jj) of the prediction, failing only if the kernel's whole 95% CI is outside
#   that band (|kernel − pred| − 1.96 SE > 0.1·scale). The theory ignores drift's
#   contribution to the lag (about (G/N)/(2β) ≈ 2% of the red-noise lag here) and treats G
#   as constant.
# - SLiM: two-sample KS tests, kernel against SLiM, on L11, L12, L22, G11, G12, G22 of each
#   scenario; exits 1 if any p < 0.001 (as the other validate/slim scripts).
# Writes valid_traits_fluctuating.log.
#
# Run from the package root (R runs per backend and scenario, default 40; at most 3 threads here):
#   julia -t 3 --project=validate/slim validate/slim/compare_traits_fluctuating.jl [R]

using QTLSim, Statistics, HypothesisTests, Printf, Random
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const R = isempty(ARGS) ? 40 : parse(Int, ARGS[1])
const LOG = open(joinpath(@__DIR__, "valid_traits_fluctuating.log"), "w")
say(s) = (print(stdout, s); print(LOG, s); flush(LOG))
const LO, HI = 1500, 2499

base = Params(PopSize = 500, NoQTL = 10, Ploidy = 2, NoGenerations = 2500, PhenotypeVar = 1.0, Heritability = 0.5,
              MutationScalar = 0.01, MutationRatePerGamete = 0.02, RecombPoissonMean = 2.0, InitFitnessSD = 1.5,
              NoTraits = 2, MutCorr = 0.4, SelCorr = -0.3, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.0)
const Σ = [0.2 0.05; 0.05 0.1]
pF1 = Params(base; MoveType = :red_noise, MoveNoiseCov = Σ, MoveAutocorr = 0.7)
pF2 = Params(base; MoveType = :periodic, MoveAmplitude = 0.6, MovePeriod = 25.0, FluctDirection = [1.0, -1.0])

"A run's columns, from the kernel or from SLiM: gen, the trait means, G's upper triangle and the optimum."
function kernel_run(p, seed)
    o = simulate(p; seed)
    (; gen = o.gen, z = o.trait_mean, G = [[g[1, 1], g[1, 2], g[2, 2]] for g in o.gmatrix],
       theta = [o.trait_mean[i] + o.trait_lag[i] for i in eachindex(o.gen)], extinct = o.extinct)
end
function slim_run(path, seed)
    out = read(`$SLIM -s $seed $path`, String)
    rows = [split(l, '\t') for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    f(r, i) = parse(Float64, r[i])
    (; gen = [parse(Int, r[1]) for r in rows], z = [[f(r, 7), f(r, 8)] for r in rows],
       G = [[f(r, 9), f(r, 10), f(r, 11)] for r in rows], theta = [[f(r, 12), f(r, 13)] for r in rows],
       extinct = occursin("Pop all died.", out))
end

window(r) = findall(g -> LO <= g <= HI, r.gen)
# per run: the window means of λλᵀ's entries and of G's
function stats(r)
    w = window(r)
    λ = [r.theta[j] - r.z[j] for j in w]
    (L = [mean(x[1]^2 for x in λ), mean(x[1] * x[2] for x in λ), mean(x[2]^2 for x in λ)], G = mean(r.G[w]))
end

function both(p, seed0)
    path = tempname() * ".slim"; write(path, slim_script(lifecycle().generation, p))
    k = Vector{Any}(undef, R); s = Vector{Any}(undef, R)
    Threads.@threads :greedy for j in 2R:-1:1          # the slower SLiM runs first
        j <= R ? (k[j] = kernel_run(p, seed0 + j)) : (s[j - R] = slim_run(path, seed0 + j))
    end
    rm(path; force = true)
    any(r -> r.extinct, k) || any(r -> r.extinct, s) ? @warn("extinctions in a scenario meant to survive") : nothing
    stats.(k), stats.(s)
end

const PMIN = Ref(1.0)
const FAIL = String[]
function report(name, a, b)
    pv = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
    PMIN[] = min(PMIN[], pv)
    say(@sprintf("  %-10s kernel %10.4g ± %-9.2g SLiM %10.4g ± %-9.2g KS p = %.3f\n", name, mean(a), std(a), mean(b), std(b), pv))
end

function theory(name, k, Λ)
    say("  kernel against theory (with the kernel's mean G), tolerance 10% of sqrt(Λ_ii Λ_jj):\n")
    for (i, (lab, (a, b))) in enumerate(zip(("L11", "L12", "L22"), ((1, 1), (1, 2), (2, 2))))
        x = [s.L[i] for s in k]
        m, se = mean(x), std(x) / sqrt(length(x))
        scale = sqrt(Λ[a, a] * Λ[b, b])
        ok = abs(m - Λ[a, b]) - 1.96se <= 0.1 * scale
        ok || push!(FAIL, "$name $lab against theory")
        say(@sprintf("    %-4s theory %9.4g  kernel %9.4g ± %-8.2g (dev %+6.1f%% of scale)  %s\n", lab, Λ[a, b], m, 1.96se,
                     100 * (m - Λ[a, b]) / scale, ok ? "ok" : "FAIL"))
    end
end

say("QTLSim kernel vs theory vs generated SLiM model: fluctuating optima with two traits, R = $R runs each, $(Threads.nthreads()) threads\n")
say(@sprintf("N = %d, %d loci, V_E = %.2g, M = %s, Ω = %s; windows %d-%d\n\n", base.PopSize, base.NoQTL, envVar(base),
             string(round.(mutational_covariance(base); digits = 5)), string(round.(selection_covariance(base); digits = 4)), LO, HI))
A = selection_variance_matrix(base)
for (name, p, seed0, pred) in (("F1 red noise, Σ = [0.2 0.05; 0.05 0.1], ρ = 0.7", pF1, 30_000,
                                G -> red_noise_lag_cov(Σ, 0.7; G, A)),
                               ("F2 periodic, amplitude 0.6, period 25, along (1, -1)/√2", pF2, 40_000,
                                G -> 0.6^2 .* periodic_lag_cov(25.0, fluct_direction(pF2); G, A)))
    say("── $name ──\n")
    t0 = time(); k, s = both(p, seed0); t = time() - t0
    for (i, nm) in enumerate(("L11", "L12", "L22"))
        report(nm, [x.L[i] for x in k], [x.L[i] for x in s])
    end
    for (i, nm) in enumerate(("G11", "G12", "G22"))
        report(nm, [x.G[i] for x in k], [x.G[i] for x in s])
    end
    Gk = mean(x.G for x in k)
    Gm = [Gk[1] Gk[2]; Gk[2] Gk[3]]
    theory(name[1:2], k, pred(Gm))
    say(@sprintf("  (%.0f s; the kernel's mean G = [%.4g %.4g; %.4g %.4g])\n\n", t, Gm[1, 1], Gm[1, 2], Gm[2, 1], Gm[2, 2]))
end

say(@sprintf("smallest KS p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY between kernel and SLiM" : "no discrepancy between kernel and SLiM"))
say(isempty(FAIL) ? "theory: every pre-specified criterion holds\n" : "theory: FAIL " * join(FAIL, "; ") * "\n")
close(LOG)
(PMIN[] < 0.001 || !isempty(FAIL)) && exit(1)
