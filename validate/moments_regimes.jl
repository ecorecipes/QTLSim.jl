# The moments lowering across regimes of mutation, selection and drift: the multilocus
# Gaussian model (`loci = :multilocus`), the rare-alleles closure (`loci = :rare_alleles`,
# src/rare_alleles.jl; since this script's runs, the default for one diploid trait) and the closed forms (`house_of_cards_vg`, Bürger's
# stochastic house of cards; `house_of_gauss_vg`, Hermisson and Wagner's 2004 stochastic
# house of Gauss) against the kernel: equilibrium V_G over a grid, the steady lag behind
# a continuous move, the time to extinction under a fast move, and the pre-move V_G of
# the paper's published populations.
#
# Usage (from the package root; it needs Catlab and XLSX, which the vignettes environment has):
#     julia -t 3 --project=vignettes validate/moments_regimes.jl [R]
# R replicates per grid cell (default 16). Writes validate/moments_regimes.log
# (validate/moments_regimes_seed<SEED_OFFSET>.log with ENV SEED_OFFSET ≠ 0). SMOKE=1 runs
# short versions of every path and writes no log. The paper-regime part needs
# Repository/Simulations and Repository/Results (read only); it is skipped without them.
#
# GRID. Diploid, V_E = 0.5 (V_P = 1, h² = 0.5), mutation scalar m = 0.01 (V_m = 0.005),
# :qtsim crossover with RecombPoissonMean = 2, fixed size N. Cells: N ∈ {100, 1000} ×
# 4Nμ ∈ {0.01, 0.1, 1, 4, 10} (μ per locus per gamete) × L ∈ {10, 100, 1000} × two
# widths of selection, V_s = ω² + V_E ∈ {50.5 (the paper's analytical ω = 7.07), 5.5
# (ω = √5)}: 60 cells. The effect of one mutation is α² = V_m/(2Lμ), and its selection
# coefficient as a heterozygote x = α²/(2V_s) (from 2e-8 to 4.5 over the grid). Genomes
# start at zero; each run is 30,000 generations and its V_G (`geno_var`) is averaged over
# generations 10,000:10:29,990; SEs are over replicates.
#
# PREDICTIONS (no constant fitted to the kernel): `gaussian_vg_loci(p; drift = true)`
# (multilocus Gaussian with drift), `rare_alleles_vg(p; drift = true)` (the rare-alleles
# closure with drift, the fixed point of `run_moments(p; loci = :rare_alleles, erosion =
# true)`), `house_of_cards_vg(p; nsurv)` and `house_of_gauss_vg(p; nsurv)` with the
# kernel's mean number of survivors.
#
# PRE-SPECIFIED CRITERIA (written before any kernel run of this script). As in
# validate/moments_loci.jl, a comparison FAILS only if the kernel's whole 95% CI lies
# outside the band: |kernel − pred| − 1.96 SE > tol·pred. "Weak" cells are those with
# x = α²/(2V_s) ≤ 0.05 (the closures assume weak selection per mutation).
# - C0: in every cell, run_moments(p; loci = :rare_alleles, erosion = true) converges to
#   rare_alleles_vg(p; drift = true) (relative difference < 1e-6).
# - R1: in every weak cell, the kernel's V_G is within 15% of rare_alleles_vg(p; drift = true).
# - R2: in every weak cell, the rare-alleles prediction is not farther from the kernel than
#   the multilocus Gaussian: |log(k/rare)| ≤ |log(k/gauss)| + 0.02.
# - R3 (the Gaussian regime: weak cells with N = 1000 and 4Nμ ≥ 4): within 10% of rare_alleles_vg.
# - R4 (the closed form, linkage equilibrium: weak cells with L ≤ 100 and 4Nμ ≤ 1): the kernel
#   within 20% of house_of_gauss_vg(p; nsurv).
# - S (steady lag; L = 100, N = 1000, μ = 1e-4 (4Nμ = 0.4), V_s = 50.5, a continuous move
#   of 0.01 per generation from generation 10,000; lag θ − z̄ and V_G averaged over
#   15,000:10:19,990; 3R replicates): the kernel's lag within 15% of steady_lag(p; VG =
#   rare_alleles_vg(p; drift = true)), within 15% of the moments lowering's (loci =
#   :rare_alleles, erosion = true, last generation), and within 10% of steady_lag with the
#   kernel's own V_G in the window.
# - E (extinction; the same population with a move of k per generation from 10,000, k
#   chosen so that the rare-alleles moments lowering predicts the deterministic lag to
#   approach 1.15 times the lag at which N·W̄ = 2; 3R replicates, 2,000 generations after the
#   move): the median generation of extinction after the move (survivors counted as
#   later than any death) within 20% of the median of `extinction_distribution(p, m)` with
#   m the rare-alleles moments trajectory (erosion = true).
# - P (the paper's populations; needs Repository): 40 published inputs drawn at random
#   from each of the 8 strata (NoQTL ∈ {100, 500, 1000, 5000} × u ∈ {3, 5}), matched to
#   their rows of the published table; the observed V_G is the table's "Genotype Var at
#   10000" (QTSim's generation 9999; one realisation per input). Prediction: the moments
#   lowering run from zero genomes for 10,000 generations (erosion = true), recorded at
#   generation 9,999, with loci = :rare_alleles ("rare") and :multilocus ("gauss").
#   P1: the median over inputs of table/rare within [0.85, 1.15]. P2: the median of
#   |log(table/pred)| smaller for rare than for gauss. P3: the kernel, 4 replicates on each
#   of 24 of these inputs (4 per stratum with NoQTL ≤ 1000): Σ kernel / Σ rare over the 24
#   inputs (replicate means) within [0.85, 1.15].
# Pre-specified expectations, reported but not criteria: in the strong cells (x > 0.05) the
# weak-selection closures misjudge V_G; the multilocus Gaussian fails (kernel < 85% of it)
# in the cells with 4Nμ ≤ 0.1 and N = 1000; under the move (S), the kernel's V_G exceeds its
# stationary value (directional selection on rare alleles, Barton and Turelli 1987), which
# no model here captures.
# Disclosure: the model was built knowing validate/moments_loci.log (and its fresh-seed
# confirmation), whose μ = 1e-4 and μ = 0.01 cells (L ≤ 500, V_s = 50.5) are near some of
# the cells here: the first version of the closure (without the LD factor β², see
# src/rare_alleles.jl) predicted 0.81–1.06 of those kernels' V_G in the house-of-cards
# cells; with β² it predicts 0.95–1.17 (1.17 at L = 1). No kernel run of this script had
# been made when these criteria were written. The script exits 1 if any criterion fails.
#
# CONFIRMATION (2026-11-07): an independent run with fresh seeds (SEED_OFFSET = 7000000,
# moments_regimes_seed7000000.log) gave the same result, 122 PASS and the same 4 FAIL (R1 at
# N = 1000, 4Nμ = 0.1, L = 1000, V_s = 5.5, ratio 0.716; both R4 at V_s = 5.5, L = 100; S),
# so those failures are systematic limits of the model, not noise.
# HISTORY. Run 1 (2026-11-07, R = 16, SEED_OFFSET = 0, moments_regimes.log, about 45 min
# on 3 threads) failed 4 criteria, kept as failures: R1 in one of 47 weak cells (N = 1000,
# 4Nμ = 0.1, L = 1000, V_s = 5.5: kernel/rare 0.72, strong selection with strong LD, between
# the regimes); R4 in two cells with V_s = 5.5 and L = 100 (the closed form ignores LD:
# 1.42 and 0.75); S's third check (the kernel's lag is 1.21 times steady_lag with its own
# V_G under the move, which rose from about 1.8 to 2.5, the Barton–Turelli effect, so the
# Gaussian response V_G λ/s does not hold there). Everything else passed: C0, R1 in 46
# of 47 weak cells (0.85–1.04, most within 5%), R2 in all 47, R3 in all 12, R4 in 10 of 12,
# S's lag within 12% of the rare-alleles prediction, E (kernel median 45 generations to
# extinction against 41 predicted; the Gaussian predicts 210), P1–P3 (paper: median
# table/rare 0.991, gauss 0.966). The model was not changed after the run.

using QTLSim, Statistics, Printf, Random
import Catlab, XLSX

const SMOKE = get(ENV, "SMOKE", "") == "1"
R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : (SMOKE ? 2 : 16)
const SEED_OFFSET = parse(Int, get(ENV, "SEED_OFFSET", "0"))
logpath = SMOKE ? tempname() : joinpath(@__DIR__, SEED_OFFSET == 0 ? "moments_regimes.log" : "moments_regimes_seed$(SEED_OFFSET).log")
io = open(logpath, "w")
say(s...) = (println(stdout, s...); println(io, s...); flush(io))
se(x) = std(x) / sqrt(length(x))
fails = String[]
check(ok, what) = (say(ok ? "  PASS  " : "  FAIL  ", what); ok || push!(fails, what))
outside(k, s, pred, tol) = abs(k - pred) - 1.96s > tol * abs(pred)
lr(a, b) = abs(log(a / b))

const NGEN, WINDOW = SMOKE ? (400, 200:10:390) : (30_000, 10_000:10:29_990)
const MOVE, LW = SMOKE ? (200, 300:10:390) : (10_000, 15_000:10:19_990)
const MGEN = SMOKE ? 20_000 : 200_000        # generations of the moments lowering for C0
base = Params(Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01, RecombPoissonMean = 2.0,
              NoGenerations = NGEN, Sample = 10)

function kernel_run(p, seed; window = WINDOW, lagwindow = nothing)
    o = simulate(p; seed)
    if o.extinct
        return (vg = NaN, nsurv = NaN, lag = NaN, extinct_gen = o.extinct_gen)
    end
    ix = findall(in(window), o.gen)
    (vg = mean(o.geno_var[ix]), nsurv = mean(o.nsurv[ix]),
     lag = lagwindow === nothing ? NaN : mean((o.theta .- o.pheno_mean)[findall(in(lagwindow), o.gen)]),
     extinct_gen = -1)
end

function replicates(p, R, seed0; kw...)
    out = Vector{Any}(undef, R)
    Threads.@threads :greedy for r in 1:R
        out[r] = kernel_run(p, SEED_OFFSET + seed0 + r; kw...)
    end
    out
end

say("QTLSim validate/moments_regimes.jl: the moments lowering across regimes (Gaussian, house of cards, drift)")
say("Julia threads: ", Threads.nthreads(), "; replicates per cell R = $R; SEED_OFFSET = $SEED_OFFSET")
say(@sprintf("V_E = %.2f, V_m = %.4f; %d generations, window %d:%d:%d", envVar(base), muteVar(base), NGEN,
             first(WINDOW), step(WINDOW), last(WINDOW)))
say()

# --- the grid ------------------------------------------------------------------------
const VSS = (50.5, 5.5)
cells = [(N, th, L, Vs) for Vs in VSS for N in (1000, 100) for th in (0.01, 0.1, 1.0, 4.0, 10.0) for L in (10, 100, 1000)]
cellparams(N, th, L, Vs) = Params(base; NoQTL = L, PopSize = N, MutationRatePerGamete = th / (4N),
                                  InitFitnessSD = Vs == 50.5 ? -1.0 : sqrt(Vs - envVar(base)))
rows = []
t0 = time()
for (i, (N, th, L, Vs)) in enumerate(cells)
    p = cellparams(N, th, L, Vs)
    @assert abs(selection_variance(p) - Vs) < 1e-6
    t = @elapsed rs = replicates(p, R, 10_000 * i)
    any(r -> r.extinct_gen >= 0, rs) && error("a replicate went extinct in cell $i")
    vk = [r.vg for r in rs]
    ns = mean(r.nsurv for r in rs)
    α2 = QTLSim.mutation_effect_variance(p)
    gauss, rare = gaussian_vg_loci(p; drift = true), rare_alleles_vg(p; drift = true)
    mo = run_moments(Params(p; NoGenerations = MGEN, Sample = 1000); loci = :rare_alleles, erosion = true).geno_var[end]
    push!(rows, (; N, th, L, Vs, μ = p.MutationRatePerGamete, α2, x = α2 / (2Vs), vk = mean(vk), sk = se(vk), gauss, rare, mo,
                 shc = house_of_cards_vg(p; nsurv = ns), hog = house_of_gauss_vg(p; nsurv = ns), t))
    @printf("cell %2d/%d  N = %4d  4Nμ = %-5g L = %4d  V_s = %4.1f  %6.1f s  (%.0f s so far)\n", i, length(cells), N, th, L, Vs, t, time() - t0)
end
say(@sprintf("%d cells in %.0f s", length(cells), time() - t0))
say()
say("Equilibrium V_G: kernel (± SE) and kernel/prediction; x = α²/(2V_s), weak cells x ≤ 0.05")
say("   N   4Nμ     L   V_s      α²        x       kernel      ±SE   gauss  ratio    rare  ratio   SHC ratio   HoG ratio")
for r in rows
    say(@sprintf("%4d %5g %5d %5.1f %8.2e %8.2e %10.4g %8.2g %7.4g %6.3f %7.4g %6.3f %6.3f %6.3f",
                 r.N, r.th, r.L, r.Vs, r.α2, r.x, r.vk, r.sk, r.gauss, r.vk / r.gauss, r.rare, r.vk / r.rare,
                 r.vk / r.shc, r.vk / r.hog))
end
say()

say("Criteria, grid")
check(all(abs(r.mo / r.rare - 1) < 1e-6 for r in rows),
      "C0: run_moments(p; loci = :rare_alleles, erosion = true) converges to rare_alleles_vg(p; drift = true) in every cell")
weak = [r for r in rows if r.x <= 0.05]
say(@sprintf("  %d weak cells of %d", length(weak), length(rows)))
cell(r) = @sprintf("N = %d, 4Nμ = %g, L = %d, V_s = %.1f", r.N, r.th, r.L, r.Vs)
for r in weak
    check(!outside(r.vk, r.sk, r.rare, 0.15), @sprintf("R1: %s: kernel %.4g within 15%% of rare %.4g (ratio %.3f)", cell(r), r.vk, r.rare, r.vk / r.rare))
end
for r in weak
    check(lr(r.vk, r.rare) <= lr(r.vk, r.gauss) + 0.02,
          @sprintf("R2: %s: rare not farther than gauss (ratios %.3f, %.3f)", cell(r), r.vk / r.rare, r.vk / r.gauss))
end
for r in [r for r in weak if r.N == 1000 && r.th >= 4]
    check(!outside(r.vk, r.sk, r.rare, 0.10), @sprintf("R3: %s: kernel within 10%% of rare (ratio %.3f)", cell(r), r.vk / r.rare))
end
for r in [r for r in weak if r.L <= 100 && r.th <= 1]
    check(!outside(r.vk, r.sk, r.hog, 0.20), @sprintf("R4: %s: kernel within 20%% of house_of_gauss_vg (ratio %.3f)", cell(r), r.vk / r.hog))
end
strong = [r for r in rows if r.x > 0.05]
say(@sprintf("  expectation (not a criterion): strong cells (x > 0.05), kernel/rare: %s",
             join([@sprintf("%.2f", r.vk / r.rare) for r in strong], ", ")))
gf = [r for r in rows if r.N == 1000 && r.th <= 0.1]
ng = count(r -> r.vk / r.gauss < 0.85, gf)
say(@sprintf("  expectation (not a criterion): multilocus Gaussian fails (kernel < 85%%) in %d of %d cells with N = 1000, 4Nμ ≤ 0.1: %s",
             ng, length(gf), ng == length(gf) ? "confirmed" : "refuted"))
say()

# --- S: steady lag in the house-of-cards regime ----------------------------------------
pH = Params(base; NoQTL = 100, PopSize = 1000, MutationRatePerGamete = 1e-4, NoGenerations = last(LW) + 10)
pS = Params(pH; MoveFitnessSDGen = MOVE, MoveType = :continuous, MoveFitnessSD = 0.01)
t = @elapsed rS = replicates(pS, 3R, 900_000; window = LW, lagwindow = LW)
any(r -> r.extinct_gen >= 0, rS) && error("a replicate of S went extinct")
lk, ls = mean(r.lag for r in rS), se([r.lag for r in rS])
vkS = mean(r.vg for r in rS)
Vr, Vg = rare_alleles_vg(pS; drift = true), gaussian_vg_loci(pS; drift = true)
λr, λg = steady_lag(pS; VG = Vr), steady_lag(pS; VG = Vg)
λmo = (o = run_moments(pS; loci = :rare_alleles, erosion = true); o.theta[end] - o.pheno_mean[end])
λk = steady_lag(pS; VG = vkS)
say(@sprintf("S, steady lag behind a continuous move of 0.01 per generation, house of cards (L = 100, 4Nμ = 0.4) (%.0f s)", t))
say(@sprintf("  kernel lag %.4f ± %.4f; kernel V_G under the move %.4f (stationary predictions: rare %.4f, gauss %.4f)", lk, ls, vkS, Vr, Vg))
say(@sprintf("  steady_lag with rare %.4f (ratio %.3f), moments lowering (rare) %.4f (ratio %.3f), with gauss %.4f (ratio %.3f), with the kernel's V_G %.4f (ratio %.3f)",
             λr, lk / λr, λmo, lk / λmo, λg, lk / λg, λk, lk / λk))
check(!outside(lk, ls, λr, 0.15), "S: kernel lag within 15% of steady_lag with rare_alleles_vg(p; drift = true)")
check(!outside(lk, ls, λmo, 0.15), "S: kernel lag within 15% of the rare-alleles moments lowering's")
check(!outside(lk, ls, λk, 0.10), "S: kernel lag within 10% of steady_lag with the kernel's V_G")
say(@sprintf("  expectation (not a criterion): kernel V_G under the move above the stationary rare prediction: %s",
             vkS > Vr ? "confirmed" : "refuted"))
say()

# --- E: extinction under a fast move ---------------------------------------------------
# the lag at which N W̄ = 2: W̄ = c ω/√s exp(−λ²/2s), s = ω² + V_E + V_G
c = 1 / (0.4 * sqrt(2π))
sE = initFitSD(pH)^2 + envVar(pH) + Vr
λdie = sqrt(2sE * log(pH.PopSize * c * initFitSD(pH) / sqrt(sE) / 2))
kE = 1.15 * λdie * Vr / sE
AFTER = SMOKE ? 200 : 2000
pE = Params(pH; NoGenerations = MOVE + AFTER, Sample = 1, MoveFitnessSDGen = MOVE, MoveType = :continuous, MoveFitnessSD = kE)
t = @elapsed rE = replicates(pE, 3R, 950_000; window = 0:0)
tk = sort([r.extinct_gen >= 0 ? r.extinct_gen - MOVE : typemax(Int) for r in rE])
medk = tk[cld(length(tk), 2)]
function predicted_median(loci)
    m = run_moments(pE; loci, erosion = true)
    d = extinction_distribution(pE, m)
    d.median - MOVE
end
medr, medg = predicted_median(:rare_alleles), predicted_median(:multilocus)
say(@sprintf("E, extinction under a move of %.4f per generation (1.15 × the dying lag %.2f × V_G/s) (%.0f s)", kE, λdie, t))
say(@sprintf("  kernel: %d of %d died within %d generations; median generation of extinction after the move %s",
             count(<(typemax(Int)), tk), length(tk), AFTER, medk == typemax(Int) ? "none (survived)" : string(medk)))
say(@sprintf("  predicted median (extinction_distribution of the moments trajectory): rare %s, gauss %s",
             medr, medg))
check(medk != typemax(Int) && medr isa Number && abs(medk - medr) <= 0.20 * medr,
      "E: the kernel's median time to extinction within 20% of the rare-alleles moments prediction")
say()

# --- P: the paper's populations --------------------------------------------------------
const ROOTDIR = normpath(joinpath(@__DIR__, "..", ".."))
const SIMDIR = joinpath(ROOTDIR, "Repository", "Simulations")
const MASTER = joinpath(ROOTDIR, "Repository", "Results", "Master Output file.xlsx")
if isdir(SIMDIR) && isfile(MASTER)
    dt = XLSX.readtable(MASTER, "Master output file")
    lab = string.(dt.column_labels)
    col(name) = dt.data[findfirst(==(name), lab)]
    num(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)
    tab = (L = round.(Int, num.(col("No. Loci"))), u = round.(Int, num.(col("recombination u=3,5"))),
           H = num.(col("Heritability")), PV = num.(col("PhenotypeVar")), MS = num.(col("MutationScalar Vm/Ve")),
           d = num.(col("d_env")), gv = num.(col("Genotype Var at 10000")))
    function match_row(L, u, q)
        hits = findall(i -> tab.L[i] == L && tab.u[i] == u && abs(tab.H[i] - q.Heritability) < 1e-9 &&
                            abs(tab.PV[i] - q.PhenotypeVar) < 1e-9 && abs(tab.d[i] - q.MoveFitnessSD) <= 5e-3 * q.MoveFitnessSD &&
                            abs(tab.MS[i] - q.MutationScalar) <= 5e-3 * q.MutationScalar, eachindex(tab.L))
        length(hits) == 1 ? hits[1] : 0
    end
    strata = [("QT0100-Possion3.0", 100, 3), ("QT0100-Possion5.0", 100, 5), ("QT0500-Possion3.0", 500, 3),
              ("QT0500-Possion5.0", 500, 5), ("QT1000-Possion3.0", 1000, 3), ("QT1000-Possion5.0", 1000, 5),
              ("QT5000-Possion3.0", 5000, 3), ("QT5000-Possion5.0", 5000, 5)]
    NPER = SMOKE ? 2 : 40
    rng = MersenneTwister(20261107 + SEED_OFFSET)
    prow = []
    t0 = time()
    for (dir, L, u) in strata
        files = sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, dir))))
        for f in sort(shuffle(rng, files)[1:NPER])
            q = read_qtsim_options(joinpath(SIMDIR, dir, f))
            q = Params(q; Seed = nothing, NoGenerations = 10_000, Sample = 1, MoveFitnessSDGen = -1)
            r = match_row(L, u, q)
            r == 0 && (say("  no unique table row for $dir/$f; skipped"); continue)
            pred(loci) = (o = run_moments(Params(q; Sample = 9_999); loci, erosion = true); o.geno_var[findfirst(==(9_999), o.gen)])
            α2 = QTLSim.mutation_effect_variance(q)
            push!(prow, (; dir, file = f, L, u, q, obs = tab.gv[r], rare = pred(:rare_alleles), gauss = pred(:multilocus),
                         rare_eq = rare_alleles_vg(q; drift = true), shc = house_of_cards_vg(q), hog = house_of_gauss_vg(q),
                         th = 4 * q.PopSize * q.MutationRatePerGamete, x = α2 / (2selection_variance(q)),
                         se = q.PopSize * α2 / selection_variance(q)))
        end
        @printf("  %s: %d inputs (%.0f s so far)\n", dir, count(r -> r.dir == dir, prow), time() - t0)
    end
    say(@sprintf("P, the paper's populations: %d published inputs matched to the table (%.0f s)", length(prow), time() - t0))
    say("  per stratum: median table/pred (moments, 10,000 generations from zero) and median |log ratio|; 4Nμ, x = α²/(2V_s), N α²/V_s medians")
    say("  stratum               n   rare  |log|   gauss |log|   rare_eq   SHC    HoG     4Nμ      x        Nα²/V_s")
    med(v) = median(v)
    for (dir, L, u) in vcat(strata, [("all", 0, 0)])
        rs = dir == "all" ? prow : [r for r in prow if r.dir == dir]
        isempty(rs) && continue
        say(@sprintf("  %-18s %3d  %6.3f %6.3f %6.3f %6.3f %7.3f %7.3f %6.3f %7.3g %9.2e %9.3g", dir, length(rs),
                     med([r.obs / r.rare for r in rs]), med([lr(r.obs, r.rare) for r in rs]),
                     med([r.obs / r.gauss for r in rs]), med([lr(r.obs, r.gauss) for r in rs]),
                     med([r.obs / r.rare_eq for r in rs]), med([r.obs / r.shc for r in rs]), med([r.obs / r.hog for r in rs]),
                     med([r.th for r in rs]), med([r.x for r in rs]), med([r.se for r in rs])))
    end
    mr, mg = med([lr(r.obs, r.rare) for r in prow]), med([lr(r.obs, r.gauss) for r in prow])
    check(0.85 <= med([r.obs / r.rare for r in prow]) <= 1.15,
          @sprintf("P1: median table/rare %.3f within [0.85, 1.15]", med([r.obs / r.rare for r in prow])))
    check(mr < mg, @sprintf("P2: median |log(table/pred)| smaller for rare (%.3f) than for gauss (%.3f)", mr, mg))
    # P3: the kernel on 24 of these inputs
    kin = vcat([[r for r in prow if r.dir == dir][1:min(SMOKE ? 1 : 4, count(r -> r.dir == dir, prow))] for (dir, L, u) in strata if L <= 1000]...)
    t = @elapsed kres = map(enumerate(kin)) do (j, r)
        rs = replicates(Params(r.q; Sample = 9_999), SMOKE ? 2 : 4, 2_000_000 + 100j; window = 9_999:9_999)
        mean(x.vg for x in rs)
    end
    say(@sprintf("  P3: kernel on %d inputs (%.0f s): Σ kernel/Σ rare %.3f, Σ kernel/Σ gauss %.3f, Σ table/Σ rare (same inputs) %.3f",
                 length(kin), t, sum(kres) / sum(r.rare for r in kin), sum(kres) / sum(r.gauss for r in kin),
                 sum(r.obs for r in kin) / sum(r.rare for r in kin)))
    for (r, k) in zip(kin, kres)
        say(@sprintf("    %s/%s  kernel %.4g  table %.4g  rare %.4g  gauss %.4g", r.dir, r.file, k, r.obs, r.rare, r.gauss))
    end
    check(0.85 <= sum(kres) / sum(r.rare for r in kin) <= 1.15, "P3: Σ kernel / Σ rare within [0.85, 1.15]")
else
    say("P skipped: Repository/Simulations or the published table not found")
end
say()

say(isempty(fails) ? "PASS: every pre-specified criterion holds" : "FAILED ($(length(fails))): " * join(fails, "; "))
close(io)
exit(isempty(fails) ? 0 : 1)
