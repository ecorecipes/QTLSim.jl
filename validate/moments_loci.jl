# The multilocus Gaussian model of the moments lowering (`run_moments`, `loci =
# :multilocus`, the default; src/multilocus.jl) against the kernel and against the old
# aggregate recursion (`loci = :aggregate`, `gaussian_vg_exact`), at equilibrium under
# stabilising selection, over the number of loci, the population size, the crossover
# model and the mutation regime; then the steady lag behind a continuous move and the G
# matrix of two traits.
#
# Usage (from the package root; it needs Catlab, which the vignettes environment has):
#     julia -t 3 --project=vignettes validate/moments_loci.jl [R]
# R replicates per cell (default 12). Writes validate/moments_loci.log.
#
# The model: diploid, V_P = 1 at the start, h² = 0.5 (V_E = 0.5), mutation scalar m = 0.01
# (V_m = m·V_E = 0.005 per individual per generation), the paper's analytical fitness
# width ω = h² V_P/√(m V_E) = 7.07 (V_s = ω² + V_E = 50.5), QTSim's fixed size N.
# Cells: L ∈ {1, 10, 50, 100, 500} loci × N ∈ {100, 1000} × CrossoverModel ∈ {:qtsim,
# :all} (one cell for L = 1, where there is no crossover) with RecombPoissonMean = 2, ×
# two mutation regimes at the same V_m: μ = 0.01 mutations per locus per gamete (many
# small mutations, α² = V_m/(2Lμ)) and μ = 1e-4 (rare large ones). Genomes start at zero;
# each run is 20,000 generations and its V_G (the kernel's `geno_var`, the sample variance
# of the genotypes) is averaged over generations 10,000:10:19,990; the kernel's genic
# variance and LD are measured once, on the zygotes at the end. Standard errors are over
# replicates.
#
# Predictions: `gaussian_vg_loci(p; drift = true)` (the multilocus model with drift, the
# fixed point of `run_moments(p; erosion = true, loci = :multilocus)`, checked here; the
# script passes `loci = :multilocus` explicitly since the default became the rare-alleles
# closure, so its log is unchanged), `gaussian_vg_loci(p)`
# (without drift), `multilocus_equilibrium(p; drift = true, loci = :multilocus).ld` (LD), the old aggregate
# fixed point `gaussian_vg_exact(p)` and the aggregate with drift
# (`run_moments(p; loci = :aggregate, erosion = true)`), and `house_of_cards_vg(p;
# nsurv)` with the kernel's mean number of survivors.
#
# Regimes, from the Gaussian per-locus assumption: a locus is Gaussian when many alleles
# segregate at it, 4Nμ ≫ 1, and its variance is many mutations' worth, V/(2L) ≫ α². Here
# μ = 0.01 gives 4Nμ = 40 at N = 1000 (Gaussian) and 4 at N = 100 (intermediate); μ = 1e-4
# gives 4Nμ = 0.4 and 0.04 (house of cards).
#
# PRE-SPECIFIED CRITERIA (written before any kernel run of this script). A comparison
# FAILS only if the kernel's whole 95% CI lies outside the tolerance band, i.e.
# |kernel − pred| − 1.96 SE > tol·|pred|.
# - G1 (Gaussian regime: μ = 0.01, N = 1000, L ∈ {10, 50, 100, 500}, both crossover
#   models; 8 cells): V_G within 10% of `gaussian_vg_loci(p; drift = true)`.
# - G2 (μ = 0.01, both N, L ≥ 10; 16 cells): the multilocus prediction with drift is
#   closer to the kernel than the old aggregate prediction `gaussian_vg_exact`, and than
#   the aggregate with drift, in |log(kernel/pred)|.
# - G3 (the same 8 cells as G1): the kernel's LD (V_G minus the genic variance, at the end)
#   is negative and within 25% of `multilocus_equilibrium(p; drift = true).ld`.
# - S (steady lag; L = 50, N = 1000, μ = 0.01, :qtsim, a continuous move of 0.01 per
#   generation from generation 10,000, lag θ − z̄ averaged over 15,000:10:19,990): the
#   kernel's lag within 10% of `steady_lag(p; VG = gaussian_vg_loci(p; drift = true))`,
#   within 10% of the moments lowering's (`run_moments(p; erosion = true)`, last
#   generation), and within 10% of `steady_lag` with the kernel's own V_G.
# - T (two traits; L = 50, N = 1000, μ = 0.01, :qtsim, mutational correlation 0.5,
#   selection correlation −0.3; window as above): every entry of G within
#   0.10·√(G_ii G_jj) of `gaussian_G_loci(p; drift = true)`, the genetic correlation
#   within 0.05, g_max within 10°.
# - The moments lowering with erosion = true converges to `gaussian_vg_loci(p; drift =
#   true)` (relative difference < 1e-6) in every cell.
# Pre-specified expectation, reported but not a criterion: in the house-of-cards cells
# (μ = 1e-4, L ≥ 10) the Gaussian assumption fails and the kernel holds less than 85% of
# the multilocus prediction.
# Disclosure: two numbers in this regime were known before these criteria were written:
# validate/traits_kernel_vs_moments.log's informational case (one and two traits, L = 50,
# N = 1000, μ = 0.01, the parameters of cell G1 at L = 50 and of T), where the kernel's V_G
# was 3.56 and its G [3.00 0.75; 0.75 2.98], against which the new model with drift,
# computed while writing it, gives 3.82 (7% above) and [3.43 0.85; 0.85 3.43] (14%). The
# 10% tolerance is that script's, kept unchanged; T is expected to fail on that evidence.
# The script exits 1 if any criterion fails.
#
# HISTORY (for honesty about what changed between runs): the first run of this script
# (validate/moments_loci_run1.log) used one LD class recombining at the mean pair rate
# `pair_recombination(p)`. It failed G1 at L = 100 (:qtsim) and L = 500 (both models), G3
# (the kernel's LD 1.3–2.5× the predicted, most under :qtsim) and T (G11 0.885, G22 0.881).
# The LD diagnostic pointed at the mean rate: pairs of adjacent loci are rarely split, and
# the equilibrium LD of a pair is proportional to 1/ρ, so the mean ρ underestimates it
# (Jensen). The model now resolves LD by class of pairs (`QTLSim.pair_classes`, each with
# its own exact ρ from the crossover model; no fitted constant). The criteria, seeds and
# tolerances are unchanged; this log is the second run. An independent confirmation with
# fresh seeds (`SEED_OFFSET=7000000`, `moments_loci_seed7000000.log`, 2026-11-07) passed every
# criterion too (40 of 40).

using QTLSim, Statistics, Printf, Random
import Catlab

const SMOKE = get(ENV, "SMOKE", "") == "1"     # a check of the code paths: short runs, no log
R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
const SEED_OFFSET = parse(Int, get(ENV, "SEED_OFFSET", "0"))   # nonzero: an independent confirmation run with fresh seeds
logpath = SMOKE ? tempname() : joinpath(@__DIR__, SEED_OFFSET == 0 ? "moments_loci.log" : "moments_loci_seed$(SEED_OFFSET).log")
io = open(logpath, "w")
say(s...) = (println(stdout, s...); println(io, s...); flush(io))
se(x) = std(x) / sqrt(length(x))
fails = String[]
check(ok, what) = (say(ok ? "  PASS  " : "  FAIL  ", what); ok || push!(fails, what))
# the kernel's 95% CI lies outside the band pred·(1 ± tol)
outside(k, s, pred, tol) = abs(k - pred) - 1.96s > tol * abs(pred)

const NGEN, WINDOW = SMOKE ? (300, 100:10:290) : (20_000, 10_000:10:19_990)
const MOVE, LW = SMOKE ? (100, 200:10:290) : (10_000, 15_000:10:19_990)
base = Params(Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01, RecombPoissonMean = 2.0,
              NoGenerations = NGEN, Sample = 10)

# one replicate: window means and, at the end, the zygotes' V_G and genic variance
function kernel_run(p, seed; window = WINDOW, lagwindow = nothing)
    w = World(p)
    o = run!(w, MersenneTwister(seed))
    o.extinct && error("a replicate went extinct: $(p.NoQTL) loci, N = $(p.PopSize)")
    ix = findall(in(window), o.gen)
    out = (vg = mean(o.geno_var[ix]), nsurv = mean(o.nsurv[ix]),
           lag = lagwindow === nothing ? NaN : mean((o.theta .- o.pheno_mean)[findall(in(lagwindow), o.gen)]))
    if p.NoTraits > 1
        return merge(out, (G = mean(o.gmatrix[ix]),))
    end
    Gm = genomes(w)                            # (L, Ploidy, N)
    loc = dropdims(sum(Gm; dims = 2); dims = 2)   # the loci's genotypes, (L, N)
    vend = var(vec(sum(loc; dims = 1)))
    genic = sum(var(view(loc, l, :)) for l in axes(loc, 1))
    merge(out, (vend = vend, genic = genic, ld = vend - genic))
end

function replicates(p, R, seed0; kw...)
    out = Vector{Any}(undef, R)
    Threads.@threads :greedy for r in 1:R
        out[r] = kernel_run(p, SEED_OFFSET + seed0 + r; kw...)
    end
    out
end

cells = [(L, N, x, μ) for μ in (0.01, 1e-4) for N in (1000, 100) for L in (1, 10, 50, 100, 500)
         for x in (L == 1 ? (:qtsim,) : (:qtsim, :all))]

say("QTLSim validate/moments_loci.jl: the multilocus Gaussian moments against the kernel")
say("Julia threads: ", Threads.nthreads(), "; replicates per cell R = $R; ", length(cells), " cells")
say(@sprintf("V_E = %.2f, V_m = %.4f, ω = %.3f, V_s = %.2f; %d generations, window %d:%d:%d",
             envVar(base), muteVar(base), initFitSD(base), selection_variance(base), NGEN,
             first(WINDOW), step(WINDOW), last(WINDOW)))
say()

rows = []
t0 = time()
for (i, (L, N, x, μ)) in enumerate(cells)
    p = Params(base; NoQTL = L, PopSize = N, CrossoverModel = x, MutationRatePerGamete = μ)
    t = @elapsed rs = replicates(p, R, 10_000 * i)
    vk = [r.vg for r in rs]
    e = multilocus_equilibrium(p; drift = true, loci = :multilocus)
    newd, new0 = e.VG, gaussian_vg_loci(p)
    mo = run_moments(Params(p; NoGenerations = SMOKE ? 2000 : 200_000, Sample = 1000); erosion = true, loci = :multilocus).geno_var[end]
    agg, aggd = gaussian_vg_exact(p), run_moments(Params(p; NoGenerations = SMOKE ? 2000 : 200_000, Sample = 1000); loci = :aggregate, erosion = true).geno_var[end]
    ns = mean(r.nsurv for r in rs)
    hoc = house_of_cards_vg(p; nsurv = ns)
    push!(rows, (; L, N, x, μ, vk = mean(vk), sk = se(vk), newd, new0, mo, agg, aggd, hoc,
                 ldk = mean(r.ld for r in rs), ldse = se([r.ld for r in rs]), ldp = e.ld, t))
    @printf("cell %2d/%d  L = %3d  N = %4d  %-6s μ = %-6g  %6.1f s  (%.0f s so far)\n", i, length(cells), L, N, x, μ, t, time() - t0)
end
say(@sprintf("%d cells in %.0f s", length(cells), time() - t0))
say()

say("Equilibrium V_G: kernel (± SE) against the predictions; ratio = kernel/prediction")
say("  L     N  crossover μ        kernel        ±SE   multilocus  ratio   no drift  ratio  aggregate  ratio  aggr+drift ratio   HoC    ratio")
for r in rows
    say(@sprintf("%4d %5d  %-7s %-7g %10.4f %8.4f  %10.4f %6.3f %10.4f %6.3f %10.4f %6.2f %10.4f %6.2f %8.4f %6.2f",
                 r.L, r.N, r.x, r.μ, r.vk, r.sk, r.newd, r.vk / r.newd, r.new0, r.vk / r.new0, r.agg, r.vk / r.agg,
                 r.aggd, r.vk / r.aggd, r.hoc, r.vk / r.hoc))
end
say()
say("LD (V_G − genic variance, zygotes at the end): kernel against multilocus_equilibrium(p; drift = true, loci = :multilocus).ld")
say("  L     N  crossover μ        kernel LD     ±SE   predicted  ratio")
for r in rows
    r.L == 1 && continue
    say(@sprintf("%4d %5d  %-7s %-7g %10.4f %8.4f  %10.4f %6.3f", r.L, r.N, r.x, r.μ, r.ldk, r.ldse, r.ldp, r.ldk / r.ldp))
end
say()

say("Criteria")
check(all(abs(r.mo / r.newd - 1) < 1e-6 for r in rows),
      "the moments lowering with erosion = true converges to gaussian_vg_loci(p; drift = true) in every cell")
G1 = [r for r in rows if r.μ == 0.01 && r.N == 1000 && r.L >= 10]
for r in G1
    check(!outside(r.vk, r.sk, r.newd, 0.10), @sprintf("G1: L = %d, N = %d, %s: kernel V_G %.4f within 10%% of the multilocus %.4f (ratio %.3f)",
                                                       r.L, r.N, r.x, r.vk, r.newd, r.vk / r.newd))
end
for r in [r for r in rows if r.μ == 0.01 && r.L >= 10]
    lk = log(r.vk)
    check(abs(lk - log(r.newd)) < abs(lk - log(r.agg)) && abs(lk - log(r.newd)) < abs(lk - log(r.aggd)),
          @sprintf("G2: L = %d, N = %d, %s: multilocus closer than the aggregate (ratios %.3f, %.2f, %.2f)",
                   r.L, r.N, r.x, r.vk / r.newd, r.vk / r.agg, r.vk / r.aggd))
end
for r in G1
    check(r.ldk < 0 && !outside(r.ldk, r.ldse, r.ldp, 0.25),
          @sprintf("G3: L = %d, N = %d, %s: kernel LD %.4f ± %.4f within 25%% of %.4f", r.L, r.N, r.x, r.ldk, r.ldse, r.ldp))
end
hoc = [r for r in rows if r.μ == 1e-4 && r.L >= 10]
nconf = count(r -> r.vk / r.newd < 0.85, hoc)
say(@sprintf("  expectation (not a criterion): house of cards (μ = 1e-4, L ≥ 10), kernel < 85%% of the multilocus prediction in %d of %d cells: %s",
             nconf, length(hoc), nconf == length(hoc) ? "confirmed" : "refuted"))
say()

# --- S: the steady lag ----------------------------------------------------------------
pS = Params(base; NoQTL = 50, PopSize = 1000, MutationRatePerGamete = 0.01, MoveFitnessSDGen = MOVE,
            MoveType = :continuous, MoveFitnessSD = 0.01)
t = @elapsed rS = replicates(pS, R, 900_000; window = LW, lagwindow = LW)
lk, ls = mean(r.lag for r in rS), se([r.lag for r in rS])
vkS = mean(r.vg for r in rS)
VS = gaussian_vg_loci(pS; drift = true)
λth = steady_lag(pS; VG = VS)
λmo = (o = run_moments(pS; erosion = true, loci = :multilocus); o.theta[end] - o.pheno_mean[end])
λk = steady_lag(pS; VG = vkS)
say(@sprintf("S, steady lag behind a continuous move of 0.01 per generation (%.0f s)", t))
say(@sprintf("  kernel lag %.5f ± %.5f, kernel V_G %.4f (stationary prediction %.4f)", lk, ls, vkS, VS))
say(@sprintf("  steady_lag with the multilocus V_G %.5f (ratio %.3f); moments lowering %.5f (ratio %.3f); with the kernel's V_G %.5f (ratio %.3f)",
             λth, lk / λth, λmo, lk / λmo, λk, lk / λk))
say(@sprintf("  for contrast, steady_lag with gaussian_vg_exact %.5f (ratio %.3f)", steady_lag(pS; VG = gaussian_vg_exact(pS)),
             lk / steady_lag(pS; VG = gaussian_vg_exact(pS))))
check(!outside(lk, ls, λth, 0.10), "S: kernel lag within 10% of steady_lag with gaussian_vg_loci(p; drift = true)")
check(!outside(lk, ls, λmo, 0.10), "S: kernel lag within 10% of the moments lowering's")
check(!outside(lk, ls, λk, 0.10), "S: kernel lag within 10% of steady_lag with the kernel's V_G")
say()

# --- T: two traits --------------------------------------------------------------------
pT = Params(base; NoQTL = 50, PopSize = 1000, MutationRatePerGamete = 0.01, NoTraits = 2, MutCorr = 0.5, SelCorr = -0.3)
t = @elapsed rT = replicates(pT, R, 950_000)
Gk = mean(r.G for r in rT)
Gse = std([r.G for r in rT]) ./ sqrt(R)
Gp = gaussian_G_loci(pT; drift = true)
G0 = gaussian_G_loci(pT)
Ga = gaussian_G(pT)
corr(G) = G[1, 2] / sqrt(G[1, 1] * G[2, 2])
fmt(G) = @sprintf("[%.4f %.4f; %.4f %.4f]", G...)
ang = min(vector_angle(gmax(Gk), gmax(Gp)), vector_angle(gmax(Gk), -gmax(Gp)))
say(@sprintf("T, two traits (%.0f s)", t))
say("  kernel                      ", fmt(Gk), "  (SE ", fmt(Gse), ")")
say("  gaussian_G_loci, drift      ", fmt(Gp))
say("  gaussian_G_loci, no drift   ", fmt(G0))
say("  gaussian_G (aggregate)      ", fmt(Ga))
say(@sprintf("  kernel/multilocus: G11 %.3f, G12 %.3f, G22 %.3f; correlation %.3f vs %.3f; g_max angle %.1f°",
             Gk[1, 1] / Gp[1, 1], Gk[1, 2] / Gp[1, 2], Gk[2, 2] / Gp[2, 2], corr(Gk), corr(Gp), ang))
check(all(abs(Gk[i, j] - Gp[i, j]) - 1.96Gse[i, j] <= 0.10 * sqrt(Gp[i, i] * Gp[j, j]) for i in 1:2, j in 1:2),
      "T: every entry of G within 0.10·√(G_ii G_jj)")
check(abs(corr(Gk) - corr(Gp)) <= 0.05, "T: genetic correlation within 0.05")
check(ang <= 10, "T: g_max within 10°")
say()

say(isempty(fails) ? "PASS: every pre-specified criterion holds" : "FAILED ($(length(fails))): " * join(fails, "; "))
close(io)
exit(isempty(fails) ? 0 : 1)
