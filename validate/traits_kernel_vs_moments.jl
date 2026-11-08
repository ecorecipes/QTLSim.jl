# Several traits: the kernel against the trait moments lowering and its theory.
#
#     julia -t 3 --project=vignettes validate/traits_kernel_vs_moments.jl [R]
#
# (the vignettes environment has QTLSim and Catlab). R replicates per scenario, 12 by
# default. Writes validate/traits_kernel_vs_moments.log and exits 1 if a pre-specified
# criterion fails.
#
# The moments recursion G' = G − G (A + G)⁻¹ G + M treats selection's reduction of the
# variance as lasting, which is exact for one haploid locus with a continuum of alleles
# (Kimura 1965; Lande 1975), so the criteria are tested there, with 2000 individuals, a
# mutation rate of 0.5 per gamete and μV_s/α² ≈ 125 (the Gaussian regime). M and Ω are
# chosen not to share eigenvectors.
#
# Pre-specified criteria (fixed before the first run):
#   E, equilibrium: the kernel's G, averaged over generations 2000–5999 and the replicates,
#      against the moments lowering with drift's erosion: every entry within
#      0.10·√(G_ii G_jj) of it, the genetic correlation within 0.05, g_max within 10°.
#   L, steady lag (continuous move by 0.001 per generation along (1, 2) from generation
#      2000): the lag θ − z̄ averaged over generations 4000–6999 within 10% (in norm) of
#      steady_lag_vector with the kernel's own G over the same window, and of the moments
#      lowering's lag, and within 5° of both in direction.
#   C, one trait: the scalar kernel in scenario E's regime, as the calibration of what
#      the Gaussian approximation misses (no criterion; the trait model should miss by as
#      much).
# Informational, no criterion: 50 diploid loci, where the aggregate recursion does not
# apply (for one trait as for two), against Lande's per-locus sum gaussian_G_loci.

using QTLSim, Statistics, Printf, Random
import Catlab
using QTLSim: symsqrt

R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
logpath = joinpath(@__DIR__, "traits_kernel_vs_moments.log")
io = open(logpath, "w")
say(s...) = (println(stdout, s...); println(io, s...); flush(io))

base = Params(PopSize = 2000, NoQTL = 1, Ploidy = 1, NoGenerations = 6000, PhenotypeVar = 0.02,
              Heritability = 0.5, MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)
pE = Params(base; NoTraits = 2, MutCov = muteVar(base) .* [1.0 0.4; 0.4 0.5], SelCorr = -0.3)
pL = Params(pE; NoGenerations = 7000, MoveFitnessSDGen = 2000, MoveType = :continuous, MoveFitnessSD = 0.001,
            MoveDirection = [1.0, 2.0])

function runs(p, R, window; seed0)
    out = Vector{Any}(undef, R)
    Threads.@threads :greedy for r in 1:R
        o = simulate(p; seed = seed0 + r)
        o.extinct && error("replicate $r went extinct")
        idx = findall(g -> g in window, o.gen)
        out[r] = haskey(o, :gmatrix) ?
            (G = mean(o.gmatrix[idx]), lag = mean(o.trait_lag[idx]), vg = mean(o.geno_var[idx])) :
            (G = fill(mean(o.geno_var[idx]), 1, 1), lag = [mean(o.theta[idx] - o.pheno_mean[idx])], vg = mean(o.geno_var[idx]))
    end
    out
end

corr(G) = G[1, 2] / sqrt(G[1, 1] * G[2, 2])
se(xs) = std(xs) / sqrt(length(xs))
fmt(G) = @sprintf("[%.5f %.5f; %.5f %.5f]", G[1, 1], G[1, 2], G[2, 1], G[2, 2])
fails = String[]
check(ok, what) = (say(ok ? "  PASS  " : "  FAIL  ", what); ok || push!(fails, what))

say("Several traits: kernel against the trait moments lowering, $R replicates per scenario, ",
    Threads.nthreads(), " threads")
say("M = ", fmt(mutational_covariance(pE)), ", Ω = ", fmt(selection_covariance(pE)), ", V_E = ", envVar(pE))
say()

# --- E: equilibrium G ---------------------------------------------------------
t = @elapsed rE = runs(pE, R, 2000:5999; seed0 = 1000)
Gk = mean(r.G for r in rE)
Gm = run_moments(pE; erosion = true).gmatrix[end]
G0 = gaussian_G(pE)
say("E, equilibrium G ($(round(t; digits = 1)) s)")
say("  kernel        ", fmt(Gk), "   (SE of entries ", fmt(std([r.G for r in rE]) ./ sqrt(R)), ")")
say("  moments+drift ", fmt(Gm))
say("  gaussian_G    ", fmt(G0), "   (no drift)")
say(@sprintf("  kernel/moments: G11 %.3f, G12 %.3f, G22 %.3f; correlation %.3f vs %.3f; g_max angle %.2f°",
             Gk[1, 1] / Gm[1, 1], Gk[1, 2] / Gm[1, 2], Gk[2, 2] / Gm[2, 2], corr(Gk), corr(Gm),
             min(vector_angle(gmax(Gk), gmax(Gm)), vector_angle(gmax(Gk), -gmax(Gm)))))
check(all(abs(Gk[i, j] - Gm[i, j]) <= 0.10 * sqrt(Gm[i, i] * Gm[j, j]) for i in 1:2, j in 1:2),
      "E: every entry of G within 0.10·√(G_ii G_jj)")
check(abs(corr(Gk) - corr(Gm)) <= 0.05, "E: genetic correlation within 0.05")
check(min(vector_angle(gmax(Gk), gmax(Gm)), vector_angle(gmax(Gk), -gmax(Gm))) <= 10, "E: g_max within 10°")

# --- C: the one-trait calibration -------------------------------------------
pC = base
rC = runs(pC, R, 2000:5999; seed0 = 3000)
vk = mean(r.vg for r in rC)
vm = run_moments(pC; erosion = true).geno_var[end]
say(@sprintf("C, one trait in the same regime: kernel V_G %.6f ± %.6f, moments+drift %.6f, ratio %.3f (two traits: %.3f, %.3f)",
             vk, se([r.vg for r in rC]), vm, vk / vm, Gk[1, 1] / Gm[1, 1], Gk[2, 2] / Gm[2, 2]))
say()

# --- L: steady lag ------------------------------------------------------------
t = @elapsed rL = runs(pL, R, 4000:6999; seed0 = 2000)
λk = mean(r.lag for r in rL)
GkL = mean(r.G for r in rL)
λG = steady_lag_vector(pL; G = GkL)
oL = run_moments(pL; erosion = true)
λm = oL.trait_lag[end]
say("L, steady lag behind an optimum moving by 0.001 along (1, 2)/√5 ($(round(t; digits = 1)) s)")
say(@sprintf("  kernel lag              (%.5f, %.5f) ± (%.5f, %.5f)", λk..., se([r.lag[1] for r in rL]), se([r.lag[2] for r in rL])))
say(@sprintf("  steady_lag_vector(G_k)  (%.5f, %.5f)   G_k = %s", λG..., fmt(GkL)))
say(@sprintf("  moments+drift           (%.5f, %.5f)", λm...))
say(@sprintf("  the move's direction is %.1f° from trait 1; the lag's %.1f° (deflected away from g_max, %.1f°)",
             vector_angle(move_direction(pL), [1.0, 0.0]), vector_angle(λk, [1.0, 0.0]),
             vector_angle(gmax(GkL), [1.0, 0.0])))
check(QTLSim.norm(λk - λG) <= 0.10 * QTLSim.norm(λG), "L: lag within 10% of steady_lag_vector with the kernel's G")
check(QTLSim.norm(λk - λm) <= 0.10 * QTLSim.norm(λm), "L: lag within 10% of the moments lowering's")
check(vector_angle(λk, λG) <= 5 && vector_angle(λk, λm) <= 5, "L: lag direction within 5° of both")
say()

# --- informational: many diploid loci -----------------------------------------
pM = Params(PopSize = 1000, NoQTL = 50, Ploidy = 2, NoGenerations = 6000, PhenotypeVar = 1.0, Heritability = 0.5,
            MutationScalar = 0.01, MutationRatePerGamete = 0.01, NoTraits = 2, MutCorr = 0.5, SelCorr = -0.3)
pM1 = Params(pM; NoTraits = 1, MutCorr = 0.0, SelCorr = 0.0)
Rm = max(3, R ÷ 4)
t = @elapsed (rM = runs(pM, Rm, 3000:5999; seed0 = 4000); rM1 = runs(pM1, Rm, 3000:5999; seed0 = 5000))
GM = mean(r.G for r in rM)
say("Informational: 50 diploid loci, N = 1000, $Rm replicates ($(round(t; digits = 1)) s)")
say("  two traits: kernel ", fmt(GM), ", aggregate gaussian_G ", fmt(gaussian_G(pM)),
    ", per-locus gaussian_G_loci ", fmt(gaussian_G_loci(pM)))
say(@sprintf("  one trait:  kernel V_G %.4f, gaussian_vg_exact %.4f, per-locus %.4f",
             mean(r.vg for r in rM1), gaussian_vg_exact(pM1), gaussian_G_loci(pM1)[1, 1]))
say(@sprintf("  kernel/aggregate: one trait %.2f, two traits %.2f and %.2f; genetic correlation kernel %.3f, theory %.3f",
             mean(r.vg for r in rM1) / gaussian_vg_exact(pM1), GM[1, 1] / gaussian_G(pM)[1, 1],
             GM[2, 2] / gaussian_G(pM)[2, 2], corr(GM), corr(gaussian_G(pM))))
say()

say(isempty(fails) ? "All pre-specified criteria pass." : "FAILED: " * join(fails, "; "))
close(io)
exit(isempty(fails) ? 0 : 1)
