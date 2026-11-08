# The rare-alleles closure for several traits and for haploids (src/rare_alleles_traits.jl)
# against the kernel and against the multilocus Gaussian model (src/multilocus.jl).
#
# Usage (from the package root; Catlab comes from the vignettes environment):
#     julia -t 3 --project=vignettes validate/rare_alleles_traits.jl [R]
# R replicates per cell (default 24; the deme scenarios use 2R). SMOKE=1 runs 2
# replicates of short runs into a temporary log. With ENV SEED_OFFSET ≠ 0 the kernel seeds
# are shifted (an independent confirmation, validate/rare_alleles_traits_seed<OFFSET>.log).
# Writes validate/rare_alleles_traits.log; exits 1 if a pre-specified criterion fails.
#
# DESIGN
# T, two traits, diploid, N = 1000, RecombPoissonMean = 3, V_E = 0.01 per trait,
#   V_m = 2e-4 (MutationScalar 0.02), ω = 2 (A ≈ 4 I + Ω's correlation), mutational
#   covariance M = V_m·[1 0.5; 0.5 0.6] and fitness correlation SelCorr = −0.4 (Ω's
#   eigenvectors (1, ±1)/√2 are not M's), L ∈ {20, 100} loci, μ ∈ {1e-5, 1e-4, 1e-3, 1e-2}
#   per locus per gamete, i.e. 4Nμ = 0.04, 0.4, 4, 40 (house of cards with drift to
#   Gaussian). 8000 generations from zero genomes; statistics averaged over generations
#   3000–7990 (every 10) and the replicates: the kernel's G (gmatrix), against the same
#   averages of run_moments(p; erosion = true) with loci = :rare_alleles (RA) and
#   loci = :multilocus (ML).
# H, one trait, haploid (clonal reproduction), N = 1000, the same V_E, V_m, ω, L and μ;
#   the kernel's geno_var against RA (the clonal closure) and ML.
# D, k = 4 demes of N = 100 diploids, 50 loci, two traits (MutCorr 0.3, SelCorr 0.2),
#   V_E = 1, V_m = 1e-3, ω = 5, island migration 0.01, local optima
#   θ_a = 2(cos φ_a, sin φ_a), φ = 0°, 60°, 120°, 180° (validate/trait_moments_demes.jl's
#   T1), at μ = 2e-4 (D1, 4Nμ = 0.08, μV_s/α² ≈ 0.1: house of cards) and μ = 0.02 (D2,
#   Gaussian). 4000 generations, census 2000:10:3990: mean within-deme G_W (kernel sample
#   covariance; moments `between = true`, metapop.within).
#
# Error measures. For a 2 × 2 prediction G and the kernel's mean Ĝ (SE per entry over the
#   replicates): diagonal entries relative to G_ii, the off-diagonal relative to
#   √(G_11 G_22); the genetic correlation r = G_12/√(G_11 G_22) (difference); the angle
#   between the leading eigenvectors (g_max). A statistic FAILS only if the kernel's whole
#   95% CI lies outside the tolerance band (|Ĝ − G| − 1.96 SE > tol·scale), as in
#   validate/moments_demes.jl; the correlation and angle are tested on the point values.
#   e = ‖Ĝ − G‖_F/‖Ĝ‖_F (one trait: |V̂ − V|/V̂) compares the models.
#
# PRE-SPECIFIED CRITERIA (written before any kernel run of this script):
# T-A (accuracy, every T cell): RA's diagonal entries within 20%, off-diagonal within
#     10% of √(G_11 G_22), |Δr| ≤ 0.10, and the g_max angle ≤ 10° when RA's eigenvalue
#     ratio λ1/λ2 ≥ 1.5 (otherwise g_max is not tested).
# T-B (never worse, every T cell): e_RA ≤ e_ML + 0.03.
# T-C (much better, T cells with 4Nμ ≤ 0.4): e_RA ≤ 0.5·e_ML.
# H-A (every H cell): RA's V_G within 20% (CI rule).
# H-B (every H cell): e_RA ≤ e_ML + 0.03.   H-C (H cells with 4Nμ ≤ 0.4): e_RA ≤ 0.5·e_ML.
# D-A (D1, D2): RA's G_W entries within 15% (CI rule; off-diagonal relative to
#     √(G_11 G_22)).   D-B: e_RA ≤ e_ML + 0.03.
# Decision rule for the default (fixed in advance): `loci = :auto` may resolve to
# :rare_alleles for several traits if T-A, T-B, T-C, D-A and D-B hold, and for haploids
# if H-A, H-B and H-C hold.
#
# HISTORY. A SMOKE run (2 replicates, 1500 generations, not equilibrated) preceded run 1;
# nothing was changed after it. Run 1 (2026-11-08, R = 24, SEED_OFFSET = 0,
# rare_alleles_traits.log, 347 s on 3 threads): every T and H criterion holds (T: e_RA
# 0.02–0.10 against e_ML 0.02–77, ML/RA 1.25–1150; H: |error| ≤ 3.5%, ML/RA 0.16–795 —
# in the Gaussian H cells both are within 2.3%), D2 holds, and D1 FAILS (D-A G11, G22 and
# D-B): with local optima and migration in the house-of-cards regime the kernel's
# within-deme G (0.84, 0.55) is the multilocus model's (0.89, 0.63), and RA's is less than
# half (0.34, 0.29). The closure takes every mutant to be deleterious wherever it is,
# removed at a²/(2s); under local adaptation an allele is favoured in some demes and
# migration keeps it common (the same kind of failure as moments_regimes.jl's S, a moving
# optimum). The one-trait closure (the current default) behaves the same way in the
# moments (one-trait analogue of D1: RA 0.50, ML 0.87; no kernel run). Forcing β = 1
# raises RA only to 0.44 (diagnostic, moments only). No model was changed.
# Confirmation (SEED_OFFSET = 7000000, rare_alleles_traits_seed7000000.log, 295 s): the same
# verdicts, the same failures (D1 G11, G22, D-B) and nothing else.
# Decision: haploids supported; several traits NOT supported by the rule (D1).

using QTLSim, Random, Statistics, Printf, LinearAlgebra
import Catlab
using Catlab: ⋅, ⊗

const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 24
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const SEED_OFFSET = parse(Int, get(ENV, "SEED_OFFSET", "0"))
const LOG = SMOKE ? tempname() : joinpath(@__DIR__, SEED_OFFSET == 0 ? "rare_alleles_traits.log" : "rare_alleles_traits_seed$(SEED_OFFSET).log")
const NREP = SMOKE ? 2 : R
const GENS = SMOKE ? 1500 : 8000
const WINDOW = SMOKE ? (500:10:1490) : (3000:10:7990)

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s); flush(stdout))
const failures = String[]
fail!(x) = push!(failures, x)

const BASE = Params(PopSize = 1000, Ploidy = 2, NoQTL = 20, PhenotypeVar = 0.02, Heritability = 0.5,
                    MutationScalar = 0.02, MutationRatePerGamete = 1e-3, InitFitnessSD = 2.0,
                    RecombPoissonMean = 3.0, NoGenerations = GENS, Sample = 10)
const MUS = [1e-5, 1e-4, 1e-3, 1e-2]
const LS = [20, 100]

tcell(L, μ) = Params(BASE; NoQTL = L, MutationRatePerGamete = μ, NoTraits = 2,
                     MutCov = muteVar(BASE) .* [1.0 0.5; 0.5 0.6], SelCorr = -0.4)
hcell(L, μ) = Params(BASE; NoQTL = L, MutationRatePerGamete = μ, Ploidy = 1)

relF(K, G) = norm(K - G) / norm(K)
gcorr(G) = G[1, 2] / sqrt(G[1, 1] * G[2, 2])
function angle(G1, G2)
    a, b = gmax(G1), gmax(G2)
    acosd(clamp(abs(dot(a, b)), 0.0, 1.0))
end
eigratio(G) = (e = eigvals(Symmetric(G)); maximum(e) / max(minimum(e), 1e-300))

# the kernel's mean over the window, one replicate
function kernel_G(p, seed)
    o = simulate(p; seed)
    o.extinct && error("extinct (seed $seed)")
    idx = findall(in(WINDOW), o.gen)
    p.NoTraits > 1 ? mean(o.gmatrix[idx]) : fill(mean(o.geno_var[idx]), 1, 1)
end

function moments_G(p, loci)
    mo = run_moments(p; erosion = true, loci)
    idx = findall(in(WINDOW), mo.gen)
    p.NoTraits > 1 ? mean(mo.gmatrix[idx]) : fill(mean(mo.geno_var[idx]), 1, 1)
end

function runcells(cells, seed0)
    jobs = [(i, r) for i in eachindex(cells) for r in 1:NREP]
    res = Vector{Matrix{Float64}}(undef, length(jobs))
    Threads.@threads :greedy for j in eachindex(jobs)
        i, r = jobs[j]
        res[j] = kernel_G(cells[i], SEED_OFFSET + seed0 + 1000i + r)
    end
    [[res[j] for j in eachindex(jobs) if jobs[j][1] == i] for i in eachindex(cells)]
end

# the CI rule: FAIL only when the kernel's whole 95% CI is outside the band
cifail(k, se, pred, scale, tol) = abs(k - pred) - 1.96se > tol * scale

function two_trait_check(tag, Ks, RA, ML; tol_d, tol_o, test_r_angle)
    n = length(Ks)
    K = mean(Ks)
    SE = [std(getindex.(Ks, a, b)) / sqrt(n) for a in 1:2, b in 1:2]
    ok = true
    for (a, b) in ((1, 1), (1, 2), (2, 2))
        scale = sqrt(RA[a, a] * RA[b, b]); tol = a == b ? tol_d : tol_o
        f = cifail(K[a, b], SE[a, b], RA[a, b], scale, tol)
        f && (ok = false; fail!("$tag G$a$b"))
        say(@sprintf("    G%d%d  kernel %.4g ± %.2g   RA %.4g (%+.1f%%)   ML %.4g (%+.1f%%)  %s", a, b, K[a, b], SE[a, b],
                     RA[a, b], 100(RA[a, b] - K[a, b]) / scale, ML[a, b], 100(ML[a, b] - K[a, b]) / sqrt(ML[a, a] * ML[b, b]),
                     f ? "FAIL" : "ok"))
    end
    if test_r_angle
        dr = gcorr(K) - gcorr(RA)
        fr = abs(dr) > 0.10
        fr && (ok = false; fail!("$tag r"))
        ang = angle(K, RA); er = eigratio(RA)
        fa = er >= 1.5 && ang > 10
        fa && (ok = false; fail!("$tag g_max"))
        say(@sprintf("    r     kernel %.3f   RA %.3f   ML %.3f  %s", gcorr(K), gcorr(RA), gcorr(ML), fr ? "FAIL" : "ok"))
        say(@sprintf("    g_max angle to kernel: RA %.1f°, ML %.1f° (RA λ1/λ2 = %.2f%s)  %s", ang, angle(K, ML), er,
                     er >= 1.5 ? "" : ", not tested", fa ? "FAIL" : "ok"))
    end
    eRA, eML = relF(K, RA), relF(K, ML)
    (; ok, eRA, eML, K)
end

function main()
    say("QTLSim validate/rare_alleles_traits.jl")
    say("threads ", Threads.nthreads(), ", replicates per cell ", NREP, ", SEED_OFFSET ", SEED_OFFSET, SMOKE ? " (SMOKE)" : "")
    t0 = time()
    summary = Tuple{String,Float64,Float64,Float64}[]

    # --- T: two traits ------------------------------------------------------------
    say("\nT: two traits, diploid, N = 1000; M = V_m[1 0.5; 0.5 0.6], SelCorr = −0.4, ω = 2; window $(first(WINDOW))–$(last(WINDOW))")
    tc = [(L, μ) for L in LS for μ in MUS]
    tps = [tcell(L, μ) for (L, μ) in tc]
    tk = runcells(tps, 1_000_000)
    for (i, (L, μ)) in enumerate(tc)
        p = tps[i]
        RA, ML = moments_G(p, :rare_alleles), moments_G(p, :multilocus)
        θ4 = 4 * p.PopSize * μ
        α2 = QTLSim.mutation_effect_variance(Params(p; NoTraits = 1, MutCov = zeros(0, 0), SelCorr = 0.0))
        say(@sprintf("\n  T L=%d μ=%g (4Nμ = %g, α² = %.3g, μV_s/α² = %.3g)", L, μ, θ4, α2, μ * selection_variance(p) / α2))
        r = two_trait_check("T L=$L μ=$μ", tk[i], RA, ML; tol_d = 0.20, tol_o = 0.10, test_r_angle = true)
        fB = r.eRA > r.eML + 0.03
        fC = θ4 <= 0.4 && r.eRA > 0.5r.eML
        fB && fail!("T-B L=$L μ=$μ"); fC && fail!("T-C L=$L μ=$μ")
        say(@sprintf("    e_RA = %.3f, e_ML = %.3f   T-B %s%s", r.eRA, r.eML, fB ? "FAIL" : "ok",
                     θ4 <= 0.4 ? "   T-C " * (fC ? "FAIL" : "ok") : ""))
        push!(summary, ("T L=$L μ=$μ", θ4, r.eRA, r.eML))
    end

    # --- H: haploid, one trait --------------------------------------------------------
    say("\nH: one trait, haploid (clones), N = 1000, ω = 2")
    hps = [hcell(L, μ) for (L, μ) in tc]
    hk = runcells(hps, 2_000_000)
    for (i, (L, μ)) in enumerate(tc)
        p = hps[i]
        Ks = [k[1, 1] for k in hk[i]]
        K, se = mean(Ks), std(Ks) / sqrt(length(Ks))
        RA, ML = moments_G(p, :rare_alleles)[1, 1], moments_G(p, :multilocus)[1, 1]
        θ4 = 4 * p.PopSize * μ
        f = cifail(K, se, RA, RA, 0.20)
        eRA, eML = abs(K - RA) / K, abs(K - ML) / K
        fB = eRA > eML + 0.03
        fC = θ4 <= 0.4 && eRA > 0.5eML
        f && fail!("H-A L=$L μ=$μ"); fB && fail!("H-B L=$L μ=$μ"); fC && fail!("H-C L=$L μ=$μ")
        say(@sprintf("  H L=%3d μ=%-6g kernel %.4g ± %.2g  RA %.4g (%+.1f%%) %s  ML %.4g (%+.1f%%)  HoG %.4g   H-B %s%s",
                     L, μ, K, se, RA, 100(RA - K) / K, f ? "FAIL" : "ok", ML, 100(ML - K) / K,
                     house_of_gauss_vg(p), fB ? "FAIL" : "ok", θ4 <= 0.4 ? "  H-C " * (fC ? "FAIL" : "ok") : ""))
        push!(summary, ("H L=$L μ=$μ", θ4, eRA, eML))
    end

    # --- D: four demes, two traits --------------------------------------------------
    say("\nD: k = 4 demes of N = 100, 50 loci, two traits, local optima, island migration 0.01")
    LC = lifecycle()
    K4, N4, RATE = 4, 100, 0.01
    PHI = [0.0, 60.0, 120.0, 180.0]
    DG = SMOKE ? 1000 : 4000
    CENSUS = collect((SMOKE ? 500 : 2000):10:DG-10)
    dbase = Params(Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                   MutationRatePerGamete = 0.02, PopSize = N4, NoGenerations = DG, Sample = 10, InitFitnessSD = 5.0,
                   NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)
    dcycle(p) = foldl(⊗, [withparams(LC.generation, p; θ = 2.0 .* [cosd(φ), sind(φ)]) for φ in PHI]) ⋅ LC.migrate(K4)
    function dkernel(p, seed)
        st = p ⊗ demes(K4)
        idx = IdDict(w => i for (i, w) in enumerate(st.demes))
        gw = [zeros(2, 2) for _ in 1:K4, _ in CENSUS]
        rec = operator(:record, function (w, rng)
            j = searchsortedfirst(CENSUS, w.gen)
            if j <= length(CENSUS) && CENSUS[j] == w.gen
                gw[idx[w], j] = cov(traitvalues(w; what = :genotype); dims = 2)
            end
            w
        end)
        o = run_composed(st, Xoshiro(seed); cycle = dcycle(p), record = false, images = Dict(:record => rec),
                         migration_rate = RATE)
        any(o.extinct_gen .>= 0) && error("a deme died")
        mean(gw)
    end
    function dmoments(p, loci)
        mo = run_moments(p; cycle = dcycle(p), migration_rate = RATE, between = true, loci)
        mean(mo.metapop.within[findall(in(CENSUS), mo.metapop.gen)])
    end
    for (tag, μ) in (("D1", 2e-4), ("D2", 0.02))
        p = Params(dbase; MutationRatePerGamete = μ)
        jobs = 1:2NREP
        Ks = Vector{Matrix{Float64}}(undef, length(jobs))
        Threads.@threads :greedy for r in jobs
            Ks[r] = dkernel(p, SEED_OFFSET + 3_000_000 + (tag == "D1" ? 0 : 10_000) + r)
        end
        RA, ML = dmoments(p, :rare_alleles), dmoments(p, :multilocus)
        say(@sprintf("\n  %s μ = %g (4Nμ = %g per deme)", tag, μ, 4N4 * μ))
        r = two_trait_check(tag, Ks, RA, ML; tol_d = 0.15, tol_o = 0.15, test_r_angle = false)
        fB = r.eRA > r.eML + 0.03
        fB && fail!("D-B $tag")
        say(@sprintf("    e_RA = %.3f, e_ML = %.3f   D-B %s", r.eRA, r.eML, fB ? "FAIL" : "ok"))
        push!(summary, (tag, 4N4 * μ, r.eRA, r.eML))
    end

    say("\nSummary (relative error of the prediction against the kernel):")
    say(@sprintf("  %-18s %8s %8s %8s %8s", "cell", "4Nμ", "e_RA", "e_ML", "ML/RA"))
    for (c, θ, a, b) in summary
        say(@sprintf("  %-18s %8g %8.3f %8.3f %8.2f", c, θ, a, b, b / a))
    end
    say(@sprintf("\nkernel runs and predictions in %.0f s", time() - t0))
    tfail = any(f -> startswith(f, "T") || startswith(f, "D"), failures)
    hfail = any(f -> startswith(f, "H"), failures)
    say("decision rule: several traits ", tfail ? "NOT supported" : "supported", "; haploids ", hfail ? "NOT supported" : "supported")
    say(isempty(failures) ? "PASS: every pre-specified criterion holds" : "FAIL: " * join(failures, "; "))
    write(LOG, String(take!(out)))
    isempty(failures) || exit(1)
end

main()
