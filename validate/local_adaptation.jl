# The rare-alleles closure with locally favoured alleles (`loci = :rare_local`,
# src/local_alleles.jl) against the kernel, the current rare-alleles closure
# (`loci = :rare_alleles`, src/rare_alleles.jl and src/rare_alleles_traits.jl) and the
# multilocus Gaussian model (`loci = :multilocus`, src/multilocus.jl), for metapopulations
# with local optima and migration.
#
# Usage (from the package root; Catlab comes from the vignettes environment):
#     julia -t 3 --project=vignettes validate/local_adaptation.jl [R]
# R replicates per one-trait cell (default 30; the two-trait cells use 36). SMOKE=1 runs
# 2 replicates of short runs into a temporary log. With ENV SEED_OFFSET ≠ 0 the kernel
# seeds are shifted (validate/local_adaptation_seed<OFFSET>.log). Writes
# validate/local_adaptation.log; exits 1 if a pre-specified criterion fails.
#
# DESIGN
# B and G, one trait (the baseline for the current default, `:auto` → `:rare_alleles`):
#   k = 4 demes of N = 100 diploids, 50 loci, RecombPoissonMean = 1 (default), V_E = 1,
#   V_m = 1e-3 (MutationScalar 1e-3), ω = 5 (V_s ≈ 26), island migration at forward rate
#   m (m/3 toward each other deme), optima θ = (2, 1, −1, −2) (the one-trait analogue
#   of validate/rare_alleles_traits.jl's D1), at μ = 2e-4 per locus per gamete (B,
#   house of cards with drift: 4Nμ = 0.08 per deme, α² = 0.05, μV_s/α² ≈ 0.1) and
#   μ = 0.02 (G, Gaussian: 4Nμ = 8), for m = 0.001, 0.01, 0.1 (cells 1–3), and a
#   common optimum θ = 0 at m = 0.01 (cell 4). 6000 generations from zero genomes;
#   census 3000:10:5990. Statistic: the within-deme genetic variance V_W, the kernel's
#   sample variance of the adults' genotypes at QTSim's census (after migration), averaged
#   over demes, census generations and replicates; the predictions are the same averages
#   of run_moments(...; between = true).metapop.within with each `loci`.
# D1, D2, two traits: validate/rare_alleles_traits.jl's D cells unchanged (k = 4,
#   N = 100, 50 loci, MutCorr 0.3, SelCorr 0.2, V_E = 1, V_m = 1e-3, ω = 5, island
#   migration 0.01, optima 2(cos φ, sin φ), φ = 0°, 60°, 120°, 180°; μ = 2e-4 (D1,
#   house of cards) and 0.02 (D2); 4000 generations, census 2000:10:3990; mean
#   within-deme G_W), on fresh seeds.
# I, identity cells (moments only): one deme (one trait and two traits, μ = 2e-4 and
#   0.02, with `erosion = true` and with `between = true`), four demes without migration
#   (local optima), and four identical demes with a common optimum (m = 0.01).
#
# Error measures as in validate/rare_alleles_traits.jl: for V_W (one trait) the relative
#   error e = |V̂ − V|/V̂ against the kernel's mean V̂; for G_W (two traits) the diagonal
#   entries relative to G_ii, the off-diagonal relative to √(G_11 G_22) of the
#   prediction, and e = ‖Ĝ − G‖_F/‖Ĝ‖_F. A tolerance on an entry FAILS only if the
#   kernel's whole 95% CI lies outside the band (|Ĝ − G| − 1.96 SE > tol·scale).
#
# PRE-SPECIFIED CRITERIA (written before any kernel run of this script):
# I    (every identity cell): `:rare_local`'s recorded predictions are identical (==) to
#      `:rare_alleles`'s.
# A    (accuracy of `:rare_local`, every kernel cell): V_W within 15% (CI rule); in D1
#      and D2 every entry of G_W within 15% (CI rule).
# R1   (decision rule, no worse, every kernel cell): e_RL ≤ e_ML + 0.03.
# R2   (decision rule, better, every house-of-cards cell, B1–B4 and D1): e_RL < e_ML.
# Baseline (reported, does not set the exit code): the current closure `:rare_alleles`
#      FAILS a cell if V_W (or an entry of G_W) misses the kernel by more than 15% (CI rule)
#      or e_RA > e_ML + 0.03.
# Decision rule for the defaults (fixed in advance, from the request):
# - one trait: if R1 and R2 hold in every one-trait cell (B, G), `loci = :auto` resolves to
#   `:rare_local` for metapopulations of one trait (it is `:rare_alleles` for one deme);
#   if the baseline fails a one-trait local-optimum cell and `:rare_local` does not fix
#   it (A, R1 or R2 fails there), `:auto` falls back to `:multilocus` for metapopulations
#   with differing optima;
# - several traits: if A, R1 and R2 hold in D1 and D2, `:auto` resolves to `:rare_local`
#   for several traits (validate/rare_alleles_traits.jl's T cells, one population, already
#   pass with `:rare_alleles`, which `:rare_local` equals for one deme).
# Either change goes in a separate commit for the user to decide on.
#
# DISCLOSURE (what was looked at before these criteria were written): the kernel's D1
# and D2 results of validate/rare_alleles_traits.jl (D1 within-deme G 0.84, 0.13, 0.55;
# D2 0.81, 0.15, 0.62), on other seeds; the moments' predictions of every cell above,
# for all three closures (one trait: B1–B4 RA 0.445, 0.510, 0.425, 0.409; RL 0.644,
# 0.916, 0.425, 0.409; ML 0.644, 0.916, 0.663, 0.582; D1 RL 0.855, 0.088, 0.488).
# A first candidate (a second genic class for the variance that migration brings with
# the differences of the deme means, left out of the closure) was tried in the moments
# only (D1 G11 0.37, one-trait analogue 0.54: little change) and dropped before this
# script was written. No kernel run of the cells B and G had been made.
#
# HISTORY. A SMOKE run (2 replicates, 1500 generations, not equilibrated; seeds 5_000_001…,
# 6_000_001…) preceded run 1 and is disclosed here: RL failed A in B1, B2 and D1, R1/R2 in
# B2, B3, B4; the baseline RA failed B2, B3, B4 and D1 (B4, a common optimum, by −16%).
# Nothing was changed after it but the replicate counts (40 → 30 one-trait, 48 → 36
# two-trait, for run time); the criteria are as written above.
# Run 1 (2026-11-08, R = 30/36, SEED_OFFSET = 0, local_adaptation.log, 3893 s on 3 threads):
# every identity cell identical. Baseline: the current closure FAILS B1 (−18%), B2 (−45%)
# and D1 (G11 −148%, G22 −85% of scale; the D1 failure of rare_alleles_traits.jl
# reproduced on fresh seeds); it passes B3, B4 and every Gaussian cell. `:rare_local`
# passes A in every cell (B1 +19%, B2 −0.8%, B3 −17%, B4 +4.3%, G within 4.5%, D1 entries
# +0.8%, −4.6%, −9.1%, D2 within 7.7%) and R1 in every cell; R2 FAILS in B2 only, on a tie:
# there the gain switches the rare loss off entirely, so `:rare_local` equals the
# multilocus model to 0.03% (e_RL = 0.008, e_ML = 0.007). Elsewhere e_RL ≤ e_ML (B3 0.170
# against 0.293, B4 0.043 against 0.485, D1 0.061 against 0.123). The rule as written
# gives: several traits supported (`:auto` → `:rare_local`); one trait not supported, and
# since the baseline fails B2 and R2 fails there, `:auto` → `:multilocus` for
# metapopulations with differing optima. The criteria were not changed; no model change
# followed, so there was no rerun.
using QTLSim, Random, Statistics, Printf, LinearAlgebra
import Catlab
using Catlab: ⋅, ⊗

const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 30
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const SEED_OFFSET = parse(Int, get(ENV, "SEED_OFFSET", "0"))
const LOG = SMOKE ? tempname() : joinpath(@__DIR__, SEED_OFFSET == 0 ? "local_adaptation.log" : "local_adaptation_seed$(SEED_OFFSET).log")
const NREP = SMOKE ? 2 : R
const NREP2 = SMOKE ? 2 : 36
const LC = lifecycle()

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s); flush(stdout))
const failures = String[]
const baseline = String[]
fail!(x) = push!(failures, x)
cifail(k, se, pred, scale, tol) = abs(k - pred) - 1.96se > tol * scale
const LOCI = (:rare_alleles, :rare_local, :multilocus)

# --- one trait ----------------------------------------------------------------------
const GENS1 = SMOKE ? 1500 : 6000
const CENSUS1 = collect((SMOKE ? 500 : 3000):10:GENS1-10)
const base1 = Params(Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                     MutationRatePerGamete = 2e-4, PopSize = 100, NoGenerations = GENS1, Sample = 10, InitFitnessSD = 5.0)
const LOCAL = [2.0, 1.0, -1.0, -2.0]
cycle1(p, θs) = foldl(⊗, [withparams(LC.generation, p; θ = t) for t in θs]) ⋅ LC.migrate(length(θs))

const CELLS1 = [(tag = "$(r)$(i)", μ, m, θs, hoc = r == "B")
                for (r, μ) in (("B", 2e-4), ("G", 0.02))
                for (i, (m, θs)) in enumerate(((0.001, LOCAL), (0.01, LOCAL), (0.1, LOCAL), (0.01, zeros(4))))]

function kernel1(p, θs, m, seed)
    k = length(θs)
    st = p ⊗ demes(k)
    idx = IdDict(w => i for (i, w) in enumerate(st.demes))
    vg = zeros(k, length(CENSUS1))
    rec = operator(:record, function (w, rng)
        j = searchsortedfirst(CENSUS1, w.gen)
        if j <= length(CENSUS1) && CENSUS1[j] == w.gen
            vg[idx[w], j] = var(view(w.genotype, 1:w.N))
        end
        w
    end)
    o = run_composed(st, Xoshiro(seed); cycle = cycle1(p, θs), record = false, images = Dict(:record => rec),
                     migration_rate = m)
    any(o.extinct_gen .>= 0) && error("a deme died (seed $seed)")
    mean(vg)
end

function moments1(p, θs, m, loci; census = CENSUS1)
    mo = run_moments(p; cycle = cycle1(p, θs), migration_rate = m, between = true, loci)
    mean(mo.metapop.within[findall(in(census), mo.metapop.gen)])
end

# --- two traits ---------------------------------------------------------------------
const GENS2 = SMOKE ? 1000 : 4000
const CENSUS2 = collect((SMOKE ? 500 : 2000):10:GENS2-10)
const PHI = [0.0, 60.0, 120.0, 180.0]
const base2 = Params(Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                     MutationRatePerGamete = 0.02, PopSize = 100, NoGenerations = GENS2, Sample = 10, InitFitnessSD = 5.0,
                     NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)
cycle2(p) = foldl(⊗, [withparams(LC.generation, p; θ = 2.0 .* [cosd(φ), sind(φ)]) for φ in PHI]) ⋅ LC.migrate(4)

function kernel2(p, seed)
    st = p ⊗ demes(4)
    idx = IdDict(w => i for (i, w) in enumerate(st.demes))
    gw = [zeros(2, 2) for _ in 1:4, _ in CENSUS2]
    rec = operator(:record, function (w, rng)
        j = searchsortedfirst(CENSUS2, w.gen)
        if j <= length(CENSUS2) && CENSUS2[j] == w.gen
            gw[idx[w], j] = cov(traitvalues(w; what = :genotype); dims = 2)
        end
        w
    end)
    o = run_composed(st, Xoshiro(seed); cycle = cycle2(p), record = false, images = Dict(:record => rec),
                     migration_rate = 0.01)
    any(o.extinct_gen .>= 0) && error("a deme died (seed $seed)")
    mean(gw)
end

function moments2(p, loci)
    mo = run_moments(p; cycle = cycle2(p), migration_rate = 0.01, between = true, loci)
    mean(mo.metapop.within[findall(in(CENSUS2), mo.metapop.gen)])
end

# --- identity cells -------------------------------------------------------------------
function identity_cells()
    say("\nI: identity of :rare_local and :rare_alleles (moments only)")
    ok = true
    p1 = Params(base1; NoGenerations = 1500)
    p2 = Params(base2; NoGenerations = 1500)
    function same(tag, f)
        a, b = f(:rare_alleles), f(:rare_local)
        eq = a == b
        eq || (ok = false; fail!("I $tag"))
        say(@sprintf("  %-46s %s", tag, eq ? "identical" : "DIFFERENT"))
    end
    for μ in (2e-4, 0.02)
        q1, q2 = Params(p1; MutationRatePerGamete = μ), Params(p2; MutationRatePerGamete = μ)
        same("one deme, one trait, μ=$μ, erosion", l -> run_moments(q1; erosion = true, loci = l).geno_var)
        same("one deme, one trait, μ=$μ, between", l -> run_moments(q1; between = true, loci = l).metapop.within)
        same("one deme, two traits, μ=$μ, erosion", l -> run_moments(q2; erosion = true, loci = l).gmatrix)
        same("one deme, two traits, μ=$μ, between", l -> run_moments(q2; between = true, loci = l).metapop.within)
        same("4 demes, local optima, no migration, μ=$μ", l -> run_moments(q1; cycle = cycle1(q1, LOCAL), migration_rate = 0.0,
                                                                            between = true, loci = l).metapop.within)
        same("4 demes, common optimum, m=0.01, μ=$μ", l -> run_moments(q1; cycle = cycle1(q1, zeros(4)), migration_rate = 0.01,
                                                                        between = true, loci = l).metapop.within)
        same("4 demes, two traits, no migration, μ=$μ", l -> run_moments(q2; cycle = cycle2(q2), migration_rate = 0.0,
                                                                          between = true, loci = l).metapop.within)
    end
    ok
end

function main()
    say("QTLSim validate/local_adaptation.jl")
    say("threads ", Threads.nthreads(), ", replicates per one-trait cell ", NREP, ", two-trait ", NREP2,
        ", SEED_OFFSET ", SEED_OFFSET, SMOKE ? " (SMOKE)" : "")
    t0 = time()
    identity_cells()
    summary = Tuple{String,Float64,Float64,Float64}[]
    verdict = Dict{String,Bool}()        # per cell: :rare_local passes A, R1 (and R2)

    # --- B and G: one trait ---------------------------------------------------------------
    say("\nB, G: one trait, k = 4 demes of N = 100, 50 loci, island migration m; census $(first(CENSUS1))–$(last(CENSUS1))")
    jobs = [(i, r) for i in eachindex(CELLS1) for r in 1:NREP]
    res = Vector{Float64}(undef, length(jobs))
    Threads.@threads :greedy for j in eachindex(jobs)
        i, r = jobs[j]
        c = CELLS1[i]
        p = Params(base1; MutationRatePerGamete = c.μ)
        res[j] = kernel1(p, c.θs, c.m, SEED_OFFSET + 5_000_000 + 1000i + r)
    end
    for (i, c) in enumerate(CELLS1)
        p = Params(base1; MutationRatePerGamete = c.μ)
        Ks = [res[j] for j in eachindex(jobs) if jobs[j][1] == i]
        K, se = mean(Ks), std(Ks) / sqrt(length(Ks))
        RA, RL, ML = (moments1(p, c.θs, c.m, l) for l in LOCI)
        eRA, eRL, eML = abs(K - RA) / K, abs(K - RL) / K, abs(K - ML) / K
        fA = cifail(K, se, RL, RL, 0.15)
        f1 = eRL > eML + 0.03
        f2 = c.hoc && !(eRL < eML)
        fA && fail!("A $(c.tag)"); f1 && fail!("R1 $(c.tag)"); f2 && fail!("R2 $(c.tag)")
        bA = cifail(K, se, RA, RA, 0.15) || eRA > eML + 0.03
        bA && push!(baseline, c.tag)
        verdict[c.tag] = !(fA || f1 || f2)
        say(@sprintf("\n  %s μ=%g m=%g optima %s", c.tag, c.μ, c.m, c.θs == LOCAL ? "(2, 1, −1, −2)" : "all 0"))
        say(@sprintf("    kernel V_W %.4g ± %.2g   RA %.4g (%+.1f%%)   RL %.4g (%+.1f%%)   ML %.4g (%+.1f%%)",
                     K, se, RA, 100(RA - K) / K, RL, 100(RL - K) / K, ML, 100(ML - K) / K))
        say(@sprintf("    RL: A %s, R1 %s%s;   baseline RA %s", fA ? "FAIL" : "ok", f1 ? "FAIL" : "ok",
                     c.hoc ? ", R2 " * (f2 ? "FAIL" : "ok") : "", bA ? "FAILS" : "ok"))
        push!(summary, (c.tag, eRA, eRL, eML))
    end

    # --- D1, D2: two traits ---------------------------------------------------------------
    say("\nD: two traits, k = 4 demes of N = 100, 50 loci, local optima, island migration 0.01; census $(first(CENSUS2))–$(last(CENSUS2))")
    for (tag, μ) in (("D1", 2e-4), ("D2", 0.02))
        p = Params(base2; MutationRatePerGamete = μ)
        Ks = Vector{Matrix{Float64}}(undef, NREP2)
        Threads.@threads :greedy for r in 1:NREP2
            Ks[r] = kernel2(p, SEED_OFFSET + 6_000_000 + (tag == "D1" ? 0 : 10_000) + r)
        end
        K = mean(Ks)
        SE = [std(getindex.(Ks, a, b)) / sqrt(NREP2) for a in 1:2, b in 1:2]
        RA, RL, ML = (moments2(p, l) for l in LOCI)
        say(@sprintf("\n  %s μ = %g", tag, μ))
        okA, bA = true, false
        for (a, b) in ((1, 1), (1, 2), (2, 2))
            f = cifail(K[a, b], SE[a, b], RL[a, b], sqrt(RL[a, a] * RL[b, b]), 0.15)
            f && (okA = false; fail!("A $tag G$a$b"))
            cifail(K[a, b], SE[a, b], RA[a, b], sqrt(RA[a, a] * RA[b, b]), 0.15) && (bA = true)
            sc(G) = sqrt(G[a, a] * G[b, b])
            say(@sprintf("    G%d%d  kernel %.4g ± %.2g   RA %.4g (%+.1f%%)   RL %.4g (%+.1f%%)   ML %.4g (%+.1f%%)  RL %s",
                         a, b, K[a, b], SE[a, b], RA[a, b], 100(RA[a, b] - K[a, b]) / sc(RA), RL[a, b],
                         100(RL[a, b] - K[a, b]) / sc(RL), ML[a, b], 100(ML[a, b] - K[a, b]) / sc(ML), f ? "FAIL" : "ok"))
        end
        eRA, eRL, eML = norm(K - RA) / norm(K), norm(K - RL) / norm(K), norm(K - ML) / norm(K)
        f1 = eRL > eML + 0.03
        f2 = tag == "D1" && !(eRL < eML)
        f1 && fail!("R1 $tag"); f2 && fail!("R2 $tag")
        (bA || eRA > eML + 0.03) && (bA = true; push!(baseline, tag))
        verdict[tag] = okA && !(f1 || f2)
        say(@sprintf("    e_RA = %.3f, e_RL = %.3f, e_ML = %.3f   R1 %s%s;   baseline RA %s", eRA, eRL, eML,
                     f1 ? "FAIL" : "ok", tag == "D1" ? ", R2 " * (f2 ? "FAIL" : "ok") : "", bA ? "FAILS" : "ok"))
        push!(summary, (tag, eRA, eRL, eML))
    end

    say("\nSummary (relative error of the prediction against the kernel):")
    say(@sprintf("  %-6s %8s %8s %8s", "cell", "e_RA", "e_RL", "e_ML"))
    for (c, a, b, d) in summary
        say(@sprintf("  %-6s %8.3f %8.3f %8.3f", c, a, b, d))
    end
    say(@sprintf("\nkernel runs and predictions in %.0f s", time() - t0))
    say("baseline (current :rare_alleles) fails: ", isempty(baseline) ? "none" : join(baseline, ", "))
    one = all(verdict[c.tag] for c in CELLS1)
    two = verdict["D1"] && verdict["D2"]
    localfail = [c.tag for c in CELLS1 if c.tag in baseline && c.θs == LOCAL && !verdict[c.tag]]
    say("decision rule, one trait: ", one ? ":auto → :rare_local for metapopulations" :
        isempty(localfail) ? "no change" : ":auto → :multilocus for metapopulations with differing optima (baseline fails, not fixed: " *
        join(localfail, ", ") * ")")
    say("decision rule, several traits: ", two ? ":auto → :rare_local" : "no change (:multilocus)")
    say(isempty(failures) ? "PASS: every pre-specified criterion holds" : "FAIL: " * join(failures, "; "))
    write(LOG, String(take!(out)))
    isempty(failures) || exit(1)
end

main()
