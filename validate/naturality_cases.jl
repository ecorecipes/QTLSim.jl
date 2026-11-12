# The table of naturality cases (roadmap C2): validation as data. Each row is a
# `NaturalityCase` — a life cycle, its parameters, two backends, the statistics with their
# tolerances, the criterion, the replicates and the horizon — ported from one cell of a
# validation script, with that script's tolerance and criterion and its recorded verdict on
# the cell (`expected`). `validate/naturality.jl` runs them with `check_naturality` and
# reports whether each verdict was reproduced. The seeds are not the scripts' (the harness
# draws replicate r from `Xoshiro(seed + r)`, the scripts from MersenneTwister or Xoshiro
# with their own offsets), so a reproduced verdict is an independent replication of the
# cell, not a rerun of it.
#
# A new step, model or closure gets its validation by adding a row here. The scripts and
# their logs stay the record; these rows do not replace them.
#
# Include after `using QTLSim; import Catlab; using Catlab: ⋅, ⊗` (validate/naturality.jl).

const NL = lifecycle()
const NG = NL.generation

# validate/moments_loci.jl and validate/moments_regimes.jl: diploid, V_P = 1, h² = 0.5,
# m = 0.01, RecombPoissonMean = 2, Sample = 10
const LOCI_BASE = Params(Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01,
                         RecombPoissonMean = 2.0, Sample = 10)
regimes_params(N, th, L, Vs) = Params(LOCI_BASE; NoQTL = L, PopSize = N, MutationRatePerGamete = th / (4N),
                                      InitFitnessSD = Vs == 50.5 ? -1.0 : sqrt(Vs - envVar(LOCI_BASE)))

# validate/moments_demes.jl: k demes of N = 100, 50 loci, island migration at forward rate 0.01
demes_base(N) = Params(Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                       MutationRatePerGamete = 0.02, PopSize = N, NoGenerations = 4000, Sample = 10, InitFitnessSD = 5.0)
demes_cycle(θs) = foldl(⊗, [withparams(NG; θ = t) for t in θs]) ⋅ NL.migrate(length(θs))
const CENSUS = 2000:10:3990

# validate/rare_alleles_traits.jl: two traits
const TRAITS_BASE = Params(PopSize = 1000, Ploidy = 2, NoQTL = 20, PhenotypeVar = 0.02, Heritability = 0.5,
                           MutationScalar = 0.02, MutationRatePerGamete = 1e-3, InitFitnessSD = 2.0,
                           RecombPoissonMean = 3.0, NoGenerations = 8000, Sample = 10)
tcell(L, μ) = Params(TRAITS_BASE; NoQTL = L, MutationRatePerGamete = μ, NoTraits = 2,
                     MutCov = muteVar(TRAITS_BASE) .* [1.0 0.5; 0.5 0.6], SelCorr = -0.4)
# its D cells (validate/local_adaptation.jl's too): 4 demes of 100, 50 loci, two traits, local optima
dcell(μ) = Params(PopSize = 100, Ploidy = 2, NoQTL = 50, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
                  MutationRatePerGamete = μ, InitFitnessSD = 5.0, NoGenerations = 4000, Sample = 10, NoTraits = 2,
                  MutCorr = 0.3, SelCorr = 0.2)
dcycle(p) = foldl(⊗, [withparams(NG, p; θ = 2.0 .* [cosd(φ), sind(φ)]) for φ in (0.0, 60.0, 120.0, 180.0)]) ⋅ NL.migrate(4)
G_entries(window; within = false) =
    [QTLSim.gmatrix_mean(1, 1, window; within) => 0.0, QTLSim.gmatrix_mean(1, 2, window; within, scale = :geometric) => 0.0,
     QTLSim.gmatrix_mean(2, 2, window; within) => 0.0]
with_tol(stats, tols) = [s => t for ((s, _), t) in zip(stats, tols)]

# the neutral metapopulation of test/runtests.jl ("between = true"): 4 demes of 30, island migration
const NEUTRAL = Params(Ploidy = 2, NoQTL = 20, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.002,
                       MutationRatePerGamete = 0.05, Drift = true, PopSize = 30, NoGenerations = 300, Sample = 1)
function neutral_theory(c)
    ts = 1:first(c.params).NoGenerations-1
    th = neutral_deme_variances(backward_migration(c.options.migration_rate, 30), 30, muteVar(first(c.params)); t = collect(ts),
                                origins = :exchange)
    (gen = collect(ts), within = [sum(x.within) / length(x.within) for x in th], among = [x.among for x in th])
end
const NEUTRAL_GENS = 1:299
neutral_stat(what) = QTLSim.Statistic(what === :within ? "V_W, every generation" : "var of deme means, every generation",
    out -> error("theory or moments only");
    moments = out -> (mp = out.metapop; ix = findall(in(NEUTRAL_GENS), mp.gen); what === :within ? mp.within[ix] : mp.var_among[ix]),
    theory = out -> what === :within ? out.within : out.among)

# validate/slim/compare_demography.jl: Beverton–Holt, 100 loci, the equilibrium scenario
const BH = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 5000, PhenotypeVar = 0.022,
                  Heritability = 0.32, MutationScalar = 0.004, MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0,
                  InitFitnessSD = 5.0, MoveFitnessSDGen = -1, MoveType = :continuous, Demography = :beverton_holt,
                  GrowthRate = 1.5)

const NATURALITY_CASES = NaturalityCase[
    NaturalityCase(name = "moments_loci G1: L = 10, N = 1000, :qtsim, μ = 0.01",
        source = "validate/moments_loci.jl, criterion G1 (log: PASS, ratio 0.960)",
        params = Params(LOCI_BASE; NoQTL = 10, PopSize = 1000, CrossoverModel = :qtsim, MutationRatePerGamete = 0.01,
                        NoGenerations = 20_000),
        statistics = [QTLSim.window_mean(:geno_var, 10_000:10:19_990) => 0.10],
        options = (moments = (loci = :multilocus, erosion = true),), replicates = 12, seed = 10_000),
    NaturalityCase(name = "moments_loci, house of cards: L = 10, N = 1000, :qtsim, μ = 1e-4",
        source = "validate/moments_loci.jl, expectation (not a criterion): kernel < 85% of the multilocus prediction " *
                 "(log: ratio 0.094, confirmed); the G1 tolerance applied to it fails",
        params = Params(LOCI_BASE; NoQTL = 10, PopSize = 1000, CrossoverModel = :qtsim, MutationRatePerGamete = 1e-4,
                        NoGenerations = 20_000),
        statistics = [QTLSim.window_mean(:geno_var, 10_000:10:19_990) => 0.10],
        options = (moments = (loci = :multilocus, erosion = true),), replicates = 12, expected = :fail, seed = 20_000),
    NaturalityCase(name = "moments_regimes R1: N = 100, 4Nμ = 1, L = 10, V_s = 50.5",
        source = "validate/moments_regimes.jl, criterion R1 (log: PASS, ratio 1.015)",
        params = Params(regimes_params(100, 1.0, 10, 50.5); NoGenerations = 30_000),
        statistics = [QTLSim.window_mean(:geno_var, 10_000:10:29_990) => 0.15],
        options = (moments = (loci = :rare_alleles, erosion = true),), replicates = 16, seed = 30_000),
    NaturalityCase(name = "moments_regimes R1: N = 1000, 4Nμ = 0.1, L = 1000, V_s = 5.5",
        source = "validate/moments_regimes.jl, criterion R1 (log: FAIL, ratio 0.721; the same in the fresh-seed run)",
        params = Params(regimes_params(1000, 0.1, 1000, 5.5); NoGenerations = 30_000),
        statistics = [QTLSim.window_mean(:geno_var, 10_000:10:29_990) => 0.15],
        options = (moments = (loci = :rare_alleles, erosion = true),), replicates = 16, expected = :fail, seed = 40_000),
    NaturalityCase(name = "moments_demes C: local optima, k = 5",
        source = "validate/moments_demes.jl, scenario C (log: PASS; −3.5%, +6.3%, +7.1%; 300 replicates)",
        cycle = demes_cycle(collect(range(-2.0, 2.0, 5))), params = fill(demes_base(100), 5),
        statistics = [QTLSim.metapop_mean(:within, CENSUS) => 0.10, QTLSim.metapop_mean(:among, CENSUS) => 0.15,
                      QTLSim.qst_stat(CENSUS, 5, 100) => 0.15],
        options = (migration_rate = 0.01, moments = (between = true, loci = :multilocus)), replicates = 100, seed = 50_000),
    NaturalityCase(name = "moments_demes A: neutral, k = 5",
        source = "validate/moments_demes.jl, scenario A (log: PASS; z = −1.36, −0.30, 0.91; 300 replicates)",
        cycle = demes_cycle(zeros(5)), params = fill(Params(demes_base(100); Drift = true), 5), criterion = :z,
        statistics = [QTLSim.metapop_mean(:within, CENSUS) => 4.0, QTLSim.metapop_mean(:among, CENSUS) => 4.0,
                      QTLSim.qst_stat(CENSUS, 5, 100) => 4.0],
        options = (migration_rate = 0.01, moments = (between = true, loci = :multilocus)), replicates = 100, seed = 60_000),
    NaturalityCase(name = "neutral metapopulation: moments = neutral_deme_variances, every generation",
        source = "test/runtests.jl, \"between = true\" (exact to 1e-10); Lean Differentiation.lean",
        cycle = foldl(⊗, fill(NG, 4)) ⋅ NL.migrate(4), params = fill(NEUTRAL, 4), from = :moments, to = :theory,
        criterion = :exact, statistics = [neutral_stat(:within) => 1e-10, neutral_stat(:among) => 1e-10],
        options = (migration_rate = island_migration(4, 0.02), moments = (between = true, origins = :exchange),
                   theory = neutral_theory), replicates = 1),
    NaturalityCase(name = "rare_alleles_traits T: two traits, L = 20, μ = 1e-3",
        source = "validate/rare_alleles_traits.jl, criterion T-A (log: ok; −9.6%, −5.9%, −5.7%)",
        params = tcell(20, 1e-3), statistics = with_tol(G_entries(3000:10:7990), (0.20, 0.10, 0.20)),
        options = (moments = (loci = :rare_alleles, erosion = true),), replicates = 24, seed = 70_000),
    NaturalityCase(name = "rare_alleles_traits D1: two traits, 4 demes, local optima, :rare_alleles",
        source = "validate/rare_alleles_traits.jl, criterion D-A (log: FAIL on G11 and G22; the same in the fresh-seed run)",
        cycle = dcycle(dcell(2e-4)), params = fill(dcell(2e-4), 4),
        statistics = with_tol(G_entries(CENSUS; within = true), (0.15, 0.15, 0.15)),
        options = (migration_rate = 0.01, moments = (loci = :rare_alleles, between = true)), replicates = 48,
        expected = :fail, seed = 80_000),
    NaturalityCase(name = "local_adaptation D1: two traits, 4 demes, local optima, :rare_local",
        source = "validate/local_adaptation.jl, criterion A on D1 (log: ok; +0.8%, −4.6%, −9.1%)",
        cycle = dcycle(dcell(2e-4)), params = fill(dcell(2e-4), 4),
        statistics = with_tol(G_entries(CENSUS; within = true), (0.15, 0.15, 0.15)),
        options = (migration_rate = 0.01, moments = (loci = :rare_local, between = true)), replicates = 36, seed = 90_000),
    NaturalityCase(name = "compare_demography: Beverton–Holt equilibrium, kernel vs SLiM",
        source = "validate/slim/compare_demography.jl, equilibrium scenario (log: KS p = 0.913, 0.913; no discrepancy)",
        params = BH, from = :kernel, to = :slim, criterion = :ks,
        statistics = [QTLSim.window_mean(:popsize, 4000:4999) => 0.001, QTLSim.window_mean(:geno_var, 4000:4999) => 0.001],
        replicates = 40, seed = 100_000),
]
