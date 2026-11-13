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

# the neutral metapopulation of test/demes.jl ("between = true"): 4 demes of 30, island migration
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

# validate/compare_to_qtsim.jl: kernel vs QTSim, Welch t, |t| ≥ 4 fails; extinct runs left out
const QT_EQUILIBRIUM = Params(Sample = 10, PopSize = 1000, NoQTL = 1, Ploidy = 2, NoGenerations = 4000, PhenotypeVar = 0.02,
                              Heritability = 0.5, MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)
const QT_DIRECTIONAL = Params(Sample = 1, PopSize = 1000, NoQTL = 1, Ploidy = 2, NoGenerations = 6500, PhenotypeVar = 0.02,
                              Heritability = 0.5, MutationScalar = 0.001, ChangeFitnessGen = 5000, ChangeFitnessSD = 0.30,
                              MoveFitnessSDGen = 5000, MoveFitnessSD = 1.0, MoveType = :phenotype)
inw(out, lo, hi) = findall(g -> lo <= g <= hi, out.gen)
inwindow(out, w) = findall(in(w), out.gen)
# the advance of the mean phenotype after the move and the survivors' share while it moves, as the script defines them
const QT_ADVANCE = QTLSim.surviving(QTLSim.Statistic("phenotypic advance",
    out -> mean(out.pheno_mean[inw(out, 6400, 6500)]) - mean(out.pheno_mean[inw(out, 5000, 5010)])))
const QT_SURVIVAL = QTLSim.surviving(QTLSim.Statistic("survival fraction", out -> mean(out.nsurv[inw(out, 5000, 6500)]) / 1000))

# validate/slim/compare_fluctuating.jl: one locus, 1000 generations of stabilising selection, then 2000 of a
# fluctuating optimum; statistics per run over 1500-2999, on surviving runs; KS p < 0.001 fails
const FLUCT_BASE = Params(PopSize = 1000, Ploidy = 2, NoQTL = 1, NoGenerations = 3000, PhenotypeVar = 0.02,
                          Heritability = 0.5, MutationScalar = 0.02, MutationRatePerGamete = 0.5,
                          InitFitnessSD = 0.3, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.0)
const FLUCT_W = 1500:2999
fluct_stats(red) = map(QTLSim.surviving, vcat(
    [QTLSim.lag_mean(FLUCT_W; squared = true), QTLSim.window_mean(:geno_var, FLUCT_W), QTLSim.window_mean(:nsurv, FLUCT_W)],
    red ? [QTLSim.Statistic("SD of θ over $FLUCT_W", out -> std(out.theta[inwindow(out, FLUCT_W)])),
           QTLSim.Statistic("lag-1 autocorrelation of θ", out -> (θ = out.theta[inwindow(out, FLUCT_W)]; cor(θ[1:end-1], θ[2:end])))] :
          QTLSim.Statistic[]))

# test/fluctuating.jl, "fluctuating optima": the moments lowering of a periodic optimum against expected_lag2
const FLUCT_TEST = Params(PopSize = 1000, NoQTL = 1, Ploidy = 2, NoGenerations = 4000, PhenotypeVar = 0.02, Heritability = 0.5,
                          MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)
const FLUCT_PM = Params(FLUCT_TEST; MoveType = :periodic, MoveFitnessSDGen = 500, MoveFitnessSD = 0.0, MoveAmplitude = 0.1,
                        MovePeriod = 200.0)
const FLUCT_V = multilocus_equilibrium(FLUCT_TEST).VG
const FLUCT_LAG2 = QTLSim.Statistic("mean (θ − z̄)² over 2000:3999",
    out -> mean((out.theta[i] - out.pheno_mean[i])^2 for i in inwindow(out, 2000:3999)); theory = out -> out.lag2)

# validate/slim/compare_demography.jl: vignette 10's parameters, Beverton–Holt, a continuous move from 5000
const BH_MOVE = Params(BH; NoGenerations = 7000, MoveFitnessSDGen = 5000, MoveFitnessSD = 0.02, GrowthRate = 2.0)

# validate/slim/compare_plasticity.jl, scenario 1: b/B = 0.5, a perfect cue, the optimum moving 0.01 from 2000
const PLASTIC_1 = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 2500, PhenotypeVar = 0.022,
                         Heritability = 0.32, MutationScalar = 0.004, MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0,
                         InitFitnessSD = 5.0, MoveType = :continuous, MoveFitnessSDGen = 2000, MoveFitnessSD = 0.01,
                         Plasticity = 0.5)

# validate/reaction_norm.jl, scenario A: a heritable reaction norm, no cost; criterion C1, |kernel − moments| ≤ 0.1
# for the mean slope (record row T0 + t is generation T0 + t − 1)
const RN_BASE = Params(NoGenerations = 2500, PopSize = 1000, Ploidy = 2, NoQTL = 50, MutationRatePerGamete = 0.01,
                       MoveType = :continuous, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.02, RecombPoissonMean = 3.0,
                       InitFitnessSD = 1.0, PhenotypeVar = 1.0, Heritability = 0.5)
const RN_A = Params(RN_BASE; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2,
                    MutCov = [muteVar(RN_BASE) 0.0; 0.0 muteVar(RN_BASE) / 4], PlasticityCostSD = Inf)
rn_slope(t) = QTLSim.Statistic("mean slope at +$t", out -> out.trait_mean[only(findall(==(1000 + t - 1), out.gen))][2];
                               scale = _ -> 1.0)

# validate/fst.jl, section 2, island k = 4, Nm = 1: neutral demes of 50, 20 loci, the forward rate whose
# backward migrant fraction is 0.02, against neutral_deme_variances fed migrate!'s own law of pairs' origins
const FST_P0 = Params(Ploidy = 2, NoQTL = 20, RecombPoissonMean = 1.0, PhenotypeVar = 1.0, Heritability = 0.5,
                      MutationScalar = 0.002, MutationRatePerGamete = 0.05, Drift = true)
function fst_forward_rate(k, N, m)              # fst.jl's bisection
    lo, hi = 0.0, 0.99
    f(r) = 1 - backward_migration(island_migration(k, r), N)[1, 1] - m
    for _ in 1:60
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    (lo + hi) / 2
end
const FST_M = island_migration(4, fst_forward_rate(4, 50, 0.02))
const FST_CENSUS = 1250:10:2490
# the law of the deme origins of two distinct adults after one call of migrate!, by Monte Carlo of the
# migration step alone (fst.jl's `migration_check`, its estimator over all calls)
function fst_pair_law(M, k, N, ncalls, seed)
    s = Params(PopSize = N, NoQTL = 1, Ploidy = 2, NoGenerations = 1, Drift = true) ⊗ demes(k)
    rng = Xoshiro(seed)
    Q, nsum, n = zeros(k^2, k^2), zeros(k, k), zeros(k, k)
    for _ in 1:ncalls
        for (a, w) in enumerate(s.demes); fill!(view(w.env, 1:w.N), a); end
        migrate!(s, rng; rate = M)
        fill!(n, 0)
        for (a, w) in enumerate(s.demes), i in 1:w.N
            n[a, Int(w.env[i])] += 1
        end
        v = vec(n); Q .+= v .* v'; nsum .+= n
    end
    J = zeros(k^2, k^2)
    for a in 1:k, b in 1:k, c in 1:k, d in 1:k
        x = Q[a + (c - 1) * k, b + (d - 1) * k]
        J[a + (b - 1) * k, c + (d - 1) * k] = a == b ? (x - (c == d) * nsum[a, c]) / (ncalls * N * (N - 1)) : x / (ncalls * N^2)
    end
    J
end
function fst_theory(c)
    J = fst_pair_law(FST_M, 4, 50, 40_000, 7002)     # fst.jl's call for this scenario (its seed 7000 + 2)
    th = neutral_deme_variances(backward_migration(FST_M, 50), 50, muteVar(FST_P0); t = collect(FST_CENSUS), origins = J)
    (within = mean(mean(x.within) for x in th), among = mean(x.among for x in th))
end
const FST_W = QTLSim.Statistic("V_W (within demes) over $FST_CENSUS",
    QTLSim.metapop_mean(:within, FST_CENSUS).of[:kernel]; theory = out -> out.within)
const FST_A = QTLSim.Statistic("var of deme means over $FST_CENSUS",
    QTLSim.metapop_mean(:among, FST_CENSUS).of[:kernel]; theory = out -> out.among)
const FST_Q = QTLSim.Statistic("Q_ST over $FST_CENSUS", [FST_W, FST_A], (w, a) -> (σ2 = 3 / 4 * (a - w / 50); σ2 / (σ2 + 2w)))

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
        source = "test/demes.jl, \"between = true\" (exact to 1e-10); Lean Differentiation.lean",
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
    NaturalityCase(name = "compare_to_qtsim A: stabilising equilibrium, kernel vs QTSim",
        source = "validate/compare_to_qtsim.jl, scenario A (no committed log; CI's QTSim parity job passes it; its fixed " *
                 "seeds rerun here 2026-11-13: t = −0.22, consistent)",
        params = QT_EQUILIBRIUM, from = :kernel, to = :qtsim, criterion = :welch,
        statistics = [QTLSim.surviving(QTLSim.window_mean(:geno_var, 3000:3999)) => 4.0], replicates = 40, seed = 200_000),
    NaturalityCase(name = "compare_to_qtsim B: directional adaptation, kernel vs QTSim",
        source = "validate/compare_to_qtsim.jl, scenario B (no committed log; CI's QTSim parity job passes it; its fixed " *
                 "seeds rerun here 2026-11-13: t = −1.66, +0.27, −0.33; no extinctions)",
        params = QT_DIRECTIONAL, from = :kernel, to = :qtsim, criterion = :welch,
        statistics = [QTLSim.surviving(QTLSim.window_mean(:geno_var, 4000:4999)) => 4.0, QT_ADVANCE => 4.0, QT_SURVIVAL => 4.0],
        replicates = 24, seed = 210_000),
    NaturalityCase(name = "compare_fluctuating: periodic, A = 0.1, P = 200, kernel vs SLiM",
        source = "validate/slim/compare_fluctuating.jl, first scenario (valid_fluctuating.log: KS p = 0.029, 0.573, 0.759; " *
                 "no discrepancy)",
        params = Params(FLUCT_BASE; MoveType = :periodic, MoveAmplitude = 0.1, MovePeriod = 200.0), from = :kernel, to = :slim,
        criterion = :ks, statistics = [s => 0.001 for s in fluct_stats(false)], replicates = 40, seed = 220_000),
    NaturalityCase(name = "compare_fluctuating: red noise, σ = 0.1, ρ = 0.9, kernel vs SLiM",
        source = "validate/slim/compare_fluctuating.jl, third scenario (valid_fluctuating.log: KS p = 0.988, 0.759, 0.913, " *
                 "0.913, 0.913; no discrepancy)",
        params = Params(FLUCT_BASE; MoveType = :red_noise, MoveNoiseSD = 0.1, MoveAutocorr = 0.9), from = :kernel, to = :slim,
        criterion = :ks, statistics = [s => 0.001 for s in fluct_stats(true)], replicates = 40, seed = 230_000),
    NaturalityCase(name = "periodic optimum: moments = expected_lag2",
        source = "test/fluctuating.jl, \"fluctuating optima\" (the moments lowering's mean squared lag over 2000-3999 " *
                 "≈ expected_lag2 at the same V_G, rtol = 1e-3; the test passes)",
        params = FLUCT_PM, from = :moments, to = :theory, criterion = :exact, statistics = [FLUCT_LAG2 => 1e-3],
        options = (moments = (VG = FLUCT_V,), theory = c -> (lag2 = expected_lag2(FLUCT_PM; VG = FLUCT_V),)), replicates = 1),
    NaturalityCase(name = "compare_demography: Beverton–Holt R = 2, optimum moving 0.02, kernel vs SLiM",
        source = "validate/slim/compare_demography.jl, second scenario (valid_demography.log: 40/40 extinct in both; " *
                 "KS p = 0.263, 0.573; no discrepancy)",
        params = BH_MOVE, from = :kernel, to = :slim, criterion = :ks,
        statistics = [QTLSim.persistence(5000) => 0.001, QTLSim.window_mean(:popsize, 4000:4999) => 0.001],
        replicates = 20, seed = 240_000, note = "20 replicates per backend, the script's 40 (cost: ~25 CPU-s per SLiM run)"),
    NaturalityCase(name = "compare_plasticity 1: b/B = 0.5, optimum moving 0.01, kernel vs SLiM",
        source = "validate/slim/compare_plasticity.jl, scenario 1 (valid_plasticity.log: KS p = 0.759, 0.759, 0.988, 0.759; " *
                 "no discrepancy)",
        params = PLASTIC_1, from = :kernel, to = :slim, criterion = :ks,
        statistics = [QTLSim.window_mean(:geno_var, 1500:1999) => 0.001, QTLSim.lag_mean(2250:2499) => 0.001,
                      QTLSim.window_mean(:geno_var, 2250:2499) => 0.001, QTLSim.window_mean(:nsurv, 2250:2499) => 0.001],
        replicates = 20, seed = 250_000, note = "20 replicates per backend, the script's 40"),
    NaturalityCase(name = "reaction_norm A: heritable reaction norm, no cost, kernel vs moments (C1)",
        source = "validate/reaction_norm.jl, scenario A, criterion C1 at +250 and +500 (reaction_norm.log: PASS; slope " *
                 "kernel 0.781 vs moments 0.772, 0.892 vs 0.883)",
        params = RN_A, criterion = :diff, statistics = [rn_slope(250) => 0.1, rn_slope(500) => 0.1],
        options = (moments = (erosion = true, loci = :multilocus),), replicates = 40, seed = 260_000),
    NaturalityCase(name = "fst: neutral island k = 4, Nm = 1, kernel vs deme_coalescence",
        source = "validate/fst.jl, section 2, island k=4 Nm=1.0 (fst.log: PASS; z = −2.57, −2.11, −0.57)",
        cycle = foldl(⊗, fill(NG, 4)) ⋅ NL.migrate(4), params = fill(Params(FST_P0; PopSize = 50, NoGenerations = 2500,
                                                                               Sample = 10), 4),
        from = :kernel, to = :theory, criterion = :z, statistics = [FST_W => 4.0, FST_A => 4.0, FST_Q => 4.0],
        options = (migration_rate = FST_M, theory = fst_theory), replicates = 400, seed = 270_000,
        note = "z uses the simulation's SE only; fst.jl also folded in the prediction's Monte Carlo SE (for the variance " *
               "among deme means 0.00053 against the simulation's 0.0017), so this z is slightly larger"),
]
