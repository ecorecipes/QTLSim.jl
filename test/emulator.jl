# Tests: emulator. Parameters as an object (QTLSim.Para, sensitivities) and the emulator of
# the published sweep (emulate, emulated_half_survival_rate, emulate_sweep). Included by
# runtests.jl; run alone with `julia --project=vignettes test/runtests.jl emulator`.

@testset "Para and sensitivities" begin
    base = Params(PhenotypeVar = 0.8, Heritability = 0.3, MutationScalar = 0.004)
    # y = V_P²·h²: elasticities 2 and 1 exactly, 0 for a parameter it ignores
    P = QTLSim.Para(p -> p.PhenotypeVar^2 * p.Heritability, base, (:PhenotypeVar, :Heritability, :MutationScalar))
    @test P() ≈ 0.8^2 * 0.3
    @test P([2.0, 0.5, 0.004]) ≈ 4 * 0.5
    @test QTLSim.para_params(P, [2.0, 0.5, 0.01]).MutationScalar == 0.01
    e = sensitivities(P)
    @test keys(e) == (:PhenotypeVar, :Heritability, :MutationScalar)
    @test e.PhenotypeVar ≈ 2 atol = 1e-9
    @test e.Heritability ≈ 1 atol = 1e-9
    @test abs(e.MutationScalar) < 1e-12
    # plain derivatives: ∂y/∂V_P = 2V_P h²
    d = sensitivities(P; log = false)
    @test d.PhenotypeVar ≈ 2 * 0.8 * 0.3 rtol = 1e-9
    @test d.Heritability ≈ 0.8^2 rtol = 1e-9
    # against an analytic derivative in theory.jl: the exact Gaussian (gaussian_vg_exact) V_G = (V_m + √(V_m² + 4V_mV_s))/2
    # with V_m = m·V_E, V_s = ω² + V_E; its elasticity in m, by hand
    p = Params(PhenotypeVar = 1.0, Heritability = 0.4, MutationScalar = 0.01, InitFitnessSD = 3.0)
    G = QTLSim.Para(gaussian_vg_exact, p, (:MutationScalar, :InitFitnessSD))
    Vm, Vs = QTLSim.muteVar(p), selection_variance(p)
    r = sqrt(Vm^2 + 4Vm * Vs)
    dVdVm = (1 + (Vm + 2Vs) / r) / 2                     # ∂V_G/∂V_m
    dVdω = (Vm / r) * 2 * 3.0                             # ∂V_G/∂V_s · ∂V_s/∂ω
    @test gaussian_vg_exact(p) ≈ (Vm + r) / 2
    s = sensitivities(G)
    @test s.MutationScalar ≈ dVdVm * Vm / gaussian_vg_exact(p) rtol = 1e-8
    @test s.InitFitnessSD ≈ dVdω * 3.0 / gaussian_vg_exact(p) rtol = 1e-8
    # a larger step is still accurate after Richardson extrapolation (error O(h⁴))
    @test sensitivities(G; rel = 0.05).MutationScalar ≈ s.MutationScalar rtol = 1e-5
    # integer, Symbol and Bool fields have no derivative; unknown names, repeats
    @test_throws ArgumentError QTLSim.Para(gaussian_vg, p, (:NoQTL,))
    @test_throws ArgumentError QTLSim.Para(gaussian_vg, p, (:MoveType,))
    @test_throws ArgumentError QTLSim.Para(gaussian_vg, p, (:Drift,))
    @test_throws ArgumentError QTLSim.Para(gaussian_vg, p, (:NoSuchField,))
    @test_throws ArgumentError QTLSim.Para(gaussian_vg, p, (:Heritability, :Heritability))
    @test_throws ArgumentError QTLSim.Para(gaussian_vg, p, ())
    @test_throws ArgumentError P([1.0])
    # an elasticity needs a positive parameter (InitFitnessSD = -1 is the analytical width)
    @test_throws ArgumentError sensitivities(QTLSim.Para(gaussian_vg, Params(), (:InitFitnessSD,)))
    @test sensitivities(QTLSim.Para(p -> -1.0, p, (:Heritability,))).Heritability |> isnan
end

@testset "the emulator of the published sweep" begin
    # a small published-protocol input: 20 loci, a short burn-in and horizon
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 20, RecombPoissonMean = 3.0, PhenotypeVar = 0.03,
               Heritability = 0.4, MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 3.0,
               MoveType = :continuous, MoveFitnessSD = 1e-3, MoveFitnessSDGen = 10000, NoGenerations = 250000)
    kw = (burnin = 300, horizon = 1500, tail = 20_000)
    # the half-survival rate first: P2000 = 1/2 there, monotone around it
    h = emulated_half_survival_rate(p; kw.burnin, kw.horizon)
    @test 0 < h.k50 < Inf && h.kc ≈ drift_critical_rate(p).k
    at(k) = emulate(Params(p; MoveFitnessSD = k); kw...)
    @test at(h.k50).P2000 ≈ 0.5 atol = 1e-6
    @test at(0.99h.k50).P2000 > 0.5 > at(1.01h.k50).P2000
    # a narrow bracket gives the same root with fewer evaluations; one that misses is moved
    h2 = emulated_half_survival_rate(p; kw.burnin, kw.horizon, bracket = (0.995h.k50, 1.005h.k50))
    @test h2.k50 ≈ h.k50 rtol = 1e-8
    @test h2.evals < h.evals && isnan(h2.kc)
    @test emulated_half_survival_rate(p; kw.burnin, kw.horizon, bracket = (1.01h.k50, 1.02h.k50)).k50 ≈ h.k50 rtol = 1e-8
    @test_throws ArgumentError emulated_half_survival_rate(p; bracket = (2.0, 1.0))
    # below it: alive, with a hazard; the pipeline by hand
    k = 0.9h.k50
    r = at(k)
    @test 0 <= r.P240k <= r.P2000 <= 1 && 0 <= r.hazard <= 1 && !r.extinct
    @test r.VG_move > 0 && 0 < r.Wbar_move <= 1 && r.lag_end > 0
    q = Params(p; MoveFitnessSD = k, Sample = 1, NoGenerations = 1800, MoveFitnessSDGen = 300)
    m = run_moments(q; VG = house_of_gauss_vg(p), erosion = true, loci = :rare_alleles)
    d = extinction_distribution(q, m; from = 300)
    @test r.P2000 == d.psurvive
    @test r.VG_move == m.geno_var[301]
    @test r.P240k ≈ r.P2000 * (1 - r.hazard)^(20_000 - 1500)
    # above it: faster moves die earlier; T is counted from 1 at the first generation of the move
    f1, f2 = at(2h.k50), at(10h.k50)
    @test f1.P2000 < 0.5 && f2.P2000 < 0.5
    @test 1 <= f2.T < f1.T <= 1500
    d2 = extinction_distribution(Params(q; MoveFitnessSD = 2h.k50), run_moments(Params(q; MoveFitnessSD = 2h.k50);
                                 VG = house_of_gauss_vg(p), erosion = true, loci = :rare_alleles); from = 300)
    @test f1.T == d2.median - 299
    @test at(0.0).P2000 > 0.99
    @test_throws ArgumentError emulate(Params(p; MoveType = :phenotype))
    @test_throws ArgumentError emulate(p; horizon = 1)
    # the batch helper: the same numbers, cached and resumed
    ps = [Params(p; MoveFitnessSD = x) for x in (0.5h.k50, 2h.k50, 10h.k50)]
    rows = emulate_sweep(ps; kw...)
    @test [x.P2000 for x in rows] == [emulate(q; kw...).P2000 for q in ps]
    @test isnan(rows[1].T) && rows[2].T == f1.T
    mktempdir() do dir
        f = joinpath(dir, "emu.csv")
        emulate_sweep(ps[1:2]; ids = ["a", "b"], cache = f, kw...)
        b = emulate_sweep(ps; ids = ["a", "b", "c"], cache = f, kw...)
        @test countlines(f) == 4
        @test [x.P2000 for x in b] == [x.P2000 for x in rows]
        @test isequal([x.T for x in b], [x.T for x in rows])
    end
    # loci = :spectrum and its options reach the moments lowering: the default unchanged bit for
    # bit, the dynamic LD (spectrum = (; ld = :dynamic)) is the pipeline by hand with that option
    ks = Params(p; MoveFitnessSD = 2h.k50)
    qs = Params(q; MoveFitnessSD = 2h.k50)
    es = emulate(ks; loci = :spectrum, kw...)
    @test isequal(es, emulate(ks; loci = :spectrum, spectrum = (ld = :none,), kw...))
    el = emulate(ks; loci = :spectrum, spectrum = (ld = :dynamic,), kw...)
    ml = run_moments(qs; VG = house_of_gauss_vg(p), erosion = true, loci = :spectrum, spectrum = (ld = :dynamic,))
    @test el.P2000 == extinction_distribution(qs, ml; from = 300).psurvive
    @test el.VG_move == ml.geno_var[301] && el.VG_move < es.VG_move          # negative LD lowers V_G at the move
    @test emulated_half_survival_rate(p; loci = :spectrum, spectrum = (ld = :dynamic,), kw.burnin, kw.horizon,
                                      bracket = (0.5h.k50, 2h.k50)).k50 > 0
    # the beta closure (spectrum = (; solver = :beta)): the grid's result unchanged, `fallback` added
    # only with the option, the pipeline by hand, deterministic; the half-survival rate counts fallbacks
    @test isequal(es, emulate(ks; loci = :spectrum, spectrum = (solver = :grid,), kw...)) && !haskey(es, :fallback)
    eb = emulate(ks; loci = :spectrum, spectrum = (solver = :beta,), kw...)
    mb = run_moments(qs; VG = house_of_gauss_vg(p), erosion = true, loci = :spectrum, spectrum = (solver = :beta,))
    @test eb.P2000 == extinction_distribution(qs, mb; from = 300).psurvive && eb.fallback == mb.spectrum_fallback
    @test isequal(eb, emulate(ks; loci = :spectrum, spectrum = (solver = :beta,), kw...))
    @test eb.VG_move ≈ es.VG_move rtol = 0.01
    @test (eb.P2000 >= 0.5) == (es.P2000 >= 0.5)
    hb = emulated_half_survival_rate(p; loci = :spectrum, spectrum = (solver = :beta,), kw.burnin, kw.horizon,
                                     bracket = (0.5h.k50, 2h.k50))
    @test hb.k50 > 0 && 0 <= hb.fallbacks <= hb.evals
    @test haskey(emulate_sweep([ks]; loci = :spectrum, spectrum = (solver = :beta,), kw...)[1], :fallback)
end

@testset "memoised pair classes and recombination" begin
    p = Params(NoQTL = 300, Ploidy = 2, RecombPoissonMean = 3.0)
    ω, ρ = QTLSim.pair_classes(p)
    @test (ω, ρ) == QTLSim._pair_classes(p)
    ω[1] = -1.0                                            # a copy: the cache is untouched
    @test QTLSim.pair_classes(p)[1][1] == 1 / 300
    @test QTLSim.pair_classes(Params(p; CrossoverModel = :all)) == QTLSim._pair_classes(Params(p; CrossoverModel = :all))
    @test QTLSim.pair_recombination(p) == QTLSim._pair_recombination(300, 3.0)
    @test QTLSim.pair_recombination(300, 3.0; model = :all) == QTLSim._pair_recombination(300, 3.0; model = :all)
end

@testset "the emulator's between-run variance" begin
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 20, RecombPoissonMean = 3.0, PhenotypeVar = 0.03,
               Heritability = 0.4, MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 3.0,
               MoveType = :continuous, MoveFitnessSD = 1e-3, MoveFitnessSDGen = 10000, NoGenerations = 250000)
    kw = (burnin = 300, horizon = 1500, tail = 20_000)
    q = Params(p; Sample = 1, NoGenerations = 1800, MoveFitnessSDGen = 300)
    VG = house_of_gauss_vg(p)
    # replicates = 0 is the deterministic emulator, bitwise
    @test emulate(p; kw..., replicates = 0) === emulate(p; kw...)
    m = run_moments(q; VG, erosion = true, loci = :rare_alleles)
    @test emulate(p; kw..., replicates = 0).P2000 == extinction_distribution(q, m; from = 300).psurvive
    # with the noise scaled to zero every replicate is the deterministic trajectory
    vg = QTLSim.vg_fluctuation(q)
    r0 = QTLSim.moments_replicates(q; replicates = 3, rng = Xoshiro(1), from = 300, VG, loci = :rare_alleles, scale = 0, vg)
    @test r0.det.gen == 0:299
    for rep in r0.replicates
        @test [r0.det.mean_fitness; rep.mean_fitness] == m.mean_fitness
        @test [r0.det.pheno_mean; rep.pheno_mean] == m.pheno_mean
        @test [r0.det.geno_var; rep.geno_var] == m.geno_var
    end
    # a fixed seed (or the same generator) gives the same result; another seed another one
    a = emulate(p; kw..., replicates = 8, seed = 3)
    @test isequal(a, emulate(p; kw..., replicates = 8, seed = 3))
    @test isequal(a, emulate(p; kw..., replicates = 8, rng = Xoshiro(3)))
    h = emulated_half_survival_rate(p; kw.burnin, kw.horizon)
    near = Params(p; MoveFitnessSD = 1.05h.k50)
    b = emulate(near; kw..., replicates = 16, seed = 3)
    @test b.P2000 != emulate(near; kw..., replicates = 16, seed = 4).P2000
    @test 0 <= b.P240k <= b.P2000 <= 1 && 0 <= b.hazard <= 1 && b.VG_move > 0
    # the noise turns the step into a graded probability around the deterministic k50
    @test emulate(Params(p; MoveFitnessSD = 0.9h.k50); kw...).P2000 > 0.999
    @test emulate(Params(p; MoveFitnessSD = 1.1h.k50); kw...).P2000 < 0.001
    pn = [emulate(Params(p; MoveFitnessSD = f * h.k50); kw..., replicates = 32, seed = 5).P2000 for f in (0.9, 1.1)]
    @test 0.05 < pn[2] < pn[1] < 0.995
    # the mean-only noise (vgnoise = false) also runs, and an rng or a seed is required
    @test 0 <= emulate(near; kw..., replicates = 4, seed = 1, vgnoise = false).P2000 <= 1
    @test_throws ArgumentError emulate(p; kw..., replicates = 2)
    @test_throws ArgumentError emulate(p; kw..., replicates = -1)
    # emulate_sweep: per-input seeds from the base seed and the id, whatever the threads and the cache
    ps = [Params(p; MoveFitnessSD = x * h.k50) for x in (0.95, 1.05, 1.2)]
    rows = emulate_sweep(ps; ids = ["a", "b", "c"], kw..., replicates = 4, seed = 11)
    @test [x.P2000 for x in rows] == [emulate(ps[i]; kw..., replicates = 4, seed = QTLSim.job_seed(11, id)).P2000
                                      for (i, id) in enumerate(["a", "b", "c"])]
    @test isequal(rows, emulate_sweep(ps; ids = ["a", "b", "c"], kw..., replicates = 4, seed = 11, ntasks = 1))
    mktempdir() do dir
        f = joinpath(dir, "emu.csv")
        emulate_sweep(ps[[3]]; ids = ["c"], cache = f, kw..., replicates = 4, seed = 11)
        @test [x.P2000 for x in emulate_sweep(ps; ids = ["a", "b", "c"], cache = f, kw..., replicates = 4, seed = 11)] ==
              [x.P2000 for x in rows]
    end
    @test_throws ArgumentError emulate_sweep(ps; kw..., replicates = 4)
    @test QTLSim.job_seed(11, "a") === QTLSim.job_seed(11, "a") != QTLSim.job_seed(12, "a")
end

@testset "drift noise of the moments replicates" begin
    # neutral: from a start with no spread, the replicates' variance of the mean genotype after
    # t generations is the sum of the innovations, Σ_s (V_G,s + V_m)/N to leading order
    pn = Params(PopSize = 500, Ploidy = 2, NoQTL = 10, RecombPoissonMean = 3.0, PhenotypeVar = 1.0,
                Heritability = 0.4, MutationScalar = 0.002, MutationRatePerGamete = 1e-3, Drift = true,
                Sample = 1, NoGenerations = 41)
    rn = QTLSim.moments_replicates(pn; replicates = 4000, rng = Xoshiro(2), from = 0, VG = 0.4, loci = :multilocus)
    vgs = rn.replicates[1].geno_var .* (500 - 1) / 500         # the recorded V_G is the sample variance
    @test all(rep.geno_var == rn.replicates[1].geno_var for rep in rn.replicates)   # no V_G noise asked for
    t = 40
    z = [rep.pheno_mean[t+1] for rep in rn.replicates]
    expect = sum(vgs[1:t]) / 500 + t * QTLSim.muteVar(pn) / 500
    @test mean(z) ≈ 0 atol = 4sqrt(expect / 4000)
    @test var(z) ≈ expect rtol = 0.08                          # SE of a variance from 4000 draws ≈ 2.2%
    # a fixed optimum: selection shrinks a deviation by 1 − V_G/s, so the spread settles
    ps = Params(pn; Drift = false, InitFitnessSD = 1.0, NoGenerations = 400)
    rs = QTLSim.moments_replicates(ps; replicates = 2000, rng = Xoshiro(3), from = 200, VG = 0.4, loci = :multilocus)
    zs = [rep.pheno_mean[end] for rep in rs.replicates]
    VGs, s = rs.replicates[1].geno_var[end], 1.0 + rs.replicates[1].geno_var[end] + QTLSim.envVar(ps)
    # the stationary variance of the AR(1) mean, about s/(2n) (Lande 1976), n = N·W̄
    n = 500 * mean(rs.replicates[1].mean_fitness)
    @test var(zs) ≈ s / (n * (2 - VGs / s)) rtol = 0.15
    # V_G noise: the innovation accumulates in V_G, relaxed only slowly over 30 generations
    vg = (stationary = 0.0, innovation = 1e-4)
    rv = QTLSim.moments_replicates(Params(pn; NoGenerations = 31); replicates = 3000, rng = Xoshiro(4), from = 0, VG = 0.4,
                                   loci = :multilocus, vg)
    v30 = [rep.geno_var[end] * 499 / 500 for rep in rv.replicates]
    @test var(v30) ≈ 30 * 1e-4 rtol = 0.15
    # vg_fluctuation, neutral: stationary 4nLμα⁴ (cv² = 1/(4nLμ)), innovation 4Lμα⁴
    f = QTLSim.vg_fluctuation(pn)
    α2 = QTLSim.mutation_effect_variance(pn)
    @test f.stationary ≈ 4 * 500 * 10 * 1e-3 * α2^2 rtol = 1e-8
    @test f.innovation ≈ 4 * 10 * 1e-3 * α2^2 rtol = 1e-8
    @test f.cv ≈ sqrt(f.stationary) / f.VG
    # under selection both fall below the neutral values
    fs = QTLSim.vg_fluctuation(Params(pn; Drift = false, InitFitnessSD = 0.5))
    @test 0 < fs.innovation < f.innovation && 0 < fs.stationary < f.stationary
end
