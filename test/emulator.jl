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
