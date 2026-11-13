# Tests: critical. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl critical`.

@testset "the critical rate under drift (drift_critical_rate, half_survival_rate)" begin
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0, PhenotypeVar = 0.03, Heritability = 0.4,
               MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous)
    r = drift_critical_rate(p)
    Vs = selection_variance(p)
    s = initFitSD(p)^2 + r.VG + QTLSim.envVar(p)
    @test r.k ≈ r.VG * r.λ / (r.VG + Vs)                                 # the response at the critical lag
    @test r.n ≈ p.PopSize * r.Wbar && r.N == p.PopSize
    @test r.M ≈ r.λ^2 / s
    rate(λ) = (st = QTLSim._feedback_state(p, λ, 0.0); st.VG * λ / (st.VG + Vs))
    @test all(rate(l) <= r.k * (1 + 1e-9) for l in range(0, 3r.λ, length = 301))       # the maximum
    @test r.VG < house_of_gauss_vg(p)                                    # fewer breeders, less variance
    @test QTLSim._house_of_gauss(p, 2effective_size(1000, 1000)) ≈ house_of_gauss_vg(p)
    # a cost: linear to first order with scale M (the envelope theorem), monotone
    h = 1e-3
    @test (drift_critical_rate(p; cost = h).k / r.k - 1) / h ≈ -1 / r.M rtol = 0.05
    @test drift_critical_rate(p; cost = 0.5).k < drift_critical_rate(p; cost = 0.25).k < r.k
    @test drift_critical_rate(Params(p; Plasticity = 0.5, PlasticityCostSD = 0.5)).k ≈ drift_critical_rate(p; cost = 0.5).k / 0.5
    @test drift_critical_rate(Params(p; Plasticity = 0.5)).k ≈ 2r.k
    # Beverton–Holt: the breeders are the equilibrium K(RW̄ − 1)/(R − 1), adults n/W̄; a lower R, a lower rate
    bh = Params(p; Demography = :beverton_holt, GrowthRate = 2.0)
    rb = drift_critical_rate(bh)
    @test rb.n ≈ 1000 * (2rb.Wbar - 1) && rb.N ≈ rb.n / rb.Wbar
    @test drift_critical_rate(Params(bh; GrowthRate = 1.5)).k < rb.k < drift_critical_rate(Params(bh; GrowthRate = 5.0)).k
    @test drift_critical_rate(bh; cost = log(2.0)).k == 0                # no growth left
    # sustainable_rate is the fixed-V_G threshold, several times higher
    @test sustainable_rate(p; VG = house_of_gauss_vg(p)) > 3r.k
    @test_throws ArgumentError drift_critical_rate(p; cost = -0.1)
    # the half-survival rate from the moments lowering with drift (the rare-alleles closure: fast)
    hs = half_survival_rate(p; loci = :rare_alleles, steps = 6, burnin = 500)
    @test 0 < hs.k50 < 2hs.kc && hs.kc ≈ r.k
    @test_throws ArgumentError half_survival_rate(Params(p; Plasticity = 0.5, PlasticityCostSD = 1.0); cost = 0.1)
end
