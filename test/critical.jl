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
    # ne = :breeders (validate/effective_size.md): N_e = the breeders, a lower critical rate; the default is unchanged
    @test drift_critical_rate(p; ne = :falconer) === r
    rn = drift_critical_rate(p; ne = :breeders)
    @test 0.7r.k < rn.k < 0.85r.k && rn.VG < r.VG                       # 0.78 for this input
    @test QTLSim._house_of_gauss(p, 2rn.n) ≈ rn.VG
    @test drift_critical_rate(bh; ne = :breeders).k < rb.k
    @test_throws ArgumentError drift_critical_rate(p; ne = :census)
end

@testset "drift_critical_rate(…; solver = :bisection) finds the steady state where the iteration does not (validate/effective_size.md, 1b)" begin
    # Xu–Shen-type inputs of validate/xu_shen.jl (diploid, 50 loci, μ = 0.1, V_E = ω²/50): the variance of a population of size K
    # carries a load near or above ln R, so the iteration from that variance overshoots to "no population"
    function xs(; r = 0.2, ω2 = 50.0, VM = 0.1, K = 1000.0)
        VE = ω2 / 50; ω = sqrt(ω2 - VE)
        R = exp(r - log(QTLSim.QTSIM_FITNESS_PEAK) + 0.5 * log(1 + VE / ω^2))
        Params(PopSize = round(Int, K * log(R) / 2), Ploidy = 2, NoQTL = 50, RecombPoissonMean = 3.0, PhenotypeVar = 2VE, Heritability = 0.5,
               MutationScalar = VM / VE, MutationRatePerGamete = 0.1, InitFitnessSD = ω, MoveType = :continuous,
               Demography = :beverton_holt, GrowthRate = R, CrossoverModel = :all)
    end
    # before the fix: 0, 0.2957, 0.4114; the maximum over all steady states (validate/xu_shen.log, "rung 0"): 0.2967, 0.4342, 0.4211
    for (q, k) in ((xs(ω2 = 20.0), 0.2967), (xs(VM = 0.2), 0.4342), (xs(K = 2000.0), 0.4211), (xs(), 0.3195), (xs(r = 0.5), 1.007))
        @test drift_critical_rate(q; solver = :bisection).k ≈ k rtol = 2e-3
    end
    # the default is the iteration, unchanged: 0, 0.2957, 0.4114 on the first three, and equal to the bisection where it converges
    @test drift_critical_rate(xs(ω2 = 20.0)).k == 0
    @test drift_critical_rate(xs(VM = 0.2)).k ≈ 0.2957 rtol = 2e-3
    @test drift_critical_rate(xs(K = 2000.0)).k ≈ 0.4114 rtol = 2e-3
    @test drift_critical_rate(xs()) === drift_critical_rate(xs(); solver = :iteration)
    @test drift_critical_rate(xs()).k ≈ drift_critical_rate(xs(); solver = :bisection).k rtol = 1e-12
    @test_throws ArgumentError drift_critical_rate(xs(); solver = :newton)
    q = xs(ω2 = 20.0)
    VK = QTLSim._house_of_gauss(q, 2.0 * q.PopSize)
    @test q.GrowthRate * QTLSim._feedback_size(q, VK, 0.0, 0.0).Wbar < 1          # size K with its own variance cannot replace itself
    @test QTLSim._feedback_state(q, 0.0, 0.0).VG == 0                             # the iteration: "no population"
    st = QTLSim._feedback_state(q, 0.0, 0.0; solver = :bisection)                 # but a smaller population with less variance can
    @test 0 < st.VG < VK && 2 < st.n < q.PopSize
    @test QTLSim._feedback_map(q, st.VG, 0.0, 0.0, :falconer) ≈ st.VG rtol = 1e-8  # a steady state: V = f(V)
    @test st.n ≈ q.PopSize * (q.GrowthRate * st.Wbar - 1) / (q.GrowthRate - 1) && st.N ≈ st.n / st.Wbar
    # f falls with V, so the root is unique
    fs = [QTLSim._feedback_map(q, V, 0.0, 0.0, :falconer) for V in range(0, VK, length = 400)]
    @test issorted(fs; rev = true) && count(i -> (fs[i] - (i - 1) * VK / 399) * (fs[i+1] - i * VK / 399) < 0, 1:399) == 1
    # no population at all: still 0
    @test drift_critical_rate(q; cost = 0.3, solver = :bisection).k == 0
    @test QTLSim._feedback_state(q, 0.0, 0.3; solver = :bisection).VG == 0
end
