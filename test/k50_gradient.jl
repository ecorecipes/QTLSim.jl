# Tests: the gradient of the emulator's half-survival rate on a target that does not saturate
# (src/k50_gradient.jl: QTLSim.log_cumulative_hazard, QTLSim.emulated_log_hazard, k50_gradient).
# Included by runtests.jl; run alone with `julia --project=vignettes test/runtests.jl k50_gradient`.

@testset "log cumulative hazard" begin
    p = Params(Ploidy = 2, PopSize = 1000)
    # where the product of survival factors is representable, exp(−H) is extinction_distribution's psurvive
    W = [0.02 * exp(-0.01t) for t in 0:199]
    d = extinction_distribution(p, W)
    c = QTLSim.log_cumulative_hazard(p, W)
    @test c.continued == 0
    @test exp(-exp(c.logH)) ≈ d.psurvive rtol = 1e-10
    @test 0 < d.psurvive < 1
    # the default path is the old product, bit for bit (it was not touched)
    let alive = 1.0
        for w in W
            alive *= 1 - QTLSim.cdf(QTLSim.Binomial(1000, w), 1)
        end
        @test d.psurvive === 1 - (1 - alive)
    end
    # where it saturates at 1 in Float64, log H is finite and close to the log of the summed hazards
    W1 = fill(0.2, 100)
    @test extinction_distribution(p, W1).psurvive === 1.0
    c1 = QTLSim.log_cumulative_hazard(p, W1)
    @test isfinite(c1.logH) && c1.logH < -200
    @test c1.logH ≈ log(100) + 999 * log(0.8) + log1p(999 * 0.2) rtol = 1e-12
    @test isfinite(QTLSim.log_cumulative_hazard(p, fill(0.9, 10)).logH)          # hazards of about 1e-1000
    # the closed form of one generation's hazard, against the binomial
    for w in (1e-4, 1e-3, 0.004, 0.02), (K, nmin) in ((1000, 2), (200, 1))
        h = QTLSim._too_few(K, w, nmin)
        @test QTLSim._log_too_few(K, w, nmin) ≈ log(h) rtol = 1e-10
        @test QTLSim._log_generation_hazard(K, w, nmin) ≈ log(-log1p(-h)) atol = 1e-8
    end
    @test QTLSim._log_generation_hazard(1000, 1e-12, 2) ≈ log(-log(1000 * 999 / 2 * 1e-24)) rtol = 1e-6   # h → 1
    @test QTLSim._log_generation_hazard(1000, 0.0, 2) == Inf && QTLSim._log_generation_hazard(1000, 1.0, 2) == -Inf
    # monotone and finite through the whole range of W̄
    ls = [QTLSim._log_generation_hazard(1000, w, 2) for w in exp.(range(log(1e-12), log(0.5); length = 400))]
    @test all(isfinite, ls) && issorted(ls; rev = true)
    # a trajectory that stops early is continued: exactly, when log W̄ is quadratic in time
    lw(t) = -3.0 - 0.2t - 0.05t^2
    Wq = [exp(lw(t)) for t in 0:29]
    full = QTLSim.log_cumulative_hazard(p, Wq)
    cut = QTLSim.log_cumulative_hazard(p, Wq[1:20]; generations = 30)
    @test cut.continued == 10
    @test cut.logH ≈ full.logH rtol = 1e-9
    @test QTLSim.log_cumulative_hazard(p, Wq[1:20]; generations = 25).logH < cut.logH
    @test QTLSim.log_cumulative_hazard(p, Wq; from = 5, generations = 10).logH < full.logH
    # a rising end is not extrapolated upward
    up = QTLSim.log_cumulative_hazard(p, [0.001, 0.002, 0.004]; generations = 6)
    @test up.logH ≈ QTLSim.log_cumulative_hazard(p, [0.001, 0.002, 0.004, 0.004, 0.004, 0.004]).logH
    @test_throws ArgumentError QTLSim.log_cumulative_hazard(Params(p; Demography = :beverton_holt, GrowthRate = 2.0), W)
    @test_throws ArgumentError QTLSim.log_cumulative_hazard(p, W; from = 200)
end

@testset "k50_gradient" begin
    base = Params(Ploidy = 2, PopSize = 200, NoQTL = 50, NoGenerations = 1000, MoveFitnessSDGen = 500, MoveType = :continuous,
                  MoveFitnessSD = 0.02, PhenotypeVar = 1.0, Heritability = 0.3, MutationScalar = 0.005,
                  MutationRatePerGamete = 0.002, RecombPoissonMean = 3.0, InitFitnessSD = 3.0)
    kw = (; burnin = 100, horizon = 200)
    scr = QTLSim.TWO_STAGE_SCREEN

    # the cheap closure: against a brute-force re-solve of emulate's own half-survival rate
    for ω in (3.0, 0.6)
        p = Params(base; InitFitnessSD = ω)
        g = k50_gradient(p; loci = :rare_alleles, kw...)
        @test g.defined && g.status === :ok && g.solver === :none && g.beta === nothing
        @test keys(g.elasticities) == (:PhenotypeVar, :Heritability, :MutationScalar, :InitFitnessSD)
        @test abs(g.h0) <= 0.02 && g.dh_dlogk < 0 && 0 < g.step <= 1e-2 && g.asymmetry < 0.2
        @test !g.straddle && !g.shift && g.fallback.base == -1
        h = emulated_half_survival_rate(p; kw...)
        @test g.k50 ≈ h.k50 rtol = 1e-3
        P = QTLSim.Para(q -> emulated_half_survival_rate(q; kw...).k50, p, (:PhenotypeVar, :InitFitnessSD))
        s = sensitivities(P; rel = 0.02)
        @test g.elasticities.PhenotypeVar ≈ s.PhenotypeVar atol = 3e-3
        @test g.elasticities.InitFitnessSD ≈ s.InitFitnessSD atol = 3e-3
        # the emulator reads V_P and h² only through V_E = (1 − h²)V_P (and m·V_E)
        @test g.elasticities.Heritability ≈ -0.3 / 0.7 * g.elasticities.PhenotypeVar rtol = 2e-3
        # a given bracket, other names, a smaller span
        g2 = k50_gradient(p; loci = :rare_alleles, names = (:InitFitnessSD,), bracket = (0.9g.k50, 1.1g.k50), dh = 0.125, kw...)
        @test g2.k50 ≈ g.k50 rtol = 1e-3
        @test g2.step < g.step
        @test g2.elasticities.InitFitnessSD ≈ g.elasticities.InitFitnessSD atol = 3e-3
        # a bracket that misses the root is moved toward it
        g3 = k50_gradient(p; loci = :rare_alleles, names = (:InitFitnessSD,), bracket = (1.02g.k50, 1.03g.k50), kw...)
        @test g3.defined
        @test g3.k50 ≈ g.k50 rtol = 1e-3
    end

    # the screen under strong selection: emulate's P2000 jumps to exactly 0 next to k50 (the
    # trajectory dies inside the horizon), where the target is finite and the gradient defined
    p = Params(base; InitFitnessSD = 0.6)
    g = k50_gradient(p; kw...)
    @test g.defined && g.solver === :beta && g.continued
    @test all(isfinite, values(g.elasticities))
    up = Params(p; MoveFitnessSD = g.k50 * exp(g.step))
    e = emulate(up; loci = :spectrum, spectrum = scr, kw...)
    l = QTLSim.emulated_log_hazard(up; kw...)
    @test e.P2000 == 0 && e.extinct                       # the old target: h = −Inf here
    @test l.extinct && l.continued > 0 && isfinite(l.h) && 0.05 < l.P2000 < 0.5
    @test l.logP ≈ log(l.P2000) && exp(l.log1mP) ≈ 1 - l.P2000
    dn = Params(p; MoveFitnessSD = g.k50 * exp(-g.step))
    e, l = emulate(dn; loci = :spectrum, spectrum = scr, kw...), QTLSim.emulated_log_hazard(dn; kw...)
    @test !l.extinct && l.continued == 0 && l.fallback == e.fallback
    @test l.P2000 ≈ e.P2000 rtol = 1e-9
    @test g.elasticities.Heritability ≈ -0.3 / 0.7 * g.elasticities.PhenotypeVar rtol = 5e-3
    @test isequal(k50_gradient(p; onswitch = :flag, kw...), g)        # no fallback switch here: the same either way
    # far below k50 emulate's P2000 is exactly 1 and the target still finite
    lo = Params(p; MoveFitnessSD = 0.5g.k50)
    @test emulate(lo; loci = :spectrum, spectrum = scr, kw...).P2000 == 1
    l = QTLSim.emulated_log_hazard(lo; kw...)
    @test isfinite(l.h) && l.h > 30 && l.logP < 0 && l.log1mP < -30

    # no root, and bad arguments
    dead = Params(base; InitFitnessSD = 0.005)             # dies in the burn-in
    @test QTLSim.emulated_log_hazard(dead; kw...).logH == Inf
    gd = k50_gradient(dead; loci = :rare_alleles, kw...)
    @test !gd.defined && gd.status === :no_k50 && all(isnan, values(gd.elasticities))
    @test_throws ArgumentError k50_gradient(base; onswitch = :other, kw...)
    @test_throws ArgumentError k50_gradient(base; names = (:NoQTL,), kw...)
    @test_throws ArgumentError k50_gradient(base; loci = :blocks, kw...)
    @test_throws ArgumentError k50_gradient(base; rel = 0, kw...)
    @test_throws ArgumentError QTLSim.emulated_log_hazard(Params(base; MoveType = :phenotype); kw...)
end
