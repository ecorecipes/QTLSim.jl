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
    # the defaults (validate/effective_size.md): N_e = the breeders, and the steady state by bisection where the iteration fails
    @test drift_critical_rate(p; ne = :breeders, solver = :bisection) === r
    @test QTLSim._house_of_gauss(p, 2r.n) ≈ r.VG                         # the house of Gauss of 2n genes, n the breeders
    @test drift_critical_rate(p; solver = :iteration) === r             # the iteration converges here: the same computation
    # ne = :falconer: the earlier effective size (about the harmonic mean of adults and breeders) and its values, which
    # validate/critical_rate_posthoc.log recorded for this input before the default changed (k_c 0.02909, M 1.751)
    rf = drift_critical_rate(p; ne = :falconer, solver = :iteration)
    @test rf.k ≈ 0.02909 rtol = 2e-4
    @test rf.M ≈ 1.751 rtol = 3e-4
    @test QTLSim._house_of_gauss(p, 2effective_size(rf.N, rf.n; scheme = :falconer)) ≈ rf.VG
    @test 0.7rf.k < r.k < 0.85rf.k && r.VG < rf.VG                      # fewer effective genes: 0.78 of it for this input
    @test rb.k < drift_critical_rate(bh; ne = :falconer).k
    @test half_survival_rate(p; loci = :rare_alleles, steps = 6, burnin = 500, ne = :falconer, solver = :iteration).kc ≈ rf.k
    @test_throws ArgumentError drift_critical_rate(p; ne = :census)
end

# The maximum over all steady states without solving V = f(V): parametrised by the variance V (as eq. 6 of Xu and Shen is).
# A steady state with variance V has the number of breeders n at which the house of Gauss of 2n genes is V (Hermisson and
# Wagner's A4 solved for 1/(2n)), so the mean fitness that leaves n breeders, so the lag λ with that mean fitness, and the
# response Vλ/(V + V_s). Nothing here calls the package's solver.
function all_states_maximum(p::Params; cost = 0.0, points = 200_001)
    L, μ, K = p.NoQTL, p.MutationRatePerGamete, p.PopSize
    ω, VE = initFitSD(p), QTLSim.envVar(p)
    Vs = ω^2 + VE
    α2 = muteVar(p) / (2L * μ)
    best = 0.0
    for V in range(0, 2L * 2μ * α2 / (α2 / (2Vs) + sqrt((α2 / (2Vs))^2 + 4μ * α2 / Vs)), length = points)[2:end-1]
        c = V / (2L)                                    # per-gene variance: c²/V_s + c(α²/(2V_s) + 1/(2n)) = μα²
        inv2n = (μ * α2 - c^2 / Vs) / c - α2 / (2Vs)
        inv2n > 0 || continue
        n = 1 / (2inv2n)
        n >= p.Ploidy || continue
        W = p.Demography === :fixed ? n / K : (1 + n * (p.GrowthRate - 1) / K) / p.GrowthRate
        s = Vs + V
        W0 = QTLSim.QTSIM_FITNESS_PEAK * ω / sqrt(s) * exp(-cost)
        W <= W0 || continue                              # no lag gives so many breeders
        λ = sqrt(2s * log(W0 / W))
        best = max(best, V * λ / (V + Vs))
    end
    best
end

@testset "drift_critical_rate is the maximum over all steady states (an independent parametrisation by V)" begin
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0, PhenotypeVar = 0.03, Heritability = 0.4,
               MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous)
    q = Params(p; NoQTL = 500, PhenotypeVar = 0.05, Heritability = 0.3, MutationScalar = 0.008, InitFitnessSD = 3.0)
    for base in (p, q), dem in ((Demography = :fixed,), (Demography = :beverton_holt, GrowthRate = 5.0),
                                (Demography = :beverton_holt, GrowthRate = 2.0), (Demography = :beverton_holt, GrowthRate = 1.5)),
        cost in (0.0, 0.1)
        x = Params(base; dem...)
        @test drift_critical_rate(x; cost).k ≈ all_states_maximum(x; cost) rtol = 1e-5
    end
end

@testset "drift_critical_rate finds the steady state where the iteration alone does not (validate/effective_size.md, 1b)" begin
    # Xu–Shen-type inputs of validate/xu_shen.jl (diploid, 50 loci, μ = 0.1, V_E = ω²/50): the variance of a population of size K
    # carries a load near or above ln R, so the iteration from that variance overshoots to "no population"
    function xs(; r = 0.2, ω2 = 50.0, VM = 0.1, K = 1000.0)
        VE = ω2 / 50; ω = sqrt(ω2 - VE)
        R = exp(r - log(QTLSim.QTSIM_FITNESS_PEAK) + 0.5 * log(1 + VE / ω^2))
        Params(PopSize = round(Int, K * log(R) / 2), Ploidy = 2, NoQTL = 50, RecombPoissonMean = 3.0, PhenotypeVar = 2VE, Heritability = 0.5,
               MutationScalar = VM / VE, MutationRatePerGamete = 0.1, InitFitnessSD = ω, MoveType = :continuous,
               Demography = :beverton_holt, GrowthRate = R, CrossoverModel = :all)
    end
    cases = (xs(ω2 = 20.0), xs(VM = 0.2), xs(K = 2000.0), xs(), xs(r = 0.5))
    # the default (breeders, bisection) is the maximum over all steady states
    for q in cases
        @test drift_critical_rate(q).k ≈ all_states_maximum(q) rtol = 1e-5
    end
    # with the earlier effective size, the values recorded in validate/xu_shen.log ("rung 0") and validate/effective_size_branch.log:
    # the iteration alone 0, 0.2957, 0.4114 on the first three; the maximum over all steady states 0.2967, 0.4342, 0.4211
    for (q, k) in zip(cases, (0.2967, 0.4342, 0.4211, 0.3195, 1.007))
        @test drift_critical_rate(q; ne = :falconer).k ≈ k rtol = 2e-3
    end
    old(q) = drift_critical_rate(q; ne = :falconer, solver = :iteration)
    @test old(xs(ω2 = 20.0)).k == 0
    @test old(xs(VM = 0.2)).k ≈ 0.2957 rtol = 2e-3
    @test old(xs(K = 2000.0)).k ≈ 0.4114 rtol = 2e-3
    # where the iteration converges the two solvers are the same computation
    @test drift_critical_rate(xs()) === drift_critical_rate(xs(); solver = :iteration)
    @test old(xs()).k ≈ drift_critical_rate(xs(); ne = :falconer).k rtol = 1e-12
    # the iteration alone misses steady states with the breeders too
    @test drift_critical_rate(xs(ω2 = 20.0); solver = :iteration).k == 0 < drift_critical_rate(xs(ω2 = 20.0)).k
    @test_throws ArgumentError drift_critical_rate(xs(); solver = :newton)
    q = xs(ω2 = 20.0)
    VK = QTLSim._house_of_gauss(q, 2.0 * q.PopSize)
    @test q.GrowthRate * QTLSim._feedback_size(q, VK, 0.0, 0.0).Wbar < 1          # size K with its own variance cannot replace itself
    for ne in (:breeders, :falconer)
        @test QTLSim._feedback_state(q, 0.0, 0.0; ne, solver = :iteration).VG == 0    # the iteration: "no population"
        st = QTLSim._feedback_state(q, 0.0, 0.0; ne)                                  # but a smaller population with less variance can
        @test st === QTLSim._feedback_state(q, 0.0, 0.0; ne, solver = :bisection)
        @test 0 < st.VG < VK && 2 < st.n < q.PopSize
        @test QTLSim._feedback_map(q, st.VG, 0.0, 0.0, ne) ≈ st.VG rtol = 1e-8        # a steady state: V = f(V)
        @test st.n ≈ q.PopSize * (q.GrowthRate * st.Wbar - 1) / (q.GrowthRate - 1) && st.N ≈ st.n / st.Wbar
        # f falls with V, so the root is unique
        fs = [QTLSim._feedback_map(q, V, 0.0, 0.0, ne) for V in range(0, VK, length = 400)]
        @test issorted(fs; rev = true) && count(i -> (fs[i] - (i - 1) * VK / 399) * (fs[i+1] - i * VK / 399) < 0, 1:399) == 1
    end
    # no population at all: still 0
    @test drift_critical_rate(q; cost = 0.3).k == 0
    @test QTLSim._feedback_state(q, 0.0, 0.3).VG == 0
end
