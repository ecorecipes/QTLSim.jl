# Tests: demography. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl demography`.

@testset "density-regulated demography (Beverton-Holt)" begin
    p = Params(PopSize=500, NoQTL=10, Ploidy=2, NoGenerations=300, PhenotypeVar=0.02, Heritability=0.4,
               MutationScalar=0.004, MutationRatePerGamete=1e-3, RecombPoissonMean=2.0)
    o = simulate(p; seed = 1)
    @test all(==(500), o.popsize)                                   # QTSim: refilled to PopSize
    b = Params(p; Demography = :beverton_holt, GrowthRate = 2.0)
    ob = simulate(b; seed = 1)
    @test !ob.extinct && all(<=(500), ob.popsize) && 300 < mean(ob.popsize) < 500
    # the population settles near K (R W̄ − 1)/((R − 1) W̄)
    W = mean(ob.mean_fitness[100:end])
    @test mean(ob.popsize[100:end]) ≈ 500 * (2W - 1) / W rtol = 0.05
    # an optimum moving faster than it can follow: the population shrinks to extinction
    fast = Params(b; GrowthRate = 1.5, MoveFitnessSDGen = 50, MoveFitnessSD = 0.05, MoveType = :continuous)
    of = simulate(fast; seed = 2)
    @test of.extinct && of.popsize[end] < of.popsize[50] / 2          # it shrinks before it dies
    # the kernel and run_composed still agree, bit for bit
    oc = run_composed(World(b), QTLSim.Random.MersenneTwister(1))
    @test oc.popsize == ob.popsize && oc.geno_var == ob.geno_var
    # the ACSet holds the current population
    w = World(b); r = QTLSim.Random.MersenneTwister(3); init_env!(w, r)
    for _ in 1:5; generation!(w, r); end
    pop = lift(w)
    @test Catlab.nparts(pop, :Individual) == w.N
    w2 = World(b); lower!(pop, w2)
    @test w2.N == w.N && genomes(w2) == genomes(w)
    # the moments lowering predicts the equilibrium size, and extinction for a fast move
    m = run_moments(Params(b; NoGenerations = 2000))
    @test m.popsize[end] ≈ 500 * (2 * m.mean_fitness[end] - 1) / m.mean_fitness[end] atol = 1
    @test run_moments(Params(fast; NoGenerations = 2000)).extinct
    @test_throws ArgumentError World(Params(p; Demography = :logistic))
    @test_throws ArgumentError World(Params(p; GrowthRate = 0.0))
end

@testset "fecundity selection (Selection = :fecundity)" begin
    p = Params(PopSize=300, NoQTL=10, Ploidy=2, NoGenerations=300, PhenotypeVar=0.02, Heritability=0.4,
               MutationScalar=0.004, MutationRatePerGamete=1e-3, RecombPoissonMean=2.0)
    @test Params().Selection === :viability
    @test simulate(Params(p; Selection = :viability); seed = 4).geno_var == simulate(p; seed = 4).geno_var
    @test_throws ArgumentError World(Params(p; Selection = :hard))
    f = Params(p; Selection = :fecundity)
    o = simulate(f; seed = 1)
    @test !o.extinct && all(==(300), o.nsurv)                     # everyone survives selection
    @test o.geno_var != simulate(p; seed = 1).geno_var
    # parents are drawn in proportion to fitness; a parent of fitness 0 never
    w = World(f); r = MersenneTwister(2); init_env!(w, r); develop!(w); QTLSim.select!(w, r)
    @test w.weighted && w.wsum ≈ sum(w.fitness[1:w.N]) && w.survivors == 1:w.N
    w.fitness[1:w.N] .= 0.0; w.fitness[1:4] .= (0.1, 0.2, 0.3, 0.4)
    empty!(w.survivors); QTLSim.weigh!(w)
    draws = [QTLSim.weighted_parent(w, r) for _ in 1:40_000]
    @test all(<=(4), draws)
    @test [count(==(i), draws) / 40_000 for i in 1:4] ≈ [0.1, 0.2, 0.3, 0.4] atol = 0.01
    # fewer than two individuals of positive fitness: no two distinct parents
    w.fitness[2:w.N] .= 0.0; w.extinct = false; empty!(w.survivors); QTLSim.weigh!(w)
    @test w.extinct
    # nearly all the weight on one: the second parent is still another one, without hanging
    w.fitness[2] = 1e-300; w.fitness[3] = 3e-300; w.extinct = false; empty!(w.survivors); QTLSim.weigh!(w)
    @test !w.extinct
    ex = [QTLSim.weighted_parent_except(w, r, 1) for _ in 1:20_000]
    @test all(in((2, 3)), ex)
    @test count(==(3), ex) / 20_000 ≈ 0.75 atol = 0.02
    QTLSim.reproduce!(w, r)
    @test w.Noff == 300
    # the selection differential recorded is the mating weights' Σ W z/ΣW − z̄
    w = World(f); r = MersenneTwister(3); init_env!(w, r); rec = QTLSim.records()
    develop!(w); QTLSim.select!(w, r); QTLSim.reproduce!(w, r); QTLSim.record!(rec, w)
    z = w.phenotype[1:w.N]; W = w.fitness[1:w.N]
    @test rec.sel_diff[1] ≈ sum(W .* z) / sum(W) - mean(z)
    # with the same phenotypes, viability's mean survivors' differential is fecundity's
    # (its expectation, to O(1/N)); viability's is noisier
    w = World(Params(p; InitFitnessSD = 0.1)); r = MersenneTwister(5); init_env!(w, r)
    w.env .+= range(-0.1, 0.1, length = w.N); develop!(w)
    w.θ = 0.05
    sv = map(1:2000) do _
        QTLSim.select!(w, r); mean(w.phenotype[i] for i in w.survivors) - mean(w.phenotype[1:w.N])
    end
    w.p = Params(w.p; Selection = :fecundity); QTLSim.select!(w, r)
    sf = sum(w.fitness[1:w.N] .* w.phenotype[1:w.N]) / w.wsum - mean(w.phenotype[1:w.N])
    @test mean(sv) ≈ sf rtol = 0.02
    # Beverton–Holt: the expected offspring number is the viability formula at n = ΣW
    b = Params(f; Demography = :beverton_holt, GrowthRate = 2.0)
    w = World(b); r = MersenneTwister(6); init_env!(w, r); develop!(w); QTLSim.select!(w, r)
    w.fitness[1:w.N] .*= 0.2; empty!(w.survivors); QTLSim.weigh!(w)    # far below the cap K
    S = w.wsum
    @test mean(QTLSim.noffspring(w, r) for _ in 1:20_000) ≈ min(300, 2S / (1 + S / 300)) rtol = 0.01
    ob = simulate(b; seed = 7)
    @test !ob.extinct && 150 < mean(ob.popsize[100:end]) < 300
    # the composed run and a select step with its own Selection agree with the kernel, bit for bit
    L = lifecycle()
    for q in (f, b, Params(f; Ploidy = 1))
        ok = simulate(q; seed = 8)
        @test run_composed(World(q), MersenneTwister(8)).geno_var == ok.geno_var
        cyc = L.develop ⋅ withparams(L.select; Selection = :fecundity) ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅
              L.mutate ⋅ L.swap
        oc = run_composed(World(Params(q; Selection = :viability)), MersenneTwister(8); cycle = cyc)
        @test oc.geno_var == ok.geno_var && oc.popsize == ok.popsize
    end
    # an optimum moving too fast: at a fixed size fecundity selection is soft and the population
    # lives on; under density regulation it dies, as under viability selection
    fast = Params(f; MoveFitnessSDGen = 50, MoveFitnessSD = 0.05, MoveType = :continuous, NoGenerations = 400)
    @test !simulate(fast; seed = 9).extinct
    @test simulate(Params(fast; Demography = :beverton_holt, GrowthRate = 1.5); seed = 9).extinct
    # moments: the same mean and variance trajectory; with drift, fecundity erodes less
    mf = run_moments(Params(fast; Demography = :beverton_holt, GrowthRate = 1.5, NoGenerations = 100))
    mv = run_moments(Params(fast; Selection = :viability, Demography = :beverton_holt, GrowthRate = 1.5, NoGenerations = 100))
    n = length(mv.gen)          # viability dies when N W̄ < 2, fecundity when fewer than one offspring are expected
    @test length(mf.gen) >= n && mf.pheno_mean[1:n] == mv.pheno_mean && mf.geno_var[1:n] == mv.geno_var && mf.popsize[1:n] == mv.popsize
    @test mf.nsurv == mf.popsize && mv.nsurv != mv.popsize
    ev = run_moments(Params(p; NoGenerations = 500); erosion = true)
    ef = run_moments(Params(f; NoGenerations = 500); erosion = true)
    @test ef.geno_var[end] > ev.geno_var[end]
    @test run_moments(Params(f; NoGenerations = 500)).geno_var == run_moments(Params(p; NoGenerations = 500)).geno_var
end

@testset "demographic stochasticity along the moments trajectory" begin
    Dist = QTLSim.Distributions
    # a Galton–Watson process with Poisson(λ) offspring: no selection (W̄ = 1) and a capacity
    # so large that the per-capita rate is R = λ
    gw(λ; ploidy = 1) = Params(PopSize = 10^15, Ploidy = ploidy, Demography = :beverton_holt, GrowthRate = λ)
    for λ in (0.8, 1.0)                                     # λ ≤ 1: extinct surely in the limit
        d = extinction_distribution(gw(λ), ones(10_000); N0 = 1, method = :branching)
        @test issorted(d.pext)
        @test d.pext[end] ≈ 1 atol = (λ < 1 ? 1e-12 : 3e-4)
        @test d.pext[1] ≈ exp(-λ)                           # extinct in the first generation
    end
    λ = 1.5                                                 # λ > 1: q = exp(λ(q − 1))
    q = 0.0; for _ in 1:10_000; q = exp(λ * (q - 1)); end
    for N0 in (1, 3)
        d = extinction_distribution(gw(λ), ones(60); N0, method = :branching)
        @test d.pext[end] ≈ q^N0 rtol = 1e-10
        @test d.psurvive ≈ 1 - q^N0 rtol = 1e-9
        @test d.median === missing
    end
    # survival thins the offspring: q = 1 − w + w exp(R(q − 1))
    w, R = 0.6, 3.0
    qt = 0.0; for _ in 1:10_000; qt = 1 - w + w * exp(R * (qt - 1)); end
    @test extinction_distribution(gw(R), fill(w, 200); N0 = 1, method = :branching).pext[end] ≈ qt rtol = 1e-10
    @test_throws ArgumentError extinction_distribution(gw(R), fill(w, 10); N0 = 1)
    # the exact chain on sizes agrees with the branching process far below the capacity
    # (haploid, so that one survivor can mate)
    small = Params(gw(R); PopSize = 10^5)
    Wt = range(0.4, 0.2, length = 60)
    br = extinction_distribution(small, Wt; N0 = 5, method = :branching)
    ch = extinction_distribution(small, Wt; N0 = 5, method = :chain)
    @test maximum(abs.(br.pext .- ch.pext)) < 2e-3
    # the cap at K: the chain's sizes never exceed K, and it agrees with a Monte Carlo of
    # the kernel's demography, min(K, Poisson(R n/(1 + (R − 1) n/K))) from Binomial(N, W̄) survivors
    b = Params(PopSize = 20, Ploidy = 2, Demography = :beverton_holt, GrowthRate = 4.0)
    Wb = vcat(fill(0.9, 10), fill(0.3, 40))
    @test QTLSim.expected_offspring(b, 20) == 20 && QTLSim.expected_offspring(b, 5) ≈ 4 * 5 / (1 + 3 * 5 / 20)
    ch = extinction_distribution(b, Wb; method = :chain)
    rng = MersenneTwister(1); nrep = 40_000; died = zeros(Int, length(Wb))
    for _ in 1:nrep
        N = 20
        for t in eachindex(Wb)
            n = rand(rng, Dist.Binomial(N, Wb[t]))
            N = n < 2 ? 0 : min(20, rand(rng, Dist.Poisson(4n / (1 + 3n / 20))))
            N == 0 && (died[t] += 1; break)
        end
    end
    @test maximum(abs.(ch.pext .- cumsum(died) ./ nrep)) < 0.01
    # the fixed size: an independent binomial trial each generation
    f = Params(PopSize = 50, Ploidy = 2)
    Wf = [0.5, 0.1, 0.05, 0.02]
    hz = [Dist.cdf(Dist.Binomial(50, x), 1) for x in Wf]
    @test extinction_distribution(f, Wf).pext ≈ 1 .- cumprod(1 .- hz)
    @test extinction_distribution(f, Wf; from = 2).pext ≈ 1 .- cumprod(1 .- hz[3:4])
    @test_throws ArgumentError extinction_distribution(f, Wf; from = 4)
    @test_throws ArgumentError extinction_distribution(b, Wb; method = :exact)
    # over a run_moments trajectory: deterministic extinction is certain extinction
    pm = Params(PopSize = 500, NoQTL = 10, Ploidy = 2, NoGenerations = 600, PhenotypeVar = 0.02,
                Heritability = 0.4, MutationScalar = 0.004, MutationRatePerGamete = 1e-3,
                Demography = :beverton_holt, GrowthRate = 1.5, MoveFitnessSDGen = 300,
                MoveFitnessSD = 0.05, MoveType = :continuous)
    m = run_moments(pm)
    d = extinction_distribution(pm, m; from = 300, method = :chain)
    @test m.extinct && d.gen == 300:m.extinct_gen && d.pext[end] == 1 && d.median <= m.extinct_gen
    @test_throws ArgumentError extinction_distribution(Params(pm; Sample = 10), run_moments(Params(pm; Sample = 10)))
    # drift in the moments recursion is off by default, and erodes V_G when on
    @test run_moments(pm) == run_moments(pm; erosion = false)
    pf = Params(pm; Demography = :fixed, MoveFitnessSDGen = -1, NoGenerations = 20_000)
    me, m0 = run_moments(pf; erosion = true), run_moments(pf)
    @test me.geno_var[end] < m0.geno_var[end]
    @test me.geno_var[end] ≈ multilocus_equilibrium(pf; drift = true).VG rtol = 1e-6      # the default, rare alleles
    @test multilocus_equilibrium(pf; drift = true).VG ≈ rare_alleles_vg(pf; drift = true)
    @test run_moments(pf; erosion = true, loci = :multilocus).geno_var[end] ≈ gaussian_vg_loci(pf; drift = true) rtol = 1e-6
    # the legacy aggregate recursion's fixed point: V = (V − V²/s)(1 − 1/(2 N W̄)) + V_m, s = V_s + V
    me = run_moments(pf; erosion = true, loci = :aggregate)
    V, s = me.geno_var[end], selection_variance(pf) + me.geno_var[end]
    @test (V - V^2 / s) * (1 - 1 / (2 * 500 * me.mean_fitness[end])) + muteVar(pf) ≈ V rtol = 1e-6
end

@testset "evolutionary rescue (Gomulkiewicz and Holt 1995)" begin
    # no lag: a supercritical branching process from the start, q = 1 − W + W exp(R(q − 1))
    W, R = 0.9, 2.0
    q = 0.0; for _ in 1:10_000; q = 1 - W + W * exp(R * (q - 1)); end
    r = rescue_probability(; N0 = 1, R, d0 = 0.0, VG = 0.0, VE = 0.0, ω = 1.0, Wmax = W)
    @test r.rescue ≈ 1 - q rtol = 1e-10
    # rescue grows with the founding size and the genetic variance, falls with the shift
    kw = (R = 1.5, VE = 0.5, ω = 2.0)
    @test rescue_probability(; N0 = 100, d0 = 3.0, VG = 0.5, kw...).rescue <
          rescue_probability(; N0 = 1000, d0 = 3.0, VG = 0.5, kw...).rescue
    @test rescue_probability(; N0 = 100, d0 = 3.0, VG = 0.2, kw...).rescue <
          rescue_probability(; N0 = 100, d0 = 3.0, VG = 0.5, kw...).rescue
    @test rescue_probability(; N0 = 100, d0 = 4.0, VG = 0.5, kw...).rescue <
          rescue_probability(; N0 = 100, d0 = 3.0, VG = 0.5, kw...).rescue
    r = rescue_probability(; N0 = 100, d0 = 3.0, VG = 0.5, kw...)
    @test 0 < r.Nmin < 100 && r.tmin > 0
    # no growth at the optimum: no rescue
    @test rescue_probability(N0 = 100, R = 1.0, d0 = 1.0, VG = 0.5, VE = 0.5, ω = 2.0).rescue == 0
end
