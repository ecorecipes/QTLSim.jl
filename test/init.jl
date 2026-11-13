# Tests: init. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl init`.

@testset "populations initialised from theory (src/init.jl)" begin
    p = Params(Ploidy = 2, PopSize = 200, NoQTL = 20, MutationRatePerGamete = 1e-3, MutationScalar = 0.01,
               NoGenerations = 60, MoveFitnessSDGen = 40, MoveType = :continuous, MoveFitnessSD = 0.01)
    # the default is unchanged: simulate is run! on a World of zeros, which is init_env!
    # followed by NoGenerations generations, with the same draws
    o = simulate(p; seed = 11)
    @test !haskey(o, :start) && o.gen == 0:59
    w = World(p); rng = MersenneTwister(11); init_env!(w, rng)
    rec = QTLSim.records()
    for it in 0:p.NoGenerations-1
        w.gen = it; generation!(w, rng; rec)
    end
    @test all(k -> getfield(o, k) == getfield(rec, k), keys(rec))            # bitwise
    @test simulate(p; seed = 11, init = :zero) == o
    @test replicate(p, 2; seed = 11)[1] == o
    o2 = run!(World(p), MersenneTwister(11); resume = false)
    @test o2.geno_var == o.geno_var
    # the functions of the stationary spectrum
    @test QTLSim.sojourn_mass(0.0) == 1.0
    @test QTLSim.sojourn_mass(1e-3) ≈ 1 - 1e-3 / 6 rtol = 1e-7
    @test QTLSim.sojourn_mass(1e4) ≈ 2e-4 rtol = 1e-3
    @test QTLSim.sojourn_mass(79.9) ≈ QTLSim.sojourn_mass(80.1) rtol = 1e-2      # across the split of the quadrature
    @test QTLSim.sojourn_ratio(0.0, 0.4) == 1.0
    @test all(0 <= QTLSim.sojourn_ratio(γ, x) <= 1 for γ in (1e-6, 0.1, 5.0, 300.0, 1e5), x in (1e-5, 0.1, 0.5, 0.9, 1 - 1e-9))
    @test QTLSim.sojourn_ratio(50.0, 0.5) < QTLSim.sojourn_ratio(5.0, 0.5) < QTLSim.sojourn_ratio(0.5, 0.5)
    st = stationary_spectrum(p)
    @test 0 < st.VG < 2 * p.PopSize * muteVar(p)
    @test st.S ≈ initFitSD(p)^2 + st.VG + envVar(p)
    @test st.nparents ≈ p.PopSize * st.Wbar && st.gene ≈ st.VG / (2 * p.NoQTL)
    @test st.VG ≈ spectrum_vg(p) * (p.PopSize - 1) / p.PopSize rtol = 0.01    # the same model on a grid
    pn = Params(p; Drift = true)
    sn = stationary_spectrum(pn)
    @test sn.VG ≈ 2 * p.PopSize * muteVar(p) && isinf(sn.S)                    # neutral: 2N V_m
    @test QTLSim.stationary_sites(sn) ≈ 4 * p.PopSize * p.NoQTL * p.MutationRatePerGamete * sum(1 / k for k in 1:2*p.PopSize-1)
    nsites = QTLSim.stationary_sites(st)
    @test nsites < QTLSim.stationary_sites(sn)
    @test st.theta ≈ 4 * st.nparents * p.NoQTL * p.MutationRatePerGamete && st.L == 20 && st.chromosomes == 400
    @test QTLSim.stationary_density(st, 0.1, 0.2) > QTLSim.stationary_density(st, 0.1, 0.6) > 0
    @test QTLSim.stationary_mean_var(p, st) ≈ st.S / (st.nparents * (2 - st.VG / st.S))
    @test QTLSim.stationary_mean_var(pn, sn) == 0.0
    # the sampler: reproducible, of the right shape and mean, and with the expected variance
    G = equilibrium_genomes(p, MersenneTwister(5))
    @test size(G) == (20, 2, 200) && G == equilibrium_genomes(p, MersenneTwister(5))
    @test sum(equilibrium_genomes(p, MersenneTwister(5); mean = 0.25)) / 200 ≈ 0.25
    genic(G) = sum(var(view(G, l, 1, :) .+ view(G, l, 2, :)) for l in axes(G, 1))
    rng = MersenneTwister(6)
    gs = [genic(equilibrium_genomes(p, rng; target = :spectrum, mean = 0.0, spectrum = st)) for _ in 1:600]
    @test abs(mean(gs) - st.VG) < 4 * std(gs) / sqrt(600)
    me = multilocus_equilibrium(p; drift = true)
    s = QTLSim.sample_equilibrium(p, MersenneTwister(7))
    @test s.scale ≈ sqrt(me.genic / st.VG)                                     # target = :genic by default
    @test QTLSim.sample_equilibrium(p, MersenneTwister(7); target = 4 * st.VG).scale ≈ 2
    @test QTLSim.sample_equilibrium(p, MersenneTwister(7); target = :spectrum).scale == 1
    ns = [QTLSim.sample_equilibrium(p, rng; target = :spectrum, spectrum = st).sites for _ in 1:600]
    @test abs(mean(ns) - nsites) < 4 * sqrt(nsites / 600)
    @test_throws ArgumentError equilibrium_genomes(Params(p; NoQTL = 21), rng; spectrum = st)
    @test all(iszero, equilibrium_genomes(Params(p; MutationScalar = 0.0, InitFitnessSD = 5.0), rng))
    @test size(equilibrium_genomes(pn, rng)) == (20, 2, 200)                   # neutral: no scaling, no shift
    # the World: generation numbering, block bookkeeping, resumption
    @test QTLSim.first_event(p) == 40
    @test QTLSim.first_event(Params(p; ChangeFitnessGen = 25, ChangeFitnessSD = 2.0)) == 25
    @test QTLSim.first_event(Params(p; MoveFitnessSDGen = -1)) == 60
    @test QTLSim.protocol_start(p; burnin = 15) == 25 && QTLSim.protocol_start(p; burnin = 500) == 0
    @test QTLSim.protocol_start(p; burnin = 5, at = 12) == 7
    w = equilibrium_world(p, MersenneTwister(8); burnin = 15)
    @test w.gen == 25 && w.N == 200 && w.θ == 0.0 && !w.extinct
    @test QTLSim.check_blocks(w) > 0
    @test any(!iszero, w.env)
    develop!(w)
    @test var(w.genotype) > 0
    rng = MersenneTwister(8)
    w = equilibrium_world(p, rng; burnin = 15)
    o = run!(w, rng; resume = true)
    @test o.gen == 25:59 && !o.extinct
    @test all(iszero, o.theta[1:15]) && o.theta[16] ≈ 0.01 && o.theta[end] ≈ 0.2   # the move starts at generation 40
    @test QTLSim.check_blocks(w) > 0
    alloc(w, rng) = @allocated generation!(w, rng)
    w = equilibrium_world(Params(p; NoGenerations = 10^6, MoveFitnessSDGen = -1), rng; burnin = 0, at = 0)
    for _ in 1:300
        generation!(w, rng)
    end
    @test alloc(w, rng) == 0
    # simulate and replicate
    o = simulate(p; seed = 9, init = :equilibrium, burnin = 15)
    @test o.start == 25 && o.gen == 25:59 && o.seed == 9
    @test o.geno_var[1] > 0                                                    # variation from the first generation
    @test o == simulate(p; seed = 9, init = :equilibrium, burnin = 15)
    @test simulate(p; seed = 9, init = :equilibrium, burnin = 100).gen == 0:59 # burn-in longer than the protocol's
    @test simulate(p; seed = 9, init = :equilibrium, burnin = 5, at = 10).start == 5
    @test simulate(Params(p; Sample = 10); seed = 9, init = :equilibrium, burnin = 15).gen == [30, 40, 50]
    r = replicate(p, 3; seed = 9, init = :equilibrium, burnin = 15, summary = o -> (o.start, o.geno_var[1]))
    @test r[1] == (25, o.geno_var[1]) && allunique(last.(r))
    v0 = mean(replicate(p, 200; init = :equilibrium, burnin = 0, target = :spectrum, summary = o -> o.geno_var[1]))
    @test v0 ≈ st.VG rtol = 0.1
    # what is not supported, and bad arguments
    @test_throws ArgumentError stationary_spectrum(Params(p; Ploidy = 1))
    @test_throws ArgumentError equilibrium_genomes(Params(p; NoTraits = 2), rng)
    @test_throws ArgumentError equilibrium_world(Params(p; Demography = :beverton_holt), rng; burnin = 0)
    @test_throws ArgumentError stationary_spectrum(Params(p; Selection = :fecundity))
    @test_throws ArgumentError stationary_spectrum(Params(p; StandardisePhenotype = true))
    @test_throws ArgumentError simulate(p; seed = 1, init = :equilibrium)
    @test_throws ArgumentError simulate(p; seed = 1, init = :theory, burnin = 3)
    @test_throws ArgumentError equilibrium_genomes(p, rng; target = :total)
    @test_throws ArgumentError equilibrium_genomes(p, rng; target = -1.0)
    @test_throws ArgumentError equilibrium_genomes(p, rng; mean = :optimum)
    @test_throws ArgumentError equilibrium_world(p, rng; burnin = -1)
    # plasticity (src/plasticity.jl): a fixed reaction norm without a cost is the same model
    # while the optimum is at 0, so the draw is the same; a cost is not covered
    pp = Params(p; Plasticity = 0.5, CueLag = 2)
    @test stationary_spectrum(pp) == stationary_spectrum(p)
    @test equilibrium_genomes(pp, MersenneTwister(4)) == equilibrium_genomes(p, MersenneTwister(4))
    op = simulate(pp; seed = 6, init = :equilibrium, burnin = 10)
    @test op.start == 30 && !op.extinct && op.gen == 30:59
    @test_throws ArgumentError stationary_spectrum(Params(pp; PlasticityCostSD = 2.0))
    @test_throws ArgumentError equilibrium_world(Params(p; NoTraits = 2, ReactionNorm = true), rng; burnin = 5)
end
