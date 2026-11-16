# Tests: the haplotype-block lineage model (src/block_lineages.jl: qtsim_switch_law,
# block_map, run_blocks, block_stats, emulate(p; loci = :blocks),
# emulated_half_survival_rate(p; loci = :blocks)). Included by runtests.jl; run alone with
# `julia --project=vignettes test/runtests.jl block_lineages`.

@testset "QTSim's switch law and the block map" begin
    # the law of the lowest of K distinct points among M gaps, by enumeration (M = 6)
    L, u = 7, 2.5
    law = QTLSim.qtsim_switch_law(L, u)
    M = L - 1
    ref = zeros(M); pn = exp(-u)
    for K in 1:M
        pk = exp(-u) * u^K / factorial(K)
        cnt = zeros(Int, M)
        for s in 0:(1<<M)-1
            count_ones(s) == K && (cnt[trailing_zeros(s)+1] += 1)
        end
        ref .+= pk .* cnt ./ binomial(M, K)
    end
    pn += 1 - sum(exp(-u) * u^K / factorial(K) for K in 0:M)
    @test law.P ≈ ref atol = 1e-14
    @test law.pnone ≈ pn atol = 1e-14
    @test law.pnone + sum(law.P) ≈ 1 atol = 1e-12
    @test QTLSim.qtsim_switch_law(1, 3.0).pnone == 1.0

    p = Params(NoQTL = 500, RecombPoissonMean = 3.0)
    for B in (1, 7, 40, 160, 500, 2000)
        bm = QTLSim.block_map(p, B)
        @test bm.nb == min(B, 500)
        @test sum(bm.nloci) == 500 && all(>(0), bm.nloci) && issorted(bm.hi) && bm.hi[end] == 500
        @test bm.pnone + sum(bm.bmass) ≈ 1 atol = 1e-12
        @test length(bm.cum) == bm.nb - 1
    end
    # one block per locus: the kernel's map exactly
    bm = QTLSim.block_map(p, 500)
    @test all(==(1), bm.nloci)
    @test bm.bmass ≈ QTLSim.qtsim_switch_law(500, 3.0).P atol = 1e-15
    @test bm.pnone ≈ QTLSim.qtsim_switch_law(500, 3.0).pnone
    # :map gives each block about the same share of the switches; under QTSim's lowest-point
    # rule the blocks at the end of the chromosome hold many loci
    bm = QTLSim.block_map(Params(NoQTL = 1000, RecombPoissonMean = 3.0), 40)
    @test maximum(bm.bmass) / minimum(bm.bmass) < 1.5
    @test bm.nloci[end] > 5bm.nloci[1]
    bp = QTLSim.block_map(Params(NoQTL = 1000, RecombPoissonMean = 3.0), 40; partition = :position)
    @test all(==(25), bp.nloci)
    @test_throws ArgumentError QTLSim.block_map(p, 10; partition = :nope)
    @test_throws ArgumentError QTLSim.block_map(p, 0)
    # :all: the expected number of switches between two loci is kept up to the end blocks
    bm = QTLSim.block_map(Params(NoQTL = 200, RecombPoissonMean = 3.0, CrossoverModel = :all), 200)
    @test all(==(1), bm.nloci) && bm.snap == collect(1:199)
end

@testset "run_blocks" begin
    p = Params(Ploidy = 2, PopSize = 200, NoQTL = 50, NoGenerations = 400, MoveFitnessSDGen = 300, MoveType = :continuous,
               MoveFitnessSD = 0.01, PhenotypeVar = 1.0, Heritability = 0.3, MutationScalar = 0.005,
               MutationRatePerGamete = 0.002, RecombPoissonMean = 3.0, InitFitnessSD = 3.0, Sample = 1)
    r1 = run_blocks(p, Xoshiro(3); blocks = 10, burnin = 100)
    r2 = run_blocks(p, Xoshiro(3); blocks = 10, burnin = 100)
    @test r1.VG == r2.VG && r1.gen == r2.gen
    @test r1.start == 200 && r1.gen[1] == 200 && r1.gen[end] == 399
    @test !r1.extinct && r1.extinct_gen == -1
    @test all(isfinite, r1.VG) && all(>(0), r1.VG)
    # θ moves by MoveFitnessSD per generation from MoveFitnessSDGen on
    @test r1.theta[r1.gen .== 299][1] == 0.0
    @test r1.theta[end] ≈ 100 * 0.01
    r3 = run_blocks(p, Xoshiro(3); blocks = 10, burnin = 100, record_from = 350)
    @test r3.gen[1] == 350
    # unsupported models are ArgumentErrors
    for bad in (Params(p; Ploidy = 1), Params(p; MoveType = :phenotype), Params(p; Selection = :fecundity),
                Params(p; NoTraits = 2))
        @test_throws ArgumentError run_blocks(bad, Xoshiro(1); blocks = 5, burnin = 10)
    end
    # extinction: a population far from a fast-moving optimum dies, at a recorded generation
    q = Params(p; MoveFitnessSD = 1.0)
    r = run_blocks(q, Xoshiro(4); blocks = 10, burnin = 50)
    @test r.extinct && 300 <= r.extinct_gen <= 399 && r.gen[end] == r.extinct_gen
    # the equilibrium start is for a fixed size
    @test_throws ArgumentError run_blocks(Params(p; Demography = :beverton_holt, GrowthRate = 2.0), Xoshiro(5); blocks = 10, burnin = 50)
end

@testset "block model against the kernel (neutral, one block per locus and coarse)" begin
    # neutral: the block partition cannot change V_G's law (no selection, so linkage is
    # irrelevant to the variance of the sum); with one block per locus the model is the
    # kernel in law. Mean V_G after 300 generations from zero genomes, 40 runs each.
    p = Params(Ploidy = 2, PopSize = 100, NoQTL = 20, NoGenerations = 300, Drift = true, PhenotypeVar = 1.0, Heritability = 0.5,
               MutationScalar = 0.01, MutationRatePerGamete = 0.01, RecombPoissonMean = 3.0, Sample = 1)
    z = zeros(20, 2, 100)
    vk = [simulate(p; seed = 900 + r).geno_var[end] for r in 1:40]
    vb = map(1:40) do r
        s = QTLSim.BlockPop(20, 100, initFitSD(p))
        bm = QTLSim.block_map(p, 20)
        QTLSim._run_blocks!(s, p, bm, Xoshiro(r)).VG[end]
    end
    vc = map(1:40) do r
        bm = QTLSim.block_map(p, 3)
        s = QTLSim.BlockPop(bm.nb, 100, initFitSD(p))
        QTLSim._run_blocks!(s, p, bm, Xoshiro(100 + r)).VG[end]
    end
    se(v) = std(v) / sqrt(length(v))
    @test abs(mean(vb) - mean(vk)) < 4 * hypot(se(vb), se(vk))
    @test abs(mean(vc) - mean(vk)) < 4 * hypot(se(vc), se(vk))
    # one generation's mutational input from zero genomes: m·V_E per individual in expectation
    q = Params(p; PopSize = 2000, NoGenerations = 2)
    bm = QTLSim.block_map(q, 4)
    s = QTLSim.BlockPop(bm.nb, 2000, initFitSD(q))
    r = QTLSim._run_blocks!(s, q, bm, Xoshiro(7))
    @test r.VG[1] == 0.0
    @test r.VG[2] ≈ q.MutationScalar * envVar(q) rtol = 0.15
end

@testset "block_stats" begin
    rng = Xoshiro(11)
    Y = randn(rng, 6, 30)
    s = QTLSim.block_stats(Y, 6, 30)
    ybar = vec(mean(Y; dims = 2))
    @test s.blockgenic ≈ 2 * sum(mean((Y[b, :] .- ybar[b]) .^ 2) for b in 1:6)
    G = vec(sum(Y; dims = 1))
    a = sum(begin
                d = Y[b, :] .- ybar[b]; r = G .- d
                mean((r .- mean(r)) .* (d .^ 2 .- mean(d .^ 2)))
            end for b in 1:6)
    @test s.assoc ≈ a
end

@testset "emulate(p; loci = :blocks)" begin
    p = Params(Ploidy = 2, PopSize = 200, NoQTL = 50, NoGenerations = 1000, MoveFitnessSDGen = 500, MoveType = :continuous,
               MoveFitnessSD = 0.02, PhenotypeVar = 1.0, Heritability = 0.3, MutationScalar = 0.005,
               MutationRatePerGamete = 0.002, RecombPoissonMean = 3.0, InitFitnessSD = 3.0)
    e = emulate(p; loci = :blocks, blocks = 10, replicates = 6, seed = 42, burnin = 100, horizon = 200)
    @test isequal(e, emulate(p; loci = :blocks, blocks = 10, replicates = 6, seed = 42, burnin = 100, horizon = 200))
    @test e.replicates == 6 && 0 <= e.alive <= 6 && e.P2000 == e.alive / 6
    @test (e.T === missing) == (e.alive >= 3)
    @test haskey(e, :VG_move) && haskey(e, :hazard) && 0 <= e.P240k <= e.P2000
    @test_throws ArgumentError emulate(p; loci = :blocks, blocks = 10, replicates = 2, burnin = 10, horizon = 10)
    # a fast move kills every replicate, a slow one none
    ef = emulate(Params(p; MoveFitnessSD = 2.0); loci = :blocks, blocks = 10, replicates = 4, seed = 1, burnin = 50, horizon = 100)
    @test ef.P2000 == 0 && ef.T isa Real && ef.extinct
    es = emulate(Params(p; MoveFitnessSD = 1e-4); loci = :blocks, blocks = 10, replicates = 4, seed = 1, burnin = 50, horizon = 100)
    @test es.P2000 == 1 && es.T === missing
    # with an rng instead of a seed
    @test emulate(p; loci = :blocks, blocks = 10, replicates = 3, rng = Xoshiro(5), burnin = 50, horizon = 50).replicates == 3
    # emulate_sweep: one seed per input from job_seed; needs a seed
    ps = [p, Params(p; MoveFitnessSD = 0.05)]
    r = emulate_sweep(ps; ids = ["a", "b"], loci = :blocks, blocks = 10, replicates = 3, seed = 9, burnin = 50, horizon = 100)
    @test length(r) == 2 && r[1].P2000 == emulate(p; loci = :blocks, blocks = 10, replicates = 3, seed = QTLSim.job_seed(9, "a"),
                                                   burnin = 50, horizon = 100).P2000
    @test_throws ArgumentError emulate_sweep(ps; loci = :blocks, blocks = 10, replicates = 3, burnin = 50, horizon = 100)
    # the half-survival rate: within the bracket, P2000 ≥ 1/2 just below it and < 1/2 just above
    h = emulated_half_survival_rate(p; loci = :blocks, blocks = 10, replicates = 8, seed = 3, burnin = 100, horizon = 200)
    @test 0 < h.k50 < Inf && h.evals > 2
    @test_throws ArgumentError emulated_half_survival_rate(p; loci = :blocks, replicates = 2, burnin = 10, horizon = 10)
end

@testset "Gaussian segregation within blocks" begin
    p = Params(Ploidy = 2, PopSize = 200, NoQTL = 50, NoGenerations = 400, MoveFitnessSDGen = 300, MoveType = :continuous,
               MoveFitnessSD = 0.01, PhenotypeVar = 1.0, Heritability = 0.3, MutationScalar = 0.005,
               MutationRatePerGamete = 0.002, RecombPoissonMean = 3.0, InitFitnessSD = 3.0, Sample = 1)
    bm = QTLSim.block_map(p, 10; within = :gaussian)
    @test bm.within === :gaussian && length(bm.gf) == 49 && all(0 .< bm.gf .<= 1)
    @test count(==(1.0), bm.gf) == 9                      # the switch points on the 9 inner boundaries
    @test bm.gcum[end] == 1 && issorted(bm.gcum)
    @test all(==(1.0), QTLSim.block_map(p, 50; within = :gaussian).gf)   # one block per locus: every switch on a boundary
    r1 = run_blocks(p, Xoshiro(3); blocks = 10, within = :gaussian, burnin = 100)
    @test r1.VG == run_blocks(p, Xoshiro(3); blocks = 10, within = :gaussian, burnin = 100).VG
    @test all(isfinite, r1.VG)
    @test_throws ArgumentError QTLSim.block_map(Params(p; CrossoverModel = :all), 10; within = :gaussian)
    @test_throws ArgumentError QTLSim.block_map(p, 10; within = :nope)
    # the segregation law keeps the variance of a block: neutral, no mutation, 2 blocks of 25
    # loci, from a population with block variance v; one generation of gametes
    q = Params(p; Drift = true, MutationScalar = 0.0, PopSize = 4000, NoGenerations = 2, MoveFitnessSDGen = -1)
    bm2 = QTLSim.block_map(q, 2; within = :gaussian)
    s = QTLSim.BlockPop(2, 4000, initFitSD(q))
    s.Y .= randn(Xoshiro(1), 2, 8000)
    r = QTLSim._run_blocks!(s, q, bm2, Xoshiro(2))
    @test r.VG[2] ≈ r.VG[1] rtol = 0.08
    e = emulate(Params(p; NoGenerations = 1000, MoveFitnessSDGen = 500); loci = :blocks, blocks = 10, within = :gaussian,
                replicates = 3, seed = 4, burnin = 50, horizon = 100)
    @test e.replicates == 3
end
