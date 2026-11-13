# Tests: blocks. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl blocks`.

# ─── genomes as shared blocks, against the dense kernel ─────────────────

# The dense kernel as it was before blocks (test/dense_reference.jl), run beside the block
# kernel on the same seed, generation by generation: alleles bit for bit, genotypes to
# rounding, survivors, sizes and extinction exactly; and the pool's reference counts match
# the slots that name its blocks after every generation.
include(joinpath(@__DIR__, "dense_reference.jl"))
const Dense = DenseReference

function same_as_dense(p, seed; blocksize = QTLSim.default_blocksize(p.NoQTL))
    w = World(p; blocksize); d = Dense.World(p)
    rb = MersenneTwister(seed); rd = MersenneTwister(seed)
    init_env!(w, rb); Dense.init_env!(d, rd)
    ok = true; gerr = 0.0
    for t in 0:p.NoGenerations-1
        w.gen = d.gen = t
        generation!(w, rb); Dense.generation!(d, rd)
        ok &= w.extinct == d.extinct && w.nsurv == d.nsurv && w.N == d.N && w.θ == d.θ
        (w.extinct || !ok) && return (; ok, gerr, extinct = w.extinct, gen = t)
        ok &= genomes(w) == d.genome[:, :, 1:d.N] && w.env[1:w.N] == d.env[1:d.N]
        gerr = max(gerr, maximum(abs, w.genotype[1:w.N] .- d.genotype[1:d.N]))
        QTLSim.check_blocks(w)                 # throws if a count is off
    end
    (; ok, gerr, extinct = false, gen = p.NoGenerations)
end

@testset "block genomes are the dense kernel's alleles" begin
    base = Params(PopSize = 120, NoGenerations = 40, RecombPoissonMean = 3.0, PhenotypeVar = 0.03,
                  Heritability = 0.3, MutationScalar = 0.003, InitFitnessSD = 1.0,
                  MoveFitnessSDGen = 20, MoveFitnessSD = 0.02, MoveType = :continuous)
    n = 0
    for L in (1, 2, 7, 100, 1000), P in (1, 2), cm in (:qtsim, :all), seed in 1:2
        p = Params(base; NoQTL = L, Ploidy = P, CrossoverModel = cm,
                   MutationRatePerGamete = min(0.05, 2 / L),
                   PopSize = L == 1000 ? 60 : 120, NoGenerations = L == 1000 ? 20 : 40)
        r = same_as_dense(p, seed)
        @test r.ok && r.gerr < 1e-12 && !r.extinct
        n += 1
    end
    @test n == 40
    # block sizes at the extremes: one locus per block, and one block per chromosome
    for B in (1, 7), cm in (:qtsim, :all)
        r = same_as_dense(Params(base; NoQTL = 7, Ploidy = 2, CrossoverModel = cm, MutationRatePerGamete = 0.1), 3; blocksize = B)
        @test r.ok && r.gerr < 1e-12
    end
    # density regulation, with populations that persist and populations that die
    bh = Params(base; NoQTL = 50, Ploidy = 2, Demography = :beverton_holt, GrowthRate = 1.5,
                MutationRatePerGamete = 0.02, NoGenerations = 80)
    for seed in 1:3
        @test same_as_dense(bh, seed).ok
        r = same_as_dense(Params(bh; MoveFitnessSD = 0.3, InitFitnessSD = 0.3), seed)
        @test r.ok && r.extinct
    end
    # QTSim's extinction, haploid and diploid
    for P in (1, 2)
        r = same_as_dense(Params(base; NoQTL = 20, Ploidy = P, InitFitnessSD = 1e-3, MoveFitnessSDGen = 0,
                                 MoveFitnessSD = 5.0, MoveType = :phenotype, MutationRatePerGamete = 0.05), 1)
        @test r.ok && r.extinct
    end
end

@testset "block genomes: migration between demes as in the dense kernel" begin
    p = Params(PopSize = 80, NoQTL = 30, Ploidy = 2, NoGenerations = 25, RecombPoissonMean = 3.0,
               PhenotypeVar = 0.03, Heritability = 0.3, MutationScalar = 0.01,
               MutationRatePerGamete = 0.05, InitFitnessSD = 1.0)
    # demes sharing a pool (⊗), with separate pools, and with pools of different block sizes
    for build in (() -> p ⊗ demes(3), () -> Ext.Stratified([World(p) for _ in 1:3]),
                  () -> Ext.Stratified([World(p; blocksize = b) for b in (1, 6, 30)]))
        s = build(); ds = [Dense.World(p) for _ in 1:3]
        rb = MersenneTwister(4); rd = MersenneTwister(4)
        for (w, d) in zip(s.demes, ds); init_env!(w, rb); Dense.init_env!(d, rd); end
        ok = true
        for t in 0:p.NoGenerations-1
            for (w, d) in zip(s.demes, ds)
                w.gen = d.gen = t
                generation!(w, rb); Dense.generation!(d, rd)
            end
            # migrate the block demes, then move the dense individuals the same way: an
            # individual is identified by its environmental deviation, which moves with it
            where = Dict(w.env[i] => (k, i) for (k, w) in enumerate(s.demes) for i in 1:w.N)
            migrate!(s, rb; rate = 0.1)
            copy!(rd, rb)
            old = [copy(d.genome) for d in ds]; oldenv = [copy(d.env) for d in ds]
            for (k, w) in enumerate(s.demes), i in 1:w.N
                (k0, i0) = where[w.env[i]]
                ds[k].genome[:, :, i] .= old[k0][:, :, i0]; ds[k].env[i] = oldenv[k0][i0]
            end
            ok &= all(genomes(w) == d.genome[:, :, 1:d.N] for (w, d) in zip(s.demes, ds))
            ok &= all(maximum(abs, w.genotype[1:w.N] .- d.genotype[1:d.N]) < 1e-12 for (w, d) in zip(s.demes, ds))
            if s.demes[1].pool === s.demes[2].pool
                QTLSim.check_blocks(s.demes...)
            else
                foreach(QTLSim.check_blocks, s.demes)
            end
        end
        @test ok
    end
end

@testset "block pool: reference counts, no leaks, accessors" begin
    p = Params(PopSize = 200, NoQTL = 100, Ploidy = 2, NoGenerations = 400, RecombPoissonMean = 3.0,
               PhenotypeVar = 0.03, Heritability = 0.3, MutationScalar = 0.003,
               MutationRatePerGamete = 0.02, InitFitnessSD = 1.0)
    w = World(p); rng = MersenneTwister(1); init_env!(w, rng)
    sizes = Int[]; live = Int[]
    for t in 0:p.NoGenerations-1
        w.gen = t; generation!(w, rng)
        push!(live, QTLSim.check_blocks(w)); push!(sizes, length(w.pool.refs))
    end
    # live blocks stay far below the slots, and the pool stops growing
    @test maximum(live) < 2 * p.Ploidy * p.PopSize * w.pool.nb
    @test sizes[end] == sizes[100]
    # steady state: a generation allocates nothing
    allocs(w, rng) = @allocated generation!(w, rng)
    allocs(w, rng)
    @test allocs(w, rng) == 0
    # a block shared by parent and offspring is copied before a mutation, never changed in place
    q = Params(p; PopSize = 2, NoQTL = 9, MutationRatePerGamete = 0.5)
    v = World(q); G = randn(MersenneTwister(2), 9, 2, 2); setgenomes!(v, G)
    v.survivors = [1, 2]; v.nsurv = 2
    r = MersenneTwister(3); reproduce!(v, r)
    off = genomes(v; offspring = true)
    mutate!(v, r)
    @test genomes(v) == G                     # the parents are untouched
    @test genomes(v; offspring = true) != off # the offspring mutated
    QTLSim.check_blocks(v)
    # cached block sums are recomputed from the values: a genotype is the sum of the alleles
    develop!(v)
    @test v.genotype ≈ vec(sum(G; dims = (1, 2))) atol = 1e-14
    # accessors: shapes, errors, and a shared copy
    @test size(genomes(w)) == (p.NoQTL, p.Ploidy, w.N)
    @test_throws ArgumentError setgenomes!(v, zeros(8, 2, 2))
    @test_throws ArgumentError setgenomes!(v, zeros(9, 2, 3))
    @test_throws ArgumentError World(p; blocksize = 101)
    u = QTLSim.shared_copy(w)
    @test u.pool === w.pool && genomes(u) == genomes(w)
    QTLSim.check_blocks(w, u)
    @test QTLSim.default_blocksize.((1, 2, 7, 100, 5000)) == (1, 2, 3, 10, 71)
end

@testset "crossover model: :qtsim switches at most once, :all at every point" begin
    L = 50
    chrA = zeros(L); chrB = ones(L); dest = zeros(L)
    nswitch(g) = count(i -> g[i] != g[i+1], 1:L-1)
    rng = MersenneTwister(7)
    sq = [nswitch(gamete!(dest, chrA, chrB, rng, 5.0, L)) for _ in 1:4000]
    sa = [nswitch(gamete!(dest, chrA, chrB, rng, 5.0, L, true)) for _ in 1:4000]
    @test maximum(sq) == 1                         # QTSim: only the lowest point switches
    @test isapprox(mean(sq), 1 - exp(-5.0); atol=0.01)   # P(1 ≤ Poisson(5) < L)
    @test maximum(sa) >= 2
    @test isapprox(mean(sa), 5.0; atol=0.2)        # every one of Poisson(5) points switches
    # with one locus both models copy a whole parental chromosome
    d1 = zeros(1)
    @test all(_ -> gamete!(d1, [0.0], [1.0], rng, 5.0, 1)[1] in (0.0, 1.0), 1:100)
    @test_throws ArgumentError World(Params(CrossoverModel=:bogus))
end

@testset "multi-locus: run_composed == run! bitwise for both crossover models" begin
    for cm in (:qtsim, :all)
        p = Params(PopSize=200, NoQTL=20, Ploidy=2, NoGenerations=40, RecombPoissonMean=3.0,
                   PhenotypeVar=0.03, Heritability=0.3, MutationScalar=0.003,
                   MutationRatePerGamete=0.01, InitFitnessSD=1.0,
                   MoveFitnessSDGen=20, MoveFitnessSD=0.01, MoveType=:continuous, CrossoverModel=cm)
        a = run!(World(p), MersenneTwister(5))
        b = run_composed(World(p), MersenneTwister(5))
        @test a.geno_var == b.geno_var && a.nsurv == b.nsurv && a.theta == b.theta
    end
end

@testset "gamete carries half the parental genotype on average (mirror Lean gamete_expected_sum)" begin
    L = 30; rng = MersenneTwister(11)
    a = randn(rng, L); b = randn(rng, L); dest = zeros(L)
    for allpoints in (false, true)
        m = mean(sum(gamete!(dest, a, b, rng, 3.0, L, allpoints)) for _ in 1:40000)
        @test isapprox(m, (sum(a) + sum(b)) / 2; atol=0.1)    # ≈5 SE over 40000 gametes
    end
end
