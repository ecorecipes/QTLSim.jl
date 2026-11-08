using QTLSim
using Test, Statistics, Random
using QTLSim: eigvals, Symmetric   # LinearAlgebra, through the package

# The categorical layer is the package extension QTLSimCatlabExt, loaded with Catlab.
import Catlab
using Catlab: ⋅, ⊗, dom, codom, Hom
const Ext = Base.get_extension(QTLSim, :QTLSimCatlabExt)
(; develop, select, reproduce, mutate, swap, generation) = operators()

# Package hygiene (ambiguities, piracy, stale dependencies, compat bounds), when Aqua is
# available, as it is under Pkg.test.
if Base.find_package("Aqua") !== nothing
    @eval using Aqua
    @testset "Aqua" begin
        Aqua.test_all(QTLSim; persistent_tasks = false)
    end
end

@testset "derived quantities (mirror Lean Heritability.lean)" begin
    p = Params(PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001)
    @test geneticVar(p) ≈ 0.5 * 0.02
    @test envVar(p) ≈ 0.02 - geneticVar(p)
    @test geneticVar(p) + envVar(p) ≈ p.PhenotypeVar            # V_P = V_G + V_E
    @test envVar(p) ≈ (1 - p.Heritability) * p.PhenotypeVar
    # ω = V_G / √(m·V_E)   (analyticalFitSD_eq)
    @test analyticalFitSD(p) ≈ geneticVar(p) / sqrt(p.MutationScalar * envVar(p))
    @test muteVar(p) ≈ p.MutationScalar * envVar(p)
end

@testset "mutational variance input ≈ m·V_E (mirror Lean mutationalVarianceInput_eq)" begin
    # Empirically the Poisson mutation operator injects ≈ muteVar per generation.
    p = Params(PopSize=20000, NoQTL=1, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.05, MutationRatePerGamete=0.01)
    w = World(p)
    rng = MersenneTwister(1)
    # the change per individual is a compound Poisson sum, very heavy-tailed (kurtosis ~150
    # here): one round's variance has a relative SE of about 9%, so average 25 rounds
    # (SE ≈ 2%) and allow 4 SE
    v = map(1:25) do _
        setgenomes!(w, zeros(p.NoQTL, p.Ploidy, p.PopSize); offspring = true)
        mutate!(w, rng)
        G = genomes(w; offspring = true)
        var(sum(@view G[:, :, j]) for j in 1:p.PopSize)
    end
    @test isapprox(mean(v), muteVar(p); rtol=0.08)
end

@testset "run smoke + extinction handling" begin
    p = Params(PopSize=300, NoQTL=1, Ploidy=2, NoGenerations=50,
               PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001, Seed=42)
    out = simulate(p)
    @test out.extinct == false
    @test length(out.gen) == 50
    @test all(out.geno_var .>= 0)

    # Guaranteed extinction: impossibly narrow fitness with a moving optimum.
    pe = Params(PopSize=200, NoQTL=1, Ploidy=2, NoGenerations=500,
                PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001,
                InitFitnessSD=1e-4, MoveFitnessSDGen=0, MoveFitnessSD=5.0, Seed=1)
    oe = simulate(pe)
    @test oe.extinct == true
end

@testset "reproducibility (same seed ⇒ identical output)" begin
    p = Params(PopSize=300, NoQTL=1, Ploidy=2, NoGenerations=40,
               PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001, Seed=7)
    a = simulate(p); b = simulate(p)
    @test a.geno_var == b.geno_var
    @test a.pheno_mean == b.pheno_mean
end

@testset "generation! composite advances one generation" begin
    p = Params(PopSize=100, NoQTL=1, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    w = World(p); rng = MersenneTwister(3)
    init_env!(w, rng)
    @test var(w.env) > 0
    generation!(w, rng)
    @test w.gen == 1
    @test w.extinct == false
end

@testset "run! is init_env! then NoGenerations of generation!" begin
    p = Params(PopSize=150, NoQTL=5, Ploidy=2, NoGenerations=60, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.004, MoveFitnessSDGen=20, MoveFitnessSD=1.0,
               ChangeFitnessGen=30, ChangeFitnessSD=0.5, Sample=4)
    a = run!(World(p), MersenneTwister(11))
    w = World(p); rng = MersenneTwister(11); rec = QTLSim.records()
    init_env!(w, rng)
    for it in 0:p.NoGenerations-1
        w.gen = it
        generation!(w, rng; rec)
    end
    @test a.gen == rec.gen == 0:4:56
    @test a.geno_var == rec.geno_var && a.theta == rec.theta && a.omega == rec.omega
    @test a.omega[end] == 0.5 && a.omega[1] == initFitSD(p)     # the width change is part of generation!
end

@testset "haploid runs, and extinction when no one survives" begin
    # a wide fitness function (all but certain survival) and frequent mutation, so the
    # outcome does not depend on the random stream, which differs between Julia versions
    p = Params(PopSize=200, NoQTL=3, Ploidy=1, NoGenerations=50, PhenotypeVar=0.02, MutationScalar=0.004,
               MutationRatePerGamete=0.05, InitFitnessSD=1.0)
    o = simulate(p; seed=3)
    @test !o.extinct && length(o.gen) == 50 && o.geno_var[end] > 0
    # one survivor is enough for a haploid, not for a diploid
    w1 = World(Params(p; PopSize=1, Drift=true)); select!(w1, MersenneTwister(1))
    @test w1.nsurv == 1 && !w1.extinct
    w0 = World(Params(p; PopSize=5, InitFitnessSD=1e-6)); w0.θ = 100.0
    develop!(w0); select!(w0, MersenneTwister(1))
    @test w0.nsurv == 0 && w0.extinct
    @test simulate(Params(p; InitFitnessSD=1e-6, MoveFitnessSDGen=0, MoveFitnessSD=50.0); seed=1).extinct_gen == 0
end

@testset "the optimum moves as QTSim moves it (SetFitnessMean)" begin
    base = Params(PopSize=50, Ploidy=2, PhenotypeVar=0.02, MoveFitnessSDGen=0, MoveFitnessSD=0.5,
                  InitFitnessSD=2.0)
    for (type, expected) in ((:phenotype, (z, w) -> z + sqrt(0.02) * 0.5),
                             (:fitness, (z, w) -> z + w.ω * 0.5),
                             (:continuous, (z, w) -> 0.3 + 0.5))
        w = World(Params(base; MoveType = type)); init_env!(w, MersenneTwister(2)); develop!(w)
        w.θ = 0.3; moveoptimum!(w)
        @test w.θ ≈ expected(mean(w.phenotype), w)
    end
    w = World(Params(base; MoveFitnessSDGen = 5)); develop!(w); w.θ = 0.3; moveoptimum!(w)
    @test w.θ == 0.3                                           # not before its generation
end

@testset "Drift and StandardisePhenotype in selection" begin
    base = Params(PopSize=200, Ploidy=2, PhenotypeVar=0.02, InitFitnessSD=0.1)
    rng = MersenneTwister(4)
    w = World(Params(base; Drift = true)); init_env!(w, rng); develop!(w); select!(w, rng)
    @test all(==(1.0), w.fitness) && w.nsurv == 200
    # QTSim's Gaussian with its 0.4: pdf(Normal(0, ω), x) / (0.4/ω) = exp(-x²/2ω²) / (0.4·√(2π))
    W(x, ω) = exp(-x^2 / (2ω^2)) / (0.4 * sqrt(2π))
    for std_ in (false, true)
        w = World(Params(base; StandardisePhenotype = std_)); init_env!(w, rng); develop!(w); w.θ = 0.05
        select!(w, rng)
        sd = std_ ? std(w.phenotype) : 1.0
        @test w.fitness ≈ [W((z - 0.05) / sd, w.ω) for z in w.phenotype]
    end
end

@testset "Sample and record = false" begin
    p = Params(PopSize=100, Ploidy=2, NoGenerations=35, PhenotypeVar=0.02, Sample=10)
    o = simulate(p; seed=1)
    @test o.gen == [0, 10, 20, 30]
    @test isempty(run!(World(p), MersenneTwister(1); record=false).gen)
    oc = run_composed(World(p), MersenneTwister(1))
    @test oc.gen == o.gen && oc.geno_var == o.geno_var
    @test isempty(run_composed(World(p), MersenneTwister(1); record=false).gen)
end

@testset "checkparams rejects what QTSim would not run" begin
    ok = Params(PopSize=10, Ploidy=2, PhenotypeVar=0.02)
    @test checkparams(ok) === ok
    for kw in ((Ploidy = 3,), (Ploidy = 0,), (PopSize = 1,), (NoQTL = 0,), (Sample = 0,),
               (NoGenerations = -1,), (PhenotypeVar = 0.0,), (Heritability = 1.5,),
               (Heritability = -0.1,), (MutationScalar = -1.0,), (MutationRatePerGamete = 0.0,),
               (RecombPoissonMean = -1.0,), (CrossoverModel = :some,), (MoveType = :Phenotype,),
               (MoveFitnessSDGen = -2,), (ChangeFitnessGen = -5,), (ChangeFitnessGen = 10,),
               (InitFitnessSD = -2.0,), (InitFitnessSD = 0.0,), (Heritability = 1.0,),
               (MutationScalar = 0.0,))
        @test_throws ArgumentError World(Params(ok; kw...))
    end
    @test World(Params(ok; Ploidy = 1, PopSize = 1)) isa World
    @test World(Params(ok; Heritability = 1.0, InitFitnessSD = 2.0)) isa World   # no V_E needs a width
    @test World(Params(ok; ChangeFitnessGen = 10, ChangeFitnessSD = 0.3)) isa World
end

# ─── the categorical layer (package extension) ───────────────────────────

@testset "Catlab loads the categorical layer" begin
    @test Ext isa Module
    @test generation isa Ext.Operator && length(operators()) == 8
end


@testset "Catlab ACSet schema + lift/lower round-trip" begin
    p = Params(PopSize=8, NoQTL=3, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    w = World(p); rng = MersenneTwister(11)
    G = zeros(p.NoQTL, p.Ploidy, p.PopSize)
    for k in eachindex(G); G[k] = randn(rng); end
    setgenomes!(w, G)
    for i in 1:p.PopSize; w.env[i] = randn(rng); end
    develop!(w)

    pop = lift(w)
    @test pop isa Ext.Population
    @test Catlab.nparts(pop, :Individual) == p.PopSize
    @test Catlab.nparts(pop, :Allele) == p.PopSize * p.NoQTL * p.Ploidy
    # genotype computed at the ACSet level via `incident(...,:owner)`
    @test genotypes(pop) ≈ w.genotype
    @test phenotypes(pop) ≈ w.phenotype
    # lower back: genome + env recovered exactly
    @test genome_from_acset(pop, p.NoQTL, p.Ploidy) == genomes(w) == G
    w2 = World(p); lower!(pop, w2)
    @test genomes(w2) == genomes(w)
    @test w2.env == w.env
    @test w2.genotype == w.genotype && w2.phenotype == w.phenotype   # recomputed by lower!
    @test_throws ArgumentError lower!(pop, World(Params(p; PopSize = 9)))
    @test_throws ArgumentError genome_from_acset(pop, p.NoQTL + 1, p.Ploidy)
    # haploid
    ph = Params(p; Ploidy = 1); wh = World(ph); setgenomes!(wh, randn(rng, ph.NoQTL, 1, ph.PopSize)); develop!(wh)
    wh2 = World(ph); lower!(lift(wh), wh2)
    @test genomes(wh2) == genomes(wh)
end

@testset "⋅ composition: run_composed == run! bitwise" begin
    p = Params(PopSize=250, NoQTL=1, Ploidy=2, NoGenerations=35,
               PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001)
    a = run!(World(p), MersenneTwister(99))
    b = run_composed(World(p), MersenneTwister(99))
    @test a.geno_var == b.geno_var
    @test a.pheno_mean == b.pheno_mean
    @test a.theta == b.theta
    @test a.nsurv == b.nsurv
    # the composite operator is really a composition of the base operators
    @test generation isa Ext.Operator
end

@testset "⊗ stratification conserves population and evolves each deme" begin
    p = Params(PopSize=150, NoQTL=1, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    s = p ⊗ demes(5)
    @test ndemes(s) == 5
    @test totalpop(s) == 5 * p.PopSize
    rng = MersenneTwister(21)
    for w in s.demes, i in 1:w.p.PopSize; w.env[i] = QTLSim.envsd(p) * randn(rng); end
    generation(s, rng)                      # one operator, applied per deme
    @test all(w -> w.gen == 1, s.demes)
    @test totalpop(s) == 5 * p.PopSize
end

@testset "migrate! moves individuals and conserves total population" begin
    p = Params(PopSize=200, NoQTL=1, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    s = p ⊗ demes(3)
    rng = MersenneTwister(31)
    for w in s.demes
        setgenomes!(w, randn(rng, p.NoQTL, p.Ploidy, p.PopSize))
    end
    before = totalpop(s)
    m = migrate!(s, rng; rate=0.2)
    @test m.count > 0
    @test totalpop(s) == before
    # migration is composable as an operator
    @test migration(0.05) isa Ext.Operator
end

@testset "⋅ migration composes on a metapopulation (regression)" begin
    # `generation ⋅ migration(rate)` must dispatch correctly on a Stratified:
    # the generation acts per deme, migration acts on the whole metapopulation.
    p = Params(PopSize=150, NoQTL=1, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    s = p ⊗ demes(3)
    rng = MersenneTwister(5)
    for w in s.demes, i in 1:w.p.PopSize; w.env[i] = QTLSim.envsd(p) * randn(rng); end
    step = generation ⋅ migration(0.1)
    @test step isa Ext.Operator
    step(s, rng)                              # must not error (was a MethodError before the fix)
    @test all(w -> w.gen == 1, s.demes)       # the generation advanced every deme
    @test totalpop(s) == 3 * p.PopSize        # migration conserved total population
end

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

@testset "QTSim options files: read_qtsim_options and write_qtsim_options" begin
    text = """
    # a comment line, as QTSim allows
    Sample 2
    popsize 300
    NoQTL 10
    Generations 1500
    Ploidy 2
    MoveFitnessSD 1000 0.05
    ChangeFitnessSD 500 0.4
    OptimumMoveType Continuous
    PhenotypeVar 0.03
    Heritability 0.4
    MutationScalar 0.005
    MutationRatePerGamete 0.0001
    RecombinationPoissonMean 3
    InitFitnessSD 2.5
    StandardisePhenotype
    Seed 1620324082
    """
    p = read_qtsim_options(IOBuffer(text))
    @test (p.Sample, p.PopSize, p.NoQTL, p.NoGenerations, p.Ploidy) == (2, 300, 10, 1500, 2)
    @test (p.MoveFitnessSDGen, p.MoveFitnessSD, p.ChangeFitnessGen, p.ChangeFitnessSD) == (1000, 0.05, 500, 0.4)
    @test p.MoveType == :continuous && p.RecombPoissonMean == 3.0 && p.InitFitnessSD == 2.5
    @test p.StandardisePhenotype && p.Seed == 1620324082
    io = IOBuffer(); write_qtsim_options(io, p)
    @test read_qtsim_options(IOBuffer(String(take!(io)))) == p          # a round trip
    @test read_qtsim_options(IOBuffer("StandardisePhenotype\nStandardisePhenotype\n")).StandardisePhenotype == false
    @test_throws ArgumentError read_qtsim_options(IOBuffer("NoSuchOption 3\n"))
    @test_throws ArgumentError read_qtsim_options(IOBuffer("OptimumMoveType sideways\n"))
    @test_throws ArgumentError read_qtsim_options(IOBuffer("MoveFitnessSD 1000\n"))
    @test (@test_logs (:warn,) read_qtsim_options(IOBuffer("Drift\n"))).Drift == false
    @test simulate(Params(p; NoGenerations = 5)).seed == p.Seed
    # the input files of the paper, when they are next to the package
    sims = joinpath(@__DIR__, "..", "..", "Repository", "Simulations")
    if isdir(sims)
        for d in filter(isdir, readdir(sims; join = true))
            q = read_qtsim_options(first(readdir(d; join = true)))
            @test checkparams(q) === q && q.NoGenerations == 250_000 && q.MoveType == :continuous
        end
    end
end

@testset "theory: equilibrium variance, lag, critical rate" begin
    p = Params(PopSize=1000, NoQTL=100, Ploidy=2, PhenotypeVar=0.02, Heritability=0.4,
               MutationScalar=0.004, MutationRatePerGamete=1e-4, InitFitnessSD=2.0)
    @test selection_variance(p) ≈ 2.0^2 + 0.012
    @test round(gaussian_vg(p), digits = 4) == 0.0139           # vignette 02's Gaussian prediction
    V = gaussian_vg_exact(p)                                    # the fixed point of the Gaussian recursion
    @test V^2 ≈ muteVar(p) * (selection_variance(p) + V) && V > gaussian_vg(p)
    @test effective_size(1000, 1000) ≈ 4000 / (2 * 999 / 1000 + 2) && effective_size(1000, 100) < 1000
    @test house_of_cards_vg(p) < house_of_cards_vg(Params(p; NoQTL = 1000)) < 4 * 1000 * 1e-4 * selection_variance(p)
    @test steady_lag(0.01; VG = 0.5, Vs = 1.5) ≈ 0.04
    @test steady_lag(Params(p; MoveFitnessSD = 0.01); VG = 0.5) ≈ 0.01 * (0.5 + selection_variance(p)) / 0.5
    kc(VG) = critical_rate(; n = 500, VG, ω = 5.0, VE = 0.01, sd50 = 0.1, N = 1000)
    @test 0 < kc(0.01) < kc(0.02) && isnan(critical_rate(; n = 2, VG = 0.01, ω = 5.0, VE = 0.01, sd50 = 0.1, N = 1000))
    Γ = QTLSim.Distributions.Gamma(1.201, 0.004635)
    m = [0.001, 0.004, 0.01]; r = [0.02, 0.03, 0.04]
    a = mutation_adjusted(r, m, Γ)
    @test mean(a ./ r) ≈ 1                                      # weights relative to their mean
    @test a ./ r ≈ QTLSim.Distributions.pdf.(Γ, m) ./ mean(QTLSim.Distributions.pdf.(Γ, m))
end

@testset "the moments lowering: the theory as a third functor" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=3000, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    # the default, loci = :auto, is the rare-alleles closure for one diploid trait with
    # mutations, and :multilocus (the multilocus Gaussian model) on request
    o = run_moments(p)
    @test QTLSim.moments_loci(p) === :rare_alleles && QTLSim.moments_loci(p, :multilocus) === :multilocus
    @test o.geno_var == run_moments(p; loci = :rare_alleles).geno_var == run_moments(p; loci = :auto).geno_var
    @test o.geno_var[end] ≈ rare_alleles_vg(p) rtol = 1e-6           # converges to the rare-alleles fixed point
    om = run_moments(p; loci = :multilocus)
    @test om.geno_var[end] ≈ gaussian_vg_loci(p) rtol = 1e-6         # and the multilocus model to its own
    # one diploid locus: half of selection's reduction (between the haplotypes) is not inherited
    @test gaussian_vg_loci(p) ≈ QTLSim.gaussian_vg_closed(2muteVar(p), selection_variance(p)) > gaussian_vg_exact(p)
    # the legacy aggregate recursion converges to gaussian_vg_exact, which the multilocus
    # model is for one haploid locus (or haploid clones of any number of loci)
    @test run_moments(p; loci = :aggregate).geno_var[end] ≈ gaussian_vg_exact(p) rtol = 1e-6
    ph = Params(p; Ploidy = 1, NoQTL = 20)
    @test QTLSim.moments_loci(ph) === :rare_alleles                   # haploids default to the closure for clones
    @test run_moments(ph).geno_var == run_moments(ph; loci = :rare_alleles).geno_var
    @test run_moments(ph; loci = :multilocus).geno_var ≈ run_moments(ph; loci = :aggregate).geno_var rtol = 1e-9
    @test gaussian_vg_loci(ph) ≈ gaussian_vg_exact(ph)
    @test_throws ArgumentError run_moments(p; loci = :other)
    @test_throws ArgumentError run_moments(p; loci = :aggregate, between = true)
    @test o.pheno_var[end] ≈ o.geno_var[end] + envVar(p)
    @test all(0 .< o.mean_fitness .<= 1) && o.nsurv[1] == round(Int, 1000 * o.mean_fitness[1])
    # a continuous move: the lag converges to Lande's steady lag
    pc = Params(p; MoveFitnessSDGen = 1000, MoveFitnessSD = 0.001, MoveType = :continuous)
    oc = run_moments(pc)
    @test oc.theta[end] - oc.pheno_mean[end] ≈ steady_lag(pc; VG = oc.geno_var[end]) rtol = 1e-6
    @test !oc.extinct
    # faster than the population can follow: deterministic extinction after the move
    ox = run_moments(Params(pc; MoveFitnessSD = 0.1))
    @test ox.extinct && ox.extinct_gen > 1000
    # drift: no selection, the variance grows by m·V_E per generation
    od = run_moments(Params(p; Drift = true, NoGenerations = 10))
    @test od.geno_var ≈ muteVar(p) .* (0:9)
    # the lowering is a functor: the generation's image is the composite of the steps' images
    steps = (:develop, :select, :reproduce, :record, :rescale, :mutate, :swap)
    Ext_ = Ext
    for cov in (nothing, Ext_.DemeCov(p, [1000.0], [0.003]; drift = false), Ext_.DemeCov(p, [1000.0], [0.003]))
        m1 = Ext_.Moments(p; VG = 0.003, cov); m2 = deepcopy(m1)
        lower(g, :moments)(m1, nothing)
        foldl(Catlab.compose, (lower(getfield(L, s), :moments) for s in steps))(m2, nothing)
        @test (m1.zbar, m1.VG, m1.θ, m1.gen) == (m2.zbar, m2.VG, m2.θ, m2.gen)
    end
    @test lower(g ⊗ g ⋅ L.migrate(2), :moments) isa Ext_.DemeOp        # several demes (see below)
    @test_throws ArgumentError run_moments(p; cycle = L.develop)
end

@testset "the multilocus Gaussian model: genic variance, LD, recombination, drift" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize = 1000, NoQTL = 50, Ploidy = 2, NoGenerations = 20_000, PhenotypeVar = 1.0, Heritability = 0.5,
               MutationScalar = 0.01, MutationRatePerGamete = 0.01)
    Vm, Vs = muteVar(p), selection_variance(p)
    # without drift: V²/(V_s + V) = 2L·V_m, whatever the crossover model
    V = gaussian_vg_loci(p)
    @test V^2 / (Vs + V) ≈ 2 * 50 * Vm
    @test V ≈ gaussian_vg_exact(Params(p; MutationScalar = 100 * p.MutationScalar, InitFitnessSD = initFitSD(p))) rtol = 1e-12
    for q in (p, Params(p; CrossoverModel = :all), Params(p; RecombPoissonMean = 5.0))
        @test run_moments(Params(q; NoGenerations = 100_000, Sample = 100); loci = :multilocus).geno_var[end] ≈ V rtol = 1e-6
        e = multilocus_equilibrium(q; loci = :multilocus)
        @test e.VG ≈ V rtol = 1e-6
        @test e.ld < 0 && e.genic ≈ e.VG - e.ld
        # the LD at equilibrium: −ω_c b/ρ_c per haplotype in each pair class, b = L V_m/2
        ω, ρ = QTLSim.pair_classes(q)
        @test sum(ω) ≈ 1 && ω[1] == 1 / 50 && ρ[1] == 0
        @test sum(ω[2:end] .* ρ[2:end]) / (1 - 1 / 50) ≈ pair_recombination(q) rtol = 1e-10
        @test e.ld ≈ -50 * Vm * sum(ω[2:end] ./ ρ[2:end]) rtol = 1e-6
    end
    # less recombination, more LD, the same V_G
    @test multilocus_equilibrium(Params(p; RecombPoissonMean = 0.5); loci = :multilocus).ld < multilocus_equilibrium(p; loci = :multilocus).ld
    # drift: erosion = true and between = true (one deme) are the same recursion, and its
    # fixed point is gaussian_vg_loci(p; drift = true), below the deterministic one
    a, b = run_moments(p; erosion = true, loci = :multilocus), run_moments(p; between = true, loci = :multilocus)
    @test a.geno_var == b.geno_var
    @test a.geno_var[end] ≈ gaussian_vg_loci(p; drift = true) rtol = 1e-6
    @test run_moments(p; erosion = true).geno_var == run_moments(p; between = true).geno_var    # also by default
    @test gaussian_vg_loci(p; drift = true) < V
    # selection: the genic reduction is the sum of the loci's shares, (H/L)²/s each, and the
    # survivors' variance is V − V²/s
    H, s, Lq = 0.7, 3.0, 50
    st = QTLSim.hap_init(2H, Inf, 2, 2)
    st2, β = QTLSim.hap_select(st, 1 / s, Inf, [1 / Lq, 1 - 1 / Lq], 2)
    @test st[1][1] - st2[1][1] ≈ sum((H / Lq)^2 / s for _ in 1:Lq) && β ≈ 2H / s
    @test QTLSim.hap_vpop(st2, Inf, 2) ≈ 2H - (2H)^2 / s
    # recombination turns the LD back into variance: one generation of mating at ρ = 1
    # restores the haplotypes' LE, losing only the between-haplotype part
    st3 = QTLSim.hap_reproduce(st2, Inf, [0.0, 1.0], 2)
    @test st3[1][2] ≈ st2[2][2] && st3[1][1] == st2[1][1]
    # pair classes: one locus has only the genic class; two loci are split by one point
    @test QTLSim.pair_classes(Params(p; NoQTL = 1)) == ([1.0], [0.0])
    ωa, ρa = QTLSim.pair_classes(Params(p; NoQTL = 2, RecombPoissonMean = 0.7))
    @test ωa ≈ [0.5, 0.5] && ρa ≈ [0.0, 0.7 * exp(-0.7)]
    # a continuous move: Lande's steady lag with the multilocus V_G
    pc = Params(p; NoGenerations = 30_000, MoveFitnessSDGen = 20_000, MoveFitnessSD = 0.01, MoveType = :continuous)
    oc = run_moments(pc; loci = :multilocus)
    @test oc.theta[end] - oc.pheno_mean[end] ≈ steady_lag(pc; VG = oc.geno_var[end]) rtol = 1e-6
    # several demes without migration are independent copies of one, drift or not
    for erosion in (false, true)
        @test run_moments(Params(p; NoGenerations = 300); cycle = g ⊗ g, erosion).geno_var[1:2:end] ==
              run_moments(Params(p; NoGenerations = 300); erosion).geno_var
    end
    @test_throws ArgumentError run_moments([p, Params(p; Ploidy = 1)]; cycle = g ⊗ g)
    # haploid drift is clonal: V_G loses 1/n per generation, not 1/(2n)
    ph = Params(p; Ploidy = 1, NoQTL = 1, Drift = true, NoGenerations = 30, MutationScalar = 0.0, InitFitnessSD = 1.0)
    oh = run_moments(ph; VG = 1.0, erosion = true)
    @test oh.geno_var ≈ [(1 - 1 / 1000)^t * 1000 / 999 for t in 0:29]      # expected sample variances
end

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

@testset "rare alleles: house of Gauss, house of cards, Gaussian and neutral limits" begin
    p = Params(PopSize = 1000, NoQTL = 10, Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5,
               MutationScalar = 0.01, MutationRatePerGamete = 1e-4, RecombPoissonMean = 2.0)
    Vm, Vs, L = muteVar(p), selection_variance(p), p.NoQTL
    @test QTLSim.mutation_effect_variance(p) ≈ Vm / (2L * 1e-4)
    @test_throws ArgumentError QTLSim.mutation_effect_variance(Params(p; MutationRatePerGamete = 0.0))
    # the closed form: rare mutations → Bürger's stochastic house of cards, and without drift
    # Turelli's 4Lμ V_s; many small ones → the Gaussian √(2L V_m V_s); V_s → ∞ → 2N_e V_m
    rare = Params(p; MutationRatePerGamete = 1e-7, MutationScalar = 0.01 * 1e-7 / 1e-4,     # the same α² and V_s, rarer
                  InitFitnessSD = initFitSD(p))
    @test house_of_gauss_vg(rare; nsurv = 950) ≈ house_of_cards_vg(rare; nsurv = 950) rtol = 1e-3
    @test house_of_gauss_vg(rare; drift = false) ≈ 4L * 1e-7 * Vs rtol = 1e-3
    many = Params(p; MutationRatePerGamete = 1.0)
    @test house_of_gauss_vg(many; drift = false) ≈ sqrt(2L * Vm * Vs) rtol = 1e-3
    neutral = Params(p; InitFitnessSD = 1e5)
    @test house_of_gauss_vg(neutral; nsurv = 1000) ≈ 2effective_size(1000, 1000) * Vm rtol = 1e-3
    @test house_of_gauss_vg(p) < house_of_cards_vg(p) && house_of_gauss_vg(p) < gaussian_vg_loci(p; drift = true)
    # the multilocus model with the rare-alleles closure: below the Gaussian, equal to it
    # when mutations are many and small, Turelli's house of cards when they are rare
    @test rare_alleles_vg(p) < gaussian_vg_loci(p)
    @test rare_alleles_vg(many) ≈ gaussian_vg_loci(many) rtol = 2e-3
    @test rare_alleles_vg(rare) ≈ 4L * 1e-7 * Vs rtol = 0.02
    @test rare_alleles_vg(rare; drift = true) ≈ house_of_cards_vg(rare; nsurv = 1000 * multilocus_equilibrium(rare; drift = true).Wbar) rtol = 0.03
    # V_s → ∞: both models are the neutral recursion, 2N V_m
    @test rare_alleles_vg(neutral; drift = true) ≈ gaussian_vg_loci(neutral; drift = true) rtol = 1e-3
    @test rare_alleles_vg(neutral; drift = true) ≈ 2 * 1000 * Vm rtol = 0.01
    # the moments lowering converges to it, without and with drift; it is the default for
    # one diploid trait with mutations, and the multilocus model stays available
    q = Params(p; NoGenerations = 100_000, Sample = 1000)
    @test run_moments(q; loci = :rare_alleles).geno_var[end] ≈ rare_alleles_vg(p) rtol = 1e-6
    @test run_moments(q; loci = :rare_alleles, erosion = true).geno_var[end] ≈ rare_alleles_vg(p; drift = true) rtol = 1e-6
    @test run_moments(q).geno_var[end] ≈ rare_alleles_vg(p) rtol = 1e-6
    @test run_moments(q; erosion = true).geno_var == run_moments(q; loci = :rare_alleles, erosion = true).geno_var
    @test run_moments(q; loci = :multilocus).geno_var[end] ≈ gaussian_vg_loci(p) rtol = 1e-6
    @test multilocus_equilibrium(p).VG == multilocus_equilibrium(p; loci = :rare_alleles).VG == rare_alleles_vg(p)
    @test multilocus_equilibrium(p; loci = :multilocus).VG ≈ gaussian_vg_loci(p) rtol = 1e-6
    # how the default resolves (QTLSim.moments_loci): the closure with locally favoured alleles
    # for several traits (the closure itself for one deme), the closure for one trait, haploid too
    @test QTLSim.moments_loci(p) === QTLSim.moments_loci(p, :auto) === :rare_alleles
    @test QTLSim.moments_loci(p, :aggregate) === :aggregate
    p2t = Params(p; NoTraits = 2, NoGenerations = 200)
    @test QTLSim.moments_loci(p2t) === :rare_local
    @test run_moments(p2t).gmatrix == run_moments(p2t; loci = :rare_alleles).gmatrix
    @test multilocus_equilibrium(p2t).VG == multilocus_equilibrium(p2t; loci = :rare_alleles).VG
    g2_ = lifecycle().generation ⊗ lifecycle().generation
    @test run_moments([p2t, p2t]; cycle = g2_).gmatrix == run_moments([p2t, p2t]; cycle = g2_, loci = :rare_local).gmatrix
    @test run_moments([p2t, p2t]; cycle = g2_, loci = :rare_alleles).gmatrix[1:2:end] ≈ run_moments(p2t; loci = :rare_alleles).gmatrix rtol = 1e-12
    @test QTLSim.moments_loci(Params(p; Ploidy = 1)) === :rare_alleles
    ph = Params(p; Ploidy = 1, NoGenerations = 200)
    @test run_moments(ph).geno_var == run_moments(ph; loci = :rare_alleles).geno_var
    @test multilocus_equilibrium(ph).VG == multilocus_equilibrium(ph; loci = :rare_alleles).VG == rare_alleles_vg(ph)
    @test QTLSim.moments_loci(Params(p; Ploidy = 1, NoTraits = 2)) === :rare_local
    @test_throws ArgumentError run_moments(q; loci = :bogus)
    @test run_moments(Params(p2t; NoGenerations = 50_000); loci = :rare_alleles).gmatrix[end] ≈ rare_alleles_G(p2t) rtol = 1e-6
    # demes: the closure acts in each deme
    L_ = lifecycle()
    g2 = L_.generation ⊗ L_.generation
    m2 = run_moments(Params(p; NoGenerations = 2000, Sample = 100); cycle = g2, loci = :rare_alleles)
    m1 = run_moments(Params(p; NoGenerations = 2000, Sample = 100); loci = :rare_alleles)
    @test m2.geno_var[m2.deme .== 1] ≈ m1.geno_var && m2.geno_var[m2.deme .== 2] ≈ m1.geno_var
    @test run_moments(Params(p; NoGenerations = 2000, Sample = 100); cycle = g2).geno_var == m2.geno_var    # the default
    # each deme has its own α²: demes of different mutation rates are independent copies of one
    pa, pb = Params(p; NoGenerations = 2000, Sample = 100), Params(p; NoGenerations = 2000, Sample = 100, MutationRatePerGamete = 1e-3)
    md = run_moments([pa, pb]; cycle = g2, between = true)
    @test md.geno_var[md.deme .== 1] ≈ run_moments(pa; between = true).geno_var rtol = 1e-12
    @test md.geno_var[md.deme .== 2] ≈ run_moments(pb; between = true).geno_var rtol = 1e-12

    # mutate steps that carry mutation fields: their α² (one law of effects), or the multilocus
    # model when they differ (an explicit :rare_alleles is then an error)
    gm = withparams(L_.generation; MutationRatePerGamete = 1e-3)
    @test run_moments(pa; cycle = gm).geno_var == run_moments(pa; cycle = gm, loci = :rare_alleles).geno_var
    @test run_moments(pa; cycle = gm).geno_var ≈ run_moments(pb).geno_var rtol = 1e-12
    gmix = gm ⋅ L_.generation
    @test run_moments(pa; cycle = gmix).geno_var == run_moments(pa; cycle = gmix, loci = :multilocus).geno_var
    @test_throws ArgumentError run_moments(pa; cycle = gmix, loci = :rare_alleles)
end
@testset "rare alleles: several traits and haploids (src/rare_alleles_traits.jl)" begin
    L_ = lifecycle(); g = L_.generation
    p = Params(PopSize = 1000, NoQTL = 10, Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, NoGenerations = 4000,
               MutationScalar = 0.01, MutationRatePerGamete = 1e-4, RecombPoissonMean = 2.0, Sample = 10)
    Vm, Vs, L, μ = muteVar(p), selection_variance(p), p.NoQTL, p.MutationRatePerGamete
    # n = 1: the matrix closure with 1×1 matrices is the one-trait closure
    @test QTLSim.mutation_effect_covariance(p) ≈ fill(QTLSim.mutation_effect_variance(p), 1, 1) rtol = 1e-14
    @test QTLSim.rare_alleles_rate(fill(0.3, 1, 1), fill(1 / 2.5, 1, 1))[1] ≈ 0.3 / 5 rtol = 1e-14
    for erosion in (false, true)
        a = run_moments(p; loci = :rare_alleles, erosion).geno_var
        b = Ext.run_trait_moments(p; loci = :rare_alleles, erosion).geno_var
        @test maximum(abs.(a .- b) ./ a[end]) < 1e-12
    end
    @test rare_alleles_G(p)[1] ≈ rare_alleles_vg(p) rtol = 1e-12
    @test rare_alleles_G(p; drift = true)[1] ≈ rare_alleles_vg(p; drift = true) rtol = 1e-12
    @test house_of_gauss_G(p)[1] ≈ house_of_gauss_vg(p) rtol = 1e-12
    @test house_of_gauss_G(p; drift = false)[1] ≈ house_of_gauss_vg(p; drift = false) rtol = 1e-12
    @test house_of_cards_G(p)[1] ≈ 4L * μ * Vs rtol = 1e-12
    # X = E[a aᵀ/(aᵀ s⁻¹ a)]: tr(s⁻¹X) = 1, X = s for one trait, the closed form for two
    # against Monte Carlo, and spherical pleiotropy dilutes the house of cards by 1/n
    C2, s2 = [1.0 0.6; 0.6 0.5], [2.0 -0.5; -0.5 1.0]
    X2 = QTLSim.rare_hc_matrix(C2, inv(s2))
    @test QTLSim.tr(inv(s2) * X2) ≈ 1 rtol = 1e-12
    @test QTLSim.rare_hc_matrix(fill(0.2, 1, 1), fill(1 / 3, 1, 1))[1] ≈ 3 rtol = 1e-12
    rng = MersenneTwister(11); Lc = QTLSim.cholesky(C2).L
    mc = sum(begin a = Lc * randn(rng, 2); a * a' ./ QTLSim.dot(a, s2 \ a) end for _ in 1:200_000) ./ 200_000
    @test mc ≈ X2 rtol = 0.01
    C3, s3 = [1.0 0.3 0.1; 0.3 0.8 -0.2; 0.1 -0.2 0.5], [1.5 0.2 0; 0.2 1 0.3; 0 0.3 2]
    @test QTLSim.tr(s3 \ QTLSim.rare_hc_matrix(C3, inv(s3))) ≈ 1 rtol = 1e-9
    @test QTLSim.rare_hc_rates([0.4, 0.4, 0.4]) ≈ fill(3 * 0.4 / 2, 3) rtol = 1e-9          # e = 1/(nλ)
    @test QTLSim.rare_hc_rates([0.2, 0.5]) ≈ [sqrt(0.2) * (sqrt(0.2) + sqrt(0.5)), sqrt(0.5) * (sqrt(0.2) + sqrt(0.5))] ./ 2
    sph = Params(p; NoTraits = 2)                                                          # MutCorr = SelCorr = 0
    @test house_of_cards_G(sph) ≈ (4L * μ * Vs / 2) .* [1 0; 0 1] rtol = 1e-12
    # two traits whose M and Ω do not share eigenvectors: the limits of the closed form and
    # of the closure in the multilocus model
    p2 = Params(p; NoTraits = 2, MutCov = Vm .* [1.0 0.5; 0.5 0.6], SelCorr = -0.4)
    A, M = selection_variance_matrix(p2), mutational_covariance(p2)
    rare = Params(p2; MutationRatePerGamete = 1e-7, MutationScalar = 0.01 * 1e-7 / 1e-4, MutCov = 1e-3 .* M,
                  InitFitnessSD = initFitSD(p2))                  # the same C and A, rarer
    @test house_of_gauss_G(rare; drift = false) ≈ house_of_cards_G(rare) rtol = 1e-3
    @test rare_alleles_G(rare) ≈ house_of_cards_G(rare) rtol = 0.02
    many = Params(p2; MutationRatePerGamete = 1.0)
    @test rare_alleles_G(many) ≈ gaussian_G_loci(many) rtol = 2e-3
    weak = Params(many; InitFitnessSD = 100.0)                   # G ≪ A: G A⁻¹ G = 2L M
    Gw = house_of_gauss_G(weak; drift = false)
    @test Gw * (selection_variance_matrix(weak) \ Gw) ≈ 2L .* mutational_covariance(weak) rtol = 1e-3
    neutral = Params(p2; InitFitnessSD = 1e5)
    @test house_of_gauss_G(neutral; nsurv = 1000) ≈ 2effective_size(1000, 1000) .* M rtol = 1e-3
    @test rare_alleles_G(neutral; drift = true) ≈ gaussian_G_loci(neutral; drift = true) rtol = 1e-3
    G2 = rare_alleles_G(p2)
    @test all(QTLSim.eigvals(G2) .> 0) && G2[1, 1] < gaussian_G_loci(p2)[1, 1]
    # the trait moments lowering converges to it, with and without drift; explicit
    # loci = :rare_alleles with traits, on several demes and through mutate steps
    q2 = Params(p2; NoGenerations = 60_000, Sample = 1000)
    @test run_moments(q2; loci = :rare_alleles).gmatrix[end] ≈ G2 rtol = 1e-6
    @test run_moments(q2; loci = :rare_alleles, erosion = true).gmatrix[end] ≈ rare_alleles_G(p2; drift = true) rtol = 1e-6
    gm = withparams(g, p2; MutationRatePerGamete = 1e-3)
    pb = Params(p2; MutationRatePerGamete = 1e-3, NoGenerations = 2000)
    @test run_moments(Params(p2; NoGenerations = 2000); cycle = gm, loci = :rare_alleles).gmatrix ≈
          run_moments(pb; loci = :rare_alleles).gmatrix rtol = 1e-12
    @test_throws ArgumentError run_moments(Params(p2; NoGenerations = 200); cycle = gm ⋅ g, loci = :rare_alleles)
    @test run_moments(Params(p2; NoGenerations = 200); cycle = gm ⋅ g).gmatrix ==
          run_moments(Params(p2; NoGenerations = 200); cycle = gm ⋅ g, loci = :multilocus).gmatrix
    md = run_moments([pb, pb]; cycle = g ⊗ g, between = true, loci = :rare_alleles)
    @test md.gmatrix[md.deme .== 1] ≈ run_moments(pb; between = true, loci = :rare_alleles).gmatrix rtol = 1e-12
    # haploids: clones are one unit of inheritance (QTLSim.hap_clonal), with U = Lμ
    ph = Params(p; Ploidy = 1)
    @test QTLSim.hap_clonal(1, [0.0, 0.3]) && QTLSim.hap_clonal(2, [0.0]) && !QTLSim.hap_clonal(2, [0.0, 0.3])
    phr = Params(rare; Ploidy = 1, NoTraits = 1, MutCov = zeros(0, 0), SelCorr = 0.0)
    Vsh = selection_variance(phr)
    @test house_of_gauss_vg(phr; drift = false) ≈ 2L * 1e-7 * Vsh rtol = 1e-3        # half the diploid's 4Lμ V_s
    @test rare_alleles_vg(phr) ≈ 2L * 1e-7 * Vsh rtol = 0.02
    phm = Params(ph; MutationRatePerGamete = 1.0)
    @test rare_alleles_vg(phm) ≈ gaussian_vg_loci(phm) rtol = 2e-3                     # √(V_m V_s): one clone
    @test house_of_gauss_vg(Params(ph; InitFitnessSD = 1e5); nsurv = 1000) ≈ effective_size(1000, 1000) * Vm rtol = 1e-3
    qh = Params(ph; NoGenerations = 100_000, Sample = 1000)
    oh = run_moments(qh; loci = :rare_alleles)
    @test oh.geno_var[end] ≈ rare_alleles_vg(ph) rtol = 1e-6
    @test run_moments(qh; loci = :rare_alleles, erosion = true).geno_var[end] ≈ rare_alleles_vg(ph; drift = true) rtol = 1e-6
    @test rare_alleles_vg(ph) ≈ house_of_gauss_vg(ph; drift = false) rtol = 0.1          # LD carries little for clones
    # several traits, haploid
    p2h = Params(p2; Ploidy = 1)
    @test house_of_gauss_G(p2h)[1, 1] < house_of_gauss_G(p2)[1, 1]
    @test run_moments(Params(p2h; NoGenerations = 60_000, Sample = 1000); loci = :rare_alleles).gmatrix[end] ≈ rare_alleles_G(p2h) rtol = 1e-6
end

@testset "allele-frequency spectrum (loci = :spectrum): neutral, equilibrium, moving-optimum limits" begin
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
               MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0)
    Vm = muteVar(p)
    # no selection: the neutral spectrum holds 2N V_m (Lynch and Hill 1986), and a move changes nothing
    e0 = spectrum_equilibrium(Params(p; Drift = true))
    @test e0.VG ≈ 2 * 1000 * Vm rtol = 0.01
    sp = e0.spectrum
    @test QTLSim.spectrum_select(sp, 3.0, Inf) == 0.0 && QTLSim.spectrum_linked(sp, 1000) == 0.0
    # no move: the lowering converges to the stationary spectrum, which for rare mutations in
    # linkage equilibrium is close to Hermisson and Wagner's stochastic house of Gauss
    e = spectrum_equilibrium(p)
    @test e.VG == spectrum_vg(p)
    @test e.VG ≈ house_of_gauss_vg(p; nsurv = e.nparents) rtol = 0.05
    o = run_moments(Params(p; NoGenerations = 30_000, Sample = 1000); loci = :spectrum)
    @test o.geno_var[end] ≈ e.VG rtol = 2e-3
    @test !o.extinct && all(o.theta .== 0)
    # the start of a move: the stationary spectrum is symmetric in the effects, so selection
    # moves the mean by V_G λ/S exactly (no skew), Lande's Gaussian response
    S = initFitSD(p)^2 + envVar(p) + QTLSim.spectrum_vg(e.spectrum)
    @test QTLSim.spectrum_select(e.spectrum, 0.7, S) ≈ QTLSim.spectrum_vg(e.spectrum) * 0.7 / S rtol = 1e-8
    # a population started at V_G takes the stationary shape scaled to V_G
    @test run_moments(Params(p; NoGenerations = 1, Sample = 1); loci = :spectrum, VG = 0.05).geno_var[1] ≈ 0.05 * 1000 / 999 rtol = 1e-10
    # a fast move: the lag decides, within a few hundred generations, and the time to extinction
    # is the Gaussian closure's (rare alleles, with drift) to within 10%
    pm = Params(p; NoGenerations = 7000, Sample = 1, MoveFitnessSDGen = 5000, MoveType = :continuous, MoveFitnessSD = 0.2)
    ms, mr = run_moments(pm; loci = :spectrum), run_moments(pm; loci = :rare_alleles, erosion = true)
    @test ms.extinct && mr.extinct
    @test abs(ms.extinct_gen - mr.extinct_gen) <= 0.1 * (mr.extinct_gen - 5000)
    # a slow move in the house-of-cards regime: directional selection drives the rare alleles of
    # large effect to intermediate frequency and V_G rises (Barton and Turelli 1987), with or
    # without the linked-selection term; the rare-alleles closure keeps its stationary value
    pH = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01,
                MutationRatePerGamete = 1e-4, RecombPoissonMean = 2.0, NoGenerations = 16_000, Sample = 10,
                MoveFitnessSDGen = 10_000, MoveType = :continuous, MoveFitnessSD = 0.01)
    for lk in (:barton, :none)
        h = run_moments(pH; loci = :spectrum, spectrum = (linkage = lk,))
        @test h.geno_var[end] > 1.3 * h.geno_var[findfirst(==(9990), h.gen)]
    end
    hr = run_moments(pH; loci = :rare_alleles, erosion = true)
    @test hr.geno_var[end] ≈ hr.geno_var[findfirst(==(9990), hr.gen)] rtol = 0.02
    # linked selection (Barton 1995): N/N_e − 1 grows as b² V_G, nothing with free recombination
    spl = QTLSim.AlleleSpectrum(p); spl.n .= sp.n
    QTLSim.spectrum_select(spl, 2.0, 25.0); c2 = QTLSim.spectrum_linked(spl, 1000)
    QTLSim.spectrum_select(spl, 4.0, 25.0); c4 = QTLSim.spectrum_linked(spl, 1000)
    @test c2 > 0 && c4 ≈ 4c2
    spf = QTLSim.AlleleSpectrum(Params(p; CrossoverModel = :all, RecombPoissonMean = 50.0)); spf.n .= sp.n
    QTLSim.spectrum_select(spf, 2.0, 25.0)
    @test QTLSim.spectrum_linked(spf, 1000) < 0.01c2
    # the slow-move collapse: under QTSim's crossover the variance falls more than with free recombination
    pv = Params(p; NoGenerations = 6000, Sample = 10, MoveFitnessSDGen = 5000, MoveType = :continuous, MoveFitnessSD = 0.02)
    vl = run_moments(pv; loci = :spectrum).geno_var[end]
    vf = run_moments(pv; loci = :spectrum, spectrum = (linkage = :none,)).geno_var[end]
    @test vl < 0.5vf
    # limits of the option
    g2 = lifecycle().generation ⊗ lifecycle().generation
    @test_throws ArgumentError run_moments([p, p]; cycle = g2, loci = :spectrum)
    @test_throws ArgumentError run_moments(Params(p; NoGenerations = 10); loci = :spectrum, between = true)
    @test_throws ArgumentError run_moments(Params(p; NoGenerations = 10, Ploidy = 1); loci = :spectrum)
    @test_throws ArgumentError run_moments(Params(p; NoGenerations = 10, NoTraits = 2); loci = :spectrum)
    @test_throws ArgumentError QTLSim.AlleleSpectrum(p; linkage = :bogus)
    # one law of effects from the mutate steps: their fields are the spectrum's
    gm = withparams(lifecycle().generation; MutationRatePerGamete = 1e-3)
    q = Params(p; NoGenerations = 300, Sample = 50)
    @test run_moments(q; cycle = gm, loci = :spectrum).geno_var ≈ run_moments(Params(q; MutationRatePerGamete = 1e-3); loci = :spectrum).geno_var
    @test_throws ArgumentError run_moments(q; cycle = gm ⋅ lifecycle().generation, loci = :spectrum)
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

@testset "migration graphs and per-deme parameters" begin
    @test island_migration(3, 0.1) == [0 0.05 0.05; 0.05 0 0.05; 0.05 0.05 0]
    @test stepping_stone(4, 0.2; circular = false) == [0 0.1 0 0; 0.1 0 0.1 0; 0 0.1 0 0.1; 0 0 0.1 0]
    @test stepping_stone(4, 0.2) == [0 0.1 0 0.1; 0.1 0 0.1 0; 0 0.1 0 0.1; 0.1 0 0.1 0]
    @test_throws ArgumentError QTLSim.check_migration([0.0 1.5; 0.0 0.0], 2)
    @test_throws ArgumentError QTLSim.check_migration(ones(2, 2), 2)
    @test_throws ArgumentError QTLSim.check_migration(zeros(3, 3), 2)
    # exchanges only along the edges of the graph: here deme 1 with deme 2, and nothing else
    p = Params(PopSize=200, NoQTL=1, Ploidy=2, NoGenerations=30, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    s = p ⊗ demes(3)
    for (d, w) in enumerate(s.demes); setgenomes!(w, fill(Float64(d), p.NoQTL, p.Ploidy, p.PopSize)); end
    M = zeros(3, 3); M[1, 2] = 0.5
    c = migrate!(s, MersenneTwister(1); rate = M)
    @test c.count > 0 && totalpop(s) == 600
    G1, G2, G3 = genomes.(s.demes)
    @test all(==(3.0), G3) && 2.0 in G1 && 3.0 ∉ G1
    @test count(==(2.0), G1) == count(==(1.0), G2) > 0   # exchanges, not moves
    # a run along a graph, and demes with their own parameters
    L = lifecycle(); g = L.generation
    cycle = g ⊗ g ⊗ g ⋅ L.migrate(3)
    o = run_composed(p ⊗ demes(3), MersenneTwister(2); cycle, migration_rate = stepping_stone(3, 0.1))
    @test length(o.gen) == 90 && !o.extinct
    s2 = [p, Params(p; InitFitnessSD = 0.2), Params(p; PopSize = 100)] ⊗ demes(3)
    @test [w.ω for w in s2.demes] == [0.3, 0.2, 0.3] && totalpop(s2) == 500
    o2 = run_composed(s2, MersenneTwister(3); cycle, migration_rate = island_migration(3, 0.05))
    @test count(==(3), o2.deme) == 30 && totalpop(s2) == 500
    @test_throws ArgumentError [p, p] ⊗ demes(3)
    @test_throws ArgumentError [p, Params(p; NoGenerations = 10)] ⊗ demes(2)
    # Eidos: the matrix is a constant, read row by row
    sc = slim_script(cycle, p; migration_rate = stepping_stone(3, 0.02; circular = false))
    @test occursin("defineConstant(\"MIG\", c(0.0, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0))", sc)
    @test !occursin("migRate", sc)
    @test_throws ArgumentError slim_script(cycle, p; migration_rate = ones(3, 3))
end

@testset "migration commutes with relabelling the demes, in distribution" begin
    # Lean `randomOrderSweep_exchangeable`: with the demes visited in a random order,
    # migrating a relabelled metapopulation is migrating the original and relabelling.
    # Two demes of two individuals at rate 1: with the demes always visited in the order
    # 1, 2, an arrangement has probability 1/8 one way and 0 the other (frequencies differ
    # by about 0.13 over these seeds). Individuals are labelled by their `env`.
    p = Params(PopSize = 2, NoQTL = 1, Ploidy = 2)
    function outcomes(start, rate; n = 20_000)
        tally = Dict{NTuple{4,Float64},Int}()
        for seed in 1:n
            s = p ⊗ demes(2)
            s.demes[1].env .= start[1]; s.demes[2].env .= start[2]
            migrate!(s, MersenneTwister(seed); rate)
            key = (s.demes[1].env..., s.demes[2].env...)
            tally[key] = get(tally, key, 0) + 1
        end
        tally
    end
    for rate in (1.0, [0.0 1.0; 1.0 0.0])
        a = outcomes(([1.0, 2.0], [3.0, 4.0]), rate)
        b = outcomes(([3.0, 4.0], [1.0, 2.0]), rate)
        b = Dict((x[3], x[4], x[1], x[2]) => v for (x, v) in b)   # relabel back
        @test sum(values(a)) == sum(values(b)) == 20_000
        @test maximum(abs(get(a, x, 0) - get(b, x, 0)) for x in union(keys(a), keys(b))) / 20_000 < 0.03
    end
end

@testset "neutral differentiation among demes: F_ST, Q_ST and variance components" begin
    N, Vm = 50, 0.001
    # the finite-island closed form: limits
    @test qst_island(N, 0.0; k = 4) == 1.0                              # no migration: complete differentiation
    @test qst_island(N, 3 / 4; k = 4) ≈ 0 atol = 1e-15                  # m = (k−1)/k is panmixia
    @test qst_island(N, 3 / 4; k = 4, among = :demes) ≈ 0 atol = 1e-15
    wright(N, m) = (1 - m)^2 / (2N - (2N - 1) * (1 - m)^2)              # discrete infinite island
    @test qst_island(N, 0.02) ≈ wright(N, 0.02)
    @test qst_island(N, 0.02; k = 10^7) ≈ wright(N, 0.02) rtol = 1e-5   # k → ∞
    @test qst_island(N, 0.02; k = 10^7, among = :demes) ≈ wright(N, 0.02) rtol = 1e-5
    for k in (2, 4, 10), m in (1e-4, 1e-3)                              # small m: Takahata 1983, Slatkin 1991 eq 14
        @test qst_island(N, m; k) ≈ 1 / (1 + 4N * m * (k / (k - 1))^2) rtol = 20m
        @test qst_island(N, m; k, among = :demes) ≈ 1 / (1 + 4N * m * k / (k - 1)) rtol = 20m
    end
    @test issorted([qst_island(N, m; k = 4) for m in 0:0.05:0.75]; rev = true)
    @test_throws ArgumentError qst_island(N, 0.1; k = 1)
    @test_throws ArgumentError qst_island(N, 0.1; among = :all)

    # the backward matrix of the kernel's migration: island model from one visiting order
    P = backward_migration(4, 0.05, N)
    @test P ≈ QTLSim.migration_mixing(island_migration(4, 0.05), fill(50.0, 4), trues(4)) rtol = 1e-12
    @test all(sum(P; dims = 2) .≈ 1) && P[1, 2] ≈ (1 - P[1, 1]) / 3
    @test 1 - backward_migration(4, 0.001, N)[1, 1] ≈ 0.002 rtol = 0.01  # m ≈ 2·rate: exchanges
    m = 1 - P[1, 1]

    # the exact recursion reproduces the closed forms (pool law)
    v = neutral_deme_variances(P, N, Vm)
    iv = island_variances(Vm, N, m, 4)
    @test v.qst ≈ qst_island(N, m; k = 4) rtol = 1e-9
    @test v.qst_demes ≈ qst_island(N, m; k = 4, among = :demes) rtol = 1e-9
    @test all(v.within .≈ iv.within) && iv.within ≈ 2N * 4 * Vm       # Lande 1992: V_W = 2NkV_m whatever m
    @test v.between ≈ iv.between rtol = 1e-9
    @test v.among ≈ iv.between + iv.within / N rtol = 1e-9
    @test iv.total_among ≈ 4N * 4 * iv.qst * Vm / (1 - iv.qst) rtol = 1e-9   # Whitlock 1999: 4NdF V_m/(1 − F)
    c = deme_coalescence(P, N)
    @test all(c.W .≈ 1 .+ [c.S[i, i] for i in 1:4])                                  # two distinct parents, no selfing
    # one deme: Lynch & Hill's 2N V_m
    @test neutral_deme_variances(ones(1, 1), N, Vm).within[1] ≈ 2N * Vm
    # transient: zero at t = 0 (genomes start equal), equilibrium as t → ∞, one pass for many t
    @test all(iszero, neutral_deme_variances(P, N, Vm; t = 0).within)
    @test neutral_deme_variances(P, N, Vm; t = 20_000).among ≈ v.among rtol = 1e-6
    vs = neutral_deme_variances(P, N, Vm; t = [10, 100, 1000])
    @test vs[2].among == neutral_deme_variances(P, N, Vm; t = 100).among && vs[1].among < vs[2].among < vs[3].among
    # panmixia (P = 1/k everywhere): no between-deme variance
    @test neutral_deme_variances(fill(0.25, 4, 4), N, Vm).between ≈ 0 atol = 1e-9

    # exchanges keep deme sizes, so two distinct adults share a deme with probability
    # (N − 1)/(kN − 1); the first-order exchange law gets V_W = 2NV_m(kN − 1)/(N − 1) nearly
    J = exchange_origins(P, N)
    @test J * ones(16) ≈ ones(16) && all(>=(0), J)
    @test neutral_deme_variances(P, N, Vm; origins = :exchange).within[1] ≈ 2N * Vm * (4N - 1) / (N - 1) rtol = 0.005
    @test_throws ArgumentError deme_coalescence(P, N; origins = :other)

    # the ring: differentiation grows with distance; Slatkin's approximation at small m
    Pr = backward_migration(stepping_stone(8, 0.005), N)
    vr = neutral_deme_variances(Pr, N, Vm)
    @test issorted([vr.pairs[1, 1 + i] for i in 1:4])
    @test vr.pairs[1, 2] ≈ vr.pairs[3, 4] rtol = 1e-9
    @test vr.within[1] ≈ 2N * 8 * Vm rtol = 1e-6
    mr = 1 - Pr[1, 1]
    for i in 1:4
        σ2 = (vr.pairs[1, 1 + i] - vr.within[1] / N) / 2
        @test σ2 / (σ2 + 2vr.within[1]) ≈ fst_stepping_stone(N, mr, 8, i) rtol = 0.1
    end
    @test fst_stepping_stone(N, 0.01, 8, 0) == 0

    # a fast simulation: island, k = 4, N = 20, against the exact transient expectation
    p = Params(PopSize = 20, NoQTL = 10, Ploidy = 2, NoGenerations = 400, MutationScalar = 0.002,
               MutationRatePerGamete = 0.1, Drift = true)
    L = lifecycle(); g = L.generation
    cycle = g ⊗ g ⊗ g ⊗ g ⋅ L.migrate(4)
    M = island_migration(4, 0.02)
    census = 200:10:399
    stats = map(1:40) do seed
        o = run_composed(p ⊗ demes(4), MersenneTwister(seed); cycle, migration_rate = M)
        # phenotypic records: subtract the environmental part, V_E per individual, V_E/N per deme mean
        within = mean(o.pheno_var[o.gen .∈ Ref(census)]) - envVar(p)
        among = mean(var(o.pheno_mean[o.gen .== t]) for t in census) - envVar(p) / 20
        (within, among)
    end
    th = neutral_deme_variances(backward_migration(M, 20), 20, muteVar(p); t = collect(census), origins = :exchange)
    for (j, pred) in ((1, mean(mean(x.within) for x in th)), (2, mean(x.among for x in th)))
        x = getindex.(stats, j)
        @test abs(mean(x) - pred) < 5 * std(x) / sqrt(length(x))
    end
end

@testset "fluctuating optima: red noise and periodic moves" begin
    p = Params(PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=4000, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    # parameters, defaults and validation
    @test (p.MoveNoiseSD, p.MoveAutocorr, p.MoveAmplitude, p.MovePeriod) == (0.0, 0.0, 0.0, 100.0)
    @test checkparams(Params(p; MoveType = :red_noise, MoveNoiseSD = 0.1, MoveAutocorr = 0.9)) isa Params
    @test checkparams(Params(p; MoveType = :periodic, MoveAmplitude = -0.1, MovePeriod = 3.5)) isa Params
    for bad in ((MoveAutocorr = 1.0,), (MoveAutocorr = -0.1,), (MoveNoiseSD = -1.0,), (MovePeriod = 0.0,),
                (MovePeriod = -5.0,), (MoveAmplitude = Inf,), (MoveType = :sinusoidal,))
        @test_throws ArgumentError World(Params(p; bad...))
    end
    # QTSim has neither move: writing one is an error, not a silent change of model
    for t in (:red_noise, :periodic)
        @test_throws ArgumentError write_qtsim_options(IOBuffer(), Params(p; MoveType = t))
    end
    # the new fields leave QTSim's moves alone: same seed, same run
    pc = Params(p; NoGenerations = 600, MoveFitnessSDGen = 300, MoveFitnessSD = 0.001, MoveType = :continuous)
    @test simulate(pc; seed = 4).pheno_mean ==
          simulate(Params(pc; MoveNoiseSD = 0.3, MoveAutocorr = 0.5, MoveAmplitude = 1.0, MovePeriod = 7.0); seed = 4).pheno_mean

    # θ trajectories in the kernel
    t0 = 300
    pp = Params(pc; MoveType = :periodic, MoveFitnessSD = 0.002, MoveAmplitude = 0.05, MovePeriod = 40.0)
    r = simulate(pp; seed = 1)
    t = r.gen .- t0
    @test all(r.theta[t .< 0] .== 0)
    @test r.theta[t .>= 0] ≈ 0.002 .* t[t .>= 0] .+ 0.05 .* sin.(2π .* t[t .>= 0] ./ 40)
    @test r.theta[t .>= 0] == fluctuating_optimum.(Ref(pp), 0.0, t[t .>= 0])
    @test simulate(pp; seed = 1) == r                                # no randomness in the move
    pr = Params(pc; MoveType = :red_noise, MoveFitnessSD = 0.001, MoveNoiseSD = 0.05, MoveAutocorr = 0.7)
    r1 = simulate(pr; seed = 2); r2 = simulate(pr; seed = 2)
    @test r1 == r2 && r1.theta != simulate(pr; seed = 3).theta     # deterministic by seed
    @test all(r1.theta[r1.gen .< t0] .== 0)
    rz = simulate(Params(pr; MoveNoiseSD = 0.0); seed = 2)          # σ = 0: the trend alone
    @test rz.theta[rz.gen .>= t0] ≈ 0.001 .* (rz.gen[rz.gen .>= t0] .- t0)
    # the noise is a stationary AR(1) process with SD σ and autocorrelation ρ
    for ρ in (0.0, 0.8)
        q = Params(p; MoveType = :red_noise, MoveFitnessSDGen = 0, MoveFitnessSD = 0.0, MoveNoiseSD = 2.0, MoveAutocorr = ρ)
        w = World(q); rng = Xoshiro(11)
        ε = map(0:199_999) do g
            w.gen = g; moveoptimum!(w, rng); w.θ
        end
        @test std(ε) ≈ 2.0 rtol = 0.03
        @test cor(ε[1:end-1], ε[2:end]) ≈ ρ atol = 0.01
        @test abs(mean(ε)) < 0.1
    end
    @test_throws ArgumentError moveoptimum!(World(Params(pr; MoveFitnessSDGen = 0)), nothing)

    # theory: the linear filter's response to a cycle and to red noise
    Vs = selection_variance(p)
    @test periodic_response(1e7; VG = 0.5, Vs = 1.5).gain ≈ 1 rtol = 1e-6
    @test periodic_response(1e7; VG = 0.5, Vs = 1.5).delay ≈ 1 / 0.25 rtol = 1e-4   # 1/β, the trend's delay
    @test red_noise_lag_var(1.0, 0.0; VG = 0.5, Vs = 1.5) ≈ 2 / (2 - 0.25)
    @test red_noise_lag_var(1.0, 0.999999; VG = 0.5, Vs = 1.5) < 1e-5
    @test expected_lag2(pc; VG = 0.004) ≈ steady_lag(pc; VG = 0.004)^2
    @test lag_load(pc; VG = 0.004) ≈ steady_lag(pc; VG = 0.004)^2 / (2 * (0.3^2 + 0.004 + envVar(p)))
    @test_throws ArgumentError expected_lag2(p; VG = 0.004)

    # the moments lowering: the periodic move lowers exactly, and the mean responds as predicted
    V = multilocus_equilibrium(p).VG          # the default's fixed point, which the move leaves alone
    P, A = 200.0, 0.1
    pm = Params(p; MoveType = :periodic, MoveFitnessSDGen = 500, MoveFitnessSD = 0.0, MoveAmplitude = A, MovePeriod = P)
    o = run_moments(pm; VG = V)
    @test o.theta[o.gen .>= 500] ≈ A .* sin.(2π .* (o.gen[o.gen .>= 500] .- 500) ./ P)
    pr_ = periodic_response(P; VG = V, Vs)
    late = o.gen .>= 2000
    @test o.pheno_mean[late] ≈ pr_.gain * A .* sin.(2π .* (o.gen[late] .- 500) ./ P .- pr_.phase) atol = 1e-8
    @test mean((o.theta[late] .- o.pheno_mean[late]) .^ 2) ≈ expected_lag2(pm; VG = V) rtol = 1e-3
    # red noise lowers to its expected optimum, the trend ...
    pn = Params(p; MoveType = :red_noise, MoveFitnessSDGen = 500, MoveFitnessSD = 1e-4, MoveNoiseSD = 0.05, MoveAutocorr = 0.8)
    on = run_moments(pn; VG = V)
    @test on.theta[on.gen .>= 500] ≈ 1e-4 .* (on.gen[on.gen .>= 500] .- 500)
    # ... or follows a realised optimum, whose mean squared lag is the theory's
    n = 100_000
    θs = zeros(n); rng = Xoshiro(5); ρ, σ = 0.8, 0.05
    θs[1] = σ * randn(rng)
    for i in 2:n; θs[i] = ρ * θs[i-1] + σ * sqrt(1 - ρ^2) * randn(rng); end
    orn = run_moments(Params(pn; NoGenerations = n, MoveFitnessSD = 0.0); VG = V, optimum = θs)
    @test orn.theta == θs
    @test mean((orn.theta[1000:end] .- orn.pheno_mean[1000:end]) .^ 2) ≈
          red_noise_lag_var(σ, ρ; VG = V, Vs) rtol = 0.05
    @test_throws ArgumentError run_moments(pn; optimum = zeros(10))

    # the kernel follows the cycle as the theory says, up to drift
    rk = simulate(Params(pm; NoGenerations = 5000); seed = 8)
    k = rk.gen .>= 1500
    VGk = mean(rk.geno_var[k]); ν = 2π / P; tk = rk.gen[k] .- 500
    X = hcat(sin.(ν .* tk), cos.(ν .* tk), ones(length(tk)))
    a, b, _ = X \ rk.pheno_mean[k]
    th = periodic_response(P; VG = VGk, Vs)
    @test hypot(a, b) / A ≈ th.gain rtol = 0.15
    @test atan(-b, a) ≈ th.phase atol = 0.15

    # the Eidos lowering: the new moves are in the generated model, and only there
    g = lifecycle().generation
    for (t, marker) in ((:red_noise, "moveAutocorr * p1.getValue(\"EPS\")"), (:periodic, "sin(2 * PI * t / movePeriod)"))
        s = slim_script(g, Params(pm; MoveType = t))
        @test occursin(marker, s) && occursin("p1.setValue(\"THETA0\", THETA)", s)
        @test occursin("defineConstant(\"moveType\", \"$t\")", s)
        for c in ("moveNoiseSD", "moveAutocorr", "moveAmplitude", "movePeriod")
            @test occursin("defineConstant(\"$c\"", s)
        end
    end
    sc = slim_script(g, pc)
    @test !occursin("THETA0", sc) && !occursin("moveNoiseSD", sc)
    @test sc == slim_script(g, Params(pc; MoveNoiseSD = 0.3, MovePeriod = 7.0))
    @test !occursin("THETA0", lower(g, :eidos).code)
end

@testset "Params(p; kw...) copies with overrides" begin
    p = Params(NoQTL = 10, Seed = 3, MoveType = :continuous)
    q = Params(p; NoGenerations = 5, MoveFitnessSD = 0.2)
    @test q.NoGenerations == 5 && q.MoveFitnessSD == 0.2
    @test q.NoQTL == 10 && q.Seed == 3 && q.MoveType === :continuous
    @test Params(p) == p
    @test_throws ArgumentError Params(p; Nonsense = 1)
end

# ─── the life cycle as syntax, lowered to the kernel and to Eidos ─────────

@testset "the life cycle is a typed morphism" begin
    L = lifecycle()
    Z = L.presentation[:Zygotes]
    @test dom(L.generation) == Z && codom(L.generation) == Z
    @test L.generation == L.develop ⋅ L.select ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅ L.mutate ⋅ L.swap
    @test codom(L.develop) == dom(L.select)                  # Adults
    @test_throws Exception L.select ⋅ L.develop              # stages do not match
end

@testset "the kernel lowering is a functor and drives run_composed" begin
    L = lifecycle()
    p = Params(PopSize=200, NoQTL=5, Ploidy=2, NoGenerations=60, RecombPoissonMean=2.0,
               PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.002, MutationRatePerGamete=0.01,
               MoveFitnessSDGen=20, MoveFitnessSD=0.5, ChangeFitnessGen=30, ChangeFitnessSD=0.5)
    # image of a composite = composite of the images, on the same state and random numbers
    f, g = L.develop ⋅ L.select, L.reproduce ⋅ L.record ⋅ L.rescale ⋅ L.mutate ⋅ L.swap
    w1 = World(p); w2 = World(p)
    for w in (w1, w2), i in 1:p.PopSize; w.env[i] = 0.1 * i / p.PopSize; end
    lower(f ⋅ g, :kernel)(w1, MersenneTwister(4))
    r = MersenneTwister(4); lower(f, :kernel)(w2, r); lower(g, :kernel)(w2, r)
    @test genomes(w1) == genomes(w2) && w1.env == w2.env && w1.θ == w2.θ && w1.gen == w2.gen == 1
    # the lowered life cycle is run!, with the width schedule and through extinction
    a = run!(World(p), MersenneTwister(8)); b = run_composed(World(p), MersenneTwister(8))
    @test all(getfield(a, k) == getfield(b, k) for k in keys(a))
    @test b.omega[end] == 0.5
    pe = Params(PopSize=200, NoQTL=1, Ploidy=2, NoGenerations=500, PhenotypeVar=0.02, Heritability=0.5,
                MutationScalar=0.001, InitFitnessSD=1e-4, MoveFitnessSDGen=0, MoveFitnessSD=5.0)
    a = run!(World(pe), MersenneTwister(1)); b = run_composed(World(pe), MersenneTwister(1))
    @test a.extinct && b.extinct && a.extinct_gen == b.extinct_gen && a.gen == b.gen
    @test_throws ArgumentError run_composed(World(p), MersenneTwister(1); cycle = L.develop)
end

@testset "the Eidos lowering is a functor and slim_script builds a SLiM model" begin
    L = lifecycle()
    steps = (:develop, :select, :reproduce, :record, :rescale, :mutate, :swap)
    code(f) = lower(f, :eidos).code
    @test code(L.generation) == join(code(getfield(L, s)) for s in steps)
    @test code(L.develop ⋅ L.select) == code(L.develop) * code(L.select)
    p = Params(PopSize=500, NoQTL=1, Ploidy=2, NoGenerations=123, PhenotypeVar=0.03, MoveType=:continuous)
    s = slim_script(L.generation, p)
    @test occursin("initializeSLiMModelType(\"nonWF\")", s) && occursin("reproduction() {", s)
    @test occursin("defineConstant(\"N\", 500)", s) && occursin("defineConstant(\"NGEN\", 123)", s)
    @test occursin("defineConstant(\"moveType\", \"continuous\")", s)
    @test issorted([findfirst("// $st:", s).start for st in steps])
    @test slim_script(L.generation, p) == s
    @test_throws ArgumentError slim_script(L.generation, Params(p; NoQTL = 10, CrossoverModel = :all))
    @test_throws ArgumentError slim_script(L.generation, Params(p; Ploidy = 1))
    @test_throws ArgumentError slim_script(L.develop ⋅ L.select, p)
end

@testset "slim_script still writes the committed SLiM model (validate/slim)" begin
    # validate/slim/check_generated.jl runs this model against the hand-written port
    # qtsim.slim, run for run, with SLiM; this golden-file test catches any change to the
    # Eidos images without SLiM. After an intended change, regenerate the file with
    # `check_generated.jl --write` and rerun the check.
    # The committed models: one locus (scenario B), 100 loci, and three demes with migration.
    dir = joinpath(@__DIR__, "..", "validate", "slim")
    include(joinpath(dir, "golden.jl"))
    for (name, script) in golden_models()
        @test script == read(joinpath(dir, name), String)
    end
end

@testset "a new step lowers through its images" begin
    L = lifecycle()
    A = codom(L.develop)
    fluctuate = Hom(:fluctuate, A, A)
    cycle = L.develop ⋅ fluctuate ⋅ L.select ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅ L.mutate ⋅ L.swap
    @test dom(cycle) == codom(cycle)
    @test_throws ArgumentError lower(cycle, :kernel)
    kick = (w, rng) -> (w.θ = 0.05 * randn(rng); w)
    p = Params(PopSize=200, NoQTL=1, Ploidy=2, NoGenerations=40, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    o = run_composed(World(p), MersenneTwister(2); cycle, images = Dict(:fluctuate => kick))
    @test length(o.gen) == 40 && length(unique(o.theta)) == 40
    s = slim_script(cycle, p; images = Dict(:fluctuate => "    @POP@.setValue(\"THETA\", rnorm(1, 0.0, 0.05));\n"))
    @test occursin("rnorm(1, 0.0, 0.05)", s)
    @test findfirst("rnorm(1, 0.0, 0.05)", s).start < findfirst("// select:", s).start
end

@testset "the life cycle on several demes: tensor products, migration and braidings" begin
    L = lifecycle(); g = L.generation; Z = dom(g)
    meta = g ⊗ g ⊗ g ⋅ L.migrate(3)
    @test dom(meta) == codom(meta) == Z ⊗ Z ⊗ Z
    @test_throws ArgumentError L.migrate(1)
    p = Params(PopSize=100, NoQTL=1, Ploidy=2, NoGenerations=30, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    init!(s, rng) = foreach(w -> (w.env .= QTLSim.envsd(p) .* randn(rng, p.PopSize)), s.demes)
    # the image of g ⊗ g ⊗ g acts deme by deme, as the generation operator does on a metapopulation
    s1 = p ⊗ demes(3); init!(s1, MersenneTwister(1)); s2 = deepcopy(s1)
    lower(g ⊗ g ⊗ g, :kernel).act(s1, [1, 2, 3], MersenneTwister(9))
    operators().generation(s2, MersenneTwister(9))
    @test all(genomes(s1.demes[i]) == genomes(s2.demes[i]) for i in 1:3)
    # a braiding relabels the demes and moves no data
    b = lower((g ⊗ g) ⋅ Catlab.braid(Z, Z), :kernel)
    @test b.perm == [2, 1]
    s3 = p ⊗ demes(2); init!(s3, MersenneTwister(2)); s4 = deepcopy(s3)
    b.act(s3, [1, 2], MersenneTwister(4)); lower(g ⊗ g, :kernel).act(s4, [1, 2], MersenneTwister(4))
    @test genomes(s3.demes[1]) == genomes(s4.demes[1]) && genomes(s3.demes[2]) == genomes(s4.demes[2])
    # runs: a row per deme and sampled generation; migration conserves every deme's size
    o = run_composed(p ⊗ demes(3), MersenneTwister(5); cycle = meta, migration_rate = 0.2)
    @test length(o.gen) == 3 * p.NoGenerations && sort(unique(o.deme)) == [1, 2, 3]
    @test !o.extinct && o.extinct_gen == [-1, -1, -1]
    @test_throws ArgumentError run_composed(p ⊗ demes(2), MersenneTwister(5); cycle = meta)
    o2 = run_composed(p ⊗ demes(2), MersenneTwister(6); cycle = (g ⊗ g) ⋅ Catlab.braid(Z, Z))
    @test length(o2.gen) == 2 * p.NoGenerations
end

@testset "metapopulations in which a deme dies, and error paths" begin
    L = lifecycle(); g = L.generation; Z = dom(g)
    p = Params(PopSize=100, NoQTL=1, Ploidy=2, NoGenerations=40, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    doomed = Params(p; InitFitnessSD=1e-3, MoveFitnessSDGen=0, MoveFitnessSD=20.0)
    for cycle in (g ⊗ g, g ⊗ g ⋅ L.migrate(2))
        s = Ext.Stratified([World(p), World(doomed)])
        o = run_composed(s, MersenneTwister(2); cycle, migration_rate = 0.1)
        @test !o.extinct && o.extinct_gen[1] == -1 && 0 <= o.extinct_gen[2] < p.NoGenerations
        @test count(==(1), o.deme) == p.NoGenerations && count(==(2), o.deme) == o.extinct_gen[2]
        @test totalpop(s) == 200
    end
    # a World ⊗ demes copies the World into independent demes
    w = World(p); s = w ⊗ demes(2)
    @test ndemes(s) == 2 && s.demes[1] !== s.demes[2] && genomes(s.demes[1]) == genomes(w)
    # migration needs demes with the same loci and ploidy, and checks before moving anyone
    s = Ext.Stratified([World(p), World(Params(p; NoQTL = 2))])
    g1 = genomes(s.demes[1])
    @test_throws ArgumentError migrate!(s, MersenneTwister(1); rate = 0.5)
    @test genomes(s.demes[1]) == g1
    @test_throws ArgumentError lower(g, :bogus)
    @test_throws ArgumentError Catlab.compose(lower(g ⊗ g, :kernel), lower(g ⊗ g ⊗ g, :kernel))
    @test_throws ArgumentError Catlab.compose(lower(g ⊗ g, :eidos), lower(g ⊗ g ⊗ g, :eidos))
    @test_throws ArgumentError run_composed(World(p), MersenneTwister(1); cycle = L.develop)
    # a user operator on the whole metapopulation
    seen = Ref(0)
    op = operator(:census, (s, rng) -> (seen[] += totalpop(s); s); scope = :meta)
    op(p ⊗ demes(3), MersenneTwister(1))
    @test seen[] == 300
end

@testset "without Catlab, the categorical functions say how to get it" begin
    code = "using QTLSim; try lifecycle() catch e; print(sprint(showerror, e)) end"
    out = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`, String)
    @test occursin("Catlab", out)
end

@testset "Eidos for several demes and many loci" begin
    L = lifecycle(); g = L.generation
    Ext_ = Base.get_extension(QTLSim, :QTLSimCatlabExt)
    p = Params(PopSize=100, NoQTL=1, Ploidy=2, NoGenerations=30)
    # the image of a tensor product runs each factor's code in its own deme block
    @test lower(g ⊗ g, :eidos).emit(["p1", "p2"]) ==
          Ext_.deme_block(lower(g, :eidos).code, "p1") * Ext_.deme_block(lower(g, :eidos).code, "p2")
    s = slim_script(g ⊗ g ⊗ g ⋅ L.migrate(3), p; migration_rate = 0.05)
    @test occursin("defineConstant(\"DEMES\", 3)", s) && occursin("defineConstant(\"migRate\", 0.05)", s)
    @test all(occursin("do {   // deme p$i", s) for i in 1:3) && occursin("// migrate:", s)
    @test findfirst("// migrate:", s).start > findfirst("do {   // deme p3", s).start
    @test_throws ArgumentError slim_script((g ⊗ g) ⋅ Catlab.braid(dom(g), dom(g)), p)
    sm = slim_script(g, Params(p; NoQTL = 20, RecombPoissonMean = 3.0))
    @test occursin("defineConstant(\"L\", 20)", sm) && occursin("defineConstant(\"RECOMB\", 3.0)", sm)
end

@testset "steps with parameters: syntax" begin
    L = lifecycle(); g = L.generation; A, S = codom(L.develop), codom(L.select)
    s3 = withparams(L.select; ω = 0.3)
    @test dom(s3) == A && codom(s3) == S
    @test s3 == withparams(L.select; ω = 0.3) && s3 != L.select && s3 != withparams(L.select; ω = 0.4)
    @test withparams(L.select; ω = 1) == withparams(L.select; ω = 1.0)              # types normalised
    @test sprint(show, s3) == "select(ω = 0.3)"
    @test withparams(s3; θ = 0.1) == withparams(L.select; ω = 0.3, θ = 0.1)        # parameters merge
    @test withparams(g) === g
    # on a composite, each step takes the parameters it has
    c = withparams(g; ω = 0.3, MutationScalar = 0.01)
    @test dom(c) == codom(c) == dom(g)
    @test c == L.develop ⋅ s3 ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅ withparams(L.mutate; MutationScalar = 0.01) ⋅ L.swap
    @test_throws Exception s3 ⋅ L.develop                                 # still typed
    @test_throws ArgumentError withparams(L.develop; ω = 0.3)             # develop takes no ω
    @test_throws ArgumentError withparams(g; omega = 0.3)                 # no step takes it
    @test_throws ArgumentError withparams(L.select; ω = -1.0)
    @test_throws ArgumentError withparams(L.select; MoveType = :bogus)    # Params rejects it
    @test_throws ArgumentError withparams(L.select; θ = 0.1, MoveFitnessSD = 0.1)
    @test_throws ArgumentError withparams(L.reproduce; PopSize = 10)
    m = L.migrate(3; rate = stepping_stone(3, 0.1))
    @test dom(m) == dom(g ⊗ g ⊗ g) && m == withparams(L.migrate(3); rate = stepping_stone(3, 0.1))
    @test_throws ArgumentError L.migrate(2; rate = 1.5)
    @test_throws ArgumentError L.migrate(2; rate = ones(2, 2))
    @test_throws ArgumentError L.migrate(3; rate = zeros(2, 2))
end

@testset "steps with parameters: kernel and moments lowerings" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize=200, NoQTL=5, Ploidy=2, NoGenerations=120, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.002, MutationRatePerGamete=0.01, InitFitnessSD=0.3, RecombPoissonMean=2.0)
    same(a, b; but = ()) = all(getfield(a, k) == getfield(b, k) for k in keys(a) if !(k in but))
    # without parameters, and with the deme's own values, the run is run!'s, bit for bit
    a = run!(World(p), MersenneTwister(1))
    @test same(a, run_composed(World(p), MersenneTwister(1)))
    @test same(a, run_composed(World(p), MersenneTwister(1); cycle = withparams(g; ω = 0.3)))
    # a step's fields of Params give the run of Params with those values
    over = (MoveType = :continuous, MoveFitnessSDGen = 40, MoveFitnessSD = 0.002, RecombPoissonMean = 0.5,
            Sample = 3, ChangeFitnessGen = 60, ChangeFitnessSD = 0.5, MutationScalar = 0.004)
    @test same(run!(World(Params(p; over...)), MersenneTwister(2)),
               run_composed(World(p), MersenneTwister(2); cycle = withparams(g; over...)))
    # a step's ω is the width for that step: the run of a deme of that width, but the
    # recorded width is the deme's
    b = run_composed(World(Params(p; InitFitnessSD = 0.1)), MersenneTwister(3); cycle = withparams(g; ω = 0.3))
    @test same(run!(World(p), MersenneTwister(3)), b; but = (:omega,)) && all(==(0.1), b.omega)
    # a step's θ: the optimum it selects toward, the deme's recorded optimum unchanged
    pθ = Params(p; NoQTL = 1, PopSize = 500, NoGenerations = 600, MutationRatePerGamete = 0.5, MutationScalar = 0.02)
    o = run_composed(World(pθ), MersenneTwister(4); cycle = withparams(g; θ = 0.1))
    @test abs(mean(o.pheno_mean[301:end]) - 0.1) < 0.02 && all(iszero, o.theta)
    # a cycle of two generations, with two selection steps of different widths
    two = withparams(g; ω = 0.2) ⋅ withparams(g; ω = 1.0)
    o2 = run_composed(World(Params(p; NoGenerations = 10)), MersenneTwister(5); cycle = two)
    @test o2.gen == 0:9
    @test same(run!(World(Params(p; NoGenerations = 10)), MersenneTwister(5)),
               run_composed(World(Params(p; NoGenerations = 10)), MersenneTwister(5); cycle = g ⋅ g))
    # migration at the step's rate is migration at the lowering's rate
    cyc = g ⊗ g ⊗ g
    r1 = run_composed(p ⊗ demes(3), MersenneTwister(6); cycle = cyc ⋅ L.migrate(3; rate = 0.2))
    r2 = run_composed(p ⊗ demes(3), MersenneTwister(6); cycle = cyc ⋅ L.migrate(3), migration_rate = 0.2)
    @test same(r1, r2)
    M = stepping_stone(3, 0.2)
    @test same(run_composed(p ⊗ demes(3), MersenneTwister(7); cycle = cyc ⋅ L.migrate(3; rate = M)),
               run_composed(p ⊗ demes(3), MersenneTwister(7); cycle = cyc ⋅ L.migrate(3), migration_rate = M))
    # the moments lowering: a step's parameters, as Params with those values
    pm = Params(pθ; NoGenerations = 300)
    mo = run_moments(pm; cycle = withparams(g; over...))
    @test same(run_moments(Params(pm; over...)), mo)
    @test same(run_moments(Params(pm; InitFitnessSD = 0.2)), run_moments(pm; cycle = withparams(g; ω = 0.2)); but = (:omega,))
    @test run_moments(pm; cycle = withparams(g; θ = 0.1)).pheno_mean[end] ≈ 0.1 rtol = 1e-4
    # a user's image for one parameterised step, keyed by the generator's name
    kick = withparams(L.select; ω = 0.5)
    seen = Ref(0)
    img = Dict(first(Catlab.args(kick)) => (w, rng) -> (seen[] += 1; QTLSim.select!(w, rng)))
    run_composed(World(Params(p; NoGenerations = 7)), MersenneTwister(1);
                 cycle = L.develop ⋅ kick ⋅ L.reproduce ⋅ L.mutate ⋅ L.swap, images = img)
    @test seen[] == 7
end

@testset "steps with parameters: the Eidos lowering" begin
    L = lifecycle(); g = L.generation
    code(f) = lower(f, :eidos).code
    plain = code(L.select)
    s = code(withparams(L.select; ω = 0.3))
    @test occursin("OMEGA = 0.3;", s) && !occursin("OMEGA = @POP@.getValue(\"OMEGA\");", s)
    s = code(withparams(L.select; θ = 0.1))
    @test occursin("THETA = 0.1;", s) && occursin("if (-1 >= 0 & gen >= -1)", s)
    s = code(withparams(L.select; MoveType = :continuous, MoveFitnessSD = 0.01, MoveFitnessSDGen = 100))
    @test occursin("if (100 >= 0 & gen >= 100)", s) && occursin("THETA = THETA + 0.01;", s)
    @test !occursin("red_noise", s)
    # the bare image of a step with a fluctuating move has the fluctuating moves
    s = code(withparams(L.select; MoveType = :periodic, MoveAmplitude = 0.1, MovePeriod = 50.0))
    @test occursin("\"periodic\" == \"periodic\"", s) && occursin("0.1 * sin(2 * PI * t / 50.0)", s)
    @test !occursin("red_noise", plain)                                   # the plain step is QTSim's
    s = code(withparams(L.mutate; MutationScalar = 0.01))
    @test occursin("((0.01 / (2 * L * rate)) * VE)", s) && !occursin("k0 * MUTVAR", s)
    @test occursin("// mutate: Poisson(L·rate) mutations", s)              # comments untouched
    @test occursin("RECOMB", code(L.reproduce)) && occursin("rpois(2 * NOFF, 3.0)", code(withparams(L.reproduce; RecombPoissonMean = 3.0)))
    @test occursin("gen % 10 == 0", code(withparams(L.record; Sample = 10)))
    @test occursin("if (gen == 50) @POP@.setValue(\"OMEGA\", 0.5);",
                   code(withparams(L.rescale; ChangeFitnessGen = 50, ChangeFitnessSD = 0.5)))
    @test_throws ArgumentError lower(withparams(L.reproduce; CrossoverModel = :all), :eidos)
    # migration at the step's rate: a literal, or a local copy of the matrix
    mig(f) = lower(f, :eidos).emit(["p1", "p2", "p3"])
    @test occursin("runif(sizes[d]) < 0.05", mig(L.migrate(3; rate = 0.05)))
    m = mig(L.migrate(3; rate = stepping_stone(3, 0.02; circular = false)))
    @test occursin("MIGstep = c(0.0, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0);", m) && occursin("cumSum(MIGstep[", m)
    # slim_script: a step's fluctuating move defines the constants it reads
    p = Params(PopSize=100, NoQTL=1, Ploidy=2, NoGenerations=30, InitFitnessSD = 0.3)
    sc = slim_script(withparams(g; MoveType = :periodic, MoveFitnessSDGen = 10, MoveAmplitude = 0.1), p)
    @test occursin("defineConstant(\"movePeriod\", 100.0)", sc) && occursin("0.1 * sin(2 * PI * t / movePeriod)", sc)
    @test !occursin("movePeriod", slim_script(g, p))
    # fecundity selection: spliced into select and reproduce only when a deme or a step uses it
    @test !occursin("selection", slim_script(g, p)) && !occursin("weighted", slim_script(g, p))
    sc = slim_script(g, Params(p; Selection = :fecundity))
    @test occursin("defineConstant(\"selection\", \"fecundity\")", sc)
    @test occursin("i1 = sample(0:(ns - 1), NOFF, T, sw);", sc) && occursin("growthRate * nw /", sc)
    @test occursin("i1 = rdunif(NOFF, 0, ns - 1);", sc)                 # viability, as before
    s = code(withparams(L.select; Selection = :fecundity))
    @test occursin("weighted = (\"fecundity\" == \"fecundity\");", s)
    @test occursin("sample(0:(ns - 1), NOFF, T, sw)", code(L.develop ⋅ withparams(L.select; Selection = :fecundity) ⋅ L.reproduce))
    sc = slim_script(g ⊗ g ⋅ L.migrate(2), [p, Params(p; Selection = :fecundity)])
    @test occursin("defineConstant(\"selection\", c(\"viability\", \"fecundity\"))", sc) &&
          occursin("p1.getValue(\"selection\") == \"fecundity\"", sc)
    # one generation of each deme per tick
    @test_throws ArgumentError slim_script(g ⋅ g, p)
    @test_throws ArgumentError slim_script((g ⋅ g) ⊗ g ⋅ L.migrate(2), p)
end

@testset "per-deme parameters in slim_script" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize=100, NoQTL=1, Ploidy=2, NoGenerations=30, InitFitnessSD = 0.3)
    cycle = g ⊗ g ⊗ g ⋅ L.migrate(3)
    # the same Params in every deme: the model of one Params
    @test slim_script(cycle, [p, p, p]; migration_rate = 0.01) == slim_script(cycle, p; migration_rate = 0.01)
    @test slim_script(g, [p]) == slim_script(g, p)
    ps = [p, Params(p; InitFitnessSD = 0.6, PopSize = 200, MoveType = :periodic, MoveFitnessSDGen = 5),
          Params(p; Heritability = 0.4, MutationScalar = 0.002)]
    sc = slim_script(cycle, ps)
    @test occursin("defineConstant(\"N\", c(100, 200, 100))", sc)
    @test occursin("defineConstant(\"initOmega\", c(0.3, 0.6, 0.3))", sc)
    @test occursin("defineConstant(\"moveType\", c(\"phenotype\", \"periodic\", \"phenotype\"))", sc)
    @test occursin("defineConstant(\"h2\", c(0.5, 0.5, 0.4))", sc) && occursin("defineConstant(\"VP\", 1.0)", sc)
    @test occursin("function (*)perDeme(* x)", sc) && occursin("pop = sim.addSubpop(i + 1, perDeme(N)[i]);", sc)
    @test occursin("pop.setValue(\"OMEGA\", perDeme(OMEGA0)[i]);", sc)
    @test occursin("ifelse(perDeme(initOmega) < 0", sc)
    # each deme's code reads its own values from its subpopulation, the derived ones too
    for (pop, name) in (("p1", "N"), ("p2", "moveType"), ("p3", "ESD"), ("p3", "MUTVAR"))
        @test occursin("$pop.getValue(\"$name\")", sc)
    end
    for name in ("N", "moveType", "moveGen", "h2", "m", "VG", "VE", "ESD", "MUTVAR")
        @test occursin("pop.setValue(\"$name\", perDeme($name)[i]);", sc)
    end
    @test !occursin("p1.getValue(\"RECOMB\")", sc) && !occursin("setValue(\"initOmega\"", sc)
    # the fluctuating moves are in the plain select, since a deme uses one
    @test occursin("p1.getValue(\"moveType\") == \"red_noise\"", sc)
    # a step's parameter wins over the deme's
    sc2 = slim_script(withparams(g; MoveFitnessSDGen = 7) ⊗ g ⊗ g ⋅ L.migrate(3), ps)
    @test occursin("if (7 >= 0 & gen >= 7)", sc2) && occursin("p2.getValue(\"moveGen\") >= 0", sc2)
    @test_throws ArgumentError slim_script(cycle, [p, p])
    @test_throws ArgumentError slim_script(cycle, [p, p, Params(p; NoQTL = 2)])
    @test_throws ArgumentError slim_script(cycle, [p, p, Params(p; Demography = :beverton_holt)])
end

@testset "the moments of several demes" begin
    L = lifecycle(); g = L.generation
    Ext_ = Ext
    # the expected composition after the kernel's migration: rows are distributions, and
    # with two equal demes at a small rate about 2·rate of each is replaced
    P = Ext_.migration_mixing([0 0.01; 0.01 0], [1000.0, 1000.0], [true, true])
    @test all(sum(P; dims = 2) .≈ 1) && all(P .>= 0)
    @test P[1, 2] ≈ 0.02 rtol = 0.02
    P3 = Ext_.migration_mixing([0 0.2 0; 0.1 0 0; 0 0 0], [100.0, 300.0, 50.0], [true, true, true])
    @test all(sum(P3; dims = 2) .≈ 1) && P3[3, :] == [0, 0, 1]
    @test Ext_.migration_mixing([0 0.2; 0.2 0], [100.0, 100.0], [true, false]) == [1 0; 0 1]
    # no migration: each deme is its own one-deme run
    p = Params(PopSize=500, NoQTL=1, Ploidy=2, NoGenerations=400, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    ps = [p, Params(p; InitFitnessSD = 0.6)]
    o = run_moments(ps; cycle = g ⊗ g)
    for d in 1:2
        o1 = run_moments(ps[d])
        @test o.geno_var[o.deme .== d] == o1.geno_var && o.pheno_mean[o.deme .== d] == o1.pheno_mean
    end
    @test run_moments(p; cycle = g ⊗ g).geno_var[1:2:end] == run_moments(p).geno_var
    @test_throws ArgumentError run_moments([p, p]; cycle = g ⊗ g ⊗ g)
    @test_throws ArgumentError run_moments([p, p]; cycle = g ⊗ g, zbar = [0.0, 0.1, 0.2])
    # neutral demes with different means: the moments right after one migration are the
    # kernel's, for unequal sizes and a migration matrix, and the means stay exact after
    pn = Params(p; PopSize = 400, NoGenerations = 4, MutationScalar = 0.0, Drift = true)
    psn = [pn, Params(pn; PopSize = 200), Params(pn; PopSize = 600)]
    z0 = [0.0, 0.2, -0.1]
    cyc = g ⊗ g ⊗ g ⋅ L.migrate(3; rate = [0 0.1 0.05; 0.2 0 0; 0.02 0.3 0])
    mo = run_moments(psn; cycle = cyc, zbar = z0)
    R = 300
    zk = zeros(4, 3); vk = zeros(4, 3)
    for r in 1:R
        s = psn ⊗ demes(3)
        for (d, w) in enumerate(s.demes); setgenomes!(w, fill(z0[d] / 2, w.p.NoQTL, w.p.Ploidy, w.p.PopSize)); end
        ok = run_composed(s, MersenneTwister(r); cycle = cyc)
        for j in eachindex(ok.gen)
            zk[ok.gen[j] + 1, ok.deme[j]] += ok.pheno_mean[j] / R; vk[ok.gen[j] + 1, ok.deme[j]] += ok.geno_var[j] / R
        end
    end
    at(t, d) = findfirst(i -> mo.gen[i] == t && mo.deme[i] == d, eachindex(mo.gen))
    for d in 1:3
        @test zk[2, d] ≈ mo.pheno_mean[at(1, d)] atol = 0.002
        @test vk[2, d] ≈ mo.geno_var[at(1, d)] rtol = 0.03
        @test all(abs(zk[t + 1, d] - mo.pheno_mean[at(t, d)]) < 0.004 for t in 2:3)
    end
    # selection toward different optima: the demes diverge, and migration holds them together,
    # in the kernel (averaged over replicates) and in the moments
    pd = Params(p; NoGenerations = 1000)
    late(o, d) = mean(o.pheno_mean[i] for i in eachindex(o.gen) if o.deme[i] == d && o.gen[i] >= 600)
    function divergence(rate)
        cyc = withparams(g; θ = 0.1) ⊗ withparams(g; θ = -0.1, ω = 0.2) ⋅ L.migrate(2; rate)
        mo = run_moments(pd; cycle = cyc)
        ks = [(o = run_composed(pd ⊗ demes(2), MersenneTwister(s); cycle = cyc); late(o, 1) - late(o, 2)) for s in 1:8]
        late(mo, 1) - late(mo, 2), mean(ks)
    end
    m0, k0 = divergence(0.0)
    m1, k1 = divergence(0.02)
    @test m0 ≈ 0.2 rtol = 1e-4
    @test k0 ≈ 0.2 rtol = 0.05
    @test m1 < 0.6m0 && k1 < 0.6k0
    @test m1 ≈ k1 rtol = 0.25          # the one-locus kernel has more variance than the Gaussian lowering
end

@testset "the moments of several demes with drift among them (between = true)" begin
    L = lifecycle(); g = L.generation
    p = Params(Ploidy = 2, NoQTL = 20, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.002,
               MutationRatePerGamete = 0.05, Drift = true, PopSize = 30, NoGenerations = 300, Sample = 1)
    Vm = muteVar(p)
    relerr(a, b) = abs(a - b) / max(abs(b), 1e-300)
    # neutral: the recursion is deme_coalescence's read forward, for both laws of origins,
    # the island model, the ring and unequal demes with a migration matrix, at every generation
    ts = [1, 2, 7, 50, 299]
    for (ps, M) in ((fill(p, 4), island_migration(4, 0.02)), (fill(p, 5), stepping_stone(5, 0.05)),
                    ([p, Params(p; PopSize = 20), Params(p; PopSize = 45)], [0 0.1 0.05; 0.2 0 0; 0.02 0.3 0]))
        k = length(ps)
        Nv = [q.PopSize for q in ps]
        cyc = foldl(⊗, fill(g, k)) ⋅ L.migrate(k)
        P = backward_migration(M, Nv)
        for origins in (:exchange, :pool)
            mo = run_moments(ps; cycle = cyc, migration_rate = M, between = true, origins)
            th = neutral_deme_variances(P, Nv, Vm; t = ts, origins)
            for (j, t) in enumerate(ts)
                i = findfirst(==(t), mo.metapop.gen)
                rows = findall(==(t), mo.gen)
                @test all(relerr.(mo.geno_var[rows], th[j].within) .< 1e-10)
                @test relerr(mo.metapop.var_among[i], th[j].among) < 1e-10
                @test abs(mo.metapop.between[i] - th[j].between) < 1e-10 * th[j].among    # a small difference early on
                @test abs(mo.metapop.qst[i] - th[j].qst) < 1e-9 && abs(mo.metapop.qst_demes[i] - th[j].qst_demes) < 1e-9
                @test maximum(abs.(mo.metapop.pairs[i] .- th[j].pairs)) < 1e-10 * maximum(th[j].pairs)
            end
        end
        # a law of origins as a matrix: kron(P, P) is :pool
        a = run_moments(ps; cycle = cyc, migration_rate = M, between = true, origins = :pool).metapop.var_among
        b = run_moments(ps; cycle = cyc, migration_rate = M, between = true, origins = kron(P, P)).metapop.var_among
        @test a ≈ b rtol = 1e-12
    end
    # ... and its equilibrium
    pe = Params(p; PopSize = 10, NoGenerations = 6000, Sample = 1000)
    cyc = g ⊗ g ⊗ g ⋅ L.migrate(3)
    mo = run_moments(pe; cycle = cyc, migration_rate = 0.05, between = true)
    eq = neutral_deme_variances(backward_migration(3, 0.05, 10), 10, Vm; origins = :exchange)
    @test mo.metapop.var_among[end] ≈ eq.among rtol = 1e-6
    @test mo.metapop.within[end] ≈ sum(eq.within) / 3 rtol = 1e-6
    @test mo.metapop.qst[end] ≈ eq.qst rtol = 1e-6
    # one deme: within-deme drift exactly; no migration: the demes are independent copies
    m1 = run_moments(p; between = true)
    @test all(relerr(m1.geno_var[t + 1], neutral_deme_variances(ones(1, 1), 30, Vm; t).within[1]) < 1e-10 for t in (5, 299))
    @test m1.metapop.var_among == zeros(300)
    m0 = run_moments(p; cycle = g ⊗ g, between = true)
    @test m0.geno_var[1:2:end] == m1.geno_var
    @test m0.metapop.var_among[end] ≈ neutral_deme_variances([1.0 0.0; 0.0 1.0], 30, Vm; t = 299).among rtol = 1e-10
    # off by default: results and columns unchanged
    ps = fill(Params(p; Drift = false, InitFitnessSD = 2.0), 3)
    cyc = g ⊗ g ⊗ g ⋅ L.migrate(3)
    @test run_moments(ps; cycle = cyc, migration_rate = 0.02) == run_moments(ps; cycle = cyc, migration_rate = 0.02, between = false)
    @test !haskey(run_moments(ps; cycle = cyc, migration_rate = 0.02), :metapop)
    @test_throws ArgumentError run_moments(fill(Params(p; Ploidy = 1), 3); cycle = cyc, between = true)
    @test_throws ArgumentError run_moments(ps; cycle = cyc, between = true, origins = :other)
    # selection: a common optimum holds the demes together (Q_ST below the neutral F_ST),
    # different optima pull them apart (above it); selection lowers the within-deme variance
    pn = Params(p; NoGenerations = 1500, Sample = 100)
    sel = Params(pn; Drift = false, InitFitnessSD = 3.0)
    run3(q, θ) = run_moments(q; cycle = foldl(⊗, [withparams(g; θ = x) for x in θ]) ⋅ L.migrate(3), migration_rate = 0.01,
                             between = true).metapop
    neu, com, loc = run3(pn, zeros(3)), run3(sel, zeros(3)), run3(sel, [-1.0, 0.0, 1.0])
    @test com.qst[end] < neu.qst[end] < loc.qst[end]
    @test com.within[end] < neu.within[end]
    # pair_recombination: none for one locus; for two loci only one crossover point can split them
    @test pair_recombination(1, 2.0) == 0.0
    @test pair_recombination(2, 0.7) ≈ 0.7 * exp(-0.7)
    @test pair_recombination(2, 0.7; model = :all) ≈ 0.7 * exp(-0.7)
    @test pair_recombination(Params(NoQTL = 50, RecombPoissonMean = 1.0)) ≈ 0.2157 atol = 1e-3
end

@testset "sweep and replicate: threads, cache, resumption, errors" begin
    p = Params(PopSize = 100, NoQTL = 1, NoGenerations = 50)
    summ(o) = (VG = o.geno_var[end], zbar = o.pheno_mean[end], ext = o.extinct)
    # results are in job order and do not depend on the number of tasks or the job order
    r1 = replicate(p, 6; seed = 10, summary = summ, ntasks = 1)
    r3 = replicate(p, 6; seed = 10, summary = summ, ntasks = 3, cost = s -> -s)
    @test r1 == r3 == [summ(simulate(p; seed = s)) for s in 10:15]
    @test sweep(x -> x^2, 1:10; ntasks = 4, cost = x -> x % 3) == (1:10) .^ 2
    @test isempty(sweep(identity, Int[]))
    mktempdir() do dir
        csv = joinpath(dir, "runs.csv")
        # a cache round-trip returns the rows, Float64s exactly; strings with commas are quoted
        f(j) = (; j, x = j / 7, flag = isodd(j), name = "a,\"b\" $j")
        rows = sweep(f, 1:5; cache = csv, ntasks = 2)
        @test rows == f.(1:5)
        @test readline(csv) == "id,j,x,flag,name"
        # resumption: done jobs are not rerun, and the rows equal an uninterrupted sweep's
        ran = Int[]; lk = ReentrantLock()
        g(j) = (lock(lk) do; push!(ran, j); end; f(j))
        @test sweep(g, 1:8; cache = csv, ntasks = 3) == f.(1:8)
        @test sort(ran) == 6:8
        # a line cut short by an interruption is dropped and its job rerun
        last_id = parse(Int, split(readlines(csv)[end], ',')[1])
        write(csv, read(csv, String)[1:end-3]); empty!(ran)
        @test sweep(g, 1:8; cache = csv) == f.(1:8) && ran == [last_id]
        @test countlines(csv) == 9
        # the cache is keyed by id; replicate keys by seed
        c2 = joinpath(dir, "rep.csv")
        @test replicate(p, 3; seed = 4, summary = summ, cache = c2) == [summ(simulate(p; seed = s)) for s in 4:6]
        @test first.(split.(readlines(c2)[2:end], ',')) ⊆ ["4", "5", "6"]
        # errors: duplicate ids, rows that do not fit the cache, results that are not rows
        @test_throws ArgumentError sweep(f, [1, 1]; cache = joinpath(dir, "dup.csv"))
        @test_logs (:error,) match_mode = :any @test_throws ArgumentError sweep(j -> (; y = j), 9:9; cache = csv)
        @test_logs (:error,) match_mode = :any @test_throws ArgumentError sweep(j -> j, 1:2; cache = joinpath(dir, "n.csv"))
        # a failing job stops the sweep, rows finished before it stay cached, and its error is rethrown
        c3 = joinpath(dir, "fail.csv")
        h(j) = j == 3 ? error("job 3 fails") : (; j)
        @test_logs (:error,) match_mode = :any @test_throws ErrorException sweep(h, 1:5; cache = c3, ntasks = 1)
        @test readlines(c3) == ["id,j", "1,1", "2,2"]
        @test sweep(j -> (; j), 1:5; cache = c3) == [(; j) for j in 1:5]
    end
    @test_throws ArgumentError sweep(identity, 1:2; ntasks = 0)
end

@testset "several traits: parameters and the kernel" begin
    base = Params(PopSize = 300, NoQTL = 20, Ploidy = 2, NoGenerations = 200, PhenotypeVar = 1.0,
                  Heritability = 0.5, MutationScalar = 0.01, MutationRatePerGamete = 0.01, RecombPoissonMean = 2.0)
    p = Params(base; NoTraits = 2, MutCorr = 0.6, SelCorr = -0.4)
    # derived matrices; one trait is the scalar model
    @test mutational_covariance(p) ≈ muteVar(p) .* [1 0.6; 0.6 1]
    @test selection_covariance(p) ≈ initFitSD(p)^2 .* [1 -0.4; -0.4 1]
    @test selection_covariance(Params(p; SelCorr = 0.0, SelCov = [2.0 0.5; 0.5 1.0])) ≈ [2.0 0.5; 0.5 1.0]
    @test mutational_covariance(base) ≈ fill(muteVar(base), 1, 1)
    @test selection_covariance(base) ≈ fill(initFitSD(base)^2, 1, 1)
    @test move_direction(Params(p; MoveDirection = [3.0, 4.0])) ≈ [0.6, 0.8] && move_direction(p) == [1.0, 0.0]
    @test Params(p) == p && hash(Params(p)) == hash(p) && Params(base; NoTraits = 2) == Params(base; NoTraits = 2)
    # checkparams
    bad(; kw...) = @test_throws ArgumentError World(Params(base; kw...))
    bad(NoTraits = 0); bad(MutCorr = 0.5)                     # trait fields need NoTraits > 1
    bad(NoTraits = 2, SelCorr = 1.0); bad(NoTraits = 3, MutCorr = -0.6)   # not a correlation of 3 traits
    bad(NoTraits = 2, MutCov = [1.0 2.0; 2.0 1.0]); bad(NoTraits = 2, SelCov = ones(3, 3))
    bad(NoTraits = 2, MutCov = [1.0 0.0; 0.0 1.0], MutCorr = 0.1); bad(NoTraits = 2, MoveDirection = [1.0])
    bad(NoTraits = 2, MoveDirection = [0.0, 0.0]); bad(NoTraits = 2, MoveNoiseCov = ones(3, 3))
    bad(NoTraits = 2, StandardisePhenotype = true)
    @test World(Params(base; NoTraits = 3, MutCorr = -0.5)) isa World       # -1/(n-1) is singular but allowed
    @test_throws ArgumentError write_qtsim_options(IOBuffer(), p)
    # genotypes are the sums of the effects, trait by trait; the block pool stays consistent;
    # a generation allocates nothing in the steady state
    w = World(p); rng = MersenneTwister(4); init_env!(w, rng)
    for t in 0:149
        w.gen = t; generation!(w, rng)
    end
    QTLSim.check_blocks(w)
    develop!(w)
    G = genomes(w)
    @test size(G) == (2, p.NoQTL, 2, w.N)
    @test traitvalues(w; what = :genotype) ≈ dropdims(sum(G; dims = (2, 3)); dims = (2, 3)) atol = 1e-12
    @test traitvalues(w) ≈ traitvalues(w; what = :genotype) .+ traitvalues(w; what = :env)
    @test w.genotype[1:w.N] == traitvalues(w; what = :genotype)[1, :]      # the scalar fields hold trait 1
    alloc(w, rng) = @allocated generation!(w, rng)
    alloc(w, rng); @test alloc(w, rng) == 0
    v = World(p); setgenomes!(v, G); develop!(v)
    @test genomes(v) == G
    @test traitvalues(v; what = :genotype) ≈ dropdims(sum(G; dims = (2, 3)); dims = (2, 3)) atol = 1e-12
    @test_throws ArgumentError setgenomes!(v, zeros(p.NoQTL, 2, 3))
    # fitness: QTSim's 0.4 at the optimum; with Ω diagonal, the product of one-trait fitnesses
    q = Params(base; NoTraits = 2, InitFitnessSD = 0.7)
    u = World(q); init_env!(u, MersenneTwister(1)); develop!(u)
    Z = traitvalues(u)
    select!(u, MersenneTwister(2))
    c = 1 / (0.4 * sqrt(2π))
    @test u.fitness[1:u.N] ≈ [c * exp(-(Z[1, i]^2 + Z[2, i]^2) / (2 * 0.7^2)) for i in 1:u.N]
    one = World(Params(base; InitFitnessSD = 0.7)); one.phenotype .= Z[1, :]; select!(one, MersenneTwister(2))
    @test u.fitness[1:u.N] ≈ one.fitness .* [exp(-Z[2, i]^2 / (2 * 0.7^2)) for i in 1:u.N]
    # mutation injects M per individual per generation (Lean mutationalVarianceInput_eq, for
    # every linear combination of the traits)
    pm = Params(base; NoTraits = 2, PopSize = 40000, NoQTL = 10, MutationRatePerGamete = 0.05,
                MutCov = [0.02 -0.012; -0.012 0.01])
    m = World(pm); m.Noff = m.N; r = MersenneTwister(9)
    mutate!(m, r); swapbuffers!(m); develop!(m)
    C = cov(traitvalues(m; what = :genotype); dims = 2)
    @test C ≈ mutational_covariance(pm) rtol = 0.05
    # the records: per-trait means, G and P matrices, lag and selection differential
    o = simulate(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.01,
                        MoveDirection = [1.0, 1.0]); seed = 3)
    @test !o.extinct && length(o.trait_mean) == length(o.gen) == length(o.gmatrix)
    @test o.pheno_mean ≈ first.(o.trait_mean) && o.geno_var ≈ [g[1, 1] for g in o.gmatrix]
    @test all(g -> g ≈ g' && all(>(-1e-12), eigvals(Symmetric(g))), o.gmatrix)
    @test o.theta ≈ [o.trait_mean[i][1] + o.trait_lag[i][1] for i in eachindex(o.gen)]
    @test o.gen[101] == 100 && o.trait_lag[101] ≈ -o.trait_mean[101] .+ 0.01 / sqrt(2) .* [1, 1]   # moved once along (1, 1)/√2
    # the same seed, the same run; the one-trait records have no trait columns
    @test simulate(p; seed = 5).gmatrix == simulate(p; seed = 5).gmatrix
    @test !haskey(simulate(Params(base; NoGenerations = 5); seed = 1), :gmatrix)
end

@testset "several traits: the kernel lowering and metapopulations" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize = 200, NoQTL = 10, Ploidy = 2, NoGenerations = 100, MutationScalar = 0.01,
               MutationRatePerGamete = 0.02, NoTraits = 2, MutCorr = 0.5, SelCorr = 0.3)
    # the life cycle is unchanged: run_composed drives the same morphism over traits
    a = simulate(p; seed = 11)
    b = run_composed(World(p), MersenneTwister(11); cycle = g)
    @test a.gmatrix == b.gmatrix && a.trait_mean == b.trait_mean
    # two demes with migration: the traits move with the migrants
    s = p ⊗ demes(2)
    o = run_composed(s, MersenneTwister(1); cycle = g ⊗ g ⋅ L.migrate(2), migration_rate = 0.05)
    @test !o.extinct && length(o.gmatrix) == length(o.gen) == 2 * p.NoGenerations
    QTLSim.check_blocks(s.demes...)
    w, v = World(p), World(p); init_env!(w, MersenneTwister(1)); init_env!(v, MersenneTwister(2))
    e1, e2 = traitvalues(w; what = :env)[:, 3], traitvalues(v; what = :env)[:, 5]
    Ext.swap_individuals!(w, 3, v, 5)
    @test traitvalues(w; what = :env)[:, 3] == e2 && traitvalues(v; what = :env)[:, 5] == e1
    @test_throws ArgumentError migrate!(Ext.Stratified([World(p), World(Params(p; NoTraits = 3))]), MersenneTwister(1))
end

@testset "several traits: the same ACSet schema, vector-valued attributes" begin
    p = Params(PopSize = 30, NoQTL = 4, Ploidy = 2, NoGenerations = 20, MutationScalar = 0.05,
               MutationRatePerGamete = 0.2, NoTraits = 3, MutCorr = 0.2)
    w = World(p); rng = MersenneTwister(2); init_env!(w, rng)
    for t in 0:19; w.gen = t; generation!(w, rng); end
    develop!(w)
    pop = lift(w)
    @test pop isa Ext.TraitPopulation && Catlab.acset_schema(pop) == Catlab.acset_schema(lift(World(Params(p; NoTraits = 1, MutCorr = 0.0))))
    @test reduce(hcat, genotypes(pop)) ≈ traitvalues(w; what = :genotype)
    @test reduce(hcat, phenotypes(pop)) ≈ traitvalues(w)
    v = World(p); lower!(pop, v)
    @test genomes(v) == genomes(w) && traitvalues(v) ≈ traitvalues(w)
end

@testset "several traits: theory (Lande 1980; Schluter 1996)" begin
    # the closed-form equilibrium solves G (A + G)⁻¹ G = M for M and A that do not commute
    M = [0.004 0.0025; 0.0025 0.003]; A = [40.0 -12.0; -12.0 25.0]
    G = gaussian_G(M, A)
    @test G * ((A + G) \ G) ≈ M rtol = 1e-10
    @test all(>(0), eigvals(Symmetric(G)))
    @test G ≈ gaussian_G_iterate(M, A) rtol = 1e-8
    # commuting M and A: the one-trait formula on each shared axis
    Q = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    Mc, Ac = Q * [0.004 0; 0 0.001] * Q', Q * [30.0 0; 0 60.0] * Q'
    one(m, a) = (m + sqrt(m^2 + 4m * a)) / 2
    @test gaussian_G(Mc, Ac) ≈ Q * [one(0.004, 30.0) 0; 0 one(0.001, 60.0)] * Q' rtol = 1e-10
    # one trait: gaussian_vg_exact
    p1 = Params(PhenotypeVar = 0.02, Heritability = 0.4, MutationScalar = 0.004, InitFitnessSD = 2.0)
    @test gaussian_G(p1)[1, 1] ≈ gaussian_vg_exact(p1) rtol = 1e-12
    @test selection_variance_matrix(p1)[1, 1] ≈ selection_variance(p1)
    # weak mutation: the matrix geometric mean A # M
    Mw = M .* 1e-6
    Ah = QTLSim.symsqrt(A); Aih = inv(Ah)
    @test gaussian_G(Mw, A) ≈ Ah * QTLSim.symsqrt(Aih * Mw * Aih) * Ah rtol = 2e-3
    # the steady lag: the response to it is the rate of change; one trait is steady_lag
    k = [0.01, 0.0]
    λ = steady_lag_vector(k; G, A)
    @test trait_response(λ; G, A) ≈ k
    @test steady_lag_vector([0.01]; G = fill(0.5, 1, 1), A = fill(1.5, 1, 1))[1] ≈ steady_lag(0.01; VG = 0.5, Vs = 1.5)
    pc = Params(p1; NoTraits = 2, MutCorr = 0.5, MoveFitnessSDGen = 10, MoveType = :continuous, MoveFitnessSD = 0.001,
                MoveDirection = [1.0, 1.0])
    @test steady_lag_vector(pc; G) ≈ steady_lag_vector(0.001 .* [1, 1] ./ sqrt(2); G, A = selection_variance_matrix(pc))
    @test_throws ArgumentError steady_lag_vector(Params(pc; MoveType = :phenotype); G)
    # the line of least resistance: the response turns from β toward g_max
    Gd = [1.0 0.8; 0.8 1.0]                 # g_max along (1, 1)/√2
    @test gmax(Gd) ≈ [1, 1] ./ sqrt(2)
    d = response_deflection([1.0, 0.0]; G = Gd, A = [10.0 0; 0 10.0])
    @test d.angle_response_gmax < d.angle_gradient_gmax && d.angle_response_gradient > 0
    @test d.response ≈ Gd * d.β && vector_angle([1.0, 0.0], [0.0, 2.0]) ≈ 90
    for _ in 1:50                            # never away from the g_max axis
        B = randn(3, 3); Gr = B * B'; β = randn(3)
        g = gmax(Gr)
        @test abs(QTLSim.dot(Gr * β, g)) / QTLSim.norm(Gr * β) >= abs(QTLSim.dot(β, g)) / QTLSim.norm(β) - 1e-12
    end
end

@testset "several traits: the moments lowering, a second table for the same syntax" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize = 1000, NoQTL = 50, Ploidy = 2, NoGenerations = 20_000, PhenotypeVar = 1.0, Heritability = 0.5,
               MutationScalar = 0.01, MutationRatePerGamete = 0.01, NoTraits = 2, MutCorr = 0.5, SelCorr = -0.3)
    o = run_moments(p; loci = :multilocus)   # dispatches to the trait table
    @test o.gmatrix[end] ≈ gaussian_G_loci(p) rtol = 1e-6       # the multilocus fixed point
    @test gaussian_G_loci(p) ≈ gaussian_G(100 .* mutational_covariance(p), selection_variance_matrix(p))
    @test run_moments(p; loci = :aggregate).gmatrix[end] ≈ gaussian_G(p) rtol = 1e-6    # the legacy recursion
    @test run_moments(p; erosion = true, loci = :multilocus).gmatrix[end] ≈ gaussian_G_loci(p; drift = true) rtol = 1e-6
    @test run_moments(p).gmatrix[end] ≈ rare_alleles_G(p) rtol = 1e-6         # the default, :rare_local, for one deme
    # one trait of the two-trait recursion is the one-trait recursion
    p1t = Params(p; NoTraits = 1, MutCorr = 0.0, SelCorr = 0.0, NoGenerations = 2000)
    for erosion in (false, true)
        @test Ext.run_trait_moments(p1t; erosion, loci = :multilocus).geno_var[2:end] ≈ run_moments(p1t; erosion, loci = :multilocus).geno_var[2:end] rtol = 1e-10
    end
    @test o.pmatrix[end] ≈ o.gmatrix[end] + envVar(p) * [1 0; 0 1]
    @test all(0 .< o.mean_fitness .<= 1)
    # a continuous move along (1, 2): the lag converges to the steady lag vector
    pc = Params(p; NoGenerations = 30_000, MoveFitnessSDGen = 20_000, MoveType = :continuous, MoveFitnessSD = 0.002,
                MoveDirection = [1.0, 2.0])
    oc = run_moments(pc)
    @test oc.trait_lag[end] ≈ steady_lag_vector(pc; G = oc.gmatrix[end]) rtol = 1e-6
    @test oc.trait_mean[end] - oc.trait_mean[end-1] ≈ 0.002 .* QTLSim.move_direction(pc) rtol = 1e-6
    # one trait: the multivariate recursion with 1×1 matrices is the one-trait table's
    p1 = Params(PopSize = 1000, NoGenerations = 3000, PhenotypeVar = 0.02, MutationScalar = 0.02,
                MutationRatePerGamete = 0.5, InitFitnessSD = 0.3, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.001,
                MoveType = :continuous, ChangeFitnessGen = 2000, ChangeFitnessSD = 0.4)
    a, b = run_moments(p1; loci = :multilocus), Ext.run_trait_moments(p1; loci = :multilocus)
    @test a.geno_var ≈ b.geno_var rtol = 1e-10
    @test a.pheno_mean ≈ b.pheno_mean rtol = 1e-10
    @test a.theta ≈ b.theta rtol = 1e-10
    @test a.mean_fitness ≈ b.mean_fitness rtol = 1e-10
    @test a.sel_diff ≈ b.sel_diff rtol = 1e-8
    ae, be = run_moments(p1; erosion = true, loci = :multilocus), Ext.run_trait_moments(p1; erosion = true, loci = :multilocus)
    @test ae.geno_var ≈ be.geno_var rtol = 1e-10
    # it is a functor: the generation's image is the composite of the steps' images
    steps = (:develop, :select, :reproduce, :record, :rescale, :mutate, :swap)
    imgs = Ext.trait_moment_images()
    m1 = Ext.TraitMoments(p; G = [0.1 0.02; 0.02 0.05]); m2 = deepcopy(m1)
    lower(g, :moments; images = imgs)(m1, nothing)
    foldl(Catlab.compose, (lower(getfield(L, s), :moments; images = imgs) for s in steps))(m2, nothing)
    @test (m1.zbar, m1.G, m1.θ, m1.gen) == (m2.zbar, m2.G, m2.θ, m2.gen)
    @test_throws ArgumentError run_moments(p; cycle = L.develop)
end

@testset "several traits: the Eidos lowering, a second table" begin
    L = lifecycle(); g = L.generation
    base = Params(PopSize = 100, NoQTL = 3, Ploidy = 2, NoGenerations = 50, PhenotypeVar = 0.02,
                  MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)
    p = Params(base; NoTraits = 2, MutCorr = 0.4, SelCorr = -0.3, MoveDirection = [1.0, 2.0])
    s = slim_script(g, p)
    @test occursin("rmvnorm(sum(k0), rep(0.0, NT), MUTCOV)", s) && occursin("defineConstant(\"NT\", 2)", s)
    @test occursin("matrixMult(matrix(SINV, nrow = NT, byrow = T), D)", s)
    @test occursin("\\tz1\\tz2\\tG11\\tG12\\tG22\\ttheta1\\ttheta2", s)
    @test occursin("rep(0.0, L * NT)", s) && occursin("pop.setValue(\"THETA\", rep(0.0, NT));", s)
    @test occursin("MCORR", s) && occursin("MCOV", slim_script(g, Params(p; MutCorr = 0.0, MutCov = [2e-4 1e-4; 1e-4 2e-4])))
    # the same syntax: each step's image comes from the trait table
    imgs = Ext.eidos_trait_images()
    @test Set(keys(imgs)) == Set(keys(Ext.eidos_images()))
    @test occursin("for (r in 0:(NT - 1)) gsum = c(gsum,", s) && occursin("from0 = (loc <= repEach(cut, L))", s)
    # one trait's model is untouched (the golden file test above), and has none of this
    @test !occursin("NT", slim_script(g, base))
    # fecundity selection splices into the trait table as into the one-trait table
    sf = slim_script(g, Params(p; Selection = :fecundity))
    @test occursin("weighted = (selection == \"fecundity\")", sf) && occursin("matrixMult", sf)
    # several demes share the trait parameters; migration moves whole allele vectors
    @test occursin("migRate", slim_script(g ⊗ g ⋅ L.migrate(2), p; migration_rate = 0.01))
    @test_throws ArgumentError slim_script(g ⊗ g ⋅ L.migrate(2), [p, Params(p; SelCorr = 0.2)])
    # what the Eidos lowering cannot do with several traits
    @test_throws ArgumentError slim_script(g, Params(p; MutCorr = 1.0))                 # rmvnorm needs M ≻ 0
    θstep = withparams(g; θ = 0.1)
    @test_throws ArgumentError slim_script(θstep, p)
    @test_throws ArgumentError run_composed(World(p), MersenneTwister(1); cycle = θstep)
    @test_throws ArgumentError run_moments(p; cycle = θstep)
end

@testset "models: step parameters typed by the model, a table of images per model" begin
    L = lifecycle(); g = L.generation
    base = Params(PopSize = 200, NoQTL = 10, Ploidy = 2, NoGenerations = 100, PhenotypeVar = 0.02,
                  MutationScalar = 0.004, MutationRatePerGamete = 0.02, InitFitnessSD = 0.3, RecombPoissonMean = 2.0)
    pt = Params(base; NoTraits = 2, MutCorr = 0.4, SelCorr = -0.3)
    @test model_of(base) == OneTraitModel() && model_of(pt) == TraitModel(2)
    @test_throws ArgumentError TraitModel(1)
    # each model has its own schema: θ is a number for one trait, a vector for two; Ω only for traits
    s1, s2 = step_schema(OneTraitModel()), step_schema(TraitModel(2))
    @test :θ in s1[:select] && :Ω ∉ s1[:select] && :Ω in s2[:select] && :MutCov in s2[:mutate] && :MutCov ∉ s1[:mutate]
    @test :StandardisePhenotype in s1[:select] && :StandardisePhenotype ∉ s2[:select]
    @test Ext.STEP_PARAMS == s1
    θv = withparams(L.select, pt; θ = [0.1, -0.1])
    @test first(Catlab.args(θv)).params.θ == [0.1, -0.1]
    @test withparams(L.select, pt; θ = [1, 0]) == withparams(L.select, TraitModel(2); θ = [1.0, 0.0])   # normalised
    @test_throws ArgumentError withparams(L.select; θ = [0.1, -0.1])                  # one trait: a number
    @test_throws ArgumentError withparams(L.select, pt; θ = 0.1)                      # two traits: a vector
    @test_throws ArgumentError withparams(L.select, pt; θ = [0.1, 0.2, 0.3])
    @test_throws ArgumentError withparams(L.select; Ω = [1.0 0; 0 1])
    @test_throws ArgumentError withparams(L.select, pt; Ω = [1.0 2; 2 1])               # not positive definite
    @test_throws ArgumentError withparams(L.select, pt; Ω = [1.0 0; 0 1], ω = 0.2)      # Ω is the whole covariance
    @test_throws ArgumentError withparams(L.select, pt; StandardisePhenotype = true)
    @test_throws ArgumentError withparams(L.mutate; MutCorr = 0.2)
    # a step written for one model is checked again when lowered for the other
    @test_throws ArgumentError lower(withparams(g; θ = 0.1), :kernel; model = TraitModel(2))
    @test_throws ArgumentError lower(withparams(g, pt; θ = [0.1, 0.0]), :moments)
    @test_throws ArgumentError slim_script(withparams(g, pt; θ = [0.1, 0.0]), base)
    # a table of images per model: the one-trait functor refuses a population of two traits
    @test_throws ArgumentError lower(g, :kernel)(World(pt), MersenneTwister(1))
    @test_throws ArgumentError lower(g, :kernel; model = TraitModel(2))(World(base), MersenneTwister(1))
    @test_throws ArgumentError lower(g, :kernel; model = TraitModel(3))(World(pt), MersenneTwister(1))
    ops = operators(pt)
    w = World(pt); init_env!(w, MersenneTwister(1)); ops.generation(w, MersenneTwister(2))
    @test w.gen == 1
    k1, k2 = Ext.kernel_images(OneTraitModel()), Ext.kernel_images(TraitModel(2))
    @test Set(keys(k1)) == Set(keys(k2))
    @test_throws ArgumentError run_composed(Ext.Stratified([World(base), World(pt)]), MersenneTwister(1); cycle = g ⊗ g)
    # n traits: a step's fields of Params give the run of Params with those values, bit for bit
    same(a, b) = all(getfield(a, k) == getfield(b, k) for k in keys(a) if k != :seed)
    over = (MutationScalar = 0.02, MutationRatePerGamete = 0.05, MoveType = :continuous, MoveFitnessSDGen = 30,
            MoveFitnessSD = 0.004, MoveDirection = [1.0, 2.0], SelCorr = 0.2)
    @test same(run!(World(Params(pt; over...)), MersenneTwister(3)),
               run_composed(World(pt), MersenneTwister(3); cycle = withparams(g, pt; over...)))
    p0 = Params(pt; MutCorr = 0.0)
    C = [0.002 -0.0005; -0.0005 0.001]
    @test same(run!(World(Params(p0; MutCov = C)), MersenneTwister(4)),
               run_composed(World(p0), MersenneTwister(4); cycle = withparams(g, p0; MutCov = C)))
    @test same(run_moments(Params(pt; over...)), run_moments(pt; cycle = withparams(g, pt; over...)))
    # a step's Ω equal to the deme's ω² S is the deme's selection (to rounding)
    Ωd = selection_covariance(pt)
    a, b = run_moments(pt), run_moments(pt; cycle = withparams(g, pt; Ω = Ωd))
    @test a.gmatrix[end] ≈ b.gmatrix[end] rtol = 1e-10
    # a step's θ, a vector: the optimum selected toward, the deme's recorded optimum unchanged
    pθ = Params(pt; NoQTL = 1, PopSize = 500, NoGenerations = 600, MutationRatePerGamete = 0.5, MutationScalar = 0.02)
    o = run_composed(World(pθ), MersenneTwister(5); cycle = withparams(g, pθ; θ = [0.1, -0.05]))
    zbar = sum(o.trait_mean[301:end]) / 300
    @test abs(zbar[1] - 0.1) < 0.02 && abs(zbar[2] + 0.05) < 0.02 && all(iszero, o.theta)
    mo = run_moments(Params(pθ; NoGenerations = 3000); cycle = withparams(g, pθ; θ = [0.1, -0.05]))
    @test mo.trait_mean[end] ≈ [0.1, -0.05] rtol = 1e-4
    # demes with optima along different directions, in the kernel and in the moments
    cyc = withparams(g, pθ; θ = [0.1, 0.0]) ⊗ withparams(g, pθ; θ = [0.0, 0.1]) ⋅ L.migrate(2; rate = 0.01)
    r = run_composed(pθ ⊗ demes(2), MersenneTwister(6); cycle = cyc)
    @test !r.extinct && length(r.gmatrix) == 2 * pθ.NoGenerations
    # Eidos: the step's θ, Ω, direction and mutation are literals in its code
    e = lower(withparams(L.select, pt; θ = [0.1, -0.1]), :eidos; model = TraitModel(2)).code
    @test occursin("THETA = c(0.1, -0.1);", e)
    e = lower(withparams(L.select, pt; Ω = [0.1 0.0; 0.0 0.05], MoveDirection = [3.0, 4.0]), :eidos; model = TraitModel(2)).code
    @test occursin("c(10.0, 0.0, 0.0, 20.0)", e) && occursin("c(0.6, 0.8)", e) && occursin("OMEGA = 1.0;", e)
    e = lower(withparams(L.mutate, pt; MutCov = C), :eidos; model = TraitModel(2)).code
    @test occursin("matrix(c(0.002, -0.0005, -0.0005, 0.001), nrow = 2, byrow = T) / (2 * L * rate)", e)
    # a step's m alone rescales the deme's MUTCOV, unless the deme gives M itself
    e = lower(withparams(L.mutate, pt; MutationScalar = 0.04), :eidos; model = TraitModel(2)).code
    @test occursin("(MUTCOV * (1.0 * (exists(\"MCOV\") ? 1.0 else 0.04 / m)))", e)
end

@testset "several traits: fluctuating optima (kernel, moments, Eidos, theory)" begin
    L = lifecycle(); g = L.generation
    base = Params(PopSize = 200, NoQTL = 10, Ploidy = 2, NoGenerations = 300, PhenotypeVar = 1.0, Heritability = 0.5,
                  MutationScalar = 0.01, MutationRatePerGamete = 0.02, InitFitnessSD = 1.5, NoTraits = 2, MutCorr = 0.4,
                  SelCorr = -0.3, MoveFitnessSDGen = 100, MoveFitnessSD = 0.001, MoveDirection = [1.0, 2.0])
    pr = Params(base; MoveType = :red_noise, MoveNoiseCov = [0.2 0.05; 0.05 0.1], MoveAutocorr = 0.6)
    pp = Params(base; MoveType = :periodic, MoveAmplitude = 0.5, MovePeriod = 20.0, FluctDirection = [1.0, -1.0])
    # parameters: defaults and checks
    @test fluct_direction(base) ≈ [1, 2] ./ sqrt(5) && fluct_direction(pp) ≈ [1, -1] ./ sqrt(2)
    @test move_noise_covariance(Params(base; MoveNoiseSD = 0.3)) ≈ 0.09 .* (fluct_direction(base) * fluct_direction(base)')
    @test move_noise_covariance(pr) == [0.2 0.05; 0.05 0.1]
    bad(; kw...) = @test_throws ArgumentError World(Params(base; kw...))
    bad(MoveNoiseCov = [1.0 2.0; 2.0 1.0]); bad(MoveNoiseCov = [0.1 0; 0 0.1], MoveNoiseSD = 0.1)
    bad(FluctDirection = [0.0, 0.0]); bad(FluctDirection = [1.0])
    @test_throws ArgumentError World(Params(PopSize = 10, FluctDirection = [1.0]))      # one trait: leave unset
    # the kernel: the periodic optimum is exactly the trend plus the cycle along f
    o = simulate(pp; seed = 1)
    θ = [o.trait_mean[i] + o.trait_lag[i] for i in eachindex(o.gen)]
    t = 150 - 100
    @test θ[151] ≈ 0.001 * t .* QTLSim.move_direction(pp) .+ 0.5 * sin(2π * t / 20) .* fluct_direction(pp)
    # red noise: reproducible, nonzero fluctuation, no allocation per generation
    @test simulate(pr; seed = 2).trait_lag == simulate(pr; seed = 2).trait_lag
    w = World(pr); rng = MersenneTwister(3); init_env!(w, rng)
    for s in 0:149; w.gen = s; generation!(w, rng); end
    alloc(w, rng) = @allocated generation!(w, rng)
    alloc(w, rng); @test alloc(w, rng) == 0
    # the moments: the periodic optimum exactly, red noise as its trend, a realised sequence
    mp = run_moments(pp)
    @test mp.trait_mean[151] + mp.trait_lag[151] ≈ θ[151]
    mr = run_moments(pr)
    @test mr.trait_mean[201] + mr.trait_lag[201] ≈ 0.001 * 100 .* QTLSim.move_direction(pr) atol = 1e-12
    seq = [o.trait_mean[i] + o.trait_lag[i] for i in eachindex(o.gen)]
    ms = run_moments(Params(pp; MoveType = :phenotype, MoveFitnessSDGen = -1); optimum = seq)
    @test ms.trait_mean[200] + ms.trait_lag[200] ≈ seq[200]
    @test_throws ArgumentError run_moments(pp; optimum = seq[1:10])
    # a step's own fluctuating move, in the kernel: the run of Params with that move
    q = Params(base; MoveType = :continuous)
    over = (MoveType = :periodic, MoveAmplitude = 0.5, MovePeriod = 20.0, FluctDirection = [1.0, -1.0])
    a = run!(World(Params(q; over...)), MersenneTwister(4))
    b = run_composed(World(q), MersenneTwister(4); cycle = withparams(g, q; over...))
    @test a.trait_lag == b.trait_lag && a.gmatrix == b.gmatrix
    # Eidos: the fluctuating code, its constants, and a step's own noise
    s = slim_script(g, pr)
    @test occursin("xi = drop(matrixMult(NOISEROOT", s) && occursin("defineConstant(\"NOISEL\"", s) &&
          occursin("defineConstant(\"NOISEROOT\", matrix(NOISEL", s)
    s = slim_script(g, pp)
    @test occursin("* FDIR;", s) && occursin("defineConstant(\"NOISEROOT\", (moveNoiseSD *", s)
    @test !occursin("NOISEROOT", slim_script(g, base))          # no fluctuating move, no fluctuating code
    e = lower(withparams(L.select, q; MoveType = :red_noise, MoveNoiseCov = [0.1 0.0; 0.0 0.05]), :eidos;
              model = TraitModel(2)).code
    @test occursin("drop(matrixMult(matrix(c(", e)
    # theory: one trait is red_noise_lag_var and periodic_response
    @test red_noise_lag_cov(fill(0.49, 1, 1), 0.6; G = fill(0.3, 1, 1), A = fill(0.7, 1, 1))[1] ≈
          red_noise_lag_var(0.7, 0.6; VG = 0.3, Vs = 0.7) rtol = 1e-12
    @test periodic_lag_cov(15.0, [1.0]; G = fill(0.3, 1, 1), A = fill(0.7, 1, 1))[1] ≈
          periodic_response(15.0; VG = 0.3, Vs = 0.7).lag_amplitude^2 / 2 rtol = 1e-12
    # two traits: against the linear recursion iterated (periodic, exactly) and simulated (red noise)
    G = [0.3 0.1; 0.1 0.2]; A = [1.0 -0.2; -0.2 0.8]; f = [0.6, 0.8]
    B = G / (A + G)
    z = zeros(2); acc = zeros(2, 2)
    for t in 0:20_000
        λ = sin(2π * t / 15) .* f - z
        t >= 10_000 && (acc .+= λ * λ')
        z = z + B * λ
    end
    @test acc ./ 10_001 ≈ periodic_lag_cov(15.0, f; G, A) rtol = 1e-3
    Σ = [0.5 0.2; 0.2 0.3]; Lc = QTLSim.cholesky(Σ).L; r = MersenneTwister(5)
    z = zeros(2); ε = Lc * randn(r, 2); acc = zeros(2, 2); T = 400_000
    for t in 1:T
        λ = ε - z; acc .+= λ * λ'
        z = z + B * λ
        ε = 0.6 .* ε + sqrt(1 - 0.36) .* (Lc * randn(r, 2))
    end
    @test acc ./ T ≈ red_noise_lag_cov(Σ, 0.6; G, A) rtol = 0.03
    # expected_lag_cov: the trend's outer product plus the fluctuation; the load
    Gp = gaussian_G_loci(pr)
    λs = steady_lag_vector(0.001 .* QTLSim.move_direction(pr); G = Gp, A = selection_variance_matrix(pr))
    @test expected_lag_cov(pr; G = Gp) ≈ λs * λs' + red_noise_lag_cov(move_noise_covariance(pr), 0.6; G = Gp,
                                                                     A = selection_variance_matrix(pr))
    @test trait_lag_load(pr; G = Gp) > 0
    @test_throws ArgumentError expected_lag_cov(Params(base; MoveType = :phenotype); G = Gp)
end

@testset "several traits: the moments of several demes" begin
    L = lifecycle(); g = L.generation
    p = Params(Ploidy = 2, NoQTL = 10, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
               MutationRatePerGamete = 0.02, PopSize = 100, NoGenerations = 400, Sample = 10, InitFitnessSD = 5.0)
    # one trait through the trait table of k demes is the one-trait k-deme moments
    θs = [-1.0, 0.0, 0.5, 2.0]
    cyc = foldl(⊗, [withparams(g; θ = t) for t in θs]) ⋅ L.migrate(4)
    imgs = Dict{Any,Any}(first(Catlab.args(withparams(L.select; θ = t))) =>
                Ext.Operator(:sel, (m, rng) -> (θ0 = m.θ; Ext.trait_moments_select!(m; θ = [t]); m.θ = θ0; m)) for t in θs)
    for (between, loci) in ((true, :multilocus), (false, :multilocus), (false, :aggregate))
        a = run_moments(p; cycle = cyc, migration_rate = 0.01, between, loci, erosion = !between)
        b = Ext.run_trait_moments(fill(p, 4); cycle = cyc, images = imgs, migration_rate = 0.01, between, loci,
                                  erosion = !between)
        @test a.geno_var ≈ b.geno_var rtol = 1e-10
        @test a.pheno_mean ≈ b.pheno_mean rtol = 1e-10 atol = 1e-14
        between && @test a.metapop.var_among ≈ [x[1, 1] for x in b.metapop.var_among] rtol = 1e-10
    end
    # two traits, local optima along different directions
    p2 = Params(p; NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)
    θ2 = [[2.0, 0.0], [0.0, 2.0], [-2.0, 0.0]]
    c2 = foldl(⊗, [withparams(g, p2; θ = t) for t in θ2]) ⋅ L.migrate(3)
    r = run_moments(p2; cycle = c2, migration_rate = 0.01, between = true)
    @test !r.extinct && length(r.trait_mean) == 3 * p2.NoGenerations ÷ p2.Sample && r.deme[1:3] == 1:3
    @test all(X -> X ≈ X' && all(>(0), eigvals(Symmetric(X))), r.metapop.within[2:end])
    @test all(q -> all(0 .<= q .< 1), r.metapop.qst[2:end])
    # deme 1 is pulled toward its optimum along trait 1, deme 2 along trait 2
    @test r.trait_mean[end-2][1] > 0.1 && r.trait_mean[end-1][2] > 0.1
    # without drift the demes are infinite: no Q, and the same means from :aggregate's start
    r0 = run_moments(p2; cycle = c2, migration_rate = 0.01)
    @test length(r0.gmatrix) == length(r.gmatrix) && !haskey(r0, :metapop)
    @test_throws ArgumentError run_moments(p2; cycle = c2, between = true, loci = :aggregate)
end

@testset "rare alleles with locally favoured alleles (loci = :rare_local, src/local_alleles.jl)" begin
    L_ = lifecycle()
    p = Params(Ploidy = 2, NoQTL = 20, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
               MutationRatePerGamete = 2e-4, PopSize = 100, NoGenerations = 600, Sample = 10, InitFitnessSD = 5.0)
    p2 = Params(p; NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)
    cyc1(θs) = foldl(⊗, [withparams(L_.generation, p; θ = t) for t in θs]) ⋅ L_.migrate(length(θs))
    cyc2(φs) = foldl(⊗, [withparams(L_.generation, p2; θ = 2.0 .* [cosd(φ), sind(φ)]) for φ in φs]) ⋅ L_.migrate(length(φs))
    flat(x) = reduce(vcat, [y isa Number ? [y] : vec(y) for y in x])
    rel(a, b) = maximum(abs.(flat(a) .- flat(b))) / maximum(abs.(flat(b)))
    # reductions to the current closure (to 1e-12; they are in fact identical)
    for (f, tag) in ((l -> run_moments(p; erosion = true, loci = l).geno_var, "one deme"),
                     (l -> run_moments(p; between = true, loci = l).metapop.within, "one deme, between"),
                     (l -> run_moments(p2; erosion = true, loci = l).gmatrix, "one deme, two traits"),
                     (l -> run_moments(p; cycle = cyc1(zeros(4)), migration_rate = 0.01, between = true, loci = l).metapop.within,
                      "common optimum"),
                     (l -> run_moments(p; cycle = cyc1([2.0, 1, -1, -2]), migration_rate = 0.0, between = true, loci = l).metapop.within,
                      "no migration"),
                     (l -> run_moments(p2; cycle = cyc2([0.0, 90.0]), migration_rate = 0.0, between = true, loci = l).metapop.within,
                      "no migration, two traits"))
        a, b = f(:rare_local), f(:rare_alleles)
        @test rel(a, b) <= 1e-12
    end
    # local optima with migration: the gain raises the within-deme variance above the current
    # closure (Lean local_fixed_point_ge) and not above the multilocus Gaussian model
    w(l) = run_moments(p; cycle = cyc1([2.0, 1, -1, -2]), migration_rate = 0.01, between = true, loci = l).metapop.within[end]
    @test w(:rare_alleles) < w(:rare_local) <= w(:multilocus) * (1 + 1e-9)
    w2(l) = QTLSim.tr(run_moments(p2; cycle = cyc2([0.0, 60, 120, 180]), migration_rate = 0.01, between = true, loci = l).metapop.within[end])
    @test w2(:rare_alleles) < w2(:rare_local) <= w2(:multilocus) * (1 + 1e-9)
    # the gain: the eigenvalue expansion, the island form, the reductions and the limits
    rng = Xoshiro(11)
    for k in (2, 3, 5)
        M = rand(rng, k, k) .* 0.05; for i in 1:k; M[i, i] = 0; end
        P = QTLSim.Diagonal(ones(k)) - QTLSim.Diagonal(vec(sum(M; dims = 2))) + M
        σ = 0.01 .* randn(rng, k)
        lead(e) = maximum(real, eigvals(P * (QTLSim.Diagonal(ones(k)) + e * QTLSim.Diagonal(σ))))
        h = 1e-3
        @test (lead(h) - 2lead(0.0) + lead(-h)) / (2h^2) ≈ σ' * QTLSim.local_weights(P) * σ rtol = 1e-4
    end
    P4(m) = QTLSim.migration_mixing(island_migration(4, m), fill(100.0, 4), trues(4))
    ℓ = [0.01, 0.004, -0.004, -0.01]
    varℓ = sum(abs2, ℓ .- sum(ℓ) / 4) / 4
    for m in (0.001, 0.01, 0.1)
        P = P4(m); mp = 1 - (P[1, 1] - P[1, 2])
        @test all(QTLSim.local_gain(P, ℓ) .≈ (1 - mp) / mp * varℓ)
    end
    @test QTLSim.local_gain(P4(0.01), ℓ)[1] > QTLSim.local_gain(P4(0.1), ℓ)[1] > 0
    @test abs(QTLSim.local_gain(fill(0.25, 4, 4), ℓ)[1]) < 1e-15        # complete mixing, m′ = 1
    @test all(iszero, QTLSim.local_gain(P4(0.01), fill(0.3, 4)))
    @test all(iszero, QTLSim.local_gain(P4(0.01), fill([0.3, -0.1], 4)))
    @test all(iszero, QTLSim.local_gain(Matrix(1.0 * QTLSim.I, 4, 4), ℓ))
    @test QTLSim.local_gain(fill(1.0, 1, 1), [0.5]) == [0.0]
    @test QTLSim.migration_classes([0.9 0.1 0; 0.1 0.9 0; 0 0 1.0]) == [[1, 2], [3]]
    @test QTLSim.local_sinv(0.04, 0.0) === 0.04
    @test QTLSim.local_sinv(0.04, 0.01) ≈ 0.02
    @test QTLSim.local_sinv(0.04, 0.05) == 0.0
    S = [0.05 0.01; 0.01 0.03]
    @test QTLSim.local_sinv(S, zeros(2, 2)) === S
    @test minimum(eigvals(Symmetric(QTLSim.local_sinv(S, [0.04 0.0; 0.0 -0.01])))) >= -1e-15
    ω, _ = QTLSim.pair_classes(p)
    st = QTLSim.hap_init(0.3, 100.0, 2, length(ω))
    α2 = QTLSim.mutation_effect_variance(p)
    @test QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2) == QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2, gain = 0.0)
    @test QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2, gain = 1.0)[1] ==
          QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2 = 0.0)[1]           # the rare loss switched off
    @test_throws ArgumentError run_moments(p; loci = :rare_locals)
end
