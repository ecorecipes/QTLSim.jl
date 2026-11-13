# Tests: catlab. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl catlab`.

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
