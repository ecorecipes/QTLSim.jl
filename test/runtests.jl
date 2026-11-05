using QTLSim
using Test, Statistics, Random

# The categorical layer is the package extension QTLSimCatlabExt, loaded with Catlab.
import Catlab
using Catlab: ⋅, ⊗
const Ext = Base.get_extension(QTLSim, :QTLSimCatlabExt)
const (develop, select, reproduce, mutate, swap, generation) = values(operators())

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
    mutate!(w, rng)
    # genotype change per individual = Σ mutations added to its two chromosomes
    Δ = [sum(@view w.sgenome[:, :, j]) for j in 1:p.PopSize]
    @test isapprox(var(Δ), muteVar(p); rtol=0.15)
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
    for i in 1:p.PopSize; w.env[i] = QTLSim.envsd(p) * randn(rng); end
    generation!(w, rng)
    @test w.gen == 1
    @test w.extinct == false
end

# ─── the categorical layer (package extension) ───────────────────────────

@testset "Catlab loads the categorical layer" begin
    @test Ext isa Module
    @test generation isa Ext.Operator && length(operators()) == 6
end


@testset "Catlab ACSet schema + lift/lower round-trip" begin
    p = Params(PopSize=8, NoQTL=3, Ploidy=2, PhenotypeVar=0.02,
               Heritability=0.5, MutationScalar=0.001)
    w = World(p); rng = MersenneTwister(11)
    for k in eachindex(w.genome); w.genome[k] = randn(rng); end
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
    @test genome_from_acset(pop, p.NoQTL, p.Ploidy) == w.genome
    w2 = World(p); lower!(pop, w2)
    @test w2.genome == w.genome
    @test w2.env == w.env
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
    for w in s.demes, i in 1:w.p.PopSize
        for k in eachindex(w.genome); w.genome[k] = randn(rng); end
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

@testset "Params(p; kw...) copies with overrides" begin
    p = Params(NoQTL = 10, Seed = 3, MoveType = :continuous)
    q = Params(p; NoGenerations = 5, MoveFitnessSD = 0.2)
    @test q.NoGenerations == 5 && q.MoveFitnessSD == 0.2
    @test q.NoQTL == 10 && q.Seed == 3 && q.MoveType === :continuous
    @test Params(p) == p
    @test_throws ArgumentError Params(p; Nonsense = 1)
end
