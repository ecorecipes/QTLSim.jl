using QTLSim
using Test, Statistics, Random

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
        fill!(w.sgenome, 0.0)
        mutate!(w, rng)
        var(sum(@view w.sgenome[:, :, j]) for j in 1:p.PopSize)
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
    p = Params(PopSize=200, NoQTL=3, Ploidy=1, NoGenerations=50, PhenotypeVar=0.02, MutationScalar=0.004)
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
    @test w2.genotype == w.genotype && w2.phenotype == w.phenotype   # recomputed by lower!
    @test_throws ArgumentError lower!(pop, World(Params(p; PopSize = 9)))
    @test_throws ArgumentError genome_from_acset(pop, p.NoQTL + 1, p.Ploidy)
    # haploid
    ph = Params(p; Ploidy = 1); wh = World(ph); wh.genome .= randn(rng, size(wh.genome)); develop!(wh)
    wh2 = World(ph); lower!(lift(wh), wh2)
    @test wh2.genome == wh.genome
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
    @test w1.genome == w2.genome && w1.env == w2.env && w1.θ == w2.θ && w1.gen == w2.gen == 1
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
    @test all(s1.demes[i].genome == s2.demes[i].genome for i in 1:3)
    # a braiding relabels the demes and moves no data
    b = lower((g ⊗ g) ⋅ Catlab.braid(Z, Z), :kernel)
    @test b.perm == [2, 1]
    s3 = p ⊗ demes(2); init!(s3, MersenneTwister(2)); s4 = deepcopy(s3)
    b.act(s3, [1, 2], MersenneTwister(4)); lower(g ⊗ g, :kernel).act(s4, [1, 2], MersenneTwister(4))
    @test s3.demes[1].genome == s4.demes[1].genome && s3.demes[2].genome == s4.demes[2].genome
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
    @test ndemes(s) == 2 && s.demes[1] !== s.demes[2] && s.demes[1].genome == w.genome
    # migration needs demes of one shape, and checks before moving anyone
    s = Ext.Stratified([World(p), World(Params(p; PopSize = 50))])
    g1 = copy(s.demes[1].genome)
    @test_throws ArgumentError migrate!(s, MersenneTwister(1); rate = 0.5)
    @test s.demes[1].genome == g1
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
