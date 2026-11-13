# Tests: core. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl core`.

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

@testset "Params(p; kw...) copies with overrides" begin
    p = Params(NoQTL = 10, Seed = 3, MoveType = :continuous)
    q = Params(p; NoGenerations = 5, MoveFitnessSD = 0.2)
    @test q.NoGenerations == 5 && q.MoveFitnessSD == 0.2
    @test q.NoQTL == 10 && q.Seed == 3 && q.MoveType === :continuous
    @test Params(p) == p
    @test_throws ArgumentError Params(p; Nonsense = 1)
end
