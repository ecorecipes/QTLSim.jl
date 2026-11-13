# Tests: io. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl io`.

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
