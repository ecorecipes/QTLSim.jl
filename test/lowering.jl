# Tests: lowering. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl lowering`.

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
