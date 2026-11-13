# Tests: traits. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl traits`.

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
    # regression: a fixed-slope reaction norm with several traits (refused by checkparams) failed
    # inside the Eidos splice of select with an internal error; it is refused before any code
    pp = Params(p; Plasticity = 0.2, MoveType = :continuous)
    for q in (pp, Params(pp; Selection = :fecundity), Params(pp; PlasticityCostSD = 1.0), Params(pp; CueLag = 2))
        @test_throws ArgumentError QTLSim.checkparams(q)
        @test_throws ArgumentError slim_script(g, q)
    end
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
