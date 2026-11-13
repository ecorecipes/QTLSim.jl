# Tests: plasticity. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl plasticity`.

@testset "phenotypic plasticity: a reaction norm of fixed slope (kernel, moments, Eidos, theory)" begin
    L = lifecycle(); g = L.generation
    same(a, b; but = ()) = all(getfield(a, k) == getfield(b, k) for k in keys(a) if !(k in but))
    base = Params(NoGenerations = 1500, PopSize = 500, Ploidy = 2, NoQTL = 20, MutationRatePerGamete = 0.01,
                  MoveType = :continuous, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.02, RecombPoissonMean = 3.0,
                  InitFitnessSD = 1.0, PhenotypeVar = 1.0, Heritability = 0.5)
    short = Params(base; NoGenerations = 150, MoveFitnessSDGen = 50, PopSize = 200)

    # --- off by default, and nothing changes without it
    @test Params().Plasticity == 0 && Params().OptimumSlope == 1 && Params().CueLag == 0 &&
          Params().PlasticityCostSD == Inf && !Params().ReactionNorm
    @test !QTLSim.isplastic(Params()) && plasticity_ratio(Params()) == 0 && plasticity_cost(Params()) == 1
    # a slope with the optimum at the origin throughout adds exactly zero: the plastic path is
    # the default path, bit for bit (the default path itself is checked against the dense
    # reference kernel and the golden SLiM files above)
    still = Params(short; MoveFitnessSDGen = -1)
    @test same(simulate(still; seed = 5), simulate(Params(still; Plasticity = 0.7, CueLag = 2); seed = 5))
    # the cost alone thins the survivors
    @test mean(simulate(Params(still; Plasticity = 1.0, PlasticityCostSD = 1.0); seed = 5).mean_fitness) <
          0.7 * mean(simulate(still; seed = 5).mean_fitness)
    # no allocation per generation with plasticity, as without
    alloc(w, rng) = @allocated generation!(w, rng)
    for q in (short, Params(short; Plasticity = 0.5, CueLag = 3, PlasticityCostSD = 2.0))
        w = World(q); rng = MersenneTwister(1); init_env!(w, rng)
        for _ in 1:80; generation!(w, rng); end
        alloc(w, rng)
        @test alloc(w, rng) == 0
    end

    # --- the kernel expresses z = g + e + (b/B) θ_{t−τ}, exactly
    q = Params(short; Plasticity = 1.2, OptimumSlope = 2.0, CueLag = 3)
    w = World(q); rng = MersenneTwister(2); init_env!(w, rng)
    θs = Float64[]
    ok = true
    for t in 0:99
        w.gen = t
        develop!(w)
        @test w.plastic == 0
        z0 = copy(w.phenotype)
        select!(w, rng)
        push!(θs, w.θ)
        cue = θs[max(t + 1 - 3, 1)]              # before the run: the first environment seen
        ok &= w.phenotype ≈ z0 .+ 0.6 * cue && w.plastic ≈ 0.6 * cue
        z1 = copy(w.phenotype)
        QTLSim.express_plasticity!(w)            # a second expression in a generation changes nothing
        ok &= w.phenotype == z1
        w.extinct && break
        reproduce!(w, rng); mutate!(w, rng); swapbuffers!(w)
    end
    @test ok && θs[end] ≈ 0.02 * 50

    # --- with a perfect cue and no cost, z − θ = g + e − (1 − b/B)θ: the run at rate k with slope
    # b is the run without plasticity at rate (1 − b/B)k (the same survivors and genotypes,
    # the fitnesses equal to rounding; the optimum and the phenotypes differ)
    ra = simulate(Params(short; MoveFitnessSD = 0.04, Plasticity = 1.5, OptimumSlope = 2.0); seed = 6)
    rb = simulate(Params(short; MoveFitnessSD = 0.01); seed = 6)
    @test ra.nsurv == rb.nsurv && ra.geno_var == rb.geno_var && ra.mean_fitness ≈ rb.mean_fitness
    @test ra.theta .- ra.pheno_mean ≈ rb.theta .- rb.pheno_mean && ra.theta[end] ≈ 4 * rb.theta[end]

    # --- checkparams
    @test_throws ArgumentError World(Params(short; Plasticity = 0.5, MoveType = :phenotype))
    @test_throws ArgumentError World(Params(short; Plasticity = 0.5, MoveType = :fitness))
    @test World(Params(short; Plasticity = 0.5, MoveType = :phenotype, MoveFitnessSDGen = -1)) isa World
    @test_throws ArgumentError World(Params(short; Plasticity = Inf))
    @test_throws ArgumentError World(Params(short; OptimumSlope = 0.0))
    @test_throws ArgumentError World(Params(short; CueLag = -1))
    @test_throws ArgumentError World(Params(short; PlasticityCostSD = 0.0))
    @test_throws ArgumentError World(Params(short; NoTraits = 2, Plasticity = 0.5))
    @test_throws ArgumentError write_qtsim_options(IOBuffer(), Params(short; Plasticity = 0.5))

    # --- theory: exact reductions at b = 0
    VG, Vs = 0.11, 1.5
    @test steady_lag(0.02; VG, Vs) === 0.02 * (VG + Vs) / VG
    @test steady_lag(0.02; VG, Vs, b = 0.0, B = 3.0) === steady_lag(0.02; VG, Vs)
    @test steady_lag(base; VG) === steady_lag(0.02; VG, Vs = selection_variance(base))
    @test red_noise_lag_var(0.3, 0.6; VG, Vs, b = 0, τ = 4) === red_noise_lag_var(0.3, 0.6; VG, Vs)
    kw = (n = 400, VG = 0.11, ω = 1.0, VE = 0.5, sd50 = 0.9, N = 1000)
    @test critical_rate(; kw..., b = 0.0, B = 2.0, ωb = 0.5) === critical_rate(; kw...)
    fl = Params(base; MoveType = :red_noise, MoveNoiseSD = 0.3, MoveAutocorr = 0.6)
    pe = Params(base; MoveType = :periodic, MoveAmplitude = 0.4, MovePeriod = 30.0)
    for q0 in (base, fl, pe)
        @test expected_lag2(Params(q0; OptimumSlope = 2.0, CueLag = 3); VG) === expected_lag2(q0; VG)
        @test lag_load(Params(q0; OptimumSlope = 2.0, CueLag = 3, PlasticityCostSD = 1.0); VG) === lag_load(q0; VG)
    end
    # with plasticity: the lag and the critical rate scale by 1 − b/B and its inverse
    @test steady_lag(0.02; VG, Vs, b = 0.5, B = 2.0) ≈ 0.75 * steady_lag(0.02; VG, Vs)
    @test critical_rate(; kw..., b = 0.5, B = 2.0) ≈ critical_rate(; kw...) / 0.75
    @test critical_rate(; kw..., b = 0.5, B = 2.0, ωb = 1.0) < critical_rate(; kw..., b = 0.5, B = 2.0)
    @test critical_rate(; kw..., b = 2.0, B = 2.0) == Inf
    @test isnan(critical_rate(; kw..., b = 2.0, B = 1.0, ωb = 0.1))          # the cost leaves no growth
    # the autocovariance of the lag behind red noise, against the lag's moving-average form
    # L_t = θ_t − β Σ_{j ≥ 1} (1 − β)^{j−1} θ_{t−j}, and the variance of L_t − φ L_{t−τ}
    σ, ρ = 0.3, 0.6
    β = VG / (VG + Vs)
    c = [1.0; [-β * (1 - β)^(j - 1) for j in 1:2000]]
    Γ(h) = σ^2 * sum(c[i] * sum(c[j] * ρ^abs(i - j + h) for j in eachindex(c)) for i in eachindex(c))
    for h in (0, 1, 3)
        @test red_noise_lag_autocov(σ, ρ, h; VG, Vs) ≈ Γ(h) rtol = 1e-8
    end
    @test red_noise_lag_autocov(σ, 1 - β, 3; VG, Vs) ≈ red_noise_lag_autocov(σ, 1 - β + 1e-9, 3; VG, Vs) rtol = 1e-6
    @test red_noise_lag_var(σ, ρ; VG, Vs, b = 0.5, B = 2.0) ≈ 0.75^2 * red_noise_lag_var(σ, ρ; VG, Vs)
    @test red_noise_lag_var(σ, ρ; VG, Vs, b = 0.5, B = 2.0, τ = 3) ≈ (1 + 0.25^2) * Γ(0) - 2 * 0.25 * Γ(3) rtol = 1e-8
    # an old cue of white noise makes plasticity add noise
    @test red_noise_lag_var(σ, 0.0; VG, Vs, b = 0.5, τ = 5) > red_noise_lag_var(σ, 0.0; VG, Vs)
    @test expected_lag2(Params(pe; Plasticity = 0.5); VG) ≈
          0.25 * expected_lag2(pe; VG)                                       # a perfect cue: (1 − φ)² throughout
    @test lag_load(Params(base; Plasticity = 0.5, PlasticityCostSD = 1.0); VG) ≈
          lag_load(Params(base; Plasticity = 0.5); VG) + 0.125

    # --- the sustainable rate is the moments recursion's threshold (density regulation)
    bh = Params(base; Demography = :beverton_holt, GrowthRate = 2.0, PopSize = 1000, NoGenerations = 3000, MoveFitnessSDGen = 10)
    V0 = gaussian_vg_exact(bh)               # the aggregate recursion's fixed point: V_G then stays put
    k0 = sustainable_rate(bh; VG = V0)
    @test sustainable_rate(Params(bh; Plasticity = 0.5); VG = V0) ≈ 2k0
    @test sustainable_rate(Params(bh; Plasticity = 0.5, OptimumSlope = 2.0); VG = V0) ≈ k0 / 0.75
    @test k0 / 0.75 > sustainable_rate(Params(bh; Plasticity = 0.5, OptimumSlope = 2.0, PlasticityCostSD = 1.0); VG = V0) > 0
    @test sustainable_rate(Params(bh; Plasticity = 3.0, PlasticityCostSD = 0.5); VG = V0) == 0
    @test sustainable_rate(Params(bh; Demography = :fixed); VG = V0) > k0       # a fixed size dies only at 2/N
    for q in (bh, Params(bh; Plasticity = 0.5), Params(bh; Plasticity = 0.6, OptimumSlope = 2.0, PlasticityCostSD = 1.0, CueLag = 2))
        kc = sustainable_rate(q; VG = V0)
        @test !run_moments(Params(q; MoveFitnessSD = 0.98kc); VG = V0, loci = :aggregate).extinct
        @test run_moments(Params(q; MoveFitnessSD = 1.02kc); VG = V0, loci = :aggregate).extinct
    end

    # --- kernel against the moments lowering and the theory: the steady lag under a moving
    # optimum, mean of 8 runs, generations 1401-1500. Tolerances: the kernel's lag within 8% of
    # the steady lag at the kernel's own V_G, and within 40% of the moments lowering's (whose
    # V_G does not fall under the move as the kernel's does: 36% without plasticity, 4-11% with
    # it, where the lag is smaller); mean fitness within 3%
    w = 1401:1500
    lags = Dict{Any,Float64}()
    for (b, B, τ, cost) in ((0.0, 1.0, 0, Inf), (0.5, 1.0, 0, Inf), (1.0, 2.0, 4, Inf), (0.5, 1.0, 0, 1.0))
        q = Params(base; Plasticity = b, OptimumSlope = B, CueLag = τ, PlasticityCostSD = cost)
        rs = [simulate(q; seed) for seed in 1:8]
        m = run_moments(q; erosion = true)
        lagk = mean(mean(r.theta[w] .- r.pheno_mean[w]) for r in rs)
        lagm = mean(m.theta[w] .- m.pheno_mean[w])
        VGk = mean(mean(r.geno_var[w]) for r in rs)
        lags[(b, B, τ, cost)] = lagk
        @test lagk ≈ steady_lag(q; VG = VGk) rtol = 0.08
        @test lagk ≈ lagm rtol = 0.4
        @test lagm ≈ steady_lag(q; VG = m.geno_var[end]) rtol = 0.01
        @test mean(mean(r.mean_fitness[w]) for r in rs) ≈ mean(m.mean_fitness[w]) rtol = 0.03
        # the moments' mean phenotype is the mean genotype plus the plastic part
        @test m.pheno_mean[end] - m.theta[end] ≈ -lagm rtol = 0.05
    end
    # plasticity at least halves the lag here (b/B = 1/2; more, as V_G falls less), and a
    # delayed cue of a steadily moving optimum makes no difference
    @test lags[(0.5, 1.0, 0, Inf)] < 0.5 * lags[(0.0, 1.0, 0, Inf)]
    @test lags[(1.0, 2.0, 4, Inf)] ≈ lags[(0.5, 1.0, 0, Inf)] rtol = 0.1

    # --- a step's own reaction norm: the run of Params with it, in the kernel and the moments
    over = (Plasticity = 0.8, OptimumSlope = 2.0, CueLag = 2, PlasticityCostSD = 1.5)
    @test all(k -> k in step_schema(OneTraitModel())[:select], keys(over))
    @test !(:Plasticity in step_schema(OneTraitModel())[:develop]) && !(:Plasticity in step_schema(TraitModel(2))[:select])
    @test same(run!(World(Params(short; over...)), MersenneTwister(3)),
               run_composed(World(short), MersenneTwister(3); cycle = withparams(g; over...)))
    @test same(run_moments(Params(short; over...)), run_moments(short; cycle = withparams(g; over...)))
    @test same(run!(World(Params(short; over...)), MersenneTwister(3)),
               run_composed(World(Params(short; over...)), MersenneTwister(3)))
    @test_throws ArgumentError withparams(L.develop; Plasticity = 0.5)
    @test_throws ArgumentError withparams(g, TraitModel(2); Plasticity = 0.5)
    # two demes, one plastic: each deme keeps its own cue
    pl = Params(short; Plasticity = 0.5, CueLag = 1)
    s = run_moments([pl, short]; cycle = g ⊗ g, loci = :multilocus)
    @test s.pheno_mean[s.deme .== 1] ≈ run_moments(pl; loci = :multilocus).pheno_mean
    @test s.pheno_mean[s.deme .== 2] ≈ run_moments(short; loci = :multilocus).pheno_mean
    sk = run_composed([pl, short] ⊗ demes(2), MersenneTwister(4); cycle = g ⊗ g)
    lagof(d) = mean((sk.theta .- sk.pheno_mean)[(sk.deme .== d) .& (sk.gen .>= 100)])
    @test lagof(1) < 0.75 * lagof(2)

    # --- Eidos: spliced in only when used
    plain = slim_script(g, short)
    @test !occursin("plasticity", plain) && !occursin("CUES", plain)
    sp = slim_script(g, Params(short; over...))
    @test occursin("defineConstant(\"plasticity\", 0.8)", sp) && occursin("defineConstant(\"cueLag\", 2)", sp) &&
          occursin("defineConstant(\"optSlope\", 2.0)", sp) && occursin("defineConstant(\"plastCostSD\", 1.5)", sp)
    @test occursin("pheno = pheno + (plasticity / optSlope) * cues[cueLag];", sp) && occursin("p1.setValue(\"CUES\", cues);", sp)
    @test occursin("defineConstant(\"plastCostSD\", INF)", slim_script(g, Params(short; Plasticity = 0.5)))
    # the rest of the model is the plain one's
    strip(s) = filter(l -> !occursin(r"plasticity|optSlope|cueLag|plastCostSD|cues|CUES|// cueLag", l), split(s, '\n'))
    @test strip(sp) == strip(plain)
    # a step's slope is a literal; per-deme slopes are deme state
    st = slim_script(withparams(g; Plasticity = 0.3), short)
    @test occursin("pheno = pheno + (0.3 / optSlope) * cues[cueLag];", st) && occursin("defineConstant(\"plasticity\", 0.0)", st)
    @test occursin("(0.3 / optSlope)", Ext.lower(withparams(L.select; Plasticity = 0.3), :eidos).code)
    sd = slim_script(g ⊗ g, [Params(short; Plasticity = 0.5), short])
    @test occursin("defineConstant(\"plasticity\", c(0.5, 0.0))", sd) && occursin("p2.getValue(\"plasticity\")", sd)
    @test Ext.eidos_literal(Inf) == "INF" && Ext.eidos_literal(-Inf) == "-INF"
end

@testset "phenotypic plasticity: a heritable reaction norm (ReactionNormModel)" begin
    L = lifecycle(); g = L.generation
    same(a, b; but = ()) = all(getfield(a, k) == getfield(b, k) for k in keys(a) if !(k in but))
    base = Params(NoGenerations = 1500, PopSize = 500, Ploidy = 2, NoQTL = 20, MutationRatePerGamete = 0.01,
                  MoveType = :continuous, MoveFitnessSDGen = 1000, MoveFitnessSD = 0.02, RecombPoissonMean = 3.0,
                  InitFitnessSD = 1.0, PhenotypeVar = 1.0, Heritability = 0.5)
    Vm = muteVar(base)
    rn = Params(base; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2, MutCov = [Vm 0.0; 0.0 Vm / 4], PlasticityCostSD = 3.0)
    @test model_of(rn) == ReactionNormModel() && QTLSim.ntraits(ReactionNormModel()) == 2 && QTLSim.isplastic(rn)
    @test sprint(show, ReactionNormModel()) == "ReactionNormModel()"
    @test_throws ArgumentError World(Params(base; ReactionNorm = true))                        # needs two heritable values
    @test_throws ArgumentError World(Params(rn; NoTraits = 3, MutCov = zeros(0, 0)))
    @test_throws ArgumentError World(Params(rn; SelCorr = 0.3))
    @test_throws ArgumentError World(Params(rn; MoveDirection = [0.0, 1.0]))
    @test_throws ArgumentError World(Params(rn; MoveType = :phenotype))
    @test :Plasticity in step_schema(ReactionNormModel())[:select] && :MutCov in step_schema(ReactionNormModel())[:mutate]

    # --- the kernel: z = g₁ + (b₀ + g₂) ε + e₁, fitness on z and the cost on the slope
    short = Params(rn; NoGenerations = 120, MoveFitnessSDGen = 40, PopSize = 200, CueLag = 1, OptimumSlope = 2.0)
    w = World(short); rng = MersenneTwister(1); init_env!(w, rng)
    prev = 0.0
    for t in 0:99
        w.gen = t
        develop!(w); select!(w, rng)
        ε = (t == 0 ? w.θ : prev) / 2
        G, E, Z = traitvalues(w; what = :genotype), traitvalues(w; what = :env), traitvalues(w)
        @test Z[2, :] ≈ 0.2 .+ G[2, :] && Z[1, :] ≈ G[1, :] .+ Z[2, :] .* ε .+ E[1, :]
        @test w.fitness[1:w.N] ≈ QTLSim.TRAIT_FITNESS_PEAK .* exp.(-(Z[1, :] .- w.θ) .^ 2 ./ 2 .- Z[2, :] .^ 2 ./ 18)
        prev = w.θ
        reproduce!(w, rng); mutate!(w, rng); swapbuffers!(w)
    end
    a = run!(World(short), MersenneTwister(2))
    @test same(a, run_composed(World(short), MersenneTwister(2)))
    @test same(a, run_composed(World(Params(short; CueLag = 0)), MersenneTwister(2); cycle = withparams(g, short; CueLag = 1)))
    @test same(run!(World(Params(short; MutCov = [Vm 0.0; 0.0 Vm])), MersenneTwister(2)),
               run_composed(World(short), MersenneTwister(2); cycle = withparams(g, short; MutCov = [Vm 0.0; 0.0 Vm])))
    @test length(a.trait_mean[1]) == 2 && size(a.gmatrix[1]) == (2, 2)
    @test a.pheno_mean ≈ first.(a.trait_mean) && a.theta .- a.pheno_mean ≈ first.(a.trait_lag)
    alloc(w, rng) = @allocated generation!(w, rng)
    alloc(w, rng)
    @test alloc(w, rng) == 0
    # the wrong model's images refuse the population
    @test_throws ArgumentError lower(g, :kernel; model = TraitModel(2))(World(short), MersenneTwister(1))
    @test_throws ArgumentError lower(g, :kernel; model = ReactionNormModel())(World(Params(base; NoTraits = 2)), MersenneTwister(1))
    @test_throws ArgumentError slim_script(g, short)
    @test_throws ArgumentError run_moments(short; loci = :rare_alleles)

    # --- theory: one generation of selection on the norm
    G = [0.2 0.03; 0.03 0.05]
    r = reaction_norm_response(G, 0.7; ω = 1.0, VE = 0.5, θ = 1.3, abar = 0.4, bbar = 0.6)
    v = [1.0, 0.7]
    s = 1.0 + 0.5 + v' * G * v
    @test r.Δ ≈ G * v .* ((1.3 - 0.4 - 0.42) / s) && r.Gs ≈ G - (G * v) * (G * v)' ./ s
    @test r.Wbar ≈ sqrt(1 / s) * exp(-(1.3 - 0.82)^2 / (2s)) && r.zbar ≈ 0.82 && r.Vz ≈ s - 1
    # without genetic variance in the slope it is one trait's selection, and the lag is steady_lag's
    r1 = reaction_norm_response([0.2 0.0; 0.0 0.0], 0.7; ω = 1.0, VE = 0.5, θ = 1.3, abar = 0.4, bbar = 0.6)
    @test r1.Δ ≈ [0.2 * (1.3 - 0.82) / 1.7, 0.0]
    @test reaction_norm_lag([0.2 0.0; 0.0 0.0], 0.7; k = 0.02, B = 2.0, bbar = 0.6, ω = 1.0, VE = 0.5) ≈
          steady_lag(0.02; VG = 0.2, Vs = 1.5, b = 0.6, B = 2.0)
    # with a cost: Monte Carlo of the selected parents' genotypes
    rc = reaction_norm_response(G, 0.7; ω = 1.0, VE = 0.5, θ = 1.3, abar = 0.4, bbar = 0.6, ωb = 0.8)
    rng = MersenneTwister(7)
    C = QTLSim.psd_sqrt(G)
    n = 400_000
    X = C * randn(rng, 2, n) .+ [0.4, 0.6]
    z = X[1, :] .+ 0.7 .* X[2, :] .+ sqrt(0.5) .* randn(rng, n)
    W = exp.(-(z .- 1.3) .^ 2 ./ 2 .- X[2, :] .^ 2 ./ (2 * 0.64))
    @test mean(W) ≈ rc.Wbar rtol = 0.01
    @test mean(W .^ 2) ≈ rc.W2bar rtol = 0.01
    μs = X * W ./ sum(W)
    @test μs ≈ [0.4, 0.6] .+ rc.Δ atol = 0.003
    D = X .- μs
    @test (D .* W') * D' ./ sum(W) ≈ rc.Gs atol = 0.003
    @test rc.Gs ≈ G - G * rc.K * G

    # --- kernel against the moments lowering (8 runs): under a sustained move the slope
    # evolves toward B = 1 as the environment leaves the reference, and its variance is used
    # up. Tolerances: mean slope within 0.1, lag within 0.1 phenotypic SD (absolute), G_aa 25%
    q = Params(rn; PlasticityCostSD = Inf)
    rs = [simulate(q; seed) for seed in 1:8]
    m = run_moments(q; erosion = true)
    @test m.trait_mean[1000] ≈ [0.0, 0.2] atol = 1e-12                     # nothing moves the means before the move
    for t in (1200, 1500)
        @test mean(r.trait_mean[t][2] - r.trait_mean[1000][2] for r in rs) ≈ m.trait_mean[t][2] - 0.2 atol = 0.1
        @test mean(mean(x[1] for x in r.trait_lag[t-49:t]) for r in rs) ≈ mean(x[1] for x in m.trait_lag[t-49:t]) atol = 0.1
    end
    @test m.trait_mean[1500][2] > 0.7                                      # plasticity has evolved
    @test m.gmatrix[1500][2, 2] < 0.5 * m.gmatrix[1000][2, 2]
    @test mean(r.gmatrix[1500][2, 2] for r in rs) < 0.5 * mean(r.gmatrix[1000][2, 2] for r in rs)
    @test mean(r.gmatrix[1000][1, 1] for r in rs) ≈ m.gmatrix[1000][1, 1] rtol = 0.25
    # the quasi-steady lag of the moments, at their own G, slope and environment
    t = 1300
    @test m.trait_lag[t][1] ≈ reaction_norm_lag(m.gmatrix[t] .* (499 / 500), m.theta[t]; k = 0.02, bbar = m.trait_mean[t][2],
                                                ω = 1.0, VE = 0.5) rtol = 0.15
    # a cost holds the slope below B
    mc = run_moments(rn; erosion = true)
    @test 0.2 < mc.trait_mean[1500][2] < m.trait_mean[1500][2]
    # a step's own optimum and width, in the moments
    mθ = run_moments(Params(rn; MoveFitnessSDGen = -1, NoGenerations = 50); cycle = withparams(g, rn; θ = 0.5, ω = 0.8))
    @test all(==(0.0), mθ.theta) && mθ.trait_mean[end][1] > 0
    mo = run_moments(Params(rn; NoGenerations = 30); optimum = fill(0.4, 30))
    @test all(==(0.4), mo.theta)
end
