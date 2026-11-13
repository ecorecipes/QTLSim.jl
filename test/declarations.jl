# Tests: declarations. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl declarations`.

@testset "closures declared as data: :auto read from the declarations is the hand-written rule" begin
    L = lifecycle(); g = L.generation
    # the rule as written before the declarations (moments.jl's and trait_moments.jl's
    # resolve_loci, QTLSim.moments_loci at 45280fb), kept here as the reference
    legacy_loci(p) = p.NoTraits > 1 ? :rare_local : :rare_alleles
    function legacy_one(ps, cycle)
        mf = Ext.mutate_fields(cycle)
        mf === nothing && return :multilocus
        qs = isempty(mf) ? ps : [Params(q; mf...) for q in ps]
        rs = unique([legacy_loci(q) for q in qs])
        length(rs) == 1 && rs[1] in (:rare_alleles, :rare_local) || return :multilocus
        l = rs[1]
        l === :rare_alleles && length(ps) > 1 && (l = :rare_local)
        l
    end
    function legacy_traits(ps, cycle)
        ps[1].ReactionNorm && return :multilocus
        fields = (:MutationScalar, :MutationRatePerGamete, :MutCorr, :MutCov)
        fs = unique([NamedTuple(k => v for (k, v) in pairs(Ext.stepparams(x)) if k in fields)
                     for x in Ext.stepsof(cycle) if !(x isa Catlab.ObExpr) && Ext.stepname(x) === :mutate])
        length(fs) > 1 && return :multilocus
        qs = isempty(fs) ? ps : [Params(q; only(fs)...) for q in ps]
        rs = unique([legacy_loci(q) for q in qs])
        length(rs) == 1 && rs[1] in (:rare_alleles, :rare_local) || return :multilocus
        rs[1]
    end
    base = Params(PopSize = 30, NoQTL = 5, NoGenerations = 4, MutationScalar = 0.01, MutationRatePerGamete = 0.01)
    models = [base, Params(base; NoTraits = 2, MutCorr = 0.3), Params(base; NoTraits = 3),
              Params(base; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2)]
    n = 0
    # run_moments resolves :auto without `between` (a haploid metapopulation with between = true
    # is rejected afterwards, by check_between, as before)
    for p0 in models, ploidy in (1, 2), k in (1, 2, 4), mut in (:none, :same, :different)
        p = Params(p0; Ploidy = ploidy)
        same = withparams(g, p; MutationScalar = 0.02)
        other = withparams(g, p; MutationScalar = 0.03)
        cycle = k == 1 ? (mut === :none ? g : mut === :same ? same ⋅ same : same ⋅ other) :
                (mut === :none ? foldl(⊗, fill(g, k)) : mut === :same ? foldl(⊗, fill(same, k)) :
                 foldl(⊗, [same; fill(other, k - 1)])) ⋅ L.migrate(k)
        ps = fill(p, k)
        new = p.NoTraits == 1 ? Ext.resolve_loci(ps, cycle, :auto)[1] : Ext.resolve_trait_loci(ps, cycle, :auto)[1]
        old = p.NoTraits == 1 ? legacy_one(ps, cycle) : legacy_traits(ps, cycle)
        @test new === old
        n += 1
        # the Params-only rule (QTLSim.moments_loci, multilocus_equilibrium's default); the
        # legacy rule gave :rare_local for a heritable reaction norm, fixed below
        @test QTLSim.moments_loci(p) === (p.ReactionNorm ? :multilocus : legacy_loci(p))
    end
    @test n == 4 * 2 * 3 * 3
    # regression: a heritable reaction norm, whose moments run with :multilocus, had
    # moments_loci = :rare_local (support matrix, 2026-11-11); now the declarations decide
    rn = models[4]
    @test QTLSim.auto_closure(QTLSim.closure_context(rn)) === :multilocus
    @test QTLSim.moments_loci(rn) === :multilocus
    @test Ext.resolve_trait_loci([rn], g, :auto)[1] === :multilocus
    # so multilocus_equilibrium's default is the multilocus model; without drift its elevation
    # entry is the reaction norm's run_moments at ε = 0 (a fixed optimum at 0)
    rq = Params(rn; PopSize = 1000, NoQTL = 20, NoGenerations = 3000, MoveType = :continuous, MoveFitnessSD = 0.0)
    eq = multilocus_equilibrium(rq)
    @test eq.VG == multilocus_equilibrium(rq; loci = :multilocus).VG
    @test eq.VG[1, 1] ≈ run_moments(rq).gmatrix[end][1, 1] rtol = 1e-9
    # the declarations: ranks increase along the order, every arrow points to a coarser closure
    cs = QTLSim.closures()
    @test [c.name for c in cs] == [:aggregate, :multilocus, :rare_alleles, :rare_local, :spectrum]
    @test all(QTLSim.closure(r).rank < c.rank for c in cs for r in c.refines)
    @test_throws ArgumentError QTLSim.closure(:gaussian)
    @test QTLSim.auto_closure(QTLSim.closure_context(base; ndemes = 3)) === :rare_local
    @test QTLSim.auto_closure(QTLSim.closure_context(base; uniform_mutation = false)) === :multilocus
end

@testset "closures declared as data: the domains are the lowering's guards" begin
    L = lifecycle(); g = L.generation
    base = Params(PopSize = 30, NoQTL = 5, NoGenerations = 4, MutationScalar = 0.01, MutationRatePerGamete = 0.01)
    models = [base, Params(base; NoTraits = 2, MutCorr = 0.3), Params(base; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2)]
    # regression: several traits in several demes accepted loci = :spectrum (its guard was only on
    # the one-deme path) and ran the rare-alleles closure; every probe now runs exactly in the domain
    for p0 in models, ploidy in (1, 2), k in (1, 3), between in (false, true), uniform in (true, false)
        p = Params(p0; Ploidy = ploidy)
        gm = withparams(g, p; MutationScalar = 0.02)
        cycle = k == 1 ? (uniform ? g : g ⋅ gm) : (uniform ? foldl(⊗, fill(g, k)) : gm ⊗ g ⊗ g) ⋅ L.migrate(k)
        ctx = QTLSim.closure_context(p; ndemes = k, uniform_mutation = uniform, between)
        for c in QTLSim.closures()
            runs = try
                k == 1 ? run_moments(p; cycle, loci = c.name, between) :
                         run_moments(fill(p, k); cycle, loci = c.name, between, migration_rate = 0.01)
                true
            catch e
                e isa ArgumentError || rethrow()
                false
            end
            @test runs == c.domain(ctx)
            c.auto(ctx) && c.name !== :multilocus && @test c.domain(ctx)     # the automatic domain is in the domain
        end
    end
    # every path into the spectrum with matrices is refused: several traits in several demes,
    # and the trait recursion called directly, with one trait or several, in one deme
    p2 = models[2]
    @test_throws ArgumentError run_moments(fill(p2, 3); cycle = foldl(⊗, fill(g, 3)) ⋅ L.migrate(3), loci = :spectrum,
                                           migration_rate = 0.01)
    @test_throws ArgumentError run_moments(p2; cycle = foldl(⊗, fill(g, 3)) ⋅ L.migrate(3), loci = :spectrum)
    @test_throws ArgumentError Ext.run_trait_moments(p2; loci = :spectrum)
    @test_throws ArgumentError Ext.run_trait_moments(base; loci = :spectrum)
end

@testset "declared model morphisms: generated reduction tests" begin
    ms = vcat(model_morphisms(), QTLSim.closure_order())
    @test length(ms) == 8 && all(m -> !isempty(m.reductions), ms)
    @test all(m -> m.kind in (:embedding, :projection, :restriction, :conjugacy, :refinement), ms)
    for m in ms, r in m.reductions
        res = QTLSim.check_reduction(r)
        @test res.pass
        res.pass || @info "reduction failed" m.name r.name res.difference
    end
    # every closure but the coarsest reduces to another along the declared order
    arrows = [(m.from, m.to) for m in QTLSim.closure_order()]
    @test all(c -> c.name === :aggregate || any(a -> a[1] == ":$(c.name)", arrows), QTLSim.closures())
end

@testset "naturality harness: comparisons, the square of each step, a case" begin
    S = QTLSim.Statistic("x", out -> out.x)
    outs = [(x = v,) for v in (1.0, 1.2, 0.8, 1.1, 0.9)]
    r = QTLSim.compare_stat(:band, S, 0.05, :kernel, :moments, outs, [(x = 1.0,)])
    @test r.pass && r.from ≈ 1.0 && abs(r.measure) < 1e-12
    @test !QTLSim.compare_stat(:band, S, 0.05, :kernel, :moments, outs, [(x = 2.0,)]).pass
    @test !QTLSim.compare_stat(:band, S, 0.05, :kernel, :moments, [(x = NaN,), (x = 1.0,)], [(x = 1.0,)]).pass  # NaN fails
    @test QTLSim.compare_stat(:z, S, 4.0, :kernel, :moments, outs, [(x = 1.05,)]).pass
    @test QTLSim.compare_stat(:exact, S, 1e-12, :moments, :theory, [(x = [1.0, 2.0],)], [(x = [1.0, 2.0],)]).pass
    @test !QTLSim.compare_stat(:exact, S, 1e-12, :moments, :theory, [(x = [1.0, 2.0],)], [(x = [1.0, 2.1],)]).pass
    # the KS p-value against HypothesisTests' asymptotic value for D = 0.5, n = m = 10: P(K > √5·0.5)
    @test QTLSim.ks_pvalue(1:10, 6:15) ≈ QTLSim.ccdf(QTLSim.Kolmogorov(), sqrt(5) * 0.5) rtol = 1e-12
    @test QTLSim.ks_pvalue(randn(Xoshiro(1), 200), randn(Xoshiro(2), 200)) > 0.01
    # a derived statistic: combine the parts' means, SE by bootstrap
    D = QTLSim.Statistic("ratio", [QTLSim.Statistic("a", o -> o.a), QTLSim.Statistic("b", o -> o.b)], (a, b) -> a / b)
    e = QTLSim.replicate_estimate(D, :kernel, [(a = 1.0, b = 2.0), (a = 3.0, b = 2.0)])
    @test e.value == 1.0 && e.se > 0
    @test_throws ArgumentError QTLSim.check_case(NaturalityCase(name = "x", params = Params(), from = :moments, to = :kernel,
                                                                statistics = [S => 0.1]))
    @test_throws ArgumentError QTLSim.check_case(NaturalityCase(name = "x", params = Params(), criterion = :ks,
                                                                statistics = [S => 0.1]))
    # the exact square: the neutral moments are deme_coalescence's recursion read forward
    L = lifecycle(); g = L.generation
    pn = Params(Ploidy = 2, NoQTL = 20, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.002,
                MutationRatePerGamete = 0.05, Drift = true, PopSize = 30, NoGenerations = 60, Sample = 1)
    M = island_migration(3, 0.02)
    th(c) = (x = neutral_deme_variances(backward_migration(M, 30), 30, muteVar(pn); t = collect(1:59), origins = :exchange);
             (gen = collect(1:59), within = [sum(y.within) / 3 for y in x]))
    W = QTLSim.Statistic("V_W", o -> error(); moments = o -> o.metapop.within[2:60], theory = o -> o.within)
    rc = check_naturality(foldl(⊗, fill(g, 3)) ⋅ L.migrate(3), fill(pn, 3); from = :moments, to = :theory, criterion = :exact,
                          statistic = W, tolerance = 1e-10, replicates = 1,
                          options = (migration_rate = M, moments = (between = true,), theory = th))
    @test rc.pass && rc.reproduced && rc.rows[1].measure < 1e-10
    @test occursin("reproduced", sprint(show, MIME"text/plain"(), rc))
    # a stochastic square: the kernel's V_G against the multilocus moments (a short, loose case)
    pk = Params(PopSize = 200, NoQTL = 10, Ploidy = 2, NoGenerations = 1500, PhenotypeVar = 1.0, Heritability = 0.5,
                MutationScalar = 0.01, MutationRatePerGamete = 0.01, Sample = 10)
    rk = check_naturality(nothing, pk; statistic = QTLSim.window_mean(:geno_var, 1000:10:1490), tolerance = 0.25, replicates = 6,
                          options = (moments = (loci = :multilocus, erosion = true),), seed = 3)
    @test rk.pass && length(rk.values[:kernel]) == 6 && haskey(rk.values[:kernel][1], :geno_mean)
    # the square of each step: the aggregate recursion inherits all of selection's reduction,
    # so the square of reproduce fails, and develop's and select's hold; the multilocus model's all hold
    ps = Params(PopSize = 1000, NoQTL = 100, Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.02,
                MutationRatePerGamete = 0.01, InitFitnessSD = 4.0, RecombPoissonMean = 3.0)
    agg = naturality_steps(ps; loci = :aggregate, replicates = 60, burnin = 2000, seed = 1, tolerance = 0.1)
    f = QTLSim.first_failure(agg)
    @test f !== nothing && f.step === :reproduce && f.observable === :VG
    @test all(r -> r.pass, filter(r -> r.step in (:develop, :select), agg))
    ml = naturality_steps(ps; loci = :multilocus, replicates = 60, burnin = 2000, seed = 1, tolerance = 0.1)
    @test QTLSim.first_failure(ml) === nothing
    @test [r.step for r in ml if r.observable === :VG] == [:develop, :select, :reproduce, :record, :rescale, :mutate, :swap]
    @test_throws ArgumentError naturality_steps(Params(ps; NoTraits = 2))
end

@testset "support matrix and explain (generated from the probes and the declarations)" begin
    rows = support_matrix()
    @test length(rows) == 20 * 3 * 7
    @test all(r -> r.status in (:supported, :untested, :unsupported, :na, :error), rows)
    @test only(r for r in rows if r.feature == "QTSim's generation" && r.model == "one trait" && r.column === :kernel).status === :supported
    @test only(r for r in rows if r.feature == "QTSim's generation" && r.model == "reaction norm" && r.column === :slim).status === :unsupported
    # the committed page is the generated one (scripts/support_matrix.jl rewrites it)
    page = joinpath(pkgdir(QTLSim), "docs", "src", "support.md")
    @test isfile(page) && read(page, String) == support_markdown(rows)
    p = Params(PopSize = 50, NoQTL = 5, Ploidy = 2, NoGenerations = 10)
    e = explain(p; io = nothing)
    @test e.closure === :rare_alleles && e.model == OneTraitModel() && all(b -> b.status === :ok, e.backends)
    @test explain(fill(p, 3); migration_rate = 0.01, io = nothing).closure === :rare_local
    er = explain(Params(p; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2); io = nothing)
    @test er.closure === :multilocus && er.backends[2].status === :unsupported
    @test occursin("closure: :rare_alleles", sprint(io -> explain(p; io)))
end
