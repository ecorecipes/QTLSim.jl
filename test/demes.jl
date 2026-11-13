# Tests: demes. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl demes`.

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
