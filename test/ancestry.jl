# Tests: ancestry. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl ancestry`.

@testset "ancestry: block identity, pedigrees, coalescence times, tree-sequence tables (src/ancestry.jl)" begin
    L = lifecycle()
    base = Params(PopSize = 24, NoQTL = 12, Ploidy = 2, NoGenerations = 30, RecombPoissonMean = 1.0,
                  MutationRatePerGamete = 0.03, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01)
    samealleles(ped, w) = (s = samples(ped); G = genomes(w); P = size(G, 2);
        all(alleles(ped, l)[s.node[P*(i-1)+c]] === G[l, c, i] for l in axes(G, 1), c in 1:P, i in 1:w.N))

    @testset "recording changes nothing; off by default and allocation-free" begin
        alloc(w, rng) = @allocated generation!(w, rng)
        for p in (base, Params(base; Ploidy = 1), Params(base; CrossoverModel = :all, RecombPoissonMean = 3.0),
                  Params(base; Selection = :fecundity), Params(base; Demography = :beverton_holt, GrowthRate = 2.0),
                  Params(base; NoTraits = 2, MutCorr = 0.5))
            plain = World(p); a = run!(plain, MersenneTwister(4))
            @test plain.anc === nothing
            for every in (0, 7)
                w = World(p); ped = record_pedigree!(w; simplify_every = every)
                @test run!(w, MersenneTwister(4)) == a               # the same run, bit for bit
                @test genomes(w) == genomes(plain)
                @test check_tables(ped)                              # complete cover when never simplified
                @test (every == 0) == (ped.nsimplified == 0)
                @test pedigree_time(ped) == p.NoGenerations
                p.NoTraits == 1 ? (@test samealleles(ped, w)) :      # every allele rebuilt from the pedigree, bit for bit
                    (@test all(alleles(ped, l; trait = r)[samples(ped).node[2(i-1)+c]] === genomes(w)[r, l, c, i]
                               for r in 1:2, l in 1:p.NoQTL, c in 1:2, i in 1:w.N))
                QTLSim.check_blocks(w)
            end
            # without a pedigree a generation still allocates nothing
            w = World(Params(p; NoGenerations = 10^6)); rng = MersenneTwister(1); init_env!(w, rng)
            for _ in 1:300; generation!(w, rng); end
            @test minimum(alloc(w, rng) for _ in 1:50) == 0
        end
    end

    @testset "a pedigree checked by hand" begin
        # two diploids of four loci; founders' chromosomes are nodes 1, 2 (individual 1) and 3, 4
        w = World(Params(base; PopSize = 2, NoQTL = 4)); ped = record_pedigree!(w)
        a = w.anc
        @test a.nodes == [1 3; 2 4] && ped.node_time == zeros(4)
        # offspring 1: a gamete of parent 1 that starts on its chromosome 1 and switches after the
        # second locus, and chromosome 2 of parent 2; offspring 2: chromosome 1 of each parent
        QTLSim.newborn!(a, 1, 1, 2)
        QTLSim.inherit!(a, 1, 1, 1, 1, [1], 1); QTLSim.inherit!(a, 2, 1, 2, 2)
        QTLSim.newborn!(a, 2, 1, 2)
        QTLSim.inherit!(a, 1, 2, 1, 1); QTLSim.inherit!(a, 2, 2, 2, 1)
        swapbuffers!(w)
        @test collect(zip(ped.left, ped.right, ped.parent, ped.child)) ==
              [(0, 2, 1, 5), (2, 4, 2, 5), (0, 4, 4, 6), (0, 4, 1, 7), (0, 4, 3, 8)]
        @test samples(ped).node == [5, 6, 7, 8] && samples(ped).ind == [3, 3, 4, 4]
        @test (ped.ind_p1[3], ped.ind_p2[3]) == (1, 2)
        @test marginal_tree(ped, 2) == [0, 0, 0, 0, 1, 4, 1, 3]
        @test marginal_tree(ped, 3) == [0, 0, 0, 0, 2, 4, 1, 3]
        @test tmrca(ped, 1, 5, 7) == 1 && tmrca(ped, 2, 5, 7) == 1      # copies of node 1
        @test tmrca(ped, 3, 5, 7) === nothing && tmrca(ped, 1, 5, 6) === nothing
        pc = pair_coalescence(ped)
        @test pc.t == 1 && pc.nloci == 4
        @test pc.W == reshape([0, 8], 2, 1)               # the two copies in an individual: never the same founder
        @test pc.S == reshape([2, 14], 2, 1, 1)           # of 4 pairs × 4 loci, nodes 5 and 7 at loci 1 and 2
        @test coalescence_cdf(pc, 1).S == [2 / 16;;] && coalescence_mean(pc).S == [1.0;;]
        @test check_tables(ped)
        simplify!(ped)                                    # leaves the samples and node 1, over loci 1:2
        @test ped.node_time == [0, 1, 1, 1, 1] && samples(ped).node == [2, 3, 4, 5]
        @test collect(zip(ped.left, ped.right, ped.parent, ped.child)) == [(0, 2, 1, 2), (0, 2, 1, 4)]
        @test pair_coalescence(ped).S == pc.S && pair_coalescence(ped).W == pc.W
        @test check_tables(ped)
        @test_throws ErrorException check_tables(ped; complete = true)
        @test_throws ArgumentError record_pedigree!(w)    # already recording
        stop_recording!(ped); @test w.anc === nothing
    end

    @testset "simplification keeps every marginal tree" begin
        for p in (base, Params(base; CrossoverModel = :all, RecombPoissonMean = 3.0), Params(base; Ploidy = 1))
            w0 = World(p); full = record_pedigree!(w0); run!(w0, MersenneTwister(8))
            w1 = World(p); simp = record_pedigree!(w1; simplify_every = 4); run!(w1, MersenneTwister(8))
            a, b = pair_coalescence(full), pair_coalescence(simp)
            @test a.W == b.W && a.S == b.S
            @test sum(a.S) == p.NoQTL * (p.Ploidy == 2 ? 48 * 47 ÷ 2 - 24 : 24 * 23 ÷ 2)   # every pair once
            s0, s1 = samples(full).node, samples(simp).node
            @test all(tmrca(full, l, s0[i], s0[j]) == tmrca(simp, l, s1[i], s1[j])
                      for l in (1, 5, 12), i in 1:3:length(s0), j in 2:5:length(s0) if i != j)
            before = pedigree_bytes(simp)
            simplify!(simp)
            @test pedigree_bytes(simp) <= before < pedigree_bytes(full)
            c = pair_coalescence(simp)
            @test a.W == c.W && a.S == c.S && samealleles(simp, w1)
            @test length(simp.node_time) <= length(s1) + p.NoQTL * (length(s1) - 1)
        end
    end

    @testset "identity of blocks against the pedigree" begin
        w = World(base); @test ibd_sharing(w).deme == [1.0]             # one founding block of zeros
        @test tag_founders!(w) == 2 * 24 * w.pool.nb && ibd_sharing(w).deme == [0.0]
        QTLSim.check_blocks(w)
        # blocks of one locus, tagged founders: the same block iff identical by descent without mutation
        w = World(base; blocksize = 1); ped = record_pedigree!(w); run!(w, MersenneTwister(2))
        s = samples(ped)
        chr(n) = (i = (n - 1) ÷ 2 + 1; c = (n - 1) % 2 + 1; (i, c))
        nsame = 0
        for l in 1:base.NoQTL
            lab = pedigree_ibd(ped, l)
            for x in 1:length(lab), y in x+1:length(lab)
                same = shared_blocks(w, chr(x)..., chr(y)...)[l]
                nsame += same
                same == (lab[x] == lab[y]) || error("locus $l, chromosomes $x and $y: block $same, pedigree $(lab[x] == lab[y])")
            end
        end
        @test nsame > 100
        sh = ibd_sharing(w)
        @test sh.deme[1] ≈ nsame / (base.NoQTL * (48 * 47 ÷ 2)) rtol = 0.05   # ind pairs are 24 of 1128
        pc = pair_coalescence(ped)
        @test coalescence_cdf(pc, pc.t).S[1] >= sh.deme[1]                   # mutation only lowers block identity
        # default blocks: the same block implies identity by descent without mutation at every locus of it
        w = World(base); ped = record_pedigree!(w); run!(w, MersenneTwister(3))
        pool = w.pool
        nsame = nibd = 0
        for k in 1:pool.nb, l in (k-1)*pool.B+1:min(k * pool.B, pool.L)
            lab = pedigree_ibd(ped, l; block = (k-1)*pool.B+1:min(k * pool.B, pool.L))
            for x in 1:length(lab), y in x+1:length(lab)
                same = shared_blocks(w, chr(x)..., chr(y)...)[k]
                nsame += same; nibd += lab[x] == lab[y]
                same && lab[x] != lab[y] && error("block $k, locus $l: the same block without identity by descent")
            end
        end
        @test 0 < nsame <= nibd
        # after simplification a mutation is kept on the lineage of its own locus only
        simplify!(ped)
        @test_throws ArgumentError pedigree_ibd(ped, 1; block = 1:pool.B)
        @test pedigree_ibd(ped, 1) isa Vector{Int32}
        @test ibd_fraction(w, 1, 1, 1, 1) == 1.0
        @test_throws ArgumentError shared_blocks(w, 1, 1, World(base), 1, 1)
    end

    @testset "metapopulations: migration, demes, the categorical layer" begin
        k = 3
        p = Params(base; PopSize = 10, Drift = true)
        cycle = foldl(⊗, fill(L.generation, k)) ⋅ L.migrate(k)
        plain = p ⊗ demes(k)
        a = run_composed(plain, MersenneTwister(6); cycle, migration_rate = 0.2)
        s = p ⊗ demes(k); ped = record_pedigree!(s; simplify_every = 0)
        @test run_composed(s, MersenneTwister(6); cycle, migration_rate = 0.2) == a
        @test all(genomes(s.demes[d]) == genomes(plain.demes[d]) for d in 1:k)
        @test length(ped.mig_ind) > 0 && iseven(length(ped.mig_ind))
        smp = samples(ped)
        @test ped.ind_deme[smp.ind] == smp.deme                    # everyone is where its migrations took it
        @test all(ped.ind_time[smp.ind] .== p.NoGenerations)
        # a migrant's deme history: born in ind_birth, moved by its migrations in order
        where = copy(ped.ind_birth)
        for e in eachindex(ped.mig_ind)
            @test where[ped.mig_ind[e]] == ped.mig_from[e]
            where[ped.mig_ind[e]] = ped.mig_to[e]
        end
        @test where == ped.ind_deme
        @test check_tables(ped)
        pc = pair_coalescence(ped)
        @test size(pc.S) == (p.NoGenerations + 1, k, k) && pc.S[:, 1, 2] == pc.S[:, 2, 1]
        @test sum(pc.S[:, 1, 2]) == p.NoQTL * 20 * 20 && sum(pc.W[:, 3]) == p.NoQTL * 10
        @test all(samealleles_d -> samealleles_d, [begin
                      G = genomes(s.demes[d]); n = smp.node[smp.deme .== d]
                      all(alleles(ped, l)[n[2(i-1)+c]] === G[l, c, i] for l in 1:p.NoQTL, c in 1:2, i in 1:10)
                  end for d in 1:k])
        # the same run by hand, generation! on each deme and then migrate!, with simplification
        s2 = p ⊗ demes(k); ped2 = record_pedigree!(s2; simplify_every = 5)
        rng = MersenneTwister(6); foreach(w -> init_env!(w, rng), s2.demes)
        for _ in 1:p.NoGenerations
            foreach(w -> generation!(w, rng), s2.demes)
            migrate!(s2, rng; rate = 0.2)
        end
        pc2 = pair_coalescence(ped2)
        @test ped2.nsimplified == 6 && pc2.W == pc.W && pc2.S == pc.S
        @test samples(ped2).deme == smp.deme
        # identity by descent among demes, from blocks and from the pedigree
        f = ibd_fst(s); cdf = coalescence_cdf(pc, pc.t)
        g = ibd_fst([cdf.S[d, d] for d in 1:k], cdf.S)
        @test 0 < f.fst < 1 && 0 < g.fst < 1 && f.f_within <= g.f_within && f.f_between <= g.f_between
        @test ibd_sharing(s).between == ibd_sharing(s.demes).between
        @test ibd_fst([0.5, 0.5], [0.0 0.1; 0.1 0.0]).fst ≈ (0.5 - 0.3) / (1 - 0.3)
        # migration between a recorded and an unrecorded population is refused
        v = World(p; pool = s.demes[1].pool)
        @test_throws ArgumentError QTLSim.exchange_genomes!(s.demes[1], 1, v, 1)
        # a copy of a recording World records nothing
        @test all(w -> w.anc === nothing, (s.demes[1] ⊗ demes(2)).demes)
    end

    @testset "tree-sequence tables: text round trip, the same trees" begin
        w = World(base); ped = record_pedigree!(w; simplify_every = 10); run!(w, MersenneTwister(9))
        for simplified in (false, true)
            simplified && simplify!(ped)
            tb = tskit_tables(ped)
            @test check_tables(tb)
            dir = mktempdir()
            write_tskit(dir, ped)
            back = read_tskit(dir)
            @test back == tb && check_tables(back)
            @test readline(joinpath(dir, "nodes.txt")) == "id\tis_sample\ttime\tpopulation\tindividual"
            @test readline(joinpath(dir, "edges.txt")) == "id\tleft\tright\tparent\tchild"
            @test readlines(joinpath(dir, "populations.txt"))[2] == "0\teyJuYW1lIjoiZGVtZV8xIn0="   # {"name":"deme_1"}
            s = samples(ped).node
            @test sum(tb.nodes.is_sample) == length(s) && all(tb.nodes.time[s] .== 0)
            @test all(table_tmrca(back, l, s[i] - 1, s[j] - 1) == tmrca(ped, l, s[i], s[j])
                      for l in (1, 6, 12), i in 1:5:length(s), j in 3:7:length(s) if i != j)
            # a mutation's derived state is the allele of its node after it
            x = alleles(ped, Int(tb.sites.position[1]) + 1)
            last = Dict(tb.mutations.node[m] => tb.mutations.derived_state[m] for m in eachindex(tb.mutations.site) if tb.mutations.site[m] == 0)
            @test all(parse(Float64, v) === x[n+1] for (n, v) in last)
        end
        # broken tables are caught
        tb = tskit_tables(ped)
        bad = merge(tb, (edges = merge(tb.edges, (parent = reverse(tb.edges.parent),)),))
        @test_throws ErrorException check_tables(bad)
        @test QTLSim._unbase64(QTLSim._base64("{\"name\":\"deme_12\"}")) == "{\"name\":\"deme_12\"}"
        @test QTLSim._base64("Man") == "TWFu" && QTLSim._base64("Ma") == "TWE=" && QTLSim._base64("M") == "TQ=="
    end
end
