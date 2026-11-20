# Tests: laws. The targets of the lowerings as symmetric monoidal categories (roadmap C1,
# ext/QTLSimCatlabExt/typed.jl): images typed and checked when lowered, the GATlab
# instances, and the laws on random composites — exactly where an image is deterministic
# (the moments lowering, Eidos text, and the kernel when the order of random draws is the
# same), in law where it is not (the kernel under interchange and the braiding's
# naturality: Lean's `kernelProd_interchange`). Included by runtests.jl; run alone with
# `julia --project=vignettes test/runtests.jl laws`.

const LAWS = lifecycle()
const LG = LAWS.generation
const LP = Params(PopSize = 60, NoQTL = 10, Ploidy = 2, NoGenerations = 30, PhenotypeVar = 1.0, Heritability = 0.5,
                  MutationScalar = 0.05, MutationRatePerGamete = 0.02, InitFitnessSD = 2.0, Sample = 1)
const LZ = dom(LG)

# the images of the steps of `n` generations, as typed morphisms, in order
step_morphs(target; n = 1) = [lower(LAWS[s], target; typed = true) for _ in 1:n for s in Ext.STEPS]

# a random bracketing of a sequence of morphisms, with identities inserted at random
function random_composite(ms, rng)
    length(ms) == 1 && return rand(rng) < 0.3 ? Ext.morph_compose(Ext.morph_id(ms[1].dom), ms[1]) : ms[1]
    k = rand(rng, 1:length(ms)-1)
    Ext.morph_compose(random_composite(ms[1:k], rng), random_composite(ms[k+1:end], rng))
end

# a kernel population after its environment is set, and the state a test compares
fresh_world(seed) = (w = World(LP); init_env!(w, Xoshiro(seed)); w)
kstate(w) = (genomes(w), w.genotype[1:w.N], w.env[1:w.N], w.θ, w.N, w.nsurv, w.extinct)
function run_image(op, w, seed, ngen)
    rng = Xoshiro(seed)
    for t in 0:ngen-1
        w.extinct && break
        w.gen = t
        op(w, rng)
    end
    kstate(w)
end
# every field of a moments state but its parameters
mstate(m) = Tuple(getfield(m, f) for f in fieldnames(typeof(m)) if f ∉ (:p, :cov))

@testset "images are typed and checked when the life cycle is lowered" begin
    for t in (:kernel, :eidos, :moments)
        m = lower(LG, t; typed = true)
        @test m isa Ext.Morph{t}
        @test m.dom.names == m.codom.names == [:Zygotes]
        # the default return is the image, as before
        @test typeof(lower(LG, t)) == typeof(m.image)
    end
    k2 = (LG ⊗ LG) ⋅ LAWS.migrate(2)
    @test lower(k2, :kernel; typed = true).dom.names == [:Zygotes, :Zygotes]
    # wrong calling convention, wrong number of demes, wrong image type: caught at lowering
    @test_throws ArgumentError lower(LG, :kernel; images = Dict(:mutate => w -> w))
    @test_throws ArgumentError lower(k2, :kernel; images = Dict(:migrate => Ext.DemeOp(3, [1, 2, 3], (s, d, r) -> s)))
    @test_throws ArgumentError lower(LG, :kernel; images = Dict(:mutate => Ext.EidosCode(:mutate, "")))
    @test_throws ArgumentError lower(LG, :eidos; images = Dict(:mutate => Ext.Operator(:mutate, (w, r) -> w)))
    # composites of morphisms whose stages do not meet
    dev, sel = lower(LAWS.develop, :kernel; typed = true), lower(LAWS.select, :kernel; typed = true)
    @test Ext.morph_compose(dev, sel).codom.names == [:Survivors]
    @test_throws ArgumentError Ext.morph_compose(sel, dev)
end

@testset "the targets are models of the theory of symmetric monoidal categories" begin
    for t in (:kernel, :eidos, :moments)
        M = Ext.LoweringCategory{t}()
        @test Catlab.GATlab.implements(M, Catlab.ThSymmetricMonoidalCategory)
        a, b = Ext.Stage{t}(:Zygotes), Ext.Stage{t}(:Adults)
        f = lower(LAWS.develop, t; typed = true)
        @test Catlab.compose[M](Catlab.id[M](a), f).codom.names == [:Adults]
        @test Catlab.otimes[M](a, b).names == [:Zygotes, :Adults]
        @test Catlab.braid[M](a, b).codom.names == [:Adults, :Zygotes]
        @test Catlab.munit[M]().names == Symbol[]
        @test Catlab.Hom[M](f, a, b) === f
        @test_throws ArgumentError Catlab.Hom[M](f, b, a)
    end
end

@testset "associativity, units and functoriality, exactly, on random composites" begin
    rng = MersenneTwister(20261113)
    for trial in 1:8
        n = rand(rng, 1:3)
        whole_k = lower(foldl(⋅, fill(LG, n)), :kernel)
        # the kernel: a random bracketing of the steps' images, with identities, acts on a
        # population exactly as the image of the whole composite (the same draws in the same order)
        comp = random_composite(step_morphs(:kernel; n), rng)
        @test comp.dom.names == comp.codom.names == [:Zygotes]
        seed = 1000trial
        @test run_image(comp.image, fresh_world(seed), seed + 1, 4) == run_image(whole_k, fresh_world(seed), seed + 1, 4)
        # Eidos: the same text
        ce = random_composite(step_morphs(:eidos; n), rng)
        @test ce.image.code == lower(foldl(⋅, fill(LG, n)), :eidos).code
        # moments: the same state
        cm = random_composite(step_morphs(:moments; n), rng)
        m1, m2 = Ext.Moments(LP; VG = 0.3), Ext.Moments(LP; VG = 0.3)
        whole_m = lower(foldl(⋅, fill(LG, n)), :moments)
        for t in 0:5
            m1.gen = m2.gen = t
            cm.image(m1, nothing); whole_m(m2, nothing)
        end
        @test isequal(mstate(m1), mstate(m2))
    end
end

# two demes with their own optima, so that they differ and a law cannot hold by symmetry
const LG1, LG2, LG3, LG4 = (withparams(LG; θ = x) for x in (0.5, -0.5, 1.0, -1.0))
two(p = LP) = [p, Params(p; InitFitnessSD = 3.0)]

# the records of a run, sorted by generation and deme (a braiding changes the order of rows, not their content)
function sorted_rows(o)
    ix = sortperm(collect(zip(o.gen, o.deme)))
    (o.gen[ix], o.deme[ix], o.pheno_mean[ix], o.geno_var[ix], o.nsurv[ix])
end

@testset "braiding: an involution, and natural (exactly in the moments lowering)" begin
    for t in (:kernel, :eidos, :moments)
        a = Ext.Stage{t}([:Zygotes]); b = Ext.Stage{t}([:Zygotes, :Zygotes])
        σσ = Ext.morph_compose(Ext.morph_braid(a, b), Ext.morph_braid(b, a))
        @test σσ.image.perm == 1:3
        @test σσ.dom.names == σσ.codom.names
    end
    s = Catlab.braid(LZ, LZ)
    # (G1 ⊗ G2) ⋅ σ and σ ⋅ (G2 ⊗ G1): deme 1 meets G1 and deme 2 meets G2 either way
    c1, c2 = (LG1 ⊗ LG2) ⋅ s, s ⋅ (LG2 ⊗ LG1)
    for loci in (:aggregate, :multilocus)
        @test sorted_rows(run_moments(two(); cycle = c1, loci)) == sorted_rows(run_moments(two(); cycle = c2, loci))
    end
end

@testset "interchange, exactly in the moments lowering" begin
    c1 = (LG1 ⊗ LG2) ⋅ (LG3 ⊗ LG4)
    c2 = (LG1 ⋅ LG3) ⊗ (LG2 ⋅ LG4)
    for loci in (:aggregate, :multilocus)
        @test sorted_rows(run_moments(two(); cycle = c1, loci)) == sorted_rows(run_moments(two(); cycle = c2, loci))
    end
    # with migration after the generation, too
    m = LAWS.migrate(2)
    @test sorted_rows(run_moments(two(); cycle = c1 ⋅ m, migration_rate = 0.05, loci = :multilocus)) ==
          sorted_rows(run_moments(two(); cycle = c2 ⋅ m, migration_rate = 0.05, loci = :multilocus))
end

# In the kernel, interchange and the braiding's naturality change the order in which the
# demes draw their random numbers, so they hold in law (Lean: `kernelProd_interchange`),
# not run for run. Each statistic of each deme at the last census is compared between the
# two composites over independent replicates: Welch's z, |z| < 4 (fixed seeds, so
# deterministic).
function law_z(c1, c2; R = 150, ngen = 20, stats = (:geno_var, :pheno_mean))
    p = Params(LP; NoGenerations = ngen)
    vals(c, off) = [begin
            o = run_composed(two(p) ⊗ demes(2), Xoshiro(off + r); cycle = c)
            last_gen = maximum(o.gen)
            [getproperty(o, st)[findfirst(i -> o.gen[i] == last_gen && o.deme[i] == d, eachindex(o.gen))]
             for st in stats for d in 1:2]
        end for r in 1:R]
    a, b = reduce(hcat, vals(c1, 0)), reduce(hcat, vals(c2, 10^6))
    [(mean(a[i, :]) - mean(b[i, :])) / sqrt(var(a[i, :]) / R + var(b[i, :]) / R) for i in axes(a, 1)]
end

@testset "interchange and the braiding's naturality hold in law in the kernel" begin
    s = Catlab.braid(LZ, LZ)
    z1 = law_z((LG1 ⊗ LG2) ⋅ (LG3 ⊗ LG4), (LG1 ⋅ LG3) ⊗ (LG2 ⋅ LG4))
    @test all(abs.(z1) .< 4)
    z2 = law_z((LG1 ⊗ LG2) ⋅ s, s ⋅ (LG2 ⊗ LG1))
    @test all(abs.(z2) .< 4)
    # and the test has power: two composites that differ in law are told apart. Sized so
    # that the expected |z| is about 7-12, not just over 4: at 20 generations and 150
    # replicates it was 4.6 on Julia 1.12 and fell below 4 on 1.10's random streams.
    z3 = law_z((LG1 ⊗ LG2) ⋅ (LG3 ⊗ LG4), (LG3 ⊗ LG4) ⋅ (LG3 ⊗ LG4); stats = (:pheno_mean,), R = 300, ngen = 60)
    @test any(abs.(z3) .> 4)
end
