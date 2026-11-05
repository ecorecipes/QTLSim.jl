# The lowering of the life cycle to SLiM's Eidos language, and a complete SLiM model built
# around the image of a life cycle.
#
# Each step becomes a fragment of Eidos. The fragments implement QTSim's single-locus
# diploid model with continuum-of-alleles genetics, statement for statement as the
# hand-written port slim_port/qtsim.slim does: each individual carries its two allele
# values (`a0`, `a1`) and SLiM supplies only the nonWF life cycle. Composition is
# concatenation inside one reproduction() callback per tick. A step that ends the run
# (`select`, when fewer than two survive) returns from that callback, which skips the rest
# of the composite as `⋅` does in the kernel.

"""
    EidosCode(name, code)

A fragment of Eidos, the image of a life-cycle morphism under the lowering to SLiM. The
image of a composite is the concatenation of the images of its steps.
"""
struct EidosCode
    name::Symbol
    code::String
end

compose(a::EidosCode, b::EidosCode) = EidosCode(Symbol(a.name, :⋅, b.name), a.code * b.code)
compose(a::EidosCode, b::EidosCode, c::EidosCode, rest::EidosCode...) =
    foldl(compose, (a, b, c, rest...))
⋅(a::EidosCode, b::EidosCode) = compose(a, b)

Base.show(io::IO, c::EidosCode) = print(io, "EidosCode(:", c.name, ")")
Base.show(io::IO, ::MIME"text/plain", c::EidosCode) = print(io, c.code)

"What Eidos code at each stage of the life cycle can rely on."
const EIDOS_SCOPE = Dict(
    :Zygotes => "the individuals of p1 carry allele values a0 and a1; globals THETA and OMEGA",
    :Adults => "adults, a0, a1, genotype, env, pheno",
    :Survivors => "W, survived, sidx and ns ≥ 2",
    :Offspring => "c0 and c1, the alleles of the N offspring")

"The Eidos image of each step of the life cycle, by name."
eidos_images() = Dict{Symbol,Any}(
    :develop => EidosCode(:develop, """
    // develop: genotype a0 + a1, phenotype genotype + e with e ~ N(0, V_E)   [Zygotes → Adults]
    adults = p1.individuals;
    a0 = adults.getValue("a0");
    a1 = adults.getValue("a1");
    genotype = a0 + a1;
    env = rnorm(size(adults), 0.0, ESD);
    pheno = genotype + env;
"""),
    :select => EidosCode(:select, """
    // select: move the optimum, Gaussian fitness with QTSim's 0.4, viability   [Adults → Survivors]
    if (moveGen >= 0 & gen >= moveGen) {
        if (moveType == "phenotype") defineGlobal("THETA", mean(pheno) + sqrt(VP) * moveSD);
        else if (moveType == "fitness") defineGlobal("THETA", mean(pheno) + OMEGA * moveSD);
        else defineGlobal("THETA", THETA + moveSD);
    }
    x = standardise ? (pheno - THETA) / sd(pheno) else pheno - THETA;
    W = drift ? rep(1.0, size(adults)) else dnorm(x, 0.0, OMEGA) / (0.4 / OMEGA);
    survived = runif(size(adults)) < W;
    sidx = which(survived);
    ns = size(sidx);
    if (ns < 2) { catn("Pop all died."); sim.simulationFinished(); self.active = 0; return; }
"""),
    :reproduce => EidosCode(:reproduce, """
    // reproduce: N offspring of two distinct survivors, a random allele from each   [Survivors → Offspring]
    sa0 = a0[sidx]; sa1 = a1[sidx];
    i1 = rdunif(N, 0, ns - 1);
    i2 = rdunif(N, 0, ns - 1);
    dup = which(i1 == i2);
    while (size(dup) > 0) { i2[dup] = rdunif(size(dup), 0, ns - 1); dup = which(i1 == i2); }
    c0 = ifelse(runif(N) < 0.5, sa1[i1], sa0[i1]);
    c1 = ifelse(runif(N) < 0.5, sa1[i2], sa0[i2]);
"""),
    :record => EidosCode(:record, """
    // record: the parent generation's statistics, every `sampleEvery` generations   [Offspring → Offspring]
    if (gen == 0) catn("gen\\tphenoMean\\tgenoVar\\tnSurv\\ttheta\\tomega");
    if (gen % sampleEvery == 0)
        catn(gen + "\\t" + mean(pheno) + "\\t" + var(genotype) + "\\t" + ns + "\\t" + THETA + "\\t" + OMEGA);
"""),
    :rescale => EidosCode(:rescale, """
    // rescale: the scheduled change of the fitness width, from the next generation   [Offspring → Offspring]
    if (gen == changeGen) defineGlobal("OMEGA", changeSD);
"""),
    :mutate => EidosCode(:mutate, """
    // mutate: Poisson(rate) increments of N(0, MUTVAR) per offspring chromosome   [Offspring → Offspring]
    k0 = rpois(N, rate);
    k1 = rpois(N, rate);
    c0 = c0 + rnorm(N, 0.0, sqrt(k0 * MUTVAR));
    c1 = c1 + rnorm(N, 0.0, sqrt(k1 * MUTVAR));
"""),
    :swap => EidosCode(:swap, """
    // swap: the offspring replace their parents   [Offspring → Zygotes]
    for (k in seqLen(N)) {
        child = p1.addEmpty();
        child.setValue("a0", c0[k]);
        child.setValue("a1", c1[k]);
    }
    adults.fitnessScaling = 0.0;
"""))

eidos_literal(x::Bool) = x ? "T" : "F"
eidos_literal(x::Integer) = string(x)
eidos_literal(x::AbstractFloat) = repr(Float64(x))
eidos_literal(x::Symbol) = string('"', x, '"')

function slim_script(f, p::Params; images = Dict{Symbol,Any}())
    zygotes = LifeCycle[:Zygotes]
    dom(f) == zygotes && codom(f) == zygotes ||
        throw(ArgumentError("a SLiM model iterates a life cycle Zygotes → Zygotes; got $(dom(f)) → $(codom(f))"))
    p.NoQTL == 1 && p.Ploidy == 2 ||
        throw(ArgumentError("the lowering to Eidos covers one locus in diploids (NoQTL = 1, Ploidy = 2); " *
                            "got NoQTL = $(p.NoQTL), Ploidy = $(p.Ploidy)"))
    body = lower(f, :eidos; images)
    cycle = sprint(show, MIME"text/plain"(), f)
    defaults = ("N" => p.PopSize, "VP" => p.PhenotypeVar, "h2" => p.Heritability, "m" => p.MutationScalar,
                "rate" => p.MutationRatePerGamete, "NGEN" => p.NoGenerations, "sampleEvery" => p.Sample,
                "moveGen" => p.MoveFitnessSDGen, "moveSD" => p.MoveFitnessSD, "moveType" => p.MoveType,
                "changeGen" => p.ChangeFitnessGen, "changeSD" => p.ChangeFitnessSD,
                "initOmega" => p.InitFitnessSD, "drift" => p.Drift, "standardise" => p.StandardisePhenotype)
    w = maximum(length ∘ first, defaults)
    io = IOBuffer()
    print(io, """
    // QTSim as a nonWF SLiM model, generated by QTLSim.jl from the life cycle
    //
    //     $cycle
    //
    // by the lowering to Eidos: each step of the life cycle is a fragment of the
    // reproduction() callback below, in the order of the composite. Regenerate it with
    // QTLSim's `slim_script` rather than editing it.
    //
    // Parameters can be overridden with  slim -d "NAME=value"  and the seed set with -s.

    initialize() {
        initializeSLiMModelType("nonWF");

    """)
    for (k, v) in defaults
        println(io, "    if (!exists(\"", k, "\"))", " "^(w - length(k)), " defineConstant(\"", k, "\", ",
                eidos_literal(v), ");")
    end
    print(io, """

        defineConstant("VG", h2 * VP);
        defineConstant("VE", VP - VG);
        defineConstant("ESD", sqrt(VE));
        defineConstant("MUTVAR", (m / (2 * rate)) * VE);           // m/(Ploidy·NoQTL·rate)·V_E
        defineConstant("OMEGA0", initOmega < 0 ? VG / sqrt(m * VE) else initOmega);

        // vestigial genetics: SLiM supplies only the life cycle
        initializeMutationRate(0.0);
        initializeMutationType("m1", 0.5, "f", 0.0);
        initializeGenomicElementType("g1", m1, 1.0);
        initializeGenomicElement(g1, 0, 0);
        initializeRecombinationRate(0.0);
    }

    1 first() {
        sim.addSubpop("p1", N);
        p1.individuals.setValue("a0", 0.0);      // founders: allele values 0
        p1.individuals.setValue("a1", 0.0);
        defineGlobal("THETA", 0.0);
        defineGlobal("OMEGA", OMEGA0);
    }

    // one generation per tick, done in the first call of the callback
    reproduction() {
        gen = sim.cycle - 1;                       // QTSim's generation index
    """)
    print(io, body.code)
    print(io, """
        if (gen >= NGEN - 1) sim.simulationFinished();
        self.active = 0;
    }

    1000000 late() { sim.simulationFinished(); }   // backstop
    """)
    String(take!(io))
end
