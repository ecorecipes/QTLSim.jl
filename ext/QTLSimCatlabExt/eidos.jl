# The lowering of the life cycle to SLiM's Eidos language, and a complete SLiM model built
# around the image of a life cycle.
#
# Each step becomes a fragment of Eidos for one deme, a template in which `@POP@` stands
# for the deme's subpopulation. The fragments implement QTSim's diploid model with
# continuum-of-alleles genetics: each individual carries its two chromosomes as vectors
# of allele values (`a0`, `a1`, one value per locus) and SLiM supplies only the nonWF life
# cycle. With one locus the fragments are, statement for statement, the hand-written port
# slim_port/qtsim.slim; with more they recombine as QTSim does, switching strands at most
# once per gamete. Composition is concatenation within a deme. A morphism on several
# demes runs each deme's code in its own block, which a step that ends the deme
# (`select`, when fewer than two survive) leaves with `break`, skipping the rest of that
# deme's composite as `⋅` does in the kernel; migration acts on the offspring of all the
# demes, in one reproduction() callback per tick.

"""
    EidosCode(name, code)

A fragment of Eidos, the image of a life-cycle morphism within one deme under the
lowering to SLiM; `@POP@` in `code` stands for the deme's subpopulation. The image of a
composite is the concatenation of the images of its steps.
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

"The code of a deme-level fragment for the deme `pop`, in a block that it leaves with `break`."
function deme_block(code::AbstractString, pop::AbstractString)
    body = replace(code, "@POP@" => pop)
    "    do {   // deme $pop\n        if ($pop.getValue(\"extinct\")) break;\n" *
        replace(body, r"^"m => "    ") * "    } while (F);\n"
end

"""
    EidosMeta(n, perm, emit, name)

The image of a life-cycle morphism between tensor products of `n` stages, one per deme,
under the lowering to SLiM: `emit(pops)` is its Eidos for the subpopulations `pops`, and
output position `i` holds the deme at input position `perm[i]`.
"""
struct EidosMeta
    n::Int
    perm::Vector{Int}
    emit::Function
    name::Symbol
end

Base.show(io::IO, m::EidosMeta) = print(io, "EidosMeta(:", m.name, ", ", m.n, " demes)")
Base.show(io::IO, ::MIME"text/plain", m::EidosMeta) = print(io, m.emit(["p$i" for i in 1:m.n]))

eidosmeta(c::EidosCode) = EidosMeta(1, [1], pops -> deme_block(c.code, pops[1]), c.name)
eidosmeta(m::EidosMeta) = m

const AnyEidos = Union{EidosCode,EidosMeta}

function compose(a::EidosMeta, b::EidosMeta)
    a.n == b.n || throw(ArgumentError("composing Eidos on $(a.n) and $(b.n) demes"))
    EidosMeta(a.n, a.perm[b.perm], pops -> a.emit(pops) * b.emit(pops[a.perm]), Symbol(a.name, :⋅, b.name))
end
compose(a::AnyEidos, b::AnyEidos) = compose(eidosmeta(a), eidosmeta(b))
compose(a::AnyEidos, b::AnyEidos, c::AnyEidos, rest::AnyEidos...) = foldl(compose, (a, b, c, rest...))
⋅(a::EidosMeta, b::AnyEidos) = compose(a, b)
⋅(a::EidosCode, b::EidosMeta) = compose(a, b)

otimes(a::EidosMeta, b::EidosMeta) =
    EidosMeta(a.n + b.n, vcat(a.perm, a.n .+ b.perm),
              pops -> a.emit(pops[1:a.n]) * b.emit(pops[a.n+1:end]), Symbol(a.name, :⊗, b.name))
otimes(a::AnyEidos, b::AnyEidos) = otimes(eidosmeta(a), eidosmeta(b))
otimes(a::AnyEidos, b::AnyEidos, c::AnyEidos, rest::AnyEidos...) = foldl(otimes, (a, b, c, rest...))

"What Eidos code at each stage of the life cycle can rely on, for the deme `@POP@`."
const EIDOS_SCOPE = Dict(
    :Zygotes => "the individuals of @POP@ carry allele vectors a0 and a1; @POP@ holds THETA, OMEGA and extinct",
    :Adults => "adults, a0, a1 (concatenated, L values per individual), genotype, env, pheno",
    :Survivors => "W, survived, sidx and ns ≥ 2",
    :Offspring => "c0 and c1, the offspring's chromosomes (concatenated, L values per offspring)")

"The Eidos image of each step of the life cycle, by name."
eidos_images() = Dict{Symbol,Any}(
    :develop => EidosCode(:develop, """
    // develop: genotype the sum of allele values, phenotype genotype + e with e ~ N(0, V_E)   [Zygotes → Adults]
    adults = @POP@.individuals;
    a0 = adults.getValue("a0");
    a1 = adults.getValue("a1");
    genotype = (L == 1) ? a0 + a1 else colSums(matrix(a0, nrow = L)) + colSums(matrix(a1, nrow = L));
    env = rnorm(size(adults), 0.0, ESD);
    pheno = genotype + env;
"""),
    :select => EidosCode(:select, """
    // select: move the optimum, Gaussian fitness with QTSim's 0.4, viability   [Adults → Survivors]
    THETA = @POP@.getValue("THETA");
    OMEGA = @POP@.getValue("OMEGA");
    if (moveGen >= 0 & gen >= moveGen) {
        if (moveType == "phenotype") THETA = mean(pheno) + sqrt(VP) * moveSD;
        else if (moveType == "fitness") THETA = mean(pheno) + OMEGA * moveSD;
        else THETA = THETA + moveSD;
        @POP@.setValue("THETA", THETA);
    }
    x = standardise ? (pheno - THETA) / sd(pheno) else pheno - THETA;
    W = drift ? rep(1.0, size(adults)) else dnorm(x, 0.0, OMEGA) / (0.4 / OMEGA);
    survived = runif(size(adults)) < W;
    sidx = which(survived);
    ns = size(sidx);
    if (ns < 2) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
"""),
    :reproduce => EidosCode(:reproduce, """
    // reproduce: N offspring of two distinct survivors; each gamete copies a parental
    // chromosome, switching strands at most once, as QTSim does   [Survivors → Offspring]
    i1 = rdunif(N, 0, ns - 1);
    i2 = rdunif(N, 0, ns - 1);
    dup = which(i1 == i2);
    while (size(dup) > 0) { i2[dup] = rdunif(size(dup), 0, ns - 1); dup = which(i1 == i2); }
    if (L == 1) {
        sa0 = a0[sidx]; sa1 = a1[sidx];
        c0 = ifelse(runif(N) < 0.5, sa1[i1], sa0[i1]);
        c1 = ifelse(runif(N) < 0.5, sa1[i2], sa0[i2]);
    } else {
        par = c(sidx[i1], sidx[i2]);                        // the parent of each of the 2N gametes
        npts = rpois(2 * N, RECOMB);                        // crossover points drawn for each gamete
        first = runif(2 * N) < 0.5;                         // T: the gamete starts on the parent's a0
        cut = rep(L - 1, 2 * N);                            // the last locus taken from the starting strand
        for (j in which(npts > 0 & npts < L)) cut[j] = min(sample(0:(L - 2), npts[j]));
        loc = rep(0:(L - 1), 2 * N);
        at = repEach(par * L, L) + loc;
        g = ifelse((loc <= repEach(cut, L)) == repEach(first, L), a0[at], a1[at]);
        c0 = g[0:(N * L - 1)];
        c1 = g[(N * L):(2 * N * L - 1)];
    }
"""),
    :record => EidosCode(:record, """
    // record: the parent generation's statistics, every `sampleEvery` generations   [Offspring → Offspring]
    if (gen % sampleEvery == 0)
        catn(gen + "\\t" + mean(pheno) + "\\t" + var(genotype) + "\\t" + ns + "\\t" + @POP@.getValue("THETA") + "\\t" + @POP@.getValue("OMEGA") + ((DEMES > 1) ? "\\t" + @POP@.id else ""));
"""),
    :rescale => EidosCode(:rescale, """
    // rescale: the scheduled change of the fitness width, from the next generation   [Offspring → Offspring]
    if (gen == changeGen) @POP@.setValue("OMEGA", changeSD);
"""),
    :mutate => EidosCode(:mutate, """
    // mutate: Poisson(L·rate) mutations per offspring chromosome, each adding N(0, MUTVAR)
    // to a random locus   [Offspring → Offspring]
    k0 = rpois(N, L * rate);
    k1 = rpois(N, L * rate);
    if (L == 1) {
        c0 = c0 + rnorm(N, 0.0, sqrt(k0 * MUTVAR));
        c1 = c1 + rnorm(N, 0.0, sqrt(k1 * MUTVAR));
    } else {
        at = repEach(0:(N - 1), k0) * L + rdunif(sum(k0), 0, L - 1);
        eff = rnorm(sum(k0), 0.0, sqrt(MUTVAR));
        for (j in seqAlong(at)) c0[at[j]] = c0[at[j]] + eff[j];
        at = repEach(0:(N - 1), k1) * L + rdunif(sum(k1), 0, L - 1);
        eff = rnorm(sum(k1), 0.0, sqrt(MUTVAR));
        for (j in seqAlong(at)) c1[at[j]] = c1[at[j]] + eff[j];
    }
"""),
    :swap => EidosCode(:swap, """
    // swap: the offspring replace their parents   [Offspring → Zygotes]
    kids_@POP@ = @POP@.addEmpty(count = N);
    for (k in seqLen(N)) {
        kids_@POP@[k].setValue("a0", c0[(k * L):(k * L + L - 1)]);
        kids_@POP@[k].setValue("a1", c1[(k * L):(k * L + L - 1)]);
    }
    adults.fitnessScaling = 0.0;
"""))

"The Eidos image of migration among the `k` demes `pops`: offspring swap places between demes."
function migration_code(pops)
    k = length(pops)
    """
    // migrate: each offspring swaps places, with probability migRate, with a random
    // offspring of a random other deme   [Zygotes^$k → Zygotes^$k]
    kids = c($(join(("kids_" * p for p in pops), ", ")));
    sizes = c($(join(("size(kids_" * p * ")" for p in pops), ", ")));
    starts = c(0, cumSum(sizes))[0:$(k - 1)];
    for (d in 0:$(k - 1)) {
        for (i in which(runif(sizes[d]) < migRate)) {
            e = rdunif(1, 0, $(k - 2)); if (e >= d) e = e + 1;
            if (sizes[e] == 0) next;
            a = kids[starts[d] + i]; b = kids[starts[e] + rdunif(1, 0, sizes[e] - 1)];
            t0 = a.getValue("a0"); t1 = a.getValue("a1");
            a.setValue("a0", b.getValue("a0")); a.setValue("a1", b.getValue("a1"));
            b.setValue("a0", t0); b.setValue("a1", t1);
        }
    }
"""
end

eidos_migrate(k) = EidosMeta(k, collect(1:k), migration_code, :migrate)

eidos_literal(x::Bool) = x ? "T" : "F"
eidos_literal(x::Integer) = string(x)
eidos_literal(x::AbstractFloat) = repr(Float64(x))
eidos_literal(x::Symbol) = string('"', x, '"')

function slim_script(f, p::Params; images = Dict{Symbol,Any}(), migration_rate::Real = 0.0)
    k = nfactors(dom(f))
    zk = zygotes_word(max(k, 1))
    dom(f) == zk && codom(f) == zk ||
        throw(ArgumentError("a SLiM model iterates a life cycle on zygotes, Zygotes → Zygotes or " *
                            "Zygotes⊗…⊗Zygotes → Zygotes⊗…⊗Zygotes; got $(dom(f)) → $(codom(f))"))
    p.Ploidy == 2 || throw(ArgumentError("the lowering to Eidos covers diploids (Ploidy = 2); got Ploidy = $(p.Ploidy)"))
    p.NoQTL == 1 || p.CrossoverModel === :qtsim ||
        throw(ArgumentError("with more than one locus the lowering to Eidos recombines as QTSim does " *
                            "(CrossoverModel = :qtsim); got :$(p.CrossoverModel)"))
    img = lower(f, :eidos; images)
    pops = ["p$i" for i in 1:k]
    if img isa EidosCode
        body = deme_block(img.code, "p1")
    else
        img.perm == 1:k ||
            throw(ArgumentError("a SLiM model needs a life cycle that returns every deme to its place; " *
                                "its braidings do not cancel"))
        body = img.emit(pops)
    end
    cycle = sprint(show, MIME"text/plain"(), f)
    defaults = ("N" => p.PopSize, "VP" => p.PhenotypeVar, "h2" => p.Heritability, "m" => p.MutationScalar,
                "rate" => p.MutationRatePerGamete, "L" => p.NoQTL, "RECOMB" => p.RecombPoissonMean,
                "NGEN" => p.NoGenerations, "sampleEvery" => p.Sample,
                "moveGen" => p.MoveFitnessSDGen, "moveSD" => p.MoveFitnessSD, "moveType" => p.MoveType,
                "changeGen" => p.ChangeFitnessGen, "changeSD" => p.ChangeFitnessSD,
                "initOmega" => p.InitFitnessSD, "drift" => p.Drift, "standardise" => p.StandardisePhenotype,
                "migRate" => Float64(migration_rate))
    w = maximum(length ∘ first, defaults)
    io = IOBuffer()
    print(io, """
    // QTSim as a nonWF SLiM model, generated by QTLSim.jl from the life cycle
    //
    //     $cycle
    //
    // by the lowering to Eidos: each step of the life cycle is a fragment of the
    // reproduction() callback below, in the order of the composite, for each of its
    // $k deme$(k == 1 ? "" : "s"). Regenerate it with QTLSim's `slim_script` rather than editing it.
    //
    // Parameters can be overridden with  slim -d "NAME=value"  and the seed set with -s.

    initialize() {
        initializeSLiMModelType("nonWF");

    """)
    for (key, v) in defaults
        println(io, "    if (!exists(\"", key, "\"))", " "^(w - length(key)), " defineConstant(\"", key, "\", ",
                eidos_literal(v), ");")
    end
    print(io, """
        defineConstant("DEMES", $k);                             // fixed by the life cycle

        defineConstant("VG", h2 * VP);
        defineConstant("VE", VP - VG);
        defineConstant("ESD", sqrt(VE));
        defineConstant("MUTVAR", (m / (2 * L * rate)) * VE);       // m/(Ploidy·NoQTL·rate)·V_E
        defineConstant("OMEGA0", initOmega < 0 ? VG / sqrt(m * VE) else initOmega);

        // vestigial genetics: SLiM supplies only the life cycle
        initializeMutationRate(0.0);
        initializeMutationType("m1", 0.5, "f", 0.0);
        initializeGenomicElementType("g1", m1, 1.0);
        initializeGenomicElement(g1, 0, 0);
        initializeRecombinationRate(0.0);
    }

    1 first() {
        for (i in seqLen(DEMES)) {
            pop = sim.addSubpop(i + 1, N);
            pop.individuals.setValue("a0", rep(0.0, L));   // founders: allele values 0
            pop.individuals.setValue("a1", rep(0.0, L));
            pop.setValue("THETA", 0.0);
            pop.setValue("OMEGA", OMEGA0);
            pop.setValue("extinct", F);
        }
    }

    // one generation per tick, done in the first call of the callback
    reproduction() {
        gen = sim.cycle - 1;                       // QTSim's generation index
        if (gen == 0) catn("gen\\tphenoMean\\tgenoVar\\tnSurv\\ttheta\\tomega" + ((DEMES > 1) ? "\\tdeme" else ""));
    """)
    for pop in pops
        println(io, "    kids_", pop, " = ", pop, ".individuals[integer(0)];")
    end
    print(io, body)
    print(io, """
        dead = 0;
        for (pop in sim.subpopulations) dead = dead + asInteger(pop.getValue("extinct"));
        if (dead == DEMES) { catn("Pop all died."); sim.simulationFinished(); }
        else if (gen >= NGEN - 1) sim.simulationFinished();
        self.active = 0;
    }

    1000000 late() { sim.simulationFinished(); }   // backstop
    """)
    String(take!(io))
end
