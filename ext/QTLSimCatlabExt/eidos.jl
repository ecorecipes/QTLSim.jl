# The lowering of the life cycle to SLiM's Eidos language, and a complete SLiM model built
# around the image of a life cycle.
#
# Each step becomes a fragment of Eidos for one deme, a template in which `@POP@` stands
# for the deme's subpopulation. The fragments implement QTSim's diploid model with
# continuum-of-alleles genetics: each individual carries its two chromosomes as vectors
# of allele values (`a0`, `a1`, one value per locus) and SLiM supplies only the nonWF life
# cycle. With one locus the fragments are, statement for statement, the hand-written port
# validate/slim/qtsim.slim; with more they recombine as QTSim does, switching strands at most
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
    :Survivors => "W, survived, sidx and ns ≥ 2 (and, when a model has fecundity selection, weighted, sw and nw)",
    :Offspring => "c0 and c1, the offspring's chromosomes (concatenated, L values per offspring)")

"The Eidos image of each step of the life cycle in `model` (by default QTSim's one trait), by name."
eidos_images() = eidos_images(OneTraitModel())
eidos_images(::OneTraitModel) = Dict{Symbol,Any}(
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
    // reproduce: NOFF offspring of two distinct survivors; each gamete copies a parental
    // chromosome, switching strands at most once, as QTSim does   [Survivors → Offspring]
    // NOFF is N, as in QTSim, or density-regulated (QTLSim `noffspring`): Beverton-Holt
    NOFF = N;
    if (demography == "beverton_holt")
        NOFF = min(N, rpois(1, growthRate * ns / (1 + (growthRate - 1) * ns / N)));
    if (NOFF == 0) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
    i1 = rdunif(NOFF, 0, ns - 1);
    i2 = rdunif(NOFF, 0, ns - 1);
    dup = which(i1 == i2);
    while (size(dup) > 0) { i2[dup] = rdunif(size(dup), 0, ns - 1); dup = which(i1 == i2); }
    if (L == 1) {
        sa0 = a0[sidx]; sa1 = a1[sidx];
        c0 = ifelse(runif(NOFF) < 0.5, sa1[i1], sa0[i1]);
        c1 = ifelse(runif(NOFF) < 0.5, sa1[i2], sa0[i2]);
    } else {
        par = c(sidx[i1], sidx[i2]);                        // the parent of each of the 2N gametes
        npts = rpois(2 * NOFF, RECOMB);                        // crossover points drawn for each gamete
        first = runif(2 * NOFF) < 0.5;                         // T: the gamete starts on the parent's a0
        cut = rep(L - 1, 2 * NOFF);                            // the last locus taken from the starting strand
        for (j in which(npts > 0 & npts < L)) cut[j] = min(sample(0:(L - 2), npts[j]));
        loc = rep(0:(L - 1), 2 * NOFF);
        at = repEach(par * L, L) + loc;
        g = ifelse((loc <= repEach(cut, L)) == repEach(first, L), a0[at], a1[at]);
        c0 = g[0:(NOFF * L - 1)];
        c1 = g[(NOFF * L):(2 * NOFF * L - 1)];
    }
"""),
    :record => EidosCode(:record, """
    // record: the parent generation's statistics, every `sampleEvery` generations   [Offspring → Offspring]
    if (gen % sampleEvery == 0)
        catn(gen + "\\t" + mean(pheno) + "\\t" + var(genotype) + "\\t" + ns + "\\t" + @POP@.getValue("THETA") + "\\t" + @POP@.getValue("OMEGA") + ((DEMES > 1) ? "\\t" + @POP@.id else "") + ((demography != "fixed") ? "\\t" + size(adults) else ""));
"""),
    :rescale => EidosCode(:rescale, """
    // rescale: the scheduled change of the fitness width, from the next generation   [Offspring → Offspring]
    if (gen == changeGen) @POP@.setValue("OMEGA", changeSD);
"""),
    :mutate => EidosCode(:mutate, """
    // mutate: Poisson(L·rate) mutations per offspring chromosome, each adding N(0, MUTVAR)
    // to a random locus   [Offspring → Offspring]
    k0 = rpois(NOFF, L * rate);
    k1 = rpois(NOFF, L * rate);
    if (L == 1) {
        c0 = c0 + rnorm(NOFF, 0.0, sqrt(k0 * MUTVAR));
        c1 = c1 + rnorm(NOFF, 0.0, sqrt(k1 * MUTVAR));
    } else {
        at = repEach(0:(NOFF - 1), k0) * L + rdunif(sum(k0), 0, L - 1);
        eff = rnorm(sum(k0), 0.0, sqrt(MUTVAR));
        for (j in seqAlong(at)) c0[at[j]] = c0[at[j]] + eff[j];
        at = repEach(0:(NOFF - 1), k1) * L + rdunif(sum(k1), 0, L - 1);
        eff = rnorm(sum(k1), 0.0, sqrt(MUTVAR));
        for (j in seqAlong(at)) c1[at[j]] = c1[at[j]] + eff[j];
    }
"""),
    :swap => EidosCode(:swap, """
    // swap: the offspring replace their parents   [Offspring → Zygotes]
    kids_@POP@ = @POP@.addEmpty(count = NOFF);
    for (k in seqLen(NOFF)) {
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
    // offspring of a random other deme, the demes in a random order   [Zygotes^$k → Zygotes^$k]
    kids = c($(join(("kids_" * p for p in pops), ", ")));
    sizes = c($(join(("size(kids_" * p * ")" for p in pops), ", ")));
    starts = c(0, cumSum(sizes))[0:$(k - 1)];
    for (d in sample(0:$(k - 1), $k)) {
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

# The fluctuating moves (QTLSim `moveoptimum!`, MoveType :red_noise and :periodic), which
# QTSim lacks. `slim_script` splices them into the image of `select` only when `p` uses one,
# so that the models of QTSim's own moves are unchanged. The optimum before the move
# (THETA0) and the AR(1) noise (EPS) are deme state on @POP@, as THETA is.
const EIDOS_CONTINUOUS_MOVE = "        else THETA = THETA + moveSD;\n"
const EIDOS_FLUCTUATING_MOVE = """
        else if (moveType == "red_noise" | moveType == "periodic") {
            // θ = θ0 + d t + ε, ε AR(1) with stationary SD moveNoiseSD; or + A sin(2π t/P)
            t = gen - moveGen;
            if (t == 0) @POP@.setValue("THETA0", THETA);
            THETA = @POP@.getValue("THETA0") + moveSD * t;
            if (moveType == "red_noise") {
                EPS = (t == 0) ? rnorm(1, 0.0, moveNoiseSD) else
                    moveAutocorr * @POP@.getValue("EPS") + moveNoiseSD * sqrt(1 - moveAutocorr^2) * rnorm(1);
                @POP@.setValue("EPS", EPS);
                THETA = THETA + EPS;
            } else THETA = THETA + moveAmplitude * sin(2 * PI * t / movePeriod);
        }
"""

"The Eidos image of `select` (by default QTSim's) with the fluctuating moves of `moveoptimum!` added."
function eidos_fluctuating_select(code::AbstractString = eidos_images()[:select].code)
    occursin(EIDOS_FLUCTUATING_MOVE, code) && return EidosCode(:select, code)
    occursin(EIDOS_CONTINUOUS_MOVE, code) ||
        error("the Eidos image of select has changed; update EIDOS_CONTINUOUS_MOVE")
    EidosCode(:select, replace(code, EIDOS_CONTINUOUS_MOVE => EIDOS_FLUCTUATING_MOVE * EIDOS_CONTINUOUS_MOVE))
end

# Fecundity selection (QTLSim `Selection = :fecundity`), which QTSim lacks. `slim_script`
# splices it into the images of `select` and `reproduce` only when a deme or a step uses
# it, so that the models of viability selection are unchanged. In such a model `select`
# leaves at the stage Survivors whether parents are weighted (`weighted`), their mating
# weights `sw` (one per survivor) and their sum `nw`, the count that density regulation
# uses; under viability `weighted` is F, `sw` all 1 and `nw` = ns, and `reproduce` draws
# the parents as before. The selection mode is the constant `selection`, so that a step
# (`withparams(L.select; Selection = :fecundity)`) or a deme can set its own.
const EIDOS_VIABILITY = """
    survived = runif(size(adults)) < W;
    sidx = which(survived);
    ns = size(sidx);
    if (ns < 2) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
"""
const EIDOS_SELECTION = """
    weighted = (selection == "fecundity");
    if (weighted) {
        // fecundity: everyone survives, and parents are drawn in proportion to W
        survived = rep(T, size(adults));
        sidx = seqLen(size(adults));
        ns = size(sidx);
        sw = W;
        if (sum(W > 0) < 2) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
    } else {
        survived = runif(size(adults)) < W;
        sidx = which(survived);
        ns = size(sidx);
        sw = rep(1.0, ns);
        if (ns < 2) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
    }
    nw = sum(sw);
"""
const EIDOS_UNIFORM_COUNT = "growthRate * ns / (1 + (growthRate - 1) * ns / N)"
const EIDOS_WEIGHTED_COUNT = "growthRate * nw / (1 + (growthRate - 1) * nw / N)"
const EIDOS_UNIFORM_PARENTS = """
    i1 = rdunif(NOFF, 0, ns - 1);
    i2 = rdunif(NOFF, 0, ns - 1);
    dup = which(i1 == i2);
    while (size(dup) > 0) { i2[dup] = rdunif(size(dup), 0, ns - 1); dup = which(i1 == i2); }
"""
const EIDOS_WEIGHTED_PARENTS = """
    if (weighted) {
        // two distinct parents, each drawn with probability proportional to its weight
        i1 = sample(0:(ns - 1), NOFF, T, sw);
        i2 = sample(0:(ns - 1), NOFF, T, sw);
        dup = which(i1 == i2);
        for (t in 1:100) { if (size(dup) == 0) break; i2[dup] = sample(0:(ns - 1), size(dup), T, sw); dup = which(i1 == i2); }
        for (j in dup) { w2 = sw; w2[i1[j]] = 0.0; i2[j] = sample(0:(ns - 1), 1, T, w2); }   // nearly all weight on one
    } else {
        i1 = rdunif(NOFF, 0, ns - 1);
        i2 = rdunif(NOFF, 0, ns - 1);
        dup = which(i1 == i2);
        while (size(dup) > 0) { i2[dup] = rdunif(size(dup), 0, ns - 1); dup = which(i1 == i2); }
    }
"""

"Splice `new` for `old` into an Eidos image, which must contain `old` once."
function eidos_splice(code::AbstractString, old::AbstractString, new::AbstractString, what)
    occursin(new, code) && return String(code)
    count(old, code) == 1 || error("the Eidos image of $what has changed; update the fecundity splice")
    replace(code, old => new)
end

"The Eidos images of `select` and `reproduce` (by default QTSim's) with fecundity selection added."
eidos_fecundity_select(code::AbstractString = eidos_images()[:select].code) =
    EidosCode(:select, eidos_splice(code, EIDOS_VIABILITY, EIDOS_SELECTION, "select"))
function eidos_fecundity_reproduce(code::AbstractString = eidos_images()[:reproduce].code)
    code = eidos_splice(code, EIDOS_UNIFORM_COUNT, EIDOS_WEIGHTED_COUNT, "reproduce")
    EidosCode(:reproduce, eidos_splice(code, EIDOS_UNIFORM_PARENTS, EIDOS_WEIGHTED_PARENTS, "reproduce"))
end

# --- steps with parameters, and demes with parameters of their own ---------

"The Eidos constant of `slim_script` that holds each field of `Params` a step can override."
const EIDOS_CONSTANT = Dict{Symbol,String}(
    :MoveType => "moveType", :MoveFitnessSD => "moveSD", :MoveFitnessSDGen => "moveGen",
    :MoveNoiseSD => "moveNoiseSD", :MoveAutocorr => "moveAutocorr", :MoveAmplitude => "moveAmplitude",
    :MovePeriod => "movePeriod", :StandardisePhenotype => "standardise", :Drift => "drift",
    :RecombPoissonMean => "RECOMB", :Demography => "demography", :GrowthRate => "growthRate",
    :Sample => "sampleEvery", :ChangeFitnessGen => "changeGen", :ChangeFitnessSD => "changeSD",
    :MutationScalar => "m", :MutationRatePerGamete => "rate", :Selection => "selection",
    :PopSize => "N", :PhenotypeVar => "VP", :Heritability => "h2", :InitFitnessSD => "initOmega",
    :NoQTL => "L", :NoGenerations => "NGEN")

"""
Replace whole-word identifiers in the code of an Eidos fragment, leaving its comments
alone: `subs` maps an identifier to its replacement. One pass, so a replacement is not
itself rewritten.
"""
function substitute_identifiers(code::AbstractString, subs)
    isempty(subs) && return String(code)
    d = Dict{String,String}(String(k) => String(v) for (k, v) in subs)
    re = Regex("\\b(" * join(sort!(collect(keys(d)); by = length, rev = true), "|") * ")\\b")
    lines = map(split(code, '\n')) do line
        i = findfirst("//", line)
        body, comment = i === nothing ? (line, "") : (line[1:prevind(line, first(i))], line[first(i):end])
        replace(body, re => s -> d[s]) * comment
    end
    join(lines, '\n')
end

"""
    eidos_with_params(model, s::Step, c::EidosCode)

The Eidos image of a step with parameters, from the image `c` of the plain step: each
constant of `slim_script` the step overrides (`EIDOS_CONSTANT`) becomes the step's value,
a literal, and a step's `ω` and `θ` replace the deme's `OMEGA` and `THETA` where the
step reads them. A `select` whose own `MoveType` is a fluctuating move gets the code of
the fluctuating moves. As in the kernel, the deme's state is left alone.

For `TraitModel(n)` (the Eidos of eidos_traits.jl) a step's `θ` is a vector literal in
place of the deme's `THETA`, `Ω` replaces the width (`OMEGA`, 1.0) and the inverse fitness
shape (`SINV`), `SelCorr` the shape, `MoveDirection` the direction `UDIR`, and the mutation
fields the covariance of one mutation `MUTCOV`.
"""
function eidos_with_params(::OneTraitModel, s, c::EidosCode)
    pr = s.params
    get(pr, :CrossoverModel, :qtsim) === :qtsim ||
        throw(ArgumentError("the lowering to Eidos recombines as QTSim does (CrossoverModel = :qtsim); the step $s does not"))
    code = c.code
    s.name === :select && get(pr, :MoveType, :none) in (:red_noise, :periodic) &&
        (code = eidos_fluctuating_select(code).code)
    haskey(pr, :ω) && (code = replace(code, "@POP@.getValue(\"OMEGA\")" => eidos_literal(pr.ω)))
    haskey(pr, :θ) && (code = replace(code, "@POP@.getValue(\"THETA\")" => eidos_literal(pr.θ)))
    # the mutational variance per mutation, from the step's mutation scalar or rate
    (haskey(pr, :MutationScalar) || haskey(pr, :MutationRatePerGamete)) &&
        (code = substitute_identifiers(code, ["MUTVAR" => "((m / (2 * L * rate)) * VE)"]))
    subs = Dict{String,String}(EIDOS_CONSTANT[k] => eidos_literal(v) for (k, v) in pairs(pr) if haskey(EIDOS_CONSTANT, k))
    haskey(pr, :θ) && (subs["moveGen"] = "-1")
    EidosCode(Symbol(string(s)), substitute_identifiers(code, subs))
end

"""
The Eidos constants of `slim_script` that become deme state when the demes have
parameters of their own: the code of a deme reads `@POP@.getValue("NAME")` in place of
the constant `NAME`.
"""
deme_state(code::AbstractString, names) =
    substitute_identifiers(code, [n => "@POP@.getValue(\"$n\")" for n in names])

"Migration among `k` demes at the step's own `rate`, a number or a matrix, written into its Eidos."
function eidos_migrate_at(k, rate::Real)
    lit = eidos_literal(Float64(rate))
    EidosMeta(k, collect(1:k), pops -> replace(migration_code(pops), r"\bmigRate\b" => lit), :migrate)
end
function eidos_migrate_at(k, M::AbstractMatrix)
    QTLSim.check_migration(M, k)
    lit = eidos_literal(vec(permutedims(Float64.(M))))
    EidosMeta(k, collect(1:k), pops -> "    MIGstep = $lit;   // this step's migration matrix, row by row\n" *
                                       replace(migration_graph_code(pops), r"\bMIG\b" => "MIGstep"), :migrate)
end

eidos_migrate(k) = EidosMeta(k, collect(1:k), migration_code, :migrate)

function migration_graph_code(pops)
    k = length(pops)
    """
    // migrate: each offspring of deme d swaps places, with probability MIG[d, e], with a
    // random offspring of deme e (MIG row by row), the demes in a random order
    // [Zygotes^$k → Zygotes^$k]
    kids = c($(join(("kids_" * p for p in pops), ", ")));
    sizes = c($(join(("size(kids_" * p * ")" for p in pops), ", ")));
    starts = c(0, cumSum(sizes))[0:$(k - 1)];
    for (d in sample(0:$(k - 1), $k)) {
        cum = cumSum(MIG[(d * $k):(d * $k + $(k - 1))]);
        u = runif(sizes[d]);
        for (i in which(u < cum[$(k - 1)])) {
            e = min(which(cum > u[i]));
            if (sizes[e] == 0) next;
            a = kids[starts[d] + i]; b = kids[starts[e] + rdunif(1, 0, sizes[e] - 1)];
            t0 = a.getValue("a0"); t1 = a.getValue("a1");
            a.setValue("a0", b.getValue("a0")); a.setValue("a1", b.getValue("a1"));
            b.setValue("a0", t0); b.setValue("a1", t1);
        }
    }
"""
end

eidos_migrate(k, graph::Bool) = graph ? EidosMeta(k, collect(1:k), migration_graph_code, :migrate) : eidos_migrate(k)

eidos_literal(x::Bool) = x ? "T" : "F"
eidos_literal(x::Integer) = string(x)
eidos_literal(x::AbstractFloat) = repr(Float64(x))
eidos_literal(x::Symbol) = string('"', x, '"')
eidos_literal(x::AbstractVector) = "c(" * join(eidos_literal.(x), ", ") * ")"

const FLUCTUATING_MOVES = (:red_noise, :periodic)

# The fields that must be the same in every deme of a SLiM model.
const SHARED_FIELDS = (:Ploidy, :NoQTL, :NoGenerations, :CrossoverModel, :Demography)

slim_script(f, p::Params; images = Dict{Symbol,Any}(), migration_rate::Union{Real,AbstractMatrix} = 0.0) =
    slim_model(f, fill(p, max(nfactors(dom(f)), 1)), images, migration_rate)

function slim_script(f, ps::AbstractVector{Params}; images = Dict{Symbol,Any}(),
                     migration_rate::Union{Real,AbstractMatrix} = 0.0)
    k = max(nfactors(dom(f)), 1)
    length(ps) == k || throw(ArgumentError("$(length(ps)) parameter sets for a life cycle on $k deme$(k == 1 ? "" : "s")"))
    slim_model(f, collect(Params, ps), images, migration_rate)
end

function slim_model(f, ps::Vector{Params}, images, migration_rate)
    k = nfactors(dom(f))
    zk = zygotes_word(max(k, 1))
    dom(f) == zk && codom(f) == zk ||
        throw(ArgumentError("a SLiM model iterates a life cycle on zygotes, Zygotes → Zygotes or " *
                            "Zygotes⊗…⊗Zygotes → Zygotes⊗…⊗Zygotes; got $(dom(f)) → $(codom(f))"))
    # one tick of SLiM is one generation of each deme: a step into Zygotes from another stage
    # (swap) ends a generation, and the life cycle has one per deme
    Z = LifeCycle[:Zygotes]
    ngen = count(g -> !(g isa ObExpr) && codom(g) == Z && dom(g) != Z, stepsof(f))
    ngen == max(k, 1) ||
        throw(ArgumentError("a SLiM model runs one generation of each deme per tick, so the life cycle has one " *
                            "step into Zygotes (swap) per deme; it has $ngen for $(max(k, 1)) deme$(k > 1 ? "s" : "")"))
    p = ps[1]
    p.Ploidy == 2 || throw(ArgumentError("the lowering to Eidos covers diploids (Ploidy = 2); got Ploidy = $(p.Ploidy)"))
    p.NoQTL == 1 || p.CrossoverModel === :qtsim ||
        throw(ArgumentError("with more than one locus the lowering to Eidos recombines as QTSim does " *
                            "(CrossoverModel = :qtsim); got :$(p.CrossoverModel)"))
    for fld in SHARED_FIELDS
        all(q -> getfield(q, fld) == getfield(p, fld), ps) ||
            throw(ArgumentError("in a SLiM model the demes have the same $fld; got $([getfield(q, fld) for q in ps])"))
    end
    migration_rate isa AbstractMatrix && QTLSim.check_migration(migration_rate, k)
    # several traits: a second table of images (eidos_traits.jl), shared trait parameters
    traits = p.NoTraits > 1
    if traits
        for fld in TRAIT_FIELDS
            all(q -> getfield(q, fld) == getfield(p, fld), ps) ||
                throw(ArgumentError("in a SLiM model of several traits the demes have the same $fld"))
        end
        QTLSim.mineigen(QTLSim.mutational_covariance(p)) > 0 ||
            throw(ArgumentError("the lowering to Eidos draws mutations with rmvnorm, which needs a positive definite mutational covariance"))
    end
    model = model_of(p)
    table = eidos_images(model)
    # QTLSim's fluctuating moves: in the image of the plain select when a deme uses one, and
    # their constants when a deme or a step does
    fluctuating = any(q -> q.MoveType in FLUCTUATING_MOVES, ps)
    # fecundity selection: in the images of the plain select and reproduce when a deme or a
    # step uses it, with the constant `selection`
    fecund = any(q -> q.Selection === :fecundity, ps) ||
             any(g -> !(g isa ObExpr) && haskey(stepparams(g), :Selection), stepsof(f))
    builtin = Dict{Any,Any}()
    fluctuating && (builtin[:select] = eidos_fluctuating_select(table[:select].code))
    if fecund
        builtin[:select] = eidos_fecundity_select(get(builtin, :select, table[:select]).code)
        builtin[:reproduce] = eidos_fecundity_reproduce(table[:reproduce].code)
    end
    isempty(builtin) || (images = merge(builtin, Dict{Any,Any}(images)))
    stepmoves = any(g -> !(g isa ObExpr) && get(stepparams(g), :MoveType, :none) in FLUCTUATING_MOVES, stepsof(f))
    # a constant whose value differs between demes is a vector, one value per deme, and deme state
    value(fld) = (vs = [getfield(q, fld) for q in ps]; all(==(vs[1]), vs) ? vs[1] : vs)
    defaults = ("N" => value(:PopSize), "VP" => value(:PhenotypeVar), "h2" => value(:Heritability),
                "m" => value(:MutationScalar), "rate" => value(:MutationRatePerGamete), "L" => p.NoQTL,
                "RECOMB" => value(:RecombPoissonMean), "NGEN" => p.NoGenerations, "sampleEvery" => value(:Sample),
                "moveGen" => value(:MoveFitnessSDGen), "moveSD" => value(:MoveFitnessSD), "moveType" => value(:MoveType),
                "changeGen" => value(:ChangeFitnessGen), "changeSD" => value(:ChangeFitnessSD),
                "initOmega" => value(:InitFitnessSD), "drift" => value(:Drift), "standardise" => value(:StandardisePhenotype),
                "demography" => p.Demography, "growthRate" => value(:GrowthRate),
                (migration_rate isa AbstractMatrix ? "MIG" => vec(permutedims(Float64.(migration_rate))) :
                 "migRate" => Float64(migration_rate)))
    fecund && (defaults = (defaults..., "selection" => value(:Selection)))
    traits && (defaults = (defaults..., eidos_trait_constants(p)...))
    (fluctuating || stepmoves) &&
        (defaults = (defaults..., "moveNoiseSD" => value(:MoveNoiseSD), "moveAutocorr" => value(:MoveAutocorr),
                     "moveAmplitude" => value(:MoveAmplitude), "movePeriod" => value(:MovePeriod)))
    perdeme = [key for (key, v) in defaults if v isa AbstractVector && !(key in VECTOR_CONSTANTS)]
    vary(names...) = any(in(perdeme), names)
    derived = [n for (n, on) in ("VG" => vary("VP", "h2"), "VE" => vary("VP", "h2"), "ESD" => vary("VP", "h2"),
                                 "MUTVAR" => vary("m", "rate", "VP", "h2"),
                                 "OMEGA0" => vary("initOmega", "VP", "h2", "m")) if on]
    # what the steps read from the deme's subpopulation (initOmega and OMEGA0 found the deme)
    state = filter(n -> !(n in ("initOmega", "OMEGA0")), vcat(perdeme, derived))
    img = lower(f, :eidos; images, migration_rate, deme_constants = state, model)
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
    """)
    isempty(perdeme) || print(io, """
    // The demes have parameters of their own: $(join(perdeme, ", ")) hold one value per deme
    // (a single value applies to all), and each deme keeps its values on its subpopulation.
    """)
    print(io, """

    initialize() {
        initializeSLiMModelType("nonWF");

    """)
    for (key, v) in defaults
        println(io, "    if (!exists(\"", key, "\"))", " "^(w - length(key)), " defineConstant(\"", key, "\", ",
                eidos_literal(v), ");")
    end
    omega0 = "OMEGA0" in derived ? "ifelse(perDeme(initOmega) < 0, perDeme(VG / sqrt(m * VE)), perDeme(initOmega))" :
                                 "initOmega < 0 ? VG / sqrt(m * VE) else initOmega"
    print(io, """
        defineConstant("DEMES", $k);                             // fixed by the life cycle

        defineConstant("VG", h2 * VP);
        defineConstant("VE", VP - VG);
        defineConstant("ESD", sqrt(VE));
        defineConstant("MUTVAR", (m / (2 * L * rate)) * VE);       // m/(Ploidy·NoQTL·rate)·V_E
        defineConstant("OMEGA0", $omega0);
    $(traits ? eidos_mutcov(p) : "")
        // vestigial genetics: SLiM supplies only the life cycle
        initializeMutationRate(0.0);
        initializeMutationType("m1", 0.5, "f", 0.0);
        initializeGenomicElementType("g1", m1, 1.0);
        initializeGenomicElement(g1, 0, 0);
        initializeRecombinationRate(0.0);
    }

    """)
    at(n) = n in perdeme || n in derived ? "perDeme($n)[i]" : n
    isempty(perdeme) || print(io, """
    // a constant with one value for all demes, or one per deme, as one per deme
    function (*)perDeme(* x) { return (size(x) == 1) ? rep(x, DEMES) else x; }

    """)
    print(io, """
    1 first() {
        for (i in seqLen(DEMES)) {
            pop = sim.addSubpop(i + 1, $(at("N")));
            pop.individuals.setValue("a0", rep(0.0, $(traits ? "L * NT" : "L")));   // founders: allele values 0
            pop.individuals.setValue("a1", rep(0.0, $(traits ? "L * NT" : "L")));
            pop.setValue("THETA", $(traits ? "rep(0.0, NT)" : "0.0"));
            pop.setValue("OMEGA", $(at("OMEGA0")));
            pop.setValue("extinct", F);
    """)
    for n in state
        println(io, "        pop.setValue(\"", n, "\", perDeme(", n, ")[i]);")
    end
    print(io, """
        }
    }

    // one generation per tick, done in the first call of the callback
    reproduction() {
        gen = sim.cycle - 1;                       // QTSim's generation index
        if (gen == 0) catn("gen\\tphenoMean\\tgenoVar\\tnSurv\\ttheta\\tomega" + ((DEMES > 1) ? "\\tdeme" else "") +
                       ((demography != "fixed") ? "\\tpopsize" else "")$(traits ? trait_header(p.NoTraits) : ""));
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
