# The lowering to Eidos for several traits: a second table of Eidos images for the same
# steps of the same life cycle. An individual carries its two chromosomes as `a0` and
# `a1`, L·NT values each, locus by locus with the NT effects of a locus together (the
# kernel's layout, blocks.jl); the optimum THETA is a vector of NT on the deme's
# subpopulation. `slim_script` uses this table when `NoTraits > 1`; the one-trait table
# (eidos.jl) and its output are untouched. The constants are NT, SINV (the inverse of the
# fitness function's shape S, Ω = OMEGA² S, row by row), MUTCOV (the covariance of one
# mutation's effects, M/(2·L·rate)) and UDIR (the unit direction the optimum moves in).
# The segments of `select` and `reproduce` that the fecundity splice replaces are the
# one-trait table's, verbatim, so `Selection = :fecundity` lowers here too.

"What Eidos code at each stage of the life cycle of several traits can rely on, for the deme `@POP@`."
const EIDOS_TRAIT_SCOPE = Dict(
    :Zygotes => "the individuals of @POP@ carry a0 and a1, L·NT values each; @POP@ holds THETA (NT values), OMEGA and extinct",
    :Adults => "adults, a0, a1, na, G, env and pheno (NT × na matrices: trait by individual)",
    :Survivors => "W, survived, sidx and ns ≥ 2",
    :Offspring => "c0 and c1, the offspring's chromosomes (L·NT values per offspring)")

"The Eidos image of each step of the life cycle of several traits, by name."
eidos_trait_images() = Dict{Symbol,Any}(
    :develop => EidosCode(:develop, """
    // develop: genotype the sum of the effect vectors, phenotype genotype + e with
    // e ~ N(0, V_E I); G, env and pheno are NT × na matrices   [Zygotes → Adults]
    adults = @POP@.individuals;
    na = size(adults);
    a0 = adults.getValue("a0");
    a1 = adults.getValue("a1");
    gm = matrix(a0 + a1, nrow = NT * L);
    gsum = NULL;
    for (r in 0:(NT - 1)) gsum = c(gsum, colSums(gm[seq(r, NT * L - 1, by = NT), ]));
    G = matrix(gsum, nrow = NT, byrow = T);
    env = matrix(rnorm(NT * na, 0.0, ESD), nrow = NT);
    pheno = G + env;
    genotype = drop(G[0, ]);
"""),
    :select => EidosCode(:select, """
    // select: move the optimum along UDIR, multivariate Gaussian fitness with QTSim's 0.4,
    // W = exp(-(z - θ)' S⁻¹ (z - θ) / 2 OMEGA²) / (0.4 √(2π)), viability   [Adults → Survivors]
    THETA = @POP@.getValue("THETA");
    OMEGA = @POP@.getValue("OMEGA");
    zbar = rowSums(pheno) / na;
    if (moveGen >= 0 & gen >= moveGen) {
        if (moveType == "phenotype") THETA = zbar + sqrt(VP) * moveSD * UDIR;
        else if (moveType == "fitness") THETA = zbar + OMEGA * moveSD * UDIR;
        else THETA = THETA + moveSD * UDIR;
        @POP@.setValue("THETA", THETA);
    }
    D = pheno - matrix(rep(THETA, na), nrow = NT);
    Q = colSums(D * matrixMult(matrix(SINV, nrow = NT, byrow = T), D));
    W = drift ? rep(1.0, size(adults)) else exp(-0.5 * Q / OMEGA^2) / (0.4 * sqrt(2 * PI));
    survived = runif(size(adults)) < W;
    sidx = which(survived);
    ns = size(sidx);
    if (ns < 2) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
"""),
    :reproduce => EidosCode(:reproduce, """
    // reproduce: NOFF offspring of two distinct survivors; each gamete copies a parental
    // chromosome, switching strands at most once, as QTSim does; a locus carries its NT
    // effects together   [Survivors → Offspring]
    NOFF = N;
    if (demography == "beverton_holt")
        NOFF = min(N, rpois(1, growthRate * ns / (1 + (growthRate - 1) * ns / N)));
    if (NOFF == 0) { @POP@.setValue("extinct", T); adults.fitnessScaling = 0.0; break; }
    i1 = rdunif(NOFF, 0, ns - 1);
    i2 = rdunif(NOFF, 0, ns - 1);
    dup = which(i1 == i2);
    while (size(dup) > 0) { i2[dup] = rdunif(size(dup), 0, ns - 1); dup = which(i1 == i2); }
    par = c(sidx[i1], sidx[i2]);                        // the parent of each of the 2N gametes
    npts = rpois(2 * NOFF, RECOMB);                        // crossover points drawn for each gamete
    first = runif(2 * NOFF) < 0.5;                         // T: the gamete starts on the parent's a0
    cut = rep(L - 1, 2 * NOFF);                            // the last locus taken from the starting strand
    for (j in which(npts > 0 & npts < L)) cut[j] = min(sample(0:(L - 2), npts[j]));
    loc = rep(0:(L - 1), 2 * NOFF);
    from0 = (loc <= repEach(cut, L)) == repEach(first, L); // per gamete and locus: T takes a0
    at = repEach(repEach(par * L, L) + loc, NT) * NT + rep(0:(NT - 1), 2 * NOFF * L);
    g = ifelse(repEach(from0, NT), a0[at], a1[at]);
    c0 = g[0:(NOFF * L * NT - 1)];
    c1 = g[(NOFF * L * NT):(2 * NOFF * L * NT - 1)];
"""),
    :record => EidosCode(:record, """
    // record: the parent generation's statistics, every `sampleEvery` generations: QTSim's
    // columns for trait 1, then the trait means, the G matrix (upper triangle, row by row)
    // and the optimum   [Offspring → Offspring]
    if (gen % sampleEvery == 0) {
        zb = rowSums(pheno) / na;
        gv = NULL;
        for (r in 0:(NT - 1)) for (s in r:(NT - 1)) gv = c(gv, cov(drop(G[r, ]), drop(G[s, ])));
        catn(gen + "\\t" + zb[0] + "\\t" + gv[0] + "\\t" + ns + "\\t" + @POP@.getValue("THETA")[0] + "\\t" + @POP@.getValue("OMEGA") + ((DEMES > 1) ? "\\t" + @POP@.id else "") + ((demography != "fixed") ? "\\t" + size(adults) else "") + "\\t" + paste(zb, sep = "\\t") + "\\t" + paste(gv, sep = "\\t") + "\\t" + paste(@POP@.getValue("THETA"), sep = "\\t"));
    }
"""),
    :rescale => EidosCode(:rescale, """
    // rescale: the scheduled change of the fitness width, from the next generation; Ω = OMEGA² S   [Offspring → Offspring]
    if (gen == changeGen) @POP@.setValue("OMEGA", changeSD);
"""),
    :mutate => EidosCode(:mutate, """
    // mutate: Poisson(L·rate) mutations per offspring chromosome, each adding a draw from
    // N(0, MUTCOV) to the NT effects of a random locus   [Offspring → Offspring]
    k0 = rpois(NOFF, L * rate);
    k1 = rpois(NOFF, L * rate);
    if (sum(k0) > 0) {
        at = repEach(0:(NOFF - 1), k0) * L + rdunif(sum(k0), 0, L - 1);
        eff = rmvnorm(sum(k0), rep(0.0, NT), MUTCOV);
        for (j in seqAlong(at)) { x = at[j] * NT + 0:(NT - 1); c0[x] = c0[x] + drop(eff[j, ]); }
    }
    if (sum(k1) > 0) {
        at = repEach(0:(NOFF - 1), k1) * L + rdunif(sum(k1), 0, L - 1);
        eff = rmvnorm(sum(k1), rep(0.0, NT), MUTCOV);
        for (j in seqAlong(at)) { x = at[j] * NT + 0:(NT - 1); c1[x] = c1[x] + drop(eff[j, ]); }
    }
"""),
    :swap => EidosCode(:swap, """
    // swap: the offspring replace their parents   [Offspring → Zygotes]
    kids_@POP@ = @POP@.addEmpty(count = NOFF);
    for (k in seqLen(NOFF)) {
        kids_@POP@[k].setValue("a0", c0[(k * L * NT):(k * L * NT + L * NT - 1)]);
        kids_@POP@[k].setValue("a1", c1[(k * L * NT):(k * L * NT + L * NT - 1)]);
    }
    adults.fitnessScaling = 0.0;
"""))

"The fields of `Params` for several traits, which the demes of a SLiM model share."
const TRAIT_FIELDS = (:NoTraits, :MutCorr, :SelCorr, :MutCov, :SelCov, :MoveDirection)

"""
The constants a SLiM model of several traits adds, as `slim_script`'s `defaults` pairs:
NT, the fitness shape's inverse SINV and the move direction UDIR, and either the
mutational correlation MCORR (M = m·V_E·MCORR) or the mutational covariance MCOV, from
which `initialize()` derives MUTCOV = M/(2·L·rate).
"""
function eidos_trait_constants(p::Params)
    n = p.NoTraits
    Mshape = isempty(p.MutCov) ? ("MCORR" => rowwise(QTLSim.equicorrelation(n, p.MutCorr))) : ("MCOV" => rowwise(p.MutCov))
    ("NT" => n, "SINV" => rowwise(inv(QTLSim.Symmetric(QTLSim.selection_shape(p)))), "UDIR" => QTLSim.move_direction(p), Mshape)
end

"The Eidos that derives MUTCOV in `initialize()`."
eidos_mutcov(p::Params) = isempty(p.MutCov) ?
    "    defineConstant(\"MUTCOV\", matrix(MCORR, nrow = NT) * ((m / (2 * L * rate)) * VE));   // M/(2·L·rate), M = m·V_E·MCORR\n" :
    "    defineConstant(\"MUTCOV\", matrix(MCOV, nrow = NT) / (2 * L * rate));               // M/(2·L·rate)\n"

"The columns a model of several traits adds to the SLiM output's header, as Eidos: trait means, G's upper triangle, the optimum."
trait_header(n) = " + \"" * join(vcat(["\\tz$r" for r in 1:n], ["\\tG$r$s" for r in 1:n for s in r:n],
                                      ["\\ttheta$r" for r in 1:n])) * "\""

"The constants of `slim_script` whose values are vectors in every deme, not one value per deme."
const VECTOR_CONSTANTS = ("MIG", "SINV", "UDIR", "MCORR", "MCOV")

eidos_images(::TraitModel) = eidos_trait_images()

rowwise(A) = vec(permutedims(Matrix{Float64}(A)))
eidos_matrix(A) = "matrix(" * eidos_literal(rowwise(A)) * ", nrow = $(size(A, 1)), byrow = T)"

function eidos_with_params(m::TraitModel, s, c::EidosCode)
    pr = s.params
    get(pr, :CrossoverModel, :qtsim) === :qtsim ||
        throw(ArgumentError("the lowering to Eidos recombines as QTSim does (CrossoverModel = :qtsim); the step $s does not"))
    code = c.code
    n = m.n
    θ, Ω = get(pr, :θ, nothing), get(pr, :Ω, nothing)
    ω = Ω === nothing ? get(pr, :ω, nothing) : 1.0
    ω === nothing || (code = replace(code, "@POP@.getValue(\"OMEGA\")" => eidos_literal(ω)))
    θ === nothing || (code = replace(code, "@POP@.getValue(\"THETA\")" => eidos_literal(θ)))
    # the derived constants, from the step's parameters: the inverse fitness shape, the
    # move's direction, and one mutation's covariance (before the constants m and rate are
    # replaced, so that a step's m and rate enter it)
    subs1 = Dict{String,String}()
    Ω === nothing || (subs1["SINV"] = eidos_literal(rowwise(inv(QTLSim.Symmetric(Ω)))))
    haskey(pr, :SelCorr) && (subs1["SINV"] = eidos_literal(rowwise(inv(QTLSim.Symmetric(QTLSim.equicorrelation(n, pr.SelCorr))))))
    haskey(pr, :MoveDirection) && (subs1["UDIR"] = eidos_literal(pr.MoveDirection ./ sqrt(sum(abs2, pr.MoveDirection))))
    if haskey(pr, :MutCov)
        subs1["MUTCOV"] = "(" * eidos_matrix(pr.MutCov) * " / (2 * L * rate))"
    elseif haskey(pr, :MutCorr)
        subs1["MUTCOV"] = "(" * eidos_matrix(QTLSim.equicorrelation(n, pr.MutCorr)) * " * ((m / (2 * L * rate)) * VE))"
    end
    code = substitute_identifiers(code, subs1)
    subs = Dict{String,String}(EIDOS_CONSTANT[k] => eidos_literal(v) for (k, v) in pairs(pr) if haskey(EIDOS_CONSTANT, k))
    θ === nothing || (subs["moveGen"] = "-1")
    code = substitute_identifiers(code, subs)
    # a step's m or rate alone rescales the deme's MUTCOV = M/(2·L·rate): M ∝ m unless the
    # deme gives M itself (MCOV); the deme's own m and rate are the constants left in place
    if !haskey(pr, :MutCov) && !haskey(pr, :MutCorr) && (haskey(pr, :MutationScalar) || haskey(pr, :MutationRatePerGamete))
        f = "1.0"
        haskey(pr, :MutationRatePerGamete) && (f *= " * (rate / " * eidos_literal(pr.MutationRatePerGamete) * ")")
        haskey(pr, :MutationScalar) &&
            (f *= " * (exists(\"MCOV\") ? 1.0 else " * eidos_literal(pr.MutationScalar) * " / m)")
        code = substitute_identifiers(code, ["MUTCOV" => "(MUTCOV * ($f))"])
    end
    EidosCode(Symbol(string(s)), code)
end
