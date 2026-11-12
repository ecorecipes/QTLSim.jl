# The support matrix (roadmap C7, F11): features × models × backends × closures, generated
# by probing the lowerings' own guards (a combination is unsupported where it throws an
# ArgumentError) and read against the declared validation evidence (`SUPPORT_EVIDENCE`):
# supported = it runs and a validation script compares it with another backend; untested =
# it runs and nothing compares it; unsupported = a guard rejects it. `explain(p; cycle)`
# answers the same question for one parameter set: which model, closure and backends
# apply, and why.

import QTLSim: support_matrix, explain, support_markdown

# the probes are tiny: a few generations of a few individuals
const PROBE_BASE = Params(PopSize = 30, NoQTL = 5, Ploidy = 2, NoGenerations = 4, PhenotypeVar = 1.0, Heritability = 0.5,
                          MutationScalar = 0.01, MutationRatePerGamete = 0.01, InitFitnessSD = 2.0, Sample = 1)

"The models of the support matrix, by name, each with its base parameters."
const PROBE_MODELS = [
    ("one trait", PROBE_BASE),
    ("two traits", Params(PROBE_BASE; NoTraits = 2, MutCorr = 0.3)),
    ("reaction norm", Params(PROBE_BASE; NoTraits = 2, ReactionNorm = true, Plasticity = 0.2))]

const MOVE = (MoveFitnessSDGen = 2, MoveFitnessSD = 0.1)

# A feature: a name, the change to a model's parameters, the number of demes and how the
# demes differ (`perdeme`), whether a step carries its own parameters (`step`), and options
# of the moments lowering (`between`).
const PROBE_FEATURES = [
    (name = "QTSim's generation", f = p -> p),
    (name = "optimum jump (MoveType = :phenotype)", f = p -> Params(p; MOVE..., MoveType = :phenotype)),
    (name = "continuous move", f = p -> Params(p; MOVE..., MoveType = :continuous)),
    (name = "red-noise optimum", f = p -> Params(p; MOVE..., MoveType = :red_noise, MoveNoiseSD = 0.1, MoveAutocorr = 0.5)),
    (name = "periodic optimum", f = p -> Params(p; MOVE..., MoveType = :periodic, MoveAmplitude = 0.5, MovePeriod = 10.0)),
    (name = "width change (ChangeFitnessGen)", f = p -> Params(p; ChangeFitnessGen = 2, ChangeFitnessSD = 3.0)),
    (name = "StandardisePhenotype", f = p -> Params(p; StandardisePhenotype = true)),
    (name = "neutral (Drift)", f = p -> Params(p; Drift = true)),
    (name = "haploid", f = p -> Params(p; Ploidy = 1)),
    (name = "every crossover point (:all)", f = p -> Params(p; CrossoverModel = :all)),
    (name = "Beverton–Holt demography", f = p -> Params(p; Demography = :beverton_holt, GrowthRate = 2.0)),
    (name = "fecundity selection", f = p -> Params(p; Selection = :fecundity)),
    (name = "plasticity, fixed slope", f = p -> Params(p; MOVE..., MoveType = :continuous, Plasticity = 0.5, OptimumSlope = 2.0)),
    (name = "plasticity with a cost", f = p -> Params(p; MOVE..., MoveType = :continuous, Plasticity = 0.5, PlasticityCostSD = 1.0)),
    (name = "plasticity with a delayed cue", f = p -> Params(p; MOVE..., MoveType = :continuous, Plasticity = 0.5, CueLag = 2)),
    (name = "a step's own parameters (withparams θ)", f = p -> p, step = true),
    (name = "3 demes, island migration", f = p -> p, demes = 3),
    (name = "3 demes, stepping stone (matrix)", f = p -> p, demes = 3, matrix = true),
    (name = "3 demes, their own Params", f = p -> p, demes = 3, perdeme = true),
    (name = "3 demes, drift among demes (between)", f = p -> p, demes = 3, between = true)]

const PROBE_CLOSURES = (:aggregate, :multilocus, :rare_alleles, :rare_local, :spectrum)
const SUPPORT_COLUMNS = (:kernel, :slim, PROBE_CLOSURES...)

# Validation evidence: (feature, model, column) => the script that compares it with
# another backend (QTSim, SLiM, the kernel or theory) and found agreement. From CLAUDE.md
# and the scripts' logs; a column missing here is `untested` when it runs.
const SUPPORT_EVIDENCE = Dict{Tuple{String,String,Symbol},String}(
    ("QTSim's generation", "one trait", :kernel) => "gsl_locked.jl, compare_to_qtsim.jl",
    ("optimum jump (MoveType = :phenotype)", "one trait", :kernel) => "compare_to_qtsim.jl",
    ("continuous move", "one trait", :kernel) => "published_sweep.jl",
    ("red-noise optimum", "one trait", :kernel) => "slim/compare_fluctuating.jl",
    ("periodic optimum", "one trait", :kernel) => "slim/compare_fluctuating.jl",
    ("neutral (Drift)", "one trait", :kernel) => "fst.jl, coalescence.jl",
    ("haploid", "one trait", :kernel) => "rare_alleles_traits.jl (H)",
    ("every crossover point (:all)", "one trait", :kernel) => "moments_loci.jl",
    ("Beverton–Holt demography", "one trait", :kernel) => "slim/compare_demography.jl",
    ("fecundity selection", "one trait", :kernel) => "slim/compare_selection.jl",
    ("plasticity, fixed slope", "one trait", :kernel) => "slim/compare_plasticity.jl, plasticity_rates.jl",
    ("plasticity with a cost", "one trait", :kernel) => "plasticity_rates.jl",
    ("a step's own parameters (withparams θ)", "one trait", :kernel) => "slim/compare_demes.jl",
    ("3 demes, island migration", "one trait", :kernel) => "slim/compare_migration.jl, fst.jl, moments_demes.jl",
    ("3 demes, stepping stone (matrix)", "one trait", :kernel) => "slim/compare_migration.jl, fst.jl",
    ("3 demes, their own Params", "one trait", :kernel) => "slim/compare_demes.jl",
    ("QTSim's generation", "one trait", :slim) => "slim/compare_to_qtsim.jl, slim/three_way_confirm.jl",
    ("optimum jump (MoveType = :phenotype)", "one trait", :slim) => "slim/compare_to_qtsim.jl",
    ("red-noise optimum", "one trait", :slim) => "slim/compare_fluctuating.jl",
    ("periodic optimum", "one trait", :slim) => "slim/compare_fluctuating.jl",
    ("Beverton–Holt demography", "one trait", :slim) => "slim/compare_demography.jl",
    ("fecundity selection", "one trait", :slim) => "slim/compare_selection.jl",
    ("plasticity, fixed slope", "one trait", :slim) => "slim/compare_plasticity.jl",
    ("a step's own parameters (withparams θ)", "one trait", :slim) => "slim/compare_demes.jl",
    ("3 demes, island migration", "one trait", :slim) => "slim/compare_migration.jl",
    ("3 demes, stepping stone (matrix)", "one trait", :slim) => "slim/compare_migration.jl",
    ("3 demes, their own Params", "one trait", :slim) => "slim/compare_demes.jl",
    ("QTSim's generation", "one trait", :multilocus) => "moments_loci.jl",
    ("continuous move", "one trait", :multilocus) => "moments_loci.jl (S)",
    ("every crossover point (:all)", "one trait", :multilocus) => "moments_loci.jl",
    ("QTSim's generation", "one trait", :rare_alleles) => "moments_regimes.jl",
    ("haploid", "one trait", :rare_alleles) => "rare_alleles_traits.jl (H)",
    ("QTSim's generation", "one trait", :spectrum) => "vg_move.jl",
    ("continuous move", "one trait", :spectrum) => "vg_move.jl",
    ("3 demes, island migration", "one trait", :multilocus) => "moments_demes.jl",
    ("3 demes, drift among demes (between)", "one trait", :multilocus) => "moments_demes.jl",
    ("3 demes, drift among demes (between)", "one trait", :rare_local) => "local_adaptation.jl",
    ("QTSim's generation", "two traits", :kernel) => "rare_alleles_traits.jl, slim/compare_traits.jl",
    ("continuous move", "two traits", :kernel) => "slim/compare_traits.jl",
    ("red-noise optimum", "two traits", :kernel) => "slim/compare_traits_fluctuating.jl",
    ("periodic optimum", "two traits", :kernel) => "slim/compare_traits_fluctuating.jl",
    ("3 demes, island migration", "two traits", :kernel) => "trait_moments_demes.jl, rare_alleles_traits.jl (D)",
    ("QTSim's generation", "two traits", :slim) => "slim/compare_traits.jl",
    ("continuous move", "two traits", :slim) => "slim/compare_traits.jl",
    ("red-noise optimum", "two traits", :slim) => "slim/compare_traits_fluctuating.jl",
    ("periodic optimum", "two traits", :slim) => "slim/compare_traits_fluctuating.jl",
    ("QTSim's generation", "two traits", :multilocus) => "moments_loci.jl (T), traits_kernel_vs_moments.jl",
    ("QTSim's generation", "two traits", :rare_alleles) => "rare_alleles_traits.jl (T)",
    ("QTSim's generation", "two traits", :rare_local) => "rare_alleles_traits.jl (T; equals :rare_alleles for one deme)",
    ("3 demes, island migration", "two traits", :multilocus) => "trait_moments_demes.jl",
    ("3 demes, drift among demes (between)", "two traits", :multilocus) => "trait_moments_demes.jl",
    ("3 demes, drift among demes (between)", "two traits", :rare_local) => "local_adaptation.jl (D1, D2)",
    ("continuous move", "reaction norm", :kernel) => "reaction_norm.jl",
    ("continuous move", "reaction norm", :multilocus) => "reaction_norm.jl")

# the parameters, life cycle and options of a feature in a model
function probe_case(feat, p0)
    L = lifecycle(); g = L.generation
    p = feat.f(p0)
    k = get(feat, :demes, 1)
    θ = p.NoTraits == 1 || p.ReactionNorm ? 0.5 : [0.5, -0.5]
    cycle = get(feat, :step, false) ? withparams(g, p; θ) :
            k == 1 ? g : foldl(⊗, fill(g, k)) ⋅ L.migrate(k)
    ps = get(feat, :perdeme, false) ? [Params(p; InitFitnessSD = 1.0 + a) for a in 1:k] : p
    rate = k == 1 ? nothing : get(feat, :matrix, false) ? stepping_stone(k, 0.05) : 0.01
    (; p, ps, k, cycle, rate, between = get(feat, :between, false))
end

# run a column on a case: :ok, or (:unsupported, message) for an ArgumentError, or (:error, message)
function probe(column::Symbol, c)
    try
        if column === :kernel
            c.between && return (:na, "between = true is an option of the moments lowering")
            if c.k == 1
                run_composed(World(c.p), Xoshiro(1); cycle = c.cycle)
            else
                run_composed((c.ps isa Params ? c.ps : collect(Params, c.ps)) ⊗ demes(c.k), Xoshiro(1); cycle = c.cycle,
                             migration_rate = c.rate)
            end
        elseif column === :slim
            c.between && return (:na, "between = true is an option of the moments lowering")
            c.k == 1 ? slim_script(c.cycle, c.p) : slim_script(c.cycle, c.ps isa Params ? fill(c.ps, c.k) : c.ps; migration_rate = c.rate)
        else
            kw = c.k == 1 ? (; loci = column) : (; loci = column, migration_rate = c.rate, between = c.between)
            c.k == 1 ? run_moments(c.p; cycle = c.cycle, kw...) :
                       run_moments(c.ps isa Params ? fill(c.ps, c.k) : c.ps; cycle = c.cycle, kw...)
        end
        (:ok, "")
    catch e
        e isa ArgumentError && return (:unsupported, e.msg)
        e isa ErrorException && occursin("ArgumentError", e.msg) && return (:unsupported, e.msg)
        (:error, sprint(showerror, e))
    end
end

function support_matrix(; features = PROBE_FEATURES, models = PROBE_MODELS)
    rows = NamedTuple[]
    for (mname, p0) in models, feat in features
        c = try
            probe_case(feat, p0)
        catch e
            e isa ArgumentError || rethrow()
            nothing
        end
        for col in SUPPORT_COLUMNS
            st, msg = c === nothing ? (:unsupported, "the parameters are rejected") : probe(col, c)
            ev = get(SUPPORT_EVIDENCE, (feat.name, mname, col), "")
            status = st === :ok ? (isempty(ev) ? :untested : :supported) : st
            push!(rows, (feature = feat.name, model = mname, column = col, status, evidence = ev, message = msg))
        end
    end
    rows
end

const STATUS_MARK = Dict(:supported => "✓", :untested => "○", :unsupported => "✗", :na => "–", :error => "!")

function support_markdown(rows = support_matrix())
    io = IOBuffer()
    println(io, """
    # Support matrix

    What QTLSim can run, and what has been checked against another backend, by feature,
    model, backend and closure of the moments lowering. **Generated** by
    `scripts/support_matrix.jl` from `support_matrix()` (ext/QTLSimCatlabExt/support.jl):
    every cell is a probe of the code (four generations of thirty individuals), so a cell
    is unsupported exactly where the lowering's own guard throws an `ArgumentError`, and
    a cell that runs is *supported* when the declared evidence names a validation script
    that compares it with another backend (QTSim, SLiM, the kernel or theory), *untested*
    otherwise. A test checks that this page is up to date.

    ✓ supported and validated · ○ runs, not validated · ✗ unsupported (a guard; see the
    note) · – not applicable. Columns: the kernel, the generated SLiM model (the probe
    generates the script; it does not run SLiM), and the moments lowering with each closure
    (`loci`; `:auto` picks one of them, see [`QTLSim.auto_closure`](@ref) and
    [`explain`](@ref)).
    """)
    notes = String[]
    note(s) = (i = findfirst(==(s), notes); i === nothing ? (push!(notes, s); length(notes)) : i)
    for mname in unique(r.model for r in rows)
        println(io, "## ", uppercasefirst(mname), "\n")
        println(io, "| feature | kernel | SLiM | aggregate | multilocus | rare alleles | rare local | spectrum |")
        println(io, "|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|")
        for feat in unique(r.feature for r in rows if r.model == mname)
            cells = String[]
            for col in SUPPORT_COLUMNS
                r = only(x for x in rows if x.model == mname && x.feature == feat && x.column === col)
                m = STATUS_MARK[r.status]
                r.status === :supported && (m *= " [$(note("validated by " * r.evidence))]")
                r.status in (:unsupported, :error) && (m *= " [$(note(first(replace(r.message, '\n' => ' ', '|' => '/'), 160)))]")
                push!(cells, m)
            end
            println(io, "| ", feat, " | ", join(cells, " | "), " |")
        end
        println(io)
    end
    println(io, "## Notes\n")
    for (i, s) in enumerate(notes)
        println(io, i, ". ", s)
    end
    String(take!(io))
end

# --- explain ----------------------------------------------------------------------------

function explain(p; cycle = nothing, loci::Symbol = :auto, migration_rate = 0.0, between::Bool = false, io = stdout)
    ps = p isa Params ? nothing : collect(Params, p)
    p1 = p isa Params ? p : first(ps)
    L = lifecycle()
    k = ps !== nothing ? length(ps) : cycle === nothing ? 1 : nfactors(dom(cycle))
    cyc = cycle !== nothing ? cycle : k == 1 ? L.generation : foldl(⊗, fill(L.generation, k)) ⋅ L.migrate(k)
    model = model_of(p1)
    mf = p1.NoTraits == 1 ? mutate_fields(cyc) !== nothing :
         length(unique([NamedTuple(x => v for (x, v) in pairs(stepparams(g)) if x in (:MutationScalar, :MutationRatePerGamete, :MutCorr, :MutCov))
                        for g in stepsof(cyc) if !(g isa ObExpr) && stepname(g) === :mutate])) <= 1
    ctx = QTLSim.closure_context(p1; ndemes = k, uniform_mutation = mf, between)
    reasons = QTLSim.closure_reasons(ctx)
    chosen = loci === :auto ? QTLSim.auto_closure(ctx) : loci
    short(q) = Params(q; NoGenerations = min(q.NoGenerations, 3))
    c = (p = short(p1), ps = ps === nothing ? short(p1) : map(short, ps), k, cycle = cyc, rate = migration_rate, between)
    backends = [(backend = b, status = first(s), message = last(s))
                for (b, s) in ((:kernel, probe(:kernel, c)), (:slim, probe(:slim, c)), (Symbol("moments ", chosen), probe(chosen, c)))]
    out = (model, ndemes = k, context = ctx, closures = reasons, closure = chosen,
           closure_reason = loci === :auto ? QTLSim.closure(chosen).auto_text : "asked for (loci = :$loci)", backends)
    if io !== nothing
        println(io, "model: ", model, "; ", k, k == 1 ? " deme" : " demes", "; ploidy ", p1.Ploidy,
                "; one law of mutation effects: ", mf, between ? "; between = true" : "")
        println(io, "closures of the moments lowering (domain: where run_moments accepts it; auto: where :auto may pick it):")
        for r in reasons
            println(io, "  ", rpad(":" * string(r.name), 15), r.domain ? "in domain" : "not in domain", ", ",
                    r.auto ? "auto" : "not auto", r.domain ? "" : " — domain: " * r.domain_text)
        end
        println(io, "closure: :", chosen, loci === :auto ? " (loci = :auto: the finest closure whose automatic domain contains this; " *
                    QTLSim.closure(chosen).auto_text * ")" : " (asked for)")
        println(io, "backends (probed on 3 generations):")
        for b in backends
            println(io, "  ", rpad(string(b.backend), 22), b.status === :ok ? "runs" : string(b.status) * ": " * first(b.message, 200))
        end
    end
    out
end
