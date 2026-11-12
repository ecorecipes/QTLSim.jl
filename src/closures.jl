# The closures of the moments lowering, declared as data (roadmap C7). A closure is the
# model of the loci with which the moments lowering (`run_moments(…; loci)`) carries a
# population: the aggregate recursion, the multilocus Gaussian model, the rare-alleles
# closure, the rare-alleles closure with locally favoured alleles, and the allele-frequency
# spectrum. They form an order of refinements (`refines`), each arrow a reduction that
# holds exactly or in a limit under a stated condition (the reductions themselves, with the
# tests they generate, are declared in ext/QTLSimCatlabExt/morphisms.jl, `closure_order`).
#
# Each closure declares its domain, the contexts (model, number of demes, ploidy, one law of
# mutation effects, drift among demes) in which `run_moments` accepts it, and its automatic
# domain, where `loci = :auto` may choose it. `:auto` is then a rule read from the data:
# the finest closure (highest `rank`) whose automatic domain contains the context
# (`auto_closure`). The tests check the domains against the lowering's own `ArgumentError`
# guards and the rule against the rule written by hand before these declarations.

"""
    QTLSim.ClosureContext(model, ndemes, ploidy, uniform_mutation, between)

What decides which closures of the moments lowering apply: the `model` (`model_of(p)`),
the number of demes of the life cycle, the ploidy, whether the `mutate` steps of the life
cycle carry one law of mutation effects (`uniform_mutation`: none or all of them carry the
same `MutationScalar`, `MutationRatePerGamete` and, for several traits, `MutCorr`,
`MutCov`) and whether drift among demes is tracked (`between = true`). See
[`QTLSim.closure_context`](@ref).
"""
struct ClosureContext
    model::Model
    ndemes::Int
    ploidy::Int
    uniform_mutation::Bool
    between::Bool
end

"""
    QTLSim.closure_context(p; ndemes = 1, uniform_mutation = true, between = false)

The [`QTLSim.ClosureContext`](@ref) of a `Params` `p` (every deme's, or the first deme's of
a metapopulation, whose demes share a model and a ploidy).
"""
closure_context(p::Params; ndemes::Integer = 1, uniform_mutation::Bool = true, between::Bool = false) =
    ClosureContext(model_of(p), ndemes, p.Ploidy, uniform_mutation, between)

ntraits(c::ClosureContext) = ntraits(c.model)

"""
    QTLSim.ClosureDecl

The declaration of a closure of the moments lowering: its `name` (the value of `loci`),
its `rank` in the order of refinements (finer is higher), the closures it `refines`, its
`domain` (a predicate on a [`QTLSim.ClosureContext`](@ref): where `run_moments` accepts
it), its automatic domain `auto` (where `loci = :auto` may choose it; a subset of the
domain, but for `:multilocus`, the fallback everywhere), each with a sentence saying it in words, the source files, and the validation
that supports it. [`QTLSim.closures`](@ref) lists them.
"""
struct ClosureDecl
    name::Symbol
    rank::Int
    refines::Vector{Symbol}
    domain::Function
    domain_text::String
    auto::Function
    auto_text::String
    source::String
    validated::String
end

Base.show(io::IO, c::ClosureDecl) = print(io, "ClosureDecl(:", c.name, ", rank ", c.rank, ")")

one_or_traits(c::ClosureContext) = c.model isa OneTraitModel || c.model isa TraitModel
# drift among demes (`between = true`) tracks the haplotype moments of diploids
between_ok(c::ClosureContext) = !c.between || c.ploidy == 2

const CLOSURES = ClosureDecl[
    ClosureDecl(:aggregate, 1, Symbol[],
        c -> !c.between,
        "every model, without drift among demes (`between = true` needs haplotype moments)",
        c -> false,
        "never: it is exact only for one haploid locus (and haploid clones); ask for it",
        "ext/QTLSimCatlabExt/moments.jl (`cov = nothing`)",
        "theory: gaussian_vg_exact, gaussian_G (test/runtests.jl)"),
    ClosureDecl(:multilocus, 2, [:aggregate],
        c -> between_ok(c),
        "every model, every number of demes, both ploidies (`between = true` only for diploids)",
        c -> true,
        "always (the fallback: the finest closure that needs no law of mutation effects)",
        "src/multilocus.jl",
        "validate/moments_loci.jl, validate/moments_demes.jl, validate/trait_moments_demes.jl"),
    ClosureDecl(:rare_alleles, 3, [:multilocus],
        c -> one_or_traits(c) && c.uniform_mutation && between_ok(c),
        "one trait or several (not a heritable reaction norm), when the mutate steps carry one law of mutation effects",
        c -> one_or_traits(c) && c.uniform_mutation && between_ok(c),
        "its whole domain",
        "src/rare_alleles.jl, src/rare_alleles_traits.jl",
        "validate/moments_regimes.jl, validate/rare_alleles_traits.jl"),
    ClosureDecl(:rare_local, 4, [:rare_alleles],
        c -> one_or_traits(c) && c.uniform_mutation && between_ok(c),
        "as `:rare_alleles` (it equals `:rare_alleles` for one deme, a common optimum or no migration)",
        c -> one_or_traits(c) && c.uniform_mutation && between_ok(c) && (ntraits(c) > 1 || c.ndemes > 1),
        "several traits, or one trait in several demes; for one trait in one deme it equals `:rare_alleles`, " *
        "the coarser name is kept",
        "src/local_alleles.jl",
        "validate/local_adaptation.jl"),
    ClosureDecl(:spectrum, 5, [:rare_alleles],
        c -> c.model isa OneTraitModel && c.ploidy == 2 && c.ndemes == 1 && c.uniform_mutation && !c.between,
        "one diploid trait, one deme, one law of mutation effects, no `between`",
        c -> false,
        "never: its V_G is too high in the Gaussian regime (validate/vg_move.jl); ask for it under a moving optimum " *
        "when 4Nμ ≲ 1 per locus",
        "src/spectrum.jl, ext/QTLSimCatlabExt/spectrum.jl",
        "validate/vg_move.jl")]

"""
    QTLSim.closures() -> Vector{ClosureDecl}

The declared closures of the moments lowering ([`QTLSim.ClosureDecl`](@ref)), from the
coarsest to the finest: `:aggregate`, `:multilocus`, `:rare_alleles`, `:rare_local`,
`:spectrum`.
"""
closures() = CLOSURES

"The declaration of the closure `name`."
function closure(name::Symbol)
    i = findfirst(c -> c.name === name, CLOSURES)
    i === nothing && throw(ArgumentError("no closure :$name; the closures are $(join((":$(c.name)" for c in CLOSURES), ", "))"))
    CLOSURES[i]
end

"""
    QTLSim.auto_closure(ctx::ClosureContext) -> Symbol

What `loci = :auto` resolves to in the context `ctx`, read from the declarations
([`QTLSim.closures`](@ref)): the finest closure whose automatic domain contains `ctx`.
"""
function auto_closure(ctx::ClosureContext)
    best = nothing
    for c in CLOSURES
        c.auto(ctx) && (best === nothing || c.rank > best.rank) && (best = c)
    end
    best === nothing && throw(ArgumentError("no closure applies in $ctx"))   # unreachable: :multilocus always applies
    best.name
end

"""
    QTLSim.closure_reasons(ctx) -> Vector{NamedTuple}

For each declared closure, whether it applies in `ctx` (`domain`), whether `:auto` may
choose it (`auto`), and the declaration's sentences: the reasons behind
[`QTLSim.auto_closure`](@ref), as [`explain`](@ref) prints them.
"""
closure_reasons(ctx::ClosureContext) =
    [(name = c.name, rank = c.rank, domain = c.domain(ctx), auto = c.auto(ctx), domain_text = c.domain_text,
      auto_text = c.auto_text) for c in CLOSURES]

# --- morphisms between models, declared as data (roadmap C7) --------------------------
#
# The declarations themselves, whose checks run the lowerings, are built by the categorical
# layer (ext/QTLSimCatlabExt/morphisms.jl): `model_morphisms()` and `closure_order()`.

"""
    QTLSim.Reduction(name, backend, exactness, tolerance, check; covers = "")

A reduction that a morphism between models (or an arrow of the order of closures) must
satisfy: `check()` returns two values, the two sides of the reduction, computed by
`backend` (`:kernel`, `:moments`, `:theory`); `exactness` is `:bitwise` (`isequal`) or
`:tolerance` (largest relative difference at most `tolerance`). `covers` names the hand
test of test/runtests.jl it repeats, if any.
"""
struct Reduction
    name::String
    backend::Symbol
    exactness::Symbol
    tolerance::Float64
    check::Function
    covers::String
end
Reduction(name, backend, exactness, tolerance, check; covers = "") =
    Reduction(String(name), backend, exactness, Float64(tolerance), check, String(covers))

"""
    QTLSim.ModelMorphism(name, from, to, kind, statement, condition, lean, reductions)

A map between models of QTLSim (or between closures of the moments lowering) declared as
data: `kind` is `:embedding`, `:projection`, `:restriction`, `:conjugacy` or `:refinement`;
`statement` says what it maps, `condition` when its reductions hold, `lean` the Lean
theorem behind it (or what is missing), and `reductions` the [`QTLSim.Reduction`](@ref)s
that must hold, from which the tests are generated ([`QTLSim.check_reduction`](@ref)).
"""
struct ModelMorphism
    name::String
    from::String
    to::String
    kind::Symbol
    statement::String
    condition::String
    lean::String
    reductions::Vector{Reduction}
end

Base.show(io::IO, m::ModelMorphism) = print(io, "ModelMorphism(\"", m.name, "\": ", m.from, " → ", m.to, ", ",
                                            length(m.reductions), " reduction", length(m.reductions) == 1 ? "" : "s", ")")

"""
    QTLSim.check_reduction(r::Reduction) -> NamedTuple

Run a declared reduction: `(pass, difference)`, the largest relative difference of its two
sides (0 when bitwise equal).
"""
function check_reduction(r::Reduction)
    a, b = r.check()
    d = isequal(a, b) ? 0.0 : relerr_max(reduce(vcat, (vec(collect(Float64, x)) for x in flat_values(a)); init = Float64[]),
                                          reduce(vcat, (vec(collect(Float64, x)) for x in flat_values(b)); init = Float64[]))
    pass = r.exactness === :bitwise ? isequal(a, b) : d <= r.tolerance
    (pass = pass, difference = d)
end
flat_values(x::Number) = [[x]]
flat_values(x::AbstractArray{<:Number}) = [x]
flat_values(x::AbstractArray) = reduce(vcat, (flat_values(y) for y in x); init = Any[])
flat_values(x::Tuple) = reduce(vcat, (flat_values(y) for y in x); init = Any[])

"""
    model_morphisms() -> Vector{ModelMorphism}

The declared maps between QTLSim's models ([`QTLSim.ModelMorphism`](@ref)): one trait into
two independent traits, the projection of `n` traits to trait 1 (QTSim's columns), `k`
demes to one, a plastic population to a plain one at rate `(1 − φ)k`; each with the
reductions that must hold. [`QTLSim.closure_order`](@ref) gives the arrows between the
closures of the moments lowering. Requires Catlab.
"""
function model_morphisms end

"""
    QTLSim.closure_order() -> Vector{ModelMorphism}

The order of refinements of the closures of the moments lowering ([`QTLSim.closures`](@ref)),
as arrows from the finer closure to the coarser one it reduces to, each with its condition
and reductions: `:rare_local` → `:rare_alleles` (one deme, a common optimum, no
migration), `:rare_alleles` → `:multilocus` (many small mutations, a limit), `:multilocus`
→ `:aggregate` (haploid, no drift), `:spectrum` → `:rare_alleles` (fixed optimum, rare
mutations; approximate). Requires Catlab.
"""
function closure_order end

# --- the support matrix and explain (implemented in ext/QTLSimCatlabExt/support.jl) ---

"""
    support_matrix(; features = all, models = all) -> Vector{NamedTuple}

The support matrix of QTLSim: for each feature (a change of the parameters or of the life
cycle), model (one trait, two traits, a heritable reaction norm) and column (the kernel,
the generated SLiM model, and the moments lowering with each closure), the `status`:
`:supported` (it runs, and `evidence` names the validation script that compares it with
another backend), `:untested` (it runs; nothing compares it), `:unsupported` (a guard
rejects it with an `ArgumentError`, whose `message` is kept), `:na` (not a question for
that column) or `:error` (any other error: a bug). Generated by probing the lowerings on
tiny populations (4 generations of 30), so the guards are the code's own, and read
against the declared evidence. [`support_markdown`](@ref) writes it as the page
docs/src/support.md. Requires Catlab.
"""
function support_matrix end

"""
    support_markdown(rows = support_matrix()) -> String

The support matrix as a Markdown page (docs/src/support.md): a table per model, a mark per
cell (✓ supported and validated, ○ runs but untested, ✗ unsupported, – not applicable),
with the validation scripts and the guards' messages as numbered notes.
"""
function support_markdown end

"""
    explain(p; cycle = nothing, loci = :auto, migration_rate = 0.0, between = false, io = stdout) -> NamedTuple

Which model, closure and backends apply to `p` (one `Params`, or one per deme) with the
life cycle `cycle` (QTSim's generation by default; for several `Params`, their demes side
by side and coupled by migration), and why: the model (`model_of(p)`), the context of the
closures ([`QTLSim.closure_context`](@ref)), each declared closure's domain and automatic
domain with the declaration's reasons, the closure `loci = :auto` resolves to
([`QTLSim.auto_closure`](@ref)), and a probe of each backend (the kernel, the generated
SLiM model, the moments lowering with that closure) on a short copy of `p`, with the
guard's message where one rejects it. Prints a summary to `io` (`nothing` to print
nothing) and returns it as a named tuple. Requires Catlab.
"""
function explain end

