# The categorical layer: declared here, implemented by the package extension
# QTLSimCatlabExt (ext/), which Julia loads when Catlab is loaded (`import Catlab`).

"""
    lifecycle()

The life cycle as syntax: a named tuple of the Catlab presentation (`presentation`) of
the free category whose objects are the stages of a generation (`Zygotes`, `Adults`,
`Survivors`, `Offspring`), its generating steps (`develop`, `select`, `reproduce`,
`record`, `rescale`, `mutate`, `swap`) and QTSim's generation (`generation`), their
composite `develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅ mutate ⋅ swap`, a morphism
`Zygotes → Zygotes`. Composites type-check: a step can follow another only where the
stages match. [`lower`](@ref) gives a life cycle its meaning in a backend. Requires
Catlab.
"""
function lifecycle end

"""
    lower(f, target; images = Dict(), record = nothing)

The image of the life-cycle morphism `f` under the lowering functor to `target`:
`:kernel` gives an operator on a [`World`](@ref), and `:eidos` a fragment of Eidos code
for SLiM. A functor out of a free category is fixed by the images of the generators, so
the image of a composite is the composite of the images. `record` steps in the kernel
append the parent generation's statistics to `record`, a named tuple of vectors like the
one [`run!`](@ref) returns, or discard them when it is `nothing`. Steps that are not in
[`lifecycle`](@ref), made with Catlab's `Hom(name, dom, codom)`, need an image in
`images`, keyed by name: an operator or a function `(state, rng) -> state` for
`:kernel`, a string of Eidos for `:eidos`. Requires Catlab.
"""
function lower end

"""
    slim_script(f, p::Params; images = Dict())

A complete SLiM model, as Eidos source, that iterates the life cycle `f`, a morphism
`Zygotes → Zygotes`: the image of `f` under [`lower`](@ref) to Eidos runs once per tick
in a `reproduction()` callback. The parameters of `p` are the defaults, which
`slim -d "NAME=value"` overrides (`N`, `VP`, `h2`, `m`, `rate`, `NGEN`, `sampleEvery`,
`moveGen`, `moveSD`, `moveType`, `changeGen`, `changeSD`, `initOmega`, `drift`,
`standardise`). The Eidos images of the steps cover one locus in diploids (`NoQTL = 1`,
`Ploidy = 2`). Requires Catlab.
"""
function slim_script end

"""
    operators()

The kernel images of the steps of the life cycle ([`lifecycle`](@ref)), as the named
tuple `(; develop, select, reproduce, record, rescale, mutate, swap, generation)`. An
operator is called as `op(world, rng)`, composes left to right with Catlab's `⋅`
(`(f ⋅ g)(w, rng)` applies `f` then `g`), and acts deme by deme on a metapopulation
built with `⊗`. `generation` is the image of QTSim's generation,
`develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅ mutate ⋅ swap`, in the order of its
generation loop; here `record` discards what it observes (see [`run_composed`](@ref)).
Requires Catlab.
"""
function operators end

"""
    operator(name::Symbol, act; scope = :deme)

A new operator that applies `act(state, rng)` and returns the state. With
`scope = :deme` it acts on one population (a [`World`](@ref)) and is applied to each
deme of a metapopulation; with `scope = :meta` it acts on the whole metapopulation, as
[`migration`](@ref) does. Requires Catlab.
"""
function operator end

"""
    demes(k)

Stratification into `k` demes: `p ⊗ demes(k)` with `p::Params` gives `k` independent
populations, and `w ⊗ demes(k)` with `w::World` gives `k` copies of `w`. Requires
Catlab.
"""
function demes end

"""
    run_composed(w::World, rng; record = true, cycle = lifecycle().generation, images = Dict())

Run a simulation by applying the kernel image of the life cycle `cycle`, a morphism
`Zygotes → Zygotes`, once per generation ([`lower`](@ref); `images` gives the kernel
images of steps that are not in [`lifecycle`](@ref)). Returns the statistics its
`record` steps collect, as [`run!`](@ref) does; with the default life cycle, for the same
seed, it is the same run as `run!`, bit for bit. Requires Catlab.
"""
function run_composed end

"""
    migrate!(s, rng; rate = 0.01)

Move individuals between the demes of a metapopulation `s`: each individual swaps
places with a random individual of a random other deme with probability `rate`, so the
deme sizes are unchanged. Returns `(; count)`, the number of exchanges. Requires Catlab.
"""
function migrate! end

"""
    migration(rate)

[`migrate!`](@ref) as an operator on the whole metapopulation, for composing into a
generation: `generation ⋅ migration(0.01)`. Requires Catlab.
"""
function migration end

"""
    totalpop(s)

The total number of individuals across the demes of a metapopulation. Requires Catlab.
"""
function totalpop end

"""
    ndemes(s)

The number of demes of a metapopulation. Requires Catlab.
"""
function ndemes end

"""
    population_schema()

The schema of the population ACSet: individuals and alleles, a map `owner` from each
allele to the individual carrying it, allele attribute `effect`, and individual
attributes `env`, `genotype`, `phenotype` and `fitness`. Requires Catlab.
"""
function population_schema end

"""
    lift(w::World)

The parents of `w` as an ACSet on [`population_schema`](@ref). Run
[`develop!`](@ref) first for the genotype and phenotype attributes to be current.
Requires Catlab.
"""
function lift end

"""
    lower!(pop, w::World)

Write the genomes and environmental deviations of the ACSet `pop` into `w`, inverting
[`lift`](@ref). Returns `w`. Requires Catlab.
"""
function lower! end

"""
    genotypes(pop)

Each individual's genotype, computed on the ACSet by summing the effects of the alleles
in the preimage of `owner` (QTSim `CaclGenotype`). Requires Catlab.
"""
function genotypes end

"""
    phenotypes(pop)

Each individual's phenotype on the ACSet, genotype plus environmental deviation.
Requires Catlab.
"""
function phenotypes end

"""
    genome_from_acset(pop, L, P)

The `(L, P, N)` genome array of a population ACSet, inverting [`lift`](@ref). Requires
Catlab.
"""
function genome_from_acset end

const CATEGORICAL = (lifecycle, lower, slim_script, operators, operator, demes,
                     run_composed, migrate!, migration, totalpop, ndemes, population_schema,
                     lift, lower!, genotypes, phenotypes, genome_from_acset)

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if any(f -> exc.f === f, CATEGORICAL) &&
           Base.get_extension(@__MODULE__, :QTLSimCatlabExt) === nothing
            print(io, "\n`", nameof(exc.f), "` is part of QTLSim's categorical layer, ",
                  "which is loaded with Catlab: run `import Catlab` first.")
        end
    end
end
