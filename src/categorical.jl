# The categorical layer: declared here, implemented by the package extension
# QTLSimCatlabExt (ext/), which Julia loads when Catlab is loaded (`import Catlab`).

"""
    lifecycle()

The life cycle as syntax: a named tuple of the Catlab presentation (`presentation`) of
the free symmetric monoidal category whose objects are the stages of a generation
(`Zygotes`, `Adults`, `Survivors`, `Offspring`), its generating steps (`develop`,
`select`, `reproduce`, `record`, `rescale`, `mutate`, `swap`), QTSim's generation
(`generation`), their composite `develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅
mutate ⋅ swap`, a morphism `Zygotes → Zygotes`, and `migrate`, where `migrate(k)` is the
step that exchanges offspring among `k` demes, `Zygotes^⊗k → Zygotes^⊗k`. Composites
type-check: a step can follow another only where the stages match. Catlab's `⊗` places
life cycles side by side, one per deme: `generation ⊗ generation ⋅ migrate(2)` is a
generation of two demes coupled by migration. [`lower`](@ref) gives a life cycle its
meaning in a backend. Requires Catlab.
"""
function lifecycle end

"""
    lower(f, target; images = Dict(), record = nothing, migration_rate = 0.0)

The image of the life-cycle morphism `f` under the lowering functor to `target`:
`:kernel` gives an operator on a [`World`](@ref), `:eidos` a fragment of Eidos code for
SLiM, in which `@POP@` stands for the deme's subpopulation, and `:moments` an update of
the moments of a Gaussian population, the predictions of theory (see
[`run_moments`](@ref); one deme). A morphism between
products of stages, one per deme, has an image on several demes: an operator on demes of
a metapopulation (`DemeOp`), or Eidos for several subpopulations (`EidosMeta`). A functor
out of a free category is fixed by the images of the generators, so the image of a
composite or a product is the composite or product of the images. `record` steps in the
kernel append the parent generation's statistics to `record`, a named tuple of vectors
like the one [`run!`](@ref) returns, or discard them when it is `nothing`; `migrate` steps
exchange individuals among demes at `migration_rate` (in Eidos, the constant `migRate`), or
along a migration graph when it is a matrix ([`island_migration`](@ref),
[`stepping_stone`](@ref); in Eidos, the constant `MIG`).
Steps that are not in [`lifecycle`](@ref), made with Catlab's `Hom(name, dom, codom)`,
need an image in `images`, keyed by name: an operator or a function `(state, rng) -> state`
for `:kernel`, a string of Eidos for `:eidos`. Eidos steps keep a deme's state on its
subpopulation: `@POP@.getValue("THETA")` and `@POP@.getValue("OMEGA")` are its optimum
and fitness width, set with `setValue`, and the steps before `select` leave `adults`,
`genotype` and `pheno` defined. Requires Catlab.
"""
function lower end

"""
    slim_script(f, p::Params; images = Dict(), migration_rate = 0.0)

A complete SLiM model, as Eidos source, that iterates the life cycle `f`, a morphism
`Zygotes → Zygotes`, or `Zygotes^⊗k → Zygotes^⊗k` for `k` demes, each a subpopulation
of `p.PopSize` individuals: the image of `f` under [`lower`](@ref) to Eidos runs once per
tick in a `reproduction()` callback. The parameters of `p` and `migration_rate` are the
defaults, which `slim -d "NAME=value"` overrides (`N`, `VP`, `h2`, `m`, `rate`, `L`,
`RECOMB`, `NGEN`, `sampleEvery`, `moveGen`, `moveSD`, `moveType`, `changeGen`,
`changeSD`, `initOmega`, `drift`, `standardise`, `migRate`). The Eidos images of the
steps cover diploids (`Ploidy = 2`) with any number of loci, recombining as QTSim does
(`CrossoverModel = :qtsim`); with several demes, each row of output ends with the deme.
The life cycle must return each deme to its place: braidings must cancel. Requires
Catlab.
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
populations, `[p1, …, pk] ⊗ demes(k)` gives demes with their own parameters (for example
the width of selection, or the size; the same `NoGenerations`), and `w ⊗ demes(k)` with
`w::World` gives `k` copies of `w`. Requires Catlab.
"""
function demes end

"""
    run_composed(w::World, rng; record = true, cycle = lifecycle().generation, images = Dict())
    run_composed(s, rng; cycle, record = true, images = Dict(), migration_rate = 0.0)

Run a simulation by applying the kernel image of the life cycle `cycle`, a morphism
`Zygotes → Zygotes`, once per generation ([`lower`](@ref); `images` gives the kernel
images of steps that are not in [`lifecycle`](@ref)). Returns the statistics its
`record` steps collect, as [`run!`](@ref) does; with the default life cycle, for the same
seed, it is the same run as `run!`, bit for bit. For a metapopulation `s = p ⊗ demes(k)`,
`cycle` is a life cycle on `k` demes, such as `g ⊗ g ⊗ g ⋅ lifecycle().migrate(3)` for
`k = 3`, and the statistics have a row per deme and sampled generation (`deme`), with
`extinct_gen` per deme. Requires Catlab.
"""
function run_composed end

"""
    run_moments(p::Params; cycle = lifecycle().generation, images = Dict(), zbar = 0.0, VG = 0.0)

The life cycle `cycle` lowered to the moments of a Gaussian population (`lower(cycle, :moments)`)
and iterated for `NoGenerations`: the deterministic, infinite-population predictions of
quantitative-genetic theory for the same morphism the kernel and SLiM simulate. Each
step's image updates the means and variances:
- `develop`: `V_P = V_G + V_E`;
- `select`: moves the optimum as `MoveType` says, then applies Gaussian selection of
  width ω to `z ~ N(z̄, V_P)`. This gives the mean fitness, the expected survivors
  `N·W̄`, the selection differential, and the survivors' genotypic mean and variance
  `z̄ + V_G λ/s` and `V_G − V_G²/s`, with `λ = θ − z̄` and `s = ω² + V_P`;
- `reproduce`: passes those to the offspring (the continuum-of-alleles approximation, no
  drift);
- `mutate`: adds `m·V_E`.

The expected population size is `PopSize`, or under `Demography = :beverton_holt` the
mean of [`noffspring`](@ref) for the expected survivors, `R n/(1 + (R − 1) n/K)`.

Returns the columns of [`run!`](@ref), predicted (`nsurv` rounded). The genetic variance
converges to [`gaussian_vg_exact`](@ref); under a continuous move the lag `θ − z̄`
converges to [`steady_lag`](@ref). The population is extinct when fewer than two
individuals (one, if haploid) are expected to survive. One deme only. Requires Catlab.
"""
function run_moments end

"""
    migrate!(s, rng; rate = 0.01)

Move individuals between the demes of a metapopulation `s`: each individual starts a
swap with probability `rate`, exchanging places with a random individual of a random
other deme, so the deme sizes are unchanged. Returns `(; count)`, the number of
exchanges. Requires Catlab.

`rate` is not the island model's migration rate m: an individual also moves when another
deme's swap picks it, so with equal deme sizes about `2·rate` of each deme is replaced
per generation (each exchange moves two individuals). Use `m/2` for an island-model
rate m. The SLiM lowering (`migRate`) does the same.

When `rate` is a matrix `M` ([`island_migration`](@ref), [`stepping_stone`](@ref)), an
individual of the deme in position `a` starts an exchange with probability `M[a, b]` with
a random individual of the deme in position `b`: migration along a graph. One uniform
draw decides whether and where, so a matrix draws differently from a number, even
`island_migration(k, r)`, which is the same model as `rate = r`. Demes must have the same
loci and ploidy, and may differ in size.
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
[`lift`](@ref), and recompute genotypes and phenotypes. `pop` must have `PopSize`
individuals (at most `PopSize`, which becomes the population size `N`, when the population
is density-regulated) with `NoQTL × Ploidy` alleles each, taken in the order of their parts (locus
fastest, then chromosome), as `lift` adds them; otherwise an `ArgumentError`. Returns
`w`. Requires Catlab.
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
                     lift, lower!, genotypes, phenotypes, genome_from_acset, run_moments)

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if any(f -> exc.f === f, CATEGORICAL) &&
           Base.get_extension(@__MODULE__, :QTLSimCatlabExt) === nothing
            print(io, "\n`", nameof(exc.f), "` is part of QTLSim's categorical layer, ",
                  "which is loaded with Catlab: run `import Catlab` first.")
        end
    end
end
