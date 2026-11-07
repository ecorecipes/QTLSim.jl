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
step that exchanges offspring among `k` demes, `Zygotes^⊗k → Zygotes^⊗k`, and
`migrate(k; rate)` the same step at its own rate, a number or a migration matrix.
Composites type-check: a step can follow another only where the stages match. Catlab's
`⊗` places life cycles side by side, one per deme: `generation ⊗ generation ⋅ migrate(2)`
is a generation of two demes coupled by migration. [`withparams`](@ref) gives steps
parameters of their own, and [`lower`](@ref) gives a life cycle its meaning in a backend.
Requires Catlab.
"""
function lifecycle end

"""
    lower(f, target; images = Dict(), record = nothing, migration_rate = 0.0)

The image of the life-cycle morphism `f` under the lowering functor to `target`:
`:kernel` gives an operator on a [`World`](@ref), `:eidos` a fragment of Eidos code for
SLiM, in which `@POP@` stands for the deme's subpopulation, and `:moments` an update of
the moments of a Gaussian population, the predictions of theory (see
[`run_moments`](@ref)). A morphism between products of stages, one per deme, has an
image on several demes: an operator on demes of a metapopulation (`DemeOp`, in the kernel
and the moments lowering), or Eidos for several subpopulations (`EidosMeta`). A functor
out of a free category is fixed by the images of the generators, so the image of a
composite or a product is the composite or product of the images. `record` steps in the
kernel append the parent generation's statistics to `record`, a named tuple of vectors
like the one [`run!`](@ref) returns, or discard them when it is `nothing`; `migrate` steps
exchange individuals among demes at `migration_rate` (in Eidos, the constant `migRate`), or
along a migration graph when it is a matrix ([`island_migration`](@ref),
[`stepping_stone`](@ref); in Eidos, the constant `MIG`), unless the step has a rate of its
own (`lifecycle().migrate(k; rate)`).

A step with parameters ([`withparams`](@ref)) is a generator of its own. Its image is the
plain step's, with the step's parameters in place of the deme's `Params`; in Eidos they
are literals in the step's code. Steps that are not in [`lifecycle`](@ref), made with
Catlab's `Hom(name, dom, codom)`, need an image in `images`, keyed by name: an operator or
a function `(state, rng) -> state` for `:kernel` and `:moments`, a string of Eidos for
`:eidos`. A key that is the name of a step with parameters, `first(Catlab.args(s))` for
the generator `s`, gives that step its own image. Eidos steps keep a deme's state on its
subpopulation: `@POP@.getValue("THETA")` and `@POP@.getValue("OMEGA")` are its optimum
and fitness width, set with `setValue`, and the steps before `select` leave `adults`,
`genotype` and `pheno` defined. Requires Catlab.
"""
function lower end

"""
    slim_script(f, p::Params; images = Dict(), migration_rate = 0.0)
    slim_script(f, ps::Vector{Params}; images = Dict(), migration_rate = 0.0)

A complete SLiM model, as Eidos source, that iterates the life cycle `f`, a morphism
`Zygotes → Zygotes`, or `Zygotes^⊗k → Zygotes^⊗k` for `k` demes, each a subpopulation
of `p.PopSize` individuals: the image of `f` under [`lower`](@ref) to Eidos runs once per
tick in a `reproduction()` callback, so `f` is one generation of each deme. The
parameters of `p` and `migration_rate` are the defaults, which `slim -d "NAME=value"`
overrides (`N`, `VP`, `h2`, `m`, `rate`, `L`, `RECOMB`, `NGEN`, `sampleEvery`, `moveGen`,
`moveSD`, `moveType`, `changeGen`, `changeSD`, `initOmega`, `drift`, `standardise`,
`demography`, `growthRate`, `migRate` or `MIG`, and with a fluctuating move
`moveNoiseSD`, `moveAutocorr`, `moveAmplitude` and `movePeriod`). The Eidos images of
the steps cover diploids (`Ploidy = 2`) with any number of loci, recombining as QTSim
does (`CrossoverModel = :qtsim`); with several demes, each row of output ends with the
deme. The life cycle must return each deme to its place: braidings must cancel.

With `ps`, one `Params` per deme, the demes have parameters of their own, as
`ps ⊗ demes(k)` gives them in the kernel. A constant whose value differs between demes
is a vector, one value per deme, which each deme keeps on its subpopulation, and the
deme's code reads `@POP@.getValue("NAME")`; so do the derived constants (`VE`, `ESD`,
`MUTVAR`). Every field of the list above can differ, except `NoQTL`, `NoGenerations`,
`Ploidy`, `CrossoverModel` and `Demography`, which must be the same in every deme. With
equal `ps` the model is the one of a single `Params`. Steps with parameters
([`withparams`](@ref)) override the deme's values in their own code. Requires Catlab.
"""
function slim_script end

"""
    withparams(f; kw...)

The life-cycle morphism `f` with parameters given to its steps. Each generator of `f`
that takes one of the parameters in `kw` becomes a new generator of the same type, named
by the step and its parameters. Its image in each lowering (`:kernel`, `:eidos`,
`:moments`) is the image of the plain step with those parameters in place of the deme's
`Params`, for that step only. The parameters each step takes:

- `select`: `ω`, a fixed width; `θ`, a fixed optimum, in place of the move; the move,
  `MoveType`, `MoveFitnessSD`, `MoveFitnessSDGen`, `MoveNoiseSD`, `MoveAutocorr`,
  `MoveAmplitude` and `MovePeriod`; `StandardisePhenotype` and `Drift`;
- `reproduce`: `RecombPoissonMean`, `CrossoverModel` (`:qtsim` only in Eidos),
  `Demography` and `GrowthRate`;
- `record`: `Sample`;
- `rescale`: `ChangeFitnessGen` and `ChangeFitnessSD`;
- `mutate`: `MutationScalar` and `MutationRatePerGamete`;
- `migrate`: `rate`, a number or a migration matrix, in place of the lowering's
  `migration_rate`. `lifecycle().migrate(k; rate)` is the same.

So `withparams(L.select; ω = 0.3)` is a selection step of width 0.3, and
`withparams(g; θ = 0.1) ⊗ g` a generation of two demes whose first selects toward 0.1.
`L.develop ⋅ withparams(L.select; ω = 0.3) ⋅ withparams(L.select; ω = 1.0) ⋅ …` is a life
cycle with two selection steps of different widths. A parameter that no step of `f`
takes is an `ArgumentError`, as is a value that `Params` would reject. Giving parameters
to a step that already has some merges them, the new ones winning. Two generators are
equal when their steps and parameters are.

A step's parameters leave the deme's state alone. Its optimum `θ` and width `ω` are
restored after a step with its own, so the `theta` and `omega` columns of the records
are the deme's; but a move made with a step's own move parameters moves the deme's
optimum, as a move does. In Eidos the overridden constants become literals in that
step's code. Steps without parameters lower exactly as before. Requires Catlab.
"""
function withparams end

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
`extinct_gen` per deme. A life cycle of several generations, such as
`withparams(g; ω = 0.2) ⋅ withparams(g; ω = 1.0)` (seasons of different selection), runs
whole cycles until `NoGenerations` generations are done; each `swap` advances the
generation. Requires Catlab.
"""
function run_composed end

"""
    run_moments(p::Params; cycle = lifecycle().generation, images = Dict(), zbar = 0.0, VG = 0.0,
                optimum = nothing, erosion = false, migration_rate = 0.0)
    run_moments(ps::Vector{Params}; cycle, images = Dict(), zbar = 0.0, VG = 0.0, erosion = false,
                migration_rate = 0.0)

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
  drift); with `erosion = true`, drift erodes the offspring's genetic variance by the
  factor `1 − 1/(2 N_e)`, with `N_e` the expected number of survivors, which mate at
  random;
- `mutate`: adds `m·V_E`.

The expected population size is `PopSize`, or under `Demography = :beverton_holt` the
mean of [`noffspring`](@ref) for the expected survivors, `R n/(1 + (R − 1) n/K)`.

Returns the columns of [`run!`](@ref), predicted (`nsurv` rounded). The genetic variance
converges to [`gaussian_vg_exact`](@ref); under a continuous move the lag `θ − z̄`
converges to [`steady_lag`](@ref). The population is extinct when fewer than two
individuals (one, if haploid) are expected to survive; [`extinction_distribution`](@ref)
adds the demographic stochasticity along the trajectory. Requires Catlab.

The fluctuating moves lower to their expected optimum ([`fluctuating_optimum`](@ref)
without noise): `MoveType = :periodic` exactly, since it is deterministic, and
`:red_noise` as its linear trend. To follow one realisation of the noise, pass the
optimum of every generation as `optimum` (`optimum[t + 1]` for generation `t`, at least
`NoGenerations` values; for example `theta` from a kernel run with `Sample = 1`). It
replaces the move whatever `MoveType` is, and the prediction is then the Gaussian
population's response to that environment ([`expected_lag2`](@ref) gives its average).
One deme only.

**Several demes.** For a life cycle on `k` demes, such as `g ⊗ g ⋅ lifecycle().migrate(2)`,
`p` is every deme's `Params`; `ps`, a vector, gives each deme its own, and `zbar` and `VG`
may be vectors, one value per deme. The statistics have a row per deme and generation
(`deme`), and `extinct_gen` is per deme, as from [`run_composed`](@ref). The image of
`migrate` is the kernel's migration ([`migrate!`](@ref)) at the level of first and second
moments. With `M[a, b]` the probability that a visit to an individual of deme `a` starts
an exchange with deme `b` (`rate/(k − 1)` for a number) and `N` the deme sizes, the
expected fraction `P[a, c]` of deme `a` that comes from deme `c` follows from visiting the
demes in order: while deme `a` is visited, deme `b` moves toward `a`'s composition by a
factor `f = (1 − M[a, b]/N_b)^{N_a}`, and `a` receives `N_b (1 − f)` times the difference
on top of what it sends. A deme then has mean `Σ_c P[a, c] z̄_c` and genetic variance
`Σ_c P[a, c] V_c + Σ_c P[a, c] (z̄_c − z̄_a)²`, within-deme plus between-deme variance, and
migrants keep their environmental variance. The means are exact in expectation, and so are
the variances right after migration, up to the sampling of the numbers of migrants. The
lowering carries the between-deme part as genic variance afterwards, as it does the
reduction by selection. In the simulations, part of it is association within individuals,
which segregation and recombination break up (at one locus, half of it goes at the next
reproduction), so the kernel's variance falls below the moments' after mixing.
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

The swaps are made one after the other, so an individual that has just arrived in a deme
can move again when that deme takes its turn. The order of the turns therefore matters:
with the demes always visited in the order 1, …, k, migration would not commute with
relabelling the demes, not even for the island model (k = 2 demes of 2 individuals are a
counterexample). Each call therefore visits the demes in a uniformly random order (a
`randperm(rng, k)`), and the individuals of each deme in order. Then migration commutes
with every relabelling that preserves the matrix, in distribution: every relabelling for
the island model (`island_migrate_after_generation_equivariant` in
`lean/CategoricalQuantGen/Migration.lean`), the rotations and reflections of the ring for
[`stepping_stone`](@ref). Extinct demes are skipped. The SLiM lowering draws its order
the same way.
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

const CATEGORICAL = (lifecycle, withparams, lower, slim_script, operators, operator, demes,
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
