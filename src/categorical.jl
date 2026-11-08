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
    lower(f, target; images = Dict(), record = nothing, migration_rate = 0.0, model = OneTraitModel())

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

`model` ([`QTLSim.Model`](@ref)) chooses the functor: each model has its own table of images
of the steps in each target, so the same `f` lowers to QTSim's one-trait operators or, with
`model = TraitModel(n)`, to those of `n` traits; an image refuses a population of the other
model. [`run_composed`](@ref), [`run_moments`](@ref) and [`slim_script`](@ref) pass the model
of their `Params`.

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
    withparams(f, model; kw...)
    withparams(f, p::Params; kw...)

The life-cycle morphism `f` with parameters given to its steps. Each generator of `f`
that takes one of the parameters in `kw` becomes a new generator of the same type, named
by the step and its parameters. Its image in each lowering (`:kernel`, `:eidos`,
`:moments`) is the image of the plain step with those parameters in place of the deme's
`Params`, for that step only. The parameters belong to a model ([`QTLSim.Model`](@ref)):
they are checked against the schema of `model` ([`step_schema`](@ref)), of `model_of(p)`
with a `Params`, or of QTSim's [`OneTraitModel`](@ref) without either, and again against the
model a lowering targets, so a step written for one model is rejected by another whose
schema it does not fit. The parameters each step takes in the one-trait model:

- `select`: `ω`, a fixed width; `θ`, a fixed optimum, in place of the move; the move,
  `MoveType`, `MoveFitnessSD`, `MoveFitnessSDGen`, `MoveNoiseSD`, `MoveAutocorr`,
  `MoveAmplitude` and `MovePeriod`; `StandardisePhenotype` and `Drift`; `Selection`
  (under `:fecundity` the survivors, all of them, carry their fitness as a mating weight,
  which the plain `reproduce` samples by);
- `reproduce`: `RecombPoissonMean`, `CrossoverModel` (`:qtsim` only in Eidos),
  `Demography` and `GrowthRate`;
- `record`: `Sample`;
- `rescale`: `ChangeFitnessGen` and `ChangeFitnessSD`;
- `mutate`: `MutationScalar` and `MutationRatePerGamete`;
- `migrate`: `rate`, a number or a migration matrix, in place of the lowering's
  `migration_rate`. `lifecycle().migrate(k; rate)` is the same.

In the model of `n` traits ([`TraitModel`](@ref)`(n)`), `select`'s `θ` is a vector of `n`, and
`select` also takes `Ω`, a fitness covariance matrix in place of `ω² S` (so not with `ω` or
`SelCorr`), `SelCorr`, `MoveDirection`, `MoveNoiseCov` and `FluctDirection`, but not
`StandardisePhenotype`; `mutate` also takes `MutCorr` and `MutCov`. `withparams(g, p; θ = [0.1, 0.0]) ⊗ withparams(g, p; θ = [0.0,
0.1])` is two demes with optima along different traits.

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
    operators(model = OneTraitModel())

The kernel images of the steps of the life cycle ([`lifecycle`](@ref)) in `model` (a
[`QTLSim.Model`](@ref) or a `Params`, whose model is `model_of(p)`), as the named
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
                optimum = nothing, erosion = false, migration_rate = 0.0, between = false, origins = :exchange,
                loci = :auto)
    run_moments(ps::Vector{Params}; cycle, images = Dict(), zbar = 0.0, VG = 0.0, erosion = false,
                migration_rate = 0.0, between = false, origins = :exchange, loci = :auto)

The life cycle `cycle` lowered to the moments of a Gaussian population (`lower(cycle, :moments)`)
and iterated for `NoGenerations`: the predictions of quantitative-genetic theory for the
same morphism the kernel and SLiM simulate. Each step's image updates the means and
variances:
- `develop`: `V_P = V_G + V_E`;
- `select`: moves the optimum as `MoveType` says, then applies Gaussian selection of
  width ω to `z ~ N(z̄, V_P)`. This gives the mean fitness, the expected survivors
  `N·W̄`, the selection differential, and the survivors' genotypic mean and variance
  `z̄ + V_G λ/s` and `V_G − V_G²/s`, with `λ = θ − z̄` and `s = ω² + V_P`;
- `reproduce`: the offspring inherit the survivors' mean; their variance is what
  recombination makes of the survivors' (below);
- `mutate`: adds `V_m = m·V_E`.

**Which model (`loci`).** The default, `loci = :auto`, is the rare-alleles closure
(`:rare_alleles`, below) for one diploid trait, and the multilocus Gaussian model
(`:multilocus`) for several traits and for haploids (the closure has not been compared
with the kernel for haploids), and when the `mutate` steps of `cycle` carry different
`MutationScalar` or `MutationRatePerGamete` (no one mutation effect variance);
[`QTLSim.moments_loci`](@ref) gives the resolution for a `Params`. Ask for
`loci = :multilocus`, `:rare_alleles` or `:aggregate` to use those models anywhere (the
closure for several traits and for haploid clones is in src/rare_alleles_traits.jl). The rare-alleles closure is within 15% of the kernel in 46 of
47 weak-selection cells of `validate/moments_regimes.jl`, never farther from it than the
multilocus Gaussian, and up to 10× closer in the house-of-cards and drift regimes; when
mutations are many and small the two agree.

**The multilocus model (`loci = :multilocus`).** The population is
`NoQTL` additive loci whose allelic values are jointly Gaussian (Lande 1975; Bulmer
1980; Turelli and Barton 1990). The lowering carries the expected second moments of the
haplotypes (one gamete's alleles summed over the loci): the variance `T` of one
haplotype, the covariance `R` of an adult's two and the covariance `Q` of haplotypes of
two distinct adults, each in full and restricted to alleles at one locus (the genic
parts; the difference is linkage disequilibrium, LD). Selection lowers a haplotype's
variance, and the covariance of an adult's two, by `D = (V_G/2)²/s`. Each locus carries
`1/NoQTL` of the variance, so the genic parts fall by `D/NoQTL`, and the rest is
negative LD. At reproduction an offspring's two haplotypes come from distinct parents
(`R ← Q`). Recombination turns LD back into variance, the LD of each class of pairs of
loci at its own rate (`QTLSim.pair_classes`, from the kernel's crossover; their mean is
[`pair_recombination`](@ref)). Mutation adds `V_m/2` to each haplotype. A haploid's
offspring is a clone of one parent. The equations are in src/multilocus.jl
(`QTLSim.hap_select`, `hap_reproduce`, `hap_mutate`). Without drift (the default) the
population is infinite, and one deme converges to [`gaussian_vg_loci`](@ref),
`V²/(V_s + V) = 2L·V_m` for a diploid of `L` loci that recombine. With `erosion = true`
the deme has its expected `N` adults and `n` parents (the expected survivors, or under
fecundity selection `N W̄²/E[W²]`), and drift enters through `Q`: two offspring share a
parent with probability `1/n`. Then `geno_var` is the expected within-population sample
variance (divisor `N − 1`, as the kernel records it), converging to
`gaussian_vg_loci(p; drift = true)`. Under a continuous move the lag `θ − z̄` converges
to [`steady_lag`](@ref) with that variance. The model is good when many small mutations
hold each locus's variance; when few large ones do (the house-of-cards regime,
[`house_of_cards_vg`](@ref)) the kernel holds less (`validate/moments_loci.jl`).

**The rare-alleles closure (`loci = :rare_alleles`; the default for one diploid
trait).** The multilocus model
with the fourth cumulant of a locus's alleles set to `α²` times their variance (`α²` =
`QTLSim.mutation_effect_variance(p)`, the variance of one mutation's effect): selection
also removes `β²α²/(2s)` of the genic variance, each rare mutant at its own rate, with
`β = H/H_genic` the regression of the genotype on a locus's alleles (`β = 1` without
LD). It equals the multilocus Gaussian when mutations are many and small, and reduces to
the (stochastic) house of cards when they are rare and large and to `2N V_m` when
selection is weak against drift; its equilibrium is [`rare_alleles_vg`](@ref), and for
loci in linkage equilibrium [`house_of_gauss_vg`](@ref). Equations in
src/rare_alleles.jl; the kernel comparison across regimes, including the paper's
populations, is `validate/moments_regimes.jl`. With drift (`erosion = true`) one deme
converges to `rare_alleles_vg(p; drift = true)`, and under a continuous move the lag to
[`steady_lag`](@ref) with that variance. Each deme of a metapopulation has its own `α²`;
when the `mutate` steps of `cycle` carry `MutationScalar` or `MutationRatePerGamete`
(all the same), `α²` is computed with them.

With several traits (`loci = :rare_alleles`) a locus's genic covariance `c` also loses
`½(𝓡 c + c 𝓡ᵀ)`, with `𝓡 = ½ C X⁻¹`, `C` the covariance of one mutation's effect vector and
`X = E[a aᵀ/(aᵀ s̃⁻¹ a)]` (`s̃⁻¹ = βᵀ s⁻¹ β`): exact at the multivariate house-of-cards
balance, `α²/(2s)` for one trait ([`QTLSim.rare_alleles_rate`](@ref)). Its equilibrium is
[`rare_alleles_G`](@ref), and for loci in linkage equilibrium [`house_of_gauss_G`](@ref), the
stabilising solution of an algebraic Riccati equation. Haploid offspring are clones, so a
haploid genome is one unit of inheritance: the closure then acts on a haplotype's whole
variance (`β = 1`, [`QTLSim.hap_clonal`](@ref)), with `L μ` mutations per genome.

**The aggregate recursion (`loci = :aggregate`).** QTLSim's moments lowering before the
multilocus model: the offspring inherit the survivors' variance, `V ← V − V²/(V_s + V) +
V_m`, converging to [`gaussian_vg_exact`](@ref). It treats all of selection's reduction
as inherited, which is exact for one haploid locus (and for haploid clones), but not for
diploids: for one diploid locus half of it (between the haplotypes) is lost at the next
reproduction, and for many loci most of it is LD that recombination breaks up. With
`erosion = true` drift multiplies the offspring's variance by `1 − 1/(2n)`.

The expected population size is `PopSize`, or under `Demography = :beverton_holt` the
mean of [`noffspring`](@ref) for the expected survivors, `R n/(1 + (R − 1) n/K)`.

Returns the columns of [`run!`](@ref), predicted (`nsurv` rounded). The population is
extinct when fewer than two individuals (one, if haploid) are expected to survive;
[`extinction_distribution`](@ref) adds the demographic stochasticity along the
trajectory. Requires Catlab.

The fluctuating moves lower to their expected optimum ([`fluctuating_optimum`](@ref)
without noise): `MoveType = :periodic` exactly, since it is deterministic, and
`:red_noise` as its linear trend. To follow one realisation of the noise, pass the
optimum of every generation as `optimum` (`optimum[t + 1]` for generation `t`, at least
`NoGenerations` values; for example `theta` from a kernel run with `Sample = 1`). It
replaces the move whatever `MoveType` is, and the prediction is then the Gaussian
population's response to that environment ([`expected_lag2`](@ref) gives its average).
One deme only.

**Several demes.** For a life cycle on `k` demes, such as `g ⊗ g ⋅ lifecycle().migrate(2)`,
`p` is every deme's `Params`; `ps`, a vector, gives each deme its own (with one ploidy),
and `zbar` and `VG` may be vectors, one value per deme. The statistics have a row per
deme and generation (`deme`), and `extinct_gen` is per deme, as from
[`run_composed`](@ref). The image of `migrate` is the kernel's migration
([`migrate!`](@ref)) at the level of first and second moments. With `M[a, b]` the
probability that a visit to an individual of deme `a` starts an exchange with deme `b`
(`rate/(k − 1)` for a number) and `N` the deme sizes, the expected fraction `P[a, c]` of
deme `a` that comes from deme `c` follows from visiting the demes in order: while deme
`a` is visited, deme `b` moves toward `a`'s composition by a factor
`f = (1 − M[a, b]/N_b)^{N_a}`, and `a` receives `N_b (1 − f)` times the difference on top
of what it sends. A deme then has mean `Σ_c P[a, c] z̄_c` and genetic variance
`Σ_c P[a, c] V_c + Σ_c P[a, c] (z̄_c − z̄_a)²`, within-deme plus between-deme variance, and
migrants keep their environmental variance. The means are exact in expectation, and so
are the variances right after migration, up to the sampling of the numbers of migrants.
In the multilocus model the difference of the deme means enters the haplotype moments as
LD, spread over the loci, which segregation and recombination then break up (with
`loci = :aggregate` it stays as inherited variance, above the kernel's).

**Drift among demes (`between = true`).** Without drift, deme means that start equal stay
equal. With `between = true` drift is on (`erosion` is ignored) and the lowering also
carries the covariances `Q[a, b]` of haplotypes of different demes, about the demes'
deterministic means. It so predicts the within-deme genetic variance, the variance among
deme means and Q_ST for any life cycle on diploid demes, with drift, mutation, migration
and selection. Selection shrinks a deme's mean deviation by `1 − β`, `β = V_G/s`, as well
as acting within the deme as above. Migration, with backward matrix `P`
([`backward_migration`](@ref)) and the law `J` of the origins of two distinct adults,
gives `T′_a = Σ_c P_ac (T_c + Δ_ac²)`, likewise `R`, and
`Q′_ab = Σ_cd J[ab, cd] (Q_cd + Δ_ac Δ_bd)`, with `Δ_ac` half the difference between
deme `c`'s mean and deme `a`'s new mean (weighted by `1/NoQTL` in the genic parts).
`origins` is `:exchange` (the default; [`exchange_origins`](@ref), the kernel's
`migrate!`, which exchanges individuals), `:pool` (independent origins, `kron(P, P)`,
the textbook island model) or a `k² × k²` matrix, as for [`deme_coalescence`](@ref).
For a neutral trait (`Drift = true`) the composite is [`deme_coalescence`](@ref)'s
recursion read forward, so the predictions are [`neutral_deme_variances`](@ref) at every
generation, exactly: in the island model the within-deme variance settles at the
metapopulation's 2NV_m(kN − 1)/(N − 1), not at an isolated deme's 2NV_m. With selection
they are the Gaussian approximation; `validate/moments_demes.jl` compares them with the
kernel (k = 5 and 10, a common optimum and an optimum per deme), with `loci = :multilocus`;
in those scenarios (many small mutations) the rare-alleles default differs from it by less
than 1%. For one deme,
`between = true` is `erosion = true` with the extra column. The result has one more
column, `metapop`, a named tuple with one entry per census generation (every `Sample`th,
at the start of a generation, before any step): `gen`; `within`, the mean over demes of
the within-deme sample variance V_W; `var_among`, the expected sample variance (k − 1) of
the deme mean genotypes; `between`, `var_among` less the sampling variance of a deme
mean (σ²_B); `qst`, `σ²/(σ² + 2V_W)` with `σ² = (k − 1)/k·between` (Wright's convention);
`qst_demes`, the same with `σ² = between`; and `pairs`, the k × k matrices of half the
expected squared differences of deme means (the deterministic differences included).
These are the statistics of [`neutral_deme_variances`](@ref), with the same names. One
deme: `var_among`, `between` and the Q_STs are 0. Generation 0's Q_ST is 0/0, NaN.

**Several traits** (`NoTraits > 1`): the multilocus (the default) and aggregate models with matrices,
`run_trait_moments` (ext/QTLSimCatlabExt/trait_moments.jl); the fixed point is
[`gaussian_G_loci`](@ref) (with `loci = :aggregate`, [`gaussian_G`](@ref)). On a life cycle of
`k` demes too, with `migration_rate`, `between` and `origins` as above: `metapop` then holds
matrices (`within`, `var_among`, `between`), `qst` a vector per trait, and `pairs` a matrix of
matrices. `optimum` (one deme) is a vector of optimum vectors. A fluctuating move lowers to
its expected optimum, as for one trait.
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
