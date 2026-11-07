# The categorical layer
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [The population schema](#the-population-schema)
- [From arrays to an ACSet and back](#from-arrays-to-an-acset-and-back)
- [A generation as a composite of
  operators](#a-generation-as-a-composite-of-operators)
- [Stratifying with ⊗](#stratifying-with-)

## Overview

QTSim runs one generation as a single loop over the population. This
package splits that loop into separate operators and keeps two
descriptions of the population:

- an **abstract** one, in which the population is an attributed C-set
  (ACSet), the categorical data structure of AlgebraicJulia ([Patterson
  et al. 2022](#ref-Patterson2022)), and each evolutionary process is an
  `Operator` that can be composed with `⋅` and stratified with `⊗`;
- a **lowered** one, the array kernel that does the numerical work and
  reproduces QTSim’s output.

This vignette shows the population as an ACSet, moves it between the two
levels, and builds a generation out of operators. The facts that make
the decomposition legitimate, that composition is associative and that a
run is the iterate of one composed generation, are proved in Lean in
`lean/CategoricalQuantGen/Composition.lean`.

## Setup

The categorical layer is a package extension of QTLSim, which Julia
loads as soon as Catlab is loaded. Catlab exports names of its own, such
as `swap` (the symmetric braiding), that would clash with the package’s
operators, so we load it with `import` and bring in only the functions
we need. `operators()` returns the operators of one generation.

``` julia
using QTLSim
import Catlab                       # loads QTLSim's categorical layer
using Catlab: ⋅, ⊗, nparts, incident, subpart, set_subpart!, generators
(; develop, select, reproduce, mutate, swap, generation) = operators()
using Printf, Random
```

## The population schema

The schema, `population_schema()`, has two kinds of part, individuals
and alleles, and a map `owner` sending each allele to the individual
that carries it. Allele effects and the individuals’ environmental
deviation, genotype, phenotype and fitness are attributes.

    Allele ──owner──▶ Individual        effect : Allele → ℝ
                                        env, genotype, phenotype, fitness : Individual → ℝ

``` julia
(objects = nameof.(generators(population_schema(), :Ob)),
 maps = nameof.(generators(population_schema(), :Hom)),
 attributes = nameof.(generators(population_schema(), :Attr)))
```

    (objects = [:Individual, :Allele], maps = [:owner], attributes = [:effect, :env, :genotype, :phenotype, :fitness])

## From arrays to an ACSet and back

The kernel keeps the population in a `World`: genomes as blocks of
alleles shared between individuals, which `genomes` and `setgenomes!`
read and write as a `(NoQTL, Ploidy, PopSize)` array, plus
per-individual vectors. `lift` turns a `World` into a `Population`
ACSet. For a population small enough to print, with effects rounded so
they are easy to read:

``` julia
p = Params(PopSize = 3, Ploidy = 2, NoQTL = 2, PhenotypeVar = 0.02, Heritability = 0.4)
w = World(p)
rng = Xoshiro(1)
setgenomes!(w, round.(0.1 .* randn(rng, p.NoQTL, p.Ploidy, p.PopSize), digits = 3))
w.env .= round.(sqrt(envVar(p)) .* randn(rng, p.PopSize), digits = 3)
develop!(w)                    # genotype = sum of effects; phenotype = genotype + env
pop = lift(w)
```

```@raw html
<div class="c-set">
```
QTLSimCatlabExt.PopulationACSet{Float64, Float64} {Individual:3, Allele:12, Effect:0, Trait:0}

| Individual |    env | genotype | phenotype | fitness |
|-----------:|-------:|---------:|----------:|--------:|
|          1 | -0.043 |   -0.021 |    -0.064 |     0.0 |
|          2 |  0.064 |    0.056 |      0.12 |     0.0 |
|          3 | -0.003 |    0.241 |     0.238 |     0.0 |

| Allele | owner | effect |
|-------:|------:|-------:|
|      1 |     1 |  0.006 |
|      2 |     1 |  0.028 |
|      3 |     1 |  -0.06 |
|      4 |     1 |  0.005 |
|      5 |     2 |  0.109 |
|      6 |     2 | -0.158 |
|      7 |     2 |  0.018 |
|      8 |     2 |  0.087 |
|      9 |     3 |  0.097 |
|     10 |     3 |  0.155 |
|     11 |     3 | -0.058 |
|     12 |     3 |  0.047 |

```@raw html
</div>
```

Each individual owns `Ploidy × NoQTL` = 4 alleles. Its genotype is the
sum of the effects of the alleles in the preimage of `owner`, which
`incident` returns. `genotypes` computes all of them this way, at the
abstract level, and they agree with the kernel’s.

``` julia
alleles = incident(pop, 1, :owner)
(alleles = alleles, effects = subpart(pop, alleles, :effect),
 genotype_of_individual_1 = sum(subpart(pop, alleles, :effect)),
 all_genotypes_agree = genotypes(pop) ≈ w.genotype)
```

    (alleles = [1, 2, 3, 4], effects = [0.006, 0.028, -0.06, 0.005], genotype_of_individual_1 = -0.020999999999999994, all_genotypes_agree = true)

The translation goes both ways. A change made to the ACSet, here giving
one allele of individual 1 a large effect, reaches the kernel through
`lower!`:

``` julia
set_subpart!(pop, alleles[1], :effect, 1.0)
w2 = World(p)
lower!(pop, w2)
develop!(w2)
(genotype_before = w.genotype[1], genotype_after = w2.genotype[1],
 others_unchanged = w2.genotype[2:3] == w.genotype[2:3])
```

    (genotype_before = -0.020999999999999998, genotype_after = 0.973, others_unchanged = true)

The round trip is exact for a full-sized population too:

``` julia
big = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 0.02, Heritability = 0.4)
wb = World(big)
setgenomes!(wb, 0.05 .* randn(Xoshiro(2), big.NoQTL, big.Ploidy, big.PopSize))
develop!(wb)
t = @elapsed popb = lift(wb)
wb2 = World(big); lower!(popb, wb2)
@printf("%d alleles; lift took %.2f s; genome recovered exactly: %s\n",
        nparts(popb, :Allele), t, genomes(wb2) == genomes(wb))
```

    200000 alleles; lift took 0.05 s; genome recovered exactly: true

## A generation as a composite of operators

Each step of the life cycle is an operator: `develop` computes genotypes
and phenotypes, `select` moves the optimum if it is due to move and
draws the survivors, `reproduce` makes the offspring, `record` notes the
parents’ statistics (here it discards them), `rescale` changes the width
of the fitness function if a change is due, `mutate` adds mutations to
the offspring and `swap` makes the offspring the next parents. One
generation is their composite, written with Catlab’s diagrammatic
composition `⋅`, which reads left to right:

``` julia
generation
```

    Operator(:develop⋅select⋅reproduce⋅record⋅rescale⋅mutate⋅swap)

`(f ⋅ g)(w, rng)` applies `f` and then `g`, skipping the rest of a
composite once the population is extinct. Because composition is
associative, how a composite is bracketed does not matter. Lean proves
this for composition in general; here it is for one case, run for 50
generations from the same state and seed:

``` julia
q = Params(PopSize = 200, Ploidy = 2, NoQTL = 20, PhenotypeVar = 0.02, Heritability = 0.4,
           MutationScalar = 0.004, MutationRatePerGamete = 0.01, RecombPoissonMean = 3.0)

function iterate_op(op, q; gens = 50, seed = 1)
    w = World(q); rng = Xoshiro(seed)
    w.env .= sqrt(envVar(q)) .* randn(rng, q.PopSize)
    for _ in 1:gens; op(w, rng); end
    w
end

left = iterate_op((develop ⋅ select) ⋅ (reproduce ⋅ mutate ⋅ swap), q)
right = iterate_op(develop ⋅ (select ⋅ (reproduce ⋅ (mutate ⋅ swap))), q)
genomes(left) == genomes(right)
```

    true

The second fact Lean proves (`loop_eq_iterate`) is that running the
generation loop $n$ times is the $n$-fold iterate of one generation.
`run_composed` drives a run by applying the composite `generation`,
recording statistics as it goes; `run!` is the hand-written loop of the
kernel, in QTSim’s order. Under the same seed they produce the same run,
bit for bit:

``` julia
qr = Params(q; NoGenerations = 300, MoveFitnessSDGen = 150, MoveFitnessSD = 0.01,
            MoveType = :continuous)
a = run!(World(qr), Xoshiro(7))
b = run_composed(World(qr), Xoshiro(7))
(genetic_variance = a.geno_var == b.geno_var, mean_phenotype = a.pheno_mean == b.pheno_mean,
 survivors = a.nsurv == b.nsurv, optimum = a.theta == b.theta)
```

    (genetic_variance = true, mean_phenotype = true, survivors = true, optimum = true)

## Stratifying with ⊗

The monoidal product `⊗` turns one population’s parameters into a
metapopulation of independent demes. An operator built for one deme acts
on every deme of a `Stratified` state, and a `migration` operator,
composed in with `⋅`, couples them.

``` julia
meta = p ⊗ demes(4)
rng = Xoshiro(3)
for d in meta.demes; d.env .= sqrt(envVar(p)) .* randn(rng, p.PopSize); end
step = generation ⋅ migration(0.25)
step(meta, rng)
(demes = ndemes(meta), total = totalpop(meta), generation_of_each_deme = [d.gen for d in meta.demes])
```

    (demes = 4, total = 12, generation_of_each_deme = [1, 1, 1, 1])

[Recomposing the
model](05_recomposing_the_model.md) uses
these pieces to change the model’s biology without touching the kernel.
The operators are themselves the kernel’s image of a life cycle written
once as syntax: [One life cycle, two
backends](06_one_life_cycle_two_backends.md)
derives them, and a SLiM model, from it.

```@raw html
<div id="refs" class="references csl-bib-body hanging-indent">
```

```@raw html
<div id="ref-Patterson2022" class="csl-entry">
```

Patterson, Evan, Owen Lynch, and James Fairbanks. 2022. “Categorical
Data Structures for Technical Computing.” *Compositionality* 4: 5.
<https://doi.org/10.32408/compositionality-4-5>.

```@raw html
</div>
```

```@raw html
</div>
```
