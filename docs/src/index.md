# QTLSim.jl

QTLSim.jl is a Julia port of QTSim, the simulation behind Pagel, Gardner and Meade,
*Quantitative-genetic analysis of directional adaptation suggests low maximum
sustainable rates of change in agreement with data from field populations*
([Scientific Reports, 2025](https://doi.org/10.1038/s41598-025-24445-2)). QTSim is a C
program by Andrew Meade (University of Reading); its source calls it QTLSim 1.0.

The model follows a population of diploid individuals with one additive quantitative
trait, under Gaussian stabilising selection on viability toward an optimum that can
move. Survivors mate at random, gametes recombine at Poisson crossover points, and
offspring receive Poisson numbers of Gaussian mutations. The paper used it to estimate
how fast a population can track a steadily changing environment before it goes extinct:
about 3% of a phenotypic standard deviation per generation.

## Installation

QTLSim.jl is not registered in the Julia General registry; install it from GitHub
(Julia 1.10 or later):

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/QTLSim.jl")
```

The categorical layer also needs Catlab (`Pkg.add("Catlab")`), loaded as a package
extension.

## Quick start

A population at mutation–selection balance whose optimum then moves by a fixed amount
every generation:

```julia
using QTLSim

p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 12_000,
           PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
           MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0,
           MoveFitnessSDGen = 10_000, MoveFitnessSD = 0.05, MoveType = :continuous)
out = simulate(p; seed = 1)
out.extinct, out.extinct_gen
```

`Params` has QTSim's defaults and follows its options, with a few renamings (see
[Parameters](parameters.md)), and `simulate` returns the mean and variance of the
phenotype, the genetic variance, the number of survivors and the optimum for every
sampled generation ([Output](output.md)).

Beyond QTSim, and off by default, the kernel has Beverton–Holt density regulation
(`Demography = :beverton_holt`), fecundity selection (`Selection = :fecundity`),
periodic and red-noise optima, metapopulations with island and stepping-stone migration,
and several traits with pleiotropy (`NoTraits`; see [Several traits](traits_design.md)).
The package also has the theory the simulations are compared with: Gaussian, multilocus
and rare-alleles equilibria of the genetic variance, house of cards, the lag behind a
moving or fluctuating optimum, extinction and rescue, and neutral F_ST/Q_ST among demes
(see the [API reference](api.md)).

## The categorical layer

Loading Catlab adds a categorical description of the same model, as the package
extension `QTLSimCatlabExt`. The life cycle is written once, as a morphism of a free
symmetric monoidal category whose objects are the stages of a generation and whose
generators are its steps, and functors lower it to the simulation kernel and to a model
for the population-genetics simulator SLiM, in its language Eidos, and to the
deterministic moments of the population (`run_moments`):

```julia
import Catlab
using Catlab: ⋅, ⊗
using Random
L = lifecycle()
qtsim = L.develop ⋅ L.select ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅ L.mutate ⋅ L.swap
# develop⋅select⋅reproduce⋅record⋅rescale⋅mutate⋅swap: Zygotes → Zygotes

run_composed(World(p), Xoshiro(1); cycle = qtsim)   # the same run as run!, bit for bit
print(slim_script(qtsim, p))                        # the same life cycle as a SLiM model

meta = qtsim ⊗ qtsim ⊗ qtsim ⋅ L.migrate(3)         # three demes coupled by migration
run_composed(p ⊗ demes(3), Xoshiro(1); cycle = meta, migration_rate = 0.01)
print(slim_script(meta, p; migration_rate = 0.01))  # the same metapopulation in SLiM
```

The SLiM model generated this way agrees with QTSim in distribution, with one locus or
many, and on a model QTSim does not have, demes coupled by migration, the kernel and the
SLiM model agree with each other. In the kernel the
steps are operators, which can be rearranged, extended or stratified into demes, and the
population is also an attributed C-set (ACSet):

```julia
(; develop, select, reproduce, mutate, swap, generation) = operators()
double = develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ mutate ⋅ swap   # twice the mutation
meta = p ⊗ demes(4)                                               # four demes
(generation ⋅ migration(0.01))(meta, Xoshiro(1))                  # with migration
```

## Validation

The kernel reproduces QTSim at three strengths:

- **Statistically**: run many times, the Julia model and QTSim give the same
  distributions of genetic variance, survival and adaptation.
- **Exactly**: a replica that draws its random numbers from QTSim's GSL routines in
  QTSim's order matches QTSim's output generation by generation.
- **Against the published results**: on about 590 of the paper's input files, in two
  runs with independent seeds, the Julia model's genetic variance before the move and its
  time to extinction are within about 1% of the published values (95% intervals inside
  0.98–1.03).

The SLiM model generated from the life cycle agrees with QTSim in distribution (a
pre-specified test of 40,000 runs per program found mean genetic variances within 0.6%),
and with the kernel on the models QTSim does not have.

The [validation tutorial](tutorials/08_validation_against_qtsim.md) shows all three,
and the [replication](tutorials/09_replicating_pagel_2025.md) reproduces the paper's
results. QTSim's C source comes with the package, unmodified, in `reference/QTSim`, and
the package's continuous integration builds it and runs the exact and statistical checks
on every change.

## Tutorials

1. [Getting started](tutorials/01_getting_started.md)
2. [Mutation–selection balance](tutorials/02_mutation_selection_balance.md)
3. [Adaptation to a moving optimum](tutorials/03_moving_optimum.md)
4. [Loci and recombination](tutorials/04_loci_and_recombination.md)
5. [The categorical layer](tutorials/05_categorical_layer.md)
6. [Recomposing the model](tutorials/06_recomposing_the_model.md)
7. [One life cycle, two backends](tutorials/07_one_life_cycle_two_backends.md)
8. [Validation against QTSim](tutorials/08_validation_against_qtsim.md)
9. [Replicating Pagel, Gardner and Meade (2025)](tutorials/09_replicating_pagel_2025.md)
10. [Demography and evolutionary rescue](tutorials/10_demography_and_rescue.md)
11. [Density regulation in the paper's populations](tutorials/11_demography_in_the_papers_populations.md)
12. [Metapopulations and migration](tutorials/12_metapopulations.md)
13. [Several traits and pleiotropy](tutorials/13_multiple_traits.md)

## Proofs in Lean

The facts that justify the decomposition and the lowerings are proved in Lean 4 with
Mathlib, in the library `CategoricalQuantGen` of the project this package comes from (not
part of the package's repository); docstrings cite it as `lean/CategoricalQuantGen/…`.

## Licence

QTLSim.jl is free software under the GNU General Public License, version 3 or (at your
option) any later version. It reimplements the algorithm of QTSim (copyright 2023 Andrew
Meade, University of Reading), which is distributed under the same licence, and includes
QTSim's source, unmodified, for testing and verification.

## Citing

Please cite the paper whose model QTLSim.jl implements, Pagel, Gardner and Meade,
*Scientific Reports* 15 (2025), [doi:10.1038/s41598-025-24445-2](https://doi.org/10.1038/s41598-025-24445-2),
and the software, as its
[`CITATION.cff`](https://github.com/ecorecipes/QTLSim.jl/blob/main/CITATION.cff) describes.
