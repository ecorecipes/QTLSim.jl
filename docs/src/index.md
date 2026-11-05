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

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/QTLSim.jl")
```

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

`Params` takes QTSim's options under the same names, and `simulate` returns the mean
and variance of the phenotype, the genetic variance, the number of survivors and the
optimum for every sampled generation.

## The categorical layer

Loading Catlab adds a categorical description of the same model, as the package
extension `QTLSimCatlabExt`. The life cycle is written once, as a morphism of a free
category whose objects are the stages of a generation and whose generators are its
steps, and functors lower it to the simulation kernel and to a model for the
population-genetics simulator SLiM, in its language Eidos:

```julia
import Catlab
using Catlab: ⋅, ⊗
using Random
L = lifecycle()
qtsim = L.develop ⋅ L.select ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅ L.mutate ⋅ L.swap
# develop⋅select⋅reproduce⋅record⋅rescale⋅mutate⋅swap: Zygotes → Zygotes

run_composed(World(p), Xoshiro(1); cycle = qtsim)   # the same run as run!, bit for bit
print(slim_script(qtsim, Params(p; NoQTL = 1)))     # the same life cycle as a SLiM model
```

The SLiM model generated this way agrees with QTSim in distribution. In the kernel the
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
- **Against the published results**: on 589 of the paper's input files, the Julia
  model's genetic variance before the move and its time to extinction are within about
  1% of the published values.

The [validation tutorial](tutorials/08_validation_against_qtsim.md) shows all three,
and the [replication](tutorials/09_replicating_pagel_2025.md) reproduces the paper's
results. QTSim's C source comes with the package, unmodified, in `reference/QTSim`, and
the package's continuous integration builds it and runs the exact and statistical checks
on every change.

## Tutorials

1. [Getting started](tutorials/01_getting_started.md)
2. [Mutation–selection balance](tutorials/02_mutation_selection_balance.md)
3. [Adaptation to a moving optimum](tutorials/03_moving_optimum.md)
4. [The categorical layer](tutorials/04_categorical_layer.md)
5. [Recomposing the model](tutorials/05_recomposing_the_model.md)
6. [One life cycle, two backends](tutorials/06_one_life_cycle_two_backends.md)
7. [Loci and recombination](tutorials/07_loci_and_recombination.md)
8. [Validation against QTSim](tutorials/08_validation_against_qtsim.md)
9. [Replicating Pagel, Gardner and Meade (2025)](tutorials/09_replicating_pagel_2025.md)

## Licence

QTLSim.jl is free software under the GNU General Public License, version 3 or (at your
option) any later version. It reimplements the algorithm of QTSim (copyright 2023 Andrew
Meade, University of Reading), which is distributed under the same licence, and includes
QTSim's source, unmodified, for testing and verification.
