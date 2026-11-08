# QTLSim.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/QTLSim.jl/dev/)

QTLSim.jl is a Julia port of QTSim, the simulation behind Pagel, Gardner and Meade,
*Quantitative-genetic analysis of directional adaptation suggests low maximum
sustainable rates of change in agreement with data from field populations*
([Scientific Reports, 2025](https://doi.org/10.1038/s41598-025-24445-2)). QTSim is a C
program by Andrew Meade (University of Reading); its source calls it QTLSim 1.0.

The model follows a population of diploid individuals with an additive quantitative
trait under Gaussian stabilising selection on viability, toward an optimum that can
move. Survivors mate at random, gametes recombine at Poisson crossover points, and
offspring receive Poisson numbers of Gaussian mutations. The paper used it to estimate
how fast a population can track a steadily changing environment before it goes extinct:
about 3% of a phenotypic standard deviation per generation.

QTLSim reproduces QTSim (bit for bit, in a replica that draws QTSim's random numbers),
writes the same life cycle once as a categorical morphism that lowers to the Julia
kernel, to a SLiM model and to moment equations, and adds models QTSim does not have:
density regulation, fecundity selection, fluctuating optima, metapopulations with
migration, and several traits with pleiotropy.

## Installation

QTLSim.jl is not registered in the Julia General registry. Install it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/QTLSim.jl")
```

It needs Julia 1.10 or later. The categorical layer also needs Catlab
(`Pkg.add("Catlab")`), which it loads as a package extension.

## Quick start

A population at mutation–selection balance whose optimum then moves by a fixed amount
every generation (about 3 s):

```julia
using QTLSim

p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 12_000,
           PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
           MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0,
           MoveFitnessSDGen = 10_000, MoveFitnessSD = 0.05, MoveType = :continuous)
out = simulate(p; seed = 1)        # per-generation statistics, as QTSim prints them
out.extinct, out.extinct_gen
```

`Params` has QTSim's defaults, and its fields follow QTSim's options with a few renamings
(`Generations` is `NoGenerations`, `RecombinationPoissonMean` is `RecombPoissonMean`,
`OptimumMoveType` is `MoveType`; the two-value `MoveFitnessSD` and `ChangeFitnessSD`
options are each split into a generation and a value). `read_qtsim_options` reads a
QTSim options file into a `Params`.

## The categorical layer

Loading Catlab adds the package extension `QTLSimCatlabExt`. The life cycle is written
once, as a morphism of a free symmetric monoidal category whose objects are the stages
of a generation and whose generators are its steps, and functors lower it to the kernel,
to a model for the simulator SLiM (in its language Eidos), and to the deterministic
moments of the population:

```julia
import Catlab
using Catlab: ⋅, ⊗
using Random
L = lifecycle()
qtsim = L.develop ⋅ L.select ⋅ L.reproduce ⋅ L.record ⋅ L.rescale ⋅ L.mutate ⋅ L.swap

run_composed(World(p), Xoshiro(1); cycle = qtsim)   # the same run as run!, bit for bit
print(slim_script(qtsim, p))                        # the same life cycle as a SLiM model
run_moments(p; cycle = qtsim)                       # its moment equations

meta = qtsim ⊗ qtsim ⊗ qtsim ⋅ L.migrate(3)         # three demes coupled by migration
run_composed(p ⊗ demes(3), Xoshiro(1); cycle = meta, migration_rate = 0.01)
print(slim_script(meta, p; migration_rate = 0.01))  # the same metapopulation in SLiM
```

Steps can carry parameters (`withparams`), demes can have their own `Params`
(`[p1, p2, p3] ⊗ demes(3)`), and `operators()` gives the kernel's operators to compose
directly. The population is also an attributed C-set (ACSet, `lift`, `lower!`).

## Features

- **Kernel**: QTSim's generation in QTSim's order, including its single-crossover
  recombination (`CrossoverModel = :qtsim`, the default; `:all` applies every crossover)
  and its `0.4/ω` fitness normalisation. Genomes are shared blocks of about √`NoQTL` loci,
  so a generation of 1000 diploids takes 0.25 ms at 100 loci and 0.78 ms at 5000, 1.9 to
  17 times faster than QTSim (one thread, Apple M3 Pro; `validate/benchmark.log`).
- **Models**: `OneTraitModel()` and `TraitModel(n)`, two families of lowerings of the
  one life cycle, with step parameters typed by the model (`step_schema`).
- **Beyond QTSim** (off by default): Beverton–Holt density regulation
  (`Demography = :beverton_holt`), fecundity selection (`Selection = :fecundity`),
  periodic and red-noise optima, metapopulations with island and stepping-stone
  migration, and `NoTraits = n` traits with mutational and selection covariances.
- **Theory**: Gaussian, multilocus and rare-alleles equilibria of the genetic variance,
  house of cards, the steady lag and critical rate under a moving optimum, fluctuating
  optima, extinction-time distributions and rescue probability, neutral F_ST/Q_ST among
  demes, and the G matrix of several traits.
- **Sweeps**: `sweep` and `replicate` run jobs on all threads with a resumable CSV cache,
  deterministic by seed.

See [`NEWS.md`](NEWS.md) for the full list and for behaviour to know about.

## Validation

| check | script | result |
|:--|:--|:--|
| exact parity | [`validate/gsl_locked.jl`](validate/gsl_locked.jl) | a replica drawing QTSim's GSL random numbers in QTSim's order is bit-identical to QTSim every generation |
| statistical parity | [`validate/compare_to_qtsim.jl`](validate/compare_to_qtsim.jl) | Welch t-tests, kernel vs QTSim |
| published results | [`validate/published_sweep.log`](validate/published_sweep.log) | on about 590 of the paper's inputs, time to extinction and pre-move V_G within 1% of the published values (95% CIs inside 0.98–1.03); `CrossoverModel = :all` is rejected |
| SLiM vs QTSim | [`validate/slim/valid_generated.log`](validate/slim/valid_generated.log), [`three_way_confirm.log`](validate/slim/three_way_confirm.log) | no discrepancy; 40,000 runs per program, mean V_G within 0.6% |
| kernel vs SLiM, Julia-only models | [`validate/slim/`](validate/slim) `valid_*.log` | demography, fecundity, fluctuating optima, demes, migration, several traits: no discrepancy |
| F_ST/Q_ST theory | [`validate/fst.log`](validate/fst.log) | 41 comparisons, largest \|z\| 2.57 |
| moments lowering | [`validate/moments_regimes.log`](validate/moments_regimes.log), [`moments_loci.log`](validate/moments_loci.log), [`moments_demes.log`](validate/moments_demes.log) | rare-alleles closure across regimes (122 criteria pass, 4 fail), multilocus model, demes |

From the package root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'                    # tests
make -C reference/QTSim                                          # build QTSim (needs a C compiler and GSL)
julia --project=. validate/gsl_locked.jl                         # exact parity with QTSim
julia --project=. validate/compare_to_qtsim.jl                   # statistical parity with QTSim
julia -t 10 --project=vignettes validate/published_sweep.jl      # against the published results
julia --project=validate/slim validate/slim/check_generated.jl  # generated SLiM model vs the port (needs SLiM)
julia -t auto --project=validate/slim validate/slim/compare_to_qtsim.jl   # SLiM model vs QTSim
```

CI builds QTSim and SLiM 5.2 and runs the exact, statistical and SLiM checks on every
push. The published-results check and vignettes 07, 08, 09 and 11 need the paper's code
and data, which are not part of this repository: they look for them in a `Repository/`
folder next to this package's folder (the paper's
[Repository.zip](http://www.evolution.reading.ac.uk/QuantitativeGenetic/Repository.zip)
unpacked).

## Documentation and vignettes

The [documentation](https://ecorecipes.github.io/QTLSim.jl/dev/) has reference pages on
the parameters, the output, running QTSim's input files and several traits, the API, and
thirteen tutorials rendered from the Quarto vignettes in [`vignettes/`](vignettes/):

1. [Getting started](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/01_getting_started/)
2. [Mutation–selection balance](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/02_mutation_selection_balance/)
3. [Adaptation to a moving optimum](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/03_moving_optimum/)
4. [Loci and recombination](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/04_loci_and_recombination/)
5. [The categorical layer](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/05_categorical_layer/)
6. [Recomposing the model](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/06_recomposing_the_model/)
7. [One life cycle, two backends](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/07_one_life_cycle_two_backends/)
8. [Validation against QTSim](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/08_validation_against_qtsim/)
9. [Replicating Pagel, Gardner and Meade (2025)](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/09_replicating_pagel_2025/)
10. [Demography and evolutionary rescue](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/10_demography_and_rescue/)
11. [Density regulation in the paper's populations](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/11_demography_in_the_papers_populations/)
12. [Metapopulations and migration](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/12_metapopulations/)
13. [Several traits and pleiotropy](https://ecorecipes.github.io/QTLSim.jl/dev/tutorials/13_multiple_traits/)

## Proofs in Lean

The facts that justify the decomposition and the lowerings (the heritability algebra,
stabilising selection, the mutational input, unbiased recombination, a generation as a
composite and a run as its iterate, a lowering fixed by its images on the steps, demes,
migration and the monoidal structure, the Gaussian moment recursion) are proved in Lean 4
with Mathlib, in the library `CategoricalQuantGen` of the project this package comes from.
It is not part of this repository; docstrings cite it as `lean/CategoricalQuantGen/…`.

## Licence

Copyright (C) 2026 Simon Frost.

QTLSim.jl is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version (GPL-3.0-or-later). See
[`LICENSE`](LICENSE).

It reimplements the algorithm of QTSim (copyright 2023 Andrew Meade, University of
Reading), written from QTSim's source, which is distributed under the same licence, and
includes QTSim's source, unmodified, in [`reference/QTSim`](reference/QTSim) for testing
and verification (its README gives the provenance).

## Citing

If you use QTLSim.jl, please cite the paper whose model it implements:

> Pagel, M., Gardner, J. D. and Meade, A. (2025). Quantitative-genetic analysis of
> directional adaptation suggests low maximum sustainable rates of change in agreement
> with data from field populations. *Scientific Reports* 15.
> https://doi.org/10.1038/s41598-025-24445-2

and the software itself (see [`CITATION.cff`](CITATION.cff)):

> Frost, S. (2026). QTLSim.jl: a Julia port and categorical reformulation of QTSim
> (version 0.1.0). https://github.com/ecorecipes/QTLSim.jl

<!-- TODO(citation): software DOI. None has been minted yet; add it here and in CITATION.cff if one is. -->
