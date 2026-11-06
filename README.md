# QTLSim.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/QTLSim.jl/dev/)

QTLSim.jl is a Julia port of QTSim, the simulation behind Pagel, Gardner and Meade,
*Quantitative-genetic analysis of directional adaptation suggests low maximum
sustainable rates of change in agreement with data from field populations*
([Scientific Reports, 2025](https://doi.org/10.1038/s41598-025-24445-2)). QTSim is a C
program by Andrew Meade (University of Reading); its source calls it QTLSim 1.0.

The model follows a population of diploid individuals with one additive quantitative
trait under Gaussian stabilising selection on viability, toward an optimum that can
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

```julia
using QTLSim

p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, NoGenerations = 12_000,
           PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
           MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0,
           MoveFitnessSDGen = 10_000, MoveFitnessSD = 0.05, MoveType = :continuous)
out = simulate(p; seed = 1)        # per-generation statistics
```

`Params` takes QTSim's options under the same names. With more than one locus, gametes
recombine as QTSim recombines them by default (`CrossoverModel = :qtsim`): QTSim draws
a Poisson number of crossover points but its copying loop never advances past the first,
so every gamete has at most one crossover. The published simulations were run that way.
`CrossoverModel = :all` applies every point instead.

## The categorical layer

Loading Catlab adds a categorical description of the same model, the package extension
`QTLSimCatlabExt`. The life cycle is written once, as a morphism of a free symmetric
monoidal category whose objects are the stages of a generation and whose generators are
its steps, and functors lower it to the simulation kernel and to a model for the
population-genetics simulator SLiM, in its language Eidos:

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
steps are operators that can be rearranged and stratified into demes, and the population
is also an attributed C-set (ACSet):

```julia
(; develop, select, reproduce, mutate, swap, generation) = operators()
p ⊗ demes(4)                       # a metapopulation of four demes
```

The facts that justify the decomposition and the lowering are proved in Lean 4 in the
project this package comes from (cited in comments as `lean/CategoricalQuantGen/…`).

## Validation

QTLSim reproduces QTSim statistically (the same distributions of genetic variance,
survival and adaptation over many runs), exactly (a replica that draws QTSim's GSL random
numbers in QTSim's order matches its output generation by generation,
`validate/gsl_locked.jl`), and against the paper's published per-run results (on 589 of
the paper's inputs, genetic variance before the move and time to extinction within about
1% of the published values, `validate/published_sweep.jl`). Applying every crossover point
instead raises both measurably, so the published results depend on QTSim's recombination.

The SLiM model generated from the life cycle (`slim_script`) is validated in
[`validate/slim`](validate/slim): run for run against a SLiM port of QTSim written by hand,
and in distribution against QTSim, including a test of 40,000 runs per program that found
no difference in genetic variance between QTSim, SLiM and the kernel (95% intervals about
±2%). CI builds SLiM 5.2 and runs the first two on every push.

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

QTSim's C source comes with the package, unmodified, in
[`reference/QTSim`](reference/QTSim) (its README gives its provenance and how to build
it), and CI builds it and runs both parity checks on every push. The published-results
check and vignette 09 also need the paper's data, which is not part of this repository:
they look for it in a `Repository/` folder next to this package's folder, the paper's code
and data unpacked, with the input files under `Repository/Simulations/` and the results
tables under `Repository/Results/`. Vignettes 06 and 08 use QTSim, and 06 also SLiM
(version 4 or later).

## Documentation and vignettes

The [documentation](https://ecorecipes.github.io/QTLSim.jl/dev/) includes nine
tutorials, rendered from the Quarto vignettes in [`vignettes/`](vignettes/), which go
from a first simulation to a replication of the paper's results.

## Licence

Copyright (C) 2026 Simon Frost.

QTLSim.jl is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version (GPL-3.0-or-later). See
[`LICENSE`](LICENSE).

It reimplements the algorithm of QTSim (copyright 2023 Andrew Meade, University of
Reading), which is distributed under the same licence, and includes QTSim's source,
unmodified, in `reference/QTSim` for testing and verification.
