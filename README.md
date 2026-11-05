# CategoricalQuantGen.jl

A compositional, AlgebraicJulia-flavoured reformulation of the quantitative-genetics
model in Pagel et al., *Quantitative genetic analysis of directional adaptation
suggests low maximum sustainable rates of change in agreement with data from field
populations* ([Scientific Reports, 2025](https://www.nature.com/articles/s41598-025-24445-2)),
and of QTSim, the C simulation used in that paper (Andrew Meade, University of
Reading). Evolutionary processes are Catlab operators over a fast array kernel
that reproduces QTSim's numerics, validated against QTSim and against the paper's
published simulations.

This package is one part of a larger project that also holds a Lean 4
formalisation of the facts the model relies on (cited below as `lean/…`), a SLiM
port, and the paper's code and data; none of those are in this repository. The
tests run on their own. The validation scripts and the vignettes' parity sections
look for the paper's code and data in a `Repository/` folder next to this
package's folder: QTSim built at `../Repository/Code/QTSim/src/QTSim`, the input
files under `../Repository/Simulations/` and the results table at
`../Repository/Results/Master Output file.xlsx`.

## Two-level design

| level | representation | where |
|-------|----------------|-------|
| **abstract** | population as a **Catlab/ACSets** attributed C-set; processes as first-class `Operator`s composed with **`⋅`** and stratified with **`⊗`** | `src/schema.jl`, `src/compose.jl` |
| **lowered** | fast array kernel (`(NoQTL, Ploidy, PopSize)` genome buffers) replicating QTSim's numerics | `src/params.jl`, `src/model.jl` |

`lift(::World) → Population` and `lower!(::Population, ::World)` translate between
the two; the round-trip preserves genome and genotype exactly. This mirrors the
project's north-star "abstract semantics → lowered backend" architecture.

## ACSet schema (`SchPopulation`)

```
Allele ──owner──▶ Individual        effect : Allele → ℝ
                                    env, genotype, phenotype, fitness : Individual → ℝ
```

An `Individual` owns `Ploidy·NoQTL` `Allele`s; its genotype is the sum of its
alleles' `effect`s, computed via the categorical preimage `incident(pop, i, :owner)`
(QTSim `CaclGenotype`). `Population = PopulationACSet{Float64,Float64}`.

## Categorical operators

Evolutionary processes are `Operator`s wrapping the validated kernel steps.
Composition uses Catlab's **diagrammatic** (left-to-right) `⋅` (`\cdot`):

```julia
generation = develop ⋅ select ⋅ reproduce ⋅ mutate ⋅ swap
```

`(f ⋅ g)(w, rng) = g(f(w, rng), rng)`, short-circuiting on extinction. This is the
reverse reading order of the mathematical `∘` in `lean/CategoricalQuantGen/Composition.lean`
(`generation = mutate ∘ reproduce ∘ select`). `run_composed(w, rng)` drives a run
through this composite and is **bit-for-bit identical** to the kernel `run!`
under the same seed.

Stratification uses Catlab's monoidal product **`⊗`** (`\otimes`):

```julia
meta = params ⊗ Demes(4)        # a metapopulation of 4 demes (Stratified)
generation(meta, rng)           # an operator acts on every deme pointwise
migrate!(meta, rng; rate=0.01)  # couples demes; conserves total population
migration(0.01)                 # the composable Operator form of migrate!
```

## Usage

```julia
using Pkg; Pkg.add(url="https://github.com/ecorecipes/CategoricalQuantGen.jl")
using CategoricalQuantGen

p   = Params(PopSize=1000, Ploidy=2, NoQTL=1, NoGenerations=6500,
             PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001,
             ChangeFitnessGen=5000, ChangeFitnessSD=0.30,
             MoveFitnessSDGen=5000, MoveFitnessSD=1.0, Seed=1)
out = simulate(p)               # NamedTuple of per-generation statistics
```

Multi-locus runs recombine by Poisson crossover. `CrossoverModel = :qtsim` (the
default) reproduces QTSim exactly: it draws `Poisson(RecombPoissonMean)` crossover
points but only the lowest one switches strands, because its point index never
advances, so every gamete has at most one crossover. The published simulations
were run that way. `CrossoverModel = :all` applies every drawn point instead.

## Tests & validation

From the package root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'           # unit + ACSet + composition tests
julia --project=. validate/compare_to_qtsim.jl          # statistical consistency vs QTSim (needs the C binary)
julia --project=. validate/gsl_locked.jl                # bit-identity vs QTSim, incl. shipped multi-locus inputs
julia -t 10 --project=docs validate/published_sweep.jl  # vs the paper's published sweep
```

The validation runs the Julia model and QTSim side by side over many replicates
and compares distributions (equilibrium genetic variance, adaptation lag,
survival, extinction rate) with Welch t-tests. Different RNGs mean runs are not
identical; the claim tested is that the two produce the **same distributions**.
To build the C reference, see `docs/README.md`.

`gsl_locked.jl` is a replica that draws all randomness from QTSim's GSL routines
in QTSim's order. Seeded like QTSim it reproduces every generation exactly, for
the single-locus protocol and for shipped multi-locus input files. A negative
control that applies every crossover point diverges within a few generations.

`published_sweep.jl` compares the model with `Repository/Results/Master Output
file.xlsx`, the per-run results of the paper's 119,984-run sweep. It first checks
that the local QTSim reproduces published rows exactly, then compares the Julia
model on a random sample of inputs from every stratum, paired by input:
pre-move genetic variance and survivors, and time to extinction after the
optimum starts moving. It needs the `docs` environment (set up as in `docs/README.md`)
for XLSX and HypothesisTests. Output: `published_sweep.log` and `published_sweep_results.csv`.

In two runs of about 590 inputs with independent seeds, the default model's ratio
to the published values was 1.00–1.01 for both time to extinction and pre-move
genetic variance, with every 95% interval inside 0.98–1.03. Applying every
crossover point instead lengthens survival by about 12% and raises genetic
variance by about 6%, so it is clearly rejected.

## Files

- `src/params.jl` — `Params` and derived quantities (mirror `lean/…/Heritability.lean`)
- `src/model.jl` — `World` state and the kernel operators / `run!` (replicate QTSim)
- `src/schema.jl` — Catlab ACSet schema, `lift`/`lower!`/`genotypes`
- `src/compose.jl` — `Operator`, `⋅`, `⊗`, `generation`, `run_composed`, `Demes`, `migrate!`
- `test/runtests.jl` — unit tests (kernel, ACSet round-trip, composition, stratification)
- `validate/compare_to_qtsim.jl` — statistical-consistency harness against QTSim (single locus)
- `validate/gsl_locked.jl` — RNG-locked bit-identity replica of QTSim (any number of loci)
- `validate/published_sweep.jl` — comparison with the paper's published sweep results
- `validate/published_sweep.log`, `validate/published_sweep_results.csv` — output of the last published-sweep run
- `docs/` — two Quarto vignettes with their rendered HTML (see `docs/README.md`)

## License

Copyright (C) 2026 Simon Frost.

CategoricalQuantGen.jl is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version (GPL-3.0-or-later). See [`LICENSE`](LICENSE).

It reimplements the algorithm of QTSim (copyright 2023 Andrew Meade, University
of Reading), which is distributed under the same licence.
