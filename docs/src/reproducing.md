# Running QTSim's input files

QTSim reads its options from a text file: one `Keyword value…` per line. The paper's
simulations are 119,984 such files. [`read_qtsim_options`](@ref) turns one into a
[`Params`](@ref), so QTLSim can run any QTSim input:

```julia
using QTLSim
p = read_qtsim_options("Repository/Simulations/QT0100-Possion3.0/Input-00001.txt")
out = simulate(p)                         # uses the file's seed, p.Seed
```

and [`write_qtsim_options`](@ref) writes a `Params` as a file QTSim can run:

```julia
write_qtsim_options("input.txt", Params(p; NoGenerations = 20_000); seed = 1)
run(`QTSim input.txt`)                     # QTSim prints its table to standard output
```

The same seed does not give QTSim's run, because the random numbers differ (see
[Output](@ref)). Compare QTSim and QTLSim over many runs, as `validate/compare_to_qtsim.jl`
does.

## QTSim itself

QTSim's C source comes with the package, unmodified, in `reference/QTSim`, with a
Makefile. With a C compiler and GSL, `make -C reference/QTSim` builds
`reference/QTSim/build/QTSim`. The validation scripts find it there, or take the program
named by `ENV["QTSIM"]`.

## The paper's data

The paper's simulations and results are not part of the package. Its data-availability
link (http://www.evolution.reading.ac.uk/QuantitativeGenetic/Repository.zip) unpacks to
a `Repository/` folder. Placed next to the package's folder, it gives:
- the input files in `Repository/Simulations/QT<loci>-Possion<u>/`;
- the results table in `Repository/Results/Master Output file.xlsx`, one row per input.
  Match rows to files by their parameters, not by position.

The table numbers generations from 1. Its "at 10000" columns are QTSim's (and QTLSim's)
output for generation 9999, and "Generations to extinction" is the dying generation
minus 9999 (240,000 for a survivor).

With that folder in place:

- `validate/published_sweep.jl` reruns a stratified sample of the published inputs in
  QTLSim and compares them with the table, input by input.
- The vignette [Replicating Pagel, Gardner and Meade (2025)](tutorials/09_replicating_pagel_2025.md)
  rebuilds the paper's results: survival, rates of change, lag, mortality, genetic
  variance and the critical rate k_c.
- `validate/gsl_locked.jl` includes six of the published inputs in its bit-identity check
  against QTSim.
