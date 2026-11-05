# Vignettes

Each vignette lives in its own directory `NN_name/NN_name.qmd` and is rendered with
[Quarto](https://quarto.org) using the native Julia engine to two formats: HTML
(`NN_name.html`, self-contained) and GitHub-flavoured Markdown (`NN_name.md` plus
`NN_name_files/`). Both outputs are committed.

| vignette | what it shows |
|:--|:--|
| `01_getting_started` | parameters, a first run, reproducibility, a moving optimum |
| `02_mutation_selection_balance` | the mutational input, and equilibrium genetic variance against the Gaussian and house-of-cards approximations |
| `03_moving_optimum` | the three ways the optimum can move, an extinction, the paper's measures, survival against the rate of change |
| `04_categorical_layer` | the population as an ACSet, and a generation as a composite of operators |
| `05_recomposing_the_model` | changing the model by changing the composite: selection, mutation, a fluctuating optimum, demes and migration |
| `06_loci_and_recombination` | how gametes are made, QTSim's single-crossover recombination, and the cost of more loci |
| `07_validation_against_qtsim` | statistical, exact and published-results checks against QTSim |
| `08_replicating_pagel_2025` | the paper's results, recomputed from its table and replicated by rerunning 1000 of its populations |

## Rendering

```sh
cd vignettes
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
../scripts/render_vignette.sh                     # every vignette
../scripts/render_vignette.sh 03_moving_optimum   # one vignette
```

The script renders each vignette twice, HTML first: in a single render of both formats,
the self-contained HTML step can delete the figure files that the Markdown needs. Make
sure the package is precompiled before rendering, otherwise Julia's precompilation output
is captured into the first cell of the rendered file. Vignettes 02, 03, 06, 07 and 08 run
their simulations on all threads (`julia: exeflags: ["--threads=auto"]`).

`julia ../scripts/sync_vignettes.jl` copies the rendered Markdown and figures into
`docs/src/tutorials/` for the documentation site, and CI runs
`sync_vignettes.jl --check`, which fails if the copies are stale.

## The paper's code and data

Vignettes 07 and 08 use the code and data published with the paper, which are not part
of this repository. They look for them in a `Repository/` folder next to the package's
folder: QTSim built at `Repository/Code/QTSim/src/QTSim`, the input files under
`Repository/Simulations/`, and the results tables under `Repository/Results/`. Vignette
07 skips the sections that need QTSim when it is missing; vignette 08 needs the results
tables. Vignette 08's simulations are stored in
`08_replicating_pagel_2025/replication_runs.csv` and read back on later renders;
regenerating them takes about an hour on ten threads.
