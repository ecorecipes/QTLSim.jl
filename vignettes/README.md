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
| `04_loci_and_recombination` | how gametes are made, QTSim's single-crossover recombination, and the cost of more loci |
| `05_categorical_layer` | the population as an ACSet, and a generation as a composite of operators |
| `06_recomposing_the_model` | changing the model by changing the composite: selection, mutation, a fluctuating optimum, demes and migration |
| `07_one_life_cycle_two_backends` | the life cycle written once as a morphism of a free symmetric monoidal category and lowered by functors to the kernel and to a SLiM model, which is the hand-written SLiM port run for run and agrees with QTSim with one locus or a hundred; demes coupled by migration in both backends |
| `08_validation_against_qtsim` | statistical, exact and published-results checks against QTSim |
| `09_replicating_pagel_2025` | the paper's results, recomputed from its table and replicated by rerunning 710 of its populations, the survivors to the full 240,000 generations |
| `10_demography_and_rescue` | a density-regulated population instead of QTSim's fixed size: how much faster populations die, and how the sustainable rate falls, with the moments lowering's predictions |
| `11_demography_in_the_papers_populations` | vignette 09's 710 populations rerun with density regulation (R = 5, 2, 1.5): how survival, time to extinction and the sustainable rate move from the paper's fixed-size estimates |
| `12_metapopulations` | demes and migration: building metapopulations, how the kernel's migration exchanges individuals, neutral Q_ST among demes against exact island and stepping-stone theory, local adaptation (Q_ST > F_ST), and the moments lowering's prediction of the variance among demes |
| `13_multiple_traits` | several traits with pleiotropy: the model as a second family of lowerings of the same life cycle, the G matrix at equilibrium against theory, correlated response and the line of least resistance, the mutational correlation, local adaptation of demes in two traits, a fluctuating optimum vector, and the kernel against SLiM |
| `14_plasticity_and_sustainable_rates` | phenotypic plasticity: a reaction norm of fixed slope expressed in `select`, the exact correspondence of a plastic population at rate k with a non-plastic one at (1 − b/B)k, survival against the rate of change from the 12,320-run experiment (and where the prediction of a cost fails), the cost trade-off, a heritable slope (kernel against moments), and a delayed cue under red noise |

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
is captured into the first cell of the rendered file. Vignettes 02, 03, 04, 07, 08, 09, 10,
11, 12 and 13 run their simulations on all threads (`julia: exeflags: ["--threads=auto"]`).

`scripts/run_vignettes.jl` runs the vignettes' code without Quarto, each vignette in its
own Julia process, and fails if any chunk throws; it does not write the rendered files.
From the package root:

```sh
julia -t 3 --project=vignettes scripts/run_vignettes.jl          # 01-06, 10, 12 and 13
julia -t 3 --project=vignettes scripts/run_vignettes.jl 03 10    # by number or name
```

Its default set is the vignettes that need neither QTSim, SLiM nor the paper's data. CI
runs 01, 02, 03, 05, 06, 10, 12 and 13 on every push (about 7 minutes on three threads locally) and
04, which takes about 6 minutes on its own, on a weekly schedule and on manual runs.

`julia ../scripts/sync_vignettes.jl` copies the rendered Markdown and figures into
`docs/src/tutorials/` for the documentation site, and CI runs
`sync_vignettes.jl --check`, which fails if the copies are stale.

## QTSim, the paper's data, and SLiM

Vignettes 07 and 08 run QTSim, which they find as the validation scripts do
(`validate/qtsim.jl`): built from the source that comes with the package
(`make -C reference/QTSim`, which needs a C compiler and GSL), or in the paper's code and
data if they are unpacked in a `Repository/` folder next to the package's folder.
Vignette 08 skips the sections that need QTSim when there is none, and uses one of the
paper's input files (`Repository/Simulations/`) when they are there. Vignettes 09 and 11 need
the paper's results tables, `Repository/Results/`. Vignette 07 also needs SLiM, version 5.2
(the version tested) (`slim/build/slim` next to the package's folder, `ENV["SLIM"]`, or `slim` on the
path); the hand-written SLiM port of QTSim it compares with comes with the package,
in `validate/slim/qtsim.slim`.

Vignette 09's simulations are stored next to it and read back on later renders:
`replication_runs.csv` holds the 710 reruns of 2000 generations after the move (about 45
minutes of single-thread CPU time to regenerate, roughly ten minutes on ten threads) and
`survivor_runs.csv` the 63 of them still alive then, run on to 240,000 generations (about
an hour of CPU time, roughly a quarter of an hour on eight threads; the longest run, at 5000
loci, takes 3 to 6 minutes).
Vignette 11 stores its reruns the same way, one file per growth rate
(`demography_runs_R5.0.csv`, `demography_runs_R2.0.csv`, `demography_runs_R1.5.csv`;
about 40 minutes of CPU time each, roughly ten minutes on ten threads), and reads
vignette 09's `replication_runs.csv` for the same populations at a fixed size.
These times are estimates from the single-thread speed of the kernel measured by
`validate/benchmark.jl` (12-core Apple M3 Pro, Julia 1.12, 2026-11-07); several threads
take up to about twice the CPU time of one.
