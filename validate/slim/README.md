# The SLiM backend: validation

QTLSim writes QTSim's life cycle once, as a morphism of a free symmetric monoidal category,

```
develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅ mutate ⋅ swap : Zygotes → Zygotes
```

and lowers it by functors to its simulation kernel and to Eidos, the scripting language of
[SLiM](https://messerlab.org/slim/) (`lower`, `slim_script`; see vignette
`06_one_life_cycle_two_backends`). This folder validates the Eidos lowering: against a SLiM
port of QTSim written by hand, run for run, and against QTSim itself, in distribution.

## Files

- `qtsim_generated.slim` — the SLiM model generated from QTSim's life cycle by
  `slim_script`. Do not edit it: regenerate it with `check_generated.jl --write`. The package
  tests check that `slim_script` still writes it unchanged.
- `qtsim.slim` — a nonWF SLiM port of QTSim written by hand (single locus, one deme). It is a
  **frozen reference**: an independent statement of QTSim in Eidos, written without the
  life-cycle machinery, that the generated model must reproduce. New features go through
  the life cycle and `slim_script`, not here.
- `check_generated.jl` — checks that `slim_script` writes `qtsim_generated.slim`, then
  checks that model against `qtsim.slim`, run for run (30 runs, including extinctions).
  Exits 1 on any difference.
- `compare_to_qtsim.jl` — statistical validation of a SLiM model (by default the generated
  one) against QTSim, scenarios A–E, KS tests, 20 runs per program. Exits 1 if any p < 0.001.
  Last results in `valid_generated.log` (generated model) and `valid.log` (`qtsim.slim`,
  scenarios A and B, run 2026-07-24 with an earlier version of the harness).
- `three_way_confirm.jl` — the large pre-specified test of QTSim vs the generated SLiM model
  vs the kernel (40,000 runs each, about 2.5 h on 8 threads); result in `three_way_confirm.log`.
- `slim.jl` — `find_slim()`: `ENV["SLIM"]`, else `slim/build/slim` next to the package's
  folder, else `slim` on the path. QTSim is found by `../qtsim.jl`.
- `Project.toml` — the environment for these scripts: QTLSim from this repository, Catlab
  and HypothesisTests.

CI's `SLiM backend` job builds SLiM 5.2 and QTSim and runs `check_generated.jl` and
`compare_to_qtsim.jl`.

## Run

SLiM 5.2: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --target slim`
in a clone of https://github.com/MesserLab/SLiM (it bundles its own GSL). QTSim:
`make -C reference/QTSim`. From the package root:

```sh
julia --project=validate/slim -e 'using Pkg; Pkg.instantiate()'           # once
julia --project=validate/slim validate/slim/check_generated.jl [--write]   # ~15 s
julia -t auto --project=validate/slim validate/slim/compare_to_qtsim.jl [model.slim] [R]   # ~3 min on 8 threads
julia -t 8 --project=validate/slim validate/slim/three_way_confirm.jl [R]  # ~2.5 h

# one run of a model (TSV to stdout): gen  phenoMean  genoVar  nSurv  theta  omega
slim -s 42 -d "N=1000" -d "NGEN=6500" -d "moveGen=5000" -d "moveSD=1.0" validate/slim/qtsim_generated.slim
```

Parameters can be overridden with `-d "NAME=value"`: `N VP h2 m rate NGEN moveGen moveSD
changeGen changeSD initOmega` in both models, and `moveType` (`"phenotype"`, `"fitness"`,
`"continuous"`), `L` (loci), `RECOMB` and `migRate` in the generated one. `initOmega < 0`
uses the analytical width `ω = V_G/√(m·V_E)`. In `qtsim.slim` the width change is tied to
`moveGen`.

## Design

QTSim represents the trait as a **continuum of alleles**: real-valued allelic effects,
mutation *adds* Gaussian increments, and the genotype is their sum. SLiM's
discrete-mutation machinery models that poorly, so both models carry the allelic values
per individual explicitly (`setValue`/`getValue`) and use SLiM only for the **nonWF life
cycle**; genetics initialisation is vestigial (`mutationRate 0`). Each generation follows
QTSim's `RunSim` order inside one `reproduction()` call:

1. **develop** — `z = G + env`, `env ~ N(0, V_E)` drawn fresh per generation.
2. **select** — the optimum moves (from `moveGen`) to `mean(z) + √V_P·moveSD`; Gaussian
   fitness `W = dnorm(z−θ, 0, ω)/(0.4/ω)` (QTSim's exact form, including the 0.4);
   **viability**: survive iff `runif() < W`.
3. **reproduce** — `N` offspring by uniform random mating of two **distinct** survivors;
   with many loci each gamete draws `Poisson(RECOMB)` crossover points and switches strands
   at the lowest one only, as QTSim does.
4. **record** — one TSV line for the parent generation (QTSim's `Output(CPop)`).
5. **mutate** — each offspring chromosome gets `Poisson(L·rate)` increments of
   `N(0, MUTVAR)` at uniformly chosen loci, `MUTVAR = m/(2·L·rate)·V_E`.

Parents are then removed, so generations do not overlap. A life cycle on several demes,
such as `g ⊗ g ⊗ g ⋅ migrate(3)`, generates a model with one subpopulation per deme and
migration at `migRate`.

SLiM's random numbers differ from QTSim's GSL stream, so the agreement with QTSim is
**statistical** (equality of distribution), not the bit-identical match that
`../gsl_locked.jl` achieves for the Julia kernel. Between the two SLiM models it is exact:
with one locus the steps of the generated model are the Eidos statements of `qtsim.slim`,
in the same order, so for the same seed both print the same table, number for number. The
one difference is intended: at extinction `qtsim.slim` also prints the line of the
generation that died, which QTSim does not.

## Results

**Run for run (2026-11-06):** 30 of 30 runs identical, 10 of them going extinct.

**Against QTSim (2026-11-06, `valid_generated.log`, 20 runs per program):**
- Scenario A (stabilising equilibrium) and scenario B (single-locus directional adaptation)
  give the hand-written port's numbers exactly.
- Scenario C (QTSim's continuous move of the optimum) is consistent: KS p ≥ 0.33.
- Scenario D (100 loci, stabilising selection): p = 0.17 for the genetic variance and
  p = 0.035 for survival.
- Scenario E (100 loci, the optimum moving continuously; every run goes extinct in both
  programs): p = 0.56 for the genetic variance before the move and p = 0.17 for the time to
  extinction.

Vignette 06 repeats D and E with 50 runs per program, adding the kernel; there the SLiM
model's p-values against QTSim are 0.18 to 0.71, survival in D included.

**Is the agreement real? (2026-11-05)** Low p-values like the 0.035 prompted a fuller
investigation:
- **Code audit:** QTSim's live code paths (`RunSim`, `CreatePop`, `SetIndividualEnv`,
  `MutateIndividualChrNo`, `CreateCombMap`/`RecombinationPoissonChr`, `NewPop`, `Fitness`,
  `SetFitnessMean`, `Output`) are algorithmically the same as the kernel and the Eidos.
- **Random numbers:** SLiM's and GSL's generators check out against theory (10⁸ draws each).
- **SLiM invariants:** population size, generation turnover and mutation counts are exact.
- **Micro-tests:** 10,000 short runs each of initial survival, drift and stabilising
  selection.
- **Trajectory study:** every generation of 1,000 runs of scenario B and 500 of D; 0 of 40
  SLiM-vs-QTSim tests below 0.05.

All agreed, except that in the rare-mutation regime QTSim's mean genetic variance came out
1–5% above both SLiM and the kernel in each batch. A pre-specified confirmatory test settled
it: `three_way_confirm.jl`, 40,000 runs per program, random seeds. Differences in mean V_G
at generation 499, with Welch p:
- SLiM − QTSim: −0.56% (p = 0.57; KS 0.61, AD 0.36)
- kernel − QTSim: +0.26% (p = 0.79)
- SLiM − kernel: −0.82% (p = 0.41)

The 95% intervals are about ±2%, and the upper quantiles agree too. The earlier gap was a
selected chance result; SLiM and the kernel appeared to deviate "together" only because both
were compared with the same QTSim sample.
