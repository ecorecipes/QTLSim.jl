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
  tests check that `slim_script` still writes it unchanged, and the same for three more
  generated models (`golden.jl`): `qtsim_generated_loci.slim` (100 loci),
  `qtsim_generated_demes.slim` (three demes, island migration) and
  `qtsim_generated_stepping.slim` (three demes on a line of stepping stones, a migration
  matrix), which `check_generated.jl` also runs in SLiM.
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
- `compare_demography.jl` — the kernel against the generated SLiM model under density
  regulation (`Demography = :beverton_holt`), which QTSim lacks; KS tests, 40 runs each,
  about 15 minutes on 10 threads; exits 1 if any p < 0.001; last result in
  `valid_demography.log`.
- `compare_fluctuating.jl` — the kernel against the generated SLiM model with a
  fluctuating optimum (`MoveType = :periodic` and `:red_noise`), which QTSim lacks; KS
  tests, 40 runs each, about 4 minutes on 3 threads; exits 1 if any p < 0.001; last
  result in `valid_fluctuating.log`.
- `compare_demes.jl` — the kernel against the generated SLiM model on demes with parameters
  of their own (`slim_script(cycle, [p1, p2, p3])`) and on steps with parameters
  (`withparams`), which QTSim lacks; KS tests, 50 runs each, about 4 minutes on 3 threads; and
  a run-for-run check that steps with parameters are the model of `Params` with those values;
  exits 1 if any p < 0.001 or any run differs; last result in `valid_demes.log`.
- `compare_migration.jl` — the kernel against the generated SLiM model on three neutral
  demes coupled by migration, which QTSim lacks: a scalar rate, the island matrix and a
  line of stepping stones; KS tests, 300 runs each, about an hour on 3 threads; exits 1
  if any p < 0.001; last result in `valid_migration.log`.
- `three_way_confirm.jl` — the large pre-specified test of QTSim vs the generated SLiM model
  vs the kernel (40,000 runs each, about 2.5 h on 8 threads); result in `three_way_confirm.log`.
- `slim.jl` — `find_slim()`: `ENV["SLIM"]`, else `slim/build/slim` next to the package's
  folder, else `slim` on the path. QTSim is found by `../qtsim.jl`.
- `Project.toml` — the environment for these scripts: QTLSim from this repository, Catlab
  and HypothesisTests.

CI's `SLiM backend` job builds SLiM 5.2 and QTSim and runs `check_generated.jl`,
`compare_to_qtsim.jl` and `compare_demography.jl`, the last with 10 runs per backend on
every push (a smoke test: so few runs flag only gross differences) and with its default
of 40 on the weekly scheduled run.

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
`"continuous"`, and `"red_noise"` or `"periodic"` when the model was generated for one of
them, with `moveNoiseSD`, `moveAutocorr`, `moveAmplitude` and `movePeriod`), `L` (loci),
`RECOMB` and `migRate` in the generated one. `initOmega < 0`
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
migration at `migRate` (or along the matrix `MIG`). Migration visits the subpopulations in
a new random order each generation (`sample`), as the kernel does (`randperm`): the swaps
are sequential, so with a fixed order migration would not commute with relabelling the
demes (`lean/CategoricalQuantGen/Migration.lean`).

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

**Migration (2026-11-07, `compare_migration.jl`, `valid_migration.log`, 300 runs each),
with the demes visited in a random order in both backends:** three neutral demes of 200
individuals, 300 generations, at a scalar rate 0.05, along `island_migration(3, 0.05)` and
along a line of stepping stones (`stepping_stone(3, 0.05; circular = false)`). On the line
the squared difference of the two end demes' means is 2.14 (kernel) and 1.98 (SLiM) times
that of neighbours, against 0.87 to 1.07 for the island model. The KS p-values between the
backends, on the end demes' and the neighbours' squared differences and the genetic
variance within a deme, are 0.10 or more (nine tests). (Before the random order, a run on
the line with the fixed order gave 1.82 and 1.73, every p 0.29 or more.)

**Density regulation (2026-11-07, `compare_demography.jl`, `valid_demography.log`, 40 runs each):**
the kernel and the generated SLiM model with `Demography = :beverton_holt`, on vignette 10's
parameter set (100 loci). Equilibrium with R = 1.5: population size 971 in both (KS
p = 0.91), genetic variance p = 0.91. Moving 0.02 per generation with R = 2: all 40 die in
both, 583 and 562 generations after the move on average (p = 0.26). Moving 0.01 with
R = 1.5, a rate every population survives at a fixed size: all 40 die in both, after
1341 and 1326 generations (p = 0.57). QTSim has a fixed size, so the kernel is the
reference here.

**Fluctuating optima (2026-11-07, `compare_fluctuating.jl`, `valid_fluctuating.log`, 40 runs
each):** the kernel and the generated SLiM model with a periodic optimum (A = 0.1, P = 200;
A = 0.15, P = 50) and red noise (σ = 0.1, ρ = 0.9 and 0), one locus, ω = 0.3. No run died.
Across the 18 KS tests on the mean squared lag, genetic variance, survivors and the
optimum's SD and autocorrelation, the smallest p was 0.029 (mean squared lag, P = 200,
kernel 0.001193 against SLiM 0.001148). A rerun of that scenario alone with 200 runs each
and new seeds, its threshold fixed in advance, gave 0.001164 against 0.001183 (p = 0.18).
The mean squared lags are 1–8% above the Gaussian prediction (`expected_lag2`), which
ignores drift. QTSim has no fluctuating optimum, so the kernel is the reference here.

**Demes with parameters of their own (2026-11-07, `compare_demes.jl`, `valid_demes.log`, 50
runs each):** three demes with ω = 0.15, 0.3, 0.6 and 400, 500, 600 individuals, migration
0.01; and two demes whose steps carry optima ±0.1 and widths 0.2 and 0.4, migration 0.02, so
that their means differ by 0.0998 (kernel) and 0.0996 (SLiM). One locus, 1000 generations.
Across the 17 KS tests on each deme's genetic variance and mean phenotype and on the
difference of the means, the smallest p was 0.112. A life cycle whose steps carry a
continuous move and a mutation scalar gives the same output as `Params` with those values in
10 of 10 seeds. QTSim has no demes, so the kernel is the reference here.

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
