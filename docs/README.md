# Vignettes

Two Quarto vignettes illustrate the `CategoricalQuantGen` package:

- **`vignette.qmd`** → `vignette.html` — the two-level (abstract Catlab/ACSets ⟷
  lowered array kernel) design, the `⋅`/`⊗` operator sugar, and a final section
  demonstrating statistical **parity with the original C simulation QTSim**
  (empirical-CDF / KS comparison across scenarios, plus the RNG-locked
  bit-identity check in `../validate/gsl_locked.jl`).
- **`compositional.qmd`** → `compositional.html` — how the operator syntax lets
  you **change features of the model by recomposition**: swapping the selection
  operator (stabilizing / neutral drift / directional truncation), changing the
  mutation input (drop or duplicate `mutate`), and stratifying with `⊗` +
  `migration` — all reusing the validated kernel unchanged.

## Prerequisites

- [Quarto](https://quarto.org) ≥ 1.5 (provides the native `julia` engine).
- Julia (the repo is developed against 1.12).
- The QTSim reference binary, for the parity section:
  ```sh
  cd ../../Repository/Code/QTSim/src
  G=$(brew --prefix gsl)
  cc *.c -O2 -o QTSim -I"$G/include" -L"$G/lib" -lgsl -lgslcblas -lm
  ```

## Set up the docs environment (once)

From this directory (`docs/`):

```sh
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
```

This dev-tracks the local `CategoricalQuantGen` package and installs CairoMakie for
the figures.

## Render

```sh
quarto render vignette.qmd          # → vignette.html
quarto render compositional.qmd     # → compositional.html
```

Quarto's `julia` engine runs the chunks in this project. The parity section
shells out to the QTSim binary; if it isn't built, that section errors with a
build hint. The computational chunks (everything except the CairoMakie figures)
can also be exercised without Quarto by running them against the main package
project — see `../validate/compare_to_qtsim.jl` for the full parity study.
