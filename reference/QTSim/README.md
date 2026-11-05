# QTSim, the reference C program

`src/` holds the C source of QTSim, the simulation that QTLSim.jl ports, as the authors
published it with Pagel, Gardner and Meade, *Quantitative genetic analysis of directional
adaptation suggests low maximum sustainable rates of change in agreement with data from
field populations*, Scientific Reports (2025),
[doi:10.1038/s41598-025-24445-2](https://doi.org/10.1038/s41598-025-24445-2). QTSim was
written by Andrew Meade (University of Reading); its source calls it QTLSim 1.0.

It is bundled for testing and verification, not used by the package itself:

- `validate/gsl_locked.jl` checks that a replica of QTLSim's kernel that draws QTSim's
  random numbers, in QTSim's order, reproduces QTSim's output bit for bit;
- `validate/compare_to_qtsim.jl` checks that QTLSim and QTSim give the same
  distributions;
- the `QTSim parity` job of the package's CI builds QTSim from this source and runs both.

## Licence and provenance

QTSim is free software: copyright 2023 Andrew Meade, distributed under the GNU General
Public License, version 3 or (at your option) any later version, as the header of every
file in `src/` states. The licence text is the package's `LICENSE`.

The files in `src/` are unmodified copies of `Code/QTSim/src/` in the paper's code and
data, http://www.evolution.reading.ac.uk/QuantitativeGenetic/Repository.zip (the paper's
data availability statement). `SHA256SUMS` lists their checksums, which `make check`
verifies. This README, the `Makefile` and `SHA256SUMS` are not part of QTSim. The rest of
the paper's code and data (QTSim's Python scripts, manual and example input, the input
files and results of its simulations) is not included.

## Building

QTSim needs a C99 compiler and GSL, the GNU Scientific Library, with `gsl-config` on the
`PATH` (Debian and Ubuntu: `libgsl-dev`; Homebrew: `gsl`). From the package root:

```sh
make -C reference/QTSim check     # the source is unmodified
make -C reference/QTSim           # builds reference/QTSim/build/QTSim
```

The validation scripts and vignettes use, in this order, the program named by the
`QTSIM` environment variable, this build, and the binary in the paper's code if it is
unpacked next to the package (`Repository/Code/QTSim/src/QTSim`).

QTSim reads an options file and prints a table to standard output:

```sh
reference/QTSim/build/QTSim options.txt > out.tsv
```

`validate/gsl_locked.jl` (`write_opts`) writes options files from QTLSim's `Params`,
whose fields have QTSim's option names.
