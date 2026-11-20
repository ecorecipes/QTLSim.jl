# QTSimFast: a faster QTSim with the same output, bit for bit

**This package contains a modified version of QTSim.** QTSim (its source calls itself
QTLSim 1.0) is by Andrew Meade, School of Biological Sciences, University of Reading,
copyright 2023, and is free software under the GNU General Public License, version 3 or
(at your option) any later version. It is the simulation of Pagel, Gardner and Meade,
*Quantitative genetic analysis of directional adaptation suggests low maximum sustainable
rates of change in agreement with data from field populations*, Scientific Reports (2025),
doi:10.1038/s41598-025-24445-2.

- `original/` is QTSim's C source, **unmodified**, as published with the paper.
- `fast/` is **QTSimFast, a version modified in November 2026 for speed**, prepared with
  QTLSim.jl (https://github.com/ecorecipes/QTLSim.jl), with Claude (Anthropic). It is
  distributed under the same licence (`COPYING`). Every modified or new file says at its
  top that it was changed, how, and when, above the original notice (GPL section 5).

For the same options file, and so the same seed, QTSimFast writes the same bytes to
standard output as QTSim and exits with the same code, for any number of threads. In the
benchmark (`bench_results.md`, quiet machine, 20 November 2026) it was 2.7 to 11.8
times as fast on one thread (100 to 5000 loci), 4.6 to 19.4 times with two threads and up
to 27 times with four. The
original remains the reference: QTSimFast is tested against it.

Revised on 19 November 2026 (`REPORT.md`, sections 3.9 to 3.13): with two threads or more
one thread now makes the random draws and another builds the genomes, which roughly
halves the wall time of a run at 100 to 1000 loci with two to four threads (benchmark in
`REPORT.md`, section 1; threads cost processor time, and one thread remains the default).

## Quick start

```sh
make                      # builds original/build/QTSim and fast/build/QTSimFast (needs GSL: gsl-config)
./demo_parity.sh          # runs both on six examples; prints checksums, sizes, times, IDENTICAL/DIFFERENT
make parity               # a small parity matrix, about a minute
fast/build/QTSimFast examples/design_1000_loci.txt > out.tsv
QTSIM_THREADS=4 fast/build/QTSimFast examples/design_5000_loci.txt > out.tsv
```

Requirements: a C99 compiler with C11 atomics and POSIX threads, `make`, and GSL (the GNU
Scientific Library) with `gsl-config` on the `PATH`. No Julia, no Python. Tested on macOS
arm64 only; written for Linux and x86-64 too, but not yet run there (run `make check` and
`make parity` first).

## What to read

- `REPORT.md`: the technical report. What was changed and why the output cannot change,
  measured speed, what bit parity rules out, threads, the start-up self-checks, how it
  was verified, limits.
- `PROOFS.md` (also `PROOFS.html`, `PROOFS.pdf`): the machine-checked proofs (Lean 4) that the algorithmic changes preserve
  the result. **The proofs are about a mathematical model of the algorithms, not about
  the C source**; the note says exactly what is proved, on which hypotheses, what links
  the model to the code, and what is not proved.
- `formal/QTSimFast.pdf` (or `.html`): the proofs themselves, as a document with the
  checked code. `formal/QTSimFastShadows.pdf`: checks that the theorems say what is
  claimed.

## Contents

| Path | What |
|---|---|
| `original/src/` | QTSim's source, byte-identical to the published files (`original/SHA256SUMS`; `make check`) |
| `fast/src/` | QTSimFast: `QTSim.c` (simulation core rewritten), `QTSim.h`, new `FastRNG.*`, `FastThreads.*`; `GenLib.*`, `PassOptionFile.*`, `RandLib.*` unmodified |
| `Makefile` | builds both programs and the test variants; targets `check`, `demo`, `parity`, `parity-full`, `model-check`, `bench`, `docs` |
| `examples/` | six option files written for this package (the paper's own input files are not included) |
| `demo_parity.sh`, `PARITY_DEMO.txt` | the parity demonstration and a transcript of it |
| `test_parity.sh`, `PARITY_FULL.txt` | the parity matrix (`quick`, `full`) and the result of the full run: 13,915 comparisons, all identical |
| `bench.sh`, `bench_results.md`, `bench_results.csv`, `bench_runs.csv` | the benchmark and its results (quiet machine, 20 November 2026, after the revision; the first version's is kept in `REPORT.md`, section 1) |
| `REPORT.md`, `PROOFS.md` | see above; `PROOFS.html` and `PROOFS.pdf` are `PROOFS.md` rendered (`formal/render.sh PROOFS`) |
| `formal/QTSimFast.lean` | the Lean 4 model and proofs (needs only Mathlib) |
| `formal/QTSimFastShadows.lean` | alignment checks: consequences, second formalizations, concrete instances, counterexamples for non-exact algorithms |
| `formal/*.md`, `*.html`, `*.pdf` | the two Lean files rendered as documents (`formal/render.sh`); each states the SHA-256 of the Lean file it was generated from |
| `formal/model_check.c`, `ModelCheck.lean`, `model_check_expected.txt` | the two C programs' own routines and the Lean model on one small instance: the same bit patterns (`make model-check`) |
| `formal/check.sh`, `render.sh`, `lean2md.awk` | re-check the Lean files in a Mathlib project; regenerate the documents |
| `COPYING` | the GNU General Public License, version 3 |

## What was verified, and how strongly

- **Tested:** the two programs' outputs are byte-identical in every one of 13,915
  comparisons (457 option files; 1 to 8 threads; the pipeline on, off and in each of its
  modes; block sizes; sanitizer builds), on one machine. See `PARITY_FULL.txt` and `REPORT.md`, section 7.
- **Proved, for a model:** the algorithmic changes (summation order kept while the loops
  are regrouped; shared blocks with copy on write; the lowest crossover point; drawing
  first and copying afterwards; their composition over generations) give the same results
  as the original algorithms, with floating-point addition modelled as an arbitrary
  operation with no algebraic law. See `PROOFS.md`.
- **Argued by inspection:** that the C code is what the model describes (supported by
  `make model-check` on one instance), the threading (including the pipeline, whose
  threads are not modelled; ThreadSanitizer finds no race in the runs of the matrix), the
  output formatting, and the statistics computed once.

## Not included

The paper's simulation input files, results, manual and Python scripts are not in this
package: they carry no licence statement. They are available from the paper's data
availability statement.

## Licence

GNU General Public License, version 3 or (at your option) any later version; see
`COPYING` and the notice at the top of every source file. QTSim: copyright 2023 Andrew
Meade. This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.
