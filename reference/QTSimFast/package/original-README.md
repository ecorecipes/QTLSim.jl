# QTSim, the original program

`src/` holds the C source of QTSim (its source calls itself QTLSim 1.0), written by Andrew
Meade, University of Reading, copyright 2023, as published with Pagel, Gardner and Meade,
Scientific Reports (2025), doi:10.1038/s41598-025-24445-2. The files are unmodified copies
of `Code/QTSim/src/` in the paper's code and data
(http://www.evolution.reading.ac.uk/QuantitativeGenetic/Repository.zip); `SHA256SUMS`
lists their checksums, which `make check` (in the directory above) verifies.

QTSim is free software under the GNU General Public License, version 3 or (at your
option) any later version, as the header of every file in `src/` states (`../COPYING`).
This README and `SHA256SUMS` are not part of QTSim. The original has no Makefile; the
Makefile of the directory above builds it with
`cc -O2 -ffp-contract=off $(gsl-config --cflags) src/*.c $(gsl-config --libs) -lm`.
