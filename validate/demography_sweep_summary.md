# Density regulation in the published sweep

1200 inputs, scale 0.01 of each stratum, 2000 generations after the move, seeds 2026110700000 + row. Strata weighted equally.

## Each demography

| | fixed (QTSim) | R = 5.0 | R = 2.0 | R = 1.5 | published (QTSim) |
|:--|--:|--:|--:|--:|--:|
| population size before the move (median) | 1000 | 985 | 977 | 969 | 1000 |
| died before the move (of 1200) | 0 | 2 | 5 | 15 | 0 |
| alive 2000 generations after the move | 8.2% | 6.7% | 5.9% | 5.0% | 8.2% |
| generations to extinction (median) | 134 | 79 | 58 | 50 | 135 |
| generations to 50% mortality (median) | 45 | 45 | 48 | 45 | – |
| rate of change to extinction (median) | 0.054 | 0.049 | 0.044 | 0.039 | – |
| standardised lag (median) | 0.345 | 0.376 | 0.419 | 0.450 | – |
| rate of change, populations alive at 2000 (median) | 0.027 | 0.021 | 0.020 | 0.016 | – |
| standardised rate at which half are alive (logistic x50) | 0.0734 | 0.0584 | 0.0599 | 0.0533 | – |

## Paired by input, against the same input at a fixed size

| | R = 5.0 | R = 2.0 | R = 1.5 |
|:--|--:|--:|--:|
| time to extinction / fixed (median, both died after the move) | 0.60 | 0.43 | 0.37 |
| alive at a fixed size, dead here | 19 | 28 | 39 |
| dead at a fixed size, alive here | 0 | 0 | 0 |

(The fixed-size column `published (QTSim)` uses the table's own runs of the same inputs, censored at 2000.)

## Cost

The seconds per run this sweep recorded were measured while the kernel was changing
(partly on the dense kernel, before block genomes) and with other heavy jobs on the
machine, and overstated the cost several times over.
Re-estimated (2026-11-07) from the generations each of these runs lasted and the
single-thread time per generation of the current kernel on the paper's inputs
(`validate/benchmark.jl`: 0.25, 0.33, 0.39 and 0.78 ms per generation at 100, 500, 1000
and 5000 loci; 12-core Apple M3 Pro, Julia 1.12), single-thread CPU time:

- NoQTL  100: 10 s per input over 4 demographies
- NoQTL  500: 13 s per input over 4 demographies
- NoQTL 1000: 16 s per input over 4 demographies
- NoQTL 5000: 32 s per input over 4 demographies
- scale 0.01: about 6 CPU-hours
- scale 0.05: about 30 CPU-hours
- scale 1.00: about 600 CPU-hours, 270 of them in the 5000-locus strata

The time per generation was measured at u = 3; u matters little, since a gamete has at
most one crossover. Several threads cost more CPU time per run than one: a 58-input slice
of `validate/published_sweep.jl` (runs of the same kind) took 733 s of CPU time on 4
threads, with other jobs running, against 360 s estimated from the single-thread times
(the 733 s include compilation and six QTSim runs). So scale 0.01 takes roughly an hour
on 10 threads (an estimate, not measured).
