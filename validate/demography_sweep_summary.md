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

Mean seconds per run (on this machine, with the threads it ran on), and the CPU-hours
they imply for every demography run here (fixed size and each R) at three scales:

- NoQTL  100, u 3: 13.3 s per input over 4 demographies
- NoQTL  100, u 5: 15.3 s per input over 4 demographies
- NoQTL  500, u 3: 18.9 s per input over 4 demographies
- NoQTL  500, u 5: 21.4 s per input over 4 demographies
- NoQTL 1000, u 3: 69.1 s per input over 4 demographies
- NoQTL 1000, u 5: 29.0 s per input over 4 demographies
- NoQTL 5000, u 3: 402.9 s per input over 4 demographies
- NoQTL 5000, u 5: 274.0 s per input over 4 demographies
- scale 0.01: about 35 CPU-hours
- scale 0.05: about 176 CPU-hours
- scale 1.00: about 3516 CPU-hours
