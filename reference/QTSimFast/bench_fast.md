# QTSimFast fast modes: timings

- Date: 2026-11-20 16:49 PST; machine: Apple M3 Pro; load before: 2.71 2.10 2.25; after: 1.92 1.91 2.14.
- Compiler: Apple clang version 21.0.0 (clang-2100.4.1.1); GSL 2.8; commit 17ada76.
- Inputs as in bench.sh: per stratum the published input with the slowest move among the first 400, run to 2000 generations after the move; median of 3 seeds.
- Wall seconds; in brackets the speed-up in wall time against the original. Exact and level 1 outputs compared with the original's: identical in 72 of 72 runs.
- Level 2 draws other random numbers: its runs are different runs (the same statistics in distribution), so a run can last a different number of generations.

| Stratum | Original | Exact, 1 thread | Exact, 4 threads | Level 1 | Level 2, 1 thread | Level 2, 2 threads | Level 2, 4 threads |
|---|---|---|---|---|---|---|---|
| QT0100-3.0 | 5.66 | 2.20 (2.6x) | 1.20 (4.7x) | 2.05 (2.8x) | 1.57 (3.6x) | 0.96 (5.9x) | 0.81 (7.0x) |
| QT0100-5.0 | 6.50 | 2.48 (2.6x) | 1.50 (4.3x) | 2.36 (2.8x) | 1.72 (3.8x) | 1.04 (6.2x) | 0.86 (7.6x) |
| QT0500-3.0 | 16.93 | 3.35 (5.1x) | 1.41 (12.0x) | 2.52 (6.7x) | 2.01 (8.4x) | 1.27 (13.3x) | 1.05 (16.1x) |
| QT0500-5.0 | 17.80 | 3.62 (4.9x) | 1.64 (10.9x) | 2.81 (6.3x) | 2.17 (8.2x) | 1.35 (13.2x) | 1.10 (16.2x) |
| QT1000-3.0 | 33.61 | 4.46 (7.5x) | 1.94 (17.3x) | 2.81 (12.0x) | 2.31 (14.5x) | 1.37 (24.5x) | 1.12 (30.0x) |
| QT1000-5.0 | 34.56 | 4.75 (7.3x) | 2.10 (16.5x) | 3.11 (11.1x) | 2.48 (13.9x) | 1.45 (23.8x) | 1.15 (30.1x) |
| QT5000-3.0 | 153.34 | 12.85 (11.9x) | 5.54 (27.7x) | 5.05 (30.4x) | 4.25 (36.1x) | 2.67 (57.4x) | 1.85 (82.9x) |
| QT5000-5.0 | 151.30 | 13.06 (11.6x) | 5.52 (27.4x) | 5.36 (28.2x) | 4.44 (34.1x) | 2.76 (54.8x) | 1.89 (80.1x) |

CPU seconds (user + system), medians:

| Stratum | Original | Exact, 1 | Exact, 4 | Level 1 | Level 2, 1 | Level 2, 2 | Level 2, 4 |
|---|---|---|---|---|---|---|---|
| QT0100-3.0 | 5.65 | 2.17 | 2.86 | 2.03 | 1.55 | 2.02 | 3.29 |
| QT0100-5.0 | 6.49 | 2.47 | 3.46 | 2.33 | 1.70 | 2.20 | 3.44 |
| QT0500-3.0 | 16.91 | 3.33 | 6.07 | 2.49 | 1.99 | 2.65 | 4.26 |
| QT0500-5.0 | 17.79 | 3.60 | 7.00 | 2.79 | 2.15 | 2.83 | 4.44 |
| QT1000-3.0 | 33.59 | 4.43 | 8.02 | 2.78 | 2.28 | 2.84 | 4.36 |
| QT1000-5.0 | 34.54 | 4.72 | 8.69 | 3.08 | 2.46 | 3.03 | 4.48 |
| QT5000-3.0 | 153.29 | 12.82 | 25.02 | 5.02 | 4.23 | 5.48 | 7.46 |
| QT5000-5.0 | 151.28 | 13.03 | 24.85 | 5.34 | 4.41 | 5.66 | 7.63 |
