# QTSimFast benchmark

Written by `bench.sh` on 2026-11-20 06:57 PST.

- Machine: Apple M3 Pro, 12 processors; Darwin arm64.
- Load averages (1, 5, 15 min) before: 3.45 3.58 4.90; after: 2.92 3.27 3.46. Runs were made one at a time.
- Compiler: Apple clang version 21.0.0 (clang-2100.4.1.1); GSL 2.8; Julia: julia version 1.12.6; commit ab95295.
- Protocol: the published inputs' options with `Generations` = the move's first generation + 2000; in each stratum the input with the slowest-moving optimum among the first 400; 3 seeds (1001 ...); medians over seeds.
- Each cell: wall seconds / CPU seconds (user + system), and for all but the original the speed-up in wall time against the original.
- QTSimFast's output was compared with the original's on every run; the Julia kernel has its own random numbers (other realisations of the same model, nothing compared).

| Stratum | Generations | Original C | QTSimFast, 1 thread(s) | QTSimFast, 2 thread(s) | QTSimFast, 4 thread(s) | QTSimFast, 8 thread(s) | Julia kernel, 1 thread | Output identical |
|---|---|---|---|---|---|---|---|---|
| QT0100-Possion3.0 | 12000 | 5.32 / 5.29 | 1.98 / 1.96 (2.7x) | 1.10 / 2.68 (4.8x) | 1.10 / 2.69 (4.8x) | 1.10 / 2.69 (4.8x) | 3.06 / 3.06 (1.7x) | yes |
| QT0100-Possion5.0 | 12000 | 6.22 / 6.21 | 2.22 / 2.20 (2.8x) | 1.35 / 3.19 (4.6x) | 1.35 / 3.19 (4.6x) | 1.34 / 3.16 (4.6x) | 3.61 / 3.61 (1.7x) | yes |
| QT0500-Possion3.0 | 12000 | 16.15 / 16.12 | 3.10 / 3.08 (5.2x) | 1.74 / 3.98 (9.3x) | 1.40 / 5.98 (11.5x) | 1.42 / 8.67 (11.4x) | 3.95 / 3.95 (4.1x) | yes |
| QT0500-Possion5.0 | 12000 | 16.94 / 16.93 | 3.32 / 3.30 (5.1x) | 1.85 / 4.20 (9.2x) | 1.60 / 6.82 (10.6x) | 1.61 / 9.83 (10.5x) | 4.55 / 4.55 (3.7x) | yes |
| QT1000-Possion3.0 | 12000 | 32.30 / 32.29 | 4.22 / 4.20 (7.7x) | 2.53 / 5.55 (12.8x) | 1.83 / 7.81 (17.7x) | 1.69 / 10.80 (19.1x) | 4.72 / 4.71 (6.8x) | yes |
| QT1000-Possion5.0 | 12000 | 33.52 / 33.49 | 4.49 / 4.47 (7.5x) | 2.64 / 5.78 (12.7x) | 1.98 / 8.41 (16.9x) | 1.79 / 11.45 (18.7x) | 5.3 / 5.3 (6.3x) | yes |
| QT5000-Possion3.0 | 12000 | 150.35 / 150.32 | 12.72 / 12.70 (11.8x) | 7.75 / 18.58 (19.4x) | 5.48 / 24.92 (27.4x) | 6.75 / 34.31 (22.3x) | 9.38 / 9.38 (16.0x) | yes |
| QT5000-Possion5.0 | 12000 | 152.44 / 152.37 | 13.27 / 13.24 (11.5x) | 7.88 / 18.89 (19.3x) | 5.69 / 25.50 (26.8x) | 6.72 / 34.59 (22.7x) | 9.83 / 9.83 (15.5x) | yes |
