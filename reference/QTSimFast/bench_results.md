# QTSimFast benchmark

Written by `bench.sh` on 2026-11-19 11:42 PST.

- Machine: Apple M3 Pro, 12 processors; Darwin arm64.
- Load averages (1, 5, 15 min) before: 2.40 3.67 5.12; after: 2.12 2.62 3.38. Runs were made one at a time.
- Compiler: Apple clang version 21.0.0 (clang-2100.4.1.1); GSL 2.8; Julia: julia version 1.12.6; commit e35f5e3.
- Protocol: the published inputs' options with `Generations` = the move's first generation + 2000; in each stratum the input with the slowest-moving optimum among the first 400; 3 seeds (1001 ...); medians over seeds.
- Each cell: wall seconds / CPU seconds (user + system), and for all but the original the speed-up in wall time against the original.
- QTSimFast's output was compared with the original's on every run; the Julia kernel has its own random numbers (other realisations of the same model, nothing compared).

| Stratum | Generations | Original C | QTSimFast, 1 thread(s) | QTSimFast, 2 thread(s) | QTSimFast, 4 thread(s) | QTSimFast, 8 thread(s) | Julia kernel, 1 thread | Output identical |
|---|---|---|---|---|---|---|---|---|
| QT0100-Possion3.0 | 12000 | 5.61 / 5.58 | 2.38 / 2.36 (2.4x) | 1.92 / 2.54 (2.9x) | 1.91 / 2.54 (2.9x) | 1.92 / 2.54 (2.9x) | 3.01 / 3.01 (1.9x) | yes |
| QT0100-Possion5.0 | 12000 | 6.49 / 6.47 | 2.65 / 2.62 (2.4x) | 2.18 / 2.81 (3.0x) | 2.18 / 2.81 (3.0x) | 2.19 / 2.81 (3.0x) | 3.56 / 3.56 (1.8x) | yes |
| QT0500-Possion3.0 | 12000 | 16.85 / 16.82 | 3.43 / 3.41 (4.9x) | 2.96 / 3.59 (5.7x) | 2.98 / 3.60 (5.7x) | 2.96 / 3.58 (5.7x) | 3.93 / 3.93 (4.3x) | yes |
| QT0500-Possion5.0 | 12000 | 18.41 / 18.39 | 3.79 / 3.77 (4.9x) | 3.25 / 3.86 (5.7x) | 3.25 / 3.87 (5.7x) | 3.24 / 3.85 (5.7x) | 4.56 / 4.56 (4.0x) | yes |
| QT1000-Possion3.0 | 12000 | 33.51 / 33.49 | 4.59 / 4.57 (7.3x) | 3.23 / 5.34 (10.4x) | 2.99 / 6.68 (11.2x) | 3.53 / 8.49 (9.5x) | 4.66 / 4.66 (7.2x) | yes |
| QT1000-Possion5.0 | 12000 | 35.84 / 35.82 | 5.03 / 5.00 (7.1x) | 4.20 / 6.25 (8.5x) | 4.49 / 7.68 (8.0x) | 4.57 / 10.50 (7.8x) | 5.24 / 5.24 (6.8x) | yes |
| QT5000-Possion3.0 | 12000 | 158.48 / 158.42 | 13.60 / 13.58 (11.7x) | 8.63 / 15.36 (18.4x) | 6.36 / 17.36 (24.9x) | 7.23 / 23.88 (21.9x) | 8.98 / 8.98 (17.6x) | yes |
| QT5000-Possion5.0 | 12000 | 155.93 / 155.87 | 12.98 / 12.95 (12.0x) | 8.58 / 14.99 (18.2x) | 6.35 / 16.86 (24.6x) | 7.11 / 21.49 (21.9x) | 9.55 / 9.55 (16.3x) | yes |
