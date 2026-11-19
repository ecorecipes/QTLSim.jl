# QTSimFast benchmark

**To be measured on a quiet machine.** No timings have been recorded yet.

Run `make -C reference/QTSimFast bench` (about half an hour, one run at a time). It
replaces this file with one table (the original C program, QTSimFast with 1, 2, 4 and 8
threads, the Julia kernel on one thread; medians of three seeds on one published input
per stratum; wall and CPU seconds; speed-ups against the original; the machine and its
load), and writes `bench_results.csv` and `bench_runs.csv` beside it.
`BENCH_ABLATION=1` adds a table with one optimisation switched off at a time.
