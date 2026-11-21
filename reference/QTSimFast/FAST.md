# QTSimFast's fast modes (`QTSIM_FAST`): faster, not bit-identical

By default QTSimFast writes the same bytes as QTSim for the same options file. The
environment variable `QTSIM_FAST` turns on optimisations that give that up:

| `QTSIM_FAST` | what changes | output |
|:--|:--|:--|
| unset or `0` | nothing: the exact program | byte-identical to QTSim (the parity matrix) |
| `1` | genotypes are sums of cached block sums | QTSim's draws, in QTSim's order; genotypes differ in the last bits; the printed output is almost always still identical |
| `2` | level 1, and another generator and other sampling methods | equal to QTSim's in distribution, not run for run |

The program says so on standard error at start-up whenever `QTSIM_FAST` is set. Nothing
changes in the exact mode: every fast path is behind the flag, and the parity matrix
(`./test_parity.sh`) checks the exact mode as before.

## Level 1: genotypes from cached block sums

QTSim sums an individual's genotype locus by locus, `G = (…((a₀ + c₀) + (a₁ + c₁)) + …)`.
Bit parity forces QTSimFast to do the same, `2·L` additions per individual per
generation, in a chain each addition of which waits on the one before. With many loci
that is most of the run (73% at 5000 loci, 51% at 1000, 13% at 100, profiled).

At level 1 every block of loci (16 to 64 alleles) carries the sum of its alleles, set
when the block is made (a crossover inside it, or a mutation in it). A genotype is then
the sum over block positions of `(sum of the block on chromosome 0 + sum of the block on
chromosome 1)`: `2·L/B` additions per individual, and `B` per block made. It is the same
real number rounded differently, about 1e-16 relative.

Everything else is unchanged, including every random draw and its order. A run can depart
from QTSim's only where a rounding difference of that size changes the outcome of a
comparison (a survival draw U < W, for instance), which is rare: see the validation below.
The pipeline (two threads or more) is not used at level 1; the threaded sums are.

## Level 2: another generator and other methods

On top of level 1:

- **Generator**: xoshiro256++ (Blackman and Vigna 2018) instead of MT19937, seeded from
  the options file's `Seed` through SplitMix64.
- **Uniforms**: 53 random bits instead of 32.
- **Integers** (parents, crossover points, mutated loci): Lemire's (2019) multiply and
  reject, unbiased, with no division, instead of GSL's division and reject.
- **Poisson** (crossovers, mutations; means up to 10): inversion, one uniform per draw,
  instead of GSL's product of uniforms, one per event. Larger means: GSL's algorithm on the
  same stream.
- **Normal** (environmental deviations, mutation effects): Doornik's (2005) ZIGNOR
  ziggurat, 128 layers, instead of GSL's Box–Muller (polar) method.

Each sampler was checked on 10⁷ draws against the exact distribution (moments, tails,
Kolmogorov–Smirnov, χ²; integers on five more seeds with 5·10⁷ draws each): no departure.

- **Independent streams.** Every individual's draws in a generation come from a stream of
  its own: its survival from a key (Seed, generation, 1, individual), and all the draws of
  offspring *i* (its parents, its two gametes, its environmental deviation, then the
  mutations of its chromosomes, in QTSim's order for one offspring) from a xoshiro256++
  stream seeded from (Seed, generation, 2, *i*). So the draws are made on the threads of
  `QTSIM_THREADS`, each thread a range of offspring, and recorded; the genomes are then
  built from the records in the offspring's order (mixed blocks numbered in order, filled
  on the threads). The fitness loop runs on the threads too.
- **One-pass statistics.** The output line's means, variances, standard deviation and the
  regression of fitness on (reflected) phenotype are computed in one loop over the
  population and one over the survivors, from sums about the first value, instead of
  GSL's separate passes.

The same `Seed` gives the same run, every time, **for any number of threads**: the output
is byte-identical with 1, 2, 4 or 8 threads (checked on every example, and in the
validation below). It does not give QTSim's run.

## Validation

`validate/qtsimfast_fast.jl` in QTLSim.jl (in the hand-over package: `validation/`, with
the criteria, both logs and the plots), with criteria committed before any comparison
(`validate/qtsimfast_fast_criteria.md`): seven scenarios (the single-locus protocol;
100, 1000 and 5000 loci at equilibrium; a haploid design; a move to extinction; a move
near the survival threshold), 1000 runs per program (500 at 5000 loci), 31 summary
statistics, their empirical distributions compared by Kolmogorov–Smirnov tests (Fisher's
test for proportions), Holm-corrected. Level 1 is checked by counting byte-identical
outputs. Results (2026-11-20, `validate/qtsimfast_fast.log`, plots in
`validate/qtsimfast_fast_figs/`):

- **Level 1: PASS.** All 6500 runs byte-identical to the exact program's. Also the full
  parity matrix run with `QTSIM_FAST=1`: 13,915 of 13,915 comparisons with the original
  QTSim identical (457 option files incl. 72 published inputs, 48 of them to extinction;
  1–8 threads; sanitizer builds clean). No case was found in which level 1's output
  differs from QTSim's; it is possible in principle.
- **Level 2 (run 1, the first implementation, one stream): PASS.** No Holm rejection
  among 30 tests (smallest p 0.048).
- **Level 2 as it is now (run 2, addendum 1: independent streams, parallel draws, one-pass
  statistics; fresh seeds; 4 threads; `validate/qtsimfast_fast_run2.log`): PASS.** No Holm
  rejection among 30 tests; smallest p 0.026 (V_G before the move in the extinction
  scenario, +1.1% in the mean, Anderson–Darling p 0.010; p 0.65 for the same statistic in
  run 1); KS distances 0.017–0.084; means within about 2% wherever the mean is not near
  zero. **Thread invariance: PASS**, 700 of 700 runs byte-identical with 1 and with 4
  threads. Level 1 again 6500 of 6500 identical.
- In both runs one statistic (time of death in the single-locus protocol) was not
  testable: no run died. The tests detect distribution differences of D ≥ 0.06–0.08 with
  1000 runs per sample (0.12 at 5000 loci, 500 runs); smaller differences are not
  excluded.

## Speed

`./bench_fast.sh` times the original, the exact mode and both fast modes on `bench.sh`'s
inputs (one published input per stratum, 12,000 generations, median of three seeds);
results in `bench_fast.md`, every run in `bench_fast_runs.csv`. On a quiet machine
(2026-11-20, Apple M3 Pro, load 1.9–2.7), wall time, speed-up against the original QTSim
(the range over the two crossover means of each locus count):

| loci | exact, 1 thread | exact, 4 threads | level 1 | level 2, 1 thread | level 2, 4 threads |
|:--|:--|:--|:--|:--|:--|
| 100 | 2.6× | 4.3–4.7× | 2.8× | 3.6–3.8× | 7.0–7.6× |
| 500 | 4.9–5.1× | 11–12× | 6.3–6.7× | 8.2–8.4× | 16× |
| 1000 | 7.3–7.5× | 17× | 11–12× | 14–15× | 30× |
| 5000 | 12× | 27–28× | 28–30× | 34–36× | 80–83× |

Level 2 on four threads also uses less processor time than the exact mode on four (at
5000 loci 7.5 CPU-seconds against 25; the original 153). Exact and level-1 outputs were
identical to the original's in all 72 runs of the benchmark.

## When to use which

- To reproduce a QTSim run, or anything published from QTSim: the exact mode.
- For new simulations of QTSim's model, many runs, many loci: level 1 gives most of the
  gain at no measurable cost in fidelity.
- Level 2 gives up run-for-run agreement with QTSim for more speed, mostly from threads:
  it is the only mode in which the draws themselves run in parallel, and its output does
  not depend on the number of threads. Use it for large sweeps of new simulations where
  agreement in distribution is enough.
