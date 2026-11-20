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

The same `Seed` gives the same run, every time, on the same machine; it does not give
QTSim's run.

## Validation

`validate/qtsimfast_fast.jl` in QTLSim.jl, with criteria committed before any comparison
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
- **Level 2: PASS.** No Holm rejection among 30 tests (smallest p 0.048; KS distances
  0.017–0.074, means within about 2% wherever the mean is not near zero). One statistic
  (time of death in the single-locus protocol) was not testable: no run died. The test
  detects distribution differences of D ≥ 0.06–0.08 with 1000 runs per sample (0.12 at
  5000 loci, 500 runs); smaller differences are not excluded.

## Speed

See `bench_results.md` (the `QTSIM_FAST` columns), measured on a quiet machine.

## When to use which

- To reproduce a QTSim run, or anything published from QTSim: the exact mode.
- For new simulations of QTSim's model, many runs, many loci: level 1 gives most of the
  gain at no measurable cost in fidelity.
- Level 2 adds little speed on top of level 1 (the random numbers are a small part of a
  run once the genotypes are cheap) and gives up run-for-run agreement; use it only where
  that does not matter.
