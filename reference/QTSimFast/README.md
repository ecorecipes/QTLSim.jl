# QTSimFast: a faster QTSim with the same output, bit for bit

**This is a modified version of QTSim.** QTSim (its source calls it QTLSim 1.0) is by
Andrew Meade, University of Reading, copyright 2023, and is free software under the GNU
General Public License, version 3 or (at your option) any later version. The unmodified
source is in [`../QTSim`](../QTSim). The copy in `src/` here was **modified in November
2026 for speed**, for the QTLSim.jl project, with Claude (Anthropic); it is distributed
under the same licence (the package's `LICENSE`). Every file states at its top whether
and how it was changed.

For the same options file (and so the same seed) QTSimFast writes the same bytes to
standard output as QTSim, and exits with the same code. It is not the reference: the
original stays the reference implementation, the package's validation uses the original,
and `find_qtsim()` never returns this program. Use it when many QTSim runs are needed and
the original's output is what is wanted.

```sh
make -C reference/QTSimFast              # builds reference/QTSimFast/build/QTSimFast
make -C reference/QTSimFast check        # the files said to be unmodified are unmodified
make -C reference/QTSimFast parity       # quick parity test against reference/QTSim (about a minute)
make -C reference/QTSimFast parity-full  # the whole matrix (about 1.5 CPU-hours)
make -C reference/QTSimFast bench        # timings; writes bench_results.md, .csv, bench_runs.csv
reference/QTSimFast/build/QTSimFast options.txt > out.tsv
```

It needs what QTSim needs (a C99 compiler, GSL with `gsl-config` on the `PATH`) and, for
threads, C11 atomics and POSIX threads (`-DQTSIM_NO_THREADS` builds without either).

## What was changed

| File | Status |
|---|---|
| `GenLib.c`, `GenLib.h`, `PassOptionFile.c`, `PassOptionFile.h`, `RandLib.c`, `RandLib.h` | unmodified copies (`make check` compares them with `../QTSim/src`) |
| `QTSim.h` | five function prototypes removed, notice added |
| `QTSim.c` | options code (defaults, option file, printed header, analytical fitness width, `main`) unchanged; the simulation core rewritten; functions the original never calls removed |
| `FastRNG.h`, `FastRNG.c` | new: MT19937 inlined, and the GSL draws built on it |
| `FastThreads.h`, `FastThreads.c` | new: a small pool of POSIX threads |

The original's identifiers are kept, misspellings included, where the functions are kept
(`CaclPopPhenotype`, `SetSurviving`, `RecombinationPoissonChr`, `NewPop`, `MutatePop`,
`OuputReg`'s code in `WriteOutput`, ...).

## Where the original spends its time

Sampled with macOS `sample` (1 ms) on published inputs, 1000 diploid individuals. These
are shares of one run each on a busy machine: a guide to what to optimise, not a
benchmark.

| Share of samples | 100 loci | 1000 loci | 5000 loci |
|---|---|---|---|
| Building gametes locus by locus (`RecombinationPoisson`, a branch per locus) | 33% | 62% | 66% |
| Summing genotypes over loci (`CaclGenotype`, inlined in `PopFitness`) | 9% | 24% | 30% |
| GSL statistics for the output line (means, variances, regression, an unused correlation) | 15% | 2% | <1% |
| Random numbers (MT19937 through GSL's function pointers, Poisson, Gaussian, `exp`, `log`) | 23% | 5% | 2% |
| Sorting the crossover points (`qsort`) and allocating the crossover map (`malloc`/`free`) | 13% | 4% | 1% |
| Writing the output (`write`) | 1% | <1% | <1% |

So with many loci nearly all the time goes into copying genomes and reading them again
to sum them; with few loci it is spread over the random numbers, the per-gamete
allocation and sorting, and the statistics.

## What was optimised

1. **Genomes as shared blocks.** A chromosome is a row of block numbers; a block holds
   the alleles of 16, 32 or 64 consecutive loci (`QTSIM_BLOCK`; any size gives the same
   output). A gamete takes its blocks from one of the parent's chromosomes up to the
   crossover and from the other after it, so it copies block numbers, and only the block
   the crossover falls in is built anew. A mutation copies one block (or changes it in
   place if the chromosome alone holds it). Blocks no chromosome uses are collected
   while the genotypes are summed. This is the idea of QTLSim.jl's `blocks.jl`, without
   cached block sums (see the limits below).
2. **Genotypes summed in the original order, eight individuals side by side.** A
   genotype is `((0 + (a0+b0)) + (a1+b1)) + ...` over the loci; each addition must wait
   for the one before, but the sums of different individuals are independent, so eight
   are advanced together (two per vector register where the compiler has vector
   extensions; elementwise IEEE addition gives the scalar bits). The sums advance one
   block position at a time over the whole population, so that the few hundred blocks in
   use at a position stay in the first-level cache while every individual reads them.
3. **Random numbers.** MT19937 is inlined with `gsl_rng_uniform`, `gsl_rng_uniform_int`
   and the small-mean branch of `gsl_ran_poisson` (integer arithmetic, one division by
   2^32, one multiplication, `exp(-mean)` from the same libm, computed once). Every other
   draw (`gsl_ran_gaussian`, `gsl_ran_poisson` for means above 10) is made by GSL itself
   from the same stream, through a `gsl_rng` whose type reads the inlined state. At
   start the stream is compared with `gsl_rng_mt19937` for the run's seed over three
   regenerations of the state; the program refuses to run if they differ.
4. **Crossover maps.** The original allocates a map per gamete, sorts the points and
   then switches strand at the first point only. Here the points are drawn the same way
   (the same rejection of repeated points) and the lowest is kept; nothing is allocated
   or sorted.
5. **Two passes per generation.** The draws of a generation (survival, parents,
   crossover points, strand, environmental deviation; then mutations) are made in the
   original order; none depends on a genome. The block numbers are copied afterwards,
   position by position.
6. **Statistics computed once.** The original computes the phenotypes' mean up to three
   times and their variance twice per generation, the genotypes' variance twice, and an
   SD, the variance of the environmental deviations and a correlation that it never
   uses. Each needed value is computed once, by the same GSL function on the same
   numbers (`gsl_stats_variance_m` and `gsl_stats_sd_m` with the mean already computed
   are GSL's own definition of `gsl_stats_variance` and `gsl_stats_sd`; at start the
   program checks on test vectors that the library in use gives the same bits both ways,
   and calls the plain functions if not).
7. **Threads** (`QTSIM_THREADS`, default 1). The genotype sums and the copying of block
   numbers are divided among threads by individual and by position; nothing is reduced
   across threads, so the output is the same for every thread count. With more than one
   thread the output line (statistics and formatting, which draw no random numbers) is
   computed and written by a further thread while the next generation starts.
8. Smaller things: no allocation in the generation loop; survival and strand choice
   without a branch; `sqrt` of the variances and `exp(-mean)` taken once.

Tried and dropped: forming the first eight products of the Poisson loop ahead of the
draws to avoid its unpredictable exit (slower: the multiplications wait on each other
longer than the branch costs).

### Speed

**To be measured on a quiet machine**: `make bench` writes
[`bench_results.md`](bench_results.md) (original C, QTSimFast at 1, 2, 4 and 8 threads,
the Julia kernel on one thread; medians of three seeds on one published input per
stratum, wall and CPU seconds). With `BENCH_ABLATION=1` it adds a table with one
optimisation switched off at a time.

| Optimisation switched off | 100 loci | 1000 loci | 5000 loci | Parity held |
|---|---|---|---|---|
| Genotypes one individual at a time (`QTSIM_LANES=1`) | to be measured | to be measured | to be measured | yes |
| One block per chromosome (`QTSIM_BLOCK` = number of loci) | to be measured | to be measured | to be measured | yes |
| GSL's generator through function pointers, no vector extensions (`build/QTSimFast-plain`) | to be measured | to be measured | to be measured | yes |
| Blocks collected every 8 generations (`QTSIM_GC=8`) | to be measured | to be measured | to be measured | yes |

## How parity is kept

Three rules, kept throughout `QTSim.c`:

1. Every random draw is made in the original order by the same distribution function
   with the same arguments. Draws that GSL computes in floating point beyond a single
   multiplication are left to GSL: a compiler may contract `x*x + y*y` in
   `gsl_ran_gaussian` into a fused multiply-add (Homebrew's GSL on arm64 does), so a
   copy of that code could differ from the installed library.
2. Every floating-point value that is printed or compared is computed by the same
   operations in the same order: the genotype sum locus by locus, phenotype = genotype +
   environment, the fitness expression with `gsl_ran_gaussian_pdf` and the constant 0.4,
   and the statistics by GSL's own functions on the same arrays. The build uses
   `-ffp-contract=off` and no option that relaxes IEEE arithmetic; no library `exp` or
   `log` is replaced.
3. Block size, thread count and collection interval decide where an allele is stored
   and who computes a value, never how it is computed.

`test_parity.sh` compares standard output (`cmp`) and exit codes with the original's:

- the paper's inputs: three from each of the eight strata with three seeds each, to 2000
  generations after the move starts (5000 loci: two full-length runs per stratum, the
  others with the move brought forward), populations that die and that survive;
- the single-locus protocol (6500 generations, jump at 5000, phenotype move) and QTSim's
  demonstration input;
- a synthetic matrix: haploid and diploid; no move and each `OptimumMoveType`;
  recombination means 0, below 1, 3, 12 and far above the number of loci; 1, 2, 3, 7 and
  100 loci and numbers of loci around the block boundaries; population sizes 1 to 1001;
  populations that die at once or under the move; `Sample`, `ChangeFitnessSD`,
  `StandardisePhenotype`, `Drift`; heritability 0, 1 and above 1; no mutation and
  mutation means above 10; seeds 0 and above 2^32; option-file syntax and every error
  exit;
- each of these with 1, 2, 3, 4 and 8 threads (forced on small inputs), repeated runs,
  block sizes 1, 5, 16 and 64, sums one individual at a time, and builds with
  AddressSanitizer + UndefinedBehaviorSanitizer, with ThreadSanitizer, and without the
  inlined generator, threads and vector extensions.

### Parity results

`./test_parity.sh full`, 2026-11-19, macOS arm64 (Apple M3 Pro, Apple clang 21, GSL 2.8),
on the code as committed:

- 457 option files run by the original (66 MB of output; 77 runs end with
  `Pop all died.`), 6406 comparisons of QTSimFast with it: **6406 identical, 0
  different**.
- Published inputs: 72 option files from the eight strata (48 populations die, 24
  survive), 58 of them to 12,000 generations, four of those with 5000 loci.
- By variant: one thread 498 (with the 5000-locus thread counts 2, 3, 4, 8 among them);
  threads forced on small inputs, 2, 3, 4 and 8, 1518; repeats of 4 threads 730; output
  on the main thread 365; `auto` 365; block sizes 1100; one individual at a time 369;
  collection every 50 generations 365; the plain build 366; AddressSanitizer +
  UndefinedBehaviorSanitizer 365 and ThreadSanitizer 365, both without a report.
- The original built with AddressSanitizer + UndefinedBehaviorSanitizer ran the 365
  synthetic option files without a report, so no undefined behaviour was observed in
  it on these inputs. (Not exercised: `Sample 0`, see the limits.)

An earlier version of the code passed the same matrix (6233 of 6233) before the cases
with no loci, no individuals and no generations were added.

## Environment variables

None changes the output. Invalid values are ignored with a note on standard error.

| Variable | Meaning | Default |
|---|---|---|
| `QTSIM_THREADS` | threads for the genotype sums, or `auto` (up to 4 when PopSize x Ploidy x NoQTL is at least 4 million, else 1) | 1 |
| `QTSIM_BLOCK` | loci per block | 16 below 400 loci, 32 below 2500, else 64 |
| `QTSIM_GC` | generations between collections of unused blocks | 1 |
| `QTSIM_LANES` | 1: sum genotypes one individual at a time | 8 |
| `QTSIM_PAR_MIN` | smallest PopSize x Ploidy x NoQTL summed on several threads (provisional, until the benchmark is run) | 2000000 |
| `QTSIM_ASYNC_OUTPUT` | 0: write output on the main thread even with several threads | 1 with several threads |
| `QTSIM_STATS` | 1: block statistics on standard error at the end | 0 |

The default is one thread because sweeps run many processes at once, and because threads
pay only with many loci. With several threads the workers wait by polling briefly before
they sleep, so CPU time rises more than wall time falls.

## What bit parity rules out

- **Cached partial sums.** QTLSim.jl's kernel keeps the sum of each block and adds block
  sums; that changes the order of the additions and the last bits of a genotype, which
  in the end changes which individuals survive. Here every allele of every individual is
  added again each generation, in order. This is what limits the speed with many loci:
  2 x NoQTL x PopSize additions that cannot be skipped, each individual's in sequence.
- **A different or parallel random stream.** All draws come from one MT19937 stream in a
  fixed order, with a data-dependent number of uniforms per Gaussian and per Poisson; it
  cannot be split among threads. With few loci this serial part is most of the run.
- **Faster statistics.** GSL's running mean and variance divide at every element, each
  step waiting on the last; several could be interleaved in one loop, but that would
  mean reimplementing GSL's arithmetic (which uses `long double` and, in some builds,
  fused operations) rather than calling it.
- **Vector `exp`, fused multiply-add, reassociation**: none is used.

## Limits and notes

- The original stays the reference. QTSimFast is tested against it; nothing is proved.
- `PopSize` or `NoQTL` below 0 are refused with a message on standard error (the
  original fails on them in its own ways). `Sample 0` divides by zero in both programs
  (`Itter % Opt->Sample`), which C leaves undefined; it is not tested.
- Without a `Seed` both programs seed from the time and the process number, so two runs
  cannot be compared.
- Standard error is not part of the comparison (neither program writes to it in a normal
  run; QTSimFast uses it for the notes above).
- Memory: the block tables are allocated for the worst case (twice the original's
  genomes) but touched only as blocks are used.
- Portability: written for gcc and clang on arm64 and x86-64 (vector extensions and
  `yield`/`pause` hints behind `#if`, with scalar fallbacks). It has been run on macOS
  arm64 (Apple clang 21, GSL 2.8) only; the `QTSim parity` CI job is set up to build and
  test it on Linux x86-64 with gcc, which has not run yet.
  The two start-up checks guard the assumptions about GSL (generator, statistics).
