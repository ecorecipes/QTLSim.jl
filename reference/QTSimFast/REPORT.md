# QTSimFast: a faster QTSim with the same output

Technical report, 19 November 2026; revised the same day for the changes of sections
3.9 to 3.13 (the summary, section 5 and section 7 say which results are from before and
which from after).

QTSim (its source calls itself QTLSim 1.0) is by Andrew Meade, School of Biological
Sciences, University of Reading, copyright 2023, and is free software under the GNU General
Public License, version 3 or any later version. **QTSimFast is a modified version of
QTSim**, prepared in November 2026 with QTLSim.jl (https://github.com/ecorecipes/QTLSim.jl)
and distributed under the same licence. The original is included unmodified beside it.

The original program was written to be clear and correct, and it is both: this work
changed how its computation is organised, not what it computes, and the original remains
the reference against which the fast copy is tested.

Paths below are those of the hand-over package (`original/src`, `fast/src`); in the
QTLSim.jl repository the same files are `reference/QTSim/src` and
`reference/QTSimFast/src`. Line numbers of the original refer to `original/src/QTSim.c`
(SHA-256 `0416dd0a…154aa`, as published).

## 1. Summary

QTSimFast takes the same options file as QTSim and writes the same standard output.

**The guarantee that was tested.** For the same options file (and therefore the same
seed), the bytes QTSimFast writes to standard output and its exit code are those of the
original, for every thread count, block size and block-collection interval, with the
pipeline on or off. This was found to hold in all 13,915 comparisons of the test matrix
(section 7; 6,406 for the first version of the code) and in every run of the benchmark
and of the timings. It is an empirical statement about the two programs as built on one
machine (macOS, arm64). Section 7 also describes machine-checked proofs that the
*algorithms* agree; those are about a mathematical model, not about the C source.

**Speed.** From `bench_results.md`, measured on a quiet machine on 20 November 2026,
after the revision of 19 November (sections 3.9 to 3.13):

- Machine: Apple M3 Pro, 12 processors (6 performance cores); Darwin arm64. Load
  averages (1, 5, 15 min) before: 3.45 3.58 4.90; after: 2.92 3.27 3.46 (desktop
  processes only). Runs were made one at a time.
- Compiler: Apple clang 21.0.0; GSL 2.8; Julia 1.12.6; commit ab95295.
- Protocol: the published inputs' options with `Generations` = the move's first
  generation + 2000 (12,000 generations); in each stratum the input with the
  slowest-moving optimum among the first 400; 3 seeds; medians over seeds.
- Each cell: wall seconds / CPU seconds (user + system), and the speed-up in wall time
  against the original.
- QTSimFast's output was compared with the original's on every run (96 runs, all
  identical). The Julia kernel of QTLSim.jl is shown for scale only: it has its own random
  numbers and is not bit-identical to QTSim.

| Stratum | Original C | QTSimFast, 1 thread | 2 threads | 4 threads | 8 threads | Julia kernel, 1 thread | Output identical |
|---|---|---|---|---|---|---|---|
| QT0100-Possion3.0 | 5.32 / 5.29 | 1.98 / 1.96 (2.7x) | 1.10 / 2.68 (4.8x) | 1.10 / 2.69 (4.8x) | 1.10 / 2.69 (4.8x) | 3.06 / 3.06 (1.7x) | yes |
| QT0100-Possion5.0 | 6.22 / 6.21 | 2.22 / 2.20 (2.8x) | 1.35 / 3.19 (4.6x) | 1.35 / 3.19 (4.6x) | 1.34 / 3.16 (4.6x) | 3.61 / 3.61 (1.7x) | yes |
| QT0500-Possion3.0 | 16.15 / 16.12 | 3.10 / 3.08 (5.2x) | 1.74 / 3.98 (9.3x) | 1.40 / 5.98 (11.5x) | 1.42 / 8.67 (11.4x) | 3.95 / 3.95 (4.1x) | yes |
| QT0500-Possion5.0 | 16.94 / 16.93 | 3.32 / 3.30 (5.1x) | 1.85 / 4.20 (9.2x) | 1.60 / 6.82 (10.6x) | 1.61 / 9.83 (10.5x) | 4.55 / 4.55 (3.7x) | yes |
| QT1000-Possion3.0 | 32.30 / 32.29 | 4.22 / 4.20 (7.7x) | 2.53 / 5.55 (12.8x) | 1.83 / 7.81 (17.7x) | 1.69 / 10.80 (19.1x) | 4.72 / 4.71 (6.8x) | yes |
| QT1000-Possion5.0 | 33.52 / 33.49 | 4.49 / 4.47 (7.5x) | 2.64 / 5.78 (12.7x) | 1.98 / 8.41 (16.9x) | 1.79 / 11.45 (18.7x) | 5.30 / 5.30 (6.3x) | yes |
| QT5000-Possion3.0 | 150.35 / 150.32 | 12.72 / 12.70 (11.8x) | 7.75 / 18.58 (19.4x) | 5.48 / 24.92 (27.4x) | 6.75 / 34.31 (22.3x) | 9.38 / 9.38 (16.0x) | yes |
| QT5000-Possion5.0 | 152.44 / 152.37 | 13.27 / 13.24 (11.5x) | 7.88 / 18.89 (19.3x) | 5.69 / 25.50 (26.8x) | 6.72 / 34.59 (22.7x) | 9.83 / 9.83 (15.5x) | yes |

On one thread QTSimFast is 2.7 to 2.8 times as fast as the original at 100 loci, 5.1 to
5.2 at 500, 7.5 to 7.7 at 1000 and 11.5 to 11.8 at 5000. With two threads (the pipeline)
4.6 to 4.8, 9.2 to 9.3, 12.7 to 12.8 and 19.3 to 19.4 times; with four, up to 11.5, 17.7
and 27.4 times at 500, 1000 and 5000 loci. At 100 loci two threads are as fast as more.
Eight threads gain nothing over four at 100 and 500 loci, a little at 1000 (1.69 against
1.83 s, 1.79 against 1.98), and lose at 5000 (6.75 against 5.48 s). Threads cost
processor time (they poll while they wait: at 5000 loci four threads use 25 CPU seconds
for a run that one thread does in 13), so for a sweep of many runs one thread per run is
the efficient setting, and it is the default. The Julia kernel on one thread is faster
than QTSimFast on one thread at 5000 loci (9.4 to 9.8 against 12.7 to 13.3 s; it does not
have to be bit-exact), slower than QTSimFast on two threads at every size.

The same benchmark of the first version of the code (19 November 2026, commit e35f5e3,
load averages before 2.40 3.67 5.12, after 2.12 2.62 3.38), kept for comparison:

| Stratum | Original C | QTSimFast, 1 thread | 2 threads | 4 threads | 8 threads | Julia kernel, 1 thread | Output identical |
|---|---|---|---|---|---|---|---|
| QT0100-Possion3.0 | 5.61 / 5.58 | 2.38 / 2.36 (2.4x) | 1.92 / 2.54 (2.9x) | 1.91 / 2.54 (2.9x) | 1.92 / 2.54 (2.9x) | 3.01 / 3.01 (1.9x) | yes |
| QT0100-Possion5.0 | 6.49 / 6.47 | 2.65 / 2.62 (2.4x) | 2.18 / 2.81 (3.0x) | 2.18 / 2.81 (3.0x) | 2.19 / 2.81 (3.0x) | 3.56 / 3.56 (1.8x) | yes |
| QT0500-Possion3.0 | 16.85 / 16.82 | 3.43 / 3.41 (4.9x) | 2.96 / 3.59 (5.7x) | 2.98 / 3.60 (5.7x) | 2.96 / 3.58 (5.7x) | 3.93 / 3.93 (4.3x) | yes |
| QT0500-Possion5.0 | 18.41 / 18.39 | 3.79 / 3.77 (4.9x) | 3.25 / 3.86 (5.7x) | 3.25 / 3.87 (5.7x) | 3.24 / 3.85 (5.7x) | 4.56 / 4.56 (4.0x) | yes |
| QT1000-Possion3.0 | 33.51 / 33.49 | 4.59 / 4.57 (7.3x) | 3.23 / 5.34 (10.4x) | 2.99 / 6.68 (11.2x) | 3.53 / 8.49 (9.5x) | 4.66 / 4.66 (7.2x) | yes |
| QT1000-Possion5.0 | 35.84 / 35.82 | 5.03 / 5.00 (7.1x) | 4.20 / 6.25 (8.5x) | 4.49 / 7.68 (8.0x) | 4.57 / 10.50 (7.8x) | 5.24 / 5.24 (6.8x) | yes |
| QT5000-Possion3.0 | 158.48 / 158.42 | 13.60 / 13.58 (11.7x) | 8.63 / 15.36 (18.4x) | 6.36 / 17.36 (24.9x) | 7.23 / 23.88 (21.9x) | 8.98 / 8.98 (17.6x) | yes |
| QT5000-Possion5.0 | 155.93 / 155.87 | 12.98 / 12.95 (12.0x) | 8.58 / 14.99 (18.2x) | 6.35 / 16.86 (24.6x) | 7.11 / 21.49 (21.9x) | 9.55 / 9.55 (16.3x) | yes |


**Before and after the revision, side by side** (sections 3.9 to 3.13: the means of the
output line, block collection, and the pipeline, in which one thread makes the draws and
another builds the genomes). These timings were taken before the benchmark above could be
repeated, on a busy machine, and are **indicative only**; they add three threads and the
pipeline switched off. The program before the revision
and after it were run alternately on one synthetic input of the published design per
number of loci (1000 diploid individuals, 12,000 generations; minimum of three runs;
load averages 5.5 to 10 during the runs; every output identical to the original's):

| Loci | Threads | Before: wall / CPU s | After: wall / CPU s | Wall, before / after |
|---|---|---|---|---|
| 100 | 1 | 2.39 / 2.37 | 2.05 / 2.04 | 1.17 |
| 100 | 2 | 1.99 / 2.53 | 1.09 / 2.60 | 1.83 |
| 100 | 4 | 1.99 / 2.52 | 1.09 / 2.63 | 1.83 |
| 100 | 8 | 1.94 / 2.50 | 1.09 / 2.62 | 1.78 |
| 500 | 1 | 3.48 / 3.45 | 3.24 / 3.22 | 1.07 |
| 500 | 2 | 3.03 / 3.60 | 1.75 / 3.95 | 1.73 |
| 500 | 4 | 3.07 / 3.61 | 1.36 / 6.00 | 2.26 |
| 500 | 8 | 3.04 / 3.61 | 1.57 / 9.77 | 1.94 |
| 1000 | 1 | 4.67 / 4.65 | 4.41 / 4.39 | 1.06 |
| 1000 | 2 | 3.32 / 5.47 | 2.55 / 5.53 | 1.30 |
| 1000 | 4 | 3.05 / 6.67 | 1.72 / 7.30 | 1.77 |
| 1000 | 8 | 3.47 / 8.32 | 1.98 / 12.32 | 1.75 |
| 5000 | 1 | 13.61 / 13.57 | 13.93 / 13.91 | 0.98 |
| 5000 | 2 | 8.86 / 15.62 | 8.01 / 19.10 | 1.11 |
| 5000 | 4 | 6.36 / 17.17 | 5.64 / 25.43 | 1.13 |
| 5000 | 8 | 6.96 / 22.77 | 7.44 / 36.42 | 0.94 |

With three threads the revised program took 1.09, 1.49, 2.01 and 6.50 s. With four
threads and the pipeline switched off (`QTSIM_PIPELINE=0`) it took 1.70, 2.92, 2.89 and
6.37 s, so the gain with several threads is the pipeline's. What the table also shows,
and should not be overlooked:

- **Processor time rises** with several threads (at 5000 loci and four threads 25 CPU
  seconds against 17 before): the pipeline's threads poll while they wait. For a sweep
  of many runs one thread per run remains the efficient setting, and the default.
- **Eight threads are slower than four** at 500 loci and above (the machine has six
  performance cores), and at 5000 loci slower than the earlier program with eight (7.44
  against 6.96 s). Use at most four, which is what `auto` gives.
- **On one thread at 5000 loci nothing was gained** (13.93 against 13.61 s here, 12.46
  against 12.65 s in another comparison the same day: no difference beyond the noise).
- At 5000 loci the pipeline gains about a tenth: the genotype sums, which it cannot start
  before the mutations there, are nearly all of the work.

## 2. Where the original spends its time

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

With many loci nearly all the time goes into copying genomes and reading them again to
sum them. With few loci it is spread over the random numbers, the per-gamete allocation
and sorting, and the statistics.

## 3. The optimisations

Three rules keep the output unchanged, and each optimisation below is justified by them:

- **R1, random numbers.** Every random draw is made in the original order, by the same
  distribution function, with the same arguments.
- **R2, arithmetic.** Every floating-point value that is printed or compared is computed
  by the same operations in the same order as in the original.
- **R3, output.** The same `printf` conversions are applied to the same values in the same
  order.

Block size, thread count and collection interval decide where an allele is stored and
which thread computes a value; they never decide how a value is computed.

Contributions were measured once with one optimisation switched off at a time (one
thread, one seed, three strata). **These numbers are indicative: the machine was shared
with another job when they were taken**; the table of section 1 is the benchmark.

| Stratum (12,000 generations) | Original | QTSimFast | genotypes one individual at a time (`QTSIM_LANES=1`) | one block per chromosome (`QTSIM_BLOCK` = loci) | collect every 8 generations (`QTSIM_GC=8`) | GSL's generator through function pointers, no vector extensions (`QTSimFast-plain`) |
|---|---|---|---|---|---|---|
| QT0100-Possion3.0 | 5.79 | 2.44 | 2.52 (1.0x) | 2.70 (1.1x) | 2.35 (1.0x) | 3.05 (1.2x) |
| QT1000-Possion3.0 | 36.38 | 4.79 | 5.77 (1.2x) | 17.02 (3.6x) | 4.54 (0.9x) | 5.52 (1.2x) |
| QT5000-Possion3.0 | 168.02 | 13.89 | 25.18 (1.8x) | 68.65 (4.9x) | 21.01 (1.5x) | 26.17 (1.9x) |

CPU seconds, and the ratio to QTSimFast as built. Every one of these variants also gave
the original's output. Optimisations 3.3 (crossover maps), 3.5 (two passes) and 3.6
(statistics) have no switch and were not measured separately; the gap between the
original and the "one block per chromosome" column (5.79 to 2.70, 36.4 to 17.0, 168 to
68.7) is an upper bound on what they and the generator contribute together.

### 3.1 Genomes as shared blocks

*Original.* An individual holds `Ploidy × NoQTL` doubles (`CreateIndividual`, lines
239–259). A gamete is copied locus by locus:

```c
PointIndex = 0;
for(Index=0;Index<Opt->NoQTL;Index++)
{
    Chr[Index] = CP[Index];

    if(Index == CombMap[PointIndex])
    {
        if(CP == P1) CP = P2; else CP = P1;
    }
}
```

(`RecombinationPoissonChr`, lines 682–693). With 1000 individuals and 5000 loci that is
ten million copies and ten million comparisons per generation.

*Fast copy.* A chromosome is a row of block numbers; a block holds the alleles of 16, 32
or 64 consecutive loci and may be held by many chromosomes (`SIM`, `fast/src/QTSim.c`).
A gamete takes its block numbers from one of the parent's chromosomes up to the crossover
and from the other after it (`FillGametes`); only the block the crossover falls in is
built anew, from the two parental blocks (`MixBlocks`), and not even that if the two are
the same block. A mutation changes a block in place if this chromosome alone holds it (it
was made this generation for this chromosome: `Owner`), and otherwise copies the block
first (`MutatePop`). Blocks that no chromosome of the current population holds are found
while the genotypes are summed (`Mark`) and returned to a free list (`CollectBlocks`).

*Why the output cannot change.* No arithmetic is involved (R2 holds trivially): alleles
are copied, never recomputed, and a mutation is the same `+=` of the same Gaussian to the
same allele. What has to hold is that every chromosome, read through its block numbers,
is the array the original holds. That needs a block never to be written while another
chromosome holds it, and never to be handed out while in use; section 7.3 states the
invariant and the proof for a model of this code (`fillMix_spec`, `mutStep_sim`,
`collect_inv`).

*Contribution.* The largest: without sharing below the chromosome the run takes 3.6 times
as long at 1000 loci and 4.9 times at 5000 (indicative).

### 3.2 Genotypes summed in the original order, individuals side by side

*Original.* `CaclGenotype` (lines 323–351), for a diploid:

```c
Genotype = 0;
for(Index=0;Index<Opt->NoQTL;Index++)
    Genotype += Ind->Genome[0][Index] + Ind->Genome[1][Index];
```

that is `((0 + (a0+b0)) + (a1+b1)) + …`, one individual after another.

*Fast copy.* `SumGenotypesW` keeps the partial sum of every individual in `Geno[]` and
advances all of them over the blocks at position 0, then over the blocks at position 1,
and so on, so that the few hundred blocks in use at a position stay in the first-level
cache while every individual reads them. Inside a block, eight individuals are advanced
together (their sums are independent chains of additions, so the processor can overlap
them), two per vector register where the compiler has vector extensions.

*Why the output cannot change.* Each individual's sum receives the same terms `a_l + b_l`
in the same order, each computed by one IEEE addition, and is advanced by one IEEE
addition per locus (R2). An elementwise vector addition is the scalar addition in each
lane. Regrouping the *loop* across individuals does not regroup any *sum*. Proved for the
model for an arbitrary, lawless `add` (`foldl_pointwise`, `fastGeno_eq`).

*Contribution.* Summing one individual at a time costs about 1.0x, 1.2x and 1.8x at 100,
1000 and 5000 loci (indicative).

### 3.3 Crossover maps

*Original.* `CreateCombMap` (lines 608–636) allocates a map per gamete, draws
`gsl_ran_poisson(rng, RecombPoissonMean)` points, each by
`gsl_rng_uniform_int(rng, NoQTL-1)` redrawn while it is already in the map, sorts them
with `qsort`, and returns the map; the loop quoted in 3.1 compares the locus with
`CombMap[PointIndex]`, and `PointIndex` stays 0, so the gamete switches strand once, after
the lowest point.

*Fast copy.* `LowestCrossover` makes the same draws with the same redraw test (a scan of
the points drawn so far; a table of stamps when there are more than 16) and keeps the
running minimum. Nothing is allocated or sorted.

*Why the output cannot change.* The draws are the same calls in the same order with the
same arguments and the same redraw decisions (R1). The only thing the original's loop
reads from the sorted map is its first entry, which is the minimum of the points. Proved
for the model (`origLoop_eq`, `runMin_eq_head`, `fastLow_eq`).

QTSimFast reproduces this behaviour deliberately: the published results were produced
with it, and reproducing them is the point.

### 3.4 Random numbers

*Original.* Every number comes from `gsl_rng_mt19937` through GSL's function pointers.

*Fast copy.* `FastRNG.h` inlines MT19937 (seeded as GSL's `mt_set`), `gsl_rng_uniform`
(one division by 2^32), `gsl_rng_uniform_int` (integer arithmetic) and the small-mean
branch of `gsl_ran_poisson` (a product of uniforms compared with `exp(-mean)`, the
exponential taken once from the same C library). Every other draw (`gsl_ran_gaussian`,
`gsl_ran_poisson` for means above 10) is made by GSL itself, through a `gsl_rng` whose
type reads the inlined state, so the stream is one stream.

*Why the output cannot change.* The same generator in the same state gives the same
numbers; the three inlined functions use only operations that have one IEEE result (R1).
Draws whose GSL code does more floating-point work are left to GSL, because a compiler
may contract `x*x + y*y` into a fused multiply-add (Homebrew's GSL on arm64 does) and a
copy of that code could then differ from the installed library. This part is tested, not
proved: at start the program compares the inlined stream with GSL's (section 6).

*Contribution.* Together with the vector extensions, about 1.2x at 100 and 1000 loci
(`QTSimFast-plain`, indicative).

### 3.5 Two passes per generation

*Original.* `NewPop` (lines 755–788) handles one offspring at a time: parents, the first
gamete's map, strand and copy, the second gamete's, the environmental deviation.

*Fast copy.* `NewPop` first makes all the draws for all offspring, in that order,
recording for each gamete the two source chromosomes and the crossover; then
`FillGametes` copies the block numbers position by position (possibly on several
threads), and `MixBlocks` builds the blocks in which a crossover fell.

*Why the output cannot change.* The order of the draws is the original's (R1), provided
no draw depends on a genome that the original would already have written. None does: the
arguments of the draws are `NoSPop`, `RecombPoissonMean`, `NoQTL-1` and `EnvVar`; the
redraw tests compare survivor indices and points; the genomes read are the parents', which
`NewPop` does not write. This is the hypothesis `DrawsIgnoreData` of the proof
(`interleaved_eq_two_pass`); it was checked by reading the code (lines 225–236, 479–492,
591–636, 667–697, 707–740, 755–788), with `AssortativeMating` fixed at −1 as the options
parser leaves it.

### 3.6 Statistics computed once

*Original.* Per generation the phenotypes' mean is computed up to three times
(`SetFitnessMean`, twice in `Output`), their variance twice, the genotypes' variance
twice; `SetFitnessMean` computes an SD it does not use, `Output` the variance of the
environmental deviations, and `OuputReg` a correlation, neither of which is printed.

*Fast copy.* Each value that is used is computed once, by the same GSL function on the
same numbers in the same order; the unused ones are not computed. Where the mean is
already known, `gsl_stats_variance_m` and `gsl_stats_sd_m` are called with it; these are
GSL's own definitions of `gsl_stats_variance` and `gsl_stats_sd` (section 6 for the
check).

*Why the output cannot change.* A function without side effects gives the same value each
time it is called on the same data (R2), and removing a computation whose result is
unused and which draws no random number changes nothing (R1). No theorem was written for
this: it is checked by inspection, and by the parity tests.

*Contribution.* Mostly at few loci, where the statistics were 15% of the original's time;
not measured separately.

### 3.7 Threads

See section 5.

### 3.8 Smaller things

No allocation inside the generation loop (with the pipeline the list of a generation's
mutations grows when a generation has more than any before); survival and strand choice
without a branch (`S->NoSurv += (u < Fit[Index])`); `sqrt` of the variances and
`exp(-mean)` taken once (same function, same argument, same result).

Sections 3.9 to 3.13 were added on 19 November 2026, after the benchmark of section 1.

### 3.9 The three means of the output line in one loop

*Before.* A line of output needs the means of the phenotypes, of the genotypes and of
the survivors' phenotypes (inside `gsl_stats_variance` and `gsl_stats_sd`, which compute
the mean first). Each is GSL's running mean,

```c
long double mean = 0;
for (i = 0; i < size; i++)
    mean += (data[i * stride] - mean) / (i + 1);
```

a chain of `size` divisions in which each waits for the one before: with 1000 individuals
about three times 1000 division latencies per line, a fifth of the fast program's time at
100 loci (section 3.13).

*Now.* `ThreeMeans` advances the three chains side by side in one loop, each written as
above (a `long double` accumulator, a `size_t` index), so that the processor overlaps
their divisions; the variances and the SD are then GSL's own `gsl_stats_variance_m` and
`gsl_stats_sd_m` about those means.

*Why the output cannot change.* A step of the recurrence is one subtraction, one
conversion of the index, one division and one addition: there is no product next to an
addition, so there is nothing a compiler could contract into a fused operation, in this
copy or in the library's, and each of the operations has one IEEE result (R2). The three
chains do not share a variable. `gsl_stats_variance(x)` is, in GSL,
`gsl_stats_variance_m(x, gsl_stats_mean(x))`, which was already relied on (3.6). Both
steps are tested against the library in use at start-up, on vectors of 1 to 1000 numbers
(section 6); if either differed the plain GSL calls would be made, as
`QTSIM_PLAIN_STATS=1` forces (a variant of the parity matrix). The variance itself and
the regression are still left to GSL: their code has products and sums that a compiler
may fuse.

### 3.10 Collecting blocks: less often, without branches

- **Interval.** Unused blocks were collected in every generation. They are now collected
  every fourth generation below 2500 loci (`QTSIM_GC` overrides). Measured on one thread
  (minimum of three runs, load 2.5 to 3.9, so indicative): every 1 / 2 / 4 / 8 generations
  1.90 / 1.85 / 1.84 / 1.84 s at 100 loci, 2.94 / 2.89 / 2.87 / 2.88 s at 500, 4.07 /
  3.92 / 3.96 / 3.86 s at 1000, and 12.46 / 12.40 / 12.40 s (1 / 2 / 4) at 5000. At 5000
  loci there is nothing to gain, and an interval of 8 was 1.5 times slower in the
  ablation of section 3 (blocks left uncollected spread the alleles in use over more
  memory), so the interval stays 1 from 2500 loci.
- **The free list** is built without a branch: the highest block in use is found first;
  then every block number below it is written to the list and the list's length advances
  only past the unmarked ones. The list is the same list in the same order as before.
- **The block in which a crossover falls** is copied by one loop with a select per allele
  (`ND[I] = (I <= Off) ? CD[I] : OD[I]`) instead of two loops whose lengths depend on the
  crossover.

*Why the output cannot change.* Which blocks are free, and when, decides only where an
allele is stored (rule 3; proved for the model for every schedule of collections:
`run_refines`, `collect_inv`). The mixed block has the same alleles: a select copies one
of two doubles and computes nothing.

### 3.11 The pipeline: draws on one thread, genomes on another

*Before.* With few loci most of a run is the serial stream of random numbers, so threads
helped only through the asynchronous output line (section 1's old table: 2.9x at 100
loci with any number of threads).

*Now* (with `QTSIM_THREADS` of 2 or more; `QTSIM_PIPELINE=0` turns it off; one thread is
the serial program unchanged). After `SetSurviving`, every remaining draw of a generation
depends on the survivors and the options and never on a genome: the parents
(`gsl_rng_uniform_int(rng, NoSPop)`), the crossover number and points
(`RecombPoissonMean`, `NoQTL-1`), the strand, the environmental deviation (`EnvVar`), and
then, in `MutatePop`, the number of mutations of each chromosome (`NoQTL ×
MutationRatePerGamete`), each mutation's locus (`NoQTL`) and its effect (a variance
computed from options). So:

- **the main thread** makes all these draws in the original order and only records them:
  per gamete the two source chromosomes, the number of leading blocks from the first, and
  the crossover's place inside a block if it falls inside one (`GameteDraws`,
  `NewPopDraws`); per mutation the chromosome, the locus and the effect (`MutateDraws`);
- **the genome thread** follows it (`PipeFollowJob`, `PipeGeneration`): it copies the
  block numbers of the gametes whose records are written and builds the blocks in which a
  crossover fell (`FillBatch`), then, when the mutations are recorded, applies them in
  the recorded order with the same copy-on-write rule as `MutatePop`
  (`ApplyMutations`);
- **the genotype sums** are started before the mutations are known, when an individual
  expects at most 0.7 mutations per generation (`QTSIM_SPEC`): the genome thread and, for
  populations large enough, the pool's threads claim batches of 64 built individuals and
  sum them while the main thread is still drawing. When the mutations arrive, the
  individuals a mutation changed are summed again from the start and those not yet
  reached are summed for the first time; the main thread, which has finished its draws,
  takes a share of that list if it is long enough. With more mutations (5000 loci in the
  published design: 63% of individuals mutate in a generation) the sums are made after
  the mutations, on the pool, as before; then the pipeline moves only the gametes and
  the mutations off the main thread;
- the genotypes of the next population go to an array of their own, exchanged with the
  current one at the end of the generation, because the output line still needs the
  current ones.

*Why the output cannot change.*

- R1: the main thread makes the same calls in the same order with the same arguments as
  the serial program; recording a draw instead of using it at once does not change the
  generator. This needs that no draw depends on a genome, which is the hypothesis
  `DrawsIgnoreData` of the two-pass theorem (3.5, section 7.3), now also for the draws of
  `MutatePop`: its arguments are options, and its loop bounds are the population size and
  the Poisson counts just drawn.
- Genomes: every gamete is built from its own record exactly as `FillGametes` and
  `MixBlocks` build it (the same block numbers; the mixed blocks allocated in the same
  gamete order), and every mutation is the same `+=` of the same double to the same
  allele, applied in the original order, so two mutations of one allele add in the same
  order. All allocation of blocks is done by one thread at a time.
- R2: a genotype is the in-order sum over the loci of the genome it is read from. A sum
  made before a mutation changed the individual is discarded and made again in full from
  the final genome; no partial sum is reused. Which thread sums an individual, and in
  which batch, does not enter its sum (3.2).
- Block collection: the genome thread asks whether to collect at the same point of the
  generation as the serial program (when the new population is the only one alive) and
  marks exactly the blocks that population holds, by one pass over its table; so the
  free list, and with it every block number, is the serial program's (rule 3 does not
  require that; `QTSIM_STATS=1` shows the same counts of blocks at the end of a run in
  every mode).

*Synchronisation.* The two threads share two numbers. `Prog`, set by the main thread, is
`generation × (N + 2) + k`: the records of individuals `0 … k−1` are written, and `k = N
+ 1` says that the mutations are recorded. `Done`, set by the genome thread, says that
the generation is finished (and, one step before, that the main thread may sum its
share). Within the genome thread's job two more numbers hand out work: how many
individuals are built, and how many are claimed for summing (an atomic add, so no two
threads get the same batch). A thread reads only what another wrote before the number it
has read; the main thread writes no genome and reads none between its first `Prog` of a
generation and `Done`. The numbers are sequentially consistent atomics (`FastThreads.c`);
a waiting thread polls a bounded number of times and then sleeps on a condition variable
until the number is set, so no signal is lost when the machine is busy. There is no lock
per batch. The polling budgets are set in microseconds (the cost of a poll is measured at
start; that measurement, like everything about waiting, cannot reach the output). ThreadSanitizer finds no race (section 7.2).

*What limits it.* The floor is the main thread's own work: the draws, the fitness of each
individual (one `exp` in `gsl_ran_gaussian_pdf`) and handing the line of output to its
thread. At 100 loci the genome thread keeps up and a run takes what the main thread
takes. With more loci the sums take longer than the draws and more threads shorten what
is left.

### 3.12 Threads: performance cores

`QTSIM_THREADS=auto` now gives 2 (the pipeline) for small populations and up to 4 for
large ones, never more than the performance cores where the system tells them apart
(macOS: `hw.perflevel0.physicalcpu`; elsewhere half the processors, as before), and the
pool that sums genotypes is not made larger than the performance cores either
(`QTSIM_POOL` overrides): in the benchmark of 19 November, 8 threads on a machine with 6
performance cores were slower than 4.

### 3.13 Where the fast program spends its time, and what did not help

Shares of samples (`sample`, 1 ms) of the fast program on one thread before the changes
of 3.9 to 3.11, on inputs of the published design (1000 individuals; one run each, a
guide and not a benchmark):

| Share of samples | 100 loci | 1000 loci | 5000 loci |
|---|---|---|---|
| The main loop: draws (inlined), gametes, mutations, fitness, collection; `gsl_ran_gaussian`, `log`, `exp` | 65% | 40% | 21% |
| Genotype sums (`SumGenotypes`) | 12% | 49% | 75% |
| Statistics of the output line (`gsl_stats_mean` 6%, `gsl_fit_linear` 5%, variances and SD 7%, at 100 loci) | 20% | 10% | 3% |
| Writing the output | 2% | 1% | <1% |

So with few loci the serial draws dominate (hence the pipeline and the means), and with
many the sums, which bit parity fixes at one addition per allele in order.

Measured on 19 November 2026 and not kept, because each gave no gain or a loss (the
output was identical in every case):

- wider and narrower lanes in the genotype sums (more or fewer than eight individuals
  side by side), and other sum kernels;
- `-O3`, link-time optimisation and profile-guided optimisation in place of `-O2`;
- forming the first products of the Poisson loop ahead of the draws, to avoid its
  unpredictable exit (2.36 s against 2.18 s at 100 loci: the multiplications wait on each
  other longer than the branch costs);
- with the pipeline: summing genotypes behind the draws at 5000 loci (6.25 s against
  5.62 s on four threads, 7.89 s against 7.80 s on two, indicative: most of those sums are
  thrown away), hence the switch on the expected number of mutations;
- with the pipeline: larger batches for the sums behind the draws (128 individuals
  instead of 64 was slower with four threads at 500 and 1000 loci in two runs each, on a
  busy machine, and no different with two threads; smaller batches were not clearly
  better).

## 4. What bit parity rules out, and what it costs

- **Cached partial sums.** QTLSim.jl's kernel keeps the sum of each block and adds block
  sums. That changes the order of the additions, hence the last bits of a genotype, and in
  time which individuals survive. Here every allele of every individual is added again in
  every generation, in order: `2 × NoQTL × PopSize` additions that cannot be skipped, and
  for each individual they form one chain in which every addition waits for the one
  before. This is why, **on one thread at 5000 loci, QTSimFast (13.6 s) remains slower
  than the Julia kernel (9.0 s)**, which is not bit-exact: the kernel's work per
  generation grows with the number of blocks touched, QTSimFast's with the number of
  loci.
- **Splitting the random stream.** All draws come from one MT19937 stream in a fixed
  order, and a Gaussian or a Poisson consumes a data-dependent number of uniforms, so the
  position in the stream of the n-th draw is not known in advance: the draws cannot be
  divided among threads. With few loci this serial part is most of the run. The pipeline
  (3.11) lets everything else go on beside it; the draws remain the floor of a run's
  wall time.
- **Re-implementing GSL's variance and regression.** GSL's running variance divides at
  every element, each step waiting on the last, and its code has a product next to a sum,
  which some builds of the library fuse; so it is called, not copied. The running mean
  has no product, and only it is written out (3.9), with a check against the library at
  start.
- **Vector `exp`, fused multiply-add, reassociation.** None is used: the build has
  `-ffp-contract=off` and no option that relaxes IEEE arithmetic, and the fitness is the
  original's expression calling `gsl_ran_gaussian_pdf`.

The formal development makes the first point precise: for an addition with no algebraic
law, summing cached block sums, or the loci in reverse, gives a different result on an
explicit example (`QTSimFastShadows.lean`, "negative shadows"), so the theorem proved for
the in-order sums could not have been proved for them.

## 5. Threads and determinism

With one thread (the default) the program is serial. With `QTSIM_THREADS` of 2 or more
there are three kinds of thread: the main thread, which makes every random draw; the
genome thread of the pipeline (3.11); and a pool that sums genotypes; and one more thread
writes the output line. This section describes the pool and the output thread, which
are as they were; the pipeline's own argument is in 3.11.

`QTSIM_THREADS` sets the size of a small pool of POSIX threads (`FastThreads.c`), at most
the performance cores (3.12). Two pieces of work are divided among them:

- the genotype sums (`SumGenotypesJob`): each thread sums a contiguous range of
  individuals, a multiple of eight, and writes only those individuals' entries of `Geno[]`
  and its own array of marks;
- the copying of block numbers (`FillGametesJob`; only without the pipeline, whose
  genome thread copies them): each thread fills a contiguous range of block positions of
  the next population's table, reading the current population's table and the gamete
  records, which no thread writes.

Nothing is reduced across threads: there is no sum, minimum or count whose value could
depend on how the work was divided. The marks are combined afterwards by one thread
with a bitwise OR, which does not depend on order. All random draws and the survivors
are made by the main thread alone; the mutations and the allocation of blocks by one
thread at a time (the main thread, or with the pipeline the genome thread). So each value
is computed by the same operations whichever thread computes it, and the output is the
same for every thread count.

With more than one thread, a further thread computes the statistics and writes the output
line while the next generation starts (`Output`, `WriteOutput`). It works on copies of the
phenotypes, genotypes and survivors made by the main thread before it goes on, draws no
random number, and the main thread waits for it before the next line, before `Pop all
died.` and before exit; so the lines appear in order and are complete.
`QTSIM_ASYNC_OUTPUT=0` turns this off.

This argument is by inspection. It is supported by the tests: every case of the matrix
was run with 2, 3, 4 and 8 threads forced on, twice more with 4 threads, and under
ThreadSanitizer, with no difference and no report (section 7).

## 6. The start-up self-checks

Three assumptions about the GSL library in use are tested each time the program starts.

1. **The inlined generator** (`FastRngSelfCheck`): the inlined MT19937, seeded with the
   run's seed, is compared with `gsl_rng_mt19937` seeded with the same seed over
   3 × 624 + 7 outputs (three regenerations of the state), together with GSL's range for
   the generator. If they differ the program writes a message to standard error and exits
   with code 3; it does not run. Building with `-DQTSIM_GSL_RNG` uses GSL's generator
   directly.
2. **The statistics with a given mean** (`StatsMSame`): on test vectors of 1 to 1000
   numbers, `gsl_stats_variance_m(x, mean(x))` and `gsl_stats_sd_m(x, mean(x))` are
   compared bit for bit with `gsl_stats_variance(x)` and `gsl_stats_sd(x)`. If any
   differs, the program calls the plain functions throughout, as the original does (it is
   then slightly slower and still exact).
3. **The means in one loop** (`StatsMSame`, `ThreeMeans`; 3.9): on the same vectors, and
   on a second vector of mixed magnitudes with a part of it as the third, the three means
   computed side by side are compared bit for bit with `gsl_stats_mean`. If any differs,
   `gsl_stats_mean` and the plain functions are called (`QTSIM_PLAIN_STATS=1` forces
   this path, and the parity matrix runs it).

None of the checks has failed on the machine used.

## 7. Verification

### 7.1 The parity matrix

`./test_parity.sh full` compares standard output (`cmp`) and exit codes with the
original's. Result on the revised code as packaged (run started 19 November 2026 and
finished the next morning, the machine having slept in between), macOS arm64 (Apple M3
Pro, Apple clang 21, GSL 2.8), recorded in `PARITY_FULL.txt`:

- **457 option files** run by the original (66 MB of output; 77 runs end with `Pop all
  died.`), **13,915 comparisons** of QTSimFast with it: **13,915 identical, 0
  different**. (The first version of the code, without the pipeline, passed the matrix as
  it then was: 6,406 of 6,406.)
- The paper's inputs: 72 option files from the eight strata (48 populations die, 24
  survive), 58 of them to 12,000 generations, four of those with 5000 loci. (The paper's
  input files are not in the package; without them the script uses synthetic inputs of
  the same design.)
- The single-locus protocol (6500 generations, jump at 5000) and QTSim's demonstration
  input.
- A synthetic matrix: haploid and diploid; no move and each `OptimumMoveType`;
  recombination means 0, below 1, 3, 12 and far above the number of loci; 1, 2, 3, 7 and
  100 loci and numbers of loci around the block boundaries; population sizes 1 to 1001;
  populations that die at once or under the move; `Sample`, `ChangeFitnessSD`,
  `StandardisePhenotype`, `Drift`; heritability 0, 1 and above 1; no mutation and mutation
  means above 10; seeds 0 and above 2^32; option-file syntax and every error exit; no
  loci, no individuals, no generations.
- By variant of the fast program: one thread 498, and the published inputs with 2, 3, 4
  and 8 threads as they come (36 more); threads forced on small inputs (2, 3, 4, 8) 1518;
  the same with the pipeline off 1558; the pipeline with sums only after the mutations
  1117, with sums behind the draws eight individuals at a time and the main thread always
  given a share 1490, with threads that sleep at once instead of polling 388; repeats with
  4 threads 730; output on the main thread 365; `auto` 365; block sizes 1, 5, 16, 64:
  1100; one individual at a time 369; collection every generation and every 50
  generations 730; GSL's plain statistics calls 365; the build without the inlined
  generator, threads and vector extensions 366; the sanitizer builds 2920 (7.2).

The first version of the code was also run from a clean extraction of its package
(6,277 of 6,277). **For the revised code the full matrix was not repeated from an
extracted package**; there, `make check`, `make parity` (the quick matrix), `make
model-check` and `./demo_parity.sh` were run (`PARITY_FULL.txt` says so).

`PARITY_DEMO.txt` is the transcript of `./demo_parity.sh`: six example inputs, two seeds,
1, 2, 4 and 8 threads, with the SHA-256 and size of each program's output and the wall
time of each run (48 comparisons, all identical), and a negative control (another seed
gives a different output, so "identical" is informative).

### 7.2 Sanitizers

- QTSimFast built with AddressSanitizer and UndefinedBehaviorSanitizer: 365 option files
  with 3 threads, each with the pipeline on, with it off, and with sums behind the draws
  in the smallest batches (1095 runs), no report, output identical.
- QTSimFast built with ThreadSanitizer: 365 option files, each in five modes (the
  pipeline on with 4 threads, off with 4, sums after the mutations with 3, sums behind
  the draws in the smallest batches with 4, threads that sleep at once with 3: 1825
  runs), no report, output identical. Before that, during development: four inputs
  (100, 1000 and 5000 loci of the published design and a haploid one in which every
  block mutates) with 2, 3, 4 and 8 threads in seven modes each, 112 runs under each
  sanitizer build, no report.
- ThreadSanitizer observes the runs it is given: it found no race in them, which is not
  a proof that none exists.
- The original built with AddressSanitizer and UndefinedBehaviorSanitizer ran the 365
  synthetic option files without a report: no undefined behaviour was observed in it on
  these inputs.

### 7.3 Formal proofs

`formal/QTSimFast.lean` (Lean 4 with Mathlib; rendered as `formal/QTSimFast.html` and
`.pdf`) defines a mathematical model of the two programs' algorithms and proves that they
agree. The Lean files were not changed by the revision of 19 November. `PROOFS.md` lists every theorem with its plain statement, the C code it models in
both versions, and its hypotheses.

**What is proved, in plain words.** For a model in which the addition of doubles is an
arbitrary operation with no algebraic law, and the random generator an arbitrary state
machine:

- (a) summing genotypes block position by block position over shared blocks gives every
  individual exactly the original's locus-by-locus sum (`fastGeno_eq`);
- (b) a gamete built by copying block numbers, and one block when the crossover falls
  inside it, stands for the chromosome the original's loop writes (`fillMix_spec`); a
  mutation with copy on write changes that chromosome as the original's `+=` does and
  changes no other (`mutStep_sim`); collecting unused blocks changes nothing and never
  frees a block in use (`collect_inv`);
- (c) the original's loop switches strand once, after the first entry of its sorted map
  (`origLoop_eq`); the fast program's draws consume the generator identically and return
  that entry (`fastLow_eq`);
- (d) making all of `NewPop`'s draws first and copying afterwards equals the original's
  interleaving, given that no draw reads a genome (`interleaved_eq_two_pass`);
- (f) from the common starting state, n generations of the fast model print the same
  records as n generations of the original model, the population dies in the same
  generation if it dies, and the generator ends in the same state, for every block size
  and every schedule of block collections (`gen_refines`, `run_refines`, `init_rel`).

**What is not proved.** This matters more than the list above.

- **The proofs are about a model, not about the C source.** No verified semantics of C
  was used. That the definitions describe what the two C programs do rests on reading the
  code (each definition names the function it models), on an executable cross-check on a
  tiny instance (`make model-check`: the two programs' own routines and the Lean model
  give the same bit patterns for genotypes and gametes, section 7.4), and on the parity
  matrix.
- The compiler, the C library and GSL are outside the model.
- The inlined MT19937 against GSL's is tested at start-up and by the parity matrix, not
  proved.
- `printf` formatting (the same calls on the same values) is by inspection.
- Thread scheduling is not modelled: sections 3.11 and 5 are arguments, supported by
  ThreadSanitizer. For the pipeline the model-level justification is the two-pass
  theorem (d) with its hypothesis `DrawsIgnoreData`; its use for the mutations is by the
  same generic theorem but is not instantiated in the composed theorem (f), which
  describes the program with one thread (`PROOFS.md`).
- Target (e), statistics computed once, is by inspection: in the model both programs
  apply the same uninterpreted functions to the genotypes.
- The model covers the diploid path with at least one individual and one locus. The
  haploid path and the cases `NoQTL = 0` and `PopSize = 0` are covered by tests only.
- The traversal by which `FillGametes` copies block numbers (lists of gametes by cut
  position, a working array per thread) is modelled by its result, not step by step.
- One hypothesis concerns a draw's value: a mutation's locus is below the number of loci,
  which `gsl_rng_uniform_int(rng, NoQTL)` guarantees by its own loop. Another concerns
  `qsort`: that it returns its input sorted.

**Do the theorems say what is claimed?** A second Lean file,
`formal/QTSimFastShadows.lean`, checks each main theorem against its informal claim: it
derives consequences the claim obviously has (special cases such as one block, block size
1, one locus, no crossover, one crossover point, zero generations; structural facts such
as "the other genomes are unchanged by a copy-on-write mutation" and "the generator ends
in the same state"), proves a second, independently phrased formalization equivalent to
the first for (a), (c) and (d), instantiates every hypothesis on a concrete example with
an addition that is neither associative nor commutative, and refutes the corresponding
statement for non-exact algorithms (loci in reverse order, cached block sums, switching
at the highest point or at every point, writing a shared block in place, draws in another
order). The idea comes from the SA-Pass metric of Han et al. (ShadowBench,
arXiv:2608.29270); **it is used here in an adapted form, not as the published benchmark
procedure**, which presupposes a reference formalization that does not exist for QTSim.
`PROOFS.md` has the table of which checks were achieved for which theorem.

### 7.4 The executable cross-check

`formal/model_check.c` includes each program's `QTSim.c` whole and calls its own routines
(the original's `CaclGenotype` and `RecombinationPoissonChr`; the fast copy's
`SumGenotypes`, `RecombinationPoissonChr`, `FillGametes` and `MixBlocks`, with block sizes
1, 2, 3 and 6) on 4 individuals with 6 loci, for alleles whose sums are exact and for
alleles whose sums are rounded, and for five scripted sets of random decisions (a
rejected repeated point, no point, a point at a block's end, three points, as many points
as loci). `formal/ModelCheck.lean` evaluates the Lean model of both programs on the same
instance with hardware doubles. All outputs (bit patterns of 8 genotypes and 60 alleles)
are identical: `formal/model_check_expected.txt`, checked by `make model-check`. This is
one instance, not a proof; it would catch a model that, say, switched strand before
instead of after the crossover locus, or summed the two chromosomes in the other order.

## 8. Limits, portability, use

**Limits.**

- The original stays the reference. QTSimFast has been tested against it on one machine.
- `Sample 0` divides by zero in both programs (`Itter % Opt->Sample`), which C leaves
  undefined; it is not tested.
- `PopSize` or `NoQTL` below 0 are refused by QTSimFast with a message on standard error
  and exit code 1; the original fails on them in its own ways. For all other invalid
  inputs tried (unknown options, wrong arguments, `Ploidy 3`, a missing file, and so on)
  the two programs print the same and exit with the same code.
- Without a `Seed` both programs seed from the time and the process number, so two runs
  cannot be compared.
- Run with no argument, both programs print the name of their main source file as the
  compiler was given it (`__FILE__`). The Makefile compiles each program from inside its
  own directory, as `src/QTSim.c`, so that this line too is the same; built otherwise the
  two messages differ in that path and in nothing else.
- Standard error is not part of the comparison (neither program writes to it in a normal
  run; QTSimFast uses it for notes about ignored environment variables).
- Memory: the block tables are allocated for the worst case (twice the original's
  genomes) but touched only as blocks are used.

**Portability.** Written for gcc and clang on arm64 and x86-64 (vector extensions and
`yield`/`pause` hints behind `#if`, with scalar fallbacks). It has been built and run on
macOS arm64 only (Apple clang 21, GSL 2.8). Linux and x86-64 are written for but
**untested**: run `make check` and `make parity` (a minute) after building elsewhere; the
two start-up checks guard the assumptions about GSL. Both programs must be built without
`-ffast-math` and with `-ffp-contract=off`, as the Makefile does.

**Build and run** (a C99 compiler with C11 atomics and POSIX threads, `make`, and GSL
with `gsl-config` on the `PATH`; no Julia, no Python):

```sh
make                    # builds original/build/QTSim and fast/build/QTSimFast
make check              # the original's files have their published checksums; unmodified files are unmodified
./demo_parity.sh        # six cases, checksums and times, a negative control (a few minutes)
make parity             # quick parity matrix (about a minute)
make parity-full        # the whole matrix (about 1.5 CPU-hours; builds the sanitizer variants)
make model-check        # the C routines against the recorded output of the Lean model
make bench              # timings; writes bench_results.md, .csv, bench_runs.csv
fast/build/QTSimFast options.txt > out.tsv
QTSIM_THREADS=4 fast/build/QTSimFast options.txt > out.tsv
```

Environment variables of QTSimFast (none changes the output; invalid values are ignored
with a note on standard error): `QTSIM_THREADS` (a number or `auto`, default 1; 2 or more
turn on the pipeline),
`QTSIM_BLOCK` (loci per block; default 16 below 400 loci, 32 below 2500, else 64),
`QTSIM_GC` (generations between collections, default 4 below 2500 loci, else 1),
`QTSIM_LANES` (1: one individual at a time), `QTSIM_PAR_MIN` (smallest `PopSize × Ploidy ×
NoQTL` summed on several threads, default 2,000,000; a quarter of it for sums behind the
draws), `QTSIM_ASYNC_OUTPUT` (0: output on the main thread), `QTSIM_PLAIN_STATS` (1: GSL's
plain statistics calls), `QTSIM_STATS` (1: block statistics on standard error at the
end); for the pipeline `QTSIM_PIPELINE` (0: off), `QTSIM_SPEC` (1 or 0: genotypes summed
behind the draws or after the mutations; default by the expected number of mutations per
individual, at most 0.7), `QTSIM_POOL` (threads that sum genotypes), `QTSIM_BATCH`
(individuals claimed at a time, default 64), `QTSIM_HELP_MIN` (fewest alleles left for the
main thread to take a share, default 40,000), `QTSIM_SPIN_MAIN`, `QTSIM_SPIN_GENOME`,
`QTSIM_SPIN_POOL` (microseconds a waiting thread polls before it sleeps: 300, 100, 200).

**Licence.** GNU General Public License, version 3 or (at your option) any later version
(`COPYING`). QTSim is copyright 2023 Andrew Meade. QTSimFast is a modified version (GPL
section 5): every modified or added file carries a notice saying so, with the date, at
its top, above the original notice where there is one.

| File in `fast/src` | Status |
|---|---|
| `GenLib.c`, `GenLib.h`, `PassOptionFile.c`, `PassOptionFile.h`, `RandLib.c`, `RandLib.h` | byte-identical to the original's (`make check`) |
| `QTSim.h` | five function prototypes removed, notice added |
| `QTSim.c` | options code (defaults, option file, printed header, analytical fitness width, `main`) unchanged; the simulation core rewritten; functions the original never calls removed |
| `FastRNG.h`, `FastRNG.c` | new |
| `FastThreads.h`, `FastThreads.c` | new |

`original/src` is byte-identical to the published source (`original/SHA256SUMS`).
The paper's input files, results, manual and scripts carry no licence statement and are
not included.
