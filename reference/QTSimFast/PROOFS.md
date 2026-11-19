# QTSimFast: the formal proofs

This note lists what is proved in `formal/QTSimFast.lean` about the algorithmic changes of
QTSimFast, what each definition models in the two C programs, on which hypotheses, and
how the statements were checked against what they are meant to say
(`formal/QTSimFastShadows.lean`). Rendered versions of the two Lean files, with their
explanations, are `formal/QTSimFast.html`/`.pdf` and `formal/QTSimFastShadows.html`/`.pdf`.

**Read this first.** The theorems are about a mathematical model of the two programs'
algorithms. They are **not** theorems about the C source, the compiler, the C library or
GSL: no verified semantics of C was used. The alignment between the model and the C code
is supported by (i) inspection: every definition names the C function it models, and the
correspondence is tabulated below; (ii) an executable cross-check on a small instance, in
which the two programs' own routines and the Lean model produce the same bit patterns
(section 5); (iii) the parity matrix of the two programs (6,406 comparisons, all
identical). It is not proved. Even with the alignment checks of section 4, this remains
so.

Checked with Lean 4.32.1 and Mathlib v4.32.1. No `sorry`, no added axiom, no
`native_decide`. The main theorems depend only on Lean's three standard axioms (`propext`,
`Classical.choice`, `Quot.sound`), as `formal/check.sh` prints.

## 1. How the model treats numbers and randomness

- **Floating point.** Alleles and sums live in an arbitrary type `F` with an arbitrary
  operation `add : F → F → F` and constant `zero`. Nothing is assumed of `add`: not
  associativity, not commutativity, not `add zero x = x`. Any theorem `lhs = rhs` in `F`
  therefore holds in particular for IEEE doubles with the processor's rounded addition,
  where it says that the two sides are the same bit pattern. (Real-number algebra is
  never used.)
- **Random numbers.** A draw is an arbitrary function `S → α × S` on an arbitrary type
  `S` of generator states (`Rand S α`). Nothing is assumed about the values drawn or
  about how many raw numbers a draw consumes, with one exception: `hpos` below.
- **What the two programs share.** Everything computed in the same way by both programs
  from the genotypes (phenotype, the move of the optimum, fitness, the survival test, the
  output line, the change of the fitness width, each elementary draw) is a field of the
  structure `Prims`, with no property assumed. The model thus assumes that both programs
  apply the *same* function there; that is target (e), checked by reading the code.

## 2. Informal claims (the reference text)

These were written before, and independently of, the Lean statements; the shadow checks
of section 4 are derived from them.

**(a) Genotype sums.** For each individual the original computes
G = (…((0 ⊕ (a₀ ⊕ b₀)) ⊕ (a₁ ⊕ b₁)) …) ⊕ (a_{L−1} ⊕ b_{L−1}), where a_l and b_l are the
alleles of its two chromosomes at locus l, ⊕ is the rounded addition of doubles and L is
`NoQTL` (`CaclGenotype`). The fast program keeps one partial sum per individual and
advances all of them over the loci of block position 0, then of position 1, and so on
(`SumGenotypesW`); inside a block eight individuals advance together, two loci per vector
step; ranges of individuals may be given to different threads. "Same" means: for every
individual, the fast program's final partial sum is the double G above, bit for bit,
whatever the block size, grouping or number of threads. No random number is involved.

**(b) Shared blocks.** The original holds every chromosome as `NoQTL` doubles and copies
a gamete locus by locus. The fast program holds a chromosome as block numbers into a
store of blocks. "Same" means: at the points where the original's arrays are read (when
genotypes are summed), every chromosome read through its block numbers has, at every
locus, the double the original's array has. This must survive (i) building gametes by
copying block numbers and building one mixed block where a crossover falls inside a
block; (ii) mutations, which add a Gaussian to one allele of one chromosome and must not
change any other chromosome even though blocks are shared; (iii) the reuse of blocks
that no chromosome holds any more.

**(c) Crossover.** For one gamete the original draws K ~ Poisson; if K = 0 or K ≥ `NoQTL`
there is no crossover; otherwise it draws K distinct points in 0 … `NoQTL`−2, redrawing
any point already drawn, sorts them, then draws the starting strand; its copying loop
switches strand once, after the locus equal to the smallest point (its index into the
map is never advanced). The fast program makes the same draws with the same redraw
decisions and keeps the running minimum. "Same" means: the generator is left in the same
state, and the gamete "starting strand for loci ≤ minimum, other strand after" is, locus
for locus, the gamete the original's loop writes.

**(d) Two passes.** The original makes the draws for offspring 1, copies its genome,
makes the draws for offspring 2, and so on. The fast program makes all the draws, then
copies all the genomes. "Same" means: the same decisions for every offspring and the same
generator state at the end, hence the same genomes. This requires that no draw depends on
the contents of a genome.

**(f) Composition.** Started from the same options and seed, after every generation the
two programs have: the generator in the same state; equal environmental deviations and
fitness parameters; chromosomes that are equal in the sense of (b); the same output line.
A population that dies does so in the same generation in both.

## 3. The theorems

For each: the Lean name, the statement in plain words, what the model's objects are in
the C code of both versions, and the hypotheses with their justification.

### (a) `foldl_pointwise`, `fastGeno_eq` — fully proved

| | |
|---|---|
| Statement | `foldl_pointwise`: updating every accumulator at each step of an outer loop leaves in accumulator `i` the left fold of `i`'s own steps. `fastGeno_eq`: for every layout of the loci in blocks, every store and table of block numbers and every individual `i`, `fastGenoAll add zero y data ids i = origGeno add zero y.L (fun g => absChr y.B data (ids g)) i`. |
| Original | `origGeno` is `CaclGenotype` (lines 323–351, diploid branch): `acc` is `Genotype`, the fold over `List.range L` is the `for` over `Index`, the term `add (pop (2i) l) (pop (2i+1) l)` is `Genome[0][Index] + Genome[1][Index]`, added to `Genotype` on the right of `+=`. |
| Fast | `fastGenoAll` is `SumGenotypesW`: the outer fold over `List.range y.NB` is the loop over `K`; the function `acc` is the array `Geno[]`; `blockSum` is the inner loop over `Index < Len` for one individual; `data (ids g k) o` is `Block(S, Ids[K*Row + g])[o]`. The lanes `J`, the two-loci vector step (`foldl_pairs`) and the division among threads are instances of `foldl_pointwise` (independent accumulators); the padding of a short last group with a repeated individual is not modelled (its result is not stored). |
| Hypotheses | None beyond the layout: `0 < B`, `0 < LastLen ≤ B` (as `CreateSim` sets them). |
| Assumed outside the model | An elementwise vector addition gives each lane the scalar sum (IEEE; no fused or reassociated operations: `-ffp-contract=off`). |

### (b) `fillMix_spec`, `mutStep_sim`, `collect_inv` — fully proved for the model described; one traversal not modelled

| | |
|---|---|
| Statement | `fillMix_spec`: after `FillGametes` and `MixBlocks`, for every gamete `g` and locus `l`, the gamete read through its block numbers is the starting chromosome's allele if `l ≤ low g` and the other chromosome's otherwise; the store invariant `MutInv` holds. `mutStep_sim`: one mutation keeps `MutInv`, and every chromosome read through its block numbers equals the dense population after the original's `Genome[ChrNo][Pos] += effect` (so the mutated chromosome changes at that locus only, and no other chromosome changes). `collect_inv`: after `CollectBlocks`, the allocator is sound and can hand out only blocks that no chromosome of the population holds; block contents are untouched. |
| Original | Dense chromosomes `pop g : ℕ → F` are `Genome[c]` of individual `g / 2`. `mutDense` is the `+=` in `MutateIndividualChrNo` (line 409). |
| Fast | `Store`: `data` = `S->Data`, `free` = `S->Free` (top first), `hw` = `S->HW`, `owner` = `S->Owner`. `Store.newBlock` = `NewBlock`. `absChr B data ids l = data (ids (l / B)) (l % B)`. `GDesc`/`gdesc` = what `RecombinationPoissonChr` (fast) records: `Src1`, `Src2`, `Cut`, `Mix`. `fillRow` = the result of `FillGametes`; `mixStep` = one entry of `MixBlocks` (including its `if(CId == OId) continue`). `mutStep` = the body of the `while` in `MutatePop` (`Owner[Id] != Tag` decides copy or in place). `scan`/`collect` = `CollectBlocks`; `markOf` = the marks set in `SumGenotypesW`. |
| Invariant | `MutInv`: the free list has no repeats and lies below the high-water mark; no block the allocator can hand out is held by a chromosome; a tagged block is held by exactly one cell, of the chromosome its tag names. `MixInv` adds: no block held by the parents has been written. Proved to hold initially (`init_rel`) and to be preserved by every step. |
| Hypotheses | Gametes' parents are individuals of the population; `low ≤ L`; the mutated locus is below `L` (see `hpos`). |
| Not modelled | The fast program has no reference counts (the task description suggested them); the model follows the code: tags and mark-and-collect. `FillGametes`' position-major traversal with the lists `Head`/`Next` and the array `Cur` (and its threads) is modelled by its result `fillRow`. `MixBlocks` iterates over the list `Mix`; the model iterates over all gametes and acts on those with a `Mix` entry (the same entries in the same order). The marks of several threads combined by OR are modelled by one mark function. The capacity check of `NewBlock` (an internal error exit) is not modelled. |

### (c) `origLoop_eq`, `runMin_eq_head`, `fastLow_eq` — fully proved

| | |
|---|---|
| Statement | `origLoop_eq`: the original's loop writes, at each locus `l < L`, the starting strand's allele if `l ≤ CombMap[0]` and the other strand's otherwise, whatever else the map contains. `runMin_eq_head`: the running minimum of the drawn points is the first entry of any sorted rearrangement of them. `fastLow_eq`: `LowestCrossover` leaves the generator in the state `CreateCombMap` leaves it in and returns `min L CombMap[0]`. |
| Original | `origMap` = `CreateCombMap` (lines 608–636): `count` = `gsl_ran_poisson(rng, RecombPoissonMean)`; the test `K = 0 ∨ L ≤ K` = `*NoPoints == 0 \|\| *NoPoints >= NoQTL`; `drawPoints` = the `for` with its `do … while(PointInMap(…))`; `sort` = `qsort(…, IntComp)`; the appended 0 = `Map[*NoPoints] = 0`. `recombStep`/`origLoop` = the loop of `RecombinationPoissonChr` (lines 682–693) with `PointIndex` the constant 0; `start` = the outcome of `gsl_rng_uniform(rng) < 0.5` (true: `CP = P2`); `chr0` = what the gamete's array held before. |
| Fast | `fastLow` = `LowestCrossover`; `drawPointsMin` = its loop with `S->Points` and `Low`; `stepMin` = `if(Point < Low) Low = Point`. For more than 16 points the code tests membership in a table of stamps: `seenTable_eq_mem` shows such a table marks exactly the points drawn. |
| Hypotheses | `SortSpec sort`: `qsort` returns a sorted rearrangement of its input (nothing else about it). `fresh`, the redrawn `gsl_rng_uniform_int(rng, NoQTL-1)`, is an arbitrary function of the membership test and the generator state: both programs call it with the same test. |
| Note | No assumption is made on the range of the points: the result is stated with `min L`, which is all the gamete needs (`shadow_gamete`). |

### (d) `interleaved_eq_two_pass` — fully proved, under a named hypothesis

| | |
|---|---|
| Statement | If the draws of every item do not depend on the data (`DrawsIgnoreData`), the interleaved loop (draw, apply, draw, apply, …) equals: all the draws in the same order (`drawAll`), then all the applications (`applyAll`); same data, same final generator state. |
| Original | `interleaved` with `origDecide` (the draws of `NewIn`, `RecombinationPoisson`, `SetIndividualEnv` for one offspring) and `origApply` (the two copying loops and the store of `Env`) is the loop of `NewPop` (lines 780–786). |
| Fast | `drawAll` with `fastDecide` is the first loop of the fast `NewPop`; the second pass is `fillMix`. |
| **Hypothesis `DrawsIgnoreData`** | No argument of a GSL call and no branch around one, between `SetSurviving` and the end of `NewPop`'s loop, reads `Genome[][]`. Checked in `original/src/QTSim.c`: `SetIndividualEnv` 225–236 (argument `EnvVar`); `ValidDiploidParents` 479–492 (compares the two pointers; `AssortativeMating` is −1, the options parser offers no way to set it, so no phenotype is read); `PointInMap`/`CreateCombMap` 591–636 (arguments `RecombPoissonMean`, `NoQTL-1`; the redraw test reads the map); `RecombinationPoissonChr` 667–697 (`gsl_rng_uniform(rng) < 0.5`; the genomes are only copied); `NewIn` 707–740 (argument `NoSPop`); `NewPop` 755–788. In the model the hypothesis holds by construction (`origDecide` has no genome argument); the theorem is stated for draws that *could* read the data so that the hypothesis is explicit, and `QTSimFastShadows.lean` shows it cannot be dropped. |

### (e) Statistics computed once — not a theorem

A function without side effects evaluated once equals the same function evaluated twice;
that is true by definition in Lean and a lemma would say nothing. By inspection: the fast
`PopFitness`, `SetFitnessMean`, `Output`/`WriteOutput` call the same GSL functions on the
same arrays as the original's, each once; the values the original computes and never uses
(an SD in `SetFitnessMean`, the variance of `Env` in `Output`, a correlation in
`OuputReg`) draw no random number. `gsl_stats_variance_m`/`gsl_stats_sd_m` with the mean
already computed are tested against the plain functions at start-up (`StatsMSame`).

### (f) `newPop_refines`, `mutLoop_sim`, `gen_refines`, `run_refines`, `init_rel` — fully proved for the diploid model

| | |
|---|---|
| Statement | `gen_refines`: from related states (`Rel`: same generator state, deviations and parameters; every chromosome's block numbers stand for the original's chromosome; the store invariant), one generation of the fast model and one of the original model either both end with the population dead, or give the same output record and related states. `run_refines`: hence for `n` generations: the same list of output records, death in the same generation, related final states, for every schedule `gc` of block collections. `init_rel`: the two starting states are related. |
| Original | `origGen` = the body of `RunSim`'s loop (lines 1005–1020): `origGeno` + `Prims.popFitness` (= `PopFitness`), `survLoop` (= `SetSurviving`), `Prims.dead` (= `DeadPopulation`; `none` = `Pop all died.`), `origNewPop` (= the loop of `NewPop`), `Prims.output` (= `Output`), `Prims.change` (= the `ChangeFitnessGen` test), `mutLoop … mutDense` (= `MutatePop`), swap. `OState.other` is `NPop`, overwritten below row `2N`. |
| Fast | `fastGen` = the body of the fast `RunSim`'s loop: `fastGenoAll`, `collect` when `gc` says so (`CaclPopPhenotype`'s `Marking`), the same `Prims`, `drawAll fastDecide`, `fillMix`, `mutLoop … mutStep`, the tags cleared (`Store.clearOwner`), swap. |
| Hypotheses | `SortSpec sort` (as in c). `hpos : ∀ s, (pr.mutPos s).1 < y.L`: a mutation's locus is below the number of loci; `gsl_rng_uniform_int(rng, NoQTL)` (and `FastUniformInt`) return only values below their bound, by their own loop, whatever the generator. `hN : 0 < N`: the population is not empty (`PopSize = 0` is a separate branch of the fast `RunSim`, covered by tests). |
| Scope | Diploid, `L ≥ 1`. Not modelled: the haploid path, `NoQTL = 0`, `PopSize = 0`, the creation of the populations (both programs draw `2 × PopSize` environmental deviations before the loop: by inspection), `Sample` (inside `Prims.output`), threads, the asynchronous output thread. |

### Not provable here

That the C code implements the model; the inlined MT19937 against GSL's (tested at
start-up and by the parity matrix); `printf` formatting; thread scheduling (argued:
disjoint writes, no reduction across threads; ThreadSanitizer clean on 365 option files).

## 4. Do the statements say what is claimed? (alignment checks)

A compiled theorem may still miss its intended meaning. `formal/QTSimFastShadows.lean`
checks each main theorem against the informal claims of section 2, following the idea of
SA-Pass (Han et al., *ShadowBench: Toward Reliable Automatic Evaluation of Semantic
Alignment in Autoformalization*, arXiv:2608.29270): auxiliary "shadow" statements S₁ … Sₙ
around an intended statement T, with machine-checked forward implications T ⇒ Sᵢ
(failures reveal omitted conclusions or added assumptions) and a backward implication
(failures reveal added conclusions or omitted assumptions).

**This is an adaptation, not the published benchmark procedure.** SA-Pass is an
evaluation metric that presupposes a reference formalization and a shadow set built from
it. No reference formalization of QTSim exists. Here the shadows were derived from the
informal claims; the backward check is replaced, where it was done, by a second,
independently styled formalization T′ proved equivalent to T; and non-vacuity and
negative examples were added. Passing these checks raises confidence that the Lean
statements mean what section 2 says; it does not certify it.

| Theorem | Informal claim | Forward shadows proved | Backward / second formalization | Non-vacuity instance | Negative shadows (refuted on an explicit example) | Executable cross-check with C |
|---|---|---|---|---|---|---|
| (a) `fastGeno_eq` | 2(a) | 5: one block (`shadow_one_block`), block size 1 (`shadow_unit_blocks`), one locus (`shadow_one_locus`), a genotype depends only on its own alleles (`shadow_local`), block size irrelevant (`shadow_layout_free`); plus 3 numeric instances (= 1485) | second formalization: recursion `sumRec` instead of a fold; `sums_T_iff_T'` proves T ⇔ T′ without using either | layout 6 loci in 3 blocks with shared blocks, `add a b = 2a + b` (not associative, not commutative) | 2: loci in reverse order (738 ≠ 1485); cached block sums, the Julia kernel's method (453 ≠ 1485). A third candidate, summing each chromosome separately, happens to agree for this `add` and was dropped: not every regrouping is caught by one test operation | yes: genotypes of 4 individuals, exact and rounded alleles |
| (c) `origLoop_eq`, `fastLow_eq` | 2(c) | 5: whole gamete equal (`shadow_gamete`), generator state equal (`shadow_rng_state`), no points (`shadow_no_points`), one point (`shadow_one_point`), `K ≥ L` (`shadow_many_points`); plus numeric instances incl. a rejected redraw | second formalization of the loop: recursion `gamRec` producing a list; `origLoop_eq_gamRec`. None for the draws | insertion sort satisfies `SortSpec` (`sortSpec_insertionSort`); a counter generator with redraws | 4: every crossover applied (`recombAll`), switching at the highest point, the first point drawn instead of the lowest, the strand drawn before the map | yes: gametes for 5 scripted decision sets |
| (b) `mutStep_sim`, `fillMix_spec`, `collect_inv` | 2(b) | 5: other genomes unchanged (`shadow_others_unchanged`), the mutated genome changes at one locus (`shadow_mutated`), no crossover copies one row (`shadow_no_crossover`), collection leaves contents and tags (`shadow_collect_data`), a block handed out after a collection is unused (`shadow_collect_fresh`); plus numeric instances incl. a second mutation in place | none | the starting store satisfies `MutInv` (`mutInv_init`) | 1: the same write without copying changes the other chromosome | yes for gametes (blocks of 1, 2, 3, 6 loci); mutation and collection only through the parity matrix |
| (d) `interleaved_eq_two_pass` | 2(d) | 3: generator state (`shadow_two_pass_state`), one item (`shadow_two_pass_one`), data determined by the decisions (`shadow_two_pass_data`) | second formalization: a fold over (data, state) pairs; `interleaved_eq_fold` | `drawGood` satisfies `DrawsIgnoreData`; both orders give (210, 3) | 2: a draw that reads the data violates the hypothesis and the two orders differ (the hypothesis is necessary); draws made in another item order differ | no (the order of draws is exercised only by the parity matrix) |
| (f) `gen_refines`, `run_refines` | 2(f) | 6: zero generations (`shadow_zero_generations`), n+1 generations = n then one more (`shadow_one_more`), final generator state, deviations and parameters equal (`shadow_final_rng`), same extinction (`shadow_same_extinction`), collection schedule irrelevant (`shadow_gc_free`), block size irrelevant (`shadow_layout_free_run`) | none | a complete instance (4 individuals, 6 loci, counter generator, `add a b = 2a + b`) satisfies every hypothesis: `instance_outputs`; executed for 3 generations under three collection schedules (`#guard`) | 1: a generation with cached block sums prints other records (`wrongGen`, `#guard`) | no: one whole generation was not run against C (the model's `Prims` are abstract); gametes and sums only |

Blanks are real: no backward or second-formalization check was made for (b) and (f), and
none for the draw half of (c). Examples closed by `decide +kernel` are checked by Lean's
kernel; those marked `#guard` are evaluated by Lean's compiler and are tests, not proofs.
Two lessons from doing this: four of the numeric instances first written by hand were
wrong and were corrected by the checker (so the instances are not rubber stamps), and one
intended negative example turned out not to be one.

## 5. Executable cross-check of the model against the C code

`formal/model_check.c` is compiled twice (`make model-check`), each time including one
program's `QTSim.c` whole with `main` renamed, so that it calls the program's own
routines and not copies of them:

- original: `CaclGenotype`, `RecombinationPoissonChr` (which calls `CreateCombMap`);
- fast: `SumGenotypes`, `RecombinationPoissonChr`, `FillGametes`, `MixBlocks`, with block
  sizes 1, 2, 3 and 6 (built with `-DQTSIM_GSL_RNG`, so that its draws go through a
  `gsl_rng`; the inlined generator is not part of this check).

Instance: 4 diploid individuals, 6 loci; allele of chromosome g at locus l either
`10g + l` (sums exact) or `scale[l] + 1/(g+3)` with scale = 1, 2^53, 3, −2^53, 0.5, 7
(sums rounded; `ModelCheck.lean` verifies that for these values the reversed sum and the
sum of block sums differ from the in-order sum, so a wrong order would show). The random
generator is replaced by a script of raw 32-bit outputs chosen so that GSL returns the
decisions: (1) two points 4, 4 (rejected), 1, second chromosome first; (2) no point; (3)
one point, 3 (the end of a block of 2); (4) three points 2, 4, 3; (5) six points, as many
as loci, so no crossover.

`formal/ModelCheck.lean` evaluates the Lean definitions (`origGeno`, `fastGenoAll`,
`origMap` + `origLoop`, `fastLow` + `gdesc` + `fillMix` + `absChr`) on the same instance
with `F := Float` (hardware doubles), the generator being the list of decisions. It
checks that the model of the original and the model of the fast program print the same,
and prints the bit patterns.

Result (19 November 2026, macOS arm64): the three outputs are byte-identical: 8
genotypes and 60 gamete alleles as 64-bit patterns. They are recorded in
`formal/model_check_expected.txt`, which begins

```
variant exact
geno 0 4636033603912859648
geno 1 4644513037586268160
geno 2 4648224988841639936
geno 3 4650336051166969856
gamete 1 4629137466983448576 4629418941960159232 4626885667169763328 4627167142146473984 4627448617123184640 4627730092099895296
```

(gamete 1: loci 0 and 1 from chromosome 3, values 30 and 31; loci 2 to 5 from chromosome
2, values 22 to 25: the switch comes after locus 1, the lowest point, although 4 was also
drawn). `make model-check` rebuilds the two C harnesses and compares their output with
this file; `formal/check.sh` reruns the Lean side.

What this does and does not show: on one instance, the Lean definitions compute what the
C routines compute, including the direction of the strand choice, "after, not before, the
crossover locus", the order of the two chromosomes in each term, the rejection of a
repeated point, and the treatment of `K ≥ NoQTL`. It does not cover mutation, block
collection, a whole generation, or the inlined generator, and one instance is not a
proof.

## 6. How to check the proofs

```sh
formal/check.sh /path/to/a/Lake/project/that/depends/on/Mathlib
```

type-checks `QTSimFast.lean` (which imports only Mathlib), then `QTSimFastShadows.lean`
(which imports it), prints the axioms of the main theorems, and runs `ModelCheck.lean`.
Toolchain used: `leanprover/lean4:v4.32.1`, Mathlib `v4.32.1`. By hand, in such a
project: copy `QTSimFast.lean` to its root and run `lake env lean QTSimFast.lean`.

`formal/render.sh` regenerates the `.md`, `.html` and `.pdf` documents from the Lean
files (mdgen, or the fallback converter `formal/lean2md.awk`; Quarto). Each document
states the SHA-256 of the Lean file it was generated from.
