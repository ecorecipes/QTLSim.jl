# Naturality and model morphisms

Roadmap items C2 (naturality as the validation harness) and C7 (morphisms between
models), with the usability items of F11 that fall out of them (`explain`, the support
matrix). The aim is consolidation: the category theory does work that the package did by
hand, one script or one `if` at a time.

## Agreement between backends is a naturality square

The life cycle is one morphism `f` of a free symmetric monoidal category, and each backend
is a functor out of it: the kernel `K(f)` (a random map of `World`s), SLiM (`lower(f, :eidos)`),
the moments of a Gaussian population `M(f)`, and closed forms of theory. A *statistic* σ
maps the states, or the output tables, of one backend to those of another: a kernel
population's sample mean and variance are a moments state; the kernel's records are the
columns SLiM prints. The central empirical claims of the package are squares

    σ ∘ K(f) = M(f) ∘ σ

that commute exactly (in law: the kernel and SLiM; to rounding: the neutral moments and
`deme_coalescence`) or up to a closure's error (Gaussian selection, the rare-alleles
closure). Every script of `validate/` checks some such square on some cells.

### Cases as data

A square is a [`NaturalityCase`](@ref): the life cycle (`cycle`, a morphism; `nothing` is
QTSim's generation), the parameters (`params`, one per deme for a metapopulation), the two
backends (`from`, `to` among `:kernel`, `:slim`, `:qtsim`, `:moments`, `:theory`), the statistics
(`statistic => tolerance` pairs of [`QTLSim.Statistic`](@ref)s, each a function of a
backend's output table, the same function for the backends that share columns), the
criterion, the number of replicates of the stochastic side, the horizon, options for the
backends, and the verdict (`expected`) the ported script recorded on the same cell. The
criteria are the scripts' rules, written before the runs:

| criterion | sides | fails when |
|:--|:--|:--|
| `:band` | stochastic → deterministic | the whole 95% CI of the replicates' mean lies outside `prediction·(1 ± tol)` |
| `:z` | stochastic → deterministic (exact predictions) | `|mean − prediction|/SE ≥ tol` |
| `:exact` | deterministic → deterministic | the largest relative difference (a vector: every generation) exceeds `tol` |
| `:ks` | stochastic → stochastic | the two-sample Kolmogorov–Smirnov p-value is below `tol` |
| `:welch` | stochastic → stochastic | Welch's \|t\| of the replicates' values is at least `tol` |
| `:diff` | stochastic → deterministic | the point difference \|mean − prediction\| exceeds `tol·scale` (no allowance for sampling error) |

[`check_naturality`](@ref) runs a case and returns a [`NaturalityResult`](@ref) (per
statistic the two sides with standard errors, the measure, pass/fail, and whether the
expected verdict was reproduced). The kernel side records QTSim's columns plus each deme's
mean genotype (`geno_mean`), which the variance among deme means needs, through the image
of the `record` step; the SLiM side runs the generated model (`slim_script`) and parses
its table.

`validate/naturality_cases.jl` is the table, `validate/naturality.jl` runs it (writes
`validate/naturality.log`). A new step, model or closure gets its validation by adding a
row. The scripts and their logs stay: they are the record, and the rows are independent
replications of some of their cells (the harness's seeds are its own).

### The ported rows

`validate/naturality.log` (2026-11-11 for these eleven rows, 3 threads, 826 s; rerun 2026-11-13 with the
nine rows below, 4 threads, 1052 s, 72 CPU-minutes, the same numbers for these rows): every verdict reproduced, on the harness's own seeds (replicate r draws from `Xoshiro(seed + r)`).

| case (script, cell) | backends, criterion | statistics: harness (script's log) | verdict |
|:--|:--|:--|:--|
| `moments_loci.jl` G1: L = 10, N = 1000, μ = 0.01 | kernel → moments (`:multilocus`, erosion), band 10%, R = 12 | V_G −2.8% (log −4.0%) | PASS, as logged |
| `moments_loci.jl` house of cards: L = 10, μ = 1e-4 (an expectation, not a criterion) | as above | V_G −90% (log −91%) | FAIL, as logged |
| `moments_regimes.jl` R1: N = 100, 4Nμ = 1, L = 10, V_s = 50.5 | kernel → moments (`:rare_alleles`, erosion), band 15%, R = 16 | V_G +1.5% (log +1.5%) | PASS, as logged |
| `moments_regimes.jl` R1: N = 1000, 4Nμ = 0.1, L = 1000, V_s = 5.5 | as above | V_G −29.7% (log −27.9%) | FAIL, as logged |
| `moments_demes.jl` C: local optima, k = 5 | kernel → moments (`between`), band 10/15/15%, R = 100 (log 300) | V_W −3.5%, among +6.1%, Q_ST +6.9% (log −3.5, +6.3, +7.1%) | PASS, as logged |
| `moments_demes.jl` A: neutral, k = 5 | kernel → moments (`between`), z < 4, R = 100 | z = −2.15, −2.52, −1.57 (log −1.36, −0.30, 0.91) | PASS, as logged |
| neutral metapopulation (test/demes.jl) | moments → `neutral_deme_variances`, exact 1e-10, every generation | 1.3e-14, 2.4e-14 | PASS, as tested |
| `rare_alleles_traits.jl` T: two traits, L = 20, μ = 1e-3 | kernel → moments (`:rare_alleles`), band 20/10/20%, R = 24 | G₁₁ +8.7%, G₁₂ +5.0%, G₂₂ +5.9% | PASS, as logged |
| `rare_alleles_traits.jl` D1, `:rare_alleles` | kernel → moments (`between`), band 15%, R = 48 | G₁₁ +141%, G₁₂ +12.9%, G₂₂ +89.7% | FAIL on G₁₁ and G₂₂, as logged |
| `local_adaptation.jl` D1, `:rare_local` | as above, R = 36 | G₁₁ +2.3%, G₁₂ +8.7%, G₂₂ +16.3% (within the CI rule) | PASS, as logged |
| `slim/compare_demography.jl`, Beverton–Holt equilibrium | kernel → SLiM, KS p ≥ 0.001, R = 40 | size p = 0.91, V_G p = 0.055 (log 0.91, 0.91) | PASS, as logged |

The percentages are (kernel − prediction)/prediction (the scripts report some of them the
other way round); in the `local_adaptation` row G₂₂ is 16% off and passes only because the
rule fails a cell when the whole CI is outside the band, as the script's rule does.

### QTSim as a backend, and the remaining families (2026-11-13)

`:qtsim` is a stochastic backend: each replicate writes an options file with
`write_qtsim_options` and its own seed, runs the C reference (`ENV["QTSIM"]`, the bundled
build `reference/QTSim/build/QTSim`, or the paper's binary `../Repository/Code/QTSim/src/QTSim`;
`options.qtsim.qtsim` overrides) and parses its tab-separated table into the kernel's
columns `gen`, `pheno_mean`, `geno_var`, `nsurv`, `theta`, `omega` and `popsize` (fixed),
with `extinct` set when QTSim prints `Pop all died.` QTSim prints six decimals and no line
for the generation in which the population died (the kernel records none either). It runs
QTSim's model only: a case with several demes, a life cycle other than QTSim's generation
(steps with their own parameters included), or a Julia-only setting (drift, several
traits, a reaction norm or plasticity, a fluctuating move, `CrossoverModel = :all`, density
regulation, fecundity selection) is an `ArgumentError` before anything runs, since QTSim
would silently run another model.

Two criteria were added for the scripts that use them: `:welch` (two stochastic backends,
Welch's t of the replicates' values, fails when |t| ≥ tol; `compare_to_qtsim.jl`) and
`:diff` (stochastic → deterministic, the point difference alone, fails when
|mean − prediction| > tol·scale, no allowance for sampling error; `reaction_norm.jl`). A
statistic may return `missing` for a replicate; `QTLSim.surviving(s)` does so on extinct
runs, as the scripts summarise surviving runs only, and the result reports how many
replicates were used. New statistics: `QTLSim.lag_mean(window; squared)` (θ − z̄ or its
square) and `QTLSim.persistence(start)` (generations from `start` to the last record, plus
one: the scripts' "generations after the move").

| case (script, cell) | backends, criterion | statistics: harness (script's log) | verdict |
|:--|:--|:--|:--|
| `compare_to_qtsim.jl` A, stabilising equilibrium | kernel → QTSim, Welch \|t\| < 4, R = 40 | V_G t = 0.84 (script −0.22) | PASS, as recorded |
| `compare_to_qtsim.jl` B, directional (phenotype move) | as above, R = 24 | V_G t = 0.96, advance 1.16, survival −0.76 (script −1.66, +0.27, −0.33) | PASS, as recorded |
| `slim/compare_fluctuating.jl` periodic A = 0.1, P = 200 | kernel → SLiM, KS p ≥ 0.001, R = 40 | lag², V_G, survivors p = 0.99, 0.91, 1.0 (log 0.029, 0.573, 0.759) | PASS, as logged |
| `slim/compare_fluctuating.jl` red noise σ = 0.1, ρ = 0.9 | as above | lag², V_G, survivors, SD θ, autocorr. θ p = 0.57, 0.40, 0.91, 0.57, 0.76 (log 0.99, 0.76, 0.91, 0.91, 0.91) | PASS, as logged |
| periodic optimum (test/fluctuating.jl) | moments → `expected_lag2`, exact 1e-3 | 2.1e-13 | PASS, as tested |
| `slim/compare_demography.jl` R = 2, move 0.02 | kernel → SLiM, KS p ≥ 0.001, R = 20 (script 40) | extinction time, size before the move p = 0.82, 0.56 (log 0.263, 0.573) | PASS, as logged |
| `slim/compare_plasticity.jl` 1: b/B = 0.5, move 0.01 | kernel → SLiM, KS p ≥ 0.001, R = 20 (script 40) | V_G before, lag, V_G, survivors p = 0.56, 0.98, 0.56, 0.98 (log 0.759, 0.759, 0.988, 0.759) | PASS, as logged |
| `reaction_norm.jl` A, C1 at +250, +500 | kernel → moments (`:multilocus`, erosion), \|Δ\| ≤ 0.1, R = 40 | slope 0.765 vs 0.772, 0.880 vs 0.883 (log 0.781, 0.892 vs the same) | PASS, as logged |
| `fst.jl` island k = 4, Nm = 1 | kernel → `neutral_deme_variances` with `migrate!`'s pair law, z < 4, R = 400 | V_W, among, Q_ST z = −0.59, −0.45, −0.03 (log −2.57, −2.11, −0.57) | PASS, as logged |

`compare_to_qtsim.jl` has no committed log; its verdict is that of CI's QTSim parity job
and of a rerun of its fixed seeds here. The `fst.jl` row's prediction is the script's
exactly (the same Monte Carlo estimate of `migrate!`'s law of pairs' origins, seed 7002),
but its z omits the prediction's Monte Carlo SE that the script folded in, so it is
slightly stricter. Not ported: `reaction_norm.jl`'s two recorded failures (C2, the lag
against the quasi-steady lag at each run's own G and slope, a function of the same kernel
runs rather than of a second backend); `coalescence.jl` (needs pedigree recording in the
kernel runner, and its statistics are per pair of genomes); `slim/three_way_confirm.jl`
(40,000 runs per program).


### The square of each step

[`naturality_steps`](@ref) checks one generation step by step. It starts from matched
states: a kernel population `w` (burned in, then its alleles shuffled among chromosomes
locus by locus, so it is in linkage equilibrium and Hardy–Weinberg proportions) and the
moments state σ(w) with its mean genotype and genetic variance. It applies the steps of
the life cycle in order, to many copies of `w` with their own random numbers and to σ(w),
and compares after each step the observables of the step's codomain stage (mean genotype,
genetic variance, number of individuals). The first square that fails, after squares that
hold, names the step that breaks naturality.

On 100 loci, N = 1000, V_P = 1, h² = 0.5, MutationScalar 0.02, μ = 0.01, burned in for
2000 generations and shuffled to linkage equilibrium, 200 replicates of one generation
(`naturality_steps(p; loci, replicates = 200, burnin = 2000)`); V_G at each stage:

| step | stage | kernel (± SE) | `:aggregate` | `:multilocus` |
|:--|:--|--:|--:|--:|
| develop | Adults | 7.538 | 7.538 (exact) | 7.538 (exact) |
| select | Survivors | 5.109 ± 0.013 | 5.174 (−1.3%) | 5.174 (−1.3%) |
| reproduce | Offspring | 6.127 ± 0.020 | 5.174 (+18.4%, **fails**) | 6.356 (−3.6%) |
| mutate | Offspring | 6.137 ± 0.020 | 5.184 | 6.366 |
| swap | Zygotes | 6.137 ± 0.020 | 5.184 | 6.366 |

(ω = 4.) The aggregate recursion's square fails first at `reproduce`: it inherits all of
selection's reduction, while recombination and segregation restore about half of it. The
multilocus model's squares hold at a 5% tolerance (its reproduce is 3.6% high, the same
sign as its equilibrium bias in `validate/moments_loci.jl`, kernel/prediction 0.96–0.99).
Under strong selection (ω = 1, V_G/s ≈ 0.7: V_G 4.06 → 1.09 among the survivors) select
still commutes (−0.2%) but reproduce misses by +9.5% for `:multilocus` and +9.6% for
`:rare_alleles`: the kernel restores more variance than the Gaussian closures predict. This
localisation is new: the scripts compare only whole runs at equilibrium. It is a finding
about the closures (one generation from linkage equilibrium), not a change to them.


## The lax version and Lean

Lean proves the exact version: `stochastic_lowerings_agree` (Lowering.lean) says that two
lowerings into the Kleisli category of `PMF` that agree on every step agree on every life
cycle, and `lower_conj` (Plasticity.lean) says the same for lowerings that agree after a
map of states κ (a naturality square with σ = κ). The harness tests a lax version: each
step's square commutes up to an error ε_s, measured on the laws of the statistic. If every
step's image is L_s-Lipschitz for that measure (a step does not amplify a discrepancy by
more than L_s), a composite of n steps commutes up to

    ε(f) ≤ Σ_s ε_s · Π_{s' after s} L_{s'},

by telescoping: replace one step's image at a time. For the exact steps (develop, record,
rescale, mutate's expectation, swap, neutral migration) ε_s = 0, which is the Lean case;
selection's ε is the closure's error (Gaussian truncation, the rare-alleles closure's
fourth cumulant); reproduce's is the LD the closure forgets (zero for the multilocus model
from linkage equilibrium, all of it for the aggregate recursion, as the per-step demo
shows). What the scripts and the harness measure is ε(f) on whole runs; `naturality_steps`
measures ε_s.

Claims that would need Lean lemmas:

1. the telescoping bound above, for Kleisli morphisms of `PMF` with a metric on laws (for
   example total variation, or a Wasserstein distance on a real statistic), from per-step
   bounds and Lipschitz constants;
2. that the kernel's develop, record, rescale and swap images commute exactly with the
   moments statistic (σ = sample mean and variance), and mutate in expectation
   (`Mutation.lean` has the expectation; the square is not stated);
3. the reductions declared below that Lean does not have: the n = 1 and independent-trait
   embeddings (two copies of the 1 × 1 recursion), the projection of the records, one deme
   as a metapopulation (`Demes.lean` has the products, not `p ⊗ demes(1)` = `World(p)`),
   the haploid reduction of the multilocus model to the aggregate recursion, and the
   spectrum's limits;
4. the order of closures as an order: that each arrow's reduction composes (rare local →
   rare alleles → multilocus → aggregate under the conjunction of their conditions).

## Model morphisms and the order of closures

`model_morphisms()` declares the maps between models, and `QTLSim.closure_order()` the
arrows between the closures of the moments lowering. Each [`QTLSim.ModelMorphism`](@ref)
states what it maps, under which condition its reductions hold, and which Lean theorem
(if any) proves them, and lists [`QTLSim.Reduction`](@ref)s on small cases: bitwise where
the code computes the same numbers, to a stated tolerance otherwise. A testset iterates
over them ("declared model morphisms: generated reduction tests"); the hand tests they
repeat (`covers`) stay.

| morphism | kind | condition | reductions (backend, exactness) | Lean |
|:--|:--|:--|:--|:--|
| one trait ↪ two independent traits | embedding | multilocus moments, no drift | G₁₁ = V_G at a fixed optimum and under a move along trait 1; the n-trait recursion with 1 × 1 matrices is the one-trait recursion (moments, 1e-12, 1e-10) | — |
| n traits → trait 1 | projection | always | pheno_mean, geno_var are trait 1's (kernel 1e-12; moments bitwise) | — |
| k demes → one deme | restriction | no migration, equal parameters | `p ⊗ demes(1)` is `World(p)` (kernel, bitwise); three uncoupled demes, migration at rate 0, `between` for one deme (moments, bitwise) | `demes_iterate`, `tensor_iterate` |
| plastic at k → plain at (1 − φ)k | conjugacy | τ = 0, no cost, continuous move | survivors and V_G (kernel, same seed, bitwise); mean fitness and lag (kernel 1e-10); V_G, survivors, lag (moments 1e-12) | `lower_conj`, `plastic_run_conj`, `plastic_statistic_law` |
| `:rare_local` → `:rare_alleles` | refinement | one deme, a common optimum, no migration | four cases (moments, 1e-12; in fact bitwise) | `local_fixed_point_zero`, `gain_const`, `gain_zero_weights` |
| `:rare_alleles` → `:multilocus` | refinement | α² → 0 (a limit) | equilibria at μ = 1 within 1e-3 | `rare_gaussian_limit` |
| `:multilocus` → `:aggregate` | refinement | haploid, no drift | one locus and ten (clones) (moments, 1e-12) | — |
| `:spectrum` → `:rare_alleles` | refinement | rare mutations, fixed optimum (approximate) | stationary V_G within 5% of `house_of_gauss_vg`; neutral 2N·V_m within 1% | — |

8 declarations, 21 reductions, all passing (generated testset). The hand tests they repeat
(the n = 1 recursion, the trait-1 columns, `between` for one deme, the plastic conjugacy on
one seed, the `:rare_local` identities, the spectrum's limits) are kept.


### `:auto` from the declarations

Each closure ([`QTLSim.ClosureDecl`](@ref), `QTLSim.closures()`) declares its rank in the
order, its domain (where `run_moments` accepts it) and its automatic domain (where
`loci = :auto` may choose it), as predicates on a [`QTLSim.ClosureContext`](@ref): the
model, the number of demes, the ploidy, whether the `mutate` steps carry one law of
mutation effects, and `between`. `loci = :auto` is [`QTLSim.auto_closure`](@ref): the
finest closure whose automatic domain contains the context.

| closure | rank | domain | automatic domain |
|:--|:-:|:--|:--|
| `:aggregate` | 1 | every model, no `between` | never |
| `:multilocus` | 2 | everything (`between` for diploids) | always (the fallback) |
| `:rare_alleles` | 3 | one trait or several (no reaction norm), one law of mutation effects | its domain |
| `:rare_local` | 4 | as `:rare_alleles` | several traits, or one trait in several demes |
| `:spectrum` | 5 | one diploid trait, one deme, one law of mutation effects, no `between` | never |

A test enumerates 72 contexts (four models, both ploidies, 1, 2 or 4 demes, no, equal or
different mutation fields on the `mutate` steps) and checks that the declarative rule gives
exactly what the hand-written rule gave; it does in all 72. `run_moments` resolves `:auto`
without `between`, as before: with `between` in the context the declarations would pick
`:multilocus` for haploid metapopulations, which the lowering then rejects anyway (a first
version passed `between` and differed from the old rule in exactly those 18 rejected
cases). The one difference is outside the lowering and is kept as
found: `QTLSim.moments_loci(p)` (which `multilocus_equilibrium` uses) returns
`:rare_local` for a heritable reaction norm, whose moments `run_moments` runs with
`:multilocus`, as the declarations do. Another test checks the declared domains against
`run_moments`' guards over 240 combinations: they agree but for one gap, also kept as
found: several traits in several demes accept `loci = :spectrum` (the guard is only on
the one-deme path) and silently run the rare-alleles closure.

## The support matrix and `explain`

[`support_matrix`](@ref) probes 20 features × 3 models × 7 columns (the kernel, the
generated SLiM model, the moments lowering with each closure) on tiny populations, so a
cell is *unsupported* exactly where the code's own guard throws an `ArgumentError`; a cell
that runs is *supported* when the declared evidence names a validation script that
compares it with another backend, *untested* otherwise. [`support_markdown`](@ref) writes
the [support matrix](support.md) page (`scripts/support_matrix.jl`; a test keeps it
current). [`explain`](@ref) answers the same question for one parameter set: the model,
each closure's domain and automatic domain with the declarations' reasons, the closure
`:auto` picks, and a probe of each backend.

The first generation of the matrix (54 supported, 223 untested, 134 unsupported, 6 not
applicable) found 3 cells that fail with an error other than an `ArgumentError` (marked
`!`): `slim_script` for two traits with `Plasticity` set, which `checkparams` rejects
("Plasticity with NoTraits > 1 needs ReactionNorm"), fails inside the Eidos splice of the
fecundity images ("the Eidos image of select has changed") before any check of the
parameters. Kept as found (reported, not fixed).

## What is left

- C7: diploid → haploid clones as a declared morphism; generating the support matrix's
  evidence from the naturality table rather than declaring it by hand; per-step squares
  for several traits and several demes; making the `ArgumentError` guards themselves read
  the closures' domains (they are only checked against them).
- C2: every family has at least one ported cell; not ported are `coalescence.jl`, the
  pre-specified 40,000-run `slim/three_way_confirm.jl`, `published_sweep.jl` (the paper's
  table as a backend), and recorded failures that compare a run with a function of itself
  (`reaction_norm.jl` C2).
- F11: tests split by topic; the categorical layer as its own package.
