# Why does the `:spectrum` emulator put k50 too high near the threshold?

Diagnostic study, exploratory (2026-11-14). Script: `validate/spectrum_diagnosis.jl`. Running log,
with each hypothesis written before its test and the per-input output in appendices:
`validate/spectrum_diagnosis.log`. Nothing here is a pre-specified validation, and no default
changed.

## Question

On 160 fresh inputs (E6′ of `validate/emulator_spectrum.jl`, 40 kernel runs each),
`emulate(p; loci = :spectrum)` called 25 inputs alive at +2000 that the kernel called dead, and
1 the other way round. The pooled x50 ratio was 1.66. We asked which link of the emulator's
chain is wrong:

burn-in V_G → V_G during the move → lag θ − z̄ → W̄ (Gaussian) → extinction given W̄

and, if V_G is the wrong link, why.

## Method

- **Inputs.**
  - The 26 discordant E6′ inputs were split in id order: odd positions went to discovery (D, 13
    inputs) and even positions to confirmation (C, 13).
  - Concordant controls were drawn at random, 2 per stratum (`MersenneTwister(20261114)`): one to
    Dctl and one to Cctl (8 each).
  - The list is in `validate/spectrum_diagnosis_inputs.csv`.
  - The kernel runs are E6′'s own runs r = 1…10 (same seeds, same `init = :equilibrium`, same
    1000-generation burn-in), re-run with extra statistics. The extras use no random draws, so the
    runs reproduce E6′'s.
  - Recorded every generation: V_G, V_P, z̄, ḡ, θ, W̄, survivors, and the survivors' mean genotype.
  - Recorded every 10 generations: the genic variance (2 Σ_loci var(allele)), giving LD = V_G − genic,
    plus genotype skew and kurtosis.
  - Recorded every 50 generations: the number of distinct alleles.
- **Substitutions.**
  - The Gaussian lag recursion of the moments lowering, fed:
    - the kernel's mean V_G over the runs alive;
    - the kernel's or the spectrum's pre-move V_G, held constant.
  - In each case, P2000 and the median T then come from `extinction_distribution`.
  - The response link: the realised survivors' genotype shift against V_G·λ/S.
  - The W̄ link: the realised W̄ against c·ω/√S·exp(−λ²/2S).
  - The extinction link: W̄ in a run's last generation.
- **Driven spectrum.** The spectrum's genic variance driven, generation by generation, by each kernel
  run's own lag, S and survivors (no feedback), compared with the kernel's genic variance.
- **Looser linkage.** The kernel and the spectrum with every one of Poisson(20) crossovers applied
  (`VARIANT=many`); QTSim applies one.
- **k50.** On 6 D inputs, kernel k50 from 10 runs at ×0.25/0.35/0.5/0.7/1.0 of the input's rate,
  against the emulator's k50 with and without an LD correction.
- **Prototype P1** (keyword `ldratio` of `AlleleSpectrum`/`spectrum_equilibrium`, default 1, so the
  default is unchanged). r = V_G/genic, taken from `multilocus_equilibrium(p; drift = true,
  loci = :multilocus)`, a model quantity with nothing fitted. The lowering's V_G is r × genic, and
  each mutation's directional selection is r·aλ/S. P1 was confirmed on:
  - the 139 E6′ inputs not used for discovery;
  - 64 fresh published inputs (20 kernel runs each; none in any earlier sample or pool).

## Results

### 1. The lag, W̄ and extinction links are right; the whole bias is in V_G

The measurements cover all 42 diagnostic inputs (D, Dctl, C, Cctl).

| link | measurement | result |
|---|---|---|
| response | realised survivors' shift / V_G·λ/S | 0.981–1.010 |
| W̄ | realised W̄ / Gaussian formula | 1.000 |
| extinction | median W̄ in a run's last generation | 0.0011–0.0059, so Binomial(1000, W̄) < 2 is likely |
| V_G | kernel's mean V_G in the Gaussian recursion | reproduces the kernel's median death time (examples below) |

The V_G substitution gives P2000 = 0 for every discordant input with kernel deaths. Examples of
predicted vs kernel median death times: 253 vs 244, 101 vs 99, 735 vs 720, 574 vs 544, 220 vs 218,
445 vs 443, 348 vs 347. The spectrum's own V_G in the same recursion reproduces the emulator.

The kernel's deaths are deterministic. Each input dies in a narrow window (e.g. 177–238
generations), not by stochastic escape from a steady state, so between-run noise is not the cause.

### 2. V_G has two errors, both from QTSim's tight linkage

**(a) Missing negative LD at the move.** The pre-move V_G, kernel/spectrum, is 0.34–0.98:

- median about 0.76;
- 5000 loci 0.34–0.78.

The kernel's *genic* variance at the move equals or exceeds the spectrum's V_G, so the shortfall
is negative linkage disequilibrium:

- kernel LD/genic is −0.02 to −0.81;
- for example Input-12000: V_G 0.142 = genic 0.379 + LD −0.237, against the spectrum's 0.303.

The spectrum has no LD (V_G = genic). The multilocus haplotype moments do model it (Input-12000:
0.132 = 0.416 − 0.284), but 10–40% too strongly. Their V_G/genic is below the kernel's in 39 of
42 inputs.

**(b) The genic variance falls too late during the move.** In the kernel, the genic variance of
the discordant inputs falls by 5–45% within 100 generations, and by 20–90% within 250.

Driven by each run's own lag and survivors, the spectrum's genic variance at +100 is still
0.88–1.55 of its start (it even rises, the Barton–Turelli rise, where V_G is small). The kernel
is at 0.56–0.98. The spectrum catches up only by about +500. Its decline comes entirely from
Barton's linked-selection term; with `linkage = :none` it hardly falls.

Holding V_G at the kernel's pre-move level (constant) keeps every discordant input alive. The LD
at the move alone is not enough; the decline during the move is needed too.

Changing δ does not fix the timing:

- δ = 1/(2n) (selection's part dropped) collapses far too fast, including in the controls.
- δ = 1/(2n) + d·V_G/S with d = 0.5, 0.25 or 0.1 matches at +250 for d ≈ 0.5 but overshoots by
  +500.
- So the kernel's decline is earlier and more gradual than a term in b² = (λ/S)² can give. The
  shape is wrong, not only the scale.

**Looser linkage removes the discordance.** With Poisson(20) crossovers all applied, all 78 runs
(13 D inputs × 6) survive to +2000, as the emulator predicts. With QTSim's crossover, 128 of the
130 E6′ runs of the same inputs had died. The kernel's V_G/genic rises to 0.71–0.99, and its V_G
no longer collapses.

### 3. How much of the k50 bias each part explains (6 D inputs)

| input | kernel k50 | spectrum / kernel | + LD at kernel's level | P1 (model's LD) |
|---|---|---|---|---|
| QT0500-3/04943 | 0.0406 | 2.09 | 1.17 | 0.93 |
| QT1000-3/09697 | 0.0770 | 1.59 | 1.15 | 1.06 |
| QT1000-5/14110 | 0.0365 | 1.46 | 0.96 | 0.80 |
| QT5000-3/12000 | 0.0577 | 3.42 | 1.71 | 1.50 |
| QT5000-5/12463 | 0.0439 | 1.78 | 1.25 | 1.10 |
| QT5000-5/01925 | 0.0376 | 1.23 | 1.02 | 0.99 |
| geometric mean | | 1.82 | 1.19 | 1.04 |

**Shares of the log k50 bias:**
- **Missing LD at the kernel's own level:** about 71% (mean log 0.599 → 0.173).
- **Too-late genic decline:** the remaining ~29%.

The same split shows in the calls. With the kernel's measured LD ratio, the spectrum still calls
7 of the 25 discordant D and C inputs alive. P1, with the model's stronger LD, fixes all 25, by
compensation. The kernel's k50 is resolved only to about ±17% (rate grid).

### 4. Prototype P1 on inputs not used to develop it (exploratory)

| sample | model | alive / kernel dead | dead / kernel alive | accuracy | pooled x50 ratio |
|---|---|---|---|---|---|
| 139 E6′ inputs (C, Cctl, unselected) | spectrum | 13 | 0 | 0.906 | 1.361 |
| | P1 | 0 | 13 | 0.906 | 0.798 |
| 64 fresh inputs | spectrum | 5 | 1 | 0.906 | 0.935 |
| | P1 | 0 | 1 | 0.984 | 0.954 |

- **E6′ sample:** P1 reverses the bias but overcorrects. Its new errors are in the 100- and
  500-locus strata, where LD is weak and the model's LD is relatively too strong.
- **Fresh sample:** most inputs are far from the threshold (54 of 64 with no run alive), so its
  x50 says little. It confirms the direction of the spectrum's errors (5 against 1), but not the
  predicted ratio above 1.25.

**P1 is not established and should not become a default.**

## Conclusion

The emulator's lag, W̄ and extinction links are right. The bias is entirely in V_G, and it comes
from QTSim's tight linkage (at most one crossover per gamete), which the spectrum closure
represents only through Barton's loose-linkage N_e term:

1. Stabilising selection with tight linkage builds strong negative LD. The kernel's V_G at the
   move is 0.34–0.98 of its genic variance, and the spectrum (V_G = genic) misses this. This is
   the larger part, about 70% of the log k50 bias on the inputs measured.
2. Under the move, the kernel's genic variance falls from the start and gradually. The spectrum's
   linked-selection drift, which grows as (λ/S)² and is damped by δ, starts too late. Retuning δ
   does not reproduce the shape. This is the remaining ~30%.

Both disappear when linkage is loosened, and so does the discordance.

**Confidence.**
- **High** that the error is in V_G and not in the lag, W̄ or extinction links: exact substitution
  on 42 inputs, every discordant one.
- **High** that tight linkage is the cause: the looser-linkage experiment.
- **Moderate** on the 70/30 split: 6 inputs, coarse kernel k50.
- **Low** on any fix.

## Open questions

- A dynamic LD term in the spectrum. For example, carry the multilocus haplotype moments' LD
  classes alongside the spectrum (they already exist in `src/multilocus.jl`), instead of P1's
  constant r. Also: why are those moments' LD 10–40% too strong for QTSim's crossover?
- The early genic decline. The candidates were not tested:
  - the kernel's allele spectrum at the move is richer in intermediate-frequency variance than
    the spectrum's equilibrium shape;
  - linked selection already acts under stabilising selection at λ = 0, where Barton's term is
    zero;
  - multiple alleles per locus.
- Is the 100-locus reverse case (06159: kernel alive, every emulator variant dead) the same
  mechanism in the other direction (V_G under the move underestimated)?

## Mechanics

- CPU time: about 7.5 CPU-hours (rough tally from wall times and thread counts).
  - Kernel runs: 196 + 210 trajectory runs, 78 looser-linkage runs, 300 k50 runs and 1280 fresh
    runs, plus a 382 s free-recombination pilot.
  - Driven-spectrum tests, emulations and tests.
- Caches (committed):
  - `spectrum_diagnosis_k50_kernel.csv`;
  - `spectrum_diagnosis_proto_e6p.csv`;
  - `spectrum_diagnosis_fresh_{inputs.txt,kernel.csv,emulator.csv}`.
- Raw per-run trajectories (74 MB) are in `validate/spectrum_diagnosis_raw/`. They are git-ignored
  and regenerated by the `kernel` mode.
