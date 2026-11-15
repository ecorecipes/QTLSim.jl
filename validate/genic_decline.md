# Why does the genetic variance fall as soon as the optimum starts to move?

Diagnostic study, exploratory (2026-11-14). Script: `validate/genic_decline.jl`. Running log,
with each hypothesis written before its test, and the full per-input output in appendices A–D:
`validate/genic_decline.log`. Nothing in `src/` changed and no default changed.

## Question

On the published protocol (continuous move from generation 10000, N = 1000, QTSim's single
crossover), the kernel's genic variance (2 Σ_l of the variance of the allele values at locus l)
starts to fall as soon as the optimum moves. By +100 it is at 0.56–0.98 of its pre-move value.
`:spectrum` does not show this, even with `ld = :dynamic` and even when driven by each kernel run's
own lag, S and survivors (validate/spectrum_diagnosis.md 2(b)). This early fall is about 40% of
the emulator's remaining log k50 bias. What is the mechanism?

## Method

**Exact accounting over the stages of one generation.** Genic variance G is computed with
population denominators at every stage:

    adults A → π-weighted adults → W-weighted adults → survivors S → offspring O → after mutation (= next A)

- π_i = E[W | g_i] is the survival probability given the genotype, i.e. deterministic selection
  on genotypes.
- W adds the environmental noise and S the Bernoulli survival. O adds the sampling of parents and
  gametes.
- The expected noise is computed alongside: the first-order Bernoulli(π) formula, and
  −G_S/(2N) for reproduction.

The selection step G_π − G_A is split per locus by Price's equation:

- **Squared mean response**, −2Σ(Δx̄_l)². This is split into:
  - own2LD: the locus's own response, scaled by the systematic LD;
  - hitch: its random association with the rest of the genome, i.e. Robertson-style hitchhiking
    of the mean;
  - a nonlinear rest.
- **Within-locus term**, 2Σ Cov(π/π̄, (x − x̄_l)²). This is split into:
  - own3 = 2β Σκ3_l: the locus's own third moment, the Barton–Turelli term;
  - other3 = 2β Σ[Cov(g, d_l²) − κ3_l]: the association between a locus's squared deviation and
    the rest of the genotype. Its regression across loci separates a κ3-proportional part from a
    sign-symmetric part (sym3);
  - quad: the curvature of π in g, including the locus's own kurtosis.

**Bitwise check.** The instrumentation draws nothing. V_G in every instrumented generation is
bitwise that of `simulate` with the same seed (the `check` mode, baseline and many-crossover runs).

**Inputs.** From validate/spectrum_diagnosis_inputs.csv:

- **Discovery:** D (13) + Dctl (8).
- **Confirmation:** C (13) + Cctl (8), analysed only after the decomposition was fixed.
- **Runs:** 8 runs per input, the E6′ seeds r = 1..8, `init = :equilibrium` with 1000
  generations of burn-in, instrumented from −200 to +299.
- **Census:** a census of distinct alleles every 25 generations.

**Counterfactuals.** Four runs each on six D inputs (CF: 04943, 00854, 14110, 07594, 12000, 05562):

- many crossovers: every one of Poisson(20) applied;
- `Selection = :fecundity`;
- PopSize × 4;
- `sites10`: NoQTL × 10 at a tenth of the rate per locus, so the same V_m, effects and map, but
  closer to infinite sites. Run for L ≤ 1000 only.

**Spectrum partition.** The spectrum was driven by the kernel runs, as in the diagnosis. Its
first-order selection change, Barton's linked drift −g·c_d, its drift −g/(2n) and its mutation
were compared with the kernel's terms.

## Results

### 1. Which stage carries the decline

| term (cumulative 0..99, relative to genic before the move) | discovery D + Dctl | confirmation C + Cctl |
|---|---|---|
| sym3 (sign-symmetric background term) | 76–145% of the total change in inputs falling by ≥ 5% by +100 (01192, −0.04, has 2.7×; e.g. 12000 −0.29 of −0.33, 05562 −0.24 of −0.27, 00854 −0.35 of −0.24) | 88–115% in all 13 falling inputs (04095 −0.58 of −0.51, 01803 −0.58 of −0.64, 09326 −0.41 of −0.45) |
| hitch (random linked response of the mean: Robertson, draft) | 1–9% (largest 14771 −0.08 of −0.85) | 1–11% |
| own2LD (own response squared) | ≤ 5% | ≤ 5% |
| quad (curvature, own kurtosis) | −0.03..+0.03 | similar |
| own3 + κ3-proportional other3 | net ≈ 0, or positive where rare favourable alleles carry variance | same |
| noise (environment, Bernoulli, reproduction) | at their expectations given the realised survivors (±0.01 of b²V_G per generation) | same |
| mutation | constant ≈ V_m | same |

In the falling inputs the loss per generation is
0.01–0.03 × β²V_G·G in generations 0..9, rising to 0.02–0.05 (10..24), 0.03–0.07 (25..49) and
0.07–0.18 (50..99). It grows 3–5× over the first 100 generations, so it is not a constant
coefficient on b²V_G. The association behind it, Σ_l[Cov(g, d_l²) − κ3_l], is zero before the
move and builds over 25–75 generations: it is created by the move, not merely switched on by β.

### 2. Counterfactuals (CF; genic(+100)/genic(pre) and sym3 over 0..99)

| input | base | many crossovers | fecundity | sites10 | 4N | sym3: base / many / fec / sites10 / 4N |
|---|---|---|---|---|---|---|
| 04943 | 0.56 | 1.71 | 0.65 | 0.61 | 1.11 | −0.54 / −0.03 / −0.56 / −0.58 / −0.12 |
| 00854 | 0.69 | 1.61 | 0.87 | 0.68 | 1.04 | −0.40 / −0.01 / −0.24 / −0.40 / −0.05 |
| 14110 | 0.78 | 1.00 | 0.88 | 0.87 | 0.98 | −0.18 / 0.00 / −0.10 / −0.15 / −0.00 |
| 07594 | 0.77 | 1.19 | 0.75 | 0.73 | 0.99 | −0.27 / −0.01 / −0.30 / −0.29 / +0.00 |
| 12000 | 0.71 | 0.98 | 0.70 | — | 0.97 | −0.26 / −0.01 / −0.29 / — / −0.02 |
| 05562 | 0.78 | 0.99 | 0.75 | — | 0.98 | −0.21 / −0.00 / −0.27 / — / −0.02 |

- **Many crossovers:** removes sym3 and the decline.
  - Where rare favourable alleles carry variance, genic variance then rises (own3 up to +0.96 by
    +100), the Barton–Turelli rise the spectrum predicts.
  - With QTSim's crossover the κ3-proportional part of other3 cancels 40–70% of that rise.
- **Fecundity selection:** sym3 is unchanged in 4 inputs and about halved in 2.
- **sites10:** unchanged, so the continuum of alleles is not the cause.
- **4N:** sym3 is 5–15× smaller. This fits LD generated by drift (∝ 1/n) as the seed. Caveat: 4N
  also changes the equilibrium, the lag and the survivors.

### 3. Why `:spectrum` misses it

Driven by the same runs (CF inputs, at +100):

- Barton's term ≈ the kernel's hitch, drift ≈ drift and mutation ≈ mutation.
- The selection term differs:
  - kernel: −0.17 to −0.37;
  - spectrum with `ld = :none`: +0.00 to +0.50;
  - spectrum with `ld = :dynamic`: −0.00 to +0.14.

The spectrum's sites are in linkage equilibrium with their backgrounds, apart from the random
variance of Barton's term. So it has no counterpart of the systematic association between a
site's deviation and its background.

## Mechanism

Under the move, deviant alleles at a locus, of **either sign**, increasingly sit on genomic
backgrounds that lag behind the population. Directional selection on those backgrounds then
removes them. Within-locus variance therefore falls while the locus's mean hardly moves, which
is why the squared-response hitchhiking term is small.

The association is built by the move itself, under QTSim's tight linkage (at most one crossover
per gamete):

1. **Seeding.** Drift creates random LD between each allele and its linked background. In a
   finite population this happens at rate ∝ 1/n.
2. **Persistence.** With little recombination, an association persists for many generations.
3. **Systematic sorting.** While it persists, directional selection keeps moving the allele's
   frequency the same way. Minority alleles on poor backgrounds keep declining, while those on
   good backgrounds rise. Once a good background becomes common, its alleles are no longer
   deviant.
4. **Accumulated association.** Sorting leaves minority status correlated with poor backgrounds.
5. **Loss.** Selection acting on that correlation removes within-locus variance, at first order
   in β.

This is Robertson's (1961) cumulative effect of linked selection, the process Barton's (1995)
N_e term summarises. But it is its autocorrelated part, the 1/ρ² part. In a per-generation
accounting that part appears as a deterministic first-order term (other3/sym3), not as Var(Δp)
(hitch, the 1/ρ part). Barton's term is the steady state of the two, which the spectrum applies
at once with δ = 1/(2n) + V_G/S. The kernel's version:

- starts from zero at the move;
- builds as the associations accumulate;
- with ρ ≈ 0 for most pairs under the single crossover, exceeds the δ-capped steady state within
  a few tens of generations.

Against the hypotheses stated beforehand:

- (a)/(e) linked selection: confirmed in nature, but in a different term than predicted (sym3,
  not hitch).
- (b) continuum of alleles: rejected (sites10).
- (c) sweeps of pre-existing favourable alleles via κ3: rejected as the cause of the fall. κ3 terms
  push the variance up and are cancelled.
- (d) survivors: rejected; noise is at expectation given realised n.
- (f): nothing else needed.

## Exploratory theory term

Barton's term was rewritten as the variance increment of the cumulative linked change of p,
starting at the move. Nothing was fitted.

    c(t) = Σ_c ω_c b_t [2R_c(t) + b_t Q_c(t)],   Q_c ← λ_c² Q_c + V/(4n),   R_c ← λ_c (R_c + b Q_c)

Its steady state is Barton's form. It was tested driven on the 13 D inputs (discovery, not a
confirmation) with two choices of λ_c:

- **λ_c = (1 − ρ_c)(1 − 1/(2n)):** far too strong (median log(model/kernel) −0.92 at +100,
  collapse by +200). A fully linked association keeps pushing the site for ~2n generations, while
  in the kernel the background's sweep ends sooner.
- **λ_c × (1 − V/S), the spectrum's selection damping:** identical to Barton's term (+0.32 at +100,
  too weak).

The kernel lies between these bounds. The missing quantity is how fast an association decays
under directional selection, i.e. the turnover of selected backgrounds, which saturates the
linear Δp = bD. No term is proposed, and the bracket is the result. A pre-specified test of any
future term would be a separate step.

## Coordinator's hypotheses (from a parallel closure study)

- **Σκ4/Σμ4 ≈ 0 and Σκ3 at noise level at the move: broadly confirmed.**
  - Σκ4/Σμ4 is −0.10 to +0.12: +0.04..+0.12 at 500/1000 loci, −0.06..−0.10 at 5000 loci.
  - Σκ3 in the favoured direction is ≤ 0.004 of genic.
  - But under the move the kernel's own3 (2βΣκ3) is clearly positive in the 500–1000-locus inputs
    (+0.11..+0.43 over 0..99). Under loose linkage it drives a real rise (counterfactual many).
  - So the spectrum's rise is not spurious as a process. With QTSim's crossover it is cancelled by
    linked selection.
- **Loss ≈ 0.05–0.11 b²V_G per generation, a whole-genome linked-selection coupling: scale and kind
  confirmed, constancy not.**
  - Selection is the stage carrying it; the noise stages are at expectation.
  - The coefficient grows 3–5× over 10–100 generations, as an association that has to accumulate
    does.
  - It is seeded by drift-generated LD: 4N cuts it 5–15×.

## Confidence

- **High** that the decline is deterministic selection acting through the association of
  within-locus deviations with the background (sym3). The accounting is exact, it holds on 42
  inputs, and it was confirmed on 13 independent falling inputs.
- **High** that it needs tight linkage (many crossovers remove it) and not the continuum of
  alleles (sites10).
- **Moderate** on the interpretation as the autocorrelated part of Robertson/Barton linked
  selection. This is supported by the time course and by the 4N and many-crossover results, but
  not derived quantitatively.
- **Low** on any closure term. The derived transient brackets the kernel but does not match it.

## Open questions

- The decay rate of an allele–background association under directional selection with tight
  linkage: the saturation of the linear Δp = bD by the background's own sweep. This is the
  δ problem of spectrum_ld_dev.log step 6 and of the diagnosis, now located precisely.
- Whether a haplotype-block description (block values, their sweep times) gives that rate
  without fitting.
- The share of the k50 bias was not re-measured here; the 40% figure is the diagnosis's.

## Mechanics

- **CPU:** about 9 CPU-hours.
  - Instrumented kernel runs: 2 × 168 discovery runs, 168 confirmation runs, 88 counterfactual runs.
  - Driven-spectrum variants and analyses.
- **Raw per-run accounting** goes to `validate/genic_decline_raw/` (git-ignored; regenerated by the
  `kernel` mode, deterministically from the seeds).
