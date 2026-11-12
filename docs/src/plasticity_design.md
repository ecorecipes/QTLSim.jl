# Phenotypic plasticity: a reaction norm through the same life cycle

Plasticity is the standard reason a population might keep up with an optimum that moves
faster than its genes can follow. This note records how a linear reaction norm (Chevin,
Lande and Mace 2010; Lande 2009) went through QTLSim's life cycle and its three lowerings,
what it is in the terms of [Several traits](traits_design.md), and what is not supported.

## The model

The environment `ε` is measured by the optimum it sets: `θ = B·ε`, with `B`
(`OptimumSlope`, 1 by default) the environmental sensitivity of selection. An individual's
phenotype is

    z = a + b·ε_cue + e,      ε_cue = θ_{t−τ}/B,      e ~ N(0, V_E),

with elevation `a` (its breeding value in the reference environment `ε = 0`, where the
optimum is the origin of the trait's scale), slope `b`, and a cue that is the environment
`τ = CueLag` generations before selection. `τ = 0` is a perfect cue, development and
selection in the same environment. The environment before the first generation is taken
to have been the first one seen. A cost of plasticity is stabilising selection on the
slope around 0: fitness is multiplied by `exp(−b²/(2ω_b²))`, `ω_b = PlasticityCostSD`
(`Inf` by default: no cost).

Only the ratio `b/B` matters for a fixed slope: the plastic part of the phenotype is
`(b/B)·θ_{t−τ}`, the fraction of the optimum's displacement that is tracked without any
genetic change. `B` is kept as a parameter so that slopes and costs can be given in the
units of the literature.

Two models:

- **A fixed slope** (`Plasticity = b`, one trait). Every individual has the same slope, so
  the plastic part is common to all and `a = g`, QTSim's genotype. This is a parameter of
  one step, not a new model.
- **A heritable reaction norm** (`ReactionNorm = true`, `NoTraits = 2`). The two heritable
  values of the model of several traits are the elevation, `a = g₁`, and the deviation of
  the slope from `Plasticity`, `b = Plasticity + g₂`. Mutation is the pleiotropic mutation
  of two traits (`MutCov` is the mutational covariance of elevation and slope). Selection
  acts on the one phenotype `z`, with width `ω`, and the cost on the individual's slope.

Everything is off by default: with `Plasticity = 0` and `ReactionNorm = false` the kernel's
output is bit for bit what it was, a generation allocates nothing, and `slim_script`
writes the same file (the golden-file test).

## Where plasticity sits in the life cycle

QTSim's generation is `develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅ mutate ⋅ swap`, and
the optimum of generation `t` is set in `select`, after `develop` has formed the
phenotypes (QTSim's `PopFitness`: phenotypes, then `SetFitnessMean`, then fitness). The
environment a perfect cue reports is therefore not known when `develop` runs. So:

- `develop` is unchanged: it gives the phenotype in the reference environment, `g + e`.
- `select` moves the optimum, **then expresses the norm**, then computes fitness. The
  plasticity fields (`Plasticity`, `OptimumSlope`, `CueLag`, `PlasticityCostSD`) are
  parameters of `select` (`step_schema(OneTraitModel())[:select]`), so
  `withparams(L.select; Plasticity = 0.5)` is a step with its own slope, and demes can
  have their own. `record` then records the expressed phenotypes.

The alternative, moving the optimum in `develop`, would keep "development" in `develop`
but would change which step owns the move (the move's fields are `select`'s, a `select`
step can carry its own move, and QTSim's `:phenotype` and `:fitness` moves need the
phenotypes). It was not done. A lagged cue (`τ ≥ 1`) could have been expressed in
`develop`; it is expressed in `select` too, so that there is one place.

No new generator was needed, and for a fixed slope no new object image either: the
syntax, the stages and the ACSet schema are unchanged.

## Is the heritable norm a third model?

Yes, a small one: `ReactionNormModel()` (`model_of(p)` when `ReactionNorm` is set). It is
best described against `TraitModel(2)`:

| | `TraitModel(2)` | `ReactionNormModel()` |
|---|---|---|
| object image (an individual) | `g, e, z ∈ ℝ²` | `g ∈ ℝ²`, `e, z ∈ ℝ` |
| `develop` | `z = g + e` | the same image (the pair's sums); `z` is formed in `select` |
| `select` | `W = c·exp(−½ (z−θ)ᵀΩ⁻¹(z−θ))`, `θ ∈ ℝ²` | `z = (1, ε)·(a, b) + e`, `W = c·exp(−(z−θ)²/2ω²)·exp(−b²/2ω_b²)`, `θ ∈ ℝ` |
| `reproduce`, `rescale`, `mutate`, `swap`, `migrate` | | the same images |
| `record` | trait means, `G`, `P` | the same code; the "traits" recorded are `(z, b)` |

So it is `TraitModel(2)` with **a linear map before selection**: the pair `(a, b)` is read
through `(a, b) ↦ a + b·ε`, a map that changes every generation with the environment, and
selection is one-dimensional on the image. It is not `TraitModel(2)` with a clever `Ω`:
the fitness function on the pair, `exp(−½ (g − g*)ᵀ Q (g − g*))` with
`Q = v vᵀ/(ω² + V_E) + e₂e₂ᵀ/ω_b²`, `v = (1, ε)`, is singular without a cost (rank one), has
no fixed optimum and changes with `ε`; and the environmental deviation is one number, not
two. But all the genetics is shared: the block pool with two values per locus,
recombination, pleiotropic mutation, migration. In the kernel the model's table is
`TraitModel(2)`'s with `select` replaced by `QTLSim.select_norm!`; in the moments lowering
`trait_moment_images` with `select` and `record` replaced (the haplotype moments of the
multilocus model take selection through the matrix `K` with `G_s = G − G K G`, where
`TraitModel` uses `(Ω + P)⁻¹`). Its step parameters are typed as one trait's for `select`
(`ω` and `θ` numbers, the plasticity fields) and as two traits' for `mutate` (`MutCorr`,
`MutCov`).

In the language of the north star: the object image did not change from `TraitModel(2)`'s
for the heritable state, and the new ingredient is a morphism of the environment,
`ℝ² × Env → ℝ`, composed before the image of `select`. A general framework would want
that map as a first-class part of a model (a genotype–phenotype map parametrised by the
environment), of which QTSim's identity and this linear map are instances.

## The lowerings

| | fixed slope | heritable norm |
|---|---|---|
| kernel | `select!` expresses the norm (`QTLSim.express_plasticity!`) and applies the cost | `QTLSim.select_norm!` |
| moments | the lag is `θ − ḡ − (b/B)θ_{t−τ}`; `W̄` times the cost | `QTLSim.reaction_norm_response` on `(ā, b̄)` and `G` (2×2), multilocus Gaussian |
| Eidos | spliced into `select` when a deme or a step uses it (`plasticity`, `optSlope`, `cueLag`, `plastCostSD`; the deme's cues are `CUES` on the subpopulation) | not implemented |

The Eidos splice follows the pattern of the fluctuating moves and of fecundity selection:
the default output is byte-identical, and `validate/slim/compare_plasticity.jl` compares
the kernel with the generated SLiM model.

## Theory

With a fixed slope and `φ = b/B` the mean genotype follows
`ḡ′ = ḡ + β (θ_t − ḡ − φ θ_{t−τ})`, `β = V_G/(V_G + V_s)`: the recursion without plasticity
for the effective optimum `θ_t − φ θ_{t−τ}`. All of the following are derived from that
recursion (and checked symbolically in `validate/symbolic/check_plasticity.jl`), not
transcribed from the literature; the equations of Chevin, Lande and Mace (2010) are
images in the published article and could not be read by the tools used here, so only
their statements in the text were checked against the source (that the critical rate
increases as `|B − b|` falls, that a cost lowers the maximum growth rate, and that above a
threshold cost an intermediate plasticity maximises the critical rate).

- **An exact correspondence.** With a perfect cue (`τ = 0`) and no cost,
  `z − θ = g + e − (1 − φ)θ`: fitness depends on the genotype and the optimum only through
  `(1 − φ)θ`. So a plastic population whose optimum moves `k` per generation is, in
  distribution and not only at steady state, a population without plasticity whose optimum
  moves `(1 − φ)k`: every statistic of the genotypes, the survivors and the time to
  extinction is that of the slower move (tested in the kernel: the same survivors and
  genetic variances for the same seed). The scaling of the sustainable rate by
  `1/(1 − φ)` is therefore exact in this model, whatever the genetics, and the experiment
  of `validate/plasticity_rates.jl` confirms it as a check of the implementation. What can
  depart from it is a cost, a delayed cue of a fluctuating optimum, and a heritable slope.
- **Steady lag** behind an optimum moving `k` per generation: `(1 − φ) k (V_G + V_s)/V_G`
  ([`steady_lag`](@ref)), whatever the cue's delay: the genotype has to supply only
  `(1 − φ) k` per generation.
- **Critical rate**: [`critical_rate`](@ref) (Bürger and Lynch's, in the paper's table
  form) divided by `|1 − φ|`, with `ln B₀` reduced by the cost `b²/(2ω_b²)`.
  [`sustainable_rate`](@ref) is the same threshold for the package's own demography: the
  rate at which the steady lag brings the mean fitness to `1/R` (Beverton–Holt) or to
  `Ploidy/N` (QTSim's fixed size). With a cost an intermediate slope can be best.
- **Fluctuating optima**: the lag is `L_t − φ L_{t−τ}` with `L` the lag without plasticity.
  For red noise its variance is `(1 + φ²)Γ(0) − 2φΓ(τ)` ([`red_noise_lag_var`](@ref),
  [`red_noise_lag_autocov`](@ref)): a perfect cue scales the variance by `(1 − φ)²`, an old
  cue of a fast environment makes plasticity add noise. For a cycle the squared amplitude
  is multiplied by `1 + φ² − 2φ cos(2πτ/P)` ([`expected_lag2`](@ref)). The cost adds
  `b²/(2ω_b²)` to [`lag_load`](@ref).
- **A heritable norm** ([`reaction_norm_response`](@ref), [`reaction_norm_lag`](@ref)):
  without a cost `Δ(ā, b̄) = G v λ/s` with `s = ω² + V_E + vᵀGv`. The slope responds in
  proportion to `ε`, so plasticity evolves only away from the reference environment; there
  the slope's variance inflates the phenotypic variance by `G_bb ε²` (a variance load) and
  is used up by selection. Under a sustained move the mean slope approaches `B` as
  `d b̄/dt ≈ (B − b̄)/t`, and the lag falls with it: with heritable, cost-free plasticity
  there is no critical rate in this model, only a transient. With a cost, in a constant
  environment the only equilibrium has `b̄ = 0` and the elevation at the optimum (Lande's
  genetic assimilation); under a sustained move the slope settles below `B`.

## What is rejected, and what is not covered

- `MoveType = :phenotype` and `:fitness` with plasticity: QTSim re-anchors the optimum to
  the current mean phenotype every generation (`θ = z̄ + lead`), so a plastic shift of the
  mean is carried along by the optimum and has no meaning. `checkparams` rejects them
  (unless the optimum never moves). Plasticity needs an optimum with its own trajectory:
  `:continuous`, `:red_noise`, `:periodic`, a step's fixed `θ`, or none.
- `Plasticity ≠ 0` with `NoTraits > 1` needs `ReactionNorm`; `ReactionNorm` needs
  `NoTraits = 2`, the default optimum direction and no `SelCorr`/`SelCov`.
- The heritable norm has no Eidos lowering (`slim_script` says so), no rare-alleles
  closure in the moments (the fitness function on the pair is singular without a cost, so
  the closure's rate is undefined; `loci` is `:multilocus` or `:aggregate`), and its ACSet
  is that of two traits.
- With red noise the moments lowering follows the trend only, as without plasticity; the
  lag's variance is `red_noise_lag_var`'s.
- `StandardisePhenotype`, fecundity selection, density regulation, demes and migration all
  work with a fixed slope (the shift is common to a deme; a migrant expresses the norm in
  the deme it develops in). They are tested lightly, not validated separately.
- The measure of the environment is the optimum itself. An environmental variable that
  moves the optimum nonlinearly, or cues that are noisy for reasons other than their age,
  are not modelled.

## Validation

- `test/runtests.jl`: the default path bit for bit (a slope with the optimum at the origin
  adds exactly zero), no allocation, the kernel's expression of the norm, reductions at
  `b = 0` (`===`), kernel against moments and theory, steps and demes with their own
  slopes, the Eidos splice, and the heritable norm against Monte Carlo and the moments.
- `validate/symbolic/check_plasticity.jl`: the algebra above.
- `validate/slim/compare_plasticity.jl`: kernel against the generated SLiM model.
- `validate/plasticity_rates.jl`: the pre-specified experiment on sustainable rates
  (`plasticity_rates.log`, 12,320 kernel runs: two published-like parameter sets with 100
  and 500 loci, fixed size and Beverton–Holt R = 2, b/B = 0, 0.25, 0.5, 0.75, with and
  without a cost ω_b = 0.5, 40 burn-ins of 10,000 generations per cell, censored 2000
  generations after the move). Without a cost the half-survival rate rose by 1.33–1.35,
  1.98–2.03 and 3.99–4.07 times (95% bootstrap intervals within 1.32–1.38, 1.95–2.07 and
  3.93–4.13), against the predicted 1.33, 2 and 4: all 12 cells pass the 15% criterion,
  as they must, since the correspondence above is exact. With the cost the gain shrinks:
  at fixed size 1.24–1.26, 1.46–1.55 and 2.11–2.13 times; under Beverton–Holt 1.03–1.06
  times at b/B = 0.25, below 0.63 (every run died on the grid) at 0.5, and extinction
  before the move at 0.75. The deterministic prediction of the cost's effect failed 4 of
  12 cells. The cost also lowers V_G before the move (more drift among fewer survivors,
  and under Beverton–Holt a population of about 340), but that is not why: in
  `validate/plasticity_erosion.jl` (pre-specified, fresh seeds) the pre-move loss accounts
  for 3–9% of the cost's effect on k50; the cost acts after the move, 2–7 times more
  strongly than `sustainable_rate`'s fitness-margin term predicts. `sustainable_rate`
  therefore mispredicts a cost of plasticity. See the script's header. The explanation is a
  feedback `sustainable_rate` leaves out: the lag and the cost both lower the number of
  breeders, drift then lowers `V_G`, and the lag grows. [`drift_critical_rate`](@ref) is
  the critical rate of that feedback; by the envelope theorem a log cost `c` lowers it by
  the factor `1 − c/M`, `M = λ*²/s` (linear, not the square root). In
  `validate/critical_rate.jl` (pre-specified, fresh seeds) its cost ratio is within 0.054
  of the kernel's in 24 of 24 cells and the linear form fits better in 8 of 8 groups (the
  observed `M` 9–22% above `λ*²/s`); [`half_survival_rate`](@ref) predicts the
  half-survival rate itself (26 of 32 cells within 0.75–1.33).
- `validate/reaction_norm.jl` (`reaction_norm.log`): the heritable norm, kernel against the
  moments lowering, 40 runs; 28 of 30 criteria pass. The mean slope (0.2 → 0.96 without a
  cost, → 0.56 with ω_b = 1) is within 0.07 of the moments at every time, the lag within
  0.01 phenotypic SD, G_bb within a factor 1.5. The two failures are the quasi-steady lag
  formula late in the move (1000 and 1500 generations), where the lag has fallen to 0.003
  and below, and the formula's ratio is noise.

## For Lean

Claims the Julia code relies on that should get Lean lemmas (none is in `lean/` yet):

1. *Steady lag with plasticity*: for `ḡ′ = ḡ + β (θ_t − ḡ − φ θ_{t−τ})` and `θ_t = k t`,
   the lag converges to `(1 − φ) k/β` for every `τ` (from `FluctuatingOptimum.trend_follows`
   by superposition: the effective optimum is a trend of slope `(1 − φ) k` plus a constant).
2. *Critical rate*: `W̄(λ) > W_min` iff `|k| < β √(2 s ln(W̄(0)/W_min))/|1 − φ|` at the steady
   lag, and its reduction to the rate without plasticity at `φ = 0`.
3. *Autocovariance of the lag behind AR(1) noise*: `Γ(h + 1) = a Γ(h) − (1 − ρ) ρ^h c₀`,
   its closed form, and `Var(L_t − φ L_{t−τ}) = (1 + φ²)Γ(0) − 2φΓ(τ)` (extends
   `red_noise_lag_variance`).
4. *A cycle with a delayed cue*: the amplitude factor `|1 − φ e^{−iντ}|`.
5. *Selection on a reaction norm*: `(I + GQ)⁻¹G = G − G K G` with `K = Q (I + GQ)⁻¹`, and for
   rank-one `Q` (no cost) the Sherman–Morrison form `G − G v vᵀ G/s` (a case of
   `Traits.posterior_cov` when `Q` is invertible; the singular case is new).
6. *Equilibria of the norm in a constant environment*: with a cost, `Δ = 0` iff `λ = 0` and
   `b̄ = 0`; without, iff `λ = 0`.
7. *The expressed phenotype is exact when the optimum has not moved*: `(b/B)·0 = 0`, so the
   run with a slope equals the run without one before the move (the bitwise test).
8. *The exact correspondence*: with `τ = 0` and no cost the life cycle's kernel with slope `b`
   and optimum sequence `θ_t` equals, as a Markov kernel on genotypes (in `Lowering.lean`'s
   Kleisli category of `PMF`), the kernel without plasticity and optimum `(1 − b/B)θ_t`: the
   two `select` images agree on every population, so by `stochastic_lowerings_agree` the
   generations do. This is the statement behind "rates scale by `1/(1 − b/B)`".
