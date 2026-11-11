# Parameters

A simulation is described by a [`Params`](@ref) value. Its defaults are QTSim's
(`SetDefOptions` in `reference/QTSim/src/PassOptionFile.c`), and [`checkparams`](@ref) checks
the values when a `World` is made. Copy one with changes with `Params(p; field = value)`.
[`read_qtsim_options`](@ref) and [`write_qtsim_options`](@ref) translate QTSim's options
files.

## Fields

| Field | Meaning | Default | QTSim option | QTSim `OPTIONS` member |
|:--|:--|:--|:--|:--|
| `NoGenerations` | generations to run | 20000 | `Generations` | `NoGenerations` |
| `Sample` | record every `Sample`-th generation | 1 | `Sample` | `Sample` |
| `PopSize` | individuals, N, constant from generation to generation | 1000 | `PopSize` | `PopSize` |
| `NoQTL` | loci, each with additive real-valued effects | 1 | `NoQTL` | `NoQTL` |
| `Ploidy` | 1 (haploid, clonal) or 2 (diploid, random mating of two distinct survivors) | 1 | `Ploidy` | `Ploidy` |
| `PhenotypeVar` | V_P, the phenotypic variance the model is calibrated to | 1.0 | `PhenotypeVar` | `PheotypeVar` |
| `Heritability` | h²: V_G = h²·V_P and V_E = (1 − h²)·V_P | 0.5 | `Heritability` | `Heritabitlity` |
| `MutationScalar` | m: mutation adds m·V_E of genetic variance per individual per generation | 0.001 | `MutationScalar` | `MuteVarScalar` |
| `MutationRatePerGamete` | mutations per locus per gamete; each chromosome gets Poisson(`NoQTL`·rate) | 0.0002 | `MutationRatePerGamete` | `MutationRatePerGamete` |
| `RecombPoissonMean` | mean number of crossover points drawn per gamete | 1.0 | `RecombinationPoissonMean` | `RecombinationPoissonMean` |
| `CrossoverModel` | `:qtsim`: only the lowest point switches strands, as in QTSim; `:all`: every point does | `:qtsim` | — | — |
| `InitFitnessSD` | ω, the width of the fitness function; −1 for the analytical ω = V_G/√(m·V_E) | −1 | `InitFitnessSD` | `InitFitnessSD` |
| `ChangeFitnessGen` | generation at whose end ω becomes `ChangeFitnessSD`; −1 for never | −1 | `ChangeFitnessSD g w` (g) | `ChangeFitnessGen` |
| `ChangeFitnessSD` | the new width | 0.0 | `ChangeFitnessSD g w` (w) | `ChangeFitnessSD` |
| `MoveFitnessSDGen` | generation from which the optimum moves; −1 for never | −1 | `MoveFitnessSD g d` (g) | `MoveFitnessSDGen` |
| `MoveFitnessSD` | d, how far it moves (see below) | 1.0 | `MoveFitnessSD g d` (d) | `MoveFitnessSD` |
| `MoveType` | `:phenotype`, `:fitness` or `:continuous`; or `:red_noise` or `:periodic` (Julia-only) | `:phenotype` | `OptimumMoveType` | `MoveType` |
| `MoveNoiseSD` | σ, the stationary standard deviation of the optimum's noise, for `:red_noise` | 0.0 | — | — |
| `MoveAutocorr` | ρ, in [0, 1), the noise's autocorrelation from one generation to the next, for `:red_noise` | 0.0 | — | — |
| `MoveAmplitude` | A, the amplitude of the optimum's cycle, for `:periodic` | 0.0 | — | — |
| `MovePeriod` | P > 0, the cycle's period in generations, for `:periodic` | 100.0 | — | — |
| `Drift` | every fitness 1: a neutral run (QTSim has the option but switches it off) | `false` | (`Drift`, ignored) | `Drift` |
| `StandardisePhenotype` | fitness of (z − θ)/sd(z) instead of z − θ | `false` | `StandardisePhenotype` | `StandardisePhenotype` |
| `Demography` | `:fixed`: the survivors always produce `PopSize` offspring, as in QTSim; `:beverton_holt`: density-regulated (see below) | `:fixed` | — | — |
| `GrowthRate` | R, offspring per survivor at low density, for `:beverton_holt` | 2.0 | — | — |
| `Selection` | `:viability`: survive iff U < W, then mate at random, as in QTSim; `:fecundity`: all survive and mate in proportion to W (see below) | `:viability` | — | — |
| `Plasticity` | b, the slope of a linear reaction norm z = g + b·ε + e on the environment ε (see [Phenotypic plasticity](@ref plasticity)); 0: none | 0.0 | — | — |
| `OptimumSlope` | B ≠ 0, the optimum's slope on the environment, θ = B·ε; only b/B matters for a fixed slope | 1.0 | — | — |
| `CueLag` | τ ≥ 0, generations from the cue that sets the plastic response to selection; 0: a perfect cue | 0 | — | — |
| `PlasticityCostSD` | ω_b, the width of stabilising selection on the slope: fitness × exp(−b²/(2ω_b²)); `Inf`: no cost | `Inf` | — | — |
| `ReactionNorm` | a heritable reaction norm: needs `NoTraits = 2`, trait 1 the elevation and trait 2 the slope (around `Plasticity`) | `false` | — | — |
| `NoTraits` | n, the number of traits each locus affects (pleiotropy; see [Several traits](traits_design.md)) | 1 | — | — |
| `MutCorr` | r_μ, the mutational correlation between every two traits (`NoTraits > 1`) | 0.0 | — | — |
| `SelCorr` | the correlation of the fitness function between every two traits (`NoTraits > 1`) | 0.0 | — | — |
| `MutCov` | M, the mutational covariance per generation, n×n, in place of m·V_E with `MutCorr` | empty | — | — |
| `SelCov` | Ω, the initial fitness covariance, n×n, in place of ω² with `SelCorr` | empty | — | — |
| `MoveDirection` | the direction the optimum moves in, normalised; trait 1 when empty | empty | — | — |
| `FluctDirection` | the direction a fluctuating optimum cycles in (`:periodic`) and, by default, wanders in (`:red_noise`), normalised; `MoveDirection` when empty | empty | — | — |
| `MoveNoiseCov` | Σ, the stationary covariance of a `:red_noise` optimum's noise, n×n, in place of `MoveNoiseSD²·f fᵀ` | empty | — | — |
| `Seed` | seed for [`simulate`](@ref); `nothing` draws one | `nothing` | `Seed` | `Seed` |

## What the parameters set, and what they don't

`PhenotypeVar` and `Heritability` set the environmental variance V_E = (1 − h²)·V_P, the
mutation variance and the analytical width ω. They do not set the realised variances.
The genetic variance starts at zero (every allelic effect is 0) and settles where
mutation and selection balance. That level depends on m, ω, the number of loci and N,
and can be well above or below h²·V_P. So the realised V_P and h² of a run are outputs,
not inputs.

Each generation the optimum θ moves, from `MoveFitnessSDGen` on, by `MoveType`:

- `:phenotype`: θ = z̄ + √V_P·d, a fixed lead over the current mean phenotype;
- `:fitness`: θ = z̄ + ω·d, the same in units of the fitness width;
- `:continuous`: θ = θ + d, a constant rate of environmental change; the paper's
  published runs use this.

Two more moves, which QTSim lacks, make the optimum fluctuate around a linear trend of d
per generation, as field environments do. With t₀ = `MoveFitnessSDGen` and θ₀ the optimum
before the move (0, unless something else has set it), from generation t₀ on:

- `:red_noise`: θ_t = θ₀ + d (t − t₀) + ε_t, where ε is a stationary AR(1) process with
  standard deviation σ = `MoveNoiseSD` and autocorrelation ρ = `MoveAutocorr`:
  ε_{t₀} ~ N(0, σ²), then ε_t = ρ ε_{t−1} + σ √(1 − ρ²) ξ_t with ξ_t ~ N(0, 1). σ sets
  the size of the fluctuations and ρ their time scale, about −1/ln ρ generations, without
  changing their variance. ρ = 0 is white noise. The move takes one normal draw per
  generation from the run's random numbers, so it changes the rest of the run's draws.
- `:periodic`: θ_t = θ₀ + d (t − t₀) + A sin(2π (t − t₀)/P), with A = `MoveAmplitude` and
  P = `MovePeriod`. It draws nothing.

Both are exogenous: θ does not depend on the population. d may be 0. The trend starts at
θ₀ in generation t₀, one step behind `:continuous`, which has already added d by then. In
a metapopulation each deme has its own noise. [`write_qtsim_options`](@ref) refuses these
moves. The SLiM model that `slim_script` generates has both (`moveType`, `moveNoiseSD`,
`moveAutocorr`, `moveAmplitude`, `movePeriod`), and agrees with the kernel
(`validate/slim/compare_fluctuating.jl`). The moments lowering (`run_moments`) follows a
periodic optimum exactly and red noise's trend, or a given sequence of optima. The steady
lag and load are [`periodic_response`](@ref), [`red_noise_lag_var`](@ref),
[`expected_lag2`](@ref) and [`lag_load`](@ref).

Fitness is QTSim's Gaussian `pdf(Normal(0, ω), z − θ) / (0.4/ω)`, which keeps QTSim's
0.4 in place of 1/√(2π), so the fitness at the optimum is 0.997 rather than 1. Selection
is viability selection: an individual survives if a uniform draw falls below its
fitness. The population dies when fewer than two survive (one, if haploid).

## [Phenotypic plasticity](@id plasticity)

With `Plasticity = b ≠ 0` (Julia-only) every individual has a linear reaction norm,
`z = g + b·ε + e`, after Chevin, Lande and Mace (2010). The environment `ε` is measured by
the optimum it sets, `θ = B·ε` with `B = OptimumSlope`, so the plastic part of the
phenotype is `(b/B)·θ`: a fraction `b/B` of the optimum's displacement from the origin of
the trait's scale is tracked without genetic change. The cue is the environment `CueLag =
τ` generations before selection, `z_t = g + (b/B)·θ_{t−τ} + e`; with `τ = 0` development
and selection see the same environment. `PlasticityCostSD = ω_b` puts a cost on
plasticity: every fitness is multiplied by `exp(−b²/(2ω_b²))`.

The norm is expressed in `select`, once the optimum of the generation is set, and the
recorded phenotypes include it. The mean genotype is `pheno_mean − (b/B)·θ_{t−τ}`. It needs
an optimum with its own trajectory: `:continuous`, `:red_noise`, `:periodic`, a step's
fixed `θ`, or none. QTSim's `:phenotype` and `:fitness` moves re-anchor the optimum to
the current mean phenotype every generation, which a plastic shift only carries along,
and [`checkparams`](@ref) rejects them with plasticity.

Under a steadily moving optimum the genotype has to follow only `(1 − b/B)` of the move:
the lag is [`steady_lag`](@ref) times `1 − b/B` and the rate a population can sustain
([`sustainable_rate`](@ref), [`critical_rate`](@ref)) is divided by `|1 − b/B|`, less the
cost. A delayed cue changes nothing under a steady move, and matters when the optimum
fluctuates ([`red_noise_lag_var`](@ref), [`expected_lag2`](@ref)). All three lowerings
have it: the kernel, `run_moments` and `slim_script`
(`validate/slim/compare_plasticity.jl`).

With `ReactionNorm = true` and `NoTraits = 2` the norm is heritable: trait 1 is the
elevation and trait 2 the slope, an individual's slope is `Plasticity + g₂`, its phenotype
`z = g₁ + (Plasticity + g₂)·ε + e`, selection acts on `z` and the cost on the individual's
slope. `MutCov` is the mutational covariance of elevation and slope. The records'
`trait_mean` is then the mean phenotype and the mean slope, and `gmatrix` the genetic
covariance matrix of elevation and slope. The kernel and `run_moments` have it; SLiM does
not. See [Phenotypic plasticity](plasticity_design.md) for the design.

## Population size

QTSim keeps the population at `PopSize`: however few survive selection, they produce
`PopSize` offspring. A population there dies only when fewer than two survive a single
generation, and it pays nothing for a low mean fitness until then.

With `Demography = :beverton_holt` the survivors' offspring are density-regulated. With
`n` survivors, carrying capacity `K = PopSize` and `R = GrowthRate`, the next generation
has `Poisson(R n/(1 + (R − 1) n/K))` individuals, at most `K` ([`noffspring`](@ref)). A
population with mean fitness `W̄` settles near `K (R W̄ − 1)/((R − 1) W̄)` and shrinks when
`R W̄ < 1`, so it can decline and be rescued by adaptation, or die. The SLiM model that
`slim_script` generates draws the number of offspring the same way, and the two backends
agree (`validate/slim/compare_demography.jl`). A population that leaves no offspring dies
in that generation, which, as at selection, has no output row.

QTSim's selection is on viability: an individual survives iff `U < W(z)`, and the
survivors mate at random. With `Selection = :fecundity` everyone survives and each
parent is drawn with probability `W/ΣW` (two distinct parents when diploid), and under
`:beverton_holt` the offspring number is the formula above with `n = ΣW = N W̄`, the
expected number of survivors under viability. The two modes then have the same expected
selection; viability adds the variance of who survives. At a fixed size fecundity
selection is soft: the population dies only when fewer than `Ploidy` individuals have
positive fitness, which in practice means that fitness underflows to zero. The `nsurv`
column is then `N`, and `sel_diff` the differential of the mating weights,
`Σ W z/ΣW − z̄`. The SLiM model draws parents with the same weights
(`validate/slim/compare_selection.jl`).

## Parameters of a deme, and of a step

A metapopulation can give each deme its own `Params`: `[p1, p2, p3] ⊗ demes(3)` in the
kernel, `run_moments([p1, p2, p3]; cycle)` in the moments lowering and
`slim_script(cycle, [p1, p2, p3])` in SLiM. The demes run for the same `NoGenerations`, and
migration needs the same `NoQTL` and `Ploidy`. In the SLiM model a constant whose value
differs between demes holds one value per deme, which each deme keeps on its
subpopulation; `NoQTL`, `NoGenerations`, `Ploidy`, `CrossoverModel` and `Demography` must
be the same in every deme there.

A step of the life cycle can also carry parameters, which override the deme's for that
step only ([`withparams`](@ref)): `withparams(L.select; ω = 0.3)` selects with width 0.3
whatever the deme's width is, `withparams(g; θ = 0.1)` is a generation whose selection is
toward a fixed optimum 0.1, and `lifecycle().migrate(k; rate)` migrates at its own rate.
The fields a step can override are those its image reads: the move, `ω`, `θ`,
`StandardisePhenotype`, `Drift`, `Selection` and the reaction norm (`Plasticity`, `OptimumSlope`,
`CueLag`, `PlasticityCostSD`) for `select`; `RecombPoissonMean`, `CrossoverModel`,
`Demography` and `GrowthRate` for `reproduce`; `Sample` for `record`; `ChangeFitnessGen`
and `ChangeFitnessSD` for `rescale`; `MutationScalar` and `MutationRatePerGamete` for
`mutate`. A step with the deme's own values is the plain step, bit for bit, and a step's
fields of `Params` give the same run as `Params` with those values. A step's `ω` and `θ`
leave the deme's state alone, so the `omega` and `theta` columns record the deme's.
