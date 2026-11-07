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
