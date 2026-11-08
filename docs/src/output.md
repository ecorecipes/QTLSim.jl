# Output

## `simulate` and `run!`

[`simulate`](@ref) and [`run!`](@ref) return a `NamedTuple` of columns, one row per
recorded generation, every `Sample` generations from generation 0:

| Column | Meaning |
|:--|:--|
| `gen` | the generation, counted from 0 |
| `pheno_mean` | mean phenotype z̄ |
| `pheno_var` | phenotypic variance |
| `geno_var` | additive genetic variance V_G |
| `nsurv` | individuals that survived viability selection |
| `theta` | the optimum θ that generation |
| `omega` | the fitness width ω that generation |
| `mean_fitness` | mean fitness W̄ of the parents (one minus the load); not in QTSim's output |
| `sel_diff` | selection differential: the survivors' mean phenotype minus the mean; not in QTSim's output |
| `popsize` | the population size N (`PopSize`, unless density-regulated); not in QTSim's output; in SLiM's table, the last column when density-regulated |

These columns are joined by:
- `extinct`, whether the population died;
- `extinct_gen`, the generation in which it died, or −1;
- `seed`, the seed used (from `simulate` only).

The details, which follow QTSim's `Output`:

- A row describes the **parents** of the generation: their phenotypes and genotypes, the
  optimum they were selected against, and how many survived. The offspring appear in
  the next row.
- Variances are sample variances (divided by n − 1), as GSL's `gsl_stats_variance`.
- The generation in which the population dies has no row: `extinct_gen` is the
  generation after the last row (with `Sample = 1`).
- QTSim also prints a phenotype ratio, slopes of survival on phenotype, and
  w = V_G/√MuteVar; QTLSim does not compute these.

## `run_composed` on a metapopulation

`run_composed(p ⊗ demes(k), rng; cycle)` returns the same columns with a row per deme
and recorded generation, plus a `deme` column (1 to k). For metapopulations:
- `extinct_gen` is a vector, the generation in which each deme died, or −1;
- `extinct` is true when every deme has died;
- a deme that dies stops producing rows, and the others run on.

## `run_moments`

[`run_moments`](@ref) returns the same columns, predicted by the moments lowering: the
means and variances of a Gaussian population, `nsurv` the expected survivors rounded. By
default (the rare-alleles closure for one diploid trait, the multilocus Gaussian model for
several traits and haploids; `loci`) `geno_var` is the infinite population's genetic variance,
or with `erosion = true` (or `between = true`) the expected within-population sample
variance, as the kernel records it.
[`extinction_distribution`](@ref) takes a `run_moments` result with `Sample = 1` and
returns `gen`, `pext` (the probability of extinction by the end of each generation, from
demographic stochasticity along the trajectory), `median` and `psurvive`.

## The SLiM models

A SLiM model written by [`slim_script`](@ref) prints a tab-separated table to standard
output. Its header is `gen  phenoMean  genoVar  nSurv  theta  omega`, with a last column
`deme` on several demes. The rows correspond to `run!`'s. A population that dies prints
`Pop all died.`, as QTSim does.

## Reproducibility

`simulate(p; seed)` gives the same run for the same seed, with the same versions of
Julia and Distributions.jl. The run is not QTSim's run with that seed: QTLSim's random
numbers come from Julia's `MersenneTwister`, not GSL's. Runs agree with QTSim in
distribution, and `validate/gsl_locked.jl` shows the algorithm is QTSim's by replaying
it with GSL's random numbers, bit for bit.
