# Symbolic checks of the moments closures

`check_closures.jl` checks with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)
the algebra behind the closures of the moments lowering, and evaluates each identity at
random parameter values against the package's own functions, so that the code and the
maths stay in step. Lean (`../../../lean/CategoricalQuantGen/`) proves the theorems;
this script checks the algebra the code uses.

```sh
julia --project=validate/symbolic -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=validate/symbolic validate/symbolic/check_closures.jl  # about 30 s
```

It writes `check_closures.log` and exits 1 if any check fails. It covers:

- the closure with locally favoured alleles (`loci = :rare_local`, src/local_alleles.jl):
  a mutant's marginal fitness in a lagging deme, the second-order growth rate of a rare
  lineage under selection and migration (the characteristic polynomial of `P(I + εΣ)`
  vanishes to `O(ε³)`), the closed form of `W = Π(I − P)^#P` for two demes, the gain's
  invariance under a common shift of the lags, the island form `(1 − m′)/m′·Var_π(ℓ)`,
  and the reductions to the current closure (one deme, lags alike, no migration,
  complete migration);
- the existing closures: the rare-alleles per-gene fixed point and Hermisson and Wagner's
  equation A4, the multivariate closure's exactness at the house of cards and its
  one-trait reduction, the two-trait house-of-cards weights, and the multilocus Gaussian
  fixed point `V²/(V_s + V) = 2L·V_m`;
- numerical spot checks of `house_of_gauss_vg`, `gaussian_vg_loci`, `rare_alleles_rate`,
  `rare_hc_matrix`, `rare_hc_rates`, `local_weights` (also against finite differences of
  the leading eigenvalue for k = 2–6), `local_gain`, `local_sinv` and `hap_select`.

Square roots are handled by reducing a polynomial in `r = √D` modulo `r² = D`.
Symbolics is a dependency of this environment only, not of QTLSim. CI runs the script in
the `symbolic` job.
