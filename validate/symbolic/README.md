# Symbolic checks of the theory

The scripts here check with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)
the algebra behind the package's theory and the closures of the moments lowering, and
evaluate each identity at random parameter values against the package's own functions,
so that the code and the maths stay in step. Lean (`../../../lean/CategoricalQuantGen/`)
proves the theorems; these scripts check the algebra the code uses. Where a documented
identity holds only to some order, the check says so and checks it to that order.

```sh
julia --project=validate/symbolic -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=validate/symbolic validate/symbolic/check_all.jl      # all of them, about 1 min
julia --project=validate/symbolic validate/symbolic/check_spectrum.jl # or one topic
```

Each script writes its `.log` and exits 1 if any check fails; `check_all.jl` runs them in
one process (Symbolics loads once) and exits 1 if any fails. CI runs it in the `symbolic`
job. `common.jl` holds the shared helpers.

- `check_closures.jl` (73 checks) — the closure with locally favoured alleles
  (`loci = :rare_local`, src/local_alleles.jl): a mutant's marginal fitness in a lagging
  deme, the second-order growth rate of a rare lineage under selection and migration, the
  closed form of `W = Π(I − P)^#P` for two demes, the gain's invariance under a common
  shift of the lags, the island form, and the reductions to the current closure; the
  rare-alleles per-gene fixed point and Hermisson and Wagner's equation A4, the
  multivariate closure's exactness at the house of cards, the two-trait house-of-cards
  weights and the multilocus Gaussian fixed point; spot checks of `house_of_gauss_vg`,
  `gaussian_vg_loci`, `rare_alleles_rate`, `rare_hc_matrix`, `rare_hc_rates`,
  `local_weights`, `local_gain`, `local_sinv`, `hap_select`.
- `check_spectrum.jl` (45) — the spectrum model (`loci = :spectrum`, src/spectrum.jl):
  `E[Δp]` is Wright's formula on `ln W̄ = −λ²/(2S) − ln S/2` exactly, and the exact change
  under viability selection on a diploid locus in a Gaussian background to second order in
  the effect, at any lag (Hayward and Sella eq. 7); Kopp and Hermisson's `s(α)` up to
  `a²λ²/(2S²)`; the `2c` term gives the rare-alleles closure's `c²/S` (and the `a²` term
  `α²c/(2S)` at the house of cards); the neutral random field gives `(2N − 1)V_m ≈ 2NV_m`
  (Wright–Fisher and diffusion), `(2N − 1/(2N))V_m` with the implicit Euler step; the
  symmetric-spectrum identity `Σ2aE[Δp] = V_Gλ/S`; the sweep identity
  (`2aλ/S · 2aS/λ = 4a²`, so `2N_eV_m`); Barton's linked-selection factor; spot checks of
  `spectrum_seff`, `spectrum_select`, `spectrum_vg`, `spectrum_linked` (also against the sum
  over the pairs themselves), `spectrum_mutate!`, the Gauss–Hermite effects and
  `spectrum_equilibrium` (neutral: to 2e-5; selected: close to `house_of_gauss_vg`, an
  approximation).
- `check_traits.jl` (52) — several traits (src/traits_theory.jl,
  src/rare_alleles_traits.jl): `gaussian_G` (the substitution `G = A^½XA^½`, the root on
  each eigenvalue and on shared axes; numerically for general `M`, `A`, against
  `gaussian_G_iterate`, and `gaussian_G_loci = gaussian_G(2L·M, A)`); `steady_lag_vector`;
  `red_noise_lag_cov` (the lag–optimum recursion, the one-trait Lyapunov solution, which is
  `red_noise_lag_var`, and `X = FXFᵀ + Q` for two traits); `periodic_response` and
  `periodic_lag_cov` (the transfer function, the lag amplitude, simulations of the
  recursion); the deflection inequality (symbolic 2 × 2, 300 random examples); the
  rare-alleles Riccati equation at `house_of_gauss_G`, its stabilising root, its one-trait
  reduction to `house_of_gauss_vg`, and its limits (`house_of_cards_G` for rare mutations
  when `α² ≪ A`; the weak-mutation Gaussian `G A⁻¹ G = 2L·M` for many small ones).
- `check_fst.jl` (49) — neutral differentiation (src/theory.jl, src/migration.jl): the
  fixed point of the island model's three-variable coalescence recursion, `V_W = 2NkV_m`,
  `σ²_B = 2kV_mλ/(1 − λ)`, `Q_ST = λ/(λ + 2N(1 − λ))` (both conventions), Wright's and
  Takahata's forms; under the kernel's exchanges `V_W = 2NV_m(kN − 1)/(N − 1)` to leading
  order in `m` (relative error `O(m/N)`); Slatkin's stepping stone; spot checks of
  `qst_island`, `island_variances`, `deme_coalescence`, `neutral_deme_variances`,
  `exchange_origins` and `backward_migration`.
- `check_demography.jl` (~15) — the Beverton–Holt fixed points `K(RW̄ − 1)/((R − 1)W̄)`
  (adults) and `K(RW̄ − 1)/(R − 1)` (survivors) and their stability, against
  `expected_offspring`; the thinned Poisson pgf of the Galton–Watson process, against
  `rescue_probability` and `extinction_distribution(:branching)`.
- `check_plasticity.jl` (58) — phenotypic plasticity (src/plasticity.jl and the plasticity
  keywords of src/theory.jl): the steady lag `(1 − b/B) k/β` for every delay of the cue and
  the exact reductions at `b = 0`; the sustainable and critical rates (`W̄ = W_min` at the
  steady lag, the factor `1/|1 − b/B|`, the cost in `ln B₀`); the lag behind red noise as
  `L′ = aL + θ′ − θ`, its autocovariance and the variance of `L_t − φL_{t−τ}`, and the cycle's
  factor `1 + φ² − 2φ cos ντ`; a heritable reaction norm: `G_s = (I + GQ)⁻¹G = G − GKG`
  (exactly, in rational arithmetic, with a cost), the Sherman–Morrison form without one, the
  equilibria in a constant environment, the quasi-steady lag and `ε·Δb̄ → (1 − b̄/B)k`, the
  mean fitness with the cost by quadrature; spot checks of `steady_lag`, `critical_rate`,
  `sustainable_rate`, `red_noise_lag_var`, `red_noise_lag_autocov`, `expected_lag2`,
  `lag_load`, `reaction_norm_response` and `reaction_norm_lag`.

Square roots are handled by reducing a polynomial in `r = √D` modulo `r² = D`. Large
rational identities are checked by substituting a closed form into the equations it
solves, not by solving them symbolically (Symbolics' simplifier is slow on those).
Symbolics is a dependency of this environment only, not of QTLSim.
