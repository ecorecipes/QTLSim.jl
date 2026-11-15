# Criteria: the beta closure of the `:spectrum` model (`spectrum = (; solver = :beta)`)

Written and committed before `validate/spectrum_beta.jl` is run. This is a check of
**numerical equivalence between two solvers of one model**, the beta closure against the
frequency grid of `loci = :spectrum` (`solver = :grid`, the default). It is not a comparison
with the kernel or with the published table, and passing it says nothing about how well
either solver predicts the kernel.

## The option under test

Per effect class, the pq-weighted spectrum ψ = pqφ is carried by its moments
P_k = ∫p^kψ, k = 0, 1, 2, closed by a beta law ψ/P₀ = Beta(α, β) (exact on the neutral
spectrum ψ ∝ 1 − p), advanced by one explicit Euler step per generation of the exact
ψ-moment equations, with Barton's linked-selection term and the same mutation input, and
projected onto the Hausdorff moment conditions (src/spectrum.jl). After each generation,
classes carrying at least 0.1% of V_G (and not an atom of new mutations) are checked; the
run falls back to the grid, for that generation and the rest of the run, if any such class
has (i) population-scaled directional selection Γ = 4n|r₁aλ|/S > 100, (ii) a U-shaped closed
law (α < 1 and β < 1), or (iii) a Hausdorff projection moving (μ₁, μ₂) by more than 1e-6 of
μ₁(1 − μ₁). The last state inside the domain is handed to the grid with V_G kept exactly.

**How the rule was chosen (disclosed).** The shape criteria (ii) and (iii), suggested in
advance, were tried first on four published inputs from an earlier scratch study (04943,
14110, 12000, 06159). They fired too late (04943: V_G 38% above the grid 100 generations
into the move, U-shape at generation 126) or never (06159: 18% below at +600). A tracked
third moment compared with the beta law's did not discriminate either. Criterion (i) and
its threshold of 100 were set on those four inputs. Two fallback policies were then
compared on a pilot sample of 40 published inputs (5 per stratum, `MersenneTwister(777)`;
none is reused below unless the new seed draws it): handing over in place, or rerunning the
whole run on the grid. Handing over was kept: rerunning improved only the medians and could
cost more than the grid alone. Pilot results (no verdict): fallback rate 25%, P2000 calls
40/40, median |V_G ratio − 1| 0.0005 / 0.0024 / 0.0053 (equilibrium / +100 / +500),
95th percentiles 0.0095 / 0.098 / 0.049, median speed-up about 300× per input. The
thresholds below are the ones suggested before the pilot; none was changed after it.

## Sample

A fresh uniform sample of published inputs (`Repository/Simulations`), **50 per stratum**
(8 strata: NoQTL ∈ {100, 500, 1000, 5000} × RecombinationPoissonMean ∈ {3, 5}), **400 inputs**,
drawn with `MersenneTwister(20261130)` (a seed not used before) by a random permutation of
each stratum's files in sorted order. Earlier studies have used published inputs, and some of
these may overlap with them. That is acceptable here, since nothing is compared with
simulations or with the table: the question is only whether two numerical methods agree.

- **ld subsample:** the first 25 per stratum of the same permutations (200 inputs) for the
  same comparison with `ld = :dynamic`.
- **Half-survival subsample:** the first 6 per stratum (48 inputs) for the half-survival rate.

## Protocol

Every input is run as the emulator runs it (`emulate` defaults: `burnin = 2000`,
`horizon = 2000`, `VG = house_of_gauss_vg(p)`, `erosion = true`, fixed size as published),
once with `spectrum = (;)` and once with `spectrum = (; solver = :beta)` (each with
`ld = :dynamic` added for the ld subsample). Both solvers of one input run back to back in
one job, on 8 threads.

## Metrics and thresholds (beta against grid)

Ratios are beta/grid. A V_G ratio at +t is taken over the inputs whose trajectories are both
alive at generation burnin + t.

1. **Equilibrium V_G** (`spectrum_vg(p)`): median |ratio − 1| ≤ 0.02 **and** 95th percentile ≤ 0.10.
2. **V_G at +100 and at +500** generations of the move (`run_moments` as `emulate` runs it):
   for each of the two, median |ratio − 1| ≤ 0.02 **and** 95th percentile ≤ 0.10.
3. **The emulator** (`emulate`):
   - agreement of the calls `P2000 ≥ 1/2` ≥ 98% of inputs;
   - T (median extinction time after the move): median ratio within [0.98, 1.02], over inputs
     where both give a finite T. Inputs where exactly one gives a finite T count as disagreements
     and are reported; no threshold.
4. **Half-survival rate** (`emulated_half_survival_rate`, defaults) on the 48-input subsample:
   median ratio of k50 within [0.98, 1.02] **and** 95th percentile of |ratio − 1| ≤ 0.10.
   Inputs where either k50 is 0 or Inf count as agreeing only if both are the same, and are
   reported.
5. **With `ld = :dynamic`** (200-input subsample): metrics 1–3 with the same thresholds.
6. **Fallback rate**: the fraction of inputs whose `emulate` run fell back to the grid at
   any point (the field `fallback ≥ 0`) ≤ 0.50. Also reported: the fraction falling back
   before the move, and the distribution of the fallback generation.
7. **Speed-up**: per input, the wall time of `emulate` with the grid divided by that with the
   beta closure, measured in the same job (fallbacks included): median ≥ 10×. Also reported:
   quartiles, minimum, and the median among inputs that fell back.

Each numbered item gets PASS or FAIL. All verdicts are reported, including failures. The
default (`solver = :grid`) is not changed by this check; changing defaults is the user's
decision.

## Amendments

(none)
