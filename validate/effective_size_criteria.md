# Is `effective_size(N, n)` the effective size of QTSim's life cycle? Criteria

Written and committed before any kernel run of `validate/effective_size.jl`. (Seen before
writing: the side observation in `validate/xu_shen.md`, a plain simulation outside the
kernel that gave N_e ≈ n; and the algebra below. No kernel measurement.)

## The question

`effective_size(N, n) = 4N/(V_k + 2)`, `V_k = 2N(n − 1)/n²`, is close to the harmonic mean
`2Nn/(N + n)` of the adults N and the breeders n. It is used by `house_of_gauss_vg`,
`house_of_gauss_G`, `house_of_cards_vg` (keyword `nsurv`, default `PopSize`) and by
`drift_critical_rate` (with the breeders the lag leaves). The moments recursions
(`hap_reproduce`, the spectrum's δ) do not call it: they use the number of parents n
directly.

## Derivation (first principles; algebra in `validate/symbolic/check_effective_size.jl`)

Life cycle, neutral: N adults; n of them survive, each with the same probability; N
offspring, each with two DISTINCT parents drawn uniformly from the n survivors, one gamete
from each (a random one of the parent's two copies at a locus); mutation; the offspring are
the next adults. Let x be a moment of two gene copies in two different adults and y of the
two copies of one adult. Looking back one generation, two copies in different offspring
come from the same parent with probability 1/n (the parents of different offspring are
independent draws) and then from the same copy with probability ½; the two copies of one
offspring come from two different parents. With ε = E[1/n] (n is random in the kernel,
independent of the genes when survival is neutral, so expectations use E[1/n]):

    probability of identity by descent / expected squared difference d = ½E(a − a′)²:
      d_x′ = (1 − ε) d_x + (ε/2) d_y + u,     d_y′ = d_x + u,      u = V_m/2 (mutation)

1. **Inbreeding (coalescence) size.** Two copies in different adults coalesce in the
   previous generation with probability ε/2 = 1/(2n_h), n_h = 1/E[1/n] the harmonic mean of
   the breeders. Mean coalescence times: T_x = 2n_h + 1, T_y = 2n_h + 2.
2. **Mutation–drift balance.** d_x = (V_m/2)(2n_h + 1), d_y = d_x + V_m/2, and the expected
   sample variance of the genotype among adults is E[V_G] = 2d_x + 2(d_x − d_y) = 2 n_h V_m
   exactly (no linkage disequilibrium in expectation: mutations are symmetric and
   independent). So the N_e of Lynch and Hill's 2N_eV_m is n_h.
3. **Eigenvalue size.** Without mutation both moments decay by the largest root λ of
   λ² − (1 − ε)λ − ε/2 = 0 per generation, 1 − λ = 1/(2n_h + 1) + O(n_h⁻³):
   N_e = 1/(2(1 − λ)) ≈ n_h + ½.
4. **Variance size.** Var(Δp) from adults to adults: hypergeometric survival of n of N
   diploids, then 2N gametes drawn from the survivors: pq[(1/(2n))(1 − n/N) + 1/(2N)] =
   pq/(2n) to leading order: N_e = n.
5. **The standard formula.** Counting all N adults as potential parents (the dead have 0
   offspring) their gamete number has mean 2 and variance V_k = 2(1 − 1/n) + 4N/n − 4, and
   4N/(V_k + 2) = n·4N/(4N − 2) ≈ n. Counting only the n breeders, the mean is 2N/n, not 2,
   and the general form (N_p k̄ − 1)/(k̄ − 1 + V_k/k̄) with N_p = n, k̄ = 2N/n, V_k =
   2N(n − 1)/n² gives exactly n. `effective_size` puts the breeders' V_k into the formula
   for a mean of 2: that is what gives ≈ 2Nn/(N + n). The two agree at n = N.
6. Haploids (one parent per offspring, clones): coalescence 1/n, N_e = n. Fecundity
   selection (all N adults breed, with probability ∝ W): N_e = (ΣW)²/ΣW² = N·W̄²/E[W²], which
   the moments lowering already uses; `effective_size` is not used for it.

All four definitions give n (to O(1)); none gives 2Nn/(N + n). Which one each use needs:
`house_of_gauss_vg`/`house_of_cards_vg`/`house_of_gauss_G` and `drift_critical_rate` use
`1/(Ploidy·N_e)` as the per-generation loss rate of a gene's variance in a mutation–drift
(–selection) balance: definitions 2 and 3.

Prediction to be tested: **derived** N_e = n_h (mutation–drift balance and coalescence) and
n_h + ½ (decay); **current** = `effective_size(N, n̄)`, n̄ = N·W; **n̄** itself (differs
from n_h only by the binomial spread of n, ≤ 5% in the smallest cell).

## Test in the kernel (no change to `src/`)

Neutral survival of a fraction of adults is obtained with existing `Params`:
`InitFitnessSD = 1e6` (the Gaussian fitness is flat), a cost of plasticity on a reaction
norm that expresses nothing (`Plasticity = 1`, `OptimumSlope = 1e12`, `PlasticityCostSD =
1/√(2c)`), which multiplies everyone's fitness by e^{−c}: every adult survives with the
same probability W = 0.99736·e^{−c} (0.99736 = QTSim's fitness at the optimum), so n ~
Binomial(N, W). The cell n/N = 1 uses `Drift = true` (everyone survives). The runs are the
kernel's own `generation!`. n_h is computed exactly from the binomial (conditional on
n ≥ 2).

Grid: N ∈ {200, 1000} × W ∈ {1, 0.75, 0.5, 0.25, 0.1} (10 cells). Diploid, 50 loci,
`CrossoverModel = :all`, `RecombinationPoissonMean = 10`.

Three measurements per cell, each turned into an estimate of N_e with a standard error
from the replicates (delta method on the mean over replicates):

- **(i) decay**: no mutation, alleles started as independent N(0,1); H(t) = sum over loci
  of the variance among the 2N gene copies; N̂_e = 1/(2(1 − (H̄(t₁)/H̄(t₀))^{1/T})), t₀ = 10,
  T = t₁ − t₀ = round(n̄), H̄ the mean over replicates; 300 replicates.
- **(ii) mutation–drift balance**: m = 0.001, h² = 0.5, V_P = 1 (V_m = 5·10⁻⁴), 0.01
  mutations per locus per gamete; from zero, 12n̄ generations discarded, then the mean of the
  recorded `geno_var` over the next 12n̄ generations (about 240 samples per run);
  N̂_e = mean/(2V_m); 120 replicates.
- **(iii) coalescence**: `record_pedigree!` (simplified every 100 generations), 10 loci,
  16n̄ generations from tagged founders, `pair_coalescence`; N̂_e = (mean time of pairs in
  different individuals − 1)/2 (pairs not coalesced count 16n̄: bias < 0.1%); 80 replicates.

If a cell's standard error exceeds 3% of the estimate, its replicates are doubled once
(new seeds added, none dropped).

## Verdict rule (fixed)

For each formula F ∈ {current, n̄, derived} and each measurement, a cell **agrees** if
|N̂_e − F| ≤ max(0.05·F, 2 SE). A formula **fits** if it agrees in all 30
cell × measurement combinations; it is **rejected** if it disagrees in at least one cell by
more than both 5% and 3 SE. Otherwise inconclusive. The current formula is expected to
agree at W = 1 under every hypothesis; the informative cells are W ≤ 0.5, where it is 33%,
60% and 82% above n̄.

Seeds: 77_000_000 + 1000·cell + replicate, per measurement offset 0, 10⁵, 2·10⁵.

## What follows from the verdict

If the current formula is rejected, the derived one is added as an option
(`effective_size(N, n; scheme)`, and a keyword on the functions that call it), off by
default: every default stays bitwise the same, and the choice of default is left to the
user. The impact with the option on and off is then reported for `drift_critical_rate` on
published-type inputs (fixed size and Beverton–Holt R = 5, 2, 1.5), for the pre-move V_G
against the published table, for `validate/critical_rate.jl`'s headline numbers recomputed
from its caches, and for the emulator's half-survival rate on about 30 published inputs.
If the current formula fits, nothing is changed and the earlier observation is explained.

## Amendment 1 (before any run on the grid)

A smoke test of the three run functions off the grid (N = 40, W = 0.5 and 1, one run each;
and one timing run each of (ii) and (iii) at N = 1000, W = 0.1, results not looked at)
showed 0.85–1 ms per generation at N = 1000 with these crossovers, which would put the
plan at about 3 CPU-hours, the whole budget of stage 1. Replicates are therefore
(i) 300, (ii) 60, (iii) 50; the doubling rule for cells with a standard error above 3% is
unchanged and decides whether that is enough. Nothing else changes.
