# `effective_size` and `drift_critical_rate`: two checks of existing code

2026-11-19. Files: `validate/effective_size_criteria.md` (derivation, design and verdict
rule, committed before any kernel run; two amendments, both before any result),
`validate/effective_size.jl` + `effective_size_parts.jl`, logs `effective_size.log` (part
1a, kernel), `effective_size_branch.log` (part 1b), `effective_size_impact.log`,
`validate/symbolic/check_effective_size.jl` (29 checks), cache `effective_size_runs.csv`.
Nothing in the kernel changed, and no default: both findings are implemented as options.

**Decision, 2026-11-19 (after this study; the text below is left as written).** Both
recommendations were adopted: `effective_size(N, n)` now returns the breeders `n` by
default (`scheme = :breeders`; `ne = :breeders` in `house_of_gauss_vg`, `house_of_gauss_G`,
`house_of_cards_vg`, `drift_critical_rate`, `half_survival_rate`), and
`drift_critical_rate` uses `solver = :bisection` by default. The earlier behaviours are
the options `:falconer` and `:iteration`. Where the text below says "the default" or "the
current formula" it means the Falconer-type form and the iteration, the defaults when it
was written; this study's scripts now name them explicitly, so they still reproduce its
logs. `validate/critical_rate.jl` was rerun on fresh seeds under the new defaults (its
addendum and `critical_rate.log`; run 1 is `critical_rate_falconer.log`), and the Lean
statements are in `CategoricalQuantGen/EffectiveSize.lean`. On the harmonic mean: the
functions are called at a deterministic steady state with the expected breeders n̄, so
`:breeders` uses n̄; for binomial survival n_h ≈ n̄ − (1 − n̄/N), less than one individual
lower, and the correction is not applied (see `effective_size`'s docstring).

## 1a. `effective_size(N, n)` is not the effective size of QTSim's life cycle

**What it computes.** `4N/(V_k + 2)` with `V_k = 2N(n − 1)/n²`, i.e. `2Nn²/(N(n − 1) + n²)`,
within `1/(2n²)` (in the reciprocal) of the harmonic mean `2Nn/(N + n)` of the adults and
the breeders.

**What the life cycle gives** (N adults; n survive, neutrally; N offspring, each with two
distinct parents drawn uniformly from the survivors). With x a moment of two gene copies
in different adults and y of the two copies of one adult,

    x′ = (1 − 1/n) x + y/(2n) + u,    y′ = x + u,

because two copies in different offspring share a parent with probability 1/n and then a
copy with probability ½, and the two copies of one offspring come from different parents.
From this:

| definition | value | used by |
|---|---|---|
| inbreeding / coalescence: probability 1/(2n) per generation, mean pair time 2n + 1 | n | ancestry (`deme_coalescence` has it already) |
| mutation–drift balance: E[V_G] = 2nV_m exactly | n | `house_of_gauss_vg`/`_G`, `house_of_cards_vg`, `drift_critical_rate` (the 1/(Ploidy·N_e) loss of a gene's variance) |
| eigenvalue (decay of variance without mutation): 1 − λ = 1/(2n + 1) | n + ½ | the same |
| variance: Var(Δp) = pq[(1 − n/N)/(2n) + 1/(2N)] | n | drift of the mean (the moments lowering has V_G/n already) |
| offspring-number formula over all N adults (mean 2, V_k = 2(1 − 1/n) + 4N/n − 4) | 4Nn/(4N − 2) ≈ n | — |
| the general form over the n breeders (mean 2N/n, V_k = 2N(n − 1)/n²) | n exactly | — |

N does not enter. `effective_size` puts the breeders' V_k, whose mean is 2N/n, into the
formula that holds for a mean of 2; the two agree only at n = N. It is close to the
recipe of the paper's table (N₀ = N(2n − 3)/(2n(n − 2)), V_k = 2N₀(1 − 2N₀/N), N_e =
4N/(V_k + 2), which is ≈ 2Nn/(N + n) as well), from which it was taken; whether that
recipe is right for Bürger and Lynch's own life cycle (many offspring, then culling) is not
examined here, and their equations were not consulted. Haploid clones: n. Fecundity
selection: (ΣW)²/ΣW² = N W̄²/E[W²], which the moments lowering already uses;
`drift_critical_rate` has no fecundity form. In the kernel n varies between generations, so
n is the harmonic mean n_h of the breeders.

**The kernel** (neutral survival of a fraction W through a cost of plasticity on a
reaction norm that expresses nothing; `Drift = true` for W = 1; 50 loci, 3940 runs,
about 1 CPU-hour). N̂_e with its standard error (1–2.3%; no cell needed the doubling
rule), and its ratio to each formula:

| N | W | (i) decay | (ii) E[V_G]/(2V_m) | (iii) (pair time − 1)/2 | `effective_size` | n̄ = NW | n_h |
|---|---|---|---|---|---|---|---|
| 200 | 1 | 195.9 ± 2.9 | 196.6 ± 2.3 | 193.7 ± 4.4 | 200.5 | 200 | 200 |
| 200 | 0.75 | 152.6 ± 2.2 | 151.8 ± 1.8 | 149.9 ± 2.9 | 172.1 | 150 | 149.8 |
| 200 | 0.5 | 97.2 ± 1.5 | 98.5 ± 1.1 | 101.3 ± 2.2 | 134.2 | 100 | 99.5 |
| 200 | 0.25 | 49.5 ± 0.8 | 48.4 ± 0.7 | 51.6 ± 1.0 | 81.3 | 50 | 49.2 |
| 200 | 0.1 | 19.3 ± 0.4 | 19.2 ± 0.3 | 18.9 ± 0.4 | 38.1 | 20 | 19.0 |
| 1000 | 1 | 1002.7 ± 15.7 | — | (1015.9 ± 21.6, 40 runs, outside the verdict) | 1000.5 | 1000 | 1000 |
| 1000 | 0.75 | 725.6 ± 10.2 | 735.9 ± 7.9 | — | 857.8 | 750 | 749.8 |
| 1000 | 0.5 | 494.0 ± 7.7 | 492.4 ± 5.8 | 485.0 ± 9.0 | 667.6 | 500 | 499.5 |
| 1000 | 0.25 | 250.4 ± 3.9 | 248.4 ± 2.5 | 250.4 ± 5.5 | 401.3 | 250 | 249.3 |
| 1000 | 0.1 | 102.6 ± 1.8 | 100.7 ± 1.3 | 101.0 ± 2.1 | 183.5 | 100 | 99.1 |

Verdicts by the pre-set rule (agree: within 5% or 2 SE; rejected: off by more than 5% and
3 SE in some cell), over the 27 cell × measurement combinations that were run:

- **current `effective_size(N, n̄)`: REJECTED** — agrees in 4 of 27 (the four W = 1 cells),
  rejected in 23. Measured/formula is 0.85–0.89 at W = 0.75, 0.72–0.76 at 0.5, 0.60–0.64
  at 0.25 and 0.50–0.56 at 0.1.
- **derived n_h (n_h + ½ for the decay): FITS** — 27 of 27, ratios 0.967–1.048.
- n̄ = N·W: agrees in 26 of 27, rejected in one (coalescence, N = 200, W = 0.1: 18.9 ± 0.4
  against 20; the harmonic mean is 19.05), so by the rule also rejected; it is n_h without
  the 0–5% correction for the binomial spread of n.

The side observation in `validate/xu_shen.md` (a plain simulation, 198 / 100 / 51 against
200 / 134 / 81) is confirmed in the kernel three ways: 196 / 97–101 / 48–52.

Amendments (both before any estimate was computed): replicates of (ii) and (iii) reduced
to 60 and 50 after timing; the N = 1000 cells with W = 1 (balance, coalescence) and
W = 0.75 (coalescence) dropped for cost. One row of a dropped cell was displayed by
accident while reading the cache's header; it is in the criteria file.

**What was implemented** (opt-in; every default bitwise unchanged):
`effective_size(N, n; scheme = :falconer | :breeders)` (`:breeders` returns n), and the
keyword `ne` on `house_of_gauss_vg`, `house_of_gauss_G`, `house_of_cards_vg`,
`drift_critical_rate` and `half_survival_rate`. Nothing else needed a switch:
**the moments lowering does not call `effective_size`**. `hap_reproduce` uses the parents n
directly ((T + R)/(2n) + (1 − 1/n)Q is the recursion above), and its neutral fixed point is
2nV_m to 6 digits (symbolic check (h), and a unit test). So `run_moments`,
`multilocus_equilibrium`, `gaussian_vg_loci(…; drift)`, `rare_alleles_vg`, the spectrum
closure, `emulate` and the survival part of `half_survival_rate` were right already and do
not move. What moves is `drift_critical_rate`, and the closed forms when a caller passes
`nsurv < PopSize` (inside the package only one declared reduction in
`ext/QTLSimCatlabExt/morphisms.jl` does, at W̄ ≈ 0.99).

**Impact of `ne = :breeders`** (`effective_size_impact.log`; 600 published inputs, 75 per
stratum):

| quantity | `:falconer` (default) | `:breeders` |
|---|---|---|
| `drift_critical_rate`, k_c(:breeders)/k_c(:falconer), min / median / max: fixed size | 1 | 0.72 / 0.77 / 0.89 |
| … Beverton–Holt R = 5 | 1 | 0.80 / 0.83 / 0.95 |
| … R = 2 | 1 | 0.89 / 0.90 / 0.97 |
| … R = 1.5 | 1 | 0.92 / 0.94 / 0.98 |
| pre-move V_G, published table / house of Gauss with N_e from the expected breeders (median, all strata) | 0.912 | 0.915 |
| … / `multilocus_equilibrium(drift = true)` (does not use the function) | 0.985 | 0.985 |
| `critical_rate.jl` C1, cost ratio k50(c)/k50(0) within 0.06 of k_c(c)/k_c(0) | 24 of 24 (max error 0.054) | **12 of 24** (max 0.140; observed ≥ predicted in 24 of 24) |
| `critical_rate.jl` C2, M_obs/M_theory within 0.8–1.25 | 8 of 8 (1.09–1.22) | **1 of 8** (1.18–1.77) |
| `critical_rate.jl` k50/k_c, 100 loci / 500 loci | 0.77–1.03 / 0.44–0.75 | 0.97–1.38 / 0.57–0.87 |
| `half_survival_rate` (C3, 26 of 32) | — | unchanged: two cells recomputed, ln ratio 1·10⁻⁴ and 1·10⁻⁵ (the bisection's resolution) |
| emulator k50 (`:spectrum`, beta solver), 32 published inputs | — | unchanged: max \|ln ratio\| 5·10⁻⁹ |
| emulator k50/k_c, median: 100 / 500 / 1000 / 5000 loci | 1.12 / 1.29 / 0.68 / 0.57 | 1.44 / 1.68 / 0.89 / 0.74 |

The pre-move V_G barely moves because before the move 99.5% survive (n̄/N = 0.995), where
the two schemes agree; the existing comparison with the table (median 0.99 with the
rare-alleles moments) never used the function.

The uncomfortable line is the third block. With the correct N_e the steady-state theory
predicts a stronger effect of a cost than the kernel shows: the pre-specified test of
`validate/critical_rate.jl` would have failed C1 (12 of 24) and C2 (1 of 8) instead of
passing both, and "the theory's M is about 18% low" becomes 18–77% low. The reported
agreement owed part of its quality to the overstated N_e (an N_e that falls more slowly
with the breeders makes the feedback weaker, which is the direction of the kernel's
deviation). The qualitative findings stand under either scheme (the cost ratio is far
from `sustainable_rate`'s, and the kernel's ratio is always above the prediction), but
"within 0.054 in 24 of 24 cells" is a statement about the default formula only. On the
other side, k50/k_c at 100 loci moves from 0.77–1.03 to 0.97–1.38, so the kernel's
half-survival rate is then at or above the deterministic critical rate, which a threshold
that ignores the transient (the variance takes time to erode) allows.

**Recommendation.** The derivation and the kernel agree and leave no room for the current
formula when n < N, so `:breeders` should become the default eventually; the drift of
every moments recursion is already that. It is not a small edit, though: it changes
`drift_critical_rate` by 2–28% everywhere it is quoted (CLAUDE.md, the docstrings of
`critical.jl`, `validate/critical_rate.jl`'s header and log, vignettes that print k_c, the
paper draft), and it turns the cost-ratio validation into a partial failure that would have
to be reported as such, ideally with `validate/critical_rate.jl` rerun on fresh seeds under
the new default (its rate grid was centred with the old k_c). Until that is decided the
default stays, the docstrings say what the formula is, and new work can pass
`ne = :breeders` (stage 2 does).

## 1b. `drift_critical_rate` misses steady states where its iteration fails

**Reproduced.** On the 16 Xu–Shen-type inputs of `validate/xu_shen.jl`, 3 differ from the
maximum over all steady states: 0 against 0.2967 (n_e = 100, ω² = 20), 0.2957 against
0.4342 (V_M = 0.2), 0.4114 against 0.4211 (K = 2000).

**Characterised.** It is not a second branch. At a given lag, let f(V) be the house-of-Gauss
variance of a population whose size is set by the load of variance V. f falls with V
(640 of 640 input × lag pairs with N_e = breeders; with the default formula 633 of 640,
the exceptions being haploid inputs with 1.01–1.02 breeders, where f < V), so V = f(V) has
one root. The function iterates V ← f(V) from the variance of N = K, and a decreasing map
can overshoot: (1) if a population of size K with its own variance cannot replace itself
(R·W̄ ≤ 1 at that variance; 0.997 in the first case) the iteration returns "no population"
at once, although a smaller population with less variance is a steady state; (2) if
|f′| > 1 at the root (1.08; 0.87 and 0.58 in the other two, where it fails at some lags
only) it does not converge and the last iterate is returned. Both happen when the variance
of a large population carries a load near ln R.

**Is the missed state real?** In the deterministic recursion with the same ingredients
(the P-D model of stage 2), started at N = K with the variance of N = K, the population
settles on exactly the bisection's steady state (V 5.863, 33.1 breeders for the
Xu–Shen-type input; V 0.1135, 192.8 breeders for a published one) and follows the optimum
at 0.5 and 0.9 of the bisection's k_c, dying at 1.1. In the kernel:

- published input QT5000-3 Input-12144 under Beverton–Holt R = 1.5 (old k_c = 0, "no
  population"): 12 of 12 populations persist 2000 generations at a fixed optimum with
  about 530 adults and 430 breeders, 12 of 12 follow the optimum for 2000 generations at
  half the bisection's k_c (0.038 per generation), 0 of 12 at 1.5 times it (dead within 32
  generations). So "0" was wrong for this input. (The kernel's V_G, 0.064, is below the
  house of Gauss's 0.114, so its load is lighter and its population larger than the steady
  state's 264.)
- the Xu–Shen-type input (K = 106): the steady state has 33 breeders, and no kernel
  population lasts (median time to extinction 320 generations at a fixed optimum), while
  alive with 40 breeders and V_G 4.5 against the predicted 33 and 5.9: the state exists
  but is too small to persist stochastically.

**Published-type values** (600 published inputs × 4 demographies, a fifth of them also
with costs 0.25 and 1; 3360 comparisons). Bisection against the old function, differences
above 10⁻¹² (the rest agree to 3·10⁻¹⁵):

| | no cost | cost 0.25 or 1 |
|---|---|---|
| fixed size | 0 of 600 | 0 of 240 |
| Beverton–Holt R = 5 | 0 of 600 | 0 of 240 |
| R = 2 | 0 of 600 | 2 of 240 |
| R = 1.5 | 5 of 600 | 5 of 240 |

In all 12 the bisection equals an independent maximum over all steady states (parametrised
by V, as eq. 6 of Xu and Shen is) to 4·10⁻¹⁵ and is above the old value (five from 0);
they are small populations at the critical lag (42–173 breeders).

**What was done, and why opt-in.** `drift_critical_rate(p; solver = :bisection)` (and the
internal `_feedback_state`) falls back on bisection where the iteration fails and is
otherwise the same computation. The task allowed making this the default only if every
published-type value stayed unchanged; 12 of 3360 change (none at a fixed size, which is
the published demography, but Beverton–Holt with R = 2 and 1.5 on published genetics is
what `validate/critical_rate.jl` and the density-regulated reruns use), so it is an option and the default
(`solver = :iteration`) is bitwise the old function: 3360 of 3360 comparisons identical in
every field. Recommendation: make `:bisection` the default together with the `ne`
decision; every value it changes was demonstrably not the maximum.

## Tests and checks

`julia --project=vignettes test/runtests.jl theory critical closures emulator xu_shen declarations`
passes; new tests in `test/theory.jl` (the scheme, the neutral limits, the moments
lowering's 2nV_m) and `test/critical.jl` (`ne`, the three fixed cases, uniqueness, the
default unchanged). `validate/symbolic/check_effective_size.jl`: 29 checks.

CPU time: kernel part about 57 CPU-minutes (plus about 12 for a run stopped before
amendment 2), branch part under 1, impact part 8.
