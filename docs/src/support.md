# Support matrix

What QTLSim can run, and what has been checked against another backend, by feature,
model, backend and closure of the moments lowering. **Generated** by
`scripts/support_matrix.jl` from `support_matrix()` (ext/QTLSimCatlabExt/support.jl):
every cell is a probe of the code (four generations of thirty individuals), so a cell
is unsupported exactly where the lowering's own guard throws an `ArgumentError`, and
a cell that runs is *supported* when the declared evidence names a validation script
that compares it with another backend (QTSim, SLiM, the kernel or theory), *untested*
otherwise. A test checks that this page is up to date.

✓ supported and validated · ○ runs, not validated · ✗ unsupported (a guard; see the
note) · – not applicable. Columns: the kernel, the generated SLiM model (the probe
generates the script; it does not run SLiM), and the moments lowering with each closure
(`loci`; `:auto` picks one of them, see [`QTLSim.auto_closure`](@ref) and
[`explain`](@ref)).

## One trait

| feature | kernel | SLiM | aggregate | multilocus | rare alleles | rare local | spectrum |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| QTSim's generation | ✓ [1] | ✓ [2] | ○ | ✓ [3] | ✓ [4] | ○ | ✓ [5] |
| optimum jump (MoveType = :phenotype) | ✓ [6] | ✓ [7] | ○ | ○ | ○ | ○ | ○ |
| continuous move | ✓ [8] | ○ | ○ | ✓ [9] | ○ | ○ | ✓ [5] |
| red-noise optimum | ✓ [10] | ✓ [10] | ○ | ○ | ○ | ○ | ○ |
| periodic optimum | ✓ [10] | ✓ [10] | ○ | ○ | ○ | ○ | ○ |
| width change (ChangeFitnessGen) | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| StandardisePhenotype | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| neutral (Drift) | ✓ [11] | ○ | ○ | ○ | ○ | ○ | ○ |
| haploid | ✓ [12] | ✗ [13] | ○ | ○ | ✓ [12] | ○ | ✗ [14] |
| every crossover point (:all) | ✓ [3] | ✗ [15] | ○ | ✓ [3] | ○ | ○ | ○ |
| Beverton–Holt demography | ✓ [16] | ✓ [16] | ○ | ○ | ○ | ○ | ○ |
| fecundity selection | ✓ [17] | ✓ [17] | ○ | ○ | ○ | ○ | ○ |
| plasticity, fixed slope | ✓ [18] | ✓ [19] | ○ | ○ | ○ | ○ | ○ |
| plasticity with a cost | ✓ [20] | ○ | ○ | ○ | ○ | ○ | ○ |
| plasticity with a delayed cue | ○ | ○ | ○ | ○ | ○ | ○ | ○ |
| a step's own parameters (withparams θ) | ✓ [21] | ✓ [21] | ○ | ○ | ○ | ○ | ○ |
| 3 demes, island migration | ✓ [22] | ✓ [23] | ○ | ✓ [24] | ○ | ○ | ✗ [25] |
| 3 demes, stepping stone (matrix) | ✓ [26] | ✓ [23] | ○ | ○ | ○ | ○ | ✗ [25] |
| 3 demes, their own Params | ✓ [21] | ✓ [21] | ○ | ○ | ○ | ○ | ✗ [25] |
| 3 demes, drift among demes (between) | – | – | ✗ [27] | ✓ [24] | ○ | ✓ [28] | ✗ [29] |

## Two traits

| feature | kernel | SLiM | aggregate | multilocus | rare alleles | rare local | spectrum |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| QTSim's generation | ✓ [30] | ✓ [31] | ○ | ✓ [32] | ✓ [33] | ✓ [34] | ✗ [35] |
| optimum jump (MoveType = :phenotype) | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [35] |
| continuous move | ✓ [31] | ✓ [31] | ○ | ○ | ○ | ○ | ✗ [35] |
| red-noise optimum | ✓ [36] | ✓ [36] | ○ | ○ | ○ | ○ | ✗ [35] |
| periodic optimum | ✓ [36] | ✓ [36] | ○ | ○ | ○ | ○ | ✗ [35] |
| width change (ChangeFitnessGen) | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [35] |
| StandardisePhenotype | ✗ [37] | ✗ [37] | ✗ [37] | ✗ [37] | ✗ [37] | ✗ [37] | ✗ [35] |
| neutral (Drift) | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [35] |
| haploid | ○ | ✗ [13] | ○ | ○ | ○ | ○ | ✗ [35] |
| every crossover point (:all) | ○ | ✗ [15] | ○ | ○ | ○ | ○ | ✗ [35] |
| Beverton–Holt demography | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [35] |
| fecundity selection | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [35] |
| plasticity, fixed slope | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [35] |
| plasticity with a cost | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [35] |
| plasticity with a delayed cue | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [38] | ✗ [35] |
| a step's own parameters (withparams θ) | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [35] |
| 3 demes, island migration | ✓ [39] | ○ | ○ | ✓ [40] | ○ | ○ | ✗ [41] |
| 3 demes, stepping stone (matrix) | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [41] |
| 3 demes, their own Params | ○ | ○ | ○ | ○ | ○ | ○ | ✗ [41] |
| 3 demes, drift among demes (between) | – | – | ✗ [27] | ✓ [40] | ○ | ✓ [42] | ✗ [29] |

## Reaction norm

| feature | kernel | SLiM | aggregate | multilocus | rare alleles | rare local | spectrum |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| QTSim's generation | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| optimum jump (MoveType = :phenotype) | ✗ [46] | ✗ [43] | ✗ [46] | ✗ [46] | ✗ [46] | ✗ [46] | ✗ [35] |
| continuous move | ✓ [47] | ✗ [43] | ○ | ✓ [47] | ✗ [44] | ✗ [45] | ✗ [35] |
| red-noise optimum | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| periodic optimum | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| width change (ChangeFitnessGen) | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| StandardisePhenotype | ✗ [37] | ✗ [43] | ✗ [37] | ✗ [37] | ✗ [37] | ✗ [37] | ✗ [35] |
| neutral (Drift) | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| haploid | ○ | ✗ [13] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| every crossover point (:all) | ○ | ✗ [15] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| Beverton–Holt demography | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| fecundity selection | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| plasticity, fixed slope | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| plasticity with a cost | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| plasticity with a delayed cue | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| a step's own parameters (withparams θ) | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [35] |
| 3 demes, island migration | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [48] |
| 3 demes, stepping stone (matrix) | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [48] |
| 3 demes, their own Params | ○ | ✗ [43] | ○ | ○ | ✗ [44] | ✗ [45] | ✗ [48] |
| 3 demes, drift among demes (between) | – | – | ✗ [27] | ○ | ✗ [44] | ✗ [45] | ✗ [29] |

## Notes

1. validated by gsl_locked.jl, compare_to_qtsim.jl
2. validated by slim/compare_to_qtsim.jl, slim/three_way_confirm.jl
3. validated by moments_loci.jl
4. validated by moments_regimes.jl
5. validated by vg_move.jl
6. validated by compare_to_qtsim.jl
7. validated by slim/compare_to_qtsim.jl
8. validated by published_sweep.jl
9. validated by moments_loci.jl (S)
10. validated by slim/compare_fluctuating.jl
11. validated by fst.jl, coalescence.jl
12. validated by rare_alleles_traits.jl (H)
13. the lowering to Eidos covers diploids (Ploidy = 2); got Ploidy = 1
14. loci = :spectrum is for diploids; Ploidy = 1
15. with more than one locus the lowering to Eidos recombines as QTSim does (CrossoverModel = :qtsim); got :all
16. validated by slim/compare_demography.jl
17. validated by slim/compare_selection.jl
18. validated by slim/compare_plasticity.jl, plasticity_rates.jl
19. validated by slim/compare_plasticity.jl
20. validated by plasticity_rates.jl
21. validated by slim/compare_demes.jl
22. validated by slim/compare_migration.jl, fst.jl, moments_demes.jl
23. validated by slim/compare_migration.jl
24. validated by moments_demes.jl
25. loci = :spectrum is for one deme; got 3
26. validated by slim/compare_migration.jl, fst.jl
27. between = true tracks the haplotype moments of the multilocus model; it has no loci = :aggregate
28. validated by local_adaptation.jl
29. loci = :spectrum is for one deme; it has no between = true
30. validated by rare_alleles_traits.jl, slim/compare_traits.jl
31. validated by slim/compare_traits.jl
32. validated by moments_loci.jl (T), traits_kernel_vs_moments.jl
33. validated by rare_alleles_traits.jl (T)
34. validated by rare_alleles_traits.jl (T; equals :rare_alleles for one deme)
35. loci = :spectrum is for one trait; NoTraits = 2
36. validated by slim/compare_traits_fluctuating.jl
37. Params: StandardisePhenotype is for one trait
38. Params: Plasticity with NoTraits > 1 needs ReactionNorm = true (and NoTraits = 2)
39. validated by trait_moments_demes.jl, rare_alleles_traits.jl (D)
40. validated by trait_moments_demes.jl
41. loci = :spectrum is for one trait in one deme; NoTraits = 2, 3 demes
42. validated by local_adaptation.jl (D1, D2)
43. a heritable reaction norm (ReactionNorm) has no lowering to Eidos yet; the kernel and the moments have
44. the moments of a heritable reaction norm are Gaussian (loci = :multilocus or :aggregate); got :rare_alleles
45. the moments of a heritable reaction norm are Gaussian (loci = :multilocus or :aggregate); got :rare_local
46. Params: plasticity needs an optimum with its own trajectory (MoveType :continuous, :red_noise or :periodic); :phenotype re-anchors the optimum to the mean pheno
47. validated by reaction_norm.jl
48. the moments of a heritable reaction norm are Gaussian (loci = :multilocus or :aggregate); got :spectrum
