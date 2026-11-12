# Morphisms between models, and the order of the closures, declared as data (roadmap C7;
# the types are in src/closures.jl). Each declaration says what it maps, when its
# reductions hold and which Lean theorem (if any) proves them, and lists the reductions on
# small cases; test/runtests.jl generates one test per reduction ("declared model
# morphisms"). Several of them repeat hand tests written before the declarations
# (`covers`); those stay.

using QTLSim: ModelMorphism, Reduction
import QTLSim: model_morphisms, closure_order

# small cases, fast enough for the tests
const MORPH_P = Params(PopSize = 200, NoQTL = 20, Ploidy = 2, NoGenerations = 400, PhenotypeVar = 1.0, Heritability = 0.5,
                       MutationScalar = 0.01, MutationRatePerGamete = 0.01, InitFitnessSD = 2.0, Sample = 10)
const MORPH_MOVE = Params(MORPH_P; MoveFitnessSDGen = 200, MoveType = :continuous, MoveFitnessSD = 0.01)
g11(o) = [G[1, 1] for G in o.gmatrix]

function model_morphisms()
    L = lifecycle(); g = L.generation
    p, pm = MORPH_P, MORPH_MOVE
    two(q; kw...) = Params(q; NoTraits = 2, MutCorr = 0.0, SelCorr = 0.0, kw...)
    ModelMorphism[
        ModelMorphism("one trait ↪ two independent traits", "OneTraitModel()", "TraitModel(2)", :embedding,
            "a trait is trait 1 of two traits that share nothing: uncorrelated mutations (MutCorr = 0), a fitness " *
            "function without correlation (SelCorr = 0), the optimum moving along trait 1 (MoveDirection = [1, 0])",
            "the moments lowering of the multilocus model without drift (erosion = false): the second trait's selection " *
            "changes only the expected survivors N·W̄, which drift and the population size read",
            "not in Lean (the 2 × 2 recursion with diagonal matrices is two copies of the 1 × 1 one)",
            [Reduction("G₁₁ is the one-trait V_G at a fixed optimum", :moments, :tolerance, 1e-12,
                       () -> (g11(run_moments(two(p); loci = :multilocus))[2:end], run_moments(p; loci = :multilocus).geno_var[2:end])),
             Reduction("G₁₁ and the mean of trait 1 under a continuous move along trait 1", :moments, :tolerance, 1e-12,
                       () -> (let a = run_moments(two(pm; MoveDirection = [1.0, 0.0]); loci = :multilocus)
                                  (g11(a)[2:end], first.(a.trait_mean)[2:end])
                              end,
                              let b = run_moments(pm; loci = :multilocus)
                                  (b.geno_var[2:end], b.pheno_mean[2:end])
                              end)),
             Reduction("the n-trait recursion with one trait (1 × 1 matrices) is the one-trait recursion", :moments, :tolerance, 1e-10,
                       () -> ((x = run_trait_moments(pm; loci = :multilocus); (x.geno_var, x.pheno_mean, x.theta)),
                              (y = run_moments(pm; loci = :multilocus); (y.geno_var, y.pheno_mean, y.theta)));
                       covers = "several traits: the moments lowering, a second table for the same syntax")]),
        ModelMorphism("n traits → trait 1 (QTSim's columns)", "TraitModel(n)", "OneTraitModel() records", :projection,
            "the records of n traits project to QTSim's columns: pheno_mean, geno_var and theta are those of trait 1",
            "always (a statistic map, not a model reduction)",
            "not in Lean (a definition of the records)",
            [Reduction("kernel: pheno_mean and geno_var are trait 1's mean phenotype and G₁₁", :kernel, :tolerance, 1e-12,
                       () -> (o = simulate(two(pm; MutCorr = 0.4, SelCorr = -0.2, MoveDirection = [1.0, 1.0], NoGenerations = 300); seed = 3);
                              ((o.pheno_mean, o.geno_var), (first.(o.trait_mean), g11(o))));
                       covers = "several traits: parameters and the kernel"),
             Reduction("moments: pheno_mean and geno_var are trait 1's", :moments, :bitwise, 0.0,
                       () -> (o = run_moments(two(pm; MutCorr = 0.4, SelCorr = -0.2, MoveDirection = [1.0, 1.0]));
                              ((o.pheno_mean, o.geno_var), (first.(o.trait_mean), g11(o)))))]),
        ModelMorphism("k demes → one deme", "k demes (Zygotes^⊗k)", "one deme", :restriction,
            "a metapopulation whose demes are not coupled is k copies of one deme, and one deme with the statistics " *
            "among demes is one deme",
            "no migration (no migrate step, or rate 0), the same parameters in every deme; for the kernel, one deme",
            "Demes.lean: demes_iterate, tensor_iterate (demes without migration iterate independently); " *
            "Differentiation.lean for the neutral between-deme recursion",
            [Reduction("kernel: a metapopulation of one deme (p ⊗ demes(1)) is the World", :kernel, :bitwise, 0.0,
                       () -> (a = run_composed(World(p), Xoshiro(5)); b = run_composed(p ⊗ demes(1), Xoshiro(5); cycle = g);
                              ((a.geno_var, a.pheno_mean, a.nsurv), (b.geno_var, b.pheno_mean, b.nsurv)))),
             Reduction("moments: three demes without a migrate step are three copies of one deme (erosion)", :moments, :bitwise, 0.0,
                       () -> (a = run_moments(p; erosion = true); b = run_moments(p; cycle = g ⊗ g ⊗ g, erosion = true);
                              (b.geno_var[b.deme .== 2], a.geno_var))),
             Reduction("moments: migrate at rate 0 is the identity (erosion; :auto is :rare_local for 3 demes, :rare_alleles for 1)",
                       :moments, :bitwise, 0.0,
                       () -> (a = run_moments(p; erosion = true);
                              b = run_moments(p; cycle = g ⊗ g ⊗ g ⋅ L.migrate(3), migration_rate = 0.0, erosion = true);
                              (b.geno_var[b.deme .== 3], a.geno_var))),
             Reduction("moments: between = true for one deme is erosion = true", :moments, :bitwise, 0.0,
                       () -> (run_moments(p; between = true).geno_var, run_moments(p; erosion = true).geno_var);
                       covers = "the moments of several demes with drift among them (between = true)")]),
        ModelMorphism("plastic at rate k → plain at rate (1 − φ)k", "OneTraitModel(), Plasticity = b, OptimumSlope = B",
            "OneTraitModel() without plasticity", :conjugacy,
            "with φ = b/B, a population expressing (b/B)·θ of its phenotype, under a continuous move of k per generation, " *
            "is a population without plasticity under a move of (1 − φ)k; the map of states subtracts φθ from every " *
            "expressed phenotype, so genotypes, survivors, population sizes and the lag θ − z̄ are the same",
            "CueLag = 0, no cost (PlasticityCostSD = Inf), MoveType = :continuous",
            "Plasticity.lean: lower_conj, plastic_generation_conj, plastic_run_conj, plastic_statistic_law",
            [Reduction("kernel, the same seed: survivors and V_G", :kernel, :bitwise, 0.0,
                       () -> (a = simulate(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.04,
                                                  Plasticity = 1.5, OptimumSlope = 2.0); seed = 6);
                              b = simulate(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.01); seed = 6);
                              ((a.nsurv, a.geno_var), (b.nsurv, b.geno_var)));
                       covers = "phenotypic plasticity: a reaction norm of fixed slope"),
             Reduction("kernel, the same seed: mean fitness and the lag", :kernel, :tolerance, 1e-10,
                       () -> (a = simulate(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.04,
                                                  Plasticity = 1.5, OptimumSlope = 2.0); seed = 6);
                              b = simulate(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.01); seed = 6);
                              ((a.mean_fitness, a.theta .- a.pheno_mean), (b.mean_fitness, b.theta .- b.pheno_mean)))),
             Reduction("moments: V_G, survivors and the lag", :moments, :tolerance, 1e-12,
                       () -> (a = run_moments(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.04,
                                                     Plasticity = 1.5, OptimumSlope = 2.0));
                              b = run_moments(Params(p; MoveFitnessSDGen = 100, MoveType = :continuous, MoveFitnessSD = 0.01));
                              ((a.geno_var, a.nsurv, a.theta .- a.pheno_mean), (b.geno_var, b.nsurv, b.theta .- b.pheno_mean))))])]
end

function closure_order()
    L = lifecycle(); g = L.generation
    p = Params(Ploidy = 2, NoQTL = 20, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
               MutationRatePerGamete = 2e-4, PopSize = 100, NoGenerations = 600, Sample = 10, InitFitnessSD = 5.0)
    p2 = Params(p; NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)
    cyc(θs) = foldl(⊗, [withparams(g, p; θ = t) for t in θs]) ⋅ L.migrate(length(θs))
    within(l, θs, rate) = run_moments(p; cycle = cyc(θs), migration_rate = rate, between = true, loci = l).metapop.within
    sp = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
                MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0)
    ModelMorphism[
        ModelMorphism(":rare_local → :rare_alleles", ":rare_local", ":rare_alleles", :refinement,
            "the closure with locally favoured alleles uses [s⁻¹ − 2Λ]₊ for s⁻¹; Λ = 0 reduces it to the rare-alleles closure",
            "one deme, demes that lag alike (a common optimum), or no migration (Λ = 0)",
            "LocalAdaptation.lean: local_fixed_point_zero, gain_const, gain_zero_weights",
            [Reduction("one deme (erosion)", :moments, :tolerance, 1e-12,
                       () -> (run_moments(p; erosion = true, loci = :rare_local).geno_var,
                              run_moments(p; erosion = true, loci = :rare_alleles).geno_var);
                       covers = "rare alleles with locally favoured alleles"),
             Reduction("one deme, two traits (erosion)", :moments, :tolerance, 1e-12,
                       () -> (run_moments(p2; erosion = true, loci = :rare_local).gmatrix,
                              run_moments(p2; erosion = true, loci = :rare_alleles).gmatrix)),
             Reduction("four demes, a common optimum, migration 0.01", :moments, :tolerance, 1e-12,
                       () -> (within(:rare_local, zeros(4), 0.01), within(:rare_alleles, zeros(4), 0.01))),
             Reduction("four demes, local optima, no migration", :moments, :tolerance, 1e-12,
                       () -> (within(:rare_local, [2.0, 1, -1, -2], 0.0), within(:rare_alleles, [2.0, 1, -1, -2], 0.0)))]),
        ModelMorphism(":rare_alleles → :multilocus", ":rare_alleles", ":multilocus", :refinement,
            "the rare-alleles closure adds the loss β²α²/(2s) of the genic variance; α² → 0 (many small mutations) removes it",
            "a limit: α² = V_m/(Ploidy·L·μ) → 0 as μ → ∞; at μ = 1 per locus the equilibria differ by less than 0.1%",
            "RareAlleles.lean: rare_gaussian_limit (a filter limit); rareFixedPoint_le_gauss, houseOfGauss_le_gaussian",
            [Reduction("equilibrium V_G at μ = 1 (no drift and with drift)", :theory, :tolerance, 1e-3,
                       () -> (q = Params(MORPH_P; MutationRatePerGamete = 1.0, PopSize = 1000);
                              ((rare_alleles_vg(q), rare_alleles_vg(q; drift = true)),
                               (gaussian_vg_loci(q), gaussian_vg_loci(q; drift = true)))))]),
        ModelMorphism(":multilocus → :aggregate", ":multilocus", ":aggregate", :refinement,
            "the multilocus model's haplotype moments carry LD that recombination breaks up; a haploid genome is inherited " *
            "whole (clones), so all of selection's reduction is inherited, as the aggregate recursion assumes",
            "haploid (any number of loci), without drift (the two models put drift in differently)",
            "Theory.lean (the aggregate fixed point gaussian_fixed_point); Multilocus.lean (multilocus_fixed_point_iff); " *
            "the haploid reduction itself is not in Lean",
            [Reduction("haploid, one locus", :moments, :tolerance, 1e-12,
                       () -> (q = Params(MORPH_P; Ploidy = 1, NoQTL = 1);
                              (run_moments(q; loci = :multilocus).geno_var[2:end], run_moments(q; loci = :aggregate).geno_var[2:end]))),
             Reduction("haploid clones, ten loci", :moments, :tolerance, 1e-12,
                       () -> (q = Params(MORPH_P; Ploidy = 1, NoQTL = 10);
                              (run_moments(q; loci = :multilocus).geno_var[2:end], run_moments(q; loci = :aggregate).geno_var[2:end])))]),
        ModelMorphism(":spectrum → :rare_alleles", ":spectrum", ":rare_alleles", :refinement,
            "at a fixed optimum the spectrum's stationary V_G is close to the rare-alleles closure in linkage equilibrium " *
            "(house_of_gauss_vg), and neutral it holds 2N·V_m",
            "approximate: rare mutations (4Nμ ≲ 1 per locus), a fixed optimum, one deme",
            "not in Lean (the spectrum's diffusion is not formalised)",
            [Reduction("stationary V_G within 5% of house_of_gauss_vg", :theory, :tolerance, 0.05,
                       () -> (e = spectrum_equilibrium(sp); (e.VG, house_of_gauss_vg(sp; nsurv = e.nparents)));
                       covers = "allele-frequency spectrum (loci = :spectrum)"),
             Reduction("neutral: 2N·V_m within 1%", :theory, :tolerance, 0.01,
                       () -> (spectrum_equilibrium(Params(sp; Drift = true)).VG, 2 * 1000 * muteVar(sp)))])]
end
