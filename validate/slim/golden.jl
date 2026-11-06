# The SLiM models committed in this folder, as `slim_script` writes them from QTLSim's
# life cycle. The package tests check that `slim_script` still writes each one unchanged
# (without SLiM); check_generated.jl checks the same, regenerates them with `--write`,
# runs the first against the hand-written port, run for run, and runs the others in SLiM.
# Needs QTLSim and Catlab loaded.

"The committed SLiM models: file name => Eidos source, as `slim_script` writes it now."
function golden_models()
    L = lifecycle(); g = L.generation
    # scenario B of compare_to_qtsim.jl: one locus, the optimum jumps at 5000
    B = Params(PopSize = 1000, NoQTL = 1, Ploidy = 2, NoGenerations = 6500, PhenotypeVar = 0.02,
               Heritability = 0.5, MutationScalar = 0.001, MutationRatePerGamete = 0.0002,
               MoveFitnessSDGen = 5000, MoveFitnessSD = 1.0, MoveType = :phenotype,
               ChangeFitnessGen = 5000, ChangeFitnessSD = 0.30)
    # scenario D: 100 loci with QTSim's recombination
    D = Params(PopSize = 1000, NoQTL = 100, Ploidy = 2, NoGenerations = 2000, PhenotypeVar = 0.03,
               Heritability = 0.4, MutationScalar = 0.005, MutationRatePerGamete = 0.0001,
               RecombPoissonMean = 3.0, InitFitnessSD = 4.0)
    # three demes coupled by migration
    M = Params(PopSize = 500, NoQTL = 1, Ploidy = 2, NoGenerations = 1000, PhenotypeVar = 0.02,
               Heritability = 0.5, MutationScalar = 0.02, MutationRatePerGamete = 0.5, InitFitnessSD = 0.3)
    ["qtsim_generated.slim" => slim_script(g, B),
     "qtsim_generated_loci.slim" => slim_script(g, D),
     "qtsim_generated_demes.slim" => slim_script(Catlab.compose(Catlab.otimes(g, g, g), L.migrate(3)), M;
                                                 migration_rate = 0.01)]
end
