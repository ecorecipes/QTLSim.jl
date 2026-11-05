# The population as an attributed C-set (ACSet): the abstract counterpart of the
# kernel's arrays. `lift` and `lower!` translate between the two.
#
#     Allele ──owner──▶ Individual        effect : Allele → ℝ
#                                         env, genotype, phenotype, fitness : Individual → ℝ
#
# An individual owns `Ploidy·NoQTL` alleles; its genotype is the sum of its alleles'
# effects (QTSim `CaclGenotype`) and its phenotype is genotype + env.

@present SchPopulation(FreeSchema) begin
    Individual::Ob
    Allele::Ob
    owner::Hom(Allele, Individual)

    Effect::AttrType
    effect::Attr(Allele, Effect)

    Trait::AttrType
    env::Attr(Individual, Trait)
    genotype::Attr(Individual, Trait)
    phenotype::Attr(Individual, Trait)
    fitness::Attr(Individual, Trait)
end

@acset_type PopulationACSet(SchPopulation, index=[:owner])

"The population ACSet, with `Float64` allele effects and traits."
const Population = PopulationACSet{Float64,Float64}

population_schema() = SchPopulation

function lift(w::World)
    p = w.p
    pop = Population()
    N, L, P = p.PopSize, p.NoQTL, p.Ploidy
    add_parts!(pop, :Individual, N)
    @inbounds for i in 1:N
        set_subpart!(pop, i, :env, w.env[i])
        set_subpart!(pop, i, :genotype, w.genotype[i])
        set_subpart!(pop, i, :phenotype, w.phenotype[i])
        set_subpart!(pop, i, :fitness, w.fitness[i])
        for c in 1:P, l in 1:L
            a = add_part!(pop, :Allele)
            set_subpart!(pop, a, :owner, i)
            set_subpart!(pop, a, :effect, w.genome[l, c, i])
        end
    end
    pop
end

# genotype via the preimage owner⁻¹(i), which `incident` returns
function genotypes(pop::Population)
    N = nparts(pop, :Individual)
    eff = pop[:, :effect]
    g = zeros(Float64, N)
    @inbounds for i in 1:N
        for a in incident(pop, i, :owner)
            g[i] += eff[a]
        end
    end
    g
end

phenotypes(pop::Population) = genotypes(pop) .+ pop[:, :env]

# alleles of each individual in stored order, reshaped column-major to (L, P)
function genome_from_acset(pop::Population, L::Int, P::Int)
    N = nparts(pop, :Individual)
    eff = pop[:, :effect]
    genome = zeros(Float64, L, P, N)
    @inbounds for i in 1:N
        as = incident(pop, i, :owner)
        length(as) == L * P || error("individual $i has $(length(as)) alleles, expected $(L*P)")
        for (k, a) in enumerate(as)
            genome[(k - 1) % L + 1, (k - 1) ÷ L + 1, i] = eff[a]
        end
    end
    genome
end

function lower!(pop::Population, w::World)
    p = w.p
    w.genome = genome_from_acset(pop, p.NoQTL, p.Ploidy)
    w.env = copy(pop[:, :env])
    w
end
