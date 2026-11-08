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

"""
The population ACSet of a model of several traits: the same schema, with vector-valued
allele effects and traits (one entry per trait). Only the attribute types change.
"""
const TraitPopulation = PopulationACSet{Vector{Float64},Vector{Float64}}

function lift(w::World)
    w.traits === nothing || return lift_traits(w)
    p = w.p
    pop = Population()
    N, L, P = w.N, p.NoQTL, p.Ploidy
    add_parts!(pop, :Individual, N)
    G = genomes(w)
    @inbounds for i in 1:N
        set_subpart!(pop, i, :env, w.env[i])
        set_subpart!(pop, i, :genotype, w.genotype[i])
        set_subpart!(pop, i, :phenotype, w.phenotype[i])
        set_subpart!(pop, i, :fitness, w.fitness[i])
        for c in 1:P, l in 1:L
            a = add_part!(pop, :Allele)
            set_subpart!(pop, a, :owner, i)
            set_subpart!(pop, a, :effect, G[l, c, i])
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

function lift_traits(w::World)
    p = w.p
    pop = TraitPopulation()
    N, L, P = w.N, p.NoQTL, p.Ploidy
    add_parts!(pop, :Individual, N)
    G = genomes(w)
    E, Gt, Z = (QTLSim.traitvalues(w; what) for what in (:env, :genotype, :phenotype))
    for i in 1:N
        set_subpart!(pop, i, :env, E[:, i])
        set_subpart!(pop, i, :genotype, Gt[:, i])
        set_subpart!(pop, i, :phenotype, Z[:, i])
        set_subpart!(pop, i, :fitness, [w.fitness[i]])
        for c in 1:P, l in 1:L
            a = add_part!(pop, :Allele)
            set_subpart!(pop, a, :owner, i)
            set_subpart!(pop, a, :effect, G[:, l, c, i])
        end
    end
    pop
end

# a genotype is the sum of the effect vectors of the preimage owner⁻¹(i)
genotypes(pop::TraitPopulation) =
    [reduce(+, (pop[a, :effect] for a in incident(pop, i, :owner))) for i in 1:nparts(pop, :Individual)]
phenotypes(pop::TraitPopulation) = genotypes(pop) .+ pop[:, :env]

function lower!(pop::TraitPopulation, w::World)
    p = w.p
    n, N = p.NoTraits, nparts(pop, :Individual)
    w.traits === nothing && throw(ArgumentError("a population of several traits for a World of one trait"))
    (p.Demography === :fixed ? N == p.PopSize : 1 <= N <= p.PopSize) ||
        throw(ArgumentError("the population has $N individuals; the World holds PopSize = $(p.PopSize)"))
    G = zeros(n, p.NoQTL, p.Ploidy, N)
    for i in 1:N
        as = incident(pop, i, :owner)
        length(as) == p.NoQTL * p.Ploidy || throw(ArgumentError("individual $i has $(length(as)) alleles"))
        for (k, a) in enumerate(as)
            G[:, (k - 1) % p.NoQTL + 1, (k - 1) ÷ p.NoQTL + 1, i] = pop[a, :effect]
        end
    end
    setgenomes!(w, G)
    for i in 1:N
        w.traits.E[:, i] = pop[i, :env]
    end
    w.N = N
    QTLSim.develop!(w)
end

# alleles of each individual in stored order, reshaped column-major to (L, P)
function genome_from_acset(pop::Population, L::Int, P::Int)
    N = nparts(pop, :Individual)
    eff = pop[:, :effect]
    genome = zeros(Float64, L, P, N)
    @inbounds for i in 1:N
        as = incident(pop, i, :owner)
        length(as) == L * P ||
            throw(ArgumentError("individual $i has $(length(as)) alleles, expected $(L*P)"))
        for (k, a) in enumerate(as)
            genome[(k - 1) % L + 1, (k - 1) ÷ L + 1, i] = eff[a]
        end
    end
    genome
end

function lower!(pop::Population, w::World)
    p = w.p
    n = nparts(pop, :Individual)
    (p.Demography === :fixed ? n == p.PopSize : 1 <= n <= p.PopSize) ||
        throw(ArgumentError("the population has $n individuals; the World holds PopSize = $(p.PopSize)" *
                            (p.Demography === :fixed ? "" : " at most")))
    setgenomes!(w, genome_from_acset(pop, p.NoQTL, p.Ploidy))
    w.env[1:n] .= pop[:, :env]
    w.N = n
    QTLSim.develop!(w)                 # keep genotype and phenotype in step with the new genomes
end
