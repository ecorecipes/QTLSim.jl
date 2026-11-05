# AlgebraicJulia data model: the population as an attributed C-set (ACSet).
#
# This is the *abstract* representation of population state, following notes.md §1.
# The fast array kernel in `model.jl` is the *lowered* backend; `lift`/`lower!`
# translate between the two levels, realizing the abstract→lowered architecture.
#
# Schema (a finitely-presented category with attributes):
#
#     Allele ──owner──▶ Individual        effect : Allele → ℝ
#                                         env, genotype, phenotype, fitness : Individual → ℝ
#
# An Individual owns `Ploidy·NoQTL` Alleles; its genotype is the sum of its
# alleles' effects (QTSim `CaclGenotype`), and phenotype = genotype + env.

using Catlab, Catlab.CategoricalAlgebra
using ACSets

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

"""
    Population

Concrete ACSet type for the population schema, with `Float64` allele effects and
traits. Construct an empty one with `Population()`, or build one from a `World`
with [`lift`](@ref).
"""
const Population = PopulationACSet{Float64,Float64}

# --- lift: array World  →  ACSet ----------------------------------------

"""
    lift(w::World) :: Population

Build the ACSet representation of a `World`'s current parent generation. Each
individual becomes an `Individual` part owning `Ploidy·NoQTL` `Allele` parts
whose `effect`s are the entries of `w.genome[:, :, i]`. Env/genotype/phenotype/
fitness attributes are copied from the World's scratch vectors (run `develop!`
first if you want genotype/phenotype populated).
"""
function lift(w::World)::Population
    p = w.p
    pop = Population()
    N, L, P = p.PopSize, p.NoQTL, p.Ploidy
    inds = add_parts!(pop, :Individual, N)
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

# --- genotype accessor via the categorical `incident` query -------------

"""
    genotypes(pop::Population) :: Vector{Float64}

Genotype of every individual, computed by summing the effects of its incident
alleles — `incident(pop, i, :owner)` is the categorical preimage `owner⁻¹(i)`.
This mirrors QTSim `CaclGenotype` at the ACSet level.
"""
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

"""
    phenotypes(pop::Population) :: Vector{Float64}

Phenotype = genotype + env, computed at the ACSet level.
"""
phenotypes(pop::Population) = genotypes(pop) .+ pop[:, :env]

# --- lower: ACSet  →  array genome --------------------------------------

"""
    genome_from_acset(pop::Population, L, P) :: Array{Float64,3}

Reconstruct the `(NoQTL, Ploidy, PopSize)` array genome from the ACSet. Alleles
of each individual are taken in their stored order and reshaped column-major to
`(L, P)`, inverting [`lift`](@ref).
"""
function genome_from_acset(pop::Population, L::Int, P::Int)
    N = nparts(pop, :Individual)
    eff = pop[:, :effect]
    genome = zeros(Float64, L, P, N)
    @inbounds for i in 1:N
        as = incident(pop, i, :owner)
        length(as) == L * P || error("individual $i has $(length(as)) alleles, expected $(L*P)")
        for (k, a) in enumerate(as)
            l = (k - 1) % L + 1
            c = (k - 1) ÷ L + 1
            genome[l, c, i] = eff[a]
        end
    end
    genome
end

"""
    lower!(pop::Population, w::World) :: World

Write the ACSet's genome and environment back into a `World`'s array buffers,
the inverse of [`lift`](@ref). Returns the updated World.
"""
function lower!(pop::Population, w::World)
    p = w.p
    w.genome = genome_from_acset(pop, p.NoQTL, p.Ploidy)
    w.env = copy(pop[:, :env])
    w
end
