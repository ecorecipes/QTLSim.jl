# The life cycle, written once as a morphism of a free symmetric monoidal category, and
# its lowering to the simulation kernel (compose.jl) and to SLiM's Eidos (eidos.jl).
#
# The objects of the category are the stages of a generation and its generators are the
# steps between them:
#
#     develop:   Zygotes → Adults           record:  Offspring → Offspring
#     select:    Adults → Survivors         rescale: Offspring → Offspring
#     reproduce: Survivors → Offspring      mutate:  Offspring → Offspring
#                                           swap:    Offspring → Zygotes
#
# QTSim's generation is the composite develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅
# mutate ⋅ swap: Zygotes → Zygotes, in the order of its generation loop. The monoidal
# product ⊗ puts demes side by side: Zygotes⊗Zygotes is a metapopulation of two demes at
# the start of a generation, generation ⊗ generation lets each evolve on its own, and
# `migrate(k)`: Zygotes^⊗k → Zygotes^⊗k couples k of them. A functor out of a free
# category is fixed by where it sends the generators, so a backend is a table of images,
# one per step, and the image of any composite or product is the composite or product of
# the images (lean/CategoricalQuantGen/Lowering.lean, Demes.lean). Composites whose
# stages do not match are not morphisms at all: `select ⋅ develop` is a domain error.

@present LifeCycle(FreeSymmetricMonoidalCategory) begin
    (Zygotes, Adults, Survivors, Offspring)::Ob
    develop::Hom(Zygotes, Adults)
    select::Hom(Adults, Survivors)
    reproduce::Hom(Survivors, Offspring)
    record::Hom(Offspring, Offspring)
    rescale::Hom(Offspring, Offspring)
    mutate::Hom(Offspring, Offspring)
    swap::Hom(Offspring, Zygotes)
end

const STEPS = (:develop, :select, :reproduce, :record, :rescale, :mutate, :swap)

"The number of stages, one per deme, in a product of stages."
nfactors(x) = head(x) === :otimes ? length(args(x)) : head(x) === :munit ? 0 : 1

"Zygotes⊗…⊗Zygotes with `k` factors: `k` demes at the start of a generation."
zygotes_word(k::Integer) = k == 1 ? LifeCycle[:Zygotes] : foldl(otimes, fill(LifeCycle[:Zygotes], k))

"The step that exchanges offspring among `k` demes."
function migrate_step(k::Integer)
    k >= 2 || throw(ArgumentError("migration needs at least two demes; got $k"))
    Hom(:migrate, zygotes_word(k), zygotes_word(k))
end

function lifecycle()
    steps = map(s -> LifeCycle[s], STEPS)
    NamedTuple{(:presentation, STEPS..., :generation, :migrate)}(
        (LifeCycle, steps..., foldl(⋅, steps), migrate_step))
end

"""
A product of stages, one per deme, as the image of an object under the lowering to
backend `B` (`:kernel` or `:eidos`). The kernel carries each deme's stages in one `World`;
in Eidos, a stage is the set of variables a step can rely on (`EIDOS_SCOPE`).
"""
struct Stage{B}
    names::Vector{Symbol}
end

Stage{B}(name::Symbol) where {B} = Stage{B}([name])
nfactors(s::Stage) = length(s.names)

otimes(a::Stage{B}, b::Stage{B}) where {B} = Stage{B}(vcat(a.names, b.names))
otimes(a::Stage{B}, b::Stage{B}, c::Stage{B}, rest::Stage{B}...) where {B} = foldl(otimes, (a, b, c, rest...))
munit(::Type{Stage{B}}) where {B} = Stage{B}(Symbol[])

id(s::Stage{:kernel}) = nfactors(s) == 1 ? Operator(:id, (w, rng) -> w) :
                        DemeOp(nfactors(s), collect(1:nfactors(s)), (st, d, rng) -> st)
id(s::Stage{:eidos}) = nfactors(s) == 1 ? EidosCode(Symbol(:id_, s.names[1]), "") :
                       EidosMeta(nfactors(s), collect(1:nfactors(s)), pops -> "", :id)

"Where a braiding X⊗Y → Y⊗X sends the demes of X (m of them) and Y (n of them)."
braid_perm(m, n) = vcat(collect(m+1:m+n), collect(1:m))
braid(a::Stage{:kernel}, b::Stage{:kernel}) =
    DemeOp(nfactors(a) + nfactors(b), braid_perm(nfactors(a), nfactors(b)), (st, d, rng) -> st)
braid(a::Stage{:eidos}, b::Stage{:eidos}) =
    EidosMeta(nfactors(a) + nfactors(b), braid_perm(nfactors(a), nfactors(b)), pops -> "", :braid)

"The generators occurring in a life-cycle expression."
stepsof(f) = head(f) === :generator ? Any[f] :
             reduce(vcat, (stepsof(a) for a in args(f) if a isa GATExpr); init = Any[])

asimage(::Val{:kernel}, name, x::Union{Operator,DemeOp}) = x
asimage(::Val{:kernel}, name, x) = Operator(name, x)               # a function (state, rng) -> state
asimage(::Val{:eidos}, name, x::Union{EidosCode,EidosMeta}) = x
asimage(::Val{:eidos}, name, x::AbstractString) = EidosCode(name, x)
asimage(::Val{:moments}, name, x::Operator) = x
asimage(::Val{:moments}, name, x) = Operator(name, x)               # a function (moments, rng) -> moments

function lower(f, target::Symbol; images = Dict{Symbol,Any}(), record = nothing,
               migration_rate::Union{Real,AbstractMatrix} = 0.0)
    table = target === :kernel ? kernel_images(record) :
            target === :eidos ? eidos_images() :
            target === :moments ? moment_images(record) :
            throw(ArgumentError("unknown target :$target; the targets are :kernel, :eidos and :moments"))
    merge!(table, Dict{Symbol,Any}(images))
    gens = Dict{Any,Any}()
    for g in stepsof(f)
        name = nameof(g)
        if g isa ObExpr
            gens[g] = Stage{target}(name)
        elseif haskey(table, name)
            gens[g] = asimage(Val(target), name, table[name])
        elseif name === :migrate
            k = nfactors(dom(g))
            gens[g] = target === :kernel ? migrate_op(k, migration_rate) :
                      target === :eidos ? eidos_migrate(k, migration_rate isa AbstractMatrix) :
                      throw(ArgumentError("the moments lowering covers one deme, so it has no image of migrate"))
        else
            throw(ArgumentError("no $target image for the step $name; give one with `images = Dict(:$name => …)`"))
        end
    end
    functor((Stage{target}, Any), f; generators = gens)
end

function operators()
    k = kernel_images()
    NamedTuple{(STEPS..., :generation)}((map(s -> k[s], STEPS)..., lower(lifecycle().generation, :kernel)))
end
