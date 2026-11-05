# The life cycle, written once as a morphism of a free category, and its lowering to the
# simulation kernel (compose.jl) and to SLiM's Eidos (eidos.jl).
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
# mutate ⋅ swap: Zygotes → Zygotes, in the order of its generation loop. A functor out of a
# free category is fixed by where it sends the generators, so a backend is a table of
# images, one per step, and the image of any composite is the composite of the images
# (lean/CategoricalQuantGen/Lowering.lean). Composites whose stages do not match are not
# morphisms at all: `select ⋅ develop` is a domain error.

@present LifeCycle(FreeCategory) begin
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

function lifecycle()
    steps = map(s -> LifeCycle[s], STEPS)
    NamedTuple{(:presentation, STEPS..., :generation)}((LifeCycle, steps..., foldl(⋅, steps)))
end

"""
A stage of the life cycle as the image of an object under the lowering to backend `B`
(`:kernel` or `:eidos`). The kernel carries every stage in one `World`; in Eidos, a stage
is the set of variables a step can rely on (`EIDOS_SCOPE`).
"""
struct Stage{B}
    name::Symbol
end

id(::Stage{:kernel}) = Operator(:id, (w, rng) -> w)
id(s::Stage{:eidos}) = EidosCode(Symbol(:id_, s.name), "")

"The generators occurring in a life-cycle expression."
stepsof(f) = head(f) === :generator ? Any[f] :
             reduce(vcat, (stepsof(a) for a in args(f) if a isa GATExpr); init = Any[])

asimage(::Val{:kernel}, name, x::Operator) = x
asimage(::Val{:kernel}, name, x) = Operator(name, x)               # a function (state, rng) -> state
asimage(::Val{:eidos}, name, x::EidosCode) = x
asimage(::Val{:eidos}, name, x::AbstractString) = EidosCode(name, x)

function lower(f, target::Symbol; images = Dict{Symbol,Any}(), record = nothing)
    table = target === :kernel ? kernel_images(record) :
            target === :eidos ? eidos_images() :
            throw(ArgumentError("unknown target :$target; the targets are :kernel and :eidos"))
    merge!(table, Dict{Symbol,Any}(images))
    gens = Dict{Any,Any}()
    for g in stepsof(f)
        name = nameof(g)
        if g isa ObExpr
            gens[g] = Stage{target}(name)
        else
            haskey(table, name) || throw(ArgumentError(
                "no $target image for the step $name; give one with `images = Dict(:$name => …)`"))
            gens[g] = asimage(Val(target), name, table[name])
        end
    end
    functor((Stage{target}, target === :kernel ? Operator : EidosCode), f; generators = gens)
end

function operators()
    k = kernel_images()
    NamedTuple{(STEPS..., :generation)}((map(s -> k[s], STEPS)..., lower(lifecycle().generation, :kernel)))
end
