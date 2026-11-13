# The targets of the lowerings as categories (roadmap C1).
#
# A lowering is a functor out of the free symmetric monoidal category of the life cycle
# (lifecycle.jl). Its target was a pair (Stage{B}, Any): objects were typed, morphisms
# were whatever image a table held, and a wrong image (an operator on the wrong number of
# demes, a function of the wrong arity) surfaced only when a run reached it. Here a
# morphism of the target is a `Morph{B}`: an image together with its domain and codomain
# stages. Each generator's image is checked when the life cycle is lowered (deme count,
# calling convention), composites check that the stages meet, and `LoweringCategory{B}`
# declares the target as a model of Catlab's theory of symmetric monoidal categories. The
# laws of that theory (associativity, units, interchange, braiding) are property-tested
# on random composites in test/laws.jl: exactly where the image is deterministic, and in
# law where the order of random draws differs (Lean: `kernelProd_interchange`,
# `stochastic_lowerings_agree`).

"""
    Morph{B}(dom, codom, image)

A morphism of the target of the lowering to backend `B` (`:kernel`, `:eidos` or
`:moments`): an image (an `Operator` or `DemeOp` for the kernel and the moments lowering,
an `EidosCode` or `EidosMeta` for Eidos) from the stages `dom` to the stages `codom`, one
stage per deme. `lower(f, B; typed = true)` returns one; by default `lower` returns its
image.
"""
struct Morph{B}
    dom::Stage{B}
    codom::Stage{B}
    image::Any
end

Base.show(io::IO, f::Morph{B}) where {B} =
    print(io, "Morph{:", B, "}(", join(f.dom.names, "⊗"), " → ", join(f.codom.names, "⊗"), ", ", f.image, ")")

# how many demes an image acts on (`nothing`: an operator on a whole metapopulation, which
# does not say), and whether it can be called as the backend calls it
image_arity(x::Operator) = x.scope === :deme ? 1 : nothing
image_arity(x::DemeOp) = x.n
image_arity(x::EidosCode) = 1
image_arity(x::EidosMeta) = x.n
image_arity(x) = nothing
image_callable(x::Operator) = hasmethod(x.act, Tuple{Any,Any})
image_callable(x::DemeOp) = hasmethod(x.act, Tuple{Any,Any,Any})
image_callable(x::EidosMeta) = hasmethod(x.emit, Tuple{Any})
image_callable(x) = true

const IMAGE_TYPES = Dict(:kernel => Union{Operator,DemeOp}, :moments => Union{Operator,DemeOp},
                         :eidos => Union{EidosCode,EidosMeta})

"""
    typed_image(B, dom, codom, image, name)

The image of the generator `name` as a morphism `dom → codom` of the target of the
lowering to `B`, after checking it: of an image type of that target, acting on as many
demes as `dom` has stages, and callable as the backend calls it (`(state, rng)` for an
`Operator`, `(state, demes, rng)` for a `DemeOp`, `(pops)` for an `EidosMeta`).
"""
function typed_image(::Val{B}, dom::Stage{B}, codom::Stage{B}, x, name) where {B}
    T = IMAGE_TYPES[B]
    x isa T || throw(ArgumentError("the $B image of $name is a $(typeof(x)), not a morphism of the $B target ($T)"))
    n = nfactors(dom)
    nfactors(codom) == n ||
        throw(ArgumentError("$name goes from $n to $(nfactors(codom)) stages; every step keeps the number of demes"))
    a = image_arity(x)
    a === nothing || a == n ||
        throw(ArgumentError("the $B image of $name acts on $a deme$(a == 1 ? "" : "s"), but $name goes " *
                            "$(join(dom.names, "⊗")) → $(join(codom.names, "⊗")) ($n deme$(n == 1 ? "" : "s"))"))
    image_callable(x) ||
        throw(ArgumentError("the $B image of $name cannot be called as the $B backend calls it " *
                            "(an Operator's act takes (state, rng), a DemeOp's (state, demes, rng))"))
    Morph{B}(dom, codom, x)
end

function morph_compose(f::Morph{B}, g::Morph{B}) where {B}
    f.codom.names == g.dom.names ||
        throw(ArgumentError("ill-typed composite in the $B target: $(join(f.codom.names, "⊗")) ≠ $(join(g.dom.names, "⊗"))"))
    Morph{B}(f.dom, g.codom, compose(f.image, g.image))
end
morph_otimes(f::Morph{B}, g::Morph{B}) where {B} =
    Morph{B}(otimes(f.dom, g.dom), otimes(f.codom, g.codom), otimes(f.image, g.image))
morph_id(s::Stage{B}) where {B} = Morph{B}(s, s, id_image(s))
morph_braid(a::Stage{B}, b::Stage{B}) where {B} = Morph{B}(otimes(a, b), otimes(b, a), braid_image(a, b))

# the functions GATlab's `functor` calls on the images of the generators
compose(f::Morph{B}, g::Morph{B}) where {B} = morph_compose(f, g)
compose(f::Morph{B}, g::Morph{B}, h::Morph{B}, rest::Morph{B}...) where {B} = foldl(morph_compose, (f, g, h, rest...))
otimes(f::Morph{B}, g::Morph{B}) where {B} = morph_otimes(f, g)
otimes(f::Morph{B}, g::Morph{B}, h::Morph{B}, rest::Morph{B}...) where {B} = foldl(morph_otimes, (f, g, h, rest...))
id(s::Stage) = morph_id(s)
braid(a::Stage{B}, b::Stage{B}) where {B} = morph_braid(a, b)

"""
    LoweringCategory{B}()

The target of the lowering to backend `B` as a model of Catlab's theory of symmetric
monoidal categories, with objects `Stage{B}` (words of stages, one per deme) and
morphisms `Morph{B}`. Its operations are those the lowering uses; `Hom` checks a
morphism's stages.
"""
struct LoweringCategory{B} end

for B in (:kernel, :eidos, :moments)
    S, M, L = :(Stage{$(QuoteNode(B))}), :(Morph{$(QuoteNode(B))}), :(LoweringCategory{$(QuoteNode(B))})
    @eval @instance ThSymmetricMonoidalCategory{$S, $M} [model::$L] begin
        Ob(x::$S) = x
        Hom(f::$M, a::$S, b::$S) =
            f.dom.names == a.names && f.codom.names == b.names ? f :
            throw(ArgumentError("not a morphism $(join(a.names, "⊗")) → $(join(b.names, "⊗"))"))
        dom(f::$M) = f.dom
        codom(f::$M) = f.codom
        id(a::$S) = QTLSimCatlabExt.morph_id(a)
        compose(f::$M, g::$M) = QTLSimCatlabExt.morph_compose(f, g)
        otimes(a::$S, b::$S) = QTLSimCatlabExt.otimes(a, b)
        otimes(f::$M, g::$M) = QTLSimCatlabExt.morph_otimes(f, g)
        munit() = $S(Symbol[])
        braid(a::$S, b::$S) = QTLSimCatlabExt.morph_braid(a, b)
    end
end

"The names of the stages of an object of the life cycle (a word of stages)."
obnames(ob) = head(ob) === :generator ? Symbol[nameof(ob)] :
              head(ob) === :munit ? Symbol[] : reduce(vcat, (obnames(a) for a in args(ob)); init = Symbol[])
