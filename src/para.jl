# Parameters as an object (ROADMAP C3): a prediction seen as a function of some of the
# continuous parameters of a `Params`, and its sensitivities.

"""
    QTLSim.Para(f, base::Params, names)

A prediction `f(::Params) -> Real` seen as a function of the continuous parameters
`names` (a tuple or vector of `Params` field names, e.g. `(:PhenotypeVar, :Heritability)`),
the other fields held at their values in `base`. `P(x)` evaluates `f` at `base` with the
named fields set to the numbers `x` (in the order of `names`); `P()` at `base` itself;
`P.x0` holds the base values. Only `Float64` fields can be named: an integer, `Bool`,
`Symbol` or array field (`NoQTL`, `Ploidy`, `MoveType`, …) is an `ArgumentError`, as it
has no derivative; compare such parameters across separate `Para`s instead.

See [`sensitivities`](@ref).

```julia
P = QTLSim.Para(p -> house_of_gauss_vg(p), Params(NoQTL = 100, Ploidy = 2),
                (:MutationScalar, :PhenotypeVar))
P([0.002, 1.0])
sensitivities(P)          # elasticities ∂log V_G/∂log θ
```
"""
struct Para{F}
    f::F
    base::Params
    names::Vector{Symbol}
    x0::Vector{Float64}
    function Para(f::F, base::Params, names) where {F}
        ns = Symbol[Symbol(n) for n in names]
        isempty(ns) && throw(ArgumentError("Para needs at least one parameter name"))
        allunique(ns) || throw(ArgumentError("Para: a parameter is named twice in $(ns)"))
        for n in ns
            hasfield(Params, n) || throw(ArgumentError("Params has no field $n"))
            fieldtype(Params, n) === Float64 ||
                throw(ArgumentError("Para: $n is a $(fieldtype(Params, n)) field, not a continuous parameter; " *
                                    "only Float64 fields of Params can be varied"))
        end
        new{F}(f, base, ns, Float64[getfield(base, n) for n in ns])
    end
end

Base.show(io::IO, P::Para) = print(io, "Para(", P.f, ", ", join(P.names, ", "), " at ", P.x0, ")")

"`QTLSim.para_params(P, x)`: the `Params` at which `P(x)` evaluates its prediction."
para_params(P::Para, x::AbstractVector{<:Real}) =
    (length(x) == length(P.names) ||
         throw(ArgumentError("Para has $(length(P.names)) parameters, got $(length(x)) values"));
     Params(P.base; (n => Float64(v) for (n, v) in zip(P.names, x))...))

(P::Para)(x::AbstractVector{<:Real}) = P.f(para_params(P, x))
(P::Para)() = P.f(P.base)

"""
    sensitivities(P::Para; rel = 1e-3, log = true)

The derivatives of the prediction `P` with respect to each of its parameters, at its base
values, as a `NamedTuple` keyed by the parameter names.

With `log = true` (the default) they are elasticities `∂log y/∂log θ`: the relative
change of the prediction per relative change of a parameter. Each is a central
difference in `log θ` with step `h = rel`, `D(h) = (log y(θe^h) − log y(θe^{−h}))/(2h)`,
improved by one level of Richardson extrapolation, `(4D(h/2) − D(h))/3`, whose error is
`O(h⁴)` for a smooth prediction (four evaluations per parameter). A parameter must then be
positive. With `log = false` they are plain derivatives `∂y/∂θ`, the same scheme with step
`h = rel·|θ|` (`rel` if `θ = 0`).

These are finite differences, not automatic differentiation: `Params` holds `Float64`
fields, so a dual number cannot pass through it; forward-mode AD would need `Params` (and
the code below it) generic over the number type. Their accuracy is therefore limited by
the precision of `P`: a prediction computed to a relative tolerance `ε` gives derivatives
good to about `ε/rel`, so a root found by an iterative search must be found far more
precisely than `rel` (see [`emulated_half_survival_rate`](@ref)). An elasticity is `NaN`
when the prediction is not positive and finite at one of the points (`log = true`) or not
finite (`log = false`).
"""
function sensitivities(P::Para; rel::Real = 1e-3, log::Bool = true)
    rel > 0 || throw(ArgumentError("rel must be positive, got $rel"))
    out = Float64[]
    for (i, θ) in enumerate(P.x0)
        log && !(θ > 0) &&
            throw(ArgumentError("an elasticity needs a positive parameter; $(P.names[i]) = $θ (use log = false)"))
        at(t) = (x = copy(P.x0); x[i] = t; P(x))
        g(y) = log ? (y > 0 && isfinite(y) ? Base.log(y) : NaN) : (isfinite(y) ? Float64(y) : NaN)
        D(h) = log ? (g(at(θ * exp(h))) - g(at(θ * exp(-h)))) / 2h :
                     (s = θ == 0 ? h : h * abs(θ); (g(at(θ + s)) - g(at(θ - s))) / 2s)
        d1, d2 = D(rel), D(rel / 2)
        push!(out, (4d2 - d1) / 3)
    end
    NamedTuple{Tuple(P.names)}(Tuple(out))
end
