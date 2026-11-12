# Naturality as the validation harness (roadmap C2).
#
# The life cycle is one morphism `f` of a free symmetric monoidal category, lowered by
# functors to backends: the kernel (`lower(f, :kernel)`, a random map of `World`s), SLiM
# (`lower(f, :eidos)`), the moments of a Gaussian population (`lower(f, :moments)`), and
# closed forms of theory. A statistic σ maps the states (or the output tables) of one
# backend to those of another: the sample moments of a kernel population are a moments
# state, the kernel's records are the columns SLiM prints. The package's central empirical
# claims are naturality squares
#
#     σ ∘ K(f) = M(f) ∘ σ        (in distribution, or in expectation),
#
# exact for some steps and lifecycles (the neutral recursion is `deme_coalescence`'s, to
# rounding; the kernel and SLiM agree in law) and true up to a closure's error for others
# (Gaussian selection, the rare-alleles closure). Each validation script of `validate/`
# checks one such square, by hand, on some cells. Here a square is data, a
# `NaturalityCase`: the life cycle, the parameters, the two backends, the statistics, the
# criterion and its tolerance, written before the run, and the replicates. `check_naturality`
# runs it; `naturality_steps` checks the square of every step of one generation from
# matched states, so a disagreement is localised to the step whose square fails. A new
# step, model or closure gets its validation by adding a row (validate/naturality_cases.jl).
#
# The exact version is proved in Lean: `stochastic_lowerings_agree` (Lowering.lean, two
# lowerings that agree in distribution on every step agree on every life cycle) and
# `lower_conj` (Plasticity.lean, the same for lowerings that agree after a map of states).
# The lax version this harness tests is the same statement with an error: if the square of
# each step s commutes up to ε_s in a metric on the laws of the statistic, and each step's
# image is L_s-Lipschitz in that metric, the square of a composite of n steps commutes up
# to Σ_s ε_s Π_{s' after s} L_{s'} (a telescoping sum over the steps); this is not proved
# in Lean (see docs/src/naturality_design.md for the claims that would need lemmas).

"""
    QTLSim.Statistic(name, f; kernel = f, moments = f, slim = f, theory = f, scale = nothing)

A statistic σ of a naturality square: for each backend, a function from the output of
one run (the columns of [`run!`](@ref) or [`run_moments`](@ref), SLiM's printed table, or a
theory table) to a number, or for an exact comparison a vector. `scale`, a function of the
target backend's output, sets the width of the tolerance band (default: the absolute
value of the target's statistic). Ready-made statistics: [`QTLSim.window_mean`](@ref),
[`QTLSim.metapop_mean`](@ref), [`QTLSim.gmatrix_mean`](@ref), [`QTLSim.qst_stat`](@ref).
"""
struct Statistic
    name::String
    of::Dict{Symbol,Function}
    scale::Union{Nothing,Function}
    parts::Vector{Statistic}        # a derived statistic: `combine` of the parts' means
    combine::Union{Nothing,Function}
end

Statistic(name::AbstractString, f::Function; kernel = f, moments = f, slim = f, theory = f, scale = nothing) =
    Statistic(String(name), Dict{Symbol,Function}(:kernel => kernel, :moments => moments, :slim => slim, :theory => theory),
              scale, Statistic[], nothing)

"""
    QTLSim.Statistic(name, parts::Vector{Statistic}, combine; scale = nothing)

A derived statistic: `combine` applied to the values of `parts`, for a stochastic
backend to their means over replicates (its standard error by a bootstrap over
replicates), as Q_ST is computed from the mean within-deme variance and the mean variance
among deme means.
"""
Statistic(name::AbstractString, parts::Vector{Statistic}, combine::Function; scale = nothing) =
    Statistic(String(name), Dict{Symbol,Function}(), scale, parts, combine)

Base.show(io::IO, s::Statistic) = print(io, "Statistic(\"", s.name, "\")")

isderived(s::Statistic) = s.combine !== nothing

"The value of `s` on the output `out` of one run of `backend`."
function statvalue(s::Statistic, backend::Symbol, out)
    isderived(s) && return s.combine(map(q -> statvalue(q, backend, out), s.parts)...)
    f = get(s.of, backend, nothing)
    f === nothing && throw(ArgumentError("the statistic $(s.name) has no value for the backend :$backend"))
    f(out)
end

# the rows of a table whose generation is in `window`
inwindow(out, window) = findall(in(window), out.gen)

"""
    QTLSim.window_mean(column, window; name)

The mean of the column `column` (for example `:geno_var`, `:popsize`, `:nsurv`) over the
recorded generations in `window` (and over the demes): the same function of the kernel's
records, SLiM's table and the moments lowering's prediction.
"""
window_mean(col::Symbol, window; name = "mean $col over $window") =
    Statistic(name, out -> (ix = inwindow(out, window); isempty(ix) ? NaN : mean(getproperty(out, col)[ix])))

"""
    QTLSim.metapop_mean(what, window)

A statistic of the demes of a metapopulation over the census generations in `window`:
`what = :within`, the mean within-deme genetic variance (the kernel's sample variances;
the moments' `metapop.within`), or `:among`, the variance of the deme mean genotypes (the
kernel's sample variance of the demes' mean genotypes at each census, which needs the
`geno_mean` column the harness records; the moments' `metapop.var_among`). For several
traits `:within` is the mean within-deme G matrix.
"""
function metapop_mean(what::Symbol, window; name = "$(what == :within ? "V_W (within demes)" : "var of deme means") over $window")
    what in (:within, :among) || throw(ArgumentError("what is :within or :among; got :$what"))
    kernel = what === :within ?
        (out -> (ix = inwindow(out, window); haskey(out, :gmatrix) ? mean(out.gmatrix[ix]) : mean(out.geno_var[ix]))) :
        (out -> mean(var(out.geno_mean[findall(==(t), out.gen)]) for t in unique(out.gen[inwindow(out, window)])))
    moments = out -> (mp = out.metapop; ix = inwindow(mp, window);
                      what === :within ? mean(mp.within[ix]) : mean(mp.var_among[ix]))
    Statistic(name, kernel; moments, slim = kernel)
end

"""
    QTLSim.gmatrix_mean(i, j, window; within = false, scale = nothing)

Entry `(i, j)` of the mean G matrix over the recorded generations in `window` (the
kernel's sample covariances, `gmatrix`; with `within = true` and the moments lowering of
several demes, `metapop.within`). `scale = :geometric` measures the tolerance on an
off-diagonal entry relative to `√(G_ii G_jj)` of the prediction, as
validate/rare_alleles_traits.jl does.
"""
function gmatrix_mean(i::Integer, j::Integer, window; within::Bool = false, scale = nothing, name = "G$i$j over $window")
    G(out) = (ix = inwindow(out, window); mean(out.gmatrix[ix]))
    Gm(out) = within && haskey(out, :metapop) ? (mp = out.metapop; mean(mp.within[inwindow(mp, window)])) : G(out)
    sc = scale === :geometric ? (out -> (g = Gm(out); sqrt(g[i, i] * g[j, j]))) : scale
    Statistic(name, out -> G(out)[i, j]; moments = out -> Gm(out)[i, j], scale = sc)
end

"""
    QTLSim.qst_stat(window, k, N)

Q_ST = σ²/(σ² + 2V_W), σ² = (k − 1)/k·(among − V_W/N), from the means of
[`QTLSim.metapop_mean`](@ref)`(:within)` and `(:among)` over `window`, for `k` demes of
`N` (as validate/moments_demes.jl computes it).
"""
qst_stat(window, k, N; name = "Q_ST over $window") =
    Statistic(name, [metapop_mean(:within, window), metapop_mean(:among, window)],
              (w, a) -> (σ2 = (k - 1) / k * (a - w / N); σ2 / (σ2 + 2w)))

"""
    NaturalityCase(; name, source = "", cycle = nothing, params, from = :kernel, to = :moments,
                   criterion = :band, statistics, replicates = 20, horizon = nothing,
                   options = (;), expected = :pass, seed = 1, note = "")

A naturality square as data: one row of the table of validation cases
(validate/naturality_cases.jl). `cycle` is the life-cycle morphism (`nothing`: QTSim's
generation; for `k` demes a morphism on `Zygotes^⊗k`), `params` a `Params` or one per
deme, `from` and `to` two of the backends `:kernel`, `:slim` (stochastic: replicated),
`:moments` and `:theory` (deterministic: run once). `statistics` is a vector of
`statistic => tolerance` pairs ([`QTLSim.Statistic`](@ref)), and `criterion` how each is
compared, with its tolerance written before the run:

- `:band`: a stochastic backend against a deterministic one; FAILS only if the whole 95%
  CI of the replicates' mean lies outside the band `prediction·(1 ± tol)`, i.e.
  `|mean − pred| − 1.96 SE > tol·scale` (the rule of validate/moments_loci.jl and its
  successors);
- `:z`: the same pair, `|mean − pred|/SE < tol` (for exact predictions: neutral cases);
- `:exact`: two deterministic backends, the largest relative difference over the
  statistic's values (a vector: every generation) at most `tol`;
- `:ks`: two stochastic backends, a two-sample Kolmogorov–Smirnov test of the replicates'
  values, p ≥ `tol` (the rule of validate/slim/compare_*.jl).

`horizon` replaces `NoGenerations`; `options` passes keywords to the backends:
`options.kernel` and `options.moments` (for example `(loci = :multilocus, erosion = true)`,
`between = true`), `options.slim` (`slim`, the executable), `options.migration_rate`
(both), and `options.theory`, a function of the case giving the theory's table. `expected`
is the verdict of the script the case was ported from on the same cell (`:pass` or
`:fail`), so [`check_naturality`](@ref) reports whether it was reproduced.
"""
Base.@kwdef struct NaturalityCase
    name::String
    source::String = ""
    cycle::Any = nothing
    params::Any
    from::Symbol = :kernel
    to::Symbol = :moments
    criterion::Symbol = :band
    statistics::Vector{Pair{Statistic,Float64}}
    replicates::Int = 20
    horizon::Union{Nothing,Int} = nothing
    options::NamedTuple = (;)
    expected::Symbol = :pass
    seed::Int = 1
    note::String = ""
end

Base.show(io::IO, c::NaturalityCase) = print(io, "NaturalityCase(\"", c.name, "\": :", c.from, " → :", c.to, ", ",
                                              length(c.statistics), " statistic", length(c.statistics) == 1 ? "" : "s", ")")

const STOCHASTIC_BACKENDS = (:kernel, :slim)
const BACKENDS = (:kernel, :slim, :moments, :theory)

"""
    NaturalityResult

What [`check_naturality`](@ref) returns: the `case`, one row per statistic (`statistic`,
`tolerance`, the `from` backend's value and standard error, the `to` backend's, the
`measure` the criterion uses — relative error, z, largest relative difference or KS p —
and `pass`), whether every statistic passed (`pass`), whether that is the verdict the case
expected (`reproduced`), the replicates' values (`values`) and the time taken.
"""
struct NaturalityResult
    case::NaturalityCase
    rows::Vector{NamedTuple}
    pass::Bool
    reproduced::Bool
    values::Dict{Symbol,Any}
    seconds::Float64
end

sig(x::Real, d = 5) = isfinite(x) ? string(round(x; sigdigits = d)) : string(x)

function Base.show(io::IO, ::MIME"text/plain", r::NaturalityResult)
    c = r.case
    println(io, c.name, "  (:", c.from, " → :", c.to, ", ", c.criterion, ", ", c.replicates, " replicates",
            isempty(c.source) ? "" : "; from $(c.source)", ")")
    for x in r.rows
        m = c.criterion === :band ? string(sig(100x.measure, 3), "%") :
            c.criterion === :z ? string("z = ", sig(x.measure, 3)) :
            c.criterion === :ks ? string("p = ", sig(x.measure, 3)) : string("max rel. diff. ", sig(x.measure, 2))
        println(io, "  ", rpad(x.statistic, 36), " ", rpad(string(c.from), 8), rpad(fmtvalue(x.from, x.from_se), 24),
                rpad(string(c.to), 8), rpad(fmtvalue(x.to, x.to_se), 24), rpad(m, 24), " tol ", x.tolerance, "  ",
                x.pass ? "ok" : "FAIL")
    end
    print(io, "  verdict ", r.pass ? "PASS" : "FAIL", "; expected ", c.expected === :pass ? "PASS" : "FAIL", ": ",
          r.reproduced ? "reproduced" : "NOT reproduced", " (", round(Int, r.seconds), " s)")
end
Base.show(io::IO, r::NaturalityResult) = print(io, "NaturalityResult(\"", r.case.name, "\", ", r.pass ? "pass" : "fail", ")")

fmtvalue(x::Real, se) = isnan(se) ? sig(x) : string(sig(x), " ± ", sig(se, 2))
fmtvalue(x, se) = string("[", length(x), " values]")

# --- comparing the two sides of a square (no Catlab needed) -----------------------------

# the mean and standard error of a statistic over replicates' outputs; for a derived
# statistic, `combine` of the parts' means and a bootstrap standard error
function replicate_estimate(s::Statistic, backend, outs; nboot = 2000, rng = Xoshiro(20261111))
    if isderived(s)
        P = [statvalue(q, backend, o) for o in outs, q in s.parts]
        n = size(P, 1)
        v = s.combine((mean(P[:, j]) for j in axes(P, 2))...)
        boot = [(ix = rand(rng, 1:n, n); s.combine((mean(P[ix, j]) for j in axes(P, 2))...)) for _ in 1:nboot]
        return (value = v, se = std(boot), samples = Float64[])
    end
    x = [statvalue(s, backend, o) for o in outs]
    (value = mean(x), se = length(x) > 1 ? std(x) / sqrt(length(x)) : NaN, samples = x)
end

"""
    QTLSim.ks_pvalue(x, y)

The asymptotic p-value of the two-sample Kolmogorov–Smirnov test of samples `x` and `y`
(as HypothesisTests' `ApproximateTwoSampleKSTest`): `P(K > √(nm/(n + m))·D)` with `D`
the largest distance between the empirical distribution functions and `K` Kolmogorov's
distribution.
"""
function ks_pvalue(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n, m = length(x), length(y)
    xs, ys = sort(x), sort(y)
    D = 0.0
    for v in vcat(xs, ys)
        D = max(D, abs(searchsortedlast(xs, v) / n - searchsortedlast(ys, v) / m))
    end
    ccdf(Kolmogorov(), sqrt(n * m / (n + m)) * D)
end

function relerr_max(a, b)
    a, b = vec(collect(Float64, a)), vec(collect(Float64, b))
    length(a) == length(b) || return Inf
    maximum(abs(x - y) / max(abs(y), 1e-300) for (x, y) in zip(a, b); init = 0.0)
end

# one row of a result: the statistic `s` with tolerance `tol` compared under `criterion`
function compare_stat(criterion, s::Statistic, tol, from, to, fouts, touts)
    if criterion === :exact
        a, b = statvalue(s, from, only(fouts)), statvalue(s, to, only(touts))
        d = relerr_max(a, b)
        return (statistic = s.name, tolerance = tol, from = a, from_se = NaN, to = b, to_se = NaN, measure = d, pass = isfinite(d) && d <= tol)
    elseif criterion === :ks
        x = replicate_estimate(s, from, fouts); y = replicate_estimate(s, to, touts)
        pv = ks_pvalue(x.samples, y.samples)
        return (statistic = s.name, tolerance = tol, from = x.value, from_se = x.se, to = y.value, to_se = y.se, measure = pv,
                pass = pv >= tol)
    end
    x = replicate_estimate(s, from, fouts)
    pred = statvalue(s, to, only(touts))
    scale = s.scale === nothing ? abs(pred) : s.scale(only(touts))
    d = x.value - pred
    if criterion === :band
        return (statistic = s.name, tolerance = tol, from = x.value, from_se = x.se, to = pred, to_se = NaN,
                measure = d / scale, pass = isfinite(d) && !(abs(d) - 1.96x.se > tol * scale))
    elseif criterion === :z
        z = d / x.se
        return (statistic = s.name, tolerance = tol, from = x.value, from_se = x.se, to = pred, to_se = NaN, measure = z,
                pass = isfinite(z) && abs(z) < tol)
    end
    throw(ArgumentError("criterion is :band, :z, :exact or :ks; got :$criterion"))
end

function check_case(c::NaturalityCase)
    c.from in BACKENDS && c.to in BACKENDS || throw(ArgumentError("the backends are $(join(BACKENDS, ", ")); got :$(c.from), :$(c.to)"))
    c.expected in (:pass, :fail) || throw(ArgumentError("expected is :pass or :fail; got :$(c.expected)"))
    sf, st = c.from in STOCHASTIC_BACKENDS, c.to in STOCHASTIC_BACKENDS
    want = sf && st ? (:ks,) : sf ? (:band, :z) : st ? () : (:exact,)
    isempty(want) && throw(ArgumentError("put the stochastic backend first: from = :$(c.to), to = :$(c.from)"))
    c.criterion in want ||
        throw(ArgumentError("from :$(c.from) to :$(c.to) the criterion is $(join((":$w" for w in want), " or ")); got :$(c.criterion)"))
    c
end

"""
    check_naturality(case::NaturalityCase; seed = case.seed, ntasks = Threads.nthreads(), verbose = false) -> NaturalityResult
    check_naturality(cycle, p; from = :kernel, to = :moments, statistic, tolerance, criterion = :band,
                     replicates = 20, generations = nothing, seed = 1, options = (;), name = "") -> NaturalityResult

Check one naturality square: run the life cycle `case.cycle` with `case.params` in the
backends `case.from` (replicated `case.replicates` times when stochastic, seeds `seed + r`)
and `case.to`, evaluate each statistic on both sides and compare it under the case's
criterion and tolerance (see [`NaturalityCase`](@ref)). The second form builds the case
from its arguments, one `statistic` (a [`QTLSim.Statistic`](@ref)) or several with one
`tolerance` each. Replicates run on `ntasks` tasks. Requires Catlab.

The squares are those of the lowering functors: σ ∘ K(f) = M(f) ∘ σ for the statistic σ,
exact when the two lowerings agree in law (Lean `stochastic_lowerings_agree`, and
`lower_conj` after a map of states), and up to a closure's error otherwise; see
[`naturality_steps`](@ref) for the square of every step of a generation, and
validate/naturality_cases.jl for the table of cases ported from the validation scripts.
"""
function check_naturality end

"""
    naturality_steps(p; cycle = lifecycle().generation, loci = :multilocus, erosion = false, replicates = 200,
                     world = nothing, burnin = 1000, linkage_equilibrium = true, tolerance = 0.05, seed = 1)

The naturality square of each step of one generation, from matched states, for one deme
of one trait: a kernel population `w` (`world`, or `burnin` generations of `p` from zero
genomes, its alleles then shuffled among chromosomes locus by locus when
`linkage_equilibrium`, so that it has no LD and is in Hardy–Weinberg proportions) and the
moments state σ(w) with its mean genotype and genetic variance (`loci`, `erosion` as for
[`run_moments`](@ref); the moments' initial haplotypes are in linkage equilibrium). The
steps of `cycle` (a composite of steps from Zygotes to Zygotes) are applied in order, to
`replicates` copies of `w` with their own random numbers and to σ(w), and after each step
the observables of its codomain stage are compared: the mean genotype (relative to the
predicted √V_G), the genetic variance and the number of individuals (survivors, offspring).
A square whose comparison FAILS (whole 95% CI outside `tolerance`) after a step whose
predecessors all passed localises the disagreement to that step. Returns a vector of
rows `(step, stage, observable, kernel, se, moments, measure, pass)`; `first_failure`
of the result names the step. Requires Catlab.

For example (100 loci, N = 1000, ω = 4; docs/src/naturality_design.md) the squares of
develop and select hold for every closure; that of reproduce fails by 18% for the aggregate
recursion (`loci = :aggregate`), which inherits all of selection's reduction of the
variance, and is within 4% for the multilocus model (`loci = :multilocus`). Under strong
selection (ω = 1, V_G/s ≈ 0.7) the multilocus and rare-alleles reproduce squares miss by
about 10%: the closure's error sits in reproduce, not in select.
"""
function naturality_steps end

"""
    QTLSim.first_failure(rows) -> Union{Nothing,NamedTuple}

The first row of [`naturality_steps`](@ref)' result that fails: the step whose square
does not commute, given that those before it do.
"""
first_failure(rows) = (i = findfirst(r -> !r.pass, rows); i === nothing ? nothing : rows[i])
