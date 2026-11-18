# The eco-evolutionary model of Xu and Shen (2026): a re-implementation from the equations of
#
#   Xu, K. and Shen, H. (2026). Eco-evolutionary feedbacks generate bistability in population
#   persistence under gradual environmental change. bioRxiv preprint,
#   doi:10.64898/2026.08.05.743160 (version 1, posted 9 August 2026; not peer reviewed).
#
# The preprint gives equations and figures but, as far as we know, no code; everything here
# is written from its equations (numbers in brackets are the preprint's). It is a theory
# module off to the side of the simulation: nothing else in the package uses it, and no
# default changes.
#
# The model (a haploid population with overlapping generations, continuous time, a fully
# heritable Gaussian trait with mean μ and variance V, an optimum θ(t) = θ₀ + vt, lag
# Δ = θ − μ, logistic growth; eqs. 1–5, S1):
#
#     dΔ/dt = v − VΔ/ω²                                  (4)
#     dV/dt = V_M − V²/(n_e ω²) − V/N                    (5)
#     dN/dt = N [r − (V + Δ²)/(2ω²) − N/K]               (1, 2)
#
# with r the maximum growth rate, ω² the width of the Gaussian fitness function, V_M the
# mutational input of variance, n_e the effective number of recombining factors (1 for an
# asexual population) and K the strength of density regulation (the equilibrium size is
# K × the mean growth rate, so the carrying capacity without load is rK). The standard
# parameters of the preprint's figures are r = 0.2, ω² = 50, V_M = 0.1, K = 1000, n_e = 1.
#
# The feedback (lag → fewer individuals → more drift → less variance → more lag) was
# described by Bürger and Lynch (1995, Evolution 49:151–163); the preprint's eq. 6 defines
# the critical rate as a maximum over steady states with drift-reduced variance, and its
# Supplement S2 shows that below it a viable equilibrium and extinction are both locally
# stable. QTLSim's `drift_critical_rate` (src/critical.jl) is a later discrete-generation
# counterpart of eq. 6 matched to QTSim's life cycle; validate/xu_shen.jl compares the two.
#
# Checks: test/xu_shen.jl; validate/symbolic/check_xu_shen.jl (the algebra of S1–S2, with
# notes on where the printed equations and the derivation differ); validate/xu_shen.jl
# (the preprint's figures).

"""
    XuShen(; r = 0.2, ω2 = 50, VM = 0.1, K = 1000, ne = 1)
    XuShen(m; kw...)

The parameters of Xu and Shen's (2026, bioRxiv preprint, doi:10.64898/2026.08.05.743160)
model of a population following a moving optimum with feedback among the lag `Δ`, the
genetic variance `V` and the population size `N` (their eqs. 1–5):

    dΔ/dt = v − VΔ/ω²,   dV/dt = V_M − V²/(n_e ω²) − V/N,   dN/dt = N[r − (V + Δ²)/(2ω²) − N/K],

with `r` the maximum growth rate, `ω2` = ω² the width of the fitness function, `VM` the
mutational variance per unit time, `K` the strength of density regulation and `ne` the
effective number of recombining factors. The defaults are the standard parameters of the
preprint's figures. The second form copies `m` with some fields replaced. The rate `v` of
the optimum is an argument of the functions, not a parameter.

A re-implementation from the preprint's equations (it is not peer reviewed, and we know of
no published code). See [`xu_shen_critical_rate`](@ref), [`xu_shen_equilibria`](@ref),
[`xu_shen_stability`](@ref), [`xu_shen_trajectory`](@ref), [`xu_shen_boundary`](@ref).
"""
struct XuShen
    r::Float64
    ω2::Float64
    VM::Float64
    K::Float64
    ne::Float64
    function XuShen(r::Real, ω2::Real, VM::Real, K::Real, ne::Real)
        (r > 0 && ω2 > 0 && VM > 0 && K > 0 && ne > 0) ||
            throw(ArgumentError("XuShen needs r, ω2, VM, K, ne > 0; got r = $r, ω2 = $ω2, VM = $VM, K = $K, ne = $ne"))
        (isfinite(r) && isfinite(ω2) && isfinite(VM) && isfinite(ne)) ||
            throw(ArgumentError("XuShen: only K may be infinite"))
        new(r, ω2, VM, K, ne)
    end
end
XuShen(; r::Real = 0.2, ω2::Real = 50.0, VM::Real = 0.1, K::Real = 1000.0, ne::Real = 1.0) = XuShen(r, ω2, VM, K, ne)
XuShen(m::XuShen; r::Real = m.r, ω2::Real = m.ω2, VM::Real = m.VM, K::Real = m.K, ne::Real = m.ne) = XuShen(r, ω2, VM, K, ne)

# --- the vector field and the functions of Supplement S1 ---------------------------------

"""
    xu_shen_field(m, Δ, V, N; v)

The right-hand sides `(dΔ/dt, dV/dt, dN/dt)` of Xu and Shen's (2026) system at lag `Δ`,
variance `V` and size `N`, with the optimum moving at rate `v` (eqs. 4, 5 and 2; S1).
"""
xu_shen_field(m::XuShen, Δ::Real, V::Real, N::Real; v::Real) =
    (v - V * Δ / m.ω2, m.VM - V^2 / (m.ne * m.ω2) - V / N, N * (m.r - (V + Δ^2) / (2m.ω2) - N / m.K))

"""
    xu_shen_vmax(m)

`V_c^M = ω√(n_e V_M)`, the equilibrium variance of an infinitely large population (below
eq. 6): the upper end of the range of variances a finite population can hold.
"""
xu_shen_vmax(m::XuShen) = sqrt(m.ne * m.ω2 * m.VM)

"""
    xu_shen_balance_variance(m, N)

The variance at mutation–selection–drift balance in a population held at size `N`, the
positive root of `V_M − V²/(n_e ω²) − V/N = 0` (eq. 5): the initial variance of the
preprint's first scenario (Figure 4a). It tends to [`xu_shen_vmax`](@ref) as `N → ∞` and
to `V_M N` as `N → 0`.
"""
function xu_shen_balance_variance(m::XuShen, N::Real)
    N > 0 || throw(ArgumentError("N must be positive, got $N"))
    c = m.ne * m.ω2
    2m.VM / (1 / N + sqrt(1 / N^2 + 4m.VM / c))
end

# N₁(V), the size needed to hold variance V (S3), its derivative (S6a), f(V) and h(V) (S1)
_xs_N1(m::XuShen, V) = (c = m.ne * m.ω2; V * c / (m.VM * c - V^2))
_xs_dN1(m::XuShen, V) = (c = m.ne * m.ω2; c * (m.VM * c + V^2) / (m.VM * c - V^2)^2)
_xs_f(m::XuShen, V) = m.r - V / (2m.ω2) - _xs_N1(m, V) / m.K
_xs_h(m::XuShen, V) = 2V^2 / m.ω2 * _xs_f(m, V)
# h′(V) = (2V/ω²)·φ(V) (S5), φ decreasing from 2r to −∞ on (0, V_c^M)
_xs_phi(m::XuShen, V) = 2 * _xs_f(m, V) - V * _xs_dN1(m, V) / m.K - V / (2m.ω2)

# the root of g in (a, b), with g(a) > 0 > g(b) assumed (g need not be evaluable at a or b)
function _xs_bisect(g, a::Float64, b::Float64)
    for _ in 1:200
        mid = (a + b) / 2
        (mid <= a || mid >= b) && break
        g(mid) > 0 ? (a = mid) : (b = mid)
    end
    (a + b) / 2
end

"""
    xu_shen_critical_rate(m)

The critical rate of Xu and Shen (2026), eq. 6: above it no viable equilibrium exists,

    v_c² = max over V ∈ (0, V_c^M) of h(V),   h(V) = (2V²/ω²)[r − V/(2ω²) − N₁(V)/K],

with `N₁(V) = V n_e ω²/(V_M n_e ω² − V²)` the size needed to hold the variance `V` and
`V_c^M = ω√(n_e V_M)`. By S5–S6 `h′(V) = (2V/ω²)φ(V)` with `φ` decreasing from `2r` to
`−∞`, so `h` has one maximum, found here as the root of `φ` by bisection (no grid).

Returns `(; v, V, Δ, N)`: the critical rate and the variance, lag `vω²/V` and size `N₁(V)`
of the equilibrium at it. With `K = Inf` there is no drift and the result is the classical
critical rate at `V = min(V_c^M, 4rω²/3)`, `v = (V/ω)√(2r − V/ω²)` (then `N = Inf`).
"""
function xu_shen_critical_rate(m::XuShen)
    Vm = xu_shen_vmax(m)
    V = if isinf(m.K)
        min(Vm, 4m.r * m.ω2 / 3)
    else
        _xs_bisect(x -> _xs_phi(m, x), 0.0, Vm)
    end
    h = isinf(m.K) ? 2V^2 / m.ω2 * (m.r - V / (2m.ω2)) : _xs_h(m, V)
    v = sqrt(max(h, 0.0))
    (; v, V, Δ = v * m.ω2 / V, N = isinf(m.K) ? Inf : _xs_N1(m, V))
end

function _xs_finiteK(m::XuShen)
    isfinite(m.K) || throw(ArgumentError("this needs a finite K (with K = Inf the population size has no equilibrium)"))
    nothing
end

"""
    xu_shen_equilibria(m, v)

The equilibria with `N > 0` of Xu and Shen's (2026) system at rate `v ≥ 0`, in increasing
order of `V`: the roots in `(0, V_c^M)` of their eq. S4 (eq. 8), `v² = h(V)`, each with
`Δ* = vω²/V*` and `N* = V* n_e ω²/(V_M n_e ω² − V*²)` (eq. 7, S2). A vector of
`(; Δ, V, N, stable)`:

- none for `v` above [`xu_shen_critical_rate`](@ref);
- two for `0 < v < v_c`: the one with the smaller variance is a saddle (`det J > 0`, S15),
  the one with the larger variance, `V₂` in S2, is locally stable (`stable = true`); the
  saddle's stable manifold is the basin boundary;
- one for `v = 0` (the stable one, `Δ = 0`; the other root is `V = 0`, i.e. `N = 0`).

`stable` is the Routh–Hurwitz verdict of [`xu_shen_stability`](@ref).
"""
function xu_shen_equilibria(m::XuShen, v::Real)
    _xs_finiteK(m)
    v >= 0 || throw(ArgumentError("the rate v must be ≥ 0, got $v"))
    T = @NamedTuple{Δ::Float64, V::Float64, N::Float64, stable::Bool}
    out = T[]
    c = xu_shen_critical_rate(m)
    v2 = float(v)^2
    hc = _xs_h(m, c.V)
    v2 > hc && return out
    V0 = _xs_bisect(x -> _xs_f(m, x), 0.0, xu_shen_vmax(m))          # the zero of f: h > 0 on (0, V0)
    eq(V) = (N = _xs_N1(m, V); Δ = v * m.ω2 / V; (; Δ, V, N, stable = xu_shen_stability(m, Δ, V, N).stable))
    if v > 0
        Vlo = _xs_bisect(x -> v2 - _xs_h(m, x), 0.0, c.V)             # h increasing on (0, V_c)
        push!(out, eq(Vlo))
    end
    Vhi = v > 0 ? _xs_bisect(x -> _xs_h(m, x) - v2, c.V, V0) : V0     # h decreasing on (V_c, V0)
    (isempty(out) || Vhi > out[1].V) && push!(out, eq(Vhi))
    out
end

"""
    xu_shen_equilibrium(m, v)

The locally stable viable equilibrium `(; Δ, V, N)` of Xu and Shen's (2026) system at rate
`v` (eqs. 7–8; `V₂` of their S2), or `nothing` when there is none (`v` at or above
[`xu_shen_critical_rate`](@ref)). See [`xu_shen_equilibria`](@ref).
"""
function xu_shen_equilibrium(m::XuShen, v::Real)
    for e in xu_shen_equilibria(m, v)
        e.stable && return (; e.Δ, e.V, e.N)
    end
    nothing
end

# --- the Jacobian and local stability (Supplement S2) ------------------------------------

"""
    xu_shen_jacobian(m, Δ, V, N)

The Jacobian of [`xu_shen_field`](@ref) with respect to `(Δ, V, N)` at any state (it does
not depend on the rate `v`). At an equilibrium it is Xu and Shen's (2026) matrix S8:

    J = [ −V/ω²     −Δ/ω²              0
           0        −2V/(n_e ω²) − 1/N  V/N²
          −NΔ/ω²    −N/(2ω²)           r − (V + Δ²)/(2ω²) − 2N/K ],

where `Δ/ω² = v/V` and the last entry is `−N/K`.
"""
function xu_shen_jacobian(m::XuShen, Δ::Real, V::Real, N::Real)
    [-V/m.ω2 -Δ/m.ω2 0.0
     0.0 (-2V / (m.ne * m.ω2)-1 / N) V/N^2
     -N*Δ/m.ω2 -N/(2m.ω2) (m.r - (V + Δ^2) / (2m.ω2)-2N / m.K)]
end

"""
    xu_shen_stability(m, Δ, V, N)

Local stability of a state of Xu and Shen's (2026) system from the Jacobian
([`xu_shen_jacobian`](@ref)): the coefficients of the characteristic polynomial
`λ³ + a₁λ² + a₂λ + a₃` (S16), the Routh–Hurwitz conditions `a₁ > 0`, `a₃ > 0`,
`a₁a₂ − a₃ > 0` (S9), and the eigenvalues. Returns `(; a1, a2, a3, hurwitz, stable,
eigenvalues)` with `hurwitz = a₁a₂ − a₃` and `stable` the Routh–Hurwitz verdict (meaningful
at an equilibrium).

At an equilibrium `a₃ = −det J = −h′(V*)/(2N*)` (S15), so the equilibrium on the falling
side of `h` is the candidate. `a₁ > 0` because the diagonal is negative. For the third
condition, with the signs `J₁₁, J₂₂, J₃₃, J₁₂, J₃₁, J₃₂ < 0 < J₂₃` and `J₁₃ = J₂₁ = 0`,

    a₁a₂ − a₃ = J₂₃J₃₁J₁₂ + J₂₂J₂₃J₃₂ + J₃₃J₂₃J₃₂ − (J₁₁ + J₂₂)(J₁₁J₂₂ + J₃₃(J₁₁ + J₂₂ + J₃₃)),

every term of which is positive, so it holds at every equilibrium with `v > 0`. (The
preprint's S17 carries the term `J₂₃(J₁₁J₃₂ − J₃₁J₁₂)` with the opposite sign to the one
that follows from S16; the conclusion of S2 is unaffected. See
`validate/symbolic/check_xu_shen.jl`.)
"""
function xu_shen_stability(m::XuShen, Δ::Real, V::Real, N::Real)
    J = xu_shen_jacobian(m, Δ, V, N)
    a1 = -(J[1, 1] + J[2, 2] + J[3, 3])
    a2 = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1] + J[1, 1] * J[3, 3] - J[1, 3] * J[3, 1] + J[2, 2] * J[3, 3] - J[2, 3] * J[3, 2]
    a3 = -det(J)
    hurwitz = a1 * a2 - a3
    (; a1, a2, a3, hurwitz, stable = a1 > 0 && a3 > 0 && hurwitz > 0, eigenvalues = eigvals(J))
end

# --- trajectories -------------------------------------------------------------------------
#
# The system is stiff when N is small (the variance relaxes at rate 1/N), so the integrator
# is a linearly implicit one with step-size control: the Rosenbrock pair of order 2(3) of
# Shampine and Reichelt (1997, SIAM J. Sci. Comput. 18:1–22; "ode23s"), on y = (Δ, V, ln N)
# with the exact Jacobian. ln N keeps N positive and makes the tolerance on N relative.

@inline function _xs_rhs(m::XuShen, y::NTuple{3,Float64}, v::Float64)
    Δ, V, ℓ = y
    (v - V * Δ / m.ω2, m.VM - V^2 / (m.ne * m.ω2) - V * exp(-ℓ), m.r - (V + Δ^2) / (2m.ω2) - exp(ℓ) / m.K)
end

# W = I − hd·J in (Δ, V, ln N), inverted by the adjugate; row-major 9-tuple
@inline function _xs_winv(m::XuShen, y::NTuple{3,Float64}, hd::Float64)
    Δ, V, ℓ = y
    e = exp(-ℓ)
    a11, a12, a13 = 1 + hd * V / m.ω2, hd * Δ / m.ω2, 0.0
    a21, a22, a23 = 0.0, 1 + hd * (2V / (m.ne * m.ω2) + e), -hd * V * e
    a31, a32, a33 = hd * Δ / m.ω2, hd / (2m.ω2), 1 + hd * exp(ℓ) / m.K
    c11, c12, c13 = a22 * a33 - a23 * a32, a13 * a32 - a12 * a33, a12 * a23 - a13 * a22
    c21, c22, c23 = a23 * a31 - a21 * a33, a11 * a33 - a13 * a31, a13 * a21 - a11 * a23
    c31, c32, c33 = a21 * a32 - a22 * a31, a12 * a31 - a11 * a32, a11 * a22 - a12 * a21
    d = a11 * c11 + a12 * c21 + a13 * c31
    (c11 / d, c12 / d, c13 / d, c21 / d, c22 / d, c23 / d, c31 / d, c32 / d, c33 / d)
end
@inline _xs_mul(A::NTuple{9,Float64}, x::NTuple{3,Float64}) =
    (A[1] * x[1] + A[2] * x[2] + A[3] * x[3], A[4] * x[1] + A[5] * x[2] + A[6] * x[3], A[7] * x[1] + A[8] * x[2] + A[9] * x[3])

const _XS_D = 1 / (2 + sqrt(2.0))
const _XS_E32 = 6 + sqrt(2.0)

# one phase at a constant rate v, from y at time t for at most a time T. Stops when
# ln N < lnext (:extinct), when within `radius` of `target` (:converged; relative distances,
# target = nothing for no test), or at T (:time). `rec`, if not nothing, records each accepted step.
function _xs_run(m::XuShen, y::NTuple{3,Float64}, t::Float64, v::Float64, T::Float64; rtol::Float64, lnext::Float64,
                 maxstep::Float64, target, radius::Float64, rec)
    scΔ, scV = 1e-3 * sqrt(m.ω2), 1e-6 * xu_shen_vmax(m)
    tend = t + T
    h = min(1e-3, T, maxstep)
    F0 = _xs_rhs(m, y, v)
    nrej = 0
    while t < tend
        y[3] < lnext && return (y, t, :extinct)
        if target !== nothing
            dist = max(abs(y[1] - target[1]) / max(abs(target[1]), 0.01 * sqrt(m.ω2)), abs(y[2] - target[2]) / target[2], abs(y[3] - target[3]))
            dist < radius && return (y, t, :converged)
        end
        h = min(h, tend - t, maxstep)
        Wi = _xs_winv(m, y, h * _XS_D)
        k1 = _xs_mul(Wi, F0)
        F1 = _xs_rhs(m, y .+ (0.5h) .* k1, v)
        k2 = _xs_mul(Wi, F1 .- k1) .+ k1
        yn = y .+ h .* k2
        F2 = _xs_rhs(m, yn, v)
        k3 = _xs_mul(Wi, F2 .- _XS_E32 .* (k2 .- F1) .- 2.0 .* (k1 .- F0))
        er = (h / 6) .* (k1 .- 2.0 .* k2 .+ k3)
        err = max(abs(er[1]) / (rtol * (scΔ + max(abs(y[1]), abs(yn[1])))), abs(er[2]) / (rtol * (scV + max(abs(y[2]), abs(yn[2])))),
                  abs(er[3]) / rtol)
        if isfinite(err) && err <= 1 && yn[2] > 0
            t += h
            y, F0 = yn, F2
            rec === nothing || rec(t, y)
            h *= min(5.0, max(0.2, 0.8 * err^(-1 / 3)))
            nrej = 0
        else
            h *= isfinite(err) ? max(0.1, 0.8 * err^(-1 / 3)) : 0.1
            (nrej += 1) > 200 && error("xu_shen: the step size collapsed at t = $t, state (Δ, V, N) = ($(y[1]), $(y[2]), $(exp(y[3])))")
        end
    end
    y[3] < lnext && return (y, t, :extinct)
    (y, t, :time)
end

"The default extinction threshold of the Xu–Shen trajectories: a population is extinct once `N` falls below it."
const XU_SHEN_NEXT = 1e-3

"""
    xu_shen_trajectory(m, Δ0, V0, N0; v, phases = (), tmax = 1e5, Next = QTLSim.XU_SHEN_NEXT,
                       rtol = 1e-7, maxstep = 0.5/r, save = true)

A trajectory of Xu and Shen's (2026) system (eqs. 2, 4, 5) from lag `Δ0`, variance `V0` and
size `N0`, with the optimum moving at rate `v`. `phases` is a list of `(rate, duration)`
pairs that are run first, in order, before the open-ended final phase at `v` (so a
temporary shift to `v1` for a time `T` is `phases = [(v1, T)]`).

Returns `(; t, Δ, V, N, fate, equilibrium)`: the accepted steps (only the end point with
`save = false`), the stable equilibrium at `v` (or `nothing`) and

- `fate = :extinct` if `N` falls below `Next` (in any phase);
- `fate = :persists` if, in the final phase, the state comes within a relative distance
  `1e-3` of the stable viable equilibrium (less when the saddle is close: 1% of the
  distance between the two equilibria), from where it converges;
- `fate = :undecided` if neither happens by time `tmax` after the phases.

**The extinction rule.** In the equations `N` never reaches 0, and a small `N` can grow
back while the lag is below `ω√(2r)`; extinction is the attracting state `N → 0`,
`V/N → V_M`, `Δ → ∞` (S1). The rule is therefore a threshold far below one individual:
the default `Next = 1e-3` leaves the boundary between the two fates where the equations
put it (through the saddle equilibrium, [`xu_shen_equilibria`](@ref)); `Next = 1` ("fewer
than one individual") is the rule to use for a population of whole individuals, and moves
the boundary only where it passes through sizes of a few individuals
(`validate/xu_shen.log` gives the sensitivity).

Integrated by an adaptive Rosenbrock method of order 2 (Shampine and Reichelt's `ode23s`)
on `(Δ, V, ln N)` with relative tolerance `rtol`, since the variance equation is stiff at
small `N`. `maxstep` caps the step: the method damps every mode when its steps are long, the
unstable one of the saddle included, so steps are kept below half the time scale `1/r` of
population growth (otherwise a trajectory started next to the boundary can sit on the
saddle).
"""
function xu_shen_trajectory(m::XuShen, Δ0::Real, V0::Real, N0::Real; v::Real, phases = (), tmax::Real = 1e5,
                            Next::Real = XU_SHEN_NEXT, rtol::Real = 1e-7, maxstep::Real = 0.5 / m.r, save::Bool = true)
    _xs_finiteK(m)
    (V0 > 0 && N0 > 0) || throw(ArgumentError("V0 and N0 must be positive, got V0 = $V0, N0 = $N0"))
    (v >= 0 && Next > 0 && rtol > 0 && tmax > 0) || throw(ArgumentError("need v ≥ 0 and Next, rtol, tmax > 0"))
    ts, Δs, Vs, Ns = Float64[], Float64[], Float64[], Float64[]
    rec = save ? (t, y) -> (push!(ts, t); push!(Δs, y[1]); push!(Vs, y[2]); push!(Ns, exp(y[3])); nothing) : nothing
    y = (float(Δ0), float(V0), log(float(N0)))
    t = 0.0
    save && rec(t, y)
    lnext = log(float(Next))
    status = :time
    for (v1, T) in phases
        (v1 >= 0 && T >= 0) || throw(ArgumentError("phases are (rate ≥ 0, duration ≥ 0) pairs, got ($v1, $T)"))
        T == 0 && continue
        y, t, status = _xs_run(m, y, t, float(v1), float(T); rtol = float(rtol), lnext, maxstep = float(maxstep),
                               target = nothing, radius = 0.0, rec)
        status === :extinct && break
    end
    eqs = xu_shen_equilibria(m, v)
    stable = findfirst(e -> e.stable, eqs)
    equilibrium = stable === nothing ? nothing : (; eqs[stable].Δ, eqs[stable].V, eqs[stable].N)
    if status !== :extinct
        target, radius = nothing, 0.0
        if equilibrium !== nothing
            target = (equilibrium.Δ, equilibrium.V, log(equilibrium.N))
            radius = 1e-3
            if length(eqs) == 2
                s = eqs[1]
                gap = max(abs(s.Δ - equilibrium.Δ) / max(equilibrium.Δ, 0.01 * sqrt(m.ω2)), abs(s.V - equilibrium.V) / equilibrium.V,
                          abs(log(s.N / equilibrium.N)))
                radius = min(radius, 0.01 * gap)
            end
        end
        y, t, status = _xs_run(m, y, t, float(v), float(tmax); rtol = float(rtol), lnext, maxstep = float(maxstep), target, radius, rec)
    end
    fate = status === :extinct ? :extinct : status === :converged ? :persists : :undecided
    save || (push!(ts, t); push!(Δs, y[1]); push!(Vs, y[2]); push!(Ns, exp(y[3])))
    (; t = ts, Δ = Δs, V = Vs, N = Ns, fate, equilibrium)
end

"""
    xu_shen_persists(m, Δ0, V0, N0; v, kw...)

Whether the population of Xu and Shen's (2026) model that starts at lag `Δ0`, variance `V0`
and size `N0` persists at rate `v`, i.e. reaches the stable viable equilibrium rather than
extinction (their Figure 3): `fate === :persists` of [`xu_shen_trajectory`](@ref), whose
keywords (`phases`, `Next`, `rtol`, `tmax`) it takes. Throws if the fate is undecided at
`tmax`.
"""
function xu_shen_persists(m::XuShen, Δ0::Real, V0::Real, N0::Real; v::Real, kw...)
    fate = xu_shen_trajectory(m, Δ0, V0, N0; v, save = false, kw...).fate
    fate === :undecided && error("xu_shen_persists: undecided at tmax from (Δ0, V0, N0) = ($Δ0, $V0, $N0) at v = $v; raise tmax")
    fate === :persists
end

# persistence inside the searches: a start still undecided at tmax sits on the boundary to
# rounding (the trajectory rests on the saddle) and counts as not persisting
_xs_ok(m::XuShen, Δ0, V0, N0; v, kw...) = xu_shen_trajectory(m, Δ0, V0, N0; v, save = false, kw...).fate === :persists

# the smallest x in [lo, hi] (on a log scale) at which `ok` holds, assuming it holds from some
# point upward: `lo` if it holds there, Inf if it does not hold at hi
function _xs_threshold(ok, lo::Float64, hi::Float64; steps::Int = 50)
    ok(lo) && return lo
    ok(hi) || return Inf
    a, b = log(lo), log(hi)
    for _ in 1:steps
        mid = (a + b) / 2
        ok(exp(mid)) ? (b = mid) : (a = mid)
    end
    exp((a + b) / 2)
end

"""
    xu_shen_boundary(m, N0; Δ0, v, Vmin = 1e-6·V_c^M, Vmax = 4·V_c^M, kw...)

The basin boundary of Xu and Shen's (2026) Figure 3 at initial size `N0`: the smallest
initial variance `V0` from which a population with initial lag `Δ0` persists at rate `v`,
by bisection on `log V0` between `Vmin` and `Vmax` (it assumes that more initial variance
never hurts within that range, which `validate/xu_shen.jl` checks on grids). Returns `Vmin`
if the population persists from there, `Inf` if it does not from `Vmax`, and `Inf` for `v`
at or above the critical rate. Keywords of [`xu_shen_trajectory`](@ref) (`Next`, `rtol`,
`tmax`) pass through. In this and the other searches a start still undecided at `tmax`
(on the boundary to rounding, where the trajectory rests on the saddle) counts as not
persisting.
"""
function xu_shen_boundary(m::XuShen, N0::Real; Δ0::Real, v::Real, Vmin::Real = 1e-6 * xu_shen_vmax(m),
                          Vmax::Real = 4 * xu_shen_vmax(m), kw...)
    xu_shen_equilibrium(m, v) === nothing && return Inf
    _xs_threshold(V0 -> _xs_ok(m, Δ0, V0, N0; v, kw...), float(Vmin), float(Vmax))
end

"""
    xu_shen_critical_size(m; Δ0, v, V0 = :balance, Nmin = 1.0, Nmax = 1e6, kw...)

The smallest initial size `N0` from which a population of Xu and Shen's (2026) model with
initial lag `Δ0` persists at rate `v`, by bisection on `log N0` between `Nmin` and `Nmax`.
With `V0 = :balance` the initial variance is at mutation–selection–drift balance for that
size ([`xu_shen_balance_variance`](@ref)), the preprint's first scenario (Figure 4a); a
number fixes it. Returns `Nmin` if the population persists from there (the default 1 is the
lowest value on Figure 4a's axis), and `Inf` if it does not from `Nmax` or if `v` is at or
above the critical rate.
"""
function xu_shen_critical_size(m::XuShen; Δ0::Real, v::Real, V0 = :balance, Nmin::Real = 1.0, Nmax::Real = 1e6, kw...)
    V0 === :balance || V0 isa Real || throw(ArgumentError("V0 is :balance or a number, got $V0"))
    xu_shen_equilibrium(m, v) === nothing && return Inf
    V(N) = V0 === :balance ? xu_shen_balance_variance(m, N) : float(V0)
    _xs_threshold(N0 -> _xs_ok(m, Δ0, V(N0), N0; v, kw...), float(Nmin), float(Nmax))
end

"""
    xu_shen_max_rate(m; v0, δ = 0, steps = 40, kw...)

The preprint's second scenario (Figure 4b): a population at the stable equilibrium of rate
`v0` meets an abrupt shift `δ` of the optimum (its lag becomes `Δ* + δ`) and a new rate
`v1`. Returns the largest `v1` at which it persists, by bisection between 0 and the
critical rate `v_c`: `v_c` itself if it persists just below it (at `(1 − 1e-6)v_c`), 0 if
it does not persist even with the optimum at rest afterwards.
"""
function xu_shen_max_rate(m::XuShen; v0::Real, δ::Real = 0.0, steps::Integer = 40, kw...)
    e = xu_shen_equilibrium(m, v0)
    e === nothing && throw(ArgumentError("no viable equilibrium at v0 = $v0 (critical rate $(xu_shen_critical_rate(m).v))"))
    vc = xu_shen_critical_rate(m).v
    ok(v1) = _xs_ok(m, e.Δ + δ, e.V, e.N; v = v1, kw...)
    ok((1 - 1e-6) * vc) && return vc
    ok(0.0) || return 0.0
    a, b = 0.0, vc
    for _ in 1:steps
        mid = (a + b) / 2
        ok(mid) ? (a = mid) : (b = mid)
    end
    (a + b) / 2
end

"""
    xu_shen_max_duration(m; v0, v1, Tmax = 1e5, steps = 40, kw...)

The preprint's third scenario (Figure 4d): a population at the stable equilibrium of rate
`v0` meets the rate `v1` for a time `T`, after which the rate returns to `v0`. Returns the
longest `T` after which it still persists, by bisection; `Inf` if it persists after `Tmax`
(as it does for `v1` below the critical rate).
"""
function xu_shen_max_duration(m::XuShen; v0::Real, v1::Real, Tmax::Real = 1e5, steps::Integer = 40, kw...)
    e = xu_shen_equilibrium(m, v0)
    e === nothing && throw(ArgumentError("no viable equilibrium at v0 = $v0 (critical rate $(xu_shen_critical_rate(m).v))"))
    ok(T) = _xs_ok(m, e.Δ, e.V, e.N; v = v0, phases = ((v1, T),), kw...)
    ok(float(Tmax)) && return Inf
    a, b = 0.0, float(Tmax)
    for _ in 1:steps
        mid = (a + b) / 2
        ok(mid) ? (a = mid) : (b = mid)
    end
    (a + b) / 2
end
