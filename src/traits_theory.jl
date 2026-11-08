# Theory for several traits (Lande 1979, 1980; Jones, Arnold and Bürger 2003): the
# multivariate Gaussian mutation–selection balance, the steady lag behind an optimum that
# moves in a direction, and the response deflected toward the genetic line of least
# resistance (Schluter 1996). `gaussian_G` is the fixed point of the aggregate recursion
# (`run_moments(p; loci = :aggregate)`, exact for one haploid locus); the trait moments
# lowering's default, the multilocus model, converges to `gaussian_G_loci` (multilocus.jl),
# which is `gaussian_G` with `2L·M` in place of `M`.
#
# Notation: G the additive genetic covariance matrix, P = G + V_E I the phenotypic one,
# Ω the covariance of the Gaussian fitness function, A = Ω + V_E I (the matrix analogue of
# V_s = ω² + V_E), M the mutational covariance per generation. Under Gaussian selection a
# Gaussian population's genotypes are again Gaussian, with
#
#     Δz̄ = G (A + G)⁻¹ (θ − z̄),       G_survivors = G − G (A + G)⁻¹ G,
#
# and with mutation G' = G − G (A + G)⁻¹ G + M.

"""
    selection_variance_matrix(p; ω = initFitSD(p))

`A = Ω + V_E I`, the matrix analogue of [`selection_variance`](@ref) `V_s = ω² + V_E`: the
covariance of the fitness function on genotypes, once the environmental deviation is
integrated out.
"""
selection_variance_matrix(p::Params; ω::Real = initFitSD(p)) =
    selection_covariance(p; ω) + envVar(p) * Matrix{Float64}(I, p.NoTraits, p.NoTraits)

symsqrt(A) = (F = eigen(Symmetric(Matrix{Float64}(A))); F.vectors * (sqrt.(max.(F.values, 0.0)) .* F.vectors'))

"""
    gaussian_G(M, A)
    gaussian_G(p)

The equilibrium G matrix of the aggregate Gaussian recursion `G' = G − G (A + G)⁻¹ G + M`,
the multivariate [`gaussian_vg_exact`](@ref) (exact for one haploid locus; for `L`
recombining diploid loci the equilibrium is [`gaussian_G_loci`](@ref), `gaussian_G(2L·M, A)`): the unique positive semidefinite solution of
`G (A + G)⁻¹ G = M`,

    G = A^{1/2} X A^{1/2},   X = ½ (M̃ + (M̃² + 4 M̃)^{1/2}),   M̃ = A^{-1/2} M A^{-1/2}.

(Substituting `G = A^{1/2} X A^{1/2}` turns the equation into `X (I + X)⁻¹ X = M̃`; `X`
commutes with `(I + X)⁻¹`, so `M̃ = X²(I + X)⁻¹` is a function of `X`, and as
`x ↦ x²/(1 + x)` increases on `[0, ∞)`, `X` is the function of `M̃` taking each eigenvalue
`μ` to `(μ + √(μ² + 4μ))/2`.) When `M` and `A` share eigenvectors, `G` shares them too,
with eigenvalues `(m_i + √(m_i² + 4 m_i a_i))/2`, the one-trait formula on each axis. For
weak mutation (`M̃ ≪ I`) `G ≈ A^{1/2} M̃^{1/2} A^{1/2}`, the matrix geometric mean of `A`
and `M`, Lande's (1980) `√(M V_s)`.
"""
function gaussian_G(M::AbstractMatrix, A::AbstractMatrix)
    Ah = symsqrt(A)
    Aih = inv(Symmetric(Ah))
    Mt = Symmetric(Aih * M * Aih)
    F = eigen(Mt)
    x = [(μ + sqrt(max(μ^2 + 4μ, 0.0))) / 2 for μ in max.(F.values, 0.0)]
    G = Ah * (F.vectors * (x .* F.vectors')) * Ah
    (G + G') / 2
end
gaussian_G(p::Params) = gaussian_G(mutational_covariance(p), selection_variance_matrix(p))

"""
    gaussian_G_iterate(M, A; G0 = M, tol = 1e-13, maxiter = 10^7, erosion = 0.0)

The fixed point of the Gaussian recursion found by iterating it, `G ← (G − G (A + G)⁻¹ G)
(1 − erosion) + M`, from `G0` until no entry changes by more than `tol` relative to the
largest: a numerical check on [`gaussian_G`](@ref), and with `erosion = 1/(2N_e)` the
equilibrium of the aggregate recursion with drift's loss of variance (as
`run_moments(…; loci = :aggregate, erosion = true)` with a constant `N_e`).
"""
function gaussian_G_iterate(M::AbstractMatrix, A::AbstractMatrix; G0 = Matrix{Float64}(M), tol = 1e-13,
                            maxiter::Integer = 10^7, erosion::Real = 0.0)
    G = Matrix{Float64}(G0)
    for _ in 1:maxiter
        Gs = G - G * ((A + G) \ G)
        Gn = Gs .* (1 - erosion) .+ M
        Gn = (Gn + Gn') / 2
        maximum(abs, Gn - G) <= tol * max(maximum(abs, Gn), eps()) && return Gn
        G = Gn
    end
    error("gaussian_G_iterate did not converge in $maxiter iterations")
end

"""
    selection_gradient(λ; G, A)

The selection gradient on the mean, `β = (A + G)⁻¹ λ` (Lande and Arnold 1983), with
`λ = θ − z̄` the vector from the mean to the optimum: the gradient of `ln W̄` for Gaussian
fitness of covariance `Ω = A − V_E I`.
"""
selection_gradient(λ::AbstractVector; G::AbstractMatrix, A::AbstractMatrix) = (A + G) \ λ

"""
    trait_response(λ; G, A)

Lande's response of the mean to one generation of selection, `Δz̄ = G β = G (A + G)⁻¹ λ`.
"""
trait_response(λ::AbstractVector; G::AbstractMatrix, A::AbstractMatrix) = G * selection_gradient(λ; G, A)

"""
    steady_lag_vector(k; G, A)
    steady_lag_vector(p; G)

The steady lag `λ = θ − z̄` of a population tracking an optimum that moves by the vector
`k` per generation (`MoveType = :continuous`, `k = MoveFitnessSD·u`): the mean must move by
`k` too, `G (A + G)⁻¹ λ = k`, so `λ = (A + G) G⁻¹ k`. For one trait, [`steady_lag`](@ref)
`k (V_G + V_s)/V_G`. The lag is longest along the directions of least genetic variance:
in the eigenbasis of `G`, when `A` shares it, each component of `k` is scaled by
`(a_i + g_i)/g_i`.
"""
steady_lag_vector(k::AbstractVector; G::AbstractMatrix, A::AbstractMatrix) = (A + G) * (G \ k)
function steady_lag_vector(p::Params; G::AbstractMatrix)
    p.MoveType === :continuous ||
        throw(ArgumentError("a steady lag needs a constant rate of change, MoveType = :continuous; got :$(p.MoveType)"))
    steady_lag_vector(p.MoveFitnessSD .* move_direction(p); G, A = selection_variance_matrix(p))
end

"""
    gmax(G)

The leading eigenvector of `G`, of unit length with a non-negative first nonzero entry:
`g_max`, the genetic line of least resistance (Schluter 1996).
"""
function gmax(G::AbstractMatrix)
    F = eigen(Symmetric(Matrix{Float64}(G)))
    v = F.vectors[:, argmax(F.values)]
    i = findfirst(x -> abs(x) > 1e-12, v)
    i !== nothing && v[i] < 0 ? -v : v
end

"The angle between two vectors in degrees, in [0, 180]."
vector_angle(a::AbstractVector, b::AbstractVector) =
    acosd(clamp(dot(a, b) / (norm(a) * norm(b)), -1.0, 1.0))

"""
    response_deflection(λ; G, A)

How selection toward an optimum at `λ = θ − z̄` from the mean is deflected by `G`: a
named tuple of the selection gradient `β`, the response `Δz̄ = G β` and the angles (in
degrees) between `Δz̄` and `β`, `Δz̄` and the axis `±g_max` (at most 90°), and `β` and
`±g_max`. The response turns from `β` toward `g_max` (Schluter 1996): with `G = Σ g_i e_i e_iᵀ`,
`Δz̄ = Σ g_i (e_i·β) e_i` scales the component along `e_i` by `g_i`, so the angle to the
`g_max` axis is never larger for the response than for the gradient.
"""
function response_deflection(λ::AbstractVector; G::AbstractMatrix, A::AbstractMatrix)
    β = selection_gradient(λ; G, A)
    r = G * β
    g = gmax(G)
    axis(v) = min(vector_angle(v, g), vector_angle(v, -g))
    (β = β, response = r, angle_response_gradient = vector_angle(r, β),
     angle_response_gmax = axis(r), angle_gradient_gmax = axis(β))
end


# --- fluctuating optima with several traits ----------------------------------------

"""
    red_noise_lag_cov(Σ, ρ; G, A)

The stationary covariance matrix of the lag `λ = θ − z̄` when the optimum vector is AR(1)
noise with stationary covariance `Σ` and autocorrelation `ρ`,
`θ_t = ρ θ_{t−1} + √(1 − ρ²) η_t`, `η_t ~ N(0, Σ)`, and the mean follows Lande's
`z̄_{t+1} = z̄_t + B (θ_t − z̄_t)` with `B = G (A + G)⁻¹`. The lag and the optimum form the
linear system

    λ_{t+1} = (I − B) λ_t + (ρ − 1) θ_t + √(1 − ρ²) η_{t+1},
    θ_{t+1} = ρ θ_t + √(1 − ρ²) η_{t+1},

whose stationary covariance `X` solves the discrete Lyapunov equation `X = F X Fᵀ + Q`;
the lag's is its upper-left block. For one trait it is [`red_noise_lag_var`](@ref),
`2σ²(1 − ρ)/((2 − β)(1 − ρ(1 − β)))`. When `Σ` and `B` share eigenvectors each axis has the
one-trait value with its own `β`; otherwise noise along a direction of little genetic
variance is tracked worst, and the lag's covariance turns away from `g_max`.
"""
function red_noise_lag_cov(Σ::AbstractMatrix, ρ::Real; G::AbstractMatrix, A::AbstractMatrix)
    n = size(G, 1)
    B = G / (A + G)
    Id = Matrix{Float64}(I, n, n)
    F = [Id-B (ρ - 1) .* Id; zeros(n, n) ρ .* Id]
    Q = (1 - ρ^2) .* [Σ Σ; Σ Σ]
    X = reshape((Matrix{Float64}(I, 4n^2, 4n^2) - kron(F, F)) \ vec(Q), 2n, 2n)
    Λ = X[1:n, 1:n]
    (Λ + Λ') / 2
end

"""
    periodic_lag_cov(P, f; G, A)

The covariance matrix of the lag `λ = θ − z̄`, averaged over a cycle, once the response is
steady, when the optimum cycles as `θ_t = sin(2π t/P)·f` (amplitude 1 along the vector `f`;
multiply by `A²` for amplitude `A`) and the mean follows `z̄_{t+1} = z̄_t + B (θ_t − z̄_t)`,
`B = G (A + G)⁻¹`. With `ν = 2π/P` the lag is the real part of `l e^{iνt}` with
`l = (I − e^{−iν}(I − B))⁻¹ (1 − e^{−iν}) f`, so its mean outer product is `Re(l l*)/2`.
For one trait it is `lag_amplitude²/2` of [`periodic_response`](@ref).
"""
function periodic_lag_cov(P::Real, f::AbstractVector; G::AbstractMatrix, A::AbstractMatrix)
    n = size(G, 1)
    B = G / (A + G)
    ν = 2π / P
    e = exp(-im * ν)
    l = (Matrix{ComplexF64}(I, n, n) - e .* (I - B)) \ ((1 - e) .* f)
    Λ = real(l * l') ./ 2
    (Λ + Λ') / 2
end

"""
    expected_lag_cov(p; G)

The expected outer product of the lag `E[λλᵀ]`, `λ = θ − z̄`, of a population of several
traits with genetic covariance `G` under the move of `p`, once the response is steady, with
`A` = [`selection_variance_matrix`](@ref)`(p)`: the outer product of the
[`steady_lag_vector`](@ref) of the trend `d·u` (`MoveFitnessSD` along
[`move_direction`](@ref)), plus, for `MoveType = :red_noise`, [`red_noise_lag_cov`](@ref)
with `Σ` = [`move_noise_covariance`](@ref)`(p)`, or, for `:periodic`, `A²` times
[`periodic_lag_cov`](@ref) along [`fluct_direction`](@ref). The trend and the fluctuation add
because the response is linear and the fluctuation has mean zero. The matrix version of
[`expected_lag2`](@ref).
"""
function expected_lag_cov(p::Params; G::AbstractMatrix)
    p.MoveType in (:continuous, :red_noise, :periodic) ||
        throw(ArgumentError("expected_lag_cov is for an exogenous optimum (:continuous, :red_noise, :periodic), got :$(p.MoveType)"))
    A = selection_variance_matrix(p)
    λ = steady_lag_vector(p.MoveFitnessSD .* move_direction(p); G, A)
    trend = λ * λ'
    p.MoveType === :red_noise && return trend + red_noise_lag_cov(move_noise_covariance(p), p.MoveAutocorr; G, A)
    p.MoveType === :periodic && return trend + p.MoveAmplitude^2 .* periodic_lag_cov(p.MovePeriod, fluct_direction(p); G, A)
    trend
end

"""
    trait_lag_load(p; G)

The expected lag load of several traits on the log scale, `tr((Ω + P)⁻¹ E[λλᵀ])/2`, with
`P = G + V_E I`, `Ω` the initial fitness covariance and `E[λλᵀ]` from
[`expected_lag_cov`](@ref): Gaussian selection gives `ln W̄(λ) = ln W̄(0) − λᵀ(Ω + P)⁻¹λ/2`.
For one trait it is [`lag_load`](@ref).
"""
function trait_lag_load(p::Params; G::AbstractMatrix)
    n = size(G, 1)
    s = selection_covariance(p) + G + envVar(p) * Matrix{Float64}(I, n, n)
    tr(s \ expected_lag_cov(p; G)) / 2
end
