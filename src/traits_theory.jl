# Theory for several traits (Lande 1979, 1980; Jones, Arnold and Bürger 2003): the
# multivariate Gaussian mutation–selection balance, the steady lag behind an optimum that
# moves in a direction, and the response deflected toward the genetic line of least
# resistance (Schluter 1996). These are the fixed points of the trait moments lowering
# (ext/QTLSimCatlabExt/trait_moments.jl), whose one-trait case is the scalar recursion
# of theory.jl.
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

The equilibrium G matrix of the Gaussian recursion `G' = G − G (A + G)⁻¹ G + M`, the
multivariate [`gaussian_vg_exact`](@ref): the unique positive semidefinite solution of
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
equilibrium with drift's loss of variance (as `run_moments(…; erosion = true)` with a
constant `N_e`).
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
