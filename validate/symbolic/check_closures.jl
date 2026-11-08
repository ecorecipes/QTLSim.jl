# Symbolic checks, with Symbolics.jl, of the algebra behind the moments lowering's closures,
# and numerical spot checks that tie each identity to the Julia functions that use it.
#
# Usage (from the package root):
#     julia --project=validate/symbolic -e 'using Pkg; Pkg.instantiate()'    # once
#     julia --project=validate/symbolic validate/symbolic/check_closures.jl
# Writes validate/symbolic/check_closures.log; exits 1 if any check fails.
#
# (a) the rare-alleles closure with locally favoured alleles (src/local_alleles.jl):
#     the marginal fitness of a mutant in a lagging deme (completing the square), the
#     second-order growth rate of a rare lineage under selection and migration (the
#     characteristic polynomial of P(I + εΣ) vanishes to O(ε³) at 1 + ε π·σ + ε² σᵀWσ),
#     the a² coefficient that gives the closure, the invariance of σᵀWσ under a common
#     shift of the lags, the island form (1 − m′)/m′·Var_π, and the reductions to the
#     current closure (one deme, lags alike, no migration, complete migration);
# (b) the existing closures: the rare-alleles per-gene fixed point (rationalised and root
#     forms) and Hermisson and Wagner's (2004) equation A4, the multivariate closure's
#     exactness ½(𝓡X + X𝓡ᵀ) = ½C for 𝓡 = ½CX⁻¹, its one-trait reduction α²/(2s), the
#     two-trait house-of-cards weights (Σλ_i e_i = 1), and the multilocus Gaussian fixed
#     point V²/(V_s + V) = 2L·V_m;
# (c) numerical spot checks at random parameter values of each identity against
#     house_of_gauss_vg, gaussian_vg_loci, rare_alleles_rate, rare_hc_matrix,
#     rare_hc_rates, local_weights, local_gain, local_sinv and hap_select.
#
# Square roots: an identity in r = √D is checked as a polynomial in r reduced modulo
# r² = D (both its even and odd parts must vanish).

using Symbolics, LinearAlgebra, Random, Printf
using QTLSim

const LOG = joinpath(@__DIR__, "check_closures.log")
const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s); flush(stdout))
const failures = String[]
function check(name, ok::Bool)
    ok || push!(failures, name)
    say(@sprintf("  %-72s %s", name, ok ? "ok" : "FAIL"))
    ok
end

# --- symbolic helpers ------------------------------------------------------------------

simp(e) = Symbolics.simplify(Symbolics.simplify_fractions(Symbolics.expand(e)); expand = true)
iszero_sym(e) = (z = simp(e); isequal(z, 0) || (Symbolics.value(z) isa Number && iszero(Symbolics.value(z))))
iszero_sym(e::AbstractArray) = all(iszero_sym, e)

# the coefficient of x^j in a polynomial e (by derivatives at 0)
function coeffj(e, x, j)
    d = e
    for _ in 1:j
        d = Symbolics.expand_derivatives(Differential(x)(d))
    end
    simp(Symbolics.substitute(d, Dict(x => 0)) / factorial(j))
end

# e, a polynomial in r of degree ≤ deg, reduced modulo r² = D: is it zero?
function iszero_mod_sqrt(e, r, D; deg = 6)
    e = Symbolics.expand(e)
    ev, od = 0, 0
    for j in 0:deg
        c = coeffj(e, r, j)
        if iseven(j)
            ev += c * D^(j ÷ 2)
        else
            od += c * D^(j ÷ 2)
        end
    end
    iszero_sym(Symbolics.expand(e - sum(coeffj(e, r, j) * r^j for j in 0:deg))) || return false   # degree check
    iszero_sym(ev) && iszero_sym(od)
end

num(e, d) = Float64(Symbolics.value(Symbolics.substitute(e, d)))

inv2(A) = [A[2, 2] -A[1, 2]; -A[2, 1] A[1, 1]] ./ (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
det2(A) = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]

const rng = Xoshiro(20261108)
runif(a, b) = a + (b - a) * rand(rng)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_closures.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the closure with locally favoured alleles (src/local_alleles.jl)")
    @variables y a θ zb P ω2 s λ
    # a mutant of effect a in a background y ~ N(z̄, P) under W ∝ exp(−(y + a − θ)²/(2ω²)):
    # completing the square in y gives exp(−(y − y*)²/(2Pω²/s) − (a − λ)²/(2s)), s = ω² + P, λ = θ − z̄
    sP = ω2 + P
    ystar = (zb * ω2 + (θ - a) * P) / sP
    lhs = -(y - zb)^2 / (2P) - (y + a - θ)^2 / (2ω2)
    rhs = -(y - ystar)^2 / (2P * ω2 / sP) - (a - (θ - zb))^2 / (2sP)
    check("completing the square: the background-averaged fitness of a mutant", iszero_sym(lhs - rhs))
    # so the mutant's log marginal fitness relative to the resident is aλ/s − a²/(2s), exactly
    check("log marginal fitness −((a − λ)² − λ²)/(2s) = aλ/s − a²/(2s)",
          iszero_sym(-((a - λ)^2 - λ^2) / (2s) - (a * λ / s - a^2 / (2s))))

    # the growth rate of a rare lineage, two demes, general P = [1−p p; q 1−q]
    @variables p q σ1 σ2 ε r1 r2 ℓ1 ℓ2 κ1 κ2 c m
    Pm = [1-p p; q 1-q]
    πv = [q, p] ./ (p + q)
    check("π is stationary: πP = π", iszero_sym(Pm' * πv - πv))
    one_π = [1, 1] * πv'
    Ag = inv2(Matrix{Num}(I, 2, 2) - Pm + one_π) - one_π              # the group inverse (I − P)^#
    check("(I − P)^# is the group inverse: A A^# A = A, A^# A A^# = A^#, A A^# = A^# A",
          iszero_sym((I - Pm) * Ag * (I - Pm) - (I - Pm)) && iszero_sym(Ag * (I - Pm) * Ag - Ag) &&
          iszero_sym((I - Pm) * Ag - Ag * (I - Pm)))
    W = simp.(Diagonal(πv) * Ag * Pm)
    Wclosed = p * q * (1 - p - q) / (p + q)^3 .* [1 -1; -1 1]
    check("two demes: W = Π(I − P)^#P = pq(1 − p − q)/(p + q)³ · [1 −1; −1 1]", iszero_sym(W - Wclosed))
    check("W's rows and columns sum to zero (so σᵀWσ sees only deviations; Lean gain_shift)",
          iszero_sym(W * [1, 1]) && iszero_sym([1 1] * W))
    W = Wclosed
    σv = [σ1, σ2]
    λ1 = dot(πv, σv)
    λ2 = σv' * W * σv
    A = Pm * (I + ε * Diagonal(σv))
    @variables t3
    lam = 1 + ε * λ1 + ε^2 * λ2 + ε^3 * t3
    χ = det2(A - lam * I)
    check("det(P(I + εΣ) − λ(ε)I) = O(ε³) at λ = 1 + ε π·σ + ε² σᵀWσ",
          all(j -> iszero_sym(coeffj(χ, ε, j)), 0:2))
    # with σ_d = aℓ_d − a²κ_d, the a² coefficient of r(a) = π·σ + σᵀWσ is −π·κ + ℓᵀWℓ:
    # a rare mutant is selected against at ½a²(s̄⁻¹ − 2Λ), κ_d = 1/(2s_d), Λ = ℓᵀWℓ
    ra = dot(πv, [a * ℓ1 - a^2 * κ1, a * ℓ2 - a^2 * κ2]) +
         [a * ℓ1 - a^2 * κ1, a * ℓ2 - a^2 * κ2]' * W * [a * ℓ1 - a^2 * κ1, a * ℓ2 - a^2 * κ2]
    Λ2 = [ℓ1, ℓ2]' * W * [ℓ1, ℓ2]
    check("a² coefficient of r(a): −π·κ + ℓᵀWℓ", iszero_sym(coeffj(ra, a, 2) - (-(dot(πv, [κ1, κ2])) + Λ2)))
    check("a¹ coefficient of r(a): π·ℓ, the common lag (left out)", iszero_sym(coeffj(ra, a, 1) - dot(πv, [ℓ1, ℓ2])))
    check("σᵀWσ is unchanged by a common shift of the lags",
          iszero_sym(([ℓ1 + c, ℓ2 + c]' * W * [ℓ1 + c, ℓ2 + c]) - Λ2))
    check("reduction, lags alike (ℓ₁ = ℓ₂): Λ = 0", iszero_sym(Symbolics.substitute(Λ2, Dict(ℓ2 => ℓ1))))
    # island model of two demes: p = q = m (forward rate m to the other deme), m′ = 2m
    Λisl = Symbolics.substitute(Λ2, Dict(p => m, q => m))
    varπ = ((ℓ1 - ℓ2) / 2)^2
    check("island (k = 2): Λ = (1 − m′)/m′ · Var_π(ℓ), m′ = 2m", iszero_sym(Λisl - (1 - 2m) / (2m) * varπ))
    check("reduction, complete migration (m′ = 1): Λ = 0",
          iszero_sym(Symbolics.substitute((1 - 2m) / (2m) * varπ, Dict(m => 1 // 2))))
    # the closure: the rare-alleles loss rate α²·[s⁻¹ − 2Λ]₊/2 is the current α²/(2s) at Λ = 0
    @variables α2 sinv Λs
    check("reduction, Λ = 0: α²(s⁻¹ − 2Λ)/2 = α²/(2s)",
          iszero_sym(Symbolics.substitute(α2 * (sinv - 2Λs) / 2, Dict(Λs => 0)) - α2 * sinv / 2))
    # several traits: ℓ_d vectors, Λ = Σ W_dc ℓ_d ℓ_cᵀ; Λ = 0 for ℓ_d all equal (W's rows and columns
    # weighted by the lags sum to zero against a common vector)
    @variables u1 u2 v1 v2
    L1, L2 = [u1, u2], [v1, v2]
    Λm = W[1, 1] * L1 * L1' + W[1, 2] * L1 * L2' + W[2, 1] * L2 * L1' + W[2, 2] * L2 * L2'
    check("several traits, lags alike: Λ = Σ W_dc ℓ_d ℓ_cᵀ = 0",
          iszero_sym(Symbolics.substitute.(Λm, Ref(Dict(v1 => u1, v2 => u2)))))
    check("several traits: aᵀΛa = σᵀWσ for σ_d = aᵀℓ_d (the quadratic form of the gain)",
          iszero_sym([r1, r2]' * Λm * [r1, r2] - [dot([r1, r2], L1), dot([r1, r2], L2)]' * W * [dot([r1, r2], L1), dot([r1, r2], L2)]))

    # ===================================================================================
    say("\n(b) the existing closures")
    @variables μ v Vs ng rr
    b = v / (2Vs) + 1 / ng
    D = b^2 + 4μ * v / Vs
    croot = Vs * (rr - b) / 2                  # the positive root of c²/V_s + b c − μα² = 0, rr = √D
    check("rare alleles: the root form solves c²/V_s + c(α²/(2V_s) + 1/ng) = μα²",
          iszero_mod_sqrt(croot^2 / Vs + b * croot - μ * v, rr, D))
    check("rare alleles: root form = rationalised form 2μα²/(b + √D) (house_of_gauss_vg)",
          iszero_mod_sqrt(croot * (b + rr) - 2μ * v, rr, D))
    se = ng * v / (2Vs)
    Θ = 4ng * μ
    DA4 = se^2 + 2(1 + Θ) * se + 1
    xA4 = (rr - se - 1) / (4se)                # Hermisson and Wagner (2004) eq. A4, c/α², rr = √DA4
    check("Hermisson–Wagner A4: c = α²(√(s_e² + 2(1+Θ)s_e + 1) − s_e − 1)/(4s_e) is the fixed point",
          iszero_mod_sqrt(simp(ng * ((xA4 * v)^2 / Vs + b * xA4 * v - μ * v)), rr, DA4))
    # multivariate closure: ½(𝓡X + X𝓡ᵀ) = ½C for 𝓡 = ½CX⁻¹, C and X symmetric
    @variables c11 c12 c22 x11 x12 x22
    Cm, Xm = [c11 c12; c12 c22], [x11 x12; x12 x22]
    Rm = Cm * inv2(Xm) ./ 2
    check("multivariate closure exact at the house of cards: ½(𝓡X + X𝓡ᵀ) = ½C", iszero_sym((Rm * Xm + Xm * Rm') ./ 2 - Cm ./ 2))
    # one trait: X = E[a²/(a² s⁻¹)] = s, so 𝓡 = ½α²/s
    check("one trait: 𝓡 = ½ α² X⁻¹ with X = s is α²/(2s)", iszero_sym(v * (1 / Vs) / 2 - v / (2Vs)))
    @variables ρ1 ρ2
    e1, e2 = 1 / (ρ1 * (ρ1 + ρ2)), 1 / (ρ2 * (ρ1 + ρ2))
    check("two traits: e_i = 1/(√λ_i(√λ₁+√λ₂)) gives Σ λ_i e_i = 1 (tr(s⁻¹X) = 1)", iszero_sym(ρ1^2 * e1 + ρ2^2 * e2 - 1))
    @variables U Vm
    Vroot = (U * Vm + rr) / 2                  # gaussian_vg_closed, rr = √(U²V_m² + 4U V_m V_s), U = 2L
    check("multilocus Gaussian: V = (U V_m + √(U²V_m² + 4U V_m V_s))/2 gives V²/(V_s + V) = U V_m",
          iszero_mod_sqrt(Vroot^2 - U * Vm * (Vs + Vroot), rr, U^2 * Vm^2 + 4U * Vm * Vs))
    @variables L H
    check("multilocus Gaussian: genic change −V²/(2Ls) + V_m = 0 ⇔ V²/(V_s + V) = 2L V_m",
          iszero_sym((-H^2 / (2L * (Vs + H)) + Vm) * (-2L) - (H^2 / (Vs + H) - 2L * Vm)))

    # ===================================================================================
    say("\n(c) numerical spot checks against the Julia implementation")
    for trial in 1:3
        L0 = rand(rng, 5:60); μ0 = 10^runif(-5, -2); Ms = 10^runif(-3, -1); N0 = rand(rng, 50:2000)
        ω0 = runif(1.0, 6.0)
        pp = Params(Ploidy = 2, NoQTL = L0, MutationRatePerGamete = μ0, MutationScalar = Ms, PopSize = N0,
                    PhenotypeVar = runif(0.5, 3.0), Heritability = runif(0.2, 0.8), InitFitnessSD = ω0, RecombPoissonMean = 2.0)
        α20, Vs0 = QTLSim.mutation_effect_variance(pp), selection_variance(pp)
        ng0 = 2effective_size(pp.PopSize, pp.PopSize)
        nu = 2L0
        b0 = α20 / (2Vs0) + 1 / ng0
        d = Dict(v => α20, Vs => Vs0, ng => ng0, μ => μ0, rr => sqrt(b0^2 + 4μ0 * α20 / Vs0))
        cs = num(croot, d)
        check(@sprintf("house_of_gauss_vg = 2L·(root form), L=%d μ=%.2g N=%d", L0, μ0, N0),
              isapprox(house_of_gauss_vg(pp), nu * cs; rtol = 1e-10))
        se0 = ng0 * α20 / (2Vs0)
        dA = Dict(v => α20, Vs => Vs0, ng => ng0, μ => μ0, rr => sqrt(num(DA4, Dict(v => α20, Vs => Vs0, ng => ng0, μ => μ0))))
        check("house_of_gauss_vg = 2L·α²·(eq. A4)", isapprox(house_of_gauss_vg(pp), nu * α20 * num(xA4, dA); rtol = 1e-10))
        Vm0 = muteVar(pp)
        dV = Dict(U => nu, Vm => Vm0, Vs => Vs0, rr => sqrt(nu^2 * Vm0^2 + 4nu * Vm0 * Vs0))
        Vg = gaussian_vg_loci(pp)
        check("gaussian_vg_loci = root form, and V²/(V_s+V) = 2L V_m",
              isapprox(Vg, num(Vroot, dV); rtol = 1e-12) && isapprox(Vg^2 / (Vs0 + Vg), nu * Vm0; rtol = 1e-10))
        s1 = runif(0.01, 2.0)
        check("rare_alleles_rate, one trait: α²·s⁻¹/2", isapprox(QTLSim.rare_alleles_rate(fill(α20, 1, 1), fill(s1, 1, 1))[1, 1],
                                                                min(α20 * s1 / 2, 1.0); rtol = 1e-12))
        # two traits
        Bc = randn(rng, 2, 2); C0 = Bc * Bc' .* 1e-3 + 1e-4I
        Bs = randn(rng, 2, 2); S0 = Bs * Bs' + 2I
        Si = inv(S0)
        X0 = QTLSim.rare_hc_matrix(C0, Si)
        R0 = QTLSim.rare_alleles_rate(C0, Si)
        dC = Dict(c11 => C0[1, 1], c12 => C0[1, 2], c22 => C0[2, 2], x11 => X0[1, 1], x12 => X0[1, 2], x22 => X0[2, 2])
        Rs = [num(Rm[i, j], dC) for i in 1:2, j in 1:2]
        check("rare_alleles_rate = ½ C X⁻¹ (symbolic 𝓡 at rare_hc_matrix's X), two traits", isapprox(R0, Rs; rtol = 1e-8))
        check("½(𝓡X + X𝓡ᵀ) = ½C with rare_alleles_rate and rare_hc_matrix", isapprox((R0 * X0 + X0 * R0') ./ 2, C0 ./ 2; rtol = 1e-8))
        check("tr(s⁻¹ X) = 1 for rare_hc_matrix", isapprox(tr(Si * X0), 1.0; rtol = 1e-10))
        lam = rand(rng, 2) .* 3
        ρh = QTLSim.rare_hc_rates(lam)
        dl = Dict(ρ1 => sqrt(lam[1]), ρ2 => sqrt(lam[2]))
        check("rare_hc_rates (two traits) = 1/(2e_i) with the symbolic e_i", isapprox(ρh, [1 / (2num(e1, dl)), 1 / (2num(e2, dl))]; rtol = 1e-12))
        lam3 = rand(rng, 3) .* 3
        ρ3 = QTLSim.rare_hc_rates(lam3)
        check("rare_hc_rates (three traits, quadrature): Σ λ_i e_i = 1", isapprox(sum(lam3 ./ (2 .* ρ3)), 1.0; rtol = 1e-8))
    end

    # the gain: W against the symbolic two-deme W, and against the eigenvalue
    for trial in 1:3
        p0, q0 = runif(0.01, 0.4), runif(0.01, 0.4)
        P2 = [1-p0 p0; q0 1-q0]
        Ws = [num(W[i, j], Dict(p => p0, q => q0)) for i in 1:2, j in 1:2]
        check(@sprintf("local_weights = symbolic W (two demes, p=%.3f q=%.3f)", p0, q0), isapprox(QTLSim.local_weights(P2), Ws; rtol = 1e-10))
    end
    for k in (2, 3, 4, 6)
        M = rand(rng, k, k) .* 0.05; M[diagind(M)] .= 0
        Pk = I - Diagonal(vec(sum(M; dims = 2))) + M          # a backward migration matrix, rows sum to 1
        σk = 0.01 .* randn(rng, k)          # selection weak against migration, as the expansion assumes
        Wk = QTLSim.local_weights(Pk)
        lead(e) = maximum(real, eigvals(Pk * (I + e * Diagonal(σk))))
        h = 1e-3
        d2 = (lead(h) - 2lead(0.0) + lead(-h)) / (2h^2)        # λ″(0)/2
        d1 = (lead(h) - lead(-h)) / (2h)
        πk = [I - Pk'; ones(1, k)] \ [zeros(k); 1.0]
        check("leading eigenvalue of P(I + εΣ), k = $k: λ′ = π·σ, λ″/2 = σᵀWσ (local_weights)",
              isapprox(d1, dot(πk, σk); rtol = 1e-6) && isapprox(d2, σk' * Wk * σk; rtol = 1e-4, atol = 1e-8))
    end
    for k in (2, 4, 7)
        m0 = runif(0.002, 0.5)
        Pk = QTLSim.migration_mixing(island_migration(k, m0), fill(100.0, k), trues(k))
        ℓ = randn(rng, k) .* 0.01
        mp = 1 - (Pk[1, 1] - Pk[1, 2])                       # 1 − m′, the eigenvalue on the deviations
        Λk = QTLSim.local_gain(Pk, ℓ)
        varℓ = sum((ℓ .- sum(ℓ) / k) .^ 2) / k
        check(@sprintf("local_gain, island k = %d, m = %.3f: Λ = (1 − m′)/m′·Var_π(ℓ)", k, m0),
              all(x -> isapprox(x, (1 - mp) / mp * varℓ; rtol = 1e-8), Λk))
    end
    Pk = QTLSim.migration_mixing(island_migration(4, 0.05), fill(100.0, 4), trues(4))
    check("local_gain = 0 exactly for lags alike (one trait)", all(iszero, QTLSim.local_gain(Pk, fill(0.0123, 4))))
    check("local_gain = 0 exactly for lags alike (two traits)", all(iszero, QTLSim.local_gain(Pk, fill([0.01, -0.02], 4))))
    check("local_gain = 0 exactly without migration (P = I)", all(iszero, QTLSim.local_gain(Matrix(1.0I, 4, 4), randn(rng, 4))))
    check("local_gain = 0 exactly for one deme", all(iszero, QTLSim.local_gain(fill(1.0, 1, 1), [0.3])))
    Pb = [0.9 0.1 0 0; 0.1 0.9 0 0; 0 0 0.8 0.2; 0 0 0.2 0.8]
    ℓb = [0.01, -0.01, 0.02, 0.02]
    Λb = QTLSim.local_gain(Pb, ℓb)
    check("local_gain per migration class: two classes, the second with lags alike",
          Λb[3] == 0 && Λb[4] == 0 && Λb[1] == Λb[2] && isapprox(Λb[1], QTLSim.local_gain(Pb[1:2, 1:2], ℓb[1:2])[1]; rtol = 1e-12))
    check("local_sinv(s⁻¹, 0) === s⁻¹; local_sinv = [s⁻¹ − 2Λ]₊",
          QTLSim.local_sinv(0.04, 0.0) === 0.04 && QTLSim.local_sinv(0.04, 0.01) ≈ 0.02 && QTLSim.local_sinv(0.04, 0.03) == 0.0)
    Sm = [0.05 0.01; 0.01 0.03]
    check("local_sinv, two traits: Λ = 0 gives s⁻¹ itself; the positive part is PSD",
          QTLSim.local_sinv(Sm, zeros(2, 2)) === Sm && minimum(eigvals(QTLSim.local_sinv(Sm, [0.04 0; 0 -0.01]))) >= -1e-15)
    # hap_select: gain = nothing and gain = 0 give identical moments
    pp = Params(Ploidy = 2, NoQTL = 20, MutationRatePerGamete = 1e-4, PopSize = 200, RecombPoissonMean = 2.0)
    ω, ρ = QTLSim.pair_classes(pp)
    st = QTLSim.hap_init(0.3, 200.0, 2, length(ω))
    α20 = QTLSim.mutation_effect_variance(pp)
    s1 = QTLSim.hap_select(st, 0.03, 200.0, ω, 2; α2 = α20)
    s2 = QTLSim.hap_select(st, 0.03, 200.0, ω, 2; α2 = α20, gain = 0.0)
    s3 = QTLSim.hap_select(st, 0.03, 200.0, ω, 2; α2 = α20, gain = 0.005)
    s4 = QTLSim.hap_select(st, 0.02, 200.0, ω, 2; α2 = α20)
    check("hap_select: gain = 0 is identical to the current closure", s1 == s2)
    check("hap_select: gain Λ with β = 1 equals the current closure at s⁻¹ − 2Λ in the rare loss",
          isapprox(s3[1][1][1] - s1[1][1][1], s4[1][1][1] - s1[1][1][1] - (QTLSim.hap_select(st, 0.02, 200.0, ω, 2)[1][1][1] -
                                                                             QTLSim.hap_select(st, 0.03, 200.0, ω, 2)[1][1][1]); rtol = 1e-10))
    # the marginal fitness, by quadrature with QTSim's fitness
    for trial in 1:3
        zb0, P0, ω0, θ0, a0 = randn(rng), runif(0.2, 2.0), runif(1.0, 5.0), randn(rng), randn(rng) * 0.5
        Wf(z) = exp(-(z - θ0)^2 / (2ω0^2))
        ys = range(zb0 - 12sqrt(P0), zb0 + 12sqrt(P0); length = 20001)
        dens = exp.(-(ys .- zb0) .^ 2 ./ (2P0))
        ratio = sum(dens .* Wf.(ys .+ a0)) / sum(dens .* Wf.(ys))
        s0, λ0 = ω0^2 + P0, θ0 - zb0
        check(@sprintf("marginal fitness by quadrature = exp(aλ/s − a²/(2s)), a=%.2f λ=%.2f", a0, λ0),
              isapprox(ratio, exp(a0 * λ0 / s0 - a0^2 / (2s0)); rtol = 1e-9))
    end

    say(@sprintf("\n%d checks failed; %.0f s", length(failures), time() - t0))
    say(isempty(failures) ? "PASS" : "FAIL: " * join(failures, "; "))
    write(LOG, String(take!(out)))
    isempty(failures) || exit(1)
end

main()
