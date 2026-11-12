# The critical rate of a moving optimum when drift limits the genetic variance (ROADMAP F8).
#
# `sustainable_rate` (plasticity.jl) asks when the steady lag of a population with a GIVEN
# genetic variance leaves too little fitness to persist. In QTSim's regime that threshold is
# 5–7 times the rate at which half the populations survive (validate/plasticity_rates.log),
# and a constant cost of fitness lowers the half-survival rate 2–7 times more than its
# fitness-margin term says (validate/plasticity_erosion.log). That a fixed-V_G threshold
# overestimates persistence is known (Bürger and Lynch 1995; Kopp and Matuszewski 2014), and
# so is the reason, a feedback the threshold leaves out: the lag lowers mean fitness, so
# fewer adults survive to breed, so drift is stronger and the variance that feeds the
# response is smaller, so the lag grows. Bürger and Lynch (1995, Evolution 49:151–163)
# describe it in words, from simulations, as a "synergistic interaction (positive
# feedback)"; Xu and Shen (2026, bioRxiv preprint, "Eco-evolutionary feedbacks generate
# bistability in population persistence under gradual environmental change") define a
# critical rate as the maximum over steady states with drift-reduced variance, in
# continuous time with logistic growth; Polechová and Barton (2015, PNAS 112:6401–6406) and
# Polechová (2026, PNAS 123:e2604510123) treat the analogue along spatial gradients. Below a
# critical rate the feedback has a steady lag; above it the lag runs away and the
# population dies. A cost of fitness lowers the number of breeders in the same way as the
# lag load, so it lowers the critical rate through the variance as well.
#
# `drift_critical_rate` is a discrete-generation implementation of that known feedback,
# matched to QTSim's life cycle (viability selection, so the breeders, not the census,
# carry it, and it acts at a fixed size too), and tested against the kernel
# (validate/critical_rate.jl). Per generation at a steady lag λ:
#   W̄(λ) = c₀ e^{−c} ω/√s · exp(−λ²/(2s)),  s = ω² + V_G + V_E,  c₀ = 1/(0.4√(2π)), c the log cost;
#   breeders n = N W̄ at a fixed size; under Beverton–Holt the equilibrium n = K(RW̄ − 1)/(R − 1)
#   with N = n/W̄ adults (Lean `fitness_fixed_survivors`);
#   N_e = effective_size(N, n), V_G = the house-of-Gauss variance with that N_e (Hermisson and
#   Wagner 2004 eq. A4; house_of_gauss_vg);
#   the response V_G λ/(V_G + V_s) per generation, V_s = ω² + V_E (steady_lag).
# k(λ) = V_G(λ)λ/(V_G(λ) + V_s) has a maximum, the critical rate k_c, at a lag λ*.
#
# Writing u = λ²/(2s) + c, W̄ = W̄₀e^{−u}, so n, N_e and V_G depend on λ and c only through
# u, and ln k = ln λ + H(u) (s held at its value at λ*). The first-order condition
# 1/λ + H′(u)λ/s = 0 and the envelope theorem give d ln k_c/dc = H′(u*) = −s/λ*²:
#     k_c(c)/k_c(0) = 1 − c/M + o(c) (O(c²) if k_c is twice differentiable; Lean: ratio_linear, powerLaw_linear),   M = λ*²/s = 2 × (the lag load at the critical lag).
# This first-order form is generic, not special to drift: whenever the cost enters only
# through u, the slope is −s/λ*² at the critical lag. A fixed-V_G threshold such as
# `sustainable_rate`, k ∝ √(2s·(L₀ − c)) with L₀ = ln(W̄(0)/W_min), gives
# √(1 − c/L₀) = 1 − c/(2L₀) + O(c²), the same form with M′ = 2L₀ = λ_thr²/s at its
# threshold lag. What drift changes is the scale: it makes the critical lag, and so M, much
# smaller (M = 1.75–1.79 at a fixed size and 0.31–1.09 under Beverton–Holt in
# validate/critical_rate.log, 2.5–7 times below the M′ that `sustainable_rate` implies
# there), so the same cost removes several times more of the critical rate. Over the
# tested range (observed cost ratios k50(c)/k50(0) down to 0.25) a ratio linear in the cost
# fitted better than a square root in 8 of 8 groups, and `sustainable_rate`'s predicted
# ratio was within tolerance in 0 of 24 cells. Symbolic checks:
# validate/symbolic/check_critical.jl; test: validate/critical_rate.jl.

# the house-of-Gauss V_G (house_of_gauss_vg) with `ng` = Ploidy·N_e genes per locus
function _house_of_gauss(p::Params, ng::Real)
    α2, Vs = mutation_effect_variance(p), selection_variance(p)
    nu, μ = rare_units(p)
    b = α2 / (2Vs) + 1 / ng
    nu * 2μ * α2 / (b + sqrt(b^2 + 4μ * α2 / Vs))
end

# the steady state of the feedback at lag λ and log cost c: V_G, breeders, adults, W̄
function _feedback_state(p::Params, λ::Real, c::Real)
    ω, VE, K = initFitSD(p), envVar(p), p.PopSize
    VG = _house_of_gauss(p, p.Ploidy * float(K))
    n = N = W = 0.0
    for _ in 1:100
        s = ω^2 + VG + VE
        W = QTSIM_FITNESS_PEAK * ω / sqrt(s) * exp(-λ^2 / (2s) - c)
        if p.Demography === :fixed
            N, n = float(K), K * W
        else
            R = p.GrowthRate
            R * W <= 1 && return (VG = 0.0, n = 0.0, N = 0.0, Wbar = W)
            n = K * (R * W - 1) / (R - 1)
            N = n / W
        end
        n < p.Ploidy && return (VG = 0.0, n, N, Wbar = W)
        V = _house_of_gauss(p, p.Ploidy * effective_size(N, n))
        done = abs(V - VG) <= 1e-13 * V
        VG = V
        done && break
    end
    (; VG, n, N, Wbar = W)
end

"""
    drift_critical_rate(p; cost = -log(plasticity_cost(p)))

The critical rate of a continuously moving optimum when the genetic variance is set by
mutation, stabilising selection and the drift of the breeders that the lag leaves: the
largest steady response `k(λ) = V_G(λ) λ/(V_G(λ) + V_s)` over lags `λ`, where at lag `λ`
mean fitness is `W̄ = c₀e^{−cost}ω/√s·exp(−λ²/(2s))` (`s = ω² + V_G + V_E`), the breeders are
`n = N W̄` at a fixed size and the Beverton–Holt equilibrium `n = K(RW̄ − 1)/(R − 1)` (adults
`n/W̄`) under `Demography = :beverton_holt`, `N_e` = [`effective_size`](@ref)`(N, n)`, and
`V_G` is the house-of-Gauss variance with that `N_e` ([`house_of_gauss_vg`](@ref)). Above
this rate no steady lag exists: the lag, the loss of breeders and the loss of variance feed
each other until the population dies. `cost` is a log-fitness cost paid by every
individual (by default the cost of plasticity of `p`, `b²/(2ω_b²)`). With plasticity the
rate is divided by `|1 − b/B|` as in [`sustainable_rate`](@ref).

Returns `(; k, λ, VG, n, N, Wbar, M)`: the critical rate, and at the critical lag `λ` the
variance, breeders, adults and mean fitness, and `M = λ²/s`, the cost scale: by the
envelope theorem `k(cost)/k(0) = 1 − cost/M + o(cost)` (`s` held fixed), with `M` twice
the lag load at the critical lag. A fixed-`V_G` threshold such as [`sustainable_rate`](@ref)
has the same first-order form, `√(1 − c/L₀) = 1 − c/(2L₀) + O(c²)`; what drift changes is
the scale: the critical lag, and so `M`, is much smaller (2.5–7 times in
`validate/critical_rate.log`), so a cost lowers the critical rate several times more.

The feedback is not new: Bürger and Lynch (1995, *Evolution* 49:151–163) describe it from
simulations, Xu and Shen (2026, bioRxiv preprint, "Eco-evolutionary feedbacks generate
bistability in population persistence under gradual environmental change") define a
critical rate as the maximum over steady states with drift-reduced variance in continuous
time with logistic growth, and Polechová and Barton (2015) and Polechová (2026) treat the
spatial analogue. This function is a discrete-generation implementation matched to
QTSim's life cycle (breeders `N W̄` or `K(RW̄ − 1)/(R − 1)`, house-of-Gauss `V_G`,
[`effective_size`](@ref)), tested against the kernel; the envelope result for a constant
cost is derived here (`validate/symbolic/check_critical.jl`).

This is a deterministic threshold, not the rate at which half of finite populations survive
a finite time: in QTSim's regime the kernel's half-survival rate within 2000 generations
is 0.29–1.06 times it (lower with more loci; the paper's 710 populations in vignettes 09
and 11, `validate/critical_rate_posthoc.log`), against 0.07–0.56 times
[`sustainable_rate`](@ref); [`half_survival_rate`](@ref) predicts it. What it predicts well
is the RELATIVE effect of a cost (`validate/critical_rate.jl`, pre-specified, fresh seeds:
`k50(c)/k50(0)` within 0.054 of `k(c)/k(0)` in 24 of 24 cells, fixed size and Beverton–Holt
R = 5, 2, 1.5, 100 and 500 loci; over the tested range, observed ratios down to 0.25, a
ratio linear in the cost fits better than a square root in all 8 groups, and
[`sustainable_rate`](@ref)'s ratio is within tolerance in 0 of 24 cells).
The error has a sign: the observed scale `M` is 1.09–1.22 times `λ²/s`.
"""
function drift_critical_rate(p::Params; cost::Real = -log(plasticity_cost(p)))
    p.NoTraits == 1 || throw(ArgumentError("drift_critical_rate is for one trait; NoTraits = $(p.NoTraits)"))
    p.Demography in (:fixed, :beverton_holt) ||
        throw(ArgumentError("drift_critical_rate needs Demography = :fixed or :beverton_holt, got :$(p.Demography)"))
    cost >= 0 || throw(ArgumentError("the cost is a log-fitness cost ≥ 0, got $cost"))
    Vs = selection_variance(p)
    rate(λ) = (st = _feedback_state(p, λ, cost); st.VG * λ / (st.VG + Vs))
    λs = range(0, 12 * sqrt(initFitSD(p)^2 + envVar(p) + _house_of_gauss(p, p.Ploidy * float(p.PopSize))), length = 600)
    ks = rate.(λs)
    i = argmax(ks)
    ks[i] > 0 || return (k = 0.0, λ = NaN, VG = 0.0, n = 0.0, N = 0.0, Wbar = NaN, M = NaN)
    a, b = λs[max(i - 1, 1)], λs[min(i + 1, end)]
    for _ in 1:80                                    # golden-section search on the bracket
        m1, m2 = b - 0.618034 * (b - a), a + 0.618034 * (b - a)
        rate(m1) < rate(m2) ? (a = m1) : (b = m2)
    end
    λ = (a + b) / 2
    st = _feedback_state(p, λ, cost)
    s = initFitSD(p)^2 + st.VG + envVar(p)
    (k = rate(λ) / abs(1 - plasticity_ratio(p)), λ, st.VG, st.n, st.N, st.Wbar, M = λ^2 / s)
end

"""
    half_survival_rate(p; cost = -log(plasticity_cost(p)), horizon = 2000, burnin = 2000,
                       loci = :spectrum, VG = house_of_gauss_vg(p), steps = 12)

The rate of a continuous move of the optimum at which a population of `p` survives
`horizon` generations after the move begins with probability one half, predicted by the
moments lowering with drift ([`run_moments`](@ref) with `erosion = true` and `loci`; the
allele-frequency spectrum by default, which follows the change of `V_G` under the move)
along the deterministic trajectory, plus demographic stochasticity along it
([`extinction_distribution`](@ref), the exact chain under Beverton–Holt, binomial survival
at a fixed size). The population is first run `burnin` generations at a fixed optimum from
genetic variance `VG`. `cost` is a log-fitness cost paid by every individual from the
start; when `p` has no plasticity it is applied as the cost of a reaction norm of slope 1
whose optimum slope is 10¹², whose plastic response is then negligible (1e-12 of the move;
`validate/plasticity_erosion.jl`, H4). The rate is found by bisection on `log k`, `steps`
times, between [`drift_critical_rate`](@ref)`/8` and twice it. Requires Catlab.

Returns `(; k50, kc)`, `kc` being `drift_critical_rate(p; cost).k`. Each evaluation is one
moments run of `burnin + horizon` generations (about 10 s with the spectrum).

`k50 = Inf` if half survive at twice `kc`, 0 if not at `kc/8`.

Accuracy (`validate/critical_rate.jl`, pre-specified, fresh seeds; 32 cells, fixed size and
Beverton–Holt R = 5, 2, 1.5, log costs up to 1): the kernel's half-survival rate over this
prediction is 0.89–1.11 with 100 loci (16 of 16 within 0.75–1.33) and 0.68–1.09 with 500
(10 of 16; it over-predicts under Beverton–Holt, R = 5 without a cost and R = 2 and 1.5 without or with a small cost,
0.68–0.74).
"""
function half_survival_rate(p::Params; cost::Real = -log(plasticity_cost(p)), horizon::Integer = 2000,
                            burnin::Integer = 2000, loci::Symbol = :spectrum, VG::Real = house_of_gauss_vg(p),
                            steps::Integer = 12)
    p.MoveType === :continuous || throw(ArgumentError("half_survival_rate is for MoveType = :continuous"))
    cost >= 0 || throw(ArgumentError("the cost is a log-fitness cost ≥ 0, got $cost"))
    q = Params(p; Sample = 1, NoGenerations = burnin + horizon, MoveFitnessSDGen = burnin)
    if p.Plasticity == 0 && cost > 0
        q = Params(q; Plasticity = 1.0, OptimumSlope = 1e12, PlasticityCostSD = 1 / sqrt(2cost))
    elseif cost != -log(plasticity_cost(p))
        throw(ArgumentError("p has plasticity; its cost is that of its PlasticityCostSD (leave `cost` out)"))
    end
    kc = drift_critical_rate(p; cost).k
    kc > 0 || return (k50 = 0.0, kc)
    psurv(k) = (r = Params(q; MoveFitnessSD = k);
                extinction_distribution(r, run_moments(r; VG, erosion = true, loci); from = burnin).psurvive)
    lo, hi = log(kc / 8), log(2kc)
    psurv(exp(hi)) >= 0.5 && return (k50 = Inf, kc)
    psurv(exp(lo)) < 0.5 && return (k50 = 0.0, kc)
    for _ in 1:steps
        mid = (lo + hi) / 2
        psurv(exp(mid)) >= 0.5 ? (lo = mid) : (hi = mid)
    end
    (k50 = exp((lo + hi) / 2), kc)
end
