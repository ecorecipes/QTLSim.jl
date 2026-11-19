# An individual-based test of the bistability that Xu and Shen (2026, bioRxiv preprint,
# doi:10.64898/2026.08.05.743160) predict below the critical rate of a moving optimum: does
# persistence in the kernel depend on the initial state (N₀, V₀, Δ₀), and where is the
# boundary relative to deterministic predictions?
#
# Design, cells, statistics, predictions and PASS/FAIL rules: validate/bistability_criteria.md
# (committed before any non-pilot run). Findings: validate/bistability.md.
#
#   julia -t 4 --project=vignettes validate/bistability.jl pilot    # one cell, own seeds (disclosed in the criteria)
#   julia -t 4 --project=vignettes validate/bistability.jl          # the experiment; resumable (validate/bistability_raw/)
#   julia -t 4 --project=vignettes validate/bistability.jl analyse  # verdicts and figures from the cache
#
# Writes validate/bistability.log and validate/bistability_figs/*.png.

# 2026-11-19: `ne = :falconer`, `solver = :iteration` and `QTLSim.burnin_start_vg` below keep the effective size, the solver and
# the starting variance this script was run with (the package's defaults are now `:breeders` and `:bisection`; NEWS.md), so that
# it still reproduces its log.
using QTLSim, Statistics, Random, Printf, LinearAlgebra

const HERE = @__DIR__
const RAW = joinpath(HERE, "bistability_raw")
const FIGS = joinpath(HERE, "bistability_figs")
const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s); flush(stdout))
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))
const PEAK = QTLSim.QTSIM_FITNESS_PEAK

# ── design constants (criteria file) ───────────────────────────────────────────────────────
# horizon: alive T generations after the move starts. Per cell, the smallest of 500, 1000, 2000 above the latest death
# of the deterministic prediction P-D (N_e = breeders) on the cell's main grid (criteria, "Pilot and horizon")
const HORIZONS = Dict((:G, 1000, 2.0) => 500, (:G, 1000, 1.5) => 500, (:G, 3000, 2.0) => 500, (:G, 3000, 1.5) => 1000,
                      (:Q, 1000, 2.0) => 2000, (:Q, 1000, 1.5) => 2000, (:Q, 3000, 2.0) => 1000, (:Q, 3000, 1.5) => 2000)
const HORIZON = Ref(500)                        # the horizon of the cell being run
const REPS = 30
const BURN = 500                                # generations at a fixed optimum for the pool populations
const NFRAC = (0.02, 0.05, 0.15, 0.4, 0.9)      # N₀/K
const VFRAC = (0.1, 0.2, 0.35, 0.6, 1.0)        # V₀/V_eq
const RATEF = (0.6, 0.85)                       # k/k_c
const LAGF = (0.0, 0.5, 1.0)                    # Δ₀/λ*
const NULLF = 2.0                               # the above-critical null arm, k/k_c
const SEED0 = 81_000_000_000
const PILOT_SEED0 = 80_000_000_000
# cells in the order they are run; if the budget runs out the last ones are dropped (criteria)
const CELLS = [(:Q, 1000, 2.0), (:G, 1000, 2.0), (:Q, 1000, 1.5), (:G, 1000, 1.5), (:Q, 3000, 2.0), (:G, 3000, 2.0)]
const DROPPED = [(:Q, 3000, 1.5), (:G, 3000, 1.5)]      # dropped before the run, from the pilot's cost (criteria, "Budget")

"The parameters of a cell: the Gaussian arm (:G) or the QTSim-type arm (:Q), carrying capacity K, growth rate R."
arm(a::Symbol, K, R) = a === :G ?
    Params(PopSize = K, Ploidy = 2, NoQTL = 50, RecombPoissonMean = 5.0, CrossoverModel = :all, PhenotypeVar = 1.0, Heritability = 0.5,
           MutationScalar = 0.002, MutationRatePerGamete = 0.05, InitFitnessSD = 5.0, MoveType = :continuous,
           Demography = :beverton_holt, GrowthRate = R) :
    Params(PopSize = K, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0, CrossoverModel = :qtsim, PhenotypeVar = 0.03, Heritability = 0.4,
           MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous,
           Demography = :beverton_holt, GrowthRate = R)

"What defines a cell's grid: the critical rate and lag (N_e = breeders), and the steady state without a move."
function cellinfo(a, K, R)
    p = arm(a, K, R)
    d = drift_critical_rate(p; ne = :breeders, solver = :bisection)
    st = QTLSim._feedback_state(p, 0.0, 0.0; ne = :breeders, solver = :bisection)
    (; p, kc = d.k, λc = d.λ, Veq = st.VG, Neq = st.N, kcf = drift_critical_rate(p; ne = :falconer, solver = :iteration).k)
end
cellname(a, K, R) = @sprintf("%s-K%d-R%g", a, K, R)

# ── predictions ───────────────────────────────────────────────────────────────────────────

"""
P-D: the deterministic discrete-generation recursion with QTLSim's ingredients, in the kernel's order
(the optimum moves, selection, response, variance, density regulation). Returns the generation at
which fewer than two breeders are left, or 0 if the population is alive after `T` generations.
"""
function pd_fate(p::Params, λ0, V0, N0; k, ne::Symbol, T::Integer)
    ω, VE, K, R = initFitSD(p), QTLSim.envVar(p), float(p.PopSize), p.GrowthRate
    Vs = ω^2 + VE; α2 = QTLSim.mutation_effect_variance(p); nu, μ = QTLSim.rare_units(p); Vm = nu * μ * α2
    λ, V, N = float(λ0), float(V0), float(N0)
    for t in 1:T
        λ += k
        s = ω^2 + V + VE
        W = PEAK * ω / sqrt(s) * exp(-λ^2 / (2s))
        n = N * W
        n < 2 && return t
        λ -= V * λ / s
        Ne = effective_size(N, n; scheme = ne)
        V = max(V + Vm - V^2 / (nu * Vs) - V * (α2 / (2Vs) + 1 / (2Ne)), 0.0)
        N = min(K, R * n / (1 + (R - 1) * n / K))
    end
    0
end
const PD_T = 10_000
pd_persists(p, λ0, V0, N0; k, ne) = pd_fate(p, λ0, V0, N0; k, ne, T = PD_T) == 0

"The Xu–Shen parameters of a cell (validate/xu_shen_targets.md, part B): K1 or K2."
function xs_map(p::Params; Kmap::Symbol = :K1)
    ω, VE, R = initFitSD(p), QTLSim.envVar(p), p.GrowthRate
    XuShen(r = log(R) + log(PEAK) - 0.5 * log(1 + VE / ω^2), ω2 = ω^2 + VE, VM = QTLSim.muteVar(p),
           K = p.Ploidy * p.PopSize / (Kmap === :K1 ? log(R) : R - 1), ne = QTLSim.gaussian_units(p))
end
"P-XS: Xu and Shen's model from the mapped state (N ↔ Ploidy·N₀ gene copies); above its own critical rate nothing persists."
function xs_persists(m::XuShen, p::Params, λ0, V0, N0; k)
    xu_shen_equilibrium(m, k) === nothing && return false
    QTLSim._xs_ok(m, λ0, V0, p.Ploidy * N0; v = k)
end

const PREDS = (:pd_breeders, :pd_falconer, :xs_K1, :xs_K2)
const PREDNAME = Dict(:pd_breeders => "P-D, N_e = breeders", :pd_falconer => "P-D, N_e = effective_size", :xs_K1 => "P-XS, K1", :xs_K2 => "P-XS, K2")
function predict(pred::Symbol, p::Params, λ0, V0, N0; k)
    pred === :pd_breeders ? pd_persists(p, λ0, V0, N0; k, ne = :breeders) :
    pred === :pd_falconer ? pd_persists(p, λ0, V0, N0; k, ne = :falconer) :
    pred === :xs_K1 ? xs_persists(xs_map(p), p, λ0, V0, N0; k) : xs_persists(xs_map(p; Kmap = :K2), p, λ0, V0, N0; k)
end

"""
The predicted boundary in V₀ at (N₀, Δ₀): the smallest initial variance in [Veq/1000, 3·Veq] from which the prediction
persists, by bisection on ln V₀. 0 if it persists from Veq/1000, Inf if not from 3·Veq.
"""
function boundary(pred::Symbol, p::Params, Veq, λ0, N0; k)
    ok(V) = predict(pred, p, λ0, V, N0; k)
    lo, hi = Veq / 1000, 3Veq
    ok(lo) && return 0.0
    ok(hi) || return Inf
    a, b = log(lo), log(hi)
    for _ in 1:40
        m = (a + b) / 2
        ok(exp(m)) ? (b = m) : (a = m)
    end
    exp((a + b) / 2)
end

# ── the kernel runs ───────────────────────────────────────────────────────────────────────

"Pool population `r` of a cell: sampled at the stationary spectrum, then `BURN` generations at a fixed optimum under the cell's own demography."
function pool_population(p::Params, seed)
    rng = Xoshiro(seed)
    q = Params(p; NoGenerations = BURN, MoveFitnessSDGen = -1)
    w = World(q)
    G0 = equilibrium_genomes(Params(q; Demography = :fixed), rng)
    setgenomes!(w, G0)
    init_env!(w, rng)
    for _ in 1:BURN
        generation!(w, rng)
    end
    w.extinct && error("a pool population died during the burn-in")
    develop!(w)
    (G = genomes(w), VG = var(view(w.genotype, 1:w.N)), N = w.N, lag = w.θ - mean(view(w.genotype, 1:w.N)))
end

"""
A kernel population started at (N₀, V₀, Δ₀): `N0` individuals drawn without replacement from the pool population `G`, all
allele values scaled by √(V₀/V) (V the sample's genetic variance), the optimum set `Δ0` above the mean genotype. The move
starts in the first generation. Returns the `World` and the realised (N, V_G, lag).
"""
function start_world(p::Params, G::Array{Float64,3}, N0::Integer, V0::Real, Δ0::Real, rng)
    idx = randperm(rng, size(G, 3))[1:N0]
    H = G[:, :, idx]
    g = vec(sum(H; dims = (1, 2)))
    v = var(g)
    v > 0 || return nothing
    H .*= sqrt(V0 / v)
    w = World(p)
    setgenomes!(w, H)
    w.N = N0
    init_env!(w, rng)
    develop!(w)
    gbar = mean(view(w.genotype, 1:N0))
    w.θ = gbar + Δ0
    (w, (N = w.N, VG = var(view(w.genotype, 1:N0)), lag = w.θ - gbar))
end

"One run: alive after `T` generations of a move at rate `k`? Returns (alive, generation of death or T, final N, V_G, lag)."
function one_run(p::Params, G, N0, V0, Δ0, k, seed)
    rng = Xoshiro(seed)
    T = HORIZON[]
    q = Params(p; NoGenerations = T, MoveFitnessSDGen = 0, MoveFitnessSD = k)
    sw = start_world(q, G, N0, V0, Δ0, rng)
    sw === nothing && return (alive = false, t = 0, N = 0, VG = NaN, lag = NaN, errN = 0.0, errV = NaN, errL = NaN)
    w, real = sw
    t = T
    for g in 1:T
        generation!(w, rng)
        if w.extinct
            t = g
            break
        end
    end
    if w.extinct
        return (alive = false, t, N = 0, VG = NaN, lag = NaN, errN = real.N - N0, errV = real.VG / V0 - 1, errL = real.lag - Δ0)
    end
    develop!(w)
    (alive = true, t, N = w.N, VG = var(view(w.genotype, 1:w.N)), lag = w.θ - mean(view(w.genotype, 1:w.N)),
     errN = float(real.N - N0), errV = real.VG / V0 - 1, errL = real.lag - Δ0)
end

"The states of a cell: the main grid (2 rates × 3 lags × 5 × 5) and the two null arms (no move; 2·k_c), lag 0."
function states(ci)
    v = NamedTuple[]
    for (ir, rf) in enumerate((RATEF..., 0.0, NULLF)), (il, lf) in enumerate(LAGF)
        (rf == 0.0 || rf == NULLF) && il > 1 && continue
        for (iN, nf) in enumerate(NFRAC), (iV, vf) in enumerate(VFRAC)
            push!(v, (; ir, il, iN, iV, rf, lf, k = rf * ci.kc, Δ0 = lf * ci.λc, N0 = round(Int, nf * ci.p.PopSize), V0 = vf * ci.Veq))
        end
    end
    v
end

function run_cell(a, K, R; seed0 = SEED0, tag = "")
    ci = cellinfo(a, K, R)
    name = cellname(a, K, R)
    icell = findfirst(==((a, K, R)), CELLS)
    HORIZON[] = HORIZONS[(a, K, R)]
    mkpath(RAW)
    t0 = time()
    pools = Vector{Any}(undef, REPS)
    Threads.@threads :greedy for r in 1:REPS
        pools[r] = pool_population(ci.p, seed0 + 10^8 * icell + r)
    end
    tp = time() - t0
    jobs = [(; s..., st = i) for (i, s) in enumerate(states(ci))]
    rows = sweep(jobs; cache = joinpath(RAW, "$(name)$(tag).csv"), id = j -> "$(j.ir)-$(j.il)-$(j.iN)-$(j.iV)",
                 cost = j -> j.rf == NULLF ? 0.2 : 1.0, progress = 50) do j
        rs = [one_run(ci.p, pools[r].G, min(j.N0, pools[r].N), j.V0, j.Δ0, j.k, seed0 + 10^8 * icell + 10^4 * j.st + 100 + r) for r in 1:REPS]
        al = [r for r in rs if r.alive]
        (; j.ir, j.il, j.iN, j.iV, j.rf, j.lf, j.k, j.Δ0, j.N0, j.V0, alive = length(al), reps = REPS,
         tdead = isempty(al) || length(al) < REPS ? mean(r.t for r in rs if !r.alive) : NaN,
         Nend = isempty(al) ? NaN : mean(r.N for r in al), Vend = isempty(al) ? NaN : mean(r.VG for r in al),
         lagend = isempty(al) ? NaN : mean(r.lag for r in al),
         maxerrN = maximum(abs(r.errN) for r in rs), maxerrV = maximum(x -> isnan(x) ? 0.0 : abs(x), (r.errV for r in rs)),
         maxerrL = maximum(x -> isnan(x) ? 0.0 : abs(x), (r.errL for r in rs)))
    end
    pool = (VG = mean(x.VG for x in pools), sdVG = std([x.VG for x in pools]), N = mean(x.N for x in pools), lag = mean(x.lag for x in pools), secs = tp)
    open(joinpath(RAW, "$(name)$(tag).pool"), "w") do io
        @printf(io, "%.10g %.10g %.10g %.10g\n", pool.VG, pool.sdVG, pool.N, pool.lag)
    end
    (; ci, rows, pool, secs = time() - t0)
end

include(joinpath(HERE, "bistability_analysis.jl"))

function main()
    mode = isempty(ARGS) ? "run" : ARGS[1]
    if mode == "pilot"
        return pilot()
    elseif mode == "run"
        budget = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : Inf        # wall-clock minutes for this session
        t0 = time()
        for c in CELLS
            isfile(joinpath(RAW, cellname(c...) * ".done")) && continue
            (time() - t0) / 60 > budget && (println("budget reached before ", cellname(c...)); break)
            r = run_cell(c...)
            write(joinpath(RAW, cellname(c...) * ".done"), @sprintf("%.1f\n", r.secs))
            @printf("%s done in %.1f min\n", cellname(c...), r.secs / 60)
        end
    end
    analyse()
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
