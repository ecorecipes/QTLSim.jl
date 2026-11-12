# Why does a cost of plasticity lower the sustainable rate of change more than the theory
# says? A pre-specified factorial experiment with the kernel, following
# validate/plasticity_rates.jl (ROADMAP F2).
#
# Everything from here to "END OF PRE-SPECIFICATION" was written before the first run of
# the main experiment. A pilot (`pilot` argument, 8 replicates, its own seeds, logged in
# plasticity_erosion_pilot.log) was run first, only to place each arm's grid of rates; its
# half-survival rates are the constants K50_PILOT below.
#
# ── Background ────────────────────────────────────────────────────────────────────────
# In plasticity_rates.log, criterion P2 (with a cost of plasticity ω_b = 0.5, the
# half-survival rate k50 within 25% of the prediction of sustainable_rate) failed in 4 of 12
# cells: at a fixed size with b/B = 0.75 observed/predicted was 1.39 (L100) and 1.26 (L500);
# under Beverton–Holt (R = 2) with b/B = 0.5 every run died at every rate of the grid. The
# explanation offered there, inferred from the per-cell means and not tested: the cost
# multiplies every fitness by c = exp(−b²/(2ω_b²)) (0.61 at b = 0.5, 0.32 at b = 0.75), so
# fewer individuals survive selection, drift is stronger and V_G before the move is lower
# (by 33–58% at a fixed size, about 3× under Beverton–Holt, where N at the move also falls
# to about 340).
#
# A post-hoc look at the same CSV (the `posthoc` argument; plasticity_erosion_posthoc.log),
# done before this design was written, showed:
#  - P2's prediction already used each cell's measured V_G at the move. At a fixed size the
#    prediction is too PESSIMISTIC (observed k50 above it), not too optimistic: it is
#    nearly proportional to V_G, and the observed k50 falls by less than V_G does.
#  - Across the 40 replicates of a cell without a cost, V_G at the move does not predict
#    survival at all: the elasticity of k50 with respect to a replicate's V_G at the move
#    is 0.01, 0.04, −0.02, 0.02 (95% intervals within −0.04 to 0.07).
#  - The cost lowers V_G during the move too (generations 500–2000), by about as much as or
#    more than before it.
# So "the cost erodes V_G before the move" is not by itself an account of the failures:
# the cost acts before the move (V_G, and N under Beverton–Holt) and after it (the fitness
# margin, and fewer parents during the move). This experiment separates the two.
#
# ── Design ────────────────────────────────────────────────────────────────────────────
# The regime, parameter sets (L100, L500), demographies (fixed size; Beverton–Holt R = 2,
# K = 1000), burn-in (10,000 generations at a fixed optimum), continuous move, cap (2000
# generations after the move) and cost (ω_b = 0.5, B = 1, perfect cue) are those of
# plasticity_rates.jl. Cells: the two failed slopes and one passing slope per demography:
#   fixed size: b/B = 0.5 (passed) and 0.75 (failed), for L100 and L500;
#   Beverton–Holt: b/B = 0.25 (passed) and 0.5 (failed), for L100 and L500.
# Arms, a 2 × 2 factorial in when the cost applies, plus a control:
#   0   no cost before or after the move (b = 0; the same in distribution as any b without a
#       cost at the same effective rate, plasticity_rates.jl's P1);
#   C   no cost before the move, the cost from the move on (the cost's post-move effect);
#   B   the cost before the move, removed at the move (its pre-move effect: V_G, and N under
#       Beverton–Holt, at the move);
#   A   the cost before and after (plasticity_rates.jl's cost cells);
#   A′  as A, but without the plastic correspondence: the same slope b and so the same cost
#       c, with OptimumSlope = 1e12 (so b/B ≈ 0), at the effective rate (1 − b/B)k itself.
# Arms 0 and C share their burn-ins (no cost), one set per (parameter set, demography);
# A, B and A′ share theirs (with the cost), one set per cell. Paired within each family.
# Rates: every arm's rates are in EFFECTIVE units (1 − b/B)·k (A′: k), eleven per arm,
# K50_PILOT[arm] × {0.6, 0.8, 0.88, 0.94, 0.97, 1, 1.03, 1.06, 1.12, 1.25, 1.6}, with the
# actual rate k = effective/(1 − b/B). Replicates: 40 burn-ins per family. Fresh seeds
# (SEED0 = 7_000_000_000; the pilot 8_000_000_000; run 1 of plasticity_rates.jl used
# 20_261_110-based seeds).
#
# ── Measures ──────────────────────────────────────────────────────────────────────────
# k50 per arm from a logistic regression of survival (alive 2000 generations after the
# move) on ln(effective rate), as in plasticity_rates.jl. D_X = ln(k50_0/k50_X), the log
# reduction of the half-survival rate by arm X. The factorial decomposition of D_A:
#   M_pre  = ½[D_B + (D_A − D_C)]   the main effect of the cost before the move,
#   M_post = ½[D_C + (D_A − D_B)]   the main effect of the cost after the move,
#   I      = D_A − D_B − D_C        the interaction;   M_pre + M_post = D_A,
# and share_pre = M_pre/D_A. Theory (sustainable_rate, which is ∝ V_G/s·√(2s ln(W̄(0)/W_min))):
#   T_fit = ln[S(no cost; V_G0)/S(cost; V_G0)], the fitness-margin effect at a fixed V_G;
#   T_VG  = ln[S(no cost; V_G0)/S(no cost; V_GB)], the effect of the pre-move V_G loss alone;
# V_G0 and V_GB the mean V_G at the move of the no-cost and cost burn-ins. Also reported:
# V_G and N at the move, and their means over the first 200 generations after the move.
# Intervals: 95% percentile intervals from 1000 bootstrap resamples, the two families of
# burn-ins resampled independently, each replicate with all its branches.
#
# ── Criteria ──────────────────────────────────────────────────────────────────────────
# H1 (the explanation: the cost acts mainly through the state it leaves at the move):
#     share_pre ≥ 0.5. 8 cells; the four failed cells decide the verdict.
# H2 (the cost's post-move effect is the theory's fitness margin): |M_post − T_fit| ≤ 0.10.
#     8 cells.
# H3 (the pre-move V_G loss acts as the theory says, k50 nearly ∝ V_G): D_B/T_VG within
#     0.75–1.25. 8 cells. (The post hoc look predicts a FAIL: k50 falls less than V_G.)
# H4 (no role for the plastic correspondence): k50_A′/k50_A within 0.95–1.05. 8 cells.
#     This is an identity (with a perfect cue z − θ = g + e − (1 − b/B)θ, and the cost does not
#     depend on θ), so it checks the implementation.
# Verdict on the explanation: SUPPORTED if H1 passes in at least 3 of the 4 failed cells;
# REFUTED if share_pre's interval lies wholly below 0.5 in at least 3 of them; otherwise
# MIXED. A cell whose k50 of arm A (or any arm) is not bracketed by its grid gives NaN and
# fails every criterion that uses it.
# Power: with 40 replicates and eleven rates around k50, run 1's k50 intervals were about
# ±1.5%, so each D_X has a standard error of about 0.01–0.02 and M_pre about 0.02; at a
# fixed size with b/B = 0.75, D_A ≈ 0.64, so share_pre is known to about ±0.05, and in the
# weakest cell (fixed, b/B = 0.5, D_A ≈ 0.31) to about ±0.1.
# Failures stay in the log. If the model or the analysis is changed after a failed run, the
# change is stated here and confirmed with fresh seeds.
#
# END OF PRE-SPECIFICATION
#
# Notes added after the pilot, before the main run (nothing above changed): the pilot
# (plasticity_erosion_pilot.log; 8 replicates, a coarse grid, so its k50 are rough) already
# showed arm B close to arm 0 and arm C close to arm A in every cell, i.e. a small pre-move
# effect. The criteria were written before the pilot and are kept as they are.
#
# Run from the package root, at most 3 threads (pilot about 0.4 CPU-hours, main about
# 2 CPU-hours):
#   julia -t 3 --project=vignettes validate/plasticity_erosion.jl posthoc   # from plasticity_rates.csv
#   julia -t 3 --project=vignettes validate/plasticity_erosion.jl pilot     # places the grids
#   julia -t 3 --project=vignettes validate/plasticity_erosion.jl [SEED0]   # the experiment; resumes
# Writes validate/plasticity_erosion.csv (via sweep; one row per burn-in with all its
# branches) and prints the analysis, logged in validate/plasticity_erosion.log.

using QTLSim, Statistics, Random, Printf

const BURN = 10_000
const CAP = 2_000
const NREP = 40
const NPILOT = 8
const GRID = (0.6, 0.8, 0.88, 0.94, 0.97, 1.0, 1.03, 1.06, 1.12, 1.25, 1.6)
const PILOT_GRID = Tuple(exp.(range(log(0.03), log(1.5); length = 15)))
const COST = 0.5
const BIG_SLOPE = 1e12
const SETS = (L100 = Params(Sample = 10, PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0,
                            PhenotypeVar = 0.03, Heritability = 0.4, MutationScalar = 0.005,
                            MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous),
              L500 = Params(Sample = 10, PopSize = 1000, Ploidy = 2, NoQTL = 500, RecombPoissonMean = 3.0,
                            PhenotypeVar = 0.05, Heritability = 0.3, MutationScalar = 0.008,
                            MutationRatePerGamete = 1e-4, InitFitnessSD = 3.0, MoveType = :continuous))
const DEMS = (fixed = (Demography = :fixed,), bh = (Demography = :beverton_holt, GrowthRate = 2.0))
const SLOPES = Dict(:fixed => (0.5, 0.75), :bh => (0.25, 0.5))
const FAILED = Set([(:L100, :fixed, 0.75), (:L500, :fixed, 0.75), (:L100, :bh, 0.5), (:L500, :bh, 0.5)])
# run 1's half-survival rates without plasticity (plasticity_rates.log), the pilot's reference
const K_REF = Dict((:L100, :fixed) => 0.02262, (:L100, :bh) => 0.01455, (:L500, :fixed) => 0.05113, (:L500, :bh) => 0.03751)
const ARMS = (:zero, :C, :B, :A, :Aprime)
const ARMID = Dict(:zero => 1, :C => 2, :B => 3, :A => 4, :Aprime => 5)
# effective half-survival rates of the pilot, per (set, dem, b, arm); b = 0 for arm 0. Filled in from the pilot.
const K50_PILOT = Dict{Tuple{Symbol,Symbol,Float64,Symbol},Float64}(
    (:L100, :bh, 0.00, :zero) => 0.014346,
    (:L100, :bh, 0.25, :A) => 0.012155,
    (:L100, :bh, 0.25, :B) => 0.014346,
    (:L100, :bh, 0.25, :C) => 0.012289,
    (:L100, :bh, 0.50, :A) => 0.0026841,
    (:L100, :bh, 0.50, :B) => 0.012815,
    (:L100, :bh, 0.50, :C) => 0.0028834,
    (:L100, :fixed, 0.00, :zero) => 0.022303,
    (:L100, :fixed, 0.50, :A) => 0.015062,
    (:L100, :fixed, 0.50, :B) => 0.022303,
    (:L100, :fixed, 0.50, :C) => 0.016866,
    (:L100, :fixed, 0.75, :A) => 0.011176,
    (:L100, :fixed, 0.75, :B) => 0.019923,
    (:L100, :fixed, 0.75, :C) => 0.014295,
    (:L500, :bh, 0.00, :zero) => 0.036984,
    (:L500, :bh, 0.25, :A) => 0.027968,
    (:L500, :bh, 0.25, :B) => 0.036984,
    (:L500, :bh, 0.25, :C) => 0.031681,
    (:L500, :bh, 0.50, :A) => 0.0066673,
    (:L500, :bh, 0.50, :B) => 0.033038,
    (:L500, :bh, 0.50, :C) => 0.0079571,
    (:L500, :fixed, 0.00, :zero) => 0.050413,
    (:L500, :fixed, 0.50, :A) => 0.038123,
    (:L500, :fixed, 0.50, :B) => 0.050413,
    (:L500, :fixed, 0.50, :C) => 0.038123,
    (:L500, :fixed, 0.75, :A) => 0.025464,
    (:L500, :fixed, 0.75, :B) => 0.050413,
    (:L500, :fixed, 0.75, :C) => 0.028829,)

cells() = [(set, dem, b) for set in keys(SETS) for dem in keys(DEMS) for b in SLOPES[dem]]
basep(set, dem) = Params(SETS[set]; DEMS[dem]..., NoGenerations = BURN)

"Burn in a population of `p` for BURN generations; the World, or `nothing` if it died."
function burnin(p::Params, seed)
    w = World(p); rng = MersenneTwister(seed)
    init_env!(w, rng)
    for it in 0:BURN-1
        w.gen = it
        generation!(w, rng)
        w.extinct && return nothing
    end
    w
end

"""
Run a copy of `w0` for up to CAP generations with the optimum moving at `k`, slope `b`,
cost width `ωb` (Inf: none) and optimum slope `B`. Survival, generations to extinction,
and mean V_G and N over the first 200 generations after the move.
"""
function branch(w0::World, k, b, ωb, B, seed)
    w = deepcopy(w0)
    w.p = Params(w0.p; MoveFitnessSDGen = BURN, MoveFitnessSD = k, Plasticity = b, PlasticityCostSD = ωb,
                 OptimumSlope = B, NoGenerations = BURN + CAP)
    rng = MersenneTwister(seed)
    rec = QTLSim.records()
    text = CAP
    for it in BURN:BURN+CAP-1
        w.gen = it
        generation!(w, rng; rec)
        w.extinct && (text = it - BURN; break)
    end
    early = findall(g -> g < BURN + 200, rec.gen)
    (alive = !w.extinct, text, vge = isempty(early) ? NaN : mean(rec.geno_var[early]),
     ne = isempty(early) ? NaN : mean(rec.popsize[early]))
end

genvar(w::World) = (develop!(w); var(view(w.genotype, 1:w.N)))

# the arm's effective grid: the pilot's 15 rates around run 1's k50, or the main grid around the pilot's k50
effgrid(set, dem, b, arm, pilot) = pilot ? [K_REF[(set, dem)] * g for g in PILOT_GRID] :
                                   [K50_PILOT[(set, dem, arm == :zero ? 0.0 : b, arm == :Aprime ? :A : arm)] * g for g in GRID]

"""
One job: a burn-in and all its branches. `fam = :nocost` (arms 0 and C for each slope of
the demography) or `:cost` (arms A, B, A′ for slope `b`). Returns one row, the branches
encoded as `arm:b:i:keff:alive:text:vge:ne` joined by `;`.
"""
function job(j, seed0, pilot)
    (; fam, set, dem, b, rep) = j
    cost = fam === :cost
    p = cost ? Params(basep(set, dem); Plasticity = b, PlasticityCostSD = COST) : basep(set, dem)
    famid = 10 * findfirst(==((set, dem)), [(s, d) for s in keys(SETS) for d in keys(DEMS)]) + (cost ? round(Int, 4b) : 0)
    seed = seed0 + 1_000_000 * famid + rep
    w0 = burnin(p, seed)
    vg0, n0 = w0 === nothing ? (NaN, 0) : (genvar(w0), w0.N)
    specs = cost ? [(:A, b, COST, 1.0), (:B, b, Inf, 1.0), (:Aprime, b, COST, BIG_SLOPE)] :
                   [(:zero, 0.0, Inf, 1.0); [(:C, bb, COST, 1.0) for bb in SLOPES[dem]]]
    out = String[]
    for (arm, bb, ωb, B) in specs, (i, ke) in enumerate(effgrid(set, dem, bb, arm, pilot))
        k = ke / (1 - bb / B)
        bseed = seed + 1000 * (100 * ARMID[arm] + 10 * round(Int, 4bb) + i)
        r = w0 === nothing ? (alive = false, text = 0, vge = NaN, ne = NaN) : branch(w0, k, bb, ωb, B, bseed)
        push!(out, join((arm, bb, i, ke, r.alive, r.text, r.vge, r.ne), ':'))
    end
    (; fam, set, dem, b, rep, vg0, n0, runs = join(out, ';'))
end

function jobs(nrep)
    js = NamedTuple[]
    for set in keys(SETS), dem in keys(DEMS), rep in 1:nrep
        push!(js, (fam = :nocost, set, dem, b = 0.0, rep))
        for b in SLOPES[dem]
            push!(js, (fam = :cost, set, dem, b, rep))
        end
    end
    js
end

"Rows of the sweep → one NamedTuple per branch."
function branches(rows)
    out = NamedTuple[]
    for r in rows, s in split(r.runs, ';')
        f = split(s, ':')
        push!(out, (fam = Symbol(r.fam), set = Symbol(r.set), dem = Symbol(r.dem), famb = Float64(r.b), rep = r.rep,
                    vg0 = Float64(r.vg0), n0 = r.n0, arm = Symbol(f[1]), b = parse(Float64, f[2]), i = parse(Int, f[3]),
                    keff = parse(Float64, f[4]), alive = f[5] == "true", text = parse(Int, f[6]),
                    vge = parse(Float64, f[7]), ne = parse(Float64, f[8])))
    end
    out
end

# --- analysis ---------------------------------------------------------------------------

"Maximum-likelihood logistic regression of `y` (Bool) on `x`: (a, c) of logit P = a + c x."
function logistic(x::AbstractVector, y::AbstractVector)
    a, c = 0.0, 0.0
    xm = mean(x); xc = x .- xm
    for _ in 1:100
        η = clamp.(a .+ c .* xc, -30, 30)
        p = 1 ./ (1 .+ exp.(-η))
        g1, g2 = sum(y .- p), sum((y .- p) .* xc)
        wt = p .* (1 .- p) .+ 1e-9
        h11, h12, h22 = sum(wt), sum(wt .* xc), sum(wt .* xc .^ 2)
        det = h11 * h22 - h12^2
        da, dc = (h22 * g1 - h12 * g2) / det, (h11 * g2 - h12 * g1) / det
        s = min(1.0, 5 / max(abs(da), abs(dc), 1e-12))
        a += s * da; c += s * dc
        max(abs(da), abs(dc)) < 1e-10 && break
    end
    (a - c * xm, c)
end

"Half-survival effective rate of branches `rs`; NaN when survival does not fall with the rate."
function k50(rs)
    isempty(rs) && return NaN
    y = [r.alive for r in rs]
    (all(y) || !any(y)) && return NaN
    a, c = logistic([log(r.keff) for r in rs], y)
    c < 0 ? exp(-a / c) : NaN
end

meanor(itr) = (v = collect(itr); isempty(v) ? NaN : mean(v))
quant(v, q) = (w = sort(filter(isfinite, v)); isempty(w) ? NaN : w[clamp(round(Int, q * length(w) + 0.5), 1, length(w))])
ci(v) = (quant(v, 0.025), quant(v, 0.975))
S(p, VG) = sustainable_rate(p; VG) * (1 - plasticity_ratio(p))     # in effective units

"The arms' k50 for one cell, from the branches of the given replicates of each family."
function cellk50(byrep0, byrepc, reps0, repsc, b)
    nb = reduce(vcat, (byrep0[r] for r in reps0); init = NamedTuple[])
    cb = reduce(vcat, (byrepc[r] for r in repsc); init = NamedTuple[])
    Dict(:zero => k50([x for x in nb if x.arm === :zero]),
         :C => k50([x for x in nb if x.arm === :C && x.b == b]),
         :B => k50([x for x in cb if x.arm === :B]),
         :A => k50([x for x in cb if x.arm === :A]),
         :Aprime => k50([x for x in cb if x.arm === :Aprime]))
end

function decompose(k)
    D = Dict(a => log(k[:zero] / k[a]) for a in (:A, :B, :C))
    Mpre = (D[:B] + D[:A] - D[:C]) / 2
    Mpost = (D[:C] + D[:A] - D[:B]) / 2
    (; DA = D[:A], DB = D[:B], DC = D[:C], Mpre, Mpost, I = D[:A] - D[:B] - D[:C], share = Mpre / D[:A],
     Aratio = k[:Aprime] / k[:A])
end

function analyse(rows, seed0)
    br = branches(rows)
    rng = MersenneTwister(seed0 + 99)
    pass = Dict(h => [0, 0] for h in (:H1, :H2, :H3, :H4))
    verdict(h, ok) = (pass[h][1] += ok; pass[h][2] += 1; ok ? "PASS" : "FAIL")
    failedcells = Dict{Any,Tuple{Float64,Float64,Float64}}()
    for (set, dem, b) in cells()
        p0 = basep(set, dem)
        pc = Params(p0; Plasticity = b, PlasticityCostSD = COST)
        nb = [x for x in br if x.fam === :nocost && x.set === set && x.dem === dem]
        cb = [x for x in br if x.fam === :cost && x.set === set && x.dem === dem && x.famb == b]
        (isempty(nb) || isempty(cb)) && continue
        g(v) = (d = Dict{Int,Vector{NamedTuple}}(); for x in v; push!(get!(d, x.rep, NamedTuple[]), x); end; d)
        B0, BC = g(nb), g(cb)
        r0, rc = sort(collect(keys(B0))), sort(collect(keys(BC)))
        k = cellk50(B0, BC, r0, rc, b)
        d = decompose(k)
        boot = [decompose(cellk50(B0, BC, rand(rng, r0, length(r0)), rand(rng, rc, length(rc)), b)) for _ in 1:1000]
        bci(f) = ci([f(x) for x in boot])
        vg0 = meanor(x.vg0 for x in nb if x.arm === :zero && x.i == 1 && isfinite(x.vg0))
        vgB = meanor(x.vg0 for x in cb if x.arm === :A && x.i == 1 && isfinite(x.vg0))
        n0 = meanor(x.n0 for x in nb if x.arm === :zero && x.i == 1)
        nB = meanor(x.n0 for x in cb if x.arm === :A && x.i == 1)
        dead = mean(x.n0 == 0 for x in cb if x.arm === :A && x.i == 1)
        Tfit = log(S(p0, vg0) / S(pc, vg0))
        TVG = log(S(p0, vg0) / S(p0, vgB))
        @printf("\n══ %s, %s, b/B = %.2f%s ══\n", set, dem, b, (set, dem, b) in FAILED ? " (failed P2 in plasticity_rates.log)" : "")
        @printf("  at the move: V_G %.4g without the cost, %.4g with it (×%.2f); N %.0f and %.0f; %.0f%% of cost burn-ins died\n",
                vg0, vgB, vgB / vg0, n0, nB, 100dead)
        println("  first 200 generations after the move, mean over the arm's runs (V_G, N) by arm:")
        for a in ARMS
            xs = a === :zero ? [x for x in nb if x.arm === :zero] : a === :C ? [x for x in nb if x.arm === :C && x.b == b] :
                 [x for x in cb if x.arm === a]
            ks = sort(unique(x.keff for x in xs))
            surv = [mean(x.alive for x in xs if x.keff == kk) for kk in ks]
            @printf("    %-6s V_G %.4g, N %.0f; survival by effective rate %s\n", a, meanor(x.vge for x in xs if isfinite(x.vge)),
                    meanor(x.ne for x in xs if isfinite(x.ne)), join((@sprintf("%.2f", s) for s in surv), " "))
            @printf("           effective rates %s\n", join((@sprintf("%.4g", kk) for kk in ks), " "))
        end
        for a in ARMS
            @printf("  k50 %-6s %.5g [%.5g, %.5g]\n", a, k[a], ci([cellk50(B0, BC, rand(rng, r0, length(r0)), rand(rng, rc, length(rc)), b)[a] for _ in 1:200])...)
        end
        @printf("  D_A %.3f [%.3f, %.3f], D_B %.3f [%.3f, %.3f], D_C %.3f [%.3f, %.3f]\n", d.DA, bci(x -> x.DA)..., d.DB, bci(x -> x.DB)..., d.DC, bci(x -> x.DC)...)
        @printf("  M_pre %.3f [%.3f, %.3f], M_post %.3f [%.3f, %.3f], interaction %.3f [%.3f, %.3f]\n", d.Mpre, bci(x -> x.Mpre)...,
                d.Mpost, bci(x -> x.Mpost)..., d.I, bci(x -> x.I)...)
        sh = bci(x -> x.share)
        @printf("  H1 share_pre = M_pre/D_A %.3f [%.3f, %.3f]  %s\n", d.share, sh..., verdict(:H1, isfinite(d.share) && d.share >= 0.5))
        @printf("  H2 M_post %.3f vs the theory's fitness margin T_fit %.3f: difference %.3f  %s\n", d.Mpost, Tfit, d.Mpost - Tfit,
                verdict(:H2, isfinite(d.Mpost) && abs(d.Mpost - Tfit) <= 0.10))
        @printf("  H3 D_B %.3f vs the theory's V_G effect T_VG %.3f: ratio %.3f [%.3f, %.3f]  %s\n", d.DB, TVG, d.DB / TVG,
                bci(x -> x.DB / TVG)..., verdict(:H3, isfinite(d.DB) && 0.75 <= d.DB / TVG <= 1.25))
        @printf("  H4 k50_A′/k50_A %.3f [%.3f, %.3f]  %s\n", d.Aratio, bci(x -> x.Aratio)...,
                verdict(:H4, isfinite(d.Aratio) && 0.95 <= d.Aratio <= 1.05))
        @printf("  descriptive: D_A vs the theory's total effect T_fit + T_VG = %.3f (the P2 prediction, in logs)\n", Tfit + TVG)
        (set, dem, b) in FAILED && (failedcells[(set, dem, b)] = (d.share, sh...))
    end
    println()
    for h in (:H1, :H2, :H3, :H4)
        @printf("%s: %d of %d cells pass\n", h, pass[h]...)
    end
    nsup = count(v -> isfinite(v[1]) && v[1] >= 0.5, values(failedcells))
    nref = count(v -> isfinite(v[3]) && v[3] < 0.5, values(failedcells))
    v = nsup >= 3 ? "SUPPORTED" : nref >= 3 ? "REFUTED" : "MIXED"
    @printf("\nVerdict on the explanation (the cost acts mainly through the state at the move): %s (H1 passes in %d of the %d failed cells; interval wholly below 0.5 in %d)\n",
            v, nsup, length(failedcells), nref)
    all(h -> pass[h][1] == pass[h][2], keys(pass))
end

# --- the post-hoc look at plasticity_rates.csv (no new runs) -----------------------------

function logistic2(X, y)
    β = zeros(size(X, 2))
    for _ in 1:200
        η = clamp.(X * β, -30, 30); p = 1 ./ (1 .+ exp.(-η))
        gr = X' * (y .- p); H = X' * (X .* (p .* (1 .- p) .+ 1e-9))
        d = H \ gr; s = min(1.0, 5 / maximum(abs.(d))); β .+= s .* d
        maximum(abs.(d)) < 1e-10 && break
    end
    β
end

function posthoc()
    println("POST HOC (not pre-specified): plasticity_rates.csv (run 1 of plasticity_rates.jl), no new runs")
    path = joinpath(@__DIR__, "plasticity_rates.csv")
    rows = NamedTuple[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ',')
        length(f) == 13 || continue
        push!(rows, (set = Symbol(f[1]), dem = Symbol(f[2]), cost = f[3] == "true", b = parse(Float64, f[4]),
                     rep = parse(Int, f[5]), keff = parse(Float64, f[7]), alive = f[8] != "true",
                     vg0 = parse(Float64, f[10]), n0 = parse(Int, f[11]), vgmove = parse(Float64, f[13])))
    end
    rng = MersenneTwister(1)
    for set in keys(SETS), dem in keys(DEMS)
        p0 = Params(SETS[set]; DEMS[dem]...)
        nc = [r for r in rows if r.set == set && r.dem == dem && !r.cost]
        base = [r for r in nc if r.b == 0]
        vg0 = mean(r.vg0 for r in base)
        kb = k50(base)
        fit(rr) = (β = logistic2(hcat(ones(length(rr)), [log(r.keff) for r in rr], [log(r.vg0) for r in rr]), [r.alive for r in rr]); -β[3] / β[2])
        e = fit(nc)
        reps = unique(r.rep for r in nc)
        es = [fit(reduce(vcat, ([r for r in nc if r.rep == i] for i in rand(rng, reps, length(reps))))) for _ in 1:300]
        @printf("\n%s %s: V_G at the move %.4f (replicates %.4f–%.4f), N %.0f; k50 (b = 0) %.5f\n", set, dem, vg0,
                extrema(r.vg0 for r in base)..., mean(r.n0 for r in base), kb)
        @printf("  elasticity of k50 with respect to a replicate's V_G at the move (no-cost runs, all slopes, by effective rate): %.3f [%.3f, %.3f]\n", e, ci(es)...)
        for b in (0.25, 0.5, 0.75)
            rs = [r for r in rows if r.set == set && r.dem == dem && r.cost && r.b == b]
            alive0 = [r for r in rs if r.n0 > 0]
            isempty(alive0) && (@printf("  b/B = %.2f cost: every burn-in died\n", b); continue)
            vgc = mean(r.vg0 for r in alive0)
            pc = Params(p0; Plasticity = b, PlasticityCostSD = COST)
            obs = k50(rs) / kb
            ffit, fvg, full = S(pc, vg0) / S(p0, vg0), S(p0, vgc) / S(p0, vg0), S(pc, vgc) / S(p0, vg0)
            @printf("  b/B = %.2f cost: V_G at the move ×%.2f, N %.0f; observed k50 ratio (effective) %.3f; theory: fitness margin alone %.3f, V_G alone %.3f, both (P2) %.3f; observed/P2 %.2f; margin × (V_G ratio)^e %.3f\n",
                    b, vgc / vg0, mean(r.n0 for r in alive0), obs, ffit, fvg, full, obs / full, ffit * (vgc / vg0)^e)
            lo = sort(unique(r.keff for r in rs))[1:min(2, end)]
            ks = unique(r.keff for r in nc)
            for kk in lo
                c = [r.vgmove for r in rs if r.keff == kk && isfinite(r.vgmove)]
                k0 = ks[argmin(abs.(log.(ks) .- log(kk)))]
                n = [r.vgmove for r in nc if r.keff == k0 && isfinite(r.vgmove)]
                isempty(c) || isempty(n) || @printf("      V_G in generations 500–2000 after the move at effective rate %.4g: %.4g (cost, %d alive) vs %.4g (no cost at %.4g, %d alive): ×%.2f\n",
                                                    kk, mean(c), length(c), mean(n), k0, length(n), mean(c) / mean(n))
            end
        end
    end
end

function pilot(seed0)
    println("Pilot: $(NPILOT) burn-ins per family, effective rates K_REF × 0.03–1.5 (15, geometric), seeds from $seed0, $(Threads.nthreads()) threads")
    rows = sweep(js -> job(js, seed0, true), jobs(NPILOT); cache = joinpath(@__DIR__, "plasticity_erosion_pilot.csv"),
                 id = j -> "$(j.fam)-$(j.set)-$(j.dem)-$(j.b)-$(j.rep)", progress = 20)
    br = branches(rows)
    println("const K50_PILOT entries:")
    for (set, dem, b) in cells()
        nb = [x for x in br if x.fam === :nocost && x.set === set && x.dem === dem]
        cb = [x for x in br if x.fam === :cost && x.set === set && x.dem === dem && x.famb == b]
        for (a, xs) in ((:zero, [x for x in nb if x.arm === :zero]), (:C, [x for x in nb if x.arm === :C && x.b == b]),
                        (:B, [x for x in cb if x.arm === :B]), (:A, [x for x in cb if x.arm === :A]),
                        (:Aprime, [x for x in cb if x.arm === :Aprime]))
            ks = sort(unique(x.keff for x in xs))
            surv = [mean(x.alive for x in xs if x.keff == kk) for kk in ks]
            @printf("  (%s, %s, %.2f, %s): k50 %.5g; survival %s\n", set, dem, a === :zero ? 0.0 : b, a, k50(xs),
                    join((@sprintf("%.2f", s) for s in surv), " "))
        end
    end
end

function main()
    if !isempty(ARGS) && ARGS[1] == "posthoc"
        return posthoc()
    elseif !isempty(ARGS) && ARGS[1] == "pilot"
        return pilot(8_000_000_000)
    end
    seed0 = isempty(ARGS) ? 7_000_000_000 : parse(Int, ARGS[1])
    isempty(K50_PILOT) && error("run the pilot and enter its half-survival rates in K50_PILOT first")
    path = joinpath(@__DIR__, seed0 == 7_000_000_000 ? "plasticity_erosion.csv" : "plasticity_erosion_seed$(seed0).csv")
    println("Plasticity's cost and the sustainable rate: factorial, kernel, $(NREP) burn-ins per family, seeds from $seed0, ",
            "$(Threads.nthreads()) threads")
    t0 = time()
    rows = sweep(j -> job(j, seed0, false), jobs(NREP); cache = path, id = j -> "$(j.fam)-$(j.set)-$(j.dem)-$(j.b)-$(j.rep)",
                 cost = j -> (j.set === :L500 ? 2 : 1), progress = 20)
    @printf("%d burn-ins in %s (%.1f min this session)\n", length(rows), basename(path), (time() - t0) / 60)
    ok = analyse(rows, seed0)
    println(ok ? "\nALL CRITERIA PASS" : "\nSOME CRITERIA FAIL")
    ok || exit(1)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
