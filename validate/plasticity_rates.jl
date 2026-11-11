# What phenotypic plasticity does to the rate of change a population can sustain: a
# pre-specified experiment with the kernel (ROADMAP F2, stage 3).
#
# Everything from here to "END OF PRE-SPECIFICATION" was written before the first run of
# the main experiment. A pilot WITHOUT plasticity (b = 0 only; `pilot` argument, logged in
# plasticity_rates_pilot.log) was run first, only to place the grid of rates; its
# half-survival rates are the constants K50_PILOT below.
#
# ── Design ────────────────────────────────────────────────────────────────────────────
# Regime: the published sweep's design (1000 diploids, CrossoverModel = :qtsim, 10,000
# generations of stabilising selection, then the optimum moves continuously), with two
# synthetic parameter sets within the paper's priors, those of validate/gsl_locked.jl's
# synthetic inputs (nothing is read from Repository/):
#   L100: NoQTL 100, u 3, V_P 0.03, h² 0.4, m 0.005, ω 4, μ 1e-4
#   L500: NoQTL 500, u 3, V_P 0.05, h² 0.3, m 0.008, ω 3, μ 1e-4
# Demography: QTSim's fixed size, and Beverton–Holt with R = 2 (K = 1000).
# Plasticity: b/B ∈ {0, 0.25, 0.5, 0.75} with B = 1 (OptimumSlope), a perfect cue
# (CueLag = 0); without a cost, and with the cost ω_b = 0.5 (PlasticityCostSD; fitness
# times exp(−b²/(2ω_b²)) = 0.88, 0.61, 0.32 for the three slopes).
# Rates: eleven rates per cell, k = K50_PILOT[L, demography] × {0.6, 0.8, 0.88, 0.94, 0.97, 1,
# 1.03, 1.06, 1.12, 1.25, 1.6} / (1 − b/B): the same grid of EFFECTIVE rates (1 − b/B)·k in
# every no-cost cell, so that under the theory's scaling their survival curves coincide.
# (The grid first written here was {0.4, …, 2.2} in nine steps; the pilot's second stage,
# still b = 0 only, showed that survival falls from 1 to 0 within about ±10% of k50, so the
# grid was made finer around 1 before any run with plasticity.) In the cells with a cost the
# grid is multiplied by the predicted effect of the cost, sustainable_rate(b, cost)·(1 − b/B)/
# sustainable_rate(0) at the nominal V_G = h²V_P, so that it brackets the prediction; if the
# prediction is off by more than the grid's span the half-survival rate is not bracketed,
# and the cell fails.
# Replicates: 40 burn-ins per (parameter set, demography). Without a cost the burn-in does
# not depend on the slope (the optimum is at the origin, so the plastic part is exactly 0;
# tested bitwise in test/runtests.jl), so each burn-in is shared by all four slopes and all
# nine rates (a paired design); with the cost each slope has its own 40 burn-ins. Each
# branch gets its own random numbers after the move. A run is censored 2000 generations
# after the move (as validate/published_sweep.jl's default cap).
#
# ── Measures ──────────────────────────────────────────────────────────────────────────
# survival: alive 2000 generations after the move. k50: the rate at which half survive,
# from a logistic regression of survival on ln k over a cell's eleven rates (maximum
# likelihood), k50 = exp(−a/c) for logit P(survive) = a + c ln k. Intervals: 95% percentile
# intervals from 1000 bootstrap resamples of the replicates (blocks: a replicate's whole set
# of rates; for ratios within the no-cost design the same resampled replicates in numerator
# and denominator).
#
# ── Predictions and criteria ──────────────────────────────────────────────────────────
# P1 (no cost): k50(b)/k50(0) = 1/(1 − b/B), i.e. 1.33, 2 and 4. PASS for a cell when the
#     estimate of k50(b)(1 − b/B)/k50(0) is within 15% of 1 (0.85–1.15). 12 cells
#     (2 parameter sets × 2 demographies × 3 slopes). The interval is reported; a cell whose
#     interval excludes 1 is flagged even if it passes.
# P2 (cost): k50(b, cost)/k50(0) = sustainable_rate(b, cost)/sustainable_rate(0), evaluated
#     at the cell's mean genetic variance at the move (QTLSim.sustainable_rate: the
#     deterministic threshold of the moments recursion). PASS when observed/predicted is
#     within 25% of 1. For Beverton–Holt with b = 0.75 the prediction is extinction before
#     the move (R·W̄ = 2·0.997·0.32 < 1): PASS when at least 90% of the burn-ins die out.
#     12 cells, of which 2 are of the second kind.
# P3 (lag, no cost): at equal effective rate the steady lag is the same. At the two lowest
#     rates of the grid, the ratio of the mean lag θ − z̄ over generations 500–2000 after the
#     move (runs alive at 2000) to that at b = 0 is within 15% of 1. 12 cells.
# Descriptive, no criterion: k50 against sustainable_rate itself (the deterministic
# threshold should overestimate), survival by rate, and the median time to extinction.
# The headline number is the P1 ratio: how much the half-survival rate rises with b/B.
# Failures stay in the log. If the model or the script's analysis is changed after a
# failed run, the change is stated here and confirmed with fresh seeds (SEED0).
#
# END OF PRE-SPECIFICATION
#
# Notes added after the run (the pre-specification above is unchanged):
# - The first pass of the analysis stopped on a programming error in the printing of the
#   results (`mean(...; init)`, which Statistics.mean does not take), after the runs were in
#   the CSV. It was fixed (`meanor`) and the analysis rerun on the same rows; no run was
#   repeated and no criterion changed.
# - Run 1 (seeds from 20261110) ran all 640 jobs (12,320 runs, validate/plasticity_rates.csv)
#   in 610.6 min of wall clock; the machine slept during jobs 540-640, which took most of
#   that. Its analysis then crashed at the line above (stderr kept in
#   plasticity_rates_run1.err). The analysis was rerun from the CSV on 2026-11-11 with the fix;
#   plasticity_rates.log is that rerun. Result: P1 12/12, P3 12/12, P2 8/12.
# - Why P2 fails, found after the run (the criteria are unchanged): the cost lowers the
#   genetic variance before the move. Viability selection with every fitness times
#   exp(−b²/(2ω_b²)) leaves fewer parents (more drift), and under Beverton–Holt a smaller
#   population too: V_G at the move falls 15%/33%/58% (L100, fixed size, b/B = 0.25/0.5/0.75),
#   5%/27%/54% (L500, fixed), and 22%/72% (L100) and 15%/66% (L500) under Beverton–Holt,
#   where N at the move falls to about 340 at b/B = 0.5. The grid's cost factor (costfactor)
#   used the nominal V_G for every cell, so it ignored this and sat too high: under
#   Beverton–Holt at b/B = 0.5 every run died even at the lowest rate (k50 not bracketed,
#   NaN, below 0.63·k50(0)); the prediction there, with the realised V_G, was 0.29–0.36·k50(0),
#   the same direction. At fixed size the prediction, which does use each cell's realised V_G,
#   underestimates the costly cells' k50 by 0.5–39% (more at larger slopes): the deterministic
#   threshold is 5–7 times the observed k50 at fixed size, and k50 responds less than
#   proportionally to V_G there. A fresh-seed confirmation was not run: no model or criterion
#   was changed.
# - P1 and P3 test an identity, not an approximation: with a perfect cue and no cost,
#   z − θ = g + e − (1 − b/B)θ, so a plastic population at rate k is, in distribution, a
#   population without plasticity at rate (1 − b/B)k. They check the implementation.
#
# Run from the package root, at most 3 threads (about 1.5 CPU-hours):
#   julia -t 3 --project=. validate/plasticity_rates.jl pilot [fine]   # places the grid (b = 0 only)
#   julia -t 3 --project=. validate/plasticity_rates.jl [SEED0]   # the experiment; resumes from its CSV
# Writes validate/plasticity_rates.csv (one row per run) and prints the analysis, logged in
# validate/plasticity_rates.log.

using QTLSim, Statistics, Random, Printf

const BURN = 10_000
const CAP = 2_000
const NREP = 40
const GRID = (0.6, 0.8, 0.88, 0.94, 0.97, 1.0, 1.03, 1.06, 1.12, 1.25, 1.6)
const SLOPES = (0.0, 0.25, 0.5, 0.75)
const COST = 0.5
const SETS = (L100 = Params(Sample = 10, PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0,
                            PhenotypeVar = 0.03, Heritability = 0.4, MutationScalar = 0.005,
                            MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous),
              L500 = Params(Sample = 10, PopSize = 1000, Ploidy = 2, NoQTL = 500, RecombPoissonMean = 3.0,
                            PhenotypeVar = 0.05, Heritability = 0.3, MutationScalar = 0.008,
                            MutationRatePerGamete = 1e-4, InitFitnessSD = 3.0, MoveType = :continuous))
const DEMS = (fixed = (Demography = :fixed,), bh = (Demography = :beverton_holt, GrowthRate = 2.0))
# half-survival rates of the pilot's second stage (b = 0, 12 replicates, plasticity_rates_pilot.log), per generation
const K50_PILOT = Dict((:L100, :fixed) => 0.02219, (:L100, :bh) => 0.0153, (:L500, :fixed) => 0.05212, (:L500, :bh) => 0.03774)

"The predicted effect of the cost on the half-survival effective rate, at the nominal V_G = h²V_P (1 without a cost)."
function costfactor(set, dem, b)
    p0 = Params(SETS[set]; DEMS[dem]...)
    f = sustainable_rate(Params(p0; Plasticity = b, PlasticityCostSD = COST); VG = geneticVar(p0)) * (1 - b) /
        sustainable_rate(p0; VG = geneticVar(p0))
    f > 0 ? f : 1.0
end

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

"Run a copy of the burnt-in `w0` for up to CAP generations with the optimum moving at `k` and slope `b`."
function branch(w0::World, k, b, seed)
    w = deepcopy(w0)
    w.p = Params(w0.p; MoveFitnessSDGen = BURN, MoveFitnessSD = k, Plasticity = b, NoGenerations = BURN + CAP)
    rng = MersenneTwister(seed)
    rec = QTLSim.records()
    text = CAP
    for it in BURN:BURN+CAP-1
        w.gen = it
        generation!(w, rng; rec)
        w.extinct && (text = it - BURN; break)
    end
    late = findall(g -> g >= BURN + 500, rec.gen)
    (extinct = w.extinct, text, lag = w.extinct || isempty(late) ? NaN : mean(rec.theta[late] .- rec.pheno_mean[late]),
     vgmove = w.extinct || isempty(late) ? NaN : mean(rec.geno_var[late]))
end

genvar(w::World) = (develop!(w); var(view(w.genotype, 1:w.N)))

# one job: a burn-in and all its branches. Without a cost one burn-in serves all slopes.
function job(set, dem, cost, bs, rep, rates, seed0)
    p = Params(SETS[set]; DEMS[dem]..., NoGenerations = BURN,
               (cost ? (Plasticity = only(bs), PlasticityCostSD = COST) : (;))...)
    ci = findfirst(==((set, dem)), [(s, d) for s in keys(SETS) for d in keys(DEMS)])
    seed = seed0 + 100_000 * ci + 10_000 * (cost ? round(Int, 4 * only(bs)) : 0) + rep
    w0 = burnin(p, seed)
    rows = NamedTuple[]
    vg0, n0 = w0 === nothing ? (NaN, 0) : (genvar(w0), w0.N)
    for b in bs, (i, k) in enumerate(rates(b))
        r = w0 === nothing ? (extinct = true, text = 0, lag = NaN, vgmove = NaN) :
            branch(w0, k, b, seed + 1_000_000 * (i + 10 * round(Int, 4b)))
        push!(rows, (; set, dem, cost, b, rep, rate = k, keff = (1 - b) * k, r.extinct, r.text, vg0, n0, r.lag, r.vgmove))
    end
    rows
end

const COLS = (:set, :dem, :cost, :b, :rep, :rate, :keff, :extinct, :text, :vg0, :n0, :lag, :vgmove)

function readrows(path)
    isfile(path) || return NamedTuple[]
    rows = NamedTuple[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ',')
        length(f) == length(COLS) || continue           # a line cut short by an interruption
        push!(rows, (set = Symbol(f[1]), dem = Symbol(f[2]), cost = f[3] == "true", b = parse(Float64, f[4]),
                     rep = parse(Int, f[5]), rate = parse(Float64, f[6]), keff = parse(Float64, f[7]),
                     extinct = f[8] == "true", text = parse(Int, f[9]), vg0 = parse(Float64, f[10]),
                     n0 = parse(Int, f[11]), lag = parse(Float64, f[12]), vgmove = parse(Float64, f[13])))
    end
    rows
end

function runjobs(jobs, path, seed0)
    done = readrows(path)
    key(r) = (r.set, r.dem, r.cost, r.rep, r.cost ? r.b : -1.0)
    have = Dict{Any,Int}()
    for r in done; have[key(r)] = get(have, key(r), 0) + 1; end
    todo = [j for j in jobs if get(have, (j.set, j.dem, j.cost, j.rep, j.cost ? only(j.bs) : -1.0), 0) < length(j.bs) * length(GRID)]
    # drop the partial rows of jobs to be redone
    redo = Set((j.set, j.dem, j.cost, j.rep, j.cost ? only(j.bs) : -1.0) for j in todo)
    keep = [r for r in done if !(key(r) in redo)]
    open(path, "w") do io
        println(io, join(COLS, ','))
        for r in keep; println(io, join((getfield(r, c) for c in COLS), ',')); end
    end
    lk = ReentrantLock()
    n = Threads.Atomic{Int}(0)
    t0 = time()
    Threads.@threads :greedy for j in todo
        rows = job(j.set, j.dem, j.cost, j.bs, j.rep, j.rates, seed0)
        lock(lk) do
            open(path, "a") do io
                for r in rows; println(io, join((getfield(r, c) for c in COLS), ',')); end
            end
            Threads.atomic_add!(n, 1)
            n[] % 20 == 0 && (@printf(stderr, "  %d/%d jobs, %.1f min\n", n[], length(todo), (time() - t0) / 60); flush(stderr))
        end
    end
    readrows(path)
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
        s = min(1.0, 5 / max(abs(da), abs(dc), 1e-12))          # damped Newton
        a += s * da; c += s * dc
        max(abs(da), abs(dc)) < 1e-10 && break
    end
    (a - c * xm, c)
end

"The half-survival rate of rows `rs`: exp(−a/c); NaN when survival does not fall with the rate in them."
function k50(rs)
    y = [!r.extinct for r in rs]
    (all(y) || !any(y)) && return NaN
    a, c = logistic([log(r.rate) for r in rs], y)
    c < 0 ? exp(-a / c) : NaN
end

meanor(itr) = (v = collect(itr); isempty(v) ? NaN : mean(v))
quant(v, q) = (w = sort(filter(isfinite, v)); isempty(w) ? NaN : w[clamp(round(Int, q * length(w) + 0.5), 1, length(w))])
ci(v) = (quant(v, 0.025), quant(v, 0.975))

function analyse(rows, seed0)
    rng = MersenneTwister(seed0 + 99)
    cell(set, dem, cost, b) = [r for r in rows if r.set == set && r.dem == dem && r.cost == cost && r.b == b]
    byrep(rs) = (d = Dict{Int,Vector{eltype(rs)}}(); for r in rs; push!(get!(d, r.rep, eltype(rs)[]), r); end; d)
    npass = Dict(:P1 => [0, 0], :P2 => [0, 0], :P3 => [0, 0])
    verdict(crit, ok) = (npass[crit][1] += ok; npass[crit][2] += 1; ok ? "PASS" : "FAIL")
    for set in keys(SETS), dem in keys(DEMS)
        p0 = Params(SETS[set]; DEMS[dem]...)
        sdp = sqrt(p0.PhenotypeVar)
        base = cell(set, dem, false, 0.0)
        isempty(base) && continue
        vg0 = mean(r.vg0 for r in base)
        kb = k50(base)
        B0 = byrep(base)
        reps = sort(collect(keys(B0)))
        boot0 = [k50(reduce(vcat, (B0[reps[i]] for i in rand(rng, eachindex(reps), length(reps))))) for _ in 1:1000]
        @printf("\n══ %s, %s: V_G at the move %.4g (h² realised %.2f), N at the move %.0f ══\n", set, dem, vg0,
                vg0 / (vg0 + envVar(p0)), mean(r.n0 for r in base))
        @printf("  b = 0: k50 = %.4g per generation [%.4g, %.4g] = %.4f phenotypic SD; deterministic sustainable_rate %.4g (k50/it %.2f)\n",
                kb, ci(boot0)..., kb / sdp, sustainable_rate(p0; VG = vg0), kb / sustainable_rate(p0; VG = vg0))
        println("  survival by effective rate (1 − b/B)k, as multiples of the pilot's k50 (with a cost: times the predicted effect of the cost, ",
                join((@sprintf("%.3f", costfactor(set, dem, b)) for b in SLOPES[2:end]), ", "), "): ", join(GRID, "  "))
        for cost in (false, true), b in SLOPES
            cost && b == 0 && continue
            rs = cell(set, dem, cost, b)
            isempty(rs) && continue
            ks = sort(unique(r.keff for r in rs))
            surv = [mean(!r.extinct for r in rs if isapprox(r.keff, k; rtol = 1e-9)) for k in ks]
            med = [median(r.text for r in rs if isapprox(r.keff, k; rtol = 1e-9)) for k in ks]
            @printf("  b/B = %.2f%s  survival %s\n               median generations to extinction (2000 = alive) %s\n", b, cost ? " cost" : "     ",
                    join((@sprintf("%.2f", s) for s in surv), " "), join((@sprintf("%.0f", m) for m in med), " "))
        end
        println("  half-survival rates:")
        for cost in (false, true), b in SLOPES[2:end]
            rs = cell(set, dem, cost, b)
            isempty(rs) && continue
            kc = k50(rs)
            R = byrep(rs)
            rp = sort(collect(keys(R)))
            # paired bootstrap without a cost (the same burn-ins as b = 0); with the cost, independent replicates
            boot = map(1:1000) do _
                if cost
                    k50(reduce(vcat, (R[rp[i]] for i in rand(rng, eachindex(rp), length(rp))))) /
                    k50(reduce(vcat, (B0[reps[i]] for i in rand(rng, eachindex(reps), length(reps)))))
                else
                    ix = rand(rng, eachindex(reps), length(reps))
                    k50(reduce(vcat, (get(R, reps[i], eltype(rs)[]) for i in ix))) / k50(reduce(vcat, (B0[reps[i]] for i in ix)))
                end
            end
            ratio = kc / kb
            lo, hi = ci(boot)
            if !cost
                pred = 1 / (1 - b)
                ok = isfinite(ratio) && 0.85 <= ratio / pred <= 1.15
                flag = isfinite(lo) && (lo > pred || hi < pred) ? "  (interval excludes the prediction)" : ""
                @printf("    b/B = %.2f       k50 = %.4g, ratio to b = 0: %.3f [%.3f, %.3f]; predicted %.3f; observed/predicted %.3f [%.3f, %.3f]  P1 %s%s\n",
                        b, kc, ratio, lo, hi, pred, ratio / pred, lo / pred, hi / pred, verdict(:P1, ok), flag)
            else
                pc = Params(p0; Plasticity = b, PlasticityCostSD = COST)
                vgc = meanor(r.vg0 for r in rs if isfinite(r.vg0))
                dead = mean(r.n0 == 0 for r in rs)
                pred = sustainable_rate(pc; VG = isfinite(vgc) ? vgc : vg0) / sustainable_rate(p0; VG = vg0)
                if pred == 0
                    @printf("    b/B = %.2f cost  predicted to die out before the move (R·W̄ < 1): %.0f%% of burn-ins died  P2 %s\n",
                            b, 100dead, verdict(:P2, dead >= 0.9))
                else
                    ok = isfinite(ratio) && 0.75 <= ratio / pred <= 1.25
                    @printf("    b/B = %.2f cost  k50 = %.4g, ratio to b = 0: %.3f [%.3f, %.3f]; predicted %.3f (no cost: %.3f); observed/predicted %.3f [%.3f, %.3f]; %.0f%% of burn-ins died  P2 %s\n",
                            b, kc, ratio, lo, hi, pred, 1 / (1 - b), ratio / pred, lo / pred, hi / pred, 100dead, verdict(:P2, ok))
                end
            end
        end
        println("  lag θ − z̄ at equal effective rate (generations 500–2000 after the move), ratio to b = 0:")
        for b in SLOPES[2:end]
            rs = cell(set, dem, false, b)
            isempty(rs) && continue
            ks = sort(unique(r.keff for r in base))[1:2]
            ratios = map(ks) do k
                l1 = [r.lag for r in rs if isapprox(r.keff, k; rtol = 1e-9) && isfinite(r.lag)]
                l0 = [r.lag for r in base if isapprox(r.keff, k; rtol = 1e-9) && isfinite(r.lag)]
                isempty(l1) || isempty(l0) ? NaN : mean(l1) / mean(l0)
            end
            ok = all(x -> isfinite(x) && 0.85 <= x <= 1.15, ratios)
            @printf("    b/B = %.2f       %.3f, %.3f at the two lowest rates (lag at b = 0: %.3g, %.3g phenotypic SD)  P3 %s\n", b, ratios...,
                    (meanor(r.lag for r in base if isapprox(r.keff, k; rtol = 1e-9) && isfinite(r.lag)) / sdp for k in ks)...,
                    verdict(:P3, ok))
        end
    end
    println()
    for c in (:P1, :P2, :P3)
        @printf("%s: %d of %d cells pass\n", c, npass[c]...)
    end
    all(v -> v[1] == v[2], values(npass))
end

# the pilot's second stage: the rates between which survival fell from 1 to 0 in its first
const PILOT_BRACKET = Dict((:L100, :fixed) => (0.0148, 0.0245), (:L100, :bh) => (0.00899, 0.0245),
                           (:L500, :fixed) => (0.0404, 0.0667), (:L500, :bh) => (0.0245, 0.0404))

function pilot(fine::Bool)
    nrep = fine ? 12 : 8
    println("Pilot without plasticity (b = 0), $nrep replicates, to place the grid of rates", fine ? " (second stage: within the first stage's brackets)" : "",
            "; $(Threads.nthreads()) threads")
    for set in keys(SETS), dem in keys(DEMS)
        p = Params(SETS[set]; DEMS[dem]..., NoGenerations = BURN)
        rates = fine ? exp.(range(log.(PILOT_BRACKET[(set, dem)])...; length = 12)) : exp.(range(log(0.002), log(0.3); length = 11))
        out = Vector{Any}(undef, nrep)
        Threads.@threads :greedy for rep in 1:nrep
            w0 = burnin(p, (fine ? 666_000 : 555_000) + rep + 100 * findfirst(==((set, dem)), [(s, d) for s in keys(SETS) for d in keys(DEMS)]))
            out[rep] = [(rate = k, extinct = w0 === nothing ? true : branch(w0, k, 0.0, 777_000 + 100rep + i).extinct)
                        for (i, k) in enumerate(rates)]
        end
        rs = reduce(vcat, out)
        @printf("%s %s: survival by rate %s\n   rates %s\n   k50 = %.4g per generation (%.4f phenotypic SD)\n", set, dem,
                join((@sprintf("%.2f", mean(!r.extinct for r in rs if r.rate == k)) for k in rates), " "),
                join((@sprintf("%.3g", k) for k in rates), " "), k50(rs), k50(rs) / sqrt(p.PhenotypeVar))
        flush(stdout)
    end
end

function main()
    if !isempty(ARGS) && ARGS[1] == "pilot"
        return pilot(length(ARGS) > 1 && ARGS[2] == "fine")
    end
    seed0 = isempty(ARGS) ? 20_261_110 : parse(Int, ARGS[1])
    any(isnan, values(K50_PILOT)) && error("run the pilot and enter its half-survival rates in K50_PILOT first")
    path = joinpath(@__DIR__, isempty(ARGS) ? "plasticity_rates.csv" : "plasticity_rates_seed$(seed0).csv")
    jobs = NamedTuple[]
    for set in keys(SETS), dem in keys(DEMS), rep in 1:NREP
        rates = b -> [K50_PILOT[(set, dem)] * g / (1 - b) for g in GRID]
        push!(jobs, (; set, dem, cost = false, bs = collect(SLOPES), rep, rates))
        for b in SLOPES[2:end]
            f = costfactor(set, dem, b)
            push!(jobs, (; set, dem, cost = true, bs = [b], rep, rates = bb -> [K50_PILOT[(set, dem)] * f * g / (1 - bb) for g in GRID]))
        end
    end
    println("Plasticity and the sustainable rate of change: kernel, $(NREP) burn-ins per cell, seeds from $seed0, ",
            "$(Threads.nthreads()) threads")
    t0 = time()
    rows = runjobs(jobs, path, seed0)
    @printf("%d runs in %s (%.1f min this session)\n", length(rows), basename(path), (time() - t0) / 60)
    ok = analyse(rows, seed0)
    println(ok ? "\nALL CRITERIA PASS" : "\nSOME CRITERIA FAIL")
    ok || exit(1)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
