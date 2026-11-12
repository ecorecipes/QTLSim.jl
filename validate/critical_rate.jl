# The half-survival rate of a continuously moving optimum, and how a constant cost of
# fitness lowers it: a pre-specified test of drift_critical_rate and half_survival_rate
# (src/critical.jl; ROADMAP F8 and the open problem of validate/plasticity_erosion.jl).
#
# Everything from here to "END OF PRE-SPECIFICATION" was written before the first run of
# the experiment (seeds from SEED0 below). The post hoc evaluation (`posthoc` argument,
# critical_rate_posthoc.log) was run first, on existing data only, and the tolerances below
# were chosen knowing its results; that is stated where it matters.
#
# ── Theory (src/critical.jl, validate/symbolic/check_critical.jl) ─────────────────────────
# At a steady lag λ mean fitness is W̄ = c₀e^{−c}ω/√s·e^{−λ²/(2s)}; the breeders are n = NW̄
# (fixed size) or the Beverton–Holt equilibrium K(RW̄ − 1)/(R − 1); N_e = effective_size(N, n);
# V_G is the house-of-Gauss variance with that N_e; the response is V_G λ/(V_G + V_s). The
# critical rate k_c = max over λ of the response (drift_critical_rate): above it the lag,
# the loss of breeders and the loss of variance feed each other. Since n depends on λ and
# the log cost c only through u = λ²/(2s) + c, the envelope theorem gives
#     k_c(c)/k_c(0) = 1 − c/M + O(c²),   M = λ*²/s
# (the linear-margin form, with M twice the lag load at the critical lag). The half-survival
# rate within 2000 generations is predicted by half_survival_rate: the moments lowering with
# the allele-frequency spectrum and drift along the move, plus demographic stochasticity
# along that trajectory (extinction_distribution), bisected for P(survive) = 1/2.
#
# Post hoc (critical_rate_posthoc.log; not a test): on plasticity_erosion.log's cells,
# drift_critical_rate's ratio k_c(c)/k_c(0) was within 0.05 of the observed k50 ratio of
# arm A (the cost before and after the move) in 7 of 8 cells (max 0.052 at c = 1.125, fixed
# size; 0.027 at the lowest ratio, 0.16 vs 0.19); the M it implies, 1.75–1.79 (fixed) and
# 0.52 (R = 2), against 2.0–2.5 and 0.59–0.63 fitted from the arms. half_survival_rate's
# k50 without a cost was 0.94–1.08 of the kernel's (100 loci) and 0.69–0.77 (500 loci).
# On vignettes 09 and 11's 710 populations, k50 was 0.29–1.06 × drift_critical_rate (falling
# with the number of loci) against 0.07–0.56 × sustainable_rate.
#
# ── Design ────────────────────────────────────────────────────────────────────────────
# Parameter sets L100 and L500 of plasticity_erosion.jl (published-like: N = K = 1000,
# diploid, MutationRatePerGamete 1e-4, RecombPoissonMean 3, CrossoverModel = :qtsim, the
# default, continuous move). Demographies: fixed size; Beverton–Holt R = 5, 2, 1.5.
# Log-fitness costs c, paid by every individual from generation 0 (burn-in and move), as a
# reaction norm of slope 1 with OptimumSlope 1e12 and PlasticityCostSD 1/√(2c) (fitness ×
# e^{−c}; the plastic response is 1e-12 of the move, H4 of plasticity_erosion.log):
#     fixed and R = 5: c ∈ {0, 0.25, 0.5, 1.0};  R = 2: {0, 0.1, 0.2, 0.4};  R = 1.5: {0, 0.05, 0.1, 0.2}
# (chosen so that drift_critical_rate's ratio spans about 0.2–0.86 in every demography,
# and R e^{−c} W̄(0) > 1.2, so the populations persist before the move).
# 32 cells, 30 burn-ins each (10,000 generations at a fixed optimum), each branched at 13
# rates of the move, k = center × {0.5, 0.62, 0.72, 0.8, 0.86, 0.92, 0.97, 1.03, 1.09,
# 1.17, 1.3, 1.55, 2.0}, center = drift_critical_rate(p; cost = c).k × F, F = 0.84 (L100),
# 0.52 (L500) (the post hoc k50/k_c without a cost; only places the grid), run up to 2000
# generations after the move. Fresh seeds: SEED0 = 9_100_000_000 (plasticity_rates.jl used
# 20_261_110-based seeds, plasticity_erosion.jl 7e9 and 8e9).
#
# ── Measures ──────────────────────────────────────────────────────────────────────────
# k50 per cell: logistic regression of survival (alive 2000 generations after the move) on
# ln k, as in plasticity_erosion.jl; NaN if the grid does not bracket it. ρ_obs(c) =
# k50(c)/k50(0) within a (set, demography) group. Intervals: 95% percentile intervals from
# 1000 bootstrap resamples of burn-ins (each with all its branches), cells independently.
#
# ── Criteria ──────────────────────────────────────────────────────────────────────────
# C1 (the relative effect of a cost): |ρ_obs(c) − ρ_pred(c)| ≤ 0.06, ρ_pred =
#     drift_critical_rate(p; cost = c).k / drift_critical_rate(p; cost = 0).k. 24 cells.
#     (The tolerance was set after the post hoc look, where the largest error was 0.052.)
# C2 (the linear-margin form and its scale): in each of the 8 groups, M_obs, the least-
#     squares fit of 1 − ρ_obs = c/M over its three costs, within 0.8–1.25 of M_theory =
#     drift_critical_rate(p).M. Also reported: the same fit of the square-root form
#     ρ = √(1 − c/M′) (sustainable_rate's), and which has the smaller sum of squares.
# C3 (the half-survival rate itself): k50_obs/half_survival_rate(p; cost = c).k50 within
#     0.75–1.33. 32 cells. (The post hoc look predicts failures at 500 loci, near 0.7.)
# Verdicts:
#   relative effect of a cost (drift_critical_rate): SUPPORTED if C1 passes in ≥ 20 of 24;
#     REFUTED if in < 12; otherwise MIXED.
#   linear-margin form: SUPPORTED if C2 passes in ≥ 6 of 8 groups and the linear form has
#     the smaller sum of squares in ≥ 6 of 8; REFUTED if the square-root form has it in
#     ≥ 6 of 8; otherwise MIXED.
#   half_survival_rate: SUPPORTED if C3 passes in ≥ 26 of 32; REFUTED if in < 16; else MIXED.
# Descriptive: k50/k_c, k50/sustainable_rate (V_G the cell's mean at the move), and
# sustainable_rate's ratio for the cost, which C1's tolerance is also applied to.
# Power: with 30 burn-ins and 13 rates, k50's interval should be about ±2% (±1.5% with 40
# and 11 in plasticity_erosion.log), so ρ's standard error is about 0.015–0.02 and C1's
# tolerance about 3 of them; M_obs about ±5–10%.
# Failures stay in the log. If the theory or the analysis is changed after a failed run,
# the change is stated here and confirmed on fresh seeds.
#
# END OF PRE-SPECIFICATION
#
# Run from the package root, at most 3 threads (estimate: about 960 jobs × 7–10 s, 2–3
# CPU-hours, plus about 30 min for the 32 half_survival_rate predictions):
#   julia -t 3 --project=vignettes validate/critical_rate.jl posthoc   # existing data only
#   julia -t 3 --project=vignettes validate/critical_rate.jl [SEED0]   # the experiment; resumes
# Writes validate/critical_rate.csv (sweep cache, one row per burn-in with its branches),
# validate/critical_rate_pred.csv (the predictions), and prints the analysis
# (critical_rate.log).

using QTLSim, Statistics, Random, Printf, DelimitedFiles, XLSX
import Catlab

const BURN = 10_000
const CAP = 2_000
const NREP = 30
const SEED0 = 9_100_000_000
const GRID = (0.5, 0.62, 0.72, 0.8, 0.86, 0.92, 0.97, 1.03, 1.09, 1.17, 1.3, 1.55, 2.0)
const F = Dict(:L100 => 0.84, :L500 => 0.52)
const SETS = (L100 = Params(Sample = 10, PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0,
                            PhenotypeVar = 0.03, Heritability = 0.4, MutationScalar = 0.005,
                            MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous),
              L500 = Params(Sample = 10, PopSize = 1000, Ploidy = 2, NoQTL = 500, RecombPoissonMean = 3.0,
                            PhenotypeVar = 0.05, Heritability = 0.3, MutationScalar = 0.008,
                            MutationRatePerGamete = 1e-4, InitFitnessSD = 3.0, MoveType = :continuous))
const DEMS = (:fixed, :R5, :R2, :R15)
const GROWTH = Dict(:R5 => 5.0, :R2 => 2.0, :R15 => 1.5)
const COSTS = Dict(:fixed => (0.0, 0.25, 0.5, 1.0), :R5 => (0.0, 0.25, 0.5, 1.0), :R2 => (0.0, 0.1, 0.2, 0.4),
                   :R15 => (0.0, 0.05, 0.1, 0.2))

demp(set, dem) = dem === :fixed ? Params(SETS[set]; Demography = :fixed) :
                 Params(SETS[set]; Demography = :beverton_holt, GrowthRate = GROWTH[dem])
withcost(p, c) = c == 0 ? p : Params(p; Plasticity = 1.0, OptimumSlope = 1e12, PlasticityCostSD = 1 / sqrt(2c))
cells() = [(set, dem, c) for set in keys(SETS) for dem in DEMS for c in COSTS[dem]]
center(set, dem, c) = drift_critical_rate(demp(set, dem); cost = c).k * F[set]
rates(set, dem, c) = [center(set, dem, c) * g for g in GRID]

function burnin(p::Params, seed)
    w = World(Params(p; NoGenerations = BURN)); rng = MersenneTwister(seed)
    init_env!(w, rng)
    for it in 0:BURN-1
        w.gen = it
        generation!(w, rng)
        w.extinct && return nothing
    end
    w
end

"Survival and generations to extinction of a copy of `w0` with the optimum moving at `k`."
function branch(w0::World, k, seed)
    w = deepcopy(w0)
    w.p = Params(w0.p; MoveFitnessSDGen = BURN, MoveFitnessSD = k, NoGenerations = BURN + CAP)
    rng = MersenneTwister(seed)
    for it in BURN:BURN+CAP-1
        w.gen = it
        generation!(w, rng)
        w.extinct && return (alive = false, text = it - BURN)
    end
    (alive = true, text = CAP)
end

genvar(w::World) = (develop!(w); var(view(w.genotype, 1:w.N)))

function job(j, seed0)
    (; set, dem, c, rep) = j
    id = findfirst(==((set, dem, c)), cells())
    seed = seed0 + 1_000_000 * id + rep
    p = withcost(demp(set, dem), c)
    w0 = burnin(p, seed)
    vg0, n0 = w0 === nothing ? (NaN, 0) : (genvar(w0), w0.N)
    out = String[]
    for (i, k) in enumerate(rates(set, dem, c))
        r = w0 === nothing ? (alive = false, text = 0) : branch(w0, k, seed + 1000i)
        push!(out, join((i, k, r.alive, r.text), ':'))
    end
    (; set, dem, c, rep, vg0, n0, runs = join(out, ';'))
end

jobs(nrep) = [(; set, dem, c, rep) for (set, dem, c) in cells() for rep in 1:nrep]

function branches(rows)
    out = NamedTuple[]
    for r in rows, s in split(r.runs, ';')
        f = split(s, ':')
        push!(out, (set = Symbol(r.set), dem = Symbol(r.dem), c = Float64(r.c), rep = r.rep, vg0 = Float64(r.vg0), n0 = r.n0,
                    k = parse(Float64, f[2]), alive = f[3] == "true", text = parse(Int, f[4])))
    end
    out
end

# --- analysis ---------------------------------------------------------------------------

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

"Half-survival rate of branches `rs` (fields `k`, `alive`); NaN unless bracketed."
function k50(rs)
    isempty(rs) && return NaN
    y = [r.alive for r in rs]
    (all(y) || !any(y)) && return NaN
    a, c = logistic([log(r.k) for r in rs], y)
    k = c < 0 ? exp(-a / c) : NaN
    ks = [r.k for r in rs]
    minimum(ks) <= k <= maximum(ks) ? k : NaN
end

quant(v, q) = (w = sort(filter(isfinite, v)); isempty(w) ? NaN : w[clamp(round(Int, q * length(w) + 0.5), 1, length(w))])
ci(v) = (quant(v, 0.025), quant(v, 0.975))
byrep(v) = (d = Dict{Int,Vector{NamedTuple}}(); for x in v; push!(get!(d, x.rep, NamedTuple[]), x); end; d)
resample(rng, d) = (ks = collect(keys(d)); reduce(vcat, (d[r] for r in rand(rng, ks, length(ks))); init = NamedTuple[]))
fitM(cs, ρ) = (ok = isfinite.(ρ); sum(cs[ok] .^ 2) / sum(cs[ok] .* (1 .- ρ[ok])))
function fitsqrt(cs, ρ)                       # ρ = √(1 − c/M′): least squares over M′ by a scan
    ok = isfinite.(ρ); best = (Inf, NaN)
    for M in exp.(range(log(0.05), log(50), length = 4000))
        sse = sum((sqrt.(max.(1 .- cs[ok] ./ M, 0)) .- ρ[ok]) .^ 2)
        sse < best[1] && (best = (sse, M))
    end
    best
end

"The predictions of half_survival_rate for every cell, cached (they take about a minute each)."
function predictions(path)
    cache = Dict{Tuple{Symbol,Symbol,Float64},Float64}()
    if isfile(path)
        for (i, line) in enumerate(eachline(path))
            i == 1 && continue
            f = split(line, ',')
            cache[(Symbol(f[1]), Symbol(f[2]), parse(Float64, f[3]))] = parse(Float64, f[4])
        end
    end
    todo = [x for x in cells() if !haskey(cache, x)]
    res = Vector{Float64}(undef, length(todo))
    Threads.@threads :greedy for i in eachindex(todo)
        set, dem, c = todo[i]
        res[i] = half_survival_rate(demp(set, dem); cost = c).k50
    end
    for (x, v) in zip(todo, res)
        cache[x] = v
    end
    open(path, "w") do io
        println(io, "set,dem,c,k50_pred")
        for x in cells()
            println(io, x[1], ',', x[2], ',', x[3], ',', cache[x])
        end
    end
    cache
end

function analyse(rows, predpath, seed0)
    br = branches(rows)
    rng = MersenneTwister(seed0 + 99)
    hsr = predictions(predpath)
    pass = Dict(h => [0, 0] for h in (:C1, :C2, :C3, :old))
    tally(h, ok) = (pass[h][1] += ok; pass[h][2] += 1; ok ? "PASS" : "FAIL")
    nlin = 0
    nsqrt = 0
    for set in keys(SETS), dem in DEMS
        p0 = demp(set, dem)
        dc0 = drift_critical_rate(p0)
        cs = collect(COSTS[dem])
        cellb = Dict(c => byrep([x for x in br if x.set === set && x.dem === dem && x.c == c]) for c in cs)
        any(isempty, values(cellb)) && continue
        k = Dict(c => k50(reduce(vcat, values(cellb[c]); init = NamedTuple[])) for c in cs)
        boot = [Dict(c => k50(resample(rng, cellb[c])) for c in cs) for _ in 1:1000]
        vg = Dict(c => mean(filter(isfinite, [first(v).vg0 for v in values(cellb[c])])) for c in cs)
        nmove = Dict(c => mean([first(v).n0 for v in values(cellb[c])]) for c in cs)
        @printf("\n══ %s, %s ══  drift_critical_rate (no cost): k_c %.5f at λ* %.3f, V_G* %.4f, n* %.0f; M = λ*²/s %.3f\n",
                set, dem, dc0.k, dc0.λ, dc0.VG, dc0.n, dc0.M)
        ρobs = Float64[]
        for c in cs
            pc = withcost(p0, c)
            kc = drift_critical_rate(p0; cost = c).k
            sr = sustainable_rate(pc; VG = vg[0.0])
            kp = hsr[(set, dem, c)]
            kci = ci([b[c] for b in boot])
            ks = sort(unique(x.k for v in values(cellb[c]) for x in v))
            surv = [mean(x.alive for v in values(cellb[c]) for x in v if x.k == kk) for kk in ks]
            @printf("  c = %.2f: V_G at the move %.4f, N %.0f; survival by rate %s\n", c, vg[c], nmove[c],
                    join((@sprintf("%.2f", s) for s in surv), " "))
            @printf("           rates %s\n", join((@sprintf("%.4g", kk) for kk in ks), " "))
            r3 = k[c] / kp
            @printf("    k50 %.5f [%.5f, %.5f]; half_survival_rate %.5f (obs/pred %.3f)  C3 %s;  k50/k_c %.3f, k50/sustainable_rate %.3f\n",
                    k[c], kci..., kp, r3, tally(:C3, isfinite(r3) && 0.75 <= r3 <= 1.33), k[c] / kc, k[c] / sr)
            if c > 0
                ρ = k[c] / k[0.0]
                ρci = ci([b[c] / b[0.0] for b in boot])
                ρp = kc / dc0.k
                ρold = sr / sustainable_rate(p0; VG = vg[0.0])
                push!(ρobs, ρ)
                @printf("    ratio k50(c)/k50(0) %.3f [%.3f, %.3f]; drift_critical_rate %.3f (diff %+.3f)  C1 %s; linear 1 − c/M %.3f; sustainable_rate %.3f (diff %+.3f) %s\n",
                        ρ, ρci..., ρp, ρ - ρp, tally(:C1, isfinite(ρ) && abs(ρ - ρp) <= 0.06), 1 - c / dc0.M, ρold, ρ - ρold,
                        tally(:old, isfinite(ρ) && abs(ρ - ρold) <= 0.06))
            end
        end
        cp = cs[2:end]
        Mobs = fitM(cp, ρobs)
        Mb = [fitM(cp, [b[c] / b[0.0] for c in cp]) for b in boot]
        sl = sum(((1 .- cp ./ Mobs) .- ρobs) .^ 2)
        ss, Ms = fitsqrt(cp, ρobs)
        lin = sl < ss
        nlin += lin; nsqrt += !lin
        @printf("  M_obs %.3f [%.3f, %.3f] vs M_theory %.3f (ratio %.3f)  C2 %s; SSE linear %.2e, square root %.2e (M′ %.3f): %s\n",
                Mobs, ci(Mb)..., dc0.M, Mobs / dc0.M, tally(:C2, isfinite(Mobs) && 0.8 <= Mobs / dc0.M <= 1.25), sl, ss, Ms,
                lin ? "linear" : "square root")
    end
    println()
    for h in (:C1, :C2, :C3)
        @printf("%s: %d of %d pass\n", h, pass[h]...)
    end
    @printf("(sustainable_rate's ratio within 0.06: %d of %d)\n", pass[:old]...)
    v1 = pass[:C1][1] >= 20 ? "SUPPORTED" : pass[:C1][1] < 12 ? "REFUTED" : "MIXED"
    v2 = pass[:C2][1] >= 6 && nlin >= 6 ? "SUPPORTED" : nsqrt >= 6 ? "REFUTED" : "MIXED"
    v3 = pass[:C3][1] >= 26 ? "SUPPORTED" : pass[:C3][1] < 16 ? "REFUTED" : "MIXED"
    println("\nVerdict, relative effect of a cost (drift_critical_rate): $v1")
    println("Verdict, linear-margin form (C2 $(pass[:C2][1])/8; linear form better in $nlin of 8): $v2")
    println("Verdict, half_survival_rate: $v3")
    all(h -> pass[h][1] == pass[h][2], (:C1, :C2, :C3))
end

# --- post hoc: existing data only ----------------------------------------------------------

# k50 of plasticity_erosion.log (run 1): arm 0 (no cost) and arm A (cost before and after)
# in effective units, with c = b²/(2·0.5²)
const EROSION = [(:L100, :fixed, 0.022342, [(0.5, 0.016779), (1.125, 0.011608)]),
                 (:L100, :R2, 0.014624, [(0.125, 0.01147), (0.5, 0.0027514)]),
                 (:L500, :fixed, 0.051718, [(0.5, 0.038191), (1.125, 0.025801)]),
                 (:L500, :R2, 0.03784, [(0.125, 0.029976), (0.5, 0.0069727)])]

function posthoc()
    println("POST HOC (not pre-specified): existing data, no new runs")
    println("\n1. plasticity_erosion.log's cells (arm 0 and arm A, k50 in effective units)")
    for (set, dem, k0, costs) in EROSION
        p = demp(set, dem)
        d0 = drift_critical_rate(p)
        h0 = half_survival_rate(p).k50
        @printf("  %s %s: k50 %.5f; drift_critical_rate %.5f (k50/k_c %.3f, M %.3f); half_survival_rate %.5f (k50/pred %.3f)\n",
                set, dem, k0, d0.k, k0 / d0.k, d0.M, h0, k0 / h0)
        for (c, kc) in costs
            ρ = kc / k0
            ρp = drift_critical_rate(p; cost = c).k / d0.k
            hc = half_survival_rate(p; cost = c).k50
            @printf("    c = %.3f: ratio %.3f; drift_critical_rate %.3f (diff %+.3f); linear 1 − c/M %.3f; M implied by the data %.3f; half_survival_rate ratio %.3f, k50/pred %.3f\n",
                    c, ρ, ρp, ρ - ρp, 1 - c / d0.M, c / (1 - ρ), hc / h0, kc / hc)
        end
    end
    posthoc_populations()
end

# vignettes 09 (fixed size) and 11 (Beverton–Holt) reruns of the paper's 710 populations:
# survival 2000 generations after the move against the rate in units of drift_critical_rate
function posthoc_populations()
    V = normpath(joinpath(@__DIR__, "..", "vignettes"))
    master = normpath(joinpath(@__DIR__, "..", "..", "Repository", "Results", "Master Output file.xlsx"))
    isfile(master) || (println("\n2. skipped: no ", master); return)
    dt = XLSX.readtable(master, "Master output file")
    lab = strip.(string.(dt.column_labels))
    num(v) = v isa AbstractString ? parse(Float64, v) : v === missing ? NaN : Float64(v)
    f(name) = num.(dt.data[findfirst(==(name), lab)])
    t = (; L = round.(Int, f("No. Loci")), u = f("recombination u=3,5"), h2 = f("Heritability"), VP = f("PhenotypeVar"),
         m = f("MutationScalar Vm/Ve"), d = f("d_env"), ω = f("w, fitness function sd"))
    P(i, dem) = Params(; PopSize = 1000, Ploidy = 2, NoQTL = t.L[i], PhenotypeVar = t.VP[i], Heritability = t.h2[i],
                       MutationScalar = t.m[i], MutationRatePerGamete = 1e-4, RecombPoissonMean = t.u[i],
                       InitFitnessSD = t.ω[i], MoveType = :continuous, dem...)
    x50(x, y) = begin ok = isfinite.(x); (a, c) = logistic(log.(x[ok]), y[ok]); c < 0 ? exp(-a / c) : NaN end
    println("\n2. vignettes 09 and 11: the paper's 710 populations, alive 2000 generations after the move; x50 is the",
            " rate at which half survive in units of the prediction (1 = exact)")
    sets = Any[("fixed", joinpath(V, "09_replicating_pagel_2025", "replication_runs.csv"), (Demography = :fixed,))]
    for R in (5.0, 2.0, 1.5)
        push!(sets, ("R = $R", joinpath(V, "11_demography_in_the_papers_populations", "demography_runs_R$R.csv"),
                     (Demography = :beverton_holt, GrowthRate = R)))
    end
    for (name, path, dem) in sets
        d, h = readdlm(path, ','; header = true); h = vec(h)
        col(n) = d[:, findfirst(==(n), h)]
        rows, alive = Int.(col("row")), string.(col("extinct")) .!= "true"
        x, xo, Ls = Float64[], Float64[], Int[]
        for i in rows
            p = P(i, dem)
            push!(x, t.d[i] / drift_critical_rate(p).k)
            push!(xo, t.d[i] / sustainable_rate(p; VG = house_of_gauss_vg(p)))
            push!(Ls, t.L[i])
        end
        @printf("  %-8s %d populations, %.1f%% alive: x50 %.3f (drift_critical_rate), %.3f (sustainable_rate)\n", name,
                length(rows), 100mean(alive), x50(x, alive), x50(xo, alive))
        for L in (100, 500, 1000, 5000)
            s = Ls .== L
            @printf("      L = %4d: x50 %.3f, %.3f\n", L, x50(x[s], alive[s]), x50(xo[s], alive[s]))
        end
    end
end

function main()
    !isempty(ARGS) && ARGS[1] == "posthoc" && return posthoc()
    seed0 = isempty(ARGS) ? SEED0 : parse(Int, ARGS[1])
    tag = seed0 == SEED0 ? "" : "_seed$(seed0)"
    path = joinpath(@__DIR__, "critical_rate$(tag).csv")
    println("Critical rate under drift and a cost of fitness: kernel, $(NREP) burn-ins per cell, seeds from $seed0, ",
            "$(Threads.nthreads()) threads")
    t0 = time()
    rows = sweep(j -> job(j, seed0), jobs(NREP); cache = path, id = j -> "$(j.set)-$(j.dem)-$(j.c)-$(j.rep)",
                 cost = j -> (j.set === :L500 ? 2 : 1), progress = 40)
    @printf("%d burn-ins in %s (%.1f min this session)\n", length(rows), basename(path), (time() - t0) / 60)
    ok = analyse(rows, joinpath(@__DIR__, "critical_rate_pred.csv"), seed0)
    println(ok ? "\nALL CRITERIA PASS" : "\nSOME CRITERIA FAIL")
    ok || exit(1)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
