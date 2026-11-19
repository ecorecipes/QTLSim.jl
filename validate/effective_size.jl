# Is `effective_size(N, n)` the effective size of QTSim's life cycle? And does
# `drift_critical_rate` find every steady state?
#
# Criteria: validate/effective_size_criteria.md (committed before any kernel run).
# Findings: validate/effective_size.md.
#
#   julia -t 4 --project=vignettes validate/effective_size.jl kernel   # part 1a: neutral kernel runs → effective_size.log
#   julia -t 4 --project=vignettes validate/effective_size.jl branch   # part 1b: drift_critical_rate's steady states → effective_size_branch.log
#   julia -t 4 --project=vignettes validate/effective_size.jl impact   # impact of ne = :breeders → effective_size_impact.log
#
# The kernel runs are cached in validate/effective_size_runs.csv (resumable).

using QTLSim, Statistics, Random, Printf, Distributions, LinearAlgebra
import Catlab

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "..", ".."))
const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s); flush(stdout))
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))
const PEAK = QTLSim.QTSIM_FITNESS_PEAK

# ------------------------------------------------------------------------------------------
# Part 1a: neutral kernel runs

const NS = (200, 1000)
const WS = (1.0, 0.75, 0.5, 0.25, 0.1)
const SEED0 = 77_000_000
const REPS = (decay = 300, balance = 60, coal = 50)   # amendment 1
const OFFSET = (decay = 0, balance = 100_000, coal = 200_000)
const VM = 5e-4

"A neutral population in which every adult survives with probability W (criteria file)."
function neutral(N, W; L = 50, mutation = true, gens = 1, sample = 1)
    kw = W == 1 ? (; Drift = true) :
         (; Plasticity = 1.0, OptimumSlope = 1e12, PlasticityCostSD = 1 / sqrt(2 * log(PEAK / W)))
    Params(; PopSize = N, Ploidy = 2, NoQTL = L, CrossoverModel = :all, RecombPoissonMean = 10.0, PhenotypeVar = 1.0,
           Heritability = 0.5, MutationScalar = mutation ? 1e-3 : 0.0, MutationRatePerGamete = 0.01, InitFitnessSD = 1e6,
           MoveType = :continuous, MoveFitnessSD = 0.0, NoGenerations = gens, Sample = sample, kw...)
end

"Harmonic mean of the breeders, n ~ Binomial(N, W) given n ≥ 2 (n = N when W = 1)."
function harmonic_breeders(N, W)
    W == 1 && return float(N)
    d = Binomial(N, W)
    num = sum(pdf(d, n) / n for n in 2:N); den = sum(pdf(d, n) for n in 2:N)
    den / num
end
"N_e of the decay of variance: 1/(2(1 − λ)), λ the largest root of λ² − (1 − ε)λ − ε/2."
eigen_size(nh) = (ε = 1 / nh; 1 / (2 * (1 - ((1 - ε) + sqrt(1 + ε^2)) / 2)))

genic(G) = sum(var(vec(view(G, l, :, :))) for l in axes(G, 1))

function decay_run(N, W, seed)
    nbar = N * W
    T = round(Int, nbar)
    p = neutral(N, W; mutation = false)
    rng = Xoshiro(seed)
    w = World(p)
    setgenomes!(w, randn(rng, p.NoQTL, 2, N))
    init_env!(w, rng)
    for _ in 1:10; generation!(w, rng); end
    H0 = genic(genomes(w))
    for _ in 1:T; generation!(w, rng); end
    (; a = H0, b = genic(genomes(w)), T, extinct = w.extinct)
end

function balance_run(N, W, seed)
    nbar = N * W
    burn = round(Int, 12nbar)
    p = neutral(N, W; gens = 2burn, sample = max(1, round(Int, nbar / 20)))
    o = simulate(p; seed)
    keep = o.gen .>= burn
    (; a = mean(o.geno_var[keep]), b = float(count(keep)), T = burn, extinct = o.extinct)
end

function coal_run(N, W, seed)
    nbar = N * W
    T = round(Int, 16nbar)
    p = neutral(N, W; L = 10, mutation = false)
    rng = Xoshiro(seed)
    w = World(p)
    ped = record_pedigree!(w; simplify_every = 100, mutations = false)
    init_env!(w, rng)
    for _ in 1:T; generation!(w, rng); end
    pc = pair_coalescence(ped)
    m = coalescence_mean(pc)
    (; a = m.S[1, 1], b = m.W[1], T = pc.t, extinct = w.extinct)
end

const RUN = (decay = decay_run, balance = balance_run, coal = coal_run)
cellno(N, W) = 10 * findfirst(==(N), NS) + findfirst(==(W), WS)

# amendment 2: the cells dropped for cost
dropped(what, N, W) = N == 1000 && ((what === :coal && W > 0.5) || (what === :balance && W == 1))
function kernel_jobs(reps_of)
    [(; what, N, W, rep) for what in keys(RUN) for N in NS for W in WS for rep in 1:reps_of(what, N, W) if !dropped(what, N, W)]
end

function run_kernel(jobs)
    sweep(jobs; cache = joinpath(HERE, "effective_size_runs.csv"), id = j -> "$(j.what)-$(j.N)-$(j.W)-$(j.rep)",
          cost = j -> j.N * j.N * j.W * (j.what === :decay ? 1 : j.what === :balance ? 24 : 60), progress = 500) do j
        r = RUN[j.what](j.N, j.W, SEED0 + OFFSET[j.what] + 1000 * cellno(j.N, j.W) + j.rep)
        (; what = String(j.what), j.N, j.W, j.rep, r.a, r.b, r.T, r.extinct)
    end
end

"N̂_e and its standard error for one cell and measurement from its rows."
function estimate(what, rows)
    R = length(rows)
    a = [r.a for r in rows]
    if what === :decay
        b = [r.b for r in rows]; T = rows[1].T
        r = mean(b) / mean(a)
        vln = (var(b) / mean(b)^2 + var(a) / mean(a)^2 - 2cov(a, b) / (mean(a) * mean(b))) / R
        λ = r^(1 / T)
        ne = 1 / (2 * (1 - λ))
        (ne, λ / (2 * (1 - λ)^2 * T) * sqrt(vln))
    elseif what === :balance
        (mean(a) / (2VM), std(a) / sqrt(R) / (2VM))
    else
        ((mean(a) - 1) / 2, std(a) / sqrt(R) / 2)
    end
end

function kernel_part()
    say("validate/effective_size.jl kernel — neutral kernel runs against three formulas for N_e (criteria: effective_size_criteria.md)")
    say("threads: ", Threads.nthreads())
    cpu0 = time()
    rows = run_kernel(kernel_jobs((what, N, W) -> REPS[what]))
    # rows of cells dropped by amendment 2 that had finished: read from the cache by hand
    partial = NamedTuple[]
    for ln in eachline(joinpath(HERE, "effective_size_runs.csv"))
        f = split(ln, ',')
        (length(f) == 9 && f[1] != "id") || continue
        r = (what = String(f[2]), N = parse(Int, f[3]), W = parse(Float64, f[4]), rep = parse(Int, f[5]), a = parse(Float64, f[6]),
             b = parse(Float64, f[7]), T = parse(Int, f[8]), extinct = parse(Bool, f[9]))
        dropped(Symbol(r.what), r.N, r.W) && push!(partial, r)
    end
    pick(rows, what, N, W) = [r for r in rows if r.what == String(what) && r.N == N && r.W == W]
    # the pre-specified doubling of cells whose standard error exceeds 3%
    doubled = Set{Tuple{Symbol,Int,Float64}}()
    for what in keys(RUN), N in NS, W in WS
        dropped(what, N, W) && continue
        ne, se = estimate(what, pick(rows, what, N, W))
        se / ne > 0.03 && push!(doubled, (what, N, W))
    end
    if !isempty(doubled)
        rows = run_kernel(kernel_jobs((what, N, W) -> ((what, N, W) in doubled ? 2 : 1) * REPS[what]))
    end
    say(@sprintf("%d runs (%d extinct); cells doubled by the 3%% rule: %s", length(rows), count(r -> r.extinct, rows),
                 isempty(doubled) ? "none" : join(("$(d[1]) N=$(d[2]) W=$(d[3])" for d in sort(collect(doubled))), ", ")))
    say("\nFormulas: current = effective_size(N, n̄), n̄ = N·W; derived = n_h (harmonic mean of the binomial breeders) for the balance and")
    say("coalescence measurements and the eigenvalue size ≈ n_h + ½ for the decay. agree: |N̂_e − F| ≤ max(0.05 F, 2 SE); reject: > 5% and > 3 SE.")
    verdict = Dict(f => (agree = Ref(0), reject = Ref(0)) for f in (:current, :nbar, :derived))
    names = (decay = "(i) decay of variance, no mutation", balance = "(ii) mutation–drift balance, E[V_G]/(2V_m)",
             coal = "(iii) pair coalescence, (mean time − 1)/2")
    for what in keys(RUN)
        say("\n", names[what])
        sayf("  %5s %5s %5s | %9s %7s %6s | %9s %7s %4s | %9s %7s %4s | %9s %7s %4s", "N", "W", "reps", "N̂_e", "SE", "SE %",
             "current", "ratio", "", "n̄", "ratio", "", "derived", "ratio", "")
        for N in NS, W in WS
            isdrop = dropped(what, N, W)
            rs = pick(isdrop ? partial : rows, what, N, W)
            length(rs) < 2 && continue
            ne, se = estimate(what, rs)
            nh = harmonic_breeders(N, W)
            F = (current = effective_size(N, N * W), nbar = N * W, derived = what === :decay ? eigen_size(nh) : nh)
            mark = Dict{Symbol,String}()
            for f in keys(F)
                d = abs(ne - F[f])
                ag = d <= max(0.05 * F[f], 2se); rej = d > 0.05 * F[f] && d > 3se
                isdrop || (verdict[f].agree[] += ag; verdict[f].reject[] += rej)
                mark[f] = (ag ? "ok" : rej ? "REJ" : "?") * (isdrop ? "†" : "")
            end
            sayf("  %5d %5.2f %5d | %9.2f %7.2f %6.2f | %9.2f %7.3f %4s | %9.2f %7.3f %4s | %9.2f %7.3f %4s", N, W, length(rs), ne, se,
                 100se / ne, F.current, ne / F.current, mark[:current], F.nbar, ne / F.nbar, mark[:nbar], F.derived, ne / F.derived, mark[:derived])
        end
    end
    say("  († a cell dropped by amendment 2: the runs that had finished, outside the verdict)")
    say("\nVerdicts (27 cell × measurement combinations each; amendment 2)")
    for (f, name) in ((:current, "current  effective_size(N, n̄)"), (:nbar, "n̄ = N·W"), (:derived, "derived  n_h (n_h + ½ for the decay)"))
        a, r = verdict[f].agree[], verdict[f].reject[]
        sayf("  %-40s agrees in %2d of 27, rejected in %2d: %s", name, a, r, a == 27 ? "FITS" : r > 0 ? "REJECTED" : "INCONCLUSIVE")
    end
    say("\nThe side observation of validate/xu_shen.md, in the kernel (N = 200): W = 1, 0.5, 0.25 ↔ n̄ = 200, 100, 50; effective_size = ",
        join((@sprintf("%.0f", effective_size(200, n)) for n in (200, 100, 50)), ", "), ".")
    sayf("wall time this session %.1f min", (time() - cpu0) / 60)
    write(joinpath(HERE, "effective_size.log"), String(take!(out)))
end

# parts 1b and impact
include(joinpath(HERE, "effective_size_parts.jl"))

function main()
    part = isempty(ARGS) ? "kernel" : ARGS[1]
    part == "kernel" ? kernel_part() : part == "branch" ? branch_part() : part == "impact" ? impact_part() :
    error("unknown part $part")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
