# Populations initialised from theory (src/init.jl: `equilibrium_world`, `simulate(p; init =
# :equilibrium, burnin)`) against populations burned in from zero genomes: how many
# generations of burn-in B does a population drawn from the stationary allele-frequency
# spectrum need before it is indistinguishable, before and after the optimum starts to
# move, from one that ran the protocol's 10,000 generations, and what does that save?
#
# Usage (from the package root; HypothesisTests and XLSX come from the vignettes environment):
#     julia -t 3 --project=vignettes validate/init_equilibrium.jl [R]
# R replicates per arm and cell (default 200, at most 2000). Writes validate/init_equilibrium.log and the
# per-run table validate/init_equilibrium_runs.csv, from which an interrupted run resumes
# (with ENV SEED_OFFSET ≠ 0: init_equilibrium_seed<SEED_OFFSET>.log and _runs_seed<…>.csv).
# SMOKE=1 runs a short version of every path and writes nothing. The PUB part reads
# Repository/Simulations (read only, as validate/moments_regimes.jl does); it is skipped
# without it. Exits 1 if no burn-in of the grid passes in some cell.
#
# CELLS (fixed before the first run). Diploid, N = 1000, `CrossoverModel = :qtsim`, the
# optimum fixed at 0 until generation T0 = 10,000 and then moving continuously
# (`MoveType = :continuous`), as in the published sweep.
# - P100:  the paper's regime, 100 loci. V_P = 0.0348, h² = 0.2325, m = 0.00268, μ = 1e-4
#          (4Nμ = 0.4), ω = 5.13, RecombPoissonMean 3: the values of one published input
#          (QT0100-Possion3.0, the first file), drift-dominated (V_G ≈ 2N_e V_m);
# - P1000: the same with 1000 loci and RecombPoissonMean 5;
# - G:     the Gaussian regime: V_P = 1, h² = 0.5, m = 0.01, 100 loci, μ = 0.01 (4Nμ = 40),
#          the analytical ω = 7.07, RecombPoissonMean 2 (a cell of validate/moments_regimes.jl);
# - H:     the house-of-cards regime with weak selection per mutation: the same with
#          μ = 2.5e-5 (4Nμ = 0.1; α² = 1, α²/(2V_s) = 0.01, Nα²/V_s = 20).
#
# ARMS. REF: zero genomes, generations 0:T0−1 (the usual burn-in). B = 0, 100, 300, 1000,
# 3000: `equilibrium_world(p, rng; burnin = B, at = T0)` with the defaults (`target =
# :genic`, `mean = :stationary`), then generations T0−B:T0−1. At generation T0 every
# population is measured and then continued twice from a copy, with fresh generators:
# - "fast": the optimum moves by k_fast per generation, to extinction or at most CAP =
#   2000 generations (the cap of the sweeps this is meant for);
# - "slow": it moves by k_slow per generation for 1000 generations.
# Rates (chosen from pilot runs of the REF arm only, seeds 9,400,001–9,500,009, so that the
# fast move kills within the cap and the slow one mostly does not; not tuned on any
# comparison of arms): P100 and P1000 k_fast = 0.03 (pilot: extinct after about 850–1060
# generations), k_slow = 0.00492 (the published input's own rate; all 9 pilot runs alive
# at 2000); G 0.3 (extinct after 300–350) and 0.1 (not piloted); H 0.5 (not piloted; at
# 0.4 eight of nine survived 2000 generations, at 0.7 all died within 121) and 0.1.
#
# MEASURES, per population:
#   vg0, genic0, ld0   the variance of the genotypes at generation T0 (before selection),
#                      its genic part Σ_loci Var(sum of the two alleles) and ld0 = vg0 − genic0;
#   zbar0              the mean genotype at T0 (the optimum is 0);
#   text               fast move: generations from T0 to extinction (CAP if alive: censored);
#   adv200             slow move: (mean genotype at T0+200 − zbar0)/200, the rate of adaptation;
#   vg200, vg1000      slow move: the variance of the genotypes at T0+200 and T0+1000;
#   lag1000            slow move: θ − mean genotype at T0+1000;
#   and whether it was alive at the fast move's cap and at the slow move's end.
#
# PRE-SPECIFIED CRITERIA (written before the first run of this script; R = 200 per arm).
# For each cell and burn-in B, against REF:
# - KS: two-sample Kolmogorov–Smirnov tests (asymptotic p) of vg0, genic0, ld0, zbar0, text,
#   adv200, vg200, vg1000, lag1000 (the slow measures over the populations alive then);
# - ratio: the ratio of means B/REF with its 95% interval (delta method) for vg0, genic0,
#   text, adv200, vg200, vg1000, lag1000;
# - survival: Fisher's exact test of the numbers alive at the fast cap and at the slow end.
# B FAILS in a cell if any KS or Fisher p < 0.01, or if any ratio's whole 95% interval lies
# outside [0.95, 1.05] (a difference of more than 5% that is not noise). B PASSES if it
# does not fail and every ratio's 95% interval lies inside [0.90, 1.10] (equivalence within
# 10% shown). Otherwise it is INCONCLUSIVE: no difference detected, equivalence not shown
# at this R (expected in H, where V_G varies by a third among populations).
# The cell's answer is the smallest B from which no larger B of the grid fails.
# Power: with 200 per arm the KS test at α = 0.01 detects a shift of the mean by 0.4 SD,
# or a standard deviation 0.6 times the reference's, about four times in five (the script
# prints the power by simulation); the ratio interval's half-width is 1.96·CV·√(2/R), 2.5%
# for CV = 0.13 (P1000's vg0) and 6.5% for CV = 0.33 (H). With eleven tests at 0.01 a
# (cell, B) that is truly equivalent fails by chance about one time in ten, so about two of
# the twenty are expected to: a failure that larger B do not share and that a fresh-seed run
# (SEED_OFFSET) does not repeat is reported as such.
# - PUB (the paper's populations; needs Repository): 2 published inputs drawn at random
#   from each of the 8 strata (16 inputs, their own parameters and move), 10 populations
#   each of REF, B = 300 and B = 1000, run to CAP after the move. Per input the ratio of
#   the arms' mean vg0 and mean text (censored at CAP); pooled as the geometric mean over
#   inputs with a t interval over inputs. Same FAIL/PASS/INCONCLUSIVE rule on the two
#   pooled ratios (no KS: 10 runs per input).
# - TIME: wall time of one run, single-threaded, of the whole protocol to T0 + CAP without
#   extinction, from zero and with each B (sampling included), 2 runs each: the saving for
#   a sweep capped 2000 generations after the move (runs that die earlier save more).
#
# DISCLOSURE. The sampler was written from the theory (src/init.jl) and then piloted
# (seeds 9,000,001–9,300,120; 18–60 populations per arm, the four cells above with H first
# at μ = 2.5e-6, which is outside the weak-selection regime and was replaced by the
# moments_regimes cell before any run of it): the pilots showed that the sampler's genic
# variance matches the reference's within their noise in every cell, that V_G starts above
# the reference's by the missing linkage disequilibrium (50% in G, 10% in P1000) and comes
# down with the burn-in, and that the spread of V_G among populations is too small at B = 0
# in P1000 (CV 0.07 against 0.13) and recovers within about 100–300 generations. Nothing in
# the sampler was changed after the pilots (the draw of the mean genotype from its
# stationary law was in it before the first pilot). The pilots led to adding B = 3000 to the
# suggested grid {0, 100, 300, 1000}. No run of this script had been made when the criteria
# were written.
#
# HISTORY. Run 1 (SEED_OFFSET = 0, R = 200, 3 threads) started 2026-11-10 15:30 at commit
# 33838d6 and was killed by a session restart on 2026-11-11 (last row 11:08) with 2,552 of
# the 4,800 cell rows in the cache (REF and B = 3000 of every cell, all of P1000, 152 of
# P100's B = 1000). It was resumed from the cache on 2026-11-11 with the same code, grid,
# criteria and seeds (the cache is keyed by job id and each job's seed is a function of its
# id, so the rows are those of an uninterrupted run; two cached rows, P100|1000|3 and
# P100|3000|7, were recomputed before resuming and agreed bit for bit). The only change
# was the log line counting cached rows.

using QTLSim, Statistics, Printf, Random
using HypothesisTests: ApproximateTwoSampleKSTest, FisherExactTest, pvalue
import Distributions

const SMOKE = get(ENV, "SMOKE", "") == "1"
const R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : (SMOKE ? 6 : 200)
R <= 2000 || error("at most 2000 replicates per arm (the seeds of the branches)")
const SEED_OFFSET = parse(Int, get(ENV, "SEED_OFFSET", "0"))
const TAG = SEED_OFFSET == 0 ? "" : "_seed$(SEED_OFFSET)"
const LOGPATH = SMOKE ? tempname() : joinpath(@__DIR__, "init_equilibrium$(TAG).log")
const CACHE = SMOKE ? tempname() : joinpath(@__DIR__, "init_equilibrium_runs$(TAG).csv")
const PUBCACHE = SMOKE ? tempname() : joinpath(@__DIR__, "init_equilibrium_pub_runs$(TAG).csv")
io = open(LOGPATH, "w")
say(s...) = (println(stdout, s...); println(io, s...); flush(io); flush(stdout))

const T0 = SMOKE ? 300 : 10_000
const BS = SMOKE ? (0, 50) : (0, 100, 300, 1000, 3000)
const CAP = SMOKE ? 150 : 2000
const SLOW = SMOKE ? 100 : 1000
const MID = SMOKE ? 20 : 200
const PUBBS = SMOKE ? (50,) : (300, 1000)
const NPUB, RPUB = SMOKE ? (1, 2) : (2, 10)

base(; kw...) = Params(; Ploidy = 2, PopSize = 1000, Sample = 1, NoGenerations = T0 + CAP, MoveFitnessSDGen = T0,
                       MoveType = :continuous, kw...)
paper(; kw...) = base(; PhenotypeVar = 0.0348, Heritability = 0.2325, MutationScalar = 0.00268, MutationRatePerGamete = 1e-4,
                      InitFitnessSD = 5.13, kw...)
const CELLS = [
    (name = "P100", p = paper(NoQTL = 100, RecombPoissonMean = 3.0), kfast = 0.03, kslow = 0.00492),
    (name = "P1000", p = paper(NoQTL = 1000, RecombPoissonMean = 5.0), kfast = 0.03, kslow = 0.00492),
    (name = "G", p = base(NoQTL = 100, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01, MutationRatePerGamete = 0.01,
                          RecombPoissonMean = 2.0), kfast = 0.3, kslow = 0.1),
    (name = "H", p = base(NoQTL = 100, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01, MutationRatePerGamete = 2.5e-5,
                          RecombPoissonMean = 2.0), kfast = 0.5, kslow = 0.1),
]

# --- one population ---------------------------------------------------------------------

function genic_variance(w::World)
    G = genomes(w)
    sum(var(view(G, l, 1, :) .+ view(G, l, 2, :)) for l in axes(G, 1))
end

# the World at generation T0: burned in from zero (B < 0) or from the stationary spectrum
function world_at_move(p::Params, B::Integer, rng, theory)
    if B < 0
        w = World(p)
        init_env!(w, rng)
    else
        w = equilibrium_world(p, rng; burnin = B, at = T0, target = theory.target, spectrum = theory.st)
    end
    for it in w.gen:T0-1
        w.gen = it
        generation!(w, rng)
        w.extinct && break
    end
    w
end

# generations from T0 on at move rate k, up to n of them; `at(w, t)` is called before generation T0 + t
function move!(w::World, p::Params, k::Real, n::Integer, rng, at = (w, t) -> nothing)
    w.p = Params(p; MoveFitnessSD = k)
    for t in 0:n-1
        at(w, t)
        w.gen = T0 + t
        generation!(w, rng)
        w.extinct && return t
    end
    at(w, n)
    n
end

gmean(w) = (develop!(w); mean(view(w.genotype, 1:w.N)))
gvar(w) = (develop!(w); var(view(w.genotype, 1:w.N)))

function run_one(cell, B::Integer, seed::Integer, theory)
    p = cell.p
    w = world_at_move(p, B, MersenneTwister(seed), theory)
    nan = NaN
    w.extinct && return (; vg0 = nan, genic0 = nan, ld0 = nan, zbar0 = nan, text = -1, fast_alive = false,
                         adv200 = nan, vg200 = nan, vg1000 = nan, lag1000 = nan, slow_alive = false)
    vg0, zbar0, genic0 = gvar(w), gmean(w), genic_variance(w)
    text = move!(deepcopy(w), p, cell.kfast, CAP, MersenneTwister(seed + 5000))
    adv200 = vg200 = vg1000 = lag1000 = nan
    v = deepcopy(w)
    tslow = move!(v, p, cell.kslow, SLOW, MersenneTwister(seed + 7000), (v, t) -> begin
        if t == MID
            adv200, vg200 = (gmean(v) - zbar0) / MID, gvar(v)
        elseif t == SLOW
            vg1000, lag1000 = gvar(v), v.θ - gmean(v)
        end
    end)
    (; vg0, genic0, ld0 = vg0 - genic0, zbar0, text, fast_alive = text == CAP,
     adv200, vg200, vg1000, lag1000, slow_alive = tslow == SLOW && !v.extinct)
end

# --- statistics -------------------------------------------------------------------------

sem(x) = std(x) / sqrt(length(x))
ks(x, y) = (length(x) < 2 || length(y) < 2) ? NaN : pvalue(ApproximateTwoSampleKSTest(x, y))
function ratio(x, y)              # mean(y)/mean(x) and its 95% interval (delta method)
    (length(x) < 2 || length(y) < 2) && return (NaN, NaN, NaN)
    r = mean(y) / mean(x)
    se = abs(r) * sqrt((sem(y) / mean(y))^2 + (sem(x) / mean(x))^2)
    (r, r - 1.96se, r + 1.96se)
end
fisher(a, na, b, nb) = pvalue(FisherExactTest(a, na - a, b, nb - b))
finite(v) = filter(isfinite, v)

const KSM = (:vg0, :genic0, :ld0, :zbar0, :text, :adv200, :vg200, :vg1000, :lag1000)
const RATIOM = (:vg0, :genic0, :text, :adv200, :vg200, :vg1000, :lag1000)

# FAIL / PASS / INCONCLUSIVE of one arm against the reference, with the reasons
function verdict(ref, arm; ksm = KSM, ratiom = RATIOM, survival = true)
    fails, open = String[], String[]
    col(rows, m) = finite(Float64[getfield(r, m) for r in rows])
    for m in ksm
        pv = ks(col(ref, m), col(arm, m))
        isnan(pv) && (push!(open, "$m: too few alive"); continue)
        pv < 0.01 && push!(fails, @sprintf("KS %s p = %.2g", m, pv))
    end
    for m in ratiom
        r, lo, hi = ratio(col(ref, m), col(arm, m))
        isnan(r) && continue
        (lo > 1.05 || hi < 0.95) && push!(fails, @sprintf("ratio %s %.3f [%.3f, %.3f]", m, r, lo, hi))
        (lo < 0.90 || hi > 1.10) && push!(open, @sprintf("%s [%.3f, %.3f]", m, lo, hi))
    end
    if survival
        for m in (:fast_alive, :slow_alive)
            pv = fisher(count(r -> getfield(r, m), ref), length(ref), count(r -> getfield(r, m), arm), length(arm))
            pv < 0.01 && push!(fails, @sprintf("Fisher %s p = %.2g", m, pv))
        end
    end
    (isempty(fails) ? (isempty(open) ? "PASS" : "INCONCLUSIVE") : "FAIL", fails, open)
end

# the power of the KS test at α = 0.01 with n per arm, by simulation
function ks_power(n; shift = 0.0, scale = 1.0, reps = 2000, rng = MersenneTwister(1))
    count(_ -> ks(randn(rng, n), shift .+ scale .* randn(rng, n)) < 0.01, 1:reps) / reps
end

# --- the cells --------------------------------------------------------------------------

say("QTLSim validate/init_equilibrium.jl: R = $R per arm, T0 = $T0, burn-ins $(BS), cap $CAP, slow move $SLOW generations, ",
    "SEED_OFFSET = $SEED_OFFSET, $(Threads.nthreads()) threads", SMOKE ? " (SMOKE)" : "")
say(@sprintf("power of the KS test at α = 0.01, %d per arm (2000 simulated pairs): shift 0.3 SD %.2f, 0.4 SD %.2f, 0.5 SD %.2f; SD ratio 0.75 %.2f, 0.6 %.2f, 0.5 %.2f; no difference %.3f",
             R, ks_power(R; shift = 0.3), ks_power(R; shift = 0.4), ks_power(R; shift = 0.5), ks_power(R; scale = 0.75),
             ks_power(R; scale = 0.6), ks_power(R; scale = 0.5), ks_power(R)))

theories = map(CELLS) do c
    st = stationary_spectrum(c.p)
    me = multilocus_equilibrium(c.p; drift = true)
    (; st, me, target = QTLSim.init_target(c.p, st, :genic))
end
say("\ntheory per cell: the stationary spectrum's V (linkage equilibrium), the moments equilibrium (rare alleles, drift), the sampler's target and sites")
for (c, th) in zip(CELLS, theories)
    say(@sprintf("  %-6s spectrum V %.5g  moments V_G %.5g genic %.5g ld %.4g  target %.5g  sites %.0f  n %.1f  SD(mean genotype) %.4g  k_fast %g k_slow %g",
                 c.name, th.st.VG, th.me.VG, th.me.genic, th.me.ld, th.target, QTLSim.stationary_sites(th.st), th.st.nparents,
                 sqrt(QTLSim.stationary_mean_var(c.p, th.st; VG = th.target)), c.kfast, c.kslow))
end

const ARMS = (-1, BS...)
armname(B) = B < 0 ? "REF" : "B=$B"
jobs = [(; ci, B, i) for ci in eachindex(CELLS) for B in ARMS for i in 1:R]
jobseed(j) = SEED_OFFSET + 1_000_000 * j.ci + 10_000 * (findfirst(==(j.B), ARMS) - 1) + j.i
jobcost(j) = CELLS[j.ci].p.NoQTL^0.5 * ((j.B < 0 ? T0 : j.B) + CAP + SLOW)
ncached = isfile(CACHE) ? max(0, countlines(CACHE) - 1) : 0
ncached > 0 && say("\nresuming: $ncached rows already in $(basename(CACHE)) (not rerun; same seeds, keyed by job id)")
t0 = time()
rows = sweep(jobs; cache = CACHE, id = j -> "$(CELLS[j.ci].name)|$(j.B)|$(j.i)", cost = jobcost, progress = SMOKE ? 0 : 200) do j
    run_one(CELLS[j.ci], j.B, jobseed(j), theories[j.ci])
end
say(@sprintf("\n%d populations in %.0f s (cached rows are not rerun)", length(jobs), time() - t0))

answers = Dict{String,Any}()
allfails = String[]
for (ci, c) in enumerate(CELLS)
    of(B) = [rows[k] for k in eachindex(jobs) if jobs[k].ci == ci && jobs[k].B == B]
    ref = of(-1)
    say("\n", "="^100, "\n", c.name, ": ", R, " populations per arm")
    say("  arm      dead<T0   vg0 mean ± se (CV)        genic0            ld0        SD(zbar0)   text median [q10, q90]  alive fast  alive slow   adv200      vg200      vg1000     lag1000")
    for B in ARMS
        a = of(B)
        col(m) = finite(Float64[getfield(r, m) for r in a])
        tx = col(:text); tx = tx[tx .>= 0]
        q(v, x) = isempty(v) ? NaN : quantile(v, x)
        mm(v) = isempty(v) ? NaN : mean(v)
        say(@sprintf("  %-8s %4d    %9.4g ± %-8.2g (%.3f) %9.4g ± %-8.2g %9.3g   %8.4g     %6.0f [%5.0f, %5.0f]      %4d       %4d      %9.4g  %9.4g  %9.4g  %9.4g",
                     armname(B), count(r -> r.text < 0, a), mm(col(:vg0)), sem(col(:vg0)), std(col(:vg0)) / mm(col(:vg0)),
                     mm(col(:genic0)), sem(col(:genic0)), mm(col(:ld0)), std(col(:zbar0)), q(tx, 0.5), q(tx, 0.1), q(tx, 0.9),
                     count(r -> r.fast_alive, a), count(r -> r.slow_alive, a), mm(col(:adv200)), mm(col(:vg200)), mm(col(:vg1000)), mm(col(:lag1000))))
    end
    say("  against REF: KS p of each measure; ratio of means [95% interval]")
    verdicts = String[]
    for B in BS
        a = of(B)
        col(rs, m) = finite(Float64[getfield(r, m) for r in rs])
        say("  ", armname(B), ":")
        say("    KS p   ", join((@sprintf("%s %.3g", m, ks(col(ref, m), col(a, m))) for m in KSM), "  "))
        say("    ratio  ", join((begin
                                     r, lo, hi = ratio(col(ref, m), col(a, m))
                                     @sprintf("%s %.3f [%.3f, %.3f]", m, r, lo, hi)
                                 end for m in RATIOM), "  "))
        say(@sprintf("    Fisher p  alive at the fast cap %.3g, alive at the slow end %.3g",
                     fisher(count(r -> r.fast_alive, ref), length(ref), count(r -> r.fast_alive, a), length(a)),
                     fisher(count(r -> r.slow_alive, ref), length(ref), count(r -> r.slow_alive, a), length(a))))
        v, fails, open = verdict(ref, a)
        push!(verdicts, v)
        say("    ", v, isempty(fails) ? "" : ": " * join(fails, "; "), v == "INCONCLUSIVE" ? " (equivalence within 10% not shown: " * join(open, "; ") * ")" : "")
        v == "FAIL" && push!(allfails, "$(c.name) $(armname(B)): " * join(fails, "; "))
    end
    lastfail = findlast(==("FAIL"), verdicts)
    answers[c.name] = lastfail === nothing ? (BS[1], verdicts[1]) : lastfail == length(BS) ? (nothing, "FAIL") : (BS[lastfail+1], verdicts[lastfail+1])
    b, v = answers[c.name]
    say("  ", c.name, ": ", b === nothing ? "every burn-in of the grid fails" : "smallest burn-in from which none fails: B = $b ($v)")
end

# --- PUB: the paper's populations ---------------------------------------------------------

const SIMDIR = joinpath(normpath(joinpath(@__DIR__, "..", "..")), "Repository", "Simulations")
pubverdict = nothing
if isdir(SIMDIR)
    strata = ["QT0100-Possion3.0", "QT0100-Possion5.0", "QT0500-Possion3.0", "QT0500-Possion5.0",
              "QT1000-Possion3.0", "QT1000-Possion5.0", "QT5000-Possion3.0", "QT5000-Possion5.0"]
    rng = MersenneTwister(20261110 + SEED_OFFSET)
    inputs = []
    for dir in strata
        files = sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, dir))))
        for f in sort(shuffle(rng, files)[1:NPUB])
            q = read_qtsim_options(joinpath(SIMDIR, dir, f))
            q = Params(q; Seed = nothing, Sample = 1, NoGenerations = T0 + CAP, MoveFitnessSDGen = T0)
            st = stationary_spectrum(q)
            push!(inputs, (; dir, file = f, p = q, theory = (; st, target = QTLSim.init_target(q, st, :genic))))
        end
    end
    function pub_one(inp, B, seed)
        w = world_at_move(inp.p, B, MersenneTwister(seed), inp.theory)
        w.extinct && return (; vg0 = NaN, text = -1)
        vg0 = gvar(w)
        (; vg0, text = move!(w, inp.p, inp.p.MoveFitnessSD, CAP, MersenneTwister(seed + 5000)))
    end
    parms = (-1, PUBBS...)
    pjobs = [(; ii, B, i) for ii in eachindex(inputs) for B in parms for i in 1:RPUB]
    t0 = time()
    prows = sweep(pjobs; cache = PUBCACHE, id = j -> "$(inputs[j.ii].dir)/$(inputs[j.ii].file)|$(j.B)|$(j.i)",
                  cost = j -> inputs[j.ii].p.NoQTL^0.5 * ((j.B < 0 ? T0 : j.B) + CAP)) do j
        pub_one(inputs[j.ii], j.B, SEED_OFFSET + 50_000_000 + 100_000 * j.ii + 10_000 * (findfirst(==(j.B), parms) - 1) + j.i)
    end
    say("\n", "="^100, @sprintf("\nPUB: %d published inputs × %d populations per arm (%.0f s)", length(inputs), RPUB, time() - t0))
    say("  input                               L    u   d_env      REF vg0   text (alive)  ", join(("B=$B: vg0 ratio, text ratio (alive)" for B in PUBBS), "   "))
    lv = Dict(B => Float64[] for B in PUBBS); lt = Dict(B => Float64[] for B in PUBBS)
    for (ii, inp) in enumerate(inputs)
        of(B) = [prows[k] for k in eachindex(pjobs) if pjobs[k].ii == ii && pjobs[k].B == B]
        ref = of(-1)
        mv(rs) = mean(finite([r.vg0 for r in rs])); mt(rs) = mean(Float64[r.text for r in rs if r.text >= 0])
        line = @sprintf("  %-17s %-16s %5d %3.0f %9.3g %9.4g %7.0f (%2d)", inp.dir, inp.file, inp.p.NoQTL, inp.p.RecombPoissonMean,
                        inp.p.MoveFitnessSD, mv(ref), mt(ref), count(r -> r.text == CAP, ref))
        for B in PUBBS
            a = of(B)
            push!(lv[B], log(mv(a) / mv(ref))); push!(lt[B], log(mt(a) / mt(ref)))
            line *= @sprintf("      %6.3f  %6.3f (%2d)   ", mv(a) / mv(ref), mt(a) / mt(ref), count(r -> r.text == CAP, a))
        end
        say(line)
    end
    function pooled(l)
        n = length(l); m = mean(l); h = Distributions.quantile(Distributions.TDist(n - 1), 0.975) * std(l) / sqrt(n)
        (exp(m), exp(m - h), exp(m + h))
    end
    pubverdict = Dict{Int,String}()
    for B in PUBBS
        rv, rt = pooled(lv[B]), pooled(lt[B])
        fail = any(r -> r[2] > 1.05 || r[3] < 0.95, (rv, rt))
        inside = all(r -> r[2] >= 0.90 && r[3] <= 1.10, (rv, rt))
        pubverdict[B] = fail ? "FAIL" : inside ? "PASS" : "INCONCLUSIVE"
        say(@sprintf("  B = %d: geometric mean over inputs of B/REF: vg0 %.3f [%.3f, %.3f], text %.3f [%.3f, %.3f]: %s",
                     B, rv..., rt..., pubverdict[B]))
        fail && push!(allfails, "PUB B=$B")
    end
else
    say("\nPUB skipped: Repository/Simulations not found")
end

# --- TIME ---------------------------------------------------------------------------------

say("\n", "="^100, "\nTIME: one run of the whole protocol to T0 + CAP = $(T0 + CAP) generations without a move (no extinction), single-threaded, mean of 2")
say("  cell      from zero    ", join((@sprintf("B=%-5d", B) for B in BS), "            "))
for (c, th) in zip(CELLS, theories)
    p = Params(c.p; MoveFitnessSDGen = -1, Sample = 10)
    tz = mean(@elapsed(simulate(p; seed = SEED_OFFSET + 90_000_000 + s)) for s in 1:2)
    line = @sprintf("  %-6s %8.2f s   ", c.name, tz)
    for B in BS
        tb = mean(@elapsed(simulate(p; seed = SEED_OFFSET + 91_000_000 + s, init = :equilibrium, burnin = B, at = T0,
                                    target = th.target, spectrum = th.st)) for s in 1:2)
        line *= @sprintf("%6.2f s (×%4.1f)   ", tb, tz / tb)
    end
    say(line)
end
say("  generations run: from zero $(T0 + CAP); ", join((@sprintf("B=%d: %d (×%.1f)", B, B + CAP, (T0 + CAP) / (B + CAP)) for B in BS), ", "))

# --- summary ------------------------------------------------------------------------------

say("\n", "="^100, "\nSUMMARY")
for c in CELLS
    b, v = answers[c.name]
    say(@sprintf("  %-6s %s", c.name, b === nothing ? "no burn-in of the grid passes" :
                 "smallest burn-in from which none fails: B = $b ($v); generations saved in a run capped $CAP after the move: ×$(round((T0 + CAP) / (b + CAP), digits = 1))"))
end
pubverdict === nothing || say("  PUB    ", join(("B = $B: $(pubverdict[B])" for B in PUBBS), ", "))
say(isempty(allfails) ? "no failures" : "failures (" * string(length(allfails)) * "):")
foreach(f -> say("  FAIL  ", f), allfails)
close(io)
SMOKE || println("log: ", LOGPATH)
exit(any(c -> answers[c.name][1] === nothing, CELLS) ? 1 : 0)
