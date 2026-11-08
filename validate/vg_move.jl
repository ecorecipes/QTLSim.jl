# The genetic variance under a moving optimum: the allele-frequency-spectrum model of the
# moments lowering (`run_moments(p; loci = :spectrum)`, src/spectrum.jl) against the
# kernel, and the rare-alleles closure (`loci = :rare_alleles, erosion = true`, the
# default with drift) for comparison: V_G(t), lag(t) and the time to extinction after
# the optimum starts moving, at rates from slow to fast, in the house-of-cards, Gaussian
# and drift regimes, at a fixed size and with Beverton–Holt density regulation, and on a
# sample of the paper's own inputs against the published table and the kernel.
#
# Usage (from the package root; Catlab and XLSX are in the vignettes environment):
#     julia -t 3 --project=vignettes validate/vg_move.jl
# Writes validate/vg_move.log (validate/vg_move_seed<SEED_OFFSET>.log with ENV
# SEED_OFFSET ≠ 0). SMOKE=1 runs short versions of every part and writes no log. Part P
# needs Repository/Simulations and Repository/Results (read only); it is skipped without
# them. The script exits 1 if any criterion fails.
#
# DESIGN. Each kernel replicate runs to the move once and is copied, with its random
# number generator, for every rate (as in vignette 10; the same runs as `simulate` with
# that seed). V_G (`geno_var`) and the lag θ − z̄ are averaged over windows of generations
# relative to the move: pre [−200, 0), +100 [80, 120), +250 [230, 270), +500 [480, 520),
# +1000 [980, 1020), +2000 [1960, 2000) (kernel: every 10th generation, over the
# replicates alive at the window; predictions: every generation). Extinction: the kernel's
# median generation of extinction after the move (survivors count as later than any
# death; "survive" when at least half are alive 2000 generations after the move) and the
# fraction alive at +2000; predictions: the median of `extinction_distribution(p, m; from
# = MOVE)` along the model's trajectory m (demographic noise along the deterministic
# trajectory) and its probability of survival to +2000.
# - A, regimes, base V_P = 1, h² = 0.5, m = 0.01, u = 2 (:qtsim crossover), N = 1000,
#   L = 100, move at 10,000 (from zero genomes), 20 replicates:
#   H, house of cards: μ = 1e-4 (4Nμ = 0.4), V_s = 50.5 (moments_regimes.jl's cell S);
#   G, Gaussian: μ = 0.0025 (4Nμ = 10), V_s = 50.5;
#   H2: H with Beverton–Holt R = 2.
#   Rates: 0.01 (cell S's) and f·k_die, f ∈ {0.075, 0.15, 0.3, 0.6, 1.2}, with k_die =
#   λ_die V/S the rate whose Gaussian steady lag is the lag λ_die at which N W̄ = 2, V the
#   spectrum model's V_G at the move and S = ω² + V_E + V.
# - B, drift regime with demography: vignette 10's population (N = 1000, L = 100, μ =
#   1e-4, V_P = 0.022, h² = 0.32, m = 0.004, u = 3, ω = 5, move at 5000), demography
#   fixed, R = 5, 2, 1.5; vignette 10's rates 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1;
#   40 replicates.
# - P, the paper's inputs: 30 published inputs drawn at random from each of the 8 strata
#   (NoQTL ∈ {100, 500, 1000, 5000} × u ∈ {3, 5}), matched to their table rows
#   ("Generations to extinction", 240000 = survived; dying generation − 9999). Times are
#   capped at 2001 generations after the move (survival beyond that counts as 2001).
#   Predictions from zero genomes (spectrum, and rare with erosion). The kernel: 3
#   replicates on 4 inputs per stratum with NoQTL ≤ 1000 (24 inputs), run to +2000.
#
# PRE-SPECIFIED CRITERIA (written before any kernel run of this script). A comparison of
# a mean FAILS only if the kernel's whole 95% CI lies outside the band:
# |kernel − pred| − 1.96 SE > tol·pred.
# - V1: pre-move V_G within 15% of the spectrum model's, in H, H2 and every demography of B.
# - V2: V_G at +500 and +1000 within 25% of the spectrum model's, in every cell of H, H2
#   and B where at least half the kernel's replicates and the model's population are
#   alive at that window.
# - V3: over the V2 comparisons, the median |log(kernel/spectrum)| is smaller than the
#   median |log(kernel/rare)|.
# - V4 (the rise, H at rate 0.01): the kernel's V_G at +2000 exceeds its pre-move V_G by
#   more than 10%, and the spectrum model's ratio (+2000 over pre) is within 15% of the
#   kernel's.
# - L1: the lag at +500 and +1000 within 20% of the spectrum model's, in the V2 cells.
# - L2: over the L1 comparisons, median |log(kernel/spectrum)| < median |log(kernel/rare)|.
# - E1: in every cell of H, H2 and B where at least half the kernel's replicates died
#   within 2000 generations, the spectrum model's predicted median time to extinction is
#   within 20% of the kernel's median.
# - E2: where at least 90% of the kernel's replicates survive to +2000, the spectrum
#   model's probability of survival is at least 0.5; where at most 10% survive, it is at
#   most 0.5 (H, H2, B).
# - E3: over the cells of E1, median |log(T_kernel/T_spectrum)| < median |log(T_kernel/T_rare)|.
# - P1: over the paper's inputs, the median of log(T_table/T_spectrum) (capped times)
#   within ±log 1.15. P2: the median |log(T_table/T_pred)| smaller for the spectrum than
#   for rare. P3: over the 24 kernel inputs, the median of log(T_kernel/T_spectrum)
#   (T_kernel the median of the 3 replicates, capped) within ±log 1.2.
# - N1 (numerics): in one cell of each of H, G and B, refining the spectrum's grid and
#   effect classes (480 frequencies, 40 classes) changes its V_G at every window by less
#   than 2%.
# Pre-specified expectations, reported but not criteria: in G the spectrum model's V_G is
# about 40–50% above the kernel's before the move (infinite sites: a locus that carries
# many mutations is the Gaussian regime, where the rare-alleles closure is right), so its
# extinction times there are too late; in mixed cells (10–90% survive) the deterministic
# trajectory with demographic noise cannot give the right probability of survival.
#
# DISCLOSURE. The spectrum model was built after seeing: vignette 10's cached runs
# (rescue_runs.csv: times to extinction at the slow rates), moments_regimes.log's cell S
# (V_G 1.8 → 2.5 under the move), and exploratory kernel runs of vignette 10's population
# at a fixed size and rate 0.02 (12 replicates, seeds 501–512, V_G and lag at +100 …
# +1250) with QTSim's crossover and with nearly free recombination (CrossoverModel = :all,
# RecombPoissonMean = 50). Those runs showed that the collapse of V_G is faster with
# QTSim's single crossover than with free recombination, and that the linkage-equilibrium
# spectrum follows the free-recombination runs. A genetic-draft term (hitchhiking on
# completed sweeps, Kaplan, Hudson and Langley 1989) was tried first and rejected because
# its same-locus part made the free-recombination prediction collapse, which those runs
# did not; Barton's (1995) loose-linkage N_e term replaced it. With drift alone dissolving the
# association (δ = 1/(2n)) it predicted V_G at +500 of 0.053 against the kernel's 0.069
# with QTSim's crossover (0.100 against 0.093 with free recombination), and vignette 10's
# cached times 10–15% early at a fixed size and R = 5. A smoke run of this script
# (SMOKE=1: 3 replicates, move at generation 600, not equilibrated) then showed the
# kernel's V_G rising under fast moves in H while that version collapsed it: b² V_G is
# then large and Barton's loose-linkage approximation (weak selection per locus against
# recombination) fails. δ was then given the third factor of Santiago and Caballero
# (1998), the reduction of the selected variance by selection, V_G/S per generation
# (Bulmer 1971): δ = 1/(2n) + V_G/S. With it, vignette 10's cached times at a fixed size
# and R = 5 are within 5% and at R = 2 and 1.5 up to 20% late (the first version was
# within 5% there). No constant was fitted. The full run uses seeds not used by the smoke
# run (offset 32,000,000; the smoke run used 31,000,000).
#
# HISTORY. Run 1 (2026-11-08, SEED_OFFSET = 0, vg_move.log, 35 min on 3 threads): 148
# PASS, 17 FAIL, kept as failures; the model was not changed after it. Passed: V1 (pre-move
# V_G 0.95–1.10 of the spectrum in H, H2 and B), V3 (median |log| error of V_G under the
# move 0.161 against the rare-alleles closure's 0.286, 43 comparisons), V4 (the rise in H
# at 0.01: kernel ×1.49, spectrum ×1.58; rare ×1.00), L2 (0.107 against 0.109), every E2
# (survival classified right in all 44 cells with ≤ 10% or ≥ 90% survivors), E3 (median
# |log(T_kernel/T_pred)| 0.038 against 0.106, 24 cells), N1, P1–P3 (paper's 240 inputs:
# median T_table/T_spectrum 1.000, |log| 0.036 against rare's 0.048; kernel/spectrum on 24
# inputs 1.015). E1 passed in 22 of 24 cells, including every cell of B (vignette 10:
# within 0–6% at a fixed size and R = 5, 0.5–20% late at R = 2 and 1.5, e.g. R = 1.5 at
# 0.01: 1508 against the kernel's 1261, where the rare-alleles closure predicts survival).
# Failed: in H and H2 at intermediate rates (0.12, 0.24) the spectrum over-predicts the rise
# of V_G (5.7–6.7 against the kernel's 3.3–4.3; V2 ×7) and so under-predicts the lag (L1
# ×7: kernel lag 1.34–1.75 times the spectrum's); E1 in H at 0.47 (kernel 232, spectrum
# 405, rare 220) and H2 at 0.47 (48 against 61); V2 in D, fixed, 0.02 at +1000 (kernel V_G
# 0.030, spectrum 0.019: the spectrum collapses later but deeper; its extinction time, 1424
# against 1509, passed). Expectation for G confirmed and worse than expected: the kernel's
# V_G is 0.68 of the spectrum's before the move, and under moves of 0.21 and 0.42 the
# kernel's V_G collapses (4.0 → 1.2 within 500 generations) and it dies (942 and 209
# generations) where both models predict survival.
# After run 1, without kernel runs, two other decay rates of the linked-selection
# association were tried on the model alone: δ = 1/(2n) (the first version) captures G's
# collapse (dies at +740; kernel 942) but collapses H at 0.24, where the kernel survives;
# δ = 1/(2n) + α|λ|/S (per-allele selection against recombination) is worse in G, H and
# vignette 10. Neither was adopted; no run 2 was made. The open problem is the rate at which
# selection dissolves a site's association with its linked background, which decides
# whether linked selection collapses V_G (G, D) or lets it rise (H).


using QTLSim, Statistics, Printf, Random
import Catlab, XLSX

const SMOKE = get(ENV, "SMOKE", "") == "1"
const SEED_OFFSET = parse(Int, get(ENV, "SEED_OFFSET", "0"))
const SEED0 = 32_000_000 + SEED_OFFSET
logpath = SMOKE ? tempname() : joinpath(@__DIR__, SEED_OFFSET == 0 ? "vg_move.log" : "vg_move_seed$(SEED_OFFSET).log")
io = open(logpath, "w")
say(s...) = (println(stdout, s...); println(io, s...); flush(io))
se(x) = length(x) > 1 ? std(x) / sqrt(length(x)) : NaN
fails = String[]
check(ok, what) = (say(ok ? "  PASS  " : "  FAIL  ", what); ok || push!(fails, what))
outside(k, s, pred, tol) = abs(k - pred) - 1.96 * (isnan(s) ? 0.0 : s) > tol * abs(pred)
lr(a, b) = abs(log(a / b))

const AFTER = SMOKE ? 300 : 2000
const CAP = AFTER + 1
const WINDOWS = SMOKE ? [("pre", -100, 0), ("+100", 80, 120), ("+250", 230, 270)] :
                [("pre", -200, 0), ("+100", 80, 120), ("+250", 230, 270), ("+500", 480, 520), ("+1000", 980, 1020), ("+2000", 1960, 2000)]
const CMPW = SMOKE ? ("+100", "+250") : ("+500", "+1000")
const LASTW = WINDOWS[end][1]
wmean(gen, x, move, lo, hi) = (ix = findall(g -> move + lo <= g < move + hi, gen); isempty(ix) ? NaN : mean(x[ix]))

# --- the kernel: one pre-move run per replicate, copied for every rate ----------------
function kernel_cell(p::Params, move::Int, rates, seed::Int)
    w, rng = World(Params(p; Sample = 10)), MersenneTwister(seed)
    init_env!(w, rng)
    rec = QTLSim.records()
    for it in 0:move-1
        w.gen = it
        generation!(w, rng; rec)
        w.extinct && error("a replicate died before the move (seed $seed)")
    end
    map(rates) do d
        w2, rng2 = deepcopy(w), deepcopy(rng)
        w2.p = Params(w.p; MoveFitnessSD = d, MoveFitnessSDGen = move, MoveType = :continuous)
        r2 = QTLSim.records()
        T = -1
        for it in move:move+AFTER-1
            w2.gen = it
            generation!(w2, rng2; rec = r2)
            w2.extinct && (T = it - move; break)
        end
        gen = vcat(rec.gen, r2.gen); vg = vcat(rec.geno_var, r2.geno_var); lag = vcat(rec.theta .- rec.pheno_mean, r2.theta .- r2.pheno_mean)
        alive(hi) = T < 0 || T >= hi
        (T = T, vg = [alive(hi) ? wmean(gen, vg, move, lo, hi) : NaN for (_, lo, hi) in WINDOWS],
         lag = [alive(hi) ? wmean(gen, lag, move, lo, hi) : NaN for (_, lo, hi) in WINDOWS])
    end
end

function kernel_runs(p, move, rates, R, seed0)
    out = Vector{Any}(undef, R)
    Threads.@threads :greedy for r in 1:R
        out[r] = kernel_cell(p, move, rates, SEED0 + seed0 + r)
    end
    [[out[r][i] for r in 1:R] for i in eachindex(rates)]      # by rate, then replicate
end

# --- predictions ----------------------------------------------------------------------
function predict(p::Params, move::Int, d; loci, spectrum = (;))
    q = Params(p; MoveFitnessSD = d, MoveFitnessSDGen = move, MoveType = :continuous, NoGenerations = move + AFTER, Sample = 1)
    m = loci === :spectrum ? run_moments(q; loci, spectrum) : run_moments(q; loci, erosion = true)
    e = extinction_distribution(q, m; from = move)
    lag = m.theta .- m.pheno_mean
    (T = ismissing(e.median) ? CAP : min(e.median - move, CAP), psurv = e.psurvive, med = e.median,
     det = m.extinct ? m.extinct_gen - move : -1,
     vg = [wmean(m.gen, m.geno_var, move, lo, hi) for (_, lo, hi) in WINDOWS],
     lag = [wmean(m.gen, lag, move, lo, hi) for (_, lo, hi) in WINDOWS])
end

# λ_die: N c ω_e/√S exp(−λ²/2S) = 2, and k_die = λ_die V/S
function kdie(p::Params, V)
    c = 1 / (0.4 * sqrt(2π))
    S = initFitSD(p)^2 + envVar(p) + V
    λ = sqrt(2S * log(p.PopSize * c * initFitSD(p) / sqrt(S) / 2))
    λ * V / S
end

# one cell's kernel summary
function ksummary(rs)
    nw = length(WINDOWS)
    Ts = sort([r.T < 0 ? typemax(Int) : r.T for r in rs])
    med = Ts[cld(length(Ts), 2)]
    vgs = [filter(!isnan, [r.vg[k] for r in rs]) for k in 1:nw]
    lags = [filter(!isnan, [r.lag[k] for r in rs]) for k in 1:nw]
    (alive = count(r -> r.T < 0, rs) / length(rs), med = med == typemax(Int) ? missing : med,
     nalive = [length(v) for v in vgs], vg = mn.(vgs), vgse = se.(vgs),
     lag = [isempty(l) ? NaN : mean(l) for l in lags], lagse = se.(lags))
end
mn(v) = isempty(v) ? NaN : mean(v)

rows = []          # (regime, demog, rate, kernel summary, spectrum, rare)
say("QTLSim validate/vg_move.jl: V_G, lag and extinction under a moving optimum (spectrum model vs kernel)")
say("Julia threads: ", Threads.nthreads(), "; SEED_OFFSET = $SEED_OFFSET; ", SMOKE ? "SMOKE" : "full")
say()

base = Params(Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01, RecombPoissonMean = 2.0,
              NoQTL = 100, PopSize = 1000)
const MOVEA = SMOKE ? 600 : 10_000
RA = SMOKE ? 3 : 20
regimesA = [("H", Params(base; MutationRatePerGamete = 1e-4)), ("G", Params(base; MutationRatePerGamete = 0.0025)),
            ("H2", Params(base; MutationRatePerGamete = 1e-4, Demography = :beverton_holt, GrowthRate = 2.0))]
for (ri, (name, p)) in enumerate(regimesA)
    q0 = Params(p; NoGenerations = MOVEA, Sample = 1)
    V = run_moments(q0; loci = :spectrum).geno_var[end]
    rates = vcat(0.01, [f * kdie(p, V) for f in (0.075, 0.15, 0.3, 0.6, 1.2)])
    t = @elapsed ks = kernel_runs(Params(p; NoGenerations = MOVEA + AFTER), MOVEA, rates, RA, 100_000 * ri)
    say(@sprintf("A/%s: 4Nμ = %g, V_s = %.1f, spectrum V_G at the move %.4f, k_die %.4f; %d replicates (%.0f s)",
                 name, 4 * p.PopSize * p.MutationRatePerGamete, selection_variance(p), V, kdie(p, V), RA, t))
    for (i, d) in enumerate(rates)
        push!(rows, (; regime = name, demog = p.Demography === :fixed ? "fixed" : "R = $(p.GrowthRate)", rate = d, k = ksummary(ks[i]),
                     s = predict(p, MOVEA, d; loci = :spectrum), r = predict(p, MOVEA, d; loci = :rare_alleles), p, move = MOVEA))
    end
end

pB = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
            MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0)
const MOVEB = SMOKE ? 600 : 5000
RB = SMOKE ? 3 : 40
ratesB = [0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1]
demogB = [("fixed", (Demography = :fixed,)), ("R = 5", (Demography = :beverton_holt, GrowthRate = 5.0)),
          ("R = 2", (Demography = :beverton_holt, GrowthRate = 2.0)), ("R = 1.5", (Demography = :beverton_holt, GrowthRate = 1.5))]
for (di, (dname, kw)) in enumerate(demogB)
    p = Params(pB; kw...)
    t = @elapsed ks = kernel_runs(Params(p; NoGenerations = MOVEB + AFTER), MOVEB, ratesB, RB, 1_000_000 + 100_000 * di)
    say(@sprintf("B/%s: vignette 10's population, %d replicates (%.0f s)", dname, RB, t))
    for (i, d) in enumerate(ratesB)
        push!(rows, (; regime = "D", demog = dname, rate = d, k = ksummary(ks[i]),
                     s = predict(p, MOVEB, d; loci = :spectrum), r = predict(p, MOVEB, d; loci = :rare_alleles), p, move = MOVEB))
    end
end
say()

# --- tables -----------------------------------------------------------------------------
fmtT(x) = ismissing(x) ? "surv" : x >= CAP ? "surv" : string(round(Int, x))
say("V_G by window: kernel (alive replicates) / spectrum / rare; columns ", join(first.(WINDOWS), ", "))
for r in rows
    say(@sprintf("%-3s %-8s rate %.4f  kernel %s", r.regime, r.demog, r.rate, join([@sprintf("%.4g(%d)", r.k.vg[j], r.k.nalive[j]) for j in eachindex(WINDOWS)], " ")))
    say(@sprintf("%28s spectrum %s", "", join([@sprintf("%.4g", x) for x in r.s.vg], " ")))
    say(@sprintf("%28s rare     %s", "", join([@sprintf("%.4g", x) for x in r.r.vg], " ")))
end
say()
say("Lag θ − z̄ by window: kernel / spectrum / rare")
for r in rows
    say(@sprintf("%-3s %-8s rate %.4f  kernel %s | spectrum %s | rare %s", r.regime, r.demog, r.rate,
                 join([@sprintf("%.3g", x) for x in r.k.lag], " "), join([@sprintf("%.3g", x) for x in r.s.lag], " "),
                 join([@sprintf("%.3g", x) for x in r.r.lag], " ")))
end
say()
say("Extinction: kernel alive at +$AFTER and median time after the move; predictions: median (P survive), deterministic death")
for r in rows
    say(@sprintf("%-3s %-8s rate %.4f  kernel %3.0f%% alive, median %5s | spectrum %5s (%.2f) det %5d | rare %5s (%.2f) det %5d",
                 r.regime, r.demog, r.rate, 100r.k.alive, fmtT(r.k.med), fmtT(r.s.T), r.s.psurv, r.s.det, fmtT(r.r.T), r.r.psurv, r.r.det))
end
say()

# --- criteria, A and B --------------------------------------------------------------------
say("Criteria, regimes and demography")
crit = [r for r in rows if r.regime != "G"]
cellname(r) = @sprintf("%s %s rate %.4f", r.regime, r.demog, r.rate)
for r in unique(r -> (r.regime, r.demog), crit)
    check(!outside(r.k.vg[1], r.k.vgse[1], r.s.vg[1], 0.15),
          @sprintf("V1: %s %s pre-move V_G: kernel %.4g within 15%% of spectrum %.4g (ratio %.3f)", r.regime, r.demog, r.k.vg[1], r.s.vg[1], r.k.vg[1] / r.s.vg[1]))
end
v2 = []; l1 = []
for r in crit, wname in CMPW
    j = findfirst(w -> w[1] == wname, WINDOWS)
    (r.k.nalive[j] >= 1 && r.k.nalive[j] >= 0.5 * (r.regime == "D" ? RB : RA) && !isnan(r.s.vg[j])) || continue
    push!(v2, (r, j))
    check(!outside(r.k.vg[j], r.k.vgse[j], r.s.vg[j], 0.25),
          @sprintf("V2: %s at %s: kernel V_G %.4g within 25%% of spectrum %.4g (ratio %.3f; rare %.3f)", cellname(r), wname,
                   r.k.vg[j], r.s.vg[j], r.k.vg[j] / r.s.vg[j], r.k.vg[j] / r.r.vg[j]))
end
if !isempty(v2)
    ms = median([lr(r.k.vg[j], r.s.vg[j]) for (r, j) in v2]); mr = median([lr(r.k.vg[j], r.r.vg[j]) for (r, j) in v2 if !isnan(r.r.vg[j])])
    check(ms < mr, @sprintf("V3: median |log(kernel/pred)| of V_G under the move: spectrum %.3f < rare %.3f (%d comparisons)", ms, mr, length(v2)))
end
rH = only([r for r in rows if r.regime == "H" && r.rate == 0.01])
kr, sr = rH.k.vg[end] / rH.k.vg[1], rH.s.vg[end] / rH.s.vg[1]
check(kr > 1.1 && abs(sr / kr - 1) <= 0.15,
      @sprintf("V4: the rise in H at 0.01: kernel V_G %s/pre %.3f (> 1.1), spectrum %.3f (within 15%%: ratio %.3f); rare %.3f", LASTW, kr, sr, sr / kr, rH.r.vg[end] / rH.r.vg[1]))
for (r, j) in v2
    check(!outside(r.k.lag[j], r.k.lagse[j], r.s.lag[j], 0.20),
          @sprintf("L1: %s at %s: kernel lag %.4g within 20%% of spectrum %.4g (ratio %.3f; rare %.3f)", cellname(r), WINDOWS[j][1],
                   r.k.lag[j], r.s.lag[j], r.k.lag[j] / r.s.lag[j], r.k.lag[j] / r.r.lag[j]))
end
if !isempty(v2)
    ok = [(r, j) for (r, j) in v2 if r.k.lag[j] > 0 && r.s.lag[j] > 0 && r.r.lag[j] > 0]
    ms = median([lr(r.k.lag[j], r.s.lag[j]) for (r, j) in ok]); mr = median([lr(r.k.lag[j], r.r.lag[j]) for (r, j) in ok])
    check(ms < mr, @sprintf("L2: median |log(kernel/pred)| of the lag: spectrum %.3f < rare %.3f", ms, mr))
end
e1 = [r for r in crit if !ismissing(r.k.med)]
for r in e1
    check(r.s.T < CAP && abs(r.k.med - r.s.T) <= 0.20 * r.s.T,
          @sprintf("E1: %s: kernel median %d within 20%% of spectrum %s (rare %s)", cellname(r), r.k.med, fmtT(r.s.T), fmtT(r.r.T)))
end
for r in crit
    r.k.alive >= 0.9 && check(r.s.psurv >= 0.5, @sprintf("E2: %s: kernel %.0f%% alive, spectrum P(survive) %.2f ≥ 0.5", cellname(r), 100r.k.alive, r.s.psurv))
    r.k.alive <= 0.1 && check(r.s.psurv <= 0.5, @sprintf("E2: %s: kernel %.0f%% alive, spectrum P(survive) %.2f ≤ 0.5", cellname(r), 100r.k.alive, r.s.psurv))
end
if !isempty(e1)
    ms = median([lr(r.k.med, r.s.T) for r in e1]); mr = median([lr(r.k.med, r.r.T) for r in e1])
    check(ms < mr, @sprintf("E3: median |log(T_kernel/T_pred)|: spectrum %.3f < rare %.3f (%d cells)", ms, mr, length(e1)))
end
say("  expectations (not criteria), G: pre-move kernel/spectrum V_G " *
    @sprintf("%.3f", [r for r in rows if r.regime == "G"][1].k.vg[1] / [r for r in rows if r.regime == "G"][1].s.vg[1]) *
    "; G extinction medians kernel/spectrum/rare: " *
    join([@sprintf("%s/%s/%s", fmtT(r.k.med), fmtT(r.s.T), fmtT(r.r.T)) for r in rows if r.regime == "G"], ", "))
mixed = [r for r in crit if 0.1 < r.k.alive < 0.9]
say("  expectations (not criteria), mixed cells (10–90% alive): " *
    join([@sprintf("%s: kernel %.0f%%, spectrum P %.2f", cellname(r), 100r.k.alive, r.s.psurv) for r in mixed], "; "))

# N1: numerics
for r in [rows[findfirst(r -> r.regime == "H", rows) + 2], rows[findfirst(r -> r.regime == "G", rows) + 2], rows[findfirst(r -> r.regime == "D", rows) + 3]]
    f = predict(r.p, r.move, r.rate; loci = :spectrum, spectrum = (grid = 480, classes = 40))
    dev = maximum(abs(f.vg[j] / r.s.vg[j] - 1) for j in eachindex(WINDOWS) if !isnan(f.vg[j]) && !isnan(r.s.vg[j]); init = 0.0)
    check(dev < 0.02, @sprintf("N1: %s: grid 480, 40 classes changes V_G by at most %.2g%% (< 2%%)", cellname(r), 100dev))
end
say()

# --- P: the paper's inputs ------------------------------------------------------------------
const ROOTDIR = normpath(joinpath(@__DIR__, "..", ".."))
const SIMDIR = joinpath(ROOTDIR, "Repository", "Simulations")
const MASTER = joinpath(ROOTDIR, "Repository", "Results", "Master Output file.xlsx")
if isdir(SIMDIR) && isfile(MASTER)
    dt = XLSX.readtable(MASTER, "Master output file")
    lab = string.(dt.column_labels)
    col(name) = dt.data[findfirst(==(name), lab)]
    num(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)
    tab = (L = round.(Int, num.(col("No. Loci"))), u = round.(Int, num.(col("recombination u=3,5"))),
           H = num.(col("Heritability")), PV = num.(col("PhenotypeVar")), MS = num.(col("MutationScalar Vm/Ve")),
           d = num.(col("d_env")), T = round.(Int, num.(col("Generations to extinction"))))
    function match_row(L, u, q)
        hits = findall(i -> tab.L[i] == L && tab.u[i] == u && abs(tab.H[i] - q.Heritability) < 1e-9 &&
                            abs(tab.PV[i] - q.PhenotypeVar) < 1e-9 && abs(tab.d[i] - q.MoveFitnessSD) <= 5e-3 * q.MoveFitnessSD &&
                            abs(tab.MS[i] - q.MutationScalar) <= 5e-3 * q.MutationScalar, eachindex(tab.L))
        length(hits) == 1 ? hits[1] : 0
    end
    strata = [("QT0100-Possion3.0", 100, 3), ("QT0100-Possion5.0", 100, 5), ("QT0500-Possion3.0", 500, 3),
              ("QT0500-Possion5.0", 500, 5), ("QT1000-Possion3.0", 1000, 3), ("QT1000-Possion5.0", 1000, 5),
              ("QT5000-Possion3.0", 5000, 3), ("QT5000-Possion5.0", 5000, 5)]
    NPER = SMOKE ? 1 : 30
    rng = MersenneTwister(20261108 + SEED_OFFSET)
    jobs = []
    for (dir, L, u) in strata
        files = sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, dir))))
        for f in sort(shuffle(rng, files)[1:NPER])
            q = read_qtsim_options(joinpath(SIMDIR, dir, f))
            r = match_row(L, u, q)
            r == 0 && (say("  no unique table row for $dir/$f; skipped"); continue)
            push!(jobs, (; dir, file = f, L, u, q = Params(q; Seed = nothing), row = r))
        end
    end
    t0 = time()
    preds = Vector{Any}(undef, length(jobs))
    Threads.@threads :greedy for k in eachindex(jobs)
        j = jobs[k]
        move = j.q.MoveFitnessSDGen
        q = SMOKE ? Params(j.q; MoveFitnessSDGen = 300) : j.q
        mv = q.MoveFitnessSDGen
        preds[k] = (s = predict(q, mv, q.MoveFitnessSD; loci = :spectrum), r = predict(q, mv, q.MoveFitnessSD; loci = :rare_alleles))
    end
    capT(T) = T >= 240_000 ? CAP : min(T, CAP)
    prow = [(; jobs[k]..., obs = capT(tab.T[jobs[k].row]), s = preds[k].s.T, r = preds[k].r.T) for k in eachindex(jobs)]
    say(@sprintf("P, the paper's inputs: %d matched to the table (%.0f s); times capped at %d generations after the move", length(prow), time() - t0, CAP))
    say("  stratum               n  table died  spec died  rare died  median log(table/spec)  log(table/rare)  |log| spec  |log| rare")
    for (dir, L, u) in vcat(strata, [("all", 0, 0)])
        rs = dir == "all" ? prow : [r for r in prow if r.dir == dir]
        isempty(rs) && continue
        say(@sprintf("  %-18s %3d  %10d %10d %10d %22.3f %16.3f %11.3f %10.3f", dir, length(rs), count(r -> r.obs < CAP, rs),
                     count(r -> r.s < CAP, rs), count(r -> r.r < CAP, rs), median([log(r.obs / r.s) for r in rs]),
                     median([log(r.obs / r.r) for r in rs]), median([lr(r.obs, r.s) for r in rs]), median([lr(r.obs, r.r) for r in rs])))
    end
    m1 = median([log(r.obs / r.s) for r in prow])
    check(abs(m1) <= log(1.15), @sprintf("P1: median log(T_table/T_spectrum) %.3f within ±log 1.15 (ratio %.3f)", m1, exp(m1)))
    ms, mr = median([lr(r.obs, r.s) for r in prow]), median([lr(r.obs, r.r) for r in prow])
    check(ms < mr, @sprintf("P2: median |log(T_table/T_pred)|: spectrum %.3f < rare %.3f", ms, mr))
    kin = vcat([[r for r in prow if r.dir == dir][1:min(SMOKE ? 1 : 4, count(r -> r.dir == dir, prow))] for (dir, L, u) in strata if L <= 1000]...)
    NK = SMOKE ? 1 : 3
    kj = [(i, rep) for i in eachindex(kin) for rep in 1:NK]
    kT = Vector{Int}(undef, length(kj))
    t = @elapsed Threads.@threads :greedy for k in eachindex(kj)
        i, rep = kj[k]
        q = kin[i].q
        q = SMOKE ? Params(q; MoveFitnessSDGen = 300) : q
        mv = q.MoveFitnessSDGen
        o = simulate(Params(q; NoGenerations = mv + AFTER, Sample = 1000); seed = SEED0 + 5_000_000 + 1000i + rep)
        kT[k] = o.extinct ? min(o.extinct_gen - mv + 1, CAP) : CAP
    end
    kmed = [median([kT[k] for k in eachindex(kj) if kj[k][1] == i]) for i in eachindex(kin)]
    say(@sprintf("  P3: kernel on %d inputs × %d replicates (%.0f s)", length(kin), NK, t))
    for (r, km) in zip(kin, kmed)
        say(@sprintf("    %s/%s  kernel %s  table %s  spectrum %s  rare %s", r.dir, r.file, fmtT(km), fmtT(r.obs), fmtT(r.s), fmtT(r.r)))
    end
    m3 = median([log(km / r.s) for (r, km) in zip(kin, kmed)])
    say(@sprintf("  P3 (rare, not a criterion): median log(T_kernel/T_rare) %.3f", median([log(km / r.r) for (r, km) in zip(kin, kmed)])))
    check(abs(m3) <= log(1.2), @sprintf("P3: median log(T_kernel/T_spectrum) %.3f within ±log 1.2 (ratio %.3f)", m3, exp(m3)))
else
    say("P skipped: Repository/Simulations or the published table not found")
end
say()

say(isempty(fails) ? "PASS: every pre-specified criterion holds" : "FAILED ($(length(fails))): " * join(fails, "; "))
close(io)
exit(isempty(fails) ? 0 : 1)
