# The beta closure of the :spectrum model (`spectrum = (; solver = :beta)`, src/spectrum.jl)
# against the frequency grid (the default), under the criteria of validate/spectrum_beta_criteria.md
# (committed before this was run). A check of numerical equivalence between two solvers of one
# model, not a comparison with the kernel. Writes validate/spectrum_beta.log.
#
# Run from the package root (vignettes environment: Catlab):
#   julia -t 8 --project=vignettes validate/spectrum_beta.jl          # everything; resumes from the caches
# Caches (resumable `sweep` CSVs, one row per input): validate/spectrum_beta_runs.csv (metrics 1–3,
# 6, 7), spectrum_beta_ld.csv (5), spectrum_beta_k50.csv (4).

# 2026-11-19: `ne = :falconer`, `solver = :iteration` and `QTLSim.burnin_start_vg` below keep the effective size, the solver and
# the starting variance this script was run with (the package's defaults are now `:breeders` and `:bisection`; NEWS.md), so that
# it still reproduces its log.
using QTLSim, Printf, Random, Statistics
import Catlab

const HERE = @__DIR__
const SIMDIR = normpath(joinpath(HERE, "..", "..", "Repository", "Simulations"))
const STRATA = ["QT0100-Possion3.0", "QT0100-Possion5.0", "QT0500-Possion3.0", "QT0500-Possion5.0",
                "QT1000-Possion3.0", "QT1000-Possion5.0", "QT5000-Possion3.0", "QT5000-Possion5.0"]
const SEED, NPER, NLD, NK50 = 20261130, 50, 25, 6
const BURNIN, HORIZON = 2000, 2000
const BETA = (; solver = :beta)

const LOG = Ref{IO}(stdout)
say(s...) = (println(stdout, s...); LOG[] === stdout || (println(LOG[], s...); flush(LOG[])))
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))
pf(b) = b ? "PASS" : "FAIL"

# the sample: per stratum, a random permutation of its files in sorted order; the first NPER,
# of which the first NLD (ld) and the first NK50 (half-survival rate)
function sample_inputs()
    rng = MersenneTwister(SEED)
    main, ld, k50 = String[], String[], String[]
    for d in STRATA
        fs = sort(filter(endswith(".txt"), readdir(joinpath(SIMDIR, d))))
        pick = ["$d/$f" for f in fs[randperm(rng, length(fs))[1:NPER]]]
        append!(main, pick); append!(ld, pick[1:NLD]); append!(k50, pick[1:NK50])
    end
    main, ld, k50
end

readp(id) = Params(read_qtsim_options(joinpath(SIMDIR, id)); Seed = nothing)

# one solver on one input: equilibrium V_G, the trajectory as emulate runs it, the emulator
function solve(p, spec)
    q = Params(p; Sample = 1, NoGenerations = BURNIN + HORIZON, MoveFitnessSDGen = BURNIN)
    VG = QTLSim.burnin_start_vg(p)
    eq = spectrum_vg(p; spec...)
    m = run_moments(q; VG, erosion = true, loci = :spectrum, spectrum = spec)
    v(t) = BURNIN + t + 1 <= length(m.geno_var) ? Float64(m.geno_var[BURNIN+t+1]) : NaN
    t = @elapsed e = emulate(p; loci = :spectrum, burnin = BURNIN, horizon = HORIZON, VG, spectrum = spec)
    (; eq, v0 = v(0), v100 = v(100), v500 = v(500), P = e.P2000, T = e.T === missing ? NaN : Float64(e.T), t,
     fb = haskey(e, :fallback) ? Float64(e.fallback) : NaN)
end

function pair_job(id, extra)
    p = readp(id)
    g = solve(p, extra)
    b = solve(p, merge(extra, BETA))
    (; eq_g = g.eq, eq_b = b.eq, v0_g = g.v0, v0_b = b.v0, v100_g = g.v100, v100_b = b.v100, v500_g = g.v500,
     v500_b = b.v500, P_g = g.P, P_b = b.P, T_g = g.T, T_b = b.T, t_g = g.t, t_b = b.t, fb = b.fb)
end

function k50_job(id)
    p = readp(id)
    tg = @elapsed g = emulated_half_survival_rate(p; loci = :spectrum)
    tb = @elapsed b = emulated_half_survival_rate(p; loci = :spectrum, spectrum = BETA)
    (; k_g = g.k50, k_b = b.k50, evals_g = Float64(g.evals), evals_b = Float64(b.evals), fallbacks = Float64(b.fallbacks), t_g = tg, t_b = tb)
end

q95(x) = isempty(x) ? NaN : quantile(x, 0.95)
function vgmetric(rows, g, b, label)
    r = [getfield(x, b) / getfield(x, g) for x in rows if isfinite(getfield(x, g)) && isfinite(getfield(x, b)) && getfield(x, g) > 0]
    e = abs.(r .- 1)
    ok = !isempty(e) && median(e) <= 0.02 && q95(e) <= 0.10
    sayf("  %-22s n = %3d  median ratio %.4f  median |r−1| %.4f  95th pct %.4f  max %.4f   %s", label, length(r),
         median(r), median(e), q95(e), maximum(e), pf(ok))
    ok
end

function emumetric(rows, label)
    agree = mean([(x.P_g >= 0.5) == (x.P_b >= 0.5) for x in rows])
    both = [x for x in rows if isfinite(x.T_g) && isfinite(x.T_b)]
    one = count(x -> isfinite(x.T_g) != isfinite(x.T_b), rows)
    rT = [x.T_b / x.T_g for x in both]
    okP = agree >= 0.98
    okT = !isempty(rT) && 0.98 <= median(rT) <= 1.02
    sayf("  %-22s P2000 call agreement %.4f (%d/%d)   %s", label, agree, round(Int, agree * length(rows)), length(rows), pf(okP))
    sayf("  %-22s T ratio: n = %d, median %.4f, IQR %.4f–%.4f, 95th pct of |r−1| %.4f; finite in one only: %d   %s", label,
         length(rT), median(rT), quantile(rT, 0.25), quantile(rT, 0.75), q95(abs.(rT .- 1)), one, pf(okT))
    sayf("  %-22s mean |ΔP2000| %.4f, max %.4f", label, mean(abs(x.P_b - x.P_g) for x in rows), maximum(abs(x.P_b - x.P_g) for x in rows))
    okP, okT
end

function main()
    open(joinpath(HERE, "spectrum_beta.log"), "w") do io
        LOG[] = io
        say("validate/spectrum_beta.jl — beta closure vs grid (criteria: validate/spectrum_beta_criteria.md)")
        say("threads: ", Threads.nthreads(), "   ", Libc.strftime("%Y-%m-%d %H:%M", time()))
        main_ids, ld_ids, k50_ids = sample_inputs()
        say("sample: ", length(main_ids), " inputs (", NPER, " per stratum, seed ", SEED, "); ld ", length(ld_ids), "; k50 ", length(k50_ids))
        # warm up (compilation out of the timings)
        w = readp(main_ids[1]); solve(Params(w; MoveFitnessSD = 10.0), (;)); solve(Params(w; MoveFitnessSD = 10.0), BETA)
        t0 = time()
        rows = sweep(id -> pair_job(id, (;)), main_ids; cache = joinpath(HERE, "spectrum_beta_runs.csv"), progress = 25)
        sayf("main runs done (%.1f min)", (time() - t0) / 60)
        t0 = time()
        lrows = sweep(id -> pair_job(id, (; ld = :dynamic)), ld_ids; cache = joinpath(HERE, "spectrum_beta_ld.csv"), progress = 25)
        sayf("ld runs done (%.1f min)", (time() - t0) / 60)
        t0 = time()
        krows = sweep(k50_job, k50_ids; cache = joinpath(HERE, "spectrum_beta_k50.csv"), progress = 8)
        sayf("half-survival runs done (%.1f min)", (time() - t0) / 60)

        verdicts = Pair{String,Bool}[]
        say("\n1. Equilibrium V_G")
        push!(verdicts, "1 equilibrium V_G" => vgmetric(rows, :eq_g, :eq_b, "equilibrium"))
        say("\n2. V_G along the move (both alive)")
        sayf("  (at the move, +0, no verdict:)"); vgmetric(rows, :v0_g, :v0_b, "+0")
        push!(verdicts, "2 V_G +100" => vgmetric(rows, :v100_g, :v100_b, "+100"))
        push!(verdicts, "2 V_G +500" => vgmetric(rows, :v500_g, :v500_b, "+500"))
        say("\n3. The emulator")
        okP, okT = emumetric(rows, "emulate")
        push!(verdicts, "3 P2000 calls" => okP, "3 T ratio" => okT)
        say("\n4. Half-survival rate (", length(krows), " inputs)")
        fin = [x for x in krows if isfinite(x.k_g) && isfinite(x.k_b) && x.k_g > 0 && x.k_b > 0]
        special = [x for x in krows if !(x in fin)]
        same = count(x -> x.k_g == x.k_b, special)
        rk = [x.k_b / x.k_g for x in fin]
        ok4 = !isempty(rk) && 0.98 <= median(rk) <= 1.02 && q95(abs.(rk .- 1)) <= 0.10 && same == length(special)
        sayf("  k50 ratio: n = %d, median %.4f, IQR %.4f–%.4f, 95th pct |r−1| %.4f, max %.4f; 0/Inf: %d (agreeing %d)   %s",
             length(rk), median(rk), quantile(rk, 0.25), quantile(rk, 0.75), q95(abs.(rk .- 1)), maximum(abs.(rk .- 1)),
             length(special), same, pf(ok4))
        sayf("  evaluations with a fallback: %d of %d (%.1f%%); time grid %.0f s, beta %.0f s (median per input %.1f×)",
             sum(x.fallbacks for x in krows), sum(x.evals_b for x in krows), 100sum(x.fallbacks for x in krows) / sum(x.evals_b for x in krows),
             sum(x.t_g for x in krows), sum(x.t_b for x in krows), median([x.t_g / x.t_b for x in krows]))
        push!(verdicts, "4 half-survival rate" => ok4)
        say("\n5. With ld = :dynamic (", length(lrows), " inputs)")
        push!(verdicts, "5 ld: equilibrium V_G" => vgmetric(lrows, :eq_g, :eq_b, "equilibrium"))
        push!(verdicts, "5 ld: V_G +100" => vgmetric(lrows, :v100_g, :v100_b, "+100"))
        push!(verdicts, "5 ld: V_G +500" => vgmetric(lrows, :v500_g, :v500_b, "+500"))
        lP, lT = emumetric(lrows, "emulate (ld)")
        push!(verdicts, "5 ld: P2000 calls" => lP, "5 ld: T ratio" => lT)
        say("\n6. Fallback rate")
        fb = [x.fb for x in rows]
        rate = mean(fb .>= 0)
        sayf("  emulate fell back in %d of %d inputs (%.3f)   %s", count(>=(0), fb), length(fb), rate, pf(rate <= 0.5))
        sayf("  before or at the move (generation ≤ %d): %d; generation of fallback after the move, quartiles: %s", BURNIN + 1,
             count(x -> 0 <= x <= BURNIN + 1, fb), let a = [x - BURNIN for x in fb if x > BURNIN + 1]
                 isempty(a) ? "-" : join(round.(Int, quantile(a, [0.25, 0.5, 0.75])), " / ")
             end)
        fbl = [x.fb for x in lrows]
        sayf("  (ld = :dynamic: %d of %d, %.3f)", count(>=(0), fbl), length(fbl), mean(fbl .>= 0))
        push!(verdicts, "6 fallback rate" => rate <= 0.5)
        say("\n7. Speed-up (emulate, grid time / beta time, same job)")
        sp = [x.t_g / x.t_b for x in rows]
        spf = [x.t_g / x.t_b for x in rows if x.fb >= 0]
        sayf("  median %.1f×  quartiles %.1f / %.1f  min %.2f×  max %.0f×; among fallbacks median %.1f× (n = %d)   %s", median(sp),
             quantile(sp, 0.25), quantile(sp, 0.75), minimum(sp), maximum(sp), isempty(spf) ? NaN : median(spf), length(spf), pf(median(sp) >= 10))
        sayf("  total emulate time: grid %.0f s, beta %.0f s; ld: median %.1f×", sum(x.t_g for x in rows), sum(x.t_b for x in rows),
             median([x.t_g / x.t_b for x in lrows]))
        push!(verdicts, "7 speed-up" => median(sp) >= 10)
        say("\nVerdicts")
        for (k, v) in verdicts
            sayf("  %-26s %s", k, pf(v))
        end
        sayf("  %d of %d PASS", count(last, verdicts), length(verdicts))
        LOG[] = stdout
    end
end

main()
