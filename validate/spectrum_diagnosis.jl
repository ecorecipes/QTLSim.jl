# Diagnostic study (exploratory): why does the emulator with loci = :spectrum put the
# half-survival rate too high near the threshold (validate/emulator_spectrum.log, E6′)?
# A causal decomposition of the emulator's chain
#   burn-in V_G → V_G during the move → lag θ − z̄ → W̄ (Gaussian) → extinction given W̄
# against kernel trajectories of the same inputs and seeds as E6′. Running log with the
# hypotheses stated before each test: validate/spectrum_diagnosis.log; findings:
# validate/spectrum_diagnosis.md.
#
# Run from the package root (vignettes environment: XLSX, Catlab):
#   julia -t 8 --project=vignettes validate/spectrum_diagnosis.jl inputs    # E6′ inputs in selection order, diagnostic sets
#   julia -t 8 --project=vignettes validate/spectrum_diagnosis.jl kernel SET [R]   # kernel trajectories (resumable)
#   julia -t 8 --project=vignettes validate/spectrum_diagnosis.jl analyse SET      # the decomposition
#   julia -t 8 --project=vignettes validate/spectrum_diagnosis.jl k50              # k50: kernel grid vs emulator
#   julia -t 8 --project=vignettes validate/spectrum_diagnosis.jl confirm          # prototype P1 on E6′ inputs not used for discovery
#   julia -t 8 --project=vignettes validate/spectrum_diagnosis.jl fresh            # prototype P1 on 64 fresh inputs
# SET: D, C, Dctl, Cctl joined by +, or all. ENV VARIANT=many (or free) runs the kernel and the
# spectrum with looser linkage. The driven-spectrum comparisons call driven_spectrum /
# driven_variant from small scratch drivers (their output is in the log's appendix).
# Raw per-run trajectories go to validate/spectrum_diagnosis_raw/ (not committed; one
# Serialization file per run, so a run is resumable); the per-input summaries that the
# analysis prints are in the log.

# 2026-11-19: `ne = :falconer`, `solver = :iteration` and `QTLSim.burnin_start_vg` below keep the effective size, the solver and
# the starting variance this script was run with (the package's defaults are now `:breeders` and `:bisection`; NEWS.md), so that
# it still reproduces its log.
using QTLSim, Statistics, Random, Printf, Serialization
import Catlab

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "..", ".."))
const SIMDIR = joinpath(ROOT, "Repository", "Simulations")
const RAW = get(ENV, "RAWDIR", joinpath(HERE, "spectrum_diagnosis_raw"))
const INPUTS = joinpath(HERE, "spectrum_diagnosis_inputs.csv")
const MOVE, KBURN, KSEED, RK = 10_000, 1000, 30261116000, 40
const BURNIN, HORIZON = 2000, 2000
const EVERY = 10                       # genic variance etc. every EVERY generations

kernel_params(p) = Params(p; NoGenerations = 12000, Sample = 1, Seed = nothing)
moment_params(p) = Params(p; Sample = 1, NoGenerations = BURNIN + HORIZON, MoveFitnessSDGen = BURNIN)
# VARIANT (ENV): "" the input as published; "free": free recombination (CrossoverModel = :all,
# RecombPoissonMean = (NoQTL − 1)/2, so that any two loci split with probability ≈ ½)
const VARIANT = get(ENV, "VARIANT", "")
function readp(id)
    p = Params(read_qtsim_options(joinpath(SIMDIR, id)); Seed = nothing)
    VARIANT == "" && return p
    if VARIANT == "free"
        q = Params(p; CrossoverModel = :all, RecombPoissonMean = (p.NoQTL - 1) / 2)
        # the exact pair classes of so many crossovers are slow to compute (pair_splits); every
        # pair of distinct loci then splits with probability ½ to within O(1/L), so enter that
        # in pair_classes' memo
        key = (q.NoQTL, Float64(q.RecombPoissonMean), :all)
        lock(QTLSim._PAIR_CLASSES_LOCK) do
            QTLSim._PAIR_CLASSES[key] = ([1 / q.NoQTL, 1 - 1 / q.NoQTL], [0.0, 0.5])
        end
        return q
    end
    # "many": every one of Poisson(20) crossover points applied (QTSim applies one), much looser
    # linkage at an affordable cost (free recombination costs ~400 s per 1000-locus run)
    VARIANT == "many" && return Params(p; CrossoverModel = :all, RecombPoissonMean = 20.0)
    error("unknown VARIANT $VARIANT")
end

# ---------------------------------------------------------------------------------------
# inputs: the E6′ inputs in selection order j (so that seeds KSEED + 100j + r reproduce the
# E6′ kernel runs), with the kernel's fraction alive and the emulator's P2000

# validate/emulator_spectrum.jl's functions (selection, statistics), without running its main
function e6p_module()
    src = read(joinpath(HERE, "emulator_spectrum.jl"), String)
    src = replace(src, r"\nmain\(\)\s*$" => "\n")
    m = Module(:E6P)
    Core.eval(m, :(include(f) = Base.include($m, f)))
    Base.include_string(m, "using QTLSim, Statistics, Random, Printf, LinearAlgebra; import Catlab\n" *
                           replace(src, "const MODE = isempty(ARGS) ? \"full\" : ARGS[1]" => "const MODE = \"select\""),
                        joinpath(HERE, "emulator_spectrum.jl"))
    m
end

function e6p_order()
    m = e6p_module()
    Base.invokelatest() do
        ids, ps, strat = m.read_inputs()
        ex = m.excluded_ids()
        S = m.main_sample(ids, strat, ex)
        E, = m.e6p_inputs(ids, ps, strat, union(ex, Set(ids[S])))
        ids[E]
    end
end

function load_e6p()
    kr = Dict{String,Vector{Tuple{Int,Bool,Int}}}()
    for l in readlines(joinpath(HERE, "emulator_spectrum_kernel.csv"))[2:end]
        f = split(l, ","); id, r = split(f[1], "#")
        push!(get!(kr, String(id), Tuple{Int,Bool,Int}[]), (parse(Int, r), f[2] == "true", parse(Int, f[3])))
    end
    em = Dict{String,Float64}()
    for l in readlines(joinpath(HERE, "emulator_spectrum_runs.csv"))[2:end]
        f = split(l, ","); startswith(f[1], "e6p-spectrum|") || continue
        em[String(split(f[1], "|")[2])] = parse(Float64, f[2])
    end
    kr, em
end

function make_inputs()
    E = e6p_order()
    length(E) == 160 || error("expected 160 E6′ inputs, got $(length(E))")
    kr, em = load_e6p()
    rows = [(j = j, id = id, f = count(x -> x[2], kr[id]) / RK, P = em[id]) for (j, id) in enumerate(E)]
    disc = [r for r in rows if (r.P > 0.5) != (r.f > 0.5)]
    conc = [r for r in rows if (r.P > 0.5) == (r.f > 0.5)]
    # discovery / confirmation: discordant inputs alternately in id order; controls: one
    # concordant input per stratum drawn at random (MersenneTwister(20261114)), split the same way
    sort!(disc; by = r -> r.id)
    set = Dict{String,String}()
    for (k, r) in enumerate(disc)
        set[r.id] = isodd(k) ? "D" : "C"
    end
    rng = MersenneTwister(20261114)
    strata = sort(unique(first(split(r.id, "/")) for r in rows))
    for (s, st) in enumerate(strata)
        cs = sort([r for r in conc if startswith(r.id, st)]; by = r -> r.id)
        pick = cs[randperm(rng, length(cs))[1:2]]
        set[pick[1].id] = "Dctl"; set[pick[2].id] = "Cctl"
    end
    open(INPUTS, "w") do io
        println(io, "j,id,set,kernel_f,P2000")
        for r in rows
            println(io, r.j, ",", r.id, ",", get(set, r.id, ""), ",", r.f, ",", r.P)
        end
    end
    println("wrote $INPUTS: $(length(disc)) discordant (D $(count(==("D"), values(set))), C $(count(==("C"), values(set)))), " *
            "controls Dctl $(count(==("Dctl"), values(set))), Cctl $(count(==("Cctl"), values(set)))")
end

function read_inputs_csv()
    out = NamedTuple[]
    for l in readlines(INPUTS)[2:end]
        f = split(l, ",")
        push!(out, (j = parse(Int, f[1]), id = String(f[2]), set = String(f[3]), f = parse(Float64, f[4]), P = parse(Float64, f[5])))
    end
    out
end
inset(rows, set) = set == "all" ? filter(r -> !isempty(r.set), rows) :
                   filter(r -> r.set in split(set, "+"), rows)

# ---------------------------------------------------------------------------------------
# the kernel, with the extra statistics; the same random draws as `simulate`, so the run of
# seed KSEED + 100j + r is E6′'s run r of input j

function genic_stats(w; seg = true)
    G = genomes(w)                        # (L, 2, N) of the current zygotes
    L, P, N = size(G)
    genic = 0.0; nseg = 0
    buf = Vector{Float64}(undef, P * N)
    for l in 1:L
        k = 0
        for i in 1:N, c in 1:P
            k += 1; buf[k] = G[l, c, i]
        end
        genic += P * var(buf)
        seg || continue
        sort!(buf)
        for k in 2:length(buf)
            buf[k] != buf[k-1] && (nseg += 1)
        end
    end
    g = vec(sum(G; dims = (1, 2)))
    m = mean(g); v = var(g)
    sk = v > 0 ? mean((g .- m) .^ 3) / v^1.5 : NaN
    ku = v > 0 ? mean((g .- m) .^ 4) / v^2 - 3 : NaN
    (genic, seg ? float(nseg) : NaN, sk, ku)
end

function kernel_traj(p, seed; from = MOVE - 1000, every = EVERY)
    q = kernel_params(p)
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    n = q.NoGenerations - w.gen
    gen = Int[]; VG = Float64[]; VP = Float64[]; zbar = Float64[]; gbar = Float64[]; θ = Float64[]
    Wbar = Float64[]; nsurv = Int[]; sg = Float64[]
    xg = Int[]; genic = Float64[]; nseg = Float64[]; skew = Float64[]; kurt = Float64[]
    extinct_gen = -1
    for it in w.gen:q.NoGenerations-1
        w.gen = it
        if it >= MOVE - 300 && it % every == 0
            s = genic_stats(w; seg = it % 50 == 0)
            push!(xg, it); push!(genic, s[1]); push!(nseg, s[2]); push!(skew, s[3]); push!(kurt, s[4])
        end
        # generation!, step by step, to read the survivors' mean genotype
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        N = w.N
        if it >= from
            push!(gen, it); push!(VG, var(view(w.genotype, 1:N))); push!(VP, var(view(w.phenotype, 1:N)))
            push!(zbar, mean(view(w.phenotype, 1:N))); push!(gbar, mean(view(w.genotype, 1:N))); push!(θ, w.θ)
            push!(Wbar, mean(view(w.fitness, 1:N))); push!(nsurv, w.nsurv)
            push!(sg, w.nsurv > 0 ? mean(w.genotype[i] for i in w.survivors) - gbar[end] : NaN)
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    (; seed, extinct_gen, gen, VG, VP, zbar, gbar, θ, Wbar, nsurv, sg, xg, genic, nseg, skew, kurt)
end

rawfile(id, r) = joinpath(RAW, replace(id, "/" => "__", ".txt" => "") * (VARIANT == "" ? "" : "@" * VARIANT) * "#$r.jls")

function run_kernel(set, R)
    mkpath(RAW)
    rows = inset(read_inputs_csv(), set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(rawfile(r.id, k))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) kernel runs to do for set $set (R = $R), $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0)
    t0 = time()
    Threads.@threads :greedy for (id, j, k) in jobs
        tr = kernel_traj(readp(id), KSEED + 100j + k)
        serialize(rawfile(id, k), tr)
        d = Threads.atomic_add!(done, 1) + 1
        d % 20 == 0 && @printf("  %d/%d  %.0f s\n", d, length(jobs), time() - t0)
    end
    @printf("kernel: %d runs, %.0f s wall\n", length(jobs), time() - t0)
end

# ---------------------------------------------------------------------------------------
# the emulator's trajectory (as `emulate` runs it), and the Gaussian lag recursion for the
# substitutions

function spectrum_traj(p; spectrum = (;), VG = QTLSim.burnin_start_vg(p), burnin = BURNIN)
    q = Params(p; Sample = 1, NoGenerations = burnin + HORIZON, MoveFitnessSDGen = burnin)
    m = run_moments(q; VG, erosion = true, loci = :spectrum, spectrum)
    t = m.gen .- burnin
    (; t, VG = m.geno_var, VP = m.pheno_var, lag = m.theta .- m.pheno_mean, Wbar = m.mean_fitness,
     nsurv = m.nsurv, extinct = m.extinct, extinct_gen = m.extinct ? m.extinct_gen - burnin : -1, q, m)
end

const C0 = 1 / (0.4 * sqrt(2π))

"""
Gaussian lag recursion of the moments lowering (moments_select!) at a fixed size, given a
V_G for each generation t = 0, 1, … of the move (`VGt[t+1]`, held at its last value past
its end), the lag `λ0` before the move, V_E and ω: θ moves by k, λ = θ − z̄, z̄ ← z̄ + V_G λ/S,
S = ω² + V_G + V_E, W̄ = c ω/√S exp(−λ²/2S). Optionally `resp` scales the response.
Returns (; lag, Wbar) for t = 0:T−1.
"""
function lag_recursion(p, VGt; T = HORIZON, λ0 = 0.0, resp = 1.0)
    VE = QTLSim.envVar(p); ω = QTLSim.initFitSD(p); k = p.MoveFitnessSD
    lag = zeros(T); W = zeros(T)
    λ = λ0
    for t in 1:T
        V = VGt[min(t, length(VGt))]
        isfinite(V) || (V = 0.0)
        λ += k
        S = ω^2 + V + VE
        lag[t] = λ
        W[t] = C0 * ω / sqrt(S) * exp(-λ^2 / (2S))
        λ -= resp * V * λ / S
    end
    (; lag, Wbar = W)
end

"P2000 and median T (counted from 1) of extinction_distribution along a W̄ for t = 0:T−1 of the move."
function survival(p, W)
    q = kernel_params(p)
    d = extinction_distribution(q, collect(Float64, W))
    (P = d.psurvive, T = d.median === missing ? missing : d.median + 1)
end

# ---------------------------------------------------------------------------------------
# analysis of one input

function load_runs(id, R)
    runs = []
    for k in 1:R
        f = rawfile(id, k)
        isfile(f) && push!(runs, deserialize(f))
    end
    runs
end

# the mean over the runs alive at each generation t = gen − MOVE in `ts` of the field `f`
function mean_alive(runs, f, ts; x = :gen)
    out = fill(NaN, length(ts)); nal = zeros(Int, length(ts))
    for (i, t) in enumerate(ts)
        s = 0.0; n = 0
        for r in runs
            g = getfield(r, x)
            k = searchsortedfirst(g, MOVE + t)
            (k <= length(g) && g[k] == MOVE + t) || continue
            v = getfield(r, f)[k]
            isfinite(v) || continue
            s += v; n += 1
        end
        out[i] = n > 0 ? s / n : NaN; nal[i] = n
    end
    out, nal
end

wmean(v) = (u = filter(isfinite, v); isempty(u) ? NaN : mean(u))
fmt(x) = x === missing ? "  surv" : (x isa Real && isnan(x) ? "   NaN" : @sprintf("%6.0f", x))

function analyse_input(r, R; io = stdout, verbose = true)
    p = readp(r.id)
    runs = load_runs(r.id, R)
    isempty(runs) && return nothing
    st = spectrum_traj(p)
    death = [x.extinct_gen < 0 ? missing : x.extinct_gen - MOVE + 1 for x in runs]
    nd = count(!ismissing, death)
    Tk = nd > length(runs) / 2 ? median(sort(collect(skipmissing(death)))[1:end]) : missing
    # the kernel's own fraction alive at +2000 among these R runs
    fk = 1 - nd / length(runs)
    ts = -200:HORIZON-1
    kVG, nal = mean_alive(runs, :VG, ts)
    kVP, _ = mean_alive(runs, :VP, ts)
    klag = mean_alive(runs, :θ, ts)[1] .- mean_alive(runs, :zbar, ts)[1]
    kW, _ = mean_alive(runs, :Wbar, ts)
    tx = collect(-200:EVERY:HORIZON-1)
    kgen, _ = mean_alive(runs, :genic, tx; x = :xg)
    kseg, _ = mean_alive(runs, :nseg, tx; x = :xg)
    kkurt, _ = mean_alive(runs, :kurt, tx; x = :xg)
    kskew, _ = mean_alive(runs, :skew, tx; x = :xg)
    sidx(t) = findfirst(==(t), st.t)
    sval(v, t) = (k = sidx(t); k === nothing ? NaN : v[k])
    kidx(t) = t - first(ts) + 1
    xidx(t) = findfirst(==(t), tx)
    # substitutions: the Gaussian recursion with (a) the spectrum's own V_G (a check of the
    # recursion), (b) the kernel's mean V_G over the runs alive, (c) per run, its own V_G
    nT = HORIZON
    sVG = [sval(st.VG, t) for t in 0:nT-1]
    lastfin(v) = (k = findlast(isfinite, v); k === nothing ? NaN : v[k])
    sVG = [isfinite(v) ? v : lastfin(sVG) for v in sVG]
    λ0s = sval(st.lag, -1); λ0k = klag[kidx(-1)]
    a = lag_recursion(p, sVG; λ0 = λ0s)
    sa = survival(p, a.Wbar)
    kv = kVG[kidx(0):end]
    kv = [isfinite(v) ? v : lastfin(kv) for v in kv]
    b = lag_recursion(p, kv; λ0 = λ0k)
    sb = survival(p, b.Wbar)
    # (b′) kernel's pre-move V_G held constant through the move (burn-in link alone)
    vpre = wmean(kVG[kidx(-200):kidx(-1)])
    b2 = lag_recursion(p, fill(vpre, nT); λ0 = λ0k)
    sb2 = survival(p, b2.Wbar)
    # (a′) spectrum's pre-move V_G held constant
    spre = wmean([sval(st.VG, t) for t in -200:-1])
    a2 = lag_recursion(p, fill(spre, nT); λ0 = λ0s)
    sa2 = survival(p, a2.Wbar)
    # (d) the response link: kernel's realised survivors' genotype shift vs V_G λ/S, pooled
    # over generations 0..min(death)−1 of all runs (relative bias of the mean), and the W̄ link
    num = 0.0; den = 0.0; wnum = 0.0; wden = 0.0
    ω2 = QTLSim.initFitSD(p)^2; VE = QTLSim.envVar(p)
    lastW = Float64[]
    for x in runs
        for k in eachindex(x.gen)
            t = x.gen[k] - MOVE
            0 <= t < HORIZON || continue
            λ = x.θ[k] - x.zbar[k]; S = ω2 + x.VP[k]
            isfinite(x.sg[k]) || continue
            num += x.sg[k]; den += x.VG[k] * λ / S
            wnum += x.Wbar[k]; wden += C0 * sqrt(ω2 / S) * exp(-λ^2 / (2S))
        end
        x.extinct_gen >= 0 && push!(lastW, x.Wbar[end])
    end
    res = (; id = r.id, set = r.set, L = p.NoQTL, u = p.RecombPoissonMean, rate = p.MoveFitnessSD / sqrt(p.PhenotypeVar),
           f40 = r.f, fk, Tk, P_em = r.P, T_sp = st.extinct ? st.extinct_gen + 1 : missing,
           VGpre_k = vpre, VGpre_s = spre, genpre_k = wmean(kgen[1:xidx(-10)]),
           VG = [(t, kVG[kidx(t)], sval(st.VG, t), nal[kidx(t)]) for t in (100, 250, 500, 1000, 1999)],
           genic = [(t, kgen[xidx(t)], kseg[xidx(t)], kkurt[xidx(t)], kskew[xidx(t)]) for t in (-50, 100, 250, 500, 1000)],
           lag = [(t, klag[kidx(t)], sval(st.lag, t)) for t in (100, 250, 500, 1000, 1999)],
           P_a = sa.P, T_a = sa.T, P_b = sb.P, T_b = sb.T, P_b2 = sb2.P, T_b2 = sb2.T, P_a2 = sa2.P, T_a2 = sa2.T,
           resp_ratio = num / den, W_ratio = wnum / wden, lastW = isempty(lastW) ? NaN : median(lastW),
           seg_s = nothing)
    if verbose
        @printf(io, "\n  %s [%s] L=%d u=%g rate=%.4f  E6′ f=%.3f, P2000=%.3f; these %d runs: alive %.2f, kernel median T %s; spectrum T %s\n",
                r.id, r.set, p.NoQTL, p.RecombPoissonMean, res.rate, r.f, r.P, length(runs), fk, fmt(Tk), fmt(res.T_sp))
        @printf(io, "    pre-move V_G: kernel %.4g (genic %.4g), spectrum %.4g, ratio k/s %.3f\n", vpre, res.genpre_k, spre, vpre / spre)
        @printf(io, "    V_G (kernel alive-mean / spectrum / n alive):")
        for (t, a1, b1, n) in res.VG
            @printf(io, "  +%d %.4g/%.4g/%d", t, a1, b1, n)
        end
        @printf(io, "\n    genic/LD (kernel): ")
        for (t, g, s, ku, sk) in res.genic
            vk = kVG[kidx(t)]
            @printf(io, " %+d genic %.4g LD %.3g nseg %.0f kurt %.2f skew %.2f;", t, g, vk - g, s, ku, sk)
        end
        @printf(io, "\n    lag (kernel / spectrum):")
        for (t, a1, b1) in res.lag
            @printf(io, "  +%d %.3g/%.3g", t, a1, b1)
        end
        @printf(io, "\n    substitutions (P2000, T): spectrum's V_G in the Gaussian recursion %.3f %s; kernel mean V_G %.3f %s; kernel pre-move V_G constant %.3f %s; spectrum pre-move V_G constant %.3f %s\n",
                res.P_a, fmt(res.T_a), res.P_b, fmt(res.T_b), res.P_b2, fmt(res.T_b2), res.P_a2, fmt(res.T_a2))
        @printf(io, "    response realised/Gaussian %.3f; W̄ realised/Gaussian %.3f; median W̄ in a run's last generation %.4f\n",
                res.resp_ratio, res.W_ratio, res.lastW)
    end
    res
end

"""
The spectrum's genic variance driven by one kernel run's realised drivers: from generation
MOVE − 300 (the spectrum's equilibrium shape scaled to the run's genic variance then), each
generation selection with the run's lag λ and S = ω² + V_P, drift among its realised
survivors, and 2N L μ mutations; returns the spectrum's V_G at the run's genic sampling
generations, so that its genic dynamics are compared with the kernel's for the same lag
and parents (no feedback). `linkage` is the spectrum's (:barton or :none).
"""
function driven_spectrum(p, x; linkage = :barton, t0 = MOVE - 300)
    k0 = findfirst(==(t0), x.xg)
    sp = QTLSim.spectrum_scale!(deepcopy(QTLSim.spectrum_equilibrium(p; linkage).spectrum), x.genic[k0])
    ω2 = QTLSim.initFitSD(p)^2
    out = fill(NaN, length(x.xg))
    out[k0] = QTLSim.spectrum_vg(sp)
    for k in eachindex(x.gen)
        g = x.gen[k]
        g < t0 && continue
        λ = x.θ[k] - x.zbar[k]; S = ω2 + x.VP[k]
        QTLSim.spectrum_select(sp, λ, S)
        QTLSim.spectrum_step!(sp, max(x.nsurv[k], 1))
        QTLSim.spectrum_mutate!(sp, p.PopSize)
        j = findfirst(==(g + 1), x.xg)
        j === nothing || (out[j] = QTLSim.spectrum_vg(sp))
    end
    out
end

"""
Barton's linked-selection drift `c` (Var(Δp) = pq(1/(2n) + c)) as `spectrum_linked`, with the
association decay δ = 1/(2n) + dS·V/S (dS = 1 is the closure's, 0 drops selection's part) and
gradient variance b² V with V = `V` (the closure uses its genic variance).
"""
function linked_c(sp, n, V; dS = 1.0, Vδ = V)
    (isinf(sp.S) || sp.λ == 0) && return 0.0
    b2V = (sp.λ / sp.S)^2 * V
    δ = 1 / (2n) + dS * Vδ / sp.S
    acc = 0.0
    for k in 2:length(sp.ω)
        acc += sp.ω[k] / (2 * (sp.ρ[k] + δ)^2)
    end
    b2V * acc / (2n)
end

"""
driven_spectrum with variants of the linked-selection drift, applied through the number of
parents (n_eff = 1/(1/n + 2c), with the spectrum's own term off): `mode` = :barton (as the
closure), :dS0 (δ = 1/(2n)), :total (b²V and δ with the run's total V_G instead of the genic).
"""
function driven_variant(p, x; mode = :barton, t0 = MOVE - 300)
    k0 = findfirst(==(t0), x.xg)
    sp = QTLSim.spectrum_scale!(deepcopy(QTLSim.spectrum_equilibrium(p).spectrum), x.genic[k0])
    sp.linkage = :none
    ω2 = QTLSim.initFitSD(p)^2
    out = fill(NaN, length(x.xg))
    out[k0] = QTLSim.spectrum_vg(sp)
    for k in eachindex(x.gen)
        g = x.gen[k]
        g < t0 && continue
        λ = x.θ[k] - x.zbar[k]; S = ω2 + x.VP[k]
        QTLSim.spectrum_select(sp, λ, S)
        n = max(x.nsurv[k], 1)
        Vg = QTLSim.spectrum_vg(sp)
        c = mode isa Real ? linked_c(sp, n, Vg; dS = mode) :
            mode === :barton ? linked_c(sp, n, Vg) : mode === :dS0 ? linked_c(sp, n, Vg; dS = 0.0) :
            mode === :total ? linked_c(sp, n, x.VG[k]) : error("mode")
        QTLSim.spectrum_step!(sp, 1 / (1 / n + 2c))
        QTLSim.spectrum_mutate!(sp, p.PopSize)
        j = findfirst(==(g + 1), x.xg)
        j === nothing || (out[j] = QTLSim.spectrum_vg(sp))
    end
    out
end

function analyse(set, R)
    rows = inset(read_inputs_csv(), set)
    out = []
    for r in rows
        x = analyse_input(r, R)
        x === nothing || push!(out, x)
    end
    out
end

# ---------------------------------------------------------------------------------------
# k50 decomposition on a few inputs: the kernel's k50 from runs at a grid of rates (multiples
# of the input's MoveFitnessSD), the emulator's, and the Gaussian recursion's with V_G held
# constant at the spectrum's or the kernel's pre-move value

const K50_IDS = ["QT0500-Possion3.0/Input-04943.txt", "QT1000-Possion3.0/Input-09697.txt", "QT1000-Possion5.0/Input-14110.txt",
                 "QT5000-Possion3.0/Input-12000.txt", "QT5000-Possion5.0/Input-12463.txt", "QT5000-Possion5.0/Input-01925.txt"]
const K50_MULT = [0.25, 0.35, 0.5, 0.7, 1.0]
const K50_R, K50_SEED = 10, 50261114000

function kernel_k50_runs(ids)
    jobs = [(; key = "$(id)×$(m)#$r", id, m, seed = K50_SEED + 1000i + 100j + r)
            for (i, id) in enumerate(ids) for (j, m) in enumerate(K50_MULT) for r in 1:K50_R]
    sweep(jobs; id = j -> j.key, cost = j -> readp(j.id).NoQTL, progress = 50,
          cache = joinpath(HERE, "spectrum_diagnosis_k50_kernel" * (VARIANT == "" ? "" : "_" * VARIANT) * ".csv")) do j
        p = readp(j.id)
        o = simulate(kernel_params(Params(p; MoveFitnessSD = j.m * p.MoveFitnessSD)); seed = j.seed, init = :equilibrium, burnin = KBURN)
        (; alive = !o.extinct, extinct_gen = o.extinct_gen)
    end
    jobs
end

# x50 by linear interpolation in log rate between the grid points around fraction alive ½
function interp50(x, f)
    for i in 1:length(x)-1
        if f[i] >= 0.5 > f[i+1]
            return exp(log(x[i]) + (f[i] - 0.5) / (f[i] - f[i+1]) * (log(x[i+1]) - log(x[i])))
        end
    end
    f[end] >= 0.5 ? Inf : (f[1] < 0.5 ? 0.0 : NaN)
end

# the rate at which the Gaussian recursion with V_G(t) = VGfun(k) keeps half alive at +2000
function bisect_k50(p, VGfun; lo = 1e-4, hi = 10.0)
    P(k) = survival(p, lag_recursion(Params(p; MoveFitnessSD = k), VGfun(k)).Wbar).P
    P(hi) >= 0.5 && return Inf
    P(lo) < 0.5 && return 0.0
    for _ in 1:50
        m = sqrt(lo * hi)
        P(m) >= 0.5 ? (lo = m) : (hi = m)
    end
    sqrt(lo * hi)
end

# ---------------------------------------------------------------------------------------
# prototype (exploratory): the spectrum with the LD factor r = V_G/genic of the multilocus
# haplotype moments at the fixed optimum with drift (multilocus_equilibrium, a model
# quantity, nothing fitted to the kernel), `spectrum = (; ldratio = r)`

ld_ratio(p) = (e = multilocus_equilibrium(p; drift = true, loci = :multilocus); clamp(e.VG / e.genic, 1e-3, 1.0))

"P2000 of the emulator pipeline (emulate's) with the spectrum keywords `kw`."
function proto_P2000(p; kw = (;))
    st = spectrum_traj(p; spectrum = kw)
    q = st.q
    m = st.m
    length(m.gen) <= BURNIN && return 0.0
    extinction_distribution(q, m; from = BURNIN).psurvive
end

function proto_k50(p; kw = (;), lo = nothing, hi = nothing)
    kc = drift_critical_rate(p; ne = :falconer, solver = :iteration).k
    lo = something(lo, kc / 64); hi = something(hi, 8kc)
    P(k) = proto_P2000(Params(p; MoveFitnessSD = k); kw)
    P(hi) >= 0.5 && return Inf
    P(lo) < 0.5 && return 0.0
    for _ in 1:14                       # log k to about 2^-14 of the bracket
        m = sqrt(lo * hi)
        P(m) >= 0.5 ? (lo = m) : (hi = m)
    end
    sqrt(lo * hi)
end

function k50_decomposition(R = 10; protos = Dict{String,Any}())
    jobs = kernel_k50_runs(K50_IDS)
    res = Dict{Tuple{String,Float64},Vector{Bool}}()
    cache = Dict(String(split(l, ",")[1]) => split(l, ",")[2] == "true"
                 for l in readlines(joinpath(HERE, "spectrum_diagnosis_k50_kernel" * (VARIANT == "" ? "" : "_" * VARIANT) * ".csv"))[2:end])
    for j in jobs
        push!(get!(res, (j.id, j.m), Bool[]), cache[j.key])
    end
    out = []
    for id in K50_IDS
        p = readp(id); k0 = p.MoveFitnessSD
        f = [count(res[(id, m)]) / length(res[(id, m)]) for m in K50_MULT]
        kK = interp50(K50_MULT .* k0, f)
        kE = emulated_half_survival_rate(p; loci = :spectrum).k50
        runs = load_runs(id, R)
        tsel(x) = [x.VG[k] for k in eachindex(x.gen) if MOVE - 200 <= x.gen[k] < MOVE]
        vk = mean(mean(tsel(x)) for x in runs)
        st = spectrum_traj(p)
        vs = mean(st.VG[k] for k in eachindex(st.t) if -200 <= st.t[k] < 0)
        kA2 = bisect_k50(p, k -> [vs]); kB2 = bisect_k50(p, k -> [vk])
        pr = Dict(name => fn(p) for (name, fn) in protos)
        @printf("  %-34s kernel f at ×%s = %s → k50 %.4g; emulator %.4g (ratio %.2f); const V_G spectrum %.4g, kernel %.4g",
                id, join(K50_MULT, "/"), join(f, "/"), kK, kE, kE / kK, kA2, kB2)
        for (name, v) in pr
            @printf("; %s %.4g (ratio %.2f)", name, v, v / kK)
        end
        println()
        @printf("      log(E/K) = %.3f = log(E/constS) %.3f + log(constS/constK) %.3f [V_G at the move] + log(constK/K) %.3f [decline during the move]\n",
                log(kE / kK), log(kE / kA2), log(kA2 / kB2), log(kB2 / kK))
        push!(out, (; id, kK, kE, kA2, kB2, pr))
    end
    out
end

# ---------------------------------------------------------------------------------------
# confirmation of the LD prototype on E6′ inputs not used for discovery (sets other than D
# and Dctl): P2000 with ldratio = ld_ratio(p) against the kernel's 40 E6′ runs

function e6p_compare(label, P, f, rate; m)
    M = f .> 0.5; C = P .> 0.5
    @printf("  %-28s n %d: emulator alive & kernel dead %d, emulator dead & kernel alive %d, accuracy %.3f, mean |f − P2000| %.3f\n",
            label, length(P), count(C .& .!M), count(.!C .& M), mean(C .== M), mean(abs.(f .- P)))
    rk = reduce(vcat, [fill(rate[i], RK) for i in eachindex(rate)])
    yk = reduce(vcat, [[k <= round(Int, f[i] * RK) for k in 1:RK] for i in eachindex(rate)])
    xk, hk = Base.invokelatest(() -> m.x50_rule(rk, yk))
    xe, he = Base.invokelatest(() -> m.x50_rule(rate, P))
    @printf("      pooled x50: kernel %.4f (%s), emulator %.4f (%s), ratio %.3f\n", xk, hk, xe, he, xe / xk)
end

function confirm_e6p(; sets = ("", "C", "Cctl"))
    m = e6p_module()
    rows = filter(r -> r.set in sets, read_inputs_csv())
    ps = [readp(r.id) for r in rows]
    res = sweep(collect(zip(rows, ps)); id = x -> x[1].id, cost = x -> x[2].NoQTL, progress = 20,
                cache = joinpath(HERE, "spectrum_diagnosis_proto_e6p.csv")) do (r, p)
        rr = ld_ratio(p)
        (; ldratio = rr, P = proto_P2000(p; kw = (; ldratio = rr)))
    end
    f = [r.f for r in rows]; Ps = [r.P for r in rows]; Pp = [x.P for x in res]
    rate = [p.MoveFitnessSD / sqrt(p.PhenotypeVar) for p in ps]
    println("confirmation on $(length(rows)) E6′ inputs not used for discovery (sets $(join(sets, ", ")))")
    e6p_compare("spectrum (as in E6′)", Ps, f, rate; m)
    e6p_compare("spectrum + ldratio", Pp, f, rate; m)
    strata = [first(split(r.id, "/")) for r in rows]
    for s in sort(unique(strata))
        k = strata .== s
        @printf("    %s: discordant spectrum %d → spectrum+ldratio %d (of %d)\n", s, count((Ps[k] .> 0.5) .!= (f[k] .> 0.5)),
                count((Pp[k] .> 0.5) .!= (f[k] .> 0.5)), count(k))
    end
end

# ---------------------------------------------------------------------------------------
# fresh confirmation: new published inputs (none excluded by emulator_spectrum.jl, none in its
# E1–E4 sample, its E6′ pool or E6′), NFRESH per stratum with MoveFitnessSD/kc in [0.3, 3]
# (E6′'s prefilter), each run FRESH_R times in the kernel exactly as E6′ with fresh seeds

const NFRESH, FRESH_R, FRESH_SEED, FRESH_KSEED = 8, 20, 20261119, 60261114000

function fresh_inputs()
    m = e6p_module()
    Base.invokelatest() do
        ids, ps, strat = m.read_inputs()
        ex = m.excluded_ids()
        S = m.main_sample(ids, strat, ex)
        pool = m.e6p_pool(ids, ps, strat, union(ex, Set(ids[S])))
        used = union(ex, Set(ids[S]), Set(ids[pool]))
        rng = MersenneTwister(FRESH_SEED)
        pick = Int[]
        for s in 1:8
            el = [i for i in eachindex(ids) if strat[i] == s && !(ids[i] in used)]
            got = 0
            for i in el[randperm(rng, length(el))]
                kc = drift_critical_rate(ps[i]; ne = :falconer, solver = :iteration).k
                (kc > 0 && 0.3 <= ps[i].MoveFitnessSD / kc <= 3.0) || continue
                push!(pick, i); got += 1
                got == NFRESH && break
            end
        end
        ids[pick]
    end
end

function confirm_fresh()
    file = joinpath(HERE, "spectrum_diagnosis_fresh_inputs.txt")
    if !isfile(file)
        write(file, join(fresh_inputs(), "\n") * "\n")
    end
    fids = filter(!isempty, readlines(file))
    println("fresh inputs: $(length(fids)) ($file)")
    jobs = [(; key = "$(id)#$r", id, seed = FRESH_KSEED + 100j + r) for (j, id) in enumerate(fids) for r in 1:FRESH_R]
    kr = sweep(jobs; id = j -> j.key, cost = j -> readp(j.id).NoQTL, progress = 200,
               cache = joinpath(HERE, "spectrum_diagnosis_fresh_kernel.csv")) do j
        o = simulate(kernel_params(readp(j.id)); seed = j.seed, init = :equilibrium, burnin = KBURN)
        (; alive = !o.extinct, extinct_gen = o.extinct_gen)
    end
    ps = [readp(id) for id in fids]
    em = sweep(collect(zip(fids, ps)); id = x -> x[1], cost = x -> x[2].NoQTL, progress = 20,
               cache = joinpath(HERE, "spectrum_diagnosis_fresh_emulator.csv")) do (id, p)
        rr = ld_ratio(p)
        (; P_spectrum = proto_P2000(p), ldratio = rr, P_proto = proto_P2000(p; kw = (; ldratio = rr)))
    end
    f = [count(kr[(j-1)*FRESH_R+r].alive for r in 1:FRESH_R) / FRESH_R for j in eachindex(fids)]
    rate = [p.MoveFitnessSD / sqrt(p.PhenotypeVar) for p in ps]
    m = e6p_module()
    @printf("kernel fraction alive: mean %.3f; inputs with none %d, all %d, some %d\n", mean(f), count(==(0), f), count(==(1), f),
            count(x -> 0 < x < 1, f))
    e6p_compare_R("spectrum", [x.P_spectrum for x in em], f, rate, FRESH_R; m)
    e6p_compare_R("spectrum + ldratio", [x.P_proto for x in em], f, rate, FRESH_R; m)
end

function e6p_compare_R(label, P, f, rate, R; m)
    M = f .> 0.5; C = P .> 0.5
    @printf("  %-22s n %d: emulator alive & kernel dead %d, emulator dead & kernel alive %d, accuracy %.3f, mean |f − P2000| %.3f\n",
            label, length(P), count(C .& .!M), count(.!C .& M), mean(C .== M), mean(abs.(f .- P)))
    rk = reduce(vcat, [fill(rate[i], R) for i in eachindex(rate)])
    yk = reduce(vcat, [[k <= round(Int, f[i] * R) for k in 1:R] for i in eachindex(rate)])
    xk, hk = Base.invokelatest(() -> m.x50_rule(rk, yk))
    xe, he = Base.invokelatest(() -> m.x50_rule(rate, P))
    @printf("      pooled x50: kernel %.4f (%s), emulator %.4f (%s), ratio %.3f\n", xk, hk, xe, he, xe / xk)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "inputs" : ARGS[1]
    if mode == "inputs"
        make_inputs()
    elseif mode == "kernel"
        run_kernel(ARGS[2], length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10)
    elseif mode == "analyse"
        analyse(ARGS[2], length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10)
    elseif mode == "k50"
        k50_decomposition()
    elseif mode == "confirm"
        confirm_e6p()
    elseif mode == "fresh"
        confirm_fresh()
    end
end
