# Development of the dynamic linkage disequilibrium of the `:spectrum` closure
# (`spectrum = (; ld = :dynamic)`, src/spectrum.jl), on the discovery inputs of
# validate/spectrum_diagnosis.jl (sets D and Dctl) and their kernel trajectories. Exploratory:
# what was tried and why is in validate/spectrum_ld_dev.log; the pre-specified validation is
# validate/spectrum_ld.jl under validate/spectrum_ld_criteria.md.
#
# Run from the package root (vignettes environment):
#   julia -t 8 --project=vignettes validate/spectrum_ld_dev.jl static  [SET]   # V_G, genic, LD at the move
#   julia -t 8 --project=vignettes validate/spectrum_ld_dev.jl move    [SET]   # V_G, genic, LD during the move
#   julia -t 8 --project=vignettes validate/spectrum_ld_dev.jl driven  [SET]   # genic and LD driven by each run's own lag and survivors
#   julia -t 8 --project=vignettes validate/spectrum_ld_dev.jl bulmer  [SET]   # the LD recursion given the kernel's genic variance
# SET defaults to D+Dctl. The kernel trajectories are the diagnosis's raw files (git-ignored):
# ENV RAWDIR, default ../wt-spectrum-diagnosis/validate/spectrum_diagnosis_raw next to this
# worktree, or regenerate them with `spectrum_diagnosis.jl kernel D+Dctl`.

const DEVRAW = get(ENV, "RAWDIR", normpath(joinpath(@__DIR__, "..", "..", "wt-spectrum-diagnosis", "validate",
                                                     "spectrum_diagnosis_raw")))
ENV["RAWDIR"] = DEVRAW
include(joinpath(@__DIR__, "spectrum_diagnosis.jl"))
using Distributions

# the variants compared in development (keyword tuples of `run_moments(…; spectrum)`)
const VARIANTS = [
    "spectrum" => (;),
    "ld" => (; ld = :dynamic),
    "ld-meanfield" => (; ld = :dynamic, ldbins = 1),
    "ld-noturn" => (; ld = :dynamic, turnover = false),
    "ld-nobarton" => (; ld = :dynamic, linkage = :none),
]

# kernel means over runs alive: V_G (every generation) and genic (every 10) in a window
function kernel_window(runs, field, x, lo, hi)
    v = Float64[]
    for r in runs
        g = getfield(r, x)
        for k in eachindex(g)
            lo <= g[k] - MOVE <= hi && isfinite(getfield(r, field)[k]) && push!(v, getfield(r, field)[k])
        end
    end
    isempty(v) ? NaN : mean(v)
end

function static(set)
    rows = inset(read_inputs_csv(), set)
    res = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        r = rows[i]
        p = readp(r.id)
        runs = load_runs(r.id, 10)
        kV = kernel_window(runs, :VG, :gen, -300, -1)
        kg = kernel_window(runs, :genic, :xg, -300, -10)
        out = Any[]
        for (name, kw) in VARIANTS
            e = spectrum_equilibrium(p; kw...)
            g = haskey(e, :genic) ? e.genic : e.VG
            push!(out, (name, e.VG, g))
        end
        me = multilocus_equilibrium(p; drift = true, loci = :multilocus)
        res[i] = (r, p, kV, kg, out, me)
    end
    @printf("%-34s %-4s %5s  %-28s", "input", "set", "L", "kernel V_G genic ratio")
    for (name, _) in VARIANTS
        @printf("  %-26s", name)
    end
    println("  multilocus")
    rk = Dict(name => Float64[] for (name, _) in VARIANTS); rk["ml"] = Float64[]
    gk = Dict(name => Float64[] for (name, _) in VARIANTS); gk["ml"] = Float64[]
    vk = Dict(name => Float64[] for (name, _) in VARIANTS); vk["ml"] = Float64[]
    for (r, p, kV, kg, out, me) in res
        @printf("%-34s %-4s %5d  %.4f %.4f %.3f      ", r.id, r.set, p.NoQTL, kV, kg, kV / kg)
        for (name, V, g) in out
            @printf("  %.4f %.4f %.3f      ", V, g, V / g)
            push!(rk[name], (V / g) / (kV / kg)); push!(gk[name], g / kg); push!(vk[name], V / kV)
        end
        @printf("  %.4f %.4f %.3f\n", me.VG, me.genic, me.VG / me.genic)
        push!(rk["ml"], (me.VG / me.genic) / (kV / kg)); push!(gk["ml"], me.genic / kg); push!(vk["ml"], me.VG / kV)
    end
    println("\nmedian (range) of model/kernel at the move:")
    for name in [first.(VARIANTS); "ml"]
        @printf("  %-12s V_G %.3f (%.3f–%.3f)  genic %.3f (%.3f–%.3f)  V_G/genic %.3f (%.3f–%.3f)\n", name,
                median(vk[name]), extrema(vk[name])..., median(gk[name]), extrema(gk[name])...,
                median(rk[name]), extrema(rk[name])...)
    end
end

# the moments trajectory as the emulator runs it, with the genic variance and the LD
function ld_traj(p; spectrum = (;), burnin = BURNIN)
    q = Params(p; Sample = 1, NoGenerations = burnin + HORIZON, MoveFitnessSDGen = burnin)
    gen = Float64[]; ld = Float64[]
    rec = Dict{Symbol,Any}(:swap => nothing)
    m = run_moments(q; VG = house_of_gauss_vg(p), erosion = true, loci = :spectrum, spectrum)
    (; t = m.gen .- burnin, VG = m.geno_var, extinct = m.extinct, T = m.extinct ? m.extinct_gen - burnin + 1 : missing,
     lag = m.theta .- m.pheno_mean)
end

const TS = (0, 50, 100, 250, 500, 1000)

function move(set; variants = VARIANTS)
    rows = inset(read_inputs_csv(), set)
    res = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        r = rows[i]
        p = readp(r.id)
        runs = load_runs(r.id, 10)
        death = [x.extinct_gen < 0 ? missing : x.extinct_gen - MOVE + 1 for x in runs]
        nd = count(!ismissing, death)
        Tk = nd > length(runs) / 2 ? median(collect(skipmissing(death))) : missing
        kv = [kernel_window(runs, :VG, :gen, t - 5, t + 5) for t in TS]
        tr = [(name, ld_traj(p; spectrum = kw)) for (name, kw) in variants]
        res[i] = (r, p, Tk, nd, kv, tr)
    end
    ratios = Dict((name, t) => Float64[] for (name, _) in variants for t in TS)
    for (r, p, Tk, nd, kv, tr) in res
        @printf("\n%-34s [%s] L=%d u=%g: kernel dead %d/10, median T %s\n", r.id, r.set, p.NoQTL, p.RecombPoissonMean, nd, fmt(Tk))
        @printf("    %-12s", "kernel V_G")
        for (t, v) in zip(TS, kv)
            @printf("  +%-4d %.4f", t, v)
        end
        println()
        for (name, s) in tr
            @printf("    %-12s", name)
            for (t, v) in zip(TS, kv)
                k = findfirst(==(t), s.t)
                sv = k === nothing ? NaN : s.VG[k]
                @printf("  +%-4d %.4f", t, sv)
                isfinite(v) && isfinite(sv) && push!(ratios[(name, t)], sv / v)
            end
            @printf("   T %s\n", fmt(s.T))
        end
    end
    println("\nmedian model/kernel V_G (runs alive), and the number of inputs:")
    for (name, _) in variants
        @printf("  %-12s", name)
        for t in TS
            v = ratios[(name, t)]
            @printf("  +%-4d %.3f (%d)", t, isempty(v) ? NaN : median(v), length(v))
        end
        println()
    end
end

"""
The spectrum (with `spectrum` options) driven by one kernel run's own lag, S = ω² + V_P and
survivors from generation `t0`, started from its equilibrium shape scaled to the run's genic
variance and, with `ld = :dynamic`, its LD scaled to the run's LD at `t0`. Returns the model's
genic variance and LD at the run's genic sampling generations (no feedback on the lag).
"""
function driven_ld(p, x; spectrum = (;), t0 = MOVE - 300)
    k0 = findfirst(==(t0), x.xg)
    sp = deepcopy(QTLSim.spectrum_equilibrium(p; spectrum...).spectrum)
    g0 = QTLSim.spectrum_vg(sp)
    sp.n .*= x.genic[k0] / g0
    j0 = findfirst(==(t0), x.gen)
    if sp.ld === :dynamic
        D = sum(sp.d)
        D < 0 && (sp.d .*= (x.VG[j0] - x.genic[k0]) / D)
    end
    ω2 = QTLSim.initFitSD(p)^2
    G = fill(NaN, length(x.xg)); Dv = fill(NaN, length(x.xg))
    G[k0] = QTLSim.spectrum_vg(sp); Dv[k0] = QTLSim.spectrum_ld(sp)
    for k in eachindex(x.gen)
        g = x.gen[k]
        g < t0 && continue
        λ = x.θ[k] - x.zbar[k]; S = ω2 + x.VP[k]
        QTLSim.spectrum_select(sp, λ, S)
        QTLSim.spectrum_step!(sp, max(x.nsurv[k], 1))
        QTLSim.spectrum_mutate!(sp, p.PopSize)
        j = findfirst(==(g + 1), x.xg)
        j === nothing || (G[j] = QTLSim.spectrum_vg(sp); Dv[j] = QTLSim.spectrum_ld(sp))
    end
    G, Dv
end

const DTS = (-10, 50, 100, 250, 500)

function driven(set; variants = VARIANTS, R = 4)
    rows = inset(read_inputs_csv(), set)
    res = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        r = rows[i]
        p = readp(r.id)
        runs = load_runs(r.id, R)
        # per variant, mean over runs alive of genic(t)/genic(−300) and LD(t)/genic(−300)
        kg = Dict(t => Float64[] for t in DTS); kd = Dict(t => Float64[] for t in DTS)
        mg = Dict((n, t) => Float64[] for (n, _) in variants for t in DTS)
        md = Dict((n, t) => Float64[] for (n, _) in variants for t in DTS)
        for x in runs
            k0 = findfirst(==(MOVE - 300), x.xg)
            out = [(n, driven_ld(p, x; spectrum = kw)) for (n, kw) in variants]
            for t in DTS
                j = findfirst(==(MOVE + t), x.xg)
                j === nothing && continue
                jv = findfirst(==(MOVE + t), x.gen)
                push!(kg[t], x.genic[j] / x.genic[k0])
                jv === nothing || push!(kd[t], (x.VG[jv] - x.genic[j]) / x.genic[k0])
                for (n, (G, Dv)) in out
                    push!(mg[(n, t)], G[j] / x.genic[k0]); push!(md[(n, t)], Dv[j] / x.genic[k0])
                end
            end
        end
        res[i] = (r, p, kg, kd, mg, md)
    end
    for (r, p, kg, kd, mg, md) in res
        @printf("\n%-34s [%s] L=%d u=%g   genic(t)/genic(−300) and LD(t)/genic(−300), mean over runs alive\n", r.id, r.set, p.NoQTL, p.RecombPoissonMean)
        @printf("    %-12s", "kernel")
        for t in DTS
            @printf("  %+4d %.3f %+.3f", t, wmean(kg[t]), wmean(kd[t]))
        end
        println()
        for (n, _) in variants
            @printf("    %-12s", n)
            for t in DTS
                @printf("  %+4d %.3f %+.3f", t, wmean(mg[(n, t)]), wmean(md[(n, t)]))
            end
            println()
        end
    end
end

# k50 on the diagnosis's six k50 inputs: the kernel's (10 runs at ×0.25/0.35/0.5/0.7/1.0 of
# the input's rate, the diagnosis's cache) against the emulator pipeline's with each variant
function k50(; variants = VARIANTS)
    jobs = kernel_k50_runs(K50_IDS)
    cache = Dict(String(split(l, ",")[1]) => split(l, ",")[2] == "true"
                 for l in readlines(joinpath(HERE, "spectrum_diagnosis_k50_kernel.csv"))[2:end])
    res = Dict{Tuple{String,Float64},Vector{Bool}}()
    for j in jobs
        push!(get!(res, (j.id, j.m), Bool[]), cache[j.key])
    end
    out = Vector{Any}(undef, length(K50_IDS) * length(variants))
    tasks = [(id, v) for id in K50_IDS for v in variants]
    Threads.@threads :greedy for t in eachindex(tasks)
        id, (name, kw) = tasks[t]
        out[t] = proto_k50(readp(id); kw)
    end
    lr = Dict(first(v) => Float64[] for v in variants)
    for (i, id) in enumerate(K50_IDS)
        p = readp(id)
        f = [count(res[(id, m)]) / length(res[(id, m)]) for m in K50_MULT]
        kK = interp50(K50_MULT .* p.MoveFitnessSD, f)
        @printf("  %-34s kernel k50 %.4g", id, kK)
        for (j, (name, _)) in enumerate(variants)
            v = out[(i - 1) * length(variants) + j]
            @printf("  %s %.4g (%.2f)", name, v, v / kK)
            push!(lr[name], log(v / kK))
        end
        println()
    end
    for (name, _) in variants
        @printf("  %-12s geometric mean emulator/kernel k50 %.3f (mean log %.3f)\n", name, exp(mean(lr[name])), mean(lr[name]))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "static" : ARGS[1]
    set = length(ARGS) >= 2 ? ARGS[2] : "D+Dctl"
    vs = length(ARGS) >= 3 ? filter(v -> first(v) in split(ARGS[3], ","), VARIANTS) : VARIANTS
    mode == "static" ? static(set) : mode == "move" ? move(set; variants = vs) : mode == "driven" ? driven(set; variants = vs) : mode == "k50" ? k50(; variants = vs) :
        error("unknown mode $mode")
end
