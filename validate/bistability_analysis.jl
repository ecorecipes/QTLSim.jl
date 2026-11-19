# Analysis, verdicts, figures and the pilot of validate/bistability.jl (included by it).
# Every rule here is in validate/bistability_criteria.md.

const HI, LO = 0.9, 0.1            # a state persists / dies "reliably"
const DROP = 0.25                   # B2: a fall of the survival fraction between neighbours that counts as a violation

"The cached rows of a finished cell, or nothing."
function load_cell(a, K, R)
    name = cellname(a, K, R)
    isfile(joinpath(RAW, name * ".done")) || return nothing
    ci = cellinfo(a, K, R)
    jobs = [(; s..., st = i) for (i, s) in enumerate(states(ci))]
    rows = sweep(j -> error("cache of $name incomplete"), jobs; cache = joinpath(RAW, name * ".csv"), id = j -> "$(j.ir)-$(j.il)-$(j.iN)-$(j.iV)")
    pl = parse.(Float64, split(readline(joinpath(RAW, name * ".pool"))))
    (; a, K, R, name, ci, rows, pool = (VG = pl[1], sdVG = pl[2], N = pl[3], lag = pl[4]), secs = parse(Float64, readline(joinpath(RAW, name * ".done"))))
end

"The four predictions for every main-grid state of a cell and their boundaries in V₀, cached."
function predictions(c)
    ci = c.ci
    main = [r for r in c.rows if r.ir <= length(RATEF)]
    jobs = [(; r.ir, r.il, r.iN, r.k, r.Δ0, r.N0) for r in main if r.iV == 1]
    b = sweep(jobs; cache = joinpath(RAW, c.name * "_boundaries.csv"), id = j -> "$(j.ir)-$(j.il)-$(j.iN)") do j
        (; j.ir, j.il, j.iN, (p => boundary(p, ci.p, ci.Veq, j.Δ0, j.N0; k = j.k) for p in PREDS)...)
    end
    s = sweep(main; cache = joinpath(RAW, c.name * "_predictions.csv"), id = j -> "$(j.ir)-$(j.il)-$(j.iN)-$(j.iV)") do j
        (; j.ir, j.il, j.iN, j.iV, (p => predict(p, ci.p, j.Δ0, j.V0, j.N0; k = j.k) for p in PREDS)...,
         pd_breeders_T = pd_fate(ci.p, j.Δ0, j.V0, j.N0; k = j.k, ne = :breeders, T = HORIZONS[(c.a, c.K, c.R)]) == 0)
    end
    bd = Dict((r.ir, r.il, r.iN) => r for r in b)
    sd = Dict((r.ir, r.il, r.iN, r.iV) => r for r in s)
    (; bd, sd)
end

"Logistic regression of binomial counts y of n on x: (a, b, se_b, status). status :ok, :separated, :undefined."
function logit_fit(x::Vector{Float64}, y::Vector{Int}, n::Vector{Int})
    length(x) >= 10 || return (NaN, NaN, NaN, :undefined)
    (sum(y) == 0 || sum(y) == sum(n)) && return (NaN, NaN, NaN, :undefined)
    # perfect separation: every state below some x has no survivor and every state above only survivors
    o = sortperm(x)
    if any(i -> all(y[o[1:i]] .== 0) && all(y[o[i+1:end]] .== n[o[i+1:end]]), 1:length(x)-1)
        i = findfirst(i -> all(y[o[1:i]] .== 0) && all(y[o[i+1:end]] .== n[o[i+1:end]]), 1:length(x)-1)
        return (-(x[o[i]] + x[o[i+1]]) / 2, 1.0, NaN, :separated)     # a = −x50 with b = 1, so that −a/b is the gap's midpoint
    end
    a, b = 0.0, 0.0
    H = zeros(2, 2)
    for _ in 1:200
        η = clamp.(a .+ b .* x, -30, 30)
        p = 1 ./ (1 .+ exp.(-η))
        g = [sum(y .- n .* p), sum((y .- n .* p) .* x)]
        w = n .* p .* (1 .- p) .+ 1e-12
        H = [sum(w) sum(w .* x); sum(w .* x) sum(w .* x .^ 2)]
        d = H \ g
        st = min(1.0, 5 / max(maximum(abs, d), 1e-12))
        a += st * d[1]; b += st * d[2]
        maximum(abs, d) < 1e-10 && break
    end
    abs(b) > 200 && return (a, b, NaN, :separated)
    (a, b, sqrt(inv(H)[2, 2]), :ok)
end

frac(r) = r.alive / r.reps

function analyse()
    cells = filter(!isnothing, [load_cell(c...) for c in CELLS])
    isempty(cells) && return println("no finished cell")
    say("validate/bistability.jl — individual-based test of the bistability of Xu and Shen (2026, bioRxiv preprint), criteria: bistability_criteria.md")
    sayf("cells run: %d of %d (%s); %d replicates per state; wall time of the cells %.0f min", length(cells), length(CELLS),
         join((@sprintf("%s (T = %d)", c.name, HORIZONS[(c.a, c.K, c.R)]) for c in cells), ", "), REPS, sum(c.secs for c in cells) / 60)
    say("cells dropped before the run for cost (criteria): ", join((cellname(c...) for c in DROPPED), ", "))
    notrun = [cellname(c...) for c in CELLS if !isfile(joinpath(RAW, cellname(c...) * ".done"))]
    isempty(notrun) || say("cells not finished: ", join(notrun, ", "))
    preds = Dict(c.name => predictions(c) for c in cells)

    say("\n1. Cells: critical rates, steady state, pool populations, and the helper's realised initial states")
    sayf("  %-12s %8s %8s %7s %8s %7s | %8s %8s | %-28s | %s", "cell", "k_c(br)", "k_c(es)", "λ*", "V_eq", "N_eq", "v_c K1", "v_c K2", "pool: V_G (sd), N, V_G/V_eq", "max |ΔN|, |V/V₀−1|, |Δlag|")
    for c in cells
        ci = c.ci
        sayf("  %-12s %8.4f %8.4f %7.3f %8.4f %7.0f | %8.4f %8.4f | %.4f (%.4f), %5.0f, %.2f       | %.0f, %.1e, %.1e", c.name, ci.kc, ci.kcf, ci.λc, ci.Veq, ci.Neq,
             xu_shen_critical_rate(xs_map(ci.p)).v, xu_shen_critical_rate(xs_map(ci.p; Kmap = :K2)).v, c.pool.VG, c.pool.sdVG, c.pool.N, c.pool.VG / ci.Veq,
             maximum(r.maxerrN for r in c.rows), maximum(r.maxerrV for r in c.rows), maximum(r.maxerrL for r in c.rows))
    end
    say("  k_c(br): drift_critical_rate with N_e = breeders (the rates of the design are 0.6 and 0.85 of it); k_c(es): with effective_size (the default).")

    say("\n2. Survival by state: alive of $(REPS) at the cell's horizon T. Rows V₀/V_eq = ", join(reverse(VFRAC), ", "), " (top to bottom); columns N₀/K = ", join(NFRAC, ", "))
    for c in cells
        say("  ", c.name)
        for ir in 1:length(RATEF)+2
            lags = ir <= length(RATEF) ? (1:length(LAGF)) : (1:1)
            rf = ir <= length(RATEF) ? RATEF[ir] : ir == length(RATEF) + 1 ? 0.0 : NULLF
            say(@sprintf("    k = %.2f k_c%s", rf, ir == length(RATEF) + 1 ? " (null: no move)" : ir == length(RATEF) + 2 ? " (null: above the critical rate)" : ""),
                "      ", join((@sprintf("%-22s", "Δ₀ = $(LAGF[il]) λ*") for il in lags), "  "))
            for iV in length(VFRAC):-1:1
                say("      ", join((join((@sprintf("%3d", only(r.alive for r in c.rows if r.ir == ir && r.il == il && r.iN == iN && r.iV == iV)) for iN in 1:length(NFRAC)), " ") * "    "
                                    for il in lags), "    "))
            end
        end
    end

    arms = unique(c.a for c in cells)
    verdicts = String[]

    say("\n3. Null arms")
    for a in arms
        ok0 = ok2 = true; n0 = n2 = 0
        for c in cells
            c.a === a || continue
            for r in c.rows
                if r.rf == 0.0 && r.ir == length(RATEF) + 1
                    NFRAC[r.iN] >= 0.05 && (n0 += 1; ok0 &= frac(r) >= HI)
                elseif r.ir == length(RATEF) + 2
                    n2 += 1; ok2 &= frac(r) <= LO
                end
            end
        end
        small = [frac(r) for c in cells if c.a === a for r in c.rows if r.ir == length(RATEF) + 1 && NFRAC[r.iN] < 0.05]
        v0 = ok0 ? "PASS" : "FAIL"; v2 = ok2 ? "PASS" : "FAIL"
        sayf("  arm %s: no move — ≥ %.0f%% alive in every state with N₀ ≥ 0.05K (%d states): %s (states with N₀ = 0.02K: survival %.2f–%.2f); k = %g k_c — ≤ %.0f%% alive in every state (%d): %s",
             a, 100HI, n0, v0, minimum(small), maximum(small), NULLF, 100LO, n2, v2)
        push!(verdicts, "null, no move, arm $a: $v0", "null, above critical, arm $a: $v2")
    end

    say("\n4. B1 — bistability: at one rate below the critical rate, a state with ≥ 90% survival and a state with ≤ 10% (unit: cell × rate; states = N₀ × V₀ × Δ₀)")
    for a in arms
        units = 0; yes = 0
        for c in cells
            c.a === a || continue
            for ir in 1:length(RATEF)
                rs = [r for r in c.rows if r.ir == ir]
                hi, lo = count(r -> frac(r) >= HI, rs), count(r -> frac(r) <= LO, rs)
                units += 1; yes += hi > 0 && lo > 0
                sayf("    %-12s k = %.2f k_c: %2d states ≥ 90%%, %2d ≤ 10%%, %2d between: %s", c.name, RATEF[ir], hi, lo, length(rs) - hi - lo, hi > 0 && lo > 0 ? "bistable" : "not")
                for il in 1:length(LAGF)
                    r2 = [r for r in rs if r.il == il]
                    sayf("        Δ₀ = %.1f λ*: %2d ≥ 90%%, %2d ≤ 10%%", LAGF[il], count(r -> frac(r) >= HI, r2), count(r -> frac(r) <= LO, r2))
                end
            end
        end
        v = units < 4 ? "INCONCLUSIVE (fewer than 4 units)" : yes >= 0.75 * units ? "PASS" : "FAIL"
        sayf("  arm %s: bistable in %d of %d units: B1 %s", a, yes, units, v)
        push!(verdicts, "B1, arm $a: $v ($yes/$units)")
    end

    say("\n5. B2 — monotonicity: survival does not fall by more than $(DROP) from a state to its neighbour with more N₀ or more V₀")
    for a in arms
        pairs = Dict(:N => [0, 0], :V => [0, 0])
        for c in cells
            c.a === a || continue
            f = Dict((r.ir, r.il, r.iN, r.iV) => frac(r) for r in c.rows if r.ir <= length(RATEF))
            for ir in 1:length(RATEF), il in 1:length(LAGF), i in 1:5, j in 1:4
                pairs[:N][2] += 1; pairs[:N][1] += f[(ir, il, j + 1, i)] < f[(ir, il, j, i)] - DROP
                pairs[:V][2] += 1; pairs[:V][1] += f[(ir, il, i, j + 1)] < f[(ir, il, i, j)] - DROP
            end
        end
        tot = pairs[:N][1] + pairs[:V][1]; np = pairs[:N][2] + pairs[:V][2]
        v = tot <= 0.05 * np ? "PASS" : "FAIL"
        sayf("  arm %s: violations %d of %d pairs along N₀, %d of %d along V₀ (%.1f%% in all): B2 %s", a, pairs[:N]..., pairs[:V]..., 100tot / np, v)
        push!(verdicts, "B2, arm $a: $v ($tot/$np)")
    end
    # how much does N₀ matter? the mean change of survival from the smallest to the largest N₀ at fixed V₀
    for a in arms
        d = [frac(only(r for r in c.rows if r.ir == ir && r.il == il && r.iN == 5 && r.iV == iV)) - frac(only(r for r in c.rows if r.ir == ir && r.il == il && r.iN == 1 && r.iV == iV))
             for c in cells if c.a === a for ir in 1:length(RATEF) for il in 1:length(LAGF) for iV in 1:5]
        dv = [frac(only(r for r in c.rows if r.ir == ir && r.il == il && r.iN == iN && r.iV == 5)) - frac(only(r for r in c.rows if r.ir == ir && r.il == il && r.iN == iN && r.iV == 1))
              for c in cells if c.a === a for ir in 1:length(RATEF) for il in 1:length(LAGF) for iN in 1:5]
        sayf("  arm %s (descriptive): mean change of survival from N₀ = 0.02K to 0.9K at fixed V₀ %+.3f (max %+.2f); from V₀ = 0.1 to 1.0 V_eq at fixed N₀ %+.3f", a, mean(d), maximum(d), mean(dv))
    end

    say("\n6. B3 — the deterministic boundary against the majority outcome of each state (ties excluded), and a logistic regression of survival on the")
    say("   signed log distance d = ln(V₀/V_b) to the predicted boundary in V₀ (states where the prediction has a boundary inside (V_eq/1000, 3V_eq)).")
    say("   x50 = e^{d50}: the kernel's 50% point as a multiple of the predicted boundary variance (> 1: the extinction basin is larger than predicted).")
    for a in arms
        say("  arm ", a)
        for p in (PREDS..., :pd_breeders_T)
            hit = n = ties = 0; nper = [0, 0]
            xs = Float64[]; ys = Int[]; ns = Int[]
            percell = String[]
            for c in cells
                c.a === a || continue
                pr = preds[c.name]
                h = m = 0
                x1 = Float64[]; y1 = Int[]; n1 = Int[]
                for r in c.rows
                    r.ir <= length(RATEF) || continue
                    pp = getproperty(pr.sd[(r.ir, r.il, r.iN, r.iV)], p)
                    nper[1] += pp; nper[2] += 1
                    if 2r.alive == r.reps
                        ties += 1
                    else
                        m += 1; h += (2r.alive > r.reps) == pp
                    end
                    if p !== :pd_breeders_T
                        vb = getproperty(pr.bd[(r.ir, r.il, r.iN)], p)
                        if 0 < vb < Inf
                            push!(x1, log(r.V0 / vb)); push!(y1, r.alive); push!(n1, r.reps)
                        end
                    end
                end
                hit += h; n += m
                append!(xs, x1); append!(ys, y1); append!(ns, n1)
                fa, fb, fse, fst = logit_fit(x1, y1, n1)
                push!(percell, @sprintf("%s %d/%d%s", c.name, h, m, fst === :undefined ? "" : @sprintf(", x50 %.2f", exp(-fa / fb))))
            end
            name = p === :pd_breeders_T ? "P-D, N_e = breeders, read at T (secondary)" : PREDNAME[p]
            v = n < 50 ? "INCONCLUSIVE (fewer than 50 states)" : hit >= 0.8n ? "PASS" : "FAIL"
            fa, fb, fse, fst = logit_fit(xs, ys, ns)
            reg = p === :pd_breeders_T ? "" :
                  fst === :undefined ? "; regression undefined ($(length(xs)) states with a boundary)" :
                  fst === :separated ? @sprintf("; regression: separated (slope > 0), x50 %.2f, %d states", exp(-fa / fb), length(xs)) :
                  @sprintf("; regression: slope %.2f [%.2f, %.2f] per ln V₀, x50 = %.2f (d50 %+.2f), %d states", fb, fb - 1.96fse, fb + 1.96fse, exp(-fa / fb), -fa / fb, length(xs))
            sayf("    %-42s predicts persistence in %3d of %3d states; majority outcome matched in %3d of %3d (%.0f%%; %d ties)%s%s", name, nper..., hit, n, 100hit / max(n, 1), ties,
                 p === :pd_breeders_T ? "" : ": B3 " * v, reg)
            say("        per cell: ", join(percell, "; "))
            p === :pd_breeders_T || push!(verdicts, "B3, arm $a, $(PREDNAME[p]): $v ($hit/$n)")
        end
    end

    say("\n7. Where the kernel's boundary is: for each cell, rate and lag, the V₀/V_eq at which half survive (logistic in ln V₀ over all N₀), against the predictions' V_b/V_eq at N₀ = 0.4K")
    sayf("  %-12s %5s %5s | %-22s | %8s %8s %8s %8s", "cell", "k/k_c", "Δ₀/λ*", "kernel V50/V_eq", "P-D br", "P-D es", "XS K1", "XS K2")
    for c in cells, ir in 1:length(RATEF), il in 1:length(LAGF)
        rs = [r for r in c.rows if r.ir == ir && r.il == il]
        fa, fb, fse, fst = logit_fit([log(r.V0 / c.ci.Veq) for r in rs], [r.alive for r in rs], [r.reps for r in rs])
        ker = fst === :undefined ? (all(r -> r.alive == r.reps, rs) ? "all alive" : all(r -> r.alive == 0, rs) ? "all dead" : sum(r.alive for r in rs) > sum(r.reps for r in rs) / 2 ? "< 0.1 (mostly alive)" : "> 1 (mostly dead)") :
              @sprintf("%.3f%s", exp(-fa / fb), fst === :separated ? " (sep.)" : @sprintf(" (slope %.1f)", fb))
        b = preds[c.name].bd[(ir, il, 4)]
        fmt(x) = x == 0 ? "< 0.001" : x == Inf ? "> 3" : @sprintf("%.3f", x / c.ci.Veq)
        sayf("  %-12s %5.2f %5.1f | %-22s | %8s %8s %8s %8s", c.name, RATEF[ir], LAGF[il], ker, fmt(b.pd_breeders), fmt(b.pd_falconer), fmt(b.xs_K1), fmt(b.xs_K2))
    end

    say("\n8. Survivors' state at the horizon against the steady state (states with ≥ 90% survival; mean over states): N, V_G/V_eq, lag/λ(steady, P-D)")
    for c in cells, ir in 1:length(RATEF)
        rs = [r for r in c.rows if r.ir == ir && frac(r) >= HI]
        isempty(rs) && continue
        sayf("    %-12s k = %.2f k_c: N %.0f, V_G/V_eq %.2f (pool %.2f), lag %.3f (k(V+V_s)/V at the survivors' V_G: %.3f)", c.name, RATEF[ir], mean(r.Nend for r in rs),
             mean(r.Vend for r in rs) / c.ci.Veq, c.pool.VG / c.ci.Veq, mean(r.lagend for r in rs),
             RATEF[ir] * c.ci.kc * (mean(r.Vend for r in rs) + selection_variance(c.ci.p)) / mean(r.Vend for r in rs))
    end

    say("\nVERDICTS")
    foreach(v -> say("  ", v), verdicts)
    figures(cells, preds)
    write(joinpath(HERE, "bistability.log"), String(take!(out)))
end

# ── figures ───────────────────────────────────────────────────────────────────────────────

function figures(cells, preds)
    @eval using CairoMakie
    Base.invokelatest(_figures, cells, preds)
end

function _figures(cells, preds)
    mkpath(FIGS)
    ink, muted = "#1f2328", "#6b7280"
    # boundary lines: four series, fixed order, each with its own dash pattern (identity is not colour alone)
    style = Dict(:pd_breeders => ("#b45309", :solid), :pd_falconer => ("#b45309", :dash), :xs_K1 => ("#7c3aed", :solid), :xs_K2 => ("#7c3aed", :dash))
    for c in cells
        fig = Figure(size = (1150, 720), fontsize = 13, backgroundcolor = :white)
        Label(fig[0, 1:3], @sprintf("%s arm, K = %d, R = %g: fraction of %d kernel populations alive %d generations after the move starts, by initial state",
                                   c.a === :G ? "Gaussian" : "QTSim-type", c.K, c.R, REPS, HORIZONS[(c.a, c.K, c.R)]); fontsize = 15, color = ink, halign = :left)
        hm = nothing
        for ir in 1:length(RATEF), il in 1:length(LAGF)
            ax = Axis(fig[ir, il]; title = @sprintf("k = %.2f k_c,  Δ₀ = %.1f λ*", RATEF[ir], LAGF[il]), titlecolor = ink, titlesize = 13,
                      xlabel = ir == length(RATEF) ? "initial size N₀ / K" : "", ylabel = il == 1 ? "initial variance V₀ / V_eq" : "",
                      xticks = (1:5, collect(string.(NFRAC))), yticks = (1:5, collect(string.(VFRAC))), xgridvisible = false, ygridvisible = false,
                      xticklabelcolor = muted, yticklabelcolor = muted, xlabelcolor = ink, ylabelcolor = ink)
            Z = [frac(only(r for r in c.rows if r.ir == ir && r.il == il && r.iN == iN && r.iV == iV)) for iN in 1:5, iV in 1:5]
            hm = heatmap!(ax, 1:5, 1:5, Z; colormap = :Blues, colorrange = (0, 1))
            for iN in 1:5, iV in 1:5
                text!(ax, iN, iV; text = @sprintf("%.2f", Z[iN, iV]), align = (:center, :center), fontsize = 11, color = Z[iN, iV] > 0.6 ? :white : ink)
            end
            # predicted boundaries: the position of V_b on the (logarithmically interpolated) V₀ axis
            lv = log.(collect(VFRAC))
            ypos(v) = v <= 0 ? 0.5 : v == Inf ? 5.5 : (x = log(v / c.ci.Veq); x <= lv[1] ? max(0.5, 1 + (x - lv[1]) / (lv[2] - lv[1])) : x >= lv[end] ? min(5.5, 5 + (x - lv[end]) / (lv[end] - lv[end-1])) :
                       (j = searchsortedlast(lv, x); j + (x - lv[j]) / (lv[j+1] - lv[j])))
            for p in PREDS
                ys = [ypos(getproperty(preds[c.name].bd[(ir, il, iN)], p)) for iN in 1:5]
                lines!(ax, 1:5, ys; color = style[p][1], linestyle = style[p][2], linewidth = 2)
            end
            limits!(ax, 0.5, 5.5, 0.5, 5.5)
        end
        Colorbar(fig[1:2, 4], hm; label = "fraction alive", labelcolor = ink, ticklabelcolor = muted)
        Legend(fig[3, 1:3], [LineElement(color = style[p][1], linestyle = style[p][2], linewidth = 2) for p in PREDS], [PREDNAME[p] * " (persists above the line)" for p in PREDS];
               orientation = :horizontal, framevisible = false, labelcolor = ink, nbanks = 2)
        save(joinpath(FIGS, c.name * ".png"), fig)
    end
end

# ── pilot ─────────────────────────────────────────────────────────────────────────────────

"""
The pilot (criteria file, "Pilot"): the cell Q-K1000-R2 with its own seeds and a horizon of 1500 generations, to check the helper,
measure the cost and choose T. Prints only what the criteria say the pilot may be used for: realised initial states, cost, and the
distribution of the times of death (not survival by state). Also times 20 runs of G-K1000-R2 and G-K3000-R2 (cost only).
"""
function pilot()
    HORIZON[] = 1500
    a, K, R = :Q, 1000, 2.0
    ci = cellinfo(a, K, R)
    t0 = time()
    pools = Vector{Any}(undef, REPS)
    Threads.@threads :greedy for r in 1:REPS
        pools[r] = pool_population(ci.p, PILOT_SEED0 + r)
    end
    @printf("pilot %s: pool of %d populations in %.0f s; V_G %.4f (sd %.4f) against V_eq %.4f; N %.0f\n", cellname(a, K, R), REPS, time() - t0,
            mean(x.VG for x in pools), std([x.VG for x in pools]), ci.Veq, mean(x.N for x in pools))
    sts = states(ci)
    res = Vector{Any}(undef, length(sts))
    cpu = Threads.Atomic{Float64}(0.0)
    t1 = time()
    Threads.@threads :greedy for i in eachindex(sts)
        j = sts[i]
        tt = @elapsed res[i] = [one_run(ci.p, pools[r].G, min(j.N0, pools[r].N), j.V0, j.Δ0, j.k, PILOT_SEED0 + 10^4 * i + 100 + r) for r in 1:10]
        Threads.atomic_add!(cpu, tt)
    end
    all = reduce(vcat, res)
    @printf("  %d states × 10 runs to %d generations: %.0f s wall, %.0f thread-seconds\n", length(sts), HORIZON[], time() - t1, cpu[])
    @printf("  helper: max |ΔN| %.0f, max |V/V₀ − 1| %.1e, max |Δlag| %.1e\n", maximum(abs(r.errN) for r in all), maximum(abs(r.errV) for r in all if !isnan(r.errV)),
            maximum(abs(r.errL) for r in all if !isnan(r.errL)))
    for (lab, sel) in (("main grid", i -> sts[i].ir <= length(RATEF)), ("null, no move", i -> sts[i].rf == 0.0), ("null, $(NULLF) k_c", i -> sts[i].rf == NULLF))
        ts = [r.t for i in eachindex(sts) if sel(i) for r in res[i] if !r.alive]
        n = sum(length(res[i]) for i in eachindex(sts) if sel(i))
        isempty(ts) && (@printf("  %-16s %d runs, no death\n", lab, n); continue)
        q(x) = sort(ts)[clamp(round(Int, x * length(ts)), 1, length(ts))]
        @printf("  %-16s %d runs, %d deaths; time of death: median %d, 90%% %d, 99%% %d, max %d; deaths after 300: %d, after 500: %d, after 1000: %d\n", lab, n, length(ts),
                q(0.5), q(0.9), q(0.99), maximum(ts), count(>(300), ts), count(>(500), ts), count(>(1000), ts))
    end
    HORIZON[] = 500
    for (a2, K2) in ((:G, 1000), (:G, 3000), (:Q, 3000))
        c2 = cellinfo(a2, K2, 2.0)
        tp = @elapsed pl = pool_population(c2.p, PILOT_SEED0 + 999)
        tr = @elapsed for r in 1:20
            one_run(c2.p, pl.G, round(Int, 0.9K2), c2.Veq, 0.0, 0.0, PILOT_SEED0 + 5000 + r)
        end
        @printf("  cost only, %s-K%d-R2: one pool population (%d generations) %.1f s; one run of %d generations at k = 0 from N₀ = 0.9K: %.2f s\n", a2, K2, BURN, tp, 500, tr / 20)
    end
end
