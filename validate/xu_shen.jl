# Xu and Shen (2026): the preprint's figures from src/xu_shen.jl, next to values read off the
# figures, and a descriptive comparison of its eq. 6 with QTLSim's discrete-generation
# counterpart `drift_critical_rate`.
#
#   Xu, K. and Shen, H. (2026). Eco-evolutionary feedbacks generate bistability in population
#   persistence under gradual environmental change. bioRxiv preprint,
#   doi:10.64898/2026.08.05.743160 (version 1; not peer reviewed).
#
# Usage (from the package root; the vignettes environment has CairoMakie and XLSX):
#     julia --project=vignettes validate/xu_shen.jl
# Writes validate/xu_shen.log and validate/xu_shen_figs/fig{1,2,3,4}.png. About 6 minutes of CPU
# time on one thread. Part 5 (iv) reads the published inputs and table in ../Repository when present.
#
# The comparison targets (figure read-offs) and the parameter mapping are in
# validate/xu_shen_targets.md, committed before this script computed anything they are
# compared with; the script prints that file first. There is no published data, so every
# "target" is a value read off a figure; nothing here is a pass/fail verdict.

using QTLSim, Printf, Random, Statistics
using CairoMakie
include(joinpath(@__DIR__, "published_table.jl"))      # load_master, row_index (XLSX)

const DIR = @__DIR__
const FIGS = joinpath(DIR, "xu_shen_figs")
const ROOT = normpath(joinpath(DIR, "..", ".."))
const LOGIO = open(joinpath(DIR, "xu_shen.log"), "w")
say(s...) = (println(stdout, s...); println(LOGIO, s...); flush(LOGIO); flush(stdout))
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))
const T0 = time()
section(s) = say("\n", "="^100, "\n", s, @sprintf("   [%.0f s]", time() - T0), "\n", "="^100)

const STD = XuShen()                       # r = 0.2, ω² = 50, V_M = 0.1, K = 1000, n_e = 1
om2(x) = 10.0^(-x)                         # ω² at x = log₁₀(1/ω²)

# linear interpolation of y at x0 on a grid (NaN outside or next to a NaN)
function interp(xs, ys, x0)
    i = searchsortedlast(xs, x0)
    (i < 1 || i >= length(xs)) && return i == length(xs) && xs[end] == x0 ? ys[end] : NaN
    ys[i] + (ys[i+1] - ys[i]) * (x0 - xs[i]) / (xs[i+1] - xs[i])
end
# first x at which y crosses the level c (linear interpolation), going up or down
function crossing(xs, ys, c; from = 1)
    for i in from:length(xs)-1
        (isfinite(ys[i]) && isfinite(ys[i+1])) || continue
        (ys[i] - c) * (ys[i+1] - c) <= 0 && ys[i] != ys[i+1] && return xs[i] + (c - ys[i]) * (xs[i+1] - xs[i]) / (ys[i+1] - ys[i])
    end
    NaN
end

# ------------------------------------------------------------------------------------------
function figure1()
    section("1. Figure 1: v_c, V* and N* against the mutational input V_M (standard parameters)")
    VMs = range(0.01, 0.2, length = 96)
    vs = range(0.005, 0.2, length = 118)
    vc = [xu_shen_critical_rate(XuShen(STD; VM)).v for VM in VMs]
    Vst = fill(NaN, length(VMs), length(vs)); Nst = copy(Vst)
    for (i, VM) in enumerate(VMs), (j, v) in enumerate(vs)
        e = xu_shen_equilibrium(XuShen(STD; VM), v)
        e === nothing || (Vst[i, j] = e.V; Nst[i, j] = e.N)
    end
    say("A1  boundary v_c(V_M)                 read-off          ours")
    for (VM, t) in ((0.05, "0.08–0.09"), (0.10, "0.115–0.12"), (0.15, "≈ 0.16"), (0.20, "≈ 0.19"))
        sayf("      V_M = %.2f                     %-14s    %.4f", VM, t, xu_shen_critical_rate(XuShen(STD; VM)).v)
    end
    e(VM, v) = xu_shen_equilibrium(XuShen(STD; VM), v)
    say("A2  V*: read-off < 0.6 at V_M ≈ 0.01, > 3.0 at V_M = 0.2 and low v; V* = 1.8 near V_M ≈ 0.075–0.08; nearly independent of v")
    sayf("      ours: V*(V_M = 0.01, v = 0.01) = %.3f; V*(0.2, 0.01) = %.3f; V_M with V* = 1.8 at v = 0.05: %.4f",
         e(0.01, 0.01).V, e(0.2, 0.01).V, crossing(VMs, [x === nothing ? NaN : x.V for x in e.(VMs, 0.05)], 1.8))
    sayf("      ours: V* at V_M = 0.1 for v = 0.01, 0.06, 0.12: %.3f, %.3f, %.3f (change with v: %.1f%%)", e(0.1, 0.01).V, e(0.1, 0.06).V,
         e(0.1, 0.12).V, 100 * (1 - e(0.1, 0.12).V / e(0.1, 0.01).V))
    say("A3  N*: read-off > 180 at v → 0.01 only for V_M ≲ 0.085, 160–180 on the rest of the bottom edge; < 60 just under the boundary at V_M ≥ 0.1")
    Nlow = [e(VM, 0.01) === nothing ? NaN : e(VM, 0.01).N for VM in VMs]
    sayf("      ours at v = 0.01: N* = %.1f, %.1f, %.1f, %.1f at V_M = 0.03, 0.05, 0.1, 0.2; N* = 180 at V_M = %.4f (falling); max %.1f at V_M = %.3f",
         e(0.03, 0.01).N, e(0.05, 0.01).N, e(0.1, 0.01).N, e(0.2, 0.01).N, crossing(VMs, Nlow, 180.0; from = argmax(Nlow)), maximum(Nlow), VMs[argmax(Nlow)])
    sayf("      ours at v = 0.99 v_c: N* = %.1f, %.1f, %.1f at V_M = 0.1, 0.15, 0.2 (at v_c itself: %.1f, %.1f, %.1f)",
         (e(VM, 0.99 * xu_shen_critical_rate(XuShen(STD; VM)).v).N for VM in (0.1, 0.15, 0.2))...,
         (xu_shen_critical_rate(XuShen(STD; VM)).N for VM in (0.1, 0.15, 0.2))...)
    fig = Figure(size = (1100, 450))
    for (k, (Z, ttl, lev)) in enumerate(((Vst, "(a) Genetic variance V*", 0.3:0.3:3.3), (Nst, "(b) Population size N*", 40:10:190)))
        ax = Axis(fig[1, 2k-1]; title = ttl, xlabel = "Variance input from mutation (V_M)", ylabel = "moving rate (v)")
        cf = contourf!(ax, VMs, vs, Z; levels = lev, colormap = :thermal)
        contour!(ax, VMs, vs, Z; levels = lev, color = (:black, 0.5), linewidth = 0.6)
        lines!(ax, VMs, vc; color = :black, linewidth = 2)
        limits!(ax, 0.005, 0.205, 0.0, 0.2)
        Colorbar(fig[1, 2k], cf)
    end
    Label(fig[0, :], "Xu & Shen (2026) Figure 1, recomputed from eqs. 6–8 (r = 0.2, ω² = 50, K = 1000, n_e = 1); black line: v_c"; fontsize = 13)
    save(joinpath(FIGS, "fig1.png"), fig)
end

# ------------------------------------------------------------------------------------------
function figure2()
    section("2. Figure 2: v_c, V* and N* against the strength of selection log₁₀(1/ω²)")
    xs = range(-3, 0, length = 301)
    sets = (("blue: r = 0.2, V_M = 0.1, K = 10³", STD, "0.087 at −3; peak 0.136 near −1.2; 0.093 at 0"),
            ("yellow: V_M = 0.5", XuShen(STD; VM = 0.5), "0.14 at −3; peak 0.19 near −1.5; 0.095 at 0"),
            ("green: K = 10⁴", XuShen(STD; K = 1e4), "0.157 at −3; peak 0.17 near −1.9; 0.097 at 0"),
            ("red: r = 0.5", XuShen(STD; r = 0.5), "0.19 at −3; peak 0.26 near −0.9; 0.24 at 0"))
    curves = [[xu_shen_critical_rate(XuShen(m; ω2 = om2(x))).v for x in xs] for (_, m, _) in sets]
    say("A4–A7  panel (a), v_c:")
    for ((name, _, target), c) in zip(sets, curves)
        i = argmax(c)
        sayf("  %-36s read-off: %s", name, target)
        sayf("  %-36s ours:     %.4f at −3; peak %.4f at %.2f; %.4f at 0", "", c[1], c[i], xs[i], c[end])
    end
    say("  Post hoc (after seeing that the yellow curve does not match V_M = 0.5): the same curve for other values of V_M, all else standard:")
    for VM in (0.15, 0.2, 0.25, 0.3)
        c = [xu_shen_critical_rate(XuShen(STD; VM, ω2 = om2(x))).v for x in xs]
        i = argmax(c)
        sayf("    V_M = %-5g %.4f at −3; peak %.4f at %.2f; %.4f at 0; at −1.7 (ω² = 50): %.4f", VM, c[1], c[i], xs[i], c[end], interp(xs, c, -1.7))
    end
    say("    (read-off of the yellow curve: 0.14 at −3; peak 0.19 near −1.5; 0.095 at 0)")
    function panels(m)
        vs = range(0.002, 0.14, length = 139)
        Vst = fill(NaN, length(xs), length(vs)); Nst = copy(Vst)
        for (i, x) in enumerate(xs), (j, v) in enumerate(vs)
            e = xu_shen_equilibrium(XuShen(m; ω2 = om2(x)), v)
            e === nothing || (Vst[i, j] = e.V; Nst[i, j] = e.N)
        end
        vc = [xu_shen_critical_rate(XuShen(m; ω2 = om2(x))).v for x in xs]
        (; vs, Vst, Nst, vc)
    end
    function describe(tag, m, P)
        i = argmax(P.vc)
        up = crossing(xs, P.vc, 0.12); down = crossing(xs, P.vc, 0.12; from = i)
        sayf("  %-22s boundary %.4f at −3, %.4f at 0; above 0.12 from %s to %s (peak %.4f at %.2f); max N* = %.1f; V* at v = 0.02: %.2f at −3, %.2f at −1.1, %.2f at −0.5",
             tag, P.vc[1], P.vc[end], isnan(up) ? "–" : @sprintf("%.2f", up), isnan(down) ? "–" : @sprintf("%.2f", down), P.vc[i], xs[i],
             maximum(filter(isfinite, P.Nst)), xu_shen_equilibrium(XuShen(m; ω2 = om2(-3.0)), 0.02).V,
             xu_shen_equilibrium(XuShen(m; ω2 = om2(-1.1)), 0.02).V, xu_shen_equilibrium(XuShen(m; ω2 = om2(-0.5)), 0.02).V)
    end
    say("\nA8  panels (b), (c). Read-off: boundary 0.045–0.05 at −3, above 0.12 between ≈ −1.2 and ≈ −0.8, 0.075–0.08 at 0; N* below about 100 everywhere;")
    say("    V* > 6 at the far left, < 0.5–1 right of ≈ −0.5.")
    P = panels(STD)
    describe("caption parameters:", STD, P)
    say("    With the caption's parameters the boundary is the blue line of panel (a) and N* reaches ", round(maximum(filter(isfinite, P.Nst)), digits = 1),
        ", as in Figure 1b; neither matches the read-off of panels (b), (c).")
    say("    Post hoc (after seeing that): one-parameter changes from the caption's set, looking for the read-off of (b), (c):")
    alts = (("K = 500", XuShen(STD; K = 500)), ("K = 300", XuShen(STD; K = 300)), ("r = 0.1", XuShen(STD; r = 0.1)), ("V_M = 0.05", XuShen(STD; VM = 0.05)),
            ("r = 0.1, K = 1000", XuShen(STD; r = 0.1)), ("r = 0.15", XuShen(STD; r = 0.15)), ("r = 0.1, V_M = 0.2", XuShen(STD; r = 0.1, VM = 0.2)))
    Ps = Dict{String,Any}()
    for (tag, m) in alts
        Ps[tag] = panels(m)
        describe("  " * tag * ":", m, Ps[tag])
    end
    say("    Post hoc scan: r ∈ {0.1, 0.15, …, 0.4}, K ∈ {100/r, 500, 1000}, V_M ∈ {0.05, 0.1, 0.15, 0.2}; distance to the read-off = the largest of")
    say("    |boundary(−3) − 0.0475|/0.0475, |boundary(0) − 0.0775|/0.0775, |up-crossing of 0.12 + 1.2|/0.2, |down-crossing + 0.8|/0.2, |max N* − 97|/97. The five closest:")
    scan = Tuple{Float64,String}[]
    for r in 0.1:0.05:0.4, K in unique((100 / r, 500.0, 1000.0)), VM in (0.05, 0.1, 0.15, 0.2)
        m = XuShen(STD; r, K, VM)
        vc = [xu_shen_critical_rate(XuShen(m; ω2 = om2(x))).v for x in xs]
        i = argmax(vc)
        up, down = crossing(xs, vc, 0.12), crossing(xs, vc, 0.12; from = i)
        Nmax = maximum(xu_shen_equilibrium(XuShen(m; ω2 = om2(x)), 0.002).N for x in (-3.0, -2.5, -2.0))
        d = max(abs(vc[1] - 0.0475) / 0.0475, abs(vc[end] - 0.0775) / 0.0775, isnan(up) ? 9.0 : abs(up + 1.2) / 0.2, isnan(down) ? 9.0 : abs(down + 0.8) / 0.2,
                abs(Nmax - 97) / 97)
        push!(scan, (d, @sprintf("r = %.2f, K = %-6.5g V_M = %.2f: boundary %.4f at −3, %.4f at 0; above 0.12 from %.2f to %.2f; max N* %.0f (distance %.2f)", r, K, VM, vc[1],
                                 vc[end], up, down, Nmax, d)))
    end
    for (_, line) in sort(scan; by = first)[1:5]
        say("      ", line)
    end
    fig = Figure(size = (1500, 820))
    ax = Axis(fig[1, 1]; title = "(a) Critical moving rate", xlabel = "Selection strength (log₁₀(1/ω²))", ylabel = "v_c")
    for ((name, _, _), c, col) in zip(sets, curves, (:royalblue, :goldenrod, :olivedrab, :orangered))
        lines!(ax, xs, c; label = name, color = col, linewidth = 2.5)
    end
    limits!(ax, -3, 0, 0, 0.3); axislegend(ax; position = :lb, labelsize = 10)
    function draw(row, P, lab)
        for (k, (Z, ttl, lev)) in enumerate(((P.Vst, "(b) V*, $lab", 0:0.5:10), (P.Nst, "(c) N*, $lab", 0:10:200)))
            a = Axis(fig[row, 2k]; title = ttl, xlabel = "Selection strength (log₁₀(1/ω²))", ylabel = "Moving rate (v)")
            cf = contourf!(a, xs, P.vs, Z; levels = lev, colormap = :thermal)
            lines!(a, xs, P.vc; color = :black, linewidth = 2)
            hlines!(a, [0.12]; color = :gray, linestyle = :dash)
            limits!(a, -3, 0, 0, 0.14)
            Colorbar(fig[row, 2k+1], cf)
        end
    end
    draw(1, P, "caption parameters")
    draw(2, Ps["K = 500"], "K = 500 (post hoc)")
    Label(fig[0, :], "Xu & Shen (2026) Figure 2, recomputed from eqs. 6–8; dashed line: the top of the preprint's frame (v = 0.12)"; fontsize = 13)
    save(joinpath(FIGS, "fig2.png"), fig)
end

# ------------------------------------------------------------------------------------------
const PANELS3 = (("(a)", 0.08, 0.5), ("(b)", 0.08, 3.0), ("(c)", 0.08, 5.0), ("(d)", 0.12, 0.5), ("(e)", 0.12, 3.0), ("(f)", 0.12, 5.0))
const V0GRID = (0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
const N0GRID = (10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 125.0, 150.0)

critN(m, V0; Δ0, v, Next = 1e-4, Nmin = 0.01) = xu_shen_critical_size(m; Δ0, v, V0, Nmin, Nmax = 1e5, Next)

function figure3()
    section("3. Figure 3: basins of the viable equilibrium and of extinction in (V₀, N₀) (standard parameters)")
    m = STD
    say("Extinction rule: N < Next, Next = 1e-3 by default (1e-4 where the search goes down to N₀ = 0.01); sensitivity below.")
    say("A9  viable equilibrium (red dot)      read-off (V*, N*)     ours (V*, N*, Δ*)           saddle (V, N, Δ)")
    for (v, t) in ((0.08, "(2.06, 142)"), (0.12, "(1.96, 88)"))
        s, e = xu_shen_equilibria(m, v)
        sayf("      v = %.2f                        %-18s    (%.3f, %.1f, %.3f)      (%.3f, %.2f, %.3f)", v, t, e.V, e.N, e.Δ, s.V, s.N, s.Δ)
    end
    targets = Dict("(a)" => "all plotted starts persist; boundary along the axes",
                   "(b)" => "all plotted starts persist; boundary roughly horizontal at N₀ ≈ 2–5",
                   "(c)" => "(0.3, 180), (0.5, 75), (1.0, 38), (1.5, 25), (2.0, 17), (3.0, 11)",
                   "(d)" => "all plotted starts persist; boundary along the axes",
                   "(e)" => "near-vertical at V₀ ≈ 0.02–0.05 above N₀ ≈ 50; (0.5, 15), (1.0, 10), (2.0, 7), (3.0, 6)",
                   "(f)" => "(1.7, 180), (2.0, 100), (2.5, 55), (3.0, 35)")
    fig = Figure(size = (1350, 900))
    Vb = 10 .^ range(log10(0.01), log10(3.2), length = 90)
    say("\nA10–A14  basin boundary: critical N₀ at given V₀ (bisection), and critical V₀ at given N₀")
    for (k, (lab, v, Δ0)) in enumerate(PANELS3)
        e = xu_shen_equilibrium(m, v)
        ax = Axis(fig[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1]; title = "$lab v = $v, Δ₀ = $Δ0", xlabel = "Genetic variance (V)", ylabel = "Population size (N)")
        npers = 0; nover = 0; ndip = 0; nrise = 0
        for V0 in V0GRID, N0 in N0GRID
            tr = xu_shen_trajectory(m, Δ0, V0, N0; v, maxstep = 0.25)
            ok = tr.fate === :persists
            npers += ok
            ok && N0 < e.N && maximum(tr.N) > 1.02 * e.N && (nover += 1)
            ok && minimum(tr.N) < 0.9 * N0 && minimum(tr.V) < 0.9 * V0 && (ndip += 1)
            !ok && maximum(tr.N) > 1.02 * N0 && maximum(tr.V) > 1.02 * V0 && (nrise += 1)
            lines!(ax, tr.V, tr.N; color = ok ? (:orange, 0.7) : (:gray40, 0.7), linewidth = 0.8)
            scatter!(ax, [V0], [N0]; color = :transparent, strokecolor = ok ? :goldenrod : :gray40, strokewidth = 1, markersize = 6)
        end
        Nb = [critN(m, V0; Δ0, v) for V0 in Vb]
        lines!(ax, Vb, [isfinite(x) ? x : 1e5 for x in Nb]; color = :blue, linewidth = 2.5)
        scatter!(ax, [e.V], [e.N]; color = :red, markersize = 12)
        scatter!(ax, [0.0], [0.0]; color = :black, markersize = 10)
        limits!(ax, -0.1, 3.2, -8, 190)
        sayf("  %s v = %.2f, Δ₀ = %.1f   read-off: %s", lab, v, Δ0, targets[lab])
        sayf("      ours: %d of %d plotted starts persist; critical N₀ at V₀ = 0.2, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0: %s", npers, length(V0GRID) * length(N0GRID),
             join(map(x -> x <= 0.01 ? "<0.01" : isfinite(x) ? @sprintf("%.1f", x) : "none", [critN(m, V0; Δ0, v) for V0 in (0.2, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)]), ", "))
        cv = [xu_shen_boundary(m, N0; Δ0, v) for N0 in (180.0, 100.0, 50.0, 20.0, 10.0)]
        sayf("      ours: critical V₀ at N₀ = 180, 100, 50, 20, 10: %s", join(map(x -> x < 1e-5 ? "≈0" : isfinite(x) ? @sprintf("%.3f", x) : ">8.9", cv), ", "))
        sayf("      A15: persisting starts below N* that overshoot it by > 2%%: %d; persisting starts whose N and V both first fall > 10%%: %d; extinct starts whose N and V both first rise > 2%%: %d",
             nover, ndip, nrise)
    end
    Label(fig[0, :], "Xu & Shen (2026) Figure 3, recomputed from eqs. 2, 4, 5 (standard parameters). Orange: persists; gray: extinct; blue: basin boundary; red: viable equilibrium";
          fontsize = 13)
    save(joinpath(FIGS, "fig3.png"), fig)

    say("\nIs persistence monotone in V₀ and in N₀? (30 × 30 grid, V₀ ∈ [0.05, 3], N₀ ∈ [1, 180]; the bisections assume it)")
    Vg = range(0.05, 3, length = 30); Ng = 10 .^ range(0, log10(180), length = 30)
    for (lab, v, Δ0) in PANELS3
        ok = [xu_shen_persists(m, Δ0, V0, N0; v) for V0 in Vg, N0 in Ng]
        badV = count(ok[i, j] && !ok[i+1, j] for i in 1:29, j in 1:30)
        badN = count(ok[i, j] && !ok[i, j+1] for i in 1:30, j in 1:29)
        sayf("  %s persisting %3d of 900; violations: %d in V₀, %d in N₀", lab, count(ok), badV, badN)
    end

    say("\nSensitivity to the extinction threshold Next: critical N₀ at V₀ = 0.5, 1, 2, 3 (searched down to max(Next, 0.01))")
    for (lab, v, Δ0) in PANELS3[[2, 3, 5, 6]]
        for Next in (1.0, 1e-2, 1e-4, 1e-8)
            c = [critN(m, V0; Δ0, v, Next, Nmin = max(1.0001 * Next, 0.01)) for V0 in (0.5, 1.0, 2.0, 3.0)]
            sayf("  %s Next = %-6g  %s", lab, Next, join(map(x -> isfinite(x) ? @sprintf("%9.4f", x) : "     none", c), " "))
        end
    end
    say("  Panels (a), (d) (Δ₀ = 0.5): critical N₀ at V₀ = 1 with Next = 1e-8, searched down to 1e-6:")
    for (lab, v, Δ0) in PANELS3[[1, 4]]
        sayf("  %s %.3g", lab, critN(m, 1.0; Δ0, v, Next = 1e-8, Nmin = 1e-6))
    end
    say("  Integration tolerance: critical N₀ in panels (c), (f) with rtol = 1e-5, 1e-7 (default), 1e-9:")
    for (lab, v, Δ0, V0) in (("(c), V₀ = 1.0", 0.08, 5.0, 1.0), ("(f), V₀ = 2.5", 0.12, 5.0, 2.5))
        sayf("  %s: %s", lab, join((@sprintf("%.6f", xu_shen_critical_size(m; Δ0, v, V0, Nmin = 0.01, Nmax = 1e5, Next = 1e-4, rtol)) for rtol in (1e-5, 1e-7, 1e-9)), ", "))
    end
end

# ------------------------------------------------------------------------------------------
fmtN(x) = isnan(x) ? "  (v ≥ v_c)" : isinf(x) ? "      > 1e6" : @sprintf("%11.4g", x)

function figure4()
    section("4. Figure 4: the three scenarios")
    fig = Figure(size = (1400, 1050))

    # ---- (a) V₀ at mutation–selection–drift balance -----------------------------------------
    say("A16  panel (a): critical N₀ with V₀ at balance for that N₀ (floor 1 as on the preprint's axis; ours searched from 1 to 1e6)")
    xs = collect(-3:0.05:0)
    curvesA = (("Δ₀ = 1, v = 0.05", 1.0, 0.05, :black, :solid), ("Δ₀ = 1, v = 0.1", 1.0, 0.1, :red, :solid),
               ("Δ₀ = 3, v = 0.05", 3.0, 0.05, :black, :dash), ("Δ₀ = 3, v = 0.1 (legend: Δ₀ = 1)", 3.0, 0.1, :red, :dash))
    function curveA(Δ0, v; Next = QTLSim.XU_SHEN_NEXT, Nmin = 1.0)
        [(m = XuShen(STD; ω2 = om2(x)); xu_shen_critical_rate(m).v <= v ? NaN : xu_shen_critical_size(m; Δ0, v, Nmin, Nmax = 1e6, Next)) for x in xs]
    end
    ax = Axis(fig[1, 1]; title = "(a) V₀ at mutation–selection–drift balance", xlabel = "Selection strength (log₁₀(1/ω²))", ylabel = "Critical initial population size",
              yscale = log10)
    resA = Dict{String,Vector{Float64}}()
    for (name, Δ0, v, col, sty) in curvesA
        c = curveA(Δ0, v)
        resA[name] = c
        lines!(ax, xs, [isinf(y) ? 1e6 : y for y in c]; color = col, linestyle = sty, linewidth = 2, label = name)
    end
    limits!(ax, -3, 0, 0.8, 1.2e4); axislegend(ax; position = :lt, labelsize = 10)
    at(name, x) = fmtN(interp(xs, resA[name], x))
    leaves(name) = (c = resA[name]; i = findfirst(y -> isfinite(y) && y > 1.05, c); i === nothing ? NaN : xs[i])
    drawn(name) = (c = resA[name]; i = findall(!isnan, c); isempty(i) ? "–" : @sprintf("%.2f to %.2f", xs[first(i)], xs[last(i)]))
    say("  black solid (Δ₀ = 1, v = 0.05)  read-off: 1 up to ≈ −0.65; ≈ 3 at 0")
    sayf("                                  ours:     first above 1.05 at %.2f; %s at 0", leaves(curvesA[1][1]), at(curvesA[1][1], 0.0))
    say("  red solid (Δ₀ = 1, v = 0.1)     read-off: drawn from ≈ −2.6 to ≈ −0.1; leaves 1 at ≈ −1.35; ≈ 2 at −1.1; ≈ 5 at −0.6; ≈ 18 at −0.1")
    sayf("                                  ours:     defined (v < v_c) from %s; first above 1.05 at %.2f; %s at −1.1; %s at −0.6; %s at −0.1",
         drawn(curvesA[2][1]), leaves(curvesA[2][1]), at(curvesA[2][1], -1.1), at(curvesA[2][1], -0.6), at(curvesA[2][1], -0.1))
    say("  black dashed (Δ₀ = 3, v = 0.05) read-off: leaves 1 at ≈ −1.65; ≈ 13 at −1.0; ≈ 40 at −0.7; ≈ 1e3 at −0.15; 1e4 at 0")
    sayf("                                  ours:     first above 1.05 at %.2f; %s at −1.0; %s at −0.7; %s at −0.15; %s at 0", leaves(curvesA[3][1]),
         at(curvesA[3][1], -1.0), at(curvesA[3][1], -0.7), at(curvesA[3][1], -0.15), at(curvesA[3][1], 0.0))
    say("  red dashed (taken as Δ₀ = 3, v = 0.1) read-off: leaves 1 at ≈ −2.1; ≈ 8 at −1.6; ≈ 30 at −1.1; ≈ 1e2 at −0.8; ≈ 1e3 at −0.45; off the top (1e4) at ≈ −0.35")
    c4 = resA[curvesA[4][1]]
    sayf("                                  ours:     first above 1.05 at %.2f; %s at −1.6; %s at −1.1; %s at −0.8; %s at −0.45; reaches 1e4 at %.2f; no persistence from any N₀ ≤ 1e6 from %s",
         leaves(curvesA[4][1]), at(curvesA[4][1], -1.6), at(curvesA[4][1], -1.1), at(curvesA[4][1], -0.8), at(curvesA[4][1], -0.45),
         crossing(xs, [isinf(y) ? 1e7 : y for y in c4], 1e4), (i = findfirst(isinf, c4); i === nothing ? "–" : @sprintf("%.2f", xs[i])))
    say("  The legend's fourth entry read as printed (Δ₀ = 1, v = 0.1) would repeat the red solid curve; with Δ₀ = 3 it is the curve above.")
    say("  Text, line 160 and caption: the dashed curves go to infinity (persistence impossible however large N₀, with V₀ at balance) although v < v_c.")
    say("  Ours, Δ₀ = 3: does a population persist from N₀ = 1e4, 1e6, 1e9, 1e12 (V₀ → V_c^M as N₀ → ∞)?")
    for (v, x) in ((0.05, -0.5), (0.05, 0.0), (0.1, -0.5), (0.1, -0.3), (0.1, -0.2), (0.1, -0.15))
        m = XuShen(STD; ω2 = om2(x))
        vcm = xu_shen_critical_rate(m).v
        v < vcm || (sayf("    v = %.2f, x = %.2f: v ≥ v_c = %.4f", v, x, vcm); continue)
        sayf("    v = %.2f, x = %5.2f (v_c = %.4f): %s, %s, %s, %s; critical V₀ at N₀ = 1e12: %.4f against V_c^M = %.4f", v, x, vcm,
             (xu_shen_persists(m, 3.0, xu_shen_balance_variance(m, N0), N0; v) for N0 in (1e4, 1e6, 1e9, 1e12))..., xu_shen_boundary(m, 1e12; Δ0 = 3.0, v), xu_shen_vmax(m))
    end
    say("  Sensitivity of panel (a) to the extinction rule: Next = 1 (fewer than one individual) against the default 1e-3, at x = −2, −1.5, −1, −0.5, 0:")
    for (name, Δ0, v, _, _) in curvesA
        a = [(m = XuShen(STD; ω2 = om2(x)); xu_shen_critical_rate(m).v <= v ? NaN : xu_shen_critical_size(m; Δ0, v, Nmin = 1.0001, Nmax = 1e6, Next = 1.0)) for x in (-2.0, -1.5, -1.0, -0.5, 0.0)]
        b = [interp(xs, resA[name], x) for x in (-2.0, -1.5, -1.0, -0.5, 0.0)]
        sayf("    %-34s Next = 1:    %s", name, join(fmtN.(a), " "))
        sayf("    %-34s Next = 1e-3: %s", "", join(fmtN.(b), " "))
    end
    say("  Below the floor: the critical N₀ of the equations themselves (Next = 1e-9, searched down to 1e-6), at x = −3, −2, −1:")
    for (name, Δ0, v, _, _) in curvesA
        a = [(m = XuShen(STD; ω2 = om2(x)); xu_shen_critical_rate(m).v <= v ? NaN : xu_shen_critical_size(m; Δ0, v, Nmin = 1e-6, Nmax = 1e6, Next = 1e-9)) for x in (-3.0, -2.0, -1.0)]
        sayf("    %-34s %s", name, join(fmtN.(a), " "))
    end

    # ---- (b) a change of rate with an abrupt shift --------------------------------------------
    say("\nA17  panel (b): the largest new rate v₁ after an abrupt shift δ, from the equilibrium at v₀")
    δs = collect(0:0.1:10.5)
    curvesB = ((20.0, 0.04, :black, :solid, "2.6", "5.65"), (20.0, 0.08, :red, :solid, "1.9", "4.95"), (100.0, 0.04, :black, :dash, "4.8", "9.2"),
               (100.0, 0.08, :red, :dash, "3.3", "7.7"))
    ax = Axis(fig[1, 2]; title = "(b) shift in environmental change rate", xlabel = "Magnitude of abrupt shift (δ)", ylabel = "Max v₁ for persistence")
    sayf("  plateau v_c: read-off ≈ 0.1355 (ω² = 20), ≈ 0.123 (ω² = 100); ours %.4f, %.4f", xu_shen_critical_rate(XuShen(STD; ω2 = 20)).v,
         xu_shen_critical_rate(XuShen(STD; ω2 = 100)).v)
    for (ω2, v0, col, sty, tl, tz) in curvesB
        m = XuShen(STD; ω2)
        vc = xu_shen_critical_rate(m).v
        c = [xu_shen_max_rate(m; v0, δ, steps = 30) for δ in δs]
        lines!(ax, δs, c; color = col, linestyle = sty, linewidth = 2, label = "ω² = $(Int(ω2)), v₀ = $v0")
        ileave = findfirst(<(vc * (1 - 1e-9)), c); izero = findfirst(==(0.0), c)
        # refine the two points by bisection in δ
        function refine(pred, lo, hi)
            for _ in 1:30
                mid = (lo + hi) / 2
                pred(mid) ? (hi = mid) : (lo = mid)
            end
            (lo + hi) / 2
        end
        δl = ileave === nothing ? NaN : refine(δ -> xu_shen_max_rate(m; v0, δ, steps = 2) < vc, δs[max(ileave - 1, 1)], δs[ileave])
        δz = izero === nothing ? NaN : refine(δ -> xu_shen_max_rate(m; v0, δ, steps = 2) == 0, δs[max(izero - 1, 1)], δs[izero])
        sayf("  ω² = %3d, v₀ = %.2f   read-off: leaves the plateau at δ ≈ %s, reaches 0 at ≈ %s;   ours: %.2f (1%% below the plateau at %.2f, 3%% at %.2f), %.2f;   ours at δ = 2, 4, 6, 8: %.4f, %.4f, %.4f, %.4f",
             ω2, v0, tl, tz, δl, crossing(δs, c, 0.99vc), crossing(δs, c, 0.97vc), δz, (interp(δs, c, x) for x in (2.0, 4.0, 6.0, 8.0))...)
    end
    limits!(ax, 0, 10.5, 0, 0.15); axislegend(ax; position = :rt, labelsize = 10)
    say("  Text, line 174: with δ = 0 a population at the equilibrium of any v₀ persists at any v₁ < v_c. Ours, standard parameters, v₀ and v₁ on a grid of 14 × 14")
    vc = xu_shen_critical_rate(STD).v
    grid = range(0, 0.999vc, length = 14)
    np = count(xu_shen_persists(STD, e.Δ, e.V, e.N; v = v1) for v0 in grid for e in (xu_shen_equilibrium(STD, v0),) for v1 in grid)
    sayf("  values in [0, 0.999 v_c]: %d of %d persist.", np, length(grid)^2)

    # ---- (c) the boundary surface and the curve of equilibria --------------------------------
    say("\nA18  panel (c): boundary V₀ over (Δ₀, N₀). The preprint does not state the rate of the surface; ours is drawn for v = 0.08 (and v = 0.12 in the log).")
    Δg = range(0, 7, length = 57); Ng = 10 .^ range(0, log10(300), length = 50)
    surf(v) = [xu_shen_boundary(STD, N0; Δ0, v, Vmax = 3.0) for Δ0 in Δg, N0 in Ng]
    S = surf(0.08)
    gc = fig[2, 1] = GridLayout()
    ax = Axis(gc[1, 1]; title = "(c) Basin boundary V₀ (v = 0.08) and viable equilibria for v from 0 to v_c", xlabel = "Phenotypic lag (Δ₀)", ylabel = "Population size (N₀)",
              yscale = log10)
    hm = heatmap!(ax, Δg, Ng, [isinf(x) ? 3.0 : x for x in S]; colormap = :blues, colorrange = (0, 3))
    contour!(ax, Δg, Ng, [isinf(x) ? 3.5 : x for x in S]; levels = [0.5, 1.0, 1.5, 2.0, 2.5, 2.99], color = :black, linewidth = 0.7)
    vv = range(0, vc, length = 80)
    eqs = [v < vc ? xu_shen_equilibrium(STD, v) : nothing for v in vv]
    eqs = [e === nothing ? (c = xu_shen_critical_rate(STD); (; c.Δ, c.V, c.N)) : e for e in eqs]
    lines!(ax, [e.Δ for e in eqs], [e.N for e in eqs]; color = :red, linewidth = 2.5)
    scatter!(ax, [eqs[end].Δ], [eqs[end].N]; color = :white, strokecolor = :red, strokewidth = 2, markersize = 10)
    Colorbar(gc[1, 2], hm; label = "Boundary V₀ (darkest: > 3)")
    for v in (0.08, 0.12)
        Sv = v == 0.08 ? S : surf(v)
        j50 = argmin(abs.(Ng .- 50)); j10 = argmin(abs.(Ng .- 10)); j200 = argmin(abs.(Ng .- 200))
        f(j) = join(map(x -> isinf(x) ? ">3" : x < 1e-4 ? "≈0" : @sprintf("%.2f", x), [interp(Δg, [isinf(y) ? 99.0 : y for y in Sv[:, j]], d) > 3 ? Inf : interp(Δg, Sv[:, j], d) for d in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)]), ", ")
        sayf("  v = %.2f: boundary V₀ at Δ₀ = 1, 2, 3, 4, 5, 6:  N₀ = 10: %s;  N₀ = 50: %s;  N₀ = 200: %s", v, f(j10), f(j50), f(j200))
    end
    say("  Qualitative target: ≈ 0 at small lags, rising steeply beyond a lag of a few units, more steeply at small N; the equilibria lie on the persisting side.")
    below = count(1:length(vv)-1) do i
        e = eqs[i]
        xu_shen_boundary(STD, e.N; Δ0 = e.Δ, v = 0.08, Vmax = 3.0) < e.V
    end
    sayf("  ours: of the equilibria for 79 rates v′ in [0, v_c), %d lie above the v = 0.08 boundary surface (V*(v′) > boundary V₀ at (Δ*(v′), N*(v′))).", below)
    sayf("  ours: at v → v_c the equilibrium (Δ, V, N) = (%.3f, %.3f, %.1f) is the point where saddle and node merge, on the boundary of its own rate.", eqs[end].Δ,
         eqs[end].V, eqs[end].N)

    # ---- (d) a temporary shift ----------------------------------------------------------------
    say("\nA19  panel (d): the longest duration T at rate v₁ after which the population still persists back at v₀")
    v1s = collect(0.095:0.0025:0.25)
    ax = Axis(fig[2, 2]; title = "(d) Temporary environmental shift", xlabel = "New movement rate (v₁)", ylabel = "Max duration for persistence")
    resD = Dict{Tuple{Float64,Float64},Vector{Float64}}()
    for (ω2, col) in ((100.0, :black), (10.0, :red), (1.0, :blue)), (v0, sty) in ((0.02, :solid), (0.09, :dash))
        m = XuShen(STD; ω2)
        vcm = xu_shen_critical_rate(m).v
        if v0 >= vcm
            sayf("  ω² = %g, v₀ = %.2f: no equilibrium at v₀ (v_c = %.4f); curve not computed", ω2, v0, vcm)
            continue
        end
        c = [v1 <= vcm ? Inf : xu_shen_max_duration(m; v0, v1, Tmax = 5e4, steps = 34) for v1 in v1s]
        resD[(ω2, v0)] = c
        lines!(ax, v1s, [isinf(y) ? NaN : y for y in c]; color = col, linestyle = sty, linewidth = 2.5, label = "ω² = $(Int(ω2)), v₀ = $v0")
    end
    limits!(ax, 0.09, 0.255, 0, 800); axislegend(ax; position = :rt, labelsize = 10)
    sayf("  divergence (v_c): read-off ≈ 0.095–0.10 (ω² = 1), ≈ 0.123 (ω² = 100), ≈ 0.135–0.14 (ω² = 10); ours %.4f, %.4f, %.4f",
         (xu_shen_critical_rate(XuShen(STD; ω2)).v for ω2 in (1.0, 100.0, 10.0))...)
    T(ω2, v0, v1) = haskey(resD, (ω2, v0)) ? interp(v1s, resD[(ω2, v0)], v1) : NaN
    say("  solid (v₀ = 0.02)        read-off                                     ours")
    sayf("    ω² = 100: T at v₁ = 0.15, 0.20, 0.25:  ≈ 230, ≈ 100, ≈ 60            %.0f, %.0f, %.0f", T(100.0, 0.02, 0.15), T(100.0, 0.02, 0.2), T(100.0, 0.02, 0.25))
    sayf("    ω² = 10:  T at v₁ = 0.148, 0.20, 0.25: ≈ 220, ≈ 60, ≈ 30             %.0f, %.0f, %.0f", T(10.0, 0.02, 0.148), T(10.0, 0.02, 0.2), T(10.0, 0.02, 0.25))
    sayf("    ω² = 1:   T at v₁ = 0.11, 0.15, 0.20:  ≈ 190, ≈ 60, ≈ 20             %.0f, %.0f, %.0f", T(1.0, 0.02, 0.11), T(1.0, 0.02, 0.15), T(1.0, 0.02, 0.2))
    d = [T(100.0, 0.02, x) - T(10.0, 0.02, x) for x in v1s]
    sayf("    the ω² = 10 and ω² = 100 solid curves cross at v₁ = %.4f (read-off ≈ 0.148)", crossing(v1s, d, 0.0))
    sayf("  dashed (v₀ = 0.09): ω² = 100 at v₁ = 0.15: read-off ≈ 180; ours %.0f. Dashed below solid at every v₁ > v_c: ω² = 100: %s; ω² = 10: %s; ω² = 1: %s",
         T(100.0, 0.09, 0.15), (haskey(resD, (w, 0.09)) ? string(all(a >= b for (a, b) in zip(resD[(w, 0.02)], resD[(w, 0.09)]) if isfinite(a))) : "–" for w in (100.0, 10.0, 1.0))...)
    say("\n  Re-reading after the comparison (NOT blind; the read-offs of A19 above were taken from a crop without the x-axis ticks and misplaced the v₁ axis).")
    say("  With the ticks in view (0.10, 0.15, 0.20, 0.25 at 220 px per 0.05; T = 0, 200, …, 800 at 0.685 px per unit):")
    sayf("    ω² = 100 solid: T ≈ 170, 89, 55 at v₁ = 0.15, 0.20, 0.25 (ours %.0f, %.0f, %.0f); dashed ≈ 125 at 0.15 (ours %.0f)", T(100.0, 0.02, 0.15), T(100.0, 0.02, 0.2),
         T(100.0, 0.02, 0.25), T(100.0, 0.09, 0.15))
    sayf("    ω² = 10 solid: crosses the ω² = 100 curve at v₁ ≈ 0.143, T ≈ 220 (ours %.4f, T = %.0f); ≈ 50, 30 at 0.20, 0.25 (ours %.0f, %.0f)", crossing(v1s, d, 0.0),
         T(100.0, 0.02, crossing(v1s, d, 0.0)), T(10.0, 0.02, 0.2), T(10.0, 0.02, 0.25))
    sayf("    ω² = 1 solid: T ≈ 185, 105, 40, 22 at v₁ = 0.10, 0.11, 0.15, 0.20 (ours %.0f, %.0f, %.0f, %.0f)", T(1.0, 0.02, 0.10), T(1.0, 0.02, 0.11), T(1.0, 0.02, 0.15),
         T(1.0, 0.02, 0.2))
    say("  Likewise for panel (b), read again with the ticks in view (not blind): max v₁ ≈ 0.095 (ω² = 20, v₀ = 0.04) and ≈ 0.0565 (20, 0.08) at δ = 4; ≈ 0.107 (100, 0.04) and")
    say("  ≈ 0.0605 (100, 0.08) at δ = 6; ≈ 0.0446 (100, 0.04) at δ = 8: compare \"ours at δ = 2, 4, 6, 8\" above.")
    say("  And Figure 1a's boundary at V_M = 0.10, read again along the concave edge (not blind): ≈ 0.13 (the first read-off, 0.115–0.12, was too low).")
    Label(fig[0, :], "Xu & Shen (2026) Figure 4, recomputed from eqs. 2, 4, 5 (panel c as a map: the preprint draws a surface whose rate is not stated)"; fontsize = 13)
    save(joinpath(FIGS, "fig4.png"), fig)
end

# ------------------------------------------------------------------------------------------
# 5. eq. 6 against drift_critical_rate. The mapping is validate/xu_shen_targets.md, part B.

const C0 = QTLSim.QTSIM_FITNESS_PEAK

"The Xu–Shen parameters of a QTLSim `Params` under Beverton–Holt (targets file, part B)."
function xs_map(p::Params; Kmap::Symbol = :K1)
    ω, VE, R = initFitSD(p), QTLSim.envVar(p), p.GrowthRate
    XuShen(r = log(R) + log(C0) - 0.5 * log(1 + VE / ω^2), ω2 = ω^2 + VE, VM = QTLSim.muteVar(p),
           K = p.Ploidy * p.PopSize / (Kmap === :K1 ? log(R) : R - 1), ne = QTLSim.gaussian_units(p))
end

"""
The critical rate of the steady-state feedback with each ingredient switchable (a local copy of
`QTLSim._feedback_state` and `drift_critical_rate`'s maximisation; src/critical.jl is not changed):
`variance` `:hog` (house of Gauss) or `:gauss` (the α²/(2V_s) term dropped: eq. 5); `drift` `:breeders`
(N_e = effective_size(N, n)), `:census` (the adults N) or `:n` (the breeders n themselves);
`form` `:discrete` or `:continuous` (first-order response and loads); `density` `:bh` or `:logistic`.
Unlike `drift_critical_rate`, which follows one branch (its inner iteration starts from the variance at
N = K), this takes the maximum over every steady state, parametrised by its variance as eq. 6 is.
"""
function ladder_kc(p::Params; variance::Symbol, drift::Symbol, form::Symbol, density::Symbol, Kmap::Symbol = :K1)
    ω, VE, K, R = initFitSD(p), QTLSim.envVar(p), float(p.PopSize), p.GrowthRate
    Vs = ω^2 + VE
    nu, μ = QTLSim.rare_units(p)
    α2 = QTLSim.mutation_effect_variance(p)
    r = log(R) + log(C0) - 0.5 * log(1 + VE / ω^2)
    vg(ng) = ng <= 0 ? 0.0 : (b = (variance === :hog ? α2 / (2Vs) : 0.0) + 1 / ng; nu * 2μ * α2 / (b + sqrt(b^2 + 4μ * α2 / Vs)))
    # the gene copies that drift acts on in a population with variance V at lag λ (0: no population)
    function genes(V, λ)
        ρ, W = if form === :discrete
            s = ω^2 + V + VE
            w = C0 * ω / sqrt(s) * exp(-λ^2 / (2s))
            log(R * w), w
        else
            x = r - (V + λ^2) / (2Vs)
            x, exp(x) / R
        end
        ρ <= 0 && return 0.0
        N = density === :bh ? K * (R * W - 1) / (R - 1) / W : K * ρ / (Kmap === :K1 ? log(R) : R - 1)
        n = N * W
        n < p.Ploidy && return 0.0
        p.Ploidy * (drift === :breeders ? effective_size(N, n) : drift === :census ? N : n)
    end
    # every steady state, by its variance: the lag at which the population holds exactly V
    # (vg(genes(V, λ)) falls with λ), and the response there
    function rate(V)
        vg(genes(V, 0.0)) > V || return 0.0
        lo, hi = 0.0, sqrt(Vs)
        while vg(genes(V, hi)) > V
            hi *= 2
        end
        for _ in 1:100
            mid = (lo + hi) / 2
            vg(genes(V, mid)) > V ? (lo = mid) : (hi = mid)
        end
        λ = (lo + hi) / 2
        form === :discrete ? V * λ / (V + Vs) : V * λ / Vs
    end
    Vg = range(0, vg(Inf), length = 4002)[2:end-1]
    ks = rate.(Vg)
    i = argmax(ks)
    ks[i] > 0 || return 0.0
    a, b = Vg[max(i - 1, 1)], Vg[min(i + 1, end)]
    for _ in 1:80
        m1, m2 = b - 0.618034 * (b - a), a + 0.618034 * (b - a)
        rate(m1) < rate(m2) ? (a = m1) : (b = m2)
    end
    max(rate((a + b) / 2), ks[i])
end

rungs(p; Kmap = :K1) = (ladder_kc(p; variance = :hog, drift = :breeders, form = :discrete, density = :bh),
                        ladder_kc(p; variance = :gauss, drift = :breeders, form = :discrete, density = :bh),
                        ladder_kc(p; variance = :gauss, drift = :census, form = :discrete, density = :bh),
                        ladder_kc(p; variance = :gauss, drift = :census, form = :continuous, density = :bh),
                        ladder_kc(p; variance = :gauss, drift = :census, form = :continuous, density = :logistic, Kmap))

function ladder_table(cases)
    say("  rung 0 = drift_critical_rate's steady states (house of Gauss, N_e from breeders, discrete, Beverton–Holt), maximised over all of them (\"rung 0\" column;")
    say("  equal to k_c unless marked *); 1: variance equation → eq. 5 (no α²/(2V_s) term);")
    say("  2: drift from the census of adults; 3: first-order response and loads; 4: logistic regulation (K1) = eq. 6 at the mapped parameters.")
    say("  Each ratio is rung i / rung i−1; their product is v_c(eq. 6)/rung 0. The attribution depends on the order of the rungs. \"ratio\" = v_c(eq. 6)/rung 0.")
    sayf("  %-44s %10s %10s %10s %7s | %7s %7s %7s %7s | %9s %9s %7s", "case", "k_c", "rung 0", "v_c(eq.6)", "ratio", "1/0", "2/1", "3/2", "4/3", "rung4 chk", "v_c (K2)", "K2/rung0")
    for (name, p) in cases
        kc = drift_critical_rate(p).k
        vc = xu_shen_critical_rate(xs_map(p)).v
        vc2 = xu_shen_critical_rate(xs_map(p; Kmap = :K2)).v
        g = rungs(p)
        sayf("  %-44s %10.4g %10.4g%s %9.4g %7.3f | %7.3f %7.3f %7.3f %7.3f | %9.4g %9.4g %7.3f", name, kc, g[1], abs(g[1] - kc) <= 1e-3 * g[1] ? " " : "*", vc, vc / g[1],
             g[2] / g[1], g[3] / g[2], g[4] / g[3], g[5] / g[4], g[5], vc2, vc2 / g[1])
    end
end

"A QTLSim `Params` whose mapped Xu–Shen parameters are (r, ω2, VM, K, ne): haploid clone for ne = 1, diploid with ne/2 loci otherwise."
function qtl_from(; r, ω2, VM, K, ne, μ = 0.1, L = 100)
    VE = ω2 / 50                                   # a small environmental variance (QTSim scales mutations by V_E)
    ω = sqrt(ω2 - VE)
    R = exp(r - log(C0) + 0.5 * log(1 + VE / ω^2))
    ploidy = ne == 1 ? 1 : 2
    loci = ne == 1 ? L : Int(ne ÷ 2)
    Params(PopSize = round(Int, K * log(R) / ploidy), Ploidy = ploidy, NoQTL = loci, RecombPoissonMean = 3.0, PhenotypeVar = 2VE, Heritability = 0.5,
           MutationScalar = VM / VE, MutationRatePerGamete = μ, InitFitnessSD = ω, MoveType = :continuous, Demography = :beverton_holt, GrowthRate = R,
           CrossoverModel = :all)
end

function comparison()
    section("5. Eq. 6 against QTLSim's discrete-generation counterpart, drift_critical_rate (descriptive; mapping in the targets file, part B)")
    say("Xu and Shen's eq. 6 has priority; the feedback is Bürger and Lynch (1995). drift_critical_rate is a later discrete-generation implementation")
    say("matched to QTSim's life cycle. Nothing below is a verdict on either.")

    say("\n(i) Xu–Shen-type inputs written in QTLSim's parameters (V_E = ω²/50, mutation rate 0.1 per locus so that α²/(2V_s) is small; PopSize rounded,")
    say("    and the Xu–Shen parameters recomputed from the rounded Params). n_e = 1: a haploid clone of 100 loci; n_e = 100: a diploid with 50 loci.")
    cases = Pair{String,Params}[]
    for ne in (1, 100)
        base = (; r = 0.2, ω2 = 50.0, VM = 0.1, K = 1000.0, ne)
        push!(cases, @sprintf("n_e = %d: standard", ne) => qtl_from(; base...))
        for (k, vals) in ((:ω2, (20.0, 100.0)), (:VM, (0.05, 0.2)), (:K, (500.0, 2000.0)), (:r, (0.5,)))
            for x in vals
                push!(cases, @sprintf("n_e = %d: %s = %g", ne, k, x) => qtl_from(; merge(base, NamedTuple{(k,)}((x,)))...))
            end
        end
    end
    for (name, p) in cases[[1, 9]]
        m = xs_map(p)
        sayf("    %-22s Params: Ploidy %d, NoQTL %d, PopSize %d, R = %.4f, ω = %.4f, V_E = %.3f, m = %.3f, α² = %.4g → XuShen(r = %.4f, ω² = %.3f, V_M = %.4f, K = %.1f, n_e = %g)",
             name, p.Ploidy, p.NoQTL, p.PopSize, p.GrowthRate, initFitSD(p), QTLSim.envVar(p), p.MutationScalar, QTLSim.mutation_effect_variance(p), m.r, m.ω2, m.VM, m.K, m.ne)
    end
    ladder_table(cases)
    say("  * drift_critical_rate differs from the maximum over all steady states. In these inputs the variance of a large population carries a load close to or above")
    say("    the growth rate (V_c^M/(2ω²) against r: n_e = 100 standard 0.224 against 0.2), so viable steady states exist only at sizes where drift has reduced the")
    say("    variance; drift_critical_rate's inner iteration starts from the variance at N = K and follows one branch. Published-type inputs (V ≪ ω²) are far from this.")

    say("\n(ii), (iii) Published-type QTSim inputs (diploid, N = 1000) under Beverton–Holt, with the ladder; and the fixed-size k_c, which has no mapping")
    say("    (no dynamics of N: the feedback acts through N_e(N, N W̄) alone).")
    pubs = Pair{String,Params}[]
    push!(pubs, "test/critical.jl's input, L = 100" => Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0, PhenotypeVar = 0.03, Heritability = 0.4,
                                                              MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous))
    push!(pubs, "the same with L = 500" => Params(pubs[1][2]; NoQTL = 500))
    sim = joinpath(ROOT, "Repository", "Simulations")
    if isdir(sim)
        for (dir, f) in (("QT0100-Possion3.0", "Input-00001.txt"), ("QT0500-Possion3.0", "Input-00002.txt"), ("QT1000-Possion5.0", "Input-00003.txt"), ("QT5000-Possion3.0", "Input-00004.txt"))
            q = Params(read_qtsim_options(joinpath(sim, dir, f)); Seed = nothing)
            push!(pubs, "$dir/$f" => q)
        end
    else
        say("    (", sim, " not found: published inputs skipped)")
    end
    for (name, p) in pubs
        sayf("    %-36s L = %d, ω = %.3f, V_s = %.3f, V_E = %.4g, V_M = m V_E = %.4g, μ = %g, α² = %.4g, α²/(2V_s) = %.3g vs 1/(2N) = %.3g; house-of-Gauss V_G = %.4g",
             name, p.NoQTL, initFitSD(p), selection_variance(p), QTLSim.envVar(p), QTLSim.muteVar(p), p.MutationRatePerGamete, QTLSim.mutation_effect_variance(p),
             QTLSim.mutation_effect_variance(p) / (2selection_variance(p)), 1 / (2p.PopSize), house_of_gauss_vg(p))
    end
    cases2 = Pair{String,Params}[]
    for (name, p) in pubs, R in (1.5, 2.0, 5.0)
        push!(cases2, @sprintf("%s, R = %g", name, R) => Params(p; Demography = :beverton_holt, GrowthRate = R))
    end
    ladder_table(cases2)
    say("\n    Fixed size (Demography = :fixed): k_c from drift_critical_rate; for scale, eq. 6 at the K1 mapping of the R = 2 case and with K = ∞ (no drift, the classical rate)")
    for (name, p) in pubs
        m = xs_map(Params(p; Demography = :beverton_holt, GrowthRate = 2.0))
        sayf("    %-36s k_c(fixed) = %.4g;  k_c(BH, R = 2) = %.4g;  v_c(eq. 6, R = 2 mapping) = %.4g;  v_c(K = ∞) = %.4g", name, drift_critical_rate(p).k,
             drift_critical_rate(Params(p; Demography = :beverton_holt, GrowthRate = 2.0)).k, xu_shen_critical_rate(m).v, xu_shen_critical_rate(XuShen(m; K = Inf)).v)
    end

    say("\n    Post hoc side line (not in the targets file): effective_size(N, n) is close to the harmonic mean 2Nn/(N + n) of adults and breeders. With the breeders")
    say("    themselves as N_e (drift = n; rung 0′), all else as drift_critical_rate:")
    for (name, p) in vcat(cases[[1, 9]], cases2[[2, 5, 8, 11]])
        k0 = ladder_kc(p; variance = :hog, drift = :breeders, form = :discrete, density = :bh)
        k1 = ladder_kc(p; variance = :hog, drift = :n, form = :discrete, density = :bh)
        sayf("    %-44s k_c = %.4g; with N_e = n: %.4g (ratio %.3f)", name, k0, k1, k1 / k0)
    end

    say("\n(iv) The n_e that eq. 5 would need to hold a given equilibrium variance at v = 0: n_e = V²/(V_s(V_M − V/(2N_e))), fixed size N = 1000, published-type inputs.")
    say("     A rough estimate. V from: the house of Gauss (linkage equilibrium); multilocus_equilibrium(p; drift = true) (rare alleles + LD under QTSim's single")
    say("     crossover); its loci = :multilocus form (Gaussian + LD); for sampled published inputs also the table's pre-move V_G (one QTSim run each).")
    say("     N_e = effective_size(1000, n), n the breeders at equilibrium (as QTLSim's functions use). 2L is the value for free recombination in the Gaussian limit.")
    ne_needed(V, Vs, VM, ng) = (d = VM - V / ng; d > 0 ? V^2 / (Vs * d) : NaN)
    function row(name, p)
        Vs, VM = selection_variance(p), QTLSim.muteVar(p)
        me = multilocus_equilibrium(p; drift = true)
        mg = multilocus_equilibrium(p; drift = true, loci = :multilocus)
        ng = p.Ploidy * effective_size(p.PopSize, me.nparents)
        hog = house_of_gauss_vg(p; nsurv = me.nparents)
        sayf("    %-36s 2L = %5d | V: HoG %.4g, rare+LD %.4g (LD/genic %+.2f), Gauss+LD %.4g; neutral V_M·2N_e = %.4g | n_e: HoG %s, rare+LD %s, Gauss+LD %s", name, 2p.NoQTL,
             hog, me.VG, me.ld / me.genic, mg.VG, VM * ng, map(x -> isnan(x) ? "undef" : @sprintf("%.3g", x), (ne_needed(hog, Vs, VM, ng), ne_needed(me.VG, Vs, VM, ng),
                                                                                                             ne_needed(mg.VG, Vs, VM, ng)))...)
    end
    for (name, p) in pubs
        row(name, p)
    end
    master = joinpath(ROOT, "Repository", "Results", "Master Output file.xlsx")
    if isdir(sim) && isfile(master)
        tab = load_master(master)
        find = row_index(tab)
        rng = Xoshiro(20261118)
        say("\n     Sampled published inputs (60 per stratum, seed 20261118), V from the table's \"Genotype Var at 10000\", n from its \"P Survived at 10000\":")
        sayf("     %-20s %6s %10s %10s %10s %12s %12s %10s", "stratum", "2L", "used", "undefined", "median n_e", "quartiles", "HoG median", "n_e/2L")
        for dir in sort(filter(startswith("QT"), readdir(sim)))
            files = sort(filter(endswith(".txt"), readdir(joinpath(sim, dir))))
            pick = files[randperm(rng, length(files))[1:60]]
            L = parse(Int, dir[3:6]); u = round(Int, parse(Float64, split(dir, "Possion")[2]))
            est = Float64[]; hogs = Float64[]; nund = 0; nused = 0
            for f in pick
                q = Params(read_qtsim_options(joinpath(sim, dir, f)); Seed = nothing)
                i = find(L, u, q)
                i == 0 && continue
                tab.ns[i] >= 3 || continue
                nused += 1
                Vs, VM = selection_variance(q), QTLSim.muteVar(q)
                ng = 2 * effective_size(1000, tab.ns[i])
                x = ne_needed(tab.gv[i], Vs, VM, ng)
                isnan(x) ? (nund += 1) : push!(est, x)
                push!(hogs, ne_needed(house_of_gauss_vg(q; nsurv = tab.ns[i]), Vs, VM, ng))
            end
            qs = isempty(est) ? (NaN, NaN, NaN) : quantile(est, (0.25, 0.5, 0.75))
            sayf("     %-20s %6d %10d %10d %10.3g %5.2g–%-6.2g %12.3g %10.3g", dir, 2L, nused, nund, qs[2], qs[1], qs[3], median(hogs), qs[2] / 2L)
        end
        say("     \"undefined\": the single run's V_G is at or above the neutral expectation 2N_e V_M, so no n_e > 0 satisfies eq. 5 (one realisation; V_G fluctuates).")
    else
        say("     (published inputs or table not found: the table-based estimate is skipped)")
    end
end

function main()
    mkpath(FIGS)
    say("QTLSim validate/xu_shen.jl; QTLSim at ", pathof(QTLSim), "; threads ", Threads.nthreads())
    say("\n", "="^100, "\nThe targets file (validate/xu_shen_targets.md), as committed before the comparisons below were computed\n", "="^100)
    say(read(joinpath(DIR, "xu_shen_targets.md"), String))
    figure1()
    figure2()
    figure3()
    figure4()
    comparison()
    sayf("\nDone: %.0f s of CPU time (process clock), %.0f s wall, on %d thread(s).", ccall(:clock, Clong, ()) / 1e6, time() - T0, Threads.nthreads())
    close(LOGIO)
end

main()
