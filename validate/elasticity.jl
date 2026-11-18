# Are the emulator's gradients worth having? (Phase A of the AD plan; no AD and no change to src/.)
#
# The elasticities of the half-survival rate k50 (the rate of the continuous move at which half
# the populations are alive 2000 generations into it) with respect to the inputs that the paper's
# published sweep varies, (1) as implied by the published table itself (a logistic regression of
# survival to 2000 on the logs of the inputs, per stratum: e_θ = −β_θ/β_k), and (2) of the
# emulator (implicit-function formula on the log cumulative hazard, by central differences).
# Pre-specified criteria: validate/elasticity_criteria.md (committed before any emulator elasticity
# was computed). Writes validate/elasticity.log (mode `analyse`) and validate/elasticity_table.log
# (mode `table`).
#
# Run from the package root (vignettes environment: XLSX, Catlab):
#   julia -t 8 --project=vignettes validate/elasticity.jl table      # the table's elasticities (no emulator)
#   julia -t 8 --project=vignettes validate/elasticity.jl timing     # emulator call timings only (no elasticity)
#   julia -t 8 --project=vignettes validate/elasticity.jl emulator   # the emulator's elasticities (resumable caches)
#   julia -t 8 --project=vignettes validate/elasticity.jl analyse    # comparison and verdicts
# Caches (git-ignored): validate/elasticity_raw/.

let src = read(joinpath(@__DIR__, "emulator.jl"), String)
    # emulator.jl's helpers (read_inputs, logistic, midranks, spearman, say, …) without its main
    include_string(Main, replace(src, r"\nmain\(\)\s*$" => "\n"), joinpath(@__DIR__, "emulator.jl"))
end
using DelimitedFiles
BLAS.set_num_threads(1)                         # the bootstrap is threaded over resamples

const GMODE = isempty(ARGS) ? "analyse" : ARGS[1]
const GRAW = joinpath(@__DIR__, "elasticity_raw")
gcache(name) = (mkpath(GRAW); joinpath(GRAW, name))

# the inputs the published sweep draws at random per file (established in mode `table`), in the
# order used everywhere below; the rate k = MoveFitnessSD is the first column of every design
const GIN = (:PhenotypeVar, :Heritability, :MutationScalar, :InitFitnessSD)
const GLAB = ("V_P", "h2", "m", "omega")
const GBOOT, GBOOTSEED = 1000, 20261118          # bootstrap resamples and seed (table)
const GCELLBOOT = 600                            # bootstrap resamples of the tercile cells
const GLS = (100, 500, 1000, 5000)               # the cells pool the two recombination strata of a locus number
const GSAMPLESEED = 20261119                     # the emulator's sample
const GNPER = 16                                 # inputs per stratum (screen and :rare_alleles)
const GNSPEC = 3                                 # of them, per stratum, also default :spectrum
const GSTEP = 1e-2                               # relative step of the central differences (log θ, log k)
const GSCREEN = (; ld = :dynamic, solver = :beta)

# ---------------------------------------------------------------------------------------------
# the table joined to the input files: one row per input
# ---------------------------------------------------------------------------------------------
function g_data()
    f = gcache("inputs.tsv")
    if !isfile(f)
        ids, ps, strat = read_inputs()
        tab = load_master(MASTER)
        find = row_index(tab)
        open(f, "w") do io
            println(io, join(["id", "strat", "k", "VP", "h2", "m", "omega", "mu", "N", "T", "surv", "row"], '\t'))
            for i in eachindex(ps)
                p = ps[i]
                r = find(STRATA[strat[i]][2], STRATA[strat[i]][3], p)
                println(io, join(Any[ids[i], strat[i], p.MoveFitnessSD, p.PhenotypeVar, p.Heritability, p.MutationScalar,
                                     p.InitFitnessSD, p.MutationRatePerGamete, p.PopSize, r > 0 ? tab.T[r] : -1,
                                     r > 0 ? tab.surv[r] : -1, r], '\t'))
            end
        end
    end
    A, _ = readdlm(f, '\t'; header = true)
    col(j) = Float64.(A[:, j])
    (; id = String.(A[:, 1]), strat = Int.(A[:, 2]), k = col(3), VP = col(4), h2 = col(5), m = col(6), omega = col(7),
     mu = col(8), N = Int.(A[:, 9]), T = Int.(A[:, 10]), surv = Int.(A[:, 11]), row = Int.(A[:, 12]))
end

"The design of the linear model for the rows `ix`: 1, log k, log V_P, log h², log m, log ω."
g_design(D, ix) = [ones(length(ix)) log.(D.k[ix]) log.(D.VP[ix]) log.(D.h2[ix]) log.(D.m[ix]) log.(D.omega[ix])]

"Elasticities −β_θ/β_k of a fit β (columns as `g_design`)."
g_elast(β) = -β[3:6] ./ β[2]

# a fit is usable if it converged, its rate coefficient is negative (survival falls with the rate)
# and no coefficient is absurd (separation)
g_ok(β, conv) = conv && β[2] < 0 && all(isfinite, β) && maximum(abs, β) < 1e3

"Delta-method standard errors of the elasticities from the Fisher information of the fit."
function g_delta_se(X, β)
    μ = @. 1 / (1 + exp(-$(X * β)))
    V = inv(Symmetric(X' * ((μ .* (1 .- μ)) .* X)))
    map(3:6) do j
        g = zeros(length(β)); g[j] = -1 / β[2]; g[2] = β[j] / β[2]^2
        sqrt(max(g' * V * g, 0.0))
    end
end

"The table's fit for the rows `ix`: (β, ok, e, se (delta), lo, hi (bootstrap percentile 95%), nboot_ok)."
function g_fit(D, ix; boot = GBOOT, seed = GBOOTSEED)
    X = g_design(D, ix); y = Float64.(D.T[ix] .> 2000)
    β, conv = logistic(X, y)
    ok = g_ok(β, conv)
    e = ok ? g_elast(β) : fill(NaN, 4)
    se = ok ? g_delta_se(X, β) : fill(NaN, 4)
    lo, hi, nok = fill(NaN, 4), fill(NaN, 4), 0
    if ok && boot > 0
        n = length(ix)
        E = Vector{Vector{Float64}}(undef, boot)
        Threads.@threads :greedy for b in 1:boot
            rng = MersenneTwister(seed + 7919 * b + first(ix))
            r = rand(rng, 1:n, n)
            βb, cb = logistic(X[r, :], y[r])
            E[b] = g_ok(βb, cb) ? g_elast(βb) : fill(NaN, 4)
        end
        nok = count(v -> all(isfinite, v), E)
        for j in 1:4
            v = sort(filter(isfinite, [x[j] for x in E]))
            length(v) >= 0.9boot && ((lo[j], hi[j]) = (quantile(v, 0.025), quantile(v, 0.975)))
        end
    end
    (; β, ok, e, se, lo, hi, nok, n = length(ix), alive = count(>(0.5), y), X, y)
end

# the quadratic model: the linear terms plus all squares and products of the five centred logs
function g_quad(X)
    Z = X[:, 2:6] .- mean(X[:, 2:6]; dims = 1)
    cols = [Z[:, a] .* Z[:, b] for a in 1:5 for b in a:5]
    [X reduce(hcat, cols)], Z
end

"Local elasticities of the quadratic fit γ at the rows of X (n × 4): −(∂η/∂log θ)/(∂η/∂log k)."
function g_quad_local(γ, X, Z)
    pairs = [(a, b) for a in 1:5 for b in a:5]
    G = repeat(γ[2:6]', size(X, 1))                # ∂η/∂x_a, a = k, V_P, h², m, ω
    for (c, (a, b)) in enumerate(pairs)
        G[:, a] .+= γ[6+c] .* Z[:, b]
        G[:, b] .+= γ[6+c] .* Z[:, a]
    end
    -G[:, 2:5] ./ G[:, 1], G[:, 1]
end


const GPAIRS = [(a, b) for a in 1:5 for b in a:5]

"""
The quadratic fit γ (centred at `xbar`, the column means of the logs k, V_P, h², m, ω) at the
inputs `x4` (log V_P, log h², log m, log ω): the boundary rate (the root in log k of η = 0 nearest
to the linear fit's boundary `lk0`; `lk0` itself when there is none with ∂η/∂log k < 0) and the
local elasticities −(∂η/∂log θ)/(∂η/∂log k) there. Returns (e (4), log k, usedroot).
"""
function g_quad_at(γ, xbar, x4, lk0)
    z = [0.0; x4 .- xbar[2:5]]
    Q = zeros(5, 5)
    for (c, (a, b)) in enumerate(GPAIRS)
        Q[a, b] += γ[6+c] / 2; Q[b, a] += γ[6+c] / 2          # η = γ1 + l'x + z'Qz
    end
    l = γ[2:6]
    C = Q[1, 1]
    B = l[1] + 2 * sum(Q[1, b] * z[b] for b in 2:5)
    A = γ[1] + l[1] * xbar[1] + sum(l[a] * (z[a] + xbar[a]) for a in 2:5) + sum(Q[a, b] * z[a] * z[b] for a in 2:5, b in 2:5)
    z0 = lk0 - xbar[1]
    roots = Float64[]
    if abs(C) < 1e-12
        B != 0 && push!(roots, -A / B)
    else
        Δ = B^2 - 4C * A
        Δ >= 0 && append!(roots, [(-B - sqrt(Δ)) / 2C, (-B + sqrt(Δ)) / 2C])
    end
    roots = filter(r -> B + 2C * r < 0, roots)
    used = !isempty(roots)
    z[1] = used ? roots[argmin(abs.(roots .- z0))] : z0
    G = l .+ 2 .* (Q * z)
    (-G[2:5] ./ G[1], z[1] + xbar[1], used && G[1] < 0)
end

"The emulator's sample: per stratum `GNPER` rows drawn without replacement with probability ∝ μ(1 − μ) of the linear table fit."
function g_sample(D, fits)
    S = Int[]
    for s in eachindex(STRATA)
        ix = findall(==(s), D.strat)
        μ = @. 1 / (1 + exp(-$(fits[s].X * fits[s].β)))
        rng = MersenneTwister(GSAMPLESEED + s)
        key = [rand(rng)^(1 / max(m * (1 - m), 1e-300)) for m in μ]       # Efraimidis–Spirakis
        append!(S, ix[sortperm(key; rev = true)[1:GNPER]])
    end
    S
end

function g_wquantile(x, w, q)
    o = sortperm(x); c = cumsum(w[o]); x[o[searchsortedfirst(c, q * c[end])]]
end

"""
The primary cells' table side (criteria): for each locus number (its two strata pooled, with an
indicator of u = 5) and each input θ, the linear logistic model within each tercile of log θ,
the terciles being those of the information weights μ(1 − μ) of the per-stratum linear fits
(so each holds a third of the table's information about the boundary). Returns one row per cell.
"""
function g_cells(D, fits; boot = GCELLBOOT)
    cells = NamedTuple[]
    for (li, L) in enumerate(GLS)
        s1, s2 = 2li - 1, 2li
        ix = findall(s -> s == s1 || s == s2, D.strat)
        w = zeros(length(ix))
        for s in (s1, s2)
            μ = @. 1 / (1 + exp(-$(fits[s].X * fits[s].β)))
            w[D.strat[ix].==s] .= μ .* (1 .- μ)
        end
        for j in 1:4
            v = (D.VP, D.h2, D.m, D.omega)[j][ix]
            c1, c2 = g_wquantile(v, w, 1 / 3), g_wquantile(v, w, 2 / 3)
            for t in 1:3
                sel = t == 1 ? v .<= c1 : t == 2 ? (c1 .< v .<= c2) : v .> c2
                sub = ix[sel]
                X = [g_design(D, sub) Float64.(D.strat[sub] .== s2)]; y = Float64.(D.T[sub] .> 2000)
                β, conv = logistic(X, y)
                ok = g_ok(β, conv)
                e = ok ? -β[2+j] / β[2] : NaN
                n = length(sub)
                eb = fill(NaN, boot)
                ok && Threads.@threads :greedy for b in 1:boot
                    rng = MersenneTwister(GBOOTSEED + 2 + 15485863 * b + 1000 * li + 10 * j + t)
                    r = rand(rng, 1:n, n)
                    βb, cb = logistic(X[r, :], y[r])
                    g_ok(βb, cb) && (eb[b] = -βb[2+j] / βb[2])
                end
                vb = sort(filter(isfinite, eb))
                def = ok && length(vb) >= 0.9boot
                push!(cells, (; L, li, j, t, lo_θ = t == 1 ? minimum(v) : t == 2 ? c1 : c2, hi_θ = t == 1 ? c1 : t == 2 ? c2 : maximum(v),
                              e, lo = def ? quantile(vb, 0.025) : NaN, hi = def ? quantile(vb, 0.975) : NaN, n, alive = count(>(0.5), y),
                              info = sum(w[sel]) / sum(w), nboot = length(vb), ok = def))
            end
        end
    end
    cells
end

"The tercile (1–3) of the value θ of input j at locus-number index li, from the cells' cut points."
function g_tercile(cells, li, j, θ)
    c = [x for x in cells if x.li == li && x.j == j]
    θ <= c[1].hi_θ ? 1 : θ <= c[2].hi_θ ? 2 : 3
end

function table_mode()
    LOG[] = open(joinpath(@__DIR__, "elasticity_table.log"), "w")
    say("validate/elasticity.jl table — the published table's own elasticities of the half-survival rate")
    say(Libc.strftime("%Y-%m-%d %H:%M", time()), "; $(Threads.nthreads()) threads; no emulator is called in this mode")
    D = g_data()
    n = length(D.id)
    say("inputs: $n files; matched to a unique table row: $(count(>(0), D.row))")
    all(>(0), D.row) || error("unmatched inputs")

    say("\n== which inputs vary (all strata pooled; per-file Params from read_qtsim_options) ==")
    for (nm, v) in (("MoveFitnessSD k", D.k), ("PhenotypeVar", D.VP), ("Heritability", D.h2), ("MutationScalar", D.m),
                    ("InitFitnessSD ω", D.omega), ("MutationRatePerGamete", D.mu), ("PopSize", Float64.(D.N)))
        sayf("  %-22s distinct %6d   min %-10.5g q05 %-10.5g median %-10.5g q95 %-10.5g max %-10.5g", nm, length(unique(v)),
             minimum(v), quantile(v, 0.05), median(v), quantile(v, 0.95), maximum(v))
    end
    say("  (ω is given per file as InitFitnessSD: it is drawn, not derived from h², V_P and m;")
    say("   analytical ω = h²V_P/√(m·V_E) for comparison: median ",
        @sprintf("%.3g", median(@. D.h2 * D.VP / sqrt(D.m * (1 - D.h2) * D.VP))), ")")
    say("\n  Spearman correlations of the inputs (pooled), and the largest |r| within any stratum")
    V = (D.k, D.VP, D.h2, D.m, D.omega); L = ("k", GLAB...)
    for a in 1:5, b in a+1:5
        r = spearman(V[a], V[b])
        rs = maximum(abs(spearman(V[a][D.strat.==s], V[b][D.strat.==s])) for s in eachindex(STRATA))
        sayf("    %-6s %-6s  % .4f   max |r| in a stratum %.4f", L[a], L[b], r, rs)
    end
    say("\n  alive at 2000 (table 'Generations to extinction' > 2000) by stratum, and the rate's range among them")
    for s in eachindex(STRATA)
        ix = findall(==(s), D.strat); al = ix[D.T[ix].>2000]
        sayf("    %-18s n %5d  alive %5d (%.3f)  k of the alive: q05 %.4g median %.4g q95 %.4g; all: %.4g–%.4g",
             STRATA[s][1], length(ix), length(al), length(al) / length(ix), quantile(D.k[al], 0.05), median(D.k[al]),
             quantile(D.k[al], 0.95), minimum(D.k[ix]), maximum(D.k[ix]))
    end

    say("\n== the linear logistic model per stratum: logit P(alive at 2000) = β0 + βk log k + Σ βθ log θ ==")
    say("   rate in raw units (k = MoveFitnessSD per generation); the elasticity of the standardised rate")
    say("   x50 = k50/√V_P in V_P is e_VP − 1/2, the others are unchanged.")
    say("   e_θ = −βθ/βk with [bootstrap percentile 95% CI, $(GBOOT) resamples of the rows, seed $(GBOOTSEED)] and (delta-method SE)")
    fits = Any[]
    out = open(gcache("table_elasticities.tsv"), "w")
    println(out, join(["strat", "input", "e", "se", "lo", "hi", "ok", "betak"], '\t'))
    for s in eachindex(STRATA)
        ix = findall(==(s), D.strat)
        f = g_fit(D, ix)
        push!(fits, f)
        sayf("  %-18s n %5d alive %5d  ok %s  βk %.3f  β0 %.3f  bootstrap fits usable %d/%d; log k50 at the median inputs %.4f",
             STRATA[s][1], f.n, f.alive, f.ok, f.β[2], f.β[1], f.nok, GBOOT,
             -(f.β[1] + sum(f.β[3:6] .* log.([median(D.VP[ix]), median(D.h2[ix]), median(D.m[ix]), median(D.omega[ix])]))) / f.β[2])
        for j in 1:4
            sayf("       e_%-6s % .4f  [% .4f, % .4f]  (se %.4f)%s", GLAB[j], f.e[j], f.lo[j], f.hi[j], f.se[j],
                 (f.lo[j] > 0 || f.hi[j] < 0) ? "  CI excludes 0" : "")
            println(out, join(Any[s, GLAB[j], f.e[j], f.se[j], f.lo[j], f.hi[j], f.ok, f.β[2]], '\t'))
        end
    end
    close(out)

    say("\n== adequacy 1: quadratic model (all squares and products of the five centred logs, 15 terms) ==")
    say("   likelihood-ratio statistic against the linear model (χ², 15 df; 5% point 25.0, 0.1% point 37.7), and the")
    say("   local elasticities of the quadratic fit over the inputs near the boundary (linear-fit P in [0.1, 0.9]):")
    say("   median [10%, 90% quantile]; and the linear fit's e for comparison")
    ll(X, β, y) = (η = X * β; sum(@. y * η - log1p(exp(-abs(η))) - max(η, 0)))
    for s in eachindex(STRATA)
        f = fits[s]
        Xq, Z = g_quad(f.X)
        γ, cq = logistic(Xq, f.y)
        lr = 2 * (ll(Xq, γ, f.y) - ll(f.X, f.β, f.y))
        μ = @. 1 / (1 + exp(-$(f.X * f.β)))
        near = findall(x -> 0.1 <= x <= 0.9, μ)
        E, gk = g_quad_local(γ, Xq[near, :], Z[near, :])
        sayf("  %-18s converged %s  LR %.1f  near-boundary inputs %d  (∂η/∂log k > 0 at %d of them)", STRATA[s][1], cq, lr,
             length(near), count(>(0), gk))
        for j in 1:4
            v = E[:, j]
            sayf("       e_%-6s quadratic local median % .3f [% .3f, % .3f]   linear % .3f", GLAB[j], median(v),
                 quantile(v, 0.1), quantile(v, 0.9), f.e[j])
        end
    end

    say("\n== adequacy 2: the linear model within terciles of each input (rows of the stratum split by the input's terciles) ==")
    say("   e_θ of the input θ itself in its low / middle / high tercile (delta-method SE), then the spread max − min")
    out = open(gcache("table_terciles.tsv"), "w")
    println(out, join(["strat", "input", "tercile", "lo", "hi", "e", "se", "ok"], '\t'))
    nvar = 0; ncell = 0
    for s in eachindex(STRATA)
        ix = findall(==(s), D.strat)
        say("  ", STRATA[s][1])
        for j in 1:4
            v = (D.VP, D.h2, D.m, D.omega)[j][ix]
            q = quantile(v, [0, 1 / 3, 2 / 3, 1])
            es = Float64[]; ss = Float64[]
            for t in 1:3
                sub = ix[(t == 1 ? v .<= q[2] : t == 2 ? (q[2] .< v .<= q[3]) : v .> q[3])]
                f = g_fit(D, sub; boot = 0)
                push!(es, f.e[j]); push!(ss, f.se[j])
                println(out, join(Any[s, GLAB[j], t, q[t], q[t+1], f.e[j], f.se[j], f.ok], '\t'))
            end
            spread = maximum(es) - minimum(es)
            # heterogeneity: the largest pairwise |difference| / SE of the difference
            z = maximum(abs(es[a] - es[b]) / sqrt(ss[a]^2 + ss[b]^2) for a in 1:3 for b in a+1:3)
            ncell += 1; (z > 3) && (nvar += 1)
            sayf("       e_%-6s % .3f (%.3f)  % .3f (%.3f)  % .3f (%.3f)   spread %.3f  largest |z| %.1f%s", GLAB[j], es[1], ss[1],
                 es[2], ss[2], es[3], ss[3], spread, z, z > 3 ? "  VARIES" : "")
        end
    end
    close(out)
    say("\n  cells (stratum × input) whose elasticity differs between two of its terciles by more than 3 SE: $nvar of $ncell")
    say("\n== adequacy 3: is the quadratic model's local elasticity consistent with the tercile fits? ==")
    say("   mean over the tercile's rows, weighted by μ(1 − μ) of the quadratic fit, of the quadratic local e_θ at each row's")
    say("   inputs and boundary rate, beside the tercile's own linear fit (adequacy 2); low / middle / high")
    quads = Any[]
    for s in eachindex(STRATA)
        f = fits[s]; ix = findall(==(s), D.strat)
        Xq, _ = g_quad(f.X)
        γ, _ = logistic(Xq, f.y)
        xbar = vec(mean(f.X[:, 2:6]; dims = 1))
        push!(quads, (; γ, xbar))
        μ = @. 1 / (1 + exp(-$(Xq * γ)))
        w = μ .* (1 .- μ)
        loc = [g_quad_at(γ, xbar, f.X[i, 3:6], -(f.β[1] + sum(f.β[3:6] .* f.X[i, 3:6])) / f.β[2]) for i in eachindex(ix)]
        say("  ", STRATA[s][1], "   (boundary root of the quadratic found for ", count(x -> x[3], loc), " of ", length(loc), " rows)")
        for j in 1:4
            v = (D.VP, D.h2, D.m, D.omega)[j][ix]
            q = quantile(v, [0, 1 / 3, 2 / 3, 1])
            cells = map(1:3) do t
                sub = findall(t == 1 ? v .<= q[2] : t == 2 ? (q[2] .< v .<= q[3]) : v .> q[3])
                sub = [i for i in sub if loc[i][3]]
                qm = sum(w[i] * loc[i][1][j] for i in sub) / sum(w[sub])
                ft = g_fit(D, ix[findall(t == 1 ? v .<= q[2] : t == 2 ? (q[2] .< v .<= q[3]) : v .> q[3])]; boot = 0)
                @sprintf("% .3f | % .3f (%.3f)", qm, ft.e[j], ft.se[j])
            end
            sayf("       e_%-6s quadratic | tercile fit:   %s    %s    %s", GLAB[j], cells...)
        end
    end

    say("\n== the emulator's sample and the table's local elasticities on it ==")
    say("   per stratum $(GNPER) inputs drawn without replacement with probability ∝ μ(1 − μ) of the linear fit above")
    say("   (MersenneTwister($(GSAMPLESEED) + stratum), Efraimidis–Spirakis keys): the inputs that carry the table's information")
    say("   about the boundary. Table truth per input: the quadratic fit's local elasticity at the input's (V_P, h², m, ω)")
    say("   and at the quadratic fit's boundary rate; $(GBOOT) bootstrap refits of the quadratic model (rows resampled, seed $(GBOOTSEED) + 1).")
    S = g_sample(D, fits)
    E = fill(NaN, length(S), 4); LK = fill(NaN, length(S)); used = falses(length(S))
    EB = fill(NaN, GBOOT, length(S), 4)
    for s in eachindex(STRATA)
        f = fits[s]; q = quads[s]
        Xq, _ = g_quad(f.X)
        pos = findall(i -> D.strat[S[i]] == s, eachindex(S))
        x4 = [[log(D.VP[S[i]]), log(D.h2[S[i]]), log(D.m[S[i]]), log(D.omega[S[i]])] for i in pos]
        lk0 = [-(f.β[1] + sum(f.β[3:6] .* x)) / f.β[2] for x in x4]
        for (a, i) in enumerate(pos)
            e, lk, u = g_quad_at(q.γ, q.xbar, x4[a], lk0[a])
            E[i, :] .= e; LK[i] = lk; used[i] = u
        end
        n = length(f.y)
        Threads.@threads :greedy for b in 1:GBOOT
            rng = MersenneTwister(GBOOTSEED + 1 + 104729 * b + s)
            r = rand(rng, 1:n, n)
            Xb = f.X[r, :]
            Xqb, _ = g_quad(Xb)
            γb, cb = logistic(Xqb, f.y[r])
            βb, cl = logistic(Xb, f.y[r])
            (cb && cl && g_ok(βb, cl)) || continue
            xb = vec(mean(Xb[:, 2:6]; dims = 1))
            for (a, i) in enumerate(pos)
                e, _, u = g_quad_at(γb, xb, x4[a], -(βb[1] + sum(βb[3:6] .* x4[a])) / βb[2])
                u && (EB[b, i, :] .= e)
            end
        end
    end
    open(gcache("sample.tsv"), "w") do io
        println(io, join(["id", "strat", "order", "k", "VP", "h2", "m", "omega", "T", "lk50_quad", "root", "q_VP", "q_h2", "q_m", "q_omega"], '\t'))
        ord = 0; last = 0
        for (i, r) in enumerate(S)
            ord = D.strat[r] == last ? ord + 1 : 1; last = D.strat[r]
            println(io, join(Any[D.id[r], D.strat[r], ord, D.k[r], D.VP[r], D.h2[r], D.m[r], D.omega[r], D.T[r], LK[i], used[i], E[i, :]...], '\t'))
        end
    end
    open(io -> write(io, EB), gcache("sample_boot.bin"), "w")
    say("   sample: $(length(S)) inputs; quadratic boundary root found for $(count(used)); alive at 2000 in the table: $(count(D.T[S] .> 2000))")
    say("   ids: ", join(D.id[S], " "))
    say("\n   table-local elasticities on the sample: median over the stratum's inputs [bootstrap 95% CI of that median], and the mean")
    for s in eachindex(STRATA)
        pos = findall(i -> D.strat[S[i]] == s, eachindex(S))
        cells = map(1:4) do j
            bm = sort(filter(isfinite, [median(EB[b, pos, j]) for b in 1:GBOOT]))
            @sprintf("%s % .3f [% .3f, % .3f] mean % .3f", GLAB[j], median(E[pos, j]), quantile(bm, 0.025), quantile(bm, 0.975), mean(E[pos, j]))
        end
        say("     ", rpad(STRATA[s][1], 18), " ", join(cells, "   "))
    end
    say("\n== the primary cells (criteria): tercile fits by locus number, input and information tercile ==")
    say("   per locus number (u = 3 and 5 pooled, with an indicator of u = 5) and input θ: the linear logistic model within each")
    say("   tercile of θ, terciles of the information weights μ(1 − μ) of the per-stratum linear fits; e_θ of θ itself")
    say("   [bootstrap percentile 95% CI, $(GCELLBOOT) resamples of the cell's rows]; n rows, alive, share of the information,")
    say("   and the number of the emulator's sample inputs in the cell")
    cells = g_cells(D, fits)
    open(gcache("table_cells.tsv"), "w") do io
        println(io, join(["L", "li", "j", "input", "t", "lo_theta", "hi_theta", "e", "lo", "hi", "n", "alive", "info", "nboot", "ok"], '\t'))
        for c in cells
            println(io, join(Any[c.L, c.li, c.j, GLAB[c.j], c.t, c.lo_θ, c.hi_θ, c.e, c.lo, c.hi, c.n, c.alive, c.info, c.nboot, c.ok], '\t'))
        end
    end
    Sv = (D.VP, D.h2, D.m, D.omega)
    for c in cells
        ns = count(r -> (D.strat[r] + 1) ÷ 2 == c.li && g_tercile(cells, c.li, c.j, Sv[c.j][r]) == c.t, S)
        sayf("   L %4d  %-6s tercile %d  θ in [%.4g, %.4g]  e % .3f [% .3f, % .3f]%s  n %5d alive %4d info %.2f  sample %2d%s",
             c.L, GLAB[c.j], c.t, c.lo_θ, c.hi_θ, c.e, c.lo, c.hi, c.ok ? ((c.lo > 0 || c.hi < 0) ? " *" : "  ") : " UNDEFINED",
             c.n, c.alive, c.info, ns, abs(c.e) >= 0.1 ? "" : "  (|e| < 0.1)")
    end
    say("   cells defined: $(count(c -> c.ok, cells)) of $(length(cells)); CI excludes 0 (*): $(count(c -> c.ok && (c.lo > 0 || c.hi < 0), cells)); ",
        "|e| ≥ 0.1: $(count(c -> c.ok && abs(c.e) >= 0.1, cells))")
    close(LOG[])
end

GMODE == "table" && table_mode()

# ---------------------------------------------------------------------------------------------
# timings of single emulator calls (no elasticity, no outcome looked at): for the budget
# ---------------------------------------------------------------------------------------------
function timing_mode()
    D = g_data()
    say("timing: one emulate call per closure at half the :rare_alleles k50, Input-00000 and Input-00001 of each stratum; seconds")
    jobs = [(s, j) for s in eachindex(STRATA) for j in 0:1]
    p0 = Params(read_qtsim_options(joinpath(SIMDIR, STRATA[1][1], "Input-00000.txt")); Seed = nothing)
    for sp in (GSCREEN, (; ld = :dynamic), (;))            # compile
        emulate(Params(p0; MoveFitnessSD = 1e-3); loci = :spectrum, spectrum = sp, burnin = 50, horizon = 50)
    end
    emulated_half_survival_rate(p0)
    rs = Vector{Any}(undef, length(jobs))
    Threads.@threads :greedy for i in eachindex(jobs)
        s, j = jobs[i]
        p = Params(read_qtsim_options(joinpath(SIMDIR, STRATA[s][1], @sprintf("Input-%05d.txt", j))); Seed = nothing)
        tr = @elapsed h = emulated_half_survival_rate(p; rtol = 1e-6)
        q = Params(p; MoveFitnessSD = 0.5 * (0 < h.k50 < Inf ? h.k50 : p.MoveFitnessSD))
        tb = @elapsed emulate(q; loci = :spectrum, spectrum = GSCREEN)
        tg = @elapsed emulate(q; loci = :spectrum, spectrum = (; ld = :dynamic))
        ts = @elapsed emulate(q; loci = :spectrum)
        rs[i] = (; s, j, tr, evals = h.evals, tb, tg, ts)
    end
    for r in rs
        sayf("  %-18s %d  :rare_alleles k50 %.2f s (%d evals)   screen (beta, ld) %.2f   grid, ld %.2f   default :spectrum %.2f",
             STRATA[r.s][1], r.j, r.tr, r.evals, r.tb, r.tg, r.ts)
    end
end

GMODE == "timing" && timing_mode()

# ---------------------------------------------------------------------------------------------
# the emulator's elasticities
# ---------------------------------------------------------------------------------------------
const GCLOS = (screen = (:spectrum, GSCREEN),                 # the two-stage screen: dynamic LD, beta solver
               grid = (:spectrum, (; ld = :dynamic)),         # the same closure on the frequency grid (no fallback switch)
               spectrum = (:spectrum, (;)),                   # default :spectrum
               rare = (:rare_alleles, (;)))
const GROUNDS, GMINROUNDS = 12, 8        # screen: inputs per stratum, in draw order ("rounds" of one input per stratum)
const GSCREENBUDGET = 6000.0             # seconds summed over the screen's jobs after which no new round starts (≥ GMINROUNDS)
const GBRUTE = (1, 3, 5, 7, 8)           # strata whose first sample input is re-solved by brute force
const GBRUTESTEP = 0.05
const GREGRID_STATUS, GREGRID_GEN = 10, 3

"The search's function of P2000: the log cumulative hazard relative to log 2 (0 at P2000 = 1/2; decreasing in the rate)."
g_hfun(P) = log(log(2.0)) - log(-log(P))

function g_eval(p::Params, c::Symbol)
    loci, spec = GCLOS[c]
    e = emulate(p; loci, spectrum = spec, burnin = BURNIN, horizon = HORIZON)
    (g_hfun(e.P2000), Float64(get(e, :fallback, -1)))
end

# input 0 is the rate; 1–4 are GIN. Each is one Params field: ω is InitFitnessSD, given per file
# (not derived), so the others do not move with it, as in the table's generating process.
g_setp(p::Params, j::Int, f::Real) = j == 0 ? Params(p; MoveFitnessSD = p.MoveFitnessSD * f) :
                                     Params(p; (GIN[j] => getfield(p, GIN[j]) * f,)...)

function g_derivs(p::Params, k0, c, δ)
    q = Params(p; MoveFitnessSD = k0)
    hp = zeros(5); hm = zeros(5); fp = zeros(5); fm = zeros(5)
    for j in 0:4
        hp[j+1], fp[j+1] = g_eval(g_setp(q, j, exp(δ)), c)
        hm[j+1], fm[j+1] = g_eval(g_setp(q, j, exp(-δ)), c)
    end
    (; hp, hm, fp, fm)
end

function g_k50(p::Params, c::Symbol; hint = NaN, rtol = 1e-2, width = NaN)
    loci, spec = GCLOS[c]
    c === :rare && isnan(width) && return emulated_half_survival_rate(p; burnin = BURNIN, horizon = HORIZON, rtol = 1e-6)
    bracket = !(0 < hint < Inf) ? nothing : isnan(width) ? (hint / 4, hint) : (hint * exp(-width), hint * exp(width))
    emulated_half_survival_rate(p; loci, spectrum = spec, burnin = BURNIN, horizon = HORIZON, rtol, bracket)
end

"""
One input and closure: k50 (unless `k0` is given), then the elasticities by the implicit-function
formula e_θ = −(∂h/∂log θ)/(∂h/∂log k) at k50, central differences with step δ in the logs.
status 0 ok, 1 no k50 (0 or Inf), 2 P2000 saturated at a difference point (also at δ/10), 3 ∂h/∂log k ≥ 0.
"""
function g_job(p::Params, c::Symbol; hint = NaN, k0 = NaN, δ = GSTEP)
    evals = 0; status = 0; step = δ
    e = fill(NaN, 4); hk = NaN; h0 = NaN; f0 = -1.0
    d = (; hp = fill(NaN, 5), hm = fill(NaN, 5), fp = fill(-1.0, 5), fm = fill(-1.0, 5))
    t = @elapsed begin
        if isnan(k0)
            h = g_k50(p, c; hint)
            k0 = h.k50; evals = h.evals
        end
        if !(0 < k0 < Inf)
            status = 1
        else
            h0, f0 = g_eval(Params(p; MoveFitnessSD = k0), c)
            d = g_derivs(p, k0, c, step)
            if !all(isfinite, d.hp) || !all(isfinite, d.hm)
                step = δ / 10
                d = g_derivs(p, k0, c, step)
            end
            if !all(isfinite, d.hp) || !all(isfinite, d.hm)
                status = 2
            else
                hk = (d.hp[1] - d.hm[1]) / 2step
                hk < 0 ? (e = -((d.hp[2:5] .- d.hm[2:5]) ./ 2step) ./ hk) : (status = 3)
            end
        end
    end
    sstat = count(j -> (d.fp[j] >= 0) != (d.fm[j] >= 0), 1:5)                         # pairs straddling a fallback switch
    sgen = count(j -> d.fp[j] >= 0 && d.fm[j] >= 0 && d.fp[j] != d.fm[j], 1:5)        # pairs falling back at different generations
    (; k50 = k0, evals, status, step, h0, fb0 = f0, hk, e1 = e[1], e2 = e[2], e3 = e[3], e4 = e[4], sstat, sgen,
     nfb = count(>=(0), d.fp) + count(>=(0), d.fm), seconds = t,
     hp1 = d.hp[1], hp2 = d.hp[2], hp3 = d.hp[3], hp4 = d.hp[4], hp5 = d.hp[5],
     hm1 = d.hm[1], hm2 = d.hm[2], hm3 = d.hm[3], hm4 = d.hm[4], hm5 = d.hm[5],
     fp1 = d.fp[1], fp2 = d.fp[2], fp3 = d.fp[3], fp4 = d.fp[4], fp5 = d.fp[5],
     fm1 = d.fm[1], fm2 = d.fm[2], fm3 = d.fm[3], fm4 = d.fm[4], fm5 = d.fm[5])
end

"Brute force: k50 re-solved at θ·e^{±GBRUTESTEP} for each input, from a bracket around the base k50."
function g_brute(p::Params, c::Symbol, k0)
    e = fill(NaN, 4); n = 0
    t = @elapsed for j in 1:4
        a = g_k50(g_setp(p, j, exp(GBRUTESTEP)), c; hint = k0, rtol = 1e-3, width = 0.15)
        b = g_k50(g_setp(p, j, exp(-GBRUTESTEP)), c; hint = k0, rtol = 1e-3, width = 0.15)
        n += a.evals + b.evals
        (0 < a.k50 < Inf && 0 < b.k50 < Inf) && (e[j] = (log(a.k50) - log(b.k50)) / 2GBRUTESTEP)
    end
    (; b1 = e[1], b2 = e[2], b3 = e[3], b4 = e[4], evals = n, seconds = t)
end

function g_read_sample()
    A, _ = readdlm(gcache("sample.tsv"), '\t'; header = true)
    f(j) = Float64.(A[:, j])
    (; id = String.(A[:, 1]), strat = Int.(A[:, 2]), order = Int.(A[:, 3]), k = f(4), VP = f(5), h2 = f(6), m = f(7), omega = f(8),
     T = Int.(A[:, 9]), lkq = f(10), root = Bool.(A[:, 11]), q = Float64.(A[:, 12:15]))
end

g_params(id) = Params(read_qtsim_options(joinpath(SIMDIR, id)); Seed = nothing)

function emulator_mode()
    t0 = time()
    S = g_read_sample()
    ps = Dict(id => g_params(id) for id in S.id)
    say("emulator mode: $(Threads.nthreads()) threads; QTLSim at ", pathof(QTLSim))
    # compile
    let p = ps[S.id[1]]
        for c in keys(GCLOS)
            emulate(Params(p; MoveFitnessSD = 1e-3); loci = GCLOS[c][1], spectrum = GCLOS[c][2], burnin = 20, horizon = 20)
        end
    end
    # (a) :rare_alleles on the whole sample (cheap; its k50 is the bracket hint of the others)
    rare = sweep(id -> g_job(ps[id], :rare), S.id; cache = gcache("rare.csv"), progress = 32)
    hint = Dict(S.id[i] => rare[i].k50 for i in eachindex(S.id))
    say("  :rare_alleles: $(length(rare)) inputs, ", round(Int, sum(r.seconds for r in rare)), " s summed")
    # (b) the screen, round by round, under the pre-specified budget
    used = 0.0; nround = 0
    for r in 1:GROUNDS
        (r > GMINROUNDS && used > GSCREENBUDGET) && break
        ids = S.id[S.order.==r]
        rs = sweep(id -> g_job(ps[id], :screen; hint = hint[id]), ids; cache = gcache("screen.csv"))
        used += sum(x.seconds for x in rs); nround = r
        sayf("  screen round %2d: %d inputs, summed job seconds so far %.0f (wall %.0f s)", r, length(ids), used, time() - t0)
    end
    open(io -> println(io, nround), gcache("rounds.txt"), "w")
    run = [i for i in eachindex(S.id) if S.order[i] <= nround]
    scr = sweep(id -> error("not cached"), S.id[run]; cache = gcache("screen.csv"))
    # (c) default :spectrum on round 1
    r1 = S.id[S.order.==1]
    sp = sweep(id -> g_job(ps[id], :spectrum; hint = hint[id]), r1; cache = gcache("spectrum.csv"))
    say("  default :spectrum: $(length(sp)) inputs, ", round(Int, sum(r.seconds for r in sp)), " s summed")
    # (d) numerical checks of the screen on round 1: a second step, and brute force on GBRUTE
    k1 = Dict(S.id[i] => scr[findfirst(==(i), run)].k50 for i in run if S.order[i] == 1)
    st = sweep(id -> g_job(ps[id], :screen; k0 = k1[id], δ = GSTEP / 10), [id for id in r1 if 0 < k1[id] < Inf];
               cache = gcache("screen_step.csv"))
    bid = [id for id in r1 if S.strat[findfirst(==(id), S.id)] in GBRUTE && 0 < k1[id] < Inf]
    br = sweep(id -> g_brute(ps[id], :screen, k1[id]), bid; cache = gcache("screen_brute.csv"))
    say("  step check: ", round(Int, sum(r.seconds for r in st; init = 0.0)), " s; brute force: ",
        round(Int, sum(r.seconds for r in br; init = 0.0)), " s summed")
    # (e) fallback switches: recompute with the grid at the screen's k50
    ss = [S.id[i] for (a, i) in enumerate(run) if scr[a].status == 0 && scr[a].sstat > 0]
    sg = [S.id[i] for (a, i) in enumerate(run) if scr[a].status == 0 && scr[a].sstat == 0 && scr[a].sgen > 0]
    kk = Dict(S.id[i] => scr[a].k50 for (a, i) in enumerate(run))
    todo = [ss[1:min(end, GREGRID_STATUS)]; sg[1:min(end, GREGRID_GEN)]]
    function regrid(id)
        r = g_job(ps[id], :grid; k0 = kk[id])
        r.status == 2 || return merge(r, (; resolved = false))
        h = g_k50(ps[id], :grid; hint = kk[id], width = 0.2)         # the grid saturates at the beta's k50: its own k50
        merge(g_job(ps[id], :grid; k0 = h.k50), (; resolved = true))
    end
    rg = sweep(regrid, todo; cache = gcache("screen_grid.csv"))
    say("  fallback: status straddles $(length(ss)) inputs, generation straddles $(length(sg)); recomputed with the grid: $(length(todo)), ",
        round(Int, sum(r.seconds for r in rg; init = 0.0)), " s summed")
    sayf("emulator mode done: wall %.0f s", time() - t0)
end

GMODE == "emulator" && emulator_mode()

# ---------------------------------------------------------------------------------------------
# comparison and verdicts (validate/elasticity_criteria.md)
# ---------------------------------------------------------------------------------------------
"A sweep cache as Dict id => Dict column => Float64 (Bools as 0/1)."
function g_cache(name)
    out = Dict{String,Dict{Symbol,Float64}}()
    isfile(gcache(name)) || return out
    lines = readlines(gcache(name))
    hdr = Symbol.(split(lines[1], ','))
    for l in lines[2:end]
        f = split(l, ',')
        length(f) == length(hdr) || continue
        out[String(f[1])] = Dict(hdr[c] => (f[c] == "true" ? 1.0 : f[c] == "false" ? 0.0 : parse(Float64, f[c])) for c in 2:length(hdr))
    end
    out
end
g_e(r) = r[:status] == 0 ? [r[:e1], r[:e2], r[:e3], r[:e4]] : fill(NaN, 4)
fin(v) = filter(isfinite, v)
med(v) = (w = fin(v); isempty(w) ? NaN : median(w))
avg(v) = (w = fin(v); isempty(w) ? NaN : mean(w))

function g_read_cells()
    A, _ = readdlm(gcache("table_cells.tsv"), '\t'; header = true)
    [(; L = Int(A[r, 1]), li = Int(A[r, 2]), j = Int(A[r, 3]), t = Int(A[r, 5]), lo_θ = Float64(A[r, 6]), hi_θ = Float64(A[r, 7]),
        e = Float64(A[r, 8]), lo = Float64(A[r, 9]), hi = Float64(A[r, 10]), ok = Bool(A[r, 15])) for r in axes(A, 1)]
end

"G1–G3 for the elasticities E (inputs × 4; NaN = undefined) of the sample rows `run` against the table's cells."
function g_gstats(name, E, S, run, cells; verdict::Bool, minn = 4)
    say("\n-- ", name, verdict ? " (VERDICT)" : " (comparator, no verdict)", " --")
    say("   cell: L, input, tercile; table e [CI]; emulator median (mean) over n defined of N inputs; G1 sign, G2 size")
    Sv = (S.VP, S.h2, S.m, S.omega)
    te = Float64[]; ee = Float64[]; g1 = Bool[]; g2 = Bool[]; fails = String[]; excl = String[]
    for c in cells
        c.ok || continue
        ins = [a for (a, i) in enumerate(run) if (S.strat[i] + 1) ÷ 2 == c.li && g_tercile(cells, c.li, c.j, Sv[c.j][i]) == c.t]
        v = fin(E[ins, c.j])
        lab = @sprintf("L %4d %-6s t%d", c.L, GLAB[c.j], c.t)
        if length(v) < minn
            push!(excl, lab * " (n = $(length(v)))")
            sayf("   %s  table % .3f [% .3f, % .3f]   emulator: %d defined of %d — EXCLUDED", lab, c.e, c.lo, c.hi, length(v), length(ins))
            continue
        end
        m = median(v); push!(te, c.e); push!(ee, m)
        s1 = ""; s2 = ""
        if c.lo > 0 || c.hi < 0
            ok = sign(m) == sign(c.e); push!(g1, ok); s1 = ok ? "sign ok" : "SIGN WRONG"
        end
        if abs(c.e) >= 0.1
            tol = max(0.15, 0.30 * abs(c.e))
            ok = c.lo - tol <= m <= c.hi + tol; push!(g2, ok); s2 = ok ? "size ok" : "SIZE FAIL"
            ok || push!(fails, @sprintf("%s: table % .3f [% .3f, % .3f], tol %.3f, emulator % .3f", lab, c.e, c.lo, c.hi, tol, m))
        end
        sayf("   %s  table % .3f [% .3f, % .3f]   emulator % .3f (% .3f)  n %2d/%2d   %-10s %s", lab, c.e, c.lo, c.hi, m, mean(v), length(v),
             length(ins), s1, s2)
    end
    r = length(te) >= 3 ? spearman(te, ee) : NaN
    p1 = isempty(g1) ? NaN : count(g1) / length(g1); p2 = isempty(g2) ? NaN : count(g2) / length(g2)
    v(b) = verdict ? pf(b) : (b ? "(would pass)" : "(would fail)")
    sayf("   G1 sign: %d of %d cells whose table CI excludes 0 (%.3f; ≥ 0.90): %s", count(g1), length(g1), p1, v(p1 >= 0.9))
    sayf("   G2 size: %d of %d cells with |table e| ≥ 0.1 (%.3f; ≥ 0.75): %s", count(g2), length(g2), p2, v(p2 >= 0.75))
    sayf("   G3 ranking: Spearman over %d cells = %.3f (≥ 0.8): %s", length(te), r, v(r >= 0.8))
    isempty(fails) || (say("   failing G2 cells:"); foreach(f -> say("     ", f), fails))
    isempty(excl) || say("   excluded cells (< $minn defined inputs): ", join(excl, "; "))
    verdict && say("   OVERALL GATE (G1 and G2): ", pf(p1 >= 0.9 && p2 >= 0.75))
    (; p1, p2, r)
end

function analyse_mode()
    LOG[] = open(joinpath(@__DIR__, "elasticity.log"), "w")
    say("validate/elasticity.jl analyse — the emulator's elasticities of the half-survival rate against the published table's")
    say("criteria: validate/elasticity_criteria.md; table side: validate/elasticity_table.log; ", Libc.strftime("%Y-%m-%d %H:%M", time()))
    S = g_read_sample(); cells = g_read_cells()
    nround = parse(Int, strip(read(gcache("rounds.txt"), String)))
    run = [i for i in eachindex(S.id) if S.order[i] <= nround]
    scr, rare, spc, stp, brt, grd = g_cache("screen.csv"), g_cache("rare.csv"), g_cache("spectrum.csv"), g_cache("screen_step.csv"),
                                    g_cache("screen_brute.csv"), g_cache("screen_grid.csv")
    say("\n== the runs ==")
    say("screen: $(nround) rounds of one input per stratum = $(length(run)) inputs (sample of 16 per stratum, draw order); ",
        ":rare_alleles on all $(length(S.id)); default :spectrum on $(length(spc)) (round 1)")
    tsec(d) = sum(r[:seconds] for r in values(d); init = 0.0)
    sayf("summed job seconds: screen %.0f, :rare_alleles %.0f, default :spectrum %.0f, step check %.0f, brute force %.0f, grid recomputation %.0f; total %.0f s = %.2f h",
         tsec(scr), tsec(rare), tsec(spc), tsec(stp), tsec(brt), tsec(grd), sum(tsec, (scr, rare, spc, stp, brt, grd)),
         sum(tsec, (scr, rare, spc, stp, brt, grd)) / 3600)
    isfile(gcache("emulator_run.out")) &&
        say("process times of the emulator mode (/usr/bin/time -p): ", join(filter(l -> occursin(r"^(real|user|sys) ", l), readlines(gcache("emulator_run.out"))), ", "))
    R = [scr[S.id[i]] for i in run]
    for (nm, d, ids) in (("screen", scr, S.id[run]), (":rare_alleles (screen's inputs)", rare, S.id[run]), ("default :spectrum", spc, collect(keys(spc))))
        st = [Int(d[id][:status]) for id in ids]
        sayf("  %-32s defined %3d of %3d; no k50 %d (k50 = 0: %d, ∞: %d); saturated %d; ∂h/∂log k ≥ 0: %d; step 1e-3 used after saturation at 1e-2: %d; evaluations for k50 median %d",
             nm, count(==(0), st), length(st), count(==(1), st), count(id -> d[id][:k50] == 0, ids), count(id -> d[id][:k50] == Inf, ids),
             count(==(2), st), count(==(3), st), count(id -> d[id][:status] == 0 && d[id][:step] < GSTEP / 2, ids),
             round(Int, median([d[id][:evals] for id in ids])))
    end
    ok = [r[:status] == 0 for r in R]
    sayf("  screen at k50 (defined inputs): |h| median %.3f, max %.3f; ∂h/∂log k median %.1f (range %.1f to %.1f)",
         median(abs(r[:h0]) for r in R[ok]), maximum(abs(r[:h0]) for r in R[ok]), median(r[:hk] for r in R[ok]),
         minimum(r[:hk] for r in R[ok]), maximum(r[:hk] for r in R[ok]))

    say("\n== fallback of the beta solver (screen) ==")
    nss = count(r -> r[:status] == 0 && r[:sstat] > 0, R); nsg = count(r -> r[:status] == 0 && r[:sstat] == 0 && r[:sgen] > 0, R)
    say("  inputs whose run at k50 fell back: $(count(r -> r[:fb0] >= 0, R)) of $(length(R)); with ≥ 1 of the 10 difference points falling back: ",
        count(r -> r[:nfb] > 0, R), "; all 10: ", count(r -> r[:nfb] == 10, R))
    say("  defined inputs with a pair straddling a switch (one point fell back, the other not): $nss; with only shifts (both fell back, different generations): $nsg")
    Eas = fill(NaN, length(run), 4); Eco = fill(NaN, length(run), 4); Era = fill(NaN, length(run), 4)
    nre = 0; nlost = 0
    for (a, i) in enumerate(run)
        id = S.id[i]; r = R[a]
        Eas[a, :] .= g_e(r); Era[a, :] .= g_e(rare[id])
        if r[:status] == 0 && r[:sstat] > 0
            if haskey(grd, id) && grd[id][:status] == 0
                Eco[a, :] .= g_e(grd[id]); nre += 1
            else
                nlost += 1
            end
        else
            Eco[a, :] .= Eas[a, :]
        end
    end
    say("  straddling inputs recomputed with the grid and used in the consistent screen values: $nre; straddling inputs left undefined (beyond the cap of $(GREGRID_STATUS), or the grid undefined): $nlost")
    for (a, i) in enumerate(run)
        id = S.id[i]; haskey(grd, id) || continue
        g = grd[id]; r = R[a]
        sayf("    %-34s %-8s beta e %s   grid e %s   max |diff| %.3f%s", id, r[:sstat] > 0 ? "straddle" : "shift",
             join([@sprintf("% .3f", x) for x in g_e(r)], " "), join([@sprintf("% .3f", x) for x in g_e(g)], " "),
             maximum(abs.(g_e(r) .- g_e(g))), g[:resolved] == 1 ? "  (grid at its own k50)" : "")
    end

    say("\n== numerical checks of the screen (round 1) ==")
    d1 = Float64[]
    for id in sort(collect(keys(stp)))
        a, b = g_e(scr[id]), g_e(stp[id])
        append!(d1, abs.(a .- b))
        sayf("  N1 %-34s δ = 1e-2: %s   δ = 1e-3: %s", id, join([@sprintf("% .3f", x) for x in a], " "), join([@sprintf("% .3f", x) for x in b], " "))
    end
    d1 = fin(d1)
    n1 = !isempty(d1) && median(d1) <= 0.02 && quantile(d1, 0.9) <= 0.05
    isempty(d1) || sayf("  N1 step: |e(1e-2) − e(1e-3)| over %d values: median %.4f (≤ 0.02), 90th percentile %.4f (≤ 0.05), max %.4f: %s", length(d1),
                        median(d1), quantile(d1, 0.9), maximum(d1), n1 ? "stable" : "NOT STABLE")
    d2 = Float64[]
    for id in sort(collect(keys(brt)))
        a = g_e(scr[id]); b = [brt[id][:b1], brt[id][:b2], brt[id][:b3], brt[id][:b4]]
        append!(d2, abs.(a .- b))
        sayf("  N2 %-34s implicit: %s   brute force: %s", id, join([@sprintf("% .3f", x) for x in a], " "), join([@sprintf("% .3f", x) for x in b], " "))
    end
    nd2 = length(d2); d2 = fin(d2)
    n2 = !isempty(d2) && length(d2) == nd2 && median(d2) <= 0.03 && maximum(d2) <= 0.10
    isempty(d2) || sayf("  N2 brute force: |brute − implicit| over %d values (%d undefined): median %.4f (≤ 0.03), max %.4f (≤ 0.10): %s", length(d2), nd2 - length(d2),
                        median(d2), maximum(d2), n2 ? "verified" : "NOT VERIFIED")
    say("  numerical status of the verdicts below: ", n1 && n2 ? "reliable (N1 and N2 pass)" : "NUMERICALLY UNRELIABLE (N1 or N2 failed)")

    say("\n== primary comparison: 48 cells (locus number × input × information tercile of the input) ==")
    g_gstats("the screen, consistent values (grid where a difference straddles a fallback switch)", Eco, S, run, cells; verdict = true)
    g_gstats("the screen as run (beta values for every input)", Eas, S, run, cells; verdict = false)
    g_gstats(":rare_alleles on the screen's inputs", Era, S, run, cells; verdict = false)

    say("\n== S1: the literal per-stratum reading (32 cells; no verdict; distorted by the non-constancy, see the criteria) ==")
    A, _ = readdlm(gcache("table_elasticities.tsv"), '\t'; header = true)
    for (nm, E) in (("screen (consistent)", Eco), (":rare_alleles", Era))
        te = Float64[]; em = Float64[]; g1 = Bool[]; g2 = Bool[]
        say("  ", nm)
        for s in eachindex(STRATA)
            ins = [a for (a, i) in enumerate(run) if S.strat[i] == s]
            line = ""
            for j in 1:4
                r = findfirst(q -> Int(A[q, 1]) == s && A[q, 2] == GLAB[j], axes(A, 1))
                e, lo, hi = Float64(A[r, 3]), Float64(A[r, 5]), Float64(A[r, 6])
                v = fin(E[ins, j]); isempty(v) && continue
                m = median(v); push!(te, e); push!(em, m)
                (lo > 0 || hi < 0) && push!(g1, sign(m) == sign(e))
                tol = max(0.15, 0.3abs(e))
                abs(e) >= 0.1 && push!(g2, lo - tol <= m <= hi + tol)
                line *= @sprintf("  %s table % .3f emu median % .3f mean % .3f |", GLAB[j], e, m, mean(v))
            end
            say("    ", rpad(STRATA[s][1], 18), line)
        end
        sayf("    for information: sign %d/%d, size %d/%d, Spearman %.3f over %d cells", count(g1), length(g1), count(g2), length(g2), spearman(te, em), length(te))
    end

    say("\n== S2: input by input against the table's quadratic-model local elasticity (no verdict; the quadratic is poor for ω) ==")
    for (nm, E) in (("screen (consistent)", Eco), (":rare_alleles", Era))
        for j in 1:4
            k = [a for a in eachindex(run) if isfinite(E[a, j])]
            q = S.q[run[k], j]
            sayf("  %-20s e_%-6s n %3d  Spearman %.3f  median |emu − table| %.3f  median (emu − table) % .3f   emu quartiles % .3f % .3f % .3f   table % .3f % .3f % .3f",
                 nm, GLAB[j], length(k), spearman(E[k, j], q), median(abs.(E[k, j] .- q)), median(E[k, j] .- q),
                 quantile(E[k, j], 0.25), median(E[k, j]), quantile(E[k, j], 0.75), quantile(q, 0.25), median(q), quantile(q, 0.75))
        end
    end

    say("\n== S3: default :spectrum beside the screen (round 1, one input per stratum; no verdict) ==")
    for id in S.id[S.order.==1]
        haskey(spc, id) || continue
        a = findfirst(i -> S.id[i] == id, run)
        sayf("  %-34s screen %s   default :spectrum %s   :rare_alleles %s   k50 %.4g / %.4g / %.4g", id,
             join([@sprintf("% .3f", x) for x in Eco[a, :]], " "), join([@sprintf("% .3f", x) for x in g_e(spc[id])], " "),
             join([@sprintf("% .3f", x) for x in Era[a, :]], " "), scr[id][:k50], spc[id][:k50], rare[id][:k50])
    end
    for j in 1:4
        d = fin([g_e(spc[id])[j] - Eco[findfirst(i -> S.id[i] == id, run), j] for id in keys(spc)])
        isempty(d) || sayf("  e_%-6s default :spectrum − screen: median % .3f, median |diff| %.3f, max |diff| %.3f (n %d)", GLAB[j], median(d), median(abs.(d)), maximum(abs.(d)), length(d))
    end

    say("\n== S4: levels (no verdict) ==")
    say("  k50 of the emulator / the table's boundary rate at the same inputs (quadratic fit): median [quartiles] by locus number")
    for (nm, d) in (("screen", scr), (":rare_alleles", rare))
        cellsL = map(1:4) do li
            v = [d[S.id[i]][:k50] / exp(S.lkq[i]) for i in run if (S.strat[i] + 1) ÷ 2 == li && 0 < d[S.id[i]][:k50] < Inf]
            @sprintf("L %d: %.2f [%.2f, %.2f]", GLS[li], median(v), quantile(v, 0.25), quantile(v, 0.75))
        end
        say("    ", rpad(nm, 16), join(cellsL, "   "))
    end
    say("  least squares of the screen's log k50 on the logs of the inputs (per locus number, u indicator): slope (SE); beside the")
    say("  per-stratum linear table elasticities (u = 3, u = 5). A global slope from levels, not a gradient.")
    for li in 1:4
        ins = [i for i in run if (S.strat[i] + 1) ÷ 2 == li && 0 < scr[S.id[i]][:k50] < Inf]
        X = [ones(length(ins)) Float64.(iseven.(S.strat[ins])) log.(S.VP[ins]) log.(S.h2[ins]) log.(S.m[ins]) log.(S.omega[ins])]
        y = [log(scr[S.id[i]][:k50]) for i in ins]
        b = X \ y; res = y - X * b
        se = sqrt.(diag(inv(X' * X)) .* (sum(abs2, res) / max(length(y) - size(X, 2), 1)))
        tb(j, s) = Float64(A[findfirst(q -> Int(A[q, 1]) == s && A[q, 2] == GLAB[j], axes(A, 1)), 3])
        say("    L ", GLS[li], " (n ", length(ins), "): ", join([@sprintf("%s % .3f (%.3f) | table % .3f, % .3f", GLAB[j], b[2+j], se[2+j], tb(j, 2li - 1), tb(j, 2li)) for j in 1:4], "   "))
    end
    say("\n== post hoc (added after the results were seen; not part of the criteria) ==")
    sat = [i for (a, i) in enumerate(run) if R[a][:status] == 2]; def = [i for (a, i) in enumerate(run) if R[a][:status] == 0]
    sayf("  the screen's %d saturated (undefined) inputs have ω in [%.3g, %.3g] (median %.3g); the %d defined ones ω in [%.3g, %.3g] (median %.3g);",
         length(sat), minimum(S.omega[sat]), maximum(S.omega[sat]), median(S.omega[sat]), length(def), minimum(S.omega[def]), maximum(S.omega[def]), median(S.omega[def]))
    sayf("  of the %d screen inputs with ω < 2.2, %d are undefined: the undefined inputs are not a random subset, and the four low-ω cells are all excluded.",
         count(i -> S.omega[i] < 2.2, run), count(i -> S.omega[i] < 2.2, sat))
    sayf("  cells entering G1–G3 for the screen: 31 of 48. |∂h/∂log k| of the defined inputs: median %.0f, so P2000 goes from 0.9 to 0.1 over a factor e^(%.4f) in k.",
         median(abs(r[:hk]) for r in R[ok]), (log(-log(0.1)) - log(-log(0.9))) / median(abs(r[:hk]) for r in R[ok]))
    close(LOG[])
end

GMODE == "analyse" && analyse_mode()
