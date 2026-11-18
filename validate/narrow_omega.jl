# Diagnostic study (exploratory, no verdicts): why is the emulator's elasticity of the
# half-survival rate k50 in ω = InitFitnessSD about three times the published table's in the
# lowest tercile of ω (validate/elasticity_b.log: emulator −0.43 to −0.63, table −0.14 to −0.19)?
# Which side is off: the table's summary (a leverage-weighted logistic slope over ω from 0.003
# up; hypothesis T), or the emulator under strong selection (hypothesis E)?
# Running log with the hypotheses and predictions stated before each test:
# validate/narrow_omega.log; findings: validate/narrow_omega.md. No change to src/.
#
# Run from the package root (vignettes environment: XLSX, Catlab), at most 6 threads:
#   julia -t 6 --project=vignettes validate/narrow_omega.jl table     # the table alone (T), phase B's elasticities re-binned
#   julia -t 6 --project=vignettes validate/narrow_omega.jl curve     # the emulator's k50(ω) at reference inputs
#   julia -t 6 --project=vignettes validate/narrow_omega.jl vgmove    # emulator V_G and W̄ at the move against the table's columns
#   julia -t 6 --project=vignettes validate/narrow_omega.jl select    # the inputs of the kernel's local elasticities
#   julia -t 6 --project=vignettes validate/narrow_omega.jl pilot     # kernel: coarse rate grid at each input's own ω
#   julia -t 6 --project=vignettes validate/narrow_omega.jl kernel    # kernel: rate grid at ω·e^{±δ}
#   julia -t 6 --project=vignettes validate/narrow_omega.jl emfd      # emulators: k50 at ω·e^{±δ} (finite difference) and local
#   julia -t 6 --project=vignettes validate/narrow_omega.jl local     # kernel against emulators, per input
#   julia -t 6 --project=vignettes validate/narrow_omega.jl traj      # kernel and emulator trajectories at ω·e^{±δ}, one rate
#   julia -t 6 --project=vignettes validate/narrow_omega.jl links     # which link carries the difference
#   julia -t 6 --project=vignettes validate/narrow_omega.jl blocks    # block-lineage model at ω·e^{±δ} (some inputs)
# Each mode writes validate/narrow_omega_raw/out_<mode>.txt (copied into the log by hand, after the
# hypotheses of its step) and resumable caches in validate/narrow_omega_raw/ (git-ignored).
# Needs validate/elasticity_raw/inputs.tsv, table_cells.tsv, sample.tsv, b_screen.csv
# (validate/elasticity.jl table; validate/elasticity_b.jl run).

let src = read(joinpath(@__DIR__, "emulator.jl"), String)
    # emulator.jl's helpers (logistic, say, STRATA, SIMDIR, MASTER, load_master, …) without its main
    include_string(Main, replace(src, r"\nmain\(\)\s*$" => "\n"), joinpath(@__DIR__, "emulator.jl"))
end
using DelimitedFiles, Serialization
BLAS.set_num_threads(1)

const NMODE = isempty(ARGS) ? "table" : ARGS[1]
const NRAW = joinpath(@__DIR__, "narrow_omega_raw")
const ERAW = joinpath(@__DIR__, "elasticity_raw")
ncache(name) = (mkpath(NRAW); joinpath(NRAW, name))
const NLS = (100, 500, 1000, 5000)
# bins of ω, fixed before any fit of this study (finer where ω is small)
const NEDGES = [0.0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.4, 1.8, 2.2, 2.7, 3.3, 4.0, 5.0, 6.5, 8.0, 10.01]
const NREF = (VP = 0.0223, h2 = 0.315, m = 0.004)         # reference inputs: about the table's medians
const NSCREEN = (; ld = :dynamic, solver = :beta)
const NDELTA = 0.2                                          # arms ω·e^{±δ}
const MOVE, KBURN, HOR = 10_000, 1000, 2000
const KSEED = 50261120000
const NSELSEED = 20261120
cpu_seconds() = ccall(:clock, Clong, ()) / 1e6
function finish_mode(t0)
    open(ncache("cpu.txt"), "a") do io
        @printf(io, "%s\t%s\t%.0f\t%.0f\t%d\n", Libc.strftime("%Y-%m-%d %H:%M", time()), NMODE, cpu_seconds(), time() - t0, Threads.nthreads())
    end
    sayf("\n[mode %s: process CPU %.0f s (clock(): user + system, all threads, with compilation), wall %.0f s, %d threads]", NMODE,
         cpu_seconds(), time() - t0, Threads.nthreads())
end

# ---------------------------------------------------------------------------------------------
# data: one row per input (validate/elasticity_raw/inputs.tsv), plus the table's columns at 10000
# ---------------------------------------------------------------------------------------------
function n_data()
    f = joinpath(ERAW, "inputs.tsv")
    isfile(f) || error("run validate/elasticity.jl table first (needs $f)")
    A, _ = readdlm(f, '\t'; header = true)
    col(j) = Float64.(A[:, j])
    D = (; id = String.(A[:, 1]), strat = Int.(A[:, 2]), k = col(3), VP = col(4), h2 = col(5), m = col(6), omega = col(7),
         T = Int.(A[:, 10]), row = Int.(A[:, 12]))
    g = ncache("table_extra.tsv")
    if !isfile(g)
        tab = load_master(MASTER)
        open(g, "w") do io
            println(io, "gv\tns\tpm")
            for r in D.row
                println(io, tab.gv[r], '\t', tab.ns[r], '\t', tab.pm[r])
            end
        end
    end
    B, _ = readdlm(g, '\t'; header = true)
    merge(D, (; gv = Float64.(B[:, 1]), ns = Float64.(B[:, 2]), pm = Float64.(B[:, 3]), L = [STRATA[s][2] for s in D.strat],
              u5 = Float64.([STRATA[s][3] == 5 for s in D.strat]), y = Float64.(D.T .> 2000)))
end
nbin(ω) = searchsortedlast(NEDGES, ω)
binlab(b) = @sprintf("[%.4g, %.4g)", max(NEDGES[b], 0.003), min(NEDGES[b+1], 10.0))

# ---- model A: per locus number, linear in log k, log V_P, log h², log m, u = 5, and a
#      piecewise-linear function of log ω with knots at the bin edges
const KNOTS = log.(NEDGES[2:end-1])
spline_cols(x) = hcat(x, [max.(x .- κ, 0.0) for κ in KNOTS]...)
designA(D, ix) = [ones(length(ix)) log.(D.k[ix]) log.(D.VP[ix]) log.(D.h2[ix]) log.(D.m[ix]) D.u5[ix] spline_cols(log.(D.omega[ix]))]
function fitA(D, L)
    ix = findall(==(L), D.L)
    X = designA(D, ix); y = D.y[ix]
    β, conv = logistic(X, y)
    μ = @. 1 / (1 + exp(-$(X * β)))
    V = inv(Symmetric(X' * ((μ .* (1 .- μ)) .* X)))
    (; L, ix, β, conv, V, μ)
end
"Model A's local e_ω in bin b (and its delta-method SE)."
function slopeA(f, b)
    g = zeros(length(f.β)); g[7:6+b] .= 1.0
    s = dot(g, f.β); βk = f.β[2]
    e = -s / βk
    gr = -g ./ βk; gr[2] = s / βk^2
    e, sqrt(max(gr' * f.V * gr, 0.0))
end
"Model A's log k50 at inputs (VP, h2, m, ω), u = 3."
function lk50A(f, VP, h2, m, ω; u5 = 0.0)
    x = [1.0, 0.0, log(VP), log(h2), log(m), u5, vec(spline_cols([log(ω)]))...]
    -dot(x, f.β) / f.β[2]
end
etaA(f, k, VP, h2, m, ω, u5) = dot([1.0, log(k), log(VP), log(h2), log(m), u5, vec(spline_cols([log(ω)]))...], f.β)

"A linear logistic fit on rows ix with ω as column 7: (β, ok, e_ω, se, n, alive)."
function fit_lin(D, ix)
    X = [ones(length(ix)) log.(D.k[ix]) log.(D.VP[ix]) log.(D.h2[ix]) log.(D.m[ix]) D.u5[ix] log.(D.omega[ix])]
    y = D.y[ix]
    (length(ix) < 60 || sum(y) < 12) && return (; β = fill(NaN, 7), ok = false, e = NaN, se = NaN, n = length(ix), alive = Int(sum(y)), X, μ = Float64[])
    β, conv = logistic(X, y)
    ok = conv && β[2] < 0 && all(isfinite, β) && maximum(abs, β) < 1e3
    ok || return (; β, ok, e = NaN, se = NaN, n = length(ix), alive = Int(sum(y)), X, μ = Float64[])
    μ = @. 1 / (1 + exp(-$(X * β)))
    V = inv(Symmetric(X' * ((μ .* (1 .- μ)) .* X)))
    g = zeros(7); g[7] = -1 / β[2]; g[2] = β[7] / β[2]^2
    (; β, ok, e = -β[7] / β[2], se = sqrt(max(g' * V * g, 0.0)), n = length(ix), alive = Int(sum(y)), X, μ)
end

function read_csv(f)
    L = readlines(f)
    h = split(L[1], ",")
    [NamedTuple{Tuple(Symbol.(h))}(Tuple(split(l, ","))) for l in L[2:end] if !isempty(l)]
end
pf64(s) = something(tryparse(Float64, s), NaN)
med(v) = (u = filter(isfinite, v); isempty(u) ? NaN : median(u))

function table_mode()
    D = n_data()
    say("== table mode: the published table alone (119,984 rows), and phase B's emulator elasticities re-binned ==")
    say(Libc.strftime("%Y-%m-%d %H:%M", time()), "; alive at 2000 = 'Generations to extinction' > 2000")

    say("\n-- T0: deaths that do not depend on the rate --")
    sayf("  rows with T ≤ 0 (death before the move): %d; T = 1: %d; T ≤ 5: %d; smallest T: %d", count(<=(0), D.T), count(==(1), D.T),
         count(<=(5), D.T), minimum(D.T))
    say("  by ω bin (all strata pooled): n, alive at 2000, T ≤ 5, T ≤ 5 among rates k < 0.01, alive among k < 0.003, median")
    say("  'P Survived at 10000' / 1000 (survivors of viability selection in the last generation before the move) and its 5% quantile,")
    say("  median Genotype Var at 10000 / (h²V_P)")
    for b in 1:length(NEDGES)-1
        ix = findall(i -> nbin(D.omega[i]) == b, eachindex(D.omega))
        lowk = filter(i -> D.k[i] < 0.01, ix); vlow = filter(i -> D.k[i] < 0.003, ix)
        sayf("    ω %-16s n %6d  alive %5d (%.3f)  T≤5 %4d  T≤5|k<0.01 %3d/%-4d  alive|k<0.003 %3d/%-3d  ns/N median %.3f q05 %.3f  gv/(h²V_P) %.3f",
             binlab(b), length(ix), count(i -> D.y[i] > 0, ix), mean(D.y[ix]), count(i -> D.T[i] <= 5, ix), count(i -> D.T[i] <= 5, lowk),
             length(lowk), count(i -> D.y[i] > 0, vlow), length(vlow), median(D.ns[ix]) / 1000, quantile(D.ns[ix], 0.05) / 1000,
             median(D.gv[ix] ./ (D.h2[ix] .* D.VP[ix])))
    end
    say("  rows with ω < 0.05, in detail (ω, k, T, ns at 10000, gv at 10000), first 30 by ω:")
    for i in sort(findall(<(0.05), D.omega); by = i -> D.omega[i])[1:min(30, count(<(0.05), D.omega))]
        sayf("    %-34s ω %.4f  k %.4f  T %6d  ns %4.0f  gv %.3g", D.id[i], D.omega[i], D.k[i], D.T[i], D.ns[i], D.gv[i])
    end

    say("\n-- T1: e_ω by ω bin --")
    say("  model A (per locus number, strata u = 3, 5 pooled with an indicator): logit = β0 + βk log k + β log(V_P, h², m) + u5 +")
    say("  piecewise-linear in log ω with knots at the bin edges; local e_ω = −(slope in the bin)/βk (delta-method SE).")
    say("  model B: a separate linear model per bin (all coefficients free): its own e_ω (se) and log k50 at the reference inputs")
    sayf("  (V_P %.4g, h² %.3g, m %.3g, u = 3) and the bin's geometric-mean ω; then e between neighbouring bins from B's k50.", NREF.VP, NREF.h2, NREF.m)
    fits = Dict(L => fitA(D, L) for L in NLS)
    open(ncache("table_bins.tsv"), "w") do io
        println(io, "L\tbin\tlo\thi\tgm\tn\talive\teA\tseA\tlk50A\teB\tseB\tlk50B\tinfo")
        for L in NLS
            f = fits[L]
            sayf("  L %d: model A converged %s, βk %.3f", L, f.conv, f.β[2])
            prev = (NaN, NaN)
            for b in 1:length(NEDGES)-1
                ix = [i for i in f.ix if nbin(D.omega[i]) == b]
                gm = exp(mean(log.(D.omega[ix])))
                eA, sA = slopeA(f, b)
                lkA = lk50A(f, NREF.VP, NREF.h2, NREF.m, gm)
                fb = fit_lin(D, ix)
                lkB = fb.ok ? -(fb.β[1] + fb.β[3] * log(NREF.VP) + fb.β[4] * log(NREF.h2) + fb.β[5] * log(NREF.m) + fb.β[7] * log(gm)) / fb.β[2] : NaN
                between = isfinite(prev[1]) && isfinite(lkB) ? (lkB - prev[1]) / (log(gm) - prev[2]) : NaN
                info = sum(f.μ[j] * (1 - f.μ[j]) for (j, i) in enumerate(f.ix) if nbin(D.omega[i]) == b)
                sayf("    ω %-16s n %5d alive %4d  A: e % .3f (%.3f) k50 %.4f | B: e % .3f (%.3f) k50 %.4f βk % .2f  e from the previous bin's k50 % .3f | info %.0f",
                     binlab(b), length(ix), fb.alive, eA, sA, exp(lkA), fb.e, fb.se, exp(lkB), fb.β[2], between, info)
                println(io, join(Any[L, b, NEDGES[b], NEDGES[b+1], gm, length(ix), fb.alive, eA, sA, lkA, fb.e, fb.se, lkB, info], '\t'))
                isfinite(lkB) && (prev = (lkB, log(gm)))
            end
        end
    end
    serialize(ncache("fitsA.jls"), Dict(L => (; β = fits[L].β, V = fits[L].V) for L in NLS))

    say("\n-- T2: the tercile-1 fit of the criteria (phase B's table side), and the same model on sub-ranges of ω --")
    say("  (delta-method SE; 'lev' = share of the fit's leverage Σ w (x − x̄_w)², w = μ(1 − μ), x = log ω, carried by ω < 0.3 and by ω < 1)")
    cells = readdlm(joinpath(ERAW, "table_cells.tsv"), '\t'; header = true)[1]
    for L in NLS
        c = [cells[r, :] for r in axes(cells, 1) if cells[r, 1] == L && cells[r, 4] == "omega" && cells[r, 5] == 1][1]
        cut = Float64(c[7])
        base = findall(i -> D.L[i] == L && D.omega[i] <= cut, eachindex(D.L))
        sayf("  L %d, tercile 1 = ω ≤ %.3f; criteria value %.3f [%.3f, %.3f]", L, cut, c[8], c[9], c[10])
        for (lab, lo, hi) in (("whole tercile", 0.0, cut), ("ω < 0.3", 0.0, 0.3), ("ω < 1", 0.0, 1.0), ("0.3 ≤ ω", 0.3, cut), ("0.5 ≤ ω", 0.5, cut),
                              ("1 ≤ ω", 1.0, cut), ("0.3 ≤ ω < 1", 0.3, 1.0), ("0.5 ≤ ω < 1.5", 0.5, 1.5), ("1.5 ≤ ω", 1.5, cut))
            ix = [i for i in base if lo <= D.omega[i] < hi || (hi == cut && D.omega[i] == cut && lo <= cut)]
            f = fit_lin(D, ix)
            lev = ""
            if f.ok && lab == "whole tercile"
                w = f.μ .* (1 .- f.μ); x = log.(D.omega[ix]); xb = sum(w .* x) / sum(w); l = w .* (x .- xb) .^ 2
                lev = @sprintf("  lev(ω<0.3) %.2f lev(ω<1) %.2f; info share ω<0.3 %.2f, ω<1 %.2f", sum(l[D.omega[ix].<0.3]) / sum(l), sum(l[D.omega[ix].<1.0]) / sum(l),
                               sum(w[D.omega[ix].<0.3]) / sum(w), sum(w[D.omega[ix].<1.0]) / sum(w))
            end
            sayf("    %-16s n %5d alive %4d  e_ω % .3f (%.3f)  βk % .2f%s", lab, f.n, f.alive, f.e, f.se, f.β[2], lev)
        end
    end

    say("\n-- T3: phase B's emulator elasticities (two-stage screen, k50_gradient; validate/elasticity_raw/b_screen.csv) by the same ω bins --")
    say("  beside model A's local e_ω of the input's locus number and bin; then per bin: n, emulator median (mean), table A mean")
    S = read_csv(joinpath(ERAW, "b_screen.csv"))
    rowof = Dict(id => i for (i, id) in enumerate(D.id))
    em = [(; id = String(s.id), i = rowof[String(s.id)], e = pf64(s.e4), k50 = pf64(s.k50), ok = s.status == "0") for s in S]
    for b in 1:length(NEDGES)-1
        v = [x for x in em if x.ok && nbin(D.omega[x.i]) == b]
        isempty(v) && continue
        tA = [slopeA(fits[D.L[x.i]], b)[1] for x in v]
        sayf("    ω %-16s n %2d  emulator e_ω median % .3f mean % .3f [min % .3f, max % .3f]   table A (same L, bin) mean % .3f", binlab(b), length(v),
             median([x.e for x in v]), mean([x.e for x in v]), minimum(x.e for x in v), maximum(x.e for x in v), mean(tA))
    end
    say("  levels: emulator k50 / model A's k50 at the same inputs (own ω, V_P, h², m, u), by ω bin: n, median ratio, quartiles")
    for b in 1:length(NEDGES)-1
        v = [x for x in em if isfinite(x.k50) && x.k50 > 0 && nbin(D.omega[x.i]) == b]
        isempty(v) && continue
        r = [x.k50 / exp(lk50A(fits[D.L[x.i]], D.VP[x.i], D.h2[x.i], D.m[x.i], D.omega[x.i]; u5 = D.u5[x.i])) for x in v]
        sayf("    ω %-16s n %2d  ratio median %.2f  [%.2f, %.2f]", binlab(b), length(v), median(r), quantile(r, 0.25), quantile(r, 0.75))
    end
    say("  the low-ω tercile cell, emulator side re-weighted: inputs of the cell (ω ≤ cut), their e_ω averaged (i) plainly, (ii) by a")
    say("  weight that makes the sample mimic the tercile fit's leverage in ω: the fit's leverage per ω bin / the number of sample inputs in the bin")
    for L in NLS
        c = [cells[r, :] for r in axes(cells, 1) if cells[r, 1] == L && cells[r, 4] == "omega" && cells[r, 5] == 1][1]
        cut = Float64(c[7])
        ix = findall(i -> D.L[i] == L && D.omega[i] <= cut, eachindex(D.L))
        f = fit_lin(D, ix)
        w = f.μ .* (1 .- f.μ); x = log.(D.omega[ix]); xb = sum(w .* x) / sum(w); l = w .* (x .- xb) .^ 2
        levb = [sum(l[[nbin(D.omega[i]) == b for i in ix]]) for b in 1:length(NEDGES)-1] ./ sum(l)
        v = [z for z in em if z.ok && D.L[z.i] == L && D.omega[z.i] <= cut]
        nb = [count(z -> nbin(D.omega[z.i]) == b, v) for b in 1:length(NEDGES)-1]
        covered = sum(levb[nb.>0])
        ew = sum(levb[nbin(D.omega[z.i])] / nb[nbin(D.omega[z.i])] * z.e for z in v) / covered
        sayf("    L %4d: table %.3f | emulator n %2d median % .3f mean % .3f leverage-weighted mean % .3f (bins with sample inputs carry %.2f of the leverage)",
             L, c[8], length(v), median([z.e for z in v]), mean([z.e for z in v]), ew, covered)
    end
end

# ---------------------------------------------------------------------------------------------
# emulator: k50(ω) curves at the reference inputs
# ---------------------------------------------------------------------------------------------
template(L, u = 3) = Params(read_qtsim_options(joinpath(SIMDIR, @sprintf("QT%04d-Possion%.1f", L, u), "Input-00001.txt")); Seed = nothing)
refparams(L, ω) = Params(template(L); PhenotypeVar = NREF.VP, Heritability = NREF.h2, MutationScalar = NREF.m, InitFitnessSD = ω)
const CURVE_W = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.75, 1.0, 1.3, 1.7, 2.2, 2.8, 3.6, 5.0, 7.0, 9.0]
const CSVLOCK = ReentrantLock()
function append_line(f, header, line)
    lock(CSVLOCK) do
        new = !isfile(f)
        open(f, "a") do io
            new && println(io, header)
            println(io, line)
        end
    end
end
done_keys(f) = isfile(f) ? Set(String(first(split(l, ","))) for l in readlines(f)[2:end]) : Set{String}()

"k50 and the local e_ω of a deterministic closure by k50_gradient; (k50, e, status, evals, continued, seconds)."
function kgrad(p; loci = :spectrum, spectrum = NSCREEN, bracket = nothing)
    t = time()
    r = try
        loci === :spectrum ? k50_gradient(p; loci, spectrum, names = (:InitFitnessSD,), bracket, onswitch = :flag) :
        k50_gradient(p; loci, names = (:InitFitnessSD,), bracket)
    catch err
        return (; k50 = NaN, e = NaN, status = "error:" * first(split(sprint(showerror, err), '\n')), evals = 0, continued = false, sec = time() - t)
    end
    (; k50 = r.k50, e = r.elasticities.InitFitnessSD, status = String(r.status), evals = r.evals, continued = r.continued, sec = time() - t)
end

function curve_mode()
    f = ncache("curve.csv")
    have = done_keys(f)
    jobs = [(L, ω) for L in (100, 500, 1000) for ω in CURVE_W if !("$L|$ω" in have)]
    say("== curve mode: the emulator's k50(ω) at the reference inputs; $(length(jobs)) points to do ==")
    Threads.@threads :greedy for (L, ω) in jobs
        p = refparams(L, ω)
        s = kgrad(p)
        ra = try
            emulated_half_survival_rate(p; rtol = 1e-6).k50
        catch
            NaN
        end
        append_line(f, "key,L,omega,k50,e,status,evals,continued,sec,k50_rare",
                    join(Any["$L|$ω", L, ω, s.k50, s.e, s.status, s.evals, s.continued, round(s.sec; digits = 1), ra], ","))
    end
    R = read_csv(f)
    fits = deserialize(ncache("fitsA.jls"))
    D = n_data()
    cells = readdlm(joinpath(ERAW, "table_cells.tsv"), '\t'; header = true)[1]
    for L in (100, 500, 1000)
        fa = (; β = fits[L].β)
        rs = sort([r for r in R if r.L == string(L)]; by = r -> pf64(r.omega))
        ω = pf64.(getfield.(rs, :omega)); k = pf64.(getfield.(rs, :k50)); e = pf64.(getfield.(rs, :e)); kr = pf64.(getfield.(rs, :k50_rare))
        sayf("  L %d (V_P %.4g, h² %.3g, m %.3g, u = 3): ω, screen k50, local e_ω (status), e from neighbouring k50s, table model A k50 at the same inputs, ratio; :rare_alleles k50, ratio, its e from neighbours",
             L, NREF.VP, NREF.h2, NREF.m)
        for i in eachindex(ω)
            kt = exp(lk50A(fa, NREF.VP, NREF.h2, NREF.m, ω[i]))
            fd = 1 < i < length(ω) ? (log(k[i+1]) - log(k[i-1])) / (log(ω[i+1]) - log(ω[i-1])) : NaN
            fr = 1 < i < length(ω) ? (log(kr[i+1]) - log(kr[i-1])) / (log(ω[i+1]) - log(ω[i-1])) : NaN
            sayf("    ω %5.2f  screen k50 %.5f  e % .3f (%s%s)  e_fd % .3f | table A k50 %.5f  screen/table %.2f | rare k50 %.5f rare/table %.2f e_fd % .3f",
                 ω[i], k[i], e[i], rs[i].status, rs[i].continued == "true" ? ", continued" : "", fd, kt, k[i] / kt, kr[i], kr[i] / kt, fr)
        end
        # the tercile fit mimicked on the emulator's curve: weighted least-squares slope of the curve's log k50 (interpolated in log ω)
        # on log ω over the rows of the table's tercile 1, with the weights μ(1 − μ) of the table's tercile fit
        c = [cells[r, :] for r in axes(cells, 1) if cells[r, 1] == L && cells[r, 4] == "omega" && cells[r, 5] == 1][1]
        cut = Float64(c[7])
        ok = findall(i -> isfinite(k[i]) && k[i] > 0, eachindex(k))
        lw, lk = log.(ω[ok]), log.(k[ok])
        interp(x) = (x <= lw[1] ? lk[1] + (x - lw[1]) * (lk[2] - lk[1]) / (lw[2] - lw[1]) :
                     x >= lw[end] ? lk[end] : (j = searchsortedlast(lw, x); lk[j] + (x - lw[j]) * (lk[j+1] - lk[j]) / (lw[j+1] - lw[j])))
        for (lab, lo) in (("whole tercile (curve extrapolated below ω = $(ω[ok][1]))", 0.0), ("ω ≥ 0.05", 0.05), ("ω ≥ 0.3", 0.3), ("ω ≥ 1", 1.0))
            ix = findall(i -> D.L[i] == L && lo <= D.omega[i] <= cut, eachindex(D.L))
            ft = fit_lin(D, ix)
            ft.ok || continue
            w = ft.μ .* (1 .- ft.μ); x = log.(D.omega[ix]); xb = sum(w .* x) / sum(w)
            yv = interp.(x); yb = sum(w .* yv) / sum(w)
            slope = sum(w .* (x .- xb) .* (yv .- yb)) / sum(w .* (x .- xb) .^ 2)
            sayf("    tercile-style slope over %-58s: emulator curve % .3f | table's own fit % .3f (%.3f)", lab, slope, ft.e, ft.se)
        end
    end
end

# ---------------------------------------------------------------------------------------------
# emulator at the move against the table's columns at 10000 (one published run per row)
# ---------------------------------------------------------------------------------------------
"The emulator's state in the last generation before the move: V_G (as recorded), W̄, nsurv."
function emul_premove(p; loci = :spectrum, spectrum = NSCREEN, burnin = 2000)
    q = Params(p; Sample = 1, NoGenerations = burnin + 2, MoveFitnessSDGen = burnin, MoveFitnessSD = 1e-9)
    m = loci === :spectrum ? run_moments(q; VG = house_of_gauss_vg(p), erosion = true, loci, spectrum) :
        run_moments(q; VG = house_of_gauss_vg(p), erosion = true, loci)
    length(m.gen) < burnin && return (; VG = NaN, W = 0.0, VP = NaN)
    (; VG = m.geno_var[burnin], W = m.mean_fitness[burnin], VP = m.pheno_var[burnin])
end

function vgmove_mode()
    D = n_data()
    f = ncache("vgmove.csv")
    have = done_keys(f)
    rng = MersenneTwister(NSELSEED + 7)
    pick = Int[]
    for L in (100, 500, 1000), b in 1:11                      # ω < 4
        ix = findall(i -> D.L[i] == L && nbin(D.omega[i]) == b, eachindex(D.L))
        append!(pick, ix[randperm(rng, length(ix))[1:min(20, length(ix))]])
    end
    jobs = [i for i in pick if !(D.id[i] in have)]
    say("== vgmove mode: emulator state before the move vs the table's 'Genotype Var at 10000' and 'P Survived at 10000' ==")
    say("  20 rows per (L ∈ {100, 500, 1000}, ω bin below 4), drawn uniformly (MersenneTwister($(NSELSEED + 7))): $(length(pick)) rows; to do $(length(jobs))")
    Threads.@threads :greedy for i in jobs
        p = Params(read_qtsim_options(joinpath(SIMDIR, D.id[i])); Seed = nothing)
        s = emul_premove(p)
        n = emul_premove(p; spectrum = (;))
        r = emul_premove(p; loci = :rare_alleles)
        append_line(f, "id,L,omega,gv,ns,VG_screen,W_screen,VG_noLD,W_noLD,VG_rare,W_rare",
                    join(Any[D.id[i], D.L[i], D.omega[i], D.gv[i], D.ns[i], s.VG, s.W, n.VG, n.W, r.VG, r.W], ","))
    end
    R = read_csv(f)
    say("  per ω bin (three locus numbers pooled, then each): median of emulator V_G / table V_G (one run per row, so the table side is noisy;")
    say("  medians over 20 rows), for the screen (dynamic LD, beta), default :spectrum (no LD, grid), :rare_alleles; and median emulator W̄ / (ns/1000)")
    for Lsel in (0, 100, 500, 1000)
        say(Lsel == 0 ? "  all three L" : "  L $Lsel")
        for b in 1:11
            v = [r for r in R if nbin(pf64(r.omega)) == b && (Lsel == 0 || r.L == string(Lsel))]
            isempty(v) && continue
            g(a) = med([pf64(getfield(r, a)) / pf64(r.gv) for r in v])
            w(a) = med([pf64(getfield(r, a)) / (pf64(r.ns) / 1000) for r in v])
            sayf("    ω %-16s n %3d  V_G ratio: screen %.3f  no-LD %.3f  rare %.3f | W̄ ratio: screen %.3f no-LD %.3f rare %.3f | table median V_G %.4g, ns/N %.3f",
                 binlab(b), length(v), g(:VG_screen), g(:VG_noLD), g(:VG_rare), w(:W_screen), w(:W_noLD), w(:W_rare), med([pf64(r.gv) for r in v]),
                 med([pf64(r.ns) / 1000 for r in v]))
        end
    end
    say("  slopes of log V_G at the move on log ω over 0.3 ≤ ω < 2.2 (least squares over the sampled rows, with log V_P, log h², log m as covariates):")
    v = [r for r in R if 0.3 <= pf64(r.omega) < 2.2]
    rowof = Dict(id => i for (i, id) in enumerate(D.id))
    X = reduce(vcat, [[1.0 log(pf64(r.omega)) log(D.VP[rowof[r.id]]) log(D.h2[rowof[r.id]]) log(D.m[rowof[r.id]]) (r.L == "500") (r.L == "1000")] for r in v])
    for a in (:gv, :VG_screen, :VG_noLD, :VG_rare)
        y = [log(pf64(getfield(r, a))) for r in v]
        ok = isfinite.(y)
        β = X[ok, :] \ y[ok]
        res = y[ok] .- X[ok, :] * β
        se = sqrt(sum(abs2, res) / (count(ok) - size(X, 2)) * inv(X[ok, :]' * X[ok, :])[2, 2])
        sayf("    %-10s d log V_G / d log ω = % .3f (se %.3f), n %d", a, β[2], se, count(ok))
    end
end

# ---------------------------------------------------------------------------------------------
# the kernel's local elasticity
# ---------------------------------------------------------------------------------------------
const SEL_RANGES = [(0.3, 0.6), (0.6, 1.0), (1.0, 1.5), (1.5, 2.2)]
function select_mode()
    D = n_data()
    fits = deserialize(ncache("fitsA.jls"))
    used = Set(String(s.id) for s in read_csv(joinpath(ERAW, "b_screen.csv")))
    union!(used, ["QT0100-Possion3.0/Input-14211.txt", "QT0100-Possion5.0/Input-01948.txt", "QT0100-Possion5.0/Input-09336.txt"])
    say("== select mode: inputs for the kernel's local elasticities ==")
    say("  one per (L ∈ {100, 500, 1000}, ω range ∈ $(SEL_RANGES)), from both recombination strata, drawn with probability ∝ μ(1 − μ) of model A")
    say("  (rows near the table's boundary; Efraimidis–Spirakis, MersenneTwister($NSELSEED + 10·(index of L) + range)); phase B's 128 inputs and")
    say("  its three development inputs excluded. The file's own rate is not used below (rates are scanned).")
    open(ncache("inputs.csv"), "w") do io
        println(io, "j,id,L,u,omega,VP,h2,m,kfile,T,k50_tab")
        j = 0
        for (li, L) in enumerate((100, 500, 1000)), (ri, (lo, hi)) in enumerate(SEL_RANGES)
            ix = findall(i -> D.L[i] == L && lo <= D.omega[i] < hi && !(D.id[i] in used), eachindex(D.L))
            fa = (; β = fits[L].β)
            μ = [1 / (1 + exp(-etaA(fa, D.k[i], D.VP[i], D.h2[i], D.m[i], D.omega[i], D.u5[i]))) for i in ix]
            rng = MersenneTwister(NSELSEED + 10li + ri)
            key = [rand(rng)^(1 / max(m * (1 - m), 1e-300)) for m in μ]
            i = ix[argmax(key)]
            j += 1
            kt = exp(lk50A(fa, D.VP[i], D.h2[i], D.m[i], D.omega[i]; u5 = D.u5[i]))
            println(io, join(Any[j, D.id[i], L, STRATA[D.strat[i]][3], D.omega[i], D.VP[i], D.h2[i], D.m[i], D.k[i], D.T[i], kt], ","))
            sayf("    %2d %-34s ω %.3f V_P %.4f h² %.3f m %.5f | file k %.4f T %6d | table model A k50 %.4f", j, D.id[i], D.omega[i], D.VP[i], D.h2[i],
                 D.m[i], D.k[i], D.T[i], kt)
        end
    end
end
kinputs() = [(; j = parse(Int, r.j), id = String(r.id), L = parse(Int, r.L), omega = pf64(r.omega), VP = pf64(r.VP), h2 = pf64(r.h2), m = pf64(r.m),
                k50_tab = pf64(r.k50_tab)) for r in read_csv(ncache("inputs.csv"))]
readp(id) = Params(read_qtsim_options(joinpath(SIMDIR, id)); Seed = nothing)
armω(ω, arm) = ω * exp(arm * NDELTA)

function genic_var(w)
    G = genomes(w)
    L, P, N = size(G)
    s = 0.0
    @inbounds for l in 1:L
        s += P * var(view(G, l, :, :))
    end
    s
end

"""
One kernel run of the published protocol from an equilibrium start (`equilibrium_world`, burn-in
KBURN generations before the move at generation 10000; alive at 2000 = the run completes generation
11999, the table's T > 2000). The loop is `generation!` step by step (as validate/spectrum_diagnosis.jl).
`traj = true` also records the trajectory.
"""
function krun(p, seed; traj = false, every = 10)
    q = Params(p; NoGenerations = MOVE + HOR, Sample = 1, Seed = nothing)
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    VGm = NaN; genm = NaN; Wpre = 0.0; npre = 0; extinct_gen = -1
    gen = Int[]; VG = Float64[]; VP = Float64[]; zbar = Float64[]; θ = Float64[]; Wb = Float64[]; nsurv = Int[]; sg = Float64[]
    xg = Int[]; genic = Float64[]
    for it in w.gen:q.NoGenerations-1
        w.gen = it
        it == MOVE - 1 && (genm = genic_var(w))
        if traj && it >= MOVE - 200 && it % every == 0
            push!(xg, it); push!(genic, genic_var(w))
        end
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        N = w.N
        it == MOVE - 1 && (VGm = var(view(w.genotype, 1:N)))
        if MOVE - 100 <= it < MOVE
            Wpre += w.nsurv / N; npre += 1
        end
        if traj && it >= MOVE - 200
            g = mean(view(w.genotype, 1:N))
            push!(gen, it); push!(VG, var(view(w.genotype, 1:N))); push!(VP, var(view(w.phenotype, 1:N)))
            push!(zbar, mean(view(w.phenotype, 1:N))); push!(θ, w.θ); push!(Wb, mean(view(w.fitness, 1:N))); push!(nsurv, w.nsurv)
            push!(sg, w.nsurv > 0 ? mean(w.genotype[i] for i in w.survivors) - g : NaN)
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    base = (; seed, extinct_gen, VGm, genm, Wpre = npre > 0 ? Wpre / npre : NaN)
    traj ? merge(base, (; gen, VG, VP, zbar, θ, Wbar = Wb, nsurv, sg, xg, genic)) : base
end

const KHEAD = "key,j,stage,arm,ri,r,k,omega,seed,extinct_gen,VGm,genm,Wpre"
function run_cells(cells)                # cells: (j, id, stage, arm, ri, r, k)
    f = ncache("kernel.csv")
    have = done_keys(f)
    jobs = [c for c in cells if !("$(c.j)|$(c.stage)|$(c.arm)|$(c.ri)|$(c.r)" in have)]
    sort!(jobs; by = c -> -c.L)
    say("  $(length(jobs)) kernel runs to do ($(length(cells) - length(jobs)) cached)")
    P = Dict(c.id => readp(c.id) for c in jobs)
    t0 = time(); done = Threads.Atomic{Int}(0)
    Threads.@threads :greedy for c in jobs
        p0 = P[c.id]
        ω = armω(p0.InitFitnessSD, c.arm)
        p = Params(p0; InitFitnessSD = ω, MoveFitnessSD = c.k)
        # the seed does not depend on the arm: the two arms of a (rate, replicate) share their random numbers
        seed = KSEED + 10_000_000 * c.j + 100_000 * (c.stage == "pilot" ? 0 : c.stage == "main" ? 1 : c.stage == "ext" ? 2 : 3) + 1000 * c.ri + c.r
        x = krun(p, seed)
        append_line(f, KHEAD, join(Any["$(c.j)|$(c.stage)|$(c.arm)|$(c.ri)|$(c.r)", c.j, c.stage, c.arm, c.ri, c.r, c.k, ω, seed, x.extinct_gen, x.VGm,
                                       x.genm, x.Wpre], ","))
        d = Threads.atomic_add!(done, 1) + 1
        d % 200 == 0 && @printf("  %d/%d  %.0f s\n", d, length(jobs), time() - t0)
    end
end
kruns() = [(; j = parse(Int, r.j), stage = String(r.stage), arm = parse(Int, r.arm), ri = parse(Int, r.ri), r = parse(Int, r.r), k = pf64(r.k),
              omega = pf64(r.omega), alive = parse(Int, r.extinct_gen) < 0, T = parse(Int, r.extinct_gen) < 0 ? 240000 : parse(Int, r.extinct_gen) - MOVE + 1,
              VGm = pf64(r.VGm), genm = pf64(r.genm), Wpre = pf64(r.Wpre)) for r in read_csv(ncache("kernel.csv"))]

"Pool-adjacent-violators for a non-increasing sequence with weights."
function pava_dec(p, n)
    v = collect(Float64, p); w = collect(Float64, n); idx = [[i] for i in eachindex(v)]
    i = 1
    while i < length(v)
        if v[i] < v[i+1] - 1e-12
            tot = w[i] + w[i+1]
            v[i] = (v[i] * w[i] + v[i+1] * w[i+1]) / tot; w[i] = tot
            append!(idx[i], idx[i+1]); deleteat!(v, i + 1); deleteat!(w, i + 1); deleteat!(idx, i + 1)
            i > 1 && (i -= 1)
        else
            i += 1
        end
    end
    out = similar(p, Float64)
    for (a, ix) in enumerate(idx), i in ix
        out[i] = v[a]
    end
    out
end
"""
Spearman–Kärber estimate of the mean of log k50 from fractions alive `f` (of `n` runs) at the
ascending log rates `x`: the fractions made non-increasing (PAVA), 1 assumed one grid step below
the lowest rate and 0 one step above the highest. Also returns whether the grid covers the
transition (fraction ≥ 0.8 at the lowest rate and ≤ 0.2 at the highest).
"""
function spearman_karber(x, f, n)
    p = pava_dec(f, n)
    d1 = x[2] - x[1]; d2 = x[end] - x[end-1]
    xs = [x[1] - d1; x; x[end] + d2]; ps = [1.0; p; 0.0]
    m = sum((ps[i] - ps[i+1]) * (xs[i] + xs[i+1]) / 2 for i in 1:length(xs)-1)
    (; m, covered = p[1] >= 0.8 && p[end] <= 0.2)
end

const PILOT_STEPS = -0.75:0.15:0.75
const PILOT_R = 6
const MAIN_STEPS = (-3.5:1.0:3.5) .* 0.08
const MAIN_R = 30
function pilot_mode()
    ins = kinputs()
    say("== pilot mode: kernel at each input's own ω; rates = table model A's k50 · exp($(PILOT_STEPS)), $PILOT_R runs each ==")
    cells = [(; j = x.j, id = x.id, L = x.L, stage = "pilot", arm = 0, ri, r, k = x.k50_tab * exp(s)) for x in ins for (ri, s) in enumerate(PILOT_STEPS)
             for r in 1:PILOT_R]
    run_cells(cells)
    pilot_report()
end
function pilot_k50(R, x)
    v = [r for r in R if r.j == x.j && r.stage in ("pilot", "pilot2") && r.arm == 0]
    ks = sort(unique(r.k for r in v))
    f = [mean(r.alive for r in v if r.k == k) for k in ks]; n = [count(r -> r.k == k, v) for k in ks]
    sk = spearman_karber(log.(ks), f, n)
    (; ks, f, n, k50 = exp(sk.m), covered = sk.covered)
end
function pilot_report()
    R = kruns()
    for x in kinputs()
        pk = pilot_k50(R, x)
        sayf("    %2d %-34s L %4d ω %.3f  table k50 %.4f  pilot k50 %.4f (ratio %.2f)%s  alive by rate: %s", x.j, x.id, x.L, x.omega, x.k50_tab, pk.k50,
             pk.k50 / x.k50_tab, pk.covered ? "" : "  NOT COVERED", join([@sprintf("%.0f", 100f) for f in pk.f], " "))
    end
end
function kernel_mode()
    ins = kinputs(); R = kruns()
    say("== kernel mode: rates = pilot k50 · exp($(collect(MAIN_STEPS))), $MAIN_R runs each, at ω·e^{−$NDELTA} and ω·e^{+$NDELTA} (same seeds in both arms) ==")
    cells = NamedTuple[]
    for x in ins
        pk = pilot_k50(R, x)
        for arm in (-1, 1), (ri, s) in enumerate(MAIN_STEPS), r in 1:MAIN_R
            push!(cells, (; j = x.j, id = x.id, L = x.L, stage = "main", arm, ri, r, k = pk.k50 * exp(s)))
        end
    end
    run_cells(cells)
end

"The kernel's k50 per arm and its elasticity over ±δ for input j from the main (and extension) runs; paired bootstrap."
function kernel_elasticity(R, j; boot = 1000)
    v = [r for r in R if r.j == j && r.stage in ("main", "ext") && r.arm != 0]
    ks = sort(unique(round(r.k; sigdigits = 12) for r in v))
    cell(arm, k) = sort([r for r in v if r.arm == arm && round(r.k; sigdigits = 12) == k]; by = r -> r.r)
    C = Dict((arm, k) => cell(arm, k) for arm in (-1, 1) for k in ks)
    est(sel) = begin
        out = Float64[]; cov = true
        for arm in (-1, 1)
            kk = [k for k in ks if !isempty(C[(arm, k)])]
            f = [mean(C[(arm, k)][i].alive for i in sel(arm, k)) for k in kk]; n = [length(C[(arm, k)]) for k in kk]
            s = spearman_karber(log.(kk), f, n); push!(out, s.m); cov &= s.covered
        end
        (out[1], out[2], cov)
    end
    m1, m2, cov = est((arm, k) -> eachindex(C[(arm, k)]))
    rng = MersenneTwister(NSELSEED + 1000 + j)
    eb = Float64[]; b1 = Float64[]; b2 = Float64[]
    for _ in 1:boot
        draw = Dict(k => rand(rng, 1:max(length(C[(-1, k)]), length(C[(1, k)]), 1), max(length(C[(-1, k)]), length(C[(1, k)]), 1)) for k in ks)
        a, b, _ = est((arm, k) -> [min(i, length(C[(arm, k)])) for i in draw[k]])     # the same replicate indices in both arms (paired)
        push!(eb, (b - a) / (2NDELTA)); push!(b1, a); push!(b2, b)
    end
    # logistic comparator: y ~ 1 + log k + arm
    X = [ones(length(v)) [log(r.k) for r in v] Float64.([r.arm for r in v])]; y = Float64.([r.alive for r in v])
    β, conv = logistic(X, y)
    elog = conv && β[2] < 0 && maximum(abs, β) < 1e4 ? (-2β[3] / β[2]) / (2NDELTA) : NaN
    fr(arm) = [isempty(C[(arm, k)]) ? NaN : mean(r.alive for r in C[(arm, k)]) for k in ks]
    (; ks, fm = fr(-1), fp = fr(1), k50m = exp(m1), k50p = exp(m2), e = (m2 - m1) / (2NDELTA), lo = quantile(eb, 0.025), hi = quantile(eb, 0.975),
     se = std(eb), covered = cov, elog, width = -1 / β[2], sem = std(b1), sep = std(b2), n = length(v))
end

# ---- emulators over the same ±δ
const EMHEAD = "key,j,closure,arm,omega,k50,e,status,evals,continued,sec"
function emfd_mode()
    ins = kinputs()
    f = ncache("emfd.csv")
    say("== emfd mode: emulators' k50 at ω·e^{−δ}, ω, ω·e^{+δ} (δ = $NDELTA) and local e_ω ==")
    # 1. the screen at the three points; :rare_alleles at the three points
    have = done_keys(f)
    jobs = [(x, arm) for x in ins for arm in (-1, 0, 1) if !("$(x.j)|screen|$arm" in have)]
    Threads.@threads :greedy for (x, arm) in jobs
        p = Params(readp(x.id); InitFitnessSD = armω(x.omega, arm))
        s = kgrad(p)
        append_line(f, EMHEAD, join(Any["$(x.j)|screen|$arm", x.j, "screen", arm, p.InitFitnessSD, s.k50, s.e, s.status, s.evals, s.continued, round(s.sec; digits = 1)], ","))
        t = time()
        ra = try
            emulated_half_survival_rate(p; rtol = 1e-6).k50
        catch
            NaN
        end
        append_line(f, EMHEAD, join(Any["$(x.j)|rare|$arm", x.j, "rare", arm, p.InitFitnessSD, ra, NaN, "search", 0, false, round(time() - t; digits = 1)], ","))
    end
    # 2. dynamic LD on the grid, and default :spectrum (no LD, grid), at the two arms, from brackets around the screen's k50
    E = read_csv(f)
    ks(j, arm) = (r = [pf64(e.k50) for e in E if e.key == "$j|screen|$arm"]; isempty(r) ? NaN : r[1])
    have = done_keys(f)
    jobs2 = [(x, arm, cl) for cl in ("gridld", "nold") for x in ins for arm in (-1, 1) if !("$(x.j)|$cl|$arm" in have)]
    Threads.@threads :greedy for (x, arm, cl) in jobs2
        p = Params(readp(x.id); InitFitnessSD = armω(x.omega, arm))
        k0 = ks(x.j, arm)
        br = !isfinite(k0) ? nothing : cl == "gridld" ? (k0 * exp(-0.03), k0 * exp(0.03)) : (k0 * exp(-0.1), k0 * exp(0.4))
        s = kgrad(p; spectrum = cl == "gridld" ? (; ld = :dynamic) : (;), bracket = br)
        append_line(f, EMHEAD, join(Any["$(x.j)|$cl|$arm", x.j, cl, arm, p.InitFitnessSD, s.k50, s.e, s.status, s.evals, s.continued, round(s.sec; digits = 1)], ","))
    end
    say("  done; see mode local")
end
function emfd()
    E = read_csv(ncache("emfd.csv"))
    Dict((parse(Int, e.j), String(e.closure), parse(Int, e.arm)) => (; k50 = pf64(e.k50), e = pf64(e.e), status = String(e.status), continued = e.continued == "true") for e in E)
end

function local_mode()
    ins = kinputs(); R = kruns(); E = isfile(ncache("emfd.csv")) ? emfd() : Dict()
    fits = deserialize(ncache("fitsA.jls"))
    say("== local mode: the kernel's elasticity of k50 in ω over ω·e^{±$NDELTA} against the emulators' over the same interval ==")
    say("  kernel: Spearman–Kärber mean of log k50 per arm (fractions alive made monotone; $MAIN_R runs per rate and arm), e = Δ log k50 / (2δ),")
    say("  paired bootstrap percentile 95% CI (1000 resamples of replicate indices within a rate, the same in both arms); logistic comparator")
    say("  y ~ 1 + log k + arm ('elog', 'width' = −1/βk in log k). Emulators: fd = finite difference of their k50 over the same arms; loc = local")
    say("  (k50_gradient) at ω. Table: model A's local e_ω in the input's bin (se).")
    fdof(j, cl) = (a = get(E, (j, cl, -1), nothing); b = get(E, (j, cl, 1), nothing); a === nothing || b === nothing ? NaN : (log(b.k50) - log(a.k50)) / (2NDELTA))
    rows = NamedTuple[]
    for x in ins
        any(r -> r.j == x.j && r.stage == "main", R) || continue
        ke = kernel_elasticity(R, x.j)
        eA, sA = slopeA((; β = fits[x.L].β, V = fits[x.L].V), nbin(x.omega))
        sc0 = get(E, (x.j, "screen", 0), (; k50 = NaN, e = NaN, status = "", continued = false))
        kmid = sqrt(ke.k50m * ke.k50p)
        r = (; x.j, x.id, x.L, x.omega, ke, eA, sA, fd_screen = fdof(x.j, "screen"), fd_grid = fdof(x.j, "gridld"), fd_nold = fdof(x.j, "nold"),
             fd_rare = fdof(x.j, "rare"), loc = sc0.e, k50_screen = sc0.k50, kmid,
             lvl_m = get(E, (x.j, "screen", -1), (; k50 = NaN)).k50 / ke.k50m, lvl_p = get(E, (x.j, "screen", 1), (; k50 = NaN)).k50 / ke.k50p)
        push!(rows, r)
        sayf("  %2d %-34s L %4d ω %.3f", x.j, x.id, x.L, x.omega)
        sayf("       kernel: rates %s", join([@sprintf("%.4f", k) for k in ke.ks], " "))
        sayf("               alive ω− %s", join([isnan(f) ? "  -" : @sprintf("%3.0f", 100f) for f in ke.fm], " "))
        sayf("               alive ω+ %s", join([isnan(f) ? "  -" : @sprintf("%3.0f", 100f) for f in ke.fp], " "))
        sayf("               k50(ω−) %.4f (se of log %.3f)  k50(ω+) %.4f (%.3f)  e_ω % .3f [% .3f, % .3f]%s  elog % .3f  width %.3f", ke.k50m, ke.sem, ke.k50p,
             ke.sep, ke.e, ke.lo, ke.hi, ke.covered ? "" : "  GRID DOES NOT COVER", ke.elog, ke.width)
        sayf("       emulator: screen fd % .3f (loc % .3f, %s)  grid+LD fd % .3f  no-LD fd % .3f  rare fd % .3f | table A % .3f (%.3f)", r.fd_screen, r.loc,
             sc0.status, r.fd_grid, r.fd_nold, r.fd_rare, eA, sA)
        sayf("       levels: screen k50 / kernel k50 at ω− %.2f, at ω+ %.2f; table model A k50 / kernel (geometric mean of arms) %.2f", r.lvl_m, r.lvl_p,
             x.k50_tab / kmid)
    end
    isempty(rows) && return
    say("\n  summary (inputs whose grid covers both arms):")
    ok = [r for r in rows if r.ke.covered]
    sayf("    n %d of %d", length(ok), length(rows))
    for (lab, g) in (("kernel", r -> r.ke.e), ("kernel (logistic)", r -> r.ke.elog), ("screen fd", r -> r.fd_screen), ("screen local", r -> r.loc),
                     ("grid+LD fd", r -> r.fd_grid), ("no-LD fd", r -> r.fd_nold), (":rare_alleles fd", r -> r.fd_rare), ("table model A (bin)", r -> r.eA))
        v = filter(isfinite, [g(r) for r in ok])
        isempty(v) && continue
        sayf("    %-22s median % .3f  mean % .3f  [min % .3f, max % .3f]  n %d", lab, median(v), mean(v), minimum(v), maximum(v), length(v))
    end
    w = [1 / r.ke.se^2 for r in ok]
    sayf("    kernel, inverse-variance weighted mean % .3f (se %.3f); screen fd − kernel: mean % .3f, median % .3f; inputs with the screen's fd outside the kernel's CI: %d of %d",
         sum(w .* [r.ke.e for r in ok]) / sum(w), sqrt(1 / sum(w)), mean(r.fd_screen - r.ke.e for r in ok if isfinite(r.fd_screen)),
         med([r.fd_screen - r.ke.e for r in ok]), count(r -> isfinite(r.fd_screen) && !(r.ke.lo <= r.fd_screen <= r.ke.hi), ok), length(ok))
    for (lab, sel) in (("L 100", r -> r.L == 100), ("L 500", r -> r.L == 500), ("L 1000", r -> r.L == 1000), ("ω < 1", r -> r.omega < 1), ("ω ≥ 1", r -> r.omega >= 1))
        v = filter(sel, ok)
        isempty(v) && continue
        ww = [1 / r.ke.se^2 for r in v]
        sayf("    %-7s n %d  kernel weighted mean % .3f (se %.3f)  screen fd mean % .3f  table A mean % .3f", lab, length(v), sum(ww .* [r.ke.e for r in v]) / sum(ww),
             sqrt(1 / sum(ww)), mean(filter(isfinite, [r.fd_screen for r in v])), mean(r.eA for r in v))
    end
    serialize(ncache("local.jls"), rows)
end

# ---------------------------------------------------------------------------------------------
dispatch() = NMODE == "table" ? table_mode() : NMODE == "curve" ? curve_mode() : NMODE == "vgmove" ? vgmove_mode() : NMODE == "select" ? select_mode() :
             NMODE == "pilot" ? pilot_mode() : NMODE == "pilotreport" ? pilot_report() : NMODE == "kernel" ? kernel_mode() : NMODE == "emfd" ? emfd_mode() :
             NMODE == "local" ? local_mode() : isdefined(Main, :extra_dispatch) ? extra_dispatch() : error("unknown mode $NMODE")
if abspath(PROGRAM_FILE) == @__FILE__
    t0 = time()
    LOG[] = open(ncache("out_$(NMODE).txt"), "w")
    println("QTLSim at ", pathof(QTLSim))
    dispatch()
    finish_mode(t0)
    close(LOG[]); LOG[] = stdout
end
