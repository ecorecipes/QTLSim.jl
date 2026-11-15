# How fast does the association between a locus's within-locus deviation and its linked,
# lagging background build up and decay under a moving optimum? (follow-up of
# validate/genic_decline.jl; criteria validate/linked_decay_criteria.md, committed before any
# comparison; running log validate/linked_decay.log).
#
# The derivation (validate/linked_decay_theory.jl, checked in validate/symbolic/check_linked_decay.jl):
# per locus l and class s of partner loci (position bin × split-probability class),
#   T_s = κ(B_s, x_l, x_l),  C_l = Σ_s T_s = Cov(g_rest, d_l²),  y_s = Cov(B_s, x_l),  K_s = κ(B_s, B, x_l, x_l),
#   T_s ← f₃(1 − ρ_s)[T_s + β(K_s + K3_s) − (c_s(κ3_l + C_l) + 2v_l T_s)/S_c],     E[K_s] ≈ −2E[y_s Y_l],
#   y_s ← f₁(1 − ρ_s)[y_s − (c_s(v_l + Y_l) + v_l y_s)/S_c] + noise(v_l Var(B_s)/(2n)),
# and the genic loss through the association is 2β Σ_l C_l (the accounting's other3).
#
# Instrumented kernel runs (reading the World only; `check` verifies bitwise identity with
# `simulate`) measure, every generation from MOVE − 60 to MOVE + 199 at the adult stage, per
# group a of loci l (16 position bins) and partner bin b (16): Σ_{l∈a} of y_{lb}, T_{lb},
# K_{lb}, K3_{lb} = κ(B_b, l, l, l) and y_{lb}Y_l, and the drivers (β, S_c, n, N, per-group
# genic variance, per-bin background variances and covariances).
#
#   level A (the decay rate): the T recursion driven by the run's measured sources K + K3 and
#           κ3, from the measured T at the move, against the run's measured T (C = Σ T);
#   level B (closed in the LD): K from the y recursion (drift-seeded allele–background LD),
#           driven by the run's measured n, N, β, S_c and variances, against the measured C
#           and against the measured sign-symmetric loss 2β C_sym.
#
# Run from the package root (vignettes environment):
#   julia --project=vignettes validate/linked_decay.jl check                  # bitwise identity with simulate
#   julia -t 6 --project=vignettes validate/linked_decay.jl kernel SET R      # instrumented runs (resumable)
#   julia --project=vignettes validate/linked_decay.jl analyse SET R          # levels A and B, per input and verdicts
# SET: D (development: genic_decline's discovery inputs, seeds as there), F (confirmation:
# fresh inputs and seeds, chosen by `fresh_rows`), or an input id. Raw records go to
# validate/linked_decay_raw/ (git-ignored).

include(joinpath(@__DIR__, "genic_decline.jl"))       # spectrum_diagnosis.jl: readp, kernel_params, MOVE, KBURN, KSEED, C0, inputs
include(joinpath(@__DIR__, "linked_decay_theory.jl"))

const LD_RAW = joinpath(@__DIR__, "linked_decay_raw")
const LD_FROM, LD_TO = MOVE - 60, MOVE + 200          # instrumented generations [LD_FROM, LD_TO)
const FSEED = 80261114000                             # confirmation seeds: FSEED + 100j + r
const NFRESH_PER_STRATUM = 2

# ---------------------------------------------------------------------------------------
# inputs

function fresh_rows()
    rows = filter(r -> r.set == "", read_inputs_csv())
    strata = unique([split(r.id, "/")[1] for r in rows])
    rng = MersenneTwister(20261121)
    out = NamedTuple[]
    for s in sort(strata)
        cand = filter(r -> startswith(r.id, s * "/"), rows)
        append!(out, cand[randperm(rng, length(cand))[1:NFRESH_PER_STRATUM]])
    end
    out
end

ld_rows(set) = set == "F" ? fresh_rows() : set == "D" ? inset(read_inputs_csv(), "D") :
               filter(r -> r.id == set, read_inputs_csv())
ld_seed(set, j, r) = (set == "F" ? FSEED : KSEED) + 100j + r
ld_file(id, r, set) = joinpath(LD_RAW, replace(id, "/" => "__", ".txt" => "") * "@" * (set == "F" ? "F" : "D") * "#$r.jls")

# ---------------------------------------------------------------------------------------
# measurement

mutable struct LDBuf
    L::Int; NB::Int
    binl::Vector{Int}
    xbar::Vector{Float64}
    S1::Matrix{Float64}; S2::Matrix{Float64}; S3::Matrix{Float64}; S4::Matrix{Float64}; S5::Matrix{Float64}   # (NB, L)
    Sd2::Vector{Float64}; Sd3::Vector{Float64}; SY::Vector{Float64}
    Bs::Matrix{Float64}; Gc::Vector{Float64}     # per chromosome bin sums (NB, M), totals (M)
    x::Vector{Float64}; e::Vector{Float64}
end

function LDBuf(L, Nmax; NB = LD_NB)
    LDBuf(L, NB, [ld_bin(l - 1, L; nb = NB) for l in 1:L], zeros(L),
          zeros(NB, L), zeros(NB, L), zeros(NB, L), zeros(NB, L), zeros(NB, L),
          zeros(L), zeros(L), zeros(L), zeros(NB, 2Nmax), zeros(2Nmax), zeros(L), zeros(NB))
end

# chromosome c of individual i into buf
ldchr!(buf, w, c, i) = loadchr!(buf, w.pool, w.chr, c, i)

"""
Per-generation association statistics of the adults (2N chromosomes, population denominators).
Returns a NamedTuple of group × bin matrices (sums over l ∈ a) and the per-locus summaries.
"""
function ld_measure!(B::LDBuf, w)
    L, NB, N = B.L, B.NB, w.N
    M = 2N
    x = B.x
    fill!(B.xbar, 0)
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        ldchr!(x, w, c, i)
        B.xbar .+= x
        g = 0.0
        @inbounds for b in 1:NB
            B.Bs[b, m] = 0.0
        end
        @inbounds for l in 1:L
            B.Bs[B.binl[l], m] += x[l]; g += x[l]
        end
        B.Gc[m] = g
    end
    B.xbar ./= M
    Bbar = vec(sum(view(B.Bs, :, 1:M); dims = 2)) ./ M
    Gbar = sum(view(B.Gc, 1:M)) / M
    for A in (B.S1, B.S2, B.S3, B.S4, B.S5)
        fill!(A, 0)
    end
    fill!(B.Sd2, 0); fill!(B.Sd3, 0); fill!(B.SY, 0)
    e = B.e
    S1, S2, S3, S4, S5 = B.S1, B.S2, B.S3, B.S4, B.S5
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        ldchr!(x, w, c, i)
        dG = B.Gc[m] - Gbar
        @inbounds for b in 1:NB
            e[b] = B.Bs[b, m] - Bbar[b]
        end
        @inbounds for l in 1:L
            d = x[l] - B.xbar[l]
            d2 = d * d; d3 = d2 * d
            dB = dG - d
            ob = B.binl[l]
            e[ob] -= d
            dB2 = dB * d2
            @simd for b in 1:NB
                eb = e[b]
                S1[b, l] += eb * d
                S2[b, l] += eb * d2
                S3[b, l] += eb * dB2
                S4[b, l] += eb * d3
                S5[b, l] += eb * dB
            end
            e[ob] += d
            B.Sd2[l] += d2; B.Sd3[l] += d3; B.SY[l] += dB * d
        end
    end
    Tab = zeros(NB, NB); Kab = zeros(NB, NB); K3ab = zeros(NB, NB); yab = zeros(NB, NB); yYab = zeros(NB, NB)
    V = zeros(NB); V2 = zeros(NB); k3 = zeros(NB); gb = zeros(NB)
    sumC = 0.0; sk3 = 0.0; sk3sq = 0.0; sk3C = 0.0; genic = 0.0
    @inbounds for l in 1:L
        a = B.binl[l]
        v = B.Sd2[l] / M; κ3 = B.Sd3[l] / M; Y = B.SY[l] / M
        V[a] += v; V2[a] += v^2; k3[a] += κ3; gb[a] += v
        Cl = 0.0
        for b in 1:NB
            y = S1[b, l] / M; T = S2[b, l] / M
            K = S3[b, l] / M - (S5[b, l] / M) * v - 2y * Y
            K3 = S4[b, l] / M - 3y * v
            Tab[a, b] += T; Kab[a, b] += K; K3ab[a, b] += K3; yab[a, b] += y; yYab[a, b] += y * Y
            Cl += T
        end
        sumC += Cl; sk3 += κ3; sk3sq += κ3^2; sk3C += κ3 * Cl; genic += 2v
    end
    # bin covariances and covariance with the chromosome total
    Cbb = zeros(NB, NB); cb = zeros(NB)
    for mm in 1:M
        dG = B.Gc[mm] - Gbar
        for b in 1:NB
            eb = B.Bs[b, mm] - Bbar[b]
            cb[b] += eb * dG
            for b2 in 1:NB
                Cbb[b, b2] += eb * (B.Bs[b2, mm] - Bbar[b2])
            end
        end
    end
    Cbb ./= M; cb ./= M
    VGchr = sum(abs2, view(B.Gc, 1:M) .- Gbar) / M
    (; Tab, Kab, K3ab, yab, yYab, V, V2, k3, g = gb, c = cb, Cbb, sumC, sk3, sk3sq, sk3C, genic, VGchr)
end

# ---------------------------------------------------------------------------------------
# the instrumented run (gd_run's structure; draws nothing)

const SCAL = [:gen, :N, :nsurv, :theta, :gbar, :VG, :VGchr, :Sc, :beta, :beta2, :genic, :sumC, :sk3, :sk3sq, :sk3C, :VE, :omega2, :Vm]
const MATS = [:Tab, :Kab, :K3ab, :yab, :yYab, :Cbb]
const VECS = [:V, :V2, :k3, :g, :c]

function ld_run(p, seed; from = LD_FROM, to = LD_TO)
    q = kernel_params(p)
    q.Selection === :viability || error("viability only")
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    L = q.NoQTL
    Bf = LDBuf(L, q.PopSize)
    VE = QTLSim.envVar(q)
    rec = Dict{Symbol,Any}(s => Float64[] for s in SCAL)
    for s in MATS
        rec[s] = Matrix{Float32}[]
    end
    for s in VECS
        rec[s] = Vector{Float32}[]
    end
    vgser = Float64[]
    extinct_gen = -1
    for it in w.gen:q.NoGenerations-1
        it >= to && break
        w.gen = it
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        N = w.N
        if it >= from
            g = view(w.genotype, 1:N)
            gb = mean(g); vg = var(g; corrected = false)
            push!(vgser, var(g))
            ω2 = w.ω^2; Sg = ω2 + VE
            π = [C0 * sqrt(ω2 / Sg) * exp(-(w.genotype[i] - w.θ)^2 / (2Sg)) for i in 1:N]
            πbar = mean(π)
            u = [w.genotype[i] - gb for i in 1:N]; vv = u .^ 2 .- vg; yy = π ./ πbar .- 1
            β = sum(yy .* u) / N / vg
            a11 = sum(u .* u); a12 = sum(u .* vv); a22 = sum(vv .* vv); r1 = sum(u .* yy); r2 = sum(vv .* yy)
            dt = a11 * a22 - a12^2
            β2 = dt > 0 ? (a22 * r1 - a12 * r2) / dt : 0.0
            ms = ld_measure!(Bf, w)
            Sc = ω2 + VE + ms.VGchr
            vals = (gen = it, N = N, nsurv = w.nsurv, theta = w.θ, gbar = gb, VG = var(g), VGchr = ms.VGchr, Sc = Sc,
                    beta = β, beta2 = β2, genic = ms.genic, sumC = ms.sumC, sk3 = ms.sk3, sk3sq = ms.sk3sq, sk3C = ms.sk3C,
                    VE = VE, omega2 = ω2, Vm = QTLSim.muteVar(q))
            for s in SCAL
                push!(rec[s], Float64(vals[s]))
            end
            for s in MATS
                push!(rec[s], Float32.(getfield(ms, s)))
            end
            for s in VECS
                push!(rec[s], Float32.(getfield(ms, s)))
            end
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    (; seed, extinct_gen, rec, vgser, L, u = q.RecombPoissonMean, PopSize = q.PopSize)
end

function ld_check(; id = "QT1000-Possion3.0/Input-09697.txt", j = 89, r = 1)
    p = readp(id)
    t = @elapsed tr = ld_run(p, KSEED + 100j + r; to = LD_FROM + 30)
    o = simulate(kernel_params(p); seed = KSEED + 100j + r, init = :equilibrium, burnin = KBURN)
    k = [findfirst(==(g), o.gen) for g in LD_FROM:LD_FROM+length(tr.vgser)-1]
    same = all(o.geno_var[k] .== tr.vgser)
    @printf("check %s: %d generations, V_G bitwise identical to simulate: %s (%.1f s, %.2f s per instrumented generation incl. burn-in share)\n",
            id, length(tr.vgser), same, t, t / length(tr.vgser))
    same || error("instrumentation changed the run")
end

function ld_kernel(set, R)
    mkpath(LD_RAW)
    rows = ld_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(ld_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) instrumented runs for set $set, R = $R, $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time()
    Threads.@threads :greedy for (id, j, k) in jobs
        tr = ld_run(readp(id), ld_seed(set, j, k))
        serialize(ld_file(id, k, set), tr)
        d = Threads.atomic_add!(done, 1) + 1
        @printf("  %d/%d %s #%d  %.0f s\n", d, length(jobs), id, k, time() - t0)
    end
    @printf("done: %d runs, %.0f s wall\n", length(jobs), time() - t0)
end

# ---------------------------------------------------------------------------------------
# theory driven by a run

function drivers(rec, k; avg = nothing)
    ks = avg === nothing ? (k:k) : avg
    mv(s) = mean(Float64.(rec[s][kk]) for kk in ks)
    ms(s) = mean(rec[s][kk] for kk in ks)
    LDDrivers(β = rec[:beta][k], Sc = mean(rec[:Sc][ks]), dS = mean(rec[:VG][ks]) / (rec[:omega2][k] + rec[:VE][k] + mean(rec[:VG][ks])),
              n = mean(rec[:nsurv][ks]), N = rec[:N][k], Vm = rec[:Vm][k],
              V = Float64.(ms(:V)), V2 = Float64.(ms(:V2)), k3 = Float64.(ms(:k3)), g = Float64.(ms(:g)),
              c = Float64.(ms(:c)), Cbb = Float64.(ms(:Cbb)))
end

const PRE_AVG = 1:20      # the stationary LD state from the drivers averaged over t = −60..−41

"""
Level B: from t = −60, y moments at their stationary state for the drivers averaged over
−60..−41, T = 0; each generation K = −2 P 1, T_step!, y_step!, driven by the run. Returns per
generation the predicted Σ C, Σ E[y Y] (= −K/2) and Σ_s y means.
"""
function levelB(tr, G::LDGroups; mode = :derived)
    rec = tr.rec
    n = length(rec[:gen])
    st = LDState(G)
    y_stationary!(st, G, drivers(rec, 1; avg = PRE_AVG))
    predC = fill(NaN, n); predyY = fill(NaN, n)
    NA = size(G.n, 1)
    for k in 1:n
        d = drivers(rec, k)
        predC[k] = sum(sum(st.T[a]) for a in 1:NA)
        predyY[k] = sum(sum(st.P[a]) for a in 1:NA)
        for a in 1:NA
            G.nb[a] == 0 && continue
            K = -2 .* vec(sum(st.P[a]; dims = 2))
            T_step!(st.T[a], G, d, a, K; mode)
            y_step!(st, G, d, a)
        end
    end
    (; predC, predyY)
end

"""
Level A: from the move (t = 0), T at the run's measured T_ab split over ρ classes, driven by the
run's measured K + K3 (split by `split_weights`) and, for :derived, κ3. Returns predicted Σ C per
generation (NaN before the move).
"""
function levelA(tr, G::LDGroups; mode = :derived)
    rec = tr.rec
    n = length(rec[:gen])
    k0 = findfirst(==(MOVE), Int.(rec[:gen]))
    NA = size(G.n, 1); ns = size(G.n, 2)
    predC = fill(NaN, n)
    k0 === nothing && return predC
    T = [zeros(ns) for _ in 1:NA]
    for a in 1:NA
        wts = split_weights(G, a, rec[:nsurv][k0])
        for s in 1:ns
            T[a][s] = rec[:Tab][k0][a, G.bofs[s]] * wts[s]
        end
    end
    for k in k0:n
        d = drivers(rec, k)
        predC[k] = sum(sum(T[a]) for a in 1:NA)
        for a in 1:NA
            G.nb[a] == 0 && continue
            wts = split_weights(G, a, d.n)
            K = [(rec[:Kab][k][a, G.bofs[s]] + rec[:K3ab][k][a, G.bofs[s]]) * wts[s] for s in 1:ns]
            T_step!(T[a], G, d, a, K; mode, k3 = mode === :derived ? d.k3[a] : 0.0)
        end
    end
    predC
end

# ---------------------------------------------------------------------------------------
# analysis

const WINS = [(-60, -1), (0, 24), (25, 49), (50, 99), (100, 199)]
const AMODES = (:derived, :nosel, :spectrum)

function load_ld(id, R, set)
    out = []
    for k in 1:R
        f = ld_file(id, k, set)
        isfile(f) && push!(out, deserialize(f))
    end
    out
end

# measured per-generation series of one run (aligned on t = gen − MOVE)
function measured(tr)
    r = tr.rec
    t = Int.(r[:gen]) .- MOVE
    Csym = [r[:sumC][k] - (r[:sk3sq][k] > 0 ? r[:sk3C][k] / r[:sk3sq][k] : 0.0) * r[:sk3][k] for k in eachindex(t)]
    (; t, C = r[:sumC], Csym, K = [sum(Float64.(m)) for m in r[:Kab]], K3 = [sum(Float64.(m)) for m in r[:K3ab]],
     yY = [sum(Float64.(m)) for m in r[:yYab]], beta = r[:beta], beta2 = r[:beta2], genic = r[:genic])
end

# mean over runs alive and generations in window of f(run series, k)
function winmean(series, wdw)
    s = 0.0; c = 0
    for (t, v) in series
        for k in eachindex(t)
            wdw[1] <= t[k] <= wdw[2] && isfinite(v[k]) && (s += v[k]; c += 1)
        end
    end
    c > 0 ? s / c : NaN
end

function analyse_input(r, R, set; io = stdout)
    runs = load_ld(r.id, R, set)
    isempty(runs) && return nothing
    G = ld_groups(runs[1].L, runs[1].u)
    res = Dict{Symbol,Any}()
    ser = Dict{Symbol,Vector{Any}}()
    cum = Dict{Symbol,Vector{Float64}}(s => Float64[] for s in (:meas_loss, :predB_loss, :meas_losssym))
    for tr in runs
        m = measured(tr)
        k0 = findfirst(==(-1), m.t)
        g0 = k0 === nothing ? NaN : mean(m.genic[max(1, k0 - 19):k0])
        push!(get!(ser, :C, []), (m.t, m.C)); push!(get!(ser, :Csym, []), (m.t, m.Csym))
        push!(get!(ser, :K, []), (m.t, m.K)); push!(get!(ser, :m2yY, []), (m.t, -2 .* m.yY))
        push!(get!(ser, :yY, []), (m.t, m.yY))
        for md in AMODES
            push!(get!(ser, Symbol(:A_, md), []), (m.t, levelA(tr, G; mode = md)))
        end
        B = levelB(tr, G)
        push!(get!(ser, :B_C, []), (m.t, B.predC)); push!(get!(ser, :B_yY, []), (m.t, B.predyY))
        # cumulative losses over 0..99 relative to the genic variance before the move (runs alive at +100)
        k100 = findfirst(==(100), m.t)
        if k0 !== nothing && k100 !== nothing
            rg = k0+1:k100-1
            push!(cum[:meas_loss], sum(2 .* m.beta[rg] .* m.C[rg]) / g0)
            push!(cum[:meas_losssym], sum(2 .* m.beta[rg] .* m.Csym[rg]) / g0)
            push!(cum[:predB_loss], sum(2 .* m.beta[rg] .* B.predC[rg]) / g0)
        end
    end
    for (s, v) in ser
        res[s] = Dict(wdw => winmean(v, wdw) for wdw in WINS)
    end
    for (s, v) in cum
        res[s] = isempty(v) ? NaN : mean(v)
    end
    res[:nalive100] = length(cum[:meas_loss])
    res[:L] = runs[1].L; res[:u] = runs[1].u; res[:nruns] = length(runs)
    @printf(io, "\n%s  L=%d u=%g runs %d (alive at +100: %d)\n", r.id, res[:L], res[:u], length(runs), res[:nalive100])
    @printf(io, "  %-12s", "window")
    cols = [:C, :Csym, :A_derived, :A_nosel, :A_spectrum, :B_C, :K, :m2yY, :yY, :B_yY]
    for c in cols
        @printf(io, " %11s", c)
    end
    println(io)
    for wdw in WINS
        @printf(io, "  %4d..%-6d", wdw...)
        for c in cols
            @printf(io, " %11.3e", res[c][wdw])
        end
        println(io)
    end
    @printf(io, "  cumulative loss 0..99 / genic(pre): measured 2βC %+.3f (sign-symmetric 2βC_sym %+.3f), level B %+.3f\n",
            res[:meas_loss], res[:meas_losssym], res[:predB_loss])
    res
end

lr(a, b) = (isfinite(a) && isfinite(b) && a / b > 0) ? log(a / b) : NaN
medf(v) = (u = filter(isfinite, v); isempty(u) ? NaN : median(u))
mabs(v) = medf(abs.(v))

function ld_analyse(set, R; io = stdout)
    rows = ld_rows(set)
    println(io, "\n== linked decay, set $set, R = $R: per input, window means over runs alive (sums over all loci; C, K in genic² units) ==")
    all = [(r, analyse_input(r, R, set; io)) for r in rows]
    all = filter(x -> x[2] !== nothing, all)
    # eligibility (criteria): measured cumulative 2βC over 0..99 ≤ −0.05 of genic before the move
    elig = filter(x -> isfinite(x[2][:meas_loss]) && x[2][:meas_loss] <= -0.05, all)
    println(io, "\n== summary over $(length(elig)) eligible inputs (of $(length(all))) ==")
    println(io, "H1 (source): measured K / (−2 y Y), median over inputs:")
    for wdw in WINS
        @printf(io, "  %4d..%-4d %.3f\n", wdw..., medf([x[2][:K][wdw] / x[2][:m2yY][wdw] for x in all]))
    end
    println(io, "H2 (allele–background LD): level B Σ E[yY] / measured Σ yY, median (log ratio median, |.| median):")
    for wdw in WINS
        v = [lr(x[2][:B_yY][wdw], x[2][:yY][wdw]) for x in all]
        @printf(io, "  %4d..%-4d ratio %.3f  log %+.3f  |log| %.3f\n", wdw..., exp(medf(v)), medf(v), mabs(v))
    end
    println(io, "P1 (decay rate, level A): log(predicted ΣC / measured ΣC) over eligible inputs, median and median |.|:")
    for md in AMODES, wdw in WINS[3:5]
        v = [lr(x[2][Symbol(:A_, md)][wdw], x[2][:C][wdw]) for x in elig]
        @printf(io, "  %-9s %4d..%-4d median %+.3f  |.| %.3f  (n %d)\n", md, wdw..., medf(v), mabs(v), count(isfinite, v))
    end
    println(io, "P2 (closed, level B): log(predicted ΣC / measured ΣC) and log(predicted loss / measured sign-symmetric loss), 0..99:")
    for wdw in WINS[3:5]
        v = [lr(x[2][:B_C][wdw], x[2][:C][wdw]) for x in elig]
        @printf(io, "  ΣC %4d..%-4d median %+.3f  |.| %.3f  (n %d)\n", wdw..., medf(v), mabs(v), count(isfinite, v))
    end
    v = [lr(x[2][:predB_loss], x[2][:meas_losssym]) for x in elig]
    @printf(io, "  cumulative loss 0..99: median %+.3f  |.| %.3f  (n %d)\n", medf(v), mabs(v), count(isfinite, v))
    all
end

# ---------------------------------------------------------------------------------------
# exploratory (criteria: not a verdict): the spectrum's linked-selection term with the derived
# persistence of the allele–background LD, λ_c = (1 − ρ_c)(1 − 1/(2n))(1 − V_chr/S_c),
# V_chr = V/2, S_c = S − V/2, in genic_decline's time-dependent form (Q_c, R_c from zero at the
# move); driven by genic_decline's instrumented runs, as its `driven_linked`.

function transient_derived!(tr::Transient, sp, b, n, V, S)
    g = V / (4n)
    c = 0.0
    sel = max(1 - (V / 2) / (S - V / 2), 0.0)
    for k in 2:length(sp.ω)
        λc = (1 - sp.ρ[k]) * (1 - 1 / (2n)) * sel
        tr.Q[k] = λc^2 * tr.Q[k] + g
        c += sp.ω[k] * b * (2tr.R[k] + b * tr.Q[k])
        tr.R[k] = λc * (tr.R[k] + b * tr.Q[k])
    end
    max(c, 0.0)
end

function driven_derived(p, x; spectrum_kw = (;))
    c = x.cols
    eq = QTLSim.spectrum_equilibrium(p; spectrum_kw...)
    sp = QTLSim.spectrum_scale!(deepcopy(eq.spectrum), c[:GA][1])
    sp.linkage = :none
    ω2 = QTLSim.initFitSD(p)^2
    K = length(sp.ω)
    tr = Transient(zeros(K), zeros(K))
    n0 = max(c[:nsurv][1], 1); V0 = QTLSim.spectrum_vg(sp); S0 = ω2 + c[:VP][1]
    for k in 2:K
        λc = (1 - sp.ρ[k]) * (1 - 1 / (2n0)) * max(1 - (V0 / 2) / (S0 - V0 / 2), 0.0)
        tr.Q[k] = V0 / (4n0) / (1 - λc^2)
    end
    out = fill(NaN, length(c[:gen]))
    for k in eachindex(c[:gen])
        out[k] = QTLSim.spectrum_vg(sp)
        λ = c[:theta][k] - c[:zbar][k]; S = ω2 + c[:VP][k]
        QTLSim.spectrum_select(sp, λ, S)
        n = max(c[:nsurv][k], 1)
        cc = transient_derived!(tr, sp, λ / S, n, out[k], S)
        QTLSim.spectrum_step!(sp, 1 / (1 / n + 2cc))
        QTLSim.spectrum_mutate!(sp, c[:N][k])
    end
    out
end

function ld_spectrum_proto(; R = 8, ts = (50, 100, 200), set = "D", io = stdout)
    println(io, "\n== exploratory: driven spectrum, genic(t)/genic(−50..−1), kernel vs Barton (default) vs derived transient (genic_decline's runs, set $set) ==")
    acc = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for r in inset(read_inputs_csv(), set)
        p = readp(r.id)
        vals = Dict((m, t) => Float64[] for m in (:kernel, :barton, :derived) for t in ts)
        for k in 1:R
            f = joinpath(GD_RAW_DIR, replace(r.id, "/" => "__", ".txt" => "") * "#$k.jls")
            isfile(f) || continue
            x = deserialize(f)
            T = terms(x)
            k0 = findfirst(==(-1), T.t); gk = mean(T.GA[k0-49:k0])
            sb = driven_linked(p, x; mode = :barton); sd = driven_derived(p, x)
            for t in ts
                kt = findfirst(==(t), T.t); kt === nothing && continue
                push!(vals[(:kernel, t)], T.GA[kt] / gk)
                push!(vals[(:barton, t)], sb[kt] / mean(sb[k0-49:k0]))
                push!(vals[(:derived, t)], sd[kt] / mean(sd[k0-49:k0]))
            end
        end
        @printf(io, "%-34s", r.id)
        for t in ts
            isempty(vals[(:kernel, t)]) && continue
            mk, mb, md = mean(vals[(:kernel, t)]), mean(vals[(:barton, t)]), mean(vals[(:derived, t)])
            @printf(io, "  +%d n%d kernel %.3f barton %.3f derived %.3f", t, length(vals[(:kernel, t)]), mk, mb, md)
            push!(get!(acc, (:barton, t), Float64[]), log(mb / mk)); push!(get!(acc, (:derived, t), Float64[]), log(md / mk))
        end
        println(io)
    end
    for t in ts
        @printf(io, "  +%d: log(model/kernel) median: barton %+.3f (|.| %.3f), derived %+.3f (|.| %.3f), n %d\n", t,
                median(acc[(:barton, t)]), median(abs.(acc[(:barton, t)])), median(acc[(:derived, t)]), median(abs.(acc[(:derived, t)])),
                length(acc[(:barton, t)]))
    end
end

# cross-check (D only): the same runs as genic_decline's runs 1..R (same seeds, nothing drawn by
# either instrumentation): cumulative 2β₂ΣC_sym (chromosome level, here) against genic_decline's
# sym3 (individual level), over 0..99, relative to the genic variance before the move. Reads
# genic_decline's raw runs from GD_RAW_DIR (read only; default the genic-decline worktree's).
const GD_RAW_DIR = get(ENV, "GD_RAW_DIR", joinpath(@__DIR__, "..", "..", "wt-genic-decline", "validate", "genic_decline_raw"))

function ld_crosscheck(R; io = stdout)
    println(io, "\n== cross-check: cumulative 0..99 / genic(pre), this study's 2β₂ΣC_sym vs genic_decline's sym3, same runs ==")
    for r in ld_rows("D")
        for k in 1:R
            f1 = ld_file(r.id, k, "D")
            f2 = joinpath(GD_RAW_DIR, replace(r.id, "/" => "__", ".txt" => "") * "#$k.jls")
            (isfile(f1) && isfile(f2)) || continue
            a = deserialize(f1); b = deserialize(f2)
            m = measured(a); T = terms(b)
            k0 = findfirst(==(-1), m.t); k100 = findfirst(==(100), m.t)
            j0 = findfirst(==(-1), T.t); j100 = findfirst(==(100), T.t)
            (k0 === nothing || k100 === nothing || j0 === nothing || j100 === nothing) && continue
            g0 = mean(m.genic[max(1, k0 - 19):k0]); h0 = mean(T.GA[max(1, j0 - 19):j0])
            mine = sum(2 .* m.beta2[k0+1:k100-1] .* m.Csym[k0+1:k100-1]) / g0
            theirs = sum(T.sym3[j0+1:j100-1]) / h0
            @printf(io, "  %-34s #%d  2β₂ΣC_sym %+.3f  sym3 %+.3f  genic(pre) %.4g vs %.4g\n", r.id, k, mine, theirs, g0, h0)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "check" : ARGS[1]
    if mode == "check"
        ld_check()
    elseif mode == "kernel"
        ld_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "analyse"
        ld_analyse(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "crosscheck"
        ld_crosscheck(parse(Int, ARGS[2]))
    elseif mode == "proto"
        ld_spectrum_proto(; set = length(ARGS) >= 2 ? ARGS[2] : "D")
    end
end
