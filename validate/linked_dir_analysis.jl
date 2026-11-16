# Analysis of validate/linked_dir.jl's runs (included by it).

load_dr(id, R, set) = [deserialize(dr_file(id, k, set)) for k in 1:R if isfile(dr_file(id, k, set))]

# linked_decay's record names (Tab, Kab, …) so that drivers/measured/levelA apply unchanged
function drview(tr)
    r = copy(tr.rec)
    for (a, b) in ((:Tab, :T), (:Kab, :K), (:K3ab, :K3), (:yab, :y), (:yYab, :yY))
        r[a] = tr.rec[b]
    end
    (; tr.seed, tr.extinct_gen, rec = r, tr.vgser, tr.L, tr.u, tr.PopSize)
end

# ---------------------------------------------------------------------------------------
# 1. the size of the directional terms (measured), per generation, against the statistics

"""
Per-generation series of one run: measured Σ yY, Σ H, Σ K; the exact directional changes
β Σ δ(yY), β Σ δH, β Σ δK; the identity check 2β Σ_l Y_l u_l; candidate decompositions of β Σ δK:
the PBM product terms −β(6 Σ T_l Y_l + 4 Σ J_l Y_l) and the PBM first-order pair term β Σ E5c.
"""
function dir_series(tr)
    r = tr.rec
    t = Int.(r[:gen]) .- MOVE
    β = r[:beta]
    (; t, yY = tot.(r[:yY]), H = tot.(r[:H]), K = tot.(r[:K]), C = r[:sumC],
     dyY = β .* tot.(r[:dYY]), dH = β .* tot.(r[:dH]), dK = β .* tot.(r[:dK]),
     dyYid = 2 .* β .* (tot.(r[:yC]) .+ tot.(r[:yJ])),
     dyYT = 2 .* β .* tot.(r[:yC]), dyYJ = 2 .* β .* tot.(r[:yJ]),
     dKprod = .-β .* (6 .* tot.(r[:TY]) .+ 4 .* tot.(r[:JY])), dKpbm = β .* tot.(r[:E5c]),
     dT = β .* tot.(r[:dT]), KK3 = β .* (tot.(r[:K]) .+ tot.(r[:K3])),
     dk3 = β .* tot.(r[:dk3]), k3 = r[:sk3], genic = r[:genic])
end

const DTWINS = [(-60, -1), (0, 24), (25, 49), (50, 99), (100, 199)]

function dir_terms(set, R; io = stdout)
    rows = dr_rows(set)
    println(io, "\n== directional terms, set $set, R = $R: window mean of the per-generation change / window mean of the statistic ==")
    keys_ = [(:dyY, :yY), (:dyYid, :yY), (:dyYT, :yY), (:dyYJ, :yY), (:dH, :H), (:dK, :K), (:dKprod, :K), (:dKpbm, :K), (:dT, :C), (:KK3, :C)]
    agg = Dict{Tuple{Symbol,Tuple{Int,Int}},Vector{Float64}}()
    for rw in rows
        runs = load_dr(rw.id, R, set)
        isempty(runs) && continue
        S = [dir_series(tr) for tr in runs]
        @printf(io, "%-34s L=%d runs %d\n", rw.id, runs[1].L, length(runs))
        for wdw in DTWINS
            @printf(io, "   %4d..%-4d", wdw...)
            for (a, b) in keys_
                den = winmean([(s.t, getfield(s, b)) for s in S], wdw)
                v = winmean([(s.t, getfield(s, a)) for s in S], wdw) / abs(den)
                push!(get!(agg, (a, wdw), Float64[]), v)
                @printf(io, " %s/|%s| %+.4f", a, b, v)
            end
            println(io)
        end
    end
    println(io, "\n== medians over inputs (per generation, ÷ |window mean of the statistic|) ==")
    for (a, b) in keys_
        @printf(io, "  %-7s/|%-3s|", a, b)
        for wdw in DTWINS
            @printf(io, "  %d..%d %+.4f", wdw..., medf(agg[(a, wdw)]))
        end
        println(io)
    end
    agg
end

# ---------------------------------------------------------------------------------------
# 2. level B′* driven additionally by the measured directional changes (diagnostic)

# class profile within bins for group a: π_s ∝ P_ss over the active classes of each bin
function bin_profile(g::LDGroups, P)
    ns = size(P, 1)
    w = [max(P[s, s], 0.0) for s in 1:ns]
    tot_ = Dict{Int,Float64}()
    for s in 1:ns
        tot_[g.bofs[s]] = get(tot_, g.bofs[s], 0.0) + w[s]
    end
    cnt = Dict{Int,Int}()
    for s in 1:ns
        cnt[g.bofs[s]] = get(cnt, g.bofs[s], 0) + 1
    end
    [tot_[g.bofs[s]] > 0 ? w[s] / tot_[g.bofs[s]] : 1 / cnt[g.bofs[s]] for s in 1:ns]
end

"""
Level B′* (linked_ld's `levelBfix` with `FIXES[:Bps]`), and optionally (`dir = true`) the measured
directional changes added at each generation's selection stage: δP_st = β π_s π_t δP_{b(s)b(t)}
(measured bin matrix, latest recorded generation), δK_st with row sums β π_s δK_{b(s)} spread over t
symmetrically in proportion to P_tt, δm_s = β π_s δy_{b(s)}; `which` selects :P, :K or both.
Returns per generation Σ K, Σ C, Σ E[yY], Σ H.
"""
function levelBdir(tr, G::LDGroups, ρ2; dir = false, which = (:P, :K), kw = FIXES[:Bps], tk3 = false)
    rec = tr.rec
    n = length(rec[:gen])
    kw0 = merge((; kw...), (hsel = :none,))
    d0 = drivers(rec, 1; avg = PRE_AVG)
    NA = size(G.n, 1)
    act = [a for a in 1:NA if G.nb[a] > 0]
    idx = Dict(a => findall(>(0), G.n[a, :]) for a in act)
    Gr = Dict(a => LDGroups(G.L, G.n[:, idx[a]], G.ρ[:, idx[a]], G.nb, G.bofs[idx[a]]) for a in act)
    sr = Dict(a => LSState(Gr[a]) for a in act)
    ρ2r = Dict(a => ρ2[:, idx[a]] for a in act)
    Tr = Dict(a => zeros(length(idx[a])) for a in act)
    for a in act
        g = Gr[a]; st = sr[a]
        c = class_c(g, d0, a); Σ = class_cov(g, d0, a)
        for it in 1:40000
            chk = it % 25 == 0
            chk && (P0 = copy(st.y.P[a]); H0 = copy(st.H[a]))
            yfix_step!(st, g, d0, a; c, Σ, ρ2 = ρ2r[a], kw0...)
            chk || continue
            sc = max(maximum(abs, st.y.P[a]), 1e-300)
            max(maximum(abs, st.y.P[a] .- P0), maximum(abs, st.H[a] .- H0)) <= 1e-9 * sc && break
        end
    end
    mg = rec[:mgen]
    out = (K = fill(NaN, n), C = fill(NaN, n), yY = fill(NaN, n), H = fill(NaN, n))
    km = 0
    for k in 1:n
        d = drivers(rec, k)
        while km < length(mg) && mg[km+1] <= rec[:gen][k]
            km += 1
        end
        β = rec[:beta][k]
        Kt = 0.0; Ct = 0.0; yt = 0.0; Ht = 0.0
        for a in act
            g = Gr[a]; st = sr[a]
            Ks = source_K(st, a)
            Kt += sum(Ks); Ct += sum(Tr[a]); yt += sum(st.y.P[a]); Ht += sum(st.H[a])
            if tk3     # measured K3 (split by the P profile) in the source and the measured Σκ3 coupling, as level A
                πk = bin_profile(g, st.y.P[a])
                Ks3 = Ks .+ [πk[s] * rec[:K3][k][a, g.bofs[s]] for s in eachindex(Ks)]
                T_step!(Tr[a], g, d, a, Ks3; mode = :derived, k3 = d.k3[a])
            else
                T_step!(Tr[a], g, d, a, Ks; mode = :derived)
            end
            if dir && km > 0
                P = st.y.P[a]
                π_ = bin_profile(g, P)
                bs = g.bofs
                ns = length(bs)
                if :P in which
                    dPb = rec[:dP][km][:, :, a]
                    δP = [β * π_[s] * π_[t] * dPb[bs[s], bs[t]] for s in 1:ns, t in 1:ns]
                    P .+= δP
                    st.H[a] .+= 2 .* δP
                    st.y.m[a] .+= [β * π_[s] * rec[:dy][k][a, bs[s]] for s in 1:ns]
                end
                if :K in which
                    rK = [β * π_[s] * rec[:dK][k][a, bs[s]] for s in 1:ns]
                    u = [max(P[s, s], 0.0) for s in 1:ns]
                    su = sum(u); su > 0 ? (u ./= su) : (u .= 1 / ns)
                    R = sum(rK)
                    st.H[a] .+= rK * u' .+ u * rK' .- R .* (u * u')
                end
            end
            yfix_step!(st, g, d, a; ρ2 = ρ2r[a], kw...)
        end
        out.K[k] = Kt; out.C[k] = Ct; out.yY[k] = yt; out.H[k] = Ht
    end
    out
end

const DRIVEN = [(:Bps, (dir = false,)), (:BdP, (dir = true, which = (:P,))), (:BdK, (dir = true, which = (:K,))),
                (:Bd, (dir = true, which = (:P, :K))), (:Bt, (dir = false, tk3 = true)), (:Bdt, (dir = true, which = (:P, :K), tk3 = true))]

function driven_input(r, R, set; io = stdout, variants = DRIVEN)
    runs = load_dr(r.id, R, set)
    isempty(runs) && return nothing
    G = ld_groups(runs[1].L, runs[1].u); ρ2 = class_rho2(runs[1].L, runs[1].u)
    ser = Dict{Symbol,Vector{Any}}()
    cum = Dict{Symbol,Vector{Float64}}()
    for tr0 in runs
        tr = drview(tr0)
        m = measured(tr)
        push!(get!(ser, :yY, []), (m.t, m.yY)); push!(get!(ser, :C, []), (m.t, m.C)); push!(get!(ser, :K, []), (m.t, m.K))
        push!(get!(ser, :H, []), (m.t, m.K .+ 2 .* m.yY))
        k0 = findfirst(==(-1), m.t); k100 = findfirst(==(100), m.t)
        g0 = k0 === nothing ? NaN : mean(m.genic[max(1, k0 - 19):k0])
        ok = k0 !== nothing && k100 !== nothing
        if ok
            rg = k0+1:k100-1
            push!(get!(cum, :meas_loss, Float64[]), sum(2 .* m.beta[rg] .* m.C[rg]) / g0)
            push!(get!(cum, :meas_losssym, Float64[]), sum(2 .* m.beta[rg] .* m.Csym[rg]) / g0)
        end
        for (v, kw) in variants
            o = levelBdir(tr, G, ρ2; kw...)
            for (s, x) in ((:yY, o.yY), (:K, o.K), (:C, o.C), (:H, o.H))
                push!(get!(ser, Symbol(v, :_, s), []), (m.t, x))
            end
            ok && push!(get!(cum, Symbol(v, :_loss), Float64[]), sum(2 .* m.beta[k0+1:k100-1] .* o.C[k0+1:k100-1]) / g0)
        end
    end
    res = Dict{Symbol,Any}()
    for (s, v) in ser
        res[s] = Dict(w => winmean(v, w) for w in WINS)
    end
    for (s, v) in cum
        res[s] = isempty(v) ? NaN : mean(v)
    end
    haskey(res, :meas_loss) || (res[:meas_loss] = NaN; res[:meas_losssym] = NaN)
    @printf(io, "%-34s L=%d runs %d", r.id, runs[1].L, length(runs))
    for (v, _) in variants
        @printf(io, "  %s K %.2f/%.2f/%.2f", v, (res[Symbol(v, :_K)][w] / res[:K][w] for w in ((-60, -1), (50, 99), (100, 199)))...)
    end
    println(io)
    res
end

function summarize(all, variants; io = stdout, stats = (:yY, :H, :K))
    elig = filter(x -> isfinite(x[2][:meas_loss]) && x[2][:meas_loss] <= -0.05, all)
    @printf(io, "medians over %d inputs (%d eligible):\n", length(all), length(elig))
    for v in variants
        for s in stats
            @printf(io, "  %-5s %-3s ratio", v, s)
            for w in WINS
                u = [x[2][Symbol(v, :_, s)][w] / x[2][s][w] for x in all]
                @printf(io, "  %d..%d %.3f (|lr| %.3f)", w..., medf(u), mabs([lr(x[2][Symbol(v, :_, s)][w], x[2][s][w]) for x in all]))
            end
            println(io)
        end
        @printf(io, "  %-5s ΣC lr (eligible)", v)
        for w in WINS[3:5]
            u = [lr(x[2][Symbol(v, :_C)][w], x[2][:C][w]) for x in elig]
            @printf(io, "  %d..%d %+.3f/%.3f", w..., medf(u), mabs(u))
        end
        u = [lr(x[2][Symbol(v, :_loss)], x[2][:meas_losssym]) for x in elig]
        @printf(io, "  loss %+.3f/%.3f\n", medf(u), mabs(u))
    end
    elig
end

function dir_driven(set, R; io = stdout, variants = DRIVEN)
    rows = dr_rows(set)
    println(io, "\n== level B′* driven by the measured directional changes, set $set, R = $R ==")
    bufs = [IOBuffer() for _ in rows]
    resv = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        resv[i] = driven_input(rows[i], R, set; io = bufs[i], variants)
    end
    foreach(b -> print(io, String(take!(b))), bufs)
    all_ = [(rows[i], resv[i]) for i in eachindex(rows) if resv[i] !== nothing]
    summarize(all_, first.(variants); io)
    all_
end

# ---------------------------------------------------------------------------------------
# 4. level B‡ (validate/linked_dir_theory.jl): B′* with the directional terms closed

include(joinpath(@__DIR__, "linked_dir_theory.jl"))

function ddrivers(rec, k; avg = nothing)
    ks = avg === nothing ? (k:k) : avg
    ms(s) = Float64.(mean(rec[s][kk] for kk in ks))
    dir_drivers(ms(:V), ms(:V2), ms(:sk4), ms(:ss2v))
end

"""
Level B‡ (validate/linked_dir_theory.jl). A zeroth-order level B′* (FIXES[:Bps]) runs alongside and
drives the odd state Q^T and the zeroth-order T̄ used in the increments; level B‡'s own y/H state gets
the increments (`incr`: :P (with H), :m) and its T recursion the K3 source and the Σκ3 coupling
(`tsrc`). Start: B′*'s stationary state for the drivers averaged over t = −60..−41, Q^T = 0, T = 0.
Returns per generation Σ K, Σ C, Σ E[yY], Σ H of level B‡ and Σ Q^T.
"""
function levelBdag(tr, G::LDGroups, ρ2; incr = (:P, :m), tsrc = true, kw = FIXES[:Bps])
    rec = tr.rec
    n = length(rec[:gen])
    kw0 = merge((; kw...), (hsel = :none,))
    d0 = drivers(rec, 1; avg = PRE_AVG)
    NA = size(G.n, 1)
    act = [a for a in 1:NA if G.nb[a] > 0]
    idx = Dict(a => findall(>(0), G.n[a, :]) for a in act)
    Gr = Dict(a => LDGroups(G.L, G.n[:, idx[a]], G.ρ[:, idx[a]], G.nb, G.bofs[idx[a]]) for a in act)
    sr = Dict(a => LSState(Gr[a]) for a in act)
    ρ2r = Dict(a => ρ2[:, idx[a]] for a in act)
    for a in act
        g = Gr[a]; st = sr[a]
        c = class_c(g, d0, a); Σ = class_cov(g, d0, a)
        for it in 1:40000
            chk = it % 25 == 0
            chk && (P0 = copy(st.y.P[a]); H0 = copy(st.H[a]))
            yfix_step!(st, g, d0, a; c, Σ, ρ2 = ρ2r[a], kw0...)
            chk || continue
            sc = max(maximum(abs, st.y.P[a]), 1e-300)
            max(maximum(abs, st.y.P[a] .- P0), maximum(abs, st.H[a] .- H0)) <= 1e-9 * sc && break
        end
    end
    s0 = Dict(a => deepcopy(sr[a]) for a in act)          # zeroth order (B′*)
    T0 = Dict(a => zeros(length(idx[a])) for a in act)
    Tr = Dict(a => zeros(length(idx[a])) for a in act)
    QT = Dict(a => zeros(length(idx[a]), length(idx[a])) for a in act)
    out = (K = fill(NaN, n), C = fill(NaN, n), yY = fill(NaN, n), H = fill(NaN, n), QT = fill(NaN, n))
    for k in 1:n
        d = drivers(rec, k); dd = ddrivers(rec, k)
        Kt = 0.0; Ct = 0.0; yt = 0.0; Ht = 0.0; qt = 0.0
        for a in act
            g = Gr[a]; st = sr[a]; z = s0[a]
            Ks = source_K(st, a); Kz = source_K(z, a)
            Kt += sum(Ks); Ct += sum(Tr[a]); yt += sum(st.y.P[a]); Ht += sum(st.H[a]); qt += sum(QT[a])
            inc = qt_increments(QT[a], T0[a], d)
            if tsrc
                T_step!(Tr[a], g, d, a, Ks .+ k3_source(st.y.m[a], dd, a); mode = :derived, k3 = d.k3[a])
                T_step!(T0[a], g, d, a, Kz .+ k3_source(z.y.m[a], dd, a); mode = :derived, k3 = d.k3[a])
            else
                T_step!(Tr[a], g, d, a, Ks; mode = :derived)
                T_step!(T0[a], g, d, a, Kz; mode = :derived)
            end
            qt_step!(QT[a], z.y.P[a], z.y.m[a], Kz, g, d, dd, a)
            yfix_step!(z, g, d, a; ρ2 = ρ2r[a], kw...)
            :P in incr && (st.y.P[a] .+= inc.δP; st.H[a] .+= inc.δH)
            :m in incr && (st.y.m[a] .+= inc.δm)
            yfix_step!(st, g, d, a; ρ2 = ρ2r[a], kw...)
        end
        out.K[k] = Kt; out.C[k] = Ct; out.yY[k] = yt; out.H[k] = Ht; out.QT[k] = qt
    end
    out
end

const DAG = [(:Bps, nothing), (:Bdag, (;)), (:BdT, (incr = (),))]

function dag_input(r, R, set; io = stdout, variants = DAG)
    runs = load_dr(r.id, R, set)
    isempty(runs) && return nothing
    G = ld_groups(runs[1].L, runs[1].u); ρ2 = class_rho2(runs[1].L, runs[1].u)
    ser = Dict{Symbol,Vector{Any}}()
    cum = Dict{Symbol,Vector{Float64}}()
    for tr0 in runs
        tr = drview(tr0)
        m = measured(tr)
        push!(get!(ser, :yY, []), (m.t, m.yY)); push!(get!(ser, :C, []), (m.t, m.C)); push!(get!(ser, :K, []), (m.t, m.K))
        push!(get!(ser, :H, []), (m.t, m.K .+ 2 .* m.yY))
        push!(get!(ser, :QT, []), (m.t, tot.(tr.rec[:yC])))
        k0 = findfirst(==(-1), m.t); k100 = findfirst(==(100), m.t)
        g0 = k0 === nothing ? NaN : mean(m.genic[max(1, k0 - 19):k0])
        ok = k0 !== nothing && k100 !== nothing
        if ok
            rg = k0+1:k100-1
            push!(get!(cum, :meas_loss, Float64[]), sum(2 .* m.beta[rg] .* m.C[rg]) / g0)
            push!(get!(cum, :meas_losssym, Float64[]), sum(2 .* m.beta[rg] .* m.Csym[rg]) / g0)
        end
        for (v, kw) in variants
            o = kw === nothing ? levelBdir(tr, G, ρ2; dir = false) : levelBdag(tr, G, ρ2; kw...)
            for (s, x) in ((:yY, o.yY), (:K, o.K), (:C, o.C), (:H, o.H))
                push!(get!(ser, Symbol(v, :_, s), []), (m.t, x))
            end
            haskey(o, :QT) && push!(get!(ser, Symbol(v, :_QT), []), (m.t, o.QT))
            ok && push!(get!(cum, Symbol(v, :_loss), Float64[]), sum(2 .* m.beta[k0+1:k100-1] .* o.C[k0+1:k100-1]) / g0)
        end
        if any(v -> v[1] === :A, variants) || true
            push!(get!(ser, :A_derived, []), (m.t, levelA(tr, G; mode = :derived)))
            push!(get!(ser, :A_nosel, []), (m.t, levelA(tr, G; mode = :nosel)))
            push!(get!(ser, :A_spectrum, []), (m.t, levelA(tr, G; mode = :spectrum)))
        end
    end
    res = Dict{Symbol,Any}()
    for (s, v) in ser
        res[s] = Dict(w => winmean(v, w) for w in WINS)
    end
    for (s, v) in cum
        res[s] = isempty(v) ? NaN : mean(v)
    end
    haskey(res, :meas_loss) || (res[:meas_loss] = NaN; res[:meas_losssym] = NaN)
    @printf(io, "%-34s L=%d runs %d", r.id, runs[1].L, length(runs))
    for (v, _) in variants
        @printf(io, "  %s yY %.2f/%.2f H %.2f/%.2f K %.2f/%.2f", v, (res[Symbol(v, :_, s)][w] / res[s][w] for s in (:yY, :H, :K) for w in ((-60, -1), (50, 99)))...)
    end
    println(io)
    res
end

function dag_run(set, R; io = stdout, variants = DAG)
    rows = dr_rows(set)
    println(io, "\n== level B‡ (closed directional terms), set $set, R = $R ==")
    bufs = [IOBuffer() for _ in rows]
    resv = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        resv[i] = dag_input(rows[i], R, set; io = bufs[i], variants)
    end
    foreach(b -> print(io, String(take!(b))), bufs)
    all_ = [(rows[i], resv[i]) for i in eachindex(rows) if resv[i] !== nothing]
    elig = summarize(all_, first.(variants); io)
    @printf(io, "  P1 (level A), eligible: log(predicted ΣC / measured ΣC) median / median |.|:\n")
    for md in AMODES, w in WINS[3:5]
        u = [lr(x[2][Symbol(:A_, md)][w], x[2][:C][w]) for x in elig]
        @printf(io, "    %-9s %d..%d %+.3f / %.3f (n %d)\n", md, w..., medf(u), mabs(u), count(isfinite, u))
    end
    all_
end

# ---------------------------------------------------------------------------------------
# exploratory: candidate decompositions of the measured directional change of K (not a verdict)

function k_candidates(tr)
    r = tr.rec
    t = Int.(r[:gen]) .- MOVE
    n = length(t)
    out = Dict{Symbol,Vector{Float64}}(k => zeros(n) for k in (:dK, :w, :E5c, :prodA, :prodB, :vchg, :K, :dH, :H, :hb))
    for k in 1:n
        β = r[:beta][k]
        V = Float64.(r[:V][k]); sk4 = Float64.(r[:sk4][k]); k3 = Float64.(r[:k3][k])
        wb = [V[a] > 0 ? sk4[a] / V[a] : 0.0 for a in eachindex(V)]
        sb = [V[a] > 0 ? k3[a] / V[a] : 0.0 for a in eachindex(V)]
        J = Float64.(r[:J][k]); T = Float64.(r[:T][k]); K = Float64.(r[:K][k])
        NB = size(J, 1)
        out[:dK][k] = β * tot(r[:dK][k]); out[:K][k] = tot(r[:K][k])
        out[:dH][k] = β * tot(r[:dH][k]); out[:H][k] = tot(r[:H][k])
        out[:w][k] = β * sum(wb[a] * J[a, b] + wb[b] * T[a, b] for a in 1:NB, b in 1:NB)
        out[:E5c][k] = β * tot(r[:E5c][k])
        out[:prodA][k] = -β * (6tot(r[:TY][k]) + 4tot(r[:JY][k]))
        out[:prodB][k] = -6β * (tot(r[:TY][k]) + tot(r[:JY][k]))
        out[:vchg][k] = β * sum(K[a, b] * (sb[a] + sb[b]) for a in 1:NB, b in 1:NB)
        out[:hb][k] = out[:w][k] + out[:prodB][k]          # biallelic projection: w-terms + product terms
    end
    (; t, out)
end

function k_explore(set, R; io = stdout)
    rows = dr_rows(set)
    println(io, "\n== exploratory: measured β Σ δK against candidate decompositions, set $set (÷ |window mean Σ K|) ==")
    cands = (:dK, :w, :E5c, :prodA, :prodB, :vchg, :hb)
    agg = Dict{Tuple{Symbol,Tuple{Int,Int}},Vector{Float64}}()
    for rw in rows
        runs = load_dr(rw.id, R, set)
        isempty(runs) && continue
        S = [k_candidates(tr) for tr in runs]
        for wdw in ((0, 24), (25, 49), (50, 99), (100, 199))
            den = abs(winmean([(s.t, s.out[:K]) for s in S], wdw))
            for c in cands
                push!(get!(agg, (c, wdw), Float64[]), winmean([(s.t, s.out[c]) for s in S], wdw) / den)
            end
        end
    end
    for wdw in ((0, 24), (25, 49), (50, 99), (100, 199))
        @printf(io, "  %d..%d:", wdw...)
        for c in cands
            @printf(io, " %s med %+.4f", c, medf(agg[(c, wdw)]))
        end
        x = agg[(:dK, wdw)]
        for c in cands[2:end]
            yv = agg[(c, wdw)]
            ok = isfinite.(x) .& isfinite.(yv)
            @printf(io, " | cor(dK,%s) %+.2f", c, sum(ok) > 3 ? cor(x[ok], yv[ok]) : NaN)
        end
        println(io)
    end
    agg
end

# exploratory: the second moments R of u = T + J against the biallelic closure ⟨s²⟩_l P (not a verdict)
function r_explore(set, R; io = stdout)
    rows = dr_rows(set)
    println(io, "\n== exploratory: R = Σ_l E[u uᵀ] (u = T + J), set $set: totals over groups and bins, window means over recorded generations ==")
    agg = Dict{Tuple{Symbol,Tuple{Int,Int}},Vector{Float64}}()
    ks = (:RTT_P, :RJJ_P, :RTJ_P, :s2P_P, :QT_P, :QJ_P)
    for rw in rows
        runs = load_dr(rw.id, R, set)
        isempty(runs) && continue
        ser = Dict(k => [] for k in (:P, :RTT, :RJJ, :RTJ, :s2P, :QT, :QJ))
        for tr in runs
            r = tr.rec
            tm = r[:mgen] .- MOVE
            gi = Dict(g => k for (k, g) in enumerate(Int.(r[:gen])))
            vals = Dict(k => Float64[] for k in keys(ser))
            for (q, g) in enumerate(r[:mgen])
                k = gi[g]
                V = Float64.(r[:V][k]); s2 = Float64.(r[:ss2v][k])
                Pm = r[:P][q]
                push!(vals[:P], sum(Float64, Pm))
                for s in (:RTT, :RJJ, :RTJ, :QT, :QJ)
                    push!(vals[s], sum(Float64, r[s][q]))
                end
                push!(vals[:s2P], sum((V[a] > 0 ? s2[a] / V[a] : 0.0) * sum(Float64, view(Pm, :, :, a)) for a in eachindex(V)))
            end
            for s in keys(ser)
                push!(ser[s], (tm, vals[s]))
            end
        end
        for wdw in DTWINS
            P = winmean(ser[:P], wdw)
            for (k, s) in zip(ks, (:RTT, :RJJ, :RTJ, :s2P, :QT, :QJ))
                push!(get!(agg, (k, wdw), Float64[]), winmean(ser[s], wdw) / P)
            end
        end
    end
    for k in ks
        @printf(io, "  %-6s", k)
        for wdw in DTWINS
            @printf(io, "  %d..%d med %+.3e", wdw..., medf(agg[(k, wdw)]))
        end
        println(io)
    end
    agg
end

function dr_main(mode, args)
    if mode == "terms"
        dir_terms(args[1], parse(Int, args[2]))
    elseif mode == "dev"
        t = @elapsed dag_run(args[1], parse(Int, args[2]))
        @printf("(analysis %.0f s wall)\n", t)
    elseif mode == "kexplore"
        k_explore(args[1], parse(Int, args[2]))
    elseif mode == "rexplore"
        r_explore(args[1], parse(Int, args[2]))
    elseif mode == "driven"
        t = @elapsed dir_driven(args[1], parse(Int, args[2]))
        @printf("(analysis %.0f s wall)\n", t)
    elseif isdefined(Main, :dr_main2)
        dr_main2(mode, args)
    else
        error("unknown mode $mode")
    end
end
