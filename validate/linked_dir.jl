# The directional (first order in the selection gradient b) terms of the allele–background
# association closure under a moving optimum (follow-up of validate/linked_ld.jl; criteria
# validate/linked_dir_criteria.md, committed before any comparison on the fresh set J; running log
# validate/linked_dir.log; theory validate/linked_dir_theory.jl; algebra
# validate/symbolic/check_linked_dir.jl).
#
# Notation (chromosome level, adults; locus l, class or bin s of partners, B_s the partner sum,
# B = G_chr − x_l, d_l = x_l − x̄_l, as linked_decay/linked_source/linked_ld):
#   y_s = Cov(B_s, x_l), T_s = κ(B_s, l, l), J_s = κ(B_s, B, l), u_s = T_s + J_s = κ(B_s, l, G_chr)
#   H_s = Cov(B_s B, d_l²), K_s = κ(B_s, B, l, l) = H_s − 2 y_s Y_l.
# Under weights 1 + ε h (selection's linear part, h = G − Ḡ, actual change = β × the derivative) the
# exact first-order changes are
#   Δy_s = β u_s,   ΔP_st = β Σ_l (u_s y_t + y_s u_t) = β (Q_st + Q_ts),   Q_st = Σ_l u_s y_t,
#   ΔK_s = β κ₅(B_s, B, l, l, G),   ΔH_s = ΔK_s + 2 Δ(y_s Y_l).
# Linked_ld's level B′* drops all three (they vanish in expectation before the move, by the symmetry
# G → −G of stabilising selection around the optimum).
#
# Instrumented runs (`dir_run`, reading the World only; `check` verifies bitwise identity with
# `simulate` and agreement with linked_source's statistics on the same run) measure every
# generation from MOVE − 60 to MOVE + 199, per group a of loci l and partner bin b (sums over l ∈ a):
# linked_source's y, T, K, K3, yY, H, J, the exact directional derivatives of y, T, J, H, K and
# Σ y_b Y, the bin matrices Σ_l y yᵀ, Σ_l T yᵀ, Σ_l J yᵀ, the directional Σ_l (δy yᵀ + y δyᵀ), the
# second moments Σ_l T Tᵀ, J Jᵀ, T Jᵀ, the spectrum-shape sums Σ κ3²/v, Σ κ4, and the pairwise
# biallelic-model (PBM) pieces Σ_{k∈b} s_k D_kl, Σ_{k∈b} (s_k² − 2v_k) D_kl (s = κ3/v).
#
# Run from the package root (vignettes environment):
#   julia --project=vignettes validate/linked_dir.jl check                # bitwise identity, agreement
#   julia -t 6 --project=vignettes validate/linked_dir.jl kernel SET R    # instrumented runs (resumable)
#   julia --project=vignettes validate/linked_dir.jl terms SET R          # size of the directional terms
#   julia -t 6 --project=vignettes validate/linked_dir.jl driven SET R    # B′* with the measured terms
#   julia -t 6 --project=vignettes validate/linked_dir.jl dev SET R       # level B‡ (development)
#   julia --project=vignettes validate/linked_dir.jl jrows                # the confirmation inputs J
#   julia -t 6 --project=vignettes validate/linked_dir.jl confirm J R     # pre-specified verdicts
# SET: D (linked_decay's development set and seeds), H (linked_ld's H′ and seeds), J (fresh).
# Raw records: validate/linked_dir_raw/ (git-ignored).

# 2026-11-19: `ne = :falconer`, `solver = :iteration` and `QTLSim.burnin_start_vg` below keep the effective size, the solver and
# the starting variance this script was run with (the package's defaults are now `:breeders` and `:bisection`; NEWS.md), so that
# it still reproduces its log.
using LinearAlgebra
include(joinpath(@__DIR__, "linked_ld.jl"))

const DR_RAW = joinpath(@__DIR__, "linked_dir_raw")
const DR_MSTRIDE = 2            # bin matrices every DR_MSTRIDE generations (from LD_FROM)

# ---------------------------------------------------------------------------------------
# measurement

mutable struct DRBuf
    L::Int; NB::Int
    binl::Vector{Int}
    xbar::Vector{Float64}; v::Vector{Float64}; k3::Vector{Float64}; m4::Vector{Float64}
    sl::Vector{Float64}; wl::Vector{Float64}
    XT::Matrix{Float64}                      # (L, M) deviations
    X::Matrix{Float64}                       # (M, L)
    E::Matrix{Float64}; Es::Matrix{Float64}; Ew::Matrix{Float64}   # (NB, M) centred bin sums: d, s·d, w·d
    Gc::Vector{Float64}; dG::Vector{Float64}; h::Vector{Float64}
    acc::Matrix{Float64}                     # (13, NB)
    x::Vector{Float64}
end

function DRBuf(L, Nmax; NB = LD_NB)
    M = 2Nmax
    DRBuf(L, NB, [ld_bin(l - 1, L; nb = NB) for l in 1:L], zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L),
          zeros(L, M), zeros(M, L), zeros(NB, M), zeros(NB, M), zeros(NB, M), zeros(M), zeros(M), zeros(M),
          zeros(13, NB), zeros(L))
end

"""
Association statistics of the adults' 2N chromosomes and their exact first-order changes under the
weights 1 + ε h, h = G − Ḡ (individual genotypes, centred), as derivatives at ε = 0 (multiply by β for
the change by the linear part of selection). `mats` adds the bin matrices.
"""
function dir_measure!(B::DRBuf, w; mats = true)
    L, NB, N = B.L, B.NB, w.N
    M = 2N
    x = B.x
    fill!(B.xbar, 0)
    gi = view(w.genotype, 1:N)
    gm = mean(gi)
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        ldchr!(x, w, c, i)
        g = 0.0
        @inbounds for l in 1:L
            g += x[l]; B.xbar[l] += x[l]; B.XT[l, m] = x[l]
        end
        B.Gc[m] = g
        B.h[m] = gi[i] - gm
    end
    B.xbar ./= M
    Gbar = sum(view(B.Gc, 1:M)) / M
    fill!(B.v, 0); fill!(B.k3, 0); fill!(B.m4, 0)
    @inbounds for m in 1:M
        for l in 1:L
            d = B.XT[l, m] - B.xbar[l]
            B.XT[l, m] = d
            d2 = d * d
            B.v[l] += d2; B.k3[l] += d2 * d; B.m4[l] += d2 * d2
        end
    end
    B.v ./= M; B.k3 ./= M; B.m4 ./= M
    @inbounds for l in 1:L
        vl = B.v[l]
        B.sl[l] = vl > 0 ? B.k3[l] / vl : 0.0
        B.wl[l] = B.sl[l]^2 - 2vl
    end
    X = view(B.X, 1:M, :)
    permutedims!(X, view(B.XT, :, 1:M), (2, 1))
    E, Es, Ew = B.E, B.Es, B.Ew
    hE = zeros(NB)
    @inbounds for m in 1:M
        B.dG[m] = B.Gc[m] - Gbar
        for b in 1:NB
            E[b, m] = 0.0; Es[b, m] = 0.0; Ew[b, m] = 0.0
        end
        for l in 1:L
            d = B.XT[l, m]; b = B.binl[l]
            E[b, m] += d; Es[b, m] += B.sl[l] * d; Ew[b, m] += B.wl[l] * d
        end
        for b in 1:NB
            hE[b] += B.h[m] * E[b, m]
        end
    end
    hE ./= M
    acc = B.acc
    ns = NB
    ks = (:y, :T, :K, :K3, :yY, :H, :J, :dy, :dT, :dJ, :dH, :dK, :dYY, :Jp, :Wp, :E5c, :sJ, :sy, :TY, :JY, :yC, :yJ, :CBB)
    A = Dict(k => zeros(NB, NB) for k in ks)
    MT = (:P, :QT, :QJ, :RTT, :RJJ, :RTJ, :dP)
    Mm = Dict(k => [zeros(NB, NB) for _ in 1:NB] for k in MT)
    V = zeros(NB); V2 = zeros(NB); k3g = zeros(NB); g = zeros(NB); ss2v = zeros(NB); sk4 = zeros(NB)
    dk3 = zeros(NB); dvg = zeros(NB); dC = zeros(NB)
    sumC = 0.0; sk3 = 0.0; sk3sq = 0.0; sk3C = 0.0; genic = 0.0
    yv = zeros(NB); Tv = zeros(NB); Jv = zeros(NB); dyv = zeros(NB)
    h = B.h; dGv = B.dG
    @inbounds for l in 1:L
        fill!(acc, 0)
        vl = B.v[l]; ob = B.binl[l]; sl = B.sl[l]; wl = B.wl[l]
        sY = 0.0; sC = 0.0; hd = 0.0; hd2 = 0.0; hdB = 0.0; hdBd = 0.0; hdBd2 = 0.0; hd3 = 0.0
        o = zeros(13)
        for m in 1:M
            d = X[m, l]; hm = h[m]
            d2 = d * d; d3 = d2 * d
            dB = dGv[m] - d
            dBd = dB * d; dBd2 = dBd * d
            hd_ = hm * d; hd2_ = hm * d2; hdBd_ = hm * dBd; hdBd2_ = hm * dBd2; hdB_ = hm * dB
            @simd for b in 1:NB
                eb = E[b, m]
                acc[1, b] += eb * d
                acc[2, b] += eb * d2
                acc[3, b] += eb * dBd2
                acc[4, b] += eb * d3
                acc[5, b] += eb * dB
                acc[6, b] += eb * dBd
                acc[7, b] += eb * hd_
                acc[8, b] += eb * hd2_
                acc[9, b] += eb * hdBd_
                acc[10, b] += eb * hdBd2_
                acc[11, b] += eb * hdB_
                acc[12, b] += Es[b, m] * d
                acc[13, b] += Ew[b, m] * d
            end
            # own bin: e_ob − d, es_ob − s_l d, ew_ob − w_l d
            o[1] -= d * d; o[2] -= d * d2; o[3] -= d * dBd2; o[4] -= d * d3; o[5] -= d * dB; o[6] -= d * dBd
            o[7] -= d * hd_; o[8] -= d * hd2_; o[9] -= d * hdBd_; o[10] -= d * hdBd2_; o[11] -= d * hdB_
            o[12] -= sl * d * d; o[13] -= wl * d * d
            sY += dBd; sC += dBd2
            hd += hd_; hd2 += hd2_; hdB += hdB_; hdBd += hdBd_; hdBd2 += hdBd2_; hd3 += hm * d3
        end
        for k in 1:13
            acc[k, ob] += o[k]
        end
        acc ./= M
        Y = sY / M; Cl = sC / M
        δx = hd / M; δv = hd2 / M; δB = hdB / M; δY = hdBd / M; δC0 = hdBd2 / M; δm3 = hd3 / M
        a = ob
        V[a] += vl; V2[a] += vl^2; k3g[a] += B.k3[l]; g[a] += vl
        vl > 0 && (ss2v[a] += B.k3[l]^2 / vl)
        sk4[a] += B.m4[l] - 3vl^2
        dk3[a] += δm3 - 3δx * vl
        dvg[a] += δv
        δYl = δY
        Jtot = 0.0
        for b in 1:NB
            Jtot += acc[6, b]
        end
        dCl = δC0 - δB * vl - 2δx * Y        # δ κ(B, l, l)
        dC[a] += dCl
        for b in 1:NB
            y = acc[1, b]; T = acc[2, b]; CBB = acc[5, b]; J = acc[6, b]
            H = acc[3, b] - CBB * vl
            K = H - 2y * Y
            K3 = acc[4, b] - 3y * vl
            δc = hE[b] - (b == ob ? δx : 0.0)
            δy = acc[7, b]
            δT = acc[8, b] - δc * vl - 2δx * y
            δJ = acc[9, b] - δc * Y - δB * y - δx * CBB
            δCBB = acc[11, b]
            δH = acc[10, b] - δc * Cl - δB * T - 2δx * J - δCBB * vl - CBB * δv
            δK = δH - 2δy * Y - 2y * δYl
            Jp = acc[12, b]; Wp = acc[13, b]
            A[:y][a, b] += y; A[:T][a, b] += T; A[:K][a, b] += K; A[:K3][a, b] += K3; A[:yY][a, b] += y * Y
            A[:H][a, b] += H; A[:J][a, b] += J; A[:CBB][a, b] += CBB
            A[:dy][a, b] += δy; A[:dT][a, b] += δT; A[:dJ][a, b] += δJ; A[:dH][a, b] += δH; A[:dK][a, b] += δK
            A[:dYY][a, b] += δy * Y + y * δYl
            A[:Jp][a, b] += Jp; A[:Wp][a, b] += Wp; A[:E5c][a, b] += sl * Wp + wl * Jp
            A[:sJ][a, b] += sl * J; A[:sy][a, b] += sl * y
            A[:TY][a, b] += T * Y; A[:JY][a, b] += J * Y; A[:yC][a, b] += y * Cl; A[:yJ][a, b] += y * Jtot
            yv[b] = y; Tv[b] = T; Jv[b] = J; dyv[b] = δy
        end
        if mats
            Pm = Mm[:P][a]; QT = Mm[:QT][a]; QJ = Mm[:QJ][a]; RTT = Mm[:RTT][a]; RJJ = Mm[:RJJ][a]; RTJ = Mm[:RTJ][a]; dP = Mm[:dP][a]
            for b2 in 1:NB, b in 1:NB
                Pm[b, b2] += yv[b] * yv[b2]; QT[b, b2] += Tv[b] * yv[b2]; QJ[b, b2] += Jv[b] * yv[b2]
                RTT[b, b2] += Tv[b] * Tv[b2]; RJJ[b, b2] += Jv[b] * Jv[b2]; RTJ[b, b2] += Tv[b] * Jv[b2]
                dP[b, b2] += dyv[b] * yv[b2] + yv[b] * dyv[b2]
            end
        end
        sumC += Cl; sk3 += B.k3[l]; sk3sq += B.k3[l]^2; sk3C += B.k3[l] * Cl; genic += 2vl
    end
    Cbb = zeros(NB, NB); cb = zeros(NB)
    VGchr = 0.0
    Bb = zeros(NB)
    @inbounds for m in 1:M
        dG = B.dG[m]
        VGchr += dG^2
        for b in 1:NB
            cb[b] += E[b, m] * dG
            for b2 in 1:NB
                Cbb[b, b2] += E[b, m] * E[b2, m]
            end
        end
    end
    Cbb ./= M; cb ./= M; VGchr /= M
    (; A..., mats = mats ? Mm : nothing, V, V2, k3 = k3g, g, ss2v, sk4, dk3, dv = dvg, dC, c = cb, Cbb,
     sumC, sk3, sk3sq, sk3C, genic, VGchr)
end

const DR_SCAL = [:gen, :N, :nsurv, :theta, :gbar, :VG, :VGchr, :Sc, :beta, :beta2, :bt, :genic, :sumC, :sk3, :sk3sq, :sk3C,
                 :VE, :omega2, :Vm]
const DR_MATS = [:y, :T, :K, :K3, :yY, :H, :J, :dy, :dT, :dJ, :dH, :dK, :dYY, :Jp, :Wp, :E5c, :sJ, :sy, :TY, :JY, :yC, :yJ, :CBB, :Cbb]
const DR_VECS = [:V, :V2, :k3, :g, :c, :ss2v, :sk4, :dk3, :dv, :dC]
const DR_BMATS = [:P, :QT, :QJ, :RTT, :RJJ, :RTJ, :dP]

function dir_run(p, seed; from = LD_FROM, to = LD_TO, mstride = DR_MSTRIDE)
    BLAS.set_num_threads(1)
    q = kernel_params(p)
    q.Selection === :viability || error("viability only")
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    L = q.NoQTL
    Bf = DRBuf(L, q.PopSize)
    VE = QTLSim.envVar(q)
    rec = Dict{Symbol,Any}(s => Float64[] for s in DR_SCAL)
    for s in DR_MATS
        rec[s] = Matrix{Float32}[]
    end
    for s in DR_VECS
        rec[s] = Vector{Float32}[]
    end
    rec[:mgen] = Int[]
    for s in DR_BMATS
        rec[s] = Array{Float32,3}[]
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
            domat = (it - from) % mstride == 0
            ms = dir_measure!(Bf, w; mats = domat)
            Sc = ω2 + VE + ms.VGchr
            vals = (gen = it, N = N, nsurv = w.nsurv, theta = w.θ, gbar = gb, VG = var(g), VGchr = ms.VGchr, Sc = Sc,
                    beta = β, beta2 = β2, bt = (w.θ - gb) / Sc, genic = ms.genic, sumC = ms.sumC, sk3 = ms.sk3, sk3sq = ms.sk3sq,
                    sk3C = ms.sk3C, VE = VE, omega2 = ω2, Vm = QTLSim.muteVar(q))
            for s in DR_SCAL
                push!(rec[s], Float64(vals[s]))
            end
            for s in DR_MATS
                push!(rec[s], Float32.(getfield(ms, s)))
            end
            for s in DR_VECS
                push!(rec[s], Float32.(getfield(ms, s)))
            end
            if domat
                push!(rec[:mgen], it)
                for s in DR_BMATS
                    push!(rec[s], Float32.(cat(ms.mats[s]...; dims = 3)))
                end
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

dr_file(id, r, set) = joinpath(DR_RAW, replace(id, "/" => "__", ".txt" => "") * "@dir" * set * "#$r.jls")

"""
Check: (i) V_G bitwise identical to `simulate` (the instrumentation draws nothing); (ii) y, T, K, H,
J, yY equal linked_source's `ls_run` on the same generations; (iii) the directional derivatives
equal finite differences of the statistics under weights 1 ± ε h (linked_source's `ls_measure!`).
"""
function dir_check(; id = "QT1000-Possion3.0/Input-09697.txt", j = 89, r = 1)
    p = readp(id)
    t = @elapsed tr = dir_run(p, KSEED + 100j + r; to = LD_FROM + 6)
    o = simulate(kernel_params(p); seed = KSEED + 100j + r, init = :equilibrium, burnin = KBURN)
    k = [findfirst(==(g), o.gen) for g in LD_FROM:LD_FROM+length(tr.vgser)-1]
    same = all(o.geno_var[k] .== tr.vgser)
    t2 = ls_run(p, KSEED + 100j + r; to = LD_FROM + 6, decompose = false)
    agree = true
    for s in (:y, :T, :K, :H, :J, :yY, :K3, :TY, :JY, :yC, :yJ), kk in 1:length(tr.vgser)
        a1 = Float64.(tr.rec[s][kk]); a2 = Float64.(t2.rec[s][kk])
        ok = maximum(abs.(a1 .- a2)) <= 1e-4 * max(maximum(abs.(a2)), 1e-30)
        ok || println("  mismatch $s gen $kk: ", maximum(abs.(a1 .- a2)), " vs ", maximum(abs.(a2)))
        agree &= ok
    end
    # finite differences of ls_measure! under weights 1 ± ε h at one generation of a fresh short run
    q = kernel_params(p)
    rng = MersenneTwister(KSEED + 100j + r)
    w = equilibrium_world(q, rng; burnin = KBURN)
    QTLSim.develop!(w); QTLSim.select!(w, rng)
    N = w.N
    hh = [w.genotype[i] - mean(view(w.genotype, 1:N)) for i in 1:N]
    ε = 1e-4 / maximum(abs, hh)
    Lb = LSBuf(q.NoQTL, q.PopSize, q.RecombPoissonMean)
    mp = ls_measure!(Lb, w, 1 .+ ε .* hh); mmn = ls_measure!(Lb, w, 1 .- ε .* hh)
    md = dir_measure!(DRBuf(q.NoQTL, q.PopSize), w; mats = true)
    fd = true
    for (s, sd) in ((:y, :dy), (:T, :dT), (:J, :dJ), (:H, :dH), (:K, :dK), (:yY, :dYY))
        num = (Float64.(getfield(mp, s)) .- Float64.(getfield(mmn, s))) ./ (2ε)
        an = Float64.(getfield(md, sd))
        err = maximum(abs.(num .- an)) / max(maximum(abs.(an)), 1e-300)
        @printf("  derivative %-3s vs finite difference: max rel. error %.2e\n", s, err)
        fd &= err <= 1e-4
    end
    # bin-matrix row sums equal the group × bin statistics
    rs = true
    for (s, sm) in ((:yY, :P), (:TY, :QT), (:JY, :QJ), (:dYY, :dP))
        for a in 1:LD_NB
            r1 = vec(sum(md.mats[sm][a]; dims = 2)); r0 = vec(getfield(md, s)[a, :])
            sm === :dP && (r1 = vec(sum(md.mats[sm][a]; dims = 2)) ./ 1)   # Σ_b' (δy_b y_b' + y_b δy_b') = δy_b Y + y_b δY
            rs &= maximum(abs.(r1 .- r0)) <= 1e-8 * max(maximum(abs.(r0)), 1e-30)
        end
    end
    @printf("check %s: %d generations, V_G bitwise identical to simulate: %s; statistics as ls_run: %s; derivatives as finite differences: %s; bin-matrix rows: %s (%.1f s)\n",
            id, length(tr.vgser), same, agree, fd, rs, t)
    (same && agree && fd && rs) || error("check failed")
end

dr_rows(set) = set == "J" ? jrows() : ll_rows(set)
dr_seed(set, j, r) = set == "J" ? JSEED + 100j + r : ll_seed(set, j, r)

function dr_kernel(set, R)
    mkpath(DR_RAW)
    rows = dr_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(dr_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) directional-instrumented runs for set $set, R = $R, $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time(); cpu = Threads.Atomic{Float64}(0.0)
    Threads.@threads :greedy for (id, j, k) in jobs
        t1 = time()
        tr = dir_run(readp(id), dr_seed(set, j, k))
        serialize(dr_file(id, k, set), tr)
        Threads.atomic_add!(cpu, time() - t1)
        d = Threads.atomic_add!(done, 1) + 1
        @printf("  %d/%d %s #%d  %.0f s\n", d, length(jobs), id, k, time() - t0)
    end
    @printf("done: %d runs, %.0f s wall, %.0f thread-s\n", length(jobs), time() - t0, cpu[])
end

# ---------------------------------------------------------------------------------------
# J: the fresh confirmation set (criteria validate/linked_dir_criteria.md)

const JSEED = 60261116000                         # J seeds: JSEED + 100j + r
const JPER = 3                                    # inputs per stratum

"""
J: inputs of the published sweep never used by the spectrum-diagnosis family of studies (not in
validate/spectrum_diagnosis_inputs.csv: none of D, F, G, C, Dctl, Cctl, E6′) and not in linked_ld's
H′, screened on parameters only (k²/V_G ≥ 0.01 with V_G = `house_of_gauss_vg`, k/k_c ≤ 1.4 with
k_c = `drift_critical_rate`). Per stratum (NoQTL × u) the files are taken in the order of a
MersenneTwister(20261118) permutation until `JPER` pass. j = 2000 + position (seeds distinct from
every earlier set's).
"""
function j_rows(; verbose = false)
    used = Set(r.id for r in read_inputs_csv())
    hset = Set(r.id for r in hrows())
    rng = MersenneTwister(20261118)
    out = NamedTuple[]
    j = 2000
    for s in sort(filter(x -> startswith(x, "QT"), readdir(SIMDIR)))
        files = sort(filter(f -> endswith(f, ".txt"), readdir(joinpath(SIMDIR, s))))
        npass = 0; ntried = 0
        for f in files[randperm(rng, length(files))]
            id = s * "/" * f
            (id in used || id in hset) && continue
            p = readp(id)
            vg = house_of_gauss_vg(p; ne = :falconer); kc = drift_critical_rate(p; ne = :falconer, solver = :iteration).k
            idx = p.MoveFitnessSD^2 / vg; ratio = p.MoveFitnessSD / kc
            ntried += 1
            if idx >= 0.01 && ratio <= 1.4
                j += 1; npass += 1
                push!(out, (j = j, id = id, set = "J", f = NaN, P = NaN, idx = idx, ratio = ratio))
                npass == JPER && break
            end
        end
        verbose && println("  $s: $npass pass of $ntried screened")
    end
    out
end
const _JROWS = Ref{Any}(nothing)
jrows() = (_JROWS[] === nothing && (_JROWS[] = j_rows()); _JROWS[])

include(joinpath(@__DIR__, "linked_dir_analysis.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "check" : ARGS[1]
    if mode == "check"
        dir_check()
    elseif mode == "kernel"
        dr_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "jrows"
        for r in j_rows(; verbose = true); @printf("%5d %-36s k²/V_G %.4f k/k_c %.3f\n", r.j, r.id, r.idx, r.ratio); end
    else
        dr_main(mode, ARGS[2:end])
    end
end
