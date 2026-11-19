# What is the source of the association between a locus's squared deviation and its linked
# background? (follow-up of validate/linked_decay.jl; criteria validate/linked_source_criteria.md,
# committed before any comparison; running log validate/linked_source.log; theory
# validate/linked_source_theory.jl; algebra validate/symbolic/check_linked_source.jl).
#
# linked_decay.jl closed the source of T_s = κ(B_s, x_l, x_l) as K_s = κ(B_s, B, x_l, x_l) ≈ −2 y_s Y_l
# (y_s = Cov(B_s, x_l), Y_l = Σ_s y_s). Exactly (centred moments, any population):
#
#   K_s = H_s − 2 y_s Y_l,     H_s = Cov(B_s B, d_l²) = μ(B_s, B, l, l) − Cov(B_s, B) v_l
#   H_s = H2_s + H3_s,         H2_s = Σ_{k∈s} Cov(d_k², d_l²)        (two loci: variance heterogeneity)
#                              H3_s = Σ_{k∈s} Σ_{j≠k,l} Cov(d_k d_j, d_l²) (three loci)
#   y_s Y_l = Σ_{k∈s} D_kl² + Σ_{k∈s} Σ_{j≠k,l} D_kl D_jl                (pairwise + cross-pair LD products)
#   so K_s = [H2_s − 2Σ D_kl²] + [H3_s − 2 Σ_{j≠k} D_kl D_jl] = Σ κ(k,k,l,l) + Σ κ(k,j,l,l) (two- and three-locus cumulants);
#   a biallelic locus l has d_l² − v_l = s_l d_l (s_l = κ3_l/v_l), hence H_s = s_l J_s,
#   J_s = κ(B_s, B, x_l), and H2_s = s_l Σ_{k∈s} κ(k,k,l) (the spectrum-shape part).
#
# Instrumented kernel runs (reading the World only; bitwise identical to `simulate`, `check`)
# measure every generation from MOVE − 60 to MOVE + 199, per group a of loci l and partner bin b,
# all of linked_decay.jl's statistics plus H, H2, J, J2 = Σ_{k∈b} κ(k,k,l), the biallelic
# projections, products for the selection terms, and every `LS_STRIDE`-th generation (a) the same
# under the selection weights π (survival probability given the genotype), and (b) per class
# s = (bin, split-probability class) of pairs the pairwise Σ D_kl² and Σ Cov(d_k², d_l²) (Gram
# matrices).
#
# Run from the package root (vignettes environment):
#   julia --project=vignettes validate/linked_source.jl check               # bitwise identity with simulate
#   julia -t 6 --project=vignettes validate/linked_source.jl kernel SET R   # instrumented runs (resumable)
#   julia --project=vignettes validate/linked_source.jl decompose SET R     # the decomposition of K
#   julia --project=vignettes validate/linked_source.jl analyse SET R       # predictions and verdicts
# SET: D (linked_decay's development inputs and seeds), F (its confirmation inputs and seeds),
# G (new confirmation set, `g_rows`, fresh seeds). Raw records: validate/linked_source_raw/ (git-ignored).

# 2026-11-19: `ne = :falconer`, `solver = :iteration` and `QTLSim.burnin_start_vg` below keep the effective size, the solver and
# the starting variance this script was run with (the package's defaults are now `:breeders` and `:bisection`; NEWS.md), so that
# it still reproduces its log.
using LinearAlgebra
include(joinpath(@__DIR__, "linked_decay.jl"))
include(joinpath(@__DIR__, "linked_source_theory.jl"))

const LS_RAW = joinpath(@__DIR__, "linked_source_raw")
const LS_STRIDE = 4                                   # selection-weighted pass and Gram every LS_STRIDE generations
const GSEED = 90261115000                             # G seeds: GSEED + 100j + r

# ---------------------------------------------------------------------------------------
# inputs

"""
G: fresh confirmation inputs. From the 102 inputs of validate/spectrum_diagnosis_inputs.csv in no
set and not in linked_decay's F, keep those passing a screen computed from the parameters only
(no run, no outcome after the move): index k²/V_G ≥ 0.01 (k the rate of the move, V_G the
house-of-Gauss variance) and k/k_c ≤ 1.4 (k_c = `drift_critical_rate`), then at most
`GPER` per stratum (NoQTL × u), in the order of MersenneTwister(20261115).
"""
const GPER = 3
function g_candidates()
    fids = Set(r.id for r in fresh_rows())
    rows = filter(r -> r.set == "" && !(r.id in fids), read_inputs_csv())
    out = NamedTuple[]
    for r in rows
        p = readp(r.id)
        vg = house_of_gauss_vg(p; ne = :falconer); kc = drift_critical_rate(p; ne = :falconer, solver = :iteration).k
        idx = p.MoveFitnessSD^2 / vg; ratio = p.MoveFitnessSD / kc
        push!(out, (r..., idx = idx, ratio = ratio, pass = idx >= 0.01 && ratio <= 1.4))
    end
    out
end
function g_rows()
    cand = filter(r -> r.pass, g_candidates())
    rng = MersenneTwister(20261115)
    out = NamedTuple[]
    for s in sort(unique(first(split(r.id, "/")) for r in cand))
        c = sort(filter(r -> startswith(r.id, s * "/"), cand); by = r -> r.id)
        append!(out, c[randperm(rng, length(c))[1:min(GPER, length(c))]])
    end
    out
end

ls_rows(set) = set == "G" ? g_rows() : ld_rows(set)
ls_seed(set, j, r) = set == "G" ? GSEED + 100j + r : ld_seed(set, j, r)
ls_file(id, r, set) = joinpath(LS_RAW, replace(id, "/" => "__", ".txt" => "") * "@" * set * "#$r.jls")

# ---------------------------------------------------------------------------------------
# measurement

mutable struct LSBuf
    L::Int; NB::Int
    binl::Vector{Int}
    xbar::Vector{Float64}; v::Vector{Float64}
    S::Array{Float64,3}              # (8, NB, L) accumulators
    acc::Matrix{Float64}             # (8, NB) accumulators of one locus
    Sd2::Vector{Float64}; Sd3::Vector{Float64}; Sd4::Vector{Float64}; SY::Vector{Float64}
    Bs::Matrix{Float64}; Gc::Vector{Float64}
    E::Matrix{Float64}; Qb::Matrix{Float64}; dG::Vector{Float64}; wv::Vector{Float64}   # (NB, M), (NB, M), (M), (M)
    x::Vector{Float64}
    XT::Matrix{Float64}; X::Matrix{Float64}   # deviations (L, M) and (M, L)
    Surv::Vector{Float64}            # QTSim's survival function of the lowest crossover (pair splits)
    rcl::Matrix{UInt8}               # split-probability class of each pair (0 on the diagonal)
    X32::Matrix{Float32}; Q32::Matrix{Float32}; Gt::Matrix{Float32}
end

function LSBuf(L, Nmax, u; NB = LD_NB)
    Sv = qtsim_survival(L, u)
    rcl = zeros(UInt8, L, L)
    for i in 1:L, j in 1:L
        i == j && continue
        rcl[i, j] = ld_rclass(abs(Sv[i] - Sv[j]))
    end
    nmax = maximum(count(==(a), [ld_bin(l - 1, L; nb = NB) for l in 1:L]) for a in 1:NB)
    M = 2Nmax
    LSBuf(L, NB, [ld_bin(l - 1, L; nb = NB) for l in 1:L], zeros(L), zeros(L), zeros(8, NB, L), zeros(8, NB),
          zeros(L), zeros(L), zeros(L), zeros(L), zeros(NB, M), zeros(M),
          zeros(NB, M), zeros(NB, M), zeros(M), zeros(M), zeros(L), zeros(L, M), zeros(M, L),
          Sv, rcl, zeros(Float32, M, L), zeros(Float32, M, L), zeros(Float32, L, nmax))
end

"""
Association statistics of the adults' 2N chromosomes, each weighted by `wt[i]` of its carrier
(all ones: the adults; π: the π-weighted adults, i.e. after deterministic selection).
Group × bin matrices are sums over l ∈ a of the per-(l, b) statistics. Loci outer, chromosomes
inner (the deviations are stored once), so the accumulators of one locus stay in cache.
"""
function ls_measure!(B::LSBuf, w, wt::AbstractVector; gram = false)
    L, NB, N = B.L, B.NB, w.N
    M = 2N
    x = B.x
    W = 0.0
    fill!(B.xbar, 0)
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        ldchr!(x, w, c, i)
        wi = wt[i]; W += wi; B.wv[m] = wi
        g = 0.0
        @inbounds for b in 1:NB
            B.Bs[b, m] = 0.0
        end
        @inbounds for l in 1:L
            B.Bs[B.binl[l], m] += x[l]; g += x[l]; B.xbar[l] += wi * x[l]
        end
        B.Gc[m] = g
        @inbounds for l in 1:L
            B.XT[l, m] = x[l]
        end
    end
    B.xbar ./= W
    Bbar = zeros(NB); Gbar = 0.0
    for m in 1:M
        Bbar .+= B.wv[m] .* view(B.Bs, :, m); Gbar += B.wv[m] * B.Gc[m]
    end
    Bbar ./= W; Gbar /= W
    fill!(B.v, 0)
    @inbounds for m in 1:M
        wm = B.wv[m]
        for l in 1:L
            d = B.XT[l, m] - B.xbar[l]
            B.XT[l, m] = d; B.v[l] += wm * d * d
        end
    end
    B.v ./= W
    v = B.v
    X = view(B.X, 1:M, :)
    permutedims!(X, view(B.XT, :, 1:M), (2, 1))
    E = B.E; Qb = B.Qb
    @inbounds for m in 1:M
        B.dG[m] = B.Gc[m] - Gbar
        for b in 1:NB
            E[b, m] = B.Bs[b, m] - Bbar[b]; Qb[b, m] = 0.0
        end
        for l in 1:L
            Qb[B.binl[l], m] += B.XT[l, m]^2 - v[l]
        end
    end
    S = B.S; acc = B.acc
    @inbounds for l in 1:L
        fill!(acc, 0)
        vl = v[l]; ob = B.binl[l]
        s2 = 0.0; s3 = 0.0; s4 = 0.0; sY = 0.0
        o1 = 0.0; o2 = 0.0; o3 = 0.0; o4 = 0.0; o5 = 0.0; o6 = 0.0; o7 = 0.0; o8 = 0.0
        for m in 1:M
            d = X[m, l]; wi = B.wv[m]
            d2 = d * d; d3 = d2 * d; d2c = d2 - vl
            dB = B.dG[m] - d
            wd = wi * d; wd2 = wi * d2; wdBd2 = wi * dB * d2; wd3 = wi * d3; wdB = wi * dB; wd2c = wi * d2c; wdBd = wi * dB * d
            @simd for b in 1:NB
                eb = E[b, m]; qb = Qb[b, m]
                acc[1, b] += eb * wd
                acc[2, b] += eb * wd2
                acc[3, b] += eb * wdBd2
                acc[4, b] += eb * wd3
                acc[5, b] += eb * wdB
                acc[6, b] += qb * wd2c
                acc[7, b] += qb * wd
                acc[8, b] += eb * wdBd
            end
            # own bin: e_ob − d and q_ob − d2c
            o1 -= d * wd; o2 -= d * wd2; o3 -= d * wdBd2; o4 -= d * wd3; o5 -= d * wdB
            o6 -= d2c * wd2c; o7 -= d2c * wd; o8 -= d * wdBd
            s2 += wd2; s3 += wd3; s4 += wd2 * d2; sY += wdBd
        end
        acc[1, ob] += o1; acc[2, ob] += o2; acc[3, ob] += o3; acc[4, ob] += o4
        acc[5, ob] += o5; acc[6, ob] += o6; acc[7, ob] += o7; acc[8, ob] += o8
        for b in 1:NB, k in 1:8
            S[k, b, l] = acc[k, b]
        end
        B.Sd2[l] = s2; B.Sd3[l] = s3; B.Sd4[l] = s4; B.SY[l] = sY
    end
    ks = (:y, :T, :K, :K3, :yY, :H, :H2, :J, :J2, :sJ, :sJ2, :sy, :TY, :yC, :JY, :yJ, :Hb, :CCov, :Jv, :Tv2Y)
    A = Dict(k => zeros(NB, NB) for k in ks)
    V = zeros(NB); V2 = zeros(NB); k3 = zeros(NB); g = zeros(NB)
    shape_res = 0.0; shape_tot = 0.0
    sumC = 0.0; sk3 = 0.0; sk3sq = 0.0; sk3C = 0.0; genic = 0.0
    Jl = zeros(NB)
    @inbounds for l in 1:L
        a = B.binl[l]
        vl = B.Sd2[l] / W; κ3 = B.Sd3[l] / W; μ4 = B.Sd4[l] / W; Y = B.SY[l] / W
        sl = vl > 0 ? κ3 / vl : 0.0
        V[a] += vl; V2[a] += vl^2; k3[a] += κ3; g[a] += vl
        if vl > 0
            shape_tot += μ4 - vl^2; shape_res += μ4 - vl^2 - κ3^2 / vl
        end
        Cl = 0.0; Jtot = 0.0
        for b in 1:NB
            Cl += S[2, b, l] / W; Jtot += S[8, b, l] / W
        end
        for b in 1:NB
            y = S[1, b, l] / W; T = S[2, b, l] / W
            H = S[3, b, l] / W - (S[5, b, l] / W) * vl
            K = H - 2y * Y
            K3 = S[4, b, l] / W - 3y * vl
            J = S[8, b, l] / W; J2 = S[7, b, l] / W
            A[:y][a, b] += y; A[:T][a, b] += T; A[:K][a, b] += K; A[:K3][a, b] += K3; A[:yY][a, b] += y * Y
            A[:H][a, b] += H; A[:H2][a, b] += S[6, b, l] / W; A[:J][a, b] += J; A[:J2][a, b] += J2
            A[:sJ][a, b] += sl * J; A[:sJ2][a, b] += sl * J2; A[:sy][a, b] += sl * y
            A[:TY][a, b] += T * Y; A[:yC][a, b] += y * Cl; A[:JY][a, b] += J * Y; A[:yJ][a, b] += y * Jtot
            A[:Hb][a, b] += (S[5, b, l] / W) * vl          # Cov(B_b, B) v_l (for the selection terms)
            A[:CCov][a, b] += Cl * (S[5, b, l] / W)          # C_l Cov(B_b, B)
            A[:Jv][a, b] += J * vl
            A[:Tv2Y][a, b] += T * (vl + 2Y)               # with VGchr ΣT: Σ T_lb Var(B_l)
        end
        sumC += Cl; sk3 += κ3; sk3sq += κ3^2; sk3C += κ3 * Cl; genic += 2vl
    end
    Cbb = zeros(NB, NB); cb = zeros(NB)
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        dG = B.Gc[m] - Gbar
        for b in 1:NB
            eb = B.Bs[b, m] - Bbar[b]
            cb[b] += wt[i] * eb * dG
            for b2 in 1:NB
                Cbb[b, b2] += wt[i] * eb * (B.Bs[b2, m] - Bbar[b2])
            end
        end
    end
    Cbb ./= W; cb ./= W
    VGchr = 0.0
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        VGchr += wt[i] * (B.Gc[m] - Gbar)^2
    end
    VGchr /= W
    ns = NB * LD_NR
    Dsq = fill(NaN, NB, ns); H2p = fill(NaN, NB, ns)
    if gram
        fill!(Dsq, 0); fill!(H2p, 0)
        Xv = view(B.X32, 1:M, :); Qv = view(B.Q32, 1:M, :)
        @inbounds for l in 1:L, m in 1:M
            d = B.X[m, l]
            Xv[m, l] = d; Qv[m, l] = d * d - v[l]
        end
        for a in 1:NB
            cols = let f = findall(==(a), B.binl); isempty(f) ? (1:0) : (first(f):last(f)) end
            Gt = view(B.Gt, :, 1:length(cols))
            # D_kl for l ∈ a, all k
            mul!(Gt, transpose(Xv), view(Xv, :, cols))
            Gt ./= M
            for (jj, l) in enumerate(cols), k in 1:L
                k == l && continue
                s = (B.rcl[k, l] - 1) * NB + B.binl[k]
                D = Float64(Gt[k, jj])
                Dsq[a, s] += D^2
            end
            # Cov(d_k², d_l²)
            mul!(Gt, transpose(Qv), view(Qv, :, cols))
            Gt ./= M
            for (jj, l) in enumerate(cols), k in 1:L
                k == l && continue
                s = (B.rcl[k, l] - 1) * NB + B.binl[k]
                H2p[a, s] += Float64(Gt[k, jj])
            end
        end
    end
    (; A..., V, V2, k3, g, c = cb, Cbb, sumC, sk3, sk3sq, sk3C, genic, VGchr, shape_res, shape_tot, Dsq, H2p, W)
end

# ---------------------------------------------------------------------------------------
# the instrumented run (linked_decay's ld_run, extended; draws nothing)

const LS_SCAL = [:gen, :N, :nsurv, :theta, :gbar, :VG, :VGchr, :Sc, :beta, :beta2, :genic, :sumC, :sk3, :sk3sq, :sk3C,
                 :VE, :omega2, :Vm, :shape_res, :shape_tot, :VGchr_pi]
const LS_MATS = [:y, :T, :K, :K3, :yY, :H, :H2, :J, :J2, :sJ, :sJ2, :sy, :TY, :yC, :JY, :yJ, :Hb, :CCov, :Jv, :Tv2Y, :Cbb]
const LS_VECS = [:V, :V2, :k3, :g, :c]
const LS_PI = [:y, :T, :K, :yY, :H, :H2, :J]          # under π weights (every LS_STRIDE)

function ls_run(p, seed; from = LD_FROM, to = LD_TO, stride = LS_STRIDE, decompose = true)
    BLAS.set_num_threads(1)
    q = kernel_params(p)
    q.Selection === :viability || error("viability only")
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    L = q.NoQTL
    Bf = LSBuf(L, q.PopSize, q.RecombPoissonMean)
    VE = QTLSim.envVar(q)
    rec = Dict{Symbol,Any}(s => Float64[] for s in LS_SCAL)
    for s in [LS_MATS; Symbol.(:pi_, LS_PI); :Dsq; :H2p]
        rec[s] = Matrix{Float32}[]
    end
    for s in LS_VECS
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
            dostride = decompose && (it - from) % stride == 0
            dogram = dostride && (it - from) % (L > 2000 ? 2stride : stride) == 0     # Gram every 8th generation above 2000 loci
            ms = ls_measure!(Bf, w, ones(N); gram = dogram)
            Sc = ω2 + VE + ms.VGchr
            mp = dostride ? ls_measure!(Bf, w, π ./ πbar) : nothing
            vals = (gen = it, N = N, nsurv = w.nsurv, theta = w.θ, gbar = gb, VG = var(g), VGchr = ms.VGchr, Sc = Sc,
                    beta = β, beta2 = β2, genic = ms.genic, sumC = ms.sumC, sk3 = ms.sk3, sk3sq = ms.sk3sq, sk3C = ms.sk3C,
                    VE = VE, omega2 = ω2, Vm = QTLSim.muteVar(q), shape_res = ms.shape_res, shape_tot = ms.shape_tot,
                    VGchr_pi = mp === nothing ? NaN : mp.VGchr)
            for s in LS_SCAL
                push!(rec[s], Float64(vals[s]))
            end
            for s in LS_MATS
                push!(rec[s], Float32.(getfield(ms, s)))
            end
            for s in LS_PI
                push!(rec[Symbol(:pi_, s)], mp === nothing ? fill(NaN32, 1, 1) : Float32.(getfield(mp, s)))
            end
            push!(rec[:Dsq], Float32.(ms.Dsq)); push!(rec[:H2p], Float32.(ms.H2p))
            for s in LS_VECS
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

function ls_check(; id = "QT1000-Possion3.0/Input-09697.txt", j = 89, r = 1)
    p = readp(id)
    t = @elapsed tr = ls_run(p, KSEED + 100j + r; to = LD_FROM + 9)
    o = simulate(kernel_params(p); seed = KSEED + 100j + r, init = :equilibrium, burnin = KBURN)
    k = [findfirst(==(g), o.gen) for g in LD_FROM:LD_FROM+length(tr.vgser)-1]
    same = all(o.geno_var[k] .== tr.vgser)
    # the new statistics reproduce linked_decay's on the same generations
    t2 = ld_run(p, KSEED + 100j + r; to = LD_FROM + 9)
    agree = true
    for (s, s2) in ((:K, :Kab), (:T, :Tab), (:yY, :yYab)), kk in 1:length(tr.vgser)
        a1 = Float64.(tr.rec[s][kk]); a2 = Float64.(t2.rec[s2][kk])
        agree &= maximum(abs.(a1 .- a2)) <= 1e-4 * max(maximum(abs.(a2)), 1e-30)
    end
    @printf("check %s: %d generations, V_G bitwise identical to simulate: %s; K, T, yY as linked_decay's: %s (%.1f s)\n",
            id, length(tr.vgser), same, agree, t)
    (same && agree) || error("check failed")
end

function ls_kernel(set, R; decompose = set != "G")
    mkpath(LS_RAW)
    rows = ls_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(ls_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) instrumented runs for set $set, R = $R, decompose = $decompose, $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time()
    cpu = Threads.Atomic{Float64}(0.0)
    Threads.@threads :greedy for (id, j, k) in jobs
        t1 = time()
        tr = ls_run(readp(id), ls_seed(set, j, k); decompose)
        serialize(ls_file(id, k, set), tr)
        Threads.atomic_add!(cpu, time() - t1)
        d = Threads.atomic_add!(done, 1) + 1
        @printf("  %d/%d %s #%d  %.0f s\n", d, length(jobs), id, k, time() - t0)
    end
    @printf("done: %d runs, %.0f s wall, %.0f thread-s\n", length(jobs), time() - t0, cpu[])
end


# ---------------------------------------------------------------------------------------
# the decomposition of K (development; reported, not a verdict)

function load_ls(id, R, set)
    out = []
    for k in 1:R
        f = ls_file(id, k, set)
        isfile(f) && push!(out, deserialize(f))
    end
    out
end

tot(x) = sum(Float64, x)
const DWINS = [(-60, -1), (0, 24), (25, 49), (50, 99), (100, 199)]

# per run: per generation series of the parts (sums over all l and b); NaN where not measured
function parts(tr)
    r = tr.rec
    t = Int.(r[:gen]) .- MOVE
    n = length(t)
    S(s) = [tot(r[s][k]) for k in 1:n]
    K = S(:K); yY = S(:yY); H = S(:H); H2 = S(:H2); J = S(:J)
    Dsq = [all(isnan, r[:Dsq][k]) ? NaN : tot(r[:Dsq][k]) for k in 1:n]
    H2p = [all(isnan, r[:H2p][k]) ? NaN : tot(r[:H2p][k]) for k in 1:n]
    piS(s) = [length(r[Symbol(:pi_, s)][k]) == 1 ? NaN : tot(r[Symbol(:pi_, s)][k]) for k in 1:n]
    β = r[:beta]
    selcand = [2β[k] * (tot(r[:TY][k]) + tot(r[:JY][k]) + tot(r[:yC][k]) + tot(r[:yJ][k])) for k in 1:n]
    # moment form, dropping the fifth-order central-moment heterogeneity term:
    # −β[C(y_s + Cov(B_s,B)) + T_s(Y + Var B) + 2J_s(v + Y)]
    selmom = haskey(r, :CCov) ? [-β[k] * (tot(r[:yC][k]) + tot(r[:CCov][k]) + tot(r[:TY][k]) + r[:VGchr][k] * tot(r[:T][k]) - tot(r[:Tv2Y][k]) +
                                        2tot(r[:Jv][k]) + 2tot(r[:JY][k])) for k in 1:n] : fill(NaN, n)
    msk = isnan.(piS(:H)); selcand[msk] .= NaN; selmom[msk] .= NaN      # compare on the generations with the π pass
    gm(x) = ifelse.(isnan.(Dsq), NaN, x)
    (; t, K, m2yY = -2 .* yY, yY, H, H2, yYg = gm(yY), Kg = gm(K), H2g = gm(H2), H3 = H .- H2, Dsq, H2p, cross = yY .- Dsq,
     K2 = H2p .- 2 .* Dsq, K3c = K .- (H2p .- 2 .* Dsq), sJ = S(:sJ), sJ2 = S(:sJ2), T = S(:T), sy = S(:sy),
     shape = r[:shape_res] ./ r[:shape_tot], dHsel = piS(:H) .- H, dKsel = piS(:K) .- K, dyYsel = piS(:yY) .- yY,
     dTsel = piS(:T) .- S(:T), selcand, selmom, beta = β, C = r[:sumC], genic = r[:genic], nsurv = r[:nsurv])
end

function ls_decompose(set, R; io = stdout)
    rows = ls_rows(set)
    cols = [:K, :m2yY, :H, :H2, :H3, :Dsq, :H2p, :yYg, :Kg, :H2g, :cross, :K2, :K3c, :sJ, :sJ2, :T, :sy, :dHsel, :selcand, :selmom, :dKsel, :dyYsel, :dTsel]
    println(io, "\n== decomposition of K, set $set, R = $R (sums over loci; window means over runs alive) ==")
    summ = Dict{Symbol,Dict{Tuple{Int,Int},Vector{Float64}}}()
    for r in rows
        runs = load_ls(r.id, R, set)
        isempty(runs) && continue
        P = [parts(tr) for tr in runs]
        @printf(io, "\n%s  L=%d u=%g runs %d\n  %-10s", r.id, runs[1].L, runs[1].u, length(runs), "window")
        for c in cols
            @printf(io, " %10s", c)
        end
        println(io)
        for wdw in DWINS
            @printf(io, "  %4d..%-4d", wdw...)
            vals = Dict(c => winmean([(p.t, getfield(p, c)) for p in P], wdw) for c in cols)
            for c in cols
                @printf(io, " %10.3e", vals[c])
            end
            println(io)
            for (nm, v) in ((:K_over_m2yY, vals[:K] / vals[:m2yY]), (:H_over_2yY, -vals[:H] / vals[:m2yY]),
                            (:H2_over_H, vals[:H2] / vals[:H]), (:Dsq_over_yY, vals[:Dsq] / vals[:yYg]), (:H2p_over_H2, vals[:H2p] / vals[:H2g]),
                            (:K2_over_K, vals[:K2] / vals[:Kg]), (:sJ_over_H, vals[:sJ] / vals[:H]), (:sJ2_over_H2, vals[:sJ2] / vals[:H2]),
                            (:sy_over_T, vals[:sy] / vals[:T]), (:selcand_over_dHsel, vals[:selcand] / vals[:dHsel]), (:selmom_over_dHsel, vals[:selmom] / vals[:dHsel]))
                push!(get!(get!(summ, nm, Dict{Tuple{Int,Int},Vector{Float64}}()), wdw, Float64[]), v)
            end
        end
        @printf(io, "  shape (residual within-locus variance of d² not explained by a biallelic locus, fraction): %.3f\n",
                mean(filter(isfinite, reduce(vcat, [p.shape for p in P]))))
    end
    println(io, "\n== medians over inputs ==")
    for (nm, d) in sort(collect(summ); by = first)
        @printf(io, "  %-20s", nm)
        for wdw in DWINS
            @printf(io, "  %4d..%-4d %7.3f", wdw..., medf(d[wdw]))
        end
        println(io)
    end
    summ
end


# ---------------------------------------------------------------------------------------
# predictions (criteria: validate/linked_source_criteria.md)

# linked_decay's record names, so that its levelA, levelB, drivers and measured apply unchanged
function ldview(tr)
    r = copy(tr.rec)
    for (a, b) in ((:Tab, :T), (:Kab, :K), (:K3ab, :K3), (:yab, :y), (:yYab, :yY))
        r[a] = tr.rec[b]
    end
    (; tr.seed, tr.extinct_gen, rec = r, tr.vgser, tr.L, tr.u, tr.PopSize)
end

"""
Level B′ (closed): from t = −60, the y and H moments at their joint stationary state for the drivers
averaged over −60..−41, T = 0; each generation the source K_s = Σ_t(H_st − 2P_st) feeds the T
recursion (linked_decay's derived decay), then `yH_step!`. Driven by the run's measured β, S_c, n,
N and variances. Returns per generation the predicted Σ K, Σ C, Σ E[yY] and Σ H.
"""
function levelBp(tr, G::LDGroups; couple = true, hsel = :none)
    rec = tr.rec
    n = length(rec[:gen])
    st = LSState(G)
    yH_stationary!(st, G, drivers(rec, 1; avg = PRE_AVG); couple, hsel = :none)
    NA = size(G.n, 1)
    T = [zeros(size(G.n, 2)) for _ in 1:NA]
    out = (K = fill(NaN, n), C = fill(NaN, n), yY = fill(NaN, n), H = fill(NaN, n))
    for k in 1:n
        d = drivers(rec, k)
        Ks = [G.nb[a] == 0 ? zeros(size(G.n, 2)) : source_K(st, a) for a in 1:NA]
        out.K[k] = sum(sum, Ks); out.C[k] = sum(sum, T)
        out.yY[k] = sum(sum(st.y.P[a]) for a in 1:NA); out.H[k] = sum(sum(st.H[a]) for a in 1:NA)
        for a in 1:NA
            G.nb[a] == 0 && continue
            T_step!(T[a], G, d, a, Ks[a]; mode = :derived)
            yH_step!(st, G, d, a; couple, hsel, Tvec = T[a])
        end
    end
    out
end

const BP_VARIANTS = [(:Bp, (couple = true, hsel = :bulmer)), (:Bd, (couple = true, hsel = :none))]

function ls_analyse_input(r, R, set; io = stdout)
    runs = load_ls(r.id, R, set)
    isempty(runs) && return nothing
    G = ld_groups(runs[1].L, runs[1].u)
    ser = Dict{Symbol,Vector{Any}}()
    cum = Dict{Symbol,Vector{Float64}}()
    for tr0 in runs
        tr = ldview(tr0)
        m = measured(tr)
        k0 = findfirst(==(-1), m.t)
        g0 = k0 === nothing ? NaN : mean(m.genic[max(1, k0 - 19):k0])
        for (s, v) in ((:C, m.C), (:Csym, m.Csym), (:K, m.K), (:yY, m.yY), (:H, m.K .+ 2 .* m.yY))
            push!(get!(ser, s, []), (m.t, v))
        end
        for md in AMODES
            push!(get!(ser, Symbol(:A_, md), []), (m.t, levelA(tr, G; mode = md)))
        end
        B = levelB(tr, G)
        push!(get!(ser, :B_C, []), (m.t, B.predC)); push!(get!(ser, :B_K, []), (m.t, -2 .* B.predyY))
        preds = Dict{Symbol,Any}(:B => B.predC)
        for (nm, kw) in BP_VARIANTS
            P = levelBp(tr, G; kw...)
            push!(get!(ser, Symbol(nm, :_C), []), (m.t, P.C)); push!(get!(ser, Symbol(nm, :_K), []), (m.t, P.K))
            push!(get!(ser, Symbol(nm, :_yY), []), (m.t, P.yY)); push!(get!(ser, Symbol(nm, :_H), []), (m.t, P.H))
            preds[nm] = P.C
        end
        k100 = findfirst(==(100), m.t)
        if k0 !== nothing && k100 !== nothing
            rg = k0+1:k100-1
            push!(get!(cum, :meas_loss, Float64[]), sum(2 .* m.beta[rg] .* m.C[rg]) / g0)
            push!(get!(cum, :meas_losssym, Float64[]), sum(2 .* m.beta[rg] .* m.Csym[rg]) / g0)
            for (nm, pc) in preds
                push!(get!(cum, Symbol(nm, :_loss), Float64[]), sum(2 .* m.beta[rg] .* pc[rg]) / g0)
            end
        end
    end
    res = Dict{Symbol,Any}()
    for (s, v) in ser
        res[s] = Dict(wdw => winmean(v, wdw) for wdw in WINS)
    end
    for (s, v) in cum
        res[s] = isempty(v) ? NaN : mean(v)
    end
    res[:nalive100] = length(get(cum, :meas_loss, Float64[]))
    haskey(res, :meas_loss) || (res[:meas_loss] = NaN; res[:meas_losssym] = NaN)
    @printf(io, "\n%s  L=%d u=%g runs %d (alive at +100: %d)\n", r.id, runs[1].L, runs[1].u, length(runs), res[:nalive100])
    cols = [:K, :B_K, :Bp_K, :H, :Bp_H, :yY, :Bp_yY, :C, :A_derived, :B_C, :Bp_C]
    @printf(io, "  %-12s", "window")
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
    @printf(io, "  cumulative loss 0..99 / genic(pre): measured 2βC %+.3f (sign-symmetric %+.3f), level B %+.3f, level B′ %+.3f\n",
            res[:meas_loss], res[:meas_losssym], get(res, :B_loss, NaN), get(res, :Bp_loss, NaN))
    res
end

function ls_analyse(set, R; io = stdout)
    rows = ls_rows(set)
    println(io, "\n== linked source, set $set, R = $R: per input, window means over runs alive (sums over all loci) ==")
    bufs = [IOBuffer() for _ in rows]
    resv = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        resv[i] = ls_analyse_input(rows[i], R, set; io = bufs[i])
    end
    foreach(b -> print(io, String(take!(b))), bufs)
    all = [(rows[i], resv[i]) for i in eachindex(rows)]
    all = filter(x -> x[2] !== nothing, all)
    elig = filter(x -> isfinite(x[2][:meas_loss]) && x[2][:meas_loss] <= -0.05, all)
    println(io, "\n== summary: $(length(all)) inputs, $(length(elig)) eligible (measured cumulative 2βΣC over 0..99 ≤ −0.05 of the genic variance) ==")
    println(io, "(a) source: predicted ΣK / measured ΣK, median over all inputs (level B′; level B = linked_decay's −2E[yY] for reference):")
    for wdw in WINS
        @printf(io, "  %4d..%-4d  B′ %.3f   B %.3f   B′ drift-only %.3f\n", wdw..., medf([x[2][:Bp_K][wdw] / x[2][:K][wdw] for x in all]),
                medf([x[2][:B_K][wdw] / x[2][:K][wdw] for x in all]), medf([x[2][:Bd_K][wdw] / x[2][:K][wdw] for x in all]))
    end
    println(io, "    level B′ Σ E[yY] / measured Σ yY and Σ H / measured Σ H (median):")
    for wdw in WINS
        @printf(io, "  %4d..%-4d  yY %.3f   H %.3f\n", wdw..., medf([x[2][:Bp_yY][wdw] / x[2][:yY][wdw] for x in all]),
                medf([x[2][:Bp_H][wdw] / x[2][:H][wdw] for x in all]))
    end
    println(io, "(b) closed ΣC, eligible inputs: log(predicted / measured), median and median |.| (level B′, level B):")
    for wdw in WINS[3:5]
        v = [lr(x[2][:Bp_C][wdw], x[2][:C][wdw]) for x in elig]; u = [lr(x[2][:B_C][wdw], x[2][:C][wdw]) for x in elig]
        @printf(io, "  ΣC %4d..%-4d  B′ median %+.3f |.| %.3f (n %d)   B median %+.3f |.| %.3f\n", wdw..., medf(v), mabs(v), count(isfinite, v), medf(u), mabs(u))
    end
    v = [lr(x[2][:Bp_loss], x[2][:meas_losssym]) for x in elig]; u = [lr(x[2][:B_loss], x[2][:meas_losssym]) for x in elig]
    @printf(io, "  cumulative loss 0..99 vs sign-symmetric: B′ median %+.3f |.| %.3f (n %d)   B median %+.3f |.| %.3f\n",
            medf(v), mabs(v), count(isfinite, v), medf(u), mabs(u))
    for wdw in WINS[3:5]
        v = [lr(x[2][:Bd_C][wdw], x[2][:C][wdw]) for x in elig]
        @printf(io, "  comparator B′ drift-only ΣC %4d..%-4d median %+.3f |.| %.3f\n", wdw..., medf(v), mabs(v))
    end
    v = [lr(x[2][:Bd_loss], x[2][:meas_losssym]) for x in elig]
    @printf(io, "  comparator B′ drift-only cumulative loss: median %+.3f |.| %.3f\n", medf(v), mabs(v))
    println(io, "P1 (decay rate, level A, eligible): log(predicted ΣC / measured ΣC), median and median |.|:")
    for md in AMODES, wdw in WINS[3:5]
        v = [lr(x[2][Symbol(:A_, md)][wdw], x[2][:C][wdw]) for x in elig]
        @printf(io, "  %-9s %4d..%-4d median %+.3f  |.| %.3f  (n %d)\n", md, wdw..., medf(v), mabs(v), count(isfinite, v))
    end
    all
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "check" : ARGS[1]
    if mode == "check"
        ls_check()
    elseif mode == "kernel"
        ls_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "decompose"
        ls_decompose(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "analyse"
        ls_analyse(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "candidates"
        for r in g_candidates()
            @printf("%-36s idx %.4f k/kc %.3f %s\n", r.id, r.idx, r.ratio, r.pass ? "pass" : "")
        end
        println("G: ", [r.id for r in g_rows()])
    end
end
