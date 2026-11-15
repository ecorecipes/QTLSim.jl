# Why does the allele–background LD recursion of linked_decay/linked_source under-predict
# E[y_s y_t]? (follow-up of validate/linked_source.jl; criteria validate/linked_ld_criteria.md,
# committed before any comparison on fresh data; running log validate/linked_ld.log; theory
# validate/linked_ld_theory.jl; algebra validate/symbolic/check_linked_ld.jl).
#
# Notation (chromosome level; locus l, background B_l = G_chr − x_l, d_l = x_l − x̄_l):
#   Y_l = Cov(B_l, x_l) = Σ_s y_{l,s},   the statistic under study is Σ_l Y_l² (level B's Σ E[yY]).
#
# Stage accounting (`stage` runs; reads the World only, `check` verifies bitwise identity with
# `simulate`; the survival Monte Carlo uses its own RNG). Per generation, summed over loci:
#   A   adults                       Σ Y², Σ μ22 (= Σ E[(B−B̄)² d²]), Σ Var(B) v, Σ Y, Σ v Y
#   Π   adults weighted by π (expected survival given the genotype): deterministic selection
#   S   the survivors; and E[Σ Ŷ² | adults] over `NMC` redraws of Bernoulli(π) survival
#   Gp  the survivors' expected gamete pool (QTSim's single crossover): exactly
#       Y_l^Gp = Y_l^S − Cov_S(R_l, x_l) + Z_l,  R_l = Σ_k ρ_kl x_k (prefix sums: ρ_kl = |S_k − S_l|),
#       Z_l = (1/2n) Σ_i Σ_c d_l^{i,c} (R_l^{i,c'} − R̄_l)  (the recombinant's cross-chromosome part)
#   O   the offspring before mutation; next generation's A is after mutation.
# Each measured stage change is set against (i) the exact first-order expectation given the
# measured state at the start of the stage and (ii) the recursion's formula for that stage
# (linked_decay's y recursion / linked_source's B′) evaluated at the same state.
#
# Run from the package root (vignettes environment):
#   julia --project=vignettes validate/linked_ld.jl check                 # bitwise identity with simulate
#   julia -t 6 --project=vignettes validate/linked_ld.jl stage SET R      # stage-instrumented runs (resumable)
#   julia --project=vignettes validate/linked_ld.jl stages SET R          # the per-stage accounting
# Raw records: validate/linked_ld_raw/ (git-ignored).

using LinearAlgebra
include(joinpath(@__DIR__, "linked_source.jl"))

const LL_RAW = joinpath(@__DIR__, "linked_ld_raw")
const ST_FROM, ST_TO = MOVE - 60, MOVE + 100        # stage-instrumented generations [ST_FROM, ST_TO)
const NMC = 16                                      # survival redraws per generation

# ---------------------------------------------------------------------------------------
# measurement

mutable struct STBuf
    L::Int
    X::Matrix{Float64}          # (L, 2Nmax) raw allele values of the current pool
    G::Vector{Float64}          # chromosome totals
    S::Vector{Float64}          # QTSim survival function: ρ_kl = |S_k − S_l|
    xbar::Vector{Float64}; v::Vector{Float64}; cGx::Vector{Float64}; mu22::Vector{Float64}
    r1::Vector{Float64}; r2::Vector{Float64}; rbar::Vector{Float64}
    a1::Vector{Float64}; a2::Vector{Float64}; a3::Vector{Float64}
end
STBuf(L, Nmax, u) = STBuf(L, zeros(L, 2Nmax), zeros(2Nmax), qtsim_survival(L, u),
                          zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L))

function load_pool!(B::STBuf, pool, chr, N)
    m = 0
    for i in 1:N, c in 1:2
        m += 1
        loadchr!(view(B.X, :, m), pool, chr, c, i)
        B.G[m] = sum(view(B.X, :, m))
    end
    B
end

"""
Per-locus moments of the chromosomes of `inds` (both chromosomes of each), weighted by `wt[i]`:
returns (Y, v, VarB, c = Cov(B_l, G), mu22, xbar, Gbar, VG, W) as vectors/scalars (copies).
"""
function moments(B::STBuf, inds, wt)
    L = B.L
    xbar = zeros(L); Gbar = 0.0; W = 0.0
    @inbounds for i in inds, c in 1:2
        m = 2(i - 1) + c; wi = wt[i]
        W += wi; Gbar += wi * B.G[m]
        @simd for l in 1:L
            xbar[l] += wi * B.X[l, m]
        end
    end
    xbar ./= W; Gbar /= W
    v = zeros(L); cGx = zeros(L); mu = zeros(L); m3 = zeros(L); m4 = zeros(L); VG = 0.0
    @inbounds for i in inds, c in 1:2
        m = 2(i - 1) + c; wi = wt[i]
        dG = B.G[m] - Gbar
        VG += wi * dG^2
        @simd for l in 1:L
            d = B.X[l, m] - xbar[l]
            eB = dG - d
            d2 = d * d
            v[l] += wi * d2; cGx[l] += wi * d * dG; mu[l] += wi * eB * eB * d2
            m3[l] += wi * eB * d * dG; m4[l] += wi * eB * d * dG * dG
        end
    end
    v ./= W; cGx ./= W; mu ./= W; m3 ./= W; m4 ./= W; VG /= W
    Y = cGx .- v
    VarB = VG .- 2 .* cGx .+ v
    cB = VG .- cGx
    k4 = m4 .- Y .* VG .- 2 .* cB .* cGx     # κ(B, x, G, G)
    (; Y, v, VarB, c = cB, mu22 = mu, k3 = m3, k4, cx = cGx, xbar, Gbar, VG, W)
end

# R_l = Σ_k ρ_kl x_k for one chromosome x (ρ_kl = |S_k − S_l|): 2P1_l − 2S_l P0_l + S_l T0 − T1
function rsum!(R, x, S)
    T0 = sum(x); T1 = dot(S, x)
    P0 = 0.0; P1 = 0.0
    @inbounds for l in eachindex(x)
        R[l] = 2P1 - 2S[l] * P0 + S[l] * T0 - T1
        P0 += x[l]; P1 += S[l] * x[l]
    end
    R
end

"""
The survivors' expected gamete pool. `ms` are the survivors' moments. Returns per locus
CovRx = Cov_S(R_l, x_l), VarR = Var_S(R_l), Z (the cross-chromosome part), Ygp = Y − CovRx + Z.
"""
function gamete_pool(B::STBuf, surv, ms)
    L = B.L; n = length(surv)
    rbar = rsum!(B.rbar, ms.xbar, B.S)
    CovRx = zeros(L); VarR = zeros(L); Z = zeros(L)
    r1, r2 = B.r1, B.r2
    @inbounds for i in surv
        m1 = 2(i - 1) + 1; m2 = m1 + 1
        rsum!(r1, view(B.X, :, m1), B.S); rsum!(r2, view(B.X, :, m2), B.S)
        @simd for l in 1:L
            d1 = B.X[l, m1] - ms.xbar[l]; d2 = B.X[l, m2] - ms.xbar[l]
            e1 = r1[l] - rbar[l]; e2 = r2[l] - rbar[l]
            CovRx[l] += d1 * e1 + d2 * e2
            VarR[l] += e1 * e1 + e2 * e2
            Z[l] += d1 * e2 + d2 * e1
        end
    end
    CovRx ./= 2n; VarR ./= 2n; Z ./= 2n
    (; CovRx, VarR, Z, Ygp = ms.Y .- CovRx .+ Z)
end

"""
The recursion's class-averaged split probabilities as piecewise-constant weights: for locus l
(group a), partner k gets ρ̄[a, s(k, l)] (s = partner bin × split class, `ld_groups`); stored
as contiguous segments of k (ρ_kl is monotone on each side of l, bins are contiguous).
"""
struct Segs
    off::Vector{Int}; lo::Vector{Int}; hi::Vector{Int}; w::Vector{Float64}; w2::Vector{Float64}
end
function class_segments(L, u)
    G = ld_groups(L, u); S = qtsim_survival(L, u)
    off = [1]; lo = Int[]; hi = Int[]; w = Float64[]; w2 = Float64[]
    # class mean of ρ² (for the within-class spread of the recombinant noise)
    ns = size(G.n, 2); r2 = zeros(LD_NB, ns); cnt = zeros(LD_NB, ns)
    for l in 1:L, k in 1:L
        k == l && continue
        a = ld_bin(l - 1, L); r = abs(S[k] - S[l]); s = (ld_rclass(r) - 1) * LD_NB + ld_bin(k - 1, L)
        r2[a, s] += r^2; cnt[a, s] += 1
    end
    r2 ./= max.(cnt, 1)
    for l in 1:L
        a = ld_bin(l - 1, L)
        cur = -1; start = 0
        for k in 1:L
            s = k == l ? 0 : (ld_rclass(abs(S[k] - S[l])) - 1) * LD_NB + ld_bin(k - 1, L)
            if s != cur
                cur > 0 && (push!(lo, start); push!(hi, k - 1); push!(w, G.ρ[a, cur]); push!(w2, r2[a, cur]))
                cur = s; start = k
            end
        end
        cur > 0 && (push!(lo, start); push!(hi, L); push!(w, G.ρ[a, cur]); push!(w2, r2[a, cur]))
        push!(off, length(lo) + 1)
    end
    Segs(off, lo, hi, w, w2)
end

# R̄_l = Σ_k ρ̄_{s(k,l)} x_k (and with ρ̄² → class mean of ρ² if `sq`), by prefix sums
function rbin!(R, x, sg::Segs, P0; sq = false)
    P0[1] = 0.0
    @inbounds for l in eachindex(x)
        P0[l+1] = P0[l] + x[l]
    end
    ww = sq ? sg.w2 : sg.w
    @inbounds for l in eachindex(x)
        acc = 0.0
        for q in sg.off[l]:sg.off[l+1]-1
            acc += ww[q] * (P0[sg.hi[q]+1] - P0[sg.lo[q]])
        end
        R[l] = acc
    end
    R
end

"""
The recursion's recombination step at the survivors' state, with class-averaged ρ̄: per locus
Cov_S(R̄_l, x_l) (so that its RPR is Σ (Y − Cov(R̄,x))²).
"""
function binned_cov(B::STBuf, surv, ms, sg::Segs)
    L = B.L; n = length(surv)
    P0 = zeros(L + 1); rb = zeros(L); rbar = rbin!(zeros(L), ms.xbar, sg, P0)
    C = zeros(L)
    @inbounds for i in surv, c in 1:2
        m = 2(i - 1) + c
        rbin!(rb, view(B.X, :, m), sg, P0)
        @simd for l in 1:L
            C[l] += (B.X[l, m] - ms.xbar[l]) * (rb[l] - rbar[l])
        end
    end
    C ./= 2n
end

"""
Between-chromosome association of `inds` (weights `wt`, means of `mw`), weighted by ρ: per locus
Z_l = Σ_i w_i Σ_c d_l^{i,c}(R_l^{i,c'} − R̄_l)/(2Σw), and CovRG = Cov(R_l, G_chr) (chromosome level).
"""
function between(B::STBuf, inds, wt, mw)
    L = B.L
    rbar = rsum!(B.rbar, mw.xbar, B.S)
    Z = zeros(L); CovRG = zeros(L); W = 0.0
    r1, r2 = B.r1, B.r2
    @inbounds for i in inds
        m1 = 2(i - 1) + 1; m2 = m1 + 1; wi = wt[i]; W += 2wi
        rsum!(r1, view(B.X, :, m1), B.S); rsum!(r2, view(B.X, :, m2), B.S)
        g1 = B.G[m1] - mw.Gbar; g2 = B.G[m2] - mw.Gbar
        @simd for l in 1:L
            d1 = B.X[l, m1] - mw.xbar[l]; d2 = B.X[l, m2] - mw.xbar[l]
            e1 = r1[l] - rbar[l]; e2 = r2[l] - rbar[l]
            Z[l] += wi * (d1 * e2 + d2 * e1)
            CovRG[l] += wi * (e1 * g1 + e2 * g2)
        end
    end
    Z ./= W; CovRG ./= W
    (; Z, CovRG)
end

"""
E[Σ_l Ŷ_l² | adults] over `nmc` independent redraws of Bernoulli(π_i) survival (own RNG).
Subset sums of per-individual statistics (centred at the adults' means).
"""
function survival_mc(B::STBuf, N, π, ma, rng; nmc = NMC)
    L = B.L
    I = [rand(rng) < π[i] for i in 1:N, r in 1:nmc]
    s1 = zeros(L, nmc); s2 = zeros(L, nmc); s3 = zeros(L, nmc)
    g1 = zeros(nmc); w = zeros(nmc)
    @inbounds for i in 1:N
        for c in 1:2
            m = 2(i - 1) + c
            dG = B.G[m] - ma.Gbar
            for r in 1:nmc
                I[i, r] || continue
                g1[r] += dG; w[r] += 1
                @simd for l in 1:L
                    d = B.X[l, m] - ma.xbar[l]
                    s1[l, r] += d; s2[l, r] += d * d; s3[l, r] += d * dG
                end
            end
        end
    end
    tot = 0.0; k = 0
    for r in 1:nmc
        w[r] < 4 && continue
        k += 1
        gm = g1[r] / w[r]
        acc = 0.0
        @inbounds for l in 1:L
            xm = s1[l, r] / w[r]
            v = s2[l, r] / w[r] - xm^2
            Y = s3[l, r] / w[r] - xm * gm - v
            acc += Y^2
        end
        tot += acc
    end
    k > 0 ? tot / k : NaN
end

const ST_KEYS = [:gen, :N, :n, :Noff, :Sc, :beta, :VGchr, :epsS, :m,
    # adults
    :Y2_A, :mu_A, :VBv_A, :Y_A, :vY_A, :v_A, :v2_A, :cY_A,
    # selection: π-weighted, and the recursion's Bulmer step at the adults' state
    :Y2_P, :mu_P, :VBv_P, :Y_P, :Y2_bul, :Y2_bulv,
    # survival
    :Y2_b1, :Y2_full, :Y2_bind, :Y2_k3, :Y2_k4, :ZA2, :ZP2, :Zth2, :ZPth, :ZXsel, :ZXth, :ZXpair,
    :Y2_RPRb,
    :Y2_S, :Y2_MC, :mu_S, :VBv_S, :Y_S,
    # gamete pool
    :Y2_G, :Y2_RPR, :Z2_pred, :Z2, :ZX, :Y_G,
    # offspring before mutation
    :Y2_O, :mu_O, :VBv_O, :Y_O, :v_O, :VB_O]

function st_run(p, seed; from = ST_FROM, to = ST_TO, nmc = NMC)
    BLAS.set_num_threads(1)
    q = kernel_params(p)
    q.Selection === :viability || error("viability only")
    rng = MersenneTwister(seed)
    mcrng = Xoshiro(seed + 7)
    w = equilibrium_world(q, rng; burnin = KBURN)
    L = q.NoQTL
    B = STBuf(L, q.PopSize, q.RecombPoissonMean)
    sg = class_segments(L, q.RecombPoissonMean)
    VE = QTLSim.envVar(q)
    mloc = q.MutationScalar * VE / (2L)            # mutational variance per locus per chromosome
    rec = Dict{Symbol,Vector{Float64}}(k => Float64[] for k in ST_KEYS)
    vgser = Float64[]
    extinct_gen = -1
    for it in w.gen:q.NoGenerations-1
        it >= to && break
        w.gen = it
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        N = w.N
        inst = it >= from
        local vals
        if inst
            g = view(w.genotype, 1:N)
            push!(vgser, var(g))
            gb = mean(g); vg = var(g; corrected = false)
            ω2 = w.ω^2; Sg = ω2 + VE
            π = [C0 * sqrt(ω2 / Sg) * exp(-(w.genotype[i] - w.θ)^2 / (2Sg)) for i in 1:N]
            πbar = mean(π)
            β = sum((π ./ πbar .- 1) .* (g .- gb)) / N / vg
            load_pool!(B, w.pool, w.chr, N)
            one = ones(N)
            ma = moments(B, 1:N, one)
            Sc = ω2 + VE + ma.VG
            mp = moments(B, 1:N, π)
            surv = copy(w.survivors); n = length(surv)
            ms = n >= 2 ? moments(B, surv, one) : nothing
            # the recursion's Bulmer step on Y at the adults' state (per-locus c_l, v_l; and with the
            # group-mean v̄ in the decay term, as y_step!)
            Ybul = ma.Y .- ((ma.c .+ ma.v) .* ma.Y .+ ma.c .* ma.v) ./ Sc
            vbar = similar(ma.v)
            for a in 1:LD_NB
                idx = [l for l in 1:L if ld_bin(l - 1, L) == a]
                isempty(idx) || (vbar[idx] .= mean(ma.v[idx]))
            end
            Ybulv = ma.Y .- ((ma.c .+ vbar) .* ma.Y .+ ma.c .* ma.v) ./ Sc
            # first-order selection with the cumulant terms (chromosome-level tilt e^{bu − u²/2S_c})
            bt = (w.θ - gb) / Sc
            dYb = -ma.c .* ma.cx ./ Sc
            dYf = dYb .+ bt .* ma.k3 .+ ((bt^2 - 1 / Sc) / 2) .* ma.k4
            Sind = ω2 + VE + vg
            # between-chromosome association created by selection (π pass) and its prediction
            bA = between(B, 1:N, one, ma); bP = between(B, 1:N, π, mp)
            Zth = -bA.CovRG .* ma.cx ./ Sc
            epsS = sum(π .* (1 .- π)) / (2 * sum(π)^2)
            Y2MC = survival_mc(B, N, π, ma, mcrng; nmc)
            gp = ms === nothing ? nothing : gamete_pool(B, surv, ms)
            vals = Dict{Symbol,Float64}(:gen => it, :N => N, :n => n, :Sc => Sc, :beta => β, :VGchr => ma.VG, :epsS => epsS, :m => mloc,
                :Y2_A => sum(abs2, ma.Y), :mu_A => sum(ma.mu22), :VBv_A => sum(ma.VarB .* ma.v), :Y_A => sum(ma.Y),
                :vY_A => sum(ma.v .* ma.Y), :v_A => sum(ma.v), :v2_A => sum(abs2, ma.v), :cY_A => sum(ma.c .* ma.Y),
                :Y2_P => sum(abs2, mp.Y), :mu_P => sum(mp.mu22), :VBv_P => sum(mp.VarB .* mp.v), :Y_P => sum(mp.Y),
                :Y2_bul => sum(abs2, Ybul), :Y2_bulv => sum(abs2, Ybulv), :Y2_MC => Y2MC,
                :Y2_b1 => sum(abs2, ma.Y .+ dYb), :Y2_full => sum(abs2, ma.Y .+ dYf), :Y2_bind => sum(abs2, ma.Y .- ma.c .* ma.cx ./ Sind),
                :Y2_k3 => 2sum(ma.Y .* bt .* ma.k3), :Y2_k4 => 2sum(ma.Y .* ((bt^2 - 1 / Sc) / 2) .* ma.k4),
                :ZA2 => sum(abs2, bA.Z), :ZP2 => sum(abs2, bP.Z), :Zth2 => sum(abs2, Zth), :ZPth => sum(bP.Z .* Zth))
            if ms !== nothing
                vals[:Y2_S] = sum(abs2, ms.Y); vals[:mu_S] = sum(ms.mu22); vals[:VBv_S] = sum(ms.VarB .* ms.v); vals[:Y_S] = sum(ms.Y)
                vals[:Y2_G] = sum(abs2, gp.Ygp); vals[:Y2_RPR] = sum(abs2, ms.Y .- gp.CovRx)
                vals[:Z2_pred] = sum((ms.v .* gp.VarR .+ gp.CovRx .^ 2) ./ (2n)); vals[:Z2] = sum(abs2, gp.Z)
                vals[:ZX] = 2sum((ms.Y .- gp.CovRx) .* gp.Z); vals[:Y_G] = sum(gp.Ygp)
                vals[:ZXsel] = 2sum((ms.Y .- gp.CovRx) .* bP.Z); vals[:ZXth] = 2sum((ms.Y .- gp.CovRx) .* Zth)
                vals[:ZXpair] = -2sum((ms.Y .- gp.CovRx) .* gp.CovRx) / (2n - 1)
                vals[:Y2_RPRb] = sum(abs2, ms.Y .- binned_cov(B, surv, ms, sg))
            end
        end
        if w.extinct
            extinct_gen = it
            inst && (for k in ST_KEYS; push!(rec[k], get(vals, k, NaN)); end)
            break
        end
        QTLSim.reproduce!(w, rng)
        if w.extinct
            extinct_gen = it
            inst && (for k in ST_KEYS; push!(rec[k], get(vals, k, NaN)); end)
            break
        end
        if inst
            load_pool!(B, w.pool, w.schr, w.Noff)
            mo = moments(B, 1:w.Noff, ones(w.Noff))
            vals[:Noff] = w.Noff
            vals[:Y2_O] = sum(abs2, mo.Y); vals[:mu_O] = sum(mo.mu22); vals[:VBv_O] = sum(mo.VarB .* mo.v)
            vals[:Y_O] = sum(mo.Y); vals[:v_O] = sum(mo.v); vals[:VB_O] = sum(mo.VarB)
            for k in ST_KEYS
                push!(rec[k], get(vals, k, NaN))
            end
        end
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    (; seed, extinct_gen, rec, vgser, L, u = q.RecombPoissonMean, PopSize = q.PopSize)
end

st_file(id, r, set) = joinpath(LL_RAW, replace(id, "/" => "__", ".txt" => "") * "@stage" * set * "#$r.jls")

function st_check(; id = "QT1000-Possion3.0/Input-09697.txt", j = 89, r = 1)
    p = readp(id)
    t = @elapsed tr = st_run(p, KSEED + 100j + r; to = ST_FROM + 8)
    o = simulate(kernel_params(p); seed = KSEED + 100j + r, init = :equilibrium, burnin = KBURN)
    k = [findfirst(==(g), o.gen) for g in ST_FROM:ST_FROM+length(tr.vgser)-1]
    same = all(o.geno_var[k] .== tr.vgser)
    # Σ Y² at the adults agrees with linked_source's Σ yY (same run, same generations)
    t2 = ls_run(p, KSEED + 100j + r; to = ST_FROM + 8, decompose = false)
    agree = all(abs(tr.rec[:Y2_A][kk] - tot(t2.rec[:yY][kk])) <= 1e-4 * abs(tr.rec[:Y2_A][kk]) for kk in eachindex(tr.vgser))
    # mutation stage closes: next adults' Σ Y² is measured from the offspring after mutation
    @printf("check %s: %d generations, V_G bitwise identical to simulate: %s; Σ Y² as linked_source's Σ yY: %s (%.1f s)\n",
            id, length(tr.vgser), same, agree, t)
    (same && agree) || error("check failed")
end

ll_rows(set) = set == "H" ? hrows() : ls_rows(set)
ll_seed(set, j, r) = set == "H" ? HSEED + 100j + r : ls_seed(set, j, r)

function st_kernel(set, R)
    mkpath(LL_RAW)
    rows = ll_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(st_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) stage-instrumented runs for set $set, R = $R, $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time(); cpu = Threads.Atomic{Float64}(0.0)
    Threads.@threads :greedy for (id, j, k) in jobs
        t1 = time()
        tr = st_run(readp(id), ll_seed(set, j, k))
        serialize(st_file(id, k, set), tr)
        Threads.atomic_add!(cpu, time() - t1)
        d = Threads.atomic_add!(done, 1) + 1
        @printf("  %d/%d %s #%d  %.0f s\n", d, length(jobs), id, k, time() - t0)
    end
    @printf("done: %d runs, %.0f s wall, %.0f thread-s\n", length(jobs), time() - t0, cpu[])
end

# ---------------------------------------------------------------------------------------
# the per-stage accounting

"""
Per-generation stage changes of Σ Y² (measured) and their predictions, for one run. Each
prediction uses the measured state at the start of its stage. Columns per stage: `meas`, `exact`
(first-order expectation, all terms), `recur` (the recursion's formula), and variants.
"""
function stage_terms(tr)
    r = tr.rec
    t = Int.(r[:gen]) .- MOVE
    n = length(t)
    L = tr.L
    N = r[:N]; nn = r[:n]; No = r[:Noff]
    εS = 1 ./ (2nn) .- 1 ./ (2N)            # realised-n factor (the recursion's)
    εB = r[:epsS]                            # Bernoulli(π) factor
    εO = 1 ./ (2No)
    nxt(k) = k < n && t[k+1] == t[k] + 1
    cols = Dict{Symbol,Vector{Float64}}(
        :sel_meas => r[:Y2_P] .- r[:Y2_A],
        :sel_recur => r[:Y2_bulv] .- r[:Y2_A],              # y_step!'s Bulmer form (group-mean v̄)
        :sel_bul => r[:Y2_b1] .- r[:Y2_A],                  # −c_B c_x/S_c, unlinearised
        :sel_full => r[:Y2_full] .- r[:Y2_A],               # + b κ(B,x,G) + (b² − 1/S_c)κ(B,x,G,G)/2
        :sel_ind => r[:Y2_bind] .- r[:Y2_A],                # −c_B c_x/(ω² + V_E + V_G): Gaussian conditioning
        :sel_k3 => r[:Y2_k3], :sel_k4 => r[:Y2_k4],
        :surv_meas => r[:Y2_S] .- r[:Y2_P],
        :surv_exact => r[:Y2_MC] .- r[:Y2_P],
        :surv_first => -3 .* εB .* r[:Y2_P] .+ εB .* r[:mu_P],
        :surv_recur => ((1 .- εS) .^ 2 .- 1) .* r[:Y2_P] .+ εS .* r[:VBv_P],
        :rec_meas => r[:Y2_G] .- r[:Y2_S],
        :rec_recur => r[:Y2_RPR] .- r[:Y2_S],               # P ← RPR with exact pairwise ρ
        :rec_recurb => haskey(r, :Y2_RPRb) ? r[:Y2_RPRb] .- r[:Y2_S] : fill(NaN, n),   # with the classes' mean ρ̄
        :rec_exact => r[:Y2_RPR] .- r[:Y2_S] .+ r[:Z2_pred] .+ r[:ZP2] .+ r[:ZXsel] .+ r[:ZXpair],
        :rec_theory => r[:Y2_RPR] .- r[:Y2_S] .+ r[:Z2_pred] .+ r[:Zth2] .+ r[:ZXth] .+ r[:ZXpair],
        :rec_Znoise => r[:Z2_pred], :rec_Zsel => r[:ZXsel], :rec_Zth => r[:ZXth], :rec_Zpair => r[:ZXpair],
        :rec_Z2 => r[:Z2], :rec_ZX => r[:ZX],
        :off_meas => r[:Y2_O] .- r[:Y2_G],
        :off_exact => -3 .* εO .* r[:Y2_G] .+ εO .* r[:mu_O],
        :off_recur => ((1 .- εO) .^ 2 .- 1) .* r[:Y2_G] .+ εO .* r[:VBv_O],
        :mut_meas => [nxt(k) ? r[:Y2_A][k+1] - r[:Y2_O][k] : NaN for k in 1:n],
        :mut_exact => [(r[:m][k] * (L - 1) * r[:v_O][k] + r[:m][k] * r[:VB_O][k]) / (2No[k]) for k in 1:n],
        :mut_recur => [r[:m][k] * r[:VB_O][k] / (2No[k]) for k in 1:n],
        :Y2 => r[:Y2_A], :ZA2 => r[:ZA2], :ZP2 => r[:ZP2], :Zth2 => r[:Zth2], :ZPth => r[:ZPth])
    (; t, cols)
end

const STAGES = (:sel, :surv, :rec, :off, :mut)
const SHOW = [:sel_meas, :sel_recur, :sel_bul, :sel_ind, :sel_full, :sel_k3, :sel_k4,
              :surv_meas, :surv_exact, :surv_first, :surv_recur,
              :rec_meas, :rec_exact, :rec_theory, :rec_recur, :rec_recurb, :rec_Znoise, :rec_Zsel, :rec_Zth, :rec_Zpair,
              :off_meas, :off_exact, :off_recur, :mut_meas, :mut_exact, :mut_recur]

function st_stages(set, R; io = stdout, wins = [(-60, -1), (0, 49), (50, 99)], verbose = true)
    rows = ll_rows(set)
    println(io, "\n== stage accounting of Σ_l Y_l², set $set, R = $R: mean change per generation / mean Σ Y² (window means over runs alive) ==")
    agg = Dict{Tuple{Symbol,Tuple{Int,Int}},Vector{Float64}}()
    for rw in rows
        runs = [deserialize(st_file(rw.id, k, set)) for k in 1:R if isfile(st_file(rw.id, k, set))]
        isempty(runs) && continue
        T = [stage_terms(tr) for tr in runs]
        verbose && @printf(io, "\n%s  L=%d u=%g runs %d\n", rw.id, runs[1].L, runs[1].u, length(runs))
        for wdw in wins
            y2 = winmean([(x.t, x.cols[:Y2]) for x in T], wdw)
            verbose && @printf(io, "  %4d..%-4d ΣY² %.3e", wdw..., y2)
            for c in SHOW
                v = winmean([(x.t, x.cols[c]) for x in T], wdw) / y2
                push!(get!(agg, (c, wdw), Float64[]), v)
                verbose && @printf(io, " %s %+.4f", c, v)
            end
            verbose && println(io)
        end
    end
    println(io, "\n== over inputs: mean (median) of change per generation / Σ Y² ==")
    for wdw in wins
        @printf(io, "  window %d..%d (%d inputs)\n", wdw..., length(agg[(:sel_meas, wdw)]))
        for c in SHOW
            v = agg[(c, wdw)]
            @printf(io, "    %-11s %+.4f (%+.4f)\n", c, mean(filter(isfinite, v)), medf(v))
        end
        for kind in (:meas, :exact, :recur)
            tot = [sum(agg[(Symbol(s, :_, kind == :exact && s == :sel ? :meas : kind), wdw)][i] for s in STAGES) for i in eachindex(agg[(:sel_meas, wdw)])]
            @printf(io, "    total %-6s %+.4f (%+.4f)\n", kind, mean(filter(isfinite, tot)), medf(tot))
        end
    end
    agg
end

# ---------------------------------------------------------------------------------------
# the corrected recursion driven by the runs (levels B and B′ of linked_decay/linked_source)

include(joinpath(@__DIR__, "linked_ld_theory.jl"))

const FIX_ALL = (between = true, pairing = true, mutfull = true, sind = true, vhet = true)
const FIXES = Dict(
    :B0 => (;),                                                        # linked_decay's level B
    :Bv => (vhet = true,),                                             # only the heterogeneity weight
    :Bx => (between = true, pairing = true, mutfull = true, sind = true),   # all but the heterogeneity weight
    :Bs => FIX_ALL,                                                    # the corrected level B (B*)
    :Bp0 => (couple = true, hsel = :bulmer),                           # linked_source's level B′
    :Bps => merge((couple = true, hsel = :bulmer), FIX_ALL))           # B′ fed by the corrected y (B′*)

"""
Level B / B′ with the corrections `kw` (`yfix_step!`): from t = −60 the stationary state for the
drivers averaged over −60..−41, T = 0; each generation the source (B: −2 P1; B′: Σ_t (H − 2P))
feeds linked_decay's T recursion, then `yfix_step!`. Each group is carried on its active classes only
(exact: the others stay zero). Returns per generation Σ K, Σ C, Σ E[yY], Σ H.
"""
function levelBfix(tr, G::LDGroups, ρ2; kw...)
    rec = tr.rec
    n = length(rec[:gen])
    couple = get(kw, :couple, false)
    kw0 = merge((; kw...), (hsel = :none,))
    d0 = drivers(rec, 1; avg = PRE_AVG)
    NA = size(G.n, 1)
    act = [a for a in 1:NA if G.nb[a] > 0]
    idx = Dict(a => findall(>(0), G.n[a, :]) for a in act)
    Gr = Dict(a => LDGroups(G.L, G.n[:, idx[a]], G.ρ[:, idx[a]], G.nb, G.bofs[idx[a]]) for a in act)
    sr = Dict(a => LSState(Gr[a]) for a in act)
    ρ2r = Dict(a => ρ2[:, idx[a]] for a in act)
    Tr = Dict(a => zeros(length(idx[a])) for a in act)
    for a in act    # the stationary state (yfix_stationary! on one group)
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
    out = (K = fill(NaN, n), C = fill(NaN, n), yY = fill(NaN, n), H = fill(NaN, n))
    for k in 1:n
        d = drivers(rec, k)
        Kt = 0.0; Ct = 0.0; yt = 0.0; Ht = 0.0
        for a in act
            g = Gr[a]; st = sr[a]
            Ks = couple ? source_K(st, a) : -2 .* vec(sum(st.y.P[a]; dims = 2))
            Kt += sum(Ks); Ct += sum(Tr[a]); yt += sum(st.y.P[a]); Ht += sum(st.H[a])
            T_step!(Tr[a], g, d, a, Ks; mode = :derived)
            yfix_step!(st, g, d, a; ρ2 = ρ2r[a], kw...)
        end
        out.K[k] = Kt; out.C[k] = Ct; out.yY[k] = yt; out.H[k] = Ht
    end
    out
end

# runs with association records: linked_source's raw (D, G) or linked_decay's (F), in linked_decay's names
as_file(id, r, set) = joinpath(LL_RAW, replace(id, "/" => "__", ".txt" => "") * "@assoc" * set * "#$r.jls")
function assoc_runs(id, R, set)
    if set == "H"
        return [ldview(deserialize(as_file(id, k, set))) for k in 1:R if isfile(as_file(id, k, set))]
    end
    if set == "F" || set == "Fd"
        return [deserialize(ld_file(id, k, "F")) for k in 1:R if isfile(ld_file(id, k, "F"))]
    end
    [ldview(deserialize(ls_file(id, k, set))) for k in 1:R if isfile(ls_file(id, k, set))]
end
has_H(tr) = haskey(tr.rec, :H)

function fix_input(r, R, set; variants = keys(FIXES), io = stdout, withA = false)
    runs = assoc_runs(r.id, R, set)
    isempty(runs) && return nothing
    G = ld_groups(runs[1].L, runs[1].u); ρ2 = class_rho2(runs[1].L, runs[1].u)
    ser = Dict{Symbol,Vector{Any}}()
    cum = Dict{Symbol,Vector{Float64}}()
    for tr in runs
        m = measured(tr)
        push!(get!(ser, :yY, []), (m.t, m.yY)); push!(get!(ser, :C, []), (m.t, m.C)); push!(get!(ser, :K, []), (m.t, m.K))
        has_H(tr) && push!(get!(ser, :H, []), (m.t, m.K .+ 2 .* m.yY))
        k0 = findfirst(==(-1), m.t); k100 = findfirst(==(100), m.t)
        g0 = k0 === nothing ? NaN : mean(m.genic[max(1, k0 - 19):k0])
        ok = k0 !== nothing && k100 !== nothing
        if ok
            rg = k0+1:k100-1
            push!(get!(cum, :meas_loss, Float64[]), sum(2 .* m.beta[rg] .* m.C[rg]) / g0)
            push!(get!(cum, :meas_losssym, Float64[]), sum(2 .* m.beta[rg] .* m.Csym[rg]) / g0)
        end
        if withA
            for md in AMODES
                push!(get!(ser, Symbol(:A_, md), []), (m.t, levelA(tr, G; mode = md)))
            end
        end
        for v in variants
            o = levelBfix(tr, G, ρ2; FIXES[v]...)
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
    @printf(io, "%-34s L=%d runs %d  yY ratio (pre / 50..99):", r.id, runs[1].L, length(runs))
    for v in variants
        @printf(io, "  %s %.3f/%.3f", v, res[Symbol(v, :_yY)][(-60, -1)] / res[:yY][(-60, -1)], res[Symbol(v, :_yY)][(50, 99)] / res[:yY][(50, 99)])
    end
    println(io)
    res
end

function fix_dev(set, R; variants = [:B0, :Bv, :Bx, :Bs, :Bp0, :Bps], io = stdout)
    rows = ll_rows(set == "Fd" ? "F" : set)
    println(io, "\n== corrected recursion driven by the runs, set $set, R = $R: Σ E[yY] predicted / measured ==")
    bufs = [IOBuffer() for _ in rows]
    resv = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        resv[i] = fix_input(rows[i], R, set; variants, io = bufs[i])
    end
    foreach(b -> print(io, String(take!(b))), bufs)
    all = [(rows[i], resv[i]) for i in eachindex(rows) if resv[i] !== nothing]
    elig = filter(x -> isfinite(x[2][:meas_loss]) && x[2][:meas_loss] <= -0.05, all)
    @printf(io, "medians over %d inputs (%d eligible):\n", length(all), length(elig))
    for v in variants
        @printf(io, "  %-6s yY", v)
        for w in WINS
            @printf(io, "  %d..%d %.3f", w..., medf([x[2][Symbol(v, :_yY)][w] / x[2][:yY][w] for x in all]))
        end
        if haskey(all[1][2], :H)
            @printf(io, " | H")
            for w in WINS
                @printf(io, "  %d..%d %.3f", w..., medf([x[2][Symbol(v, :_H)][w] / x[2][:H][w] for x in all]))
            end
        end
        @printf(io, " | K")
        for w in WINS
            @printf(io, "  %d..%d %.3f", w..., medf([x[2][Symbol(v, :_K)][w] / x[2][:K][w] for x in all]))
        end
        @printf(io, " | ΣC lr (eligible)")
        for w in WINS[3:5]
            u = [lr(x[2][Symbol(v, :_C)][w], x[2][:C][w]) for x in elig]
            @printf(io, "  %d..%d %+.3f/%.3f", w..., medf(u), mabs(u))
        end
        u = [lr(x[2][Symbol(v, :_loss)], x[2][:meas_losssym]) for x in elig]
        @printf(io, "  loss %+.3f/%.3f\n", medf(u), mabs(u))
    end
    all
end

# ---------------------------------------------------------------------------------------
# exact accounting in the recursion's own class space: per locus the class vector y_l (s = partner
# bin × split class, as `ld_groups`), at each stage; the recursion's stage map applied to the
# measured state at the start of the stage, against the measured state at its end.

"""
Per locus l the partner classes as contiguous segments of k (ρ_kl is monotone on each side of l),
each with its side (−1 left, +1 right) and local class index; `cls` maps local to global class s.
"""
struct CSeg
    off::Vector{Int}; lo::Vector{Int}; hi::Vector{Int}; side::Vector{Int8}; lq::Vector{Int}
    coff::Vector{Int}; cls::Vector{Int}
end
function class_seg(L, u)
    S = qtsim_survival(L, u)
    off = [1]; lo = Int[]; hi = Int[]; side = Int8[]; lq = Int[]; coff = [1]; cls = Int[]
    for l in 1:L
        loc = Dict{Int,Int}()
        cur = -1; start = 0
        flush(kend) = if cur > 0
            q = get!(loc, cur) do
                push!(cls, cur); length(loc) + 1
            end
            push!(lo, start); push!(hi, kend); push!(side, start > l ? Int8(1) : Int8(-1)); push!(lq, q)
        end
        for k in 1:L
            s = k == l ? 0 : (ld_rclass(abs(S[k] - S[l])) - 1) * LD_NB + ld_bin(k - 1, L)
            if s != cur
                flush(k - 1)
                cur = s; start = k
            end
        end
        flush(L)
        push!(off, length(lo) + 1); push!(coff, length(cls) + 1)
    end
    CSeg(off, lo, hi, side, lq, coff, cls)
end

# class sums of one centred chromosome d: Bv[coff[l] + q − 1] = Σ_{k ∈ class q of l} d_k, and Rv the ρ-weighted sums
function class_sums!(Bv, Rv, d, cs::CSeg, S, P0, P1)
    L = length(d)
    P0[1] = 0.0; P1[1] = 0.0
    @inbounds for l in 1:L
        P0[l+1] = P0[l] + d[l]; P1[l+1] = P1[l] + S[l] * d[l]
    end
    fill!(Bv, 0.0); Rv === nothing || fill!(Rv, 0.0)
    @inbounds for l in 1:L
        base = cs.coff[l] - 1; Sl = S[l]
        for q in cs.off[l]:cs.off[l+1]-1
            a, b = cs.lo[q], cs.hi[q]
            s0 = P0[b+1] - P0[a]
            Bv[base+cs.lq[q]] += s0
            if Rv !== nothing
                s1 = P1[b+1] - P1[a]
                Rv[base+cs.lq[q]] += cs.side[q] < 0 ? s1 - Sl * s0 : Sl * s0 - s1
            end
        end
    end
    nothing
end

mutable struct CSBuf
    cs::CSeg; T::Int
    d::Vector{Float64}; d2::Vector{Float64}; P0::Vector{Float64}; P1::Vector{Float64}
    B1::Vector{Float64}; R1::Vector{Float64}; B2::Vector{Float64}; R2::Vector{Float64}
end
CSBuf(cs::CSeg, L) = (T = length(cs.cls); CSBuf(cs, T, zeros(L), zeros(L), zeros(L + 1), zeros(L + 1), zeros(T), zeros(T), zeros(T), zeros(T)))

"""
Class vectors per locus (flattened, `cs.coff`) of the chromosomes of `inds` weighted by `wt`, centred
at the means `mw` (from `moments`): y (Cov(B_s, x_l)); with `pairs`, also CovR (Cov(R_s, x_l)) and Z
(the cross-chromosome part Σ_i (d¹ e² + d² e¹)/(2Σw)).
"""
function class_y(B::STBuf, C::CSBuf, inds, wt, mw; pairs = false)
    L = B.L; cs = C.cs
    y = zeros(C.T); cR = pairs ? zeros(C.T) : nothing; Z = pairs ? zeros(C.T) : nothing
    W = 0.0
    @inbounds for i in inds
        wi = wt[i]; W += 2wi
        m1 = 2(i - 1) + 1
        for l in 1:L
            C.d[l] = B.X[l, m1] - mw.xbar[l]; C.d2[l] = B.X[l, m1+1] - mw.xbar[l]
        end
        class_sums!(C.B1, pairs ? C.R1 : nothing, C.d, cs, B.S, C.P0, C.P1)
        class_sums!(C.B2, pairs ? C.R2 : nothing, C.d2, cs, B.S, C.P0, C.P1)
        for l in 1:L
            d1 = C.d[l]; d2 = C.d2[l]
            for j in cs.coff[l]:cs.coff[l+1]-1
                y[j] += wi * (d1 * C.B1[j] + d2 * C.B2[j])
                if pairs
                    cR[j] += wi * (d1 * C.R1[j] + d2 * C.R2[j])
                    Z[j] += wi * (d1 * C.R2[j] + d2 * C.R1[j])
                end
            end
        end
    end
    y ./= W
    pairs && (cR ./= W; Z ./= W)
    (; y, cR, Z)
end

# drivers of the recursion from the adults (as linked_decay's ld_measure!)
function cs_drivers(B::STBuf, N, ma, binl; Sc, n, Nn, Vm, β = 0.0)
    NB = LD_NB; M = 2N
    V = zeros(NB); V2 = zeros(NB)
    for l in 1:B.L
        V[binl[l]] += ma.v[l]; V2[binl[l]] += ma.v[l]^2
    end
    Bs = zeros(NB); cb = zeros(NB); Cbb = zeros(NB, NB)
    for m in 1:M
        fill!(Bs, 0.0)
        for l in 1:B.L
            Bs[binl[l]] += B.X[l, m] - ma.xbar[l]
        end
        dG = B.G[m] - ma.Gbar
        cb .+= Bs .* dG; Cbb .+= Bs * Bs'
    end
    cb ./= M; Cbb ./= M
    LDDrivers(β = β, Sc = Sc, dS = 0.0, n = n, N = Nn, Vm = Vm, V = V, V2 = V2, k3 = zeros(NB), g = copy(V), c = cb, Cbb = Cbb)
end

# summaries of per-locus class vectors: Σ_l (Σ_s y)², and the 7×7 split-class matrix Σ_l (Σ_{s∈r} y)(Σ_{s∈r'} y)
function cs_summ(vecs, cs::CSeg, L)
    tot = 0.0; rr = zeros(LD_NR, LD_NR); rv = zeros(LD_NR)
    @inbounds for l in 1:L
        fill!(rv, 0.0)
        for j in cs.coff[l]:cs.coff[l+1]-1
            rv[(cs.cls[j] - 1) ÷ LD_NB + 1] += vecs[j]
        end
        tot += sum(rv)^2
        rr .+= rv * rv'
    end
    (tot, rr)
end
# summary of a class-pair matrix M (ns × ns) for group a: total and 7×7 by split class
function cs_summ_mat(M)
    rr = zeros(LD_NR, LD_NR)
    for s in axes(M, 1), t in axes(M, 2)
        rr[(s - 1) ÷ LD_NB + 1, (t - 1) ÷ LD_NB + 1] += M[s, t]
    end
    (sum(M), rr)
end

const CS_KEYS = [:sel_meas, :sel_recur, :surv_meas, :surv_recur, :rec_meas, :rec_recur, :rec_exact_rho, :rec_fix, :rec_fixb, :rec_fixp,
                 :off_meas, :off_recur, :mut_meas, :mut_recur, :mut_fix, :A, :Pi, :S, :Gp, :O]

function cs_run(p, seed; from = ST_FROM, to = ST_TO)
    BLAS.set_num_threads(1)
    q = kernel_params(p)
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    L = q.NoQTL; u = q.RecombPoissonMean
    B = STBuf(L, q.PopSize, u)
    cs = class_seg(L, u); C = CSBuf(cs, L)
    G = ld_groups(L, u); ρ2 = class_rho2(L, u)
    binl = [ld_bin(l - 1, L) for l in 1:L]
    VE = QTLSim.envVar(q); Vm = QTLSim.muteVar(q)
    rec = Dict{Symbol,Any}(k => Tuple{Float64,Matrix{Float64}}[] for k in CS_KEYS)
    rec[:gen] = Int[]
    vgser = Float64[]
    prevO = nothing; prevgen = -10
    extinct_gen = -1
    for it in w.gen:q.NoGenerations-1
        it >= to && break
        w.gen = it
        QTLSim.develop!(w); QTLSim.select!(w, rng)
        N = w.N
        inst = it >= from
        if inst
            g = view(w.genotype, 1:N); push!(vgser, var(g))
            ω2 = w.ω^2; Sg = ω2 + VE
            π = [C0 * sqrt(ω2 / Sg) * exp(-(w.genotype[i] - w.θ)^2 / (2Sg)) for i in 1:N]
            load_pool!(B, w.pool, w.chr, N)
            one = ones(N)
            ma = moments(B, 1:N, one); mp = moments(B, 1:N, π)
            Sc = ω2 + VE + ma.VG
            surv = copy(w.survivors); n = length(surv)
            n < 2 && (extinct_gen = it; break)
            ms = moments(B, surv, one)
            d = cs_drivers(B, N, ma, binl; Sc, n, Nn = q.PopSize, Vm)
            yA = class_y(B, C, 1:N, one, ma).y
            yP = class_y(B, C, 1:N, π, mp).y
            cy = class_y(B, C, surv, one, ms; pairs = true)
            yS = cy.y
            vals = Dict{Symbol,Tuple{Float64,Matrix{Float64}}}()
            vals[:A] = cs_summ(yA, cs, L); vals[:Pi] = cs_summ(yP, cs, L); vals[:S] = cs_summ(yS, cs, L)
            prevO !== nothing && prevgen == it - 1 && (vals[:mut_meas] = (vals[:A][1] - prevO[1][1], vals[:A][2] .- prevO[1][2]))
            prevO !== nothing && prevgen == it - 1 && (for k in (:mut_recur, :mut_fix); vals[k] = prevO[2][k]; end)
            # per-locus recursion maps
            ysel = similar(yA); ygr = similar(yA); yex = similar(yA); yfx = similar(yA); yfb = similar(yA); yfp = similar(yA)
            εs = 1 / (2n) - 1 / (2q.PopSize); fS = 1 - εs
            cvec = [class_c(G, d, a) for a in 1:LD_NB]
            vbar = [d.V[a] / max(G.nb[a], 1) for a in 1:LD_NB]
            noise_fix = 0.0; rr_fix = zeros(LD_NR, LD_NR)
            for l in 1:L
                a = binl[l]; c = cvec[a]; r = cs.coff[l]:cs.coff[l+1]-1
                Ysum = sum(view(yA, r))
                for j in r
                    s = cs.cls[j]; ρ = G.ρ[a, s]; R = 1 - ρ; Rp = R - ρ / (2n - 1)
                    Δ = -(c[s] * Ysum + vbar[a] * yA[j]) / Sc - c[s] * ma.v[l] / Sc       # selection's change (within and between)
                    ysel[j] = yA[j] + Δ
                    ygr[j] = R * yS[j]
                    yex[j] = yS[j] - cy.cR[j]
                    yfx[j] = Rp * yS[j] + ρ * fS * Δ
                    yfb[j] = R * yS[j] + ρ * fS * Δ          # between only
                    yfp[j] = Rp * yS[j]                       # pairing only (mean part; noise added below)
                end
            end
            # pairing noise per group: (V_a CovR + Σ_l (ρ∘y_S)(ρ∘y_S)ᵀ)/(2n)
            pn_tot = 0.0; pn_rr = zeros(LD_NR, LD_NR)
            for a in 1:LD_NB
                G.nb[a] == 0 && continue
                Σ = class_cov(G, d, a); CR = class_covR(G, d, a, ρ2, Σ)
                t1, r1 = cs_summ_mat(CR .* (d.V[a] / (2n)))
                pn_tot += t1; pn_rr .+= r1
            end
            ρy = similar(yS)
            for l in 1:L, j in cs.coff[l]:cs.coff[l+1]-1
                ρy[j] = G.ρ[binl[l], cs.cls[j]] * yS[j]
            end
            t2, r2 = cs_summ(ρy, cs, L)
            pn_tot += t2 / (2n); pn_rr .+= r2 ./ (2n)
            # sampling seeds (Σ_a V_a Σ_B)
            sd_tot = 0.0; sd_rr = zeros(LD_NR, LD_NR); md_rr = zeros(LD_NR, LD_NR); md_tot = 0.0
            mloc = (Vm / 2) / L
            for a in 1:LD_NB
                G.nb[a] == 0 && continue
                t, rmat = cs_summ_mat(class_cov(G, d, a) .* d.V[a]); sd_tot += t; sd_rr .+= rmat
                for s in 1:size(G.n, 2)
                    rs = (s - 1) ÷ LD_NB + 1
                    md_rr[rs, rs] += d.V[a] * G.n[a, s] * mloc; md_tot += d.V[a] * G.n[a, s] * mloc
                end
            end
            vals[:sel_meas] = (vals[:Pi][1] - vals[:A][1], vals[:Pi][2] .- vals[:A][2])
            s1 = cs_summ(ysel, cs, L); vals[:sel_recur] = (s1[1] - vals[:A][1], s1[2] .- vals[:A][2])
            vals[:surv_meas] = (vals[:S][1] - vals[:Pi][1], vals[:S][2] .- vals[:Pi][2])
            vals[:surv_recur] = ((fS^2 - 1) * vals[:Pi][1] + εs * sd_tot, (fS^2 - 1) .* vals[:Pi][2] .+ εs .* sd_rr)
            # gamete pool, measured
            ygp = yS .- cy.cR .+ cy.Z
            vals[:Gp] = cs_summ(ygp, cs, L)
            vals[:rec_meas] = (vals[:Gp][1] - vals[:S][1], vals[:Gp][2] .- vals[:S][2])
            for (k, yy, extra) in ((:rec_recur, ygr, false), (:rec_exact_rho, yex, false), (:rec_fix, yfx, true), (:rec_fixb, yfb, false), (:rec_fixp, yfp, true))
                s2 = cs_summ(yy, cs, L)
                vals[k] = extra ? (s2[1] + pn_tot - vals[:S][1], s2[2] .+ pn_rr .- vals[:S][2]) : (s2[1] - vals[:S][1], s2[2] .- vals[:S][2])
            end
            QTLSim.reproduce!(w, rng)
            w.extinct && (extinct_gen = it; break)
            load_pool!(B, w.pool, w.schr, w.Noff)
            mo = moments(B, 1:w.Noff, ones(w.Noff))
            yO = class_y(B, C, 1:w.Noff, ones(w.Noff), mo).y
            vals[:O] = cs_summ(yO, cs, L)
            εo = 1 / (2w.Noff)
            vals[:off_meas] = (vals[:O][1] - vals[:Gp][1], vals[:O][2] .- vals[:Gp][2])
            vals[:off_recur] = (((1 - εo)^2 - 1) * vals[:Gp][1] + εo * sd_tot, ((1 - εo)^2 - 1) .* vals[:Gp][2] .+ εo .* sd_rr)
            mr = (mloc * sd_tot / (2w.Noff) / max(sum(d.V) / L, 1e-300) * 0, zeros(LD_NR, LD_NR))
            # mutation: recursion Σ_a Σ_B m n_a/(2N); fix adds Σ_a V_a n_as m/(2N) on the diagonal
            mrt = 0.0; mrr = zeros(LD_NR, LD_NR)
            for a in 1:LD_NB
                G.nb[a] == 0 && continue
                t, rmat = cs_summ_mat(class_cov(G, d, a) .* (mloc * G.nb[a] / (2w.Noff))); mrt += t; mrr .+= rmat
            end
            prevO = (vals[:O], Dict(:mut_recur => (mrt, mrr), :mut_fix => (mrt + md_tot / (2w.Noff), mrr .+ md_rr ./ (2w.Noff))))
            prevgen = it
            push!(rec[:gen], it)
            for k in CS_KEYS
                push!(rec[k], get(vals, k, (NaN, fill(NaN, LD_NR, LD_NR))))
            end
            QTLSim.rescale!(w); QTLSim.mutate!(w, rng); QTLSim.swapbuffers!(w)
            continue
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        QTLSim.rescale!(w); QTLSim.mutate!(w, rng); QTLSim.swapbuffers!(w)
    end
    (; seed, extinct_gen, rec, vgser, L, u, PopSize = q.PopSize)
end

cs_file(id, r, set) = joinpath(LL_RAW, replace(id, "/" => "__", ".txt" => "") * "@class" * set * "#$r.jls")

function cs_kernel(set, R)
    mkpath(LL_RAW)
    rows = ll_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(cs_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) class-space runs for set $set, R = $R, $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time(); cpu = Threads.Atomic{Float64}(0.0)
    Threads.@threads :greedy for (id, j, k) in jobs
        t1 = time()
        tr = cs_run(readp(id), ll_seed(set, j, k))
        serialize(cs_file(id, k, set), tr)
        Threads.atomic_add!(cpu, time() - t1)
        d = Threads.atomic_add!(done, 1) + 1
        @printf("  %d/%d %s #%d  %.0f s\n", d, length(jobs), id, k, time() - t0)
    end
    @printf("done: %d runs, %.0f s wall, %.0f thread-s\n", length(jobs), time() - t0, cpu[])
end

const CS_SHOW = [:sel_meas, :sel_recur, :surv_meas, :surv_recur, :rec_meas, :rec_recur, :rec_exact_rho, :rec_fixb, :rec_fixp, :rec_fix,
                 :off_meas, :off_recur, :mut_meas, :mut_recur, :mut_fix]

function cs_stages(set, R; io = stdout, wins = [(-60, -1), (0, 49), (50, 99)], tight = 1:3)
    rows = ll_rows(set)
    println(io, "\n== class-space stage accounting, set $set, R = $R: change per generation / Σ_l Y_l² at A; total, and the tightly linked block (split classes $(first(tight))..$(last(tight)), ρ < $(LD_REDGES[last(tight)])) ==")
    agg = Dict{Tuple{Symbol,Symbol,Tuple{Int,Int}},Vector{Float64}}()
    for rw in rows
        runs = [deserialize(cs_file(rw.id, k, set)) for k in 1:R if isfile(cs_file(rw.id, k, set))]
        isempty(runs) && continue
        for wdw in wins
            sel(tr) = [k for k in eachindex(tr.rec[:gen]) if wdw[1] <= tr.rec[:gen][k] - MOVE <= wdw[2]]
            den = mean(reduce(vcat, [[tr.rec[:A][k][1] for k in sel(tr)] for tr in runs]))
            dent = mean(reduce(vcat, [[sum(tr.rec[:A][k][2][tight, tight]) for k in sel(tr)] for tr in runs]))
            for c in CS_SHOW
                v = filter(isfinite, reduce(vcat, [[tr.rec[c][k][1] for k in sel(tr)] for tr in runs]))
                vt = filter(isfinite, reduce(vcat, [[sum(tr.rec[c][k][2][tight, tight]) for k in sel(tr)] for tr in runs]))
                push!(get!(agg, (c, :tot, wdw), Float64[]), mean(v) / den)
                push!(get!(agg, (c, :tight, wdw), Float64[]), mean(vt) / dent)
            end
            push!(get!(agg, (:share_tight, :tot, wdw), Float64[]), dent / den)
        end
    end
    for wdw in wins
        @printf(io, "  window %d..%d (%d inputs; tight block's share of ΣY² median %.2f)\n", wdw..., length(agg[(:sel_meas, :tot, wdw)]), medf(agg[(:share_tight, :tot, wdw)]))
        for c in CS_SHOW
            @printf(io, "    %-13s total %+.4f (median %+.4f)   tight %+.4f (median %+.4f)\n", c, mean(agg[(c, :tot, wdw)]), medf(agg[(c, :tot, wdw)]),
                    mean(agg[(c, :tight, wdw)]), medf(agg[(c, :tight, wdw)]))
        end
        for (nm, cols) in (("recursion", [:sel_recur, :surv_recur, :rec_recur, :off_recur, :mut_recur]), ("measured", [:sel_meas, :surv_meas, :rec_meas, :off_meas, :mut_meas]),
                           ("fixed rec+mut", [:sel_recur, :surv_recur, :rec_fix, :off_recur, :mut_fix]))
            for part in (:tot, :tight)
                tt = [sum(agg[(c, part, wdw)][i] for c in cols) for i in eachindex(agg[(:sel_meas, part, wdw)])]
                @printf(io, "    sum %-13s %-5s %+.4f (median %+.4f)\n", nm, part, mean(filter(isfinite, tt)), medf(tt))
            end
        end
    end
    agg
end

# ---------------------------------------------------------------------------------------
# the full measured class moments (P_st per group) before the move, to test the recursion's
# fixed point against the measured state entry by entry

function pm_run(p, seed; gens = MOVE-60:5:MOVE-1)
    BLAS.set_num_threads(1)
    q = kernel_params(p)
    rng = MersenneTwister(seed)
    w = equilibrium_world(q, rng; burnin = KBURN)
    L = q.NoQTL; u = q.RecombPoissonMean
    B = STBuf(L, q.PopSize, u); cs = class_seg(L, u); C = CSBuf(cs, L)
    G = ld_groups(L, u); ns = size(G.n, 2)
    binl = [ld_bin(l - 1, L) for l in 1:L]
    Psum = [zeros(ns, ns) for _ in 1:LD_NB]; msum = [zeros(ns) for _ in 1:LD_NB]; nm = 0
    for it in w.gen:last(gens)
        w.gen = it
        QTLSim.develop!(w); QTLSim.select!(w, rng)
        if it in gens
            N = w.N
            load_pool!(B, w.pool, w.chr, N)
            ma = moments(B, 1:N, ones(N))
            y = class_y(B, C, 1:N, ones(N), ma).y
            for l in 1:L
                a = binl[l]; r = cs.coff[l]:cs.coff[l+1]-1
                cl = cs.cls[r]; yl = y[r]
                msum[a][cl] .+= yl
                Psum[a][cl, cl] .+= yl * yl'
            end
            nm += 1
        end
        QTLSim.reproduce!(w, rng); QTLSim.rescale!(w); QTLSim.mutate!(w, rng); QTLSim.swapbuffers!(w)
    end
    (; P = Psum ./ nm, m = msum ./ nm, L, u)
end
pm_file(id, r, set) = joinpath(LL_RAW, replace(id, "/" => "__", ".txt" => "") * "@pm" * set * "#$r.jls")
function pm_kernel(set, R)
    rows = ll_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(pm_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) class-moment runs for set $set")
    t0 = time()
    Threads.@threads :greedy for (id, j, k) in jobs
        serialize(pm_file(id, k, set), pm_run(readp(id), ll_seed(set, j, k)))
    end
    @printf("done, %.0f s wall\n", time() - t0)
end

# ---------------------------------------------------------------------------------------
# the pre-specified verdicts (validate/linked_ld_criteria.md)

function ld_confirm(set, R; io = stdout, variants = [:B0, :Bs, :Bp0, :Bps])
    rows = ll_rows(set)
    println(io, "\n== linked LD, set $set, R = $R: per input (Σ E[yY] predicted/measured, pre / 50..99) ==")
    bufs = [IOBuffer() for _ in rows]
    resv = Vector{Any}(undef, length(rows))
    Threads.@threads :greedy for i in eachindex(rows)
        resv[i] = fix_input(rows[i], R, set; variants, io = bufs[i], withA = true)
    end
    foreach(b -> print(io, String(take!(b))), bufs)
    allr = [(rows[i], resv[i]) for i in eachindex(rows) if resv[i] !== nothing]
    elig = filter(x -> isfinite(x[2][:meas_loss]) && x[2][:meas_loss] <= -0.05, allr)
    println(io, "\n== summary: $(length(allr)) inputs, $(length(elig)) eligible (measured cumulative 2βΣC over 0..99 ≤ −0.05 of the genic variance) ==")
    ratio(x, num, den, w) = x[2][num][w] / x[2][den][w]
    stat(num, den, w; xs = allr) = (v = [lr(x[2][num][w], x[2][den][w]) for x in xs]; (exp(medf(v)), medf(v), mabs(v), count(isfinite, v)))
    for (lab, num, den) in (("Σ E[yY]", :Bs_yY, :yY), ("Σ E[yY] B0", :B0_yY, :yY), ("Σ E[yY] B′*", :Bps_yY, :yY), ("Σ E[yY] B′", :Bp0_yY, :yY),
                            ("Σ H B′*", :Bps_H, :H), ("Σ H B′", :Bp0_H, :H), ("Σ K B′*", :Bps_K, :K), ("Σ K B′", :Bp0_K, :K), ("Σ K B*", :Bs_K, :K), ("Σ K B0", :B0_K, :K))
        @printf(io, "  %-12s", lab)
        for w in WINS
            r, l, a, n = stat(num, den, w)
            @printf(io, "  %d..%d %.3f (|lr| %.3f, n %d)", w..., r, a, n)
        end
        println(io)
    end
    println(io, "  closed ΣC, eligible: log(predicted/measured) median / median |.|:")
    for v in variants
        @printf(io, "    %-4s", v)
        for w in WINS[3:5]
            u = [lr(x[2][Symbol(v, :_C)][w], x[2][:C][w]) for x in elig]
            @printf(io, "  %d..%d %+.3f/%.3f (n %d)", w..., medf(u), mabs(u), count(isfinite, u))
        end
        u = [lr(x[2][Symbol(v, :_loss)], x[2][:meas_losssym]) for x in elig]
        @printf(io, "  loss 0..99 vs sign-symmetric %+.3f/%.3f (n %d)\n", medf(u), mabs(u), count(isfinite, u))
    end
    println(io, "  P1 (level A), eligible: log(predicted ΣC / measured ΣC) median / median |.|:")
    P1 = Dict{Tuple{Symbol,Tuple{Int,Int}},Tuple{Float64,Float64}}()
    for md in AMODES, w in WINS[3:5]
        u = [lr(x[2][Symbol(:A_, md)][w], x[2][:C][w]) for x in elig]
        P1[(md, w)] = (medf(u), mabs(u))
        @printf(io, "    %-9s %d..%d %+.3f / %.3f (n %d)\n", md, w..., medf(u), mabs(u), count(isfinite, u))
    end
    # verdicts
    pre, mid = (-60, -1), (50, 99)
    vrd(k, n) = k == n ? "PASS" : (k >= 1 && n == 2 ? "PARTIAL" : (k == n - 1 && n == 3 ? "PARTIAL" : "FAIL"))
    println(io, "\n== pre-specified verdicts (validate/linked_ld_criteria.md) ==")
    oka = [begin r, _, a, _ = stat(:Bs_yY, :yY, w); 0.80 <= r <= 1.25 && a <= 0.30 end for w in (pre, mid)]
    @printf(io, "(a) Σ E[yY] level B*: %s -> %s\n", join([(r = stat(:Bs_yY, :yY, w); @sprintf("%d..%d ratio %.3f |lr| %.3f", w..., r[1], r[3])) for w in (pre, mid)], "; "),
            sum(oka) == 2 ? "PASS" : sum(oka) == 1 ? "PARTIAL" : "FAIL")
    okb = [begin r, _, a, _ = stat(:Bps_H, :H, w); 0.67 <= r <= 1.5 && a <= 0.40 end for w in (pre, mid)]
    @printf(io, "(b) Σ H level B′*: %s -> %s\n", join([(r = stat(:Bps_H, :H, w); @sprintf("%d..%d ratio %.3f |lr| %.3f", w..., r[1], r[3])) for w in (pre, mid)], "; "),
            sum(okb) == 2 ? "PASS" : sum(okb) == 1 ? "PARTIAL" : "FAIL")
    if length(elig) < 6
        println(io, "(c) INCONCLUSIVE (fewer than 6 eligible)"); println(io, "P1 INCONCLUSIVE")
    else
        okc = Bool[]
        for w in WINS[3:4]
            u = [lr(x[2][:Bps_C][w], x[2][:C][w]) for x in elig]
            push!(okc, abs(medf(u)) <= 0.3 && mabs(u) <= 0.4)
        end
        u = [lr(x[2][:Bps_loss], x[2][:meas_losssym]) for x in elig]
        push!(okc, abs(medf(u)) <= 0.3 && mabs(u) <= 0.4)
        @printf(io, "(c) closed level B′*: 25..49 %s, 50..99 %s, cumulative loss %s -> %s\n", okc[1] ? "ok" : "fails", okc[2] ? "ok" : "fails", okc[3] ? "ok" : "fails",
                sum(okc) == 3 ? "PASS" : sum(okc) == 2 ? "PARTIAL" : "FAIL")
        acc = [abs(P1[(:derived, w)][1]) <= 0.2 && P1[(:derived, w)][2] <= 0.3 for w in WINS[3:5]]
        comp = all(P1[(:derived, w)][2] < P1[(md, w)][2] for w in WINS[4:5] for md in (:nosel, :spectrum))
        @printf(io, "P1 (decay rate): accuracy %s, comparator condition %s -> %s\n", join([a ? "ok" : "fails" for a in acc], "/"), comp ? "holds" : "fails",
                comp && all(acc) ? "PASS" : comp && sum(acc) == 2 ? "PARTIAL" : "FAIL")
    end
    okd = [begin r = stat(:Bps_K, :K, w)[1]; 0.75 <= r <= 1.33 end for w in (pre, mid)]
    @printf(io, "(d) source ΣK level B′*: %s -> %s\n", join([@sprintf("%d..%d %.3f", w..., stat(:Bps_K, :K, w)[1]) for w in (pre, mid)], "; "), all(okd) ? "PASS" : "FAIL")
    allr
end

# ---------------------------------------------------------------------------------------
# H′: the fresh confirmation set (criteria validate/linked_ld_criteria.md)

const HSEED = 70261115000                         # H′ seeds: HSEED + 100j + r
const HPER = 3                                    # inputs per stratum

"""
H′: inputs of the published sweep never used by the spectrum-diagnosis family of studies (not in
validate/spectrum_diagnosis_inputs.csv, so in none of D, F, G, C, Dctl, Cctl, E6′), screened on
parameters only as linked_source's G: k²/V_G ≥ 0.01 (V_G = `house_of_gauss_vg`) and
k/k_c ≤ 1.4 (k_c = `drift_critical_rate`). Per stratum (NoQTL × u) the files are taken in the
order of a MersenneTwister(20261117) permutation until `HPER` pass. j = 1000 + position, so the
seeds HSEED + 100j + r are distinct from every earlier set's.
"""
function h_rows(; verbose = false)
    used = Set(r.id for r in read_inputs_csv())
    rng = MersenneTwister(20261117)
    out = NamedTuple[]
    j = 1000
    for s in sort(filter(x -> startswith(x, "QT"), readdir(SIMDIR)))
        files = sort(filter(f -> endswith(f, ".txt"), readdir(joinpath(SIMDIR, s))))
        npass = 0; ntried = 0
        for f in files[randperm(rng, length(files))]
            id = s * "/" * f
            id in used && continue
            p = readp(id)
            vg = house_of_gauss_vg(p); kc = drift_critical_rate(p).k
            idx = p.MoveFitnessSD^2 / vg; ratio = p.MoveFitnessSD / kc
            ntried += 1
            if idx >= 0.01 && ratio <= 1.4
                j += 1; npass += 1
                push!(out, (j = j, id = id, set = "H", f = NaN, P = NaN, idx = idx, ratio = ratio))
                npass == HPER && break
            end
        end
        verbose && println("  $s: $npass pass of $ntried screened")
    end
    out
end

"""
Association runs (linked_source's `ls_run`, no decomposition pass, as its G) for a new set, into
validate/linked_ld_raw/ (`as_file`).
"""
function as_kernel(set, R)
    mkpath(LL_RAW)
    rows = ll_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(as_file(r.id, k, set))]
    sort!(jobs; by = j -> -readp(j[1]).NoQTL)
    println("$(length(jobs)) association runs for set $set, R = $R, $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time(); cpu = Threads.Atomic{Float64}(0.0)
    Threads.@threads :greedy for (id, j, k) in jobs
        t1 = time()
        tr = ls_run(readp(id), ll_seed(set, j, k); decompose = false)
        serialize(as_file(id, k, set), tr)
        Threads.atomic_add!(cpu, time() - t1)
        d = Threads.atomic_add!(done, 1) + 1
        @printf("  %d/%d %s #%d  %.0f s\n", d, length(jobs), id, k, time() - t0)
    end
    @printf("done: %d runs, %.0f s wall, %.0f thread-s\n", length(jobs), time() - t0, cpu[])
end

const _HROWS = Ref{Any}(nothing)
hrows() = (_HROWS[] === nothing && (_HROWS[] = h_rows()); _HROWS[])

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "check" : ARGS[1]
    if mode == "check"
        st_check()
    elseif mode == "stage"
        st_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "stages"
        st_stages(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "class"
        cs_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "classes"
        cs_stages(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "pm"
        pm_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "assoc"
        as_kernel(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "hrows"
        for r in h_rows(; verbose = true); @printf("%5d %-36s k²/V_G %.4f k/k_c %.3f\n", r.j, r.id, r.idx, r.ratio); end
    elseif mode == "confirm"
        ld_confirm(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "fixdev"
        fix_dev(ARGS[2], parse(Int, ARGS[3]))
    end
end
