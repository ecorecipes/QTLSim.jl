# Symbolic and exact checks of the haplotype-block lineage model (src/block_lineages.jl:
# `QTLSim.qtsim_switch_law`, `QTLSim.block_map`, `run_blocks`, `QTLSim.block_stats`), with
# numerical spot checks against its code and against the kernel's `gamete!`.
#
# Usage (from the package root):
#     julia --project=validate/symbolic validate/symbolic/check_blocks.jl
# Writes validate/symbolic/check_blocks.log; exits 1 if any check fails.
#
# (a) QTSim's lowest-point rule: K distinct points uniform among the M = L − 1 gaps; for
#     1 ≤ K ≤ M, P(min ≥ c | K) = C(M−c, K)/C(M, K), hence P(min = c | K) = C(M−c−1, K−1)/C(M, K);
#     the code's ratio recurrence r_{c+1}/r_c = (M−c−K)/(M−c−1), r_0 = K/M (symbolic in M, c for
#     each K); exact enumeration of all K-subsets for small M; the hockey-stick sum
#     Σ_c C(M−c−1, K−1) = C(M, K) (the law sums to one);
# (b) `qtsim_switch_law` against an independent BigFloat sum, and against 4·10⁵ draws of the
#     kernel's `gamete!` (the switch point read off a gamete of two marked strands);
# (c) `block_map`: P(no switch) + Σ_j P(switch after block j) = 1; with one block per locus the
#     probability that two loci are split is the kernel's exactly; with fewer blocks it differs
#     by at most the largest within-block share of the map plus twice the largest single point (half of each end
#     block, and the discreteness of the nearer-boundary rule); exactly min(B, L) blocks; `:all` moves each point the same way;
# (d) mutation: Poisson(m·λ) mutations placed on uniform chromosomes have the law of independent
#     Poisson(λ) per chromosome (pgf identity, symbolic), and a uniform locus falls in block b
#     with probability n_b/L, so the expected variance input per individual is unchanged
#     (2·λ·σ²_mut = m·V_E, Lean `mutationalVarianceInput_eq`);
# (e) sufficiency: with switches at block boundaries, the block sums of a gamete built locus by
#     locus equal the gamete built from block sums;
# (f) `block_stats`: with one block per locus the block-genic variance is the genic variance
#     2Σ_l Var(x_l), and the association equals its brute-force definition;
# (g) the Gaussian segregation law: for Gaussian halves in proportion f : 1 − f, E[left | y] = f y,
#     Var(left | y) = f(1 − f)V, the recombinant keeps the block's variance and its covariances
#     with both parents (symbolic, and Monte Carlo against exact recombination); block_map's f.

include("common.jl")

const LOG = joinpath(@__DIR__, "check_blocks.log")
const rng = Xoshiro(20261116)

function main()
    t0 = time()
    say("QTLSim validate/symbolic/check_blocks.jl (Symbolics ", pkgversion(Symbolics), ")")

    # ===================================================================================
    say("\n(a) the lowest of K distinct points among M gaps")
    @variables M c
    binp(a, k) = k == 0 ? 1 : prod(a - i for i in 0:k-1) / factorial(k)     # C(a, k), polynomial in a
    for K in 1:6
        # P(min = c | K) C(M, K) = C(M−c, K) − C(M−c−1, K) = C(M−c−1, K−1)  (Pascal)
        check("P(min ≥ c) − P(min ≥ c+1) = C(M−c−1, K−1)/C(M, K), K = $K",
              iszero_sym(binp(M - c, K) - binp(M - c - 1, K) - binp(M - c - 1, K - 1)))
        # ratio recurrence of the code: C(M−c−2, K−1)(M−c−1) = C(M−c−1, K−1)(M−c−K)
        check("r_{c+1}/r_c = (M−c−K)/(M−c−1), K = $K",
              iszero_sym(binp(M - c - 2, K - 1) * (M - c - 1) - binp(M - c - 1, K - 1) * (M - c - K)))
        check("r_0 = C(M−1, K−1)/C(M, K) = K/M, K = $K", iszero_sym(binp(M - 1, K - 1) * M - binp(M, K) * K))
    end
    # exact enumeration
    ok = true; okhs = true
    for Mv in 2:9, K in 1:Mv
        cnt = zeros(Int, Mv)
        for s in 0:(1<<Mv)-1
            count_ones(s) == K || continue
            cnt[trailing_zeros(s)+1] += 1
        end
        tot = binomial(Mv, K)
        for cv in 0:Mv-1
            f = cv <= Mv - K ? binomial(Mv - cv - 1, K - 1) : 0
            ok &= cnt[cv+1] == f
        end
        okhs &= sum(binomial(Mv - cv - 1, K - 1) for cv in 0:Mv-K) == tot
    end
    check("enumeration of all K-subsets, M = 2..9: counts of the minimum = C(M−c−1, K−1)", ok)
    check("hockey stick Σ_c C(M−c−1, K−1) = C(M, K), M = 2..9", okhs)
    okhs = all(sum(big(binomial(big(Mv - cv - 1), K - 1)) for cv in 0:Mv-K) == binomial(big(Mv), K) for Mv in (99, 499, 4999) for K in (1, 2, 5, 17))
    check("hockey stick at M = 99, 499, 4999, K = 1, 2, 5, 17", okhs)

    # ===================================================================================
    say("\n(b) qtsim_switch_law against a BigFloat sum and the kernel's gamete!")
    for (L, u) in ((2, 3.0), (5, 3.0), (20, 5.0), (100, 3.0), (1000, 5.0))
        law = QTLSim.qtsim_switch_law(L, u)
        Mv = L - 1
        ref = zeros(BigFloat, Mv)
        pn = exp(-big(u))
        for K in 1:Mv
            pk = exp(-big(u)) * big(u)^K / factorial(big(K))
            pk < big(1e-40) && break
            for cv in 0:Mv-K
                ref[cv+1] += pk * binomial(big(Mv - cv - 1), K - 1) / binomial(big(Mv), K)
            end
        end
        pn += 1 - exp(-big(u)) * sum(big(u)^K / factorial(big(K)) for K in 0:min(Mv, 150))
        err = maximum(abs.(Float64.(ref) .- law.P); init = 0.0)
        check(@sprintf("L = %4d, u = %.0f: max |P(c) − BigFloat| = %.1e, pnone %.6f vs %.6f", L, u, err, law.pnone, Float64(pn)),
              err < 1e-13 && abs(law.pnone - pn) < 1e-13)
        check(@sprintf("L = %4d: pnone + Σ P(c) = 1 (%.1e)", L, abs(law.pnone + sum(law.P) - 1)), abs(law.pnone + sum(law.P) - 1) < 1e-12)
    end
    for (L, u) in ((12, 3.0), (30, 5.0))
        law = QTLSim.qtsim_switch_law(L, u)
        A = zeros(L); B = ones(L); d = zeros(L)
        n = 400_000
        cnt = zeros(Int, L - 1); none = 0
        for _ in 1:n
            QTLSim.gamete!(d, A, B, rng, u, L)
            k = findfirst(i -> d[i] != d[1], 2:L)
            k === nothing ? (none += 1) : (cnt[k] += 1)     # switch after 0-based locus k−1
        end
        z = [(cnt[i] - n * law.P[i]) / sqrt(n * law.P[i] * (1 - law.P[i]) + 1e-12) for i in 1:L-1]
        zn = (none - n * law.pnone) / sqrt(n * law.pnone * (1 - law.pnone))
        check(@sprintf("L = %d, u = %.0f: %d gamete! draws, max |z| over the switch points %.2f, no switch z %.2f", L, u, n,
                       maximum(abs.(z)), zn), maximum(abs.(z)) < 4.5 && abs(zn) < 4.5)
    end

    # ===================================================================================
    say("\n(c) block_map")
    for (L, u, cm) in ((100, 3.0, :qtsim), (500, 5.0, :qtsim), (1000, 3.0, :qtsim), (200, 3.0, :all))
        p = Params(; NoQTL = L, RecombPoissonMean = u, CrossoverModel = cm)
        law = cm === :qtsim ? QTLSim.qtsim_switch_law(L, u).P : fill(u / (L - 1), L - 1)
        for nb in (10, 40, L)
            bm = QTLSim.block_map(p, nb)
            cm === :qtsim && check(@sprintf("L %d, B %d: pnone + Σ bmass = 1 (%.1e)", L, nb, abs(bm.pnone + sum(bm.bmass) - 1)),
                                   abs(bm.pnone + sum(bm.bmass) - 1) < 1e-12)
            check("L $L, B $nb: exactly min(B, L) blocks covering the loci in order", bm.nb == min(nb, L) && sum(bm.nloci) == L && issorted(bm.hi) && all(>(0), bm.nloci))
            # split probability of loci i < j: kernel Σ_{c=i-1}^{j-2} P(c) (gaps between them,
            # 0-based gap c lies after 1-based locus c+1); block: Σ of the masses of the
            # boundaries between their blocks (each boundary's mass is that of the points
            # moved to it)
            bmass = zeros(bm.nb + 1)                       # boundaries 0:nb, index +1
            for cv in 0:L-2
                bmass[bm.snap[cv+1]+1] += law[cv+1]
            end
            cumb = cumsum(bmass)                           # cumb[j+1] = Σ_{boundaries ≤ j}
            cumk = [0.0; cumsum(law)]
            inner = [sum(law[cv+1] for cv in 0:L-2 if bm.blockof[cv+1] == b && (b == bm.nb || cv != bm.hi[b] - 1); init = 0.0)
                     for b in 1:bm.nb]
            maxerr = 0.0
            for _ in 1:4000
                i, j = minmax(rand(rng, 1:L), rand(rng, 1:L))
                i == j && continue
                ks = cumk[j] - cumk[i]                     # gaps i-1 … j-2
                bi, bj = bm.blockof[i], bm.blockof[j]
                bs = bi == bj ? 0.0 : cumb[bj] - cumb[bi]  # boundaries bi … bj-1 (index +1)
                maxerr = max(maxerr, abs(ks - bs))
            end
            bound = 1.0000001 * (maximum(inner) + 2maximum(law))   # half of each end block, and one point of discreteness each
            if nb >= L
                check(@sprintf("L %d, %s, one block per locus: split probabilities exact (max err %.1e)", L, cm, maxerr), maxerr < 1e-12)
            else
                check(@sprintf("L %d, %s, B %d: max |split error| %.4f ≤ largest within-block share + 2 max point %.4f", L, cm, nb, maxerr, bound),
                      maxerr <= bound)
            end
        end
    end
    bm = QTLSim.block_map(Params(; NoQTL = 1000, RecombPoissonMean = 3.0), 40)
    inner = bm.bmass
    check(@sprintf(":map partition, L 1000, u 3, B 40: boundary masses %.4f–%.4f (≈ (1 − pnone)/B = %.4f)", extrema(inner)...,
                   (1 - bm.pnone) / 39), maximum(inner) / minimum(inner) < 1.5)

    # ===================================================================================
    say("\n(d) mutation")
    @variables λ s1 s2 s3
    sv = [s1, s2, s3]; m = 3
    # pgf of the counts: total Poisson(mλ), each mutation on chromosome h w.p. 1/m:
    # E Π s_h^{n_h} = exp(mλ(Σ s_h/m − 1)); independent Poisson(λ): Π exp(λ(s_h − 1))
    check("log pgf: mλ(Σ s_h/m − 1) = Σ λ(s_h − 1)", iszero_sym(m * λ * (sum(sv) / m - 1) - sum(λ * (s - 1) for s in sv)))
    @variables L_ μ_ mm VE P_
    σ2 = mm / (P_ * L_ * μ_) * VE
    check("expected input per individual: Ploidy·(Lμ)·σ²_mut = m·V_E (unchanged by blocks)", iszero_sym(P_ * L_ * μ_ * σ2 - mm * VE))
    p = Params(; NoQTL = 500, RecombPoissonMean = 3.0)
    bm = QTLSim.block_map(p, 40)
    hits = zeros(Int, bm.nb)
    for _ in 1:200_000
        hits[bm.blockof[rand(rng, 1:500)]] += 1
    end
    z = maximum(abs.((hits .- 200_000 .* bm.nloci ./ 500) ./ sqrt.(200_000 .* bm.nloci ./ 500)))
    check(@sprintf("mutations fall in block b with probability n_b/L (max |z| %.2f)", z), z < 4.5)

    # ===================================================================================
    say("\n(e) sufficiency of block sums")
    ok = true
    for _ in 1:200
        L = rand(rng, 20:60); nbv = rand(rng, 3:10)
        bm = QTLSim.block_map(Params(; NoQTL = L, RecombPoissonMean = 3.0), nbv)
        A = randn(rng, L); B = randn(rng, L)
        j = rand(rng, 1:bm.nb-1)                   # switch after block j: after locus hi[j]
        g = [A[1:bm.hi[j]]; B[bm.hi[j]+1:L]]
        yA = [sum(A[bm.blockof .== b]) for b in 1:bm.nb]; yB = [sum(B[bm.blockof .== b]) for b in 1:bm.nb]
        yg = [sum(g[bm.blockof .== b]) for b in 1:bm.nb]
        ok &= isapprox(yg, [yA[1:j]; yB[j+1:end]]; atol = 1e-12)
        # a mutation of effect x at locus l: the block sum gains x
        l = rand(rng, 1:L); x = randn(rng); A2 = copy(A); A2[l] += x
        ok &= isapprox(sum(A2[bm.blockof .== bm.blockof[l]]), yA[bm.blockof[l]] + x; atol = 1e-12)
    end
    check("block sums of a gamete switching at a boundary, and of a mutated chromosome (200 random cases)", ok)

    # ===================================================================================
    say("\n(f) block_stats")
    ok = true
    for _ in 1:20
        nbv = rand(rng, 2:8); Mv = rand(rng, 10:40)
        Y = randn(rng, nbv, Mv) .^ 2
        s = QTLSim.block_stats(Y, nbv, Mv)
        ybar = vec(sum(Y; dims = 2)) / Mv
        G = vec(sum(Y; dims = 1)); Gb = sum(G) / Mv
        bg = 2 * sum(sum((Y[b, :] .- ybar[b]) .^ 2) / Mv for b in 1:nbv)
        as = sum(begin
                     d = Y[b, :] .- ybar[b]; r = G .- d
                     sum((r .- sum(r) / Mv) .* (d .^ 2 .- sum(d .^ 2) / Mv)) / Mv
                 end for b in 1:nbv)
        ok &= isapprox(s.blockgenic, bg; rtol = 1e-12) && isapprox(s.assoc, as; rtol = 1e-10, atol = 1e-14)
    end
    check("block_stats: block-genic and association against their definitions (20 random populations)", ok)

    # ===================================================================================
    say("\n(g) Gaussian segregation within a block (within = :gaussian)")
    @variables f V
    # halves L ~ N(0, fV), R ~ N(0, (1−f)V) independent, y = L + R: regression of L on y
    check("Cov(L, y)/Var(y) = f", iszero_sym(f * V / V - f))
    check("Var(L | y) = fV − (fV)²/V = f(1 − f)V", iszero_sym(f * V - (f * V)^2 / V - f * (1 - f) * V))
    # recombinant L_A + R_B given y_A, y_B: mean f y_A + (1−f) y_B, variance f(1−f)V + (1−f)fV
    check("Var(L_A | y_A) + Var(R_B | y_B) = 2f(1 − f)V", iszero_sym(f * (1 - f) * V + (1 - f) * f * V - 2f * (1 - f) * V))
    # the block's variance is kept: Var(f y_A + (1−f) y_B) + 2f(1−f)V = V when Cov(y_A, y_B) = 0,
    # (random mating; with homologue covariance c the law would add 2f(1−f)c, exact recombination would not)
    check("Var of the recombinant block = V (homologues uncorrelated)",
          iszero_sym(f^2 * V + (1 - f)^2 * V + 2f * (1 - f) * V - V))
    # Monte Carlo against exact recombination of Gaussian halves
    n = 400_000; fv = 0.3; Vv = 2.0
    LA = sqrt(fv * Vv) .* randn(rng, n); RA = sqrt((1 - fv) * Vv) .* randn(rng, n)
    LB = sqrt(fv * Vv) .* randn(rng, n); RB = sqrt((1 - fv) * Vv) .* randn(rng, n)
    yA = LA .+ RA; yB = LB .+ RB; ex = LA .+ RB
    gs = fv .* yA .+ (1 - fv) .* yB .+ sqrt(2fv * (1 - fv) * Vv) .* randn(rng, n)
    cov_(a, b) = sum((a .- sum(a) / n) .* (b .- sum(b) / n)) / n
    ok = abs(cov_(gs, gs) / cov_(ex, ex) - 1) < 0.01 && abs(cov_(gs, yA) / cov_(ex, yA) - 1) < 0.02 && abs(cov_(gs, yB) / cov_(ex, yB) - 1) < 0.02
    check(@sprintf("Monte Carlo (f = 0.3): variance and covariances with both parents match exact recombination (%.4f, %.4f, %.4f)",
                   cov_(gs, gs) / cov_(ex, ex), cov_(gs, yA) / cov_(ex, yA), cov_(gs, yB) / cov_(ex, yB)), ok)
    p = Params(; NoQTL = 300, RecombPoissonMean = 3.0)
    bm = QTLSim.block_map(p, 20; within = :gaussian)
    law = QTLSim.qtsim_switch_law(300, 3.0)
    okf = all(c -> (bm.gf[c+1] == 1.0) == (c + 1 in bm.hi), 0:298) &&
          all(c -> bm.gf[c+1] ≈ (c + 1 - (bm.blockof[c+1] == 1 ? 0 : bm.hi[bm.blockof[c+1]-1])) / bm.nloci[bm.blockof[c+1]], 0:298)
    check("block_map(within = :gaussian): f = loci of the block at or before the switch / block size; f = 1 exactly on boundaries", okf)
    check("block_map(within = :gaussian): the switch law is the kernel's (gcum)", bm.gcum ≈ cumsum(law.P) ./ (1 - law.pnone) && bm.pnk == law.pnone)

    finish(LOG, t0)
end

runmain(main, @__FILE__)
