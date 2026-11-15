# Diagnostic study (exploratory): why does the kernel's genic variance fall at once when the
# optimum starts to move (published protocol: continuous move, N = 1000, QTSim's single
# crossover), while the `:spectrum` closure's, even driven by the kernel's own lag, S and
# survivors, stays flat for ~100 generations (validate/spectrum_diagnosis.md, 2(b))?
#
# Method: an exact accounting of the change of the genic variance, G = 2 Σ_l v_l (v_l the
# variance of the allele values of the 2N genes at locus l, population denominators), over
# the stages of one generation of the kernel:
#
#   adults A  →  π-weighted adults  →  W-weighted adults  →  survivors S  →  offspring O  →  mutated M (= next A)
#
# π_i = E[W | g_i] = c ω/√(ω² + V_E) exp(−(g_i − θ)²/(2(ω² + V_E))) is the survival
# probability given the genotype (deterministic selection on genotypes); W_i the realised
# fitness given the phenotype (adds the environmental noise); S the realised survivors (adds
# the Bernoulli noise); O the offspring before mutation (sampling of parents and gametes).
# The deterministic selection term G_π − G_A is split per locus by Price's equation,
#   v_l^π − v_l = Cov(π/π̄, (x − x̄_l)²) − (Δx̄_l)²,
# the first part ("within") into β κ3_l (the locus's third central moment times the
# regression β of relative fitness on the genotype) and the rest (stabilising), and the
# response Δx̄_l into its own part o_l = β Var(y_l)/2 (y_l the individual's two genes at l),
# its linked part k_l = β (Cov(g, y_l) − Var(y_l))/2 and a nonlinear rest; k_l = φ o_l + ε_l
# by a regression across loci, so (1 + φ) is the systematic (Bulmer) LD's scaling of the
# response and ε_l the random association with the rest of the genome (hitchhiking, draft,
# Robertson's heritable fitness background).
#
# Instrumentation reads the World only (no random draws): `check` verifies that a run is
# bitwise the run of `simulate` with the same seed.
#
# Run from the package root (vignettes environment):
#   julia -t 8 --project=vignettes validate/genic_decline.jl check
#   julia -t 8 --project=vignettes validate/genic_decline.jl kernel SET R [VARIANT]   # instrumented runs (resumable)
#   julia -t 8 --project=vignettes validate/genic_decline.jl analyse SET R [VARIANT]  # the accounting, mean over runs alive
#   julia -t 8 --project=vignettes validate/genic_decline.jl driven SET R             # the spectrum's partition, driven by each run
# SET as in validate/spectrum_diagnosis.jl (D, Dctl, C, Cctl joined by +). VARIANT: "" (the
# input), many (CrossoverModel = :all, Poisson(20)), fec (Selection = :fecundity), bigN
# (PopSize × 4), sites10 (NoQTL × 10, MutationRatePerGamete / 10: the same V_m and law of
# effects, closer to infinite sites).
# Raw per-run accounting goes to validate/genic_decline_raw/ (git-ignored).

include(joinpath(@__DIR__, "spectrum_diagnosis.jl"))   # readp, kernel_params, MOVE, KBURN, KSEED, inputs, C0

const GD_RAW = joinpath(@__DIR__, "genic_decline_raw")
const GD_FROM, GD_TO = MOVE - 200, MOVE + 300          # instrumented generations
const GD_CENSUS = 25                                   # allele census every GD_CENSUS generations
const NBIN = 10                                        # position bins

# counterfactual inputs (CF): six discovery inputs with a clear early decline, 500–5000 loci, both u
const CF_IDS = ["QT0500-Possion3.0/Input-04943.txt", "QT1000-Possion3.0/Input-00854.txt", "QT1000-Possion5.0/Input-14110.txt",
                "QT1000-Possion5.0/Input-07594.txt", "QT5000-Possion3.0/Input-12000.txt", "QT5000-Possion5.0/Input-05562.txt"]
gd_rows(set) = set == "CF" ? filter(r -> r.id in CF_IDS, read_inputs_csv()) : inset(read_inputs_csv(), set)

function gd_params(id, variant)
    p = Params(read_qtsim_options(joinpath(SIMDIR, id)); Seed = nothing)
    variant == "" && return p
    variant == "many" && return Params(p; CrossoverModel = :all, RecombPoissonMean = 20.0)
    variant == "fec" && return Params(p; Selection = :fecundity)
    variant == "bigN" && return Params(p; PopSize = 4p.PopSize)
    variant == "sites10" && return Params(p; NoQTL = 10p.NoQTL, MutationRatePerGamete = p.MutationRatePerGamete / 10)
    error("unknown variant $variant")
end

# ---------------------------------------------------------------------------------------
# one pass over a population's genes with M weight vectors: per locus Σ w y, Σ w (x1² + x2²),
# Σ w y² (y = x1 + x2, values shifted by the locus's reference c_l), and Σ x³, Σ x⁴ (weight 1)

struct Pass
    S1::Matrix{Float64}; S2::Matrix{Float64}; Y2::Matrix{Float64}   # (M, L)
    S3::Vector{Float64}; S4::Vector{Float64}
    W0::Vector{Float64}
    ref::Vector{Float64}
    b1::Vector{Float64}; b2::Vector{Float64}
end
Pass(M, L) = Pass(zeros(M, L), zeros(M, L), zeros(M, L), zeros(L), zeros(L), zeros(M), zeros(L), zeros(L), zeros(L))

function loadchr!(buf, pool, chr, c, i)
    B = pool.B
    @inbounds for k in 1:pool.nb
        copyto!(buf, (k - 1) * B + 1, pool.vals, QTLSim.blockoffset(pool, chr[k, c, i]) + 1, QTLSim.blocklen(pool, k))
    end
    buf
end

function pass!(P::Pass, pool, chr, N, wts::AbstractMatrix; setref = false)
    M, L = size(P.S1)
    fill!(P.S1, 0); fill!(P.S2, 0); fill!(P.Y2, 0); fill!(P.S3, 0); fill!(P.S4, 0); fill!(P.W0, 0)
    setref && copyto!(P.ref, loadchr!(P.b1, pool, chr, 1, 1))
    b1, b2, ref = P.b1, P.b2, P.ref
    @inbounds for i in 1:N
        loadchr!(b1, pool, chr, 1, i); loadchr!(b2, pool, chr, 2, i)
        for m in 1:M
            P.W0[m] += wts[i, m]
        end
        for l in 1:L
            x1 = b1[l] - ref[l]; x2 = b2[l] - ref[l]
            y = x1 + x2; q = x1 * x1 + x2 * x2; y2 = y * y
            P.S3[l] += x1 * x1 * x1 + x2 * x2 * x2
            P.S4[l] += x1^4 + x2^4
            for m in 1:M
                w = wts[i, m]
                P.S1[m, l] += w * y; P.S2[m, l] += w * q; P.Y2[m, l] += w * y2
            end
        end
    end
    P
end

# per-locus gene mean and variance under weight m
@inline gmean(P, m, l) = P.S1[m, l] / (2P.W0[m])
@inline gvar(P, m, l) = P.S2[m, l] / (2P.W0[m]) - gmean(P, m, l)^2

# ---------------------------------------------------------------------------------------
# allele census of the current adults (dense genomes): per locus the distinct allele values

function census(w)
    G = genomes(w)
    L, Pl, N = size(G)
    buf = Vector{Float64}(undef, Pl * N)
    nal = 0.0; nloc = 0; v3 = 0.0; vtot = 0.0; v2common = 0.0
    fold = zeros(NBIN)                         # variance of non-major alleles by frequency (0, ½]
    top1 = 0.0                                 # Σ over loci of the largest single allele's share × locus var
    for l in 1:L
        k = 0
        for i in 1:N, c in 1:Pl
            k += 1; buf[k] = G[l, c, i]
        end
        sort!(buf)
        vals = Float64[]; cnt = Int[]
        for x in buf
            if isempty(vals) || x != vals[end]
                push!(vals, x); push!(cnt, 1)
            else
                cnt[end] += 1
            end
        end
        length(vals) == 1 && continue
        nloc += 1; nal += length(vals)
        f = cnt ./ length(buf)
        xbar = sum(f .* vals)
        v = sum(f .* (vals .- xbar) .^ 2)
        vtot += 2v
        length(vals) >= 3 && (v3 += 2v)
        count(>=(0.05), f) >= 3 && (v2common += 2v)
        j = argmax(cnt)
        best = 0.0
        for a in eachindex(vals)
            a == j && continue
            c = 2 * (vals[a] - vals[j])^2 * f[a] * (1 - f[a])
            fold[clamp(ceil(Int, f[a] * 2NBIN), 1, NBIN)] += c
            best = max(best, c)
        end
        top1 += min(best, 2v)
    end
    (; nseg_loci = nloc, alleles_per_seg = nloc > 0 ? nal / nloc : NaN, frac_multi = v3 / vtot,
     frac_3common = v2common / vtot, top1 = top1 / vtot, fold = fold ./ max(sum(fold), eps()))
end

# ---------------------------------------------------------------------------------------
# the instrumented kernel run

const COLS = [:gen, :N, :nsurv, :theta, :zbar, :VG, :VP, :Wbar, :beta,
              :GA, :Gpi, :GW, :GS, :GO, :Leff,
              :within, :k3, :own2, :phi, :own2LD, :hitch, :nonlin,
              :Esurv, :Erepro, :k3fav, :kurt, :lin3, :own3, :quad, :ownq, :sym3]

function gd_run(p, seed; from = GD_FROM, to = GD_TO, census_every = GD_CENSUS)
    q = kernel_params(p)
    rng = MersenneTwister(seed)
    w = if q.Selection === :viability
        equilibrium_world(q, rng; burnin = KBURN)
    else                                       # equilibrium_world's steps, genomes from the viability spectrum
        qv = Params(q; Selection = :viability)
        v = World(q)
        setgenomes!(v, equilibrium_genomes(qv, rng))
        QTLSim.init_env!(v, rng)
        v.gen = QTLSim.protocol_start(q; burnin = KBURN)
        v
    end
    L = q.NoQTL
    PA = Pass(7, L); PO = Pass(1, L)
    N0 = q.PopSize
    wts = zeros(N0, 7); one = ones(N0, 1)
    out = Dict(c => Float64[] for c in COLS)
    bins = Dict(s => zeros(0, NBIN) for s in (:own2, :hitch, :within, :genic))
    cens = []
    vgser = Float64[]                          # V_G of every generation from `from` (bitwise check)
    VE = QTLSim.envVar(q)
    binof(l) = min(NBIN, fld((l - 1) * NBIN, L) + 1)
    extinct_gen = -1
    GOprev = NaN
    for it in w.gen:q.NoGenerations-1
        it >= to && break
        w.gen = it
        rec = it >= from
        QTLSim.develop!(w)
        QTLSim.select!(w, rng)
        N = w.N
        if rec
            rec && it % census_every == 0 && push!(cens, (it, census(w)))
            g = view(w.genotype, 1:N)
            gb = mean(g); vg = var(g; corrected = false)
            push!(vgser, var(g))
            ω2 = w.ω^2; Sg = ω2 + VE
            for i in 1:N
                π = C0 * sqrt(ω2 / Sg) * exp(-(w.genotype[i] - w.θ)^2 / (2Sg))
                wts[i, 1] = 1.0; wts[i, 2] = π; wts[i, 3] = w.fitness[i]; wts[i, 4] = 0.0
                wts[i, 5] = π * (1 - π); wts[i, 6] = w.genotype[i] - gb; wts[i, 7] = (w.genotype[i] - gb)^2 - vg
            end
            if w.weighted                      # fecundity: the mating weights are the survivors' weights
                for i in 1:N
                    wts[i, 4] = w.fitness[i]
                end
            else
                for i in w.survivors
                    wts[i, 4] = 1.0
                end
            end
            pass!(PA, w.pool, w.chr, N, view(wts, 1:N, :); setref = true)
            πbar = PA.W0[2] / N
            β = sum(wts[i, 2] * wts[i, 6] for i in 1:N) / N / (πbar * vg)
            GA = 0.0; Gπ = 0.0; GW = 0.0; GS = 0.0; lin3 = 0.0; own3 = 0.0; quad = 0.0; ownq = 0.0; sk33 = 0.0; sk3o = 0.0
            # π/π̄ − 1 regressed on g′ and g′² − V jointly
            a11 = 0.0; a12 = 0.0; a22 = 0.0; r1 = 0.0; r2 = 0.0
            for i in 1:N
                u = wts[i, 6]; v = wts[i, 7]; yy = wts[i, 2] / πbar - 1
                a11 += u * u; a12 += u * v; a22 += v * v; r1 += u * yy; r2 += v * yy
            end
            dt = a11 * a22 - a12^2
            β2 = dt > 0 ? (a22 * r1 - a12 * r2) / dt : 0.0; q2 = dt > 0 ? (a11 * r2 - a12 * r1) / dt : 0.0
            within = 0.0; k3 = 0.0; k3fav = 0.0; so2 = 0.0; sk2 = 0.0; sok = 0.0; sdx2 = 0.0; sm4 = 0.0
            Es = 0.0; sv2 = 0.0
            ov = zeros(L); kv = zeros(L)
            bo = zeros(NBIN); bh = zeros(NBIN); bw = zeros(NBIN); bg = zeros(NBIN)
            dir = sign(w.θ - mean(view(w.phenotype, 1:N)))
            for l in 1:L
                xb = gmean(PA, 1, l); vA = gvar(PA, 1, l)
                vπ = gvar(PA, 2, l); vW = gvar(PA, 3, l); vS = gvar(PA, 4, l)
                GA += 2vA; Gπ += 2vπ; GW += 2vW; GS += 2vS
                m2 = PA.S2[1, l] / (2N); m3 = PA.S3[l] / (2N)
                κ3 = m3 - 3xb * m2 + 2xb^3
                m4 = PA.S4[l] / (2N)
                μ4 = m4 - 4xb * m3 + 6xb^2 * m2 - 3xb^4
                sm4 += 2μ4; sv2 += vA^2
                C1 = (PA.S2[6, l] - 2xb * PA.S1[6, l]) / (2N); C2 = (PA.S2[7, l] - 2xb * PA.S1[7, l]) / (2N)
                lin3 += 2β2 * C1; own3 += 2β2 * κ3; quad += 2q2 * C2; ownq += 2q2 * (μ4 - vA^2)
                sk33 += κ3^2; sk3o += κ3 * (C1 - κ3)
                vy = PA.Y2[1, l] / N - (PA.S1[1, l] / N)^2
                cgy = PA.S1[6, l] / N
                o = β * vy / 2; k = β * (cgy - vy) / 2
                dx = gmean(PA, 2, l) - xb
                ov[l] = o; kv[l] = k
                D2 = vπ - vA + dx^2
                within += 2D2; k3 += 2β * κ3; k3fav += 2dir * κ3
                so2 += o^2; sk2 += k^2; sok += o * k; sdx2 += dx^2
                # expected change by survival noise given π (independent Bernoulli(π_i), first order)
                c = gmean(PA, 2, l)                    # shifted coordinates, as the sums
                b = PA.W0[2]
                A1 = PA.Y2[5, l] - 4c * PA.S1[5, l] + 4c^2 * PA.W0[5]           # Σπ(1−π)(y − 2c)²
                A2 = PA.S2[5, l] - 2c * PA.S1[5, l] + 2c^2 * PA.W0[5]           # Σπ(1−π)Q′
                Es += 2 * (-A1 / (4b^2) - A2 / (2b^2) + vπ * PA.W0[5] / b^2)
                bn = binof(l)
                bo[bn] += -2o^2; bw[bn] += 2D2; bg[bn] += 2vA
            end
            φ = so2 > 0 ? sok / so2 : 0.0
            hitch = -2 * (sk2 - φ^2 * so2)
            for l in 1:L
                bh[binof(l)] += -2 * (kv[l] - φ * ov[l])^2
            end
            own2LD = -2 * (1 + φ)^2 * so2
            nonlin = (Gπ - GA) - within - (own2LD + hitch)    # = −2Σ[Δx̄² − (o + k)²]
            z = view(w.phenotype, 1:N)
            vals = (gen = it, N = N, nsurv = w.nsurv, theta = w.θ, zbar = mean(z), VG = var(g), VP = var(z),
                    Wbar = mean(view(w.fitness, 1:N)), beta = β, GA = GA, Gpi = Gπ, GW = GW, GS = GS, GO = NaN,
                    Leff = sv2 > 0 ? (GA / 2)^2 / sv2 : NaN,
                    within = within, k3 = k3, own2 = -2so2, phi = φ, own2LD = own2LD, hitch = hitch, nonlin = nonlin,
                    Esurv = Es, Erepro = NaN, k3fav = k3fav, kurt = sm4, lin3 = lin3, own3 = own3, quad = quad, ownq = ownq, sym3 = (lin3 - own3) - (sk33 > 0 ? sk3o / sk33 : 0.0) * own3)
            for c in COLS
                push!(out[c], Float64(vals[c]))
            end
            bins[:own2] = vcat(bins[:own2], bo'); bins[:hitch] = vcat(bins[:hitch], bh')
            bins[:within] = vcat(bins[:within], bw'); bins[:genic] = vcat(bins[:genic], bg')
        end
        w.extinct && (extinct_gen = it; break)
        QTLSim.reproduce!(w, rng)
        w.extinct && (extinct_gen = it; break)
        if rec
            pass!(PO, w.pool, w.schr, w.Noff, view(one, 1:w.Noff, :))
            GO = 0.0
            for l in 1:L
                GO += 2gvar(PO, 1, l)
            end
            out[:GO][end] = GO
            out[:Erepro][end] = -out[:GS][end] / (2w.Noff)
        end
        QTLSim.rescale!(w)
        QTLSim.mutate!(w, rng)
        QTLSim.swapbuffers!(w)
    end
    (; seed, extinct_gen, cols = out, bins, census = cens, vgser)
end

gd_file(id, r, variant) = joinpath(GD_RAW, replace(id, "/" => "__", ".txt" => "") * (variant == "" ? "" : "@" * variant) * "#$r.jls")

function gd_kernel(set, R, variant)
    mkpath(GD_RAW)
    rows = gd_rows(set)
    jobs = [(r.id, r.j, k) for r in rows for k in 1:R if !isfile(gd_file(r.id, k, variant))]
    variant == "sites10" && filter!(j -> gd_params(j[1], "").NoQTL <= 1000, jobs)
    sort!(jobs; by = j -> -gd_params(j[1], variant).NoQTL * gd_params(j[1], variant).PopSize)
    println("$(length(jobs)) instrumented runs for set $set, R = $R, variant \"$variant\", $(Threads.nthreads()) threads")
    done = Threads.Atomic{Int}(0); t0 = time()
    Threads.@threads :greedy for (id, j, k) in jobs
        tr = gd_run(gd_params(id, variant), KSEED + 100j + k)
        serialize(gd_file(id, k, variant), tr)
        d = Threads.atomic_add!(done, 1) + 1
        d % 10 == 0 && @printf("  %d/%d  %.0f s\n", d, length(jobs), time() - t0)
    end
    @printf("done: %d runs, %.0f s wall\n", length(jobs), time() - t0)
end

# bitwise identity with simulate (same seed): V_G of every instrumented generation
function gd_check(; id = "QT1000-Possion3.0/Input-09697.txt", j = 89, r = 1)
    for variant in ("", "many")
        p = gd_params(id, variant)
        t = @elapsed tr = gd_run(p, KSEED + 100j + r; to = MOVE + 60)
        o = simulate(kernel_params(p); seed = KSEED + 100j + r, init = :equilibrium, burnin = KBURN)
        k = [findfirst(==(g), o.gen) for g in GD_FROM:GD_FROM+length(tr.vgser)-1]
        same = all(o.geno_var[k] .== tr.vgser)
        @printf("check %s variant \"%s\": %d generations, V_G bitwise identical to simulate: %s (%.1f s)\n",
                id, variant, length(tr.vgser), same, t)
        same || error("instrumentation changed the run")
    end
end

# ---------------------------------------------------------------------------------------
# analysis: per input, the accounting in windows of the move, mean over runs alive

const WINDOWS = [(-200, -1), (0, 24), (25, 49), (50, 99), (100, 199), (200, 299)]

function load_gd(id, R, variant)
    out = []
    for k in 1:R
        f = gd_file(id, k, variant)
        isfile(f) && push!(out, deserialize(f))
    end
    out
end

# per-run per-generation accounting terms (genic units per generation), aligned on t = gen − MOVE
function terms(x)
    c = x.cols
    n = length(c[:gen])
    t = Int.(c[:gen]) .- MOVE
    GA = c[:GA]
    mut = [k < n ? GA[k+1] - c[:GO][k] : NaN for k in 1:n]
    (; t, GA, sel = c[:Gpi] .- GA, within = c[:within], k3 = c[:k3], stab = c[:within] .- c[:k3],
     own2LD = c[:own2LD], hitch = c[:hitch], nonlin = c[:nonlin],
     env = c[:GW] .- c[:Gpi], bern = c[:GS] .- c[:GW], Esurv = c[:Esurv],
     repro = c[:GO] .- c[:GS], Erepro = c[:Erepro], mut, phi = c[:phi], Leff = c[:Leff],
     VG = c[:VG], nsurv = c[:nsurv], N = c[:N], beta = c[:beta], lag = c[:theta] .- c[:zbar], k3fav = c[:k3fav],
     own3 = c[:own3], other3 = c[:lin3] .- c[:own3], quad = c[:quad], ownq = c[:ownq], bgq = c[:quad] .- c[:ownq],
     rest3 = c[:within] .- c[:lin3] .- c[:quad], sym3 = c[:sym3], prop3 = c[:lin3] .- c[:own3] .- c[:sym3])
end

const TERMS = [:sel, :within, :k3, :stab, :own2LD, :hitch, :nonlin, :env, :bern, :Esurv, :repro, :Erepro, :mut,
               :own3, :other3, :quad, :ownq, :bgq, :rest3, :sym3, :prop3]

# mean over runs alive and generations in window of term / G_A(−1) (relative to the genic variance before the move)
function window_means(runs, f; rel = true)
    res = Dict{Tuple{Int,Int},Float64}()
    for (a, b) in WINDOWS
        s = 0.0; n = 0
        for x in runs
            T = terms(x)
            k0 = findfirst(==(-1), T.t)
            k0 === nothing && continue
            g0 = mean(T.GA[max(1, k0 - 49):k0])
            v = getfield(T, f)
            for k in eachindex(T.t)
                a <= T.t[k] <= b || continue
                isfinite(v[k]) || continue
                s += rel ? v[k] / g0 : v[k]; n += 1
            end
        end
        res[(a, b)] = n > 0 ? s / n : NaN
    end
    res
end

# genic variance at t relative to its mean over −50..−1, mean over runs alive at t
function genic_ratio(runs, ts)
    [begin
        v = Float64[]
        for x in runs
            T = terms(x)
            k0 = findfirst(==(-1), T.t); k = findfirst(==(t), T.t)
            (k0 === nothing || k === nothing) && continue
            push!(v, T.GA[k] / mean(T.GA[max(1, k0 - 49):k0]))
        end
        (isempty(v) ? NaN : mean(v), length(v))
    end for t in ts]
end

# cumulative contribution of each term over 0..t−1, relative to G(−1..−50), mean over runs alive at t
function cumulative(runs, t)
    acc = Dict(f => Float64[] for f in TERMS)
    for x in runs
        T = terms(x)
        k0 = findfirst(==(-1), T.t); kt = findfirst(==(t), T.t)
        (k0 === nothing || kt === nothing) && continue
        g0 = mean(T.GA[max(1, k0 - 49):k0])
        for f in TERMS
            v = getfield(T, f)[k0+1:kt-1]
            push!(acc[f], sum(filter(isfinite, v)) / g0)
        end
    end
    Dict(f => (isempty(v) ? NaN : mean(v)) for (f, v) in acc)
end

function gd_analyse(set, R, variant; io = stdout)
    rows = gd_rows(set)
    println(io, "\n== accounting, set $set, variant \"$variant\" (terms per generation / genic(−50..−1); mean over runs alive) ==")
    for r in rows
        runs = load_gd(r.id, R, variant)
        isempty(runs) && continue
        p = gd_params(r.id, variant)
        gr = genic_ratio(runs, (50, 100, 200, 299))
        @printf(io, "\n%s [%s] L=%d u=%g  runs %d; genic/genic(pre): +50 %.3f (n%d) +100 %.3f (n%d) +200 %.3f (n%d) +299 %.3f (n%d)\n",
                r.id, r.set, p.NoQTL, p.RecombPoissonMean, length(runs), gr[1]..., gr[2]..., gr[3]..., gr[4]...)
        wm = Dict(f => window_means(runs, f) for f in TERMS)
        @printf(io, "  %-10s", "window")
        for f in TERMS
            @printf(io, " %8s", f)
        end
        println(io)
        for wdw in WINDOWS
            @printf(io, "  %4d..%-4d", wdw...)
            for f in TERMS
                @printf(io, " %8.5f", wm[f][wdw])
            end
            println(io)
        end
        for t in (100, 200)
            cu = cumulative(runs, t)
            tot = sum(cu[f] for f in (:sel, :env, :bern, :repro, :mut))
            @printf(io, "  cumulative 0..%d: total %+.3f = sel %+.3f [within %+.3f (k3 %+.3f stab %+.3f), own2LD %+.3f, hitch %+.3f, nonlin %+.3f] + env %+.3f + bern %+.3f (E %+.3f) + repro %+.3f (E %+.3f) + mut %+.3f\n",
                    t, tot, cu[:sel], cu[:within], cu[:k3], cu[:stab], cu[:own2LD], cu[:hitch], cu[:nonlin],
                    cu[:env], cu[:bern], cu[:Esurv], cu[:repro], cu[:Erepro], cu[:mut])
            @printf(io, "      other3 = κ3-proportional %+.3f + sign-symmetric %+.3f\n", cu[:prop3], cu[:sym3])
            @printf(io, "      within %+.3f = own3 %+.3f + other3 %+.3f + quad %+.3f (own kurtosis %+.3f, background %+.3f) + rest %+.3f\n",
                    cu[:within], cu[:own3], cu[:other3], cu[:quad], cu[:ownq], cu[:bgq], cu[:rest3])
        end
        # other diagnostics: φ, L_eff, census
        phi = window_means(runs, :phi; rel = false); le = window_means(runs, :Leff; rel = false)
        @printf(io, "  φ (systematic LD in the response): pre %.3f, 0..24 %.3f, 50..99 %.3f; L_eff pre %.1f, 50..99 %.1f, 100..199 %.1f\n",
                phi[(-200, -1)], phi[(0, 24)], phi[(50, 99)], le[(-200, -1)], le[(50, 99)], le[(100, 199)])
        for tc in (MOVE - 100, MOVE, MOVE + 100)
            cs = [c for x in runs for (g, c) in x.census if g == tc]
            isempty(cs) && continue
            @printf(io, "  census t=%+d (n%d): seg loci %.0f, alleles/seg locus %.2f, var at ≥3-allele loci %.3f, at loci with ≥3 alleles ≥5%% %.3f, top allele share %.3f, folded var by freq %s\n",
                    tc - MOVE, length(cs), mean(c.nseg_loci for c in cs), mean(c.alleles_per_seg for c in cs), mean(c.frac_multi for c in cs),
                    mean(c.frac_3common for c in cs), mean(c.top1 for c in cs), join([@sprintf("%.2f", v) for v in mean(c.fold for c in cs)], "/"))
        end
    end
end

# ---------------------------------------------------------------------------------------
# the spectrum's partition, driven by one instrumented kernel run: from t = −200, with the
# spectrum's equilibrium shape scaled to the run's genic variance, each generation selection
# with the run's lag λ and S = ω² + V_P, drift among its realised survivors (fecundity: as
# the lowering), and 2N L μ mutations; the terms: selection (first order, as spectrum_select:
# Σ 2a² pq (1 − 2p) E[Δp]/pq), Barton's linked drift −g c_d, drift −g/(2n), mutation, and
# the realised total change of the genic variance

function driven_partition(p, x; ld = :none)
    c = x.cols
    sp = QTLSim.spectrum_scale!(deepcopy(QTLSim.spectrum_equilibrium(p; ld).spectrum), c[:GA][1])
    ω2 = QTLSim.initFitSD(p)^2
    n = length(c[:gen])
    out = Dict(s => fill(NaN, n) for s in (:g, :sel, :barton, :drift, :mut, :total))
    for k in 1:n
        g = QTLSim.spectrum_vg(sp)
        out[:g][k] = g
        λ = c[:theta][k] - c[:zbar][k]; S = ω2 + c[:VP][k]
        QTLSim.spectrum_select(sp, λ, S)
        r1, r2 = QTLSim._ld_factors(sp, g)
        cc = g / (2sp.L); dg = 0.0
        for i in 2:length(sp.p)-1
            pi = sp.p[i]; pq = pi * (1 - pi)
            for j in eachindex(sp.a)
                e = QTLSim.spectrum_seff(sp.a[j], pi, λ, S, cc, r1, r2)
                dg += 2sp.a[j]^2 * pq * (1 - 2pi) * e * sp.n[j, i]
            end
        end
        npar = max(c[:nsurv][k], 1)
        QTLSim.spectrum_step!(sp, npar)
        out[:sel][k] = dg; out[:barton][k] = -g * sp.cd; out[:drift][k] = -g / (2npar)
        g1 = QTLSim.spectrum_vg(sp)
        QTLSim.spectrum_mutate!(sp, c[:N][k])
        g2 = QTLSim.spectrum_vg(sp)
        out[:mut][k] = g2 - g1; out[:total][k] = g2 - g
    end
    out
end

function gd_driven(set, R; io = stdout)
    rows = gd_rows(set)
    println(io, "\n== driven spectrum partition vs kernel accounting, set $set (cumulative over 0..t−1, / genic(−50..−1); mean over runs alive at t) ==")
    for r in rows
        runs = load_gd(r.id, R, "")
        isempty(runs) && continue
        p = gd_params(r.id, "")
        for ld in (:none, :dynamic)
            res = Dict(t => Dict{Symbol,Vector{Float64}}() for t in (50, 100, 200))
            for x in runs
                d = driven_partition(p, x; ld)
                T = terms(x)
                k0 = findfirst(==(-1), T.t)
                g0k = mean(T.GA[max(1, k0 - 49):k0]); g0s = mean(d[:g][max(1, k0 - 49):k0])
                for t in (50, 100, 200)
                    kt = findfirst(==(t), T.t); kt === nothing && continue
                    rg = k0+1:kt-1
                    push!(get!(res[t], :sp_ratio, Float64[]), d[:g][kt] / g0s)
                    push!(get!(res[t], :k_ratio, Float64[]), T.GA[kt] / g0k)
                    for (s, v) in ((:sp_sel, d[:sel]), (:sp_barton, d[:barton]), (:sp_drift, d[:drift]), (:sp_mut, d[:mut]))
                        push!(get!(res[t], s, Float64[]), sum(v[rg]) / g0s)
                    end
                    for (s, v) in ((:k_det, T.within .+ T.own2LD .+ T.nonlin), (:k_hitch, T.hitch),
                                   (:k_noise, T.env .+ T.bern .+ T.repro), (:k_mut, T.mut))
                        push!(get!(res[t], s, Float64[]), sum(filter(isfinite, v[rg])) / g0k)
                    end
                end
            end
            @printf(io, "%s [%s] ld=%s:", r.id, r.set, ld)
            for t in (50, 100, 200)
                isempty(res[t]) && continue
                m(s) = mean(res[t][s])
                @printf(io, "  +%d n%d genic k %.3f / sp %.3f | sel k %+.3f sp %+.3f | hitch k %+.3f barton sp %+.3f | drift k %+.3f sp %+.3f | mut k %+.3f sp %+.3f;",
                        t, length(res[t][:k_ratio]), m(:k_ratio), m(:sp_ratio), m(:k_det), m(:sp_sel), m(:k_hitch), m(:sp_barton),
                        m(:k_noise), m(:sp_drift), m(:k_mut), m(:sp_mut))
            end
            println(io)
        end
    end
end

# ---------------------------------------------------------------------------------------
# exploratory theory term: linked selection with its transient (time-dependent Barton term)
#
# A neutral site (frequency p) and the selected loci of class c (weight ω_c, split probability
# ρ_c). Drift among n parents creates LD between the site and a locus with E[D²] increment
# pq·v/(2n) per generation (v the locus's gene-level variance, Σ_k v_k = V/2, V the genic
# variance); an LD D created at generation s decays by λ_c = (1 − ρ_c)(1 − 1/(2n)) per generation
# and, while it lasts, moves p by b_t D λ_c^{t−s} in generation t (b = λ/S, the selection
# gradient). The cumulative linked change of p is a sum over creation times s of
# D_s Σ_{t′=s}^{t} b_{t′} λ^{t′−s}; the increment of its variance in generation t is pq·c(t),
#   c(t) = Σ_c ω_c b_t [2 R_c(t) + b_t Q_c(t)],
#   Q_c(t) = λ_c² Q_c(t−1) + g_t,          g_t = V_t/(4 n_t)  (current E[D²]/pq, summed weight),
#   R_c(t+1) = λ_c [R_c(t) + b_t Q_c(t)]   (persistence-weighted cumulated selection on existing LD).
# b_t Q_c is the per-generation Var(Δp) part (∝ 1/ρ: hitch of the accounting), 2 b_t R_c the
# autocorrelated part (∝ 1/ρ²: other3). At constant b, g it reaches b² g Σ ω_c/(1 − λ_c)² (to
# leading order), Barton's b² V Σ ω_c/(2(ρ_c + δ)²)/(2n) with δ = 1/(2n) and no selection term,
# which the spectrum's spectrum_linked applies at once (with δ = 1/(2n) + V/S); before the move
# b = 0, so R = 0 and the term starts from zero. `sel = true` adds the Gaussian dissolution of the
# association by selection, λ_c × (1 − V/S) (the spectrum's δ). Nothing is fitted.

mutable struct Transient
    Q::Vector{Float64}
    R::Vector{Float64}
end

function transient_c!(tr::Transient, sp, b, n, V, S; sel = false)
    g = V / (4n)
    c = 0.0
    for k in 2:length(sp.ω)
        λc = (1 - sp.ρ[k]) * (1 - 1 / (2n)) * (sel ? max(1 - V / S, 0.0) : 1.0)
        tr.Q[k] = λc^2 * tr.Q[k] + g
        c += sp.ω[k] * b * (2tr.R[k] + b * tr.Q[k])
        tr.R[k] = λc * (tr.R[k] + b * tr.Q[k])
    end
    max(c, 0.0)
end

"""
The spectrum's genic variance driven by one instrumented kernel run (its λ, S, survivors, N),
from t = −200 with the equilibrium shape scaled to the run's genic variance, with the linked
selection term `mode`: :barton (the closure's), :none, :transient (above), :transient_sel.
`ld` the spectrum's LD option. Returns the spectrum's genic variance at each recorded generation.
"""
function driven_linked(p, x; mode = :barton, ld = :none)
    c = x.cols
    eq = QTLSim.spectrum_equilibrium(p; ld)
    sp = QTLSim.spectrum_scale!(deepcopy(eq.spectrum), c[:GA][1])
    mode === :barton || (sp.linkage = :none)
    ω2 = QTLSim.initFitSD(p)^2
    K = length(sp.ω)
    tr = Transient(zeros(K), zeros(K))
    # the drift LD before the move at its steady state for the first recorded n and V
    n0 = max(c[:nsurv][1], 1); V0 = QTLSim.spectrum_vg(sp)
    for k in 2:K
        λc = (1 - sp.ρ[k]) * (1 - 1 / (2n0))
        tr.Q[k] = V0 / (4n0) / (1 - λc^2)
    end
    out = fill(NaN, length(c[:gen]))
    for k in eachindex(c[:gen])
        out[k] = QTLSim.spectrum_vg(sp)
        λ = c[:theta][k] - c[:zbar][k]; S = ω2 + c[:VP][k]
        QTLSim.spectrum_select(sp, λ, S)
        n = max(c[:nsurv][k], 1)
        if mode === :barton || mode === :none
            QTLSim.spectrum_step!(sp, n)
        else
            cc = transient_c!(tr, sp, λ / S, n, out[k], S; sel = mode === :transient_sel)
            QTLSim.spectrum_step!(sp, 1 / (1 / n + 2cc))
        end
        QTLSim.spectrum_mutate!(sp, c[:N][k])
    end
    out
end

const LINK_MODES = (:none, :barton, :transient_sel, :transient)

function gd_linked(set, R; io = stdout, ts = (50, 100, 200), ld = :none)
    rows = gd_rows(set)
    println(io, "\n== driven spectrum (ld = $ld) with linked-selection variants vs kernel, set $set: genic(t)/genic(−50..−1), mean over runs alive ==")
    allres = []
    for r in rows
        runs = load_gd(r.id, R, "")
        isempty(runs) && continue
        p = gd_params(r.id, "")
        acc = Dict((m, t) => Float64[] for m in (:kernel, LINK_MODES...) for t in ts)
        for x in runs
            T = terms(x)
            k0 = findfirst(==(-1), T.t)
            gk = mean(T.GA[k0-49:k0])
            sims = Dict(m => driven_linked(p, x; mode = m, ld) for m in LINK_MODES)
            for t in ts
                kt = findfirst(==(t), T.t); kt === nothing && continue
                push!(acc[(:kernel, t)], T.GA[kt] / gk)
                for m in LINK_MODES
                    s = sims[m]
                    push!(acc[(m, t)], s[kt] / mean(s[k0-49:k0]))
                end
            end
        end
        @printf(io, "%-34s [%-4s]", r.id, r.set)
        for t in ts
            isempty(acc[(:kernel, t)]) && continue
            @printf(io, "  +%d n%d k %.3f", t, length(acc[(:kernel, t)]), mean(acc[(:kernel, t)]))
            for m in LINK_MODES
                @printf(io, " %s %.3f", m === :transient_sel ? "trS" : m === :transient ? "tr" : m === :barton ? "bar" : "none", mean(acc[(m, t)]))
            end
        end
        println(io)
        push!(allres, (r, Dict(k => (isempty(v) ? NaN : mean(v)) for (k, v) in acc)))
    end
    # summary: log ratio model/kernel, median and mean |.| over inputs, at each t
    for t in ts
        @printf(io, "  +%d:", t)
        for m in LINK_MODES
            lr = [log(d[(m, t)] / d[(:kernel, t)]) for (r, d) in allres if isfinite(d[(m, t)]) && isfinite(d[(:kernel, t)]) && d[(:kernel, t)] > 0 && d[(m, t)] > 0]
            @printf(io, "  %s log(model/kernel) median %+.3f mean|.| %.3f (n%d);", m, median(lr), mean(abs.(lr)), length(lr))
        end
        println(io)
    end
    allres
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "check" : ARGS[1]
    variant = length(ARGS) >= 4 ? ARGS[4] : ""
    if mode == "check"
        gd_check()
    elseif mode == "kernel"
        gd_kernel(ARGS[2], parse(Int, ARGS[3]), variant)
    elseif mode == "analyse"
        gd_analyse(ARGS[2], parse(Int, ARGS[3]), variant)
    elseif mode == "driven"
        gd_driven(ARGS[2], parse(Int, ARGS[3]))
    elseif mode == "linked"
        gd_linked(ARGS[2], parse(Int, ARGS[3]); ld = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :none)
    end
end
