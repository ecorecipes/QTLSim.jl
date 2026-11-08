# Several quantitative traits with pleiotropy (Lande 1979, 1980): the kernel's operators
# for `NoTraits = n > 1`. Each is the image of the same step of the life cycle as its
# one-trait counterpart in model.jl, which calls it when the `World` carries traits; the
# one-trait code is untouched. A locus holds a vector of n effects (blocks.jl), an
# individual's genotype is the vector of their sums, its phenotype z = g + e with
# e ~ N(0, V_E I), fitness is QTSim's Gaussian in n dimensions,
# W(z) = c·exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ)) with c = 1/(0.4√(2π)), and a mutation adds a draw
# from N(0, M/(Ploidy·NoQTL·rate)). The scalar fields of the `World` (`genotype`,
# `phenotype`, `θ`) hold trait 1, so that the columns QTSim records describe trait 1.

# QTSim's fitness at the optimum, 1/(0.4√(2π)) ≈ 0.9974, for any number of traits
const TRAIT_FITNESS_PEAK = 1 / (0.4 * sqrt(2π))

# develop! with several traits: g = Σ of the blocks' sum vectors, chromosome by chromosome,
# block by block, and z = g + e
function develop_traits!(w::World)
    t = w.traits::TraitState
    n, sums, chr = t.n, w.pool.sums, w.chr
    nb, P = size(chr, 1), size(chr, 2)
    m = nb * P
    @inbounds for i in 1:w.N
        for r in 1:n
            t.G[r, i] = 0.0
        end
        for s in m*(i-1)+1:m*i
            o = (Int(chr[s]) - 1) * n
            for r in 1:n
                t.G[r, i] += sums[o + r]
            end
        end
        for r in 1:n
            t.Z[r, i] = t.G[r, i] + t.E[r, i]
        end
        w.genotype[i] = t.G[1, i]
        w.phenotype[i] = t.Z[1, i]
    end
    w
end

# the mean phenotype vector of the current population, into t.d
function mean_traits!(t::TraitState, N::Int)
    fill!(t.d, 0.0)
    @inbounds for i in 1:N, r in 1:t.n
        t.d[r] += t.Z[r, i]
    end
    t.d ./= N
    t.d
end

# moveoptimum! with several traits (called from the move's first generation on): along
# the unit direction u, θ = z̄ + √V_P·d·u (:phenotype), z̄ + ω·d·u (:fitness), or θ += d·u
function moveoptimum_traits!(w::World, rng = nothing)
    p, t = w.p, w.traits::TraitState
    d = p.MoveFitnessSD
    if p.MoveType === :continuous
        t.θ .+= d .* t.u
    else
        zbar = mean_traits!(t, w.N)
        lead = p.MoveType === :phenotype ? sqrt(p.PhenotypeVar) * d : w.ω * d
        t.θ .= zbar .+ lead .* t.u
    end
    w.θ = t.θ[1]
    w
end

# Gaussian fitness of each individual, c·exp(−½ (z−θ)ᵀ Ω⁻¹ (z−θ)) with Ω⁻¹ = S⁻¹/ω²
function fitness_traits!(w::World)
    t = w.traits::TraitState
    if w.p.Drift
        fill!(view(w.fitness, 1:w.N), 1.0)
        return w
    end
    n, Sinv, x = t.n, t.Sinv, t.ξ
    iω2 = 1 / w.ω^2
    @inbounds for i in 1:w.N
        for r in 1:n
            x[r] = t.Z[r, i] - t.θ[r]
        end
        q = 0.0
        for a in 1:n, b in 1:n
            q += x[a] * Sinv[a, b] * x[b]
        end
        w.fitness[i] = TRAIT_FITNESS_PEAK * exp(-0.5 * q * iω2)
    end
    w
end

# mutate! with several traits: Poisson(L·rate) mutations per chromosome, each adding
# Lmut·ξ, ξ ~ N(0, I_n), to the effects of a uniformly chosen locus (copy on write, as for
# one trait)
function mutate_traits!(w::World, rng)
    p, t = w.p, w.traits::TraitState
    L, n = p.NoQTL, t.n
    λ = L * p.MutationRatePerGamete
    pool, schr, Lm, ξ = w.pool, w.schr, t.Lmut, t.ξ
    B = pool.B
    @inbounds for j in 1:w.Noff
        for c in 1:p.Ploidy
            no = rand(rng, Poisson(λ))
            for _ in 1:no
                pos = rand(rng, 1:L)
                for r in 1:n
                    ξ[r] = randn(rng)
                end
                k = (pos - 1) ÷ B + 1
                len = blocklen(pool, k)
                id = schr[k, c, j]
                if pool.refs[id] > 1     # shared: copy on write
                    nid = copyblock!(pool, pool, id, len)
                    pool.refs[id] -= Int32(1)
                    schr[k, c, j] = id = nid
                end
                o = blockoffset(pool, id) + (pos - (k - 1) * B - 1) * n
                for r in 1:n
                    x = 0.0
                    for s in 1:n
                        x += Lm[r, s] * ξ[s]
                    end
                    pool.vals[o + r] += x
                end
                blocksum!(pool, id, len)
            end
        end
    end
    w
end

# record! with several traits: the trait means, the G and P matrices (sample covariances,
# n − 1), the lag θ − z̄ and the selection differential vector
function record_traits!(rec, w::World)
    t = w.traits::TraitState
    N = w.N
    G, Z = view(t.G, :, 1:N), view(t.Z, :, 1:N)
    zbar = vec(mean(Z; dims = 2))
    push!(rec.trait_mean, zbar)
    push!(rec.gmatrix, cov(G; dims = 2))
    push!(rec.pmatrix, cov(Z; dims = 2))
    push!(rec.trait_lag, t.θ .- zbar)
    zs = w.weighted ? (Z * view(w.fitness, 1:N)) ./ w.wsum : vec(mean(view(t.Z, :, w.survivors); dims = 2))
    push!(rec.trait_sel_diff, zs .- zbar)
    rec
end

"""
    exchange_env!(w, i, v, j)

Exchange the environmental deviations of individual `i` of `w` and individual `j` of `v`
(all of its traits' when there are several), for migration.
"""
function exchange_env!(w::World, i::Integer, v::World, j::Integer)
    w.env[i], v.env[j] = v.env[j], w.env[i]
    tw, tv = w.traits, v.traits
    if tw !== nothing && tv !== nothing
        @inbounds for r in 1:tw.n
            tw.E[r, i], tv.E[r, j] = tv.E[r, j], tw.E[r, i]
        end
    end
    nothing
end

"""
    traitvalues(w; what = :phenotype)

The traits of the current population as an `n × N` matrix (a copy): `:phenotype`,
`:genotype` or `:env`. With one trait, a `1 × N` matrix of the scalar fields.
"""
function traitvalues(w::World; what::Symbol = :phenotype)
    t = w.traits
    if t === nothing
        v = what === :phenotype ? w.phenotype : what === :genotype ? w.genotype :
            what === :env ? w.env : throw(ArgumentError("what is :phenotype, :genotype or :env; got :$what"))
        return reshape(v[1:w.N], 1, w.N)
    end
    M = what === :phenotype ? t.Z : what === :genotype ? t.G :
        what === :env ? t.E : throw(ArgumentError("what is :phenotype, :genotype or :env; got :$what"))
    M[:, 1:w.N]
end
