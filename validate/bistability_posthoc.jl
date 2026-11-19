# Post hoc diagnostic of validate/bistability.jl (written after its verdicts; exploratory, no criterion): why does the
# Gaussian arm die below the deterministic critical rate? The genic variance (2 × the sum over loci of the variance among
# gene copies) of six populations of the cell G-K1000-R2 under a move at 0.6 k_c, from the pool's own variance at
# N₀ = min(800, pool size), with the design's genetics and three variants. Seeds 999_100 + r.
#
#   julia -t 4 --project=vignettes validate/bistability_posthoc.jl      → validate/bistability_posthoc.log

include(joinpath(@__DIR__, "bistability.jl"))

function posthoc()
    ci = cellinfo(:G, 1000, 2.0)
    k = 0.6 * ci.kc
    say(@sprintf("validate/bistability_posthoc.jl — cell G-K1000-R2, k = 0.6 k_c = %.4f per generation, 6 populations per line, 500 generations (post hoc, exploratory)", k))
    variants = ("design (50 loci, μ = 0.05, Poisson(5) crossovers)" => (ci.p, k),
                "Poisson(20) crossovers" => (Params(ci.p; RecombPoissonMean = 20.0), k),
                "500 loci, μ = 0.005, Poisson(20) (same V_m and α²)" => (Params(ci.p; NoQTL = 500, MutationRatePerGamete = 0.005, RecombPoissonMean = 20.0), k),
                "design, no move" => (ci.p, 0.0))
    genic = w -> (G = genomes(w); 2 * sum(var(vec(G[l, :, :])) for l in 1:size(G, 1)))
    for (name, (p, kk)) in variants
        res = Vector{Any}(undef, 6)
        Threads.@threads for r in 1:6
            pl = pool_population(p, 999_100 + r)
            q = Params(p; NoGenerations = 500, MoveFitnessSDGen = 0, MoveFitnessSD = kk)
            rng = Xoshiro(100 + r)
            w, _ = start_world(q, pl.G, min(800, pl.N), pl.VG, 0.0, rng)
            g0 = genic(w); o = Float64[]; t = 0
            for g in 1:500
                generation!(w, rng); w.extinct && break; t = g
                g in (100, 200, 300, 400, 500) && push!(o, genic(w) / g0)
            end
            res[r] = (; g0, o, t, N = w.N, n = w.nsurv)
        end
        m(i) = (v = [r.o[i] for r in res if length(r.o) >= i]; isempty(v) ? NaN : mean(v))
        say(@sprintf("  %-52s genic variance at the start %.3f; relative to it at generations 100, 200, 300, 400, 500: %.2f %.2f %.2f %.2f %.2f; alive at 500: %d/6, N then %.0f",
                     name, mean(r.g0 for r in res), m(1), m(2), m(3), m(4), m(5), count(r -> r.t == 500, res), (v = [r.N for r in res if r.t == 500]; isempty(v) ? 0.0 : mean(v))))
    end
    say("  Drift alone among about 900 breeders would lose 1/(2·900) of the genic variance per generation (to 0.76 in 500 generations without mutation).")
    write(joinpath(HERE, "bistability_posthoc.log"), String(take!(out)))
end

abspath(PROGRAM_FILE) == (@__FILE__) && posthoc()
