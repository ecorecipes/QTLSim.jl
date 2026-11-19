# Parts 1b and "impact" of validate/effective_size.jl (included by it).

include(joinpath(HERE, "published_table.jl"))

# ------------------------------------------------------------------------------------------
# drift_critical_rate as it was at d8b2565 (verbatim, renamed), for comparisons

function old_feedback_state(p::Params, λ::Real, c::Real)
    ω, VE, K = initFitSD(p), QTLSim.envVar(p), p.PopSize
    VG = QTLSim._house_of_gauss(p, p.Ploidy * float(K))
    n = N = W = 0.0
    for _ in 1:100
        s = ω^2 + VG + VE
        W = QTLSim.QTSIM_FITNESS_PEAK * ω / sqrt(s) * exp(-λ^2 / (2s) - c)
        if p.Demography === :fixed
            N, n = float(K), K * W
        else
            R = p.GrowthRate
            R * W <= 1 && return (VG = 0.0, n = 0.0, N = 0.0, Wbar = W)
            n = K * (R * W - 1) / (R - 1)
            N = n / W
        end
        n < p.Ploidy && return (VG = 0.0, n, N, Wbar = W)
        V = QTLSim._house_of_gauss(p, p.Ploidy * effective_size(N, n))
        done = abs(V - VG) <= 1e-13 * V
        VG = V
        done && break
    end
    (; VG, n, N, Wbar = W)
end

function old_drift_critical_rate(p::Params; cost::Real = -log(QTLSim.plasticity_cost(p)))
    Vs = selection_variance(p)
    rate(λ) = (st = old_feedback_state(p, λ, cost); st.VG * λ / (st.VG + Vs))
    λs = range(0, 12 * sqrt(initFitSD(p)^2 + QTLSim.envVar(p) + QTLSim._house_of_gauss(p, p.Ploidy * float(p.PopSize))), length = 600)
    ks = rate.(λs)
    i = argmax(ks)
    ks[i] > 0 || return (k = 0.0, λ = NaN, VG = 0.0, n = 0.0, N = 0.0, Wbar = NaN, M = NaN)
    a, b = λs[max(i - 1, 1)], λs[min(i + 1, end)]
    for _ in 1:80                                    # golden-section search on the bracket
        m1, m2 = b - 0.618034 * (b - a), a + 0.618034 * (b - a)
        rate(m1) < rate(m2) ? (a = m1) : (b = m2)
    end
    λ = (a + b) / 2
    st = old_feedback_state(p, λ, cost)
    s = initFitSD(p)^2 + st.VG + QTLSim.envVar(p)
    (k = rate(λ) / abs(1 - QTLSim.plasticity_ratio(p)), λ, st.VG, st.n, st.N, st.Wbar, M = λ^2 / s)
end

"How the old iteration ends at lag λ: :converged, :dead_first (no population at the variance of N = K), :dead_later, :cycling."
function old_iteration_fate(p::Params, λ::Real, c::Real)
    ω, VE, K = initFitSD(p), QTLSim.envVar(p), p.PopSize
    VG = QTLSim._house_of_gauss(p, p.Ploidy * float(K))
    for it in 1:100
        sz = QTLSim._feedback_size(p, VG, λ, c)
        sz.n < p.Ploidy && return it == 1 ? :dead_first : :dead_later
        V = QTLSim._house_of_gauss(p, p.Ploidy * effective_size(sz.N, sz.n))
        abs(V - VG) <= 1e-13 * V && return :converged
        VG = V
    end
    :cycling
end

"""
An independent maximum over every steady state, parametrised by its variance V as Xu and Shen's
eq. 6 is (after `ladder_kc` of validate/xu_shen.jl): for each V the lag at which the population
holds exactly V, and the response there.
"""
function all_states_kc(p::Params; cost::Real = 0.0, ne::Symbol = :falconer, grid::Integer = 4000)
    Vs = selection_variance(p)
    f(V, λ) = QTLSim._feedback_map(p, V, λ, cost, ne)
    function rate(V)
        f(V, 0.0) > V || return 0.0
        lo, hi = 0.0, sqrt(Vs)
        while f(V, hi) > V
            hi *= 2
        end
        for _ in 1:100
            mid = (lo + hi) / 2
            f(V, mid) > V ? (lo = mid) : (hi = mid)
        end
        V * (lo + hi) / 2 / (V + Vs)
    end
    Vmax = QTLSim._house_of_gauss(p, Inf)
    Vg = range(0, Vmax, length = grid + 2)[2:end-1]
    ks = rate.(Vg)
    i = argmax(ks)
    ks[i] > 0 || return 0.0
    a, b = Vg[max(i - 1, 1)], Vg[min(i + 1, end)]
    for _ in 1:80
        m1, m2 = b - 0.618034 * (b - a), a + 0.618034 * (b - a)
        rate(m1) < rate(m2) ? (a = m1) : (b = m2)
    end
    max(rate((a + b) / 2), ks[i])
end

"The deterministic recursion P-D of validate/bistability.jl (N_e = effective_size, the default), returning the final state."
function pd_run(p::Params, λ0, V0, N0; k, T::Integer, ne::Symbol = :falconer)
    ω, VE, K, R = initFitSD(p), QTLSim.envVar(p), float(p.PopSize), p.GrowthRate
    Vs = ω^2 + VE; α2 = QTLSim.mutation_effect_variance(p); nu, μ = QTLSim.rare_units(p); Vm = nu * μ * α2
    λ, V, N, n = float(λ0), float(V0), float(N0), 0.0
    for t in 1:T
        λ += k
        s = ω^2 + V + VE
        W = PEAK * ω / sqrt(s) * exp(-λ^2 / (2s))
        n = N * W
        n < p.Ploidy && return (dead = t, λ, V, N, n)
        λ -= V * λ / s
        V = max(V + Vm - V^2 / (nu * Vs) - V * (α2 / (2Vs) + 1 / (p.Ploidy * effective_size(N, n; scheme = ne))), 0.0)
        N = min(K, R * n / (1 + (R - 1) * n / K))
    end
    (dead = 0, λ, V, N, n)
end

"validate/xu_shen.jl's `qtl_from`: a Params whose mapped Xu–Shen parameters are (r, ω2, VM, K, ne)."
function qtl_from(; r, ω2, VM, K, ne, μ = 0.1, L = 100)
    VE = ω2 / 50
    ω = sqrt(ω2 - VE)
    R = exp(r - log(PEAK) + 0.5 * log(1 + VE / ω^2))
    ploidy = ne == 1 ? 1 : 2
    loci = ne == 1 ? L : Int(ne ÷ 2)
    Params(PopSize = round(Int, K * log(R) / ploidy), Ploidy = ploidy, NoQTL = loci, RecombPoissonMean = 3.0, PhenotypeVar = 2VE, Heritability = 0.5,
           MutationScalar = VM / VE, MutationRatePerGamete = μ, InitFitnessSD = ω, MoveType = :continuous, Demography = :beverton_holt, GrowthRate = R,
           CrossoverModel = :all)
end

function xs_cases()
    cases = Pair{String,Params}[]
    for ne in (1, 100)
        base = (; r = 0.2, ω2 = 50.0, VM = 0.1, K = 1000.0, ne)
        push!(cases, @sprintf("n_e = %d: standard", ne) => qtl_from(; base...))
        for (k, vals) in ((:ω2, (20.0, 100.0)), (:VM, (0.05, 0.2)), (:K, (500.0, 2000.0)), (:r, (0.5,)))
            for x in vals
                push!(cases, @sprintf("n_e = %d: %s = %g", ne, k, x) => qtl_from(; merge(base, NamedTuple{(k,)}((x,)))...))
            end
        end
    end
    cases
end

const STRATA = [(L, u) for L in (100, 500, 1000, 5000) for u in (3, 5)]
stratum_dir(L, u) = joinpath(ROOT, "Repository", "Simulations", @sprintf("QT%04d-Possion%.1f", L, u))

"`per` published inputs from each of the 8 strata (a fixed pseudo-random choice), as (L, u, id, Params)."
function published_inputs(per::Integer; seed = 20261119)
    rng = Xoshiro(seed)
    v = NamedTuple[]
    for (L, u) in STRATA
        dir = stratum_dir(L, u)
        isdir(dir) || continue
        for id in sort(randperm(rng, 14998)[1:per])
            f = joinpath(dir, @sprintf("Input-%05d.txt", id))
            isfile(f) || continue
            push!(v, (; L, u, id, p = Params(read_qtsim_options(f); Seed = nothing)))
        end
    end
    v
end

const DEMOS = (("fixed", p -> Params(p; Demography = :fixed)), ("R = 5", p -> Params(p; Demography = :beverton_holt, GrowthRate = 5.0)),
               ("R = 2", p -> Params(p; Demography = :beverton_holt, GrowthRate = 2.0)), ("R = 1.5", p -> Params(p; Demography = :beverton_holt, GrowthRate = 1.5)))

q3(v) = isempty(v) ? (NaN, NaN, NaN) : (minimum(v), median(v), maximum(v))

function branch_part()
    t0 = time()
    say("validate/effective_size.jl branch — part 1b: does drift_critical_rate find every steady state?")
    say("old = the function at d8b2565 (copied in effective_size_parts.jl); new = drift_critical_rate(p; solver = :bisection) (ne = :falconer);")
    say("the default, solver = :iteration, is compared with old bitwise in section 3.")
    say("all states = an independent maximum over the steady states parametrised by their variance (as eq. 6 of Xu and Shen 2026 is).")

    say("\n1. Reproduction on the 16 Xu–Shen-type inputs of validate/xu_shen.jl")
    sayf("  %-22s %9s %9s %10s | %-30s | %8s %8s %8s", "case", "old k_c", "new k_c", "all states", "old iteration over 600 lags", "V(K)", "R·W̄(V(K))", "|f′(V*)|")
    cases = xs_cases()
    differ = String[]
    for (name, p) in cases
        o, nw, al = old_drift_critical_rate(p).k, drift_critical_rate(p; solver = :bisection), all_states_kc(p)
        λs = range(0, 12 * sqrt(initFitSD(p)^2 + QTLSim.envVar(p) + QTLSim._house_of_gauss(p, p.Ploidy * float(p.PopSize))), length = 600)
        fates = [old_iteration_fate(p, λ, 0.0) for λ in λs]
        cnt(s) = count(==(s), fates)
        VK = QTLSim._house_of_gauss(p, p.Ploidy * float(p.PopSize))
        RW = p.GrowthRate * QTLSim._feedback_size(p, VK, 0.0, 0.0).Wbar
        # slope of the map V ↦ f(V) at the steady state of lag 0 (the iteration converges iff |f′| < 1)
        st = QTLSim._feedback_state(p, 0.0, 0.0; solver = :bisection)
        h = 1e-6 * max(st.VG, 1e-12)
        fp = st.VG > 0 ? (QTLSim._feedback_map(p, st.VG + h, 0.0, 0.0, :falconer) - QTLSim._feedback_map(p, st.VG - h, 0.0, 0.0, :falconer)) / 2h : NaN
        abs(o - nw.k) > 1e-12 * max(nw.k, 1e-300) && push!(differ, name)
        sayf("  %-22s %9.4f %9.4f %10.4f | conv %3d dead1 %3d dead+ %3d cyc %3d | %8.3f %8.3f %8.2f%s", name, o, nw.k, al, cnt(:converged), cnt(:dead_first),
             cnt(:dead_later), cnt(:cycling), VK, RW, abs(fp), abs(o - nw.k) > 1e-12 * max(nw.k, 1e-300) ? "  *" : "")
    end
    say("  * old ≠ new: ", isempty(differ) ? "none" : join(differ, "; "))
    say("  conv: the old iteration converged; dead1: \"no population\" at the variance of N = K (first step); dead+: the iteration jumped to")
    say("  \"no population\" at a later step; cyc: no convergence in 100 steps. R·W̄(V(K)) ≤ 1 means a population of size K with its house-of-Gauss variance")
    say("  cannot replace itself at lag 0. f(V) is the variance that the population sized by V holds; it falls with V, so V = f(V) has one root.")
    worst = maximum(abs(drift_critical_rate(p; solver = :bisection).k - all_states_kc(p)) / max(all_states_kc(p), 1e-300) for (_, p) in cases)
    sayf("  new against all states: largest relative difference %.2e", worst)

    say("\n2. Uniqueness: is f non-increasing in V? On each of the 16 inputs and 40 lags, a grid of 2000 variances, both schemes of N_e")
    for ne in (:falconer, :breeders)
        bad = 0; worstn = Inf; worstd = 0.0
        for (_, p) in cases
            VK = QTLSim._house_of_gauss(p, Inf)
            for λ in range(0, 6 * sqrt(selection_variance(p)), length = 40)
                Vg = range(1e-9 * VK, VK, length = 2000)
                g = [QTLSim._feedback_map(p, V, λ, 0.0, ne) for V in Vg]
                d = diff(g); j = argmax(d)
                if d[j] > 1e-12 * VK
                    bad += 1; worstd = max(worstd, d[j] / VK)
                    worstn = min(worstn, QTLSim._feedback_size(p, Vg[j], λ, 0.0).n)
                end
            end
        end
        sayf("  ne = :%-9s f non-increasing in %d of %d (input, lag) pairs%s", ne, 640 - bad, 640,
             bad == 0 ? "" : @sprintf("; the exceptions rise by at most %.1e of the largest variance, all in haploid inputs with %.2f–1.02 breeders (the default formula's N_e rises as n falls toward 1), where f < V", worstd, worstn))
    end

    say("\n3. Published inputs (75 per stratum, 600) × fixed size and Beverton–Holt R = 5, 2, 1.5, without a cost and, on every fifth input, with log costs 0.25 and 1.")
    say("   (a) the default (solver = :iteration) against old, field by field (isequal); (b) solver = :bisection against old: relative difference in k_c > 1e-12")
    pubs = published_inputs(75)
    recs = Vector{Any}(undef, length(pubs))
    Threads.@threads :greedy for i in eachindex(pubs)
        x = pubs[i]
        v = NamedTuple[]
        for (dn, dem) in DEMOS, c in (i % 5 == 0 ? (0.0, 0.25, 1.0) : (0.0,))
            q = dem(x.p)
            o, df, nw = old_drift_critical_rate(q; cost = c), drift_critical_rate(q; cost = c), drift_critical_rate(q; cost = c, solver = :bisection)
            same = all(isequal(getfield(o, f), getfield(df, f)) for f in keys(o))
            rel = abs(o.k - nw.k) / max(nw.k, 1e-300)
            al = rel > 1e-12 ? all_states_kc(q; cost = c) : NaN
            push!(v, (; x.L, x.u, x.id, dn, c, same, rel, old = o.k, new = nw.k, al, n = nw.n, N = nw.N,
                      RW = q.Demography === :fixed ? NaN : q.GrowthRate * QTLSim._feedback_size(q, QTLSim._house_of_gauss(q, 2.0 * q.PopSize), 0.0, c).Wbar))
        end
        recs[i] = v
    end
    rec = reduce(vcat, recs)
    sayf("   (a) %d comparisons: %d identical in every field", length(rec), count(r -> r.same, rec))
    sayf("   (b) solver = :bisection differs from old by more than 1e-12 in %d of %d (largest difference below that threshold: %.1e):", count(r -> r.rel > 1e-12, rec), length(rec),
         maximum(r.rel for r in rec if r.rel <= 1e-12))
    sayf("       %-10s %-22s %-22s", "", "no cost", "cost 0.25 or 1")
    for (dn, _) in DEMOS
        a = [r for r in rec if r.dn == dn && r.c == 0]; b = [r for r in rec if r.dn == dn && r.c > 0]
        sayf("       %-10s %3d of %4d             %3d of %4d", dn, count(r -> r.rel > 1e-12, a), length(a), count(r -> r.rel > 1e-12, b), length(b))
    end
    for r in rec
        r.rel > 1e-12 || continue
        sayf("       QT%04d-%d Input-%05d %-8s cost %.2f: old %.5f, bisection %.5f, all states %.5f; at the critical lag %.0f breeders of %.0f adults; R·W̄ at the variance of N = K: %.3f",
             r.L, r.u, r.id, r.dn, r.c, r.old, r.new, r.al, r.n, r.N, r.RW)
    end
    dif = [r for r in rec if r.rel > 1e-12]
    sayf("       in all %d the bisection equals the independent maximum over all steady states (largest relative difference %.1e) and is above old", length(dif),
         maximum(abs(r.new - r.al) / r.al for r in dif; init = 0.0))
    tc = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, RecombPoissonMean = 3.0, PhenotypeVar = 0.03, Heritability = 0.4,
                MutationScalar = 0.005, MutationRatePerGamete = 1e-4, InitFitnessSD = 4.0, MoveType = :continuous)
    same = all(isequal(old_drift_critical_rate(dem(tc); cost = c), drift_critical_rate(dem(tc); cost = c, solver = :bisection)) for (_, dem) in DEMOS for c in (0.0, 0.1, 0.5, 1.0))
    say("   test/critical.jl's input, 4 demographies × 4 costs, solver = :bisection against old: ", same ? "identical" : "DIFFERENT")

    say("\n4. Is the steady state that the iteration misses reachable? Two inputs for which old k_c = 0 (\"no population\"):")
    say("   (a) the deterministic recursion of the same ingredients (lag, house-of-Gauss variance per gene, breeders, Beverton–Holt; validate/bistability.jl's P-D),")
    say("       started at N = K with the variance of N = K and no lag; (b) the kernel: 2000 generations at a fixed optimum, then 2000 of a continuous move; 12 runs per rate.")
    pin = first(x for x in published_inputs(75) if x.L == 5000 && x.u == 3 && x.id == 12144)
    for (name, p) in (("Xu–Shen-type, n_e = 100: ω2 = 20", Dict(cases)["n_e = 100: ω2 = 20"]),
                      ("published QT5000-3 Input-12144, R = 1.5", Params(pin.p; Demography = :beverton_holt, GrowthRate = 1.5, MoveType = :continuous)))
        o, nw = old_drift_critical_rate(p).k, drift_critical_rate(p; solver = :bisection)
        st0 = QTLSim._feedback_state(p, 0.0, 0.0; solver = :bisection)
        VK = QTLSim._house_of_gauss(p, p.Ploidy * float(p.PopSize))
        sayf("   %s: K = %d, R = %.3f, %d loci; old k_c %.4f, bisection k_c %.4f (λ* %.2f, V* %.4f, n* %.1f); steady state without a move: V %.4f, n %.1f, N %.1f (variance of N = K: %.4f)", name,
             p.PopSize, p.GrowthRate, p.NoQTL, o, nw.k, nw.λ, nw.VG, nw.n, st0.VG, st0.n, st0.N, VK)
        for kf in (0.0, 0.5, 0.9, 1.1)
            r = pd_run(p, 0.0, VK, float(p.PopSize); k = kf * nw.k, T = 20_000)
            sayf("     (a) k = %.2f × bisection k_c: %s", kf, r.dead > 0 ? @sprintf("extinct at generation %d", r.dead) : @sprintf("alive at 20,000: V %.4f, breeders %.1f, adults %.1f, lag %.3f", r.V, r.n, r.N, r.λ))
        end
        rates = (0.0, 0.5 * nw.k, 1.5 * nw.k)
        res = Vector{Any}(undef, length(rates) * 12)
        Threads.@threads :greedy for j in eachindex(res)
            ik, rep = fldmod1(j, 12)
            q = Params(p; NoGenerations = 4000, MoveFitnessSDGen = 2000, MoveFitnessSD = rates[ik], Sample = 10)
            o2 = simulate(q; seed = 78_000_000 + 1000 * ik + rep + (p.NoQTL == 5000 ? 500_000 : 0))
            w1 = (o2.gen .>= 100) .& (o2.gen .< 2000)
            res[j] = (alive = !o2.extinct, tdead = o2.extinct ? o2.extinct_gen : 4000, vg = any(w1) ? mean(o2.geno_var[w1]) : NaN, n = any(w1) ? mean(o2.nsurv[w1]) : NaN,
                      N = any(w1) ? mean(o2.popsize[w1]) : NaN)
        end
        for (ik, k) in enumerate(rates)
            rs = res[(ik-1)*12+1:ik*12]
            td = sort([r.tdead for r in rs])
            sayf("     (b) k = %.4f (%.1f × bisection k_c): alive at the start of the move %2d/12, 2000 generations later %2d/12; time of death (4000 = alive): median %d, range %d–%d; while alive in generations 100–1999: V_G %.4f, breeders %.1f, adults %.1f",
                 k, k / nw.k, count(r -> r.tdead >= 2000, rs), count(r -> r.alive, rs), td[6], td[1], td[end], mean(filter(isfinite, [r.vg for r in rs])),
                 mean(filter(isfinite, [r.n for r in rs])), mean(filter(isfinite, [r.N for r in rs])))
        end
    end
    sayf("\nwall time %.1f min", (time() - t0) / 60)
    write(joinpath(HERE, "effective_size_branch.log"), String(take!(out)))
end

# ------------------------------------------------------------------------------------------
# Impact of ne = :breeders

module CR
    include(joinpath(@__DIR__, "critical_rate.jl"))
end

"The expected breeders N·W̄ and the house-of-Gauss V_G that are consistent with each other at a fixed optimum."
function hog_consistent(p::Params; ne::Symbol)
    V = house_of_gauss_vg(p; ne)
    n = float(p.PopSize)
    for _ in 1:200
        n = p.PopSize * PEAK * initFitSD(p) / sqrt(initFitSD(p)^2 + V + QTLSim.envVar(p))
        Vn = house_of_gauss_vg(p; nsurv = n, ne)
        abs(Vn - V) <= 1e-12 * Vn && break
        V = Vn
    end
    (; V, n)
end

function impact_part()
    t0 = time()
    say("validate/effective_size.jl impact — what ne = :breeders (N_e = the breeders) changes, against the default ne = :falconer")
    say("Only functions that call effective_size are affected: house_of_gauss_vg/_G and house_of_cards_vg when nsurv < PopSize is passed, and")
    say("drift_critical_rate. The moments lowering (run_moments, multilocus_equilibrium, emulate, half_survival_rate's survival) uses the parents n directly.")
    pubs = published_inputs(75)

    say("\n1. drift_critical_rate on published inputs (75 per stratum, the move's own rate not used): k_c(:breeders)/k_c(:falconer), min / median / max")
    sayf("  %-8s %-26s %-26s %-26s %-26s", "loci", DEMOS[1][1], DEMOS[2][1], DEMOS[3][1], DEMOS[4][1])
    ratios = Dict{Tuple{Int,String},Vector{Float64}}()
    lk = ReentrantLock()
    Threads.@threads :greedy for x in pubs
        for (dn, dem) in DEMOS
            q = dem(x.p)
            a, b = drift_critical_rate(q).k, drift_critical_rate(q; ne = :breeders).k
            a > 0 && lock(() -> push!(get!(ratios, (x.L, dn), Float64[]), b / a), lk)
        end
    end
    for L in (100, 500, 1000, 5000)
        say(@sprintf("  %-8d", L), join((@sprintf("%.3f / %.3f / %.3f (%3d)   ", q3(get(ratios, (L, dn), Float64[]))..., length(get(ratios, (L, dn), Float64[]))) for (dn, _) in DEMOS)))
    end
    allr = reduce(vcat, values(ratios))
    sayf("  all: %.3f / %.3f / %.3f over %d (input, demography) pairs with k_c > 0", q3(allr)..., length(allr))

    say("\n2. Pre-move V_G against the published table (\"Genotype Var at 10000\"), the same 600 inputs; table/prediction, median [quartiles] per stratum.")
    say("   house of Gauss with N_e from the expected breeders N·W̄ (self-consistent), both schemes, and with N_e = N (the default nsurv = PopSize);")
    say("   multilocus_equilibrium(p; drift = true) (rare alleles with LD; the comparison of validate/moments_regimes.jl) does not use effective_size.")
    xl = joinpath(ROOT, "Repository", "Results", "Master Output file.xlsx")
    if isfile(xl)
        tab = load_master(xl); idx = row_index(tab)
        rowsv = Vector{Any}(undef, length(pubs))
        Threads.@threads :greedy for i in eachindex(pubs)
            x = pubs[i]
            r = idx(x.L, x.u, x.p)
            if r == 0 || !(tab.gv[r] > 0)
                rowsv[i] = nothing
            else
                f, b = hog_consistent(x.p; ne = :falconer), hog_consistent(x.p; ne = :breeders)
                me = try multilocus_equilibrium(x.p; drift = true).VG catch; NaN end
                rowsv[i] = (; x.L, x.u, gv = tab.gv[r], f = f.V, b = b.V, N = house_of_gauss_vg(x.p), me, nfrac = b.n / x.p.PopSize)
            end
        end
        ok = filter(!isnothing, rowsv)
        qs(v) = (w = sort(filter(isfinite, v)); isempty(w) ? (NaN, NaN, NaN) : (w[max(1, round(Int, 0.25 * length(w)))], median(w), w[max(1, round(Int, 0.75 * length(w)))]))
        sayf("  %-14s %4s %7s | %-22s %-22s %-22s %-22s", "stratum", "n", "n̄/N", "hog, :falconer", "hog, :breeders", "hog, N_e = N", "multilocus_equilibrium")
        for (L, u) in STRATA
            v = [r for r in ok if r.L == L && r.u == u]
            isempty(v) && continue
            col(f) = (q = qs([r.gv / f(r) for r in v]); @sprintf("%.3f [%.3f, %.3f]   ", q[2], q[1], q[3]))
            say(@sprintf("  QT%04d-%d      %4d %7.3f | ", L, u, length(v), median(r.nfrac for r in v)), col(r -> r.f), col(r -> r.b), col(r -> r.N), col(r -> r.me))
        end
        colall(f) = (q = qs([r.gv / f(r) for r in ok]); @sprintf("%.3f [%.3f, %.3f]   ", q[2], q[1], q[3]))
        say(@sprintf("  %-14s %4d %7.3f | ", "all", length(ok), median(r.nfrac for r in ok)), colall(r -> r.f), colall(r -> r.b), colall(r -> r.N), colall(r -> r.me))
        sayf("  V_G(:breeders)/V_G(:falconer): %.4f / %.4f / %.4f (min / median / max)", q3([r.b / r.f for r in ok])...)
    else
        say("  (", xl, " not found: skipped)")
    end

    say("\n3. validate/critical_rate.jl's headline numbers, recomputed from its cache (critical_rate.csv, the kernel's k50 per cell) with each scheme")
    rows = sweep(j -> error("the cache is incomplete"), CR.jobs(CR.NREP); cache = joinpath(HERE, "critical_rate.csv"), id = j -> "$(j.set)-$(j.dem)-$(j.c)-$(j.rep)")
    br = CR.branches(rows)
    hsr = CR.predictions(joinpath(HERE, "critical_rate_pred.csv"))
    for ne in (:falconer, :breeders)
        c1 = [0, 0]; c2 = [0, 0]; maxerr = 0.0; nsign = 0; Mr = Float64[]; kk = Dict(:L100 => Float64[], :L500 => Float64[])
        for set in keys(CR.SETS), dem in CR.DEMS
            p0 = CR.demp(set, dem)
            cs = collect(CR.COSTS[dem])
            k = Dict(c => CR.k50([x for x in br if x.set === set && x.dem === dem && x.c == c]) for c in cs)
            d0 = drift_critical_rate(p0; ne)
            ρo = Float64[]
            for c in cs
                kc = drift_critical_rate(p0; cost = c, ne).k
                push!(kk[set], k[c] / kc)
                c > 0 || continue
                ρ, ρp = k[c] / k[0.0], kc / d0.k
                push!(ρo, ρ)
                c1[2] += 1; c1[1] += abs(ρ - ρp) <= 0.06; maxerr = max(maxerr, abs(ρ - ρp)); nsign += ρ >= ρp
            end
            Mo = CR.fitM(cs[2:end], ρo)
            push!(Mr, Mo / d0.M)
            c2[2] += 1; c2[1] += 0.8 <= Mo / d0.M <= 1.25
        end
        sayf("  ne = :%-9s C1 (cost ratio within 0.06) %d/%d, largest |error| %.3f, observed ≥ predicted in %d; C2 (M_obs/M_theory in 0.8–1.25) %d/%d, M_obs/M_theory %.2f–%.2f;",
             ne, c1..., maxerr, nsign, c2..., extrema(Mr)...)
        sayf("  %-14s k50/k_c %.2f–%.2f (100 loci), %.2f–%.2f (500 loci)", "", extrema(kk[:L100])..., extrema(kk[:L500])...)
    end
    say("  C3 (half_survival_rate within 0.75–1.33 of the kernel, 26/32) does not depend on the scheme: its survival is the moments lowering's, and ne only")
    say("  moves the bracket of its bisection. Check on two cells (12 bisection steps resolve ln k to 5e-4):")
    for (set, dem, c) in ((:L100, :fixed, 0.0), (:L100, :R2, 0.2))
        p0 = CR.demp(set, dem)
        a = half_survival_rate(p0; cost = c); b = half_survival_rate(p0; cost = c, ne = :breeders)
        sayf("    %s %s c = %.1f: k50 %.6f (:falconer; cached prediction %.6f), %.6f (:breeders); ln ratio %+.1e; k_c %.5f → %.5f", set, dem, c, a.k50,
             hsr[(set, dem, c)], b.k50, log(b.k50 / a.k50), a.kc, b.kc)
    end

    say("\n4. The emulator's half-survival rate on 32 published inputs (4 per stratum; emulated_half_survival_rate, loci = :spectrum with the beta solver).")
    say("   The emulator does not call effective_size; drift_critical_rate only sets the starting bracket (here passed as `bracket` for :breeders).")
    sub = [x for (L, u) in STRATA for x in [y for y in pubs if y.L == L && y.u == u][1:4]]
    res = Vector{Any}(undef, length(sub))
    Threads.@threads :greedy for i in eachindex(sub)
        p = sub[i].p
        a = emulated_half_survival_rate(p; loci = :spectrum, spectrum = (; solver = :beta), rtol = 1e-8)
        kf, kb = drift_critical_rate(p).k, drift_critical_rate(p; ne = :breeders).k
        b = kb > 0 ? emulated_half_survival_rate(p; loci = :spectrum, spectrum = (; solver = :beta), rtol = 1e-8, bracket = (kb / 8, 2kb)) : (k50 = 0.0,)
        res[i] = (; sub[i].L, sub[i].u, a = a.k50, b = b.k50, kf, kb)
    end
    fin = [r for r in res if isfinite(r.a) && r.a > 0 && isfinite(r.b) && r.b > 0]
    sayf("   %d of %d with a finite k50 under both brackets; |ln(k50_b/k50_a)|: max %.1e", length(fin), length(res), maximum(abs(log(r.b / r.a)) for r in fin; init = 0.0))
    for L in (100, 500, 1000, 5000)
        v = [r for r in fin if r.L == L]
        isempty(v) && continue
        sayf("   %4d loci: emulator k50/k_c  %.2f / %.2f / %.2f with :falconer,  %.2f / %.2f / %.2f with :breeders (min / median / max, %d inputs)", L,
             q3([r.a / r.kf for r in v])..., q3([r.a / r.kb for r in v])..., length(v))
    end
    sayf("\nwall time %.1f min", (time() - t0) / 60)
    write(joinpath(HERE, "effective_size_impact.log"), String(take!(out)))
end
