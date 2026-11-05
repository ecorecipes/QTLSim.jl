# GSL-locked replica of QTSim (diploid, any NoQTL ≥ 1, all optimum-move types).
#
# This is a *validation instrument*, NOT the production model. It draws ALL
# randomness from the same GSL routines QTSim uses (gsl_rng_mt19937,
# gsl_ran_gaussian, gsl_ran_poisson, gsl_rng_uniform / _int), in QTSim's EXACT
# draw order, and computes fitness and summary statistics with the same GSL
# functions (gsl_ran_gaussian_pdf, gsl_stats_mean / _variance). Seeded
# identically to QTSim it should therefore reproduce QTSim's trajectory
# bit-for-bit. Comparing the two removes the RNG as a confound and tests whether
# the port's ALGORITHM is identical to the reference — a far stronger check than
# distributional agreement.
#
# Two regimes are checked:
#   * the single-locus 6500-generation protocol written by
#     Repository/Code/QTSim/Scrips/CreateCGenomeSimInputDist.py, and
#   * the published sweep (Repository/Simulations/QT<NoQTL>-Possion<u>/): 100–5000
#     loci, Poisson recombination, a continuously moving optimum from generation
#     10000. Shipped input files are run as-is except for a generation cap.
#
# QTSim's recombination (CreateCombMap + RecombinationPoissonChr) draws
# Poisson(RecombinationPoissonMean) distinct crossover points but never advances
# its point index, so only the LOWEST point switches strands: every gamete has at
# most one crossover. The replica reproduces that exactly (see `gamete_qtsim!`);
# the production kernel does the same by default (Params `CrossoverModel = :qtsim`).
#
# The production kernel (src/model.jl) is deliberately an independent, RNG-
# agnostic implementation; this replica exists only to prove algorithmic parity.
#
# Run from the package root:  julia --project=. validate/gsl_locked.jl
# Needs the QTSim binary (see docs/README.md) and a GSL shared library; set
# ENV["LIBGSL"] if libgsl is not in a Homebrew location.

using CategoricalQuantGen, Printf

const ROOT   = normpath(joinpath(@__DIR__, "..", ".."))
const QTSIM  = joinpath(ROOT, "Repository", "Code", "QTSim", "src", "QTSim")
const SIMDIR = joinpath(ROOT, "Repository", "Simulations")
const LIBGSL = let c = filter(isfile, [get(ENV, "LIBGSL", ""),
                                        "/opt/homebrew/opt/gsl/lib/libgsl.dylib",
                                        "/opt/homebrew/lib/libgsl.dylib",
                                        "/usr/local/lib/libgsl.dylib",
                                        "/usr/lib/libgsl.dylib",
                                        "/usr/lib/x86_64-linux-gnu/libgsl.so"])
    isempty(c) ? "" : first(c)
end

# --- GSL bindings (same functions QTSim links against) -------------------
gsl_rng_alloc() = begin
    T = unsafe_load(cglobal((:gsl_rng_mt19937, LIBGSL), Ptr{Cvoid}))
    ccall((:gsl_rng_alloc, LIBGSL), Ptr{Cvoid}, (Ptr{Cvoid},), T)
end
gset(r, s)        = ccall((:gsl_rng_set, LIBGSL), Cvoid, (Ptr{Cvoid}, Culong), r, s % Culong)
gfree(r)          = ccall((:gsl_rng_free, LIBGSL), Cvoid, (Ptr{Cvoid},), r)
guni(r)           = ccall((:gsl_rng_uniform, LIBGSL), Cdouble, (Ptr{Cvoid},), r)
guni_int(r, n)    = ccall((:gsl_rng_uniform_int, LIBGSL), Culong, (Ptr{Cvoid}, Culong), r, n % Culong)
ggauss(r, σ)      = ccall((:gsl_ran_gaussian, LIBGSL), Cdouble, (Ptr{Cvoid}, Cdouble), r, σ)
gpois(r, μ)       = ccall((:gsl_ran_poisson, LIBGSL), Cuint, (Ptr{Cvoid}, Cdouble), r, μ)
ggauss_pdf(x, σ)  = ccall((:gsl_ran_gaussian_pdf, LIBGSL), Cdouble, (Cdouble, Cdouble), x, σ)
gmean(v)          = ccall((:gsl_stats_mean, LIBGSL), Cdouble, (Ptr{Cdouble}, Csize_t, Csize_t), v, 1, length(v))
gvar(v)           = ccall((:gsl_stats_variance, LIBGSL), Cdouble, (Ptr{Cdouble}, Csize_t, Csize_t), v, 1, length(v))

# --- recombination, exactly as QTSim's CreateCombMap + RecombinationPoissonChr ---
# Draw NoPoints ~ Poisson(mean). If 0 or ≥ L the map is the sentinel {L} (no
# switch). Otherwise draw NoPoints DISTINCT points in 0:L-2 (rejection, in draw
# order), sort, and switch strands after the first (lowest) point only.
# `allpoints=true` switches at every point; it is used only as a negative control.
function gamete_qtsim!(dest, P1, P2, r, recmean, L, mapbuf, allpoints::Bool=false)
    k = Int(gpois(r, recmean))
    nsw = 0                                   # sentinel Map[0] = NoQTL ⇒ no switch
    if !(k == 0 || k >= L)
        for idx in 1:k
            x = 0
            while true
                x = Int(guni_int(r, L - 1))
                dup = false
                for q in 1:idx-1
                    if mapbuf[q] == x; dup = true; break; end
                end
                dup || break
            end
            mapbuf[idx] = x
        end
        sort!(@view mapbuf[1:k])              # qsort ascending
        nsw = allpoints ? k : 1               # QTSim never advances PointIndex
    end
    CP, OP = P1, P2
    if guni(r) < 0.5; CP, OP = P2, P1; end
    pi = 1
    @inbounds for i in 0:L-1
        dest[i+1] = CP[i+1]
        if pi <= nsw && i == mapbuf[pi]; CP, OP = OP, CP; pi += 1; end
    end
    dest
end

# --- the replica ---------------------------------------------------------
function run_gsl_locked(p::Params, seed::Integer; allpoints::Bool=false)
    p.Ploidy == 2 || error("gsl_locked replica covers Ploidy = 2")
    (p.StandardisePhenotype || p.Drift) && error("StandardisePhenotype/Drift not replicated")
    N, L = p.PopSize, p.NoQTL
    vE = envVar(p)
    envdraw(r) = vE != 0 ? ggauss(r, sqrt(vE)) : 0.0          # SetIndividualEnv
    # NLOptParam / GetAnalyticalFitSD, in QTSim's exact floating-point form
    ω = p.InitFitnessSD >= 0 ? p.InitFitnessSD :
        (p.PhenotypeVar * p.Heritability) / sqrt(p.MutationScalar * (p.PhenotypeVar * (1.0 - p.Heritability)))
    θ = 0.0
    rate = p.MutationRatePerGamete
    mutsd = sqrt((p.MutationScalar / (p.Ploidy * L * rate)) * vE)
    λmut = L * rate

    r = gsl_rng_alloc(); gset(r, seed)
    # CreatePop is called twice (CPop then NPop): 2N environment draws at init.
    G0 = zeros(Cdouble, L, N); G1 = zeros(Cdouble, L, N); env  = Cdouble[envdraw(r) for _ in 1:N]
    H0 = zeros(Cdouble, L, N); H1 = zeros(Cdouble, L, N); nenv = Cdouble[envdraw(r) for _ in 1:N]
    pheno = zeros(Cdouble, N); geno = zeros(Cdouble, N); fit = zeros(Cdouble, N)
    surv = Vector{Int}(undef, N); mapbuf = Vector{Int}(undef, max(L, 1))
    og = Int[]; opm = Float64[]; ogv = Float64[]; ons = Int[]; oθ = Float64[]

    for it in 0:p.NoGenerations-1
        @inbounds for i in 1:N                                  # CaclGenotype order
            s = 0.0
            for l in 1:L; s += G0[l, i] + G1[l, i]; end
            geno[i] = s
            pheno[i] = s + env[i]
        end
        if p.MoveFitnessSDGen != -1 && it >= p.MoveFitnessSDGen  # SetFitnessMean
            if p.MoveType == :fitness
                θ = gmean(pheno) + ω * p.MoveFitnessSD
            elseif p.MoveType == :phenotype
                θ = gmean(pheno) + sqrt(p.PhenotypeVar) * p.MoveFitnessSD
            else
                θ = θ + p.MoveFitnessSD
            end
        end
        @inbounds for i in 1:N
            fit[i] = ggauss_pdf(pheno[i] - θ, ω) / (0.4 / ω)
        end
        ns = 0                                                  # SetSurviving
        @inbounds for i in 1:N
            if guni(r) < fit[i]; ns += 1; surv[ns] = i; end
        end
        if ns < 2                                               # DeadPopulation
            gfree(r)
            return (extinct=true, extinct_gen=it, gen=og, pheno_mean=opm, geno_var=ogv, nsurv=ons, theta=oθ)
        end
        @inbounds for j in 1:N                                  # NewIn + SetIndividualEnv
            p1 = p2 = 0
            while true                                          # QTSim redraws BOTH parents
                p1 = surv[Int(guni_int(r, ns)) + 1]
                p2 = surv[Int(guni_int(r, ns)) + 1]
                p1 != p2 && break
            end
            gamete_qtsim!(@view(H0[:, j]), @view(G0[:, p1]), @view(G1[:, p1]), r, p.RecombPoissonMean, L, mapbuf, allpoints)
            gamete_qtsim!(@view(H1[:, j]), @view(G0[:, p2]), @view(G1[:, p2]), r, p.RecombPoissonMean, L, mapbuf, allpoints)
            nenv[j] = envdraw(r)
        end
        if it % p.Sample == 0                                   # Output(CPop)
            push!(og, it); push!(opm, gmean(pheno)); push!(ogv, gvar(geno))
            push!(ons, ns); push!(oθ, θ)
        end
        it == p.ChangeFitnessGen && (ω = p.ChangeFitnessSD)
        @inbounds for j in 1:N                                  # MutatePop(NPop): chr0 then chr1
            for H in (H0, H1)
                no = gpois(r, λmut)
                while no > 0
                    pos = Int(guni_int(r, L)) + 1
                    H[pos, j] += ggauss(r, mutsd)
                    no -= 1
                end
            end
        end
        G0, H0 = H0, G0; G1, H1 = H1, G1; env, nenv = nenv, env
    end
    gfree(r)
    (extinct=false, extinct_gen=-1, gen=og, pheno_mean=opm, geno_var=ogv, nsurv=ons, theta=oθ)
end

# --- QTSim option files ---------------------------------------------------
"""
    read_qtsim_opts(path) -> (Params, seed)

Parse a QTSim options file (one `Key value…` per line; `#` comments) into
`Params` plus the seed.
"""
function read_qtsim_opts(path)
    kw = Dict{Symbol,Any}(); seed = 0
    for ln in eachline(path)
        s = strip(first(split(ln, '#'))); isempty(s) && continue
        f = split(s); k = lowercase(f[1])
        if     k == "sample";                   kw[:Sample] = parse(Int, f[2])
        elseif k == "popsize";                  kw[:PopSize] = parse(Int, f[2])
        elseif k == "noqtl";                    kw[:NoQTL] = parse(Int, f[2])
        elseif k == "generations";              kw[:NoGenerations] = parse(Int, f[2])
        elseif k == "ploidy";                   kw[:Ploidy] = parse(Int, f[2])
        elseif k == "movefitnesssd";            kw[:MoveFitnessSDGen] = parse(Int, f[2]); kw[:MoveFitnessSD] = parse(Float64, f[3])
        elseif k == "changefitnesssd";          kw[:ChangeFitnessGen] = parse(Int, f[2]); kw[:ChangeFitnessSD] = parse(Float64, f[3])
        elseif k == "phenotypevar";             kw[:PhenotypeVar] = parse(Float64, f[2])
        elseif k == "mutationscalar";           kw[:MutationScalar] = parse(Float64, f[2])
        elseif k == "heritability";             kw[:Heritability] = parse(Float64, f[2])
        elseif k == "initfitnesssd";            kw[:InitFitnessSD] = parse(Float64, f[2])
        elseif k == "optimummovetype";          kw[:MoveType] = Symbol(lowercase(f[2]))
        elseif k == "mutationratepergamete";    kw[:MutationRatePerGamete] = parse(Float64, f[2])
        elseif k == "recombinationpoissonmean"; kw[:RecombPoissonMean] = parse(Float64, f[2])
        elseif k == "seed";                     seed = parse(Int, f[2])
        else error("read_qtsim_opts: unsupported option '$(f[1])' in $path")
        end
    end
    Params(; kw...), seed
end

function write_opts(path, q::Params, seed::Integer)
    open(path, "w") do io
        for (k, v) in (("Sample", q.Sample), ("PopSize", q.PopSize), ("NoQTL", q.NoQTL),
                       ("Generations", q.NoGenerations), ("Ploidy", q.Ploidy))
            println(io, "$k $v")
        end
        q.MoveFitnessSDGen != -1 && println(io, "MoveFitnessSD $(q.MoveFitnessSDGen) $(q.MoveFitnessSD)")
        q.ChangeFitnessGen  != -1 && println(io, "ChangeFitnessSD $(q.ChangeFitnessGen) $(q.ChangeFitnessSD)")
        println(io, "OptimumMoveType $(q.MoveType)")
        println(io, "PhenotypeVar $(q.PhenotypeVar)")
        println(io, "MutationScalar $(q.MutationScalar)")
        println(io, "Heritability $(q.Heritability)")
        println(io, "MutationRatePerGamete $(q.MutationRatePerGamete)")
        println(io, "RecombinationPoissonMean $(q.RecombPoissonMean)")
        q.InitFitnessSD >= 0 && println(io, "InitFitnessSD $(q.InitFitnessSD)")
        println(io, "Seed $seed")
    end
end

"Run the QTSim binary; parse generation, phenotype mean, genotype variance, survivors, θ."
function run_qtsim_full(q::Params, seed::Integer)
    o = tempname() * ".txt"; write_opts(o, q, seed)
    s = try; read(`$QTSIM $o`, String); finally; rm(o; force=true); end
    g = Int[]; pm = Float64[]; gv = Float64[]; ns = Int[]; th = Float64[]; b = false; ext = false
    for ln in split(s, '\n')
        startswith(ln, "Pop all died") && (ext = true; break)
        startswith(ln, "Generation\t") && (b = true; continue)
        (b && !isempty(strip(ln))) || continue
        f = split(ln, '\t'); length(f) < 6 && continue
        push!(g, parse(Int, f[1])); push!(pm, parse(Float64, f[2])); push!(gv, parse(Float64, f[4]))
        push!(ns, parse(Int, f[5])); push!(th, parse(Float64, f[6]))
    end
    # QTSim exits inside NewPop, before printing the dying generation
    (gen=g, pheno_mean=pm, geno_var=gv, nsurv=ns, theta=th, extinct=ext,
     extinct_gen=ext ? (isempty(g) ? 0 : last(g) + q.Sample) : -1)
end

# --- compare (QTSim prints 6 decimals; bit-identity ⇒ 6-decimal string match) --
f6(x) = @sprintf("%f", x)
function compare(p::Params, seed::Integer, label; allpoints::Bool=false)
    gl = run_gsl_locked(p, seed; allpoints)
    qt = run_qtsim_full(p, seed)
    idx = Dict(g => i for (i, g) in enumerate(gl.gen))
    n = bad = 0; firstbad = -1
    for (k, g) in enumerate(qt.gen)
        i = get(idx, g, 0); i == 0 && continue
        n += 1
        same = f6(gl.pheno_mean[i]) == f6(qt.pheno_mean[k]) && f6(gl.geno_var[i]) == f6(qt.geno_var[k]) &&
               gl.nsurv[i] == qt.nsurv[k] && f6(gl.theta[i]) == f6(qt.theta[k])
        same || (bad += 1; firstbad < 0 && (firstbad = g))
    end
    sameend = gl.extinct == qt.extinct && (!gl.extinct || p.Sample != 1 || gl.extinct_gen == qt.extinct_gen)
    ok = bad == 0 && n == length(qt.gen) && sameend
    endtxt = qt.extinct ? "extinct at gen $(qt.extinct_gen) (replica $(gl.extinct_gen))" : "survived (replica extinct=$(gl.extinct))"
    @printf("%-44s gens compared=%-6d mismatches=%-4d %s\n  ⇒ %s\n", label, n, bad, endtxt,
            ok ? "BIT-IDENTICAL (6 dp + survivor counts, every generation)" :
                 (firstbad >= 0 ? "DIVERGES at gen $firstbad" : "end state differs"))
    ok
end

function main()
    isfile(LIBGSL) || error("libgsl not found; set ENV[\"LIBGSL\"]")
    isfile(QTSIM)  || error("QTSim binary not found at $QTSIM — build it (see docs/README.md)")
    allok = true
    println("GSL-locked replica vs QTSim (same seed, same RNG stream & order)\n")
    println("── single-locus protocol (CreateCGenomeSimInputDist.py) ──")
    A = Params(Sample=10, PopSize=1000, Ploidy=2, NoQTL=1, NoGenerations=4000,
               PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.02,
               MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    B = Params(Sample=1, PopSize=1000, Ploidy=2, NoQTL=1, NoGenerations=6500,
               PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001,
               ChangeFitnessGen=5000, ChangeFitnessSD=0.30,
               MoveFitnessSDGen=5000, MoveFitnessSD=1.0, MoveType=:phenotype)
    for s in (12345, 67890, 111213)
        allok &= compare(B, s, "directional (B) seed $s")
    end
    for s in (2222, 3333)
        allok &= compare(A, s, "equilibrium (A) seed $s")
    end
    println("\n── published sweep: shipped input files, capped at 10600 generations ──")
    for (dir, file) in (("QT0100-Possion3.0", "Input-00000.txt"), ("QT0100-Possion3.0", "Input-00001.txt"),
                        ("QT0100-Possion5.0", "Input-00002.txt"), ("QT0500-Possion3.0", "Input-00003.txt"),
                        ("QT0500-Possion5.0", "Input-00004.txt"), ("QT1000-Possion3.0", "Input-00005.txt"))
        q, seed = read_qtsim_opts(joinpath(SIMDIR, dir, file))
        kw = Dict{Symbol,Any}(k => getfield(q, k) for k in fieldnames(Params))
        kw[:NoGenerations] = min(q.NoGenerations, 10600)
        q = Params(; kw...)
        allok &= compare(q, seed, "$dir/$file (L=$(q.NoQTL), u=$(q.RecombPoissonMean))")
    end
    println("\n── negative control: replica applying EVERY crossover point (not QTSim's behaviour) ──")
    q, seed = read_qtsim_opts(joinpath(SIMDIR, "QT0100-Possion3.0", "Input-00000.txt"))
    kw = Dict{Symbol,Any}(k => getfield(q, k) for k in fieldnames(Params)); kw[:NoGenerations] = 2000
    detected = !compare(Params(; kw...), seed, "QT0100-Possion3.0/Input-00000.txt, all points"; allpoints=true)
    println(detected ? "  ⇒ divergence detected, as expected: the check is sensitive to the crossover model" :
                       "  ⇒ NOT detected: the check cannot distinguish the crossover models")
    allok &= detected
    println(allok ? "\nALL CASES BIT-IDENTICAL" : "\nSOME CASES DIFFER")
    allok
end

if abspath(PROGRAM_FILE) == @__FILE__
    main() || exit(1)
end
