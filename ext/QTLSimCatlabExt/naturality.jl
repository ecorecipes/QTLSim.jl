# The naturality harness (roadmap C2; the data types and the comparisons are in
# src/naturality.jl): the runners of the backends of a square, `check_naturality`, and the
# per-step squares of one generation, `naturality_steps`.

using QTLSim: NaturalityCase, NaturalityResult, Statistic, statvalue, compare_stat, check_case, STOCHASTIC_BACKENDS
import QTLSim: check_naturality, naturality_steps

# --- runners: the output of one run of a backend ----------------------------------------

case_params(c::NaturalityCase) =
    c.horizon === nothing ? c.params :
    c.params isa Params ? Params(c.params; NoGenerations = c.horizon) : [Params(q; NoGenerations = c.horizon) for q in c.params]
case_ndemes(c::NaturalityCase) = c.params isa Params ? (c.cycle === nothing ? 1 : nfactors(dom(c.cycle))) : length(c.params)
case_cycle(c::NaturalityCase) = c.cycle === nothing ? lifecycle().generation : c.cycle
option(c::NaturalityCase, k, default = (;)) = get(c.options, k, default)

# The kernel: QTSim's columns (`records`) plus each deme's mean genotype at the census
# (`geno_mean`, which the variance among deme means needs) and a `deme` column, recorded
# by the image of `record`.
function run_backend(::Val{:kernel}, c::NaturalityCase, seed)
    p, k, cyc = case_params(c), case_ndemes(c), case_cycle(c)
    rate = option(c, :migration_rate, 0.0)
    p1 = p isa Params ? p : first(p)
    gm = Float64[]
    obs(rec) = operator(:record, function (w, rng)
        if w.gen % w.p.Sample == 0
            record!(rec, w)
            push!(gm, mean(view(w.genotype, 1:w.N)))
        end
        w
    end)
    if k == 1
        w = World(p isa Params ? p : only(p))
        rec = records(; traits = p1.NoTraits)
        o = run_composed(w, Xoshiro(seed); cycle = cyc, images = Dict(:record => obs(rec)), record = false, option(c, :kernel)...)
        return (; extinct = o.extinct, extinct_gen = o.extinct_gen, rec..., deme = fill(1, length(rec.gen)), geno_mean = gm)
    end
    s = (p isa Params ? p : collect(Params, p)) ⊗ demes(k)
    rec = records(; demes = s.demes, traits = p1.NoTraits)
    o = run_composed(s, Xoshiro(seed); cycle = cyc, images = Dict(:record => obs(rec)), migration_rate = rate,
                     record = false, option(c, :kernel)...)
    cols = NamedTuple(n => getproperty(rec, n) for n in keys(rec) if n !== :demeof)
    (; extinct = o.extinct, extinct_gen = o.extinct_gen, cols..., geno_mean = gm)
end

function run_backend(::Val{:moments}, c::NaturalityCase, seed)
    p, cyc = case_params(c), case_cycle(c)
    kw = option(c, :moments)
    rate = case_ndemes(c) > 1 ? (migration_rate = option(c, :migration_rate, 0.0),) : (;)
    p isa Params ? run_moments(p; cycle = cyc, rate..., kw...) : run_moments(collect(Params, p); cycle = cyc, rate..., kw...)
end

run_backend(::Val{:theory}, c::NaturalityCase, seed) = option(c, :theory, nothing)(c)

"The SLiM executable for the harness: `ENV[\"SLIM\"]`, `slim/build/slim` beside the package, or `slim` on the path."
function harness_slim()
    for path in (get(ENV, "SLIM", ""), joinpath(dirname(dirname(pathof(QTLSim))), "..", "slim", "build", "slim"),
                 something(Sys.which("slim"), ""))
        isfile(path) && return path
    end
    nothing
end

const SLIM_COLUMNS = Dict("gen" => :gen, "phenoMean" => :pheno_mean, "genoVar" => :geno_var, "nSurv" => :nsurv,
                          "theta" => :theta, "omega" => :omega, "deme" => :deme, "popsize" => :popsize)

# SLiM's table: QTSim's columns, as the generated model prints them (eidos.jl's `record`)
function run_backend(::Val{:slim}, c::NaturalityCase, seed; script)
    slim = get(option(c, :slim), :slim, harness_slim())
    slim === nothing && throw(ArgumentError("SLiM not found: set ENV[\"SLIM\"] or options.slim.slim"))
    out = read(`$slim -s $seed $script`, String)
    lines = split(out, '\n')
    h = findfirst(l -> startswith(l, "gen\t"), lines)
    h === nothing && error("no table in SLiM's output:\n", first(out, 2000))
    names = [get(SLIM_COLUMNS, n, Symbol(n)) for n in split(lines[h], '\t')]
    rows = [split(l, '\t') for l in lines[h+1:end] if occursin(r"^\d+\t", l)]
    # the last column, the population size, is printed only under density regulation
    cols = Dict{Symbol,Vector}()
    for (j, n) in enumerate(names)
        cols[n] = [n in (:gen, :nsurv, :deme, :popsize) ? parse(Int, r[j]) : parse(Float64, r[j]) for r in rows if length(r) >= j]
    end
    (; (n => cols[n] for n in names)..., extinct = occursin("Pop all died.", out))
end

function slim_file(c::NaturalityCase)
    p, cyc = case_params(c), case_cycle(c)
    path = tempname() * ".slim"
    write(path, slim_script(cyc, p isa Params ? p : collect(Params, p); migration_rate = option(c, :migration_rate, 0.0)))
    path
end

# the outputs of `backend` for the case: one per replicate, or one for a deterministic backend
function backend_outputs(backend::Symbol, c::NaturalityCase, seed, ntasks)
    backend in STOCHASTIC_BACKENDS || return Any[run_backend(Val(backend), c, seed)]
    script = backend === :slim ? slim_file(c) : nothing
    outs = Vector{Any}(undef, c.replicates)
    next = Threads.Atomic{Int}(1)
    work() = while true
        r = Threads.atomic_add!(next, 1)
        r > c.replicates && break
        outs[r] = backend === :slim ? run_backend(Val(:slim), c, seed + r; script) : run_backend(Val(backend), c, seed + r)
    end
    tasks = [Threads.@spawn work() for _ in 1:max(1, min(ntasks, c.replicates))]
    foreach(wait, tasks)
    script === nothing || rm(script; force = true)
    outs
end

function check_naturality(c::NaturalityCase; seed::Integer = c.seed, ntasks::Integer = Threads.nthreads(), verbose::Bool = false)
    check_case(c)
    t0 = time()
    fouts = backend_outputs(c.from, c, seed, ntasks)
    touts = backend_outputs(c.to, c, seed + 1_000_000, ntasks)
    rows = NamedTuple[compare_stat(c.criterion, s, tol, c.from, c.to, fouts, touts) for (s, tol) in c.statistics]
    ok = all(r -> r.pass, rows)
    r = NaturalityResult(c, rows, ok, ok == (c.expected === :pass), Dict{Symbol,Any}(c.from => fouts, c.to => touts), time() - t0)
    verbose && show(stdout, MIME"text/plain"(), r)
    r
end

function check_naturality(cycle, p; from::Symbol = :kernel, to::Symbol = :moments, statistic, tolerance,
                          criterion::Symbol = :band, replicates::Integer = 20, generations = nothing, seed::Integer = 1,
                          options::NamedTuple = (;), name::AbstractString = "", expected::Symbol = :pass, kw...)
    stats = statistic isa Statistic ? [statistic => Float64(tolerance)] :
            [s => Float64(t) for (s, t) in zip(statistic, tolerance isa Real ? fill(tolerance, length(statistic)) : tolerance)]
    c = NaturalityCase(; name = isempty(name) ? "$(from) → $(to)" : String(name), cycle, params = p, from, to, criterion,
                       statistics = stats, replicates, horizon = generations, options, expected, seed)
    check_naturality(c; kw...)
end

# --- the square of each step of one generation ------------------------------------------

"The generators of a composite, in the order they apply (`f ⋅ g`: f first)."
stepsequence(f) = head(f) === :compose ? reduce(vcat, (stepsequence(a) for a in args(f))) : Any[f]

# the genotypes of the first `n` individuals of the chromosomes `chr` (block ids) of `w`
function chr_genotypes(w::World, chr, n)
    sums = w.pool.sums
    m = size(chr, 1) * size(chr, 2)
    [sum(sums[chr[t]] for t in m*(i-1)+1:m*i) for i in 1:n]
end

# the kernel's observables at the codomain `stage` of a step: the individuals there and their genotypes
function kernel_observables(w::World, stage::Symbol)
    g = stage === :Zygotes ? chr_genotypes(w, w.chr, w.N) :
        stage === :Adults ? w.genotype[1:w.N] :
        stage === :Survivors ? w.genotype[w.survivors] :
        chr_genotypes(w, w.schr, w.Noff)            # Offspring
    (mean = isempty(g) ? NaN : mean(g), VG = length(g) > 1 ? var(g) : NaN, N = float(length(g)))
end

# the moments' observables at `stage` (Moments, one deme)
function moment_observables(m::Moments, stage::Symbol)
    vc() = cov_vpop(m, m.cov)
    stage in (:Zygotes, :Adults) ? (mean = m.zbar, VG = m.VG, N = m.N) :
    stage === :Survivors ? (mean = m.zs, VG = m.cov === nothing ? m.VGs : vc(), N = m.nsurv) :
                           (mean = m.zoff, VG = m.cov === nothing ? m.VGoff : vc(), N = m.Noff)
end

# shuffle each locus's alleles among all chromosomes: no LD, Hardy–Weinberg proportions
function linkage_equilibrium!(w::World, rng)
    G = genomes(w)
    L, P, N = size(G)
    for l in 1:L
        v = vec(G[l, :, :])
        G[l, :, :] = reshape(v[randperm(rng, length(v))], P, N)
    end
    setgenomes!(w, G)
    w
end

function naturality_steps(p::Params; cycle = lifecycle().generation, loci::Symbol = :multilocus, erosion::Bool = false,
                          replicates::Integer = 200, world = nothing, burnin::Integer = 1000, linkage_equilibrium::Bool = true,
                          tolerance::Real = 0.05, seed::Integer = 1)
    p.NoTraits == 1 || throw(ArgumentError("naturality_steps is for one trait; NoTraits = $(p.NoTraits)"))
    dom(cycle) == LifeCycle[:Zygotes] && codom(cycle) == LifeCycle[:Zygotes] ||
        throw(ArgumentError("naturality_steps checks a generation of one deme, Zygotes → Zygotes"))
    loci in (:aggregate, :multilocus, :rare_alleles) ||
        throw(ArgumentError("naturality_steps builds the moments from a mean and a variance: loci is :aggregate, :multilocus or :rare_alleles"))
    rng = Xoshiro(seed)
    w0 = world === nothing ? (w = World(Params(p; NoGenerations = burnin)); run!(w, rng; record = false); w) : deepcopy(world)
    w0.extinct && throw(ArgumentError("the starting population died"))
    linkage_equilibrium && linkage_equilibrium!(w0, rng)
    # σ(w0): the moments state with the population's mean genotype and genetic variance
    z0 = chr_genotypes(w0, w0.chr, w0.N)
    loci2, α2 = resolve_loci([p], cycle, loci)
    cov = loci2 === :aggregate ? nothing : DemeCov(p, [float(w0.N)], [var(z0)]; α2, drift = erosion)
    m = Moments(p; zbar = mean(z0), VG = var(z0), erosion, cov)
    m.θ, m.ω, m.gen, m.N = w0.θ, w0.ω, w0.gen, w0.N
    model = model_of(p)
    ws = [deepcopy(w0) for _ in 1:replicates]
    rngs = [Xoshiro(seed + r) for r in 1:replicates]
    rows = NamedTuple[]
    for g in stepsequence(cycle)
        kop = lower(g, :kernel; model)
        mop = lower(g, :moments; model)
        stage = Symbol(head(codom(g)) === :generator ? first(args(codom(g))) : codom(g))
        for (w, r) in zip(ws, rngs)
            w.extinct || kop(w, r)
        end
        m.extinct || mop(m, nothing)
        ko = [kernel_observables(w, stage) for w in ws if !w.extinct]
        mo = moment_observables(m, stage)
        for obs in (:mean, :VG, :N)
            x = [getproperty(o, obs) for o in ko]
            k, se = mean(x), std(x) / sqrt(length(x))
            pred = getproperty(mo, obs)
            scale = obs === :mean ? sqrt(max(mo.VG, 0.0)) : abs(pred)
            measure = (k - pred) / scale
            push!(rows, (step = stepname(g), stage, observable = obs, kernel = k, se, moments = pred, measure,
                         pass = isfinite(k - pred) && !(abs(k - pred) - 1.96se > tolerance * scale)))
        end
    end
    rows
end
