# Sweeps and replicates: a function mapped over jobs on all threads, with an optional CSV
# cache from which an interrupted sweep resumes. Standard library only.

"""
    sweep(f, jobs; cache = nothing, id = identity, cost = nothing,
          ntasks = Threads.nthreads(), progress = 0)

Apply `f` to every job in `jobs` and return the results in the order of `jobs`.

Jobs are handed out one at a time to `ntasks` tasks on Julia's threads, so a long job
never holds up a queue of others (the greedy scheduling of `Threads.@threads :greedy`,
which needs Julia 1.11; `:dynamic` gives each thread a fixed block of jobs). With `cost`,
a function of a job, the costliest jobs start first, which shortens a sweep whose jobs
differ in cost. Neither changes the results, provided `f` takes its randomness from the
job (a seed) rather than from a shared generator: the results then do not depend on the
number of threads or on the order in which the jobs run.

With `cache`, the path of a CSV file, `f` must return a `NamedTuple` of numbers, `Bool`s,
strings or symbols, with the same fields for every job. Each row is appended to the file
as soon as it is computed, after a first column `id` holding `string(id(job))`, and jobs
whose id is already in the file are not run again, so an interrupted sweep resumes where
it stopped. The rows returned are those of the file, for jobs run now as for jobs run
before, so a resumed sweep returns the same rows as an uninterrupted one. Each column
comes back as the narrowest of `Bool`, `Int`, `Float64` and `String` that reads all its
values (`Float64`s round-trip exactly; symbols, and strings that read as numbers, do
not). A line cut short by an interruption is dropped. Two jobs with the same id, and a
file whose columns differ from the rows, are errors; rows of ids not among `jobs` are
kept in the file and ignored.

If `f` throws, no further jobs start, the jobs already running finish (and are cached),
and the first error is logged with its job and rethrown.

`progress = k` prints a line every `k` finished jobs.

```julia
# 20 runs of each of three mutation scalars, cached, resumable
jobs = [(; m, seed) for m in (0.001, 0.005, 0.02) for seed in 1:20]
rows = sweep(jobs; cache = "runs.csv", id = j -> "\$(j.m)-\$(j.seed)") do j
    o = simulate(Params(MutationScalar = j.m, NoGenerations = 2000); seed = j.seed)
    (; j.m, j.seed, VG = o.geno_var[end])
end
```
"""
function sweep(f, jobs; cache::Union{Nothing,AbstractString} = nothing, id = identity,
               cost = nothing, ntasks::Integer = Threads.nthreads(), progress::Integer = 0)
    ntasks >= 1 || throw(ArgumentError("ntasks must be at least 1, got $ntasks"))
    js = collect(jobs)
    n = length(js)
    ids = cache === nothing ? String[] : [string(id(j)) for j in js]
    if cache !== nothing
        allunique(ids) || throw(ArgumentError("two jobs have the same id"))
        any(k -> occursin('\n', k) || occursin('\r', k), ids) &&
            throw(ArgumentError("a job id contains a line break"))
    end
    header, raw = cache === nothing ? (nothing, Dict{String,Vector{String}}()) : read_cache!(cache)
    todo = [i for i in 1:n if cache === nothing || !haskey(raw, ids[i])]
    if cost !== nothing
        todo = todo[sortperm([cost(js[i]) for i in todo]; rev = true, alg = Base.Sort.DEFAULT_STABLE)]
    end

    results = Vector{Any}(undef, n)
    lk = ReentrantLock()
    io = Ref{Union{Nothing,IO}}(nothing)
    next = Threads.Atomic{Int}(1)
    failed = Threads.Atomic{Bool}(false)
    firsterr = Ref{Any}(nothing)
    ndone = Threads.Atomic{Int}(0)
    t0 = time()

    function store!(i, r)
        if cache === nothing
            results[i] = r
            return
        end
        r isa NamedTuple || throw(ArgumentError("with a cache, f must return a NamedTuple; it returned a $(typeof(r))"))
        names = ["id"; [string(k) for k in keys(r)]]
        line = join((csvfield(v) for v in (ids[i], values(r)...)), ',')
        lock(lk) do
            if header === nothing
                io[] = open(cache, "a")
                println(io[], join(names, ','))
                header = names
            elseif header != names
                throw(ArgumentError("the cache $cache has the columns $(join(header, ", ")), " *
                                    "but a job returned $(join(names[2:end], ", "))"))
            end
            io[] === nothing && (io[] = open(cache, "a"))
            println(io[], line)
            flush(io[])
            raw[ids[i]] = csvsplit(line)[2:end]
        end
    end

    function worker()
        while !failed[]
            local k, i, bt, d
            k = Threads.atomic_add!(next, 1)
            k > length(todo) && break
            i = todo[k]
            try
                store!(i, f(js[i]))
            catch e
                bt = catch_backtrace()
                lock(lk) do
                    firsterr[] === nothing && (firsterr[] = (e, bt, i))
                end
                failed[] = true
                break
            end
            d = Threads.atomic_add!(ndone, 1) + 1
            progress > 0 && d % progress == 0 &&
                println("  $d/$(length(todo)) jobs done ($(round(Int, time() - t0)) s)")
        end
    end

    try
        nt = min(ntasks, length(todo))
        if nt <= 1
            worker()
        else
            foreach(wait, [Threads.@spawn worker() for _ in 1:nt])
        end
    finally
        io[] === nothing || close(io[])
    end
    if firsterr[] !== nothing
        err, errbt, erri = firsterr[]
        @error "sweep: job $erri of $n failed" job = js[erri] exception = (err, errbt)
        throw(err)
    end
    cache === nothing && return n == 0 ? Any[] : identity.(results)
    n == 0 && return NamedTuple[]
    parse_rows(header, [raw[k] for k in ids])
end

"""
    replicate(p::Params, n; seed = 1, summary = identity, init = :zero, burnin = nothing, kw...)

Run `n` simulations of `p` with the seeds `seed, seed + 1, …, seed + n - 1`, in parallel,
and return `summary` of each run's output (see [`simulate`](@ref)), in the order of the
seeds. The runs go through [`sweep`](@ref), which takes the other keywords. A `cache` is
keyed by the seed, and `summary` must then return a `NamedTuple`. `init`, `burnin`, `at`
and `target` are [`simulate`](@ref)'s: `init = :equilibrium` with a `burnin` starts every
run from a population drawn from the stationary spectrum ([`equilibrium_world`](@ref)).

```julia
VG = replicate(Params(NoGenerations = 2000), 100; summary = o -> o.geno_var[end])
```
"""
function replicate(p::Params, n::Integer; seed::Integer = 1, summary = identity, init::Symbol = :zero,
                   burnin::Union{Nothing,Integer} = nothing, at::Integer = first_event(p), target = :genic, kw...)
    init === :equilibrium || return sweep(s -> summary(simulate(p; seed = s, init, burnin, at, target)), seed:seed+n-1; kw...)
    spectrum = stationary_spectrum(p)                 # once for all the runs, as the target
    target = init_target(p, spectrum, target)
    sweep(s -> summary(simulate(p; seed = s, init, burnin, at, target, spectrum)), seed:seed+n-1; kw...)
end

# --- the CSV cache ---------------------------------------------------------

function csvfield(v)
    v isa Union{Real,AbstractString,Symbol} ||
        throw(ArgumentError("a cached value must be a number, Bool, string or symbol, not a $(typeof(v))"))
    s = string(v)
    (occursin('\n', s) || occursin('\r', s)) && throw(ArgumentError("a cached value contains a line break: $(repr(s))"))
    occursin(',', s) || occursin('"', s) ? '"' * replace(s, '"' => "\"\"") * '"' : s
end

function csvsplit(line::AbstractString)
    fields = String[]; buf = IOBuffer(); quoted = false; i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if quoted
            if c == '"'
                j = nextind(line, i)
                if j <= lastindex(line) && line[j] == '"'
                    write(buf, '"'); i = j
                else
                    quoted = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            quoted = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    fields
end

# Read a cache: its header (nothing if empty or absent) and its rows by id. A last line
# without its newline, cut short by an interruption, is removed from the file.
function read_cache!(path)
    isfile(path) || return nothing, Dict{String,Vector{String}}()
    s = read(path, String)
    if !isempty(s) && !endswith(s, '\n')
        s = s[1:something(findlast('\n', s), 0)]
        write(path, s)
    end
    lines = split(chomp(s), '\n'; keepempty = false)
    isempty(lines) && return nothing, Dict{String,Vector{String}}()
    header = csvsplit(lines[1])
    (isempty(header) || header[1] != "id") && throw(ArgumentError("$path is not a sweep cache: its first column is not `id`"))
    raw = Dict{String,Vector{String}}()
    for ln in lines[2:end]
        f = csvsplit(ln)
        length(f) == length(header) ||
            throw(ArgumentError("$path: a row has $(length(f)) fields, the header $(length(header))"))
        raw[f[1]] = f[2:end]
    end
    header, raw
end

function parse_column(v::Vector{String})
    all(x -> x == "true" || x == "false", v) && return v .== "true"
    all(x -> tryparse(Int, x) !== nothing, v) && return parse.(Int, v)
    all(x -> tryparse(Float64, x) !== nothing, v) && return parse.(Float64, v)
    v
end

function parse_rows(header, rows)
    names = Tuple(Symbol.(header[2:end]))
    cols = [parse_column([r[c] for r in rows]) for c in eachindex(names)]
    [NamedTuple{names}(Tuple(col[i] for col in cols)) for i in eachindex(rows)]
end
