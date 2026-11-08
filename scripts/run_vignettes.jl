# Run the code of the vignettes without Quarto, and fail if any chunk throws.
#
#   julia -t 3 --project=vignettes scripts/run_vignettes.jl              # the default set
#   julia -t 3 --project=vignettes scripts/run_vignettes.jl 03 10        # by number or name
#
# Each vignette runs in its own Julia process, with the driver's thread count, in the
# vignette's directory, as Quarto's Julia engine runs it. Its ```{julia} chunks are
# evaluated in order in one module; `@__DIR__` is the vignette's directory and error
# messages give lines of the .qmd file. A chunk's last value is shown, as in the rendered
# vignette (figures are drawn to PNG in memory), unless the chunk ends in `;`. Chunks with
# `#| eval: false` are skipped.
#
# The default set is the vignettes that need neither the paper's data nor QTSim nor SLiM.
# 07 needs SLiM and QTSim, 08 QTSim (it skips those sections without it, but its last
# section reads the paper's data), and 09 and 11 the paper's results tables.

import REPL

const ROOT = normpath(joinpath(@__DIR__, ".."))
const VIG = joinpath(ROOT, "vignettes")
const DEFAULT = ["01", "02", "03", "04", "05", "06", "10", "12", "13"]

"""
    chunks(qmd) -> Vector{(line, options, code)}

The executable chunks of a Quarto file: the line of each chunk's first code line, its
`#| key: value` options, and its code.
"""
function chunks(qmd::AbstractString)
    out = Tuple{Int,Dict{String,String},String}[]
    lines = readlines(qmd)
    i = 1
    while i <= length(lines)
        if occursin(r"^```\{julia[^}]*\}\s*$", lines[i])
            first = i + 1
            j = first
            while j <= length(lines) && !occursin(r"^```\s*$", lines[j])
                j += 1
            end
            j > length(lines) && error("$qmd:$i: chunk is not closed")
            body = lines[first:(j - 1)]
            opts = Dict{String,String}()
            for l in body
                m = match(r"^#\|\s*([\w-]+)\s*:\s*(.*?)\s*$", l)
                m === nothing && break
                opts[m.captures[1]] = m.captures[2]
            end
            push!(out, (first, opts, join(body, '\n')))
            i = j + 1
        else
            i += 1
        end
    end
    return out
end

ends_with_semicolon(code) = endswith(rstrip(replace(code, r"#[^\n]*$" => "")), ';')

function show_value(x)
    x === nothing && return
    if showable(MIME("image/png"), x)
        show(IOBuffer(), MIME("image/png"), x)
    else
        show(devnull, MIME("text/plain"), x)
    end
    return
end

"Run one vignette in this process; exit 1 at the first chunk that throws."
function run_one(qmd::AbstractString)
    qmd = abspath(qmd)
    cd(dirname(qmd))
    mod = Module(:Notebook)
    cs = chunks(qmd)
    for (k, (line, opts, code)) in enumerate(cs)
        get(opts, "eval", "true") == "false" && continue
        println(stderr, "-- chunk $k/$(length(cs)) (line $line)")
        t = time()
        try
            # pad with blank lines so that errors point at lines of the .qmd file
            # REPL.softscope: Quarto evaluates chunks with the REPL's soft scope, so a loop may assign a global
            x = Base.invokelatest(include_string, REPL.softscope, mod, "\n"^(line - 1) * code, qmd)
            ends_with_semicolon(code) || Base.invokelatest(show_value, x)
        catch err
            flush(stdout)
            println(stderr, "\nERROR in $(basename(qmd)), chunk $k starting at line $line:")
            showerror(stderr, err, catch_backtrace())
            println(stderr)
            exit(1)
        end
        println(stderr, "   $(round(time() - t; digits = 1)) s")
    end
    return
end

function vignette_path(key::AbstractString)
    dirs = filter(d -> occursin(r"^\d+_", d) && isdir(joinpath(VIG, d)), readdir(VIG))
    hits = filter(d -> d == key || startswith(d, key * "_"), dirs)
    length(hits) == 1 || error("no single vignette matches \"$key\" (found $(hits))")
    return joinpath(VIG, only(hits), only(hits) * ".qmd")
end

if length(ARGS) == 2 && ARGS[1] == "--one"
    run_one(ARGS[2])
else
    wanted = isempty(ARGS) ? DEFAULT : ARGS
    paths = vignette_path.(wanted)
    results = Tuple{String,Bool,Float64}[]
    for p in paths
        name = splitext(basename(p))[1]
        println(stderr, "\n== $name")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --threads=$(Threads.nthreads()) $(@__FILE__) --one $p`
        t = time()
        ok = success(run(ignorestatus(cmd)))
        push!(results, (name, ok, time() - t))
    end
    println(stderr, "\nvignette                                   result  time")
    for (name, ok, t) in results
        println(stderr, rpad(name, 43), rpad(ok ? "ok" : "FAILED", 8),
                round(Int, t), " s")
    end
    all(r -> r[2], results) || exit(1)
end
