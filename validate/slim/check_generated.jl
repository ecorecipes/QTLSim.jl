# The SLiM model generated from QTLSim's life cycle (`slim_script`, the lowering of
# develop ⋅ select ⋅ reproduce ⋅ record ⋅ rescale ⋅ mutate ⋅ swap to Eidos) against the
# hand-written port qtsim.slim, run for run: for the same seed and parameters both must
# print the same table, number for number. The one intended difference is at extinction,
# where the port also prints the dying generation's line and QTSim (like the generated
# model and the Julia kernel) does not.
#
# First checks that `slim_script` still writes the committed qtsim_generated.slim (the
# package tests check the same without SLiM); with `--write` it rewrites the file instead,
# after an intended change to the Eidos images. Exits 1 on any difference.
#
# Run from the package root:
#   julia --project=validate/slim -e 'using Pkg; Pkg.instantiate()'    # once
#   julia --project=validate/slim validate/slim/check_generated.jl [--write]

using QTLSim, Printf
import Catlab
include(joinpath(@__DIR__, "slim.jl"))
include(joinpath(@__DIR__, "golden.jl"))

const SLIM = find_slim()
SLIM === nothing && error(slim_missing())
const PORT = joinpath(@__DIR__, "qtsim.slim")
const GENERATED = joinpath(@__DIR__, "qtsim_generated.slim")

for (name, script) in golden_models()
    path = joinpath(@__DIR__, name)
    if "--write" in ARGS
        write(path, script)
        println("wrote ", name)
    elseif !isfile(path) || script != read(path, String)
        println("slim_script no longer writes $name; rerun with --write if the change is intended")
        exit(1)
    else
        println("slim_script writes $name unchanged")
    end
end

# the models with many loci and with demes run in SLiM and print their tables
for (name, rows) in (("qtsim_generated_loci.slim", 50), ("qtsim_generated_demes.slim", 150))
    out = read(`$SLIM -s 1 -d NGEN=50 $(joinpath(@__DIR__, name))`, String)
    n = count(l -> occursin(r"^\d+\t", l), split(out, '\n'))
    n == rows || (println("$name printed $n table rows, expected $rows"); exit(1))
    println("$name runs in SLiM ($n rows)")
end

function run_table(script, seed, d)
    args = ["-s", string(seed)]
    for (k, v) in d; append!(args, ["-d", "$k=$v"]); end
    out = read(`$SLIM $args $script`, String)
    rows = [l for l in split(out, '\n') if occursin(r"^\d+\t", l)]
    (; rows, died = occursin("Pop all died.", out))
end

scenarios = (
    "stabilising (A)" => Dict("N" => 1000, "VP" => 0.02, "h2" => 0.5, "m" => 0.02, "rate" => 0.5, "NGEN" => 600,
                              "moveGen" => -1, "moveSD" => 1.0, "changeGen" => -1, "changeSD" => -1.0, "initOmega" => 0.3),
    "moving optimum (B)" => Dict("N" => 1000, "VP" => 0.02, "h2" => 0.5, "m" => 0.001, "rate" => 0.0002, "NGEN" => 600,
                                 "moveGen" => 300, "moveSD" => 1.0, "changeGen" => 300, "changeSD" => 0.30, "initOmega" => -1.0),
    "extinction" => Dict("N" => 300, "VP" => 0.02, "h2" => 0.5, "m" => 0.001, "rate" => 0.0002, "NGEN" => 600,
                         "moveGen" => 50, "moveSD" => 6.0, "changeGen" => 50, "changeSD" => 0.05, "initOmega" => -1.0))

allsame = true
for (name, d) in scenarios
    same = 0; died = 0
    for seed in 1:10
        a = run_table(PORT, seed, d); b = run_table(GENERATED, seed, d)
        died += a.died
        ref = a.died ? a.rows[1:end-1] : a.rows      # the port also prints the dying generation
        ok = a.died == b.died && ref == b.rows
        same += ok
        global allsame &= ok
    end
    @printf("%-20s %2d of 10 seeds identical (%d went extinct)\n", name, same, died)
end
println(allsame ? "the generated model is the hand-written port, run for run" : "MISMATCH")
allsame || exit(1)
