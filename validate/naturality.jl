# Run the table of naturality cases (validate/naturality_cases.jl) with `check_naturality`
# and report, for each, the statistics, the verdict and whether it is the verdict the
# ported script recorded on the same cell.
#
# Usage (from the package root; Catlab comes from the vignettes environment):
#     julia -t 3 --project=vignettes validate/naturality.jl [pattern]
# runs the cases whose name contains `pattern` (all by default; the SLiM case is skipped
# when SLiM is not found). Writes validate/naturality.log when every case is run; exits 1
# if a verdict is not reproduced.

using QTLSim, Statistics, Random
import Catlab
using Catlab: ⋅, ⊗
include(joinpath(@__DIR__, "naturality_cases.jl"))

const PATTERN = isempty(ARGS) ? "" : ARGS[1]
const Ext = Base.get_extension(QTLSim, :QTLSimCatlabExt)
cases = [c for c in NATURALITY_CASES if occursin(PATTERN, c.name)]
const out = IOBuffer()
say(x...) = (s = string(x...); println(s); println(out, s); flush(stdout))

say("QTLSim validate/naturality.jl: the table of naturality cases (validate/naturality_cases.jl)")
say("threads ", Threads.nthreads(), "; ", length(cases), " of ", length(NATURALITY_CASES), " cases",
    isempty(PATTERN) ? "" : " matching \"$PATTERN\"")
notrep = String[]
skipped = String[]
t0 = time()
for c in cases
    if (c.from === :slim || c.to === :slim) && Ext.harness_slim() === nothing
        push!(skipped, c.name); say("\nSKIPPED (no SLiM): ", c.name); continue
    end
    r = check_naturality(c)
    say()
    say(sprint(show, MIME"text/plain"(), r))
    r.reproduced || push!(notrep, c.name)
end
say()
say("total ", round(Int, time() - t0), " s; ", length(cases) - length(skipped) - length(notrep), " verdicts reproduced, ",
    length(notrep), " not reproduced, ", length(skipped), " skipped")
isempty(notrep) || say("NOT REPRODUCED: ", join(notrep, "; "))
isempty(PATTERN) && write(joinpath(@__DIR__, "naturality.log"), take!(out))
isempty(notrep) || exit(1)
