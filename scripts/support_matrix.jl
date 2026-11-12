# Write the support matrix (docs/src/support.md) from `support_matrix()`: every cell a probe
# of the code, read against the declared validation evidence (ext/QTLSimCatlabExt/support.jl).
# test/runtests.jl checks that the committed page is the generated one.
#
# Usage (from the package root): julia --project=vignettes scripts/support_matrix.jl

using QTLSim
import Catlab

page = joinpath(@__DIR__, "..", "docs", "src", "support.md")
rows = support_matrix()
write(page, support_markdown(rows))
for st in (:supported, :untested, :unsupported, :na, :error)
    println(rpad(st, 12), count(r -> r.status === st, rows))
end
println("wrote ", normpath(page))
