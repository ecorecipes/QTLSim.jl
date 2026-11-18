using QTLSim
using Test, Statistics, Random
using QTLSim: eigvals, Symmetric   # LinearAlgebra, through the package

# The categorical layer is the package extension QTLSimCatlabExt, loaded with Catlab.
import Catlab
using Catlab: ⋅, ⊗, dom, codom, Hom
const Ext = Base.get_extension(QTLSim, :QTLSimCatlabExt)
(; develop, select, reproduce, mutate, swap, generation) = operators()

# Package hygiene (ambiguities, piracy, stale dependencies, compat bounds), when Aqua is
# available, as it is under Pkg.test.
if Base.find_package("Aqua") !== nothing && isempty(ARGS) && isempty(get(ENV, "QTLSIM_TEST", ""))
    @eval using Aqua
    @testset "Aqua" begin
        Aqua.test_all(QTLSim; persistent_tasks = false)
    end
end

# Topics, in the order they run. Run a subset by naming topics as arguments,
# `Pkg.test(test_args = ["traits", "demes"])` or `julia --project=vignettes test/runtests.jl
# traits demes`, or in ENV["QTLSIM_TEST"] (space- or comma-separated). Names are the files
# test/<topic>.jl.
const TOPICS = ["core", "catlab", "blocks", "io", "theory", "demography", "closures", "demes", "fluctuating", "lowering", "sweep", "traits", "init", "ancestry", "plasticity", "declarations", "critical", "laws", "emulator", "block_lineages", "two_stage", "k50_gradient", "xu_shen"]
const CHOSEN = let a = isempty(ARGS) ? split(get(ENV, "QTLSIM_TEST", ""), r"[\s,]+"; keepempty = false) : ARGS
    bad = setdiff(a, TOPICS)
    isempty(bad) || error("unknown test topic(s) $(join(bad, ", ")); topics are $(join(TOPICS, ", "))")
    isempty(a) ? TOPICS : String.(a)
end

for t in TOPICS
    t in CHOSEN || continue
    @testset "$t" begin
        include(joinpath(@__DIR__, "$t.jl"))
    end
end
