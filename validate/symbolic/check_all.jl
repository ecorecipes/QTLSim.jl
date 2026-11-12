# Runs every symbolic check in one process (Symbolics and QTLSim load once): each file
# writes its own log; exits 1 if any check fails. CI's `symbolic` job runs this.
#
#     julia --project=validate/symbolic validate/symbolic/check_all.jl

const FILES = ["check_closures.jl", "check_spectrum.jl", "check_traits.jl", "check_fst.jl", "check_demography.jl", "check_plasticity.jl", "check_init.jl", "check_critical.jl"]

t0 = time()
nfail = Dict{String,Int}()
for f in FILES
    m = Core.eval(Main, :(module $(Symbol(splitext(f)[1])) end))
    Base.include(m, joinpath(@__DIR__, f))
    nfail[f] = Base.invokelatest(getfield(m, :main))
end
println("\n", "="^80)
for f in FILES
    println(rpad(f, 24), nfail[f] == 0 ? "PASS" : "FAIL ($(nfail[f]))")
end
println("total ", round(Int, time() - t0), " s")
exit(sum(values(nfail)) == 0 ? 0 : 1)
