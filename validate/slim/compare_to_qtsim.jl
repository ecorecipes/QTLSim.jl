# Statistical validation of a SLiM/Eidos model of QTSim against the reference C
# simulation QTSim: by default qtsim_generated.slim, the model generated from QTLSim's life
# cycle (see check_generated.jl), or the model given as the first argument, such as the
# hand-written port qtsim.slim. Runs both R times on matched scenarios and compares the
# distributions of summary metrics with two-sample KS tests. Scenario C, QTSim's
# continuous move of the optimum, runs only for a model with a `moveType` parameter.
#
# SLiM and QTSim use different RNG streams, so the claim tested is equality of
# DISTRIBUTION, not identical runs.
#
# Exits 1 if any p-value is below 0.001: with fixed seeds the outcome is deterministic for
# a given SLiM and QTSim build, and at that level a failure is a discrepancy, not chance
# (about a dozen tests are run; at 0.05 one or two "CHECK"s are expected).
#
# Run from the package root (R replicates per program, default 20; uses all threads):
#   julia -t auto --project=validate/slim validate/slim/compare_to_qtsim.jl [model.slim] [R]

using QTLSim, Statistics, HypothesisTests, Printf
include(joinpath(@__DIR__, "slim.jl"))
include(joinpath(@__DIR__, "..", "qtsim.jl"))

const SLIM  = something(find_slim(), "")
const QTSIM = something(find_qtsim(), "")
const SLIMSCRIPT = abspath(get(ARGS, 1, joinpath(@__DIR__, "qtsim_generated.slim")))
const R = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20

# --- run the SLiM port; parse its TSV stdout ------------------------------
function run_slim(scen::Dict, seed::Int)
    args = ["-s", string(seed)]
    for (k, v) in scen; append!(args, ["-d", v isa AbstractString ? "$k=\"$v\"" : "$k=$v"]); end
    out = read(`$SLIM $args $SLIMSCRIPT`, String)
    gen = Int[]; pm = Float64[]; gv = Float64[]; ns = Int[]; th = Float64[]
    body = false
    for ln in split(out, '\n')
        startswith(ln, "Pop all died") && return (gen=gen, pheno_mean=pm, geno_var=gv, nsurv=ns, theta=th, extinct=true)
        if startswith(ln, "gen\t"); body = true; continue; end
        body || continue
        f = split(ln, '\t'); length(f) < 5 && continue
        g = tryparse(Int, f[1]); g === nothing && continue
        push!(gen, g); push!(pm, parse(Float64, f[2])); push!(gv, parse(Float64, f[3]))
        push!(ns, parse(Int, f[4])); push!(th, parse(Float64, f[5]))
    end
    (gen=gen, pheno_mean=pm, geno_var=gv, nsurv=ns, theta=th, extinct=false)
end

# --- QTSim reference ------------------------------------------------------
function write_opts(path, scen, seed)
    open(path, "w") do io
        println(io, "Sample 1")
        println(io, "PopSize $(scen["N"])")
        println(io, "NoQTL $(get(scen, "L", 1))")
        haskey(scen, "RECOMB") && println(io, "RecombinationPoissonMean $(scen["RECOMB"])")
        println(io, "Generations $(scen["NGEN"])")
        println(io, "Ploidy 2")
        scen["moveGen"] != -1 && println(io, "MoveFitnessSD $(scen["moveGen"]) $(scen["moveSD"])")
        scen["changeGen"] != -1 && println(io, "ChangeFitnessSD $(scen["changeGen"]) $(scen["changeSD"])")
        println(io, "OptimumMoveType $(get(scen, "moveType", "phenotype"))")
        println(io, "PhenotypeVar $(scen["VP"])")
        println(io, "MutationScalar $(scen["m"])")
        println(io, "Heritability $(scen["h2"])")
        println(io, "MutationRatePerGamete $(scen["rate"])")
        scen["initOmega"] >= 0 && println(io, "InitFitnessSD $(scen["initOmega"])")
        println(io, "Seed $seed")
    end
end

function run_qtsim(scen, seed)
    o = tempname() * ".txt"; write_opts(o, scen, seed)
    s = try; read(`$QTSIM $o`, String); finally; rm(o; force=true); end
    gen = Int[]; pm = Float64[]; gv = Float64[]; ns = Int[]; th = Float64[]; b = false
    for ln in split(s, '\n')
        startswith(ln, "Pop all died") && return (gen=gen, pheno_mean=pm, geno_var=gv, nsurv=ns, theta=th, extinct=true)
        startswith(ln, "Generation") && (b = true; continue)
        (b && !isempty(strip(ln))) || continue
        f = split(ln, '\t'); length(f) < 6 && continue
        push!(gen, parse(Int, f[1])); push!(pm, parse(Float64, f[2])); push!(gv, parse(Float64, f[4]))
        push!(ns, parse(Int, f[5])); push!(th, parse(Float64, f[6]))
    end
    (gen=gen, pheno_mean=pm, geno_var=gv, nsurv=ns, theta=th, extinct=false)
end

# --- metrics (identical for both) ----------------------------------------
win(r, lo, hi) = findall(g -> lo <= g <= hi, r.gen)
eq_VG(r, lo, hi) = mean(@view r.geno_var[win(r, lo, hi)])
end_VG(r, ng)    = mean(@view r.geno_var[win(r, ng - 100, ng)])
survfrac(r, mg, ng, N) = mean(@view r.nsurv[win(r, mg, ng)]) / N
adv(r, mg, ng) = mean(r.pheno_mean[win(r, ng - 100, ng)]) - mean(r.pheno_mean[win(r, mg, mg + 10)])

ksD(a, b) = (xs = sort(vcat(a, b)); maximum(abs(count(<=(x), a)/length(a) - count(<=(x), b)/length(b)) for x in xs))
ksP(a, b) = pvalue(ApproximateTwoSampleKSTest(Float64.(a), Float64.(b)))
const PMIN = Ref(1.0)
function report(name, a, b)
    PMIN[] = min(PMIN[], ksP(a, b))
    @printf("  %-22s SLiM %9.4g ± %-8.2g  QTSim %9.4g ± %-8.2g  D=%.3f p=%.3f  [%s]\n",
            name, mean(a), std(a), mean(b), std(b), ksD(a, b), ksP(a, b),
            ksP(a, b) > 0.05 ? "consistent" : "CHECK")
end

# --- scenarios ------------------------------------------------------------
# run R replicates on all threads, emitting a progress dot each (keeps stdout moving)
function collect_runs(runner, scen, base)
    out = Vector{Any}(undef, R)
    Threads.@threads :greedy for s in 1:R
        out[s] = runner(scen, base + s)
        print("."); flush(stdout)
    end
    println(); out
end
A = Dict("N"=>1000, "VP"=>0.02, "h2"=>0.5, "m"=>0.02, "rate"=>0.5, "NGEN"=>4000,
         "moveGen"=>-1, "moveSD"=>1.0, "changeGen"=>-1, "changeSD"=>-1.0, "initOmega"=>0.3)
B = Dict("N"=>1000, "VP"=>0.02, "h2"=>0.5, "m"=>0.001, "rate"=>0.0002, "NGEN"=>6500,
         "moveGen"=>5000, "moveSD"=>1.0, "changeGen"=>5000, "changeSD"=>0.30, "initOmega"=>-1.0)
# QTSim's continuous move, as in the paper's sweep: the optimum advances by moveSD per generation
C = Dict("N"=>1000, "VP"=>0.02, "h2"=>0.5, "m"=>0.001, "rate"=>0.0002, "NGEN"=>6500,
         "moveGen"=>5000, "moveSD"=>0.002, "moveType"=>"continuous", "changeGen"=>-1,
         "changeSD"=>-1.0, "initOmega"=>-1.0)

isempty(SLIM) && error(slim_missing())
isempty(QTSIM) && error(qtsim_missing())

println("SLiM model ", basename(SLIMSCRIPT), " vs QTSim — R=$R replicates each\n")

println("── Scenario A: stabilizing equilibrium ──")
sA = collect_runs(run_slim, A, 10_000)
qA = collect_runs(run_qtsim, A, 20_000)
report("equilibrium V_G", [eq_VG(r, 3000, 3999) for r in sA], [eq_VG(r, 3000, 3999) for r in qA])

println("\n── Scenario B: directional adaptation ──")
sB = collect_runs(run_slim, B, 30_000)
qB = collect_runs(run_qtsim, B, 40_000)
ks = filter(r -> !r.extinct, sB); kq = filter(r -> !r.extinct, qB)
report("end V_G",           [end_VG(r, 6500) for r in ks],            [end_VG(r, 6500) for r in kq])
report("survival fraction", [survfrac(r, 5000, 6500, 1000) for r in ks], [survfrac(r, 5000, 6500, 1000) for r in kq])
report("phenotypic advance",[adv(r, 5000, 6500) for r in ks],        [adv(r, 5000, 6500) for r in kq])
@printf("  %-22s SLiM %d/%d               QTSim %d/%d\n", "extinctions",
        count(r->r.extinct, sB), R, count(r->r.extinct, qB), R)

if occursin("moveType", read(SLIMSCRIPT, String))
    println("\n── Scenario C: continuous move of the optimum ──")
    sC = collect_runs(run_slim, C, 50_000)
    qC = collect_runs(run_qtsim, C, 60_000)
    ks = filter(r -> !r.extinct, sC); kq = filter(r -> !r.extinct, qC)
    report("end V_G",           [end_VG(r, 6500) for r in ks],               [end_VG(r, 6500) for r in kq])
    report("survival fraction", [survfrac(r, 5000, 6500, 1000) for r in ks], [survfrac(r, 5000, 6500, 1000) for r in kq])
    report("phenotypic advance",[adv(r, 5000, 6500) for r in ks],           [adv(r, 5000, 6500) for r in kq])
    @printf("  %-22s SLiM %d/%d               QTSim %d/%d\n", "extinctions",
            count(r->r.extinct, sC), R, count(r->r.extinct, qC), R)
end

# Many loci, the design of the paper's simulations: only for a model with QTSim's
# recombination (`RECOMB`), such as one generated from QTLSim's life cycle.
D = Dict("N"=>1000, "VP"=>0.03, "h2"=>0.4, "m"=>0.005, "rate"=>0.0001, "L"=>100, "RECOMB"=>3.0,
         "NGEN"=>2000, "moveGen"=>-1, "moveSD"=>1.0, "changeGen"=>-1, "changeSD"=>-1.0, "initOmega"=>4.0)
E = merge(D, Dict("RECOMB"=>5.0, "moveGen"=>1000, "moveSD"=>0.05, "moveType"=>"continuous"))
if occursin("RECOMB", read(SLIMSCRIPT, String))
    println("\n── Scenario D: 100 loci, three crossovers per gamete on average, stabilizing selection ──")
    sD = collect_runs(run_slim, D, 70_000)
    qD = collect_runs(run_qtsim, D, 80_000)
    report("V_G, gens 1500-1999", [eq_VG(r, 1500, 1999) for r in sD], [eq_VG(r, 1500, 1999) for r in qD])
    report("survival fraction", [survfrac(r, 1500, 1999, 1000) for r in sD], [survfrac(r, 1500, 1999, 1000) for r in qD])

    println("\n── Scenario E: 100 loci, five crossovers, the optimum moving continuously from generation 1000 ──")
    sE = collect_runs(run_slim, E, 90_000)
    qE = collect_runs(run_qtsim, E, 100_000)
    after_move(r) = (r.extinct ? last(r.gen) + 1 : 2000) - 1000      # neither prints the dying generation
    report("V_G before the move", [eq_VG(r, 900, 999) for r in sE], [eq_VG(r, 900, 999) for r in qE])
    report("gens to extinction", after_move.(sE), after_move.(qE))
    @printf("  %-22s SLiM %d/%d               QTSim %d/%d\n", "extinctions",
            count(r->r.extinct, sE), R, count(r->r.extinct, qE), R)
end
@printf("\nsmallest p = %.3g: %s\n", PMIN[], PMIN[] < 0.001 ? "DISCREPANCY" : "no discrepancy")
PMIN[] < 0.001 && exit(1)
