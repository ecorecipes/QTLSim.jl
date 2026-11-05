# Statistical-consistency check: compositional Julia model vs. the reference
# QTSim C simulation. For each scenario, runs both R times and compares the
# distributions of summary statistics with Welch t-tests.
#
# Usage (from the package root), after building the C binary
# (../Repository/Code/QTSim/src/QTSim):
#     julia --project=. validate/compare_to_qtsim.jl
#
# Different RNGs (Julia MersenneTwister vs GSL mt19937) mean runs are not
# identical; the claim being tested is that the two produce the SAME
# distributions. |t| ≳ 2.5 on any metric flags a discrepancy worth checking.

using CategoricalQuantGen
using Statistics, Random, Printf

const QTSIM = normpath(joinpath(@__DIR__, "..", "..",
                                "Repository", "Code", "QTSim", "src", "QTSim"))

# --- QTSim reference -----------------------------------------------------

function write_opts(path, p::Params, seed::Int)
    open(path, "w") do io
        println(io, "Sample $(p.Sample)")
        println(io, "PopSize $(p.PopSize)")
        println(io, "NoQTL $(p.NoQTL)")
        println(io, "Generations $(p.NoGenerations)")
        println(io, "Ploidy $(p.Ploidy)")
        p.MoveFitnessSDGen != -1 &&
            println(io, "MoveFitnessSD $(p.MoveFitnessSDGen) $(p.MoveFitnessSD)")
        p.ChangeFitnessGen != -1 &&
            println(io, "ChangeFitnessSD $(p.ChangeFitnessGen) $(p.ChangeFitnessSD)")
        println(io, "OptimumMoveType $(p.MoveType)")
        println(io, "PhenotypeVar $(p.PhenotypeVar)")
        println(io, "MutationScalar $(p.MutationScalar)")
        println(io, "Heritability $(p.Heritability)")
        println(io, "MutationRatePerGamete $(p.MutationRatePerGamete)")
        p.InitFitnessSD >= 0 && println(io, "InitFitnessSD $(p.InitFitnessSD)")
        println(io, "Seed $seed")
    end
end

# Parse QTSim stdout into the same shape as `simulate`'s output.
function run_qtsim(p::Params, seed::Int)
    opts = tempname() * ".txt"
    write_opts(opts, p, seed)
    out = try
        read(`$QTSIM $opts`, String)
    finally
        rm(opts; force=true)
    end
    gen = Int[]; pmean = Float64[]; gvar = Float64[]; nsurv = Int[]; θ = Float64[]
    inbody = false
    for line in split(out, '\n')
        if startswith(line, "Pop all died")
            return (extinct=true, gen=gen, pheno_mean=pmean, geno_var=gvar,
                    nsurv=nsurv, theta=θ)
        elseif startswith(line, "Generation")
            inbody = true; continue
        end
        (inbody && !isempty(strip(line))) || continue
        f = split(line, '\t')
        length(f) < 6 && continue
        push!(gen, parse(Int, f[1]))
        push!(pmean, parse(Float64, f[2]))   # Phenotype Mean
        push!(gvar, parse(Float64, f[4]))    # Genotype Var
        push!(nsurv, parse(Int, f[5]))       # P Survived (count)
        push!(θ, parse(Float64, f[6]))       # Fitness Function Mean (optimum)
    end
    (extinct=false, gen=gen, pheno_mean=pmean, geno_var=gvar, nsurv=nsurv, theta=θ)
end

# --- summaries & comparison ---------------------------------------------

pick(gen, lo, hi) = findall(g -> lo <= g <= hi, gen)

function metrics(r, p::Params, eqlo, eqhi)
    vg = mean(@view r.geno_var[pick(r.gen, eqlo, eqhi)])
    if p.MoveFitnessSDGen != -1
        dr = pick(r.gen, p.MoveFitnessSDGen, p.NoGenerations)
        en = pick(r.gen, p.NoGenerations - 100, p.NoGenerations)
        return (eqVG=vg, lag=mean(r.theta[en] .- r.pheno_mean[en]),
                surv=mean(r.nsurv[dr]) / p.PopSize)
    end
    (eqVG=vg, lag=NaN, surv=NaN)
end

function welch(a, b)
    ma, mb = mean(a), mean(b)
    va, vb = var(a), var(b)
    t = (ma - mb) / sqrt(va / length(a) + vb / length(b) + eps())
    (ma=ma, sa=std(a), mb=mb, sb=std(b), t=t)
end

function report(name, a, b)
    w = welch(a, b)
    verdict = abs(w.t) < 2.5 ? "consistent" : "CHECK"
    @printf("  %-20s  Julia %10.4g ± %-9.3g  C %10.4g ± %-9.3g   t=%+6.2f  [%s]\n",
            name, w.ma, w.sa, w.mb, w.sb, w.t, verdict)
end

function compare(p::Params, label; eqlo, eqhi, R=24)
    @printf("\n── %s ──\n", label)
    println("  PopSize=$(p.PopSize) Ploidy=$(p.Ploidy) NoQTL=$(p.NoQTL) " *
            "Gens=$(p.NoGenerations) h²=$(p.Heritability) V_P=$(p.PhenotypeVar) " *
            "m=$(p.MutationScalar) rate=$(p.MutationRatePerGamete) " *
            "ω₀=$(round(initFitSD(p), digits=4))  (R=$R each)")
    jl = NamedTuple[]; cl = NamedTuple[]; jext = 0; cext = 0
    for s in 1:R
        rj = simulate(p; seed=10_000 + s)
        rj.extinct ? (jext += 1) : push!(jl, metrics(rj, p, eqlo, eqhi))
        rc = run_qtsim(p, 20_000 + s)
        rc.extinct ? (cext += 1) : push!(cl, metrics(rc, p, eqlo, eqhi))
        print("."); flush(stdout)
    end
    println()
    report("equilibrium V_G", [x.eqVG for x in jl], [x.eqVG for x in cl])
    if p.MoveFitnessSDGen != -1
        report("adaptation lag", [x.lag for x in jl], [x.lag for x in cl])
        report("survival fraction", [x.surv for x in jl], [x.surv for x in cl])
        @printf("  %-20s  Julia %d/%d              C %d/%d\n",
                "extinctions", jext, R, cext, R)
    end
end

# --- scenarios -----------------------------------------------------------

isfile(QTSIM) || error("QTSim binary not found at $QTSIM — build it first " *
    "(see docs/README.md: cc *.c -O2 -o QTSim -I<gsl>/include -L<gsl>/lib -lgsl -lgslcblas -lm)")

# 1. Well-conditioned mutation–selection balance: finite stabilizing width and a
#    high mutation rate ⇒ a low-variance equilibrium V_G that sharply tests the
#    selection + mutation operators. No moving optimum.
equilibrium = Params(
    Sample=10, PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=4000,
    PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.02,
    MutationRatePerGamete=0.5, InitFitnessSD=0.3,
)
compare(equilibrium, "stabilizing equilibrium (mutation–selection balance)";
        eqlo=3000, eqhi=3999, R=40)

# 2. Directional adaptation to a moving optimum — the single-jump protocol of
#    Scrips/CreateCGenomeSimInputDist.py (the published sweep is checked by
#    published_sweep.jl). Tests the
#    optimum-move schedule, sustained directional selection and survival.
directional = Params(
    Sample=1, PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=6500,
    PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.001,
    ChangeFitnessGen=5000, ChangeFitnessSD=0.30,
    MoveFitnessSDGen=5000, MoveFitnessSD=1.0, MoveType=:phenotype,
)
compare(directional, "directional adaptation (moving optimum)";
        eqlo=4000, eqhi=4999, R=24)
