# Statistical-consistency check: compositional Julia model vs. the reference
# QTSim C simulation. For each scenario, runs both R times and compares the
# distributions of summary statistics with Welch t-tests.
#
# Usage (from the package root), with QTSim built from the bundled source
# (`make -C reference/QTSim`; validate/qtsim.jl says where else it looks):
#     julia --project=. validate/compare_to_qtsim.jl
#
# Different RNGs (Julia MersenneTwister vs GSL mt19937) mean runs are not
# identical; the claim being tested is that the two produce the SAME
# distributions. |t| ≳ 2.5 on any metric flags a discrepancy worth checking.

using QTLSim
using Statistics, Random, Printf

include(joinpath(@__DIR__, "qtsim.jl"))
const QTSIM = something(find_qtsim(), "")

# --- QTSim reference -----------------------------------------------------

write_opts(path, p::Params, seed::Int) = write_qtsim_options(path, p; seed)

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
        st = pick(r.gen, p.MoveFitnessSDGen, p.MoveFitnessSDGen + 10)
        # advance of the mean phenotype after the move (the lag θ − z̄ is no test: under the
        # phenotype move θ is re-anchored to z̄ + √V_P·d every generation, so it is constant)
        return (eqVG=vg, adv=mean(r.pheno_mean[en]) - mean(r.pheno_mean[st]),
                surv=mean(r.nsurv[dr]) / p.PopSize)
    end
    (eqVG=vg, adv=NaN, surv=NaN)
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
    abs(w.t)
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
    t = report("equilibrium V_G", [x.eqVG for x in jl], [x.eqVG for x in cl])
    if p.MoveFitnessSDGen != -1
        t = max(t, report("phenotypic advance", [x.adv for x in jl], [x.adv for x in cl]),
                report("survival fraction", [x.surv for x in jl], [x.surv for x in cl]))
        @printf("  %-20s  Julia %d/%d              C %d/%d\n",
                "extinctions", jext, R, cext, R)
    end
    t                                     # the largest |t| of the scenario
end

# --- scenarios -----------------------------------------------------------

isfile(QTSIM) || error(qtsim_missing())

# 1. Well-conditioned mutation–selection balance: finite stabilizing width and a
#    high mutation rate ⇒ a low-variance equilibrium V_G that sharply tests the
#    selection + mutation operators. No moving optimum.
equilibrium = Params(
    Sample=10, PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=4000,
    PhenotypeVar=0.02, Heritability=0.5, MutationScalar=0.02,
    MutationRatePerGamete=0.5, InitFitnessSD=0.3,
)
worst = compare(equilibrium, "stabilizing equilibrium (mutation–selection balance)";
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
worst = max(worst, compare(directional, "directional adaptation (moving optimum)";
                           eqlo=4000, eqhi=4999, R=24))

# |t| ≥ 2.5 is flagged for a look; |t| ≥ 4, which chance all but never produces with
# these replicate counts, fails the run (CI runs this script against the bundled QTSim).
if worst >= 4
    println("\nFAIL: largest |t| = $(round(worst, digits=2)), a discrepancy between the Julia model and QTSim")
    exit(1)
end
println("\nlargest |t| = $(round(worst, digits=2))")
