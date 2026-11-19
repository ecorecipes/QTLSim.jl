# The change of defaults of 2026-11-19 (effective_size → :breeders, drift_critical_rate → :bisection) leaves the moments
# lowering and the emulators unchanged: this script writes their outputs on 24 published inputs (3 per stratum, seed
# 20261119) as bit patterns, to be run once with the package before the change (e35f5e3) and once after, and compared:
#   julia --project=<package>/vignettes validate/ne_default_unchanged.jl out.txt
#   diff <(grep -v NOTBITWISE before.txt) <(grep -v NOTBITWISE after.txt)
# Lines marked NOTBITWISE are not expected to be identical (a search whose bracket starts from drift_critical_rate, and
# descriptive quantities). Result: validate/ne_default_unchanged.log.
using QTLSim, Random, Printf
import Catlab
const SIM = normpath(joinpath(@__DIR__, "..", "..", "Repository", "Simulations"))
bits(x::Float64) = string(reinterpret(UInt64, x), base = 16, pad = 16)
bits(x::Real) = bits(Float64(x))
bits(x::Missing) = "missing"
bits(x::Bool) = string(x)
bits(v::AbstractArray) = string(hash(reinterpret(UInt64, Float64.(collect(v)))), base = 16)
function inputs()
    rng = Xoshiro(20261119)
    ps = Pair{String,Params}[]
    for dir in sort(filter(startswith("QT"), readdir(SIM)))
        files = sort(filter(endswith(".txt"), readdir(joinpath(SIM, dir))))
        for f in files[randperm(rng, length(files))[1:3]]
            push!(ps, "$dir/$f" => Params(read_qtsim_options(joinpath(SIM, dir, f)); Seed = nothing))
        end
    end
    ps
end
open(ARGS[1], "w") do io
    out(name, key, val) = (println(io, name, "\t", key, "\t", val); flush(io))
    for (name, p) in inputs()
        L = p.NoQTL
        q = Params(p; NoGenerations = 2400, MoveFitnessSDGen = 2000, Sample = 1)
        for loci in (:rare_alleles, :multilocus, :aggregate)
            m = run_moments(q; loci, erosion = true)
            out(name, "run_moments $loci geno_var", bits(m.geno_var)); out(name, "run_moments $loci pheno_mean", bits(m.pheno_mean))
            out(name, "run_moments $loci mean_fitness", bits(m.mean_fitness))
        end
        m = run_moments(q)
        out(name, "run_moments default geno_var", bits(m.geno_var))
        if L <= 1000
            m = run_moments(q; loci = :spectrum)
            out(name, "run_moments spectrum geno_var", bits(m.geno_var)); out(name, "run_moments spectrum mean_fitness", bits(m.mean_fitness))
        end
        for kw in ((;), (; drift = true), (; drift = true, loci = :multilocus))
            e = multilocus_equilibrium(p; kw...)
            out(name, "multilocus_equilibrium $kw", join((bits(e.VG), bits(e.genic), bits(e.ld), bits(e.nparents), bits(e.Wbar)), " "))
        end
        out(name, "rare_alleles_vg drift", bits(rare_alleles_vg(p; drift = true)))
        out(name, "gaussian_vg_loci drift", bits(gaussian_vg_loci(p; drift = true)))
        e = emulate(p)
        out(name, "emulate rare_alleles", join((bits(getfield(e, k)) for k in keys(e)), " "))
        e = emulate(p; replicates = 4, seed = 7)
        out(name, "emulate replicates=4", join((bits(getfield(e, k)) for k in keys(e)), " "))
        if L <= 1000
            e = emulate(p; loci = :spectrum, spectrum = (; solver = :beta))
            out(name, "emulate spectrum beta", join((bits(getfield(e, k)) for k in keys(e)), " "))
            e = emulate(p; loci = :spectrum, spectrum = (; solver = :beta, ld = :dynamic))
            out(name, "emulate spectrum beta ld", join((bits(getfield(e, k)) for k in keys(e)), " "))
            L <= 100 && (e = emulate(p; loci = :spectrum); out(name, "emulate spectrum grid", join((bits(getfield(e, k)) for k in keys(e)), " ")))
        end
        h = emulated_half_survival_rate(p)
        out(name, "NOTBITWISE emulated_half_survival_rate k50", @sprintf("%.12e", h.k50)); out(name, "NOTBITWISE emulated_half_survival_rate kc", @sprintf("%.6e", h.kc))
        # what the emulator would give from the house of Gauss at N_e = N exactly instead of the kept start
        e2 = emulate(p; VG = house_of_gauss_vg(p; ne = :falconer)); e3 = emulate(p; VG = QTLSim._house_of_gauss(p, 2.0 * p.PopSize))
        out(name, "NOTBITWISE start N+½ vs N: P2000, VG_move rel diff", @sprintf("%.3e %.3e", e3.P2000 - e2.P2000, e3.VG_move / e2.VG_move - 1))
        out(name, "NOTBITWISE drift_critical_rate new-default/falconer+iteration", @sprintf("%.4f", drift_critical_rate(p).k / drift_critical_rate(p; ne = :falconer, solver = :iteration).k))
    end
end
