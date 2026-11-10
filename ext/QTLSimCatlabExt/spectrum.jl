# The moments lowering with `loci = :spectrum`: the images of select, reproduce and mutate
# on the expected allele-frequency spectrum of src/spectrum.jl (`QTLSim.AlleleSpectrum`),
# which then replaces the haplotype moments as the `cov` of a deme's `Moments`. The other
# steps (develop, record, rescale, swap) and the demography are the lowering's own.
#
# - select: the lag λ = θ − z̄ and S = ω_e² + V_P set each mutation's E[Δp]
#   (`spectrum_select`); the survivors' mean genotype is z̄ + Σ 2a E[Δp], which includes
#   the alleles' skew, instead of the Gaussian z̄ + V_G λ/S;
# - reproduce: the Wright–Fisher diffusion of every effect class over one generation,
#   with that selection and drift among the lowering's parents (`spectrum_step!`);
# - mutate: 2N L μ new mutations at frequency 1/(2N), N the offspring (`spectrum_mutate!`);
# - swap: V_G = Σ 2a² p(1 − p).
#
# One deme, one diploid trait; drift is always on (the spectrum of a finite population),
# so `erosion` does not matter.

# the spectrum of a deme: Params with the mutation fields the life cycle's mutate steps
# carry (one law of effects), started with genetic variance VG in the shape of the
# equilibrium spectrum (empty when VG = 0, the kernel's zero genomes)
function spectrum_state(p::Params, cycle, VG; kw...)
    mf = mutate_fields(cycle)
    mf === nothing && throw(ArgumentError("loci = :spectrum needs one law of mutation effects; the mutate steps of the " *
                                          "life cycle carry different MutationScalar or MutationRatePerGamete"))
    q = isempty(mf) ? p : Params(p; mf...)
    VG > 0 || return QTLSim.AlleleSpectrum(q; kw...)
    QTLSim.spectrum_scale!(QTLSim.spectrum_equilibrium(q; kw...).spectrum, VG)
end

function cov_select!(m::Moments, sp::QTLSim.AlleleSpectrum)
    S = m.sinv == 0 ? Inf : 1 / m.sinv
    m.zs = m.zbar + QTLSim.spectrum_select(sp, m.θ - m.zbar - m.shift, S)
    m
end

cov_reproduce!(m::Moments, sp::QTLSim.AlleleSpectrum, npar) = (QTLSim.spectrum_step!(sp, npar); m)
cov_mutate!(m::Moments, sp::QTLSim.AlleleSpectrum, Vm) = (QTLSim.spectrum_mutate!(sp, m.Noff); m)
cov_vpop(m::Moments, sp::QTLSim.AlleleSpectrum) = QTLSim.spectrum_vg(sp)
cov_drift(::QTLSim.AlleleSpectrum) = true
