# Replicate trajectories of the moments lowering with drift of the mean (the emulator's
# between-run variance, src/emulator.jl; criteria validate/emulator_variance_criteria.md).
#
# The drift-mode lowering (`erosion = true`) carries, among its haplotype moments, C: the
# expected square of the deviation of a population's realised mean haplotype from the
# deterministic mean z̄/Ploidy, so that Ploidy²·C is the variance among runs of the mean
# genotype around the deterministic trajectory. Here each replicate carries its own
# realised mean instead: every generation the offspring's mean receives a Gaussian
# innovation of variance Ploidy²(C_new − C_sel), the increase of C over reproduction and
# mutation (C_sel after selection, C_new in the offspring), and selection then acts on the
# replicate's own lag, shrinking a deviation by 1 − V_G/s, the factor by which `hap_select`
# shrinks C. In the linear Gaussian recursion the replicates' variance about the
# deterministic mean is then Ploidy²·C at every generation. The second moments (T, R, Q)
# are left as the lowering computes them; they depend on the realised lag only through the
# number of parents N·W̄.
#
# With `vg = (; stationary, innovation)` (QTLSim.vg_fluctuation, the Poisson random field
# of the stationary spectrum) the genetic variance fluctuates too: a replicate starts the
# noisy generations with δ ~ N(0, stationary) added to its V_G, and every generation its
# offspring's V_G receives δ ~ N(0, innovation), each added to the genic class's
# haplotype variance (T₁ += δ/2 for a diploid, so V_G changes by δ up to O(1/N)) and
# clamped so that the genic within-haplotype variance keeps at least a tenth of itself. The
# lowering's own recursion then relaxes the deviation (drift at 1/(2n), selection and the
# rare-alleles loss). Only with haplotype moments (DemeCov); ignored otherwise.

# C of the one deme of `m` at census size N (NaN without haplotype moments)
_noise_c(m::Moments, N) = m.cov isa DemeCov ? _cmean(m.cov, m.idx, N) : NaN

# the variance of the innovation of the offspring's mean genotype, at swap
function _innovation_var(m::Moments, Csel)
    c = m.cov
    if c isa DemeCov && c.drift
        return c.ploidy^2 * max(_cmean(c, m.idx, m.Noff) - Csel, 0.0)
    end
    # without haplotype moments (:spectrum, :aggregate): the leading order, V_Gs/n + V_m/N
    m.VGs / max(m.nsurv, 1.0) + QTLSim.muteVar(m.p) / max(m.Noff, 1.0)
end

# the variance of the mean genotype about the deterministic one at the start of a
# generation: Ploidy²·C, or without haplotype moments the stationary variance of the AR(1)
# recursion d′ = (1 − β)d + ε, σ²/(β(2 − β)), from the previous generation's selection
function _initial_var(m::Moments)
    c = m.cov
    c isa DemeCov && return c.drift ? c.ploidy^2 * _cmean(c, m.idx, m.N) : 0.0
    β = m.VG * m.sinv
    0 < β < 2 || return 0.0
    σ2 = m.VGs / max(m.nsurv, 1.0) + QTLSim.muteVar(m.p) / max(m.N, 1.0)
    σ2 / (β * (2 - β))
end

# add δ to the within-population genic variance of the deme of m (DemeCov), clamped
function _vg_kick!(m::Moments, δ, N)
    c = m.cov
    c isa DemeCov && c.drift || return m
    a = m.idx
    C1 = QTLSim.hap_cmean(c.T[a, 1], c.R[a, 1], c.Q[a, a, 1], _nd(c, N), c.ploidy)
    g = c.T[a, 1] - C1                         # the genic within-haplotype variance
    d = c.ploidy == 2 ? δ / 2 : δ
    c.T[a, 1] += max(d, -0.9 * max(g, 0.0))
    m
end

function moments_replicates(p::Params; replicates::Integer, rng::AbstractRNG, from::Integer,
                            VG::Real = 0.0, zbar::Real = 0.0, loci::Symbol = :auto, scale::Real = 1.0,
                            spectrum = (;), vg = nothing)
    QTLSim.checkparams(p)
    p.NoTraits == 1 || throw(ArgumentError("moments_replicates is for one trait; NoTraits = $(p.NoTraits)"))
    replicates >= 1 || throw(ArgumentError("replicates must be at least 1, got $replicates"))
    0 <= from < p.NoGenerations || throw(ArgumentError("from must be in 0:$(p.NoGenerations - 1), got $from"))
    scale >= 0 || throw(ArgumentError("scale must be non-negative, got $scale"))
    cycle = lifecycle().generation
    check_loci(loci, false)
    loci, α2 = resolve_loci([p], cycle, loci)
    loci === :rare_local && (loci = :rare_alleles)        # one deme: the same closure
    noisy = Ref(false)
    Csel = Ref(NaN)
    images = Dict{Symbol,Any}(
        :select => Operator(:select, (m, _) -> (moments_select!(m); noisy[] && (Csel[] = _noise_c(m, m.N)); m)),
        :swap => Operator(:swap, function (m, _)
            if noisy[] && !m.extinct
                v = scale^2 * _innovation_var(m, Csel[])
                v > 0 && (m.zoff += sqrt(v) * randn(rng))
                vg === nothing || vg.innovation <= 0 || _vg_kick!(m, scale * sqrt(vg.innovation) * randn(rng), m.Noff)
            end
            moments_swap!(m)
        end))
    rec = records()
    step = lower(cycle, :moments; images, record = rec)
    cov = loci === :aggregate ? nothing : loci === :spectrum ? spectrum_state(p, cycle, VG; spectrum...) :
          DemeCov(p, [float(p.PopSize)], [float(VG)]; α2, drift = true)
    m = Moments(p; zbar, VG, erosion = true, cov)
    function advance!(m, stop)
        while m.gen < stop
            it = m.gen
            step(m, nothing)
            m.extinct && return m
            m.gen = max(m.gen, it + 1)
        end
        m
    end
    m.gen = 0
    advance!(m, from)
    det = (extinct = m.extinct, extinct_gen = m.extinct ? m.gen : -1, map(copy, rec)...)
    m.extinct && return (; det, from, replicates = NamedTuple[])
    snap = deepcopy(m)
    noisy[] = true
    reps = map(1:replicates) do _
        foreach(empty!, rec)
        mr = deepcopy(snap)
        v = scale^2 * _initial_var(mr)
        v > 0 && (mr.zbar += sqrt(v) * randn(rng))
        if vg !== nothing && vg.stationary > 0 && mr.cov isa DemeCov
            _vg_kick!(mr, scale * sqrt(vg.stationary) * randn(rng), mr.N)
            mr.VG = cov_vpop(mr, mr.cov)
        end
        advance!(mr, p.NoGenerations)
        (extinct = mr.extinct, extinct_gen = mr.extinct ? mr.gen : -1, map(copy, rec)...)
    end
    (; det, from, replicates = reps)
end
