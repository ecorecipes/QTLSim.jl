# Tests: theory. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl theory`.

@testset "theory: equilibrium variance, lag, critical rate" begin
    p = Params(PopSize=1000, NoQTL=100, Ploidy=2, PhenotypeVar=0.02, Heritability=0.4,
               MutationScalar=0.004, MutationRatePerGamete=1e-4, InitFitnessSD=2.0)
    @test selection_variance(p) ≈ 2.0^2 + 0.012
    @test round(gaussian_vg(p), digits = 4) == 0.0139           # vignette 02's Gaussian prediction
    V = gaussian_vg_exact(p)                                    # the fixed point of the Gaussian recursion
    @test V^2 ≈ muteVar(p) * (selection_variance(p) + V) && V > gaussian_vg(p)
    # the effective size (validate/effective_size.md; Lean EffectiveSize.lean): by default the breeders themselves, whatever the
    # census; :falconer is the earlier form, 4N/(V_k + 2) with the survivors' V_k, about the harmonic mean of adults and survivors
    @test effective_size(1000, 1000) == 1000 && effective_size(1000, 100) == 100 && effective_size(50, 100) == 100
    @test effective_size(1000, 250) === effective_size(1000, 250; scheme = :breeders) === 250
    @test effective_size(1000, 250.5) == 250.5
    @test effective_size(1000, 250; scheme = :falconer) === 4000 / (2000 * 249 / 250^2 + 2)      # the old default, exactly
    @test effective_size(1000, 1000; scheme = :falconer) ≈ 2 * 1000^2 / (2 * 1000 - 1)             # N + ½, not N, when all survive
    @test effective_size(1000, 250; scheme = :falconer) ≈ 2 * 1000 * 250 / 1250 rtol = 0.01        # ≈ the harmonic mean of N and n
    @test effective_size(1000, 250; scheme = :falconer) > 1.5 * 250                                # 1.6 times the breeders
    @test_throws ArgumentError effective_size(1000, 250; scheme = :census)
    # the neutral two-stage recursion of the life cycle, x′ = (1 − 1/n)x + y/(2n) + u, y′ = x + u (x: two copies in different
    # adults, y: the two copies of one adult): identity by descent rises by 1/(2n) per generation whatever N, the mutation–drift
    # fixed point is x = (2n + 1)u, y = (2n + 2)u, and without mutation it decays by the eigenvalue 1 − 1/(2n + 1) + O(1/n³)
    let n = 250.0, step = (x, y, u) -> ((1 - 1 / n) * x + y / (2n) + u, x + u)
        x, y = (2n + 1) * 1.0, (2n + 2) * 1.0
        @test all(isapprox.(step(x, y, 1.0), (x, y); rtol = 1e-14))                               # the fixed point
        λ = ((1 - 1 / n) + sqrt((1 - 1 / n)^2 + 2 / n)) / 2                                        # root of λ² − (1 − 1/n)λ − 1/(2n)
        @test abs((1 - λ) - 1 / (2n + 1)) < 1 / n^3 && 1 - 1 / (2n + 1) < λ < 1 - 1 / (2n + 2)
        a, b = 1.0, 1.0
        for _ in 1:2000; a, b = step(a, b, 0.0); end
        a2, _ = step(a, b, 0.0)
        @test a2 / a ≈ λ rtol = 1e-12
    end
    # old values exactly under ne = :falconer; fewer genes, less variance, with the breeders
    Ne_old = 4000 / (2000 * 249 / 250^2 + 2)
    let n = p.NoQTL, μ = p.MutationRatePerGamete, Vs = selection_variance(p), α2 = muteVar(p) / (2 * p.NoQTL * p.MutationRatePerGamete)
        hoc(Ne) = 4n * μ * Vs / (1 + Vs / (Ne * α2))                                               # Bürger's stochastic house of cards
        hog(Ne) = (b = α2 / (2Vs) + 1 / (2Ne); 2n * 2μ * α2 / (b + sqrt(b^2 + 4μ * α2 / Vs)))      # Hermisson and Wagner's A4
        @test house_of_cards_vg(p; nsurv = 250, ne = :falconer) ≈ hoc(Ne_old) rtol = 1e-14
        @test house_of_cards_vg(p; nsurv = 250) === house_of_cards_vg(p; nsurv = 250, ne = :breeders) ≈ hoc(250)
        @test house_of_gauss_vg(p; nsurv = 250, ne = :falconer) ≈ hog(Ne_old) rtol = 1e-14
        @test house_of_gauss_vg(p; nsurv = 250) === house_of_gauss_vg(p; nsurv = 250, ne = :breeders) ≈ hog(250)
        @test house_of_gauss_vg(p; nsurv = 250) < house_of_gauss_vg(p; nsurv = 250, ne = :falconer)
        @test house_of_gauss_vg(p) ≈ hog(1000) rtol = 1e-14                                        # all breed: N, not N + ½
        @test QTLSim.burnin_start_vg(p) === house_of_gauss_vg(p; ne = :falconer) ≈ hog(2 * 1000^2 / 1999)   # the emulators' start, kept
    end
    let q = Params(p; InitFitnessSD = 1e5)                      # neutral: Lynch and Hill's 2 N_e V_m with N_e = the breeders
        @test house_of_gauss_vg(q; nsurv = 250) ≈ 2 * 250 * muteVar(q) rtol = 1e-3
        @test house_of_gauss_vg(q; nsurv = 250, ne = :falconer) ≈ 2 * Ne_old * muteVar(q) rtol = 1e-3
        # the moments lowering's drift is that of the n parents: its neutral fixed point is 2 n V_m, not 2·effective_size·V_m
        st = QTLSim.hap_init(0.0, 1000.0, 2, 1)
        for _ in 1:20_000
            st = QTLSim.hap_mutate(QTLSim.hap_reproduce(st, 250.0, [0.0], 2), muteVar(q), 2)
        end
        @test QTLSim.hap_vpop(st, 1000.0, 2) * 1000 / 999 ≈ 2 * 250 * muteVar(q) rtol = 1e-6
    end
    @test house_of_cards_vg(p) < house_of_cards_vg(Params(p; NoQTL = 1000)) < 4 * 1000 * 1e-4 * selection_variance(p)
    @test steady_lag(0.01; VG = 0.5, Vs = 1.5) ≈ 0.04
    @test steady_lag(Params(p; MoveFitnessSD = 0.01); VG = 0.5) ≈ 0.01 * (0.5 + selection_variance(p)) / 0.5
    kc(VG) = critical_rate(; n = 500, VG, ω = 5.0, VE = 0.01, sd50 = 0.1, N = 1000)
    @test 0 < kc(0.01) < kc(0.02) && isnan(critical_rate(; n = 2, VG = 0.01, ω = 5.0, VE = 0.01, sd50 = 0.1, N = 1000))
    Γ = QTLSim.Distributions.Gamma(1.201, 0.004635)
    m = [0.001, 0.004, 0.01]; r = [0.02, 0.03, 0.04]
    a = mutation_adjusted(r, m, Γ)
    @test mean(a ./ r) ≈ 1                                      # weights relative to their mean
    @test a ./ r ≈ QTLSim.Distributions.pdf.(Γ, m) ./ mean(QTLSim.Distributions.pdf.(Γ, m))
end

@testset "the moments lowering: the theory as a third functor" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=3000, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    # the default, loci = :auto, is the rare-alleles closure for one diploid trait with
    # mutations, and :multilocus (the multilocus Gaussian model) on request
    o = run_moments(p)
    @test QTLSim.moments_loci(p) === :rare_alleles && QTLSim.moments_loci(p, :multilocus) === :multilocus
    @test o.geno_var == run_moments(p; loci = :rare_alleles).geno_var == run_moments(p; loci = :auto).geno_var
    @test o.geno_var[end] ≈ rare_alleles_vg(p) rtol = 1e-6           # converges to the rare-alleles fixed point
    om = run_moments(p; loci = :multilocus)
    @test om.geno_var[end] ≈ gaussian_vg_loci(p) rtol = 1e-6         # and the multilocus model to its own
    # one diploid locus: half of selection's reduction (between the haplotypes) is not inherited
    @test gaussian_vg_loci(p) ≈ QTLSim.gaussian_vg_closed(2muteVar(p), selection_variance(p)) > gaussian_vg_exact(p)
    # the legacy aggregate recursion converges to gaussian_vg_exact, which the multilocus
    # model is for one haploid locus (or haploid clones of any number of loci)
    @test run_moments(p; loci = :aggregate).geno_var[end] ≈ gaussian_vg_exact(p) rtol = 1e-6
    ph = Params(p; Ploidy = 1, NoQTL = 20)
    @test QTLSim.moments_loci(ph) === :rare_alleles                   # haploids default to the closure for clones
    @test run_moments(ph).geno_var == run_moments(ph; loci = :rare_alleles).geno_var
    @test run_moments(ph; loci = :multilocus).geno_var ≈ run_moments(ph; loci = :aggregate).geno_var rtol = 1e-9
    @test gaussian_vg_loci(ph) ≈ gaussian_vg_exact(ph)
    @test_throws ArgumentError run_moments(p; loci = :other)
    @test_throws ArgumentError run_moments(p; loci = :aggregate, between = true)
    @test o.pheno_var[end] ≈ o.geno_var[end] + envVar(p)
    @test all(0 .< o.mean_fitness .<= 1) && o.nsurv[1] == round(Int, 1000 * o.mean_fitness[1])
    # a continuous move: the lag converges to Lande's steady lag
    pc = Params(p; MoveFitnessSDGen = 1000, MoveFitnessSD = 0.001, MoveType = :continuous)
    oc = run_moments(pc)
    @test oc.theta[end] - oc.pheno_mean[end] ≈ steady_lag(pc; VG = oc.geno_var[end]) rtol = 1e-6
    @test !oc.extinct
    # faster than the population can follow: deterministic extinction after the move
    ox = run_moments(Params(pc; MoveFitnessSD = 0.1))
    @test ox.extinct && ox.extinct_gen > 1000
    # drift: no selection, the variance grows by m·V_E per generation
    od = run_moments(Params(p; Drift = true, NoGenerations = 10))
    @test od.geno_var ≈ muteVar(p) .* (0:9)
    # the lowering is a functor: the generation's image is the composite of the steps' images
    steps = (:develop, :select, :reproduce, :record, :rescale, :mutate, :swap)
    Ext_ = Ext
    for cov in (nothing, Ext_.DemeCov(p, [1000.0], [0.003]; drift = false), Ext_.DemeCov(p, [1000.0], [0.003]))
        m1 = Ext_.Moments(p; VG = 0.003, cov); m2 = deepcopy(m1)
        lower(g, :moments)(m1, nothing)
        foldl(Catlab.compose, (lower(getfield(L, s), :moments) for s in steps))(m2, nothing)
        @test (m1.zbar, m1.VG, m1.θ, m1.gen) == (m2.zbar, m2.VG, m2.θ, m2.gen)
    end
    @test lower(g ⊗ g ⋅ L.migrate(2), :moments) isa Ext_.DemeOp        # several demes (see below)
    @test_throws ArgumentError run_moments(p; cycle = L.develop)
end

@testset "the multilocus Gaussian model: genic variance, LD, recombination, drift" begin
    L = lifecycle(); g = L.generation
    p = Params(PopSize = 1000, NoQTL = 50, Ploidy = 2, NoGenerations = 20_000, PhenotypeVar = 1.0, Heritability = 0.5,
               MutationScalar = 0.01, MutationRatePerGamete = 0.01)
    Vm, Vs = muteVar(p), selection_variance(p)
    # without drift: V²/(V_s + V) = 2L·V_m, whatever the crossover model
    V = gaussian_vg_loci(p)
    @test V^2 / (Vs + V) ≈ 2 * 50 * Vm
    @test V ≈ gaussian_vg_exact(Params(p; MutationScalar = 100 * p.MutationScalar, InitFitnessSD = initFitSD(p))) rtol = 1e-12
    for q in (p, Params(p; CrossoverModel = :all), Params(p; RecombPoissonMean = 5.0))
        @test run_moments(Params(q; NoGenerations = 100_000, Sample = 100); loci = :multilocus).geno_var[end] ≈ V rtol = 1e-6
        e = multilocus_equilibrium(q; loci = :multilocus)
        @test e.VG ≈ V rtol = 1e-6
        @test e.ld < 0 && e.genic ≈ e.VG - e.ld
        # the LD at equilibrium: −ω_c b/ρ_c per haplotype in each pair class, b = L V_m/2
        ω, ρ = QTLSim.pair_classes(q)
        @test sum(ω) ≈ 1 && ω[1] == 1 / 50 && ρ[1] == 0
        @test sum(ω[2:end] .* ρ[2:end]) / (1 - 1 / 50) ≈ pair_recombination(q) rtol = 1e-10
        @test e.ld ≈ -50 * Vm * sum(ω[2:end] ./ ρ[2:end]) rtol = 1e-6
    end
    # less recombination, more LD, the same V_G
    @test multilocus_equilibrium(Params(p; RecombPoissonMean = 0.5); loci = :multilocus).ld < multilocus_equilibrium(p; loci = :multilocus).ld
    # drift: erosion = true and between = true (one deme) are the same recursion, and its
    # fixed point is gaussian_vg_loci(p; drift = true), below the deterministic one
    a, b = run_moments(p; erosion = true, loci = :multilocus), run_moments(p; between = true, loci = :multilocus)
    @test a.geno_var == b.geno_var
    @test a.geno_var[end] ≈ gaussian_vg_loci(p; drift = true) rtol = 1e-6
    @test run_moments(p; erosion = true).geno_var == run_moments(p; between = true).geno_var    # also by default
    @test gaussian_vg_loci(p; drift = true) < V
    # selection: the genic reduction is the sum of the loci's shares, (H/L)²/s each, and the
    # survivors' variance is V − V²/s
    H, s, Lq = 0.7, 3.0, 50
    st = QTLSim.hap_init(2H, Inf, 2, 2)
    st2, β = QTLSim.hap_select(st, 1 / s, Inf, [1 / Lq, 1 - 1 / Lq], 2)
    @test st[1][1] - st2[1][1] ≈ sum((H / Lq)^2 / s for _ in 1:Lq) && β ≈ 2H / s
    @test QTLSim.hap_vpop(st2, Inf, 2) ≈ 2H - (2H)^2 / s
    # recombination turns the LD back into variance: one generation of mating at ρ = 1
    # restores the haplotypes' LE, losing only the between-haplotype part
    st3 = QTLSim.hap_reproduce(st2, Inf, [0.0, 1.0], 2)
    @test st3[1][2] ≈ st2[2][2] && st3[1][1] == st2[1][1]
    # pair classes: one locus has only the genic class; two loci are split by one point
    @test QTLSim.pair_classes(Params(p; NoQTL = 1)) == ([1.0], [0.0])
    ωa, ρa = QTLSim.pair_classes(Params(p; NoQTL = 2, RecombPoissonMean = 0.7))
    @test ωa ≈ [0.5, 0.5] && ρa ≈ [0.0, 0.7 * exp(-0.7)]
    # a continuous move: Lande's steady lag with the multilocus V_G
    pc = Params(p; NoGenerations = 30_000, MoveFitnessSDGen = 20_000, MoveFitnessSD = 0.01, MoveType = :continuous)
    oc = run_moments(pc; loci = :multilocus)
    @test oc.theta[end] - oc.pheno_mean[end] ≈ steady_lag(pc; VG = oc.geno_var[end]) rtol = 1e-6
    # several demes without migration are independent copies of one, drift or not
    for erosion in (false, true)
        @test run_moments(Params(p; NoGenerations = 300); cycle = g ⊗ g, erosion).geno_var[1:2:end] ==
              run_moments(Params(p; NoGenerations = 300); erosion).geno_var
    end
    @test_throws ArgumentError run_moments([p, Params(p; Ploidy = 1)]; cycle = g ⊗ g)
    # haploid drift is clonal: V_G loses 1/n per generation, not 1/(2n)
    ph = Params(p; Ploidy = 1, NoQTL = 1, Drift = true, NoGenerations = 30, MutationScalar = 0.0, InitFitnessSD = 1.0)
    oh = run_moments(ph; VG = 1.0, erosion = true)
    @test oh.geno_var ≈ [(1 - 1 / 1000)^t * 1000 / 999 for t in 0:29]      # expected sample variances
end
