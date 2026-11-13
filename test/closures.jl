# Tests: closures. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl closures`.

@testset "rare alleles: house of Gauss, house of cards, Gaussian and neutral limits" begin
    p = Params(PopSize = 1000, NoQTL = 10, Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5,
               MutationScalar = 0.01, MutationRatePerGamete = 1e-4, RecombPoissonMean = 2.0)
    Vm, Vs, L = muteVar(p), selection_variance(p), p.NoQTL
    @test QTLSim.mutation_effect_variance(p) ≈ Vm / (2L * 1e-4)
    @test_throws ArgumentError QTLSim.mutation_effect_variance(Params(p; MutationRatePerGamete = 0.0))
    # the closed form: rare mutations → Bürger's stochastic house of cards, and without drift
    # Turelli's 4Lμ V_s; many small ones → the Gaussian √(2L V_m V_s); V_s → ∞ → 2N_e V_m
    rare = Params(p; MutationRatePerGamete = 1e-7, MutationScalar = 0.01 * 1e-7 / 1e-4,     # the same α² and V_s, rarer
                  InitFitnessSD = initFitSD(p))
    @test house_of_gauss_vg(rare; nsurv = 950) ≈ house_of_cards_vg(rare; nsurv = 950) rtol = 1e-3
    @test house_of_gauss_vg(rare; drift = false) ≈ 4L * 1e-7 * Vs rtol = 1e-3
    many = Params(p; MutationRatePerGamete = 1.0)
    @test house_of_gauss_vg(many; drift = false) ≈ sqrt(2L * Vm * Vs) rtol = 1e-3
    neutral = Params(p; InitFitnessSD = 1e5)
    @test house_of_gauss_vg(neutral; nsurv = 1000) ≈ 2effective_size(1000, 1000) * Vm rtol = 1e-3
    @test house_of_gauss_vg(p) < house_of_cards_vg(p) && house_of_gauss_vg(p) < gaussian_vg_loci(p; drift = true)
    # the multilocus model with the rare-alleles closure: below the Gaussian, equal to it
    # when mutations are many and small, Turelli's house of cards when they are rare
    @test rare_alleles_vg(p) < gaussian_vg_loci(p)
    @test rare_alleles_vg(many) ≈ gaussian_vg_loci(many) rtol = 2e-3
    @test rare_alleles_vg(rare) ≈ 4L * 1e-7 * Vs rtol = 0.02
    @test rare_alleles_vg(rare; drift = true) ≈ house_of_cards_vg(rare; nsurv = 1000 * multilocus_equilibrium(rare; drift = true).Wbar) rtol = 0.03
    # V_s → ∞: both models are the neutral recursion, 2N V_m
    @test rare_alleles_vg(neutral; drift = true) ≈ gaussian_vg_loci(neutral; drift = true) rtol = 1e-3
    @test rare_alleles_vg(neutral; drift = true) ≈ 2 * 1000 * Vm rtol = 0.01
    # the moments lowering converges to it, without and with drift; it is the default for
    # one diploid trait with mutations, and the multilocus model stays available
    q = Params(p; NoGenerations = 100_000, Sample = 1000)
    @test run_moments(q; loci = :rare_alleles).geno_var[end] ≈ rare_alleles_vg(p) rtol = 1e-6
    @test run_moments(q; loci = :rare_alleles, erosion = true).geno_var[end] ≈ rare_alleles_vg(p; drift = true) rtol = 1e-6
    @test run_moments(q).geno_var[end] ≈ rare_alleles_vg(p) rtol = 1e-6
    @test run_moments(q; erosion = true).geno_var == run_moments(q; loci = :rare_alleles, erosion = true).geno_var
    @test run_moments(q; loci = :multilocus).geno_var[end] ≈ gaussian_vg_loci(p) rtol = 1e-6
    @test multilocus_equilibrium(p).VG == multilocus_equilibrium(p; loci = :rare_alleles).VG == rare_alleles_vg(p)
    @test multilocus_equilibrium(p; loci = :multilocus).VG ≈ gaussian_vg_loci(p) rtol = 1e-6
    # how the default resolves (QTLSim.moments_loci): the closure with locally favoured alleles
    # for several traits (the closure itself for one deme), the closure for one trait, haploid too
    @test QTLSim.moments_loci(p) === QTLSim.moments_loci(p, :auto) === :rare_alleles
    @test QTLSim.moments_loci(p, :aggregate) === :aggregate
    p2t = Params(p; NoTraits = 2, NoGenerations = 200)
    @test QTLSim.moments_loci(p2t) === :rare_local
    @test run_moments(p2t).gmatrix == run_moments(p2t; loci = :rare_alleles).gmatrix
    @test multilocus_equilibrium(p2t).VG == multilocus_equilibrium(p2t; loci = :rare_alleles).VG
    g2_ = lifecycle().generation ⊗ lifecycle().generation
    @test run_moments([p2t, p2t]; cycle = g2_).gmatrix == run_moments([p2t, p2t]; cycle = g2_, loci = :rare_local).gmatrix
    @test run_moments([p2t, p2t]; cycle = g2_, loci = :rare_alleles).gmatrix[1:2:end] ≈ run_moments(p2t; loci = :rare_alleles).gmatrix rtol = 1e-12
    @test QTLSim.moments_loci(Params(p; Ploidy = 1)) === :rare_alleles
    ph = Params(p; Ploidy = 1, NoGenerations = 200)
    @test run_moments(ph).geno_var == run_moments(ph; loci = :rare_alleles).geno_var
    @test multilocus_equilibrium(ph).VG == multilocus_equilibrium(ph; loci = :rare_alleles).VG == rare_alleles_vg(ph)
    @test QTLSim.moments_loci(Params(p; Ploidy = 1, NoTraits = 2)) === :rare_local
    @test_throws ArgumentError run_moments(q; loci = :bogus)
    @test run_moments(Params(p2t; NoGenerations = 50_000); loci = :rare_alleles).gmatrix[end] ≈ rare_alleles_G(p2t) rtol = 1e-6
    # demes: the closure acts in each deme
    L_ = lifecycle()
    g2 = L_.generation ⊗ L_.generation
    m2 = run_moments(Params(p; NoGenerations = 2000, Sample = 100); cycle = g2, loci = :rare_alleles)
    m1 = run_moments(Params(p; NoGenerations = 2000, Sample = 100); loci = :rare_alleles)
    @test m2.geno_var[m2.deme .== 1] ≈ m1.geno_var && m2.geno_var[m2.deme .== 2] ≈ m1.geno_var
    @test run_moments(Params(p; NoGenerations = 2000, Sample = 100); cycle = g2).geno_var == m2.geno_var    # the default
    # each deme has its own α²: demes of different mutation rates are independent copies of one
    pa, pb = Params(p; NoGenerations = 2000, Sample = 100), Params(p; NoGenerations = 2000, Sample = 100, MutationRatePerGamete = 1e-3)
    md = run_moments([pa, pb]; cycle = g2, between = true)
    @test md.geno_var[md.deme .== 1] ≈ run_moments(pa; between = true).geno_var rtol = 1e-12
    @test md.geno_var[md.deme .== 2] ≈ run_moments(pb; between = true).geno_var rtol = 1e-12

    # mutate steps that carry mutation fields: their α² (one law of effects), or the multilocus
    # model when they differ (an explicit :rare_alleles is then an error)
    gm = withparams(L_.generation; MutationRatePerGamete = 1e-3)
    @test run_moments(pa; cycle = gm).geno_var == run_moments(pa; cycle = gm, loci = :rare_alleles).geno_var
    @test run_moments(pa; cycle = gm).geno_var ≈ run_moments(pb).geno_var rtol = 1e-12
    gmix = gm ⋅ L_.generation
    @test run_moments(pa; cycle = gmix).geno_var == run_moments(pa; cycle = gmix, loci = :multilocus).geno_var
    @test_throws ArgumentError run_moments(pa; cycle = gmix, loci = :rare_alleles)
end

@testset "rare alleles: several traits and haploids (src/rare_alleles_traits.jl)" begin
    L_ = lifecycle(); g = L_.generation
    p = Params(PopSize = 1000, NoQTL = 10, Ploidy = 2, PhenotypeVar = 1.0, Heritability = 0.5, NoGenerations = 4000,
               MutationScalar = 0.01, MutationRatePerGamete = 1e-4, RecombPoissonMean = 2.0, Sample = 10)
    Vm, Vs, L, μ = muteVar(p), selection_variance(p), p.NoQTL, p.MutationRatePerGamete
    # n = 1: the matrix closure with 1×1 matrices is the one-trait closure
    @test QTLSim.mutation_effect_covariance(p) ≈ fill(QTLSim.mutation_effect_variance(p), 1, 1) rtol = 1e-14
    @test QTLSim.rare_alleles_rate(fill(0.3, 1, 1), fill(1 / 2.5, 1, 1))[1] ≈ 0.3 / 5 rtol = 1e-14
    for erosion in (false, true)
        a = run_moments(p; loci = :rare_alleles, erosion).geno_var
        b = Ext.run_trait_moments(p; loci = :rare_alleles, erosion).geno_var
        @test maximum(abs.(a .- b) ./ a[end]) < 1e-12
    end
    @test rare_alleles_G(p)[1] ≈ rare_alleles_vg(p) rtol = 1e-12
    @test rare_alleles_G(p; drift = true)[1] ≈ rare_alleles_vg(p; drift = true) rtol = 1e-12
    @test house_of_gauss_G(p)[1] ≈ house_of_gauss_vg(p) rtol = 1e-12
    @test house_of_gauss_G(p; drift = false)[1] ≈ house_of_gauss_vg(p; drift = false) rtol = 1e-12
    @test house_of_cards_G(p)[1] ≈ 4L * μ * Vs rtol = 1e-12
    # X = E[a aᵀ/(aᵀ s⁻¹ a)]: tr(s⁻¹X) = 1, X = s for one trait, the closed form for two
    # against Monte Carlo, and spherical pleiotropy dilutes the house of cards by 1/n
    C2, s2 = [1.0 0.6; 0.6 0.5], [2.0 -0.5; -0.5 1.0]
    X2 = QTLSim.rare_hc_matrix(C2, inv(s2))
    @test QTLSim.tr(inv(s2) * X2) ≈ 1 rtol = 1e-12
    @test QTLSim.rare_hc_matrix(fill(0.2, 1, 1), fill(1 / 3, 1, 1))[1] ≈ 3 rtol = 1e-12
    rng = MersenneTwister(11); Lc = QTLSim.cholesky(C2).L
    mc = sum(begin a = Lc * randn(rng, 2); a * a' ./ QTLSim.dot(a, s2 \ a) end for _ in 1:200_000) ./ 200_000
    @test mc ≈ X2 rtol = 0.01
    C3, s3 = [1.0 0.3 0.1; 0.3 0.8 -0.2; 0.1 -0.2 0.5], [1.5 0.2 0; 0.2 1 0.3; 0 0.3 2]
    @test QTLSim.tr(s3 \ QTLSim.rare_hc_matrix(C3, inv(s3))) ≈ 1 rtol = 1e-9
    @test QTLSim.rare_hc_rates([0.4, 0.4, 0.4]) ≈ fill(3 * 0.4 / 2, 3) rtol = 1e-9          # e = 1/(nλ)
    @test QTLSim.rare_hc_rates([0.2, 0.5]) ≈ [sqrt(0.2) * (sqrt(0.2) + sqrt(0.5)), sqrt(0.5) * (sqrt(0.2) + sqrt(0.5))] ./ 2
    sph = Params(p; NoTraits = 2)                                                          # MutCorr = SelCorr = 0
    @test house_of_cards_G(sph) ≈ (4L * μ * Vs / 2) .* [1 0; 0 1] rtol = 1e-12
    # two traits whose M and Ω do not share eigenvectors: the limits of the closed form and
    # of the closure in the multilocus model
    p2 = Params(p; NoTraits = 2, MutCov = Vm .* [1.0 0.5; 0.5 0.6], SelCorr = -0.4)
    A, M = selection_variance_matrix(p2), mutational_covariance(p2)
    rare = Params(p2; MutationRatePerGamete = 1e-7, MutationScalar = 0.01 * 1e-7 / 1e-4, MutCov = 1e-3 .* M,
                  InitFitnessSD = initFitSD(p2))                  # the same C and A, rarer
    @test house_of_gauss_G(rare; drift = false) ≈ house_of_cards_G(rare) rtol = 1e-3
    @test rare_alleles_G(rare) ≈ house_of_cards_G(rare) rtol = 0.02
    many = Params(p2; MutationRatePerGamete = 1.0)
    @test rare_alleles_G(many) ≈ gaussian_G_loci(many) rtol = 2e-3
    weak = Params(many; InitFitnessSD = 100.0)                   # G ≪ A: G A⁻¹ G = 2L M
    Gw = house_of_gauss_G(weak; drift = false)
    @test Gw * (selection_variance_matrix(weak) \ Gw) ≈ 2L .* mutational_covariance(weak) rtol = 1e-3
    neutral = Params(p2; InitFitnessSD = 1e5)
    @test house_of_gauss_G(neutral; nsurv = 1000) ≈ 2effective_size(1000, 1000) .* M rtol = 1e-3
    @test rare_alleles_G(neutral; drift = true) ≈ gaussian_G_loci(neutral; drift = true) rtol = 1e-3
    G2 = rare_alleles_G(p2)
    @test all(QTLSim.eigvals(G2) .> 0) && G2[1, 1] < gaussian_G_loci(p2)[1, 1]
    # the trait moments lowering converges to it, with and without drift; explicit
    # loci = :rare_alleles with traits, on several demes and through mutate steps
    q2 = Params(p2; NoGenerations = 60_000, Sample = 1000)
    @test run_moments(q2; loci = :rare_alleles).gmatrix[end] ≈ G2 rtol = 1e-6
    @test run_moments(q2; loci = :rare_alleles, erosion = true).gmatrix[end] ≈ rare_alleles_G(p2; drift = true) rtol = 1e-6
    gm = withparams(g, p2; MutationRatePerGamete = 1e-3)
    pb = Params(p2; MutationRatePerGamete = 1e-3, NoGenerations = 2000)
    @test run_moments(Params(p2; NoGenerations = 2000); cycle = gm, loci = :rare_alleles).gmatrix ≈
          run_moments(pb; loci = :rare_alleles).gmatrix rtol = 1e-12
    @test_throws ArgumentError run_moments(Params(p2; NoGenerations = 200); cycle = gm ⋅ g, loci = :rare_alleles)
    @test run_moments(Params(p2; NoGenerations = 200); cycle = gm ⋅ g).gmatrix ==
          run_moments(Params(p2; NoGenerations = 200); cycle = gm ⋅ g, loci = :multilocus).gmatrix
    md = run_moments([pb, pb]; cycle = g ⊗ g, between = true, loci = :rare_alleles)
    @test md.gmatrix[md.deme .== 1] ≈ run_moments(pb; between = true, loci = :rare_alleles).gmatrix rtol = 1e-12
    # haploids: clones are one unit of inheritance (QTLSim.hap_clonal), with U = Lμ
    ph = Params(p; Ploidy = 1)
    @test QTLSim.hap_clonal(1, [0.0, 0.3]) && QTLSim.hap_clonal(2, [0.0]) && !QTLSim.hap_clonal(2, [0.0, 0.3])
    phr = Params(rare; Ploidy = 1, NoTraits = 1, MutCov = zeros(0, 0), SelCorr = 0.0)
    Vsh = selection_variance(phr)
    @test house_of_gauss_vg(phr; drift = false) ≈ 2L * 1e-7 * Vsh rtol = 1e-3        # half the diploid's 4Lμ V_s
    @test rare_alleles_vg(phr) ≈ 2L * 1e-7 * Vsh rtol = 0.02
    phm = Params(ph; MutationRatePerGamete = 1.0)
    @test rare_alleles_vg(phm) ≈ gaussian_vg_loci(phm) rtol = 2e-3                     # √(V_m V_s): one clone
    @test house_of_gauss_vg(Params(ph; InitFitnessSD = 1e5); nsurv = 1000) ≈ effective_size(1000, 1000) * Vm rtol = 1e-3
    qh = Params(ph; NoGenerations = 100_000, Sample = 1000)
    oh = run_moments(qh; loci = :rare_alleles)
    @test oh.geno_var[end] ≈ rare_alleles_vg(ph) rtol = 1e-6
    @test run_moments(qh; loci = :rare_alleles, erosion = true).geno_var[end] ≈ rare_alleles_vg(ph; drift = true) rtol = 1e-6
    @test rare_alleles_vg(ph) ≈ house_of_gauss_vg(ph; drift = false) rtol = 0.1          # LD carries little for clones
    # several traits, haploid
    p2h = Params(p2; Ploidy = 1)
    @test house_of_gauss_G(p2h)[1, 1] < house_of_gauss_G(p2)[1, 1]
    @test run_moments(Params(p2h; NoGenerations = 60_000, Sample = 1000); loci = :rare_alleles).gmatrix[end] ≈ rare_alleles_G(p2h) rtol = 1e-6
end

@testset "allele-frequency spectrum (loci = :spectrum): neutral, equilibrium, moving-optimum limits" begin
    p = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 0.022, Heritability = 0.32, MutationScalar = 0.004,
               MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0, InitFitnessSD = 5.0)
    Vm = muteVar(p)
    # no selection: the neutral spectrum holds 2N V_m (Lynch and Hill 1986), and a move changes nothing
    e0 = spectrum_equilibrium(Params(p; Drift = true))
    @test e0.VG ≈ 2 * 1000 * Vm rtol = 0.01
    sp = e0.spectrum
    @test QTLSim.spectrum_select(sp, 3.0, Inf) == 0.0 && QTLSim.spectrum_linked(sp, 1000) == 0.0
    # no move: the lowering converges to the stationary spectrum, which for rare mutations in
    # linkage equilibrium is close to Hermisson and Wagner's stochastic house of Gauss
    e = spectrum_equilibrium(p)
    @test e.VG == spectrum_vg(p)
    @test e.VG ≈ house_of_gauss_vg(p; nsurv = e.nparents) rtol = 0.05
    o = run_moments(Params(p; NoGenerations = 30_000, Sample = 1000); loci = :spectrum)
    @test o.geno_var[end] ≈ e.VG rtol = 2e-3
    @test !o.extinct && all(o.theta .== 0)
    # the start of a move: the stationary spectrum is symmetric in the effects, so selection
    # moves the mean by V_G λ/S exactly (no skew), Lande's Gaussian response
    S = initFitSD(p)^2 + envVar(p) + QTLSim.spectrum_vg(e.spectrum)
    @test QTLSim.spectrum_select(e.spectrum, 0.7, S) ≈ QTLSim.spectrum_vg(e.spectrum) * 0.7 / S rtol = 1e-8
    # a population started at V_G takes the stationary shape scaled to V_G
    @test run_moments(Params(p; NoGenerations = 1, Sample = 1); loci = :spectrum, VG = 0.05).geno_var[1] ≈ 0.05 * 1000 / 999 rtol = 1e-10
    # a fast move: the lag decides, within a few hundred generations, and the time to extinction
    # is the Gaussian closure's (rare alleles, with drift) to within 10%
    pm = Params(p; NoGenerations = 7000, Sample = 1, MoveFitnessSDGen = 5000, MoveType = :continuous, MoveFitnessSD = 0.2)
    ms, mr = run_moments(pm; loci = :spectrum), run_moments(pm; loci = :rare_alleles, erosion = true)
    @test ms.extinct && mr.extinct
    @test abs(ms.extinct_gen - mr.extinct_gen) <= 0.1 * (mr.extinct_gen - 5000)
    # a slow move in the house-of-cards regime: directional selection drives the rare alleles of
    # large effect to intermediate frequency and V_G rises (Barton and Turelli 1987), with or
    # without the linked-selection term; the rare-alleles closure keeps its stationary value
    pH = Params(PopSize = 1000, Ploidy = 2, NoQTL = 100, PhenotypeVar = 1.0, Heritability = 0.5, MutationScalar = 0.01,
                MutationRatePerGamete = 1e-4, RecombPoissonMean = 2.0, NoGenerations = 16_000, Sample = 10,
                MoveFitnessSDGen = 10_000, MoveType = :continuous, MoveFitnessSD = 0.01)
    for lk in (:barton, :none)
        h = run_moments(pH; loci = :spectrum, spectrum = (linkage = lk,))
        @test h.geno_var[end] > 1.3 * h.geno_var[findfirst(==(9990), h.gen)]
    end
    hr = run_moments(pH; loci = :rare_alleles, erosion = true)
    @test hr.geno_var[end] ≈ hr.geno_var[findfirst(==(9990), hr.gen)] rtol = 0.02
    # linked selection (Barton 1995): N/N_e − 1 grows as b² V_G, nothing with free recombination
    spl = QTLSim.AlleleSpectrum(p); spl.n .= sp.n
    QTLSim.spectrum_select(spl, 2.0, 25.0); c2 = QTLSim.spectrum_linked(spl, 1000)
    QTLSim.spectrum_select(spl, 4.0, 25.0); c4 = QTLSim.spectrum_linked(spl, 1000)
    @test c2 > 0 && c4 ≈ 4c2
    spf = QTLSim.AlleleSpectrum(Params(p; CrossoverModel = :all, RecombPoissonMean = 50.0)); spf.n .= sp.n
    QTLSim.spectrum_select(spf, 2.0, 25.0)
    @test QTLSim.spectrum_linked(spf, 1000) < 0.01c2
    # the slow-move collapse: under QTSim's crossover the variance falls more than with free recombination
    pv = Params(p; NoGenerations = 6000, Sample = 10, MoveFitnessSDGen = 5000, MoveType = :continuous, MoveFitnessSD = 0.02)
    vl = run_moments(pv; loci = :spectrum).geno_var[end]
    vf = run_moments(pv; loci = :spectrum, spectrum = (linkage = :none,)).geno_var[end]
    @test vl < 0.5vf
    # limits of the option
    g2 = lifecycle().generation ⊗ lifecycle().generation
    @test_throws ArgumentError run_moments([p, p]; cycle = g2, loci = :spectrum)
    @test_throws ArgumentError run_moments(Params(p; NoGenerations = 10); loci = :spectrum, between = true)
    @test_throws ArgumentError run_moments(Params(p; NoGenerations = 10, Ploidy = 1); loci = :spectrum)
    @test_throws ArgumentError run_moments(Params(p; NoGenerations = 10, NoTraits = 2); loci = :spectrum)
    @test_throws ArgumentError QTLSim.AlleleSpectrum(p; linkage = :bogus)
    # one law of effects from the mutate steps: their fields are the spectrum's
    gm = withparams(lifecycle().generation; MutationRatePerGamete = 1e-3)
    q = Params(p; NoGenerations = 300, Sample = 50)
    @test run_moments(q; cycle = gm, loci = :spectrum).geno_var ≈ run_moments(Params(q; MutationRatePerGamete = 1e-3); loci = :spectrum).geno_var
    @test_throws ArgumentError run_moments(q; cycle = gm ⋅ lifecycle().generation, loci = :spectrum)
end

@testset "rare alleles with locally favoured alleles (loci = :rare_local, src/local_alleles.jl)" begin
    L_ = lifecycle()
    p = Params(Ploidy = 2, NoQTL = 20, PhenotypeVar = 2.0, Heritability = 0.5, MutationScalar = 1e-3,
               MutationRatePerGamete = 2e-4, PopSize = 100, NoGenerations = 600, Sample = 10, InitFitnessSD = 5.0)
    p2 = Params(p; NoTraits = 2, MutCorr = 0.3, SelCorr = 0.2)
    cyc1(θs) = foldl(⊗, [withparams(L_.generation, p; θ = t) for t in θs]) ⋅ L_.migrate(length(θs))
    cyc2(φs) = foldl(⊗, [withparams(L_.generation, p2; θ = 2.0 .* [cosd(φ), sind(φ)]) for φ in φs]) ⋅ L_.migrate(length(φs))
    flat(x) = reduce(vcat, [y isa Number ? [y] : vec(y) for y in x])
    rel(a, b) = maximum(abs.(flat(a) .- flat(b))) / maximum(abs.(flat(b)))
    # reductions to the current closure (to 1e-12; they are in fact identical)
    for (f, tag) in ((l -> run_moments(p; erosion = true, loci = l).geno_var, "one deme"),
                     (l -> run_moments(p; between = true, loci = l).metapop.within, "one deme, between"),
                     (l -> run_moments(p2; erosion = true, loci = l).gmatrix, "one deme, two traits"),
                     (l -> run_moments(p; cycle = cyc1(zeros(4)), migration_rate = 0.01, between = true, loci = l).metapop.within,
                      "common optimum"),
                     (l -> run_moments(p; cycle = cyc1([2.0, 1, -1, -2]), migration_rate = 0.0, between = true, loci = l).metapop.within,
                      "no migration"),
                     (l -> run_moments(p2; cycle = cyc2([0.0, 90.0]), migration_rate = 0.0, between = true, loci = l).metapop.within,
                      "no migration, two traits"))
        a, b = f(:rare_local), f(:rare_alleles)
        @test rel(a, b) <= 1e-12
    end
    # local optima with migration: the gain raises the within-deme variance above the current
    # closure (Lean local_fixed_point_ge) and not above the multilocus Gaussian model
    w(l) = run_moments(p; cycle = cyc1([2.0, 1, -1, -2]), migration_rate = 0.01, between = true, loci = l).metapop.within[end]
    @test w(:rare_alleles) < w(:rare_local) <= w(:multilocus) * (1 + 1e-9)
    w2(l) = QTLSim.tr(run_moments(p2; cycle = cyc2([0.0, 60, 120, 180]), migration_rate = 0.01, between = true, loci = l).metapop.within[end])
    @test w2(:rare_alleles) < w2(:rare_local) <= w2(:multilocus) * (1 + 1e-9)
    # the default: one trait in several demes, and several traits, resolve to :rare_local
    @test w(:auto) == w(:rare_local)
    @test w2(:auto) == w2(:rare_local)
    # the gain: the eigenvalue expansion, the island form, the reductions and the limits
    rng = Xoshiro(11)
    for k in (2, 3, 5)
        M = rand(rng, k, k) .* 0.05; for i in 1:k; M[i, i] = 0; end
        P = QTLSim.Diagonal(ones(k)) - QTLSim.Diagonal(vec(sum(M; dims = 2))) + M
        σ = 0.01 .* randn(rng, k)
        lead(e) = maximum(real, eigvals(P * (QTLSim.Diagonal(ones(k)) + e * QTLSim.Diagonal(σ))))
        h = 1e-3
        @test (lead(h) - 2lead(0.0) + lead(-h)) / (2h^2) ≈ σ' * QTLSim.local_weights(P) * σ rtol = 1e-4
    end
    P4(m) = QTLSim.migration_mixing(island_migration(4, m), fill(100.0, 4), trues(4))
    ℓ = [0.01, 0.004, -0.004, -0.01]
    varℓ = sum(abs2, ℓ .- sum(ℓ) / 4) / 4
    for m in (0.001, 0.01, 0.1)
        P = P4(m); mp = 1 - (P[1, 1] - P[1, 2])
        @test all(QTLSim.local_gain(P, ℓ) .≈ (1 - mp) / mp * varℓ)
    end
    @test QTLSim.local_gain(P4(0.01), ℓ)[1] > QTLSim.local_gain(P4(0.1), ℓ)[1] > 0
    @test abs(QTLSim.local_gain(fill(0.25, 4, 4), ℓ)[1]) < 1e-15        # complete mixing, m′ = 1
    @test all(iszero, QTLSim.local_gain(P4(0.01), fill(0.3, 4)))
    @test all(iszero, QTLSim.local_gain(P4(0.01), fill([0.3, -0.1], 4)))
    @test all(iszero, QTLSim.local_gain(Matrix(1.0 * QTLSim.I, 4, 4), ℓ))
    @test QTLSim.local_gain(fill(1.0, 1, 1), [0.5]) == [0.0]
    @test QTLSim.migration_classes([0.9 0.1 0; 0.1 0.9 0; 0 0 1.0]) == [[1, 2], [3]]
    @test QTLSim.local_sinv(0.04, 0.0) === 0.04
    @test QTLSim.local_sinv(0.04, 0.01) ≈ 0.02
    @test QTLSim.local_sinv(0.04, 0.05) == 0.0
    S = [0.05 0.01; 0.01 0.03]
    @test QTLSim.local_sinv(S, zeros(2, 2)) === S
    @test minimum(eigvals(Symmetric(QTLSim.local_sinv(S, [0.04 0.0; 0.0 -0.01])))) >= -1e-15
    ω, _ = QTLSim.pair_classes(p)
    st = QTLSim.hap_init(0.3, 100.0, 2, length(ω))
    α2 = QTLSim.mutation_effect_variance(p)
    @test QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2) == QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2, gain = 0.0)
    @test QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2, gain = 1.0)[1] ==
          QTLSim.hap_select(st, 0.03, 100.0, ω, 2; α2 = 0.0)[1]           # the rare loss switched off
    @test_throws ArgumentError run_moments(p; loci = :rare_locals)
end
