# Tests: the two-stage emulator (src/two_stage.jl: QTLSim.two_stage_screen,
# emulate(p; loci = :two_stage), emulated_half_survival_rate(p; loci = :two_stage),
# emulate_sweep(…; loci = :two_stage)). Included by runtests.jl; run alone with
# `julia --project=vignettes test/runtests.jl two_stage`.

@testset "two-stage emulator" begin
    p = Params(Ploidy = 2, PopSize = 200, NoQTL = 50, NoGenerations = 1000, MoveFitnessSDGen = 500, MoveType = :continuous,
               MoveFitnessSD = 0.02, PhenotypeVar = 1.0, Heritability = 0.3, MutationScalar = 0.005,
               MutationRatePerGamete = 0.002, RecombPoissonMean = 3.0, InitFitnessSD = 3.0)
    kw = (; burnin = 100, horizon = 200)
    bk = (; blocks = 10, replicates = 8)
    scr = QTLSim.TWO_STAGE_SCREEN
    @test scr == (; ld = :dynamic, solver = :beta) && QTLSim.TWO_STAGE_BAND > 0
    h = emulated_half_survival_rate(p; loci = :spectrum, spectrum = scr, kw...)
    @test 0 < h.k50 < Inf
    w = 0.2

    # the screen: the band decided by one more evaluation, toward the flip
    for (f, inside) in ((0.5, false), (0.9, true), (1.1, true), (2.0, false))
        q = Params(p; MoveFitnessSD = f * h.k50)
        s = QTLSim.two_stage_screen(q; band = w, kw...)
        e = emulate(q; loci = :spectrum, spectrum = scr, kw...)
        @test isequal(s.result, e)
        @test s.inband == inside && s.evals == 2
        @test s.probe ≈ q.MoveFitnessSD * exp(e.P2000 > 0.5 ? w : -w)
        @test (s.probe_P2000 > 0.5) == ((e.P2000 > 0.5) ⊻ inside)
    end
    q = Params(p; MoveFitnessSD = 0.9h.k50)
    @test !QTLSim.two_stage_screen(q; band = 0, kw...).inband && QTLSim.two_stage_screen(q; band = 0, kw...).evals == 1
    @test QTLSim.two_stage_screen(Params(p; MoveFitnessSD = 0.01h.k50); band = Inf, kw...).inband
    @test_throws ArgumentError QTLSim.two_stage_screen(q; band = -1, kw...)
    # the grid screen also works (no beta: fallback −1)
    @test QTLSim.two_stage_screen(q; band = w, screen = (; ld = :dynamic), kw...).fallback == -1

    # emulate: outside the band the screen's result, inside the block model's (same seed)
    far = Params(p; MoveFitnessSD = 0.5h.k50)
    e1 = emulate(far; loci = :two_stage, band = w, seed = 7, bk..., kw...)
    s1 = emulate(far; loci = :spectrum, spectrum = scr, kw...)
    @test e1.stage == 1 && !e1.inband && e1.alive == -1 && e1.replicates == 0
    @test e1.P2000 == s1.P2000 && isequal(e1.T, s1.T) && e1.VG_move == s1.VG_move && e1.screen_P2000 == s1.P2000
    near = Params(p; MoveFitnessSD = 0.95h.k50)
    e2 = emulate(near; loci = :two_stage, band = w, seed = 7, bk..., kw...)
    b2 = QTLSim.blocks_emulate(near; blocks = 10, replicates = 8, seed = 7, within = :gaussian, kw...)
    @test e2.stage == 2 && e2.inband && e2.replicates == 8
    @test e2.P2000 == b2.P2000 && e2.alive == b2.alive && isequal(e2.T, b2.T)
    @test e2.screen_P2000 == emulate(near; loci = :spectrum, spectrum = scr, kw...).P2000
    @test isequal(e2, emulate(near; loci = :two_stage, band = w, seed = 7, bk..., kw...))
    @test 0 <= e2.P240k <= e2.P2000 <= 1
    # band = 0: always the screen; band = Inf: always the block model
    @test emulate(near; loci = :two_stage, band = 0, seed = 7, bk..., kw...).stage == 1
    @test emulate(far; loci = :two_stage, band = Inf, seed = 7, bk..., kw...).stage == 2
    @test_throws ArgumentError emulate(near; loci = :two_stage, band = w, bk..., kw...)       # needs a seed or an rng
    @test emulate(near; loci = :two_stage, band = w, rng = Xoshiro(3), bk..., kw...).stage == 2
    # the block model's own defaults are unchanged (40 blocks, snapping) when loci = :blocks
    @test isequal(emulate(near; loci = :blocks, replicates = 3, seed = 2, kw...),
                  QTLSim.blocks_emulate(near; blocks = 40, replicates = 3, seed = 2, within = :snap, kw...))

    # emulate_sweep: per-input seeds from job_seed, the stage and the screen's values recorded
    ps = [far, near]
    rows = emulate_sweep(ps; ids = ["a", "b"], loci = :two_stage, band = w, seed = 11, bk..., kw...)
    @test [r.stage for r in rows] == [1, 2]
    @test rows[2].P2000 == emulate(near; loci = :two_stage, band = w, seed = QTLSim.job_seed(11, "b"), bk..., kw...).P2000
    @test rows[1].screen_P2000 == s1.P2000 && rows[1].screen_evals == 2 && isnan(rows[1].T) == (s1.T === missing)
    mktempdir() do dir
        f = joinpath(dir, "ts.csv")
        emulate_sweep(ps[[2]]; ids = ["b"], cache = f, loci = :two_stage, band = w, seed = 11, bk..., kw...)
        r2 = emulate_sweep(ps; ids = ["a", "b"], cache = f, loci = :two_stage, band = w, seed = 11, bk..., kw...)
        @test [x.P2000 for x in r2] == [x.P2000 for x in rows] && [x.stage for x in r2] == [1, 2]
    end
    @test_throws ArgumentError emulate_sweep(ps; loci = :two_stage, band = w, bk..., kw...)

    # the half-survival rate: the screen's k50, then the block model's bisected from ±band around it
    t = emulated_half_survival_rate(p; loci = :two_stage, band = w, seed = 5, bk..., kw...)
    @test t.stage == 2 && isapprox(t.k50_screen, h.k50; rtol = 2e-3) && 0 < t.k50 < Inf && t.evals_blocks >= 2
    # with common random numbers, the block model's P2000 is ≥ 1/2 just below k50 and < 1/2 just above
    b = QTLSim.blocks_half_survival_rate(p; blocks = 10, replicates = 8, seed = 5, within = :gaussian, kw...)
    @test 0.5 < t.k50 / b.k50 < 2
    # a bracket that misses the root is moved toward it (or falls back to the wide search)
    m = QTLSim.blocks_half_survival_rate(p; blocks = 10, replicates = 8, seed = 5, within = :gaussian, kw...,
                                         bracket = (0.2b.k50, 0.3b.k50))
    @test 0.5 < m.k50 / b.k50 < 2 && haskey(m, :bracketed)
    @test !haskey(b, :bracketed)
    @test_throws ArgumentError emulated_half_survival_rate(p; loci = :two_stage, band = w, bk..., kw...)
    @test_throws ArgumentError emulated_half_survival_rate(p; loci = :two_stage, band = 0, seed = 1, bk..., kw...)
end
