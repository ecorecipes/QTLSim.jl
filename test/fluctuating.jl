# Tests: fluctuating. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl fluctuating`.

@testset "fluctuating optima: red noise and periodic moves" begin
    p = Params(PopSize=1000, NoQTL=1, Ploidy=2, NoGenerations=4000, PhenotypeVar=0.02, Heritability=0.5,
               MutationScalar=0.02, MutationRatePerGamete=0.5, InitFitnessSD=0.3)
    # parameters, defaults and validation
    @test (p.MoveNoiseSD, p.MoveAutocorr, p.MoveAmplitude, p.MovePeriod) == (0.0, 0.0, 0.0, 100.0)
    @test checkparams(Params(p; MoveType = :red_noise, MoveNoiseSD = 0.1, MoveAutocorr = 0.9)) isa Params
    @test checkparams(Params(p; MoveType = :periodic, MoveAmplitude = -0.1, MovePeriod = 3.5)) isa Params
    for bad in ((MoveAutocorr = 1.0,), (MoveAutocorr = -0.1,), (MoveNoiseSD = -1.0,), (MovePeriod = 0.0,),
                (MovePeriod = -5.0,), (MoveAmplitude = Inf,), (MoveType = :sinusoidal,))
        @test_throws ArgumentError World(Params(p; bad...))
    end
    # QTSim has neither move: writing one is an error, not a silent change of model
    for t in (:red_noise, :periodic)
        @test_throws ArgumentError write_qtsim_options(IOBuffer(), Params(p; MoveType = t))
    end
    # the new fields leave QTSim's moves alone: same seed, same run
    pc = Params(p; NoGenerations = 600, MoveFitnessSDGen = 300, MoveFitnessSD = 0.001, MoveType = :continuous)
    @test simulate(pc; seed = 4).pheno_mean ==
          simulate(Params(pc; MoveNoiseSD = 0.3, MoveAutocorr = 0.5, MoveAmplitude = 1.0, MovePeriod = 7.0); seed = 4).pheno_mean

    # θ trajectories in the kernel
    t0 = 300
    pp = Params(pc; MoveType = :periodic, MoveFitnessSD = 0.002, MoveAmplitude = 0.05, MovePeriod = 40.0)
    r = simulate(pp; seed = 1)
    t = r.gen .- t0
    @test all(r.theta[t .< 0] .== 0)
    @test r.theta[t .>= 0] ≈ 0.002 .* t[t .>= 0] .+ 0.05 .* sin.(2π .* t[t .>= 0] ./ 40)
    @test r.theta[t .>= 0] == fluctuating_optimum.(Ref(pp), 0.0, t[t .>= 0])
    @test simulate(pp; seed = 1) == r                                # no randomness in the move
    pr = Params(pc; MoveType = :red_noise, MoveFitnessSD = 0.001, MoveNoiseSD = 0.05, MoveAutocorr = 0.7)
    r1 = simulate(pr; seed = 2); r2 = simulate(pr; seed = 2)
    @test r1 == r2 && r1.theta != simulate(pr; seed = 3).theta     # deterministic by seed
    @test all(r1.theta[r1.gen .< t0] .== 0)
    rz = simulate(Params(pr; MoveNoiseSD = 0.0); seed = 2)          # σ = 0: the trend alone
    @test rz.theta[rz.gen .>= t0] ≈ 0.001 .* (rz.gen[rz.gen .>= t0] .- t0)
    # the noise is a stationary AR(1) process with SD σ and autocorrelation ρ
    for ρ in (0.0, 0.8)
        q = Params(p; MoveType = :red_noise, MoveFitnessSDGen = 0, MoveFitnessSD = 0.0, MoveNoiseSD = 2.0, MoveAutocorr = ρ)
        w = World(q); rng = Xoshiro(11)
        ε = map(0:199_999) do g
            w.gen = g; moveoptimum!(w, rng); w.θ
        end
        @test std(ε) ≈ 2.0 rtol = 0.03
        @test cor(ε[1:end-1], ε[2:end]) ≈ ρ atol = 0.01
        @test abs(mean(ε)) < 0.1
    end
    @test_throws ArgumentError moveoptimum!(World(Params(pr; MoveFitnessSDGen = 0)), nothing)

    # theory: the linear filter's response to a cycle and to red noise
    Vs = selection_variance(p)
    @test periodic_response(1e7; VG = 0.5, Vs = 1.5).gain ≈ 1 rtol = 1e-6
    @test periodic_response(1e7; VG = 0.5, Vs = 1.5).delay ≈ 1 / 0.25 rtol = 1e-4   # 1/β, the trend's delay
    @test red_noise_lag_var(1.0, 0.0; VG = 0.5, Vs = 1.5) ≈ 2 / (2 - 0.25)
    @test red_noise_lag_var(1.0, 0.999999; VG = 0.5, Vs = 1.5) < 1e-5
    @test expected_lag2(pc; VG = 0.004) ≈ steady_lag(pc; VG = 0.004)^2
    @test lag_load(pc; VG = 0.004) ≈ steady_lag(pc; VG = 0.004)^2 / (2 * (0.3^2 + 0.004 + envVar(p)))
    @test_throws ArgumentError expected_lag2(p; VG = 0.004)

    # the moments lowering: the periodic move lowers exactly, and the mean responds as predicted
    V = multilocus_equilibrium(p).VG          # the default's fixed point, which the move leaves alone
    P, A = 200.0, 0.1
    pm = Params(p; MoveType = :periodic, MoveFitnessSDGen = 500, MoveFitnessSD = 0.0, MoveAmplitude = A, MovePeriod = P)
    o = run_moments(pm; VG = V)
    @test o.theta[o.gen .>= 500] ≈ A .* sin.(2π .* (o.gen[o.gen .>= 500] .- 500) ./ P)
    pr_ = periodic_response(P; VG = V, Vs)
    late = o.gen .>= 2000
    @test o.pheno_mean[late] ≈ pr_.gain * A .* sin.(2π .* (o.gen[late] .- 500) ./ P .- pr_.phase) atol = 1e-8
    @test mean((o.theta[late] .- o.pheno_mean[late]) .^ 2) ≈ expected_lag2(pm; VG = V) rtol = 1e-3
    # red noise lowers to its expected optimum, the trend ...
    pn = Params(p; MoveType = :red_noise, MoveFitnessSDGen = 500, MoveFitnessSD = 1e-4, MoveNoiseSD = 0.05, MoveAutocorr = 0.8)
    on = run_moments(pn; VG = V)
    @test on.theta[on.gen .>= 500] ≈ 1e-4 .* (on.gen[on.gen .>= 500] .- 500)
    # ... or follows a realised optimum, whose mean squared lag is the theory's
    n = 100_000
    θs = zeros(n); rng = Xoshiro(5); ρ, σ = 0.8, 0.05
    θs[1] = σ * randn(rng)
    for i in 2:n; θs[i] = ρ * θs[i-1] + σ * sqrt(1 - ρ^2) * randn(rng); end
    orn = run_moments(Params(pn; NoGenerations = n, MoveFitnessSD = 0.0); VG = V, optimum = θs)
    @test orn.theta == θs
    @test mean((orn.theta[1000:end] .- orn.pheno_mean[1000:end]) .^ 2) ≈
          red_noise_lag_var(σ, ρ; VG = V, Vs) rtol = 0.05
    @test_throws ArgumentError run_moments(pn; optimum = zeros(10))

    # the kernel follows the cycle as the theory says, up to drift
    rk = simulate(Params(pm; NoGenerations = 5000); seed = 8)
    k = rk.gen .>= 1500
    VGk = mean(rk.geno_var[k]); ν = 2π / P; tk = rk.gen[k] .- 500
    X = hcat(sin.(ν .* tk), cos.(ν .* tk), ones(length(tk)))
    a, b, _ = X \ rk.pheno_mean[k]
    th = periodic_response(P; VG = VGk, Vs)
    @test hypot(a, b) / A ≈ th.gain rtol = 0.15
    @test atan(-b, a) ≈ th.phase atol = 0.15

    # the Eidos lowering: the new moves are in the generated model, and only there
    g = lifecycle().generation
    for (t, marker) in ((:red_noise, "moveAutocorr * p1.getValue(\"EPS\")"), (:periodic, "sin(2 * PI * t / movePeriod)"))
        s = slim_script(g, Params(pm; MoveType = t))
        @test occursin(marker, s) && occursin("p1.setValue(\"THETA0\", THETA)", s)
        @test occursin("defineConstant(\"moveType\", \"$t\")", s)
        for c in ("moveNoiseSD", "moveAutocorr", "moveAmplitude", "movePeriod")
            @test occursin("defineConstant(\"$c\"", s)
        end
    end
    sc = slim_script(g, pc)
    @test !occursin("THETA0", sc) && !occursin("moveNoiseSD", sc)
    @test sc == slim_script(g, Params(pc; MoveNoiseSD = 0.3, MovePeriod = 7.0))
    @test !occursin("THETA0", lower(g, :eidos).code)
end
