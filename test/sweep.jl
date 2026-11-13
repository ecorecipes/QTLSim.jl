# Tests: sweep. Included by runtests.jl, which loads the packages and defines `Ext` and
# the operators; run alone with `julia --project=vignettes test/runtests.jl sweep`.

@testset "sweep and replicate: threads, cache, resumption, errors" begin
    p = Params(PopSize = 100, NoQTL = 1, NoGenerations = 50)
    summ(o) = (VG = o.geno_var[end], zbar = o.pheno_mean[end], ext = o.extinct)
    # results are in job order and do not depend on the number of tasks or the job order
    r1 = replicate(p, 6; seed = 10, summary = summ, ntasks = 1)
    r3 = replicate(p, 6; seed = 10, summary = summ, ntasks = 3, cost = s -> -s)
    @test r1 == r3 == [summ(simulate(p; seed = s)) for s in 10:15]
    @test sweep(x -> x^2, 1:10; ntasks = 4, cost = x -> x % 3) == (1:10) .^ 2
    @test isempty(sweep(identity, Int[]))
    mktempdir() do dir
        csv = joinpath(dir, "runs.csv")
        # a cache round-trip returns the rows, Float64s exactly; strings with commas are quoted
        f(j) = (; j, x = j / 7, flag = isodd(j), name = "a,\"b\" $j")
        rows = sweep(f, 1:5; cache = csv, ntasks = 2)
        @test rows == f.(1:5)
        @test readline(csv) == "id,j,x,flag,name"
        # resumption: done jobs are not rerun, and the rows equal an uninterrupted sweep's
        ran = Int[]; lk = ReentrantLock()
        g(j) = (lock(lk) do; push!(ran, j); end; f(j))
        @test sweep(g, 1:8; cache = csv, ntasks = 3) == f.(1:8)
        @test sort(ran) == 6:8
        # a line cut short by an interruption is dropped and its job rerun
        last_id = parse(Int, split(readlines(csv)[end], ',')[1])
        write(csv, read(csv, String)[1:end-3]); empty!(ran)
        @test sweep(g, 1:8; cache = csv) == f.(1:8) && ran == [last_id]
        @test countlines(csv) == 9
        # the cache is keyed by id; replicate keys by seed
        c2 = joinpath(dir, "rep.csv")
        @test replicate(p, 3; seed = 4, summary = summ, cache = c2) == [summ(simulate(p; seed = s)) for s in 4:6]
        @test first.(split.(readlines(c2)[2:end], ',')) ⊆ ["4", "5", "6"]
        # errors: duplicate ids, rows that do not fit the cache, results that are not rows
        @test_throws ArgumentError sweep(f, [1, 1]; cache = joinpath(dir, "dup.csv"))
        @test_logs (:error,) match_mode = :any @test_throws ArgumentError sweep(j -> (; y = j), 9:9; cache = csv)
        @test_logs (:error,) match_mode = :any @test_throws ArgumentError sweep(j -> j, 1:2; cache = joinpath(dir, "n.csv"))
        # a failing job stops the sweep, rows finished before it stay cached, and its error is rethrown
        c3 = joinpath(dir, "fail.csv")
        h(j) = j == 3 ? error("job 3 fails") : (; j)
        @test_logs (:error,) match_mode = :any @test_throws ErrorException sweep(h, 1:5; cache = c3, ntasks = 1)
        @test readlines(c3) == ["id,j", "1,1", "2,2"]
        @test sweep(j -> (; j), 1:5; cache = c3) == [(; j) for j in 1:5]
    end
    @test_throws ArgumentError sweep(identity, 1:2; ntasks = 0)
end
