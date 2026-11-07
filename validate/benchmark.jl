# Time per generation of the kernel, of the old dense kernel and of QTSim on the paper's
# inputs.
#
# For each number of loci L, the input is the first published input of the stratum
# QT<L>-Possion3.0 (PopSize 1000, diploid, RecombinationPoissonMean 3, the published
# mutation rate and scalar), or, without the paper's data, a synthetic input of the same
# design; L = 1 uses the 100-locus input with NoQTL = 1. Only generations before the move
# (generation 10000) are timed, the stabilising phase in which the published runs spend
# most of their time.
#
# Kernel: one thread, `generation!` with records (as `simulate`), after a warm-up run
# that compiles everything. The block genomes take a while to reach their working set,
# so the time per generation is measured over generations BURN+1 … BURN+T and, for
# comparison, over the first BURN generations; the median over REPS seeds.
# Dense: the kernel as it was before block genomes (test/dense_reference.jl), timed the
# same way over generations BURN+1 … BURN+T.
# QTSim: the same input run for G1 and G2 generations (its table to /dev/null, every
# generation printed as in the published inputs); the time per generation is the
# difference of the median wall-clock times over REPS seeds, divided by G2 − G1, so
# start-up and the table's header are not counted.
#
# The 1-minute load average is recorded before and after each L, since the machine may
# be shared: the numbers are only worth keeping when it stayed well below the number of
# performance cores. With ENV["BENCH_MAXLOAD"] set, each L waits until the load is below it.
#
# Run from the package root, on one thread:
#   julia -t 1 --project=. validate/benchmark.jl [L1,L2,…] [REPS=3]
# Writes validate/benchmark.log.

using QTLSim, Random, Printf, Statistics, Dates
using QTLSim: records
include(joinpath(@__DIR__, "qtsim.jl"))
include(joinpath(@__DIR__, "..", "test", "dense_reference.jl"))

const LS   = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], ',')) : [1, 100, 500, 1000, 5000]
const REPS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
const BURN, T = 4000, 1000                 # kernels: warm-up generations, timed generations
const G1, G2 = 500, 2500                   # QTSim: two run lengths
const SIMDIR = normpath(joinpath(@__DIR__, "..", "..", "Repository", "Simulations"))
const QTSIM = find_qtsim()
const LOG = joinpath(@__DIR__, "benchmark.log")
const MAXLOAD = parse(Float64, get(ENV, "BENCH_MAXLOAD", "Inf"))   # wait for a lower load before each L

loadavg() = parse(Float64, split(strip(read(`sysctl -n vm.loadavg`, String), ['{', '}', ' ']))[1])

function input(L)
    f = joinpath(SIMDIR, @sprintf("QT%04d-Possion3.0", max(L, 100)), "Input-00001.txt")
    p = isfile(f) ? read_qtsim_options(f) :
        Params(PopSize = 1000, Ploidy = 2, NoQTL = L, PhenotypeVar = 0.011, Heritability = 0.57,
               MutationScalar = 0.0013, MutationRatePerGamete = 1e-4, RecombPoissonMean = 3.0,
               InitFitnessSD = 7.6, MoveFitnessSDGen = 10_000, MoveFitnessSD = 0.15,
               MoveType = :continuous, Sample = 1)
    Params(p; NoQTL = L)
end

"Mean ms per generation over the first `BURN` generations and over the `T` after them."
function kernel_ms(p, seed; M = QTLSim)
    w = M.World(Params(p; NoGenerations = BURN + T))
    rng = MersenneTwister(seed)
    M.init_env!(w, rng)
    rec = records()
    t0 = time_ns()
    for it in 0:BURN-1
        w.gen = it; M.generation!(w, rng; rec)
    end
    t1 = time_ns()
    for it in BURN:BURN+T-1
        w.gen = it; M.generation!(w, rng; rec)
    end
    t2 = time_ns()
    w.extinct && error("extinct during the benchmark (L = $(p.NoQTL), seed $seed)")
    ((t1 - t0) / 1e6 / BURN, (t2 - t1) / 1e6 / T)
end

function qtsim_s(p, G, seed)
    path = tempname() * ".txt"
    write_qtsim_options(path, Params(p; NoGenerations = G); seed)
    t = @elapsed run(pipeline(`$QTSIM $path`, stdout = devnull))
    rm(path)
    t
end

function main()
    io = open(LOG, "w")
    out(s) = (print(s); print(io, s); flush(io))
    cpu = Sys.cpu_info()
    qt = QTSIM === nothing ? "none" : relpath(QTSIM, normpath(joinpath(@__DIR__, "..")))
    out(@sprintf("# %s; Julia %s, %d thread(s); %s, %d cores, %.0f GB; QTSim: %s\n",
                 Dates.format(now(), "yyyy-mm-dd HH:MM"), VERSION, Threads.nthreads(),
                 cpu[1].model, length(cpu), Sys.total_memory() / 2^30, qt))
    out(@sprintf("# ms per generation; kernels: median of %d seeds, over generations 1–%d (first) and %d–%d; QTSim: (t(%d) − t(%d))/%d, median of %d\n",
                 REPS, BURN, BURN + 1, BURN + T, G2, G1, G2 - G1, REPS))
    small = Params(input(100); NoGenerations = 50)
    kernel_ms(small, 1); kernel_ms(small, 1; M = DenseReference)      # compile
    out(@sprintf("%6s  %8s  %8s  %8s  %8s  %8s  %8s  %6s  %6s\n", "NoQTL", "first", "kernel", "dense",
                 "QTSim", "dense/k", "QTSim/k", "load0", "load1"))
    for L in LS
        p = input(L)
        while loadavg() >= MAXLOAD; sleep(30); end
        l0 = loadavg()
        ks = [kernel_ms(p, s) for s in 1:REPS]
        k1, k2 = median(first.(ks)), median(last.(ks))
        d = median(last(kernel_ms(p, s; M = DenseReference)) for s in 1:REPS)
        q = NaN
        if QTSIM !== nothing
            a = median(qtsim_s(p, G1, s) for s in 1:REPS)
            b = median(qtsim_s(p, G2, s) for s in 1:REPS)
            q = 1000 * (b - a) / (G2 - G1)
        end
        l1 = loadavg()
        out(@sprintf("%6d  %8.3f  %8.3f  %8.3f  %8.3f  %8.2f  %8.2f  %6.2f  %6.2f\n",
                     L, k1, k2, d, q, d / k2, q / k2, l0, l1))
    end
    out(@sprintf("# finished %s; BENCH_MAXLOAD = %s\n", Dates.format(now(), "yyyy-mm-dd HH:MM"), MAXLOAD))
    close(io)
end

main()
