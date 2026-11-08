# Reading and writing QTSim's options files (QTSim `PassOptFile`, the keywords of
# `CMD_LIST` in reference/QTSim/src/PassOptionFile.h).

"""
    read_qtsim_options(path_or_io) -> Params

Read a QTSim options file into [`Params`](@ref), as QTSim's `PassOptFile` reads it:
- one `Keyword value…` per line, keywords case-insensitive;
- a line whose first word starts with `#` is a comment;
- `MoveFitnessSD <generation> <size>` and `ChangeFitnessSD <generation> <width>` take two
  values;
- `StandardisePhenotype` and `Drift` take none and switch their option each time they
  appear.

The seed (`Seed`) goes into `Params`' `Seed` field. Note that QTLSim's random numbers are
not QTSim's, so the same seed reproduces a QTLSim run, not the QTSim run.

Where the files differ from `Params`' field names:

| QTSim keyword | `Params` field |
|:--|:--|
| `Generations` | `NoGenerations` |
| `MoveFitnessSD g d` | `MoveFitnessSDGen = g`, `MoveFitnessSD = d` |
| `ChangeFitnessSD g w` | `ChangeFitnessGen = g`, `ChangeFitnessSD = w` |
| `OptimumMoveType` | `MoveType` (`:phenotype`, `:fitness`, `:continuous`) |
| `RecombinationPoissonMean` | `RecombPoissonMean` |

The other keywords (`Sample`, `PopSize`, `NoQTL`, `Ploidy`, `PhenotypeVar`,
`Heritability`, `MutationScalar`, `MutationRatePerGamete`, `InitFitnessSD`,
`StandardisePhenotype`) are field names.

Two rules differ from QTSim's:
- An unknown keyword is an `ArgumentError`, as it is an error in QTSim. So is an unknown
  move type, which QTSim silently ignores, keeping the phenotype move.
- `Drift` is read but dropped with a warning: QTSim switches drift off at the start of
  every run, so the option never has an effect there.

See also [`write_qtsim_options`](@ref).
"""
read_qtsim_options(path::AbstractString) = open(read_qtsim_options, path)

function read_qtsim_options(io::IO)
    kw = Dict{Symbol,Any}()
    drift = false
    for (n, line) in enumerate(eachline(io))
        f = split(line)
        (isempty(f) || startswith(f[1], '#')) && continue
        k = lowercase(f[1])
        nargs(m) = length(f) == m + 1 ||
            throw(ArgumentError("line $n: $(f[1]) takes $m value$(m == 1 ? "" : "s"), got $(length(f) - 1)"))
        int(s) = parse(Int, s); float(s) = parse(Float64, s)
        if k == "help"
        elseif k == "sample";       nargs(1); kw[:Sample] = int(f[2])
        elseif k == "popsize";      nargs(1); kw[:PopSize] = int(f[2])
        elseif k == "noqtl";        nargs(1); kw[:NoQTL] = int(f[2])
        elseif k == "ploidy";       nargs(1); kw[:Ploidy] = int(f[2])
        elseif k == "generations";  nargs(1); kw[:NoGenerations] = int(f[2])
        elseif k == "phenotypevar"; nargs(1); kw[:PhenotypeVar] = float(f[2])
        elseif k == "heritability"; nargs(1); kw[:Heritability] = float(f[2])
        elseif k == "mutationscalar";        nargs(1); kw[:MutationScalar] = float(f[2])
        elseif k == "mutationratepergamete"; nargs(1); kw[:MutationRatePerGamete] = float(f[2])
        elseif k == "recombinationpoissonmean"; nargs(1); kw[:RecombPoissonMean] = float(f[2])
        elseif k == "initfitnesssd"; nargs(1); kw[:InitFitnessSD] = float(f[2])
        elseif k == "seed";          nargs(1); kw[:Seed] = int(f[2])
        elseif k == "movefitnesssd"
            nargs(2); kw[:MoveFitnessSDGen] = int(f[2]); kw[:MoveFitnessSD] = float(f[3])
        elseif k == "changefitnesssd"
            nargs(2); kw[:ChangeFitnessGen] = int(f[2]); kw[:ChangeFitnessSD] = float(f[3])
        elseif k == "optimummovetype"
            nargs(1); t = Symbol(lowercase(f[2]))
            t in (:phenotype, :fitness, :continuous) ||
                throw(ArgumentError("line $n: unknown OptimumMoveType $(f[2])"))
            kw[:MoveType] = t
        elseif k == "standardisephenotype"
            nargs(0); kw[:StandardisePhenotype] = !get(kw, :StandardisePhenotype, false)
        elseif k == "drift"
            nargs(0); drift = !drift
        else
            throw(ArgumentError("line $n: unknown QTSim option $(f[1])"))
        end
    end
    drift && @warn "QTSim ignores Drift (RunSim switches it off), so it is dropped"
    Params(; kw...)
end

"""
    write_qtsim_options(path_or_io, p::Params; seed = p.Seed)

Write `p` as a QTSim options file, the inverse of [`read_qtsim_options`](@ref): every
option QTSim reads, with QTSim's keywords. The seed is written when there is one.
QTSim cannot express two `Params` settings, so the file warns for them:
- `CrossoverModel = :all`: QTSim applies only the lowest crossover point;
- `Drift = true`: QTSim switches drift off.

It is an `ArgumentError` to write `MoveType = :red_noise` or `:periodic`: QTSim has no
fluctuating optimum, and would run a different model.
"""
write_qtsim_options(path::AbstractString, p::Params; kw...) =
    open(io -> write_qtsim_options(io, p; kw...), path, "w")

function write_qtsim_options(io::IO, p::Params; seed = p.Seed)
    p.NoTraits == 1 || throw(ArgumentError("QTSim has one trait; NoTraits = $(p.NoTraits) cannot be written"))
    p.MoveType in (:phenotype, :fitness, :continuous) ||
        throw(ArgumentError("QTSim has no MoveType :$(p.MoveType); only :phenotype, :fitness and :continuous can be written"))
    p.CrossoverModel === :qtsim || @warn "QTSim has only its own crossover model; CrossoverModel = :$(p.CrossoverModel) is not written"
    p.Drift && @warn "QTSim switches drift off; Drift = true is not written"
    p.Demography === :fixed || @warn "QTSim has a fixed population size; Demography = :$(p.Demography) is not written"
    println(io, "Sample ", p.Sample)
    println(io, "PopSize ", p.PopSize)
    println(io, "NoQTL ", p.NoQTL)
    println(io, "Generations ", p.NoGenerations)
    println(io, "Ploidy ", p.Ploidy)
    p.MoveFitnessSDGen != -1 && println(io, "MoveFitnessSD ", p.MoveFitnessSDGen, " ", p.MoveFitnessSD)
    p.ChangeFitnessGen != -1 && println(io, "ChangeFitnessSD ", p.ChangeFitnessGen, " ", p.ChangeFitnessSD)
    println(io, "OptimumMoveType ", p.MoveType)
    println(io, "PhenotypeVar ", p.PhenotypeVar)
    println(io, "Heritability ", p.Heritability)
    println(io, "MutationScalar ", p.MutationScalar)
    println(io, "MutationRatePerGamete ", p.MutationRatePerGamete)
    println(io, "RecombinationPoissonMean ", p.RecombPoissonMean)
    p.InitFitnessSD != -1 && println(io, "InitFitnessSD ", p.InitFitnessSD)
    p.StandardisePhenotype && println(io, "StandardisePhenotype")
    seed === nothing || println(io, "Seed ", seed)
    nothing
end
