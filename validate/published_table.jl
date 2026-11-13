# The published table, `Repository/Results/Master Output file.xlsx`, and the matching of its
# rows to the shipped input files by parameters. Shared by validate/published_sweep.jl and
# validate/emulator.jl; needs XLSX (the vignettes environment).
#
# In the table "Generations to extinction" counts from the move start (generation of
# death − 9999; 240000 = survived) and "Extinct (0,1)" is 1 for SURVIVORS.

using XLSX

function load_master(path)
    dt = XLSX.readtable(path, "Master output file")
    lab = string.(dt.column_labels)
    c(name) = dt.data[findfirst(==(name), lab)]
    num(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)     # some cells are stored as text
    I(name) = round.(Int, num.(c(name))); F(name) = num.(c(name))
    (L = I("No. Loci"), u = I("recombination u=3,5"), H = F("Heritability"),
     PV = F("PhenotypeVar"), MS = F("MutationScalar Vm/Ve"), d = F("d_env"),
     T = I("Generations to extinction"), surv = I("Extinct (0,1)"),
     pm = F("Phenotype Mean at 10000"), gv = F("Genotype Var at 10000"), ns = I("P Survived at 10000"))
end

# the conditions a table row must meet to be the row of an input file q
row_matches(tab, i, L, u, q::Params) =
    tab.L[i] == L && tab.u[i] == u &&
    abs(tab.H[i] - q.Heritability) < 1e-9 && abs(tab.PV[i] - q.PhenotypeVar) < 1e-9 &&
    # some rows round MutationScalar or d_env to 3 significant figures
    abs(tab.d[i] - q.MoveFitnessSD) <= 5e-3 * q.MoveFitnessSD &&
    abs(tab.MS[i] - q.MutationScalar) <= 5e-3 * q.MutationScalar

"The unique table row of input `q` in stratum (L, u), or 0. Scans the whole table."
function match_row(tab, L, u, q::Params)
    hits = findall(i -> row_matches(tab, i, L, u, q), eachindex(tab.L))
    length(hits) == 1 ? hits[1] : 0
end

"""
    row_index(tab) -> f(L, u, q)

The same matching as `match_row`, by binary search on Heritability within each stratum
(fast enough for all 119,984 inputs).
"""
function row_index(tab)
    by = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in eachindex(tab.L)
        push!(get!(by, (tab.L[i], tab.u[i]), Int[]), i)
    end
    for v in values(by)
        sort!(v; by = i -> tab.H[i])
    end
    function (L, u, q::Params)
        v = get(by, (L, u), Int[])
        a = searchsortedfirst(v, q.Heritability - 1e-9; by = x -> x isa Int ? tab.H[x] : x)
        hits = Int[]
        for k in a:length(v)
            i = v[k]
            tab.H[i] > q.Heritability + 1e-9 && break
            row_matches(tab, i, L, u, q) && push!(hits, i)
        end
        length(hits) == 1 ? hits[1] : 0
    end
end
