using Documenter
using QTLSim
import Catlab              # loads the categorical layer, so its docstrings are available

const QTLSimCatlabExt = Base.get_extension(QTLSim, :QTLSimCatlabExt)

DocMeta.setdocmeta!(QTLSim, :DocTestSetup, :(using QTLSim); recursive=true)

# Tutorials are rendered quarto vignettes copied into docs/src/tutorials by
# scripts/sync_vignettes.jl. The page list is built from the files on disk so
# it never has to be maintained by hand.
function tutorial_pages()
    dir = joinpath(@__DIR__, "src", "tutorials")
    isdir(dir) || return Pair{String,String}[]
    files = sort(filter(f -> endswith(f, ".md"), readdir(dir)))
    map(files) do f
        title = f
        for line in eachline(joinpath(dir, f))
            m = match(r"^#\s+(.*)", line)
            if m !== nothing
                title = String(strip(m.captures[1]))
                break
            end
        end
        return title => "tutorials/" * f
    end
end

pages = Any["Home" => "index.md", "Parameters" => "parameters.md", "Output" => "output.md",
             "Running QTSim input files" => "reproducing.md", "Several traits" => "traits_design.md"]
tutorials = tutorial_pages()
isempty(tutorials) || push!(pages, "Tutorials" => tutorials)
push!(pages, "API reference" => "api.md")

# Source links need a commit to point at. A checkout without commits (a fresh
# scaffold, a tarball) builds without them instead of failing.
has_commit = success(pipeline(`git -C $(@__DIR__) rev-parse HEAD`; stderr=devnull))
remote_kw = has_commit ? (; repo=Remotes.GitHub("ecorecipes", "QTLSim.jl")) : (; remotes=nothing)

makedocs(;
         remote_kw...,
         modules=[QTLSim, QTLSimCatlabExt],
         sitename="QTLSim.jl",
         authors="Simon Frost",
         checkdocs=:exports,
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://ecorecipes.github.io/QTLSim.jl",
                                repolink="https://github.com/ecorecipes/QTLSim.jl",
                                edit_link="main",
                                size_threshold_ignore=["api.md"]),
         pages=pages)

"--no-deploy" in ARGS || deploydocs(; repo="github.com/ecorecipes/QTLSim.jl.git", devbranch="main")
