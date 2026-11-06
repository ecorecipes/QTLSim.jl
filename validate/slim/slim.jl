# Where to find SLiM (https://messerlab.org/slim/, version 5.2) for the scripts in this
# folder: the program named by ENV["SLIM"]; else `slim/build/slim` in the folder that
# holds the package's folder; else `slim` on the path.

"The SLiM executable, or `nothing` if there is none."
function find_slim()
    pkg = normpath(joinpath(@__DIR__, "..", ".."))
    for path in (get(ENV, "SLIM", ""), joinpath(pkg, "..", "slim", "build", "slim"),
                 something(Sys.which("slim"), ""))
        isfile(path) && return path
    end
    nothing
end

slim_missing() = "SLiM not found: build SLiM 5.2 (`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && " *
                 "cmake --build build --target slim` in a clone of https://github.com/MesserLab/SLiM), " *
                 "then set ENV[\"SLIM\"] or put `slim` on the path"
