# Where to find QTSim, the C program this package ports, for the validation scripts and
# vignettes: the program named by ENV["QTSIM"]; else the build of the source bundled in
# reference/QTSim (`make -C reference/QTSim`); else the binary in the paper's code if it
# is unpacked next to the package, as Repository/Code/QTSim/src/QTSim.

"The QTSim executable to compare with, or `nothing` if there is none."
function find_qtsim()
    pkg = normpath(joinpath(@__DIR__, ".."))
    for path in (get(ENV, "QTSIM", ""),
                 joinpath(pkg, "reference", "QTSim", "build", "QTSim"),
                 joinpath(pkg, "..", "Repository", "Code", "QTSim", "src", "QTSim"))
        isfile(path) && return path
    end
    nothing
end

qtsim_missing() = "QTSim not found: build it from the package root with `make -C reference/QTSim` " *
                  "(it needs a C compiler and GSL; see reference/QTSim/README.md), or set ENV[\"QTSIM\"]"
