# Shared helpers of the symbolic checks (included by every check_*.jl): logging, the
# pass/fail record, and the symbolic tools. Each check file includes this once, runs its
# `main()`, writes its own log and returns the number of failures; check_all.jl runs them
# all in one process.

using Symbolics, LinearAlgebra, Random, Printf
using QTLSim

const out = IOBuffer()
say(args...) = (s = string(args...); println(s); println(out, s); flush(stdout))
const failures = String[]
const nchecks = Ref(0)
function check(name, ok::Bool)
    nchecks[] += 1
    ok || push!(failures, name)
    say(@sprintf("  %-72s %s", name, ok ? "ok" : "FAIL"))
    ok
end

# write the log and the summary; returns the number of failed checks
function finish(log::AbstractString, t0::Real)
    say(@sprintf("\n%d checks, %d failed; %.0f s", nchecks[], length(failures), time() - t0))
    say(isempty(failures) ? "PASS" : "FAIL: " * join(failures, "; "))
    write(log, String(take!(out)))
    length(failures)
end

# run `main` when the file is the script itself (not when check_all.jl includes it)
runmain(main, file) = abspath(PROGRAM_FILE) == abspath(file) && exit(main() == 0 ? 0 : 1)

# --- symbolic helpers ------------------------------------------------------------------

simp(e) = Symbolics.simplify(Symbolics.simplify_fractions(Symbolics.expand(e)); expand = true)
iszero_sym(e) = (z = simp(e); isequal(z, 0) || (Symbolics.value(z) isa Number && iszero(Symbolics.value(z))))
iszero_sym(e::AbstractArray) = all(iszero_sym, e)

# the coefficient of x^j in e (by derivatives at 0): exact for a polynomial in x, the
# Taylor coefficient for a rational function regular at 0
function coeffj(e, x, j)
    d = e
    for _ in 1:j
        d = Symbolics.expand_derivatives(Differential(x)(d))
    end
    simp(Symbolics.substitute(d, Dict(x => 0)) / factorial(j))
end

# e, a polynomial in r of degree ≤ deg, reduced modulo r² = D: is it zero?
function iszero_mod_sqrt(e, r, D; deg = 6)
    e = Symbolics.expand(e)
    ev, od = 0, 0
    for j in 0:deg
        c = coeffj(e, r, j)
        if iseven(j)
            ev += c * D^(j ÷ 2)
        else
            od += c * D^(j ÷ 2)
        end
    end
    iszero_sym(Symbolics.expand(e - sum(coeffj(e, r, j) * r^j for j in 0:deg))) || return false   # degree check
    iszero_sym(ev) && iszero_sym(od)
end

num(e, d) = Float64(Symbolics.value(Symbolics.substitute(e, d)))

inv2(A) = [A[2, 2] -A[1, 2]; -A[2, 1] A[1, 1]] ./ (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
det2(A) = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]

runif(rng::AbstractRNG, a, b) = a + (b - a) * rand(rng)
randspd(rng::AbstractRNG, n; scale = 1.0, ridge = 0.1) = (B = randn(rng, n, n); Symmetric(scale .* (B * B' ./ n + ridge * I)) |> Matrix)
