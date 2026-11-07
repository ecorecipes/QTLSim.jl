# Migration graphs among demes: a k×k matrix `M` whose entry `M[a, b]` is the probability
# that an individual of the deme in position a starts an exchange with a random individual
# of the deme in position b, each generation. Exchanges keep every deme's size; each moves
# two individuals, so about `M[a, b] + M[b, a]` of deme a's individuals come from deme b
# per generation (with equal deme sizes). The kernel (`migrate!`) and Eidos (`slim_script`)
# lowerings of `migrate(k)` take such a matrix as `migration_rate`; a number r is the
# island model `island_migration(k, r)`, drawn as before.

"""
    island_migration(k, m)

The island model on `k` demes: an individual starts an exchange with probability `m`, with
a deme chosen uniformly among the other `k − 1`. As a matrix, `m/(k − 1)` off the
diagonal. Equivalent in distribution to `migration_rate = m`.
"""
function island_migration(k::Integer, m::Real)
    k >= 2 || throw(ArgumentError("migration needs at least two demes; got $k"))
    [a == b ? 0.0 : m / (k - 1) for a in 1:k, b in 1:k]
end

"""
    stepping_stone(k, m; circular = true)

The stepping-stone model on `k` demes in a line (or a ring, when `circular`): an
individual starts an exchange with probability `m`, with one of its neighbours, each with
probability `m/2`. At the ends of a line, the one neighbour gets `m/2`.
"""
function stepping_stone(k::Integer, m::Real; circular::Bool = true)
    k >= 2 || throw(ArgumentError("migration needs at least two demes; got $k"))
    M = zeros(k, k)
    for a in 1:k, b in (a - 1, a + 1)
        if 1 <= b <= k
            M[a, b] += m / 2
        elseif circular && k > 2
            M[a, mod1(b, k)] += m / 2
        end
    end
    M
end

"""
    check_migration(M, k)

Check that `M` is a migration matrix on `k` demes: `k × k`, non-negative, zero on the
diagonal, every row summing to at most one. Returns `M`.
"""
function check_migration(M::AbstractMatrix, k::Integer)
    size(M) == (k, k) || throw(ArgumentError("a migration matrix on $k demes is $k × $k; got $(size(M))"))
    all(>=(0), M) || throw(ArgumentError("migration rates must be non-negative"))
    all(iszero, (M[a, a] for a in 1:k)) || throw(ArgumentError("a migration matrix has zeros on its diagonal"))
    all(<=(1), sum(M; dims = 2)) || throw(ArgumentError("each row of a migration matrix sums to at most one"))
    M
end

"""
    migration_mixing(M, N, alive)

The expected composition of the demes after the kernel's migration (`migrate!`), with
`M[a, b]` the probability that a visit to an individual of deme `a` starts an exchange
with a random individual of deme `b`, `N` the deme sizes and `alive` which demes take
part: `P[a, c]` is the expected fraction of deme `a` that came from deme `c`.

`migrate!` visits the demes in order, and each individual of a deme once. When deme `a`
is visited, its individuals all have the composition of the deme, `P[a, :]` (earlier
demes chose partners in it at random), and its unvisited ones do not change while it is
visited. Each visit exchanges with deme `b` with probability `r = M[a, b]`, sending an
individual of composition `P[a, :]` and receiving a random one of `b`, whose expected
composition moves toward `P[a, :]` by a factor `1 − r/N_b` per visit. Over the `N_a`
visits, with `f = (1 − r/N_b)^{N_a}`, deme `b` becomes `P[a, :] + f (P[b, :] − P[a, :])`
and deme `a` receives `N_a r P[a, :] + N_b (1 − f)(P[b, :] − P[a, :])`, in place of the
`N_a r` individuals it sends. These expectations are exact for any deme sizes, given the
order in which the demes are visited.

`migrate!` draws that order uniformly at random at every call (so that migration is
exchangeable, `Migration.lean`), and the expected composition is the average of the
above over all `k!` orders: exactly, by enumerating them, for `k ≤ 8`, and over 20,000
random orders from a fixed seed (so deterministically, to about 1%) for more demes.
"""
function migration_mixing(M::AbstractMatrix, N::AbstractVector, alive::AbstractVector{Bool})
    k = length(N)
    k <= 1 && return migration_mixing(M, N, alive, collect(1:k))
    orders = k <= 8 ? _permutations(k) : (rng = Xoshiro(20261107); [randperm(rng, k) for _ in 1:20_000])
    sum(migration_mixing(M, N, alive, o) for o in orders) ./ length(orders)
end

# all permutations of 1:k, by Heap's algorithm
function _permutations(k::Integer)
    a, c, out = collect(1:k), ones(Int, k), [collect(1:k)]
    i = 1
    while i <= k
        if c[i] < i
            j = isodd(i) ? 1 : c[i]
            a[j], a[i] = a[i], a[j]
            push!(out, copy(a))
            c[i] += 1
            i = 1
        else
            c[i] = 1
            i += 1
        end
    end
    out
end

# the expected composition for one visiting order of the demes
function migration_mixing(M::AbstractMatrix, N::AbstractVector, alive::AbstractVector{Bool}, order::AbstractVector{<:Integer})
    k = length(N)
    P = [Float64(a == c) for a in 1:k, c in 1:k]
    for a in order
        alive[a] || continue
        cs = P[a, :]
        received = zeros(k)
        sent = 0.0
        for b in 1:k
            r = M[a, b]
            (b == a || r == 0 || !alive[b]) && continue
            f = (1 - r / N[b])^N[a]
            received .+= N[a] * r .* cs .+ N[b] * (1 - f) .* (P[b, :] .- cs)
            P[b, :] .= cs .+ f .* (P[b, :] .- cs)
            sent += N[a] * r
        end
        P[a, :] .= ((N[a] - sent) .* cs .+ received) ./ N[a]
    end
    P
end

