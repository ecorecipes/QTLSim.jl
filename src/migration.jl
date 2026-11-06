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
