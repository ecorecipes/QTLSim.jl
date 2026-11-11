# The tables of a recorded pedigree. The functions that record into them, read them and
# write them are in ancestry.jl, after `World`; the types are here because a `World`
# carries one.

"""
    Pedigree

The recorded ancestry of a population or metapopulation ([`record_pedigree!`](@ref)), as
tables in the style of a tree sequence (Kelleher et al. 2018). Times are generations
since recording began: the founders are born at 0, their offspring at 1.

- **Nodes**, one per chromosome of every recorded individual: `node_time[n]`, the
  generation in which it was born, and `node_ind[n]`, its individual. Node ids increase
  with the time of birth, so a parent's id is below its child's.
- **Individuals**: `ind_time[i]`, `ind_birth[i]` the deme it was born in, `ind_deme[i]`
  the deme it was an adult in (after migration), and `ind_p1[i]`, `ind_p2[i]` its
  parents (0 for a founder, for a haploid's second parent, and for a parent that
  simplification removed). An individual's first chromosome came from `ind_p1`.
- **Edges** `(left[e], right[e], parent[e], child[e])`: node `child` inherited loci
  `left+1:right` (the half-open interval `[left, right)` of 0-based positions, as in
  tskit) from node `parent`. A gamete without a crossover is one edge, one with a
  crossover after locus `x` two, `[0, x)` and `[x, L)`.
- **Mutations** (when `mutations`): `mut_node[m]`, `mut_locus[m]` (1-based),
  `mut_time[m]` and the effect added, `mut_effect[m]` (with `n` traits, the `n` effects at
  `(m-1)n+1:mn`), in the order they happened. A founder's non-zero alleles are mutations at
  time 0, so that an allele is the sum, oldest first, of the mutations on its lineage.
- **Migrations**: individual `mig_ind[e]` moved from deme `mig_from[e]` to `mig_to[e]`
  after the generation `mig_time[e]` in which it was born, in the order they happened.

`worlds` are the populations recorded, in deme order, and `simplify_every` the interval
between simplifications ([`simplify!`](@ref); 0 for never), the last at `last_simplified`.
"""
mutable struct Pedigree{W}
    L::Int
    ploidy::Int
    ntraits::Int
    simplify_every::Int
    mutations::Bool
    last_simplified::Int
    nsimplified::Int
    node_time::Vector{Int32}
    node_ind::Vector{Int32}
    ind_time::Vector{Int32}
    ind_birth::Vector{Int32}
    ind_deme::Vector{Int32}
    ind_p1::Vector{Int32}
    ind_p2::Vector{Int32}
    left::Vector{Int32}
    right::Vector{Int32}
    parent::Vector{Int32}
    child::Vector{Int32}
    mut_node::Vector{Int32}
    mut_locus::Vector{Int32}
    mut_time::Vector{Int32}
    mut_effect::Vector{Float64}
    mig_ind::Vector{Int32}
    mig_time::Vector{Int32}
    mig_from::Vector{Int32}
    mig_to::Vector{Int32}
    worlds::Vector{W}
end

Pedigree{W}(L::Integer, ploidy::Integer, ntraits::Integer; simplify_every::Integer = 0,
            mutations::Bool = true) where {W} =
    Pedigree{W}(L, ploidy, ntraits, simplify_every, mutations, 0, 0,
                (Int32[] for _ in 1:11)..., Int32[], Int32[], Int32[], Float64[],
                Int32[], Int32[], Int32[], Int32[], W[])

function Base.show(io::IO, ped::Pedigree)
    print(io, "Pedigree(", length(ped.node_time), " nodes, ", length(ped.left), " edges, ",
          length(ped.ind_time), " individuals, ", length(ped.mut_node), " mutations, ",
          length(ped.mig_ind), " migrations; ", length(ped.worlds), " deme",
          length(ped.worlds) == 1 ? "" : "s", ")")
end

"""
    Ancestry

What a `World` carries while its pedigree is recorded (`w.anc`, otherwise `nothing`): the
shared [`Pedigree`](@ref) `ped`, the deme's index `deme`, its clock `gen` (generations
since recording began; the adults were born at `gen`), the node of each chromosome of
the adults and of the offspring in the buffer, `nodes[c, i]` and `snodes[c, j]`, and
their individuals, `inds[i]` and `sinds[j]`.
"""
mutable struct Ancestry{W}
    ped::Pedigree{W}
    deme::Int32
    gen::Int
    nodes::Matrix{Int32}
    snodes::Matrix{Int32}
    inds::Vector{Int32}
    sinds::Vector{Int32}
end

# --- recording (called by the kernel's operators when a World has an Ancestry) -------

# offspring j of the deme, of parents p1 and p2 (indices among the adults; p2 = 0 for a
# haploid): a new individual and its nodes
function newborn!(a::Ancestry, j::Int, p1::Int, p2::Int)
    ped = a.ped
    t = Int32(a.gen + 1)
    push!(ped.ind_time, t); push!(ped.ind_birth, a.deme); push!(ped.ind_deme, a.deme)
    @inbounds push!(ped.ind_p1, a.inds[p1])
    @inbounds push!(ped.ind_p2, p2 == 0 ? Int32(0) : a.inds[p2])
    id = Int32(length(ped.ind_time))
    @inbounds a.sinds[j] = id
    @inbounds for c in 1:ped.ploidy
        push!(ped.node_time, t); push!(ped.node_ind, id)
        a.snodes[c, j] = Int32(length(ped.node_time))
    end
    nothing
end

@inline function pushedge!(ped::Pedigree, l::Integer, r::Integer, parent::Int32, child::Int32)
    push!(ped.left, Int32(l)); push!(ped.right, Int32(r))
    push!(ped.parent, parent); push!(ped.child, child)
    nothing
end

# chromosome c of offspring j is a copy of chromosome `strand` of parent `par`
function inherit!(a::Ancestry, c::Int, j::Int, par::Int, strand::Int)
    @inbounds pushedge!(a.ped, 0, a.ped.L, a.nodes[strand, par], a.snodes[c, j])
end

# chromosome c of offspring j starts on chromosome `strand` of parent `par` and switches
# after each of the 0-based loci pts[1:npts] (increasing)
function inherit!(a::Ancestry, c::Int, j::Int, par::Int, strand::Int, pts::Vector{Int}, npts::Int)
    ped = a.ped
    @inbounds begin
        child = a.snodes[c, j]
        l = 0
        for ip in 1:npts
            r = pts[ip] + 1
            pushedge!(ped, l, r, a.nodes[strand, par], child)
            l = r
            strand = 3 - strand
        end
        pushedge!(ped, l, ped.L, a.nodes[strand, par], child)
    end
    nothing
end

# a mutation at locus `pos` of chromosome c of offspring j; its effects are pushed by
# the caller, one per trait
@inline function mutation!(a::Ancestry, c::Int, j::Int, pos::Int)
    ped = a.ped
    @inbounds push!(ped.mut_node, a.snodes[c, j])
    push!(ped.mut_locus, Int32(pos)); push!(ped.mut_time, Int32(a.gen + 1))
    nothing
end

# --- reading ---------------------------------------------------------------------------

"""
    marginal_tree(ped, locus)

The tree of the recorded chromosomes at `locus` (1-based), as a vector `par` over the
nodes: `par[n]` is the node from which node `n` inherited that locus, or 0 if `n` is a
root (a founder, or after simplification the most recent common ancestor of its samples).
"""
function marginal_tree(ped::Pedigree, locus::Integer)
    1 <= locus <= ped.L || throw(ArgumentError("locus $locus is not in 1:$(ped.L)"))
    par = zeros(Int32, length(ped.node_time))
    x = Int32(locus - 1)
    @inbounds for e in eachindex(ped.left)
        ped.left[e] <= x < ped.right[e] && (par[ped.child[e]] = ped.parent[e])
    end
    par
end

"""
    mrca(ped, tree, a, b)

The most recent common ancestor of nodes `a` and `b` in the marginal tree `tree`
([`marginal_tree`](@ref)), or 0 if they have none in the record.
"""
function mrca(ped::Pedigree, tree::AbstractVector{<:Integer}, a::Integer, b::Integer)
    time = ped.node_time
    @inbounds while a != b
        (a == 0 || b == 0) && return 0
        # move the younger one up; ids increase with the time of birth
        if time[a] > time[b] || (time[a] == time[b] && a > b)
            a = tree[a]
        else
            b = tree[b]
        end
    end
    Int(a)
end

"""
    allele(ped, node, locus; tree = marginal_tree(ped, locus), trait = 1)

The allelic effect of `node` at `locus`, rebuilt from the pedigree alone: the sum, oldest
first, of the recorded mutations at the locus on the node's lineage. It is the kernel's
value bit for bit, since the kernel adds the same numbers in the same order. Needs a
pedigree recorded with `mutations = true`.
"""
function allele(ped::Pedigree, node::Integer, locus::Integer; tree = marginal_tree(ped, locus),
                trait::Integer = 1)
    ped.mutations || throw(ArgumentError("this pedigree was recorded without mutations"))
    path = Int32[]
    n = Int32(node)
    while n != 0
        push!(path, n); n = tree[n]
    end
    x = 0.0
    nt = ped.ntraits
    for u in Iterators.reverse(path), m in eachindex(ped.mut_node)
        ped.mut_node[m] == u && ped.mut_locus[m] == locus && (x += ped.mut_effect[(m - 1) * nt + trait])
    end
    x
end

"""
    alleles(ped, locus; trait = 1)

The allelic effect of every node at `locus`, rebuilt from the pedigree ([`allele`](@ref)),
in one pass down the tree.
"""
function alleles(ped::Pedigree, locus::Integer; trait::Integer = 1)
    ped.mutations || throw(ArgumentError("this pedigree was recorded without mutations"))
    tree = marginal_tree(ped, locus)
    nt = ped.ntraits
    # mutations of each node at the locus, in order: chained through `next`
    first_m = zeros(Int, length(tree)); last_m = zeros(Int, length(tree))
    next = zeros(Int, length(ped.mut_node))
    for m in eachindex(ped.mut_node)
        ped.mut_locus[m] == locus || continue
        u = ped.mut_node[m]
        last_m[u] == 0 ? (first_m[u] = m) : (next[last_m[u]] = m)
        last_m[u] = m
    end
    x = zeros(length(tree))
    for u in eachindex(tree)             # parents before children: ids increase with time
        v = tree[u] == 0 ? 0.0 : x[tree[u]]
        m = first_m[u]
        while m != 0
            v += ped.mut_effect[(m - 1) * nt + trait]
            m = next[m]
        end
        x[u] = v
    end
    x
end
