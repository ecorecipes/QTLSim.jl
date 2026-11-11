# Ancestry: identity by descent read off the block pool, recorded pedigrees, realised
# coalescence times, and export as tree-sequence tables.
#
# Two sources of ancestry, of different strength:
#
# 1. Block ids (free, always there). Chromosomes are vectors of ids of shared blocks
#    (blocks.jl). Two chromosomes of the current population that hold the same id at a
#    block position are copies of one ancestral block, see `shared_blocks` for exactly
#    what that implies and what it does not.
# 2. A recorded pedigree (`record_pedigree!`, off by default): every offspring
#    chromosome, its parent's two chromosomes and the crossovers that made it, the
#    mutations, and the migrations, as tree-sequence tables (pedigree.jl). It gives the
#    complete ancestry of every locus whatever happens to block ids.

# --- identity by descent from block ids --------------------------------------------------

"""
    tag_founders!(w::World)
    tag_founders!(ws::AbstractVector{World})

Give every chromosome of the current population(s) its own blocks, copies of the ones it
holds, so that from now on two chromosomes hold the same block only if they inherited it
from one of these chromosomes. A `World(p)` starts with all chromosomes on one block of
zeros per position, and so with every pair identical at every block; after
`tag_founders!` no two are, and block identity ([`shared_blocks`](@ref)) counts descent
from the tagged generation. Allele values, block sums, genotypes and the random stream are
unchanged, bit for bit. Returns the number of blocks copied.
"""
function tag_founders!(w::World)
    pool, chr = w.pool, w.chr
    ncopied = 0
    for i in 1:w.N, c in axes(chr, 2), k in axes(chr, 1)
        id = chr[k, c, i]
        pool.refs[id] > 1 || continue
        nid = copyblock!(pool, pool, id, blocklen(pool, k))
        pool.refs[id] -= Int32(1)
        chr[k, c, i] = nid
        ncopied += 1
    end
    ncopied
end
tag_founders!(ws::AbstractVector{World}) = sum(tag_founders!, ws; init = 0)

"""
    shared_blocks(w, i, c, v, j, d) -> BitVector

For each block position, whether chromosome `c` of individual `i` of `w` and chromosome
`d` of individual `j` of `v` (current adults of populations that share a block pool, as
the demes of a metapopulation made by `⊗` do; `v` may be `w`) hold the same block.

# What block identity means

A block id is given to a chromosome when the block is made: at the start (one block of
zeros per position, shared by everyone, unless [`tag_founders!`](@ref) or `setgenomes!`
gave each chromosome its own), by a crossover strictly inside the block position between
two different blocks (the recombinant is a new block), or by a mutation in a block that
another slot also holds (copy on write; a block held by one slot only is changed in
place and keeps its id, which nobody else can then hold). Every other transmission copies
the id. Therefore, for two chromosomes present **at the same time**:

- **The same id implies** that both blocks are copies, through whole-block transmissions
  only, of one block of one ancestral chromosome `X` (the one in which the block was
  made, or last changed in place): at every locus of the block the two chromosomes are
  identical by descent from `X`, their most recent common ancestor at the locus is `X`
  or a descendant of it, no mutation has occurred anywhere in the block on either lineage
  below that ancestor, and the alleles are identical in state. With the initial shared
  blocks of zeros, `X` is the founding block itself, not a chromosome: then the same id
  says only that neither lineage has mutated or been recut in the block since the start
  (use `tag_founders!` to make ids count descent from individual founders).
- **Different ids imply nothing.** The two may still be identical by descent at some or
  all loci of the block: one lineage may carry a mutation elsewhere in the block, or the
  block may have been rebuilt by a crossover inside it (the loci on each side remain
  copies of their sources), or the genomes came through `setgenomes!` or through
  migration between pools.

So block identity is a lower bound on identity by descent: recent, unbroken and
unmutated. It is exact (same id if and only if identical by descent without an
intervening mutation at the locus) when blocks are single loci (`World(p; blocksize = 1)`)
and founders are tagged, since no crossover falls inside a block of one locus.

Ids are recycled when a block is freed, so ids saved from an earlier generation mean
nothing later; compare only chromosomes that are alive together. Ids of different pools
are unrelated (an `ArgumentError`).
"""
function shared_blocks(w::World, i::Integer, c::Integer, v::World, j::Integer, d::Integer)
    w.pool === v.pool || throw(ArgumentError("block ids of different pools cannot be compared"))
    (1 <= i <= w.N && 1 <= j <= v.N) || throw(BoundsError())
    [w.chr[k, c, i] == v.chr[k, d, j] for k in 1:w.pool.nb]
end
shared_blocks(w::World, i::Integer, c::Integer, j::Integer, d::Integer) = shared_blocks(w, i, c, w, j, d)

"""
    ibd_fraction(w, i, c, v, j, d)

The fraction of loci in blocks that the two chromosomes share ([`shared_blocks`](@ref)):
a lower bound on the fraction of the genome identical by descent between them.
"""
function ibd_fraction(w::World, i::Integer, c::Integer, v::World, j::Integer, d::Integer)
    s = shared_blocks(w, i, c, v, j, d)
    pool = w.pool
    sum(blocklen(pool, k) for k in 1:pool.nb if s[k]; init = 0) / pool.L
end
ibd_fraction(w::World, i::Integer, c::Integer, j::Integer, d::Integer) = ibd_fraction(w, i, c, w, j, d)

"""
    ibd_sharing(w::World)
    ibd_sharing(ws::AbstractVector{World}) -> (; ind, deme, between, block, weights)

Identity by descent from block ids ([`shared_blocks`](@ref)) over all pairs of
chromosomes of the current adults of the populations `ws` (demes sharing one pool), by
kind of pair, as probabilities that a pair holds the same block:

- `ind[a]`: the two chromosomes of one individual of deme `a` (`NaN` for haploids);
- `deme[a]`: chromosomes of two distinct individuals of deme `a`;
- `between[a, b]`: a chromosome of deme `a` and one of deme `b` (`between[a, a] = deme[a]`);
- `block`: the same three per block position, `block.ind[k, a]`, `block.deme[k, a]`,
  `block.between[k, a, b]`; the genome-wide values are their means weighted by the
  number of loci in each block, `weights[k]`.

All pairs are counted, exactly, in time proportional to the number of chromosomes. Extinct
demes and demes without a pair of that kind give `NaN`.
"""
function ibd_sharing(ws::AbstractVector{World})
    isempty(ws) && throw(ArgumentError("no populations"))
    pool = first(ws).pool
    all(w -> w.pool === pool, ws) || throw(ArgumentError("block ids of different pools cannot be compared"))
    k, nb = length(ws), pool.nb
    P = size(first(ws).chr, 2)
    N = [w.extinct ? 0 : w.N for w in ws]
    bind = fill(NaN, nb, k); bdeme = fill(NaN, nb, k); bbetween = fill(NaN, nb, k, k)
    count = zeros(Int, length(pool.refs), k)        # holders of each block, per deme
    ids = Int32[]
    for pos in 1:nb
        for a in 1:k, i in 1:N[a], c in 1:P
            id = ws[a].chr[pos, c, i]
            all(iszero, view(count, id, :)) && push!(ids, id)
            count[id, a] += 1
        end
        for a in 1:k
            same = P == 2 ? sum(ws[a].chr[pos, 1, i] == ws[a].chr[pos, 2, i] for i in 1:N[a]; init = 0) : 0
            P == 2 && N[a] > 0 && (bind[pos, a] = same / N[a])
            for b in a:k
                if a == b
                    pairs = sum(count[id, a] * (count[id, a] - 1) ÷ 2 for id in ids; init = 0) - same
                    total = (P * N[a]) * (P * N[a] - 1) ÷ 2 - (P == 2 ? N[a] : 0)
                    total > 0 && (bdeme[pos, a] = bbetween[pos, a, a] = pairs / total)
                else
                    pairs = sum(count[id, a] * count[id, b] for id in ids; init = 0)
                    total = P^2 * N[a] * N[b]
                    total > 0 && (bbetween[pos, a, b] = bbetween[pos, b, a] = pairs / total)
                end
            end
        end
        for id in ids
            count[id, :] .= 0
        end
        empty!(ids)
    end
    weights = [blocklen(pool, pos) / pool.L for pos in 1:nb]
    wmean(x) = sum(weights .* x)
    (ind = [wmean(view(bind, :, a)) for a in 1:k], deme = [wmean(view(bdeme, :, a)) for a in 1:k],
     between = [wmean(view(bbetween, :, a, b)) for a in 1:k, b in 1:k],
     block = (ind = bind, deme = bdeme, between = bbetween), weights)
end
ibd_sharing(w::World) = ibd_sharing([w])

"""
    ibd_fst(within, between) -> (; fst, f_within, f_between, f_total, pairwise)
    ibd_fst(ws::AbstractVector{World})

F_ST from probabilities of identity: `within[a]` for two genes of distinct individuals of
deme `a`, `between[a, b]` for genes of demes `a ≠ b` (the diagonal is ignored). With
`f₀` the mean of `within`, `f₁` the mean of `between` over `a ≠ b`, and
`f̄ = (f₀ + (k − 1) f₁)/k` the identity of two genes of the whole metapopulation (demes
weighted equally, pairs drawn with replacement of demes),

    F_ST = (f₀ − f̄)/(1 − f̄),

Nei's (1973) definition, the identity analogue of Slatkin's (1991) `(T̄ − T̄₀)/T̄` for
coalescence times. `pairwise[a, b]` is the same for the two demes `a` and `b` alone.

From a metapopulation `ws`, the identities are those of block ids
([`ibd_sharing`](@ref)): identity by descent without an intervening mutation or recut in
the block, which is identity in state under an infinite-alleles model of blocks. The
identities of a pedigree ([`pair_coalescence`](@ref), [`coalescence_cdf`](@ref)) can be
passed as `within` and `between` in the same way.
"""
function ibd_fst(within::AbstractVector{<:Real}, between::AbstractMatrix{<:Real})
    k = length(within)
    size(between) == (k, k) || throw(ArgumentError("between must be $k × $k"))
    k >= 2 || throw(ArgumentError("F_ST needs at least two demes"))
    f0 = sum(within) / k
    f1 = sum(between[a, b] for a in 1:k, b in 1:k if a != b) / (k * (k - 1))
    fT = (f0 + (k - 1) * f1) / k
    pairwise = [a == b ? 0.0 : (f = (within[a] + within[b]) / 2; t = (f + between[a, b]) / 2; (f - t) / (1 - t))
                for a in 1:k, b in 1:k]
    (fst = (f0 - fT) / (1 - fT), f_within = f0, f_between = f1, f_total = fT, pairwise)
end
ibd_fst(ws::AbstractVector{World}) = (s = ibd_sharing(ws); ibd_fst(s.deme, s.between))

# --- recording a pedigree ----------------------------------------------------------------

"""
    record_pedigree!(w::World; simplify_every = 0, mutations = true, tag = true) -> Pedigree
    record_pedigree!(ws::AbstractVector{World}; ...)

Start recording the pedigree of a population, or of the demes `ws` of a metapopulation
(for a metapopulation `s` made by `⊗`, `record_pedigree!(s)`), from the current adults,
the founders, on. Returns the [`Pedigree`](@ref), which the populations fill as they run
(`run!`, `generation!`, `run_composed`, the operators): for every offspring its parents;
for each of its chromosomes the parental chromosomes it copies and where it switches
between them (so the ancestor of every locus, whatever happens to block ids); the
mutations, with their effects (`mutations = true`); and each migration between demes.

Recording draws no random numbers and changes no value: a run with a pedigree is the same
run, bit for bit, as without one. It is off by default, and a `World` without a pedigree
does no work for it (the kernel still allocates nothing per generation).

# Cost

Without simplification the tables grow by every individual of every generation: 8 bytes
per chromosome, 20 per individual, 16 per edge (one per gamete, two with a crossover; more
under `CrossoverModel = :all`), 12 + 8·`NoTraits` per mutation and 16 per migration, about
140 bytes per diploid individual per generation with one crossover and one mutation per
gamete; see `validate/coalescence.log` for measured time and memory. The founders'
non-zero alleles are recorded as mutations at time 0, 20 bytes each, so starting from a
population that has already evolved costs up to `20·NoQTL` bytes per founder chromosome
(nothing for a fresh `World(p)`, whose alleles are zero). With
`simplify_every = n` the tables are simplified ([`simplify!`](@ref)) every `n`
generations, which keeps only what is ancestral to the living population: memory is then
bounded by the genealogy of the living chromosomes (at most `2Nk − 1` coalescence nodes
per distinct marginal tree) plus `n` generations of new records, at the price of the
individuals in between (their identities, parents and migrations).

`tag = true` first gives each founder chromosome its own blocks
([`tag_founders!`](@ref)), so that block identity agrees with the pedigree's founders.
Do not write genomes (`setgenomes!`, `lower!`) while recording: the pedigree would not
know about them.
"""
function record_pedigree!(ws::AbstractVector{World}; simplify_every::Integer = 0, mutations::Bool = true,
                          tag::Bool = true)
    isempty(ws) && throw(ArgumentError("no populations to record"))
    simplify_every >= 0 || throw(ArgumentError("simplify_every must be non-negative, got $simplify_every"))
    shape = genome_shape(first(ws))
    all(w -> genome_shape(w) == shape, ws) ||
        throw(ArgumentError("the demes of one pedigree have the same loci, ploidy and traits"))
    any(w -> w.anc !== nothing, ws) && throw(ArgumentError("a population is already recording a pedigree"))
    length(unique(objectid, ws)) == length(ws) || throw(ArgumentError("a population is listed twice"))
    p = first(ws).p
    ped = Pedigree{World}(p.NoQTL, p.Ploidy, p.NoTraits; simplify_every, mutations)
    append!(ped.worlds, ws)
    tag && tag_founders!(ws)
    nt = p.NoTraits
    for (a, w) in enumerate(ws)
        P, cap = p.Ploidy, size(w.chr, 3)
        anc = Ancestry{World}(ped, Int32(a), 0, zeros(Int32, P, cap), zeros(Int32, P, cap),
                              zeros(Int32, cap), zeros(Int32, cap))
        G = mutations ? genomes(w) : nothing
        for i in 1:(w.extinct ? 0 : w.N)
            push!(ped.ind_time, Int32(0)); push!(ped.ind_birth, Int32(a)); push!(ped.ind_deme, Int32(a))
            push!(ped.ind_p1, Int32(0)); push!(ped.ind_p2, Int32(0))
            anc.inds[i] = Int32(length(ped.ind_time))
            for c in 1:P
                push!(ped.node_time, Int32(0)); push!(ped.node_ind, anc.inds[i])
                n = anc.nodes[c, i] = Int32(length(ped.node_time))
                mutations || continue
                for l in 1:p.NoQTL        # a founder's non-zero alleles, as mutations at time 0
                    eff = nt == 1 ? (G[l, c, i],) : ntuple(r -> G[r, l, c, i], nt)
                    all(iszero, eff) && continue
                    push!(ped.mut_node, n); push!(ped.mut_locus, Int32(l)); push!(ped.mut_time, Int32(0))
                    append!(ped.mut_effect, eff)
                end
            end
        end
        w.anc = anc
    end
    ped
end
record_pedigree!(w::World; kw...) = record_pedigree!([w]; kw...)

"""
    stop_recording!(ped) -> ped

Stop recording: the populations of `ped` no longer carry it. The tables stay as they are.
"""
function stop_recording!(ped::Pedigree)
    for w in ped.worlds
        w.anc = nothing
    end
    ped
end

# swapbuffers! with a pedigree: the offspring are the adults, one generation on
function swap_ancestry!(w::World)
    a = w.anc::Ancestry{World}
    a.nodes, a.snodes = a.snodes, a.nodes
    a.inds, a.sinds = a.sinds, a.inds
    a.gen += 1
    ped = a.ped
    if ped.simplify_every > 0 && a.gen - ped.last_simplified >= ped.simplify_every &&
       all(v -> v.extinct || (v.anc::Ancestry{World}).gen == a.gen, ped.worlds)
        simplify!(ped)
    end
    nothing
end

# exchange_genomes! with a pedigree: the two individuals exchange demes
function exchange_ancestry!(w::World, i::Integer, v::World, j::Integer)
    a, b = w.anc, v.anc
    (a !== nothing && b !== nothing && a.ped === b.ped) ||
        throw(ArgumentError("migration between a population that records a pedigree and one that does not, " *
                            "or between two pedigrees; record the demes together, record_pedigree!([w, v])"))
    ped = a.ped
    @inbounds for c in axes(a.nodes, 1)
        a.nodes[c, i], b.nodes[c, j] = b.nodes[c, j], a.nodes[c, i]
    end
    @inbounds a.inds[i], b.inds[j] = b.inds[j], a.inds[i]
    @inbounds for (ind, from, to, t) in ((a.inds[i], b.deme, a.deme, b.gen), (b.inds[j], a.deme, b.deme, a.gen))
        push!(ped.mig_ind, ind); push!(ped.mig_time, Int32(t)); push!(ped.mig_from, from); push!(ped.mig_to, to)
        ped.ind_deme[ind] = to
    end
    nothing
end

"""
    pedigree_time(ped)

The current generation of the pedigree's clock: generations since recording began, for
the populations that are not extinct.
"""
pedigree_time(ped::Pedigree) = maximum((w.anc::Ancestry).gen for w in ped.worlds if !w.extinct; init = 0)

"""
    samples(ped) -> (; node, deme, ind)

The living chromosomes of a pedigree, the adults of its populations that are not
extinct: their nodes, the deme (index in `ped.worlds`) each is in, and their individuals.
Chromosomes of one individual are adjacent.
"""
function samples(ped::Pedigree)
    node, deme, ind = Int32[], Int32[], Int32[]
    for (d, w) in enumerate(ped.worlds)
        w.extinct && continue
        a = w.anc
        a === nothing && throw(ArgumentError("recording was stopped; the pedigree no longer knows who is alive"))
        for i in 1:w.N, c in axes(a.nodes, 1)
            push!(node, a.nodes[c, i]); push!(deme, Int32(d)); push!(ind, a.inds[i])
        end
    end
    (; node, deme, ind)
end

"""
    pedigree_bytes(ped)

The memory the pedigree's tables hold, in bytes (their contents, not their capacity).
"""
pedigree_bytes(ped::Pedigree) =
    sum(sizeof(getfield(ped, f)) for f in fieldnames(typeof(ped)) if getfield(ped, f) isa Vector{<:Number})

# --- simplification ------------------------------------------------------------------------

"""
    simplify!(ped) -> ped

Reduce the pedigree to the genealogy of the living chromosomes ([`samples`](@ref)): the
algorithm of Kelleher et al. (2018, PLoS Comput. Biol. 14: e1006581). It keeps the samples
and every node in which two or more lineages ancestral to samples join at some locus, with
edges only over the intervals where they do; nodes through which a single lineage passes,
ancestry above the most recent common ancestor of all samples at a locus, and everything
not ancestral to a sample are dropped. The marginal tree of the samples at every locus,
with the times of its nodes, is unchanged, so coalescence times are exact; mutations are
moved to the node below them on the same lineage, so [`allele`](@ref) still rebuilds every
living allele bit for bit. Individuals without a remaining node are dropped with their
migrations, and a parent that was dropped becomes 0.

`record_pedigree!(…; simplify_every = n)` calls it every `n` generations, when every
living deme has finished the generation. Called by hand, call it between generations.
"""
function simplify!(ped::Pedigree)
    L = Int32(ped.L)
    nn = length(ped.node_time)
    smp = samples(ped)
    out_of = zeros(Int32, nn)             # the output node of an input node, in order of creation
    out_in = Int32[]                      # and the input node of each output node
    newout(u) = (push!(out_in, u); out_of[u] = Int32(length(out_in)))
    # the ancestry of each input node: segments (left, right, output node), node u's at
    # astart[u] .+ (1:acount[u])
    aL, aR, aN = Int32[], Int32[], Int32[]
    astart, acount = zeros(Int, nn), zeros(Int32, nn)
    for u in smp.node
        newout(u)
        astart[u] = length(aL); acount[u] = 1
        push!(aL, Int32(0)); push!(aR, L); push!(aN, out_of[u])
    end
    eL, eR, eP, eC = Int32[], Int32[], Int32[], Int32[]     # output edges
    order = sortperm(ped.parent; rev = true)                # parents from the youngest: ids increase with time
    qL, qR, qN, bps = Int32[], Int32[], Int32[], Int32[]
    e = 1
    while e <= length(order)
        u = ped.parent[order[e]]
        empty!(qL); empty!(qR); empty!(qN)
        while e <= length(order) && ped.parent[order[e]] == u
            ed = order[e]; c = ped.child[ed]
            for s in astart[c]+1:astart[c]+acount[c]
                l, r = max(aL[s], ped.left[ed]), min(aR[s], ped.right[ed])
                l < r && (push!(qL, l); push!(qR, r); push!(qN, aN[s]))
            end
            e += 1
        end
        (isempty(qL) || acount[u] > 0) && continue          # nothing ancestral below u (or u is a sample)
        empty!(bps); append!(bps, qL); append!(bps, qR); sort!(bps); unique!(bps)
        astart[u] = length(aL)
        for ib in 1:length(bps)-1
            l, r = bps[ib], bps[ib+1]
            n, node = 0, Int32(0)
            for s in eachindex(qL)
                qL[s] <= l && qR[s] >= r && (n += 1; node = qN[s])
            end
            n == 0 && continue
            if n > 1                                         # a coalescence in u over [l, r)
                v = out_of[u] == 0 ? newout(u) : out_of[u]
                for s in eachindex(qL)
                    qL[s] <= l && qR[s] >= r && (push!(eL, l); push!(eR, r); push!(eP, v); push!(eC, qN[s]))
                end
                node = v
            end
            if acount[u] > 0 && aR[end] == l && aN[end] == node
                aR[end] = r
            else
                push!(aL, l); push!(aR, r); push!(aN, node); acount[u] += 1
            end
        end
    end
    # renumber the output nodes by time of birth, the oldest first (stable, so that the
    # chromosomes of an individual stay together)
    no = length(out_in)
    rank = invperm(sortperm([ped.node_time[u] for u in out_in]; alg = MergeSort))
    newid(u) = out_of[u] == 0 ? Int32(0) : Int32(rank[out_of[u]])
    # mutations move to the output node that carries their locus on the mutant's lineage
    keepm = Int[]
    mnode = Int32[]
    for m in eachindex(ped.mut_node)
        u, x = ped.mut_node[m], ped.mut_locus[m] - Int32(1)
        for s in astart[u]+1:astart[u]+acount[u]
            if aL[s] <= x < aR[s]
                push!(keepm, m); push!(mnode, Int32(rank[aN[s]]))
                break
            end
        end
    end
    nt = ped.ntraits
    ped.mut_effect = [ped.mut_effect[(m - 1) * nt + r] for m in keepm for r in 1:nt]
    ped.mut_locus = ped.mut_locus[keepm]; ped.mut_time = ped.mut_time[keepm]; ped.mut_node = mnode
    # nodes and individuals
    node_time = Vector{Int32}(undef, no); node_ind = Vector{Int32}(undef, no)
    for (o, u) in enumerate(out_in)
        node_time[rank[o]] = ped.node_time[u]; node_ind[rank[o]] = ped.node_ind[u]
    end
    keepi = sort!(unique(node_ind))
    newind = zeros(Int32, length(ped.ind_time))
    newind[keepi] .= 1:length(keepi)
    remap(x) = x == 0 ? Int32(0) : newind[x]
    ped.node_time = node_time; ped.node_ind = newind[node_ind]
    ped.ind_p1 = remap.(ped.ind_p1[keepi]); ped.ind_p2 = remap.(ped.ind_p2[keepi])
    ped.ind_time = ped.ind_time[keepi]; ped.ind_birth = ped.ind_birth[keepi]; ped.ind_deme = ped.ind_deme[keepi]
    keepg = findall(i -> newind[i] != 0, ped.mig_ind)
    ped.mig_ind = newind[ped.mig_ind[keepg]]; ped.mig_time = ped.mig_time[keepg]
    ped.mig_from = ped.mig_from[keepg]; ped.mig_to = ped.mig_to[keepg]
    # edges: renumbered, sorted, and adjacent intervals of one parent and child joined
    eP .= rank[eP]; eC .= rank[eC]
    eo = sortperm(collect(zip(eP, eC, eL)))
    left, right, parent, child = Int32[], Int32[], Int32[], Int32[]
    for i in eo
        if !isempty(left) && parent[end] == eP[i] && child[end] == eC[i] && right[end] == eL[i]
            right[end] = eR[i]
        else
            push!(left, eL[i]); push!(right, eR[i]); push!(parent, eP[i]); push!(child, eC[i])
        end
    end
    ped.left, ped.right, ped.parent, ped.child = left, right, parent, child
    # the populations' own references
    for w in ped.worlds
        a = w.anc::Ancestry{World}
        fill!(a.snodes, 0); fill!(a.sinds, 0)
        n = w.extinct ? 0 : w.N
        for i in 1:size(a.nodes, 2)
            for c in axes(a.nodes, 1)
                a.nodes[c, i] = i <= n ? newid(a.nodes[c, i]) : Int32(0)
            end
            a.inds[i] = i <= n ? remap(a.inds[i]) : Int32(0)
        end
    end
    ped.last_simplified = pedigree_time(ped)
    ped.nsimplified += 1
    ped
end

# --- realised coalescence times --------------------------------------------------------------

"""
    tmrca(ped, locus, a, b; tree = marginal_tree(ped, locus))

The time, in generations before the present ([`pedigree_time`](@ref)), of the most recent
common ancestor at `locus` of nodes `a` and `b`, or `nothing` if they have none in the
record (their lineages reach different founders). Two copies of one parental gene copy
have time 1.
"""
function tmrca(ped::Pedigree, locus::Integer, a::Integer, b::Integer; tree = marginal_tree(ped, locus))
    m = mrca(ped, tree, a, b)
    m == 0 ? nothing : pedigree_time(ped) - Int(ped.node_time[m])
end

"""
    PairCoalescence

Realised pairwise coalescence times of the living gene copies of a pedigree
([`pair_coalescence`](@ref)), as counts of pairs of gene copies at one locus, summed over
`nloci` loci, `t` generations after recording began:

- `W[T, a]`: the two copies in one individual of deme `a` whose most recent common
  ancestor lived `T` generations ago, `T` in `1:t`;
- `S[T, a, b]`: copies in two distinct individuals, one of deme `a` and one of deme `b`
  (`a == b`: of one deme; `S[T, a, b] == S[T, b, a]`);
- row `t + 1` of each: pairs with no common ancestor in the record, whose lineages reach
  two different founder chromosomes (coalescence time greater than `t`).

Every pair is counted. See [`coalescence_cdf`](@ref) and [`coalescence_mean`](@ref).
"""
struct PairCoalescence
    t::Int
    nloci::Int
    W::Matrix{Int}
    S::Array{Int,3}
end

"""
    pair_coalescence(ped; loci = 1:ped.L) -> PairCoalescence

The realised coalescence times of all pairs of living gene copies at each of `loci`, read
off the recorded pedigree: within individuals, between individuals of a deme and between
demes (a [`PairCoalescence`](@ref)). These are the times whose expectations
[`deme_coalescence`](@ref) gives for the neutral life cycle: `E[min(T, t)]` is
[`coalescence_mean`](@ref) and `P(T ≤ s)` [`coalescence_cdf`](@ref). The demes are those
of `ped.worlds`, by where each individual is now.

All pairs are counted by one pass up each marginal tree (the number of leaves of each
deme under each node), so the cost per locus is proportional to the number of nodes and
edges, not to the number of pairs.
"""
function pair_coalescence(ped::Pedigree; loci = 1:ped.L)
    smp = samples(ped)
    k, t = length(ped.worlds), pedigree_time(ped)
    nn = length(ped.node_time)
    W = zeros(Int, t + 1, k); S = zeros(Int, t + 1, k, k)
    slot = zeros(Int32, nn)
    cnt = Int[]                               # leaves below each visited node, per deme: k per slot
    nd = zeros(Int, k)
    for d in smp.deme
        nd[d] += 1
    end
    for x in loci
        tree = marginal_tree(ped, x)
        fill!(slot, 0); empty!(cnt)
        for (n, d) in zip(smp.node, smp.deme)
            slot[n] = length(cnt) ÷ k + 1
            append!(cnt, (j == d ? 1 : 0 for j in 1:k))
        end
        @inbounds for n in nn:-1:1            # children before parents
            s = slot[n]
            (s == 0 || tree[n] == 0) && continue
            q = tree[n]
            if slot[q] == 0
                slot[q] = length(cnt) ÷ k + 1
                append!(cnt, (0 for _ in 1:k))
            end
            so, qo = (s - 1) * k, (slot[q] - 1) * k
            T = t - ped.node_time[q]
            for a in 1:k
                ca, pa = cnt[so+a], cnt[qo+a]
                for b in a:k
                    pairs = a == b ? pa * ca : pa * cnt[so+b] + cnt[qo+b] * ca
                    pairs == 0 && continue
                    S[T, a, b] += pairs
                    a == b || (S[T, b, a] += pairs)
                end
            end
            for a in 1:k
                cnt[qo+a] += cnt[so+a]
            end
        end
        if ped.ploidy == 2                    # the two copies of one individual are not a pair of S
            for i in 1:2:length(smp.node)
                d = smp.deme[i]
                m = mrca(ped, tree, smp.node[i], smp.node[i+1])
                m == 0 && (W[t+1, d] += 1; continue)
                T = t - ped.node_time[m]
                W[T, d] += 1; S[T, d, d] -= 1
            end
        end
    end
    nl = length(loci)
    for a in 1:k, b in 1:k                    # the pairs that did not coalesce
        total = a == b ? nd[a] * (nd[a] - 1) ÷ 2 - (ped.ploidy == 2 ? nd[a] ÷ 2 : 0) : nd[a] * nd[b]
        S[t+1, a, b] = nl * total - sum(view(S, 1:t, a, b))
    end
    PairCoalescence(t, nl, W, S)
end

_frac(num, den) = den == 0 ? NaN : num / den

"""
    coalescence_cdf(pc, s) -> (; W, S)

The fraction of pairs of gene copies whose most recent common ancestor lived at most `s`
generations ago: `W[a]` within individuals of deme `a`, `S[a, b]` between distinct
individuals of demes `a` and `b`. With `s = pc.t`, the probability of identity by
descent from the founders. `NaN` where there are no pairs.
"""
function coalescence_cdf(pc::PairCoalescence, s::Integer)
    0 <= s <= pc.t || throw(ArgumentError("s must be in 0:$(pc.t), got $s"))
    k = size(pc.W, 2)
    (W = [_frac(sum(view(pc.W, 1:s, a)), sum(view(pc.W, :, a))) for a in 1:k],
     S = [_frac(sum(view(pc.S, 1:s, a, b)), sum(view(pc.S, :, a, b))) for a in 1:k, b in 1:k])
end

"""
    coalescence_mean(pc) -> (; W, S)

The mean of `min(T, t)` over pairs, `T` the coalescence time and `t = pc.t` the length of
the record (pairs without a common ancestor count `t`): the realised value of what
`deme_coalescence(P, N; t)` predicts.
"""
function coalescence_mean(pc::PairCoalescence)
    k, t = size(pc.W, 2), pc.t
    tt = [1:t; t]
    (W = [_frac(sum(tt .* view(pc.W, :, a)), sum(view(pc.W, :, a))) for a in 1:k],
     S = [_frac(sum(tt .* view(pc.S, :, a, b)), sum(view(pc.S, :, a, b))) for a in 1:k, b in 1:k])
end

"""
    coalescence_transform(pc, r) -> (; W, S)

The mean over pairs of `exp(−rT)` for pairs with a common ancestor in the record, 0 for
the others. With `r = 2μ` it is the probability that two gene copies are identical by
descent with no mutation on either lineage since their common ancestor, for a mutation
rate `μ` per locus per generation: what block identity measures for blocks of one locus.
"""
function coalescence_transform(pc::PairCoalescence, r::Real)
    k, t = size(pc.W, 2), pc.t
    e = [exp.(-r .* (1:t)); 0.0]
    (W = [_frac(sum(e .* view(pc.W, :, a)), sum(view(pc.W, :, a))) for a in 1:k],
     S = [_frac(sum(e .* view(pc.S, :, a, b)), sum(view(pc.S, :, a, b))) for a in 1:k, b in 1:k])
end

"""
    pedigree_ibd(ped, loci = locus:locus; block = loci) -> Vector{Int32}

Identity by descent without an intervening mutation, from the pedigree: a label for each
living chromosome ([`samples`](@ref), in that order) such that two chromosomes have the
same label if and only if, at locus `first(loci)`, they have a common ancestor in the
record and no mutation has occurred at any locus of `block` on either lineage below that
ancestor. The label is the youngest node on the chromosome's lineage that carries a
mutation in `block`, or the lineage's root.

This is the pedigree's definition of what block ids measure ([`shared_blocks`](@ref)): for
`World`s of one locus per block and tagged founders, two chromosomes hold the same block
at locus `l` if and only if `pedigree_ibd(ped, l)` gives them the same label; for larger
blocks, holding the same block at position `k` implies the same label at every locus of
the block with `block` the block's loci. Needs `mutations = true`.

A `block` wider than the locus itself needs a pedigree that was never simplified:
[`simplify!`](@ref) keeps each mutation on the lineage of its own locus only, so whether a
mutation at another locus of the block lies on this locus's lineage is no longer
recorded (an `ArgumentError`).
"""
function pedigree_ibd(ped::Pedigree, loci::AbstractUnitRange{<:Integer}; block::AbstractUnitRange{<:Integer} = loci)
    ped.mutations || throw(ArgumentError("this pedigree was recorded without mutations"))
    (ped.nsimplified == 0 || block == first(loci):first(loci)) ||
        throw(ArgumentError("a block wider than the locus needs a pedigree that was never simplified"))
    tree = marginal_tree(ped, first(loci))
    mutant = falses(length(tree))
    for m in eachindex(ped.mut_node)
        ped.mut_locus[m] in block && (mutant[ped.mut_node[m]] = true)
    end
    label = zeros(Int32, length(tree))
    for n in eachindex(tree)              # parents before children
        label[n] = (mutant[n] || tree[n] == 0) ? n : label[tree[n]]
    end
    label[samples(ped).node]
end
pedigree_ibd(ped::Pedigree, locus::Integer; kw...) = pedigree_ibd(ped, locus:locus; kw...)

# --- tree-sequence tables ----------------------------------------------------------------------

"""
    tskit_tables(ped) -> (; sequence_length, nodes, edges, individuals, populations, sites, mutations, migrations)

The pedigree as tree-sequence tables in tskit's data model (ids from 0, −1 for none; time
in generations before the present, [`pedigree_time`](@ref)), each a `NamedTuple` of
columns:

- `nodes`: `is_sample` (1 for the living chromosomes), `time`, `population` (the deme the
  node's individual was an adult in, after migration) and `individual`;
- `edges`: `left`, `right`, `parent`, `child`, positions being loci, `[left, right)` of
  `0:L`, sorted as tskit requires (by the parent's time, the youngest first, then parent,
  child and `left`);
- `individuals`: `flags` (0), `parents` (two per individual), and the extra columns
  `time`, `birth_deme` and `deme`;
- `populations`: one per deme, `metadata` the JSON `{"name":"deme_a"}`;
- `sites` and `mutations` (if mutations were recorded): a site per mutated locus, at
  `position` locus − 1 with `ancestral_state` `"0"`; a mutation's `derived_state` is the
  allelic effect after it, as a decimal string (the effects on several traits joined by
  commas), `effect` (an extra column) what it added, `parent` the mutation it follows on
  the same lineage, `time` the age of the chromosome it arose in;
- `migrations`: for each recorded move of an individual and each of its nodes, a row
  `left = 0, right = L, node, source, dest, time` in tskit's backward-in-time sense
  (`source` is the deme the individual moved *to*), at `time` = the node's age + ½, between
  the node's birth and its parents' adulthood, in order of time.

[`write_tskit`](@ref) writes them as text; [`check_tables`](@ref) checks them.
"""
function tskit_tables(ped::Pedigree)
    now = pedigree_time(ped)
    nn = length(ped.node_time)
    smp = samples(ped)
    is_sample = zeros(Int, nn); is_sample[smp.node] .= 1
    age(t) = now - Int(t)
    nodes = (is_sample = is_sample, time = age.(ped.node_time),
             population = [Int(ped.ind_deme[i]) - 1 for i in ped.node_ind], individual = Int.(ped.node_ind) .- 1)
    eo = sortperm([(age(ped.node_time[ped.parent[e]]), ped.parent[e], ped.child[e], ped.left[e]) for e in eachindex(ped.left)])
    edges = (left = Int.(ped.left[eo]), right = Int.(ped.right[eo]), parent = Int.(ped.parent[eo]) .- 1,
             child = Int.(ped.child[eo]) .- 1)
    individuals = (flags = zeros(Int, length(ped.ind_time)),
                   parents = [(Int(ped.ind_p1[i]) - 1, Int(ped.ind_p2[i]) - 1) for i in eachindex(ped.ind_time)],
                   time = age.(ped.ind_time), birth_deme = Int.(ped.ind_birth) .- 1, deme = Int.(ped.ind_deme) .- 1)
    populations = (metadata = ["{\"name\":\"deme_$a\"}" for a in eachindex(ped.worlds)],)
    # sites and mutations, site by site, in the order the mutations happened (oldest first)
    nt = ped.ntraits
    loci = sort!(unique(ped.mut_locus))
    siteof = Dict(l => i - 1 for (i, l) in enumerate(loci))
    bysite = [Int[] for _ in loci]
    for m in eachindex(ped.mut_node)
        push!(bysite[siteof[ped.mut_locus[m]] + 1], m)
    end
    msite, mnode, mtime, mparent = Int[], Int[], Int[], Int[]
    mstate, meffect = String[], String[]
    value = Vector{Float64}[]                     # the allele after each output mutation
    lastmut = zeros(Int, nn)                      # the last output mutation on each node at this site (1-based), 0 if none
    for (s, ms) in enumerate(bysite)
        tree = marginal_tree(ped, loci[s])
        for m in ms
            u = ped.mut_node[m]
            q = u
            while q != 0 && lastmut[q] == 0       # the nearest mutation at or above the node
                q = tree[q]
            end
            par = q == 0 ? 0 : lastmut[q]
            eff = ped.mut_effect[(m-1)*nt+1:m*nt]
            push!(value, par == 0 ? 0.0 .+ eff : value[par] .+ eff)
            push!(msite, s - 1); push!(mnode, u - 1); push!(mtime, age(ped.mut_time[m])); push!(mparent, par - 1)
            push!(mstate, join(repr.(value[end]), ",")); push!(meffect, join(repr.(eff), ","))
            lastmut[u] = length(msite)
        end
        for m in ms
            lastmut[ped.mut_node[m]] = 0
        end
    end
    sites = (position = Int.(loci) .- 1, ancestral_state = fill("0", length(loci)))
    mutations = (site = msite, node = mnode, derived_state = mstate, time = mtime, parent = mparent, effect = meffect)
    # migrations: each node of the individual, the most recent first
    gl, gr, gn, gs, gd = Int[], Int[], Int[], Int[], Int[]
    gt = Float64[]
    nodes_of = [Int[] for _ in ped.ind_time]
    for n in 1:nn
        push!(nodes_of[ped.node_ind[n]], n)
    end
    for e in length(ped.mig_ind):-1:1             # the table is in order of time: the most recent first
        for n in nodes_of[ped.mig_ind[e]]
            push!(gl, 0); push!(gr, ped.L); push!(gn, n - 1)
            push!(gs, ped.mig_to[e] - 1); push!(gd, ped.mig_from[e] - 1); push!(gt, age(ped.mig_time[e]) + 0.5)
        end
    end
    migrations = (left = gl, right = gr, node = gn, source = gs, dest = gd, time = gt)
    (; sequence_length = ped.L, nodes, edges, individuals, populations, sites, mutations, migrations)
end

const _B64 = collect("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/")

# base64 (RFC 4648), for tskit's population metadata; the Base64 standard library is not
# a dependency
function _base64(s::AbstractString)
    b = codeunits(s)
    io = IOBuffer()
    for i in 1:3:length(b)
        n = min(3, length(b) - i + 1)
        x = (UInt32(b[i]) << 16) | (n > 1 ? UInt32(b[i+1]) << 8 : UInt32(0)) | (n > 2 ? UInt32(b[i+2]) : UInt32(0))
        print(io, _B64[(x >> 18) & 63 + 1], _B64[(x >> 12) & 63 + 1],
              n > 1 ? _B64[(x >> 6) & 63 + 1] : '=', n > 2 ? _B64[x & 63 + 1] : '=')
    end
    String(take!(io))
end

function _unbase64(s::AbstractString)
    val = Dict(c => UInt32(i - 1) for (i, c) in enumerate(_B64))
    out = UInt8[]
    cs = collect(s)
    for i in 1:4:length(cs)
        q = cs[i:i+3]
        x = sum(get(val, q[j], UInt32(0)) << (6 * (4 - j)) for j in 1:4)
        push!(out, UInt8((x >> 16) & 0xff))
        q[3] == '=' || push!(out, UInt8((x >> 8) & 0xff))
        q[4] == '=' || push!(out, UInt8(x & 0xff))
    end
    String(out)
end

const _TSKIT_FILES = (:nodes, :edges, :individuals, :populations, :sites, :mutations, :migrations)

_cell(x::Tuple) = join(x, ",")
_cell(x::AbstractFloat) = isinteger(x) ? string(Int(x)) : string(x)
_cell(x) = string(x)

"""
    write_tskit(dir, ped) -> tables
    write_tskit(dir, tables)

Write the pedigree's tree-sequence tables ([`tskit_tables`](@ref)) as text files in the
directory `dir`: `nodes.txt`, `edges.txt`, `individuals.txt`, `populations.txt`,
`sites.txt`, `mutations.txt`, `migrations.txt`, and `sequence_length.txt` (the number of
loci). Each is tab-separated with a header line of column names, in the text format of
tskit's `load_text` (https://tskit.dev/tskit/docs/stable/file-formats.html#text-file-formats):
an `id` column first (tskit ignores it), `-1` for no reference, an individual's parents
as `p1,p2`, the population metadata base64-encoded. In Python,

    tskit.load_text(nodes=open("nodes.txt"), edges=open("edges.txt"), sites=open("sites.txt"),
                    mutations=open("mutations.txt"), individuals=open("individuals.txt"),
                    populations=open("populations.txt"), migrations=open("migrations.txt"),
                    sequence_length=L, strict=False)

should load them. **This has not been tried**: the files were written from tskit's
documentation, on a machine without Python or tskit. What is checked, in Julia, is the
tables' internal consistency ([`check_tables`](@ref)) and that reading the files back
([`read_tskit`](@ref)) gives the same tables and the same marginal trees as the pedigree.
"""
function write_tskit(dir::AbstractString, tables::NamedTuple)
    mkpath(dir)
    for name in _TSKIT_FILES
        tab = getfield(tables, name)
        cols = keys(tab)
        n = length(first(tab))
        open(joinpath(dir, "$name.txt"), "w") do io
            println(io, join(("id", cols...), '\t'))
            for i in 1:n
                print(io, i - 1)
                for c in cols
                    x = tab[c][i]
                    print(io, '\t', name === :populations && c === :metadata ? _base64(x) : _cell(x))
                end
                println(io)
            end
        end
    end
    write(joinpath(dir, "sequence_length.txt"), string(tables.sequence_length, "\n"))
    tables
end
write_tskit(dir::AbstractString, ped::Pedigree) = write_tskit(dir, tskit_tables(ped))

"""
    read_tskit(dir) -> tables

Read back the text files [`write_tskit`](@ref) wrote, as the tables of
[`tskit_tables`](@ref).
"""
function read_tskit(dir::AbstractString)
    strcols = (:derived_state, :ancestral_state, :effect, :metadata)
    function table(name)
        lines = readlines(joinpath(dir, "$name.txt"))
        cols = Symbol.(split(lines[1], '\t'))[2:end]
        rows = [split(l, '\t')[2:end] for l in lines[2:end]]
        NamedTuple{Tuple(cols)}(Tuple(map(enumerate(cols)) do (j, c)
            raw = String[r[j] for r in rows]
            c === :metadata ? _unbase64.(raw) :
            c in strcols ? raw :
            c === :parents ? [Tuple(parse.(Int, split(x, ','))) for x in raw] :
            (name === :migrations && c === :time) ? parse.(Float64, raw) : parse.(Int, raw)
        end))
    end
    (; sequence_length = parse(Int, strip(read(joinpath(dir, "sequence_length.txt"), String))),
     (name => table(name) for name in _TSKIT_FILES)...)
end

"""
    check_tables(tables; complete = false) -> true
    check_tables(ped)

Check tree-sequence tables ([`tskit_tables`](@ref), [`read_tskit`](@ref)) against the
requirements of tskit's data model
(https://tskit.dev/tskit/docs/stable/data-model.html#valid-tree-sequence-requirements),
throwing an `ErrorException` that names the first one broken:

- nodes: a finite time, a valid population and individual;
- individuals: valid parents, none its own parent;
- edges: `0 ≤ left < right ≤ L`, valid nodes, the parent older than the child; sorted by
  the parent's time (the youngest first), the edges of a parent together and sorted by
  child and `left` (so none is repeated); the intervals over which a node is a child
  disjoint, so that every locus of a child has at most one parent;
- with `complete = true` (the default for a pedigree that was never simplified): the
  intervals of every node born after the founders cover `[0, L)` exactly, each locus once;
- sites: positions in `[0, L)`, increasing; mutations: a valid site and node, sorted by
  site, not younger than their node, and after their parent mutation, which is at the same
  site and not younger;
- migrations: valid nodes and populations, sorted by time.
"""
function check_tables(tb::NamedTuple; complete::Bool = false)
    L = tb.sequence_length
    nd, ed, ind, mu, mg = tb.nodes, tb.edges, tb.individuals, tb.mutations, tb.migrations
    nn, ni, np = length(nd.time), length(ind.flags), length(tb.populations.metadata)
    bad(msg) = error("invalid tables: " * msg)
    for n in 1:nn
        (isfinite(nd.time[n]) && nd.time[n] >= 0) || bad("node $(n-1) has time $(nd.time[n])")
        -1 <= nd.population[n] < np || bad("node $(n-1) has population $(nd.population[n])")
        -1 <= nd.individual[n] < ni || bad("node $(n-1) has individual $(nd.individual[n])")
    end
    for i in 1:ni, q in ind.parents[i]
        -1 <= q < ni || bad("individual $(i-1) has parent $q")
        q == i - 1 && bad("individual $(i-1) is its own parent")
    end
    ne = length(ed.left)
    seen = Set{Int}()
    for e in 1:ne
        l, r, q, c = ed.left[e], ed.right[e], ed.parent[e], ed.child[e]
        0 <= l < r <= L || bad("edge $(e-1) has interval [$l, $r)")
        (0 <= q < nn && 0 <= c < nn) || bad("edge $(e-1) names nodes $q and $c")
        nd.time[q+1] > nd.time[c+1] || bad("edge $(e-1): parent $q is not older than child $c")
        e == 1 && continue
        q0, c0, l0 = ed.parent[e-1], ed.child[e-1], ed.left[e-1]
        if q != q0
            q in seen && bad("the edges of parent $q are not together")
            push!(seen, q0)
            nd.time[q+1] >= nd.time[q0+1] || bad("edges are not sorted by the parent's time at edge $(e-1)")
        else
            (c0, l0) < (c, l) || bad("the edges of parent $q are not sorted by child and left at edge $(e-1)")
        end
    end
    # every locus of a child has at most one parent (exactly one, if complete)
    bychild = sortperm([(ed.child[e], ed.left[e]) for e in 1:ne])
    covered = zeros(Int, nn)
    for (j, e) in enumerate(bychild)
        c = ed.child[e]
        if j > 1 && ed.child[bychild[j-1]] == c
            ed.right[bychild[j-1]] <= ed.left[e] || bad("node $c inherits loci from two parents (edges overlap)")
            complete && ed.right[bychild[j-1]] != ed.left[e] && bad("node $c has a gap in its ancestry")
        elseif complete
            ed.left[e] == 0 || bad("the ancestry of node $c does not start at 0")
        end
        covered[c+1] += ed.right[e] - ed.left[e]
    end
    if complete
        tmax = maximum(nd.time; init = 0)
        for n in 1:nn
            want = nd.time[n] == tmax ? 0 : L
            covered[n] == want || bad("node $(n-1), of age $(nd.time[n]), has parents over $(covered[n]) of $L loci")
        end
    end
    pos = tb.sites.position
    all(0 <= x < L for x in pos) || bad("a site is outside [0, $L)")
    issorted(pos; lt = <=) || bad("sites are not in increasing order of position")
    for m in eachindex(mu.site)
        0 <= mu.site[m] < length(pos) || bad("mutation $(m-1) has site $(mu.site[m])")
        0 <= mu.node[m] < nn || bad("mutation $(m-1) has node $(mu.node[m])")
        mu.time[m] >= nd.time[mu.node[m]+1] || bad("mutation $(m-1) is younger than its node")
        m > 1 && mu.site[m-1] > mu.site[m] && bad("mutations are not sorted by site")
        q = mu.parent[m]
        q == -1 && continue
        0 <= q < m - 1 || bad("mutation $(m-1) does not come after its parent $q")
        mu.site[q+1] == mu.site[m] || bad("mutation $(m-1) and its parent are at different sites")
        mu.time[q+1] >= mu.time[m] || bad("mutation $(m-1) is older than its parent")
    end
    for g in eachindex(mg.node)
        0 <= mg.node[g] < nn || bad("migration $(g-1) has node $(mg.node[g])")
        (0 <= mg.source[g] < np && 0 <= mg.dest[g] < np) || bad("migration $(g-1) names populations $(mg.source[g]), $(mg.dest[g])")
        (mg.left[g] == 0 && mg.right[g] == L) || bad("migration $(g-1) has interval [$(mg.left[g]), $(mg.right[g]))")
        g > 1 && mg.time[g-1] > mg.time[g] && bad("migrations are not sorted by time")
    end
    true
end
check_tables(ped::Pedigree; complete::Bool = ped.nsimplified == 0) = check_tables(tskit_tables(ped); complete)

"""
    table_tmrca(tables, locus, a, b)

The age of the most recent common ancestor of nodes `a` and `b` (ids from 0) at `locus`
(1-based) in tree-sequence tables, or `nothing`: computed from the tables alone, as a
check on them (it must equal [`tmrca`](@ref) on the pedigree).
"""
function table_tmrca(tb::NamedTuple, locus::Integer, a::Integer, b::Integer)
    ed, time = tb.edges, tb.nodes.time
    x = locus - 1
    up(n) = (e = findfirst(e -> ed.child[e] == n && ed.left[e] <= x < ed.right[e], eachindex(ed.left));
             e === nothing ? -1 : ed.parent[e])
    while a != b
        (a == -1 || b == -1) && return nothing
        if time[a+1] < time[b+1] || (time[a+1] == time[b+1] && a > b)
            a = up(a)
        else
            b = up(b)
        end
    end
    time[a+1]
end
