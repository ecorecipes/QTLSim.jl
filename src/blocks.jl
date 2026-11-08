# Chromosomes made of shared blocks.
#
# The L loci of a chromosome are split into nb consecutive blocks of B loci (the last
# block may be shorter). A chromosome is a vector of nb block ids into a pool. A block is
# immutable while it is shared: copying a chromosome, or the part of it on one side of a
# crossover, copies block ids and counts references, and a block is copied only when it
# must change (a crossover inside it, or a mutation in a block that someone else also
# holds). Each block caches the sum of its values, so a genotype is a sum of nb·Ploidy
# cached sums rather than of L·Ploidy alleles.
#
# With `n` traits (`Params` `NoTraits`) every locus holds a vector of `n` allelic
# effects, one per trait, stored locus by locus (the `n` values of a locus together), so
# that a run of loci is a contiguous run of values, and a block caches `n` sums. With
# `n = 1` the layout is the scalar one.

"""
    default_blocksize(L)

The number of loci per block for `L` loci, `⌈√L⌉`. Copying a chromosome then costs `√L`
block ids, and building the block that holds a crossover, or copying one before a
mutation, `√L` values.
"""
default_blocksize(L::Integer) = max(1, ceil(Int, sqrt(L)))

"""
    BlockPool(L, B; traits = 1)

The blocks of chromosomes of `L` loci, `B` loci per block. Block `b` holds its values in
`vals[(b-1)B+1 : (b-1)B+len]`, where `len` is the length of the block position it is
used at (`B`, or less for the last position), its reference count in `refs[b]` and the
sum of its values in `sums[b]`. With `n = traits` traits a locus holds `n` values, its
effects on each trait, so block `b` holds `n·len` values from `(b-1)nB+1`, locus `t`'s
effect on trait `r` at `(b-1)nB + (t-1)n + r`, and its `n` sums, one per trait, at
`sums[(b-1)n+1 : bn]`. Free blocks (count 0) are on the stack `free`. The pool
doubles when it runs out of blocks, so a run allocates only until it reaches its
steady size.

The sum of a block is always recomputed from its values, left to right, after they
change ([`blocksum!`](@ref)), so it is a function of the values alone, not of their
history.
"""
mutable struct BlockPool
    L::Int                  # loci per chromosome
    B::Int                  # loci per block (the last block may be shorter)
    nb::Int                 # blocks per chromosome
    n::Int                  # traits: values per locus
    vals::Vector{Float64}   # B values per block, block b at (b-1)B+1
    refs::Vector{Int32}     # reference counts
    sums::Vector{Float64}   # cached sum of each block's values
    free::Vector{Int32}     # ids of free blocks, a stack
end

function BlockPool(L::Integer, B::Integer = default_blocksize(L); capacity::Integer = 64,
                   traits::Integer = 1)
    L >= 1 || throw(ArgumentError("a chromosome needs at least one locus, got $L"))
    1 <= B <= L || throw(ArgumentError("the block size must be in 1:$L, got $B"))
    traits >= 1 || throw(ArgumentError("a locus affects at least one trait, got $traits"))
    pool = BlockPool(L, B, cld(L, B), traits, Float64[], Int32[], Float64[], Int32[])
    grow!(pool, capacity)
    pool
end

"The number of loci in block position `k`."
@inline blocklen(pool::BlockPool, k::Integer) = k < pool.nb ? pool.B : pool.L - (pool.nb - 1) * pool.B

"The offset of block `id`'s values in `pool.vals`: value `t` is at `offset + t` (with `n` traits, locus `t`'s effect on trait `r` at `offset + (t-1)n + r`)."
@inline blockoffset(pool::BlockPool, id::Integer) = (Int(id) - 1) * pool.B * pool.n

# Grow the pool to at least `n` blocks, doubling, and put the new ids on the free stack,
# lowest on top.
function grow!(pool::BlockPool, n::Integer = 2 * length(pool.refs))
    old = length(pool.refs)
    new = max(n, 2 * old, 16)
    resize!(pool.vals, new * pool.B * pool.n)
    resize!(pool.refs, new)
    resize!(pool.sums, new * pool.n)
    fill!(view(pool.refs, old+1:new), Int32(0))
    sizehint!(pool.free, new)
    for id in new:-1:old+1
        push!(pool.free, Int32(id))
    end
    pool
end

"""
    newblock!(pool)

A free block with reference count 1; its values and sum are for the caller to set.
Growing the pool may move `pool.vals`, so take offsets into it only after this call.
"""
@inline function newblock!(pool::BlockPool)
    isempty(pool.free) && grow!(pool)
    id = pop!(pool.free)
    @inbounds pool.refs[id] = Int32(1)
    id
end

"Take one more reference to block `id`."
@inline retain!(pool::BlockPool, id::Int32) = (@inbounds pool.refs[id] += Int32(1); id)

"Drop one reference to block `id`, freeing it when none is left."
@inline function release!(pool::BlockPool, id::Int32)
    @inbounds r = pool.refs[id] - Int32(1)
    r < 0 && throw(ErrorException("block $id released more often than it was taken"))
    @inbounds pool.refs[id] = r
    r == 0 && push!(pool.free, id)
    nothing
end

"""
    blocksum!(pool, id, len)

Recompute the cached sum of block `id`, of length `len` loci, from its values, left to
right (with several traits, one sum per trait, each left to right).
"""
@inline function blocksum!(pool::BlockPool, id::Int32, len::Integer)
    pool.n == 1 || return blocksums!(pool, id, len)
    o = blockoffset(pool, id)
    s = 0.0
    @inbounds for t in 1:len
        s += pool.vals[o + t]
    end
    @inbounds pool.sums[id] = s
    s
end

# blocksum! with n > 1 traits: sum r of block `id` is that of its values for trait r
function blocksums!(pool::BlockPool, id::Int32, len::Integer)
    n, o = pool.n, blockoffset(pool, id)
    so = (Int(id) - 1) * n
    @inbounds for r in 1:n
        s = 0.0
        for t in 1:len
            s += pool.vals[o + (t - 1) * n + r]
        end
        pool.sums[so + r] = s
    end
    0.0
end

"Point slot `(k, c, i)` of `chr` (or linear slot `idx`) at block `id`, whose reference the caller has taken, and release the block it named."
@inline setslot!(pool::BlockPool, chr::Array{Int32,3}, k, c, i, id::Int32) =
    setslot!(pool, chr, LinearIndices(chr)[k, c, i], id)
@inline function setslot!(pool::BlockPool, chr::Array{Int32,3}, idx::Int, id::Int32)
    @inbounds old = chr[idx]
    @inbounds chr[idx] = id
    release!(pool, old)
end

"Point linear slot `idx` of `chr` at block `id`, taking a reference to it, unless the slot already names it."
@inline function share!(pool::BlockPool, chr::Array{Int32,3}, idx::Int, id::Int32)
    @inbounds old = chr[idx]
    old == id && return nothing
    retain!(pool, id)
    @inbounds chr[idx] = id
    release!(pool, old)
end

"A new block of length `len` (loci) holding a copy of block `src` of `from` (possibly another pool with as many traits)."
@inline function copyblock!(pool::BlockPool, from::BlockPool, src::Int32, len::Integer)
    id = newblock!(pool)
    n = pool.n
    copyto!(pool.vals, blockoffset(pool, id) + 1, from.vals, blockoffset(from, src) + 1, len * n)
    if n == 1
        @inbounds pool.sums[id] = from.sums[src]   # the same values in the same order: the same sum
    else
        copyto!(pool.sums, (Int(id) - 1) * n + 1, from.sums, (Int(src) - 1) * n + 1, n)
    end
    id
end

"""
    check_blocks(chrs...)

Check a pool's reference counts against the block ids in the buffers `chrs` (arrays of
ids, as `w.chr` and `w.schr`), which must be every buffer that uses the pool: each count
is the number of slots naming the block, and the free blocks are exactly those with count
0, each once. Returns the number of live blocks. For tests.
"""
function check_blocks(pool::BlockPool, chrs::AbstractArray{Int32}...)
    n = length(pool.refs)
    count = zeros(Int, n)
    for chr in chrs, id in chr
        1 <= id <= n || error("slot names block $id, not in the pool of $n")
        count[id] += 1
    end
    for id in 1:n
        count[id] == pool.refs[id] ||
            error("block $id has $(pool.refs[id]) references counted but $(count[id]) slots name it")
    end
    free = sort(Int.(pool.free))
    free == findall(==(0), count) || error("the free list is not the set of unreferenced blocks")
    n - length(free)
end
