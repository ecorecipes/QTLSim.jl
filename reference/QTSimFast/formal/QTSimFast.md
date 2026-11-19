```lean
import Mathlib
```

# QTSimFast: the algorithmic changes preserve the output

QTSim (its source calls itself QTLSim 1.0; Andrew Meade, University of Reading, 2023) is
the C simulation behind Pagel, Gardner and Meade (Scientific Reports, 2025). QTSimFast is
a modified copy that runs faster and is meant to write the same bytes. This file proves,
for a mathematical model of the two programs, that each algorithmic change gives the same
result as the original algorithm, and that the changes compose: one generation, and then
$n$ generations, of the fast model produce the same output records as the original model
and leave the random generator in the same state.

**What is proved, and what is not.** The theorems are about the definitions in this file.
They are not about the C source, the compiler, the C library or GSL: there is no verified
C semantics here. That the definitions describe the C code is supported by inspection
(each definition names the C function it models), by the executable cross-check of
`formal/model_check.c`, and by the parity tests of the two programs; it is not proved.

**How floating point is modelled.** The type `F` of alleles and sums is arbitrary, with an
arbitrary binary operation `add : F → F → F` and an arbitrary constant `zero`. No
associativity, commutativity or identity law is assumed. A theorem that two procedures
return equal values of `F` therefore holds when `F` is the set of IEEE doubles and `add`
is the processor's rounded addition: it is a statement about bit patterns. Had an
optimisation changed the order of the additions, the theorem would be false for some
`add`, and could not have been proved (the companion file `QTSimFastShadows.lean` exhibits
such counterexamples).

**How the random generator is modelled.** A draw is an arbitrary function
$S \to \alpha \times S$ on an arbitrary type $S$ of generator states (`Rand S α`).
Nothing is assumed about the numbers a draw returns, with one exception, stated where it
is used: a mutation's locus is below the number of loci (`hpos` in `mutLoop_sim`).

**The targets.**

| | Change | Main theorems |
|---|---|---|
| (a) | genotypes summed block position by block position, all individuals side by side | `foldl_pointwise`, `fastGeno_eq` |
| (b) | genomes as shared blocks, copied when written; unused blocks collected | `fillMix_spec`, `mutStep_sim`, `collect_inv` |
| (c) | the lowest crossover point kept, no map built or sorted | `origLoop_eq`, `runMin_eq_head`, `fastLow_eq` |
| (d) | all the draws of `NewPop` first, then the copying | `interleaved_eq_two_pass` |
| (e) | statistics computed once, dead computations removed | not a theorem: by inspection |
| (f) | one generation and $n$ generations | `gen_refines`, `run_refines`, `init_rel` |

The model covers the diploid path with at least one individual and one locus. The haploid
path (a gamete is a copy of one survivor's chromosome: `Cut = NB`, no crossover) and the
special cases `NoQTL = 0` and `PopSize = 0` are not modelled; they are covered by the
parity tests only.

```lean
namespace QTSimFast
```

## 0. Computations over an arbitrary generator

Every loop of `NewPop` and `MutatePop` has the same shape: for each item, in order, make
some draws and change some data. `loopM` is that shape. `interleaved` is the special case
in which the body first draws, then applies what was drawn; `drawAll` followed by
`applyAll` is the two-pass form of the fast program.

**Informal claim (d).** In the original, `NewPop` handles one offspring at a time: it
draws the parents, the crossover points and the strand of each gamete, copies the two
gametes locus by locus, and draws the environmental deviation. In the fast program all
these draws are made first, for all offspring, in the same order and with the same
arguments, and the genomes are copied afterwards. If no draw reads a genome (neither its
arguments nor the number of numbers it consumes depend on one), the two orders leave the
generator in the same state and produce the same decisions, hence the same genomes.

```lean
/-- A computation that reads and advances a generator state `S` and returns an `α`: the
state monad over the generator, written out. Nothing is assumed about `S` or about what a
computation returns. -/
abbrev Rand (S : Type) (α : Type) : Type := S → α × S

section Loops

variable {S ι D G G₁ G₂ : Type}

/-- A loop over `items` whose body reads and changes some data `g : G` and advances the
generator (the shape of every loop of `NewPop` and `MutatePop`). -/
def loopM (act : ι → G → Rand S G) : List ι → G → Rand S G
  | [], g => fun s => (g, s)
  | i :: is, g => fun s => loopM act is (act i g s).1 (act i g s).2

/-- All the draws of a loop, made in order, with no data in hand (`mapM` in the state
monad over the generator): the first pass of the fast program. -/
def drawAll (draw : ι → Rand S D) : List ι → Rand S (List D)
  | [] => fun s => ([], s)
  | i :: is => fun s =>
      ((draw i s).1 :: (drawAll draw is (draw i s).2).1, (drawAll draw is (draw i s).2).2)

/-- The drawn decisions applied in order, with no generator: the second pass. -/
def applyAll (apply : ι → D → G → G) : List ι → List D → G → G
  | i :: is, d :: ds, g => applyAll apply is ds (apply i d g)
  | _, _, g => g

/-- The original's shape: for each item, draw (the draw is handed the data, so that it
*could* depend on it), then apply at once. -/
def interleaved (draw : ι → G → Rand S D) (apply : ι → D → G → G) : List ι → G → Rand S G :=
  loopM (fun i g s => (apply i (draw i g s).1 g, (draw i g s).2))

/-- **The hypothesis of the two-pass theorem**: what an item's draws return, and the
generator state they leave, do not depend on the data (the genomes). In the C source this
is the statement that between `SetSurviving` and the end of `NewPop`'s loop no argument of
a GSL call and no branch around one reads `Genome[][]`: the calls are
`gsl_rng_uniform_int(rng, CPop->NoSPop)` twice per attempt in `NewIn` (retried while the
two survivors are the same individual, a comparison of pointers; `AssortativeMating` is
−1, so `ValidDiploidParents` reads no phenotype), `gsl_ran_poisson(rng,
RecombPoissonMean)` and `gsl_rng_uniform_int(rng, NoQTL-1)` in `CreateCombMap` (retried
on `PointInMap`, which reads the map only), `gsl_rng_uniform(rng) < 0.5` in
`RecombinationPoissonChr`, and `RandNormalVar(rng, EnvVar)` in `SetIndividualEnv`
(`reference/QTSim/src/QTSim.c`, lines 225–236, 479–492, 591–636, 667–697, 707–740,
755–788). The survivor list and its length are fixed before the loop. -/
def DrawsIgnoreData (draw : ι → G → Rand S D) : Prop :=
  ∀ i g g' s, draw i g s = draw i g' s

/-- **Two passes (target d).** If the draws do not depend on the data, the interleaved
loop equals: all the draws first, in the same order, then the applications. The generator
ends in the same state and the data are equal. -/
theorem interleaved_eq_two_pass {draw : ι → G → Rand S D} (h : DrawsIgnoreData draw)
    (apply : ι → D → G → G) (g0 : G) (items : List ι) (g : G) (s : S) :
    interleaved draw apply items g s =
      (applyAll apply items (drawAll (fun i => draw i g0) items s).1 g,
        (drawAll (fun i => draw i g0) items s).2) := by
  induction items generalizing g s with
  | nil => rfl
  | cons i is ih =>
    have h1 : draw i g s = draw i g0 s := h i g g0 s
    show loopM _ is _ _ = _
    simp only [drawAll, applyAll, h1]
    exact ih _ _

theorem drawAll_length (draw : ι → Rand S D) (items : List ι) (s : S) :
    (drawAll draw items s).1.length = items.length := by
  induction items generalizing s with
  | nil => rfl
  | cons i is ih => simp [drawAll, ih]

/-- Mapping a function over each draw commutes with drawing them all. -/
theorem drawAll_map {D' : Type} (f : D → D') (draw : ι → Rand S D) (items : List ι) (s : S) :
    drawAll (fun i s => (f (draw i s).1, (draw i s).2)) items s =
      ((drawAll draw items s).1.map f, (drawAll draw items s).2) := by
  induction items generalizing s with
  | nil => rfl
  | cons i is ih => simp [drawAll, ih]

/-- Two loops whose bodies keep a relation between their data, and leave the generator in
the same state, do so as wholes. -/
theorem loopM_sim (Rl : G₁ → G₂ → Prop) (act₁ : ι → G₁ → Rand S G₁) (act₂ : ι → G₂ → Rand S G₂)
    (items : List ι)
    (h : ∀ i ∈ items, ∀ x y s, Rl x y →
      Rl (act₁ i x s).1 (act₂ i y s).1 ∧ (act₁ i x s).2 = (act₂ i y s).2)
    (x : G₁) (y : G₂) (s : S) (hxy : Rl x y) :
    Rl (loopM act₁ items x s).1 (loopM act₂ items y s).1 ∧
      (loopM act₁ items x s).2 = (loopM act₂ items y s).2 := by
  induction items generalizing x y s with
  | nil => exact ⟨hxy, rfl⟩
  | cons i is ih =>
    obtain ⟨h1, h2⟩ := h i (by simp) x y s hxy
    simp only [loopM]
    rw [h2]
    exact ih (fun j hj => h j (by simp [hj])) _ _ _ h1

end Loops
```

## 1. Sums in the original order, regrouped across individuals (target a)

**Informal claim (a).** The original's `CaclGenotype` computes, for each individual on
its own, $G = (\cdots((0 \oplus (a_0 \oplus b_0)) \oplus (a_1 \oplus b_1)) \cdots)
\oplus (a_{L-1} \oplus b_{L-1})$, where $a_l$, $b_l$ are the alleles of its two
chromosomes at locus $l$ and $\oplus$ is the rounded addition of doubles. The fast
program's `SumGenotypesW` keeps one accumulator per individual and advances all of them
over the first block of loci, then all of them over the second block, and so on (inside a
block, eight individuals side by side, two loci per vector step; the individuals may be
divided among threads). Each accumulator receives exactly the same additions in exactly
the same order as in the original, so every individual's genotype is the same double.

The heart of it is that loops over independent accumulators can be interchanged
(`foldl_pointwise`); the rest is the arithmetic of the layout of loci in blocks.

```lean
section Sums

variable {κ ι A T : Type}

/-- **Interchange of independent loops.** Updating the accumulators of all the
individuals `j` at each position `k`, position by position, leaves in each accumulator the
left fold of that individual's own steps. No law of the step is used: with
`step k i a = a + x k i` in IEEE arithmetic the two sides are the same bits. This is the
loop over `K` in `SumGenotypesW` of the fast program (accumulators `Geno[Ind]`), the loop
over the offsets `Index` inside a block with the lanes `J` as `j` (accumulators `Sum[J]`),
and the division of the individuals among threads in `SumGenotypesJob`. -/
theorem foldl_pointwise (step : κ → ι → A → A) (ks : List κ) (init : ι → A) (i : ι) :
    (ks.foldl (fun acc k => fun j => step k j (acc j)) init) i =
      ks.foldl (fun a k => step k i a) (init i) := by
  induction ks generalizing init with
  | nil => rfl
  | cons k ks ih => simp only [List.foldl_cons]; exact ih _

/-- The vector loop of `SumGenotypesW` takes two loci per step, each accumulator adding
its two terms in order: a fold over pairs is the fold over the terms. -/
theorem foldl_pairs (f : A → T → A) (ps : List (T × T)) (z : A) :
    ps.foldl (fun a p => f (f a p.1) p.2) z = (ps.flatMap (fun p => [p.1, p.2])).foldl f z := by
  induction ps generalizing z with
  | nil => rfl
  | cons p ps ih => simp [ih]

theorem foldl_ext_mem (f g : A → T → A) (l : List T) (z : A)
    (h : ∀ a, ∀ t ∈ l, f a t = g a t) : l.foldl f z = l.foldl g z := by
  induction l generalizing z with
  | nil => rfl
  | cons t l ih =>
    simp only [List.foldl_cons]
    rw [h z t (by simp)]
    exact ih _ (fun a u hu => h a u (by simp [hu]))

theorem foldl_range_add (f : A → ℕ → A) (m n : ℕ) (z : A) :
    (List.range (m + n)).foldl f z =
      (List.range n).foldl (fun a o => f a (m + o)) ((List.range m).foldl f z) := by
  rw [List.range_add, List.foldl_append, List.foldl_map]

theorem foldl_range_mul (f : A → ℕ → A) (B K : ℕ) (z : A) :
    (List.range (K * B)).foldl f z =
      (List.range K).foldl
        (fun a k => (List.range B).foldl (fun a o => f a (k * B + o)) a) z := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Nat.succ_mul, foldl_range_add, ih, List.range_succ, List.foldl_append]
    rfl

end Sums

/-- How the loci are laid out in blocks: `M` full blocks of `B` loci and a last block of
`last` loci, `0 < last ≤ B` (the fast program's `S->B`, `S->NB = M + 1`, `S->LastLen`,
set in `CreateSim`). -/
structure Layout where
  B : ℕ
  M : ℕ
  last : ℕ
  hB : 0 < B
  hlast : 0 < last
  hle : last ≤ B

namespace Layout

variable (y : Layout)

/-- Blocks per chromosome (`S->NB`). -/
def NB : ℕ := y.M + 1

/-- Loci per chromosome (`NoQTL`, `S->L`). -/
def L : ℕ := y.M * y.B + y.last

/-- `BlockLen`. -/
def blockLen (k : ℕ) : ℕ := if k = y.M then y.last else y.B

theorem blockLen_le (k : ℕ) : y.blockLen k ≤ y.B := by
  unfold blockLen; split
  · exact y.hle
  · exact le_rfl

theorem div_lt_NB {l : ℕ} (hl : l < y.L) : l / y.B < y.NB := by
  rw [Nat.div_lt_iff_lt_mul y.hB]
  have h1 := y.hle
  have h2 : (y.M + 1) * y.B = y.M * y.B + y.B := Nat.succ_mul _ _
  unfold L at hl; unfold NB
  omega

theorem foldl_range_L {A : Type} (f : A → ℕ → A) (z : A) :
    (List.range y.L).foldl f z =
      (List.range y.NB).foldl
        (fun a k => (List.range (y.blockLen k)).foldl (fun a o => f a (k * y.B + o)) a) z := by
  unfold L NB
  rw [foldl_range_add, foldl_range_mul, List.range_succ, List.foldl_append]
  simp only [List.foldl_cons, List.foldl_nil, blockLen, if_true]
  congr 1
  apply foldl_ext_mem
  intro a k hk
  have : k ≠ y.M := by
    have := List.mem_range.mp hk
    omega
  simp [this]

end Layout

theorem block_div {B k o : ℕ} (ho : o < B) : (k * B + o) / B = k := by
  have hB : 0 < B := by omega
  rw [Nat.add_comm, Nat.add_mul_div_right _ _ hB, Nat.div_eq_of_lt ho, Nat.zero_add]

theorem block_mod {B k o : ℕ} (ho : o < B) : (k * B + o) % B = o := by
  rw [Nat.add_comm, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt ho]

section Genotypes

variable {F : Type} (add : F → F → F) (zero : F)

/-- The dense chromosome a row of block numbers stands for: locus `l` is allele `l % B` of
block `ids (l / B)` (the fast program's `Block(S, Ids[...])[Off]`). -/
def absChr (B : ℕ) (data : ℕ → ℕ → F) (ids : ℕ → ℕ) (l : ℕ) : F := data (ids (l / B)) (l % B)

/-- `CaclGenotype` of the original, diploid: `Genotype = 0`, then
`Genotype += Genome[0][l] + Genome[1][l]` for `l = 0 .. NoQTL-1`. Individual `i` holds
chromosomes `2i` and `2i+1`. `add` and `zero` are arbitrary. -/
def origGeno (L : ℕ) (pop : ℕ → ℕ → F) (i : ℕ) : F :=
  (List.range L).foldl (fun acc l => add acc (add (pop (2 * i) l) (pop (2 * i + 1) l))) zero

/-- One individual's partial sum advanced over the block at position `k`. -/
def blockSum (y : Layout) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (k i : ℕ) (acc : F) : F :=
  (List.range (y.blockLen k)).foldl
    (fun a o => add a (add (data (ids (2 * i) k) o) (data (ids (2 * i + 1) k) o))) acc

/-- `SumGenotypesW` of the fast program: every accumulator starts at 0; for each block
position, every individual's accumulator is advanced over its two blocks there. -/
def fastGenoAll (y : Layout) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) : ℕ → F :=
  (List.range y.NB).foldl (fun acc k => fun i => blockSum add y data ids k i (acc i))
    (fun _ => zero)

/-- **Genotypes (target a).** The position-major sums over shared blocks give, for every
individual, exactly the original's locus-by-locus sum over the dense chromosomes the
blocks stand for, for an arbitrary `add`. -/
theorem fastGeno_eq (y : Layout) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (i : ℕ) :
    fastGenoAll add zero y data ids i =
      origGeno add zero y.L (fun g => absChr y.B data (ids g)) i := by
  unfold fastGenoAll origGeno
  rw [foldl_pointwise, y.foldl_range_L]
  apply foldl_ext_mem
  intro a k _
  unfold blockSum
  apply foldl_ext_mem
  intro a' o ho
  have ho' : o < y.B := lt_of_lt_of_le (List.mem_range.mp ho) (y.blockLen_le k)
  simp only [absChr, block_div ho', block_mod ho']

end Genotypes
```

## 2. The lowest crossover point (target c)

**Informal claim (c).** The original's `CreateCombMap` draws a Poisson number $K$ of
crossover points; if $K = 0$ or $K \ge$ `NoQTL` the map is `[NoQTL]`; otherwise it draws
$K$ distinct points (redrawing a point already in the map), sorts them and appends a 0.
`RecombinationPoissonChr` then copies the parent's chromosomes locus by locus, switching
strand when the locus equals `CombMap[PointIndex]`; `PointIndex` is never advanced, so
the gamete switches strand once, after locus `CombMap[0]`, the smallest point. The fast
program's `LowestCrossover` makes the same draws with the same tests and keeps the running
minimum; it builds and sorts nothing. Claim: the two consume the generator identically,
and the gamete defined by "the starting strand up to the running minimum, the other
strand after it" is the gamete the original's loop writes.

```lean
section Crossover

variable {F S : Type}

/-- One pass of the copying loop of the original's `RecombinationPoissonChr`:
`Chr[Index] = CP[Index]; if(Index == CombMap[PointIndex]) switch CP`. The state is whether
`CP` is the parent's second chromosome, and the gamete so far. `pointIndex` is passed as
the constant 0: the original never advances it. -/
def recombStep (p1 p2 : ℕ → F) (map : List ℕ) (pointIndex : ℕ) (st : Bool × (ℕ → F))
    (idx : ℕ) : Bool × (ℕ → F) :=
  (if idx = map.getD pointIndex 0 then !st.1 else st.1,
    Function.update st.2 idx ((if st.1 then p2 else p1) idx))

/-- The gamete the original's loop writes over `chr0`, from the map and the starting
strand (`start = true`: `CP = P2`). -/
def origLoop (p1 p2 : ℕ → F) (map : List ℕ) (start : Bool) (L : ℕ) (chr0 : ℕ → F) : ℕ → F :=
  ((List.range L).foldl (recombStep p1 p2 map 0) (start, chr0)).2

theorem recomb_fold (p1 p2 : ℕ → F) (map : List ℕ) (start : Bool) (chr0 : ℕ → F) (n : ℕ) :
    ((List.range n).foldl (recombStep p1 p2 map 0) (start, chr0)).1 =
        (if map.getD 0 0 < n then !start else start) ∧
      ∀ l, ((List.range n).foldl (recombStep p1 p2 map 0) (start, chr0)).2 l =
        if l < n then
          (if l ≤ map.getD 0 0 then (if start then p2 else p1) l else (if start then p1 else p2) l)
        else chr0 l := by
  induction n with
  | zero => simp
  | succ n ih =>
    obtain ⟨ih1, ih2⟩ := ih
    rw [List.range_succ, List.foldl_append]
    simp only [List.foldl_cons, List.foldl_nil, recombStep]
    generalize map.getD 0 0 = m at *
    refine ⟨?_, ?_⟩
    · rw [ih1]
      rcases Nat.lt_trichotomy m n with h | h | h
      · simp [h, show n ≠ m by omega, show m < n + 1 by omega]
      · subst h; simp
      · simp [show ¬ m < n by omega, show n ≠ m by omega, show ¬ m < n + 1 by omega]
    · intro l
      by_cases hl : l = n
      · subst hl
        rw [Function.update_self, ih1]
        by_cases h1 : m < l
        · have h2 : ¬ l ≤ m := by omega
          cases start <;> simp [h1, h2]
        · have h2 : l ≤ m := by omega
          cases start <;> simp [h1, h2]
      · rw [Function.update_of_ne hl, ih2 l]
        by_cases h3 : l < n
        · simp [h3, show l < n + 1 by omega]
        · simp [h3, show ¬ l < n + 1 by omega]

/-- **The original's loop switches strand once, after its map's first entry.** Since
`PointIndex` stays 0, locus `l` comes from the starting strand if `l ≤ CombMap[0]` and
from the other strand otherwise, whatever the rest of the map holds; loci from `L` on are
not written. -/
theorem origLoop_eq (p1 p2 : ℕ → F) (map : List ℕ) (start : Bool) (L : ℕ) (chr0 : ℕ → F)
    (l : ℕ) :
    origLoop p1 p2 map start L chr0 l =
      if l < L then
        (if l ≤ map.getD 0 0 then (if start then p2 else p1) l else (if start then p1 else p2) l)
      else chr0 l :=
  (recomb_fold p1 p2 map start chr0 L).2 l

/-- The fast program's `if(Point < Low) Low = Point;`. -/
def stepMin (lo p : ℕ) : ℕ := if p < lo then p else lo

theorem foldl_stepMin_le_init (pts : List ℕ) (lo : ℕ) : pts.foldl stepMin lo ≤ lo := by
  induction pts generalizing lo with
  | nil => exact le_rfl
  | cons p pts ih =>
    simp only [List.foldl_cons]
    refine le_trans (ih _) ?_
    unfold stepMin; split <;> omega

theorem foldl_stepMin_le_mem (pts : List ℕ) (lo : ℕ) : ∀ p ∈ pts, pts.foldl stepMin lo ≤ p := by
  induction pts generalizing lo with
  | nil => simp
  | cons q pts ih =>
    intro p hp
    simp only [List.foldl_cons]
    rcases List.mem_cons.mp hp with rfl | hp
    · refine le_trans (foldl_stepMin_le_init _ _) ?_
      unfold stepMin; split <;> omega
    · exact ih _ p hp

theorem foldl_stepMin_mem (pts : List ℕ) (lo : ℕ) :
    pts.foldl stepMin lo = lo ∨ pts.foldl stepMin lo ∈ pts := by
  induction pts generalizing lo with
  | nil => left; rfl
  | cons q pts ih =>
    simp only [List.foldl_cons]
    rcases ih (stepMin lo q) with h | h
    · rw [h]; unfold stepMin; split
      · right; simp
      · left; rfl
    · right; exact List.mem_cons_of_mem _ h

/-- **The running minimum is the head of the sorted list.** For any list `m` that is a
sorted rearrangement of the drawn points (what `qsort` with `IntComp` leaves in
`CombMap`; nothing else about `qsort` is used), the fast program's running minimum,
started at `lo`, is `min lo m[0]`. -/
theorem runMin_eq_head (m pts : List ℕ) (hperm : m.Perm pts) (hsorted : m.Pairwise (· ≤ ·))
    (hne : pts ≠ []) (lo : ℕ) : pts.foldl stepMin lo = min lo (m.getD 0 0) := by
  cases m with
  | nil => exact absurd (List.Perm.eq_nil (hperm.symm)) hne
  | cons h t =>
    have hh : (h :: t).getD 0 0 = h := rfl
    rw [hh]
    have hmem : h ∈ pts := hperm.subset (by simp)
    have hle : ∀ p ∈ pts, h ≤ p := by
      intro p hp
      rcases List.mem_cons.mp (hperm.symm.subset hp) with rfl | hp'
      · exact le_rfl
      · exact (List.pairwise_cons.mp hsorted).1 p hp'
    apply le_antisymm
    · exact Nat.le_min.mpr ⟨foldl_stepMin_le_init _ _, foldl_stepMin_le_mem _ _ h hmem⟩
    · rcases foldl_stepMin_mem pts lo with h1 | h1
      · rw [h1]; exact Nat.min_le_left _ _
      · exact le_trans (Nat.min_le_right _ _) (hle _ h1)

/-- The points of one gamete as the original's `CreateCombMap` draws them: each from
`fresh`, which stands for the loop `do Point = gsl_rng_uniform_int(rng, NoQTL-1) while
(PointInMap(...))` and is handed the membership test on the points drawn so far. `fresh`
is arbitrary: nothing is assumed about how many numbers it consumes or what it returns. -/
def drawPoints (fresh : (ℕ → Bool) → Rand S ℕ) : ℕ → List ℕ → Rand S (List ℕ)
  | 0, pts => fun s => (pts, s)
  | k + 1, pts => fun s =>
      drawPoints fresh k (pts ++ [(fresh (fun q => decide (q ∈ pts)) s).1])
        (fresh (fun q => decide (q ∈ pts)) s).2

/-- The same draws in the fast program's `LowestCrossover`: the same `fresh` with the same
membership test (its array `S->Points`), keeping the running minimum. -/
def drawPointsMin (fresh : (ℕ → Bool) → Rand S ℕ) : ℕ → List ℕ → ℕ → Rand S ℕ
  | 0, _, low => fun s => (low, s)
  | k + 1, pts, low => fun s =>
      drawPointsMin fresh k (pts ++ [(fresh (fun q => decide (q ∈ pts)) s).1])
        (stepMin low (fresh (fun q => decide (q ∈ pts)) s).1)
        (fresh (fun q => decide (q ∈ pts)) s).2

theorem drawPoints_length (fresh : (ℕ → Bool) → Rand S ℕ) (k : ℕ) (pts : List ℕ) (s : S) :
    (drawPoints fresh k pts s).1.length = pts.length + k := by
  induction k generalizing pts s with
  | zero => rfl
  | succ k ih => simp only [drawPoints]; rw [ih]; simp; omega

/-- The fast program's draws consume the generator exactly as the original's and return
the running minimum of the same points. -/
theorem drawPointsMin_eq (fresh : (ℕ → Bool) → Rand S ℕ) (lo0 : ℕ) (k : ℕ) (pts : List ℕ)
    (s : S) :
    drawPointsMin fresh k pts (pts.foldl stepMin lo0) s =
      ((drawPoints fresh k pts s).1.foldl stepMin lo0, (drawPoints fresh k pts s).2) := by
  induction k generalizing pts s with
  | zero => rfl
  | succ k ih =>
    simp only [drawPointsMin, drawPoints]
    rw [← ih]
    simp [List.foldl_append]

/-- With more than 16 points the fast program tests membership in a table of stamps
(`S->Seen`) instead of scanning `S->Points`: the table marks exactly the points drawn. -/
theorem seenTable_eq_mem (pts : List ℕ) (q : ℕ) :
    (pts.foldl (fun seen p => Function.update seen p true) (fun _ => false)) q =
      decide (q ∈ pts) := by
  suffices h : ∀ (seen : ℕ → Bool),
      (pts.foldl (fun seen p => Function.update seen p true) seen) q =
        (seen q || decide (q ∈ pts)) by
    simpa using h (fun _ => false)
  induction pts with
  | nil => intro seen; simp
  | cons p pts ih =>
    intro seen
    simp only [List.foldl_cons]
    rw [ih]
    by_cases hq : q = p
    · subst hq; simp
    · simp [hq]

/-- What is assumed of `qsort`: its result is a sorted rearrangement of its input. -/
structure SortSpec (sort : List ℕ → List ℕ) : Prop where
  perm : ∀ l, (sort l).Perm l
  sorted : ∀ l, (sort l).Pairwise (· ≤ ·)

/-- The original's `CreateCombMap`: the number of points from `count`
(`gsl_ran_poisson(rng, RecombPoissonMean)`); none, or `NoQTL` or more: the map is
`[NoQTL]`; otherwise the points, sorted, and a final 0. -/
def origMap (count : Rand S ℕ) (fresh : (ℕ → Bool) → Rand S ℕ) (sort : List ℕ → List ℕ)
    (L : ℕ) : Rand S (List ℕ) := fun s =>
  if (count s).1 = 0 ∨ L ≤ (count s).1 then ([L], (count s).2)
  else (sort (drawPoints fresh (count s).1 [] (count s).2).1 ++ [0],
    (drawPoints fresh (count s).1 [] (count s).2).2)

/-- The fast program's `LowestCrossover`: the same draws, the lowest point or `L`. -/
def fastLow (count : Rand S ℕ) (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) : Rand S ℕ := fun s =>
  if (count s).1 = 0 ∨ L ≤ (count s).1 then (L, (count s).2)
  else drawPointsMin fresh (count s).1 [] L (count s).2

/-- **Crossover maps (target c).** `LowestCrossover` leaves the generator where
`CreateCombMap` leaves it, and returns `min L CombMap[0]`: for every locus `l < L`,
`l ≤` the one is `l ≤` the other, which is all that the copying loop asks of the map
(`origLoop_eq`). -/
theorem fastLow_eq {sort : List ℕ → List ℕ} (hs : SortSpec sort) (count : Rand S ℕ)
    (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) (s : S) :
    fastLow count fresh L s =
      (min L ((origMap count fresh sort L s).1.getD 0 0), (origMap count fresh sort L s).2) := by
  unfold fastLow origMap
  by_cases h : (count s).1 = 0 ∨ L ≤ (count s).1
  · simp [h]
  · simp only [h, if_false]
    have key := drawPointsMin_eq fresh L (count s).1 [] (count s).2
    simp only [List.foldl_nil] at key
    rw [key]
    set pts := (drawPoints fresh (count s).1 [] (count s).2).1 with hpts
    have hlen : pts.length = (count s).1 := by
      rw [hpts, drawPoints_length]; simp
    have hne : pts ≠ [] := by
      intro h0
      rw [h0] at hlen
      simp at hlen
      omega
    have hne' : sort pts ≠ [] := by
      intro h0
      have := hs.perm pts
      rw [h0] at this
      exact hne (List.Perm.eq_nil this.symm)
    rw [runMin_eq_head (sort pts) pts (hs.perm pts) (hs.sorted pts) hne L]
    cases hsp : sort pts with
    | nil => exact absurd hsp hne'
    | cons a t => rfl

end Crossover
```

## 3. Shared blocks, copied when written (target b)

**Informal claim (b).** The original stores each chromosome as an array of `NoQTL`
doubles. The fast program stores a chromosome as a row of block numbers; a block holds the
alleles of `B` consecutive loci and may be held by many chromosomes. Claim: at every point
of a generation, each row of block numbers *stands for* (`absChr`) the dense chromosome
the original holds, where

* a gamete copies block numbers from the parent's two chromosomes, and only the block the
  crossover falls in is built anew (`FillGametes`, `MixBlocks`);
* a mutation changes a block in place only if this chromosome alone holds it (it is
  tagged in `Owner`), and otherwise first copies the block (`MutatePop`): no other
  chromosome is changed;
* blocks that no chromosome of the live population holds are returned to the allocator
  (`CollectBlocks`), which changes no block's contents.

The proof needs an invariant of the store (`MutInv`, `MixInv`): the allocator never hands
out a block in use, and a tagged block is held by exactly one cell of the table.

```lean
section Blocks

variable {F : Type}

/-- The block store of the fast program: `data id o` is allele `o` of block `id`
(`S->Data`), `free` the free list with its top first (`S->Free`), `hw` the high-water mark
(`S->HW`: blocks from `hw` on have never been used), `owner` the tags of the blocks made
this generation (`S->Owner`: chromosome number + 1, or 0). The fast program keeps no
reference counts: a block may be written in place only if its tag names the writer, and
blocks are reclaimed by marking the ones in use (`CollectBlocks`). -/
structure Store (F : Type) where
  data : ℕ → ℕ → F
  free : List ℕ
  hw : ℕ
  owner : ℕ → ℕ

/-- `NewBlock`: the top of the free list, else the block at the high-water mark; the block
is tagged. (The C program stops with an internal error if the capacity is exhausted; the
model has no capacity.) -/
def Store.newBlock (st : Store F) (tag : ℕ) : ℕ × Store F :=
  match st.free with
  | id :: rest => (id, ⟨st.data, rest, st.hw, Function.update st.owner id tag⟩)
  | [] => (st.hw, ⟨st.data, [], st.hw + 1, Function.update st.owner st.hw tag⟩)

/-- The store with other block contents. -/
def Store.withData (st : Store F) (D : ℕ → ℕ → F) : Store F := ⟨D, st.free, st.hw, st.owner⟩

/-- The blocks `NewBlock` may still hand out. -/
def Avail (st : Store F) (id : ℕ) : Prop := id ∈ st.free ∨ st.hw ≤ id

/-- The free list has no repeats and lies below the high-water mark. -/
structure AllocOK (st : Store F) : Prop where
  nodup : st.free.Nodup
  lt : ∀ id ∈ st.free, id < st.hw

theorem newBlock_data (st : Store F) (tag : ℕ) : (st.newBlock tag).2.data = st.data := by
  obtain ⟨data, free, hw, owner⟩ := st
  cases free <;> rfl

theorem newBlock_owner (st : Store F) (tag : ℕ) :
    (st.newBlock tag).2.owner = Function.update st.owner (st.newBlock tag).1 tag := by
  obtain ⟨data, free, hw, owner⟩ := st
  cases free <;> rfl

theorem newBlock_avail (st : Store F) (tag : ℕ) : Avail st (st.newBlock tag).1 := by
  obtain ⟨data, free, hw, owner⟩ := st
  cases free with
  | nil => right; exact le_rfl
  | cons a t => left; simp [Store.newBlock]

theorem newBlock_avail_after {st : Store F} (h : AllocOK st) (tag : ℕ) {id : ℕ}
    (hid : Avail (st.newBlock tag).2 id) : Avail st id ∧ id ≠ (st.newBlock tag).1 := by
  obtain ⟨data, free, hw, owner⟩ := st
  cases free with
  | nil =>
    simp only [Store.newBlock, Avail, List.not_mem_nil, false_or] at hid ⊢
    omega
  | cons a t =>
    have hnd := h.nodup
    have hlt := h.lt
    simp only [Store.newBlock, Avail] at hid hnd hlt ⊢
    rcases hid with hid | hid
    · refine ⟨Or.inl (List.mem_cons_of_mem _ hid), ?_⟩
      rintro rfl
      exact (List.nodup_cons.mp hnd).1 hid
    · refine ⟨Or.inr hid, ?_⟩
      have := hlt a (by simp)
      omega

theorem newBlock_ok {st : Store F} (h : AllocOK st) (tag : ℕ) : AllocOK (st.newBlock tag).2 := by
  obtain ⟨data, free, hw, owner⟩ := st
  cases free with
  | nil => exact ⟨List.nodup_nil, by simp [Store.newBlock]⟩
  | cons a t =>
    have hnd := h.nodup
    have hlt := h.lt
    simp only [Store.newBlock] at hnd hlt ⊢
    exact ⟨(List.nodup_cons.mp hnd).2, fun id hid => hlt id (List.mem_cons_of_mem _ hid)⟩

/-- Block `id` is in use: some chromosome `g < R` of the table holds it at a position
`k < NB`. -/
def Refd (R NB : ℕ) (t : ℕ → ℕ → ℕ) (id : ℕ) : Prop := ∃ g k, g < R ∧ k < NB ∧ t g k = id

/-- Writing blocks that a chromosome does not hold leaves the chromosome it stands for
unchanged. -/
theorem absChr_congr (y : Layout) {data data' : ℕ → ℕ → F} {ids : ℕ → ℕ}
    (h : ∀ k < y.NB, data' (ids k) = data (ids k)) {l : ℕ} (hl : l < y.L) :
    absChr y.B data' ids l = absChr y.B data ids l := by
  unfold absChr; rw [h _ (y.div_lt_NB hl)]

/-- Loci compare as (block, offset) pairs. -/
theorem le_iff_block {B : ℕ} (hB : 0 < B) (l low : ℕ) :
    l ≤ low ↔ l / B < low / B ∨ (l / B = low / B ∧ l % B ≤ low % B) := by
  have h1 := Nat.div_add_mod l B
  have h2 := Nat.div_add_mod low B
  have h3 := Nat.mod_lt l hB
  have h4 := Nat.mod_lt low hB
  generalize l / B = q at *
  generalize l % B = r at *
  generalize low / B = q' at *
  generalize low % B = r' at *
  rcases Nat.lt_trichotomy q q' with h | h | h
  · have h5 : B * (q + 1) ≤ B * q' := Nat.mul_le_mul_left B h
    rw [Nat.mul_add, Nat.mul_one] at h5
    constructor
    · intro _; left; exact h
    · intro _; omega
  · subst h
    constructor
    · intro hle; right; exact ⟨rfl, by omega⟩
    · rintro (h' | ⟨_, h'⟩) <;> omega
  · have h5 : B * (q' + 1) ≤ B * q := Nat.mul_le_mul_left B h
    rw [Nat.mul_add, Nat.mul_one] at h5
    constructor
    · intro hle; omega
    · rintro (h' | ⟨h', _⟩) <;> omega

theorem eq_of_block {B : ℕ} {l p : ℕ} (h1 : l / B = p / B) (h2 : l % B = p % B) : l = p := by
  rw [← Nat.div_add_mod l B, ← Nat.div_add_mod p B, h1, h2]

/-- What the fast program records of a gamete (`RecombinationPoissonChr`): the chromosome
it starts on (`Src1`), the parent's other chromosome (`Src2`), the number of leading
blocks taken whole from the first (`Cut`), and, if the crossover falls inside a block, that
block's position and the last offset taken from the first chromosome (`Mix`). -/
structure GDesc where
  src1 : ℕ
  src2 : ℕ
  cut : ℕ
  mix : Option (ℕ × ℕ)

/-- The record for parent individual `a` (chromosomes `2a`, `2a+1`), lowest crossover
point `low` (`L` for none) and starting strand `s`. -/
def gdesc (y : Layout) (a low : ℕ) (s : Bool) : GDesc :=
  if low = y.L then ⟨2 * a + (if s then 1 else 0), 2 * a + (if s then 0 else 1), y.NB, none⟩
  else if low % y.B + 1 = y.blockLen (low / y.B) then
    ⟨2 * a + (if s then 1 else 0), 2 * a + (if s then 0 else 1), low / y.B + 1, none⟩
  else
    ⟨2 * a + (if s then 1 else 0), 2 * a + (if s then 0 else 1), low / y.B,
      some (low / y.B, low % y.B)⟩

theorem gdesc_src (y : Layout) (a low : ℕ) (s : Bool) :
    (gdesc y a low s).src1 = 2 * a + (if s then 1 else 0) ∧
      (gdesc y a low s).src2 = 2 * a + (if s then 0 else 1) := by
  unfold gdesc; split_ifs <;> exact ⟨rfl, rfl⟩

/-- Without a crossover inside a block, "position below `Cut`" is "locus not after the
crossover point". -/
theorem gdesc_cut_none (y : Layout) (a low : ℕ) (s : Bool) (_hlow : low ≤ y.L)
    (hm : (gdesc y a low s).mix = none) {l : ℕ} (hl : l < y.L) :
    l / y.B < (gdesc y a low s).cut ↔ l ≤ low := by
  unfold gdesc at hm ⊢
  by_cases h1 : low = y.L
  · simp only [h1, if_true]
    exact ⟨fun _ => le_of_lt hl, fun _ => y.div_lt_NB hl⟩
  · by_cases h2 : low % y.B + 1 = y.blockLen (low / y.B)
    · simp only [h1, h2, if_true, if_false]
      rw [le_iff_block y.hB l low]
      constructor
      · intro h
        rcases Nat.lt_or_ge (l / y.B) (low / y.B) with h3 | h3
        · left; exact h3
        · right
          have h4 : l / y.B = low / y.B := by omega
          refine ⟨h4, ?_⟩
          unfold Layout.blockLen at h2
          by_cases h5 : low / y.B = y.M
          · rw [if_pos h5] at h2
            have h6 := Nat.div_add_mod l y.B
            rw [h4, h5, Nat.mul_comm] at h6
            unfold Layout.L at hl
            omega
          · rw [if_neg h5] at h2
            have := Nat.mod_lt l y.hB
            omega
      · rintro (h | ⟨h, _⟩) <;> omega
    · simp [h1, h2] at hm

/-- With a crossover inside a block. -/
theorem gdesc_cut_some (y : Layout) (a low : ℕ) (s : Bool) (_hlow : low ≤ y.L) {b off : ℕ}
    (hm : (gdesc y a low s).mix = some (b, off)) :
    b = low / y.B ∧ off = low % y.B ∧ (gdesc y a low s).cut = b ∧ b < y.NB := by
  unfold gdesc at hm ⊢
  by_cases h1 : low = y.L
  · simp [h1] at hm
  · by_cases h2 : low % y.B + 1 = y.blockLen (low / y.B)
    · simp [h1, h2] at hm
    · simp only [h1, h2, if_false, Option.some.injEq, Prod.mk.injEq] at hm ⊢
      obtain ⟨rfl, rfl⟩ := hm
      exact ⟨rfl, rfl, rfl, y.div_lt_NB (by omega)⟩

/-- `FillGametes`: block numbers copied from the parent's two chromosomes. (The C code
makes this copy position by position with the lists `Head`/`Next` and the array `Cur`,
possibly on several threads; that traversal is not modelled.) -/
def fillRow (cur : ℕ → ℕ → ℕ) (d : GDesc) : ℕ → ℕ :=
  fun k => if k < d.cut then cur d.src1 k else cur d.src2 k

/-- Set one cell of a table of block numbers. -/
def setCell (t : ℕ → ℕ → ℕ) (g b id : ℕ) : ℕ → ℕ → ℕ :=
  Function.update t g (Function.update (t g) b id)

theorem setCell_apply (t : ℕ → ℕ → ℕ) (g b id g' k : ℕ) :
    setCell t g b id g' k = if g' = g ∧ k = b then id else t g' k := by
  unfold setCell
  by_cases hg : g' = g
  · subst hg
    by_cases hk : k = b
    · subst hk; simp
    · simp [hk]
  · simp [hg]

/-- One entry of `MixBlocks`, for gamete `g`: if its crossover fell inside a block and the
parent's two blocks there differ, a new block gets the first chromosome's alleles up to
the offset and the other's after it. -/
def mixStep (cur : ℕ → ℕ → ℕ) (desc : ℕ → GDesc) (x : Store F × (ℕ → ℕ → ℕ)) (g : ℕ) :
    Store F × (ℕ → ℕ → ℕ) :=
  match (desc g).mix with
  | none => x
  | some (b, off) =>
    if cur (desc g).src1 b = cur (desc g).src2 b then x
    else
      ((x.1.newBlock (g + 1)).2.withData
          (Function.update x.1.data (x.1.newBlock (g + 1)).1
            (fun o => if o ≤ off then x.1.data (cur (desc g).src1 b) o
              else x.1.data (cur (desc g).src2 b) o)),
        setCell x.2 g b (x.1.newBlock (g + 1)).1)

/-- The second pass of the fast program's `NewPop`: `FillGametes`, then `MixBlocks`. -/
def fillMix (R : ℕ) (cur : ℕ → ℕ → ℕ) (desc : ℕ → GDesc) (st : Store F) :
    Store F × (ℕ → ℕ → ℕ) :=
  (List.range R).foldl (mixStep cur desc) (st, fun g => fillRow cur (desc g))

/-- The invariant of the store and the next population's table while it is built and
mutated: the allocator is sound, no block it may hand out is in use, and a tagged block is
held by one cell only, a cell of the chromosome its tag names. -/
structure MutInv (y : Layout) (R : ℕ) (x : Store F × (ℕ → ℕ → ℕ)) : Prop where
  ok : AllocOK x.1
  avail : ∀ id, Avail x.1 id → ¬ Refd R y.NB x.2 id
  own : ∀ id, x.1.owner id ≠ 0 → ∀ g k, g < R → k < y.NB → x.2 g k = id →
    x.1.owner id = g + 1 ∧ ∀ k' < y.NB, x.2 g k' = id → k' = k

/-- Giving a new block to one cell keeps the invariant, whatever is written in the block. -/
theorem MutInv.alloc {y : Layout} {R : ℕ} {x : Store F × (ℕ → ℕ → ℕ)} (h : MutInv y R x)
    {g b : ℕ} (hg : g < R) (hb : b < y.NB) (D : ℕ → ℕ → F) :
    MutInv y R ((x.1.newBlock (g + 1)).2.withData D,
      setCell x.2 g b (x.1.newBlock (g + 1)).1) := by
  have hnew : ¬ Refd R y.NB x.2 (x.1.newBlock (g + 1)).1 := h.avail _ (newBlock_avail _ _)
  refine ⟨?_, ?_, ?_⟩
  · have := newBlock_ok h.ok (g + 1)
    exact ⟨this.nodup, this.lt⟩
  · intro id hid
    have hid' : Avail (x.1.newBlock (g + 1)).2 id := hid
    obtain ⟨h1, h2⟩ := newBlock_avail_after h.ok (g + 1) hid'
    rintro ⟨g', k, hg', hk, hcell⟩
    replace hcell : setCell x.2 g b (x.1.newBlock (g + 1)).1 g' k = id := hcell
    rw [setCell_apply] at hcell
    split at hcell
    · exact h2 hcell.symm
    · exact h.avail id h1 ⟨g', k, hg', hk, hcell⟩
  · intro id hid g' k hg' hk hcell
    have hown : ((x.1.newBlock (g + 1)).2.withData D).owner =
        Function.update x.1.owner (x.1.newBlock (g + 1)).1 (g + 1) := newBlock_owner _ _
    replace hid : ((x.1.newBlock (g + 1)).2.withData D).owner id ≠ 0 := hid
    replace hcell : setCell x.2 g b (x.1.newBlock (g + 1)).1 g' k = id := hcell
    show ((x.1.newBlock (g + 1)).2.withData D).owner id = g' + 1 ∧
      ∀ k' < y.NB, setCell x.2 g b (x.1.newBlock (g + 1)).1 g' k' = id → k' = k
    rw [hown] at hid ⊢
    rw [setCell_apply] at hcell
    by_cases hc : g' = g ∧ k = b
    · rw [if_pos hc] at hcell
      obtain ⟨rfl, rfl⟩ := hc
      subst hcell
      refine ⟨by simp, ?_⟩
      intro k' hk' hcell'
      rw [setCell_apply] at hcell'
      by_contra hne
      rw [if_neg (by simp [hne])] at hcell'
      exact hnew ⟨g', k', hg', hk', hcell'⟩
    · rw [if_neg hc] at hcell
      have hne : id ≠ (x.1.newBlock (g + 1)).1 := by
        rintro rfl
        exact hnew ⟨g', k, hg', hk, hcell⟩
      rw [Function.update_of_ne hne] at hid ⊢
      obtain ⟨h1, h2⟩ := h.own id hid g' k hg' hk hcell
      refine ⟨h1, ?_⟩
      intro k' hk' hcell'
      rw [setCell_apply] at hcell'
      split at hcell'
      · exact absurd hcell'.symm hne
      · exact h2 k' hk' hcell'

/-- The invariant while `MixBlocks` runs, after the gametes below `n`: also, no block in
use by the parents has been written, the gametes done stand for the right chromosomes, and
the others still hold the block numbers `FillGametes` gave them. -/
structure MixInv (y : Layout) (R : ℕ) (cur : ℕ → ℕ → ℕ) (desc : ℕ → GDesc)
    (spec : ℕ → ℕ → F) (st0 : Store F) (n : ℕ) (x : Store F × (ℕ → ℕ → ℕ)) : Prop where
  base : MutInv y R x
  availCur : ∀ id, Avail x.1 id → ¬ Refd R y.NB cur id
  curData : ∀ id, Refd R y.NB cur id → x.1.data id = st0.data id
  done : ∀ g, g < n → g < R → ∀ l < y.L, absChr y.B x.1.data (x.2 g) l = spec g l
  todo : ∀ g, n ≤ g → x.2 g = fillRow cur (desc g)

/-- **Gametes over shared blocks (target b, first part).** After `FillGametes` and
`MixBlocks`, every gamete's row of block numbers stands for the chromosome the original's
locus-by-locus copy makes: the starting chromosome up to the lowest crossover point, the
other one after it (`origLoop_eq`). The store's invariant holds, with each new block
tagged for its gamete. -/
theorem fillMix_spec (y : Layout) (R : ℕ) (cur : ℕ → ℕ → ℕ) (st : Store F)
    (a low : ℕ → ℕ) (s : ℕ → Bool) (desc : ℕ → GDesc)
    (hdesc : ∀ g < R, desc g = gdesc y (a g) (low g) (s g))
    (ha : ∀ g < R, 2 * a g + 1 < R) (hlow : ∀ g < R, low g ≤ y.L)
    (hok : AllocOK st) (hav : ∀ id, Avail st id → ¬ Refd R y.NB cur id)
    (hown : ∀ id, st.owner id = 0) :
    MutInv y R (fillMix R cur desc st) ∧
      ∀ g < R, ∀ l < y.L,
        absChr y.B (fillMix R cur desc st).1.data ((fillMix R cur desc st).2 g) l =
          if l ≤ low g then absChr y.B st.data (cur (2 * a g + (if s g then 1 else 0))) l
          else absChr y.B st.data (cur (2 * a g + (if s g then 0 else 1))) l := by
  set spec : ℕ → ℕ → F := fun g l =>
    if l ≤ low g then absChr y.B st.data (cur (2 * a g + (if s g then 1 else 0))) l
    else absChr y.B st.data (cur (2 * a g + (if s g then 0 else 1))) l with hspec
  have hsrc : ∀ g < R, (desc g).src1 < R ∧ (desc g).src2 < R := by
    intro g hg
    rw [hdesc g hg]
    obtain ⟨h1, h2⟩ := gdesc_src y (a g) (low g) (s g)
    rw [h1, h2]
    have := ha g hg
    cases s g <;> simp <;> omega
  have hfillRef : ∀ g < R, ∀ k < y.NB, Refd R y.NB cur (fillRow cur (desc g) k) := by
    intro g hg k hk
    unfold fillRow
    split
    · exact ⟨_, k, (hsrc g hg).1, hk, rfl⟩
    · exact ⟨_, k, (hsrc g hg).2, hk, rfl⟩
  have key : ∀ n, n ≤ R → MixInv y R cur desc spec st n
      ((List.range n).foldl (mixStep cur desc) (st, fun g => fillRow cur (desc g))) := by
    intro n
    induction n with
    | zero =>
      intro _
      refine ⟨⟨hok, ?_, ?_⟩, hav, fun _ _ => rfl, ?_, fun _ _ => rfl⟩
      · rintro id hid ⟨g, k, hg, hk, hcell⟩
        exact hav id hid (hcell ▸ hfillRef g hg k hk)
      · intro id hid; exact absurd (hown id) hid
      · intro g hg; omega
    | succ n ih =>
      intro hn
      have hnR : n < R := by omega
      have I := ih (by omega)
      rw [List.range_succ, List.foldl_append]
      simp only [List.foldl_cons, List.foldl_nil]
      set x := (List.range n).foldl (mixStep cur desc) (st, fun g => fillRow cur (desc g))
        with hx
      have hrow : x.2 n = fillRow cur (desc n) := I.todo n le_rfl
      have hdn := hdesc n hnR
      obtain ⟨hs1, hs2⟩ := gdesc_src y (a n) (low n) (s n)
      -- data of the parents' blocks, as the gamete's own row reads them
      have hcur : ∀ c, c < R → ∀ l < y.L, absChr y.B x.1.data (cur c) l =
          absChr y.B st.data (cur c) l := by
        intro c hc l hl
        exact absChr_congr y (fun k hk => I.curData _ ⟨c, k, hc, hk, rfl⟩) hl
      -- the case in which nothing is allocated
      have hsame : (∀ l < y.L, absChr y.B x.1.data (fillRow cur (desc n)) l = spec n l) →
          MixInv y R cur desc spec st (n + 1) x := by
        intro hgoal
        refine ⟨I.base, I.availCur, I.curData, ?_, fun g hg => I.todo g (by omega)⟩
        intro g hg hgR l hl
        rcases Nat.lt_or_ge g n with h | h
        · exact I.done g h hgR l hl
        · have : g = n := by omega
          subst this
          rw [hrow]; exact hgoal l hl
      unfold mixStep
      cases hmix : (desc n).mix with
      | none =>
        simp only []
        apply hsame
        intro l hl
        have hcut := gdesc_cut_none y (a n) (low n) (s n) (hlow n hnR) (hdn ▸ hmix) hl
        rw [← hdn] at hcut hs1 hs2
        simp only [hspec, absChr, fillRow] at hcur ⊢
        by_cases h : l ≤ low n
        · rw [if_pos (hcut.mpr h), if_pos h, ← hs1]
          exact hcur _ (hsrc n hnR).1 l hl
        · rw [if_neg (fun h' => h (hcut.mp h')), if_neg h, ← hs2]
          exact hcur _ (hsrc n hnR).2 l hl
      | some p =>
        obtain ⟨b, off⟩ := p
        obtain ⟨hb, hoff, hcut, hbNB⟩ :=
          gdesc_cut_some y (a n) (low n) (s n) (hlow n hnR) (hdn ▸ hmix)
        rw [← hdn] at hcut hs1 hs2
        simp only []
        -- which chromosome each locus comes from, in block coordinates
        have hlex : ∀ l, l ≤ low n ↔ l / y.B < b ∨ (l / y.B = b ∧ l % y.B ≤ off) := by
          intro l; rw [hb, hoff]; exact le_iff_block y.hB l (low n)
        by_cases hco : cur (desc n).src1 b = cur (desc n).src2 b
        · rw [if_pos hco]
          apply hsame
          intro l hl
          simp only [hspec]
          by_cases h : l ≤ low n
          · rw [if_pos h, ← hs1, ← hcur _ (hsrc n hnR).1 l hl]
            unfold absChr fillRow
            rw [hcut]
            rcases (hlex l).mp h with h' | ⟨h', _⟩
            · rw [if_pos h']
            · rw [if_neg (by omega), h', hco]
          · rw [if_neg h, ← hs2, ← hcur _ (hsrc n hnR).2 l hl]
            unfold absChr fillRow
            rw [hcut]
            have : ¬ l / y.B < b := fun h' => h ((hlex l).mpr (Or.inl h'))
            rw [if_neg this]
        · rw [if_neg hco]
          set nid := (x.1.newBlock (n + 1)).1 with hnid
          have hnidAvail : Avail x.1 nid := newBlock_avail _ _
          have hnidCur : ¬ Refd R y.NB cur nid := I.availCur _ hnidAvail
          have hnidNxt : ¬ Refd R y.NB x.2 nid := I.base.avail _ hnidAvail
          refine ⟨I.base.alloc hnR hbNB _, ?_, ?_, ?_, ?_⟩
          · intro id hid
            have hid' : Avail (x.1.newBlock (n + 1)).2 id := hid
            exact I.availCur id (newBlock_avail_after I.base.ok (n + 1) hid').1
          · intro id hid
            have hne : id ≠ nid := by rintro rfl; exact hnidCur hid
            show Function.update x.1.data nid _ id = _
            rw [Function.update_of_ne hne]
            exact I.curData id hid
          · intro g hg hgR l hl
            have hlNB := y.div_lt_NB hl
            show absChr y.B (Function.update x.1.data nid _) (setCell x.2 n b nid g) l = _
            unfold absChr
            rw [setCell_apply]
            by_cases hc : g = n ∧ l / y.B = b
            · obtain ⟨rfl, hlb⟩ := hc
              rw [if_pos ⟨rfl, hlb⟩, Function.update_self]
              simp only [hspec]
              have hc1 := I.curData _ ⟨(desc g).src1, b, (hsrc g hnR).1, hbNB, rfl⟩
              have hc2 := I.curData _ ⟨(desc g).src2, b, (hsrc g hnR).2, hbNB, rfl⟩
              by_cases h : l ≤ low g
              · have h' : l % y.B ≤ off := by
                  rcases (hlex l).mp h with h' | ⟨_, h'⟩
                  · omega
                  · exact h'
                rw [if_pos h', if_pos h, hc1, ← hs1]
                unfold absChr; rw [hlb]
              · have h' : ¬ l % y.B ≤ off := fun h' => h ((hlex l).mpr (Or.inr ⟨hlb, h'⟩))
                rw [if_neg h', if_neg h, hc2, ← hs2]
                unfold absChr; rw [hlb]
            · rw [if_neg hc]
              have hne : x.2 g (l / y.B) ≠ nid := fun h => hnidNxt ⟨g, _, hgR, hlNB, h⟩
              rw [Function.update_of_ne hne]
              rcases Nat.lt_or_ge g n with h | h
              · exact I.done g h hgR l hl
              · have hgn : g = n := by omega
                subst hgn
                have hlb : l / y.B ≠ b := fun h' => hc ⟨rfl, h'⟩
                rw [hrow]
                simp only [hspec]
                unfold fillRow
                rw [hcut]
                by_cases h : l ≤ low g
                · have h' : l / y.B < b := by
                    rcases (hlex l).mp h with h' | ⟨h', _⟩
                    · exact h'
                    · exact absurd h' hlb
                  rw [if_pos h', if_pos h, ← hs1, ← hcur _ (hsrc g hnR).1 l hl]
                  rfl
                · have h' : ¬ l / y.B < b := fun h' => h ((hlex l).mpr (Or.inl h'))
                  rw [if_neg h', if_neg h, ← hs2, ← hcur _ (hsrc g hnR).2 l hl]
                  rfl
          · intro g hg
            show setCell x.2 n b nid g = _
            have hgn : g ≠ n := by omega
            unfold setCell
            rw [Function.update_of_ne hgn]
            exact I.todo g (by omega)
  have fin := key R le_rfl
  exact ⟨fin.base, fun g hg l hl => fin.done g hg hg l hl⟩

variable (add : F → F → F)

/-- One mutation in the original (`MutateIndividualChrNo`): `Genome[ChrNo][Pos] += eff`. -/
def mutDense (g pos : ℕ) (eff : F) (pop : ℕ → ℕ → F) : ℕ → ℕ → F :=
  Function.update pop g (Function.update (pop g) pos (add (pop g pos) eff))

/-- One mutation in the fast program (`MutatePop`): if the block is tagged for this
chromosome it is changed in place; otherwise a copy is made, tagged, put in the
chromosome's row and changed. -/
def mutStep (B : ℕ) (g pos : ℕ) (eff : F) (x : Store F × (ℕ → ℕ → ℕ)) :
    Store F × (ℕ → ℕ → ℕ) :=
  if x.1.owner (x.2 g (pos / B)) = g + 1 then
    (x.1.withData (Function.update x.1.data (x.2 g (pos / B))
          (Function.update (x.1.data (x.2 g (pos / B))) (pos % B)
            (add (x.1.data (x.2 g (pos / B)) (pos % B)) eff))), x.2)
  else
    ((x.1.newBlock (g + 1)).2.withData (Function.update x.1.data (x.1.newBlock (g + 1)).1
          (Function.update (x.1.data (x.2 g (pos / B))) (pos % B)
            (add (x.1.data (x.2 g (pos / B)) (pos % B)) eff))),
      setCell x.2 g (pos / B) (x.1.newBlock (g + 1)).1)

/-- **Mutation with copy on write (target b, second part).** One mutation of the fast
program changes the chromosome its row stands for exactly as the original's in-place
`+=` changes the dense chromosome, changes no other chromosome of the population, and
keeps the store's invariant. -/
theorem mutStep_sim (y : Layout) (R : ℕ) {g pos : ℕ} (hg : g < R) (hpos : pos < y.L) (eff : F)
    (x : Store F × (ℕ → ℕ → ℕ)) (pop : ℕ → ℕ → F) (hinv : MutInv y R x)
    (habs : ∀ g' < R, ∀ l < y.L, absChr y.B x.1.data (x.2 g') l = pop g' l) :
    MutInv y R (mutStep add y.B g pos eff x) ∧
      ∀ g' < R, ∀ l < y.L,
        absChr y.B (mutStep add y.B g pos eff x).1.data ((mutStep add y.B g pos eff x).2 g') l =
          mutDense add g pos eff pop g' l := by
  have hb : pos / y.B < y.NB := y.div_lt_NB hpos
  have hold : x.1.data (x.2 g (pos / y.B)) (pos % y.B) = pop g pos := habs g hg pos hpos
  have hdense : ∀ g' l, mutDense add g pos eff pop g' l =
      if g' = g ∧ l = pos then add (pop g pos) eff else pop g' l := by
    intro g' l
    unfold mutDense
    by_cases h1 : g' = g
    · subst h1
      by_cases h2 : l = pos
      · subst h2; simp
      · simp [h2]
    · simp [h1]
  unfold mutStep
  by_cases hown : x.1.owner (x.2 g (pos / y.B)) = g + 1
  · rw [if_pos hown]
    refine ⟨⟨⟨hinv.ok.nodup, hinv.ok.lt⟩, hinv.avail, hinv.own⟩, ?_⟩
    intro g' hg' l hl
    have hlNB := y.div_lt_NB hl
    show Function.update x.1.data _ _ (x.2 g' (l / y.B)) (l % y.B) = _
    rw [hdense]
    by_cases hcell : x.2 g' (l / y.B) = x.2 g (pos / y.B)
    · obtain ⟨h1, _⟩ := hinv.own _ (by omega) g' (l / y.B) hg' hlNB hcell
      have hgg : g' = g := by omega
      subst hgg
      obtain ⟨_, h2⟩ := hinv.own _ (by omega) g' (pos / y.B) hg' hb rfl
      have hlb : l / y.B = pos / y.B := h2 _ hlNB hcell
      rw [hcell, Function.update_self]
      by_cases hlp : l = pos
      · subst hlp
        rw [if_pos ⟨rfl, rfl⟩, Function.update_self, hold]
      · have hmod : l % y.B ≠ pos % y.B := fun h => hlp (eq_of_block hlb h)
        rw [if_neg (by simp [hlp]), Function.update_of_ne hmod, ← hlb]
        exact habs g' hg' l hl
    · rw [Function.update_of_ne hcell]
      have : ¬ (g' = g ∧ l = pos) := by
        rintro ⟨rfl, rfl⟩; exact hcell rfl
      rw [if_neg this]
      exact habs g' hg' l hl
  · rw [if_neg hown]
    set nid := (x.1.newBlock (g + 1)).1 with hnid
    have hnidNxt : ¬ Refd R y.NB x.2 nid := hinv.avail _ (newBlock_avail _ _)
    refine ⟨hinv.alloc hg hb _, ?_⟩
    intro g' hg' l hl
    have hlNB := y.div_lt_NB hl
    show Function.update x.1.data nid _ (setCell x.2 g (pos / y.B) nid g' (l / y.B)) (l % y.B) = _
    rw [hdense, setCell_apply]
    by_cases hc : g' = g ∧ l / y.B = pos / y.B
    · obtain ⟨rfl, hlb⟩ := hc
      rw [if_pos ⟨rfl, hlb⟩, Function.update_self]
      by_cases hlp : l = pos
      · subst hlp
        rw [if_pos ⟨rfl, rfl⟩, Function.update_self, hold]
      · have hmod : l % y.B ≠ pos % y.B := fun h => hlp (eq_of_block hlb h)
        rw [if_neg (by simp [hlp]), Function.update_of_ne hmod, ← hlb]
        exact habs g' hg' l hl
    · rw [if_neg hc]
      have hne : x.2 g' (l / y.B) ≠ nid := fun h => hnidNxt ⟨g', _, hg', hlNB, h⟩
      rw [Function.update_of_ne hne]
      have : ¬ (g' = g ∧ l = pos) := by
        rintro ⟨rfl, rfl⟩; exact hc ⟨rfl, rfl⟩
      rw [if_neg this]
      exact habs g' hg' l hl

/-- The scan of `CollectBlocks`, from block `n - 1` down to 0: the new high-water mark is
one above the highest marked block; unmarked blocks below it go to the free list. The
state is the free list (top first) and the new high-water mark (0: none found yet). -/
def scan (mark : ℕ → Bool) : ℕ → List ℕ × ℕ → List ℕ × ℕ
  | 0, st => st
  | n + 1, st =>
      scan mark n
        (if mark n then (st.1, if st.2 = 0 then n + 1 else st.2)
         else if st.2 = 0 then st else (n :: st.1, st.2))

theorem scan_spec (mark : ℕ → Bool) (HW : ℕ) (n : ℕ) (st : List ℕ × ℕ)
    (h1 : st.1.Nodup) (h2 : ∀ id ∈ st.1, n ≤ id ∧ id < st.2 ∧ mark id = false)
    (h3 : st.2 = 0 ∨ n < st.2)
    (h4 : ∀ id, n ≤ id → st.2 ≤ id → id < HW → mark id = false) :
    (scan mark n st).1.Nodup ∧ (∀ id ∈ (scan mark n st).1, id < (scan mark n st).2 ∧ mark id = false) ∧
      ∀ id, (scan mark n st).2 ≤ id → id < HW → mark id = false := by
  induction n generalizing st with
  | zero =>
    exact ⟨h1, fun id hid => ⟨(h2 id hid).2.1, (h2 id hid).2.2⟩, fun id h5 h6 => h4 id (Nat.zero_le _) h5 h6⟩
  | succ n ih =>
    simp only [scan]
    by_cases hm : mark n = true
    · rw [if_pos hm]
      apply ih
      · exact h1
      · intro id hid
        obtain ⟨a, b, c⟩ := h2 id hid
        refine ⟨by omega, ?_, c⟩
        show id < (if st.2 = 0 then n + 1 else st.2)
        split <;> omega
      · show (if st.2 = 0 then n + 1 else st.2) = 0 ∨ n < (if st.2 = 0 then n + 1 else st.2)
        right; split <;> omega
      · intro id h5 h6 h7
        have h6' : (if st.2 = 0 then n + 1 else st.2) ≤ id := h6
        apply h4 id _ _ h7
        · split at h6' <;> omega
        · split at h6' <;> omega
    · rw [if_neg hm]
      have hm' : mark n = false := by simpa using hm
      by_cases h0 : st.2 = 0
      · rw [if_pos h0]
        apply ih _ h1
        · intro id hid
          obtain ⟨a, b, c⟩ := h2 id hid
          omega
        · left; exact h0
        · intro id h5 h6 h7
          rcases Nat.eq_or_lt_of_le h5 with h | h
          · rw [← h]; exact hm'
          · exact h4 id h (by omega) h7
      · rw [if_neg h0]
        apply ih
        · show (n :: st.1).Nodup
          refine List.nodup_cons.mpr ⟨?_, h1⟩
          intro hn
          have := (h2 n hn).1
          omega
        · intro id hid
          rcases List.mem_cons.mp hid with rfl | hid
          · exact ⟨le_rfl, by show id < st.2; omega, hm'⟩
          · obtain ⟨a, b, c⟩ := h2 id hid
            exact ⟨by omega, b, c⟩
        · right; show n < st.2; omega
        · intro id h5 h6 h7
          have h6' : st.2 ≤ id := h6
          exact h4 id (by omega) h6' h7

/-- The marks `SumGenotypesW` sets while it reads the population: the blocks in use. -/
def markOf (R NB : ℕ) (t : ℕ → ℕ → ℕ) (id : ℕ) : Bool :=
  (List.range R).any (fun g => (List.range NB).any (fun k => decide (t g k = id)))

theorem markOf_iff (R NB : ℕ) (t : ℕ → ℕ → ℕ) (id : ℕ) :
    markOf R NB t id = true ↔ Refd R NB t id := by
  unfold markOf Refd
  simp only [List.any_eq_true, List.mem_range, decide_eq_true_eq]
  constructor
  · rintro ⟨g, hg, k, hk, h⟩; exact ⟨g, k, hg, hk, h⟩
  · rintro ⟨g, k, hg, hk, h⟩; exact ⟨g, hg, k, hk, h⟩

/-- `CollectBlocks`: a new free list and high-water mark from the marks. -/
def collect (R NB : ℕ) (st : Store F) (t : ℕ → ℕ → ℕ) : Store F :=
  ⟨st.data, (scan (markOf R NB t) st.hw ([], 0)).1, (scan (markOf R NB t) st.hw ([], 0)).2,
    st.owner⟩

/-- **Collecting unused blocks (target b, third part).** After a collection the allocator
is sound and hands out only blocks that no chromosome of the live population holds; the
contents of the blocks, and so every chromosome the population stands for, are untouched
(`collect` does not change `data`). -/
theorem collect_inv (R NB : ℕ) (st : Store F) (t : ℕ → ℕ → ℕ)
    (hav : ∀ id, Avail st id → ¬ Refd R NB t id) :
    AllocOK (collect R NB st t) ∧ ∀ id, Avail (collect R NB st t) id → ¬ Refd R NB t id := by
  obtain ⟨h1, h2, h3⟩ := scan_spec (markOf R NB t) st.hw st.hw ([], 0) List.nodup_nil
    (by simp) (Or.inl rfl) (by intro id h5 _ h7; omega)
  refine ⟨⟨h1, fun id hid => (h2 id hid).1⟩, ?_⟩
  rintro id (hid | hid) href
  · have := (h2 id hid).2
    rw [(markOf_iff R NB t id).mpr href] at this
    exact Bool.noConfusion this
  · have hid' : (scan (markOf R NB t) st.hw ([], 0)).2 ≤ id := hid
    rcases Nat.lt_or_ge id st.hw with h | h
    · have := h3 id hid' h
      rw [(markOf_iff R NB t id).mpr href] at this
      exact Bool.noConfusion this
    · exact hav id (Or.inr h) href

end Blocks
```

## 4. One generation, and n generations (target f)

**Informal claim (f).** Start the two programs from the same options and seed. Then after
every generation: the generator is in the same state in both; the environmental
deviations and the fitness parameters are equal; every chromosome of the fast program
stands for the original's chromosome; and the line of output is the same. If the
population dies, it dies in the same generation in both. This holds for every schedule of
block collections, every block size, and whatever the generator returns.

Everything the two programs compute in the same way from the genotypes (phenotypes, the
move of the optimum, fitness, survival, the output line) is left uninterpreted in
`Prims`: the model assumes only that both programs apply the *same* functions there. That
is target (e), checked by reading the code: where the original computes a mean or a
variance twice, the fast program calls the same GSL function once on the same array.

```lean
section Generation

variable {S F P Rec : Type}

/-- What the two programs share and the model leaves uninterpreted. Every field is
arbitrary; no law is assumed of any of them.

* `popFitness itter par geno env`: the part of `PopFitness` after the genotypes
  (phenotype = genotype + environment, the move of the optimum, the fitness of every
  individual): the new parameters and the fitnesses.
* `survives w`: `gsl_rng_uniform(rng) < w` in `SetSurviving`.
* `dead n`: `DeadPopulation`.
* `parents n`: the loop of `NewIn` that draws two different survivors out of `n`.
* `count`, `fresh`, `strand`: `gsl_ran_poisson(rng, RecombPoissonMean)`, the redrawn
  `gsl_rng_uniform_int(rng, NoQTL-1)`, and `gsl_rng_uniform(rng) < 0.5`.
* `envDraw`: `SetIndividualEnv`.
* `output`: the line of `Output`, from the generation number, the parameters, the
  genotypes, the environmental deviations, the fitnesses and the survivors.
* `change`: `if(Itter == ChangeFitnessGen) FitnessPar[1] = ChangeFitnessSD`.
* `mutCount`, `mutPos`, `mutEff`: the draws of `MutateIndividualChrNo`. -/
structure Prims (S F P Rec : Type) where
  popFitness : ℕ → P → List F → List F → P × (ℕ → F)
  survives : F → Rand S Bool
  dead : ℕ → Bool
  parents : ℕ → Rand S (ℕ × ℕ)
  count : Rand S ℕ
  fresh : (ℕ → Bool) → Rand S ℕ
  strand : Rand S Bool
  envDraw : Rand S F
  output : ℕ → P → List F → List F → (ℕ → F) → List ℕ → Rec
  change : ℕ → P → P
  mutCount : Rand S ℕ
  mutPos : Rand S ℕ
  mutEff : Rand S F

/-- `SetSurviving`: one draw per individual, in order; the survivors' indices. -/
def survLoop (surv : F → Rand S Bool) (fit : ℕ → F) : List ℕ → Rand S (List ℕ)
  | [] => fun s => ([], s)
  | i :: is => fun s =>
      ((if (surv (fit i) s).1 then [i] else []) ++ (survLoop surv fit is (surv (fit i) s).2).1,
        (survLoop surv fit is (surv (fit i) s).2).2)

theorem survLoop_mem (surv : F → Rand S Bool) (fit : ℕ → F) (is : List ℕ) (s : S) :
    ∀ j ∈ (survLoop surv fit is s).1, j ∈ is := by
  induction is generalizing s with
  | nil => simp [survLoop]
  | cons i is ih =>
    intro j hj
    simp only [survLoop, List.mem_append] at hj
    rcases hj with hj | hj
    · split at hj
      · simp at hj; simp [hj]
      · simp at hj
    · exact List.mem_cons_of_mem _ (ih _ j hj)

theorem getD_mem_or (l : List ℕ) (i d : ℕ) : l.getD i d = d ∨ l.getD i d ∈ l := by
  induction l generalizing i with
  | nil => left; simp
  | cons a l ih =>
    cases i with
    | zero => right; simp
    | succ i =>
      rcases ih i with h | h
      · left; simpa using h
      · right; simp only [List.getD_cons_succ]; exact List.mem_cons_of_mem _ h

/-- The draws for one offspring in the original, in its order: the parents, the first
gamete's map and strand, the second gamete's, the environmental deviation. -/
structure ODec (F : Type) where
  a : ℕ
  b : ℕ
  m1 : List ℕ
  s1 : Bool
  m2 : List ℕ
  s2 : Bool
  e : F

/-- The same in the fast program: lowest points in place of maps. -/
structure FDec (F : Type) where
  a : ℕ
  b : ℕ
  low1 : ℕ
  s1 : Bool
  low2 : ℕ
  s2 : Bool
  e : F

def ODec.toF (L : ℕ) (d : ODec F) : FDec F :=
  ⟨d.a, d.b, min L (d.m1.getD 0 0), d.s1, min L (d.m2.getD 0 0), d.s2, d.e⟩

variable (pr : Prims S F P Rec) (sort : List ℕ → List ℕ)

/-- `NewIn` + `RecombinationPoisson` + `SetIndividualEnv` of the original, draws only. -/
def origDecide (L : ℕ) (surv : List ℕ) : Rand S (ODec F) := fun s =>
  let r0 := pr.parents surv.length s
  let r1 := origMap pr.count pr.fresh sort L r0.2
  let r2 := pr.strand r1.2
  let r3 := origMap pr.count pr.fresh sort L r2.2
  let r4 := pr.strand r3.2
  let r5 := pr.envDraw r4.2
  (⟨surv.getD r0.1.1 0, surv.getD r0.1.2 0, r1.1, r2.1, r3.1, r4.1, r5.1⟩, r5.2)

/-- The draws for one offspring in the fast program's `NewPop`. -/
def fastDecide (L : ℕ) (surv : List ℕ) : Rand S (FDec F) := fun s =>
  let r0 := pr.parents surv.length s
  let r1 := fastLow pr.count pr.fresh L r0.2
  let r2 := pr.strand r1.2
  let r3 := fastLow pr.count pr.fresh L r2.2
  let r4 := pr.strand r3.2
  let r5 := pr.envDraw r4.2
  (⟨surv.getD r0.1.1 0, surv.getD r0.1.2 0, r1.1, r2.1, r3.1, r4.1, r5.1⟩, r5.2)

theorem fastDecide_eq {sort : List ℕ → List ℕ} (hs : SortSpec sort) (L : ℕ) (surv : List ℕ)
    (s : S) :
    fastDecide pr L surv s =
      ((origDecide pr sort L surv s).1.toF L, (origDecide pr sort L surv s).2) := by
  simp only [fastDecide, origDecide, fastLow_eq hs, ODec.toF]

/-- The original writes offspring `i`: its two gametes, by the copying loop, over what the
other population held, and its environmental deviation. -/
def origApply (L : ℕ) (cur : ℕ → ℕ → F) (i : ℕ) (d : ODec F) (x : (ℕ → ℕ → F) × List F) :
    (ℕ → ℕ → F) × List F :=
  (Function.update
      (Function.update x.1 (2 * i)
        (origLoop (cur (2 * d.a)) (cur (2 * d.a + 1)) d.m1 d.s1 L (x.1 (2 * i))))
      (2 * i + 1)
      (origLoop (cur (2 * d.b)) (cur (2 * d.b + 1)) d.m2 d.s2 L (x.1 (2 * i + 1))),
    x.2 ++ [d.e])

/-- The loop of the original's `NewPop`: draw, then copy, offspring by offspring. -/
def origNewPop (L N : ℕ) (cur : ℕ → ℕ → F) (surv : List ℕ) (old : ℕ → ℕ → F) :
    Rand S ((ℕ → ℕ → F) × List F) :=
  interleaved (fun _ _ => origDecide pr sort L surv) (origApply L cur) (List.range N) (old, [])

theorem applyAll_snoc {ι D G : Type} (apply : ι → D → G → G) (is : List ι) (ds : List D)
    (i : ι) (d : D) (g : G) (h : is.length = ds.length) :
    applyAll apply (is ++ [i]) (ds ++ [d]) g = apply i d (applyAll apply is ds g) := by
  induction is generalizing ds g with
  | nil =>
    cases ds with
    | nil => rfl
    | cons e ds => simp at h
  | cons j is ih =>
    cases ds with
    | nil => simp at h
    | cons e ds =>
      simp only [List.cons_append, applyAll]
      exact ih _ _ (by simpa using h)

/-- The gamete in closed form: the starting strand up to locus `m`, the other after. -/
def gamSpec (p1 p2 : ℕ → F) (m : ℕ) (s : Bool) (l : ℕ) : F :=
  if l ≤ m then (if s then p2 else p1) l else (if s then p1 else p2) l

theorem origNew_spec (L : ℕ) (cur old : ℕ → ℕ → F) (ds : List (ODec F)) :
    (applyAll (origApply L cur) (List.range ds.length) ds (old, [])).2 = ds.map (·.e) ∧
      ∀ i d, ds[i]? = some d → ∀ l < L,
        (applyAll (origApply L cur) (List.range ds.length) ds (old, [])).1 (2 * i) l =
            gamSpec (cur (2 * d.a)) (cur (2 * d.a + 1)) (d.m1.getD 0 0) d.s1 l ∧
          (applyAll (origApply L cur) (List.range ds.length) ds (old, [])).1 (2 * i + 1) l =
            gamSpec (cur (2 * d.b)) (cur (2 * d.b + 1)) (d.m2.getD 0 0) d.s2 l := by
  induction ds using List.reverseRecOn with
  | nil => exact ⟨rfl, by simp⟩
  | append_singleton ds d ih =>
    obtain ⟨ih1, ih2⟩ := ih
    have hlen : (ds ++ [d]).length = ds.length + 1 := by simp
    rw [hlen, List.range_succ, applyAll_snoc _ _ _ _ _ _ (by simp)]
    refine ⟨by simp [origApply, ih1], ?_⟩
    intro i d' hd' l hl
    simp only [origApply]
    by_cases hi : i < ds.length
    · rw [List.getElem?_append_left hi] at hd'
      rw [Function.update_of_ne (by omega), Function.update_of_ne (by omega),
        Function.update_of_ne (by omega), Function.update_of_ne (by omega)]
      exact ih2 i d' hd' l hl
    · have hi' : i = ds.length := by
        by_contra hne
        have : (ds ++ [d]).length ≤ i := by rw [hlen]; omega
        rw [List.getElem?_eq_none this] at hd'
        exact absurd hd' (by simp)
      subst hi'
      have : d' = d := by simpa using hd'.symm
      subst this
      rw [Function.update_of_ne (by omega), Function.update_self, Function.update_self]
      simp only [origLoop_eq, if_pos hl, gamSpec]
      simp

/-- The gamete records of the fast program, by chromosome of the next population. -/
def descOf (y : Layout) (decs : List (FDec F)) (g : ℕ) : GDesc :=
  match decs[g / 2]? with
  | none => ⟨0, 0, 0, none⟩
  | some d => if g % 2 = 0 then gdesc y d.a d.low1 d.s1 else gdesc y d.b d.low2 d.s2

theorem drawAll_forall {ι D : Type} (draw : ι → Rand S D) (Pd : D → Prop)
    (h : ∀ i s, Pd (draw i s).1) (items : List ι) (s : S) :
    ∀ d ∈ (drawAll draw items s).1, Pd d := by
  induction items generalizing s with
  | nil => simp [drawAll]
  | cons i is ih =>
    intro d hd
    simp only [drawAll, List.mem_cons] at hd
    rcases hd with rfl | hd
    · exact h i s
    · exact ih _ d hd

/-- **`NewPop` (targets b, c and d together).** The fast program's two passes (all the
draws, then `FillGametes` and `MixBlocks` over shared blocks) leave the generator where
the original's interleaved loop leaves it, draw the same environmental deviations, and
build a table of block numbers that stands for the population the original builds. -/
theorem newPop_refines {sort : List ℕ → List ℕ} (hs : SortSpec sort) (y : Layout) {N : ℕ}
    (hN : 0 < N) (cur : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (st : Store F) (surv : List ℕ)
    (hsurv : ∀ j ∈ surv, j < N)
    (habs : ∀ g < 2 * N, ∀ l < y.L, absChr y.B st.data (ids g) l = cur g l)
    (hok : AllocOK st) (hav : ∀ id, Avail st id → ¬ Refd (2 * N) y.NB ids id)
    (hown : ∀ id, st.owner id = 0) (old : ℕ → ℕ → F) (s : S) :
    (drawAll (fun _ => fastDecide pr y.L surv) (List.range N) s).2 =
        (origNewPop pr sort y.L N cur surv old s).2 ∧
      (drawAll (fun _ => fastDecide pr y.L surv) (List.range N) s).1.map (·.e) =
        (origNewPop pr sort y.L N cur surv old s).1.2 ∧
      MutInv y (2 * N) (fillMix (2 * N) ids
        (descOf y (drawAll (fun _ => fastDecide pr y.L surv) (List.range N) s).1) st) ∧
      ∀ g < 2 * N, ∀ l < y.L,
        absChr y.B (fillMix (2 * N) ids
            (descOf y (drawAll (fun _ => fastDecide pr y.L surv) (List.range N) s).1) st).1.data
          ((fillMix (2 * N) ids
            (descOf y (drawAll (fun _ => fastDecide pr y.L surv) (List.range N) s).1) st).2 g) l =
          (origNewPop pr sort y.L N cur surv old s).1.1 g l := by
  -- the original as two passes
  have htwo : origNewPop pr sort y.L N cur surv old s =
      (applyAll (origApply y.L cur) (List.range N)
          (drawAll (fun _ => origDecide pr sort y.L surv) (List.range N) s).1 (old, []),
        (drawAll (fun _ => origDecide pr sort y.L surv) (List.range N) s).2) :=
    interleaved_eq_two_pass (draw := fun _ _ => origDecide pr sort y.L surv)
      (fun _ _ _ _ => rfl) (origApply y.L cur) (old, []) (List.range N) (old, []) s
  -- the fast program's draws are the original's, with maps replaced by lowest points
  have hfd : (fun (_ : ℕ) => fastDecide pr y.L surv) =
      fun _ s' => (ODec.toF y.L (origDecide pr sort y.L surv s').1,
        (origDecide pr sort y.L surv s').2) := by
    funext _ s'; exact fastDecide_eq pr hs y.L surv s'
  have hdr := drawAll_map (ODec.toF y.L) (fun (_ : ℕ) => origDecide pr sort y.L surv)
    (List.range N) s
  rw [hfd, hdr, htwo]
  set ds := (drawAll (fun _ => origDecide pr sort y.L surv) (List.range N) s).1 with hds
  have hlen : ds.length = N := by rw [hds, drawAll_length]; simp
  have hpar : ∀ d ∈ ds, d.a < N ∧ d.b < N := by
    apply drawAll_forall (fun (_ : ℕ) => origDecide pr sort y.L surv)
      (fun d => d.a < N ∧ d.b < N)
    intro _ s'
    simp only [origDecide]
    constructor
    · rcases getD_mem_or surv (pr.parents surv.length s').1.1 0 with h | h
      · rw [h]; exact hN
      · exact hsurv _ h
    · rcases getD_mem_or surv (pr.parents surv.length s').1.2 0 with h | h
      · rw [h]; exact hN
      · exact hsurv _ h
  obtain ⟨hsp1, hsp2⟩ := origNew_spec y.L cur old ds
  rw [hlen] at hsp1 hsp2
  refine ⟨rfl, ?_, ?_⟩
  · rw [hsp1, List.map_map]; rfl
  -- the parent, lowest point and strand of every gamete
  let aF : ℕ → ℕ := fun g => match ds[g / 2]? with
    | none => 0
    | some d => if g % 2 = 0 then d.a else d.b
  let lowF : ℕ → ℕ := fun g => match ds[g / 2]? with
    | none => 0
    | some d => if g % 2 = 0 then min y.L (d.m1.getD 0 0) else min y.L (d.m2.getD 0 0)
  let sF : ℕ → Bool := fun g => match ds[g / 2]? with
    | none => false
    | some d => if g % 2 = 0 then d.s1 else d.s2
  have hget : ∀ g < 2 * N, ∃ d, ds[g / 2]? = some d ∧ d ∈ ds := by
    intro g hg
    have h2 : g / 2 < ds.length := by rw [hlen]; omega
    exact ⟨ds[g / 2], List.getElem?_eq_getElem h2, List.getElem_mem h2⟩
  have hfm := fillMix_spec y (2 * N) ids st aF lowF sF (descOf y (ds.map (ODec.toF y.L)))
    (by
      intro g hg
      obtain ⟨d, hd, _⟩ := hget g hg
      simp only [descOf, aF, lowF, sF, List.getElem?_map, hd, Option.map_some, ODec.toF]
      split <;> rfl)
    (by
      intro g hg
      obtain ⟨d, hd, hmem⟩ := hget g hg
      have := hpar d hmem
      simp only [aF, hd]
      split <;> omega)
    (by
      intro g hg
      obtain ⟨d, hd, _⟩ := hget g hg
      simp only [lowF, hd]
      split <;> exact Nat.min_le_left _ _)
    hok hav hown
  refine ⟨hfm.1, ?_⟩
  intro g hg l hl
  rw [hfm.2 g hg l hl]
  obtain ⟨d, hd, hmem⟩ := hget g hg
  obtain ⟨hda, hdb⟩ := hpar d hmem
  obtain ⟨hrow0, hrow1⟩ := hsp2 (g / 2) d hd l hl
  have hmin : ∀ m, (l ≤ min y.L m) = (l ≤ m) := by
    intro m; apply propext; constructor
    · intro h; exact le_trans h (Nat.min_le_right _ _)
    · intro h; exact Nat.le_min.mpr ⟨le_of_lt hl, h⟩
  simp only [aF, lowF, sF, hd]
  by_cases hpar2 : g % 2 = 0
  · have hg2 : 2 * (g / 2) = g := by omega
    rw [hg2] at hrow0
    rw [hrow0]
    simp only [hpar2, if_true, hmin, gamSpec]
    by_cases hle : l ≤ d.m1.getD 0 0
    · rw [if_pos hle, if_pos hle]
      cases d.s1
      · exact habs _ (by simp; omega) l hl
      · exact habs _ (by simp; omega) l hl
    · rw [if_neg hle, if_neg hle]
      cases d.s1
      · exact habs _ (by simp; omega) l hl
      · exact habs _ (by simp; omega) l hl
  · have hg2 : 2 * (g / 2) + 1 = g := by omega
    rw [hg2] at hrow1
    rw [hrow1]
    simp only [hpar2, if_false, hmin, gamSpec]
    by_cases hle : l ≤ d.m2.getD 0 0
    · rw [if_pos hle, if_pos hle]
      cases d.s2
      · exact habs _ (by simp; omega) l hl
      · exact habs _ (by simp; omega) l hl
    · rw [if_neg hle, if_neg hle]
      cases d.s2
      · exact habs _ (by simp; omega) l hl
      · exact habs _ (by simp; omega) l hl

/-- `MutatePop`, for either representation of the genomes: for each chromosome in order,
the number of mutations, then for each mutation its locus, its effect, and the change,
made at once (both programs interleave these). -/
def mutLoop {G : Type} (apply : ℕ → ℕ → F → G → G) (R : ℕ) : G → Rand S G :=
  loopM (fun g x s =>
    loopM (fun (_ : Unit) x s =>
        (apply g (pr.mutPos s).1 (pr.mutEff (pr.mutPos s).2).1 x, (pr.mutEff (pr.mutPos s).2).2))
      (List.replicate (pr.mutCount s).1 ()) x (pr.mutCount s).2) (List.range R)

/-- **`MutatePop` over shared blocks.** `hpos` is the one fact about a draw that the proof
uses: a mutation's locus is below the number of loci (`gsl_rng_uniform_int(rng, NoQTL)`
returns a value below its bound by construction, whatever the generator). -/
theorem mutLoop_sim (add : F → F → F) (y : Layout) (hpos : ∀ s, (pr.mutPos s).1 < y.L) (R : ℕ)
    (x : Store F × (ℕ → ℕ → ℕ)) (pop : ℕ → ℕ → F) (s : S) (hinv : MutInv y R x)
    (habs : ∀ g < R, ∀ l < y.L, absChr y.B x.1.data (x.2 g) l = pop g l) :
    (MutInv y R (mutLoop pr (mutStep add y.B) R x s).1 ∧
      ∀ g < R, ∀ l < y.L,
        absChr y.B (mutLoop pr (mutStep add y.B) R x s).1.1.data
            ((mutLoop pr (mutStep add y.B) R x s).1.2 g) l =
          (mutLoop pr (mutDense add) R pop s).1 g l) ∧
      (mutLoop pr (mutStep add y.B) R x s).2 = (mutLoop pr (mutDense add) R pop s).2 := by
  unfold mutLoop
  apply loopM_sim (fun (x : Store F × (ℕ → ℕ → ℕ)) (pop : ℕ → ℕ → F) =>
    MutInv y R x ∧ ∀ g < R, ∀ l < y.L, absChr y.B x.1.data (x.2 g) l = pop g l)
  · intro g hg x pop s hx
    have hgR : g < R := List.mem_range.mp hg
    apply loopM_sim (fun (x : Store F × (ℕ → ℕ → ℕ)) (pop : ℕ → ℕ → F) =>
      MutInv y R x ∧ ∀ g < R, ∀ l < y.L, absChr y.B x.1.data (x.2 g) l = pop g l)
    · intro _ _ x pop s hx
      exact ⟨mutStep_sim add y R hgR (hpos s) _ x pop hx.1 hx.2, rfl⟩
    · exact hx
  · exact ⟨hinv, habs⟩

/-- The original's state between generations: the generator, the current population
(`CPop`: chromosome, locus ↦ allele), the other buffer (`NPop`), the environmental
deviations of the current population, the parameters (`FitnessPar` and the options). -/
structure OState (S F P : Type) where
  rng : S
  pop : ℕ → ℕ → F
  other : ℕ → ℕ → F
  env : List F
  par : P

/-- The fast program's state: the generator, the block store, the current population's
block numbers (`S->Ids`), the environmental deviations, the parameters. -/
structure FState (S F P : Type) where
  rng : S
  st : Store F
  ids : ℕ → ℕ → ℕ
  env : List F
  par : P

variable (add : F → F → F) (zero : F) (y : Layout) (N : ℕ)

/-- **One generation of the original** (`RunSim`'s loop body, diploid): `PopFitness`,
`NewPop`, `Output`, the change of the fitness width, `MutatePop`, swap. `none`: the
population died (`Pop all died.`). -/
def origGen (itter : ℕ) (o : OState S F P) : Option (Rec × OState S F P) :=
  let geno := (List.range N).map (origGeno add zero y.L o.pop)
  let pf := pr.popFitness itter o.par geno o.env
  let sv := survLoop pr.survives pf.2 (List.range N) o.rng
  if pr.dead sv.1.length then none
  else
    let np := origNewPop pr sort y.L N o.pop sv.1 o.other sv.2
    let mu := mutLoop pr (mutDense add) (2 * N) np.1.1 np.2
    some (pr.output itter pf.1 geno o.env pf.2 sv.1,
      ⟨mu.2, mu.1, o.pop, np.1.2, pr.change itter pf.1⟩)

/-- The tags are cleared at the end of `MutatePop`. -/
def Store.clearOwner (st : Store F) : Store F := ⟨st.data, st.free, st.hw, fun _ => 0⟩

/-- **One generation of the fast program**: the genotypes position by position over the
blocks; a collection of unused blocks in the generations `gc` names (any schedule); then
as the original, with `NewPop` in two passes and `MutatePop` copying on write. -/
def fastGen (gc : ℕ → Bool) (itter : ℕ) (f : FState S F P) : Option (Rec × FState S F P) :=
  let geno := (List.range N).map (fastGenoAll add zero y f.st.data f.ids)
  let st0 := if gc itter then collect (2 * N) y.NB f.st f.ids else f.st
  let pf := pr.popFitness itter f.par geno f.env
  let sv := survLoop pr.survives pf.2 (List.range N) f.rng
  if pr.dead sv.1.length then none
  else
    let dr := drawAll (fun _ => fastDecide pr y.L sv.1) (List.range N) sv.2
    let x := fillMix (2 * N) f.ids (descOf y dr.1) st0
    let mu := mutLoop pr (mutStep add y.B) (2 * N) x dr.2
    some (pr.output itter pf.1 geno f.env pf.2 sv.1,
      ⟨mu.2, mu.1.1.clearOwner, mu.1.2, dr.1.map (·.e), pr.change itter pf.1⟩)

/-- The fast state stands for the original state: same generator state, deviations and
parameters; every chromosome's row of block numbers stands for the original's dense
chromosome; and the store is sound (no block that can be handed out is in use, no tags). -/
structure Rel (f : FState S F P) (o : OState S F P) : Prop where
  rng : f.rng = o.rng
  env : f.env = o.env
  par : f.par = o.par
  abs : ∀ g < 2 * N, ∀ l < y.L, absChr y.B f.st.data (f.ids g) l = o.pop g l
  ok : AllocOK f.st
  avail : ∀ id, Avail f.st id → ¬ Refd (2 * N) y.NB f.ids id
  owner : ∀ id, f.st.owner id = 0

/-- Two outcomes of a generation agree: both died, or both printed the same line and
their new states are related. -/
def RelOut {A B : Type} (Rl : A → B → Prop) : Option (Rec × A) → Option (Rec × B) → Prop
  | none, none => True
  | some p, some q => p.1 = q.1 ∧ Rl p.2 q.2
  | _, _ => False

theorem geno_eq (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (pop : ℕ → ℕ → F)
    (habs : ∀ g < 2 * N, ∀ l < y.L, absChr y.B data (ids g) l = pop g l) :
    (List.range N).map (fastGenoAll add zero y data ids) =
      (List.range N).map (origGeno add zero y.L pop) := by
  apply List.map_congr_left
  intro i hi
  have hiN : i < N := List.mem_range.mp hi
  rw [fastGeno_eq]
  unfold origGeno
  apply foldl_ext_mem
  intro a l hl
  have hlL : l < y.L := List.mem_range.mp hl
  simp only []
  rw [habs (2 * i) (by omega) l hlL, habs (2 * i + 1) (by omega) l hlL]

/-- **One generation (target f).** From related states, one generation of the fast
program and one of the original either both end with the population dead, or print the
same output record and reach related states: the same generator state, the same
environmental deviations and parameters, and block tables that stand for the original's
population. Hypotheses: `qsort` sorts; a mutation's locus is below the number of loci;
the population is not empty. -/
theorem gen_refines {sort : List ℕ → List ℕ} (hs : SortSpec sort)
    (hpos : ∀ s, (pr.mutPos s).1 < y.L) (hN : 0 < N) (gc : ℕ → Bool) (itter : ℕ)
    (f : FState S F P) (o : OState S F P) (h : Rel y N f o) :
    RelOut (Rel y N) (fastGen pr add zero y N gc itter f) (origGen pr sort add zero y N itter o) := by
  have hgeno := geno_eq add zero y N f.st.data f.ids o.pop h.abs
  -- the store after the collection, if any
  obtain ⟨st0, hst0, hdata0, hok0, hav0, hown0⟩ :
      ∃ st0 : Store F, st0 = (if gc itter then collect (2 * N) y.NB f.st f.ids else f.st) ∧
        st0.data = f.st.data ∧ AllocOK st0 ∧
        (∀ id, Avail st0 id → ¬ Refd (2 * N) y.NB f.ids id) ∧ ∀ id, st0.owner id = 0 := by
    refine ⟨_, rfl, ?_⟩
    by_cases hgc : gc itter = true
    · rw [if_pos hgc]
      obtain ⟨c1, c2⟩ := collect_inv (2 * N) y.NB f.st f.ids h.avail
      exact ⟨rfl, c1, c2, h.owner⟩
    · rw [if_neg hgc]
      exact ⟨rfl, h.ok, h.avail, h.owner⟩
  simp only [fastGen, origGen, hgeno, h.env, h.par, h.rng, ← hst0]
  set pf := pr.popFitness itter o.par (List.map (origGeno add zero y.L o.pop) (List.range N)) o.env
    with hpf
  set sv := survLoop pr.survives pf.2 (List.range N) o.rng with hsv
  by_cases hd : pr.dead sv.1.length = true
  · simp only [hd, if_true]; trivial
  · simp only [hd]
    have hsurv : ∀ j ∈ sv.1, j < N := fun j hj =>
      List.mem_range.mp (survLoop_mem pr.survives pf.2 (List.range N) o.rng j hj)
    obtain ⟨n1, n2, n3, n4⟩ := newPop_refines pr hs y hN o.pop f.ids st0 sv.1 hsurv
      (by intro g hg l hl; rw [hdata0]; exact h.abs g hg l hl) hok0 hav0 hown0 o.other sv.2
    obtain ⟨⟨m1, m2⟩, m3⟩ := mutLoop_sim pr add y hpos (2 * N) _ _
      (drawAll (fun _ => fastDecide pr y.L sv.1) (List.range N) sv.2).2 n3 n4
    refine ⟨rfl, ?_⟩
    refine ⟨?_, n2, rfl, ?_, ⟨m1.ok.nodup, m1.ok.lt⟩, m1.avail, fun _ => rfl⟩
    · show (mutLoop pr (mutStep add y.B) (2 * N) _ _).2 = (mutLoop pr (mutDense add) (2 * N) _ _).2
      rw [m3, n1]
    · intro g hg l hl
      show absChr y.B (mutLoop pr (mutStep add y.B) (2 * N) _ _).1.1.data _ l =
        (mutLoop pr (mutDense add) (2 * N) _ _).1 g l
      rw [m2 g hg l hl, n1]

/-- A run of `n` generations from generation `t`: the output lines, and the final state
(`none` if the population died). -/
def runGen {St : Type} (gen : ℕ → St → Option (Rec × St)) : ℕ → ℕ → St → List Rec × Option St
  | 0, _, s => ([], some s)
  | n + 1, t, s =>
    match gen t s with
    | none => ([], none)
    | some p => (p.1 :: (runGen gen n (t + 1) p.2).1, (runGen gen n (t + 1) p.2).2)

/-- Two final states agree: both runs ended with the population dead, or the states are
related. -/
def RelFin {A B : Type} (Rl : A → B → Prop) : Option A → Option B → Prop
  | none, none => True
  | some a, some b => Rl a b
  | _, _ => False

theorem run_sim {A B : Type} (Rl : A → B → Prop) (genA : ℕ → A → Option (Rec × A))
    (genB : ℕ → B → Option (Rec × B))
    (hgen : ∀ t a b, Rl a b → RelOut Rl (genA t a) (genB t b)) (n t : ℕ) (a : A) (b : B)
    (h : Rl a b) :
    (runGen genA n t a).1 = (runGen genB n t b).1 ∧
      RelFin Rl (runGen genA n t a).2 (runGen genB n t b).2 := by
  induction n generalizing t a b with
  | zero => exact ⟨rfl, h⟩
  | succ n ih =>
    have hg := hgen t a b h
    simp only [runGen]
    cases hA : genA t a with
    | none =>
      cases hB : genB t b with
      | none => exact ⟨rfl, trivial⟩
      | some q => rw [hA, hB] at hg; exact hg.elim
    | some p =>
      cases hB : genB t b with
      | none => rw [hA, hB] at hg; exact hg.elim
      | some q =>
        rw [hA, hB] at hg
        obtain ⟨hg1, hg2⟩ := hg
        obtain ⟨i1, i2⟩ := ih (t + 1) p.2 q.2 hg2
        exact ⟨by simp only [hg1, i1], i2⟩

/-- **n generations (target f).** From related states, `n` generations of the fast program
print the same output lines as `n` generations of the original, stop at the same
generation if the population dies, and otherwise end in related states (in particular
with the generator in the same state), for every schedule of block collections. -/
theorem run_refines {sort : List ℕ → List ℕ} (hs : SortSpec sort)
    (hpos : ∀ s, (pr.mutPos s).1 < y.L) (hN : 0 < N) (gc : ℕ → Bool) (n t : ℕ)
    (f : FState S F P) (o : OState S F P) (h : Rel y N f o) :
    (runGen (fastGen pr add zero y N gc) n t f).1 = (runGen (origGen pr sort add zero y N) n t o).1 ∧
      RelFin (Rel y N) (runGen (fastGen pr add zero y N gc) n t f).2
        (runGen (origGen pr sort add zero y N) n t o).2 :=
  run_sim (Rel y N) _ _ (fun t f o h => gen_refines pr add zero y N hs hpos hN gc t f o h) n t f o h

/-- The states the two programs start from (`CreatePop` twice; `CreateSim`): all alleles
`zero`; in the fast program every cell holds block 0, the only block in use. They are
related, given the same generator state, deviations and parameters. -/
theorem init_rel (rng : S) (env : List F) (par : P) (other : ℕ → ℕ → F) :
    Rel y N (⟨rng, ⟨fun _ _ => zero, [], 1, fun _ => 0⟩, fun _ _ => 0, env, par⟩ : FState S F P)
      ⟨rng, fun _ _ => zero, other, env, par⟩ := by
  refine ⟨rfl, rfl, rfl, fun _ _ _ _ => rfl, ⟨List.nodup_nil, by simp⟩, ?_, fun _ => rfl⟩
  rintro id (hid | hid) ⟨g, k, _, _, hcell⟩
  · simp at hid
  · have : (1 : ℕ) ≤ id := hid
    have : (0 : ℕ) = id := hcell
    omega

end Generation

end QTSimFast
```
