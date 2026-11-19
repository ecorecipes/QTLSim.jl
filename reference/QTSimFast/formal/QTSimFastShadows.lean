import QTSimFast

/-!
# QTSimFast: do the theorems say what is claimed?

A theorem can be proved and still not say what its author meant: a conclusion may have
been weakened, a hypothesis added that is never satisfied, a definition may not describe
the thing it is named after. This file checks the main theorems of `QTSimFast.lean`
against their informal claims (written in that file's section introductions and in
`PROOFS.md`) in the manner of the SA-Pass idea of Han et al., *ShadowBench: Toward
Reliable Automatic Evaluation of Semantic Alignment in Autoformalization*
(arXiv:2608.29270): surround an intended statement $T$ with auxiliary "shadow" statements
and machine-check the implications between them.

**This is an adaptation, not the published procedure.** SA-Pass is an evaluation metric
that presupposes a reference formalization and a shadow set built from it. There is no
reference formalization of QTSim, so here the shadows are derived from the informal
claims, and the checks are:

* **forward shadows** $T \Rightarrow S_i$: consequences that the informal claim obviously
  has (special cases, structural facts) are derived from the Lean theorem. A failure
  would reveal a conclusion that is too weak or a hypothesis that is too strong;
* **second formalization** $T \Leftrightarrow T'$: the same informal claim is formalized
  again in a different style (recursion instead of folds, and so on) and the two are
  proved equivalent without using the truth of either. This stands in for SA-Pass's
  backward check;
* **non-vacuity**: the hypotheses are satisfied by an explicit instance, on which the
  conclusion is exhibited, so no theorem holds only because a hypothesis is contradictory
  or a set empty;
* **negative shadows**: for algorithms that are *not* exact (sums in another order, cached
  block sums, switching at the highest crossover point, writing a shared block in place,
  draws in another order) the corresponding statement is refuted on an explicit example,
  with an addition that is neither associative nor commutative. The theorems are strong
  enough to tell these apart.

Even with these checks, the theorems concern a model. That the model describes the C
programs is supported by inspection, by `formal/model_check.c` and by the parity tests;
it is not proved.

Examples closed by `decide +kernel` are evaluated and checked by Lean's kernel (no
compiler, no additional axiom). Lines starting `#guard` are
evaluated by Lean's compiler and are tests, not proofs; they are used where an example is
too large for the kernel or involves hardware floating point.
-/

namespace QTSimFast.Shadows

open QTSimFast

/-!
## The test arithmetic

`nadd a b = 2a + b` on the integers is neither associative nor commutative, so any
regrouping or reordering of additions is visible. The second line shows the same for the
processor's doubles.
-/

/-- An "addition" with no algebraic law. -/
def nadd (a b : ℤ) : ℤ := 2 * a + b

example : nadd (nadd 1 2) 3 ≠ nadd 1 (nadd 2 3) := by decide +kernel
example : nadd 1 2 ≠ nadd 2 1 := by decide +kernel

#guard ((0.1 : Float) + 0.2) + 0.3 != 0.1 + (0.2 + 0.3)

/-- The layout used in the examples: 6 loci in 3 blocks of 2. -/
def lay : Layout := ⟨2, 2, 2, by decide, by decide, by decide⟩

example : lay.L = 6 ∧ lay.NB = 3 := by decide +kernel

/-- Block numbers with sharing: chromosome `g` holds block `(g + k) % 4` at position `k`. -/
def cIds (g k : ℕ) : ℕ := (g + k) % 4

/-- Block contents. -/
def cData (id o : ℕ) : ℤ := 3 * id + o + 1

/-!
## (a) Genotype sums

**Second formalization.** `sumRec` writes the original's sum as a recursion on the number
of loci, with no list and no fold: $G_0 = 0$, $G_{n+1} = G_n \oplus (a_n \oplus b_n)$.
-/

section Sums

variable {F : Type} (add : F → F → F) (zero : F)

/-- The original's genotype as a recursion on the number of loci. -/
def sumRec (a b : ℕ → F) : ℕ → F
  | 0 => zero
  | n + 1 => add (sumRec a b n) (add (a n) (b n))

/-- The fold of the main file and the recursion are the same function. -/
theorem origGeno_eq_sumRec (L : ℕ) (pop : ℕ → ℕ → F) (i : ℕ) :
    origGeno add zero L pop i = sumRec add zero (pop (2 * i)) (pop (2 * i + 1)) L := by
  unfold origGeno
  induction L with
  | zero => rfl
  | succ n ih => rw [List.range_succ, List.foldl_append, ih]; rfl

theorem sumRec_congr {a a' b b' : ℕ → F} (n : ℕ) (ha : ∀ l < n, a l = a' l)
    (hb : ∀ l < n, b l = b' l) : sumRec add zero a b n = sumRec add zero a' b' n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [sumRec]
    rw [ih (fun l hl => ha l (by omega)) (fun l hl => hb l (by omega)), ha n (by omega),
      hb n (by omega)]

/-- $T'$ for target (a): the fast sums equal the recursion. -/
def SumsT' : Prop :=
  ∀ (y : Layout) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (i : ℕ),
    fastGenoAll add zero y data ids i =
      sumRec add zero (absChr y.B data (ids (2 * i))) (absChr y.B data (ids (2 * i + 1))) y.L

/-- $T$ for target (a), as stated in the main file. -/
def SumsT : Prop :=
  ∀ (y : Layout) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (i : ℕ),
    fastGenoAll add zero y data ids i =
      origGeno add zero y.L (fun g => absChr y.B data (ids g)) i

/-- **Second-formalization check (a):** $T \Leftrightarrow T'$. The proof uses only that
the fold and the recursion agree, not that either statement is true. -/
theorem sums_T_iff_T' : SumsT add zero ↔ SumsT' add zero := by
  constructor
  · intro h y data ids i
    rw [h, origGeno_eq_sumRec]
  · intro h y data ids i
    rw [h, origGeno_eq_sumRec]

/-- Both hold. -/
theorem sums_T : SumsT add zero := fun y data ids i => fastGeno_eq add zero y data ids i

theorem sums_T' : SumsT' add zero := (sums_T_iff_T' add zero).mp (sums_T add zero)

/-- A layout with one block for the whole chromosome (`QTSIM_BLOCK` ≥ `NoQTL`). -/
def oneBlock (L : ℕ) (hL : 0 < L) : Layout := ⟨L, 0, L, hL, hL, le_rfl⟩

/-- A layout with one locus per block (block size 1: the dense layout). -/
def unitBlocks (M : ℕ) : Layout := ⟨1, M, 1, by decide, by decide, le_rfl⟩

/-- **Forward shadow (a1), one block:** the sum runs straight through the individual's two
blocks. -/
theorem shadow_one_block (L : ℕ) (hL : 0 < L) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (i : ℕ) :
    fastGenoAll add zero (oneBlock L hL) data ids i =
      sumRec add zero (data (ids (2 * i) 0)) (data (ids (2 * i + 1) 0)) L := by
  rw [sums_T' add zero]
  show sumRec add zero _ _ (0 * L + L) = _
  rw [Nat.zero_mul, Nat.zero_add]
  apply sumRec_congr <;> intro l hl <;>
    simp [absChr, oneBlock, Nat.div_eq_of_lt hl, Nat.mod_eq_of_lt hl]

/-- **Forward shadow (a2), block size 1:** locus `l` is the only allele of the block at
position `l`. -/
theorem shadow_unit_blocks (M : ℕ) (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (i : ℕ) :
    fastGenoAll add zero (unitBlocks M) data ids i =
      sumRec add zero (fun l => data (ids (2 * i) l) 0) (fun l => data (ids (2 * i + 1) l) 0)
        (M + 1) := by
  rw [sums_T' add zero]
  show sumRec add zero _ _ (M * 1 + 1) = _
  rw [Nat.mul_one]
  apply sumRec_congr <;> intro l _ <;> simp [absChr, unitBlocks, Nat.mod_one]

/-- **Forward shadow (a3), one locus:** $G = 0 \oplus (a_0 \oplus b_0)$. -/
theorem shadow_one_locus (data : ℕ → ℕ → F) (ids : ℕ → ℕ → ℕ) (i : ℕ) :
    fastGenoAll add zero (oneBlock 1 (by decide)) data ids i =
      add zero (add (data (ids (2 * i) 0) 0) (data (ids (2 * i + 1) 0) 0)) := by
  rw [shadow_one_block]; rfl

/-- **Forward shadow (a4), locality:** an individual's genotype depends only on its own two
chromosomes: two stores and tables that stand for the same two chromosomes of individual
`i` give it the same genotype, whatever the other individuals hold. (This is what allows
the individuals to be divided among lanes and threads.) -/
theorem shadow_local (y : Layout) (data data' : ℕ → ℕ → F) (ids ids' : ℕ → ℕ → ℕ) (i : ℕ)
    (h0 : ∀ l < y.L, absChr y.B data (ids (2 * i)) l = absChr y.B data' (ids' (2 * i)) l)
    (h1 : ∀ l < y.L,
      absChr y.B data (ids (2 * i + 1)) l = absChr y.B data' (ids' (2 * i + 1)) l) :
    fastGenoAll add zero y data ids i = fastGenoAll add zero y data' ids' i := by
  rw [sums_T' add zero, sums_T' add zero]
  exact sumRec_congr add zero _ h0 h1

/-- **Forward shadow (a5), the block size does not matter:** two layouts of the same
number of loci, over stores that stand for the same chromosomes, give the same genotypes. -/
theorem shadow_layout_free (y y' : Layout) (hL : y.L = y'.L) (data data' : ℕ → ℕ → F)
    (ids ids' : ℕ → ℕ → ℕ) (i : ℕ)
    (h : ∀ g, ∀ l < y.L, absChr y.B data (ids g) l = absChr y'.B data' (ids' g) l) :
    fastGenoAll add zero y data ids i = fastGenoAll add zero y' data' ids' i := by
  rw [sums_T' add zero, sums_T' add zero, ← hL]
  exact sumRec_congr add zero _ (h _) (h _)

end Sums

/-!
**Numeric instances (a)**, with the lawless `nadd`, 6 loci in 3 blocks of 2 and shared
blocks. Individual 1 holds chromosomes 2 and 3, that is blocks 2, 3, 0 and 3, 0, 1, with
alleles 7, 8, 10, 11, 1, 2 and 10, 11, 1, 2, 4, 5. By hand:
$(((((0 \oplus 17) \oplus 19) \oplus 11) \oplus 13) \oplus 5) \oplus 7$ with
$x \oplus y = 2x + y$ is $17, 53, 117, 247, 499, 1485$.
-/

example : fastGenoAll nadd 0 lay cData cIds 1 = 1485 := by decide +kernel
example : origGeno nadd 0 6 (fun g => absChr 2 cData (cIds g)) 1 = 1485 := by decide +kernel
example : sumRec nadd 0 (absChr 2 cData (cIds 2)) (absChr 2 cData (cIds 3)) 6 = 1485 := by
  decide +kernel

/-!
**Negative shadows (a).** Two sums that a faster program might compute, each different
from the original's on this example: the loci in reverse order (738), and one cached
partial sum per block, the block sums then added (453; this is what QTLSim.jl's kernel
does, and why it is not bit-identical).
-/

/-- The loci in reverse order. -/
def revGeno (pop : ℕ → ℕ → ℤ) (i : ℕ) : ℤ :=
  (List.range 6).reverse.foldl (fun acc l => nadd acc (nadd (pop (2 * i) l) (pop (2 * i + 1) l))) 0

/-- A cached sum per block, then the sum of the block sums. -/
def cachedGeno (y : Layout) (data : ℕ → ℕ → ℤ) (ids : ℕ → ℕ → ℕ) (i : ℕ) : ℤ :=
  (List.range y.NB).foldl (fun acc k => nadd acc (blockSum nadd y data ids k i 0)) 0

example : revGeno (fun g => absChr 2 cData (cIds g)) 1 = 738 := by decide +kernel
example : cachedGeno lay cData cIds 1 = 453 := by decide +kernel
example : cachedGeno lay cData cIds 1 ≠ fastGenoAll nadd 0 lay cData cIds 1 := by decide +kernel

/-!
## (c) The lowest crossover point

**Second formalization.** `gamRec` restates the original's copying loop as a recursion
that produces the gamete's alleles as a list, carrying the current strand and the locus
index; like the C loop it compares the index with one fixed entry `m` of the map.
-/

section Crossover

variable {F S : Type}

/-- The copying loop of `RecombinationPoissonChr` as a recursion: `cp` is the current
strand (`true`: the second chromosome), `idx` the locus, the last argument the number of
loci left. -/
def gamRec (p1 p2 : ℕ → F) (m : ℕ) : Bool → ℕ → ℕ → List F
  | _, _, 0 => []
  | cp, idx, n + 1 =>
      (if cp then p2 else p1) idx :: gamRec p1 p2 m (if idx = m then !cp else cp) (idx + 1) n

theorem gamRec_eq (p1 p2 : ℕ → F) (m : ℕ) (start : Bool) (n : ℕ) (idx : ℕ) :
    gamRec p1 p2 m (if m < idx then !start else start) idx n =
      (List.range n).map (fun j => gamSpec p1 p2 m start (idx + j)) := by
  induction n generalizing idx with
  | zero => rfl
  | succ n ih =>
    rw [List.range_succ_eq_map, List.map_cons, List.map_map, gamRec]
    congr 1
    · unfold gamSpec
      by_cases h : m < idx
      · have h' : ¬ idx ≤ m := by omega
        cases start <;> simp [h, h']
      · have h' : idx ≤ m := by omega
        cases start <;> simp [h, h']
    · have hcp : (if idx = m then !(if m < idx then !start else start)
          else (if m < idx then !start else start)) = (if m < idx + 1 then !start else start) := by
        rcases Nat.lt_trichotomy m idx with h | h | h
        · simp [h, show idx ≠ m by omega, show m < idx + 1 by omega]
        · subst h; simp
        · simp [show ¬ m < idx by omega, show idx ≠ m by omega, show ¬ m < idx + 1 by omega]
      rw [hcp, ih (idx + 1)]
      apply List.map_congr_left
      intro j _
      simp only [Function.comp, Nat.succ_eq_add_one]
      congr 1; omega

/-- **Second-formalization check (c):** the fold of the main file and the recursion write
the same gamete. -/
theorem origLoop_eq_gamRec (p1 p2 : ℕ → F) (map : List ℕ) (start : Bool) (L : ℕ)
    (chr0 : ℕ → F) :
    (List.range L).map (origLoop p1 p2 map start L chr0) =
      gamRec p1 p2 (map.getD 0 0) start 0 L := by
  have h := gamRec_eq p1 p2 (map.getD 0 0) start L 0
  simp only [Nat.not_lt_zero, if_false, Nat.zero_add] at h
  rw [h]
  apply List.map_congr_left
  intro l hl
  rw [origLoop_eq, if_pos (List.mem_range.mp hl)]
  rfl

/-- **Forward shadow (c1), the whole gamete:** for the same generator state, the gamete
the original's loop writes from its sorted map is, locus for locus, "the starting strand
up to the fast program's lowest point, the other strand after it". -/
theorem shadow_gamete {sort : List ℕ → List ℕ} (hs : SortSpec sort) (count : Rand S ℕ)
    (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) (s : S) (p1 p2 chr0 : ℕ → F) (start : Bool)
    (l : ℕ) (hl : l < L) :
    origLoop p1 p2 (origMap count fresh sort L s).1 start L chr0 l =
      gamSpec p1 p2 (fastLow count fresh L s).1 start l := by
  rw [origLoop_eq, if_pos hl, fastLow_eq hs]
  unfold gamSpec
  have : (l ≤ min L ((origMap count fresh sort L s).1.getD 0 0)) =
      (l ≤ (origMap count fresh sort L s).1.getD 0 0) := by
    apply propext
    constructor
    · intro h; exact le_trans h (Nat.min_le_right _ _)
    · intro h; exact Nat.le_min.mpr ⟨le_of_lt hl, h⟩
  simp only [this]

/-- **Forward shadow (c2), the generator:** the fast draw leaves the generator in the
state the original's leaves it in. -/
theorem shadow_rng_state {sort : List ℕ → List ℕ} (hs : SortSpec sort) (count : Rand S ℕ)
    (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) (s : S) :
    (fastLow count fresh L s).2 = (origMap count fresh sort L s).2 := by
  rw [fastLow_eq hs]

/-- **Forward shadow (c3), no crossover:** if no point is drawn the gamete is a copy of
the starting strand, and the only draw made is the count. -/
theorem shadow_no_points {sort : List ℕ → List ℕ} (hs : SortSpec sort) (count : Rand S ℕ)
    (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) (s : S) (h0 : (count s).1 = 0)
    (p1 p2 chr0 : ℕ → F) (start : Bool) (l : ℕ) (hl : l < L) :
    origLoop p1 p2 (origMap count fresh sort L s).1 start L chr0 l =
        (if start then p2 else p1) l ∧
      fastLow count fresh L s = (L, (count s).2) := by
  have h2 : fastLow count fresh L s = (L, (count s).2) := by simp [fastLow, h0]
  refine ⟨?_, h2⟩
  rw [shadow_gamete hs count fresh L s p1 p2 chr0 start l hl, h2]
  simp [gamSpec, le_of_lt hl]

/-- **Forward shadow (c4), one point:** with one point `p` drawn, the lowest point is
`min L p`, and the generator has made the count's draw and one point's. -/
theorem shadow_one_point (count : Rand S ℕ) (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) (s : S)
    (h1 : (count s).1 = 1) (hL : 1 < L) :
    fastLow count fresh L s =
      (min L (fresh (fun q => decide (q ∈ ([] : List ℕ))) (count s).2).1,
        (fresh (fun q => decide (q ∈ ([] : List ℕ))) (count s).2).2) := by
  have hc : ¬ ((1 : ℕ) = 0 ∨ L ≤ 1) := by omega
  simp only [fastLow, h1, hc, if_false, drawPointsMin, stepMin]
  congr 1
  generalize (fresh (fun q => decide (q ∈ ([] : List ℕ))) (count s).2).1 = p
  rw [Nat.min_def]
  by_cases h : p < L
  · simp [h, show ¬ L ≤ p by omega]
  · simp [h, show L ≤ p by omega]

/-- **Forward shadow (c5):** too many points (`K ≥ NoQTL`) is no crossover. -/
theorem shadow_many_points (count : Rand S ℕ) (fresh : (ℕ → Bool) → Rand S ℕ) (L : ℕ) (s : S)
    (h : L ≤ (count s).1) : fastLow count fresh L s = (L, (count s).2) := by
  simp [fastLow, h]

end Crossover

/-!
**Non-vacuity (c).** The hypothesis on `qsort` is satisfied by insertion sort. The
generator of the examples is a counter: state `s`, raw output `h s`; a point is redrawn
(at most twice) while it is in the map.
-/

/-- Insertion sort satisfies what is assumed of `qsort`. -/
theorem sortSpec_insertionSort : SortSpec (List.insertionSort (· ≤ ·)) :=
  ⟨fun l => List.perm_insertionSort _ l, fun l => List.pairwise_insertionSort _ l⟩

/-- Raw outputs of the example generator. -/
def hRaw (s : ℕ) : ℕ := (s * s + 3 * s + 1) % 5

/-- The number of crossover points: 2 at state 0. -/
def countC : Rand ℕ ℕ := fun s => (2 - s % 2, s + 1)

/-- A point not yet in the map, redrawn at most twice. -/
def freshC (mem : ℕ → Bool) : Rand ℕ ℕ := fun s =>
  if mem (hRaw s) then
    (if mem (hRaw (s + 1)) then (hRaw (s + 2), s + 3) else (hRaw (s + 1), s + 2))
  else (hRaw s, s + 1)

/-!
The raw outputs from state 0 on are 1, 0, 1, 4, 4, 1, … From state 0: the count is 2
(state 1); the first point is 0 (state 2); the second, 1, is new (state 3). The
original's map is `[0, 1, 0]`; the fast program returns 0, in the same state.
-/

example : origMap countC freshC (List.insertionSort (· ≤ ·)) 6 0 = ([0, 1, 0], 3) := by decide +kernel
example : fastLow countC freshC 6 0 = (0, 3) := by decide +kernel

/-- From state 2 a point is rejected and redrawn: the points are 4 and (after a rejected
second 4) 1; the map is sorted to `[1, 4, 0]`, the lowest point is 1, and both programs
have consumed the rejected draw (state 6). -/
example : origMap countC freshC (List.insertionSort (· ≤ ·)) 6 2 = ([1, 4, 0], 6) := by
  decide +kernel
example : fastLow countC freshC 6 2 = (1, 6) := by decide +kernel

/-- Two parental chromosomes with distinguishable alleles. -/
def pA (l : ℕ) : ℤ := 10 + l
def pB (l : ℕ) : ℤ := 20 + l

/-- The gamete for the map `[1, 4, 0]`, starting on the first chromosome: loci 0 and 1
from it, the rest from the second, **although 4 is also in the map** (the original
switches once only). -/
example : (List.range 6).map (origLoop pA pB [1, 4, 0] false 6 (fun _ => 0)) =
    [10, 11, 22, 23, 24, 25] := by decide +kernel
example : gamRec pA pB 1 false 0 6 = [10, 11, 22, 23, 24, 25] := by decide +kernel
example : (List.range 6).map (gamSpec pA pB 1 false) = [10, 11, 22, 23, 24, 25] := by decide +kernel

/-!
**Negative shadows (c).** Switching at the highest point, or at every point (what a
"corrected" loop that advances `PointIndex` would do), or taking the first point drawn
instead of the lowest, each gives a different gamete or point on these examples.
-/

/-- A loop that advances `PointIndex`: every crossover applied. -/
def recombAll (p1 p2 : ℕ → ℤ) (map : List ℕ) (start : Bool) (L : ℕ) : List ℤ :=
  ((List.range L).foldl (fun (st : (Bool × ℕ) × List ℤ) idx =>
      ((if idx = map.getD st.1.2 0 then (!st.1.1, st.1.2 + 1) else st.1),
        st.2 ++ [(if st.1.1 then p2 else p1) idx])) ((start, 0), [])).2

example : recombAll pA pB [1, 4, 0] false 6 = [10, 11, 22, 23, 24, 15] := by decide +kernel
example : recombAll pA pB [1, 4, 0] false 6 ≠
    (List.range 6).map (origLoop pA pB [1, 4, 0] false 6 (fun _ => 0)) := by decide +kernel

/-- Switching after the highest point. -/
example : (List.range 6).map (gamSpec pA pB 4 false) ≠
    (List.range 6).map (origLoop pA pB [1, 4, 0] false 6 (fun _ => 0)) := by decide +kernel

/-- The first point drawn (4, at state 3) is not the lowest (1). -/
example : (freshC (fun _ => false) 3).1 = 4 ∧ (fastLow countC freshC 6 2).1 = 1 := by
  decide +kernel

/-- Drawing the strand before the map instead of after it shifts the generator by one
draw, and the map changes: from state 3 instead of state 2 it is `[4, 0]`. -/
example : (origMap countC freshC (List.insertionSort (· ≤ ·)) 6 3).1 = [4, 0] ∧
    (origMap countC freshC (List.insertionSort (· ≤ ·)) 6 3).1 ≠
      (origMap countC freshC (List.insertionSort (· ≤ ·)) 6 2).1 := by decide +kernel

/-!
## (b) Shared blocks and copy on write
-/

section Blocks

variable {F : Type} (add : F → F → F)

/-- **Forward shadow (b1), other genomes untouched:** after a mutation of chromosome `g`,
every other chromosome of the population stands for the same dense chromosome as before. -/
theorem shadow_others_unchanged (y : Layout) (R : ℕ) {g pos : ℕ} (hg : g < R)
    (hpos : pos < y.L) (eff : F) (x : Store F × (ℕ → ℕ → ℕ)) (hinv : MutInv y R x)
    (g' : ℕ) (hg' : g' < R) (hne : g' ≠ g) (l : ℕ) (hl : l < y.L) :
    absChr y.B (mutStep add y.B g pos eff x).1.data ((mutStep add y.B g pos eff x).2 g') l =
      absChr y.B x.1.data (x.2 g') l := by
  have h := (mutStep_sim add y R hg hpos eff x (fun g l => absChr y.B x.1.data (x.2 g) l) hinv
    (fun _ _ _ _ => rfl)).2 g' hg' l hl
  rw [h]
  simp [mutDense, Function.update_of_ne hne]

/-- **Forward shadow (b2), the mutated genome:** the mutated chromosome changes at the
mutated locus only, to `old ⊕ effect` (the effect on the right, as in `+=`). -/
theorem shadow_mutated (y : Layout) (R : ℕ) {g pos : ℕ} (hg : g < R) (hpos : pos < y.L)
    (eff : F) (x : Store F × (ℕ → ℕ → ℕ)) (hinv : MutInv y R x) (l : ℕ) (hl : l < y.L) :
    absChr y.B (mutStep add y.B g pos eff x).1.data ((mutStep add y.B g pos eff x).2 g) l =
      if l = pos then add (absChr y.B x.1.data (x.2 g) pos) eff
      else absChr y.B x.1.data (x.2 g) l := by
  have h := (mutStep_sim add y R hg hpos eff x (fun g l => absChr y.B x.1.data (x.2 g) l) hinv
    (fun _ _ _ _ => rfl)).2 g hg l hl
  rw [h]
  by_cases hlp : l = pos
  · subst hlp; simp [mutDense]
  · simp [mutDense, hlp]

/-- **Forward shadow (b3), no crossover:** a gamete without a crossover takes every block
number from the chromosome it starts on, and no block is built. -/
theorem shadow_no_crossover (y : Layout) (a : ℕ) (s : Bool) (cur : ℕ → ℕ → ℕ) (k : ℕ)
    (hk : k < y.NB) :
    (gdesc y a y.L s).mix = none ∧
      fillRow cur (gdesc y a y.L s) k = cur (2 * a + (if s then 1 else 0)) k := by
  simp [gdesc, fillRow, hk]

/-- **Forward shadow (b4), collection:** collecting unused blocks changes no block's
contents and no tag. -/
theorem shadow_collect_data (R NB : ℕ) (st : Store F) (t : ℕ → ℕ → ℕ) :
    (collect R NB st t).data = st.data ∧ (collect R NB st t).owner = st.owner := ⟨rfl, rfl⟩

/-- **Forward shadow (b5), a collected block is not in use:** after a collection, the
block `NewBlock` hands out is held by no chromosome of the population. -/
theorem shadow_collect_fresh (R NB : ℕ) (st : Store F) (t : ℕ → ℕ → ℕ)
    (hav : ∀ id, Avail st id → ¬ Refd R NB t id) (tag : ℕ) :
    ¬ Refd R NB t ((collect R NB st t).newBlock tag).1 :=
  (collect_inv R NB st t hav).2 _ (newBlock_avail _ _)

/-- **Non-vacuity (b):** the store both programs start from (every cell holds block 0)
satisfies the invariant. -/
theorem mutInv_init (y : Layout) (R : ℕ) (zero : F) :
    MutInv y R ((⟨fun _ _ => zero, [], 1, fun _ => 0⟩ : Store F), fun _ _ => 0) := by
  refine ⟨⟨List.nodup_nil, by simp⟩, ?_, ?_⟩
  · rintro id (hid | hid) ⟨g, k, _, _, hcell⟩
    · simp at hid
    · have h1 : (1 : ℕ) ≤ id := hid
      have h2 : (0 : ℕ) = id := hcell
      omega
  · intro id hid
    exact absurd rfl hid

end Blocks

/-!
**Numeric instance and negative shadow (b).** Two chromosomes that share all three of
their blocks (the start of every run). A mutation of chromosome 0 at locus 3 with copy on
write changes chromosome 0 only; the same write made in place, without the copy, changes
chromosome 1 as well.
-/

/-- The starting store over the integers. -/
def st0 : Store ℤ := ⟨fun _ _ => 0, [], 1, fun _ => 0⟩

/-- Chromosome `g` after a mutation of chromosome 0 at locus 3 by 5, with copy on write. -/
def afterCow (g : ℕ) : List ℤ :=
  (List.range 6).map (absChr 2 (mutStep nadd 2 0 3 5 (st0, fun _ _ => 0)).1.data
    ((mutStep nadd 2 0 3 5 (st0, fun _ _ => 0)).2 g))

/-- The same write without the copy. -/
def afterInPlace (g : ℕ) : List ℤ :=
  (List.range 6).map (absChr 2
    (Function.update st0.data 0 (Function.update (st0.data 0) 1 (nadd (st0.data 0 1) 5)))
    ((fun _ _ => 0 : ℕ → ℕ → ℕ) g))

example : afterCow 0 = [0, 0, 0, 5, 0, 0] := by decide +kernel
example : afterCow 1 = [0, 0, 0, 0, 0, 0] := by decide +kernel
example : afterInPlace 1 ≠ [0, 0, 0, 0, 0, 0] := by decide +kernel
example : afterInPlace 0 ≠ afterCow 0 := by decide +kernel

/-- A second mutation of the same chromosome in the same block finds the block tagged and
writes in place: no further block is taken (the high-water mark stays at 2). -/
example : (mutStep nadd 2 0 2 7 (mutStep nadd 2 0 3 5 (st0, fun _ _ => 0))).1.hw = 2 := by
  decide +kernel
example : (List.range 6).map (absChr 2
    (mutStep nadd 2 0 2 7 (mutStep nadd 2 0 3 5 (st0, fun _ _ => 0))).1.data
    ((mutStep nadd 2 0 2 7 (mutStep nadd 2 0 3 5 (st0, fun _ _ => 0))).2 0)) =
    [0, 0, 7, 5, 0, 0] := by decide +kernel

/-!
## (d) Two passes

**Second formalization.** `interleavedFold` writes the original's loop as a left fold
over the pair (data, generator state).
-/

section TwoPass

variable {S ι D G : Type}

/-- The interleaved loop as a fold over the pair of data and generator state. -/
def interleavedFold (draw : ι → G → Rand S D) (apply : ι → D → G → G) (items : List ι) (g : G)
    (s : S) : G × S :=
  items.foldl (fun (p : G × S) i => (apply i (draw i p.1 p.2).1 p.1, (draw i p.1 p.2).2)) (g, s)

/-- **Second-formalization check (d):** the recursion of the main file is this fold. -/
theorem interleaved_eq_fold (draw : ι → G → Rand S D) (apply : ι → D → G → G) (items : List ι)
    (g : G) (s : S) : interleaved draw apply items g s = interleavedFold draw apply items g s := by
  induction items generalizing g s with
  | nil => rfl
  | cons i is ih => exact ih _ _

/-- **Forward shadow (d1), the generator:** the interleaved loop and the first pass alone
leave the generator in the same state. -/
theorem shadow_two_pass_state {draw : ι → G → Rand S D} (h : DrawsIgnoreData draw)
    (apply : ι → D → G → G) (g0 : G) (items : List ι) (g : G) (s : S) :
    (interleaved draw apply items g s).2 = (drawAll (fun i => draw i g0) items s).2 := by
  rw [interleaved_eq_two_pass h apply g0]

/-- **Forward shadow (d2), one item:** draw, then apply. -/
theorem shadow_two_pass_one {draw : ι → G → Rand S D} (h : DrawsIgnoreData draw)
    (apply : ι → D → G → G) (g0 : G) (i : ι) (g : G) (s : S) :
    interleaved draw apply [i] g s = (apply i (draw i g0 s).1 g, (draw i g0 s).2) := by
  rw [interleaved_eq_two_pass h apply g0]; rfl

/-- **Forward shadow (d3):** the data that results depends on the generator only through
the list of decisions drawn. -/
theorem shadow_two_pass_data {draw : ι → G → Rand S D} (h : DrawsIgnoreData draw)
    (apply : ι → D → G → G) (g0 : G) (items : List ι) (g : G) (s : S) :
    (interleaved draw apply items g s).1 =
      applyAll apply items (drawAll (fun i => draw i g0) items s).1 g := by
  rw [interleaved_eq_two_pass h apply g0]

end TwoPass

/-!
**Non-vacuity and negative shadows (d).** A draw that ignores the data satisfies the
hypothesis and the two orders agree; a draw that reads the data does not, and the two
orders then differ: the hypothesis cannot be dropped. Making the same draws for the items
in another order also changes the result.
-/

/-- A draw that ignores the data. -/
def drawGood : ℕ → ℕ → Rand ℕ ℕ := fun i _ s => (hRaw (s + i), s + 1)

/-- A draw that reads the data. -/
def drawBad : ℕ → ℕ → Rand ℕ ℕ := fun i g s => (hRaw (s + i) + g, s + 1 + g % 2)

/-- Applying a decision. -/
def applyC : ℕ → ℕ → ℕ → ℕ := fun i d g => 3 * g + d + i

theorem drawGood_ignores : DrawsIgnoreData drawGood := fun _ _ _ _ => rfl

example : interleaved drawGood applyC [0, 1, 2] 7 0 = (210, 3) := by decide +kernel
example : (applyAll applyC [0, 1, 2] (drawAll (fun i => drawGood i 0) [0, 1, 2] 0).1 7,
    (drawAll (fun i => drawGood i 0) [0, 1, 2] 0).2) = (210, 3) := by decide +kernel

example : ¬ DrawsIgnoreData drawBad := by
  intro h
  have := h 0 0 1 0
  revert this
  decide +kernel

example : interleaved drawBad applyC [0, 1, 2] 7 0 ≠
    (applyAll applyC [0, 1, 2] (drawAll (fun i => drawBad i 7) [0, 1, 2] 0).1 7,
      (drawAll (fun i => drawBad i 7) [0, 1, 2] 0).2) := by decide +kernel

/-- The draws made for the items in reverse order. -/
example : (applyAll applyC [0, 1, 2] (drawAll (fun i => drawGood i 0) [2, 1, 0] 0).1.reverse 7) ≠
    (interleaved drawGood applyC [0, 1, 2] 7 0).1 := by decide +kernel

/-!
## (f) One generation and n generations

A complete small instance: 4 diploid individuals, 6 loci in blocks of 2, the lawless
`nadd`, a counter as generator, and arbitrary definitions for everything the model leaves
uninterpreted. The output record of a generation is the list of genotypes followed by the
parameter and the survivors.
-/

/-- The uninterpreted parts, made concrete. -/
def cPr : Prims ℕ ℤ ℤ (List ℤ) where
  popFitness := fun itter par geno env =>
    (par + itter, fun i => geno.getD i 0 + env.getD i 0 + par)
  survives := fun w s => (decide ((w + s) % 5 ≠ 0), s + 1)
  dead := fun n => decide (n < 2)
  parents := fun n s => ((s % n, (s + 1) % n), s + 2)
  count := fun s => (s % 3, s + 1)
  fresh := freshC
  strand := fun s => (decide (s % 2 = 0), s + 1)
  envDraw := fun s => ((s % 7 : ℕ) - 3, s + 1)
  output := fun itter par geno _ _ surv => geno ++ [par, itter] ++ surv.map (fun j => (j : ℤ))
  change := fun itter par => if itter = 1 then par + 100 else par
  mutCount := fun s => (s % 2, s + 1)
  mutPos := fun s => (s % 6, s + 1)
  mutEff := fun s => ((s % 3 : ℕ) - 1, s + 1)

/-- The fast program's starting state. -/
def f0 : FState ℕ ℤ ℤ := ⟨11, ⟨fun _ _ => 0, [], 1, fun _ => 0⟩, fun _ _ => 0, [1, -2, 0, 3], 0⟩

/-- The original's starting state. -/
def o0 : OState ℕ ℤ ℤ := ⟨11, fun _ _ => 0, fun _ _ => 0, [1, -2, 0, 3], 0⟩

/-- Collect unused blocks in every generation, or never, or in odd generations. -/
def gcAlways : ℕ → Bool := fun _ => true
def gcNever : ℕ → Bool := fun _ => false
def gcOdd : ℕ → Bool := fun t => t % 2 == 1

/-- The fast model on the instance. -/
def fastRun (gc : ℕ → Bool) (n : ℕ) : List (List ℤ) × Option (FState ℕ ℤ ℤ) :=
  runGen (fastGen cPr nadd 0 lay 4 gc) n 0 f0

/-- The original model on the instance. -/
def origRun (n : ℕ) : List (List ℤ) × Option (OState ℕ ℤ ℤ) :=
  runGen (origGen cPr (List.insertionSort (· ≤ ·)) nadd 0 lay 4) n 0 o0

theorem cPr_pos : ∀ s, (cPr.mutPos s).1 < lay.L := by
  intro s
  show s % 6 < 6
  exact Nat.mod_lt _ (by decide)

/-- **Non-vacuity (f):** every hypothesis of `run_refines` holds on the instance, so its
conclusion does: for every number of generations and every collection schedule the fast
model prints what the original model prints. -/
theorem instance_outputs (gc : ℕ → Bool) (n : ℕ) : (fastRun gc n).1 = (origRun n).1 :=
  (run_refines cPr nadd 0 lay 4 sortSpec_insertionSort cPr_pos (by decide) gc n 0 f0 o0
    (init_rel 0 lay 4 11 [1, -2, 0, 3] 0 (fun _ _ => 0))).1

/-!
The same, executed (a test of the definitions, independent of the proof): 3 generations
under three collection schedules; the run is not trivial (genotypes become non-zero and
all three schedules print the same 3 lines). Only a few generations are run because the
model represents genomes as functions, which the evaluator recomputes at every use: the
cost grows geometrically with the number of generations. That is a property of executing
the model, not of the C programs.
-/

#guard (fastRun gcAlways 3).1 == (origRun 3).1
#guard (fastRun gcNever 3).1 == (origRun 3).1
#guard (fastRun gcOdd 3).1 == (origRun 3).1
#guard (origRun 3).1.length == 3
#guard (origRun 3).1.any (fun r => (r.take 4).any (· != 0))

#eval (origRun 3).1
#eval (fastRun gcAlways 3).1

/-- The store really is exercised: blocks are allocated, and collection returns some. -/
def storeSize (r : List (List ℤ) × Option (FState ℕ ℤ ℤ)) : ℕ × ℕ :=
  match r.2 with
  | none => (0, 0)
  | some f => (f.st.hw, f.st.free.length)

#eval (storeSize (fastRun gcNever 3), storeSize (fastRun gcAlways 3))
#guard (storeSize (fastRun gcNever 3)).1 > (storeSize (fastRun gcAlways 3)).1

section Run

variable {S F P Rec : Type}

/-- **Forward shadow (f1), zero generations:** nothing is printed and nothing changes. -/
theorem shadow_zero_generations {St : Type} (gen : ℕ → St → Option (Rec × St)) (t : ℕ) (s : St) :
    runGen gen 0 t s = ([], some s) := rfl

/-- **Forward shadow (f2), one more generation at the end:** a run of `n + 1` generations
is the run of `n` generations followed, if the population is alive, by one more
generation (numbered `t + n`) from the state reached. -/
theorem shadow_one_more {St : Type} (gen : ℕ → St → Option (Rec × St)) (n t : ℕ) (s : St) :
    runGen gen (n + 1) t s =
      match (runGen gen n t s).2 with
      | none => ((runGen gen n t s).1, none)
      | some s' =>
        match gen (t + n) s' with
        | none => ((runGen gen n t s).1, none)
        | some p => ((runGen gen n t s).1 ++ [p.1], some p.2) := by
  induction n generalizing t s with
  | zero =>
    simp only [runGen, Nat.add_zero]
    cases gen t s <;> rfl
  | succ n ih =>
    rw [runGen]
    cases hg : gen t s with
    | none => simp [runGen, hg]
    | some p =>
      simp only []
      rw [ih (t + 1) p.2]
      have ht : t + 1 + n = t + (n + 1) := by omega
      have h2 : runGen gen (n + 1) t s =
          (p.1 :: (runGen gen n (t + 1) p.2).1, (runGen gen n (t + 1) p.2).2) := by
        simp [runGen, hg]
      rw [h2, ht]
      cases (runGen gen n (t + 1) p.2).2 with
      | none => rfl
      | some s' =>
        simp only []
        cases gen (t + (n + 1)) s' <;> rfl

variable (pr : Prims S F P Rec) (add : F → F → F) (zero : F) (y : Layout) (N : ℕ)

/-- **Forward shadow (f3), the generator after n generations:** if the population is
alive after `n` generations, the two programs' generators are in the same state. -/
theorem shadow_final_rng {sort : List ℕ → List ℕ} (hs : SortSpec sort)
    (hpos : ∀ s, (pr.mutPos s).1 < y.L) (hN : 0 < N) (gc : ℕ → Bool) (n t : ℕ)
    (f : FState S F P) (o : OState S F P) (h : Rel y N f o) (f' : FState S F P)
    (o' : OState S F P) (hf : (runGen (fastGen pr add zero y N gc) n t f).2 = some f')
    (ho : (runGen (origGen pr sort add zero y N) n t o).2 = some o') :
    f'.rng = o'.rng ∧ f'.env = o'.env ∧ f'.par = o'.par := by
  have h2 := (run_refines pr add zero y N hs hpos hN gc n t f o h).2
  rw [hf, ho] at h2
  exact ⟨h2.rng, h2.env, h2.par⟩

/-- **Forward shadow (f4), extinction:** the population dies within `n` generations in the
fast program if and only if it does in the original. -/
theorem shadow_same_extinction {sort : List ℕ → List ℕ} (hs : SortSpec sort)
    (hpos : ∀ s, (pr.mutPos s).1 < y.L) (hN : 0 < N) (gc : ℕ → Bool) (n t : ℕ)
    (f : FState S F P) (o : OState S F P) (h : Rel y N f o) :
    (runGen (fastGen pr add zero y N gc) n t f).2 = none ↔
      (runGen (origGen pr sort add zero y N) n t o).2 = none := by
  have h2 := (run_refines pr add zero y N hs hpos hN gc n t f o h).2
  cases hf : (runGen (fastGen pr add zero y N gc) n t f).2 <;>
    cases ho : (runGen (origGen pr sort add zero y N) n t o).2 <;>
    rw [hf, ho] at h2 <;> first | exact h2.elim | simp

/-- **Forward shadow (f5), the collection schedule does not matter:** two fast runs that
collect unused blocks in different generations print the same lines. -/
theorem shadow_gc_free {sort : List ℕ → List ℕ} (hs : SortSpec sort)
    (hpos : ∀ s, (pr.mutPos s).1 < y.L) (hN : 0 < N) (gc gc' : ℕ → Bool) (n t : ℕ)
    (f : FState S F P) (o : OState S F P) (h : Rel y N f o) :
    (runGen (fastGen pr add zero y N gc) n t f).1 =
      (runGen (fastGen pr add zero y N gc') n t f).1 := by
  rw [(run_refines pr add zero y N hs hpos hN gc n t f o h).1,
    (run_refines pr add zero y N hs hpos hN gc' n t f o h).1]

/-- **Forward shadow (f6), the block size does not matter:** two fast runs with different
layouts of the same number of loci, each started from a state that stands for the same
original state, print the same lines. -/
theorem shadow_layout_free_run {sort : List ℕ → List ℕ} (hs : SortSpec sort) (y' : Layout)
    (hL : y.L = y'.L) (hpos : ∀ s, (pr.mutPos s).1 < y.L) (hN : 0 < N) (gc gc' : ℕ → Bool)
    (n t : ℕ) (f f' : FState S F P) (o : OState S F P) (h : Rel y N f o) (h' : Rel y' N f' o) :
    (runGen (fastGen pr add zero y N gc) n t f).1 =
      (runGen (fastGen pr add zero y' N gc') n t f').1 := by
  have horig : origGen pr sort add zero y N = origGen pr sort add zero y' N := by
    funext t o; simp only [origGen, hL]
  rw [(run_refines pr add zero y N hs hpos hN gc n t f o h).1,
    (run_refines pr add zero y' N hs (hL ▸ hpos) hN gc' n t f' o h').1, horig]

end Run

/-!
**Negative shadow (f).** A "fast" generation that differs from `fastGen` only in summing
cached block sums (as QTLSim.jl's kernel does) prints different lines on the instance:
`gen_refines` could not have been proved for it.
-/

/-- `fastGen` with cached block sums in place of the in-order sums. -/
def wrongGen (gc : ℕ → Bool) (itter : ℕ) (f : FState ℕ ℤ ℤ) :
    Option (List ℤ × FState ℕ ℤ ℤ) :=
  let geno := (List.range 4).map (cachedGeno lay f.st.data f.ids)
  let st0 := if gc itter then collect 8 lay.NB f.st f.ids else f.st
  let pf := cPr.popFitness itter f.par geno f.env
  let sv := survLoop cPr.survives pf.2 (List.range 4) f.rng
  if cPr.dead sv.1.length then none
  else
    let dr := drawAll (fun _ => fastDecide cPr lay.L sv.1) (List.range 4) sv.2
    let x := fillMix 8 f.ids (descOf lay dr.1) st0
    let mu := mutLoop cPr (mutStep nadd lay.B) 8 x dr.2
    some (cPr.output itter pf.1 geno f.env pf.2 sv.1,
      ⟨mu.2, mu.1.1.clearOwner, mu.1.2, dr.1.map (·.e), cPr.change itter pf.1⟩)

#guard (runGen (wrongGen gcAlways) 3 0 f0).1 != (origRun 3).1

end QTSimFast.Shadows
