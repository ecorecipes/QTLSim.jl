import QTSimFast

/-!
# The Lean model, executed on the instance of `model_check.c`

`lean --run ModelCheck.lean` evaluates the definitions of `QTSimFast.lean` (the model of
the original and the model of the fast program) on the tiny instance that
`formal/model_check.c` runs through the two C programs' own routines: 4 diploid
individuals, 6 loci in blocks of 2, the same alleles, the same five sets of random
decisions. `F` is Lean's `Float` (the processor's IEEE doubles) and `add` its addition.
It prints the genotypes and gametes as the 64-bit patterns of the doubles, in the format
of `model_check.c`; `make model-check` compares the outputs byte for byte.

This is a test, not a proof: it narrows the gap between the model and the C code on one
instance. The random generator here is a list of decisions (number of points, points as
drawn including a rejected repeat, strand); in C it is a list of raw 32-bit outputs chosen
so that GSL's routines return those decisions.
-/

open QTSimFast

/-- Magnitudes that differ by 2^53 from locus to locus, so that rounding depends on the
order of the additions. -/
def scale (l : Nat) : Float :=
  match l with
  | 0 => 1.0
  | 1 => Float.ofNat 9007199254740992
  | 2 => 3.0
  | 3 => -(Float.ofNat 9007199254740992)
  | 4 => 0.5
  | _ => 7.0

/-- Allele of chromosome `g` at locus `l`: variant 0 exact, variant 1 rounded. -/
def val (variant g l : Nat) : Float :=
  if variant = 0 then Float.ofNat (g * 10 + l)
  else scale l + 1.0 / Float.ofNat (g + 3)

def lay : Layout := ⟨2, 2, 2, by decide, by decide, by decide⟩

/-- The generator as a script of decisions. -/
def countS : Rand (List Nat) Nat := fun s => (s.headD 0, s.tail)

def freshS (mem : Nat → Bool) : Rand (List Nat) Nat
  | [] => (0, [])
  | p :: r => if mem p then freshS mem r else (p, r)

def strandS : Rand (List Nat) Bool := fun s => (s.headD 0 != 0, s.tail)

/-- The five cases of `model_check.c`: the count, the points as drawn, the strand
(1: the parent's second chromosome first). -/
def cases : List (List Nat) :=
  [[2, 4, 4, 1, 1], [0, 0], [1, 3, 0], [3, 2, 4, 3, 1], [6, 1]]

def sortN (l : List Nat) : List Nat := l.insertionSort (· ≤ ·)

/-- The dense population of the original's model. -/
def pop (variant : Nat) (g l : Nat) : Float := val variant g l

/-- The block store of the fast model: chromosome `g` holds blocks `3g`, `3g+1`, `3g+2`. -/
def ids (g k : Nat) : Nat := 3 * g + k
def data (variant : Nat) (id o : Nat) : Float := val variant (id / 3) ((id % 3) * 2 + o)
def store (variant : Nat) : Store Float := ⟨data variant, [], 24, fun _ => 0⟩

/-- The gamete of individual 1 in the original's model: `CreateCombMap`, the strand, the
copying loop. Also returns whether the script was used up. -/
def gameteOrig (variant : Nat) (script : List Nat) : List Float × Bool :=
  let m := origMap countS freshS sortN 6 script
  let st := strandS m.2
  ((List.range 6).map (origLoop (pop variant 2) (pop variant 3) m.1 st.1 6 (fun _ => 0)),
    st.2.isEmpty)

/-- The same gamete in the fast model: `LowestCrossover`, the strand, `FillGametes` and
`MixBlocks` over the block store. -/
def gameteFast (variant : Nat) (script : List Nat) : List Float × Bool :=
  let lo := fastLow countS freshS 6 script
  let st := strandS lo.2
  let desc : Nat → GDesc := fun g => if g = 0 then gdesc lay 1 lo.1 st.1 else ⟨0, 0, 3, none⟩
  let x := fillMix 8 ids desc (store variant)
  ((List.range 6).map (absChr 2 x.1.data (x.2 0)), st.2.isEmpty)

def bits (x : Float) : String := toString x.toBits.toNat

def report (variant : Nat) (geno : Nat → Float) (gam : List Nat → List Float × Bool) : String :=
  let head := s!"variant {if variant = 0 then "exact" else "inexact"}\n"
  let g := (List.range 4).map (fun i => s!"geno {i} {bits (geno i)}\n")
  let c := (List.zip (List.range 5) cases).map (fun p =>
    s!"gamete {p.1 + 1} " ++ " ".intercalate ((gam p.2).1.map bits) ++ "\n")
  head ++ String.join g ++ String.join c

def main : IO UInt32 := do
  let mut ok := true
  for variant in [0, 1] do
    let o := report variant (origGeno (· + ·) 0.0 6 (pop variant)) (gameteOrig variant)
    let f := report variant (fastGenoAll (· + ·) 0.0 lay (data variant) ids) (gameteFast variant)
    IO.print o
    if o != f then
      IO.eprintln s!"ModelCheck: the fast model and the original's model differ (variant {variant})"
      ok := false
    for c in cases do
      if !(gameteOrig variant c).2 || !(gameteFast variant c).2 then
        IO.eprintln "ModelCheck: a case left decisions unused"
        ok := false
  -- the rounded variant must be one in which the order of the additions matters:
  -- compare the in-order sum with the reversed sum and with the sum of block sums
  let mut differs := 0
  for i in [0, 1, 2, 3] do
    let fwd := origGeno (· + ·) 0.0 6 (pop 1) i
    let rev := (List.range 6).reverse.foldl
      (fun acc l => acc + (pop 1 (2 * i) l + pop 1 (2 * i + 1) l)) 0.0
    let blk := (List.range 3).foldl (fun acc k =>
      acc + ((0.0 + (pop 1 (2 * i) (2 * k) + pop 1 (2 * i + 1) (2 * k)))
        + (pop 1 (2 * i) (2 * k + 1) + pop 1 (2 * i + 1) (2 * k + 1)))) 0.0
    if fwd.toBits != rev.toBits then differs := differs + 1
    if fwd.toBits != blk.toBits then differs := differs + 1
    IO.eprintln s!"note: individual {i}, rounded alleles: in-order {bits fwd}, reversed {bits rev}, sum of block sums {bits blk}"
  if differs == 0 then
    IO.eprintln "ModelCheck: the rounded variant does not distinguish the order of the additions"
    ok := false
  return if ok then 0 else 1
