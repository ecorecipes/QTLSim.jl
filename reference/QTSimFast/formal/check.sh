#!/usr/bin/env bash
# Type-check the Lean files of this directory, standalone, and run the model check.
#
#     formal/check.sh /path/to/a/Lake/project/that/depends/on/Mathlib
#
# The files need only Mathlib (they were checked with Lean 4.32.1 and Mathlib v4.32.1;
# see README.md here). The project is used for its toolchain and its built Mathlib
# only: nothing is written into it. Steps:
#   1. QTSimFast.lean is checked and compiled to an .olean in a temporary directory;
#   2. QTSimFastShadows.lean (which imports it) is checked;
#   3. the axioms the main theorems depend on are printed (expected: propext,
#      Classical.choice, Quot.sound, the three standard ones; no sorryAx);
#   4. ModelCheck.lean is run and its output written to model_check_lean.txt in the
#      temporary directory (or to the file named by MODEL_CHECK_OUT).
# Exit code 0 only if every step succeeds.

set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
PROJ=${1:-${MATHLIB_PROJECT:-}}
[ -n "$PROJ" ] && [ -d "$PROJ" ] || { echo "usage: $0 <Lake project with Mathlib>"; exit 2; }
OUT=$(mktemp -d "${TMPDIR:-/tmp}/qtsimfast-lean.XXXXXX")
trap 'rm -rf "$OUT"' EXIT

cd "$PROJ"
echo "Lean: $(lake env lean --version)"
echo "1. QTSimFast.lean"
lake env lean --root="$HERE" -o "$OUT/QTSimFast.olean" "$HERE/QTSimFast.lean"
echo "2. QTSimFastShadows.lean"
OUT="$OUT" HERE="$HERE" lake env sh -c 'LEAN_PATH="$LEAN_PATH:$OUT" lean "$HERE/QTSimFastShadows.lean"'
echo "3. axioms"
cat > "$OUT/Axioms.lean" <<'EOT'
import QTSimFast
#print axioms QTSimFast.foldl_pointwise
#print axioms QTSimFast.fastGeno_eq
#print axioms QTSimFast.origLoop_eq
#print axioms QTSimFast.runMin_eq_head
#print axioms QTSimFast.fastLow_eq
#print axioms QTSimFast.interleaved_eq_two_pass
#print axioms QTSimFast.fillMix_spec
#print axioms QTSimFast.mutStep_sim
#print axioms QTSimFast.collect_inv
#print axioms QTSimFast.gen_refines
#print axioms QTSimFast.run_refines
#print axioms QTSimFast.init_rel
EOT
OUT="$OUT" lake env sh -c 'LEAN_PATH="$LEAN_PATH:$OUT" lean "$OUT/Axioms.lean"' | tee "$OUT/axioms.txt"
if grep -q "sorryAx" "$OUT/axioms.txt"; then echo "sorry found"; exit 1; fi
echo "4. ModelCheck.lean"
MC=${MODEL_CHECK_OUT:-$OUT/model_check_lean.txt}
OUT="$OUT" HERE="$HERE" lake env sh -c 'LEAN_PATH="$LEAN_PATH:$OUT" lean --run "$HERE/ModelCheck.lean"' > "$MC"
if [ -f "$HERE/model_check_expected.txt" ]; then
	cmp "$MC" "$HERE/model_check_expected.txt" && echo "model check: the Lean model's output is model_check_expected.txt"
fi
echo "LEAN CHECK OK"
