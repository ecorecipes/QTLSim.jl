#!/usr/bin/env bash
# A demonstration of bit parity: QTSimFast against the original QTSim.
#
#     ./demo_parity.sh            # builds both programs if needed, a few minutes
#
# For each option file in examples/ and each of DEMO_SEEDS seeds, the original is
# run once and QTSimFast is run with 1, 2, 4 and 8 threads. For every run the
# script prints the SHA-256 and the size of each program's standard output, the
# wall time of each, and IDENTICAL or DIFFERENT (bytes compared with cmp, exit
# codes compared too). It ends with a negative control: the same input with
# another seed must give DIFFERENT, which shows that the comparison can fail.
# The exit code is 0 only if every comparison is IDENTICAL and the control differs.
#
# QTSIM_PAR_MIN=0 is set for the fast program so that the threads are really used
# on these small inputs (by default small populations are summed by one thread).
#
# Environment: DEMO_SEEDS (default "1001 1002"), DEMO_THREADS (default "1 2 4 8"),
# DEMO_TMP (work directory, default a new temporary directory, removed at the end).
# The times are wall-clock seconds on this machine at this moment: an indication,
# not a benchmark (see bench_results.md for that).

set -u
cd "$(dirname "$0")"
HERE=$(pwd)

# Two layouts: inside QTLSim.jl (reference/QTSimFast beside reference/QTSim), or
# the hand-over package (original/ and fast/ under this directory).
if [ -d "$HERE/original/src" ]; then
	ORIGDIR=$HERE/original; FASTDIR=$HERE/fast
	make -s -C "$HERE" all || { echo "build failed"; exit 2; }
else
	ORIGDIR=$HERE/../QTSim; FASTDIR=$HERE
	make -s -C "$ORIGDIR" || { echo "build of the original failed"; exit 2; }
	make -s -C "$FASTDIR" all || { echo "build of QTSimFast failed"; exit 2; }
fi
ORIG=$ORIGDIR/build/QTSim
FAST=$FASTDIR/build/QTSimFast
[ -x "$ORIG" ] && [ -x "$FAST" ] || { echo "programs not built"; exit 2; }

SEEDS=${DEMO_SEEDS:-"1001 1002"}
THREADS=${DEMO_THREADS:-"1 2 4 8"}
if [ -n "${DEMO_TMP:-}" ]; then WORK=$DEMO_TMP; mkdir -p "$WORK"; KEEPWORK=1
else WORK=$(mktemp -d "${TMPDIR:-/tmp}/qtsim-demo.XXXXXX") || exit 2; KEEPWORK=0; fi

sha() {
	if command -v shasum > /dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1
	elif command -v sha256sum > /dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
	elif command -v openssl > /dev/null 2>&1; then openssl dgst -sha256 "$1" | sed 's/.*= *//'
	else echo "no-sha256-tool"; fi
}
bytes() { wc -c < "$1" | tr -d ' '; }

# run <outfile> <command...>: sets RC and SECS (wall seconds)
TIMEFORMAT=%R
run() {
	local out=$1; shift
	{ time "$@" > "$out" 2> /dev/null; } 2> "$WORK/time.txt"
	RC=$?
	SECS=$(tail -1 "$WORK/time.txt")
}

echo "QTSimFast parity demonstration"
echo "date:      $(date '+%Y-%m-%d %H:%M %Z')"
echo "machine:   $(uname -sm)"
echo "compiler:  $(${CC:-cc} --version 2> /dev/null | head -1)"
echo "GSL:       $(gsl-config --version 2> /dev/null)"
echo "original:  ${ORIG#$HERE/}"
echo "fast:      ${FAST#$HERE/}"
echo "sources:   original QTSim.c sha256 $(sha "$ORIGDIR/src/QTSim.c")"
echo "           fast     QTSim.c sha256 $(sha "$FASTDIR/src/QTSim.c")"
echo

nsame=0; ndiff=0
for ex in "$HERE"/examples/*.txt; do
	name=$(basename "$ex" .txt)
	echo "== $name: $(sed -n '1s/^# *//p' "$ex")"
	for seed in $SEEDS; do
		f=$WORK/$name-$seed.txt
		grep -v '^[Ss]eed' "$ex" > "$f"; echo "Seed $seed" >> "$f"
		run "$WORK/orig.out" "$ORIG" "$f"; orc=$RC; osecs=$SECS
		osha=$(sha "$WORK/orig.out"); obytes=$(bytes "$WORK/orig.out")
		last=$(tail -1 "$WORK/orig.out" | cut -f1)
		echo "   seed $seed  original          sha256 $osha  $obytes bytes  exit $orc  ${osecs}s  (last line: $last)"
		for t in $THREADS; do
			run "$WORK/fast.out" env QTSIM_THREADS=$t QTSIM_PAR_MIN=0 "$FAST" "$f"
			fsha=$(sha "$WORK/fast.out"); fbytes=$(bytes "$WORK/fast.out")
			if cmp -s "$WORK/orig.out" "$WORK/fast.out" && [ "$orc" = "$RC" ]; then verdict=IDENTICAL; nsame=$((nsame + 1))
			else verdict=DIFFERENT; ndiff=$((ndiff + 1)); fi
			echo "   seed $seed  fast, $t thread(s)  sha256 $fsha  $fbytes bytes  exit $RC  ${SECS}s  $verdict"
		done
	done
done

echo
echo "== negative control: the same input with another seed must differ"
ctl=0
ex=$HERE/examples/design_0100_loci.txt
f1=$WORK/control-1.txt; f2=$WORK/control-2.txt
grep -v '^[Ss]eed' "$ex" > "$f1"; echo "Seed 1001" >> "$f1"
grep -v '^[Ss]eed' "$ex" > "$f2"; echo "Seed 1003" >> "$f2"
run "$WORK/c1.out" "$ORIG" "$f1"
run "$WORK/c2.out" env QTSIM_THREADS=1 "$FAST" "$f2"
# the two option files differ in the seed, which both programs print in their
# header: compare the simulation's lines only (from the line "Generation ...")
sed -n '/^Generation/,$p' "$WORK/c1.out" > "$WORK/c1.body"
sed -n '/^Generation/,$p' "$WORK/c2.out" > "$WORK/c2.body"
if cmp -s "$WORK/c1.body" "$WORK/c2.body"; then cv=IDENTICAL; else cv=DIFFERENT; ctl=1; fi
echo "   original, seed 1001   sha256 $(sha "$WORK/c1.body")  $(bytes "$WORK/c1.body") bytes (lines from 'Generation')"
echo "   fast,     seed 1003   sha256 $(sha "$WORK/c2.body")  $(bytes "$WORK/c2.body") bytes (lines from 'Generation')"
echo "   $cv (expected: DIFFERENT)"
run "$WORK/c3.out" env QTSIM_THREADS=1 "$FAST" "$f1"
sed -n '/^Generation/,$p' "$WORK/c3.out" > "$WORK/c3.body"
if cmp -s "$WORK/c1.body" "$WORK/c3.body"; then cv2=IDENTICAL; else cv2=DIFFERENT; ctl=0; fi
echo "   fast,     seed 1001   sha256 $(sha "$WORK/c3.body")  $(bytes "$WORK/c3.body") bytes (lines from 'Generation')"
echo "   $cv2 with the original at seed 1001 (expected: IDENTICAL)"

[ "$KEEPWORK" = 0 ] && rm -rf "$WORK"

echo
echo "summary: $nsame comparisons IDENTICAL, $ndiff DIFFERENT; negative control $([ $ctl = 1 ] && echo passed || echo FAILED)"
if [ "$ndiff" != 0 ] || [ "$ctl" != 1 ] || [ "$nsame" = 0 ]; then echo "PARITY DEMONSTRATION FAILED"; exit 1; fi
echo "PARITY DEMONSTRATION OK"
exit 0
