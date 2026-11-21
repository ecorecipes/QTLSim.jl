#!/usr/bin/env bash
# Timings of QTSimFast's fast modes (FAST.md) against the original QTSim and the exact
# mode, on bench.sh's inputs: per stratum the published input with the slowest move among
# the first 400 (or a synthetic input of the same design), run to 2000 generations after
# the move, BENCH_SEEDS seeds (default 3). Run it on a quiet machine:
#
#     ./bench_fast.sh
#
# Writes bench_fast.md (medians) and bench_fast_runs.csv (every run). Programs: the
# original; QTSimFast exact on 1 and 4 threads; QTSIM_FAST=1 on 1 thread; QTSIM_FAST=2 on
# 1, 2 and 4 threads. For the exact mode and level 1, whether the output is identical to
# the original's is recorded (level 2's is not expected to be).

set -u
cd "$(dirname "$0")"
HERE=$(pwd)
PKG=$(cd "$HERE/../.." && pwd)
if [ -d "$HERE/original/src" ]; then ORIGDIR=$HERE/original; FASTDIR=$HERE/fast
else ORIGDIR=$HERE/../QTSim; FASTDIR=$HERE; fi
ORIG=${QTSIM:-$ORIGDIR/build/QTSim}
FAST=$FASTDIR/build/QTSimFast
SIMS=${QTSIM_SIMULATIONS:-$PKG/../Repository/Simulations}
NSEED=${BENCH_SEEDS:-3}
STRATA=${BENCH_STRATA:-"QT0100-Possion3.0 QT0100-Possion5.0 QT0500-Possion3.0 QT0500-Possion5.0 QT1000-Possion3.0 QT1000-Possion5.0 QT5000-Possion3.0 QT5000-Possion5.0"}
AFTER=${BENCH_AFTER:-2000}
WORK=$FASTDIR/build/bench_fast
OUT=${BENCH_OUT:-$HERE}

[ -x "$ORIG" ] || { echo "no original at $ORIG (run make)"; exit 2; }
[ -x "$FAST" ] || { echo "no $FAST (make)"; exit 2; }
mkdir -p "$WORK"
RUNS=$OUT/bench_fast_runs.csv
echo "stratum,input,program,seed,wall_s,cpu_s,generations,identical" > "$RUNS"
load_before=$(uptime | sed 's/.*load average[s]*: *//')

timed() {
	out=$1; shift
	/usr/bin/time -p "$@" > "$out" 2> "$WORK/time.txt"
	WALL=$(awk '/^real/{print $2}' "$WORK/time.txt")
	CPU=$(awk '/^user/{u=$2} /^sys/{s=$2} END{printf "%.2f", u+s}' "$WORK/time.txt")
}
gens() { grep -c '^[0-9]' "$1"; }

# label : environment
PROGS="exact1:QTSIM_FAST=0,QTSIM_THREADS=1 exact4:QTSIM_FAST=0,QTSIM_THREADS=4 fast1:QTSIM_FAST=1,QTSIM_THREADS=1 fast2t1:QTSIM_FAST=2,QTSIM_THREADS=1 fast2t2:QTSIM_FAST=2,QTSIM_THREADS=2 fast2t4:QTSIM_FAST=2,QTSIM_THREADS=4"

for st in $STRATA; do
	d=$SIMS/$st
	loci=$(echo "$st" | sed -e 's/^QT0*//' -e 's/-.*//')
	rec=$(echo "$st" | sed -e 's/.*Possion//' -e 's/\.0$//')
	if [ -d "$d" ]; then
		src=$(for i in $(ls "$d" | head -400); do echo "$(awk '/^MoveFitnessSD/{print $3}' "$d/$i") $d/$i"; done | sort -g | head -1 | cut -d' ' -f2)
		input=$(basename "$src")
	else
		src=$WORK/synthetic-$st.txt; input=synthetic
		printf 'Sample 1\nPopSize 1000\nNoQTL %s\nGenerations 250000\nPloidy 2\nMoveFitnessSD 10000 0.001\nPhenotypeVar 0.02\nMutationScalar 0.001\nHeritability 0.6\nSeed 1\nInitFitnessSD 1.0\nOptimumMoveType continuous\nMutationRatePerGamete 0.0001\nRecombinationPoissonMean %s\n' "$loci" "$rec" > "$src"
	fi
	move=$(awk '/^MoveFitnessSD/{print $2}' "$src")
	gen=$((move + AFTER))
	for seed in $(seq 1001 $((1000 + NSEED))); do
		f=$WORK/$st-$seed.txt
		sed -e "s/^Generations .*/Generations $gen/" -e "s/^Seed .*/Seed $seed/" "$src" > "$f"
		echo "$st $input seed $seed" >&2
		timed "$WORK/orig.out" "$ORIG" "$f"
		echo "$st,$input,original,$seed,$WALL,$CPU,$(gens "$WORK/orig.out")," >> "$RUNS"
		for p in $PROGS; do
			label=${p%%:*}; envs=$(echo "${p#*:}" | tr ',' ' ')
			timed "$WORK/fast.out" env $envs "$FAST" "$f"
			same=
			case $label in exact*|fast1) cmp -s "$WORK/orig.out" "$WORK/fast.out" && same=yes || same=NO ;; esac
			echo "$st,$input,$label,$seed,$WALL,$CPU,$(gens "$WORK/fast.out"),$same" >> "$RUNS"
		done
	done
done

load_after=$(uptime | sed 's/.*load average[s]*: *//')
med() { awk -F, -v s="$1" -v p="$2" -v c="$3" '$1==s && $3==p {print $c}' "$RUNS" | sort -g | awk '{a[NR]=$1} END{if(NR%2) print a[(NR+1)/2]; else printf "%.2f\n", (a[NR/2]+a[NR/2+1])/2}'; }
ratio() { awk -v a="$1" -v b="$2" 'BEGIN{ if(b > 0) printf "%.1f", a/b; else print "-" }'; }

commit=$(git -C "$HERE" rev-parse --short HEAD 2> /dev/null)
{
	echo "# QTSimFast fast modes: timings"
	echo
	echo "- Date: $(date '+%Y-%m-%d %H:%M %Z'); machine: $(sysctl -n machdep.cpu.brand_string 2> /dev/null || uname -m); load before: $load_before; after: $load_after."
	echo "- Compiler: $(${CC:-cc} --version | head -1); GSL $(gsl-config --version); commit ${commit:-unknown}."
	echo "- Inputs as in bench.sh: per stratum the published input with the slowest move among the first 400, run to $AFTER generations after the move; median of $NSEED seeds."
	echo "- Wall seconds; in brackets the speed-up in wall time against the original. Exact and level 1 outputs compared with the original's: identical in $(awk -F, '$8=="yes"' "$RUNS" | wc -l | tr -d ' ') of $(awk -F, '$8=="yes" || $8=="NO"' "$RUNS" | wc -l | tr -d ' ') runs."
	echo "- Level 2 draws other random numbers: its runs are different runs (the same statistics in distribution), so a run can last a different number of generations."
	echo
	echo "| Stratum | Original | Exact, 1 thread | Exact, 4 threads | Level 1 | Level 2, 1 thread | Level 2, 2 threads | Level 2, 4 threads |"
	echo "|---|---|---|---|---|---|---|---|"
	for st in $STRATA; do
		o=$(med $st original 5)
		line="| ${st/-Possion/-} | $o |"
		for p in exact1 exact4 fast1 fast2t1 fast2t2 fast2t4; do
			w=$(med $st $p 5); line="$line $w ($(ratio $o $w)x) |"
		done
		echo "$line"
	done
	echo
	echo "CPU seconds (user + system), medians:"
	echo
	echo "| Stratum | Original | Exact, 1 | Exact, 4 | Level 1 | Level 2, 1 | Level 2, 2 | Level 2, 4 |"
	echo "|---|---|---|---|---|---|---|---|"
	for st in $STRATA; do
		line="| ${st/-Possion/-} | $(med $st original 6) |"
		for p in exact1 exact4 fast1 fast2t1 fast2t2 fast2t4; do line="$line $(med $st $p 6) |"; done
		echo "$line"
	done
} > "$OUT/bench_fast.md"
echo "wrote $OUT/bench_fast.md and $RUNS"
grep -q ',NO$' "$RUNS" && { echo "note: some exact or level-1 outputs differ from the original's"; exit 1; }
exit 0
