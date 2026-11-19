#!/usr/bin/env bash
# Timings of the original QTSim, QTSimFast (1, 2, 4 and 8 threads) and the Julia
# kernel of QTLSim.jl (one thread) on the same inputs: one published input per
# stratum (number of loci x recombination mean), the full published protocol up
# to 2000 generations after the optimum starts to move, BENCH_SEEDS seeds each.
#
#     ./bench.sh            # about half an hour; run it on a quiet machine
#
# Writes bench_results.md (one table: medians, wall and CPU seconds, speed-ups
# against the original, the machine and its load), bench_results.csv (the same
# medians) and bench_runs.csv (every run). Runs are made one at a time.
#
# The inputs: in each stratum, the input with the slowest-moving optimum among
# the first 400 of the paper's (Repository/Simulations, or QTSIM_SIMULATIONS), so
# that every program runs all the generations; without the paper's inputs, a
# synthetic input of the same design. QTSimFast's output is compared with the
# original's on every run (it must be identical). The Julia kernel runs the same
# options with the same number of seeds; it has its own random numbers, so its
# runs are other realisations of the same model, and nothing is compared.
#
# Environment: BENCH_SEEDS (default 3), BENCH_THREADS (default "1 2 4 8"),
# BENCH_STRATA (default all eight), BENCH_AFTER (generations after the move
# starts, default 2000), BENCH_JULIA=0 to leave Julia out, BENCH_ABLATION=1 to
# add a second table (QTSimFast with one optimisation switched off at a time),
# JULIA (default julia), QTSIM (the original, default ../QTSim/build/QTSim),
# BENCH_OUT (where the three files are written, default this directory),
# QTLSIM_PROJECT (a Julia environment with QTLSim.jl; without Julia or without it
# the Julia column is left out).

set -u
cd "$(dirname "$0")"
HERE=$(pwd)
PKG=$(cd "$HERE/../.." && pwd)

# Two layouts: inside QTLSim.jl, or the hand-over package (original/ and fast/ here).
if [ -d "$HERE/original/src" ]; then ORIGDIR=$HERE/original; FASTDIR=$HERE/fast
else ORIGDIR=$HERE/../QTSim; FASTDIR=$HERE; fi
ORIG=${QTSIM:-$ORIGDIR/build/QTSim}
FAST=$FASTDIR/build/QTSimFast
# The Julia column needs QTLSim.jl: its vignettes environment (QTLSIM_PROJECT overrides).
JPROJ=${QTLSIM_PROJECT:-$PKG/vignettes}
SIMS=${QTSIM_SIMULATIONS:-$PKG/../Repository/Simulations}
NSEED=${BENCH_SEEDS:-3}
THREADS=${BENCH_THREADS:-"1 2 4 8"}
STRATA=${BENCH_STRATA:-"QT0100-Possion3.0 QT0100-Possion5.0 QT0500-Possion3.0 QT0500-Possion5.0 QT1000-Possion3.0 QT1000-Possion5.0 QT5000-Possion3.0 QT5000-Possion5.0"}
AFTER=${BENCH_AFTER:-2000}
JULIA=${JULIA:-julia}
WORK=$FASTDIR/build/bench

[ -x "$ORIG" ] || { echo "no original at $ORIG (run make)"; exit 2; }
[ -x "$FAST" ] || { echo "no $FAST (make)"; exit 2; }
mkdir -p "$WORK"
OUT=${BENCH_OUT:-$HERE}
RUNS=$OUT/bench_runs.csv
echo "stratum,input,program,threads,seed,wall_s,cpu_s,generations,identical" > "$RUNS"

load_before=$(uptime | sed 's/.*load average[s]*: *//')

# timed <outfile> <command...>: sets WALL and CPU (user + system) in seconds
timed() {
	out=$1; shift
	/usr/bin/time -p "$@" > "$out" 2> "$WORK/time.txt"
	WALL=$(awk '/^real/{print $2}' "$WORK/time.txt")
	CPU=$(awk '/^user/{u=$2} /^sys/{s=$2} END{printf "%.2f", u+s}' "$WORK/time.txt")
}
# the generations a run's output has (lines that start with a number)
gens() { grep -c '^[0-9]' "$1"; }

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
	seed=1001
	while [ $seed -lt $((1001 + NSEED)) ]; do
		f=$WORK/$st-$seed.txt
		sed -e "s/^Generations .*/Generations $gen/" -e "s/^Seed .*/Seed $seed/" "$src" > "$f"
		echo "$st $input seed $seed: original" >&2
		timed "$WORK/orig.out" "$ORIG" "$f"
		echo "$st,$input,original,1,$seed,$WALL,$CPU,$(gens "$WORK/orig.out")," >> "$RUNS"
		for t in $THREADS; do
			timed "$WORK/fast.out" env QTSIM_THREADS=$t "$FAST" "$f"
			cmp -s "$WORK/orig.out" "$WORK/fast.out" && same=yes || same=NO
			echo "$st,$input,fast,$t,$seed,$WALL,$CPU,$(gens "$WORK/fast.out"),$same" >> "$RUNS"
		done
		if [ "${BENCH_ABLATION:-0}" = 1 ]; then
			# one optimisation off at a time (one thread)
			for a in "lanes1:QTSIM_LANES=1" "dense:QTSIM_BLOCK=1000000" "gc8:QTSIM_GC=8"; do
				timed "$WORK/fast.out" env QTSIM_THREADS=1 "${a#*:}" "$FAST" "$f"
				cmp -s "$WORK/orig.out" "$WORK/fast.out" && same=yes || same=NO
				echo "$st,$input,fast-${a%%:*},1,$seed,$WALL,$CPU,$(gens "$WORK/fast.out"),$same" >> "$RUNS"
			done
			if [ -x "$FAST-plain" ]; then
				timed "$WORK/fast.out" "$FAST-plain" "$f"
				cmp -s "$WORK/orig.out" "$WORK/fast.out" && same=yes || same=NO
				echo "$st,$input,fast-plain,1,$seed,$WALL,$CPU,$(gens "$WORK/fast.out"),$same" >> "$RUNS"
			fi
		fi
		seed=$((seed + 1))
	done
done

# ---------------------------------------------------------------- the Julia kernel
julia_version="not run"
if [ "${BENCH_JULIA:-1}" = 1 ] && command -v "$JULIA" > /dev/null 2>&1 && [ -f "$JPROJ/Project.toml" ]; then
	cat > "$WORK/bench.jl" <<'EOF'
# Times QTLSim's kernel (simulate) on the option files of bench.sh; one thread.
using QTLSim
work, runs, strata, nseed = ARGS[1], ARGS[2], split(ARGS[3]), parse(Int, ARGS[4])
cpu() = ccall(:clock, Clong, ()) / 1e6
first_file = joinpath(work, "$(strata[1])-1001.txt")
let p = read_qtsim_options(first_file)      # compile on a short run that is not timed
    simulate(Params(p; NoGenerations = 200); seed = 1)
end
open(runs, "a") do io
    for st in strata, seed in 1001:(1000 + nseed)
        p = read_qtsim_options(joinpath(work, "$st-$seed.txt"))
        c0 = cpu(); t0 = time_ns()
        r = simulate(p; seed)
        wall = (time_ns() - t0) / 1e9; c = cpu() - c0
        gens = length(r.gen)
        println(io, "$st,,julia,1,$seed,$(round(wall, digits = 2)),$(round(c, digits = 2)),$gens,")
        println(stderr, "$st seed $seed: Julia kernel $(round(wall, digits = 2)) s")
    end
end
EOF
	julia_version=$("$JULIA" --version 2> /dev/null)
	"$JULIA" -t 1 --project="$JPROJ" "$WORK/bench.jl" "$WORK" "$RUNS" "$STRATA" "$NSEED" || echo "note: the Julia runs failed; the table will lack them" >&2
fi

load_after=$(uptime | sed 's/.*load average[s]*: *//')

# ---------------------------------------------------------------- the tables
# med <stratum> <program> <threads> <column>: the median over seeds
med() {
	awk -F, -v s="$1" -v p="$2" -v t="$3" -v c="$4" '$1==s && $3==p && $4==t {print $c}' "$RUNS" | sort -g | awk '{a[NR]=$1} END{if(NR==0){print "NA"} else if(NR%2){print a[(NR+1)/2]} else {printf "%.2f\n", (a[NR/2]+a[NR/2+1])/2}}'
}
ratio() { awk -v a="$1" -v b="$2" 'BEGIN{if(a=="NA"||b=="NA"||b+0==0) print "NA"; else printf "%.1f", a/b}'; }

CSV=$OUT/bench_results.csv
MD=$OUT/bench_results.md
hdr="stratum,loci,recombination_mean,generations,original_wall,original_cpu"
for t in $THREADS; do hdr="$hdr,fast${t}_wall,fast${t}_cpu,fast${t}_speedup_wall"; done
hdr="$hdr,fast1_speedup_cpu,julia_wall,julia_cpu,julia_generations,julia_speedup_wall,all_identical"
echo "$hdr" > "$CSV"

cpu_name=$(sysctl -n machdep.cpu.brand_string 2> /dev/null || awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo 2> /dev/null)
ncpu=$(getconf _NPROCESSORS_ONLN 2> /dev/null)
commit=$(cd "$PKG" && git rev-parse --short HEAD 2> /dev/null)
{
	echo "# QTSimFast benchmark"
	echo
	echo "Written by \`bench.sh\` on $(date '+%Y-%m-%d %H:%M %Z')."
	echo
	echo "- Machine: ${cpu_name:-unknown}, $ncpu processors; $(uname -sm)."
	echo "- Load averages (1, 5, 15 min) before: $load_before; after: $load_after. Runs were made one at a time."
	echo "- Compiler: $(${CC:-cc} --version | head -1); GSL $(gsl-config --version); Julia: $julia_version; commit ${commit:-unknown}."
	echo "- Protocol: the published inputs' options with \`Generations\` = the move's first generation + $AFTER; in each stratum the input with the slowest-moving optimum among the first 400; $NSEED seeds (1001 ...); medians over seeds."
	echo "- Each cell: wall seconds / CPU seconds (user + system), and for all but the original the speed-up in wall time against the original."
	echo "- QTSimFast's output was compared with the original's on every run; the Julia kernel has its own random numbers (other realisations of the same model, nothing compared)."
	echo
	h="| Stratum | Generations | Original C |"; l="|---|---|---|"
	for t in $THREADS; do h="$h QTSimFast, $t thread(s) |"; l="$l---|"; done
	echo "$h Julia kernel, 1 thread | Output identical |"; echo "$l---|---|"
} > "$MD"

for st in $STRATA; do
	loci=$(echo "$st" | sed -e 's/^QT0*//' -e 's/-.*//'); rec=$(echo "$st" | sed -e 's/.*Possion//')
	ow=$(med $st original 1 6); oc=$(med $st original 1 7); g=$(med $st original 1 8)
	row="$st,$loci,$rec,$g,$ow,$oc"; cells="| $st | $g | $ow / $oc |"
	for t in $THREADS; do
		w=$(med $st fast $t 6); c=$(med $st fast $t 7); r=$(ratio $ow $w)
		row="$row,$w,$c,$r"; cells="$cells $w / $c (${r}x) |"
	done
	jw=$(med $st julia 1 6); jc=$(med $st julia 1 7); jg=$(med $st julia 1 8); jr=$(ratio $ow $jw)
	nno=$(awk -F, -v s="$st" '$1==s && $9=="NO"' "$RUNS" | wc -l | tr -d ' ')
	[ "$nno" = 0 ] && same=yes || same="NO ($nno runs)"
	row="$row,$(ratio $oc $(med $st fast 1 7)),$jw,$jc,$jg,$jr,$same"
	echo "$row" >> "$CSV"
	echo "$cells $jw / $jc (${jr}x) | $same |" >> "$MD"
done

if [ "${BENCH_ABLATION:-0}" = 1 ]; then
	{
		echo
		echo "## One optimisation off at a time (one thread)"
		echo
		echo "CPU seconds, median over seeds, and the ratio to QTSimFast as built. \`lanes1\`: genotypes summed one individual at a time; \`dense\`: one block per chromosome (no sharing below a chromosome); \`gc8\`: unused blocks collected every 8 generations instead of every generation; \`plain\`: GSL's generator through its function pointers, no vector extensions, no threads."
		echo
		echo "| Stratum | Original | QTSimFast | lanes1 | dense | gc8 | plain |"
		echo "|---|---|---|---|---|---|---|"
		for st in $STRATA; do
			fc=$(med $st fast 1 7); line="| $st | $(med $st original 1 7) | $fc |"
			for a in lanes1 dense gc8 plain; do c=$(med $st fast-$a 1 7); line="$line $c ($(ratio $c $fc)x) |"; done
			echo "$line"
		done
	} >> "$MD"
fi

cat "$MD"
if grep -q ',NO$' "$RUNS"; then echo "SOME OUTPUTS DIFFER FROM THE ORIGINAL'S" >&2; exit 1; fi
exit 0
