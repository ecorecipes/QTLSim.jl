#!/usr/bin/env bash
# Parity test of QTSimFast against the original QTSim: the two programs' standard
# output must be the same bytes, and their exit codes equal, for the same options
# file. Any difference is a failure.
#
#     ./test_parity.sh quick     a small matrix, about a minute (CI)
#     ./test_parity.sh full      the whole matrix, roughly 1.5 CPU-hours
#
# full also needs the variants (make variants) and, for the published inputs, the
# paper's Simulations directory (QTSIM_SIMULATIONS, default ../../../Repository/Simulations);
# without it, synthetic inputs of the same design are used.
#
# Environment: JOBS (processes at once for single-thread runs, default 4), QTSIM
# (the original, default ../QTSim/build/QTSim), QTSIM_PARITY_TMP (work directory,
# default build/parity; its inputs and outputs are kept if anything fails or KEEP=1).
#
# Variants of a run of the fast program (letters in the result lines):
#   t<n>   QTSIM_THREADS=n                 p<n>  the same, threads forced on small inputs
#   b<n>   QTSIM_BLOCK=n (block size)      g<n>  QTSIM_GC=n (collect every n generations)
#   l1     QTSIM_LANES=1                   s<n>  n threads, output on the main thread
#   auto   QTSIM_THREADS=auto              r<k>  a repeat of p4
#   plain  build without the built-in generator, threads and vector extensions
#   asan   build with AddressSanitizer and UndefinedBehaviorSanitizer (3 threads)
#   tsan   build with ThreadSanitizer (4 threads)
#   osan   the ORIGINAL built with the sanitizers: reported, never a failure

set -u
cd "$(dirname "$0")"
HERE=$(pwd)

MODE=${1:-quick}
ORIG=${QTSIM:-$HERE/../QTSim/build/QTSim}
FAST=$HERE/build/QTSimFast
WORK=${QTSIM_PARITY_TMP:-$HERE/build/parity}
SIMS=${QTSIM_SIMULATIONS:-$HERE/../../../Repository/Simulations}
JOBS=${JOBS:-4}

# ---------------------------------------------------------------- one options file
if [ "$MODE" = "__one" ]; then
	f=$2; variants=$3
	ref=$f.ref
	if [ ! -f "$ref.rc" ]; then
		"$ORIG" "$f" > "$ref" 2> /dev/null
		echo $? > "$ref.rc"
	fi
	rc=$(cat "$ref.rc")
	for v in $(echo "$variants" | tr ',' ' '); do
		bin=$FAST; e="QTSIM_THREADS=1"
		case $v in
			t*)    e="QTSIM_THREADS=${v#t}" ;;
			p*)    e="QTSIM_THREADS=${v#p} QTSIM_PAR_MIN=0" ;;
			r*)    e="QTSIM_THREADS=4 QTSIM_PAR_MIN=0" ;;
			s*)    e="QTSIM_THREADS=${v#s} QTSIM_PAR_MIN=0 QTSIM_ASYNC_OUTPUT=0" ;;
			b*)    e="QTSIM_BLOCK=${v#b}" ;;
			g*)    e="QTSIM_GC=${v#g}" ;;
			l1)    e="QTSIM_LANES=1" ;;
			auto)  e="QTSIM_THREADS=auto" ;;
			plain) bin=$FAST-plain ;;
			asan)  bin=$FAST-asan; e="QTSIM_THREADS=3 QTSIM_PAR_MIN=0" ;;
			tsan)  bin=$FAST-tsan; e="QTSIM_THREADS=4 QTSIM_PAR_MIN=0" ;;
			osan)  bin=$HERE/build/QTSim-asan ;;
			*)     echo "FAIL $f $v unknown-variant" >> "$WORK/results"; continue ;;
		esac
		out=$f.$v.out
		env $e "$bin" "$f" > "$out" 2> "$out.err"
		rc2=$?
		san=0
		grep -q "Sanitizer\|runtime error" "$out.err" && san=1
		if [ "$v" = "osan" ]; then
			if [ $san = 1 ]; then
				echo "OSAN $f $(grep -m1 "Sanitizer\|runtime error" "$out.err")" >> "$WORK/results"
			else
				echo "OSANCLEAN $f" >> "$WORK/results"
			fi
			continue
		fi
		if cmp -s "$ref" "$out" && [ "$rc" = "$rc2" ] && [ $san = 0 ]; then
			echo "OK $f $v" >> "$WORK/results"
		else
			echo "FAIL $f $v (exit $rc vs $rc2, sanitizer $san)" >> "$WORK/results"
		fi
	done
	exit 0
fi

case $MODE in quick|full) ;; *) echo "usage: $0 quick|full"; exit 2 ;; esac
[ -x "$ORIG" ] || { echo "no original at $ORIG (make -C ../QTSim)"; exit 2; }
[ -x "$FAST" ] || { echo "no $FAST (make)"; exit 2; }
if [ "$MODE" = full ]; then
	for b in "$FAST-plain" "$FAST-asan" "$FAST-tsan" "$HERE/build/QTSim-asan"; do
		[ -x "$b" ] || { echo "no $b (make variants)"; exit 2; }
	done
fi

# A fresh work directory (only ever this script's own).
case $WORK in */parity|*/parity-*) ;; *) echo "QTSIM_PARITY_TMP must end in /parity or /parity-<name>"; exit 2 ;; esac
[ -d "$WORK" ] && find "$WORK" -type f -delete
mkdir -p "$WORK/syn" "$WORK/pub"
: > "$WORK/results"; : > "$WORK/list.serial"; : > "$WORK/list.threads"

# add <file> <single-thread variants> <threaded variants>
add() {
	[ -n "$2" ] && echo "$1 $2" >> "$WORK/list.serial"
	[ -n "${3:-}" ] && echo "$1 $3" >> "$WORK/list.threads"
	return 0
}

if [ "$MODE" = full ]; then
	SV="t1,b1,b5,b64,l1,g50,plain,asan,osan"; TV="p2,p3,p4,p8,s2,r1,r2,auto,tsan"
else
	SV="t1,b5"; TV="p3"
fi

# ------------------------------------------------- (c) the synthetic option matrix
# syn <name> <option lines...>: defaults first, then the case's own lines (QTSim
# takes the last value given for an option).
syn() {
	name=$1; shift
	f=$WORK/syn/$name.txt
	{
		printf 'Sample 1\nPopSize 200\nGenerations 400\nPhenotypeVar 0.02\nHeritability 0.5\nMutationScalar 0.001\nSeed 4711\n'
		for l in "$@"; do printf '%s\n' "$l"; done
	} > "$f"
	add "$f" "$SV" "$TV"
}

if [ "$MODE" = full ]; then LS="1 2 3 7 100"; RS="0 0.7 3 12 400"; else LS="1 3 100"; RS="0 3 400"; fi
for P in 1 2; do for L in $LS; do for R in $RS; do
	syn "m-none-P$P-L$L-R$R" "Ploidy $P" "NoQTL $L" "RecombinationPoissonMean $R"
	syn "m-fit-P$P-L$L-R$R" "Ploidy $P" "NoQTL $L" "RecombinationPoissonMean $R" "MoveFitnessSD 200 0.3" "OptimumMoveType fitness"
	syn "m-phen-P$P-L$L-R$R" "Ploidy $P" "NoQTL $L" "RecombinationPoissonMean $R" "MoveFitnessSD 200 0.5" "OptimumMoveType phenotype"
	syn "m-cont-P$P-L$L-R$R" "Ploidy $P" "NoQTL $L" "RecombinationPoissonMean $R" "MoveFitnessSD 200 0.004" "OptimumMoveType continuous"
done; done; done

for P in 1 2; do
	q="Ploidy $P"
	# tiny populations, populations that die
	for N in 1 2 3 5 10; do
		syn "die-P$P-N$N" "$q" "PopSize $N" "NoQTL 5" "InitFitnessSD 0.05"
		syn "tiny-P$P-N$N" "$q" "PopSize $N" "NoQTL 5" "InitFitnessSD 5"
	done
	syn "die-move-P$P" "$q" "NoQTL 20" "MoveFitnessSD 100 0.05" "OptimumMoveType continuous"
	syn "die-jump-P$P" "$q" "NoQTL 20" "MoveFitnessSD 100 6" "OptimumMoveType phenotype"
	# odd sizes; numbers of loci around the block boundaries
	for N in 7 8 9 15 16 17 33 1001; do syn "odd-P$P-N$N" "$q" "PopSize $N" "NoQTL 37" "Generations 150" "InitFitnessSD 3"; done
	for L in 15 16 17 31 32 33 48 63 64 65 399 400 401; do syn "loci-P$P-L$L" "$q" "NoQTL $L" "PopSize 60" "Generations 200" "RecombinationPoissonMean 3" "MutationRatePerGamete 0.01"; done
	# Sample, ChangeFitnessSD, StandardisePhenotype, Drift
	for S in 2 7 1000 100000; do syn "sample-P$P-S$S" "$q" "Sample $S" "NoQTL 4" "Generations 2500" "PopSize 50" "MoveFitnessSD 2000 0.5"; done
	syn "change-before-P$P" "$q" "NoQTL 9" "ChangeFitnessSD 50 0.3" "MoveFitnessSD 200 0.5" "Sample 3"
	syn "change-at-P$P" "$q" "NoQTL 9" "ChangeFitnessSD 200 0.17" "MoveFitnessSD 200 2.0"
	syn "change-after-P$P" "$q" "NoQTL 9" "ChangeFitnessSD 300 2.5" "MoveFitnessSD 200 0.002" "OptimumMoveType continuous"
	syn "change-0-P$P" "$q" "NoQTL 9" "ChangeFitnessSD 0 0.5"
	syn "change-only-P$P" "$q" "NoQTL 9" "ChangeFitnessSD 100 0.4"
	syn "std-none-P$P" "$q" "NoQTL 9" "StandardisePhenotype" "InitFitnessSD 3"
	syn "std-fit-P$P" "$q" "NoQTL 9" "StandardisePhenotype" "InitFitnessSD 3" "MoveFitnessSD 200 0.3" "OptimumMoveType fitness"
	syn "std-phen-P$P" "$q" "NoQTL 9" "StandardisePhenotype" "InitFitnessSD 3" "MoveFitnessSD 200 0.5"
	syn "std-cont-P$P" "$q" "NoQTL 9" "StandardisePhenotype" "InitFitnessSD 3" "MoveFitnessSD 200 0.01" "OptimumMoveType continuous"
	syn "std-twice-P$P" "$q" "NoQTL 9" "StandardisePhenotype" "StandardisePhenotype"
	syn "drift-P$P" "$q" "NoQTL 9" "Drift"
	# variances: no environmental variance, no genetic target, impossible values
	syn "h1-P$P" "$q" "NoQTL 9" "Heritability 1" "InitFitnessSD 1"
	syn "h1-analytic-P$P" "$q" "NoQTL 9" "Heritability 1"
	syn "h0-P$P" "$q" "NoQTL 9" "Heritability 0"
	syn "h0-sd-P$P" "$q" "NoQTL 9" "Heritability 0" "InitFitnessSD 1"
	syn "h12-P$P" "$q" "NoQTL 9" "Heritability 1.2" "InitFitnessSD 1"
	syn "sd0-P$P" "$q" "NoQTL 9" "InitFitnessSD 0"
	syn "sdbig-P$P" "$q" "NoQTL 9" "InitFitnessSD 1e6" "MoveFitnessSD 100 0.5" "OptimumMoveType fitness"
	syn "pv0-P$P" "$q" "NoQTL 9" "PhenotypeVar 0" "InitFitnessSD 1"
	# mutation: none, rare, common, and means above 10 (GSL's other Poisson method)
	syn "mut0-P$P" "$q" "NoQTL 9" "MutationRatePerGamete 0"
	syn "mscal0-P$P" "$q" "NoQTL 9" "MutationScalar 0" "InitFitnessSD 1"
	syn "mut5-P$P" "$q" "NoQTL 100" "MutationRatePerGamete 0.05" "PopSize 80"
	syn "mut10-P$P" "$q" "NoQTL 100" "MutationRatePerGamete 0.1" "PopSize 80"
	syn "mut20-P$P" "$q" "NoQTL 100" "MutationRatePerGamete 0.2" "PopSize 80" "Generations 200"
	syn "mut50-P$P" "$q" "NoQTL 1000" "MutationRatePerGamete 0.05" "PopSize 40" "Generations 100"
	syn "mut1locus-P$P" "$q" "NoQTL 1" "MutationRatePerGamete 3"
	# recombination: the mean 10 boundary, more points than loci, many points
	syn "rec10-P$P" "$q" "NoQTL 50" "RecombinationPoissonMean 10"
	syn "rec10.5-P$P" "$q" "NoQTL 50" "RecombinationPoissonMean 10.5"
	syn "rec25-P$P" "$q" "NoQTL 30" "RecombinationPoissonMean 25"
	syn "rec99-P$P" "$q" "NoQTL 100" "RecombinationPoissonMean 99" "PopSize 60"
	syn "rec40-P$P" "$q" "NoQTL 1000" "RecombinationPoissonMean 40" "PopSize 40" "Generations 100" "MutationRatePerGamete 0.002"
	syn "rec0.01-P$P" "$q" "NoQTL 50" "RecombinationPoissonMean 0.01"
	# seeds
	for s in 0 1 4294967295 4294967296 4294967297 25106888622.8; do syn "seed-P$P-$s" "$q" "NoQTL 6" "Generations 60" "Seed $s"; done
	# nothing to simulate: no loci, no generations; a negative Sample
	syn "zero-loci-P$P" "$q" "NoQTL 0" "MoveFitnessSD 200 0.5"
	syn "zero-loci-rec0-P$P" "$q" "NoQTL 0" "RecombinationPoissonMean 0" "Heritability 1"
	syn "zero-gens-P$P" "$q" "Generations 0"
	syn "one-gen-P$P" "$q" "Generations 1"
	syn "sample-neg-P$P" "$q" "Sample -3" "NoQTL 4"
	# many loci, small population; block collection under heavy mutation
	syn "big-P$P" "$q" "NoQTL 5000" "PopSize 64" "Generations 150" "RecombinationPoissonMean 5" "MutationRatePerGamete 0.001" "MoveFitnessSD 100 0.002" "OptimumMoveType continuous"
	syn "churn-P$P" "$q" "NoQTL 257" "PopSize 48" "Generations 300" "RecombinationPoissonMean 5" "MutationRatePerGamete 0.3"
done

# option-file syntax and errors (exit codes are compared too)
f=$WORK/syn/syntax.txt
printf '# a comment\r\n\r\nsample\t2\r\nPOPSIZE 90\r\nnoqtl   11\r\n  ploidy 2\r\nGENERATIONS 120\r\nphenotypevar 0.03\r\nHeritability 0.4\r\nMutationScalar 0.002\r\nseed 77\r\nmovefitnesssd 60 0.4\r\noptimummovetype PHENOTYPE\r\n#Seed 78\r\n' > "$f"
add "$f" "$SV" "$TV"
err() { name=$1; shift; f=$WORK/syn/err-$name.txt; printf 'PopSize 20\nGenerations 5\nSeed 3\n' > "$f"; for l in "$@"; do printf '%s\n' "$l" >> "$f"; done; add "$f" "t1" ""; }
err unknown "Bogus 1"
err ploidy3 "Ploidy 3"
err args "PopSize 20 30"
err noint "NoQTL abc"
err nofloat "Heritability x"
err negrate "MutationRatePerGamete -1"
err negrec "RecombinationPoissonMean -1"
err negsd "InitFitnessSD -1"
err help "Help"
err movetype "OptimumMoveType sideways" "MoveFitnessSD 2 0.1"
err assort "AssortativeMating 1"
err zeropop1 "PopSize 0" "Ploidy 1"
err zeropop2 "PopSize 0" "Ploidy 2" "StandardisePhenotype" "MoveFitnessSD 0 0.5"
err zeropop0 "PopSize 0" "Generations 0"
add "$WORK/syn/does-not-exist.txt" "t1" ""

# ------------------------------------------------- (b) the single-locus protocol
one_locus() {
	f=$WORK/syn/onelocus-$1.txt
	printf 'Sample 1\nPopSize %s\nNoQTL 1\nGenerations 6500\nPloidy 2\nChangeFitnessSD 5000\t%s\nMoveFitnessSD 5000 2.0\nOptimumMoveType phenotype\nPhenotypeVar 0.01\nMutationScalar 0.001\nHeritability 0.5\nseed %s\n' "$2" "$3" "$4" > "$f"
	add "$f" "$5" "$6"
}
if [ "$MODE" = full ]; then
	one_locus a 1000 0.170000 101 "t1,b1,plain" "p2,p4,p8"
	one_locus b 1000 0.100000 102 "t1" "p3"
	one_locus c 1000 0.250000 103 "t1" "t2"
	one_locus d 1000 0.170000 25106888622.8 "t1" "t4"
	# QTSim's demonstration input, if the paper's code is at hand
	ex=$SIMS/../Code/QTSim/Demo/Example1.txt
	if [ -f "$ex" ]; then cp "$ex" "$WORK/syn/Example1.txt"; add "$WORK/syn/Example1.txt" "t1" "t2,t8"; fi
else
	one_locus a 300 0.170000 101 "t1" "p3"
fi

# ------------------------------------------------- (a) the published inputs
# pub <stratum> <source file> <seed or "own"> <generations> <move generation or "own"> <serial> <threaded>
pub() {
	name=$(basename "$2" .txt)-s$3-g$4
	f=$WORK/pub/$1-$name.txt
	sed -e "s/^Generations .*/Generations $4/" "$2" > "$f"
	if [ "$3" != own ]; then sed -e "s/^Seed .*/Seed $3/" "$f" > "$f.tmp"; mv "$f.tmp" "$f"; fi
	if [ "$5" != own ]; then sed -e "s/^MoveFitnessSD [0-9]* /MoveFitnessSD $5 /" "$f" > "$f.tmp"; mv "$f.tmp" "$f"; fi
	add "$f" "$6" "$7"
}
# An input of the published design (CreateCGenomeSimInputDist's sweep) when the paper's are not at hand.
# synthetic_pub <NoQTL> <rate> <seed> <recombination> <generations> <move generation> <file>
synthetic_pub() {
	printf 'Sample 1\nPopSize 1000\nNoQTL %s\nGenerations %s\nPloidy 2\nMoveFitnessSD %s %s\nPhenotypeVar 0.02\nMutationScalar 0.001\nHeritability 0.6\nSeed %s\nInitFitnessSD 1.0\nOptimumMoveType continuous\nMutationRatePerGamete 0.0001\nRecombinationPoissonMean %s\n' "$1" "$5" "$6" "$2" "$3" "$4" > "$7"
}

if [ "$MODE" = full ] && [ -d "$SIMS/QT0100-Possion3.0" ]; then
	for st in QT0100-Possion3.0 QT0100-Possion5.0 QT0500-Possion3.0 QT0500-Possion5.0 QT1000-Possion3.0 QT1000-Possion5.0 QT5000-Possion3.0 QT5000-Possion5.0; do
		d=$SIMS/$st
		a=$d/Input-00001.txt; b=$d/Input-00002.txt
		# the slowest-moving optimum among the first 400 inputs: a population that survives
		c=$(for i in $(ls "$d" | head -400); do echo "$(awk '/^MoveFitnessSD/{print $3}' "$d/$i") $d/$i"; done | sort -g | head -1 | cut -d' ' -f2)
		case $st in
		QT5000*)
			# two full-length runs per stratum, the rest with the move brought forward
			pub $st "$a" own 12000 own "t1" "t4"
			pub $st "$c" own 12000 own "t1" "t3"
			for s in 1001 1002; do
				pub $st "$a" $s 1600 1000 "t1,b16" "t2"
				pub $st "$b" $s 1600 1000 "t1" "t8"
				pub $st "$c" $s 1600 1000 "t1,l1" "t4"
			done
			pub $st "$b" own 1600 1000 "t1" "t3"
			;;
		*)
			for s in own 1001 1002; do
				pub $st "$a" $s 12000 own "t1" "p2,t4"
				pub $st "$b" $s 12000 own "t1" "p3"
				pub $st "$c" $s 12000 own "t1" "p8"
			done
			;;
		esac
	done
elif [ "$MODE" = full ]; then
	echo "note: no published inputs at $SIMS; using synthetic inputs of the same design"
	for L in 100 500 1000; do for R in 3 5; do for k in 0.2 0.002; do
		f=$WORK/pub/synpub-L$L-R$R-k$k.txt; synthetic_pub $L $k 1001 $R 12000 10000 "$f"; add "$f" "t1" "t2,t4"
	done; done; done
	for R in 3 5; do for k in 0.2 0.002; do
		f=$WORK/pub/synpub-L5000-R$R-k$k.txt; synthetic_pub 5000 $k 1001 $R 1600 1000 "$f"; add "$f" "t1" "t3"
	done; done
else
	for L in 100 500; do
		f=$WORK/pub/synpub-L$L.txt; synthetic_pub $L 0.05 1001 3 700 500 "$f"; add "$f" "t1" "t2"
	done
fi

# ------------------------------------------------- run
echo "parity ($MODE): $(wc -l < "$WORK/list.serial" | tr -d ' ') option files, original $ORIG"
# single-thread variants, JOBS at once; then the threaded variants one at a time
tr ' ' '\n' < "$WORK/list.serial" | xargs -n 2 -P "$JOBS" "$HERE/test_parity.sh" __one
tr ' ' '\n' < "$WORK/list.threads" | xargs -n 2 -P 1 "$HERE/test_parity.sh" __one

# with no argument both programs print the name of their main source file
"$ORIG" > "$WORK/noarg.ref" 2> /dev/null; rc=$?
"$FAST" > "$WORK/noarg.out" 2> /dev/null; rc2=$?
if cmp -s "$WORK/noarg.ref" "$WORK/noarg.out" && [ $rc = $rc2 ]; then echo "OK noarg t1" >> "$WORK/results"; else echo "FAIL noarg t1" >> "$WORK/results"; fi

# ------------------------------------------------- report
ok=$(grep -c '^OK ' "$WORK/results"); fail=$(grep -c '^FAIL ' "$WORK/results")
died=$(cat "$WORK"/syn/*.ref "$WORK"/pub/*.ref 2> /dev/null | grep -c '^Pop all died')
nref=$(ls "$WORK"/syn/*.ref "$WORK"/pub/*.ref 2> /dev/null | wc -l | tr -d ' ')
bytes=$(cat "$WORK"/syn/*.ref "$WORK"/pub/*.ref 2> /dev/null | wc -c | tr -d ' ')
echo "option files run by the original: $nref ($died ended with 'Pop all died.'), $bytes bytes of output"
echo "comparisons: $ok identical, $fail different"
echo "by variant (trailing digits dropped):"
grep '^OK \|^FAIL ' "$WORK/results" | awk '{v=$3; sub(/[0-9]+$/,"",v); n[v" "$1]++} END{for(k in n) print "  " k, n[k]}' | sort
if [ "$MODE" = full ]; then
	echo "published inputs: $(ls "$WORK"/pub/*.ref 2> /dev/null | wc -l | tr -d ' ') files, $(cat "$WORK"/pub/*.ref 2> /dev/null | grep -c '^Pop all died') died"
	echo "original under ASan/UBSan: $(grep -c '^OSANCLEAN ' "$WORK/results") clean, $(grep -c '^OSAN ' "$WORK/results") with reports"
	grep '^OSAN ' "$WORK/results" | sed -e 's/^OSAN [^ ]* //' | sort | uniq -c | sed 's/^/  /'
fi
if [ "$fail" != 0 ]; then
	grep '^FAIL ' "$WORK/results" | head -50
	echo "PARITY FAILED (files kept in $WORK)"
	exit 1
fi
echo "PARITY OK"
if [ "${KEEP:-0}" != 1 ]; then find "$WORK/syn" "$WORK/pub" -type f -delete; fi
exit 0
