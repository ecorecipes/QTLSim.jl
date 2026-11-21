#!/usr/bin/env bash
# Assemble the hand-over package from this directory and ../QTSim:
#
#     ./make_package.sh [destination directory]
#
# Default destination: dist/QTSimFast-<date> beside the QTLSim.jl repository's parent
# (PACKAGE_DATE overrides the date). Writes the directory and <directory>.tar.gz.
# Nothing here is built or run; the package's own Makefile does that. The paper's input
# files, results, manual and scripts are not copied (no licence statement).

set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
DATE=${PACKAGE_DATE:-$(date +%Y-%m-%d)}
DEST=${1:-$HERE/../../../dist/QTSimFast-$DATE}
mkdir -p "$DEST"
DEST=$(cd "$DEST" && pwd)
NAME=$(basename "$DEST")

mkdir -p "$DEST/original/src" "$DEST/fast/src" "$DEST/examples" "$DEST/formal"
cp "$HERE"/../QTSim/src/*.c "$HERE"/../QTSim/src/*.h "$DEST/original/src/"
cp "$HERE/../QTSim/SHA256SUMS" "$DEST/original/SHA256SUMS"
cp "$HERE/package/original-README.md" "$DEST/original/README.md"
cp "$HERE"/src/*.c "$HERE"/src/*.h "$DEST/fast/src/"
cp "$HERE"/examples/*.txt "$DEST/examples/"
cp "$HERE/package/Makefile" "$DEST/Makefile"
cp "$HERE/package/README.md" "$DEST/README.md"
cp "$HERE/../../LICENSE" "$DEST/COPYING"
for f in demo_parity.sh test_parity.sh bench.sh bench_results.md bench_results.csv bench_runs.csv bench_fast.sh bench_fast.md bench_fast_runs.csv FAST.md REPORT.md PROOFS.md PROOFS.html PROOFS.pdf PARITY_DEMO.txt PARITY_FULL.txt; do
	[ -f "$HERE/$f" ] && cp "$HERE/$f" "$DEST/$f" || echo "note: $f is missing"
done
for f in README.md QTSimFast.lean QTSimFastShadows.lean ModelCheck.lean model_check.c model_check_expected.txt check.sh render.sh lean2md.awk QTSimFast.md QTSimFast.html QTSimFast.pdf QTSimFastShadows.md QTSimFastShadows.html QTSimFastShadows.pdf; do
	[ -f "$HERE/formal/$f" ] && cp "$HERE/formal/$f" "$DEST/formal/$f" || echo "note: formal/$f is missing"
done
# the fast modes' statistical validation (FAST.md), from QTLSim.jl's validate/
mkdir -p "$DEST/validation"
for f in qtsimfast_fast_criteria.md qtsimfast_fast.jl qtsimfast_fast.log qtsimfast_fast_run2.log; do
	[ -f "$HERE/../../validate/$f" ] && cp "$HERE/../../validate/$f" "$DEST/validation/$f" || echo "note: validate/$f is missing"
done
# the logs name the binary by its absolute path on the machine that ran them: keep only
# the part inside QTLSim.jl
for f in "$DEST"/validation/*.log; do
	[ -f "$f" ] && sed -i.bak -e 's|/[^ ;]*/reference/QTSimFast/|<QTLSim.jl>/reference/QTSimFast/|g' "$f" && rm -f "$f.bak"
done
for d in qtsimfast_fast_figs qtsimfast_fast_figs_run2; do
	[ -d "$HERE/../../validate/$d" ] && { mkdir -p "$DEST/validation/$d"; cp "$HERE/../../validate/$d"/*.png "$DEST/validation/$d/"; } || echo "note: validate/$d is missing"
done
chmod +x "$DEST"/*.sh "$DEST"/formal/*.sh

# the archive, without anything built
( cd "$DEST/.." && COPYFILE_DISABLE=1 tar --exclude='build' --exclude='.DS_Store' -czf "$NAME.tar.gz" "$NAME" )
echo "package: $DEST"
echo "archive: $DEST.tar.gz"
