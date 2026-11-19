#!/usr/bin/env bash
# Render the Lean files of this directory as documents: Markdown (mdgen's format),
# a self-contained HTML page and a PDF.
#
#     formal/render.sh
#
# For QTSimFast.lean and QTSimFastShadows.lean it writes <name>.md, <name>.html and
# <name>.pdf beside them. Every rendered document states the SHA-256 of the Lean file
# it was made from; the documents are generated, never edited.
#
# Tools (each optional; what cannot be made is skipped with a note):
#   mdgen   https://github.com/Seasawher/mdgen, the executable named by MDGEN or
#           `mdgen` on the PATH. Without it, lean2md.awk here does the same
#           conversion (block comments become prose, code becomes fenced blocks)
#           and the documents say so.
#   quarto  for HTML (one self-contained file) and PDF (through typst, which Quarto
#           bundles; no TeX is needed).
# The PDF's code font is the first available of Menlo, DejaVu Sans Mono, with STIX
# Two Math and Apple Symbols for symbols those lack.

set -u
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/qtsimfast-render.XXXXXX") || exit 2
trap 'rm -rf "$WORK"' EXIT

sha() {
	if command -v shasum > /dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1
	else sha256sum "$1" | cut -d' ' -f1; fi
}

MDGEN=${MDGEN:-$(command -v mdgen 2> /dev/null || true)}
mkdir -p "$WORK/in" "$WORK/out"
cp QTSimFast.lean QTSimFastShadows.lean "$WORK/in/"
if [ -n "$MDGEN" ] && [ -x "$MDGEN" ] && "$MDGEN" "$WORK/in" "$WORK/out" > /dev/null 2>&1 && [ -s "$WORK/out/QTSimFast.md" ]; then
	converter="mdgen $("$MDGEN" --version 2> /dev/null | head -1)"
else
	echo "note: mdgen not available; using the fallback converter lean2md.awk"
	for n in QTSimFast QTSimFastShadows; do awk -f lean2md.awk "$n.lean" > "$WORK/out/$n.md"; done
	converter="lean2md.awk (a fallback that writes mdgen's format; not mdgen itself)"
fi

status=0
for n in QTSimFast QTSimFastShadows; do
	cp "$WORK/out/$n.md" "$n.md"
	h=$(sha "$n.lean")
	echo "$n.md written ($converter); $n.lean sha256 $h"
	command -v quarto > /dev/null 2>&1 || { echo "note: quarto not found; $n.html and $n.pdf not made"; continue; }
	d=$WORK/q-$n; mkdir -p "$d"
	{
		cat <<EOT
---
pagetitle: "$n.lean"
format:
  html:
    embed-resources: true
    toc: true
    toc-depth: 2
    html-math-method: mathml
    code-overflow: scroll
    theme: cosmo
    fontsize: 1em
    monofont: "Menlo, DejaVu Sans Mono, monospace"
  typst:
    toc: true
    toc-depth: 2
    papersize: a4
    margin:
      x: 1.7cm
      y: 2cm
    fontsize: 10pt
    include-in-header:
      text: |
        #show raw: set text(font: ("Menlo", "DejaVu Sans Mono", "STIX Two Math", "Apple Symbols"))
        #show raw.where(block: true): set text(size: 7.3pt)
        #show raw.where(block: true): set par(leading: 0.5em)
        #set table(inset: 5pt)
        #show table: set text(size: 8.5pt)
---

> Generated from \`$n.lean\` (SHA-256 \`$h\`) by $converter and Quarto; not edited by
> hand. The Lean file was type-checked with Lean 4.32.1 and Mathlib v4.32.1
> (\`formal/check.sh\`); the code blocks below are its code, in order and complete, and
> the text is its comments.

EOT
		cat "$n.md"
	} > "$d/$n.md"
	( cd "$d" && quarto render "$n.md" --to html > html.log 2>&1 ) && [ -s "$d/$n.html" ] && cp "$d/$n.html" "$n.html" && echo "$n.html written" || { echo "note: HTML rendering of $n failed"; tail -5 "$d/html.log"; status=1; }
	( cd "$d" && quarto render "$n.md" --to typst > pdf.log 2>&1 ) && [ -s "$d/$n.pdf" ] && cp "$d/$n.pdf" "$n.pdf" && echo "$n.pdf written" || { echo "note: PDF rendering of $n failed"; tail -5 "$d/pdf.log"; status=1; }
done
exit $status
