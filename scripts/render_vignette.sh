#!/bin/sh
# Render vignettes to self-contained HTML and then to GitHub-flavoured Markdown.
#
#   scripts/render_vignette.sh 03_moving_optimum     # one vignette
#   scripts/render_vignette.sh                       # all of them
#
# Each vignette is rendered in two passes, HTML first: when both formats are rendered in
# one call, the self-contained HTML step can delete the figure files that the Markdown
# needs.
set -e
cd "$(dirname "$0")/../vignettes"
names="$*"
[ -z "$names" ] && names=$(ls -d [0-9][0-9]_*/ | tr -d /)
for n in $names; do
    f="$n/$n.qmd"
    echo "== $n"
    quarto render "$f" --to html 2>&1 | grep -E 'Output created|ERROR|Error' || true
    quarto render "$f" --to gfm 2>&1 | grep -E 'Output created|ERROR|Error' || true
done
