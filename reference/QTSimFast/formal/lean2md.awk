# A fallback for mdgen (https://github.com/Seasawher/mdgen): converts a Lean file to
# Markdown in the same format. Block comments that start a line with "/-!" or "/-"
# (not the documentation comments "/--") and end a line with "-/" become prose; the
# code between them becomes fenced lean blocks. Used by render.sh only when mdgen is
# not available.
function flush(   i, first, last) {
	first = 1; last = n
	while(first <= last && code[first] ~ /^[ \t]*$/) first++
	while(last >= first && code[last] ~ /^[ \t]*$/) last--
	if(first <= last) {
		print "```lean"
		for(i = first; i <= last; i++) print code[i]
		print "```"
		print ""
	}
	n = 0
}
BEGIN { n = 0; prose = 0 }
{
	line = $0
	if(!prose && (line ~ /^\/-!/ || (line ~ /^\/-/ && line !~ /^\/--/))) {
		flush()
		prose = 1
		sub(/^\/-!?[ ]?/, "", line)
		if(line ~ /-\/[ \t]*$/) { sub(/[ ]?-\/[ \t]*$/, "", line); prose = 0; if(line != "") print line; print ""; next }
		if(line != "") print line
		next
	}
	if(prose) {
		if(line ~ /-\/[ \t]*$/) { sub(/[ ]?-\/[ \t]*$/, "", line); prose = 0; if(line != "") print line; print ""; next }
		print line
		next
	}
	code[++n] = line
}
END { flush() }
