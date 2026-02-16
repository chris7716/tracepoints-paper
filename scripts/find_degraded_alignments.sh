#!/bin/bash
# Usage: ./find_degraded_alignments.sh <original.paf> <reconstructed.paf> [output.tsv]
#
# Finds alignments where the reconstructed score is worse than the original.
# Outputs: query, target, coords, original_score, reconstructed_score, score_diff

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <original.paf> <reconstructed.paf> [output.tsv]" >&2
    exit 1
fi

orig_paf="$1"
recon_paf="$2"
output="${3:-/dev/stdout}"

# Header
echo -e "query\tqlen\tqstart\tqend\tstrand\ttarget\ttlen\ttstart\ttend\torig_score\trecon_score\tscore_diff" > "$output"

# Join sorted files, compare scores, output degraded
paste \
    <(sort -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$orig_paf" | awk -F'\t' -v OFS='\t' '{
        score=""
        for (i = 13; i <= NF; i++) if ($i ~ /^sc:i:/) { score = substr($i, 6) + 0; break }
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, score
    }') \
    <(sort -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$recon_paf" | awk -F'\t' '{
        for (i = 13; i <= NF; i++) if ($i ~ /^sc:i:/) { print substr($i, 6) + 0; next }
        print ""
    }') | \
awk -F'\t' -v OFS='\t' '
    $10 != "" && $11 != "" && $11 < $10 {
        # Degraded: recon_score < orig_score (scores are negative, lower = worse)
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $11 - $10
    }
' >> "$output"

count=$(tail -n +2 "$output" | wc -l)
echo "Found $count degraded alignments" >&2
