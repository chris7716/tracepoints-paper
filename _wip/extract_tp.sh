#!/bin/bash

input="$1"

# Combine continuation lines into full records
# (Because your PAF lines are split over many physical lines)
awk '
    {
        if ($0 ~ /^pattern/) {       # new PAF record starts
            if (line != "") {
                print line
            }
            line = $0
        } else {
            line = line $0
        }
    }
    END { if (line != "") print line }
' "$input" |
while read -r full; do
    # Extract substring after tp:Z:
    tp=$(echo "$full" | sed -n 's/.*tp:Z:\(.*\)$/\1/p')

    # Split by ";" into array
    IFS=';' read -ra pairs <<< "$tp"

    T_line="X"
    X_line="T"

    for p in "${pairs[@]}"; do
        IFS=',' read a b <<< "$p"
        T_line+=" $a"
        X_line+=" $b"
    done

    
    echo "$X_line"
    echo "$T_line"
    echo
done
