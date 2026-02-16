#!/usr/bin/env python3
"""
Convert alignment output to PAF format.
Assumes global (end-to-end) alignment.
D = deletion in query (gap in query, consumes target)
I = insertion in query (gap in target, consumes query)
M = match (converted to '=' in output for standard CIGAR)
X = mismatch
"""

import re
import sys

def parse_cigar(cigar):
    """
    Parse CIGAR string and return query_len, target_len, matches, aln_block_len.

    CIGAR operations:
    - M: match (consumes both query and target)
    - X: mismatch (consumes both query and target)
    - I: insertion in query (consumes query only)
    - D: deletion in query (consumes target only)
    """
    # Find all operations: number followed by operation letter
    ops = re.findall(r'(\d+)([MIDX])', cigar)

    query_len = 0
    target_len = 0
    matches = 0
    aln_block_len = 0

    for count, op in ops:
        count = int(count)
        aln_block_len += count

        if op in ['M', 'X']:  # Match/mismatch consumes both query and target
            query_len += count
            target_len += count
            if op == 'M':
                matches += count
        elif op == 'I':  # Insertion in query (consumes query)
            query_len += count
        elif op == 'D':  # Deletion in query (consumes target)
            target_len += count

    return query_len, target_len, matches, aln_block_len


def convert_to_paf(input_file, output_file):
    """Convert alignment file to PAF format."""
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        query_id = 1  # Counter starting from 1

        for line in f_in:
            line = line.strip()
            if not line:
                continue

            # Parse the line: alignment_score CIGAR
            parts = line.split('\t')
            if len(parts) != 2:
                continue

            alignment_score = int(parts[0].strip())
            cigar = parts[1].strip()

            # Replace M with = for standard CIGAR (M -> = for sequence match)
            cigar_standard = cigar.replace('M', '=')

            # Each line is a separate alignment: query_N to target_N
            target_id = query_id

            # Parse CIGAR (using original before replacement)
            query_len, target_len, matches, aln_block_len = parse_cigar(cigar)

            # Global alignment: end-to-end
            query_start = 0
            query_end = query_len
            target_start = 0
            target_end = target_len
            strand = '+'
            mapq = 255  # Mapping quality not available

            # Generate sequence names (counter starting from 1, already in the file)
            query_name = f"query_{query_id}"
            target_name = f"target_{target_id}"

            # PAF format (12 mandatory columns):
            # 1. Query name
            # 2. Query length
            # 3. Query start (0-based)
            # 4. Query end (0-based)
            # 5. Strand (+ or -)
            # 6. Target name
            # 7. Target length
            # 8. Target start (0-based)
            # 9. Target end (0-based)
            # 10. Number of matching bases
            # 11. Alignment block length
            # 12. Mapping quality
            # Optional tags: AS:i:alignment_score, cg:Z:CIGAR
            paf_line = (f"{query_name}\t{query_len}\t{query_start}\t{query_end}\t"
                       f"{strand}\t{target_name}\t{target_len}\t{target_start}\t"
                       f"{target_end}\t{matches}\t{aln_block_len}\t{mapq}\t"
                       f"AS:i:{alignment_score}\tcg:Z:{cigar_standard}\n")
            f_out.write(paf_line)

            query_id += 1  # Increment query counter


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 out_to_paf.py <input_file> <output_file>")
        print("Example: python3 out_to_paf.py set_100_001.out set_100_001.paf")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    convert_to_paf(input_file, output_file)
    print(f"Converted {input_file} to {output_file}")

