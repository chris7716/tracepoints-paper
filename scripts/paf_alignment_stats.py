#!/usr/bin/env python3
"""
Compute alignment statistics from PAF files:
- Total alignment length (in GB)
- Min/average/std/max alignment length
- Min/average/std/max identity (gi:f: field)
"""

import argparse
import gzip
import math
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path


def parse_paf_stats(args):
    """Extract alignment lengths and identities from a PAF file."""
    paf_file, min_identity = args
    lengths = []
    identities = []

    opener = gzip.open if str(paf_file).endswith('.gz') else open
    mode = 'rt' if str(paf_file).endswith('.gz') else 'r'

    print(f"Processing: {Path(paf_file).name}", file=sys.stderr)

    with opener(paf_file, mode) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 11:
                aln_len = int(fields[10])
                lengths.append(aln_len)

                # Parse optional fields for gi:f: (identity)
                for field in fields[12:]:
                    if field.startswith('gi:f:'):
                        identity = float(field[5:])
                        if identity >= min_identity:
                            identities.append(identity)
                        break

    return paf_file, lengths, identities


def compute_stats(values):
    """Compute min, average, standard deviation, and max."""
    if not values:
        return None

    n = len(values)
    min_val = min(values)
    max_val = max(values)
    mean = sum(values) / n

    # Standard deviation (population)
    variance = sum((x - mean) ** 2 for x in values) / n
    std_dev = math.sqrt(variance)

    return {
        'count': n,
        'min': min_val,
        'avg': mean,
        'std': std_dev,
        'max': max_val,
        'sum': sum(values)
    }


def print_header():
    cols = [
        'file',
        'num_alignments',
        'total_length_GB',
        'len_min',
        'len_avg',
        'len_std',
        'len_max',
        'num_identities',
        'id_min',
        'id_avg',
        'id_std',
        'id_max'
    ]
    print('\t'.join(cols))


def print_row(name, len_stats, id_stats):
    total_gb = len_stats['sum'] / 1e9 if len_stats else 0

    row = [
        name,
        f"{len_stats['count']}" if len_stats else '0',
        f"{total_gb:.3f}",
        f"{len_stats['min']}" if len_stats else 'NA',
        f"{len_stats['avg']:.2f}" if len_stats else 'NA',
        f"{len_stats['std']:.2f}" if len_stats else 'NA',
        f"{len_stats['max']}" if len_stats else 'NA',
        f"{id_stats['count']}" if id_stats else '0',
        f"{id_stats['min']:.6f}" if id_stats else 'NA',
        f"{id_stats['avg']:.6f}" if id_stats else 'NA',
        f"{id_stats['std']:.6f}" if id_stats else 'NA',
        f"{id_stats['max']:.6f}" if id_stats else 'NA'
    ]
    print('\t'.join(row))


def main():
    parser = argparse.ArgumentParser(
        description='Compute alignment length and identity statistics from PAF files'
    )
    parser.add_argument(
        'paf_files',
        nargs='+',
        help='PAF files to process (can be gzipped)'
    )
    parser.add_argument(
        '--min-identity',
        type=float,
        default=0.1,
        help='Minimum identity threshold for identity stats (default: 0.1)'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=cpu_count(),
        help=f'Number of parallel threads (default: {cpu_count()})'
    )
    args = parser.parse_args()

    # Filter valid files
    valid_files = []
    for paf_file in args.paf_files:
        path = Path(paf_file)
        if not path.exists():
            print(f"Warning: {paf_file} not found, skipping", file=sys.stderr)
        else:
            valid_files.append(str(path))

    if not valid_files:
        print("No valid files to process", file=sys.stderr)
        sys.exit(1)

    # Process files in parallel
    work_args = [(f, args.min_identity) for f in valid_files]

    with Pool(processes=min(args.threads, len(valid_files))) as pool:
        results = pool.map(parse_paf_stats, work_args)

    # Sort results by original file order
    file_order = {f: i for i, f in enumerate(valid_files)}
    results.sort(key=lambda x: file_order[x[0]])

    all_lengths = []
    all_identities = []

    print_header()

    for paf_file, lengths, identities in results:
        all_lengths.extend(lengths)
        all_identities.extend(identities)

        len_stats = compute_stats(lengths)
        id_stats = compute_stats(identities)
        print_row(Path(paf_file).name, len_stats, id_stats)

    # Combined summary
    if all_lengths and len(valid_files) > 1:
        len_stats = compute_stats(all_lengths)
        id_stats = compute_stats(all_identities)
        print_row('TOTAL', len_stats, id_stats)


if __name__ == '__main__':
    main()
