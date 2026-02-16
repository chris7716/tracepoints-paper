#!/usr/bin/env python3
"""
Compute statistics for the 1st and 2nd values in the tp:Z: tag of PAF files.
The tp:Z: tag format is: val1,val2;val1,val2;...

Reports: min/average/std/max for both value columns.
Output: TAB-separated table with one row per file and a final row for all files.

Processes files in parallel using multiprocessing.
"""

import argparse
import gzip
import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path


def parse_tp_values(paf_file):
    """Extract 1st and 2nd values from tp:Z: tags in a PAF file."""
    values1 = []
    values2 = []

    opener = gzip.open if str(paf_file).endswith('.gz') else open
    mode = 'rt' if str(paf_file).endswith('.gz') else 'r'

    with opener(paf_file, mode) as f:
        for line in f:
            fields = line.strip().split('\t')
            # Parse optional fields for tp:Z:
            for field in fields[12:]:
                if field.startswith('tp:Z:'):
                    tp_data = field[5:]  # Remove 'tp:Z:' prefix
                    if tp_data:
                        pairs = tp_data.split(';')
                        for pair in pairs:
                            if ',' in pair:
                                parts = pair.split(',')
                                if len(parts) >= 2:
                                    try:
                                        values1.append(int(parts[0]))
                                        values2.append(int(parts[1]))
                                    except ValueError:
                                        pass
                    break

    return values1, values2


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
        'max': max_val
    }


def compute_partial_stats(values):
    """Compute partial statistics that can be combined later."""
    if not values:
        return None

    return {
        'count': len(values),
        'sum': sum(values),
        'sum_sq': sum(x * x for x in values),
        'min': min(values),
        'max': max(values)
    }


def combine_partial_stats(partials):
    """Combine partial statistics from multiple files."""
    partials = [p for p in partials if p is not None]
    if not partials:
        return None

    total_count = sum(p['count'] for p in partials)
    total_sum = sum(p['sum'] for p in partials)
    total_sum_sq = sum(p['sum_sq'] for p in partials)
    total_min = min(p['min'] for p in partials)
    total_max = max(p['max'] for p in partials)

    mean = total_sum / total_count
    # Variance = E[X^2] - E[X]^2
    variance = (total_sum_sq / total_count) - (mean * mean)
    std_dev = math.sqrt(max(0, variance))  # max(0, ...) to handle floating point errors

    return {
        'count': total_count,
        'min': total_min,
        'avg': mean,
        'std': std_dev,
        'max': total_max
    }


def process_file(paf_file):
    """Process a single PAF file and return results."""
    path = Path(paf_file)
    if not path.exists():
        return (path.name, None, None, None, None)

    print(f"Processing: {path.name}", file=sys.stderr)
    values1, values2 = parse_tp_values(path)

    stats1 = compute_stats(values1)
    stats2 = compute_stats(values2)
    partial1 = compute_partial_stats(values1)
    partial2 = compute_partial_stats(values2)

    return (path.name, stats1, stats2, partial1, partial2)


def format_row(name, stats1, stats2):
    """Format a row for output."""
    row = [
        name,
        f"{stats1['count']}" if stats1 else '0',
        f"{stats1['min']}" if stats1 else 'NA',
        f"{stats1['avg']:.2f}" if stats1 else 'NA',
        f"{stats1['std']:.2f}" if stats1 else 'NA',
        f"{stats1['max']}" if stats1 else 'NA',
        f"{stats2['min']}" if stats2 else 'NA',
        f"{stats2['avg']:.2f}" if stats2 else 'NA',
        f"{stats2['std']:.2f}" if stats2 else 'NA',
        f"{stats2['max']}" if stats2 else 'NA'
    ]
    return '\t'.join(row)


def print_header():
    cols = [
        'file',
        'num_pairs',
        'val1_min',
        'val1_avg',
        'val1_std',
        'val1_max',
        'val2_min',
        'val2_avg',
        'val2_std',
        'val2_max'
    ]
    print('\t'.join(cols))


def main():
    parser = argparse.ArgumentParser(
        description='Compute statistics for tp:Z: tag values in PAF files'
    )
    parser.add_argument(
        'paf_files',
        nargs='+',
        help='PAF files to process (can be gzipped)'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=1,
        help='Number of parallel threads (default: 1)'
    )
    args = parser.parse_args()

    print_header()

    # Process files in parallel
    with Pool(processes=args.threads) as pool:
        results = pool.map(process_file, args.paf_files)

    # Collect partial stats for combined summary
    all_partial1 = []
    all_partial2 = []

    # Print results in order
    for name, stats1, stats2, partial1, partial2 in results:
        if stats1 is None and stats2 is None and partial1 is None and partial2 is None:
            print(f"Warning: {name} not found, skipping", file=sys.stderr)
            continue
        print(format_row(name, stats1, stats2))
        all_partial1.append(partial1)
        all_partial2.append(partial2)

    # Combined summary
    if len(all_partial1) > 1:
        combined1 = combine_partial_stats(all_partial1)
        combined2 = combine_partial_stats(all_partial2)
        print(format_row('TOTAL', combined1, combined2))


if __name__ == '__main__':
    main()
