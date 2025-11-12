#!/usr/bin/env python3
"""
Aggregate contamination stats from PhiX and UniVec flagging

Parses BBDuk stats files to extract contamination percentages and
creates a summary table for identifying outlier samples.
"""

import pandas as pd
import re
from pathlib import Path


def parse_bbduk_stats(stats_file):
    """
    Parse BBDuk stats file to extract contamination metrics

    BBDuk stats format:
    #File	reads	bases	...
    Input:	10000	1500000	...
    Match:	150	22500	...

    Returns:
        dict: {'input_reads': int, 'match_reads': int, 'match_pct': float}
    """
    input_reads = 0
    match_reads = 0

    with open(stats_file) as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Parse Input line
            if line.startswith('Input:'):
                parts = line.split('\t')
                try:
                    input_reads = int(parts[1])
                except (ValueError, IndexError):
                    pass

            # Parse Match line (reads matching contaminants)
            elif line.startswith('Match:') or line.startswith('Matched:'):
                parts = line.split('\t')
                try:
                    match_reads = int(parts[1])
                except (ValueError, IndexError):
                    pass

    # Calculate percentage
    match_pct = (match_reads / input_reads * 100) if input_reads > 0 else 0.0

    return {
        'input_reads': input_reads,
        'match_reads': match_reads,
        'match_pct': match_pct
    }


def main():
    # Get inputs from snakemake
    phix_stats_files = snakemake.input.phix_stats
    univec_stats_files = snakemake.input.univec_stats
    output_file = snakemake.output[0]

    # Extract sample names and process
    results = []

    for phix_file, univec_file in zip(phix_stats_files, univec_stats_files):
        # Extract sample name from file path
        sample = Path(phix_file).stem.replace('_phix_stats', '')

        # Parse PhiX stats
        phix_stats = parse_bbduk_stats(phix_file)

        # Parse UniVec stats
        univec_stats = parse_bbduk_stats(univec_file)

        # Compile results
        results.append({
            'sample': sample,
            'input_reads': phix_stats['input_reads'],
            'phix_reads': phix_stats['match_reads'],
            'phix_percent': phix_stats['match_pct'],
            'vector_reads': univec_stats['match_reads'],
            'vector_percent': univec_stats['match_pct'],
            'total_contamination_percent': phix_stats['match_pct'] + univec_stats['match_pct']
        })

    # Create DataFrame
    results_df = pd.DataFrame(results)

    # Sort by total contamination (highest first)
    results_df = results_df.sort_values('total_contamination_percent', ascending=False)

    # Write to file
    results_df.to_csv(output_file, sep='\t', index=False, float_format='%.3f')

    # Print summary
    print("\n" + "="*80)
    print("CONTAMINATION FLAGGING SUMMARY")
    print("="*80)
    print(f"\nTotal samples: {len(results_df)}")
    print(f"\nPhiX contamination:")
    print(f"  Mean: {results_df['phix_percent'].mean():.3f}%")
    print(f"  Median: {results_df['phix_percent'].median():.3f}%")
    print(f"  Max: {results_df['phix_percent'].max():.3f}% ({results_df.loc[results_df['phix_percent'].idxmax(), 'sample']})")
    print(f"\nVector contamination:")
    print(f"  Mean: {results_df['vector_percent'].mean():.3f}%")
    print(f"  Median: {results_df['vector_percent'].median():.3f}%")
    print(f"  Max: {results_df['vector_percent'].max():.3f}% ({results_df.loc[results_df['vector_percent'].idxmax(), 'sample']})")
    print(f"\nTotal contamination:")
    print(f"  Mean: {results_df['total_contamination_percent'].mean():.3f}%")
    print(f"  Median: {results_df['total_contamination_percent'].median():.3f}%")
    print(f"  Max: {results_df['total_contamination_percent'].max():.3f}%")

    # Flag high contamination samples (>1% for either category)
    high_phix = results_df[results_df['phix_percent'] > 1.0]
    high_vector = results_df[results_df['vector_percent'] > 1.0]

    if len(high_phix) > 0:
        print(f"\n⚠️  Samples with >1% PhiX contamination ({len(high_phix)}):")
        for _, row in high_phix.iterrows():
            print(f"  {row['sample']}: {row['phix_percent']:.3f}%")

    if len(high_vector) > 0:
        print(f"\n⚠️  Samples with >1% vector contamination ({len(high_vector)}):")
        for _, row in high_vector.iterrows():
            print(f"  {row['sample']}: {row['vector_percent']:.3f}%")

    print("="*80 + "\n")


if __name__ == "__main__":
    main()
