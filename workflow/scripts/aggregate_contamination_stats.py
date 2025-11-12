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

    # Calculate statistics for outlier detection
    import numpy as np

    def identify_outliers(series, method='iqr'):
        """
        Identify statistical outliers using IQR method
        Returns boolean mask of outliers
        """
        if len(series) < 4:
            # Not enough samples for meaningful outlier detection
            return pd.Series([False] * len(series), index=series.index)

        Q1 = series.quantile(0.25)
        Q3 = series.quantile(0.75)
        IQR = Q3 - Q1

        # Outliers are values beyond 1.5 * IQR from Q1/Q3
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        return (series < lower_bound) | (series > upper_bound)

    # Identify outliers for each contamination type
    results_df['phix_outlier'] = identify_outliers(results_df['phix_percent'])
    results_df['vector_outlier'] = identify_outliers(results_df['vector_percent'])
    results_df['total_outlier'] = identify_outliers(results_df['total_contamination_percent'])

    # Print summary statistics
    print("\n" + "="*80)
    print("CONTAMINATION SUMMARY - Identifying Outliers Within This Run")
    print("="*80)
    print(f"\nTotal samples: {len(results_df)}")

    print(f"\n{'PhiX contamination:':<30}")
    print(f"  {'Median:':<15} {results_df['phix_percent'].median():.3f}%")
    print(f"  {'Mean ± StdDev:':<15} {results_df['phix_percent'].mean():.3f}% ± {results_df['phix_percent'].std():.3f}%")
    print(f"  {'Range:':<15} {results_df['phix_percent'].min():.3f}% - {results_df['phix_percent'].max():.3f}%")
    print(f"  {'25-75th %ile:':<15} {results_df['phix_percent'].quantile(0.25):.3f}% - {results_df['phix_percent'].quantile(0.75):.3f}%")

    print(f"\n{'Vector contamination:':<30}")
    print(f"  {'Median:':<15} {results_df['vector_percent'].median():.3f}%")
    print(f"  {'Mean ± StdDev:':<15} {results_df['vector_percent'].mean():.3f}% ± {results_df['vector_percent'].std():.3f}%")
    print(f"  {'Range:':<15} {results_df['vector_percent'].min():.3f}% - {results_df['vector_percent'].max():.3f}%")
    print(f"  {'25-75th %ile:':<15} {results_df['vector_percent'].quantile(0.25):.3f}% - {results_df['vector_percent'].quantile(0.75):.3f}%")

    # Report statistical outliers
    phix_outliers = results_df[results_df['phix_outlier']]
    vector_outliers = results_df[results_df['vector_outlier']]
    total_outliers = results_df[results_df['total_outlier']]

    print(f"\n{'STATISTICAL OUTLIERS (>1.5 IQR from median):':<50}")
    print("-" * 80)

    if len(phix_outliers) > 0:
        print(f"\nPhiX outliers ({len(phix_outliers)} samples):")
        for _, row in phix_outliers.iterrows():
            fold_vs_median = row['phix_percent'] / results_df['phix_percent'].median() if results_df['phix_percent'].median() > 0 else float('inf')
            print(f"  {row['sample']:<30} {row['phix_percent']:>7.3f}%  ({fold_vs_median:.1f}x median)")
    else:
        print(f"\nNo PhiX outliers detected (all samples within expected range)")

    if len(vector_outliers) > 0:
        print(f"\nVector outliers ({len(vector_outliers)} samples):")
        for _, row in vector_outliers.iterrows():
            fold_vs_median = row['vector_percent'] / results_df['vector_percent'].median() if results_df['vector_percent'].median() > 0 else float('inf')
            print(f"  {row['sample']:<30} {row['vector_percent']:>7.3f}%  ({fold_vs_median:.1f}x median)")
    else:
        print(f"\nNo vector outliers detected (all samples within expected range)")

    if len(total_outliers) > 0 and not (set(total_outliers['sample']) <= set(phix_outliers['sample']) | set(vector_outliers['sample'])):
        print(f"\nTotal contamination outliers ({len(total_outliers)} samples):")
        for _, row in total_outliers.iterrows():
            print(f"  {row['sample']:<30} {row['total_contamination_percent']:>7.3f}%")

    print("\n" + "="*80)
    print("TIP: Check visualizations to identify samples that deviate from the group")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
