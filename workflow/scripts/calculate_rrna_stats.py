#!/usr/bin/env python3
"""
Calculate rRNA removal statistics across all samples
"""

import pandas as pd
import sys
from pathlib import Path


def calculate_rrna_stats(read_counts_file):
    """
    Calculate rRNA removal percentages from read counts file

    Parameters
    ----------
    read_counts_file : str
        Path to read_counts.tsv file

    Returns
    -------
    DataFrame with rRNA statistics per sample
    """
    # Load read counts
    df = pd.read_csv(read_counts_file, sep='\t')

    # Get unique samples
    samples = df['sample'].unique()

    results = []

    for sample in samples:
        sample_data = df[df['sample'] == sample]

        # Get read counts at key steps
        host_depleted = sample_data[sample_data['step'] == 'host_depleted']['reads'].values
        clean = sample_data[sample_data['step'] == 'clean']['reads'].values

        if len(host_depleted) > 0 and len(clean) > 0:
            host_depleted_count = host_depleted[0]
            clean_count = clean[0]

            # Calculate rRNA removed
            rrna_removed = host_depleted_count - clean_count
            rrna_pct = (rrna_removed / host_depleted_count) * 100 if host_depleted_count > 0 else 0

            results.append({
                'sample': sample,
                'host_depleted_reads': int(host_depleted_count),
                'clean_reads': int(clean_count),
                'rrna_removed': int(rrna_removed),
                'rrna_percent': round(rrna_pct, 2)
            })

    return pd.DataFrame(results)


def main():
    if len(sys.argv) < 2:
        print("Usage: python calculate_rrna_stats.py <read_counts.tsv>")
        sys.exit(1)

    read_counts_file = sys.argv[1]

    if not Path(read_counts_file).exists():
        print(f"Error: File not found: {read_counts_file}")
        sys.exit(1)

    # Calculate stats
    stats_df = calculate_rrna_stats(read_counts_file)

    if len(stats_df) == 0:
        print("No samples found in read counts file")
        sys.exit(1)

    # Print summary statistics
    print("\n" + "=" * 80)
    print("rRNA REMOVAL STATISTICS")
    print("=" * 80)
    print(f"\nTotal samples: {len(stats_df)}")
    print(f"\nrRNA Removal Percentage:")
    print(f"  Mean:   {stats_df['rrna_percent'].mean():.2f}%")
    print(f"  Median: {stats_df['rrna_percent'].median():.2f}%")
    print(f"  Min:    {stats_df['rrna_percent'].min():.2f}%")
    print(f"  Max:    {stats_df['rrna_percent'].max():.2f}%")
    print(f"  Std:    {stats_df['rrna_percent'].std():.2f}%")

    # Show distribution
    print(f"\nDistribution:")
    bins = [0, 20, 40, 60, 80, 100]
    labels = ['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']
    stats_df['rrna_bin'] = pd.cut(stats_df['rrna_percent'], bins=bins, labels=labels, include_lowest=True)
    distribution = stats_df['rrna_bin'].value_counts().sort_index()
    for label, count in distribution.items():
        print(f"  {label}: {count} samples ({count/len(stats_df)*100:.1f}%)")

    # Show samples with highest and lowest rRNA
    print(f"\nHighest rRNA contamination (top 5):")
    top5 = stats_df.nlargest(5, 'rrna_percent')
    for _, row in top5.iterrows():
        print(f"  {row['sample']}: {row['rrna_percent']:.2f}%")

    print(f"\nLowest rRNA contamination (bottom 5):")
    bottom5 = stats_df.nsmallest(5, 'rrna_percent')
    for _, row in bottom5.iterrows():
        print(f"  {row['sample']}: {row['rrna_percent']:.2f}%")

    print("\n" + "=" * 80)

    # Optionally save detailed results
    output_file = read_counts_file.replace('read_counts.tsv', 'rrna_statistics.tsv')
    stats_df.drop('rrna_bin', axis=1).to_csv(output_file, sep='\t', index=False)
    print(f"\nDetailed statistics saved to: {output_file}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
