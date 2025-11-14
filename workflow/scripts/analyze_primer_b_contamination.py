"""
Analyze Primer B Cross-Contamination

Parses BBDuk stats files from primer B removal steps to detect cross-contamination.
Each sample should have primarily one primer B variant; multiple variants indicate
cross-contamination between samples during library prep or sequencing.

Outputs:
- Summary table with dominant primer B and contamination metrics
- Per-sample primer B distribution for heatmap visualization
"""

import pandas as pd
import re
from pathlib import Path

# Snakemake inputs
step1_stats = snakemake.input.step1  # List of forward primer B stats files
step2_stats = snakemake.input.step2  # List of RC primer B stats files
sample_assignments_file = snakemake.params.get("sample_assignments", None)
contamination_threshold = snakemake.params.contamination_threshold

# Snakemake outputs
summary_output = snakemake.output.summary
distribution_output = snakemake.output.distribution


def parse_bbduk_stats(stats_file):
    """
    Parse BBDuk stats file to extract primer B read counts per variant

    Returns:
        dict: {sample: str, total_reads: int, primer_b_counts: {variant: count}}
    """
    sample = Path(stats_file).stem.replace('_pb_fwd_stats', '').replace('_pb_rc_stats', '')

    data = {
        'sample': sample,
        'total_reads': 0,
        'primer_b_counts': {}
    }

    with open(stats_file) as f:
        for line in f:
            line = line.strip()

            # Get total read count
            if line.startswith("#Total"):
                fields = line.split()
                if len(fields) >= 2:
                    data['total_reads'] = int(fields[1])

            # Parse per-reference counts (skip header lines)
            elif not line.startswith("#") and line:
                fields = line.split()
                if len(fields) >= 2:
                    variant = fields[0]  # e.g., "3GB-1"
                    count = int(fields[1])
                    data['primer_b_counts'][variant] = count

    return data


def load_sample_assignments(assignments_file):
    """
    Load expected primer B assignments from TSV file

    Format: sample_name\tprimer_b_id

    Returns:
        dict: {sample: expected_primer_b}
    """
    if assignments_file is None or assignments_file == "null":
        return None

    assignments = {}
    with open(assignments_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                parts = line.split('\t')
                if len(parts) >= 2:
                    sample, primer_b = parts[0], parts[1]
                    assignments[sample] = primer_b

    return assignments


def combine_forward_and_rc_counts(fwd_data, rc_data):
    """
    Combine forward and RC primer B counts for the same sample

    Returns:
        dict: Combined counts for all primer B variants
    """
    combined = {}

    # Add forward counts
    for variant, count in fwd_data['primer_b_counts'].items():
        combined[variant] = combined.get(variant, 0) + count

    # Add RC counts
    for variant, count in rc_data['primer_b_counts'].items():
        combined[variant] = combined.get(variant, 0) + count

    return combined


# Parse all stats files
print("Parsing BBDuk stats files...")
fwd_results = [parse_bbduk_stats(f) for f in step1_stats]
rc_results = [parse_bbduk_stats(f) for f in step2_stats]

# Create lookup dictionaries by sample name
fwd_lookup = {r['sample']: r for r in fwd_results}
rc_lookup = {r['sample']: r for r in rc_results}

# Load expected primer B assignments (if provided)
expected_assignments = load_sample_assignments(sample_assignments_file)

# Analyze each sample
summary_data = []
distribution_data = []

samples = set(fwd_lookup.keys()) | set(rc_lookup.keys())

for sample in sorted(samples):
    fwd_data = fwd_lookup.get(sample, {'sample': sample, 'total_reads': 0, 'primer_b_counts': {}})
    rc_data = rc_lookup.get(sample, {'sample': sample, 'total_reads': 0, 'primer_b_counts': {}})

    # Combine forward and RC counts
    combined_counts = combine_forward_and_rc_counts(fwd_data, rc_data)
    total_pb_reads = sum(combined_counts.values())

    # Determine dominant primer B
    if combined_counts:
        dominant_pb = max(combined_counts, key=combined_counts.get)
        dominant_pb_reads = combined_counts[dominant_pb]
        dominant_pb_pct = (dominant_pb_reads / total_pb_reads * 100) if total_pb_reads > 0 else 0
    else:
        dominant_pb = "None"
        dominant_pb_reads = 0
        dominant_pb_pct = 0

    # Use expected primer B if provided, otherwise use auto-detected dominant
    expected_pb = expected_assignments.get(sample, dominant_pb) if expected_assignments else dominant_pb

    # Calculate contamination metrics
    expected_pb_reads = combined_counts.get(expected_pb, 0)
    unexpected_pb_reads = total_pb_reads - expected_pb_reads
    unexpected_pb_pct = (unexpected_pb_reads / total_pb_reads * 100) if total_pb_reads > 0 else 0

    # Flag if contamination exceeds threshold
    is_contaminated = unexpected_pb_pct > (contamination_threshold * 100)

    # Add to summary
    summary_data.append({
        'sample': sample,
        'total_input_reads': fwd_data['total_reads'],
        'total_pb_reads': total_pb_reads,
        'total_pb_pct': (total_pb_reads / fwd_data['total_reads'] * 100) if fwd_data['total_reads'] > 0 else 0,
        'dominant_pb': dominant_pb,
        'dominant_pb_reads': dominant_pb_reads,
        'dominant_pb_pct': dominant_pb_pct,
        'expected_pb': expected_pb,
        'expected_pb_reads': expected_pb_reads,
        'unexpected_pb_reads': unexpected_pb_reads,
        'unexpected_pb_pct': unexpected_pb_pct,
        'contamination_flag': is_contaminated
    })

    # Add distribution data for heatmap
    for variant, count in combined_counts.items():
        pct = (count / total_pb_reads * 100) if total_pb_reads > 0 else 0
        distribution_data.append({
            'sample': sample,
            'primer_b_variant': variant,
            'reads': count,
            'percent_of_pb_reads': pct
        })

# Create DataFrames
summary_df = pd.DataFrame(summary_data)
distribution_df = pd.DataFrame(distribution_data)

# Save outputs
print(f"Saving summary table to {summary_output}")
summary_df.to_csv(summary_output, sep='\t', index=False)

print(f"Saving distribution data to {distribution_output}")
distribution_df.to_csv(distribution_output, sep='\t', index=False)

# Print summary statistics
print("\n=== Primer B Cross-Contamination Analysis ===")
print(f"Total samples analyzed: {len(summary_df)}")
print(f"Samples with contamination >{contamination_threshold*100}%: {summary_df['contamination_flag'].sum()}")
print(f"Median unexpected primer B %: {summary_df['unexpected_pb_pct'].median():.2f}%")
print(f"Max unexpected primer B %: {summary_df['unexpected_pb_pct'].max():.2f}%")

if summary_df['contamination_flag'].sum() > 0:
    print("\nContaminated samples:")
    contaminated = summary_df[summary_df['contamination_flag']]
    for _, row in contaminated.iterrows():
        print(f"  {row['sample']}: {row['unexpected_pb_pct']:.2f}% unexpected primer B reads")
