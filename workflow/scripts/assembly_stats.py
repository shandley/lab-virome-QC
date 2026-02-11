#!/usr/bin/env python3
"""
Calculate assembly statistics from FASTA contigs

Two-Stage Assembly Workflow Statistics:
Computes metrics for both the final Flye meta-assembly and the intermediate
MEGAHIT assemblies (per-sample or per-group).

Metrics include:
- Number of contigs
- Total assembly size (bp)
- N50, L50 (standard assembly quality metrics)
- Longest contig
- Mean contig length
- GC content (%)
"""

import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np
import sys


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return list of sequence records

    Args:
        fasta_file: Path to FASTA file

    Returns:
        list: List of SeqRecord objects
    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def calculate_n50_l50(lengths):
    """
    Calculate N50 and L50 metrics

    N50: Length of the shortest contig for which longer and equal contigs
         cover at least 50% of the assembly
    L50: Number of contigs whose summed length is N50

    Args:
        lengths: List of contig lengths

    Returns:
        tuple: (N50, L50)
    """
    if not lengths:
        return 0, 0

    # Sort lengths in descending order
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    half_length = total_length / 2

    cumsum = 0
    n50 = 0
    l50 = 0

    for length in sorted_lengths:
        cumsum += length
        l50 += 1
        if cumsum >= half_length:
            n50 = length
            break

    return n50, l50


def calculate_gc_content(sequences):
    """
    Calculate overall GC content across all sequences

    Args:
        sequences: List of SeqRecord objects

    Returns:
        float: GC content as percentage
    """
    total_gc = 0
    total_bases = 0

    for seq in sequences:
        seq_str = str(seq.seq).upper()
        gc_count = seq_str.count('G') + seq_str.count('C')
        total_gc += gc_count
        total_bases += len(seq_str)

    if total_bases == 0:
        return 0.0

    return (total_gc / total_bases) * 100


def calculate_assembly_stats(fasta_file, assembly_name=None, assembly_type="megahit"):
    """
    Calculate comprehensive assembly statistics for a single assembly

    Args:
        fasta_file: Path to FASTA file
        assembly_name: Name for this assembly (sample, group, or 'flye_final')
        assembly_type: Type of assembly ('megahit' or 'flye')

    Returns:
        dict: Assembly statistics
    """
    # Parse sequences
    try:
        sequences = parse_fasta(fasta_file)
    except Exception as e:
        print(f"Warning: Could not parse {fasta_file}: {e}", file=sys.stderr)
        sequences = []

    # Basic counts
    num_contigs = len(sequences)

    # Handle empty assemblies
    if num_contigs == 0:
        return {
            'assembly_name': assembly_name or 'unknown',
            'assembly_type': assembly_type,
            'num_contigs': 0,
            'total_size_bp': 0,
            'mean_length_bp': 0,
            'longest_contig_bp': 0,
            'n50': 0,
            'l50': 0,
            'gc_percent': 0.0
        }

    # Get all contig lengths
    lengths = [len(seq.seq) for seq in sequences]

    # Calculate statistics
    total_size = sum(lengths)
    mean_length = np.mean(lengths)
    longest_contig = max(lengths)
    n50, l50 = calculate_n50_l50(lengths)
    gc_content = calculate_gc_content(sequences)

    return {
        'assembly_name': assembly_name or 'unknown',
        'assembly_type': assembly_type,
        'num_contigs': num_contigs,
        'total_size_bp': total_size,
        'mean_length_bp': int(mean_length),
        'longest_contig_bp': longest_contig,
        'n50': n50,
        'l50': l50,
        'gc_percent': round(gc_content, 2)
    }


def parse_flye_info(info_file):
    """
    Parse Flye assembly_info.txt for additional metrics

    Args:
        info_file: Path to Flye assembly_info.txt

    Returns:
        dict: Additional Flye-specific metrics
    """
    metrics = {
        'circular_contigs': 0,
        'repeat_contigs': 0
    }

    try:
        with open(info_file, 'r') as f:
            # Skip header
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    # Check circular status (column 4)
                    if parts[3] == 'Y':
                        metrics['circular_contigs'] += 1
                    # Check repeat status (column 5, if exists)
                    if len(parts) >= 5 and parts[4] == 'Y':
                        metrics['repeat_contigs'] += 1
    except Exception as e:
        print(f"Warning: Could not parse Flye info file: {e}", file=sys.stderr)

    return metrics


def main():
    """Main function called by Snakemake"""
    # Get inputs from snakemake
    flye_assembly = snakemake.input.flye_assembly
    flye_info = snakemake.input.flye_info
    megahit_contigs = snakemake.input.megahit_contigs
    output_file = snakemake.output.stats
    log_file = snakemake.log[0] if snakemake.log else None

    # Get assembly strategy
    assembly_strategy = snakemake.params.get("strategy", "by_sample")

    results = []

    # Log start
    if log_file:
        with open(log_file, 'w') as log:
            log.write(f"Calculating assembly statistics\n")
            log.write(f"Assembly strategy: {assembly_strategy}\n")
            log.write(f"Flye assembly: {flye_assembly}\n")
            log.write(f"MEGAHIT contigs: {megahit_contigs}\n\n")

    # 1. Calculate stats for final Flye assembly
    flye_stats = calculate_assembly_stats(
        flye_assembly,
        assembly_name="FINAL_FLYE",
        assembly_type="flye"
    )

    # Add Flye-specific metrics
    flye_metrics = parse_flye_info(flye_info)
    flye_stats.update(flye_metrics)

    results.append(flye_stats)

    # 2. Calculate stats for each MEGAHIT assembly (per-sample or per-group)
    if isinstance(megahit_contigs, str):
        megahit_contigs = [megahit_contigs]

    for contig_file in megahit_contigs:
        # Extract name from path
        # For individual: .../assembly/per_sample/{sample}/final.contigs.fa
        # For groups: .../assembly/per_group/{group}/final.contigs.fa
        contig_path = Path(contig_file)
        assembly_name = contig_path.parent.name

        stats = calculate_assembly_stats(
            contig_file,
            assembly_name=assembly_name,
            assembly_type="megahit"
        )
        results.append(stats)

    # Create DataFrame
    df = pd.DataFrame(results)

    # Reorder columns for better readability
    column_order = [
        'assembly_name', 'assembly_type', 'num_contigs', 'total_size_bp',
        'mean_length_bp', 'longest_contig_bp', 'n50', 'l50', 'gc_percent'
    ]

    # Add Flye-specific columns if present
    if 'circular_contigs' in df.columns:
        column_order.extend(['circular_contigs', 'repeat_contigs'])

    # Ensure all columns exist
    for col in column_order:
        if col not in df.columns:
            df[col] = 0

    df = df[column_order]

    # Sort: FINAL_FLYE first, then alphabetically
    df['sort_key'] = df['assembly_name'].apply(lambda x: '0' if x == 'FINAL_FLYE' else x)
    df = df.sort_values('sort_key').drop('sort_key', axis=1)

    # Write output
    df.to_csv(output_file, sep='\t', index=False)

    # Log completion
    if log_file:
        with open(log_file, 'a') as log:
            log.write(f"\nAssembly statistics written to {output_file}\n")
            log.write(f"Processed {len(results)} assemblies:\n")
            log.write(f"  - 1 final Flye meta-assembly\n")
            log.write(f"  - {len(results) - 1} MEGAHIT assemblies\n")

    print(f"Assembly statistics written to {output_file}")
    print(f"Processed {len(results)} assemblies (1 Flye + {len(results) - 1} MEGAHIT)")


if __name__ == '__main__':
    main()
