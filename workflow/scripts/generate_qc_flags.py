#!/usr/bin/env python3
"""
Generate QC pass/fail flags for virome samples

Evaluates samples based on:
- ViromeQC enrichment score
- % host reads
- % rRNA reads
- Final read count after QC
- Duplication rate
"""

import pandas as pd
import sys
from pathlib import Path

def parse_viromeqc(viromeqc_file):
    """Parse ViromeQC output to extract enrichment score"""
    enrichment_score = None
    with open(viromeqc_file) as f:
        for line in f:
            if "Enrichment score" in line or "enrichment" in line.lower():
                parts = line.strip().split()
                try:
                    enrichment_score = float(parts[-1])
                except (ValueError, IndexError):
                    pass
    return enrichment_score

def main():
    # Get inputs from snakemake
    viromeqc_files = snakemake.input.viromeqc
    read_counts_file = snakemake.input.read_counts
    output_file = snakemake.output[0]

    # Load thresholds from config
    config = snakemake.config
    thresholds = config.get('qc_thresholds', {})
    min_enrichment = thresholds.get('min_enrichment_score', 10)
    max_host_pct = thresholds.get('max_host_percent', 10)
    max_rrna_pct = thresholds.get('max_rrna_percent', 20)
    min_final = thresholds.get('min_final_reads', 100000)

    # Load read counts
    read_df = pd.read_csv(read_counts_file, sep='\t')

    # Initialize results
    results = []

    # Process each sample - extract from read counts file
    samples = read_df['sample'].unique()

    for sample in samples:
        flags = {
            'sample': sample,
            'enrichment_score': 'NA',
            'pass_enrichment': 'UNKNOWN',
            'pass_host': 'UNKNOWN',
            'pass_rrna': 'UNKNOWN',
            'pass_final_count': 'UNKNOWN',
            'overall_pass': 'UNKNOWN',
            'notes': []
        }

        # Parse ViromeQC enrichment score
        vqc_file = [f for f in viromeqc_files if sample in f][0]
        enrichment = parse_viromeqc(vqc_file)

        if enrichment is not None:
            flags['enrichment_score'] = enrichment
            if enrichment >= min_enrichment:
                flags['pass_enrichment'] = 'PASS'
            else:
                flags['pass_enrichment'] = 'FAIL'
                flags['notes'].append(f'Low_enrichment({enrichment:.2f})')

        # Calculate % reads retained at each step
        sample_reads = read_df[read_df['sample'] == sample]

        if len(sample_reads) > 0:
            fastp_count = sample_reads[sample_reads['step'] == 'fastp']['reads'].values[0]
            host_depleted = sample_reads[sample_reads['step'] == 'host_depleted']['reads'].values[0]
            clean_count = sample_reads[sample_reads['step'] == 'clean']['reads'].values[0]

            # Calculate % host removed (using fastp count as denominator)
            # fastp is the step immediately before host depletion, so this gives
            # the true host contamination % without including QC filtering losses
            host_removed = fastp_count - host_depleted
            host_pct = (host_removed / fastp_count) * 100 if fastp_count > 0 else 0

            if host_pct <= max_host_pct:
                flags['pass_host'] = 'PASS'
            else:
                flags['pass_host'] = 'FAIL'
                flags['notes'].append(f'High_host({host_pct:.1f}%)')

            # Calculate % rRNA removed
            rrna_removed = host_depleted - clean_count
            rrna_pct = (rrna_removed / host_depleted) * 100 if host_depleted > 0 else 0

            if rrna_pct <= max_rrna_pct:
                flags['pass_rrna'] = 'PASS'
            else:
                flags['pass_rrna'] = 'FAIL'
                flags['notes'].append(f'High_rRNA({rrna_pct:.1f}%)')

            # Check final read count
            if clean_count >= min_final:
                flags['pass_final_count'] = 'PASS'
            else:
                flags['pass_final_count'] = 'FAIL'
                flags['notes'].append(f'Low_final_reads({int(clean_count)})')

        # Overall pass/fail
        all_checks = [
            flags['pass_enrichment'],
            flags['pass_host'],
            flags['pass_rrna'],
            flags['pass_final_count']
        ]

        if 'FAIL' in all_checks:
            flags['overall_pass'] = 'FAIL'
        elif 'UNKNOWN' in all_checks:
            flags['overall_pass'] = 'WARNING'
        else:
            flags['overall_pass'] = 'PASS'

        flags['notes'] = ';'.join(flags['notes']) if flags['notes'] else 'None'
        results.append(flags)

    # Write results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)

    # Print summary
    print("\n" + "="*60)
    print("QC FLAGS SUMMARY")
    print("="*60)
    print(results_df.to_string(index=False))
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
