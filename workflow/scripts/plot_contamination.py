#!/usr/bin/env python3
"""
Generate contamination visualization plots

Creates publication-quality plots showing:
1. Bar plot of PhiX and vector contamination per sample
2. Scatter plot showing correlation between contamination types
3. Box plots showing distribution of contamination levels
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def plot_contamination_bars(df, output_prefix):
    """
    Create stacked bar plot of contamination by sample
    """
    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.5), 6))

    # Sort by total contamination
    df_sorted = df.sort_values('total_contamination_percent', ascending=True)

    # Create bars
    x = range(len(df_sorted))
    width = 0.8

    ax.bar(x, df_sorted['phix_percent'], width, label='PhiX', color='#E69F00', alpha=0.8)
    ax.bar(x, df_sorted['vector_percent'], width,
           bottom=df_sorted['phix_percent'], label='Vector/Plasmid',
           color='#56B4E9', alpha=0.8)

    # Add threshold line at 1%
    ax.axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5,
               label='1% threshold')

    # Formatting
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    ax.set_title('Contamination Levels by Sample', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(df_sorted['sample'], rotation=45, ha='right')
    ax.legend(loc='upper left', frameon=True, fancybox=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Adjust layout
    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_bars.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved bar plot: {output_file}")


def plot_contamination_boxes(df, output_prefix):
    """
    Create box plots showing distribution of contamination types
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # PhiX distribution
    axes[0].boxplot([df['phix_percent']], labels=['PhiX'],
                     patch_artist=True,
                     boxprops=dict(facecolor='#E69F00', alpha=0.7),
                     medianprops=dict(color='black', linewidth=2))
    axes[0].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    axes[0].set_title('PhiX Contamination Distribution', fontsize=12, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3, linestyle='--')
    axes[0].axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5)

    # Vector distribution
    axes[1].boxplot([df['vector_percent']], labels=['Vector/Plasmid'],
                     patch_artist=True,
                     boxprops=dict(facecolor='#56B4E9', alpha=0.7),
                     medianprops=dict(color='black', linewidth=2))
    axes[1].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    axes[1].set_title('Vector/Plasmid Contamination Distribution', fontsize=12, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3, linestyle='--')
    axes[1].axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5)

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_boxes.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved box plot: {output_file}")


def plot_contamination_scatter(df, output_prefix):
    """
    Create scatter plot showing correlation between PhiX and vector contamination
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create scatter plot
    ax.scatter(df['phix_percent'], df['vector_percent'],
               s=100, alpha=0.6, c=df['total_contamination_percent'],
               cmap='YlOrRd', edgecolors='black', linewidth=0.5)

    # Add colorbar
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('Total Contamination (%)', fontsize=11, fontweight='bold')

    # Add reference lines
    ax.axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5)

    # Label outliers (samples with >1% in either category)
    for _, row in df.iterrows():
        if row['phix_percent'] > 1.0 or row['vector_percent'] > 1.0:
            ax.annotate(row['sample'],
                       (row['phix_percent'], row['vector_percent']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)

    # Formatting
    ax.set_xlabel('PhiX Contamination (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Vector/Plasmid Contamination (%)', fontsize=12, fontweight='bold')
    ax.set_title('Contamination Correlation Analysis', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_scatter.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved scatter plot: {output_file}")


def plot_contamination_heatmap(df, output_prefix):
    """
    Create heatmap showing contamination levels across samples
    """
    # Prepare data for heatmap
    heatmap_data = df[['sample', 'phix_percent', 'vector_percent']].set_index('sample').T

    # Sort columns by total contamination
    col_order = df.sort_values('total_contamination_percent', ascending=False)['sample']
    heatmap_data = heatmap_data[col_order]

    # Create figure
    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.4), 4))

    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd',
                cbar_kws={'label': 'Contamination (%)'},
                linewidths=0.5, linecolor='gray', ax=ax)

    # Formatting
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Contamination Type', fontsize=12, fontweight='bold')
    ax.set_yticklabels(['PhiX', 'Vector/Plasmid'], rotation=0)
    ax.set_title('Contamination Heatmap', fontsize=14, fontweight='bold')

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved heatmap: {output_file}")


def main():
    # Get inputs from snakemake
    contamination_table = snakemake.input[0]
    output_prefix = snakemake.params.output_prefix

    # Load data
    df = pd.read_csv(contamination_table, sep='\t')

    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 10

    # Generate plots
    print("\n" + "="*60)
    print("GENERATING CONTAMINATION PLOTS")
    print("="*60 + "\n")

    plot_contamination_bars(df, output_prefix)
    plot_contamination_boxes(df, output_prefix)
    plot_contamination_scatter(df, output_prefix)
    plot_contamination_heatmap(df, output_prefix)

    print("\n" + "="*60)
    print("PLOTTING COMPLETE")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
