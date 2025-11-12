#!/usr/bin/env python3
"""
Generate contamination visualization plots

Creates publication-quality plots emphasizing outlier detection:
1. Bar plot highlighting samples that deviate from the group
2. Box plots showing distribution with outliers marked
3. Scatter plot showing correlation with outliers labeled
4. Heatmap showing patterns across samples
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path


def identify_outliers(series):
    """
    Identify statistical outliers using IQR method
    Returns boolean mask of outliers
    """
    if len(series) < 4:
        return pd.Series([False] * len(series), index=series.index)

    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1

    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    return (series < lower_bound) | (series > upper_bound)


def plot_contamination_bars(df, output_prefix):
    """
    Create stacked bar plot with outliers highlighted
    """
    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.5), 6))

    # Identify outliers
    df['is_outlier'] = identify_outliers(df['total_contamination_percent'])

    # Sort by total contamination
    df_sorted = df.sort_values('total_contamination_percent', ascending=True).copy()

    # Create bars with different colors for outliers
    x = range(len(df_sorted))
    width = 0.8

    # Color outliers differently
    phix_colors = ['#D55E00' if outlier else '#E69F00' for outlier in df_sorted['is_outlier']]
    vector_colors = ['#CC79A7' if outlier else '#56B4E9' for outlier in df_sorted['is_outlier']]

    # Plot bars (without label - we'll add custom legend)
    ax.bar(x, df_sorted['phix_percent'], width, color=phix_colors, alpha=0.8,
           edgecolor='black', linewidth=0.5)
    ax.bar(x, df_sorted['vector_percent'], width,
           bottom=df_sorted['phix_percent'],
           color=vector_colors, alpha=0.8, edgecolor='black', linewidth=0.5)

    # Add median reference line
    median_total = df['total_contamination_percent'].median()
    ax.axhline(y=median_total, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)

    # Add IQR shaded region
    Q1 = df['total_contamination_percent'].quantile(0.25)
    Q3 = df['total_contamination_percent'].quantile(0.75)
    ax.axhspan(Q1, Q3, alpha=0.1, color='gray')

    # Create custom legend with all color variations
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E69F00', edgecolor='black', alpha=0.8, label='PhiX (normal)'),
        Patch(facecolor='#D55E00', edgecolor='black', alpha=0.8, label='PhiX (outlier)'),
        Patch(facecolor='#56B4E9', edgecolor='black', alpha=0.8, label='Vector (normal)'),
        Patch(facecolor='#CC79A7', edgecolor='black', alpha=0.8, label='Vector (outlier)'),
        plt.Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, alpha=0.7,
                   label=f'Median ({median_total:.2f}%)'),
        Patch(facecolor='gray', alpha=0.1, label='25-75th percentile')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True, fontsize=9)

    # Formatting
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    ax.set_title('Contamination Levels by Sample (Outliers Highlighted)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(df_sorted['sample'], rotation=45, ha='right', fontsize=9)

    # Highlight outlier sample names in red
    for i, (idx, row) in enumerate(df_sorted.iterrows()):
        if row['is_outlier']:
            ax.get_xticklabels()[i].set_color('red')
            ax.get_xticklabels()[i].set_fontweight('bold')

    ax.legend(loc='upper left', frameon=True, fancybox=True, fontsize=9)
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
    Create box plots showing distribution with outliers clearly marked
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Identify outliers for each type
    phix_outliers = identify_outliers(df['phix_percent'])
    vector_outliers = identify_outliers(df['vector_percent'])
    total_outliers = identify_outliers(df['total_contamination_percent'])

    # PhiX distribution
    bp1 = axes[0].boxplot([df['phix_percent']], labels=['PhiX'],
                          patch_artist=True,
                          boxprops=dict(facecolor='#E69F00', alpha=0.7, linewidth=1.5),
                          medianprops=dict(color='black', linewidth=2),
                          whiskerprops=dict(linewidth=1.5),
                          capprops=dict(linewidth=1.5),
                          flierprops=dict(marker='D', markerfacecolor='red', markersize=8,
                                        markeredgecolor='darkred', alpha=0.8))
    axes[0].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    axes[0].set_title(f'PhiX Distribution\n({phix_outliers.sum()} outliers)', fontsize=12, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3, linestyle='--')

    # Vector distribution
    bp2 = axes[1].boxplot([df['vector_percent']], labels=['Vector/Plasmid'],
                          patch_artist=True,
                          boxprops=dict(facecolor='#56B4E9', alpha=0.7, linewidth=1.5),
                          medianprops=dict(color='black', linewidth=2),
                          whiskerprops=dict(linewidth=1.5),
                          capprops=dict(linewidth=1.5),
                          flierprops=dict(marker='D', markerfacecolor='red', markersize=8,
                                        markeredgecolor='darkred', alpha=0.8))
    axes[1].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    axes[1].set_title(f'Vector/Plasmid Distribution\n({vector_outliers.sum()} outliers)', fontsize=12, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3, linestyle='--')

    # Total contamination distribution
    bp3 = axes[2].boxplot([df['total_contamination_percent']], labels=['Total'],
                          patch_artist=True,
                          boxprops=dict(facecolor='#999999', alpha=0.7, linewidth=1.5),
                          medianprops=dict(color='black', linewidth=2),
                          whiskerprops=dict(linewidth=1.5),
                          capprops=dict(linewidth=1.5),
                          flierprops=dict(marker='D', markerfacecolor='red', markersize=8,
                                        markeredgecolor='darkred', alpha=0.8))
    axes[2].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold')
    axes[2].set_title(f'Total Contamination Distribution\n({total_outliers.sum()} outliers)', fontsize=12, fontweight='bold')
    axes[2].grid(axis='y', alpha=0.3, linestyle='--')

    # Add statistics text
    for ax, data, name in zip(axes,
                               [df['phix_percent'], df['vector_percent'], df['total_contamination_percent']],
                               ['PhiX', 'Vector', 'Total']):
        stats_text = f"Median: {data.median():.3f}%\nIQR: {data.quantile(0.75) - data.quantile(0.25):.3f}%"
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
               fontsize=9)

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_boxes.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved box plot: {output_file}")


def plot_contamination_scatter(df, output_prefix):
    """
    Create scatter plot showing correlation with statistical outliers labeled
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Identify outliers
    phix_outliers = identify_outliers(df['phix_percent'])
    vector_outliers = identify_outliers(df['vector_percent'])
    any_outlier = phix_outliers | vector_outliers

    # Separate outliers from normal samples
    df_normal = df[~any_outlier]
    df_outliers = df[any_outlier]

    # Plot normal samples
    if len(df_normal) > 0:
        ax.scatter(df_normal['phix_percent'], df_normal['vector_percent'],
                  s=100, alpha=0.6, c=df_normal['total_contamination_percent'],
                  cmap='YlOrRd', edgecolors='gray', linewidth=0.5,
                  label='Within expected range')

    # Plot outliers with different style
    if len(df_outliers) > 0:
        scatter = ax.scatter(df_outliers['phix_percent'], df_outliers['vector_percent'],
                           s=150, alpha=0.8, c=df_outliers['total_contamination_percent'],
                           cmap='YlOrRd', edgecolors='red', linewidth=2,
                           marker='D', label='Statistical outliers')

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap='YlOrRd',
                               norm=plt.Normalize(vmin=df['total_contamination_percent'].min(),
                                                vmax=df['total_contamination_percent'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Total Contamination (%)', fontsize=11, fontweight='bold')

    # Add median reference lines
    phix_median = df['phix_percent'].median()
    vector_median = df['vector_percent'].median()
    ax.axhline(y=vector_median, color='gray', linestyle='--', linewidth=1, alpha=0.5,
              label=f'Median vector ({vector_median:.2f}%)')
    ax.axvline(x=phix_median, color='gray', linestyle='--', linewidth=1, alpha=0.5,
              label=f'Median PhiX ({phix_median:.2f}%)')

    # Label outliers
    for _, row in df_outliers.iterrows():
        # Calculate fold-change vs median for annotation
        phix_fold = row['phix_percent'] / phix_median if phix_median > 0 else float('inf')
        vector_fold = row['vector_percent'] / vector_median if vector_median > 0 else float('inf')
        max_fold = max(phix_fold, vector_fold)

        ax.annotate(f"{row['sample']}\n({max_fold:.1f}x)",
                   (row['phix_percent'], row['vector_percent']),
                   xytext=(8, 8), textcoords='offset points',
                   fontsize=8, fontweight='bold', color='red',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.6, edgecolor='red'),
                   arrowprops=dict(arrowstyle='->', color='red', lw=1))

    # Formatting
    ax.set_xlabel('PhiX Contamination (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Vector/Plasmid Contamination (%)', fontsize=12, fontweight='bold')
    ax.set_title('Contamination Correlation - Outliers Highlighted', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best', frameon=True, fancybox=True, fontsize=9)

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_scatter.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved scatter plot: {output_file}")


def plot_contamination_heatmap(df, output_prefix):
    """
    Create heatmap showing contamination patterns with outliers marked
    """
    # Identify outliers
    total_outliers = identify_outliers(df['total_contamination_percent'])

    # Prepare data for heatmap
    heatmap_data = df[['sample', 'phix_percent', 'vector_percent']].set_index('sample').T

    # Sort columns by total contamination
    col_order = df.sort_values('total_contamination_percent', ascending=False)['sample']
    heatmap_data = heatmap_data[col_order]

    # Create figure
    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.4), 5))

    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd',
                cbar_kws={'label': 'Contamination (%)'},
                linewidths=0.5, linecolor='gray', ax=ax,
                vmin=0, vmax=df[['phix_percent', 'vector_percent']].max().max())

    # Highlight outlier columns with colored background
    outlier_samples = df[total_outliers]['sample'].values
    for i, sample in enumerate(col_order):
        if sample in outlier_samples:
            # Add red box around outlier columns
            col_idx = i
            ax.add_patch(plt.Rectangle((col_idx, 0), 1, 2, fill=False,
                                      edgecolor='red', lw=3, zorder=10))

    # Formatting
    ax.set_xlabel('Sample (sorted by total contamination)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Contamination Type', fontsize=12, fontweight='bold')
    ax.set_yticklabels(['PhiX', 'Vector/Plasmid'], rotation=0, fontsize=11)

    # Color outlier sample names in red
    xticklabels = ax.get_xticklabels()
    for i, label in enumerate(xticklabels):
        if label.get_text() in outlier_samples:
            label.set_color('red')
            label.set_fontweight('bold')

    ax.set_title(f'Contamination Heatmap (Outliers Boxed in Red)', fontsize=14, fontweight='bold')

    # Add median reference line annotation
    median_total = df['total_contamination_percent'].median()
    ax.text(0.02, 0.98, f"Median total: {median_total:.3f}%",
           transform=ax.transAxes, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
           fontsize=10, fontweight='bold')

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
