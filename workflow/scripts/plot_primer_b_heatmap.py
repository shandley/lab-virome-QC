"""
Visualize Primer B Cross-Contamination

Creates a heatmap showing primer B distribution across samples to identify
cross-contamination patterns.

The heatmap shows:
- Rows: Samples
- Columns: Primer B variants (3GB-1 through 3GB-24)
- Color intensity: % of primer B reads for that variant
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Snakemake inputs
distribution_file = snakemake.input.distribution
summary_file = snakemake.input.summary

# Snakemake outputs
heatmap_output = snakemake.output.heatmap

# Load data
print("Loading primer B distribution data...")
dist_df = pd.read_csv(distribution_file, sep='\t')
summary_df = pd.read_csv(summary_file, sep='\t')

# Create pivot table for heatmap (samples × primer B variants)
print("Creating heatmap matrix...")
heatmap_data = dist_df.pivot(
    index='sample',
    columns='primer_b_variant',
    values='percent_of_pb_reads'
).fillna(0)

# Sort columns naturally (3GB-1, 3GB-2, ..., 3GB-24)
def natural_sort_key(s):
    """Sort primer B names naturally (3GB-1, 3GB-2, ..., 3GB-10, 3GB-11, ...)"""
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]

heatmap_data = heatmap_data[sorted(heatmap_data.columns, key=natural_sort_key)]

# Sort rows by dominant primer B for better visualization
sample_order = summary_df.sort_values(['dominant_pb', 'sample'])['sample'].tolist()
heatmap_data = heatmap_data.reindex(sample_order)

# Create figure
fig_height = max(8, len(heatmap_data) * 0.3)  # Scale height with number of samples
fig_width = max(12, len(heatmap_data.columns) * 0.4)  # Scale width with number of primer B variants

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Create heatmap
sns.heatmap(
    heatmap_data,
    cmap='YlOrRd',
    cbar_kws={'label': '% of Primer B Reads'},
    linewidths=0.5,
    linecolor='lightgray',
    ax=ax,
    vmin=0,
    vmax=100,
    square=False
)

# Formatting
ax.set_xlabel('Primer B Variant', fontsize=12, fontweight='bold')
ax.set_ylabel('Sample', fontsize=12, fontweight='bold')
ax.set_title('Primer B Cross-Contamination Heatmap\n(% of Primer B Reads per Variant)',
             fontsize=14, fontweight='bold', pad=20)

# Rotate x-axis labels for readability
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, fontsize=8)

# Add contamination flags as colored markers
contaminated_samples = summary_df[summary_df['contamination_flag']]['sample'].tolist()
for i, sample in enumerate(heatmap_data.index):
    if sample in contaminated_samples:
        # Add red marker on y-axis for contaminated samples
        ax.text(-0.5, i + 0.5, '⚠', ha='right', va='center',
                fontsize=12, color='red', fontweight='bold')

# Adjust layout
plt.tight_layout()

# Save figure
print(f"Saving heatmap to {heatmap_output}")
plt.savefig(heatmap_output, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== Heatmap Statistics ===")
print(f"Samples: {len(heatmap_data)}")
print(f"Primer B variants: {len(heatmap_data.columns)}")
print(f"Contaminated samples (marked with ⚠): {len(contaminated_samples)}")

# Print dominant primer B distribution
print("\nDominant Primer B Distribution:")
dominant_counts = summary_df['dominant_pb'].value_counts()
for pb, count in dominant_counts.items():
    print(f"  {pb}: {count} samples")
