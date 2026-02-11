#!/usr/bin/env python3
"""
Generate user-friendly virome QC report (HTML + PDF)

This script consolidates all QC data sources into a unified structure,
implements enhanced outlier detection, and generates both interactive HTML
dashboard and publication-ready PDF reports.

Replaces MultiQC with virome-specific reporting focused on:
- Batch-level overview and outlier identification
- Quality metric visualization and assessment
- Sample ranking by practical metrics
- Interactive exploration and static sharing
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
import re
import sys
from typing import Dict, List, Tuple, Optional, Any
import warnings
from datetime import datetime
import base64

# HTML template generation
try:
    from jinja2 import Environment, FileSystemLoader, Template
    JINJA2_AVAILABLE = True
except ImportError:
    JINJA2_AVAILABLE = False
    print("Warning: jinja2 not available. HTML generation will be skipped.")

# Suppress pandas warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)

class ViromePipelineDataLoader:
    """
    Unified data loader for all virome QC data sources
    """

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.qc_thresholds = config.get('qc_thresholds', {})
        self.samples = []
        self.data = {}

    def load_read_counts(self, read_counts_file: str) -> pd.DataFrame:
        """Load and process read count data"""
        read_df = pd.read_csv(read_counts_file, sep='\t')

        # Pivot to have samples as rows and steps as columns
        read_pivot = read_df.pivot(index='sample', columns='step', values='reads')

        # Calculate retention percentages
        if 'raw' in read_pivot.columns:
            for col in read_pivot.columns:
                if col != 'raw':
                    read_pivot[f'{col}_retention'] = (read_pivot[col] / read_pivot['raw'] * 100)

        return read_pivot

    def load_contamination_data(self, contamination_file: str) -> pd.DataFrame:
        """Load contamination summary data"""
        return pd.read_csv(contamination_file, sep='\t')

    def load_qc_metrics(self, qc_metrics_file: str) -> pd.DataFrame:
        """Load QC quality metrics"""
        return pd.read_csv(qc_metrics_file, sep='\t')

    def load_primer_b_data(self, primer_b_file: str) -> pd.DataFrame:
        """Load primer B contamination data"""
        return pd.read_csv(primer_b_file, sep='\t')


    def parse_fastqc_summary(self, fastqc_files: List[str]) -> pd.DataFrame:
        """Parse FastQC summary statistics from ZIP files"""
        fastqc_data = []

        for fastqc_file in fastqc_files:
            # Extract sample name and stage (raw/trimmed/final)
            file_path = Path(fastqc_file)
            filename = file_path.stem

            # Parse sample name and stage from filename
            # Format: {sample}_R1_fastqc.zip in raw/trimmed/final directories
            stage = file_path.parent.name  # raw, trimmed, or final
            sample = filename.replace('_R1_fastqc', '').replace('_R2_fastqc', '')

            # For now, store basic metadata - can extend to parse ZIP contents if needed
            fastqc_data.append({
                'sample': sample,
                'stage': stage,
                'file_path': str(fastqc_file)
            })

        return pd.DataFrame(fastqc_data)

class ViromePipelineOutlierDetector:
    """
    Enhanced outlier detection across all QC metrics
    """

    def __init__(self):
        self.outlier_methods = {
            'iqr': self._iqr_outliers,
            'zscore': self._zscore_outliers,
            'modified_zscore': self._modified_zscore_outliers
        }

    def _iqr_outliers(self, series: pd.Series) -> pd.Series:
        """IQR-based outlier detection (existing method)"""
        if len(series) < 4:
            return pd.Series([False] * len(series), index=series.index)

        Q1 = series.quantile(0.25)
        Q3 = series.quantile(0.75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        return (series < lower_bound) | (series > upper_bound)

    def _zscore_outliers(self, series: pd.Series, threshold: float = 2.0) -> pd.Series:
        """Z-score based outlier detection"""
        if len(series) < 4:
            return pd.Series([False] * len(series), index=series.index)

        z_scores = np.abs((series - series.mean()) / series.std())
        return z_scores > threshold

    def _modified_zscore_outliers(self, series: pd.Series, threshold: float = 3.5) -> pd.Series:
        """Modified Z-score using median absolute deviation"""
        if len(series) < 4:
            return pd.Series([False] * len(series), index=series.index)

        median = series.median()
        mad = np.median(np.abs(series - median))

        if mad == 0:
            return pd.Series([False] * len(series), index=series.index)

        modified_z_scores = 0.6745 * (series - median) / mad
        return np.abs(modified_z_scores) > threshold

    def detect_outliers(self, data: pd.DataFrame, columns: List[str],
                       method: str = 'iqr') -> pd.DataFrame:
        """
        Detect outliers across multiple metrics
        """
        outlier_df = data.copy()

        if method not in self.outlier_methods:
            raise ValueError(f"Unknown outlier method: {method}")

        detect_func = self.outlier_methods[method]

        for col in columns:
            if col in data.columns:
                outlier_df[f'{col}_outlier'] = detect_func(data[col])

        return outlier_df

class ViromePipelineReportGenerator:
    """
    Main class for generating virome QC reports
    """

    def __init__(self, inputs: Dict[str, Any], outputs: Dict[str, str],
                 config: Dict[str, Any]):
        self.inputs = inputs
        self.outputs = outputs
        self.config = config

        self.loader = ViromePipelineDataLoader(config)
        self.outlier_detector = ViromePipelineOutlierDetector()
        self.unified_data = None

    def load_all_data(self):
        """Load and unify all QC data sources"""
        print("Loading QC data sources...")

        # Load individual data sources
        read_counts = self.loader.load_read_counts(self.inputs['read_counts'])
        contamination = self.loader.load_contamination_data(self.inputs['contamination_summary'])
        qc_metrics = self.loader.load_qc_metrics(self.inputs['qc_metrics'])
        primer_b = self.loader.load_primer_b_data(self.inputs['primer_b_summary'])

        # Start with QC metrics as the base (has all samples)
        unified = qc_metrics.set_index('sample')

        # Merge other data sources
        if not read_counts.empty:
            unified = unified.join(read_counts, how='left', rsuffix='_reads')

        if not contamination.empty:
            contamination_indexed = contamination.set_index('sample')
            unified = unified.join(contamination_indexed, how='left', rsuffix='_contam')

        if not primer_b.empty:
            primer_b_indexed = primer_b.set_index('sample')
            unified = unified.join(primer_b_indexed, how='left', rsuffix='_pb')

        # Reset index to have sample as column
        unified = unified.reset_index()

        # QC metrics now use simplified 3-metric system

        self.unified_data = unified
        print(f"Successfully loaded data for {len(unified)} samples")

    def rank_samples_by_metrics(self):
        """Rank samples by final read count (most practical metric for prioritization)"""
        if self.unified_data is None:
            raise ValueError("Must load data first")

        # Rank by final read count - most practical for determining sample usability
        if 'clean' in self.unified_data.columns:
            self.unified_data['read_count_rank'] = self.unified_data['clean'].rank(
                method='dense', ascending=False
            ).astype(int)
        else:
            # Fallback ranking by index if clean reads not available
            self.unified_data['read_count_rank'] = range(1, len(self.unified_data) + 1)

    def detect_all_outliers(self):
        """Detect outliers across all quality metrics"""
        outlier_columns = [
            'clean_retention',
            'host_percent',
            'rrna_percent',
            'phix_percent',
            'vector_percent'
        ]

        # Filter to columns that exist in the data
        existing_columns = [col for col in outlier_columns if col in self.unified_data.columns]

        self.unified_data = self.outlier_detector.detect_outliers(
            self.unified_data, existing_columns, method='iqr'
        )

        # Add composite outlier flag (outlier in any metric)
        outlier_cols = [col for col in self.unified_data.columns if col.endswith('_outlier')]
        self.unified_data['any_outlier'] = self.unified_data[outlier_cols].any(axis=1)

    def calculate_batch_statistics(self) -> Dict[str, Any]:
        """Calculate batch-level summary statistics"""
        stats = {
            'total_samples': len(self.unified_data),
            'outlier_count': len(self.unified_data[self.unified_data['any_outlier'] == True]),
        }

        # Calculate read count statistics (convert to native Python types for JSON serialization)
        if 'raw' in self.unified_data.columns:
            stats['total_reads'] = int(self.unified_data['raw'].sum())
            stats['total_reads_millions'] = round(float(self.unified_data['raw'].sum()) / 1000000, 1)
            stats['median_raw_reads'] = int(self.unified_data['raw'].median())
            stats['median_raw_reads_millions'] = round(float(self.unified_data['raw'].median()) / 1000000, 1)

        if 'clean' in self.unified_data.columns:
            stats['total_final_reads'] = int(self.unified_data['clean'].sum())
            stats['total_final_reads_millions'] = round(float(self.unified_data['clean'].sum()) / 1000000, 1)
            stats['median_final_reads'] = int(self.unified_data['clean'].median())
            stats['median_final_reads_millions'] = round(float(self.unified_data['clean'].median()) / 1000000, 1)

            # Calculate low coverage samples (< 1M final reads)
            stats['low_coverage_count'] = int(len(self.unified_data[self.unified_data['clean'] < 1000000]))

        # Calculate processing efficiency statistics
        if 'raw' in self.unified_data.columns and 'fastp' in self.unified_data.columns:
            trimming_loss = ((self.unified_data['raw'] - self.unified_data['fastp']) / self.unified_data['raw'] * 100)
            stats['median_trimming_loss_pct'] = round(float(trimming_loss.median()), 1)

        if 'fastp' in self.unified_data.columns and 'host_depleted' in self.unified_data.columns:
            host_removal = ((self.unified_data['fastp'] - self.unified_data['host_depleted']) / self.unified_data['fastp'] * 100)
            stats['median_host_removal_pct'] = round(float(host_removal.median()), 1)

        # Use pre-calculated rrna_percent from QC flags for consistency with plot
        if 'rrna_percent' in self.unified_data.columns:
            # Filter out non-numeric values
            rrna_values = self.unified_data['rrna_percent']
            rrna_numeric = rrna_values[rrna_values != 'NA'].astype(float)
            if len(rrna_numeric) > 0:
                stats['median_rrna_removal_pct'] = round(float(rrna_numeric.median()), 1)

        if 'raw' in self.unified_data.columns and 'clean' in self.unified_data.columns:
            total_loss = ((self.unified_data['raw'] - self.unified_data['clean']) / self.unified_data['raw'] * 100)
            stats['median_total_loss_pct'] = round(float(total_loss.median()), 1)

        # Calculate median values for key metrics (convert to native Python types for JSON serialization)
        numeric_cols = ['clean_retention', 'host_percent', 'rrna_percent', 'phix_percent',
                       'vector_percent']

        for col in numeric_cols:
            if col in self.unified_data.columns:
                stats[f'{col}_median'] = float(self.unified_data[col].median()) if pd.notna(self.unified_data[col].median()) else None
                stats[f'{col}_mean'] = float(self.unified_data[col].mean()) if pd.notna(self.unified_data[col].mean()) else None
                stats[f'{col}_std'] = float(self.unified_data[col].std()) if pd.notna(self.unified_data[col].std()) else None

        return stats

    def generate_json_data(self) -> Dict[str, Any]:
        """Generate structured data for HTML/PDF rendering"""
        batch_stats = self.calculate_batch_statistics()

        # Convert DataFrame to dict for JSON serialization
        samples_data = self.unified_data.to_dict('records')

        return {
            'batch_statistics': batch_stats,
            'samples': samples_data,
            'config': self.config,
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }

    def encode_images_as_base64(self, output_dir: str) -> Dict[str, str]:
        """
        Convert PNG images to base64 data URLs for embedding in HTML
        """
        image_files = {
            'contamination_bars': 'contamination_bars.png',
            'contamination_heatmap': 'contamination_heatmap.png',
            'primer_b_heatmap': 'primer_b_heatmap.png'
        }

        base64_images = {}

        for key, filename in image_files.items():
            image_path = Path(output_dir) / filename

            if image_path.exists():
                try:
                    with open(image_path, 'rb') as img_file:
                        img_data = img_file.read()
                        img_base64 = base64.b64encode(img_data).decode('utf-8')
                        base64_images[key] = f"data:image/png;base64,{img_base64}"
                        print(f"Encoded {filename} as base64 ({len(img_data)} bytes)")
                except Exception as e:
                    print(f"Warning: Could not encode {filename} as base64: {e}")
                    base64_images[key] = None
            else:
                print(f"Warning: Image file not found: {image_path}")
                base64_images[key] = None

        return base64_images

    def link_images_externally(self, output_dir: str) -> Dict[str, str]:
        """
        Generate relative paths to image files instead of base64 encoding
        """
        image_files = {
            'contamination_bars': 'contamination_bars.png',
            'contamination_heatmap': 'contamination_heatmap.png',
            'primer_b_heatmap': 'primer_b_heatmap.png'
        }

        image_links = {}

        for key, filename in image_files.items():
            image_path = Path(output_dir) / filename

            if image_path.exists():
                image_links[key] = filename  # Relative path for HTML
                print(f"Linked external image: {filename}")
            else:
                print(f"Warning: Image file not found: {image_path}")
                image_links[key] = None

        return image_links

    def generate_html_reports(self, template_dir: str, output_file: str, output_file_standalone: str = None) -> bool:
        """
        Generate both lightweight and standalone HTML reports
        """
        if not JINJA2_AVAILABLE:
            print("Skipping HTML generation - jinja2 not available")
            return False

        try:
            # Set up Jinja2 environment
            env = Environment(loader=FileSystemLoader(template_dir))
            template = env.get_template('virome_report.html')

            # Generate report data
            report_data = self.generate_json_data()
            output_dir = Path(output_file).parent

            # Generate lightweight version (external images)
            print("Generating lightweight HTML report with external images...")
            image_links = self.link_images_externally(str(output_dir))

            html_content_light = template.render(
                batch_statistics=report_data['batch_statistics'],
                samples=report_data['samples'],
                config=report_data['config'],
                timestamp=report_data['timestamp'],
                report_data=report_data,
                image_links=image_links,  # External image references
                base64_images={},  # Empty for lightweight version
                standalone_version=False
            )

            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html_content_light)

            print(f"Lightweight HTML report generated: {output_file}")

            # Generate standalone version (embedded images) if requested
            if output_file_standalone:
                print("Generating standalone HTML report with embedded images...")
                base64_images = self.encode_images_as_base64(str(output_dir))

                html_content_standalone = template.render(
                    batch_statistics=report_data['batch_statistics'],
                    samples=report_data['samples'],
                    config=report_data['config'],
                    timestamp=report_data['timestamp'],
                    report_data=report_data,
                    image_links={},  # Empty for standalone version
                    base64_images=base64_images,  # Base64-encoded images
                    standalone_version=True
                )

                with open(output_file_standalone, 'w', encoding='utf-8') as f:
                    f.write(html_content_standalone)

                print(f"Standalone HTML report generated: {output_file_standalone}")

            return True

        except Exception as e:
            print(f"Error generating HTML reports: {e}")
            import traceback
            traceback.print_exc()
            return False

    def generate_html_report(self, template_dir: str, output_file: str) -> bool:
        """
        Backward compatibility wrapper - generates lightweight version only
        """
        return self.generate_html_reports(template_dir, output_file)

    def generate(self):
        """Main generation method"""
        print("="*60)
        print("VIROME QC REPORT GENERATOR")
        print("="*60)

        # Phase 1: Load and unify all data
        self.load_all_data()

        # Phase 2: Rank samples by practical metrics
        self.rank_samples_by_metrics()

        # Phase 3: Detect outliers
        self.detect_all_outliers()

        # Phase 4: Generate output data structure
        report_data = self.generate_json_data()

        # Phase 5: Save structured data
        output_json = self.outputs.get('json_data', 'virome_report_data.json')
        with open(output_json, 'w') as f:
            json.dump(report_data, f, indent=2, default=str)

        print(f"\nGenerated structured report data: {output_json}")

        # Phase 6: Generate HTML reports (both versions)
        html_output = self.outputs.get('html_report')
        html_standalone_output = self.outputs.get('html_report_standalone')

        if html_output:
            # Determine template directory relative to script location
            script_dir = Path(__file__).parent
            template_dir = script_dir.parent / 'templates'

            success = self.generate_html_reports(str(template_dir), html_output, html_standalone_output)
            if success:
                print(f"HTML reports generated:")
                print(f"  Lightweight version: {html_output}")
                if html_standalone_output:
                    print(f"  Standalone version: {html_standalone_output}")
            else:
                print("Failed to generate HTML reports")

        # Phase 7: Generate PDF report (placeholder for future implementation)
        pdf_output = self.outputs.get('pdf_report')
        if pdf_output:
            print(f"PDF generation not yet implemented, skipping: {pdf_output}")

        # Print summary
        stats = report_data['batch_statistics']
        print(f"\nBATCH SUMMARY:")
        print(f"  Total samples: {stats['total_samples']}")
        print(f"  Outliers detected: {stats['outlier_count']}")
        print(f"  Median host contamination: {stats.get('host_percent_median', 'N/A'):.1f}%")
        print(f"  Median rRNA contamination: {stats.get('rrna_percent_median', 'N/A'):.1f}%")

        print("\n" + "="*60)
        return report_data

def main():
    """Main entry point for Snakemake integration"""
    # Get inputs and outputs from snakemake
    inputs = {
        'read_counts': snakemake.input.read_counts,
        'contamination_summary': snakemake.input.contamination_summary,
        'qc_metrics': snakemake.input.qc_metrics,
        'primer_b_summary': snakemake.input.primer_b_summary,
        'fastqc_files': snakemake.input.get('fastqc_files', [])
    }

    outputs = {
        'json_data': snakemake.output.json_data,
        'html_report': snakemake.output.get('html_report'),
        'html_report_standalone': snakemake.output.get('html_report_standalone'),
        'pdf_report': snakemake.output.get('pdf_report')
    }

    config = snakemake.config

    # Generate report
    generator = ViromePipelineReportGenerator(inputs, outputs, config)
    report_data = generator.generate()

    return report_data

if __name__ == "__main__":
    main()