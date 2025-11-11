"""
Utility functions for sample detection and management
"""

from pathlib import Path
import sys


def auto_detect_samples(config):
    """
    Auto-detect paired-end samples from a directory based on file patterns.

    Parameters
    ----------
    config : dict
        Pipeline configuration dictionary

    Returns
    -------
    dict
        Dictionary of samples with r1 and r2 paths
        Format: {sample_name: {"r1": path, "r2": path}}
    """
    # If auto-detection is not enabled, return manual samples
    auto_config = config.get("sample_auto_detection", {})
    if not auto_config.get("enabled", False):
        return config.get("samples", {})

    # Get auto-detection settings
    input_dir = Path(auto_config["input_dir"])
    r1_pattern = auto_config.get("r1_pattern", "*_R1.fastq.gz")
    r2_pattern = auto_config.get("r2_pattern", "*_R2.fastq.gz")

    # Validate input directory
    if not input_dir.exists():
        sys.stderr.write(f"ERROR: Sample input directory does not exist: {input_dir}\n")
        sys.exit(1)

    if not input_dir.is_dir():
        sys.stderr.write(f"ERROR: Sample input path is not a directory: {input_dir}\n")
        sys.exit(1)

    # Extract the suffix from patterns (everything after the *)
    if "*" not in r1_pattern or "*" not in r2_pattern:
        sys.stderr.write("ERROR: Sample patterns must contain '*' wildcard\n")
        sys.exit(1)

    r1_suffix = r1_pattern.split("*", 1)[1]
    r2_suffix = r2_pattern.split("*", 1)[1]

    # Find all R1 files matching the pattern
    r1_files = sorted(input_dir.glob(r1_pattern))

    if not r1_files:
        sys.stderr.write(f"WARNING: No R1 files found matching pattern '{r1_pattern}' in {input_dir}\n")
        return {}

    samples = {}

    for r1_path in r1_files:
        # Extract sample name by removing the suffix
        if not r1_path.name.endswith(r1_suffix):
            continue

        sample_name = r1_path.name[:-len(r1_suffix)]

        # Look for corresponding R2 file
        r2_path = input_dir / f"{sample_name}{r2_suffix}"

        if r2_path.exists():
            samples[sample_name] = {
                "r1": str(r1_path),
                "r2": str(r2_path)
            }
            print(f"Detected sample: {sample_name}", file=sys.stderr)
        else:
            sys.stderr.write(f"WARNING: No R2 file found for sample '{sample_name}', skipping...\n")

    if not samples:
        sys.stderr.write(f"ERROR: No valid paired-end samples found in {input_dir}\n")
        sys.exit(1)

    print(f"\nAuto-detected {len(samples)} sample(s)", file=sys.stderr)

    return samples


def get_samples(config):
    """
    Get samples either from manual config or auto-detection.

    This is the main entry point for sample management.

    Parameters
    ----------
    config : dict
        Pipeline configuration dictionary

    Returns
    -------
    dict
        Dictionary of samples with r1 and r2 paths
    """
    return auto_detect_samples(config)
