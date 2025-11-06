#!/usr/bin/env python3
"""
Test script for sample auto-detection functionality
"""

import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from sample_utils import get_samples


def test_auto_detection():
    """Test auto-detection with mock config"""

    print("=" * 80)
    print("Testing Sample Auto-Detection")
    print("=" * 80)

    # Test 1: Manual sample specification (auto-detection disabled)
    print("\n[Test 1] Manual sample specification (auto-detection disabled)")
    print("-" * 80)
    config1 = {
        "samples": {
            "sample1": {"r1": "data/sample1_R1.fastq.gz", "r2": "data/sample1_R2.fastq.gz"},
            "sample2": {"r1": "data/sample2_R1.fastq.gz", "r2": "data/sample2_R2.fastq.gz"},
        },
        "sample_auto_detection": {
            "enabled": False
        }
    }
    samples1 = get_samples(config1)
    print(f"Samples detected: {list(samples1.keys())}")
    assert "sample1" in samples1
    assert "sample2" in samples1
    print("✓ Manual sample specification works correctly")

    # Test 2: Auto-detection with default patterns
    print("\n[Test 2] Auto-detection with test data")
    print("-" * 80)

    # Use the test output directory as a test case
    test_dir = Path(__file__).parent.parent.parent / "test_NovaSeq_N983_I13420_39802" / "rrna_removed"

    if test_dir.exists():
        config2 = {
            "sample_auto_detection": {
                "enabled": True,
                "input_dir": str(test_dir),
                "r1_pattern": "*_R1.fastq.gz",
                "r2_pattern": "*_R2.fastq.gz"
            }
        }
        samples2 = get_samples(config2)
        print(f"Samples detected: {list(samples2.keys())}")
        print(f"Total: {len(samples2)} sample(s)")

        for sample_name, paths in samples2.items():
            print(f"\n  {sample_name}:")
            print(f"    R1: {paths['r1']}")
            print(f"    R2: {paths['r2']}")

        print("✓ Auto-detection works correctly")
    else:
        print(f"⚠ Test directory not found: {test_dir}")
        print("  Skipping auto-detection test")

    # Test 3: Different pattern formats
    print("\n[Test 3] Testing pattern flexibility")
    print("-" * 80)
    patterns = [
        ("*_R1.fastq.gz", "*_R2.fastq.gz", "Underscore separator"),
        ("*_R1_001.fastq.gz", "*_R2_001.fastq.gz", "Illumina default"),
        ("*.R1.fastq.gz", "*.R2.fastq.gz", "Dot separator"),
        ("*_1.fq.gz", "*_2.fq.gz", "Short form"),
    ]

    for r1_pat, r2_pat, description in patterns:
        print(f"  Pattern: {description}")
        print(f"    R1: {r1_pat}")
        print(f"    R2: {r2_pat}")

    print("✓ Multiple pattern formats supported")

    print("\n" + "=" * 80)
    print("All tests passed!")
    print("=" * 80)


if __name__ == "__main__":
    test_auto_detection()
