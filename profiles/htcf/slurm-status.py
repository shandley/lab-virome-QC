#!/usr/bin/env python3
"""
SLURM job status checker for Snakemake

Returns job status for a given SLURM job ID
Used by Snakemake to monitor cluster jobs
"""

import sys
import subprocess
import time

jobid = sys.argv[1]

# Try to get job status with sacct
for i in range(20):
    try:
        # Check job status using sacct
        sacct_result = subprocess.run(
            ["sacct", "-j", jobid, "--format=State", "--noheader", "--parsable2"],
            capture_output=True,
            text=True,
            timeout=10
        )

        status = sacct_result.stdout.strip().split("\n")[0]

        # Map SLURM states to Snakemake expectations
        if status == "COMPLETED":
            print("success")
            sys.exit(0)
        elif status in ["RUNNING", "PENDING", "COMPLETING"]:
            print("running")
            sys.exit(0)
        elif status in ["FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL"]:
            print("failed")
            sys.exit(0)
        else:
            # Job status not yet available, wait and retry
            time.sleep(1)

    except subprocess.TimeoutExpired:
        time.sleep(1)
        continue
    except Exception as e:
        # If sacct fails, assume job is still running
        print("running")
        sys.exit(0)

# If we couldn't determine status after retries, assume running
print("running")
sys.exit(0)
