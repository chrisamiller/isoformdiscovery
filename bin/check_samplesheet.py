#!/usr/bin/env python3
"""
check_samplesheet.py — Validate a transcript-assembly pipeline samplesheet CSV.

Usage:
    check_samplesheet.py <samplesheet.csv> <validated.csv>

Arguments:
    samplesheet.csv  Input CSV with required columns: group, sample, bam
    validated.csv    Output path for the validated (and optionally normalised) CSV

Validation rules:
    - Required columns: group, sample, bam
    - group must be 'g1' or 'g2'
    - sample must be unique across all rows
    - bam must be an absolute path to an existing file

Exits non-zero on the first validation error, printing a descriptive message to stderr.
"""

import csv
import os
import sys

REQUIRED_COLS = {"group", "sample", "bam"}
VALID_GROUPS  = {"g1", "g2"}


def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    input_path  = sys.argv[1]
    output_path = sys.argv[2]

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh)

        if reader.fieldnames is None:
            error("Samplesheet appears to be empty.")

        missing = REQUIRED_COLS - set(reader.fieldnames)
        if missing:
            error(f"Samplesheet is missing required columns: {', '.join(sorted(missing))}")

        rows     = []
        seen     = set()

        for i, row in enumerate(reader, start=2):  # row 1 is header
            group  = row["group"].strip()
            sample = row["sample"].strip()
            bam    = row["bam"].strip()

            if group not in VALID_GROUPS:
                error(f"Row {i}: 'group' must be 'g1' or 'g2', got '{group}'")

            if sample in seen:
                error(f"Row {i}: duplicate sample name '{sample}'")
            seen.add(sample)

            if not os.path.isabs(bam):
                error(f"Row {i}: 'bam' must be an absolute path, got '{bam}'")

            if not os.path.isfile(bam):
                error(f"Row {i}: BAM file does not exist: '{bam}'")

            rows.append({"group": group, "sample": sample, "bam": bam})

    if not rows:
        error("Samplesheet contains no data rows.")

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["group", "sample", "bam"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Validated {len(rows)} samples.")


if __name__ == "__main__":
    main()
