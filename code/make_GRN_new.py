#!/usr/bin/env python3
"""Export a global z-score/BH network from a legacy tab-separated TE matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def export_network(input_path: Path, fdr: float) -> Path:
    if not np.isfinite(fdr) or not 0 < fdr <= 1:
        raise ValueError("FDR cutoff must be in (0, 1].")
    sources, targets, scores = [], [], []
    with input_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError("Expected a header with a row-label column and target IDs.")
        names = header[1:]
        if any(not name.strip() for name in names) or len(set(names)) != len(names):
            raise ValueError("Target IDs must be nonempty and unique.")
        seen = set()
        for row in reader:
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) != len(header):
                raise ValueError(f"Row {reader.line_num}: expected {len(header)} columns, got {len(row)}.")
            source = row[0]
            if not source.strip() or source in seen:
                raise ValueError(f"Row {reader.line_num}: source IDs must be nonempty and unique.")
            seen.add(source)
            try:
                values = np.asarray(row[1:], dtype=np.float64)
            except ValueError as exc:
                raise ValueError(f"Row {reader.line_num}: TE values must be numeric.") from exc
            if not np.isfinite(values).all():
                raise ValueError(f"Row {reader.line_num}: TE values must be finite.")
            for target, value in zip(names, values):
                if value > 0:
                    sources.append(source)
                    targets.append(target)
                    scores.append(float(value))
    values = np.asarray(scores, dtype=np.float64)
    selected = np.zeros(len(values), dtype=bool)
    # A global z-score is undefined without positive variance. Match the
    # Parquet exporter by producing an empty network for these cases.
    if values.size > 1 and np.std(values) > 0:
        z = (values - np.mean(values)) / np.std(values)
        # Preserve the original formula and strict threshold for valid inputs.
        p = 1 - norm.cdf(z)
        q = multipletests(p, alpha=0.05, method="fdr_bh")[1]
        selected = q < fdr
    stem = input_path.stem if input_path.suffix == ".txt" else input_path.name
    output = input_path.with_name(f"{stem}.fdr{fdr}.sif")
    with output.open("w") as handle:
        for i in np.flatnonzero(selected):
            handle.write(f"{sources[i]}\t{scores[i]}\t{targets[i]}\n")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path, help="Headered, tab-separated source-by-target TE matrix.")
    parser.add_argument("fdr", type=float, help="FDR cutoff in (0, 1].")
    args = parser.parse_args()
    print(f"Wrote network to {export_network(args.matrix, args.fdr)}")


if __name__ == "__main__":
    main()
