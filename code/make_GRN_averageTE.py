#!/usr/bin/env python
"""
Build a GRN from average TE scores (e.g., TE_result_all.parquet) using the
same statistic as make_GRN_new: global z-score + BH-FDR filtering.

The script accepts Parquet/CSV input with at least columns
  - Source
  - Target
  - TE (average transfer entropy)

Outputs can be written to a single Parquet/CSV file or to a SIF edge list for
compatibility with the legacy pipeline.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.sandbox.stats.multicomp import multipletests


def load_table(path: str, columns: Optional[list[str]] = None) -> pd.DataFrame:
    _, ext = os.path.splitext(path)
    ext = ext.lower()
    try:
        if ext == ".parquet":
            return pd.read_parquet(path, columns=columns)
        if ext in {".csv", ".txt", ".tsv"}:
            sep = "," if ext == ".csv" else "\t"
            return pd.read_csv(path, sep=sep, usecols=columns)
    except Exception as exc:
        raise RuntimeError(f"Failed to load {path}: {exc}")
    raise ValueError(f"Unsupported input format for {path}; expected parquet/csv/tsv")


def write_output(df: pd.DataFrame, path: str, sif: bool = False) -> None:
    if sif:
        with open(path, "w") as handle:
            for row in df.itertuples(index=False):
                handle.write(f"{row.Source}\t{row.TE}\t{row.Target}\n")
        return
    _, ext = os.path.splitext(path)
    ext = ext.lower()
    if ext == ".csv":
        df.to_csv(path, index=False)
    else:
        df.to_parquet(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Average-TE GRN builder (z-score + BH-FDR).")
    parser.add_argument("--input", required=True, help="Input TE table (parquet/csv) with Source, Target, TE columns.")
    parser.add_argument("--fdr", type=float, default=0.05, help="BH-FDR threshold (default 0.05).")
    parser.add_argument("--te_cutoff", type=float, default=0.0, help="Include edges with TE strictly greater than this cutoff (default 0.0).")
    parser.add_argument("--output_file", required=True, help="Output path (parquet/csv/sif).")
    args = parser.parse_args()

    df = load_table(args.input, columns=["Source", "Target", "TE"])
    if df.empty:
        print("Input table is empty; nothing to do.")
        return

    sources = []
    targets = []
    tes = []
    cutoff = float(args.te_cutoff)
    for row in df.itertuples(index=False):
        te_val = float(row.TE)
        if te_val > cutoff:
            sources.append(row.Source)
            targets.append(row.Target)
            tes.append(te_val)

    if len(tes) < 2:
        print("Not enough TE values above cutoff to compute statistics.")
        return

    te_arr = np.asarray(tes, dtype=float)
    std = te_arr.std(ddof=0)
    if np.allclose(std, 0.0):
        print("TE values have ~zero variance; no significant edges.")
        return

    mean = te_arr.mean()
    zscores = (te_arr - mean) / std
    pvals = 1 - norm.cdf(zscores)
    _, qvals, _, _ = multipletests(pvals, alpha=args.fdr, method="fdr_bh")

    rows = []
    for src, tgt, te_val, z, p, q in zip(sources, targets, te_arr, zscores, pvals, qvals):
        if q < args.fdr:
            rows.append({
                "Source": src,
                "Target": tgt,
                "TE": te_val,
                "zscore": z,
                "p_value": p,
                "q_value": q,
            })

    if not rows:
        print("No edges passed the FDR cutoff.")
        return

    df_out = pd.DataFrame(rows)

    out_path = args.output_file
    sif = out_path.lower().endswith(".sif")
    write_output(df_out, out_path, sif=sif)
    print(f"Saved {len(df_out):,} edges to {out_path}")


if __name__ == "__main__":
    main()
