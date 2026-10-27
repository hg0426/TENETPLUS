#!/usr/bin/env python3
"""Apply TENET indirect trimming independently to each LocalTE timepoint network.

Input is the aggregated output from ``make_GRN_localTE.py``:

    Source, TE, Target, Time

The script streams the input by ``Time`` and applies the same indirect-edge
rule as ``trim_indirect.py`` within each timepoint network.  It writes one
trimmed parquet table and, optionally, a per-timepoint source outdegree table.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from time import perf_counter
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

try:  # module mode: python -m code.trim_indirect_localte_times
    from .trim_indirect import trim_indirect
except Exception:  # script mode: python code/trim_indirect_localte_times.py
    from trim_indirect import trim_indirect  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Trim indirect edges separately for each timepoint in a LocalTE "
            "timepoint network parquet/csv."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Input parquet/csv from make_GRN_localTE.")
    parser.add_argument("--output", required=True, type=Path, help="Trimmed output parquet path.")
    parser.add_argument(
        "--outdegree-output",
        type=Path,
        default=None,
        help="Optional parquet path for Time, Source, Outdegree after trimming.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=None,
        help="Optional TSV summary path with per-timepoint kept/removed counts.",
    )
    parser.add_argument("--threshold", type=float, default=-0.01, help="Indirect trim threshold.")
    parser.add_argument("--source-col", default="Source")
    parser.add_argument("--target-col", default="Target")
    parser.add_argument("--te-col", default="TE")
    parser.add_argument("--time-col", default="Time")
    parser.add_argument("--batch-size", type=int, default=262144, help="Parquet/csv rows per read batch.")
    parser.add_argument("--compression", default="zstd", help="Parquet compression.")
    parser.add_argument(
        "--allow-noncontiguous-time",
        action="store_true",
        help=(
            "Do not fail if the same Time appears again after being flushed. "
            "Use only if input has already been manually partitioned; otherwise "
            "trimming will be applied to partial networks."
        ),
    )
    return parser.parse_args()


def _input_batches(path: Path, columns: list[str], batch_size: int) -> Iterable[pd.DataFrame]:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(columns=columns, batch_size=batch_size, use_threads=True):
            if batch.num_rows:
                yield batch.to_pandas()
        return
    if suffix in {".csv", ".tsv", ".txt"}:
        sep = "\t" if suffix in {".tsv", ".txt"} else ","
        for chunk in pd.read_csv(path, sep=sep, usecols=columns, chunksize=batch_size):
            if not chunk.empty:
                yield chunk
        return
    raise ValueError(f"Unsupported input suffix for {path}; use parquet/csv/tsv/txt.")


def _trim_time_df(
    df: pd.DataFrame,
    source_col: str,
    target_col: str,
    te_col: str,
    threshold: float,
) -> tuple[pd.DataFrame, dict[str, int]]:
    source_order: list[str] = []
    targets_by_source: dict[str, list[str]] = {}
    tes_by_source: dict[str, list[float]] = {}

    sources = df[source_col].astype(str).to_numpy()
    targets = df[target_col].astype(str).to_numpy()
    tes = df[te_col].astype(float).to_numpy()

    for source, target, te in zip(sources, targets, tes):
        if source not in targets_by_source:
            source_order.append(source)
            targets_by_source[source] = []
            tes_by_source[source] = []
        targets_by_source[source].append(target)
        tes_by_source[source].append(float(te))

    indirect_by_source = trim_indirect(source_order, targets_by_source, tes_by_source, threshold)

    keep_mask = np.ones(len(df), dtype=bool)
    cursor_by_source: dict[str, int] = {source: 0 for source in source_order}
    for row_idx, source in enumerate(sources):
        pos = cursor_by_source[source]
        cursor_by_source[source] = pos + 1
        if indirect_by_source[source][pos]:
            keep_mask[row_idx] = False

    trimmed = df.loc[keep_mask].copy()
    summary = {
        "sources": len(source_order),
        "edges": int(len(df)),
        "kept": int(keep_mask.sum()),
        "removed": int((~keep_mask).sum()),
    }
    return trimmed, summary


def _outdegree_table(df: pd.DataFrame, source_col: str, time_col: str) -> pa.Table:
    if df.empty:
        return pa.table(
            {
                time_col: pa.array([], type=pa.string()),
                source_col: pa.array([], type=pa.string()),
                "Outdegree": pa.array([], type=pa.int64()),
            }
        )
    time_value = str(df[time_col].iat[0])
    counts = Counter(df[source_col].astype(str).to_numpy())
    sources = sorted(counts)
    return pa.table(
        {
            time_col: pa.array([time_value] * len(sources), type=pa.string()),
            source_col: pa.array(sources, type=pa.string()),
            "Outdegree": pa.array([counts[source] for source in sources], type=pa.int64()),
        }
    )


def _normalise_for_arrow(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.loc[:, cols].copy()
    out[cols[0]] = out[cols[0]].astype(str)
    out[cols[1]] = out[cols[1]].astype(float)
    out[cols[2]] = out[cols[2]].astype(str)
    out[cols[3]] = out[cols[3]].astype(str)
    return out


def main() -> None:
    args = parse_args()
    input_path = args.input
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.outdegree_output:
        args.outdegree_output.parent.mkdir(parents=True, exist_ok=True)
    if args.summary_output:
        args.summary_output.parent.mkdir(parents=True, exist_ok=True)

    columns = [args.source_col, args.te_col, args.target_col, args.time_col]
    writer: pq.ParquetWriter | None = None
    degree_writer: pq.ParquetWriter | None = None
    pending_time: object | None = None
    pending_parts: list[pd.DataFrame] = []
    flushed_times: set[str] = set()
    summaries: list[dict[str, object]] = []
    total_in = 0
    total_kept = 0
    total_removed = 0
    t0 = perf_counter()

    def flush_pending() -> None:
        nonlocal writer, degree_writer, pending_time, pending_parts, total_kept, total_removed
        if pending_time is None or not pending_parts:
            return
        time_key = str(pending_time)
        df_time = pd.concat(pending_parts, ignore_index=True)
        trimmed, summary = _trim_time_df(
            df_time,
            args.source_col,
            args.target_col,
            args.te_col,
            float(args.threshold),
        )
        summary_row = {"Time": time_key, **summary}
        summaries.append(summary_row)
        total_kept += int(summary["kept"])
        total_removed += int(summary["removed"])

        if not trimmed.empty:
            table = pa.Table.from_pandas(
                _normalise_for_arrow(trimmed, columns),
                preserve_index=False,
            )
            if writer is None:
                writer = pq.ParquetWriter(args.output, table.schema, compression=args.compression)
            writer.write_table(table)

            if args.outdegree_output:
                degree_table = _outdegree_table(trimmed, args.source_col, args.time_col)
                if degree_writer is None:
                    degree_writer = pq.ParquetWriter(
                        args.outdegree_output,
                        degree_table.schema,
                        compression=args.compression,
                    )
                degree_writer.write_table(degree_table)

        flushed_times.add(time_key)
        if len(summaries) == 1 or len(summaries) % 100 == 0:
            elapsed = perf_counter() - t0
            print(
                f"[trim-localTE] times={len(summaries):,} "
                f"edges={total_in:,} kept={total_kept:,} removed={total_removed:,} "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )
        pending_time = None
        pending_parts = []

    try:
        for batch_df in _input_batches(input_path, columns, max(1, int(args.batch_size))):
            total_in += len(batch_df)
            for time_value, group in batch_df.groupby(args.time_col, sort=False):
                time_key = str(time_value)
                if pending_time is None:
                    pending_time = time_value
                if str(pending_time) != time_key:
                    flush_pending()
                    if time_key in flushed_times and not args.allow_noncontiguous_time:
                        raise ValueError(
                            f"Input is not contiguous by {args.time_col!r}; time {time_key!r} "
                            "appeared again after it was already processed. Sort or partition "
                            "by Time before trimming."
                        )
                    pending_time = time_value
                pending_parts.append(group.loc[:, columns].copy())
        flush_pending()
    finally:
        if writer is not None:
            writer.close()
        if degree_writer is not None:
            degree_writer.close()

    if args.summary_output:
        with args.summary_output.open("w", encoding="utf-8", newline="") as handle:
            fieldnames = ["Time", "sources", "edges", "kept", "removed"]
            out = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            out.writeheader()
            out.writerows(summaries)

    elapsed = perf_counter() - t0
    print(
        "trim_indirect_localte_times completed: "
        f"input={input_path} output={args.output} "
        f"times={len(summaries)} edges={total_in} kept={total_kept} "
        f"removed={total_removed} seconds={elapsed:.3f}"
    )


if __name__ == "__main__":
    main()
