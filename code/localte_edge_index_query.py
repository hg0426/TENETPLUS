#!/usr/bin/env python3
"""Query or materialize compact timepoint LocalTE GRNs.

The edge-index GRN format stores one row per timepoint:

    Time, EdgeIndex[], TE[]

and a separate edge catalog:

    EdgeIndex, Source, Target, TE

This helper expands only the requested timepoints / sources / targets to a
regular edge table, avoiding full materialization unless explicitly requested.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from code.path_utils import coerce_input_path, coerce_output_path


def _split_filter(value: Optional[str]) -> Optional[set[str]]:
    if value is None or str(value).strip() == "":
        return None
    return {item.strip() for item in str(value).split(",") if item.strip()}


def _parse_times(value: Optional[str]) -> Optional[set[int]]:
    if value is None or str(value).strip() == "":
        return None
    out: set[int] = set()
    for part in str(value).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            left, right = part.split("-", 1)
            out.update(range(int(left), int(right) + 1))
        else:
            out.add(int(part))
    return out


def _edge_index_filter(catalog: pd.DataFrame, sources: Optional[set[str]], targets: Optional[set[str]], size: int) -> Optional[np.ndarray]:
    mask = pd.Series(True, index=catalog.index)
    if sources is not None:
        mask &= catalog["Source"].astype(str).isin(sources)
    if targets is not None:
        mask &= catalog["Target"].astype(str).isin(targets)
    if bool(mask.all()):
        return None
    allowed = np.zeros(size, dtype=bool)
    allowed[catalog.loc[mask, "EdgeIndex"].astype(np.int64).to_numpy()] = True
    return allowed


def _catalog_arrays(catalog: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edge_index = catalog["EdgeIndex"].astype(np.int64).to_numpy()
    size = int(edge_index.max()) + 1 if edge_index.size else 0
    source_by_edge = np.empty(size, dtype=object)
    target_by_edge = np.empty(size, dtype=object)
    globalte_by_edge = np.empty(size, dtype=np.float32)
    source_by_edge[edge_index] = catalog["Source"].to_numpy(dtype=object)
    target_by_edge[edge_index] = catalog["Target"].to_numpy(dtype=object)
    if "GlobalTE" in catalog.columns:
        global_vals = catalog["GlobalTE"].astype(np.float32).to_numpy()
    elif "TE" in catalog.columns:
        global_vals = catalog["TE"].astype(np.float32).to_numpy()
    else:
        global_vals = np.full(edge_index.size, np.nan, dtype=np.float32)
    globalte_by_edge[edge_index] = global_vals
    return source_by_edge, target_by_edge, globalte_by_edge


def _expand_rows(
    edge_index_path: Path,
    source_by_edge: np.ndarray,
    target_by_edge: np.ndarray,
    globalte_by_edge: np.ndarray,
    *,
    times: Optional[set[int]],
    time_start: Optional[int],
    time_end: Optional[int],
    allowed_edges: Optional[np.ndarray],
    max_rows: Optional[int],
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    row_count = 0
    pf = pq.ParquetFile(edge_index_path)
    for batch in pf.iter_batches(columns=["Time", "EdgeIndex", "TE"]):
        data = batch.to_pydict()
        for time_value, edge_list, te_list in zip(data["Time"], data["EdgeIndex"], data["TE"]):
            time_int = int(time_value)
            if times is not None and time_int not in times:
                continue
            if time_start is not None and time_int < time_start:
                continue
            if time_end is not None and time_int > time_end:
                continue
            if not edge_list:
                continue
            edges = np.asarray(edge_list, dtype=np.int64)
            tes = np.asarray(te_list, dtype=np.float32)
            if allowed_edges is not None:
                keep = allowed_edges[edges]
                if not np.any(keep):
                    continue
                edges = edges[keep]
                tes = tes[keep]
            if max_rows is not None:
                remaining = int(max_rows) - row_count
                if remaining <= 0:
                    break
                if edges.size > remaining:
                    edges = edges[:remaining]
                    tes = tes[:remaining]
            frame = pd.DataFrame(
                {
                    "Time": np.full(edges.size, time_int, dtype=np.int32),
                    "EdgeIndex": edges.astype(np.int32, copy=False),
                    "TE": tes,
                }
            )
            frame.insert(0, "Source", source_by_edge[edges])
            frame.insert(2, "Target", target_by_edge[edges])
            frame["GlobalTE"] = globalte_by_edge[edges]
            frames.append(frame)
            row_count += int(edges.size)
            if max_rows is not None and row_count >= int(max_rows):
                break
        if max_rows is not None and row_count >= int(max_rows):
            break
    if not frames:
        return pd.DataFrame(columns=["Source", "TE", "Target", "Time", "EdgeIndex", "GlobalTE"])
    out = pd.concat(frames, ignore_index=True)
    return out.loc[:, ["Source", "TE", "Target", "Time", "EdgeIndex", "GlobalTE"]]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--edge-index", required=True, help="edge_index parquet from make_GRN_localTE --output_mode edge_index")
    p.add_argument("--catalog", required=True, help="edge catalog parquet written beside the edge_index output")
    p.add_argument("--output", help="Optional output path. If omitted, prints a small preview.")
    p.add_argument("--time", help="Comma-separated timepoints or ranges, e.g. 1,5,10-20")
    p.add_argument("--time-start", type=int)
    p.add_argument("--time-end", type=int)
    p.add_argument("--source", help="Comma-separated source names or ids")
    p.add_argument("--target", help="Comma-separated target names or ids")
    p.add_argument("--max-rows", type=int, help="Optional row cap for previews/debugging")
    p.add_argument("--columns", choices=["full", "legacy"], default="full", help="legacy writes Source,TE,Target,Time only")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    edge_index_path = coerce_input_path(args.edge_index)
    catalog_path = coerce_input_path(args.catalog)
    catalog = pd.read_parquet(catalog_path)
    required = {"EdgeIndex", "Source", "Target"}
    missing = required - set(catalog.columns)
    if missing:
        raise ValueError(f"catalog is missing required columns: {sorted(missing)}")
    keep_cols = ["EdgeIndex", "Source", "Target"]
    if "GlobalTE" in catalog.columns:
        keep_cols.append("GlobalTE")
    elif "TE" in catalog.columns:
        keep_cols.append("TE")
    catalog = catalog.loc[:, keep_cols].copy()
    catalog["EdgeIndex"] = catalog["EdgeIndex"].astype(np.int64)

    sources = _split_filter(args.source)
    targets = _split_filter(args.target)
    source_by_edge, target_by_edge, globalte_by_edge = _catalog_arrays(catalog)
    allowed_edges = _edge_index_filter(catalog, sources, targets, len(source_by_edge))

    out = _expand_rows(
        edge_index_path,
        source_by_edge,
        target_by_edge,
        globalte_by_edge,
        times=_parse_times(args.time),
        time_start=args.time_start,
        time_end=args.time_end,
        allowed_edges=allowed_edges,
        max_rows=args.max_rows,
    )
    if args.columns == "legacy":
        out["Time"] = out["Time"].astype(str)
        out = out.loc[:, ["Source", "TE", "Target", "Time"]]

    if args.output:
        output_path = coerce_output_path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix.lower() == ".csv":
            out.to_csv(output_path, index=False)
        elif output_path.suffix.lower() in {".tsv", ".txt"}:
            out.to_csv(output_path, sep="\t", index=False)
        else:
            out.to_parquet(output_path, index=False)
        print(f"Wrote {len(out):,} rows to {output_path}")
    else:
        print(out.head(20).to_string(index=False))
        print(f"[localte_edge_index_query] rows={len(out):,}")


if __name__ == "__main__":
    main()
