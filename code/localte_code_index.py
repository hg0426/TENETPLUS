#!/usr/bin/env python3
"""Build reusable source/target time-code indexes for virtual LocalTE export.

The virtual LocalTE array stores per-edge state tables.  To materialize a full
LocalTE vector, we also need the source code sequence and target code sequence
for every timepoint.  Recomputing those sequences for every edge is wasteful
when many edges reuse the same sources and targets.  This module stores the
per-feature code sequences once, like an animation/game-engine cache.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from tqdm import tqdm

from . import runTE_for_py_python_batch as te
from .localte_lazy_export import _load_recipe, _prepare_te_globals, _target_context


def _read_virtual_feature_sets(virtual_dir: Path, limit_edges: int = 0) -> tuple[set[int], set[int], int]:
    parts = sorted((virtual_dir / "parts").rglob("part_*.parquet"))
    if not parts:
        raise FileNotFoundError(f"No virtual LocalTE parts found under {virtual_dir}/parts")
    sources: set[int] = set()
    targets: set[int] = set()
    count = 0
    dataset = ds.dataset([str(p) for p in parts], format="parquet")
    for batch in dataset.to_batches(columns=["Source", "Target"]):
        data = batch.to_pydict()
        for src, tgt in zip(data["Source"], data["Target"]):
            sources.add(int(src))
            targets.add(int(tgt))
            count += 1
            if int(limit_edges) > 0 and count >= int(limit_edges):
                return sources, targets, count
    return sources, targets, count


def _write_source_codes(path: Path, rows: list[dict], compression: str) -> None:
    table = pa.table(
        {
            "Feature": pa.array([r["Feature"] for r in rows], type=pa.int64()),
            "NObs": pa.array([r["NObs"] for r in rows], type=pa.int32()),
            "NSourceStates": pa.array([r["NSourceStates"] for r in rows], type=pa.int32()),
            "SourceCodes": pa.array([r["SourceCodes"] for r in rows], type=pa.list_(pa.int32())),
        }
    )
    pq.write_table(table, path, compression=None if compression == "none" else compression)


def _write_target_codes(path: Path, rows: list[dict], compression: str) -> None:
    table = pa.table(
        {
            "Feature": pa.array([r["Feature"] for r in rows], type=pa.int64()),
            "NObs": pa.array([r["NObs"] for r in rows], type=pa.int32()),
            "NTargetStates": pa.array([r["NTargetStates"] for r in rows], type=pa.int32()),
            "TargetCodes": pa.array([r["TargetCodes"] for r in rows], type=pa.list_(pa.int32())),
        }
    )
    pq.write_table(table, path, compression=None if compression == "none" else compression)


def _source_index_row(feature: int, history: int, kernel_width: float) -> dict:
    _values, codes = te._get_grouped_source_codes_kernel(int(feature), int(history), float(kernel_width))
    return {
        "Feature": int(feature),
        "NObs": int(codes.size),
        "NSourceStates": int(_values.size),
        "SourceCodes": np.asarray(codes, dtype=np.int32),
    }


def _target_index_row(feature: int, history: int, kernel_width: float) -> dict:
    ctx = _target_context(int(feature), int(history), float(kernel_width))
    if ctx is None:
        return {
            "Feature": int(feature),
            "NObs": 0,
            "NTargetStates": 0,
            "TargetCodes": np.zeros(0, dtype=np.int32),
        }
    target_codes = np.asarray(ctx[0], dtype=np.int32)
    return {
        "Feature": int(feature),
        "NObs": int(target_codes.size),
        "NTargetStates": int(ctx[1].size),
        "TargetCodes": target_codes,
    }


def _parallel_rows(features: list[int], fn, history: int, kernel_width: float, workers: int, desc: str) -> list[dict]:
    if int(workers) <= 1:
        return [
            fn(int(feature), int(history), float(kernel_width))
            for feature in tqdm(features, desc=desc, unit="feature", dynamic_ncols=True)
        ]
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, int(workers))) as executor:
        futures = [
            executor.submit(fn, int(feature), int(history), float(kernel_width))
            for feature in features
        ]
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc=desc, unit="feature", dynamic_ncols=True):
            rows.append(future.result())
    rows.sort(key=lambda row: int(row["Feature"]))
    return rows


def build(virtual_dir: Path, output_dir: Path, limit_edges: int, compression: str, workers: int = 1) -> dict:
    metadata = json.loads((virtual_dir / "metadata.json").read_text(encoding="utf-8"))
    recipe_path = Path(metadata["recipe"])
    if recipe_path.exists():
        recipe = _load_recipe(recipe_path)
    else:
        run_dir = virtual_dir.parent
        recipe = {
            "format": "tenet_localte_lazy_recipe_v1",
            "work_dir": str(run_dir),
            "matrix_parquet": metadata.get("matrix_parquet", str(run_dir / "cell_gene_trsps.parquet")),
            "gene_names": metadata.get("gene_names", str(run_dir / "gene_names")),
            "time_index_map": metadata.get("time_index_map", str(run_dir / "time_index_map.parquet")),
            "te_result_all": metadata.get("edges", str(run_dir / "TE_result_all.parquet")),
            "virtual_array": str(virtual_dir),
            "history_length": int(metadata.get("history_length", 1)),
            "kernel_width": float(metadata.get("kernel_width", 0.5)),
        }
    _names, history, kernel_width = _prepare_te_globals(recipe)
    sources, targets, edge_count = _read_virtual_feature_sets(virtual_dir, int(limit_edges))
    output_dir.mkdir(parents=True, exist_ok=True)

    start = time.time()
    source_rows = _parallel_rows(sorted(sources), _source_index_row, int(history), float(kernel_width), int(workers), "Build source code index")
    target_rows = _parallel_rows(sorted(targets), _target_index_row, int(history), float(kernel_width), int(workers), "Build target code index")

    _write_source_codes(output_dir / "source_codes.parquet", source_rows, compression)
    _write_target_codes(output_dir / "target_codes.parquet", target_rows, compression)
    payload = {
        "format": "tenet_localte_code_index_v1",
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "virtual_dir": str(virtual_dir),
        "recipe": str(recipe_path),
        "history_length": int(history),
        "kernel_width": float(kernel_width),
        "edges_scanned": int(edge_count),
        "sources": int(len(sources)),
        "targets": int(len(targets)),
        "compression": str(compression),
        "workers": int(workers),
        "elapsed_seconds": float(time.time() - start),
    }
    (output_dir / "metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build reusable code indexes for virtual LocalTE export.")
    parser.add_argument("--virtual_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--limit_edges", type=int, default=0)
    parser.add_argument("--compression", choices=("zstd", "snappy", "gzip", "none"), default="zstd")
    parser.add_argument("--workers", type=int, default=1, help="Thread workers for feature code generation.")
    args = parser.parse_args(argv)
    manifest = build(Path(args.virtual_dir), Path(args.output_dir), int(args.limit_edges), str(args.compression), max(1, int(args.workers)))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
