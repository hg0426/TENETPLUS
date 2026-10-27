#!/usr/bin/env python3
"""Export exact LocalTE chunks on demand from a TENET+ lazy recipe.

This avoids materializing all pair x time LocalTE payloads during the main TE run.
Given a recipe and a selected edge list, it recomputes the exact grouped-state
LocalTE vectors for those edges and writes the familiar chunk parquet layout.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse
from tqdm import tqdm

from . import runTE_for_py_python_batch as te


def _load_gene_names(path: Path) -> list[str]:
    from .feature_utils import read_names, validate_names
    names = read_names(path)
    validate_names(names)
    return names


def _coerce_endpoint(value, name_to_idx: dict[str, int]) -> int:
    if isinstance(value, (np.integer, int)) and not isinstance(value, (bool, np.bool_)):
        return int(value)
    if not isinstance(value, str):
        raise ValueError("Edge endpoints must be feature-name strings or integer Parquet IDs")
    text = value
    if not text:
        raise ValueError("empty edge endpoint")
    if text not in name_to_idx:
        raise KeyError(f"Unknown feature name in edge list: {text}")
    return int(name_to_idx[text])


def _read_edges(path: Path, name_to_idx: dict[str, int], limit: int = 0) -> list[tuple[int, int]]:
    suffix = path.suffix.lower()
    rows = []
    if suffix == ".parquet":
        table = pq.read_table(path, columns=["Source", "Target"])
        data = table.to_pydict()
        rows = list(zip(data["Source"], data["Target"]))
    elif suffix in (".csv", ".tsv", ".txt") and suffix != ".sif":
        sep = "\t" if suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
        if "Source" not in df.columns or "Target" not in df.columns:
            # Fall back to first and last columns for headerless/simple tables.
            df = pd.read_csv(path, sep=sep, header=None, dtype=str, keep_default_na=False)
            rows = list(zip(df.iloc[:, 0], df.iloc[:, -1]))
        else:
            rows = list(zip(df["Source"], df["Target"]))
    else:
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                # .sif is usually: Source interaction Target.  If only two
                # columns are present, treat them as Source Target.
                src = parts[0]
                tgt = parts[-1]
                rows.append((src, tgt))
    edges = []
    seen = set()
    for src, tgt in rows:
        s = _coerce_endpoint(src, name_to_idx)
        t = _coerce_endpoint(tgt, name_to_idx)
        key = (s, t)
        if key in seen:
            continue
        seen.add(key)
        edges.append(key)
        if limit > 0 and len(edges) >= limit:
            break
    return edges


def _load_recipe(recipe_path: Path) -> dict:
    with open(recipe_path, "r", encoding="utf-8") as fh:
        recipe = json.load(fh)
    if recipe.get("format") != "tenet_localte_lazy_recipe_v1":
        raise ValueError(f"Unsupported recipe format: {recipe.get('format')}")
    return recipe


def _prepare_te_globals(recipe: dict) -> tuple[list[str], int, float]:
    matrix_path = Path(recipe["matrix_parquet"])
    gene_path = Path(recipe["gene_names"])
    dense = pd.read_parquet(matrix_path).to_numpy(dtype=np.float64)
    te.GLOBAL_CELL_GENE_DENSE = dense
    te.GLOBAL_CELL_GENE = sparse.csr_matrix(dense.shape, dtype=np.float64)
    names = _load_gene_names(gene_path)
    te.GENE_NAMES = names
    history = int(recipe.get("history_length", 1))
    kernel_width = float(recipe.get("kernel_width", 0.5))
    if history != 1:
        raise ValueError("Lazy grouped-state LocalTE export currently requires history_length=1")
    if dense.ndim != 2 or dense.shape[0] != len(names) or not np.isfinite(dense).all():
        raise ValueError("Recipe matrix must be finite feature-by-cell data matching gene_names")
    from .feature_utils import load_manifest
    previous_output = os.environ.get("TENET_OUTPUT_DIR")
    os.environ["TENET_OUTPUT_DIR"] = str(gene_path.parent)
    try:
        load_manifest(names)
    finally:
        if previous_output is None:
            os.environ.pop("TENET_OUTPUT_DIR", None)
        else:
            os.environ["TENET_OUTPUT_DIR"] = previous_output
    return names, history, kernel_width


def _target_context(target_idx: int, history: int, kernel_width: float):
    dest = te._maybe_subsample(te._get_series(int(target_idx)))
    if dest.size <= int(history):
        return None
    dest_past = dest[:-1].astype(np.float64, copy=False)
    dest_next = dest[1:].astype(np.float64, copy=False)
    epsilon = 1e-9
    kw_past = float(kernel_width) * (np.std(dest_past, ddof=1) + epsilon)
    kw_next = float(kernel_width) * (np.std(dest_next, ddof=1) + epsilon)
    scaled_past = np.divide(dest_past, kw_past, out=np.zeros_like(dest_past), where=(kw_past != 0))
    scaled_next = np.divide(dest_next, kw_next, out=np.zeros_like(dest_next), where=(kw_next != 0))
    target_codes, target_past_values, target_next_values, target_counts = te._target_grouped_state_codes(
        scaled_past,
        scaled_next,
    )
    target_factor_context = te._target_grouped_state_factor_context(
        target_past_values,
        target_next_values,
        target_counts,
    )
    return target_codes, target_past_values, target_next_values, target_factor_context, int(scaled_past.shape[0])


def _compute_localte(source_idx: int, target_ctx, history: int, kernel_width: float) -> tuple[float, np.ndarray]:
    if target_ctx is None:
        return 0.0, np.zeros(0, dtype=np.float32)
    target_codes, target_past_values, target_next_values, target_factor_context, n_obs = target_ctx
    source_values, source_codes = te._get_grouped_source_codes_kernel(int(source_idx), int(history), float(kernel_width))
    if source_codes.shape[0] != int(n_obs):
        return 0.0, np.zeros(int(n_obs), dtype=np.float32)
    result = te._grouped_state_te_local_from_coded(
        source_values,
        source_codes,
        target_codes,
        target_past_values,
        target_next_values,
        target_factor_context=target_factor_context,
        max_unique=int(te.KERNEL_GROUPED_STATE_MAX_UNIQUE),
    )
    if result is None:
        raise RuntimeError(f"LocalTE max_unique fallback is unsupported for source={source_idx}")
    score, values = result
    return float(score), np.asarray(values, dtype=np.float16)


def _flush_chunk(out_dir: Path, chunk_idx: int, rows: list[dict], compression: str) -> None:
    if not rows:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    table = pa.table(
        {
            "Source": pa.array([r["Source"] for r in rows], type=pa.int64()),
            "Target": pa.array([r["Target"] for r in rows], type=pa.int64()),
            "ChunkIndex": pa.array([r["ChunkIndex"] for r in rows], type=pa.int32()),
            "ChunkStart": pa.array([r["ChunkStart"] for r in rows], type=pa.int32()),
            "Values": pa.array([r["Values"] for r in rows], type=pa.list_(pa.float16())),
        }
    )
    pq.write_table(table, out_dir / f"chunk_{chunk_idx:04d}.parquet", compression=compression)


def export(recipe_path: Path, edges_path: Path, output_dir: Path, chunk_size: int, limit_edges: int, compression: str) -> dict:
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise FileExistsError(f"LocalTE export destination must be empty: {output_dir}; use a new output directory")
    recipe = _load_recipe(recipe_path)
    names, history, kernel_width = _prepare_te_globals(recipe)
    name_to_idx = {name: i + 1 for i, name in enumerate(names)}
    edges = _read_edges(edges_path, name_to_idx, limit_edges)
    canonical = pd.read_parquet(recipe["te_result_all"], columns=["Source", "Target", "TE"])
    if canonical.duplicated(["Source", "Target"]).any() or not np.isfinite(canonical.TE).all():
        raise ValueError("Canonical TE table must contain unique pairs and finite scores")
    canonical_scores = canonical.set_index(["Source", "Target"])["TE"]
    if any(s < 1 or t < 1 or s > len(names) or t > len(names) for s, t in edges):
        raise ValueError("Selected endpoints must be valid 1-based feature IDs")
    if any(pair not in canonical_scores.index for pair in edges):
        raise ValueError("Selected edge is absent from canonical TE results")
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped_by_target: dict[int, list[int]] = {}
    for src, tgt in edges:
        grouped_by_target.setdefault(int(tgt), []).append(int(src))

    chunk_buffers: dict[int, list[dict]] = {}
    score_rows = []
    pbar = tqdm(total=len(edges), desc="Lazy LocalTE edges", unit="edge", dynamic_ncols=True)
    for target_idx, sources in grouped_by_target.items():
        ctx = _target_context(int(target_idx), history, kernel_width)
        for source_idx in sources:
            score, values = _compute_localte(int(source_idx), ctx, history, kernel_width)
            original_score = np.float64(canonical_scores.loc[(source_idx, target_idx)])
            if np.float64(score).view(np.uint64) != original_score.view(np.uint64):
                raise ValueError(f"Recipe TE does not match canonical result for {(source_idx, target_idx)}; use matching inputs and kernel settings")
            score_rows.append({"Source": int(source_idx), "Target": int(target_idx), "TE": float(original_score)})
            n_chunks = int(math.ceil(float(values.size) / float(chunk_size))) if values.size else 0
            for chunk_idx in range(n_chunks):
                start = int(chunk_idx) * int(chunk_size)
                end = min(start + int(chunk_size), int(values.size))
                chunk_buffers.setdefault(int(chunk_idx), []).append(
                    {
                        "Source": int(source_idx),
                        "Target": int(target_idx),
                        "ChunkIndex": int(chunk_idx),
                        "ChunkStart": int(start),
                        "Values": np.asarray(values[start:end], dtype=np.float16),
                    }
                )
            pbar.update(1)
    pbar.close()

    for chunk_idx in sorted(chunk_buffers):
        _flush_chunk(output_dir, int(chunk_idx), chunk_buffers[chunk_idx], compression)

    scores = pa.table(
        {
            "Source": pa.array([r["Source"] for r in score_rows], type=pa.int64()),
            "Target": pa.array([r["Target"] for r in score_rows], type=pa.int64()),
            "TE": pa.array([r["TE"] for r in score_rows], type=pa.float64()),
        }
    )
    pq.write_table(scores, output_dir / "TE_result_scores.parquet", compression=compression)
    manifest = {
        "format": "tenet_localte_lazy_export_v1",
        "recipe": str(recipe_path),
        "edges": str(edges_path),
        "output_dir": str(output_dir),
        "edges_exported": len(edges),
        "canonical_TE_bits_verified": True,
        "score_dtype": "float64",
        "local_payload_dtype": "float16",
        "chunk_size": int(chunk_size),
        "chunks": len(chunk_buffers),
        "compression": compression,
    }
    (output_dir / "metadata.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Export exact LocalTE chunks from a lazy TENET+ recipe.")
    parser.add_argument("--recipe", required=True, help="local_te_recipe/recipe.json from a lazy run")
    parser.add_argument("--edges", required=True, help="Selected edges: .sif, parquet, csv, or tsv")
    parser.add_argument("--output_dir", required=True, help="Output LocalTE chunk directory")
    parser.add_argument("--chunk_size", type=int, default=300)
    parser.add_argument("--limit_edges", type=int, default=0, help="Debug: export only first N edges")
    parser.add_argument("--compression", default="snappy", help="Parquet compression for chunks")
    args = parser.parse_args()
    manifest = export(
        Path(args.recipe),
        Path(args.edges),
        Path(args.output_dir),
        max(1, int(args.chunk_size)),
        max(0, int(args.limit_edges)),
        str(args.compression),
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
