#!/usr/bin/env python3
"""Build and export a virtual LocalTE array.

A virtual LocalTE array stores each edge as a compact local-state table instead
of materializing one value per pseudotime observation.  The full LocalTE vector
is reconstructed exactly by combining:

  source_codes[source, time] * n_target_states + target_codes[target, time]
  -> local_state_table[edge]

This is intended for analyses that need exact LocalTE values but do not want to
write every edge x time vector during the main TE run.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from tqdm import tqdm

from . import runTE_for_py_python_batch as te
from .localte_lazy_export import (
    _flush_chunk,
    _load_gene_names,
    _load_recipe,
    _prepare_te_globals,
    _read_edges,
    _target_context,
)


_EXPORT_CHUNK_SCHEMA = pa.schema(
    [
        ("Source", pa.int64()),
        ("Target", pa.int64()),
        ("ChunkIndex", pa.int32()),
        ("ChunkStart", pa.int32()),
        ("Values", pa.list_(pa.float16())),
    ]
)

_EXPORT_SCORE_SCHEMA = pa.schema(
    [
        ("Source", pa.int64()),
        ("Target", pa.int64()),
        ("TE", pa.float64()),
    ]
)


def _value_arrow_type(dtype: str):
    dtype = str(dtype).lower()
    if dtype in ("float32", "fp32"):
        return pa.float32(), np.float32
    return pa.float16(), np.float16


def _state_code_type():
    bits = str(os.getenv("LOCAL_TE_VIRTUAL_STATE_CODE_BITS", "32")).strip().lower()
    if bits in ("32", "int32", "i32"):
        return pa.int32(), np.int32, "int32"
    return pa.int64(), np.int64, "int64"


def _coerce_state_codes(state_codes) -> np.ndarray:
    _code_type, code_np_dtype, _code_dtype_str = _state_code_type()
    arr = np.asarray(state_codes)
    if code_np_dtype is np.int32 and arr.size:
        max_code = int(np.max(arr))
        min_code = int(np.min(arr))
        if min_code < np.iinfo(np.int32).min or max_code > np.iinfo(np.int32).max:
            raise OverflowError(
                "StateCodes exceed int32 range; rerun with LOCAL_TE_VIRTUAL_STATE_CODE_BITS=64."
            )
    return np.asarray(arr, dtype=code_np_dtype)


def _compute_localte_state_table(source_idx: int, target_ctx, history: int, kernel_width: float, value_np_dtype):
    """Return TE and a compact observed-state LocalTE table for one edge."""
    if target_ctx is None:
        return 0.0, 0, 0, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_np_dtype)
    target_codes, target_past_values, target_next_values, target_factor_context, n_obs = target_ctx
    source_values, source_codes = te._get_grouped_source_codes_kernel(int(source_idx), int(history), float(kernel_width))
    if source_codes.shape[0] != int(n_obs):
        return 0.0, 0, int(n_obs), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_np_dtype)

    n_target_states = int(target_past_values.size)
    if source_values.size == 0 or n_target_states <= 0:
        return 0.0, n_target_states, int(n_obs), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_np_dtype)

    n_source_states = int(source_values.size)
    code_space = n_source_states * n_target_states
    combined_codes = np.asarray(source_codes, dtype=np.int64) * n_target_states + np.asarray(target_codes, dtype=np.int64)
    use_bincount = str(te.KERNEL_GROUPED_STATE_BINCOUNT).lower() not in ("off", "0", "false", "no")
    if use_bincount and int(te.KERNEL_GROUPED_STATE_BINCOUNT_MAX) > 0:
        use_bincount = code_space <= int(te.KERNEL_GROUPED_STATE_BINCOUNT_MAX)
    if use_bincount:
        counts = np.bincount(combined_codes, minlength=code_space)
        combined_unique = np.flatnonzero(counts)
        weights = counts[combined_unique]
    else:
        combined_unique, weights = np.unique(combined_codes, return_counts=True)
    if combined_unique.size == 0:
        return 0.0, n_target_states, int(n_obs), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_np_dtype)
    if int(te.KERNEL_GROUPED_STATE_MAX_UNIQUE) > 0 and combined_unique.size > int(te.KERNEL_GROUPED_STATE_MAX_UNIQUE):
        raise RuntimeError(f"combined unique states exceed max_unique for source={source_idx}")

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    w = weights.astype(np.int64, copy=False)

    if target_factor_context is not None:
        requested = str(te.KERNEL_GROUPED_STATE_FACTORIZED).lower()
        if (
            requested in ("off", "0", "false", "no")
            or (requested in ("auto", "") and int(source_codes.size) < int(te.KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS))
            or (int(te.KERNEL_GROUPED_STATE_FACTORIZED_MAX) > 0 and code_space > int(te.KERNEL_GROUPED_STATE_FACTORIZED_MAX))
        ):
            target_factor_context = None

    if target_factor_context is not None:
        W = np.zeros(code_space, dtype=np.int64)
        W[combined_unique] = w
        W = W.reshape(n_source_states, n_target_states)
        target_past_left = target_factor_context[0]
        target_past_right = target_factor_context[1]
        count_past_by_target = target_factor_context[3]
        count_next_past_by_target = target_factor_context[4]
        use_source_prefix = str(te.KERNEL_GROUPED_STATE_SOURCE_PREFIX).lower() not in ("off", "0", "false", "no")
        if use_source_prefix:
            source_weighted = te._source_interval_weighted_counts(source_values, W)
        else:
            source_ok = (np.abs(source_values[:, None] - source_values[None, :]) <= 1.0).astype(np.int64, copy=False)
            source_weighted = source_ok @ W
        count_past_source_matrix = te._target_past_interval_weighted_counts(
            target_past_left,
            target_past_right,
            source_weighted,
        )
        count_next_past_source_matrix = te._target_past_next_interval_weighted_counts(
            target_factor_context,
            source_weighted,
        )
        count_past_source = count_past_source_matrix[src_code, tgt_code]
        count_next_past_source = count_next_past_source_matrix[src_code, tgt_code]
        count_past = count_past_by_target[tgt_code]
        count_next_past = count_next_past_by_target[tgt_code]
    else:
        us = source_values[src_code]
        up = target_past_values[tgt_code]
        uy = target_next_values[tgt_code]
        past_ok = np.abs(up[:, None] - up[None, :]) <= 1.0
        next_ok = np.abs(uy[:, None] - uy[None, :]) <= 1.0
        source_ok = np.abs(us[:, None] - us[None, :]) <= 1.0
        count_past = past_ok @ w
        past_next_ok = past_ok & next_ok
        count_next_past = past_next_ok @ w
        count_past_source = (past_ok & source_ok) @ w
        count_next_past_source = (past_next_ok & source_ok) @ w

    local_by_unique = np.zeros(combined_unique.size, dtype=np.float64)
    valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
    if np.any(valid):
        with np.errstate(divide="ignore", invalid="ignore"):
            numerator = count_next_past_source[valid] / count_past_source[valid]
            denominator = count_next_past[valid] / count_past[valid]
            log_terms = np.log(numerator / denominator)
            log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
        local_by_unique[valid] = log_terms / np.log(2.0)
    te_score = float(np.sum(local_by_unique * w) / int(n_obs))
    return (
        te_score,
        n_target_states,
        int(n_obs),
        np.asarray(combined_unique, dtype=np.int64),
        np.asarray(local_by_unique, dtype=value_np_dtype),
    )


def _schema(value_type) -> pa.Schema:
    code_type, _code_np_dtype, _code_dtype_str = _state_code_type()
    return pa.schema(
        [
            ("Source", pa.int64()),
            ("Target", pa.int64()),
            ("TE", pa.float64()),
            ("NTargetStates", pa.int32()),
            ("NObs", pa.int32()),
            ("NLocalStates", pa.int32()),
            ("StateCodes", pa.list_(code_type)),
            ("LocalValues", pa.list_(value_type)),
        ]
    )


def _write_part(path: Path, rows: list[dict], value_type, compression: str) -> None:
    table = pa.table(
        {
            "Source": pa.array([r["Source"] for r in rows], type=pa.int64()),
            "Target": pa.array([r["Target"] for r in rows], type=pa.int64()),
            "TE": pa.array([r["TE"] for r in rows], type=pa.float64()),
            "NTargetStates": pa.array([r["NTargetStates"] for r in rows], type=pa.int32()),
            "NObs": pa.array([r["NObs"] for r in rows], type=pa.int32()),
            "NLocalStates": pa.array([r["NLocalStates"] for r in rows], type=pa.int32()),
            "StateCodes": pa.array([r["StateCodes"] for r in rows], type=pa.list_(_state_code_type()[0])),
            "LocalValues": pa.array([r["LocalValues"] for r in rows], type=pa.list_(value_type)),
        },
        schema=_schema(value_type),
    )
    pq.write_table(table, path, compression=compression)


def build(recipe_path: Path, edges_path: Path, output_dir: Path, limit_edges: int, buffer_edges: int, value_dtype: str, compression: str) -> dict:
    recipe = _load_recipe(recipe_path)
    names, history, kernel_width = _prepare_te_globals(recipe)
    name_to_idx = {name: i + 1 for i, name in enumerate(names)}
    edges = _read_edges(edges_path, name_to_idx, int(limit_edges))
    output_dir.mkdir(parents=True, exist_ok=True)
    parts_dir = output_dir / "parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    value_type, value_np_dtype = _value_arrow_type(value_dtype)

    grouped_by_target: dict[int, list[int]] = {}
    for src, tgt in edges:
        grouped_by_target.setdefault(int(tgt), []).append(int(src))

    rows: list[dict] = []
    part_idx = 0
    total_states = 0
    start = time.time()
    pbar = tqdm(total=len(edges), desc="Build virtual LocalTE", unit="edge", dynamic_ncols=True)
    for target_idx, sources in grouped_by_target.items():
        ctx = _target_context(int(target_idx), history, kernel_width)
        for source_idx in sources:
            score, n_target_states, n_obs, state_codes, local_values = _compute_localte_state_table(
                int(source_idx),
                ctx,
                history,
                kernel_width,
                value_np_dtype,
            )
            total_states += int(state_codes.size)
            rows.append(
                {
                    "Source": int(source_idx),
                    "Target": int(target_idx),
                    "TE": float(score),
                    "NTargetStates": int(n_target_states),
                    "NObs": int(n_obs),
                    "NLocalStates": int(state_codes.size),
                    "StateCodes": _coerce_state_codes(state_codes),
                    "LocalValues": local_values,
                }
            )
            if len(rows) >= int(buffer_edges):
                _write_part(parts_dir / f"part_{part_idx:06d}.parquet", rows, value_type, compression)
                part_idx += 1
                rows = []
            pbar.update(1)
    pbar.close()
    if rows:
        _write_part(parts_dir / f"part_{part_idx:06d}.parquet", rows, value_type, compression)
        part_idx += 1

    manifest = {
        "format": "tenet_localte_virtual_array_v1",
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "recipe": str(recipe_path),
        "edges": str(edges_path),
        "parts_dir": "parts",
        "edges_stored": len(edges),
        "state_values_stored": int(total_states),
        "mean_states_per_edge": float(total_states / len(edges)) if edges else 0.0,
        "value_dtype": str(value_dtype),
        "compression": str(compression),
        "buffer_edges": int(buffer_edges),
        "elapsed_seconds": time.time() - start,
        "note": "Reconstruct exact LocalTE vectors by combining StateCodes/LocalValues with source and target time codes.",
    }
    (output_dir / "metadata.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def _read_edges_grouped(path: Path, name_to_idx: dict[str, int], limit: int = 0) -> tuple[dict[int, list[int]], int]:
    """Read edges grouped by target without materializing row dicts.

    Parquet edge tables are common for full TE outputs, so keep that path
    batch-oriented. Text/SIF inputs are usually selected networks and can reuse
    the simpler reader.
    """
    grouped: dict[int, list[int]] = {}
    count = 0
    if path.suffix.lower() == ".parquet":
        parquet = pq.ParquetFile(path)
        for batch in parquet.iter_batches(columns=["Source", "Target"], batch_size=262144):
            data = batch.to_pydict()
            for src, tgt in zip(data["Source"], data["Target"]):
                grouped.setdefault(int(tgt), []).append(int(src))
                count += 1
                if int(limit) > 0 and count >= int(limit):
                    return grouped, count
        return grouped, count
    for src, tgt in _read_edges(path, name_to_idx, int(limit)):
        grouped.setdefault(int(tgt), []).append(int(src))
        count += 1
    return grouped, count


def _virtual_build_target_worker(args_tuple):
    (
        worker_id,
        target_items,
        output_dir,
        history,
        kernel_width,
        buffer_edges,
        value_dtype,
        compression,
    ) = args_tuple
    value_type, value_np_dtype = _value_arrow_type(value_dtype)
    worker_dir = Path(output_dir) / "parts" / f"worker_{int(worker_id):03d}"
    worker_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    part_idx = 0
    edge_count = 0
    state_count = 0
    for target_idx, sources in target_items:
        ctx = _target_context(int(target_idx), int(history), float(kernel_width))
        for source_idx in sources:
            score, n_target_states, n_obs, state_codes, local_values = _compute_localte_state_table(
                int(source_idx),
                ctx,
                int(history),
                float(kernel_width),
                value_np_dtype,
            )
            edge_count += 1
            state_count += int(state_codes.size)
            rows.append(
                {
                    "Source": int(source_idx),
                    "Target": int(target_idx),
                    "TE": float(score),
                    "NTargetStates": int(n_target_states),
                    "NObs": int(n_obs),
                    "NLocalStates": int(state_codes.size),
                    "StateCodes": _coerce_state_codes(state_codes),
                    "LocalValues": local_values,
                }
            )
            if len(rows) >= int(buffer_edges):
                _write_part(worker_dir / f"part_{part_idx:06d}.parquet", rows, value_type, compression)
                part_idx += 1
                rows = []
    if rows:
        _write_part(worker_dir / f"part_{part_idx:06d}.parquet", rows, value_type, compression)
        part_idx += 1
    return {
        "worker": int(worker_id),
        "edges": int(edge_count),
        "states": int(state_count),
        "parts": int(part_idx),
    }


def _split_target_items(grouped_by_target: dict[int, list[int]], workers: int) -> list[list[tuple[int, list[int]]]]:
    bins: list[list[tuple[int, list[int]]]] = [[] for _ in range(max(1, int(workers)))]
    loads = [0 for _ in bins]
    for target_idx, sources in sorted(grouped_by_target.items(), key=lambda kv: len(kv[1]), reverse=True):
        bin_idx = min(range(len(bins)), key=lambda i: loads[i])
        bins[bin_idx].append((int(target_idx), [int(s) for s in sources]))
        loads[bin_idx] += len(sources)
    return [b for b in bins if b]


def build_parallel(
    recipe_path: Path,
    edges_path: Path,
    output_dir: Path,
    limit_edges: int,
    buffer_edges: int,
    value_dtype: str,
    compression: str,
    workers: int,
) -> dict:
    recipe = _load_recipe(recipe_path)
    names, history, kernel_width = _prepare_te_globals(recipe)
    name_to_idx = {name: i + 1 for i, name in enumerate(names)}
    grouped_by_target, edge_count = _read_edges_grouped(edges_path, name_to_idx, int(limit_edges))
    output_dir.mkdir(parents=True, exist_ok=True)
    parts_dir = output_dir / "parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    n_workers = max(1, min(int(workers), max(1, len(grouped_by_target))))
    if n_workers <= 1:
        return build(recipe_path, edges_path, output_dir, limit_edges, buffer_edges, value_dtype, compression)

    target_bins = _split_target_items(grouped_by_target, n_workers)
    tasks = [
        (
            worker_id,
            target_items,
            str(output_dir),
            int(history),
            float(kernel_width),
            int(buffer_edges),
            str(value_dtype),
            str(compression),
        )
        for worker_id, target_items in enumerate(target_bins)
    ]
    results = []
    pbar = tqdm(total=int(edge_count), desc="Build virtual LocalTE", unit="edge", dynamic_ncols=True)
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(processes=len(tasks)) as pool:
        for result in pool.imap_unordered(_virtual_build_target_worker, tasks):
            results.append(result)
            pbar.update(int(result["edges"]))
    pbar.close()
    total_states = sum(int(r["states"]) for r in results)
    total_edges = sum(int(r["edges"]) for r in results)
    manifest = {
        "format": "tenet_localte_virtual_array_v1",
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "recipe": str(recipe_path),
        "edges": str(edges_path),
        "parts_dir": "parts",
        "edges_stored": int(total_edges),
        "state_values_stored": int(total_states),
        "mean_states_per_edge": float(total_states / total_edges) if total_edges else 0.0,
        "value_dtype": str(value_dtype),
        "compression": str(compression),
        "buffer_edges": int(buffer_edges),
        "workers": int(len(tasks)),
        "worker_results": results,
        "elapsed_seconds": time.time() - start,
        "note": "Reconstruct exact LocalTE vectors by combining StateCodes/LocalValues with source and target time codes.",
    }
    (output_dir / "metadata.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def _code_index_table(path: Path, wanted_features: set[int] | None):
    if wanted_features is None:
        return pq.read_table(path).to_pydict()
    wanted = sorted(int(v) for v in wanted_features if int(v) > 0)
    if not wanted:
        return {"Feature": []}
    dataset = ds.dataset(str(path), format="parquet")
    return dataset.to_table(filter=ds.field("Feature").isin(wanted)).to_pydict()


def _load_code_index(
    code_index_dir: Path | None,
    source_filter: set[int] | None = None,
    target_filter: set[int] | None = None,
):
    if code_index_dir is None:
        return None
    code_index_dir = Path(code_index_dir)
    source_path = code_index_dir / "source_codes.parquet"
    target_path = code_index_dir / "target_codes.parquet"
    if not source_path.exists() or not target_path.exists():
        return None
    source_data = _code_index_table(source_path, source_filter)
    target_data = _code_index_table(target_path, target_filter)
    return {
        "sources": {
            int(feature): np.asarray(codes, dtype=np.int64)
            for feature, codes in zip(source_data["Feature"], source_data["SourceCodes"])
        },
        "targets": {
            int(feature): (
                np.asarray(codes, dtype=np.int64),
                int(n_target_states),
                int(n_obs),
            )
            for feature, codes, n_target_states, n_obs in zip(
                target_data["Feature"],
                target_data["TargetCodes"],
                target_data["NTargetStates"],
                target_data["NObs"],
            )
        },
    }


def _collect_virtual_feature_sets(
    virtual_dir: Path,
    edges_path: Path | None,
    name_to_idx: dict[str, int],
    limit_edges: int,
) -> tuple[set[int] | None, set[int] | None]:
    if edges_path is not None:
        edges = _read_edges(edges_path, name_to_idx, int(limit_edges))
        return {int(src) for src, _tgt in edges}, {int(tgt) for _src, tgt in edges}
    if int(limit_edges) <= 0:
        return None, None
    parts = sorted((virtual_dir / "parts").rglob("part_*.parquet"))
    sources: set[int] = set()
    targets: set[int] = set()
    seen = 0
    for batch in ds.dataset([str(p) for p in parts], format="parquet").to_batches(columns=["Source", "Target"]):
        data = batch.to_pydict()
        for src, tgt in zip(data["Source"], data["Target"]):
            sources.add(int(src))
            targets.add(int(tgt))
            seen += 1
            if seen >= int(limit_edges):
                return sources, targets
    return sources, targets


def _expand_virtual_row(row: dict, history: int, kernel_width: float, code_index=None) -> np.ndarray:
    source_idx = int(row["Source"])
    target_idx = int(row["Target"])
    n_target_states = int(row["NTargetStates"])
    n_obs = int(row["NObs"])
    if n_target_states < 0:
        return np.asarray(row["LocalValues"], dtype=np.float16)
    if n_obs <= 0 or n_target_states <= 0:
        return np.zeros(max(0, n_obs), dtype=np.float16)
    source_codes = None
    target_codes = None
    if code_index is not None:
        source_codes = code_index.get("sources", {}).get(source_idx)
        target_item = code_index.get("targets", {}).get(target_idx)
        if target_item is not None:
            target_codes, indexed_target_states, _indexed_n_obs = target_item
            if int(indexed_target_states) != int(n_target_states):
                target_codes = None
    if target_codes is None:
        target_ctx = _target_context(target_idx, history, kernel_width)
        if target_ctx is None:
            return np.zeros(max(0, n_obs), dtype=np.float16)
        target_codes = np.asarray(target_ctx[0], dtype=np.int64)
    if source_codes is None:
        _source_values, source_codes = te._get_grouped_source_codes_kernel(source_idx, history, kernel_width)
    if source_codes.shape[0] != target_codes.shape[0]:
        return np.zeros(max(0, target_codes.shape[0]), dtype=np.float16)
    combined_codes = np.asarray(source_codes, dtype=np.int64) * n_target_states + target_codes
    state_codes = np.asarray(row["StateCodes"], dtype=np.int64)
    # Keep the stored precision through lookup.  When the virtual table is
    # stored as float32 this matches the legacy path, which computes LocalTE as
    # float32 and only then casts materialized chunk values to float16.
    local_values = np.asarray(row["LocalValues"])
    if state_codes.size == 0:
        return np.zeros(combined_codes.shape[0], dtype=np.float16)
    positions = np.searchsorted(state_codes, combined_codes)
    if np.any(positions >= state_codes.size) or np.any(state_codes[positions] != combined_codes):
        raise RuntimeError(f"virtual state lookup failed for {source_idx}->{target_idx}")
    return np.asarray(local_values[positions], dtype=np.float16)


def _chunk_rows_to_table(rows: list[dict]) -> pa.Table:
    return pa.table(
        {
            "Source": pa.array([r["Source"] for r in rows], type=pa.int64()),
            "Target": pa.array([r["Target"] for r in rows], type=pa.int64()),
            "ChunkIndex": pa.array([r["ChunkIndex"] for r in rows], type=pa.int32()),
            "ChunkStart": pa.array([r["ChunkStart"] for r in rows], type=pa.int32()),
            "Values": pa.array([r["Values"] for r in rows], type=pa.list_(pa.float16())),
        },
        schema=_EXPORT_CHUNK_SCHEMA,
    )


def _score_rows_to_table(rows: list[dict]) -> pa.Table:
    return pa.table(
        {
            "Source": pa.array([r["Source"] for r in rows], type=pa.int64()),
            "Target": pa.array([r["Target"] for r in rows], type=pa.int64()),
            "TE": pa.array([r["TE"] for r in rows], type=pa.float64()),
        },
        schema=_EXPORT_SCORE_SCHEMA,
    )


def _flush_worker_chunk_rows(writers: dict[int, pq.ParquetWriter], part_dir: Path, chunk_idx: int, rows: list[dict], compression: str) -> None:
    if not rows:
        return
    writer = writers.get(int(chunk_idx))
    if writer is None:
        path = part_dir / f"chunk_{int(chunk_idx):04d}.parquet"
        writer = pq.ParquetWriter(path, _EXPORT_CHUNK_SCHEMA, compression=compression)
        writers[int(chunk_idx)] = writer
    writer.write_table(_chunk_rows_to_table(rows))
    rows.clear()


def _flush_worker_score_rows(score_writer: pq.ParquetWriter | None, score_path: Path, rows: list[dict], compression: str) -> pq.ParquetWriter | None:
    if not rows:
        return score_writer
    if score_writer is None:
        score_writer = pq.ParquetWriter(score_path, _EXPORT_SCORE_SCHEMA, compression=compression)
    score_writer.write_table(_score_rows_to_table(rows))
    rows.clear()
    return score_writer


def _part_feature_sets(part_path: Path) -> tuple[set[int], set[int]]:
    sources: set[int] = set()
    targets: set[int] = set()
    for batch in pq.ParquetFile(part_path).iter_batches(columns=["Source", "Target"]):
        data = batch.to_pydict()
        sources.update(int(v) for v in data["Source"])
        targets.update(int(v) for v in data["Target"])
    return sources, targets


def _export_part_task(args):
    (
        part_idx,
        part_path_str,
        part_out_dir_str,
        code_index_dir_str,
        history,
        kernel_width,
        chunk_size,
        compression,
        buffer_edges,
    ) = args
    part_path = Path(part_path_str)
    part_out_dir = Path(part_out_dir_str)
    part_out_dir.mkdir(parents=True, exist_ok=True)
    sources, targets = _part_feature_sets(part_path)
    code_index = _load_code_index(
        Path(code_index_dir_str),
        source_filter=sources,
        target_filter=targets,
    )
    if code_index is None:
        raise RuntimeError("Parallel virtual export requires a code index.")
    chunk_buffers: dict[int, list[dict]] = {}
    chunk_writers: dict[int, pq.ParquetWriter] = {}
    score_rows: list[dict] = []
    score_path = part_out_dir / "scores.parquet"
    score_writer = None
    exported = 0
    chunk_count = 0
    try:
        for batch in pq.ParquetFile(part_path).iter_batches():
            data = batch.to_pydict()
            for i in range(batch.num_rows):
                row = {k: data[k][i] for k in data}
                values = _expand_virtual_row(row, int(history), float(kernel_width), code_index=code_index)
                score_rows.append({"Source": int(row["Source"]), "Target": int(row["Target"]), "TE": float(row["TE"])})
                n_chunks = int(math.ceil(float(values.size) / float(chunk_size))) if values.size else 0
                for chunk_idx in range(n_chunks):
                    start = int(chunk_idx) * int(chunk_size)
                    end = min(start + int(chunk_size), int(values.size))
                    rows = chunk_buffers.setdefault(int(chunk_idx), [])
                    rows.append(
                        {
                            "Source": int(row["Source"]),
                            "Target": int(row["Target"]),
                            "ChunkIndex": int(chunk_idx),
                            "ChunkStart": int(start),
                            "Values": np.asarray(values[start:end], dtype=np.float16),
                        }
                    )
                    if len(rows) >= int(buffer_edges):
                        _flush_worker_chunk_rows(chunk_writers, part_out_dir, int(chunk_idx), rows, compression)
                exported += 1
                if len(score_rows) >= int(buffer_edges):
                    score_writer = _flush_worker_score_rows(score_writer, score_path, score_rows, compression)
        for chunk_idx, rows in list(chunk_buffers.items()):
            _flush_worker_chunk_rows(chunk_writers, part_out_dir, int(chunk_idx), rows, compression)
        score_writer = _flush_worker_score_rows(score_writer, score_path, score_rows, compression)
    finally:
        if score_writer is not None:
            score_writer.close()
        for writer in chunk_writers.values():
            writer.close()
    chunk_count = len(list(part_out_dir.glob("chunk_*.parquet")))
    return {"part_idx": int(part_idx), "edges": int(exported), "chunks": int(chunk_count)}


def _merge_parquet_files(files: list[Path], output_path: Path, schema: pa.Schema, compression: str) -> int:
    writer = None
    rows = 0
    try:
        for path in files:
            pf = pq.ParquetFile(path)
            for batch in pf.iter_batches():
                if writer is None:
                    writer = pq.ParquetWriter(output_path, schema, compression=compression)
                writer.write_table(pa.Table.from_batches([batch], schema=schema))
                rows += int(batch.num_rows)
    finally:
        if writer is not None:
            writer.close()
    return rows


def _merge_chunk_task(args):
    chunk_idx, files, output_dir, compression = args
    output_path = Path(output_dir) / f"chunk_{int(chunk_idx):04d}.parquet"
    rows = _merge_parquet_files([Path(f) for f in files], output_path, _EXPORT_CHUNK_SCHEMA, compression)
    return int(chunk_idx), int(rows)


def _read_virtual_rows(
    virtual_dir: Path,
    edges_path: Path | None,
    name_to_idx: dict[str, int],
    limit_edges: int,
    source_filter: set[int] | None = None,
    target_filter: set[int] | None = None,
):
    parts = sorted((virtual_dir / "parts").rglob("part_*.parquet"))
    if not parts:
        return
    if edges_path is None:
        for batch in ds.dataset([str(p) for p in parts], format="parquet").to_batches():
            data = batch.to_pydict()
            for i in range(batch.num_rows):
                yield {k: data[k][i] for k in data}
        return
    wanted = set(_read_edges(edges_path, name_to_idx, int(limit_edges)))
    if not wanted:
        return
    seen = 0
    dataset = ds.dataset([str(p) for p in parts], format="parquet")
    scan_filter = None
    if target_filter:
        target_values = sorted(int(x) for x in target_filter)
        scan_filter = ds.field("Target").isin(target_values)
    if source_filter:
        source_values = sorted(int(x) for x in source_filter)
        source_expr = ds.field("Source").isin(source_values)
        scan_filter = source_expr if scan_filter is None else (scan_filter & source_expr)
    for batch in dataset.to_batches(filter=scan_filter):
        data = batch.to_pydict()
        for i in range(batch.num_rows):
            key = (int(data["Source"][i]), int(data["Target"][i]))
            if key not in wanted:
                continue
            yield {k: data[k][i] for k in data}
            seen += 1
            if int(limit_edges) > 0 and seen >= int(limit_edges):
                return


def _load_virtual_export_recipe(virtual_dir: Path, code_index_dir: Path | None) -> tuple[dict, Path | None]:
    metadata = json.loads((virtual_dir / "metadata.json").read_text(encoding="utf-8"))
    recipe_path = Path(metadata["recipe"])
    if recipe_path.exists():
        recipe = _load_recipe(recipe_path)
    else:
        # Direct runTE invocations may create a virtual array without going
        # through the shell wrapper that writes local_te_recipe/recipe.json.
        # Keep export self-contained by reconstructing the minimal exact recipe
        # from metadata and files next to the virtual array.
        run_dir = virtual_dir.parent
        history = metadata.get("history_length")
        kernel_width = metadata.get("kernel_width", 0.5)
        if history is None and code_index_dir is not None and (code_index_dir / "metadata.json").exists():
            code_meta = json.loads((code_index_dir / "metadata.json").read_text(encoding="utf-8"))
            history = code_meta.get("history_length", 1)
            kernel_width = code_meta.get("kernel_width", kernel_width)
        recipe = {
            "format": "tenet_localte_lazy_recipe_v1",
            "work_dir": str(run_dir),
            "matrix_parquet": metadata.get("matrix_parquet", str(run_dir / "cell_gene_trsps.parquet")),
            "gene_names": metadata.get("gene_names", str(run_dir / "gene_names")),
            "time_index_map": metadata.get("time_index_map", str(run_dir / "time_index_map.parquet")),
            "all_pairs": metadata.get("all_pairs"),
            "te_result_all": metadata.get("edges", str(run_dir / "TE_result_all.parquet")),
            "virtual_array": str(virtual_dir),
            "code_index": str(code_index_dir) if code_index_dir is not None else metadata.get("code_index"),
            "history_length": int(history or 1),
            "kernel_width": float(kernel_width),
            "note": "Auto-reconstructed from virtual array metadata because recipe.json was not present.",
        }
    if code_index_dir is None and recipe.get("code_index"):
        candidate = Path(recipe["code_index"])
        if (candidate / "source_codes.parquet").exists() and (candidate / "target_codes.parquet").exists():
            code_index_dir = candidate
    return recipe, code_index_dir


def export_parallel(
    virtual_dir: Path,
    output_dir: Path,
    chunk_size: int,
    compression: str,
    code_index_dir: Path,
    workers: int,
    merge_workers: int,
    buffer_edges: int,
) -> dict:
    recipe, code_index_dir = _load_virtual_export_recipe(virtual_dir, code_index_dir)
    if code_index_dir is None:
        raise RuntimeError("Parallel virtual export requires --code_index_dir.")
    parts = sorted((virtual_dir / "parts").rglob("part_*.parquet"))
    if not parts:
        raise FileNotFoundError(f"No virtual array part files found under {virtual_dir / 'parts'}")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir / "_parallel_parts"
    temp_dir.mkdir(parents=True, exist_ok=True)
    tasks = [
        (
            idx,
            str(part),
            str(temp_dir / f"part_{idx:04d}"),
            str(code_index_dir),
            int(recipe.get("history_length", 1)),
            float(recipe.get("kernel_width", 0.5)),
            int(chunk_size),
            str(compression),
            max(1, int(buffer_edges)),
        )
        for idx, part in enumerate(parts)
    ]
    start = time.time()
    total_edges = 0
    with multiprocessing.Pool(processes=max(1, int(workers))) as pool:
        iterator = pool.imap_unordered(_export_part_task, tasks)
        for result in tqdm(iterator, total=len(tasks), desc="Virtual export parts", unit="part", dynamic_ncols=True):
            total_edges += int(result["edges"])

    chunk_files: dict[int, list[Path]] = {}
    score_files: list[Path] = []
    for part_dir in sorted(temp_dir.glob("part_*")):
        score_path = part_dir / "scores.parquet"
        if score_path.exists():
            score_files.append(score_path)
        for chunk_path in sorted(part_dir.glob("chunk_*.parquet")):
            chunk_idx = int(chunk_path.stem.split("_")[-1])
            chunk_files.setdefault(chunk_idx, []).append(chunk_path)

    merge_tasks = [
        (chunk_idx, [str(p) for p in files], str(output_dir), str(compression))
        for chunk_idx, files in sorted(chunk_files.items())
    ]
    chunk_rows: dict[int, int] = {}
    if merge_tasks:
        with multiprocessing.Pool(processes=max(1, int(merge_workers))) as pool:
            iterator = pool.imap_unordered(_merge_chunk_task, merge_tasks)
            for chunk_idx, rows in tqdm(iterator, total=len(merge_tasks), desc="Merge virtual chunks", unit="chunk", dynamic_ncols=True):
                chunk_rows[int(chunk_idx)] = int(rows)
    _merge_parquet_files(score_files, output_dir / "TE_result_scores.parquet", _EXPORT_SCORE_SCHEMA, str(compression))
    shutil.rmtree(temp_dir, ignore_errors=True)
    elapsed = time.time() - start
    manifest = {
        "format": "tenet_localte_virtual_export_v1",
        "virtual_dir": str(virtual_dir),
        "edges": None,
        "output_dir": str(output_dir),
        "edges_exported": int(total_edges),
        "chunk_size": int(chunk_size),
        "chunks": len(chunk_files),
        "compression": str(compression),
        "code_index": str(code_index_dir),
        "workers": int(workers),
        "merge_workers": int(merge_workers),
        "elapsed_seconds": float(elapsed),
        "chunk_rows": {str(k): int(v) for k, v in sorted(chunk_rows.items())},
    }
    (output_dir / "metadata.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def export(
    virtual_dir: Path,
    output_dir: Path,
    edges_path: Path | None,
    chunk_size: int,
    limit_edges: int,
    compression: str,
    code_index_dir: Path | None = None,
    workers: int = 1,
    merge_workers: int | None = None,
    buffer_edges: int = 5000,
) -> dict:
    recipe, code_index_dir = _load_virtual_export_recipe(virtual_dir, code_index_dir)
    if int(workers) > 1 and edges_path is None and int(limit_edges) <= 0 and code_index_dir is not None:
        return export_parallel(
            virtual_dir,
            output_dir,
            int(chunk_size),
            str(compression),
            Path(code_index_dir),
            max(1, int(workers)),
            max(1, int(merge_workers if merge_workers is not None else workers)),
            max(1, int(buffer_edges)),
        )
    history = int(recipe.get("history_length", 1))
    kernel_width = float(recipe.get("kernel_width", 0.5))
    name_to_idx: dict[str, int] = {}
    source_filter = None
    target_filter = None
    if code_index_dir is not None and edges_path is not None:
        names = _load_gene_names(Path(recipe["gene_names"]))
        name_to_idx = {name: i + 1 for i, name in enumerate(names)}
        source_filter, target_filter = _collect_virtual_feature_sets(
            virtual_dir,
            edges_path,
            name_to_idx,
            int(limit_edges),
        )
    elif code_index_dir is not None and int(limit_edges) > 0:
        source_filter, target_filter = _collect_virtual_feature_sets(
            virtual_dir,
            None,
            name_to_idx,
            int(limit_edges),
        )
    code_index = _load_code_index(code_index_dir, source_filter=source_filter, target_filter=target_filter)
    if code_index is None:
        names, history, kernel_width = _prepare_te_globals(recipe)
        name_to_idx = {name: i + 1 for i, name in enumerate(names)}
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_buffers: dict[int, list[dict]] = {}
    score_rows = []
    exported = 0
    rows_iter = _read_virtual_rows(
        virtual_dir,
        edges_path,
        name_to_idx,
        int(limit_edges),
        source_filter=source_filter,
        target_filter=target_filter,
    )
    pbar = tqdm(desc="Export virtual LocalTE", unit="edge", dynamic_ncols=True)
    for row in rows_iter:
        values = _expand_virtual_row(row, history, kernel_width, code_index=code_index)
        score_rows.append({"Source": int(row["Source"]), "Target": int(row["Target"]), "TE": float(row["TE"])})
        n_chunks = int(math.ceil(float(values.size) / float(chunk_size))) if values.size else 0
        for chunk_idx in range(n_chunks):
            start = int(chunk_idx) * int(chunk_size)
            end = min(start + int(chunk_size), int(values.size))
            chunk_buffers.setdefault(int(chunk_idx), []).append(
                {
                    "Source": int(row["Source"]),
                    "Target": int(row["Target"]),
                    "ChunkIndex": int(chunk_idx),
                    "ChunkStart": int(start),
                    "Values": np.asarray(values[start:end], dtype=np.float16),
                }
            )
        exported += 1
        pbar.update(1)
        if int(limit_edges) > 0 and exported >= int(limit_edges):
            break
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
        "format": "tenet_localte_virtual_export_v1",
        "virtual_dir": str(virtual_dir),
        "edges": str(edges_path) if edges_path else None,
        "output_dir": str(output_dir),
        "edges_exported": int(exported),
        "chunk_size": int(chunk_size),
        "chunks": len(chunk_buffers),
        "compression": str(compression),
        "code_index": str(code_index_dir) if code_index is not None and code_index_dir else None,
    }
    (output_dir / "metadata.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Build or export a virtual LocalTE array.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="Build virtual LocalTE state tables from a lazy recipe and edge list.")
    b.add_argument("--recipe", required=True)
    b.add_argument("--edges", required=True)
    b.add_argument("--output_dir", required=True)
    b.add_argument("--limit_edges", type=int, default=0)
    b.add_argument("--buffer_edges", type=int, default=5000)
    b.add_argument("--value_dtype", default="float32", choices=["float16", "float32"])
    b.add_argument("--compression", default="zstd")
    b.add_argument("--workers", type=int, default=1, help="Parallel target workers for large edge sets.")

    e = sub.add_parser("export", help="Materialize LocalTE chunks from a virtual LocalTE array.")
    e.add_argument("--virtual_dir", required=True)
    e.add_argument("--output_dir", required=True)
    e.add_argument("--edges", default=None, help="Optional selected edges to export.")
    e.add_argument("--chunk_size", type=int, default=300)
    e.add_argument("--limit_edges", type=int, default=0)
    e.add_argument("--compression", default="snappy")
    e.add_argument("--code_index_dir", default=None, help="Optional precomputed source/target code index.")
    e.add_argument("--workers", type=int, default=1, help="Parallel export workers for full virtual export (requires --code_index_dir).")
    e.add_argument("--merge_workers", type=int, default=None, help="Parallel chunk merge workers (default: same as --workers).")
    e.add_argument("--buffer_edges", type=int, default=5000, help="Rows buffered per worker/chunk before writing.")

    args = parser.parse_args()
    if args.cmd == "build":
        manifest = build_parallel(
            Path(args.recipe),
            Path(args.edges),
            Path(args.output_dir),
            int(args.limit_edges),
            max(1, int(args.buffer_edges)),
            str(args.value_dtype),
            str(args.compression),
            max(1, int(args.workers)),
        )
    else:
        manifest = export(
            Path(args.virtual_dir),
            Path(args.output_dir),
            Path(args.edges) if args.edges else None,
            max(1, int(args.chunk_size)),
            max(0, int(args.limit_edges)),
            str(args.compression),
            Path(args.code_index_dir) if args.code_index_dir else None,
            max(1, int(args.workers)),
            int(args.merge_workers) if args.merge_workers is not None else None,
            max(1, int(args.buffer_edges)),
        )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
