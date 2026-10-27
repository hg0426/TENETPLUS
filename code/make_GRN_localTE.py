#!/usr/bin/env python3
"""
Build per-timepoint GRNs from compressed LocalTE outputs (float16 bytes + len/dtype).
Applies the standard TENET GRN-FDR logic independently at each time index: filter TE>cutoff, compute
z-score using population variance, derive one-tailed p-values, apply BH-FDR, and keep
edges with q < alpha.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import os
import shutil
import tempfile
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.sandbox.stats.multicomp import multipletests
import pyarrow as pa
import pyarrow.parquet as pq

try:
    from numba import njit
except Exception:  # pragma: no cover - optional acceleration
    njit = None

from code.path_utils import (
    coerce_input_path,
    coerce_output_path,
    ensure_output_subdir,
    resolve_output,
)

try:
    from code.trim_indirect import trim_indirect as _trim_indirect_edges
except Exception:
    _trim_indirect_edges = None

try:
    from tqdm import tqdm as _tqdm
except Exception:  # pragma: no cover - progress bars are optional
    _tqdm = None


if njit is not None:
    @njit(cache=True)
    def _update_hist_numba(hist, raw_u16, mat_f32, cutoff):
        n_rows, n_cols = raw_u16.shape
        for row_idx in range(n_rows):
            for col_idx in range(n_cols):
                if mat_f32[row_idx, col_idx] > cutoff:
                    hist[col_idx, raw_u16[row_idx, col_idx]] += 1
else:
    _update_hist_numba = None


def _human_int(value: int) -> str:
    return f"{int(value):,}"


def _human_bytes(value: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    amount = float(value)
    for unit in units:
        if amount < 1024.0 or unit == units[-1]:
            return f"{amount:.1f} {unit}"
        amount /= 1024.0
    return f"{amount:.1f} TiB"


def _chunk_index_from_name(path: Path) -> Optional[int]:
    stem = path.stem
    if not stem.startswith("chunk_"):
        return None
    try:
        return int(stem.split("_")[-1])
    except Exception:
        return None


def _discover_chunk_files(chunk_dir: Path) -> Dict[int, List[Path]]:
    """Find top-level chunks or direct-stream worker chunk parts.

    Direct LocalTE streaming writes:
      <chunk_dir>/parts_wXX/chunk_0000.parquet
    while older exporters write:
      <chunk_dir>/chunk_0000.parquet

    Return a mapping from chunk index to one or more parquet files.
    """
    found: Dict[int, List[Path]] = {}
    for path in sorted(chunk_dir.glob("chunk_*.parquet")):
        idx = _chunk_index_from_name(path)
        if idx is not None:
            found.setdefault(idx, []).append(path)
    if found:
        return found
    for path in sorted(chunk_dir.glob("parts_w*/chunk_*.parquet")):
        idx = _chunk_index_from_name(path)
        if idx is not None:
            found.setdefault(idx, []).append(path)
    return found


def _merge_chunk_parts_for_processing(chunk_map: Dict[int, List[Path]], temp_parent: Path) -> Tuple[Path, List[Tuple[int, Path]]]:
    """Merge worker-part chunks into temporary top-level chunks for legacy processors."""
    merged_dir = Path(tempfile.mkdtemp(prefix="localte_chunk_parts_merged_", dir=str(temp_parent)))
    merged: List[Tuple[int, Path]] = []
    for chunk_idx, files in sorted(chunk_map.items()):
        if len(files) == 1:
            merged.append((int(chunk_idx), files[0]))
            continue
        out_path = merged_dir / f"chunk_{int(chunk_idx):04d}.parquet"
        writer: Optional[pq.ParquetWriter] = None
        try:
            for part_path in files:
                pf = pq.ParquetFile(part_path)
                for batch in pf.iter_batches():
                    if writer is None:
                        writer = pq.ParquetWriter(out_path, batch.schema, compression="snappy")
                    writer.write_table(pa.Table.from_batches([batch]))
        finally:
            if writer is not None:
                writer.close()
        merged.append((int(chunk_idx), out_path))
    return merged_dir, merged


def _env_flag(name: str, default: str = "off") -> bool:
    flag = os.environ.get(name, default).strip().lower()
    return flag in {"1", "true", "yes", "on"}


def _progress_iter(iterable, total: Optional[int], desc: str, unit: str):
    if _tqdm is None:
        return iterable
    return _tqdm(iterable, total=total, desc=desc, unit=unit, dynamic_ncols=True)


def load_time_labels(path: Optional[str], history_length: int, required_len: int) -> List[str]:
    if not path:
        return [str(i) for i in range(required_len)]
    labels: List[str] = []
    labels_path = coerce_input_path(path)
    with open(labels_path, encoding="utf-8") as handle:
        for line in handle:
            labels.append(line.strip())
    if len(labels) < required_len + history_length:
        raise ValueError(
            f"Time-label file {labels_path} has {len(labels)} entries, "
            f"need >= {required_len + history_length}."
        )
    trimmed = labels[history_length: history_length + required_len]
    return [str(lab) for lab in trimmed]


def decode_local_matrix(df: pd.DataFrame, max_len: int) -> np.ndarray:
    n_edges = len(df)
    matrix = np.full((n_edges, max_len), np.nan, dtype=np.float32)
    bytes_col = df["LocalTE_bytes"].to_numpy()
    dtype_col = df["LocalTE_dtype"].to_numpy()
    len_col = df["LocalTE_len"].fillna(0).astype(int).to_numpy()
    codec_col = df["LocalTE_codec"].to_numpy() if "LocalTE_codec" in df.columns else None
    for idx, blob in enumerate(bytes_col):
        length = len_col[idx]
        if length <= 0:
            continue
        if codec_col is not None and str(codec_col[idx]).lower() == 'zlib':
            import zlib as _z
            try:
                raw = _z.decompress(blob)
            except Exception:
                raw = blob
            series = np.frombuffer(raw, dtype=dtype_col[idx], count=length).astype(np.float32, copy=False)
        else:
            series = np.frombuffer(blob, dtype=dtype_col[idx], count=length).astype(np.float32, copy=False)
        matrix[idx, :length] = series
    return matrix


def process_timepoint(
    t_idx: int,
    matrix: np.ndarray,
    sources: Sequence[object],
    targets: Sequence[object],
    labels: Sequence[str],
    te_cutoff: float,
    fdr: float,
    orig_indices: Optional[Sequence[int]] = None,
    offset: int = 0,
) -> Optional[pd.DataFrame]:
    global_idx = offset + t_idx
    if global_idx >= len(labels):
        return None
    # Upcast to float32 for numeric stability while keeping matrix memory small (float16)
    column = matrix[:, t_idx].astype(np.float32, copy=False)
    valid_idx = np.where(~np.isnan(column))[0]
    if valid_idx.size == 0:
        return None
    vals = column[valid_idx]
    gt_mask = vals > te_cutoff
    if not np.any(gt_mask):
        return None
    sel_idx = valid_idx[gt_mask]
    vals = vals[gt_mask]
    if vals.size < 2:
        return None
    std = vals.std(ddof=0)
    if np.allclose(std, 0.0):
        return None

    mean = vals.mean()
    zscores = (vals - mean) / std
    pvals = 1 - norm.cdf(zscores)
    # Fast BH-FDR (vectorized)
    m = pvals.size
    order = np.argsort(pvals)
    ranks = np.arange(1, m + 1, dtype=np.float32)
    q_tmp = pvals[order] * (m / ranks)
    q_sorted = np.minimum.accumulate(q_tmp[::-1])[::-1]
    qvals = np.empty_like(pvals)
    qvals[order] = np.clip(q_sorted, 0.0, 1.0)
    keep = qvals < fdr
    if not np.any(keep):
        return None

    data = {
        "Source": np.asarray(sources, dtype=object)[sel_idx][keep],
        "Target": np.asarray(targets, dtype=object)[sel_idx][keep],
        "TE": vals[keep],
        "zscore": zscores[keep],
        "p_value": pvals[keep],
        "q_value": qvals[keep],
        "TimeIndex": np.full(np.count_nonzero(keep), global_idx, dtype=int),
        "TimeLabel": np.full(np.count_nonzero(keep), labels[global_idx], dtype=object),
    }
    if orig_indices is not None:
        if global_idx >= len(orig_indices):
            return None
        oi = int(orig_indices[global_idx])
        data["OriginalIndex"] = np.full(np.count_nonzero(keep), oi, dtype=int)
        # If labels correspond to original indices, TimeLabel already matches; else, optionally add OriginalLabel
    return pd.DataFrame(data)


# Globals for multiprocessing workers
_MATRIX = None
_SOURCES = None
_TARGETS = None
_LABELS = None
_ORIG_IDX_SEQ = None
_CUTOFF = 0.0
_ALPHA = 0.05
_TRIPLE_DTYPE = np.dtype([("TE", np.float32), ("Source", np.int32), ("Target", np.int32)])
_COMPACT_TRIPLE_DTYPE = np.dtype([("TE", np.float16), ("Row", np.int32)])


def _worker_init(matrix, sources, targets, labels, orig_idx_seq, cutoff, alpha):
    global _MATRIX, _SOURCES, _TARGETS, _LABELS, _ORIG_IDX_SEQ, _CUTOFF, _ALPHA
    _MATRIX = matrix
    _SOURCES = sources
    _TARGETS = targets
    _LABELS = labels
    _ORIG_IDX_SEQ = orig_idx_seq
    _CUTOFF = cutoff
    _ALPHA = alpha


def _worker_task(t_idx: int):
    return process_timepoint(t_idx, _MATRIX, _SOURCES, _TARGETS, _LABELS, _CUTOFF, _ALPHA, _ORIG_IDX_SEQ)


def _load_gene_names_list(path: str) -> List[str]:
    names: List[str] = []
    resolved = coerce_input_path(path)
    with open(resolved, encoding="utf-8") as fh:
        for line in fh:
            names.append(line.strip())
    return names


def _idx_to_name(idx: int, names: Sequence[str]) -> str:
    try:
        i = int(idx)
        if 1 <= i <= len(names):
            return names[i - 1]
    except Exception:
        pass
    return str(idx)


def _chunk_batch_rows() -> int:
    try:
        batch = int(os.environ.get("LOCAL_TE_READ_BATCH_ROWS", "8192"))
        return max(256, batch)
    except Exception:
        return 8192


def _chunk_thread_workers() -> int:
    try:
        workers = int(os.environ.get("LOCAL_TE_TIME_THREADS", "1"))
        return max(1, workers)
    except Exception:
        return 1


def _use_single_pass_mode() -> bool:
    flag = os.environ.get("LOCAL_TE_SINGLE_PASS", "").strip().lower()
    return flag in {"1", "true", "yes", "on", "fast"}


def _arrow_use_threads() -> bool:
    flag = os.environ.get("LOCAL_TE_ARROW_THREADS", "off").strip().lower()
    return flag in {"1", "true", "yes", "on"}


def _use_histogram_mode() -> bool:
    flag = os.environ.get("LOCAL_TE_HISTOGRAM_FDR", "on").strip().lower()
    return flag in {"1", "true", "yes", "on", "fast"}


def _use_vectorized_edge_writer() -> bool:
    flag = os.environ.get("LOCAL_TE_EDGE_VECTORIZED", "on").strip().lower()
    return flag in {"1", "true", "yes", "on", "fast"}


def _histogram_accumulator_mode() -> str:
    mode = os.environ.get("LOCAL_TE_HIST_ACCUM", "auto").strip().lower()
    if mode in {"legacy", "bincount", "combined"}:
        return "legacy"
    if mode in {"column", "columns"}:
        return "column"
    if mode in {"numba", "jit"}:
        return "numba" if _update_hist_numba is not None else "indexed"
    if mode in {"auto", ""}:
        return "numba" if _update_hist_numba is not None else "indexed"
    return "indexed"


def _use_numeric_edge_output() -> bool:
    return _env_flag("LOCAL_TE_EDGE_NUMERIC", "off")


def _use_dataset_output() -> bool:
    return _env_flag("LOCAL_TE_OUTPUT_DATASET", "off")


def _dataset_output_dir(output_file_path: Path) -> Path:
    configured = os.environ.get("LOCAL_TE_OUTPUT_DATASET_DIR", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    if output_file_path.suffix:
        return output_file_path.with_suffix("").with_name(output_file_path.stem + ".dataset")
    return output_file_path


def _write_dataset_manifest(
    dataset_dir: Path,
    *,
    output_mode: str,
    part_pattern: str,
    part_count: int,
    row_count: int,
    fdr: float,
    te_cutoff: float,
    chunk_size: int,
    max_timepoints: int,
) -> None:
    manifest = {
        "format": "TENETPLUS_LOCALTE_GRN_DATASET",
        "version": 1,
        "output_mode": output_mode,
        "part_pattern": part_pattern,
        "part_count": int(part_count),
        "row_count": int(row_count),
        "fdr": float(fdr),
        "te_cutoff": float(te_cutoff),
        "chunk_size": int(chunk_size),
        "max_timepoints": int(max_timepoints),
    }
    (dataset_dir / "_metadata.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _warm_hist_numba() -> None:
    if _histogram_accumulator_mode() != "numba" or _update_hist_numba is None:
        return
    hist = np.zeros((1, 4), dtype=np.int64)
    raw = np.zeros((1, 1), dtype=np.uint16)
    vals = np.zeros((1, 1), dtype=np.float32)
    _update_hist_numba(hist, raw, vals, 0.0)


def _use_compact_temp() -> bool:
    return _env_flag("LOCAL_TE_COMPACT_TEMP", "off")


def _resolve_temp_parent(output_file_path: Optional[Path]) -> Path:
    configured = os.environ.get("LOCAL_TE_TEMP_DIR", "").strip()
    if configured:
        temp_parent = Path(configured).expanduser()
    elif output_file_path is not None:
        temp_parent = output_file_path.parent
    else:
        temp_parent = Path(tempfile.gettempdir())
    temp_parent.mkdir(parents=True, exist_ok=True)
    return temp_parent.resolve()


def _estimate_chunk_temp_bytes(chunk_meta: Dict[str, object], chunk_size: int, max_len: int) -> int:
    """Worst-case worker temp bytes for the single-pass per-timepoint triples."""
    try:
        total_edges = int(chunk_meta.get("total_edges", 0))
    except Exception:
        total_edges = 0
    if total_edges <= 0:
        return 0
    time_len = max(1, min(int(chunk_size), int(max_len)))
    itemsize = _COMPACT_TRIPLE_DTYPE.itemsize if _use_compact_temp() else _TRIPLE_DTYPE.itemsize
    return int(total_edges * time_len * itemsize)


def _reserve_temp_bytes() -> int:
    try:
        gib = float(os.environ.get("LOCAL_TE_TEMP_RESERVE_GIB", "100"))
    except Exception:
        gib = 100.0
    return max(0, int(gib * (1024 ** 3)))


def _cap_workers_by_temp(requested: int, temp_parent: Path, est_per_worker: int) -> int:
    requested = max(1, int(requested))
    if est_per_worker <= 0 or _env_flag("LOCAL_TE_DISABLE_TEMP_CAP", "off"):
        return requested
    try:
        free = shutil.disk_usage(temp_parent).free
    except Exception:
        return requested
    usable = max(0, free - _reserve_temp_bytes())
    # Keep a safety margin because parquet decoding and output writes also use space.
    safe_per_worker = max(est_per_worker, int(est_per_worker * 1.15))
    capped = max(1, min(requested, usable // safe_per_worker if safe_per_worker else requested))
    if capped < requested:
        print(
            "[LocalTE GRN] Reducing workers to avoid temp-disk exhaustion: "
            f"requested={requested}, capped={capped}, "
            f"temp_dir={temp_parent}, free={_human_bytes(free)}, "
            f"estimated_temp_per_worker={_human_bytes(est_per_worker)}, "
            f"reserve={_human_bytes(_reserve_temp_bytes())}. "
            "Set LOCAL_TE_TEMP_DIR to a larger filesystem or LOCAL_TE_DISABLE_TEMP_CAP=on to override.",
            flush=True,
        )
    return int(capped)


def _infer_chunk_length(pf: pq.ParquetFile) -> int:
    if pf.num_row_groups == 0:
        return 0
    first = pf.read_row_group(0, columns=["Values"], use_threads=False)
    column = first.column("Values")
    if isinstance(column, pa.ChunkedArray):
        if column.num_chunks == 0:
            return 0
        arr = column.chunk(0)
    else:
        arr = column
    if not hasattr(arr, "offsets"):
        return 0
    offsets = arr.offsets.to_numpy()
    if offsets.size < 2:
        return 0
    return int(offsets[1] - offsets[0])


def _compute_chunk_stats(
    pf: pq.ParquetFile,
    raw_chunk_len: int,
    chunk_len: int,
    te_cutoff: float,
    batch_rows: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts = np.zeros(chunk_len, dtype=np.int64)
    sum_vals = np.zeros(chunk_len, dtype=np.float64)
    sum_sq = np.zeros(chunk_len, dtype=np.float64)
    use_arrow_threads = _arrow_use_threads()
    for batch in pf.iter_batches(columns=["Values"], batch_size=batch_rows, use_threads=use_arrow_threads):
        values = batch.column(0)
        flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
        if flat.size == 0:
            continue
        mat = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)[:, :chunk_len]
        if mat.size == 0:
            continue
        mask = mat > te_cutoff
        if not np.any(mask):
            continue
        counts += mask.sum(axis=0).astype(np.int64)
        selected = mat * mask
        sum_vals += selected.sum(axis=0, dtype=np.float64)
        sum_sq += np.square(selected, dtype=np.float64).sum(axis=0, dtype=np.float64)
    return counts, sum_vals, sum_sq


def _fill_chunk_memmaps(
    pf: pq.ParquetFile,
    raw_chunk_len: int,
    chunk_len: int,
    te_cutoff: float,
    batch_rows: int,
    counts: np.ndarray,
    tmp_dir: Path,
) -> Tuple[np.ndarray, Optional[np.memmap], Optional[np.memmap], Optional[np.memmap]]:
    offsets = np.zeros(chunk_len + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    total_entries = int(offsets[-1])
    if total_entries <= 0:
        return offsets, None, None, None

    te_mem = np.memmap(tmp_dir / "te_values.bin", dtype=np.float32, mode="w+", shape=(total_entries,))
    src_mem = np.memmap(tmp_dir / "sources.bin", dtype=np.int32, mode="w+", shape=(total_entries,))
    tgt_mem = np.memmap(tmp_dir / "targets.bin", dtype=np.int32, mode="w+", shape=(total_entries,))
    write_pos = offsets[:-1].copy()

    use_arrow_threads = _arrow_use_threads()
    for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=batch_rows, use_threads=use_arrow_threads):
        values = batch.column("Values")
        flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
        if flat.size == 0:
            continue
        mat = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)[:, :chunk_len]
        if mat.size == 0:
            continue
        mask = mat > te_cutoff
        if not np.any(mask):
            continue
        src_arr = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
        tgt_arr = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
        for local_idx in range(chunk_len):
            col_mask = mask[:, local_idx]
            if not np.any(col_mask):
                continue
            idxs = np.nonzero(col_mask)[0]
            n = idxs.size
            start = write_pos[local_idx]
            end = start + n
            te_mem[start:end] = mat[idxs, local_idx]
            src_mem[start:end] = src_arr[idxs]
            tgt_mem[start:end] = tgt_arr[idxs]
            write_pos[local_idx] = end

    return offsets, te_mem, src_mem, tgt_mem


def _build_timepoint_df(
    te_vals: np.ndarray,
    src_idx: np.ndarray,
    tgt_idx: np.ndarray,
    mean: float,
    std: float,
    global_idx: int,
    labels: Sequence[str],
    orig_idx_seq: Optional[Sequence[int]],
    gene_names_arr: np.ndarray,
    fdr: float,
) -> Optional[pd.DataFrame]:
    if te_vals.size < 2 or not np.isfinite(std) or std <= 0.0:
        return None

    te_arr = te_vals.astype(np.float32, copy=False)
    zscores = (te_arr - mean) / std
    pvals = 1 - norm.cdf(zscores)
    m = pvals.size
    if m == 0:
        return None
    order = np.argsort(pvals)
    ranks = np.arange(1, m + 1, dtype=np.float32)
    q_tmp = pvals[order] * (m / ranks)
    q_sorted = np.minimum.accumulate(q_tmp[::-1])[::-1]
    qvals = np.empty_like(pvals)
    qvals[order] = np.clip(q_sorted, 0.0, 1.0)
    keep = qvals < fdr
    if not np.any(keep):
        return None

    idx_keep = np.nonzero(keep)[0]
    te_keep = te_arr[idx_keep]
    z_keep = zscores[idx_keep]
    p_keep = pvals[idx_keep]
    q_keep = qvals[idx_keep]
    sources = src_idx[idx_keep].astype(np.int64, copy=False)
    targets = tgt_idx[idx_keep].astype(np.int64, copy=False)
    # Gene indices are 1-based; clamp to valid range defensively
    src_names = gene_names_arr[np.clip(sources - 1, 0, gene_names_arr.size - 1)]
    tgt_names = gene_names_arr[np.clip(targets - 1, 0, gene_names_arr.size - 1)]
    time_label = labels[global_idx] if global_idx < len(labels) else str(global_idx)

    data = {
        "Source": src_names,
        "Target": tgt_names,
        "TE": te_keep,
        "zscore": z_keep,
        "p_value": p_keep,
        "q_value": q_keep,
        "TimeIndex": np.full(te_keep.size, global_idx, dtype=int),
        "TimeLabel": np.full(te_keep.size, time_label, dtype=object),
    }
    if orig_idx_seq is not None and global_idx < len(orig_idx_seq):
        data["OriginalIndex"] = np.full(te_keep.size, int(orig_idx_seq[global_idx]), dtype=int)
    return pd.DataFrame(data)


def _filter_timepoint_arrays(
    te_vals: np.ndarray,
    src_idx: np.ndarray,
    tgt_idx: np.ndarray,
    mean: float,
    std: float,
    fdr: float,
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    if te_vals.size < 2 or not np.isfinite(std) or std <= 0.0:
        return None
    te_arr = te_vals.astype(np.float32, copy=False)
    zscores = (te_arr - mean) / std
    pvals = 1 - norm.cdf(zscores)
    m = pvals.size
    if m == 0:
        return None
    order = np.argsort(pvals)
    ranks = np.arange(1, m + 1, dtype=np.float32)
    q_tmp = pvals[order] * (m / ranks)
    q_sorted = np.minimum.accumulate(q_tmp[::-1])[::-1]
    qvals = np.empty_like(pvals)
    qvals[order] = np.clip(q_sorted, 0.0, 1.0)
    keep = qvals < fdr
    if not np.any(keep):
        return None
    idx_keep = np.nonzero(keep)[0]
    return (
        src_idx[idx_keep].astype(np.int32, copy=False),
        tgt_idx[idx_keep].astype(np.int32, copy=False),
        te_arr[idx_keep].astype(np.float32, copy=False),
    )


def _numeric_trim_outdegree(
    src_vals: np.ndarray,
    tgt_vals: np.ndarray,
    te_vals: np.ndarray,
    threshold: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    source_order: List[int] = []
    targets_by_source: Dict[int, List[int]] = {}
    tes_by_source: Dict[int, List[float]] = {}
    for source, target, te in zip(src_vals.astype(np.int64), tgt_vals.astype(np.int64), te_vals.astype(np.float32)):
        s = int(source)
        if s not in targets_by_source:
            source_order.append(s)
            targets_by_source[s] = []
            tes_by_source[s] = []
        targets_by_source[s].append(int(target))
        tes_by_source[s].append(float(te))
    if _trim_indirect_edges is None:
        raise RuntimeError("trim_indirect helper is unavailable; cannot trim indirect edges.")
    indirect = _trim_indirect_edges(source_order, targets_by_source, tes_by_source, threshold)
    out_sources: List[int] = []
    out_counts: List[int] = []
    kept = 0
    removed = 0
    for source in source_order:
        keep_count = 0
        for is_indirect in indirect[source]:
            if is_indirect:
                removed += 1
            else:
                kept += 1
                keep_count += 1
        if keep_count:
            out_sources.append(source)
            out_counts.append(keep_count)
    return (
        np.asarray(out_sources, dtype=np.int32),
        np.asarray(out_counts, dtype=np.int64),
        {
            "sources": len(source_order),
            "edges": int(src_vals.size),
            "kept": int(kept),
            "removed": int(removed),
        },
    )


def _process_chunk_single_pass(
    pf: pq.ParquetFile,
    chunk_idx: int,
    raw_chunk_len: int,
    chunk_len: int,
    chunk_size: int,
    te_cutoff: float,
    fdr: float,
    labels: Sequence[str],
    orig_idx_seq: Optional[Sequence[int]],
    gene_names_arr: np.ndarray,
    max_len: int,
    batch_rows: int,
    temp_root: Path,
) -> List[pd.DataFrame]:
    chunk_tmp = temp_root / f"chunk_{chunk_idx:04d}"
    chunk_tmp.mkdir(parents=True, exist_ok=True)
    counts = np.zeros(chunk_len, dtype=np.int64)
    sum_vals = np.zeros(chunk_len, dtype=np.float64)
    sum_sq = np.zeros(chunk_len, dtype=np.float64)
    handles: List[Optional[object]] = []
    compact_temp = _use_compact_temp()
    total_rows = int(pf.metadata.num_rows) if pf.metadata is not None else 0
    source_lookup = np.empty(total_rows, dtype=np.int32) if compact_temp and total_rows > 0 else None
    target_lookup = np.empty(total_rows, dtype=np.int32) if compact_temp and total_rows > 0 else None
    row_offset = 0
    for local_idx in range(chunk_len):
        handles.append((chunk_tmp / f"time_{local_idx:03d}.bin").open("wb"))

    try:
        use_arrow_threads = _arrow_use_threads()
        for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=batch_rows, use_threads=use_arrow_threads):
            values = batch.column("Values")
            flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
            if flat.size == 0:
                row_offset += batch.num_rows
                continue
            mat16_full = np.asarray(flat, dtype=np.float16).reshape(-1, raw_chunk_len)
            mat16 = mat16_full[:, :chunk_len]
            if mat16.size == 0:
                row_offset += batch.num_rows
                continue
            mat = mat16.astype(np.float32, copy=False)
            mask = mat > te_cutoff
            if not np.any(mask):
                if source_lookup is not None and target_lookup is not None:
                    source_lookup[row_offset:row_offset + batch.num_rows] = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
                    target_lookup[row_offset:row_offset + batch.num_rows] = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
                row_offset += batch.num_rows
                continue
            src_arr = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
            tgt_arr = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
            if source_lookup is not None and target_lookup is not None:
                source_lookup[row_offset:row_offset + batch.num_rows] = src_arr
                target_lookup[row_offset:row_offset + batch.num_rows] = tgt_arr
            for local_idx in range(chunk_len):
                col_mask = mask[:, local_idx]
                if not np.any(col_mask):
                    continue
                vals = mat[col_mask, local_idx]
                n = vals.size
                counts[local_idx] += n
                sum_vals[local_idx] += float(np.sum(vals, dtype=np.float64))
                sum_sq[local_idx] += float(np.sum(vals * vals, dtype=np.float64))
                if compact_temp and source_lookup is not None and target_lookup is not None:
                    arr = np.empty(n, dtype=_COMPACT_TRIPLE_DTYPE)
                    arr["TE"] = mat16[col_mask, local_idx]
                    arr["Row"] = np.nonzero(col_mask)[0].astype(np.int32, copy=False) + np.int32(row_offset)
                else:
                    arr = np.empty(n, dtype=_TRIPLE_DTYPE)
                    arr["TE"] = vals.astype(np.float32, copy=False)
                    arr["Source"] = src_arr[col_mask]
                    arr["Target"] = tgt_arr[col_mask]
                handles[local_idx].write(arr.tobytes())
            row_offset += batch.num_rows
    finally:
        for handle in handles:
            try:
                handle.close()
            except Exception:
                pass

    thread_workers = _chunk_thread_workers()
    global_offset = chunk_idx * chunk_size

    def _load_timepoint(local_idx: int) -> Optional[Tuple[int, pd.DataFrame]]:
        global_idx = global_offset + local_idx
        if global_idx >= len(labels) or global_idx >= max_len:
            return None
        cnt = int(counts[local_idx])
        if cnt < 2:
            return None
        file_path = chunk_tmp / f"time_{local_idx:03d}.bin"
        if not file_path.exists() or file_path.stat().st_size == 0:
            return None
        total = float(sum_vals[local_idx])
        mean = total / cnt
        variance = (sum_sq[local_idx] / cnt) - (mean ** 2)
        if variance <= 1e-12:
            return None
        std = math.sqrt(max(variance, 0.0))
        dtype = _COMPACT_TRIPLE_DTYPE if compact_temp and source_lookup is not None and target_lookup is not None else _TRIPLE_DTYPE
        triples = np.memmap(file_path, dtype=dtype, mode="r")
        try:
            te_vals = np.asarray(triples["TE"], dtype=np.float32)
            if dtype == _COMPACT_TRIPLE_DTYPE:
                row_vals = np.asarray(triples["Row"], dtype=np.int64)
                src_vals = np.asarray(source_lookup[row_vals], dtype=np.int32)
                tgt_vals = np.asarray(target_lookup[row_vals], dtype=np.int32)
            else:
                src_vals = np.asarray(triples["Source"], dtype=np.int32)
                tgt_vals = np.asarray(triples["Target"], dtype=np.int32)
            df_time = _build_timepoint_df(
                te_vals,
                src_vals,
                tgt_vals,
                mean,
                std,
                global_idx,
                labels,
                orig_idx_seq,
                gene_names_arr,
                fdr,
            )
            if df_time is None or df_time.empty:
                return None
            return global_idx, df_time
        finally:
            del triples

    results: List[Tuple[int, pd.DataFrame]] = []
    if thread_workers > 1 and chunk_len > 1:
        with ThreadPoolExecutor(max_workers=min(thread_workers, chunk_len)) as executor:
            for res in executor.map(_load_timepoint, range(chunk_len)):
                if res is not None:
                    results.append(res)
    else:
        for local_idx in range(chunk_len):
            res = _load_timepoint(local_idx)
            if res is not None:
                results.append(res)

    for local_idx in range(chunk_len):
        (chunk_tmp / f"time_{local_idx:03d}.bin").unlink(missing_ok=True)
    shutil.rmtree(chunk_tmp, ignore_errors=True)

    if not results:
        return []
    results.sort(key=lambda x: x[0])
    return [df for _, df in results]


def _histogram_thresholds_for_chunk(
    pf: pq.ParquetFile,
    raw_chunk_len: int,
    chunk_len: int,
    te_cutoff: float,
    fdr: float,
    batch_rows: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute exact BH-FDR TE thresholds from per-timepoint float16 histograms."""
    hist_matrix = np.zeros((chunk_len, 65536), dtype=np.int64)
    hist_mode = _histogram_accumulator_mode()
    hist_flat = hist_matrix.reshape(-1)
    offsets = (np.arange(chunk_len, dtype=np.uint32) << np.uint32(16)).reshape(1, -1)
    use_arrow_threads = _arrow_use_threads()
    _accumulate_histogram_from_parquet(
        pf,
        hist_matrix=hist_matrix,
        hist_flat=hist_flat,
        offsets=offsets,
        raw_chunk_len=raw_chunk_len,
        chunk_len=chunk_len,
        te_cutoff=te_cutoff,
        batch_rows=batch_rows,
        hist_mode=hist_mode,
        use_arrow_threads=use_arrow_threads,
    )
    return _thresholds_from_histogram(hist_matrix, fdr)


def _accumulate_histogram_from_parquet(
    pf: pq.ParquetFile,
    *,
    hist_matrix: np.ndarray,
    hist_flat: np.ndarray,
    offsets: np.ndarray,
    raw_chunk_len: int,
    chunk_len: int,
    te_cutoff: float,
    batch_rows: int,
    hist_mode: str,
    use_arrow_threads: bool,
) -> None:
    for batch in pf.iter_batches(columns=["Values"], batch_size=batch_rows, use_threads=use_arrow_threads):
        values = batch.column("Values")
        flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
        if flat.size == 0:
            continue
        mat = np.asarray(flat, dtype=np.float16).reshape(-1, raw_chunk_len)[:, :chunk_len]
        if mat.size == 0:
            continue
        if hist_mode == "numba" and _update_hist_numba is not None:
            raw_u16 = np.ascontiguousarray(mat).view(np.uint16)
            _update_hist_numba(hist_matrix, raw_u16, mat.astype(np.float32, copy=False), float(te_cutoff))
            continue
        mask = mat > te_cutoff
        if not np.any(mask):
            continue
        raw_u16 = np.ascontiguousarray(mat).view(np.uint16)
        if hist_mode == "legacy":
            raw_mat = raw_u16.astype(np.uint32, copy=False)
            combined = raw_mat + offsets
            selected = combined[mask]
            counts = np.bincount(selected.astype(np.int64, copy=False), minlength=chunk_len * 65536)
            hist_matrix += counts.reshape(chunk_len, 65536)
        elif hist_mode == "column":
            for local_idx in range(chunk_len):
                col_mask = mask[:, local_idx]
                if np.any(col_mask):
                    hist_matrix[local_idx] += np.bincount(
                        raw_u16[:, local_idx][col_mask].astype(np.int64, copy=False),
                        minlength=65536,
                    )
        else:
            raw_mat = raw_u16.astype(np.uint32, copy=False)
            combined = raw_mat + offsets
            selected = combined[mask].astype(np.int64, copy=False)
            np.add.at(hist_flat, selected, 1)


def _thresholds_from_histogram(hist_matrix: np.ndarray, fdr: float) -> Tuple[np.ndarray, np.ndarray]:
    chunk_len = int(hist_matrix.shape[0])

    thresholds = np.full(chunk_len, np.nan, dtype=np.float32)
    selected_counts = np.zeros(chunk_len, dtype=np.int64)
    for local_idx in range(chunk_len):
        hist = hist_matrix[local_idx]
        nz = np.nonzero(hist)[0].astype(np.uint16, copy=False)
        if nz.size == 0:
            continue
        raw_vals = nz
        counts = hist[nz.astype(np.int64)].astype(np.int64, copy=False)
        vals = raw_vals.view(np.float16).astype(np.float32)
        m = int(counts.sum())
        if m < 2:
            continue
        total = np.sum(vals * counts, dtype=np.float64)
        sum_sq = np.sum((vals.astype(np.float64) ** 2) * counts, dtype=np.float64)
        mean = total / m
        variance = (sum_sq / m) - (mean ** 2)
        if variance <= 1e-12:
            continue
        std = math.sqrt(max(variance, 0.0))
        order = np.argsort(vals)[::-1]
        vals_sorted = vals[order]
        counts_sorted = counts[order]
        ranks_end = np.cumsum(counts_sorted, dtype=np.int64).astype(np.float64)
        pvals = 1.0 - norm.cdf((vals_sorted - mean) / std)
        q_tmp = pvals * (m / ranks_end)
        q_group = np.minimum.accumulate(q_tmp[::-1])[::-1]
        keep = q_group < fdr
        if not np.any(keep):
            continue
        threshold = float(np.min(vals_sorted[keep]))
        thresholds[local_idx] = np.float32(threshold)
        selected_counts[local_idx] = int(np.sum(counts_sorted[keep]))
    return thresholds, selected_counts


def _histogram_thresholds_for_chunk_files(
    pfs: Sequence[pq.ParquetFile],
    raw_chunk_len: int,
    chunk_len: int,
    te_cutoff: float,
    fdr: float,
    batch_rows: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute exact BH-FDR TE thresholds across worker-part parquet files."""
    hist_matrix = np.zeros((chunk_len, 65536), dtype=np.int64)
    hist_mode = _histogram_accumulator_mode()
    hist_flat = hist_matrix.reshape(-1)
    offsets = (np.arange(chunk_len, dtype=np.uint32) << np.uint32(16)).reshape(1, -1)
    use_arrow_threads = _arrow_use_threads()
    for pf in pfs:
        _accumulate_histogram_from_parquet(
            pf,
            hist_matrix=hist_matrix,
            hist_flat=hist_flat,
            offsets=offsets,
            raw_chunk_len=raw_chunk_len,
            chunk_len=chunk_len,
            te_cutoff=te_cutoff,
            batch_rows=batch_rows,
            hist_mode=hist_mode,
            use_arrow_threads=use_arrow_threads,
        )
    return _thresholds_from_histogram(hist_matrix, fdr)


def _process_chunk_histogram_file_task(args):
    (
        chunk_path,
        chunk_idx,
        chunk_size,
        gene_names,
        te_cutoff,
        fdr,
        labels,
        output_dir,
        max_len,
        batch_rows,
    ) = args
    chunk_path = Path(chunk_path)
    pf = pq.ParquetFile(chunk_path)
    raw_chunk_len = _infer_chunk_length(pf)
    if raw_chunk_len <= 0:
        return None
    global_offset = int(chunk_idx) * int(chunk_size)
    if global_offset >= int(max_len):
        return None
    remaining = int(max_len) - global_offset
    chunk_len = min(raw_chunk_len, remaining)
    thresholds, selected_counts = _histogram_thresholds_for_chunk(
        pf,
        raw_chunk_len=raw_chunk_len,
        chunk_len=chunk_len,
        te_cutoff=float(te_cutoff),
        fdr=float(fdr),
        batch_rows=int(batch_rows),
    )
    if not np.any(np.isfinite(thresholds)):
        return None

    gene_names_arr = np.asarray(gene_names, dtype=object)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"grn_chunk_{int(chunk_idx):04d}.parquet"
    numeric_edge_output = _use_numeric_edge_output()
    if numeric_edge_output:
        schema = pa.schema(
            [
                ("Source", pa.int32()),
                ("TE", pa.float32()),
                ("Target", pa.int32()),
                ("Time", pa.int32()),
            ]
        )
    else:
        schema = pa.schema(
            [
                ("Source", pa.string()),
                ("TE", pa.float32()),
                ("Target", pa.string()),
                ("Time", pa.string()),
            ]
        )
    if _use_vectorized_edge_writer() or numeric_edge_output:
        writer: Optional[pq.ParquetWriter] = None
        written_rows = 0
        use_arrow_threads = _arrow_use_threads()
        write_rows = max(50_000, int(os.environ.get("LOCAL_TE_GRN_WRITE_ROWS", "500000")))
        finite_thresholds = np.isfinite(thresholds)

        def write_selected(
            src_idx: np.ndarray,
            tgt_idx: np.ndarray,
            time_idx: np.ndarray,
            te_vals: np.ndarray,
        ) -> None:
            nonlocal writer, written_rows
            if te_vals.size == 0:
                return
            global_time_idx = global_offset + time_idx.astype(np.int64, copy=False)
            if numeric_edge_output:
                table = pa.table(
                    {
                        "Source": pa.array(src_idx.astype(np.int32, copy=False), type=pa.int32()),
                        "TE": pa.array(te_vals.astype(np.float32, copy=False), type=pa.float32()),
                        "Target": pa.array(tgt_idx.astype(np.int32, copy=False), type=pa.int32()),
                        "Time": pa.array(global_time_idx.astype(np.int32, copy=False), type=pa.int32()),
                    },
                    schema=schema,
                )
            else:
                src_names = gene_names_arr[np.clip(src_idx.astype(np.int64, copy=False) - 1, 0, gene_names_arr.size - 1)]
                tgt_names = gene_names_arr[np.clip(tgt_idx.astype(np.int64, copy=False) - 1, 0, gene_names_arr.size - 1)]
                time_labels = np.asarray(
                    [str(labels[int(i)]) if int(i) < len(labels) else str(int(i)) for i in global_time_idx],
                    dtype=object,
                )
                table = pa.table(
                    {
                        "Source": pa.array(src_names.astype(object), type=pa.string()),
                        "TE": pa.array(te_vals.astype(np.float32, copy=False), type=pa.float32()),
                        "Target": pa.array(tgt_names.astype(object), type=pa.string()),
                        "Time": pa.array(time_labels, type=pa.string()),
                    },
                    schema=schema,
                )
            if writer is None:
                writer = pq.ParquetWriter(out_path, schema, compression="snappy")
            writer.write_table(table)
            written_rows += int(te_vals.size)

        try:
            src_parts: List[np.ndarray] = []
            tgt_parts: List[np.ndarray] = []
            time_parts: List[np.ndarray] = []
            te_parts: List[np.ndarray] = []
            buffered_rows = 0

            def flush_parts() -> None:
                nonlocal src_parts, tgt_parts, time_parts, te_parts, buffered_rows
                if buffered_rows <= 0:
                    return
                write_selected(
                    np.concatenate(src_parts),
                    np.concatenate(tgt_parts),
                    np.concatenate(time_parts),
                    np.concatenate(te_parts),
                )
                src_parts = []
                tgt_parts = []
                time_parts = []
                te_parts = []
                buffered_rows = 0

            for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=int(batch_rows), use_threads=use_arrow_threads):
                values = batch.column("Values")
                flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
                if flat.size == 0:
                    continue
                mat = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)[:, :chunk_len]
                if mat.size == 0:
                    continue
                # Broadcast thresholds across rows. NaN thresholds are masked out.
                mask = (mat >= thresholds.reshape(1, -1)) & finite_thresholds.reshape(1, -1)
                if not np.any(mask):
                    continue
                row_idx, local_idx = np.nonzero(mask)
                if row_idx.size == 0:
                    continue
                src_arr = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
                tgt_arr = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
                src_parts.append(src_arr[row_idx])
                tgt_parts.append(tgt_arr[row_idx])
                time_parts.append(local_idx.astype(np.int32, copy=False))
                te_parts.append(mat[row_idx, local_idx].astype(np.float32, copy=False))
                buffered_rows += int(row_idx.size)
                if buffered_rows >= write_rows:
                    flush_parts()
            flush_parts()
        finally:
            if writer is not None:
                writer.close()

        if written_rows <= 0:
            try:
                out_path.unlink(missing_ok=True)
            except Exception:
                pass
            return None
        return str(out_path)

    writer: Optional[pq.ParquetWriter] = None
    use_arrow_threads = _arrow_use_threads()
    flush_rows = max(50_000, int(os.environ.get("LOCAL_TE_GRN_WRITE_ROWS", "500000")))
    src_buf: List[np.ndarray] = []
    tgt_buf: List[np.ndarray] = []
    te_buf: List[np.ndarray] = []
    time_buf: List[np.ndarray] = []
    buffered_rows = 0
    written_rows = 0

    def flush() -> None:
        nonlocal writer, buffered_rows, written_rows, src_buf, tgt_buf, te_buf, time_buf
        if buffered_rows <= 0:
            return
        src_idx = np.concatenate(src_buf).astype(np.int64, copy=False)
        tgt_idx = np.concatenate(tgt_buf).astype(np.int64, copy=False)
        te_vals = np.concatenate(te_buf).astype(np.float32, copy=False)
        times = np.concatenate(time_buf).astype(object, copy=False)
        src_names = gene_names_arr[np.clip(src_idx - 1, 0, gene_names_arr.size - 1)]
        tgt_names = gene_names_arr[np.clip(tgt_idx - 1, 0, gene_names_arr.size - 1)]
        table = pa.table(
            {
                "Source": pa.array(src_names.astype(object), type=pa.string()),
                "TE": pa.array(te_vals, type=pa.float32()),
                "Target": pa.array(tgt_names.astype(object), type=pa.string()),
                "Time": pa.array(times, type=pa.string()),
            },
            schema=schema,
        )
        if writer is None:
            writer = pq.ParquetWriter(out_path, schema, compression="snappy")
        writer.write_table(table)
        written_rows += buffered_rows
        src_buf = []
        tgt_buf = []
        te_buf = []
        time_buf = []
        buffered_rows = 0

    try:
        for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=int(batch_rows), use_threads=use_arrow_threads):
            values = batch.column("Values")
            flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
            if flat.size == 0:
                continue
            mat16 = np.asarray(flat, dtype=np.float16).reshape(-1, raw_chunk_len)[:, :chunk_len]
            if mat16.size == 0:
                continue
            src_arr = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
            tgt_arr = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
            for local_idx in range(chunk_len):
                threshold = thresholds[local_idx]
                if not np.isfinite(threshold):
                    continue
                col = mat16[:, local_idx].astype(np.float32, copy=False)
                mask = col >= threshold
                if not np.any(mask):
                    continue
                n = int(np.count_nonzero(mask))
                global_idx = global_offset + local_idx
                time_label = labels[global_idx] if global_idx < len(labels) else str(global_idx)
                src_buf.append(src_arr[mask])
                tgt_buf.append(tgt_arr[mask])
                te_buf.append(col[mask].astype(np.float32, copy=False))
                time_buf.append(np.full(n, str(time_label), dtype=object))
                buffered_rows += n
                if buffered_rows >= flush_rows:
                    flush()
        flush()
    finally:
        if writer is not None:
            writer.close()

    if written_rows <= 0:
        try:
            out_path.unlink(missing_ok=True)
        except Exception:
            pass
        return None
    return str(out_path)


def _process_chunk_edge_index_hist_task(args):
    (
        chunk_path,
        chunk_idx,
        chunk_size,
        te_cutoff,
        fdr,
        labels,
        output_dir,
        max_len,
        batch_rows,
    ) = args
    chunk_path = Path(chunk_path)
    chunk_idx = int(chunk_idx)
    chunk_size = int(chunk_size)
    pf = pq.ParquetFile(chunk_path)
    raw_chunk_len = _infer_chunk_length(pf)
    if raw_chunk_len <= 0:
        return None
    global_offset = chunk_idx * chunk_size
    if global_offset >= int(max_len):
        return None
    remaining = int(max_len) - global_offset
    chunk_len = min(raw_chunk_len, remaining)
    thresholds, _ = _histogram_thresholds_for_chunk(
        pf,
        raw_chunk_len=raw_chunk_len,
        chunk_len=chunk_len,
        te_cutoff=float(te_cutoff),
        fdr=float(fdr),
        batch_rows=int(batch_rows),
    )
    finite_thresholds = np.isfinite(thresholds)
    if not np.any(finite_thresholds):
        return None

    edge_parts: List[List[np.ndarray]] = [[] for _ in range(chunk_len)]
    te_parts: List[List[np.ndarray]] = [[] for _ in range(chunk_len)]
    row_offset = 0
    use_arrow_threads = _arrow_use_threads()
    for batch in pf.iter_batches(columns=["Values"], batch_size=int(batch_rows), use_threads=use_arrow_threads):
        values = batch.column("Values")
        flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
        if flat.size == 0:
            row_offset += batch.num_rows
            continue
        mat = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)[:, :chunk_len]
        if mat.size == 0:
            row_offset += batch.num_rows
            continue
        mask = (mat >= thresholds.reshape(1, -1)) & finite_thresholds.reshape(1, -1)
        if np.any(mask):
            row_idx, local_idx = np.nonzero(mask)
            if row_idx.size:
                order = np.argsort(local_idx, kind="stable")
                row_idx = row_idx[order]
                local_idx = local_idx[order]
                split_points = np.flatnonzero(np.diff(local_idx)) + 1
                row_groups = np.split(row_idx, split_points)
                local_groups = np.split(local_idx, split_points)
                for rows_for_time, locals_for_time in zip(row_groups, local_groups):
                    local = int(locals_for_time[0])
                    edge_parts[local].append((rows_for_time.astype(np.int32, copy=False) + np.int32(row_offset)))
                    te_parts[local].append(mat[rows_for_time, local].astype(np.float32, copy=False))
        row_offset += batch.num_rows

    times: List[int] = []
    edge_lists: List[np.ndarray] = []
    te_lists: List[np.ndarray] = []
    for local_idx in range(chunk_len):
        if not edge_parts[local_idx]:
            continue
        edges = np.concatenate(edge_parts[local_idx]).astype(np.int32, copy=False)
        tes = np.concatenate(te_parts[local_idx]).astype(np.float32, copy=False)
        if edges.size == 0:
            continue
        global_idx = global_offset + local_idx
        if global_idx < len(labels):
            try:
                time_value = int(labels[global_idx])
            except Exception:
                time_value = int(global_idx)
        else:
            time_value = int(global_idx)
        times.append(time_value)
        edge_lists.append(edges)
        te_lists.append(tes)

    if not times:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"grn_chunk_{chunk_idx:04d}.parquet"
    schema = pa.schema(
        [
            ("Time", pa.int32()),
            ("EdgeIndex", pa.list_(pa.int32())),
            ("TE", pa.list_(pa.float32())),
        ]
    )
    table = pa.table(
        {
            "Time": pa.array(np.asarray(times, dtype=np.int32), type=pa.int32()),
            "EdgeIndex": pa.array(edge_lists, type=pa.list_(pa.int32())),
            "TE": pa.array(te_lists, type=pa.list_(pa.float32())),
        },
        schema=schema,
    )
    pq.write_table(table, out_path, compression="snappy")
    return str(out_path)


def _write_edge_index_catalog_from_chunk(
    chunk_path: Path,
    catalog_out: Path,
    gene_names: Sequence[str],
    batch_rows: int,
) -> None:
    """Write EdgeIndex -> Source/Target catalog using chunk row order.

    EdgeIndex values in edge_index outputs are row offsets within each LocalTE
    chunk parquet. Therefore the catalog must be generated from the chunk
    Source/Target row order, not from a separate score table that may have a
    different ordering.
    """
    pf = pq.ParquetFile(chunk_path)
    gene_names_arr = np.asarray(gene_names, dtype=object)
    schema = pa.schema(
        [
            ("EdgeIndex", pa.int32()),
            ("SourceIndex", pa.int32()),
            ("TargetIndex", pa.int32()),
            ("Source", pa.string()),
            ("Target", pa.string()),
        ]
    )
    writer: Optional[pq.ParquetWriter] = None
    edge_offset = 0
    try:
        for batch in pf.iter_batches(columns=["Source", "Target"], batch_size=int(batch_rows), use_threads=_arrow_use_threads()):
            src_idx = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
            tgt_idx = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
            n = int(src_idx.size)
            edge_idx = np.arange(edge_offset, edge_offset + n, dtype=np.int32)
            src_names = gene_names_arr[np.clip(src_idx.astype(np.int64) - 1, 0, gene_names_arr.size - 1)]
            tgt_names = gene_names_arr[np.clip(tgt_idx.astype(np.int64) - 1, 0, gene_names_arr.size - 1)]
            table = pa.table(
                {
                    "EdgeIndex": pa.array(edge_idx, type=pa.int32()),
                    "SourceIndex": pa.array(src_idx, type=pa.int32()),
                    "TargetIndex": pa.array(tgt_idx, type=pa.int32()),
                    "Source": pa.array(src_names.astype(object), type=pa.string()),
                    "Target": pa.array(tgt_names.astype(object), type=pa.string()),
                },
                schema=schema,
            )
            if writer is None:
                writer = pq.ParquetWriter(str(catalog_out), schema, compression="snappy")
            writer.write_table(table)
            edge_offset += n
    finally:
        if writer is not None:
            writer.close()


def _trim_timepoint_edges_for_outdegree(
    df_time: pd.DataFrame,
    threshold: float,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    if df_time.empty:
        return df_time, {"sources": 0, "edges": 0, "kept": 0, "removed": 0}
    if _trim_indirect_edges is None:
        raise RuntimeError("trim_indirect helper is unavailable; cannot run outdegree-direct trimming.")
    source_order: List[str] = []
    targets_by_source: Dict[str, List[str]] = {}
    tes_by_source: Dict[str, List[float]] = {}
    sources = df_time["Source"].astype(str).to_numpy()
    targets = df_time["Target"].astype(str).to_numpy()
    tes = df_time["TE"].astype(float).to_numpy()
    for source, target, te in zip(sources, targets, tes):
        if source not in targets_by_source:
            source_order.append(source)
            targets_by_source[source] = []
            tes_by_source[source] = []
        targets_by_source[source].append(target)
        tes_by_source[source].append(float(te))
    indirect = _trim_indirect_edges(source_order, targets_by_source, tes_by_source, threshold)
    keep_mask = np.ones(len(df_time), dtype=bool)
    cursor_by_source: Dict[str, int] = {source: 0 for source in source_order}
    for row_idx, source in enumerate(sources):
        pos = cursor_by_source[source]
        cursor_by_source[source] = pos + 1
        if indirect[source][pos]:
            keep_mask[row_idx] = False
    trimmed = df_time.loc[keep_mask].copy()
    return trimmed, {
        "sources": len(source_order),
        "edges": int(len(df_time)),
        "kept": int(keep_mask.sum()),
        "removed": int((~keep_mask).sum()),
    }


def _process_chunk_outdegree_task(args):
    (
        chunk_path,
        chunk_idx,
        chunk_size,
        gene_names,
        te_cutoff,
        fdr,
        labels,
        orig_idx_seq,
        output_dir,
        max_len,
        batch_rows,
        temp_root,
        trim_threshold,
    ) = args
    gene_names_arr = np.asarray(gene_names, dtype=object)
    frames = _process_chunk_streaming(
        chunk_idx=int(chunk_idx),
        chunk_path=Path(chunk_path),
        chunk_size=int(chunk_size),
        te_cutoff=float(te_cutoff),
        fdr=float(fdr),
        labels=labels,
        orig_idx_seq=orig_idx_seq,
        gene_names_arr=gene_names_arr,
        max_len=int(max_len),
        batch_rows=int(batch_rows),
        temp_root=Path(temp_root),
    )
    if not frames:
        return None
    out_rows: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, object]] = []
    do_trim = _env_flag("LOCALTE_TRIM_INDIRECT", "yes")
    for df_time in frames:
        if df_time is None or df_time.empty:
            continue
        time_label = str(df_time["TimeLabel"].iat[0])
        if do_trim:
            df_use, summary = _trim_timepoint_edges_for_outdegree(df_time, float(trim_threshold))
        else:
            df_use = df_time
            summary = {
                "sources": int(df_time["Source"].nunique()),
                "edges": int(len(df_time)),
                "kept": int(len(df_time)),
                "removed": 0,
            }
        summary_rows.append({"Time": time_label, **summary})
        if df_use.empty:
            continue
        counts = df_use.groupby("Source", sort=False).size().reset_index(name="Outdegree")
        counts.insert(0, "Time", time_label)
        out_rows.append(counts)
    if not out_rows:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"outdegree_chunk_{int(chunk_idx):04d}.parquet"
    summary_path = out_dir / f"summary_chunk_{int(chunk_idx):04d}.parquet"
    pd.concat(out_rows, ignore_index=True).to_parquet(out_path, index=False)
    if summary_rows:
        pd.DataFrame(summary_rows).to_parquet(summary_path, index=False)
    return str(out_path)


def _process_chunk_outdegree_fast_task(args):
    (
        chunk_path,
        chunk_idx,
        chunk_size,
        gene_names,
        te_cutoff,
        fdr,
        labels,
        output_dir,
        max_len,
        batch_rows,
        temp_root,
        trim_threshold,
    ) = args
    chunk_idx = int(chunk_idx)
    chunk_size = int(chunk_size)
    max_len = int(max_len)
    pf = pq.ParquetFile(Path(chunk_path))
    raw_chunk_len = _infer_chunk_length(pf)
    if raw_chunk_len <= 0:
        return None
    global_offset = chunk_idx * chunk_size
    if global_offset >= max_len:
        return None
    chunk_len = min(raw_chunk_len, max_len - global_offset)
    chunk_tmp = Path(temp_root) / f"outdegree_chunk_{chunk_idx:04d}"
    chunk_tmp.mkdir(parents=True, exist_ok=True)
    counts = np.zeros(chunk_len, dtype=np.int64)
    sum_vals = np.zeros(chunk_len, dtype=np.float64)
    sum_sq = np.zeros(chunk_len, dtype=np.float64)
    handles: List[Optional[object]] = []
    for local_idx in range(chunk_len):
        handles.append((chunk_tmp / f"time_{local_idx:03d}.bin").open("wb"))
    try:
        use_arrow_threads = _arrow_use_threads()
        for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=int(batch_rows), use_threads=use_arrow_threads):
            values = batch.column("Values")
            flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
            if flat.size == 0:
                continue
            mat = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)[:, :chunk_len]
            if mat.size == 0:
                continue
            mask = mat > float(te_cutoff)
            if not np.any(mask):
                continue
            src_arr = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
            tgt_arr = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
            for local_idx in range(chunk_len):
                col_mask = mask[:, local_idx]
                if not np.any(col_mask):
                    continue
                vals = mat[col_mask, local_idx]
                n = vals.size
                counts[local_idx] += n
                sum_vals[local_idx] += float(np.sum(vals, dtype=np.float64))
                sum_sq[local_idx] += float(np.sum(vals * vals, dtype=np.float64))
                arr = np.empty(n, dtype=_TRIPLE_DTYPE)
                arr["TE"] = vals.astype(np.float32, copy=False)
                arr["Source"] = src_arr[col_mask]
                arr["Target"] = tgt_arr[col_mask]
                handles[local_idx].write(arr.tobytes())
    finally:
        for handle in handles:
            try:
                handle.close()
            except Exception:
                pass

    gene_names_arr = np.asarray(gene_names, dtype=object)
    out_rows: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, object]] = []
    try:
        for local_idx in range(chunk_len):
            global_idx = global_offset + local_idx
            if global_idx >= len(labels):
                continue
            cnt = int(counts[local_idx])
            if cnt < 2:
                continue
            variance = (sum_sq[local_idx] / cnt) - ((sum_vals[local_idx] / cnt) ** 2)
            if variance <= 1e-12:
                continue
            file_path = chunk_tmp / f"time_{local_idx:03d}.bin"
            if not file_path.exists() or file_path.stat().st_size == 0:
                continue
            triples = np.memmap(file_path, dtype=_TRIPLE_DTYPE, mode="r")
            try:
                filtered = _filter_timepoint_arrays(
                    np.asarray(triples["TE"], dtype=np.float32),
                    np.asarray(triples["Source"], dtype=np.int32),
                    np.asarray(triples["Target"], dtype=np.int32),
                    float(sum_vals[local_idx] / cnt),
                    math.sqrt(max(float(variance), 0.0)),
                    float(fdr),
                )
            finally:
                del triples
            if filtered is None:
                continue
            src_keep, tgt_keep, te_keep = filtered
            src_out, outdegree, summary = _numeric_trim_outdegree(
                src_keep,
                tgt_keep,
                te_keep,
                float(trim_threshold),
            )
            time_label = labels[global_idx] if global_idx < len(labels) else str(global_idx)
            summary_rows.append({"Time": str(time_label), **summary})
            if src_out.size:
                src_names = gene_names_arr[np.clip(src_out.astype(np.int64) - 1, 0, gene_names_arr.size - 1)]
                out_rows.append(
                    pd.DataFrame(
                        {
                            "Time": np.full(src_out.size, str(time_label), dtype=object),
                            "Source": src_names,
                            "Outdegree": outdegree,
                        }
                    )
                )
    finally:
        shutil.rmtree(chunk_tmp, ignore_errors=True)
    if not out_rows:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"outdegree_chunk_{chunk_idx:04d}.parquet"
    summary_path = out_dir / f"summary_chunk_{chunk_idx:04d}.parquet"
    pd.concat(out_rows, ignore_index=True).to_parquet(out_path, index=False)
    if summary_rows:
        pd.DataFrame(summary_rows).to_parquet(summary_path, index=False)
    return str(out_path)


def _process_chunk_outdegree_hist_task(args):
    (
        chunk_path,
        chunk_idx,
        chunk_size,
        gene_names,
        te_cutoff,
        fdr,
        labels,
        output_dir,
        max_len,
        batch_rows,
        trim_threshold,
    ) = args
    chunk_idx = int(chunk_idx)
    chunk_size = int(chunk_size)
    max_len = int(max_len)
    if isinstance(chunk_path, (list, tuple)):
        chunk_paths = [Path(p) for p in chunk_path]
    else:
        chunk_paths = [Path(chunk_path)]
    pfs = [pq.ParquetFile(path) for path in chunk_paths]
    raw_chunk_len = _infer_chunk_length(pfs[0])
    if raw_chunk_len <= 0:
        return None
    for pf_check, path in zip(pfs[1:], chunk_paths[1:]):
        part_len = _infer_chunk_length(pf_check)
        if part_len != raw_chunk_len:
            raise ValueError(
                f"Worker-part chunk length mismatch for chunk_{chunk_idx:04d}: "
                f"{chunk_paths[0]} has {raw_chunk_len}, {path} has {part_len}"
            )
    global_offset = chunk_idx * chunk_size
    if global_offset >= max_len:
        return None
    chunk_len = min(raw_chunk_len, max_len - global_offset)
    thresholds, _ = _histogram_thresholds_for_chunk_files(
        pfs,
        raw_chunk_len=raw_chunk_len,
        chunk_len=chunk_len,
        te_cutoff=float(te_cutoff),
        fdr=float(fdr),
        batch_rows=int(batch_rows),
    )
    if not np.any(np.isfinite(thresholds)):
        return None

    src_parts: List[List[np.ndarray]] = [[] for _ in range(chunk_len)]
    tgt_parts: List[List[np.ndarray]] = [[] for _ in range(chunk_len)]
    te_parts: List[List[np.ndarray]] = [[] for _ in range(chunk_len)]
    use_arrow_threads = _arrow_use_threads()
    for pf in pfs:
        for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=int(batch_rows), use_threads=use_arrow_threads):
            values = batch.column("Values")
            flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
            if flat.size == 0:
                continue
            mat = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)[:, :chunk_len]
            if mat.size == 0:
                continue
            src_arr = np.asarray(batch.column("Source").to_numpy(), dtype=np.int32)
            tgt_arr = np.asarray(batch.column("Target").to_numpy(), dtype=np.int32)
            for local_idx in range(chunk_len):
                threshold = thresholds[local_idx]
                if not np.isfinite(threshold):
                    continue
                col = mat[:, local_idx]
                mask = col >= threshold
                if not np.any(mask):
                    continue
                src_parts[local_idx].append(src_arr[mask])
                tgt_parts[local_idx].append(tgt_arr[mask])
                te_parts[local_idx].append(col[mask].astype(np.float32, copy=False))

    gene_names_arr = np.asarray(gene_names, dtype=object)
    out_rows: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, object]] = []
    do_trim = _env_flag("LOCALTE_TRIM_INDIRECT", "yes")
    for local_idx in range(chunk_len):
        if not src_parts[local_idx]:
            continue
        src_keep = np.concatenate(src_parts[local_idx]).astype(np.int32, copy=False)
        tgt_keep = np.concatenate(tgt_parts[local_idx]).astype(np.int32, copy=False)
        te_keep = np.concatenate(te_parts[local_idx]).astype(np.float32, copy=False)
        if do_trim:
            src_out, outdegree, summary = _numeric_trim_outdegree(
                src_keep,
                tgt_keep,
                te_keep,
                float(trim_threshold),
            )
        else:
            src_out, outdegree = np.unique(src_keep, return_counts=True)
            src_out = src_out.astype(np.int32, copy=False)
            outdegree = outdegree.astype(np.int64, copy=False)
            summary = {
                "sources": int(src_out.size),
                "edges": int(src_keep.size),
                "kept": int(src_keep.size),
                "removed": 0,
            }
        global_idx = global_offset + local_idx
        time_label = labels[global_idx] if global_idx < len(labels) else str(global_idx)
        summary_rows.append({"Time": str(time_label), **summary})
        if src_out.size:
            src_names = gene_names_arr[np.clip(src_out.astype(np.int64) - 1, 0, gene_names_arr.size - 1)]
            out_rows.append(
                pd.DataFrame(
                    {
                        "Time": np.full(src_out.size, str(time_label), dtype=object),
                        "Source": src_names,
                        "Outdegree": outdegree,
                    }
                )
            )
    if not out_rows:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"outdegree_chunk_{chunk_idx:04d}.parquet"
    summary_path = out_dir / f"summary_chunk_{chunk_idx:04d}.parquet"
    pd.concat(out_rows, ignore_index=True).to_parquet(out_path, index=False)
    if summary_rows:
        pd.DataFrame(summary_rows).to_parquet(summary_path, index=False)
    return str(out_path)


def _process_chunk_streaming(
    chunk_idx: int,
    chunk_path: Path,
    chunk_size: int,
    te_cutoff: float,
    fdr: float,
    labels: Sequence[str],
    orig_idx_seq: Optional[Sequence[int]],
    gene_names_arr: np.ndarray,
    max_len: int,
    batch_rows: int,
    temp_root: Path,
) -> List[pd.DataFrame]:
    pf = pq.ParquetFile(chunk_path)
    raw_chunk_len = _infer_chunk_length(pf)
    if raw_chunk_len <= 0:
        return []
    global_offset = chunk_idx * chunk_size
    if global_offset >= max_len:
        return []
    remaining = max_len - global_offset
    chunk_len = min(raw_chunk_len, remaining)
    if _use_single_pass_mode():
        return _process_chunk_single_pass(
            pf=pf,
            chunk_idx=chunk_idx,
            raw_chunk_len=raw_chunk_len,
            chunk_len=chunk_len,
            chunk_size=chunk_size,
            te_cutoff=te_cutoff,
            fdr=fdr,
            labels=labels,
            orig_idx_seq=orig_idx_seq,
            gene_names_arr=gene_names_arr,
            max_len=max_len,
            batch_rows=batch_rows,
            temp_root=temp_root,
        )

    chunk_tmp = temp_root / f"chunk_{chunk_idx:04d}"
    chunk_tmp.mkdir(parents=True, exist_ok=True)
    counts, sum_vals, sum_sq = _compute_chunk_stats(pf, raw_chunk_len, chunk_len, te_cutoff, batch_rows)
    if counts.sum() == 0:
        shutil.rmtree(chunk_tmp, ignore_errors=True)
        return []
    offsets, te_mem, src_mem, tgt_mem = _fill_chunk_memmaps(
        pf,
        raw_chunk_len,
        chunk_len,
        te_cutoff,
        batch_rows,
        counts,
        chunk_tmp,
    )
    results: List[Tuple[int, pd.DataFrame]] = []

    def _build_for_local_idx(local_idx: int) -> Optional[Tuple[int, pd.DataFrame]]:
        global_idx = global_offset + local_idx
        if global_idx >= len(labels):
            return None
        start = int(offsets[local_idx])
        end = int(offsets[local_idx + 1])
        n = end - start
        if n < 2:
            return None
        cnt = float(counts[local_idx])
        if cnt < 2:
            return None
        total = float(sum_vals[local_idx])
        mean = total / cnt if cnt else 0.0
        variance = (sum_sq[local_idx] / cnt) - (mean ** 2)
        if variance <= 1e-12:
            return None
        std = math.sqrt(max(variance, 0.0))
        te_vals = np.asarray(te_mem[start:end], dtype=np.float32)
        src_vals = np.asarray(src_mem[start:end], dtype=np.int32)
        tgt_vals = np.asarray(tgt_mem[start:end], dtype=np.int32)
        df_time = _build_timepoint_df(
            te_vals,
            src_vals,
            tgt_vals,
            mean,
            std,
            global_idx,
            labels,
            orig_idx_seq,
            gene_names_arr,
            fdr,
        )
        if df_time is None or df_time.empty:
            return None
        return global_idx, df_time

    try:
        thread_workers = _chunk_thread_workers()
        if thread_workers > 1 and chunk_len > 1:
            with ThreadPoolExecutor(max_workers=min(thread_workers, chunk_len)) as executor:
                for result in executor.map(_build_for_local_idx, range(chunk_len)):
                    if result is not None:
                        results.append(result)
        else:
            for local_idx in range(chunk_len):
                result = _build_for_local_idx(local_idx)
                if result is not None:
                    results.append(result)
    finally:
        if te_mem is not None:
            del te_mem
        if src_mem is not None:
            del src_mem
        if tgt_mem is not None:
            del tgt_mem
        shutil.rmtree(chunk_tmp, ignore_errors=True)
    if not results:
        return []
    results.sort(key=lambda x: x[0])
    return [df for _, df in results]


def _process_chunk_file_task(args):
    (
        chunk_path,
        chunk_idx,
        chunk_size,
        gene_names,
        te_cutoff,
        fdr,
        labels,
        orig_idx_seq,
        output_dir,
        max_len,
        batch_rows,
        temp_root,
    ) = args
    gene_names_arr = np.asarray(gene_names, dtype=object)
    frames = _process_chunk_streaming(
        chunk_idx=chunk_idx,
        chunk_path=Path(chunk_path),
        chunk_size=chunk_size,
        te_cutoff=te_cutoff,
        fdr=fdr,
        labels=labels,
        orig_idx_seq=orig_idx_seq,
        gene_names_arr=gene_names_arr,
        max_len=max_len,
        batch_rows=batch_rows,
        temp_root=Path(temp_root),
    )
    if not frames:
        return None
    part = pd.concat(frames, ignore_index=True)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"grn_chunk_{chunk_idx:04d}.parquet"
    part.to_parquet(out_path, index=False)
    return str(out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create timepoint GRNs from compressed LocalTE columns or chunked payloads.")
    parser.add_argument("--input", required=True, help="Input Parquet file containing TE scores (with or without LocalTE payload columns).")
    parser.add_argument("--output_prefix", default=None, help="Optional prefix for output file when --output_file is not set.")
    parser.add_argument("--output_file", help="Optional output file path (parquet/csv). Defaults to <output_prefix or localte_times>.parquet.")
    parser.add_argument("--fdr", type=float, default=0.05, help="BH-FDR alpha threshold (default 0.05).")
    parser.add_argument("--te_cutoff", type=float, default=0.0, help="Include edges with TE strictly greater than this value (default 0.0).")
    parser.add_argument("--history_length", type=int, default=1, help="History length k (for aligning time labels).")
    parser.add_argument("--time_labels", help="Optional text file with one pseudotime label per filtered cell.")
    parser.add_argument("--max_timepoints", type=int, help="Limit number of time indices to evaluate (debug speedup).")
    parser.add_argument("--time_map", help="Optional Parquet file with columns [LocalIndex, OriginalIndex] from the TE run.")
    parser.add_argument(
        "--gene_names",
        help=(
            "Optional gene_names file. If omitted, the script first tries "
            "<input parquet directory>/gene_names, then falls back to the package input directory."
        ),
    )
    parser.add_argument("--workers", type=int, default=0, help="Worker processes (0/1 sequential, >1 multiprocessing).")
    parser.add_argument("--chunk_dir", help="Directory containing chunked LocalTE parquet files (from localte_chunk_export.py).")
    parser.add_argument("--chunk_size", type=int, help="Timepoints per chunk when using --chunk_dir (default: metadata.json value).")
    parser.add_argument("--chunk_workers", type=int, default=0, help="Parallel workers for chunk_dir processing (0/1=single-process; >1 parallel by chunk). Only used when --output_file is set.")
    parser.add_argument(
        "--output_mode",
        choices=["edge", "edge_index", "outdegree"],
        default=os.environ.get("LOCALTE_OUTPUT_MODE", "edge"),
        help=(
            "edge writes materialized timepoint GRN edges; "
            "edge_index writes compact per-timepoint EdgeIndex/TE lists plus an edge catalog; "
            "outdegree writes per-timepoint source outdegree after optional indirect trimming."
        ),
    )
    args = parser.parse_args()

    def resolve_gene_names_path() -> Path:
        if args.gene_names:
            path = coerce_input_path(args.gene_names)
            if not path.exists():
                raise FileNotFoundError(f"gene_names file not found: {path}")
            return path
        input_neighbor = Path(args.input).resolve().parent / "gene_names"
        if input_neighbor.exists():
            return input_neighbor
        path = coerce_input_path("gene_names")
        if not path.exists():
            raise FileNotFoundError(
                "gene_names file not found. Provide --gene_names or run from an output directory containing gene_names."
            )
        return path

    gene_names_path = resolve_gene_names_path()

    output_file_path: Optional[Path] = None
    if args.output_file:
        output_file_path = coerce_output_path(args.output_file)
        output_file_path.parent.mkdir(parents=True, exist_ok=True)

    def load_gene_index_map() -> Dict[str, int]:
        gene_file = gene_names_path
        if not gene_file.exists():
            raise FileNotFoundError("gene_names file not found; required to map gene names to indices.")
        mapping: Dict[str, int] = {}
        with gene_file.open(encoding="utf-8") as handle:
            for idx, line in enumerate(handle, start=1):
                mapping[line.strip()] = idx
        return mapping

    def ensure_indices(series: pd.Series, gene_to_idx: Optional[Dict[str, int]] = None) -> np.ndarray:
        indices: List[int] = []
        need_map = False
        for value in series:
            if isinstance(value, (int, np.integer)):
                indices.append(int(value))
                continue
            try:
                indices.append(int(str(value)))
            except Exception:
                need_map = True
                break
        if not need_map:
            return np.asarray(indices, dtype=int)
        if gene_to_idx is None:
            gene_to_idx = load_gene_index_map()
        indices = []
        for value in series:
            if isinstance(value, (int, np.integer)):
                indices.append(int(value))
                continue
            try:
                indices.append(int(str(value)))
                continue
            except Exception:
                name = str(value)
                if name not in gene_to_idx:
                    raise ValueError(f"Gene/peak name '{name}' not found in gene_names mapping.")
                indices.append(gene_to_idx[name])
        return np.asarray(indices, dtype=int)

    try:
        pq_file = pq.ParquetFile(args.input)
    except Exception as exc:
        print(f"Failed to open {args.input}: {exc}", file=sys.stderr)
        sys.exit(1)

    available_cols = set(pq_file.schema.names)
    base_cols = ["Source", "Target", "TE"]
    cols_to_read = base_cols.copy()
    length_col = None
    if "LocalTE_total_len" in available_cols:
        cols_to_read.append("LocalTE_total_len")
        length_col = "LocalTE_total_len"
    if "LocalTE_len" in available_cols and "LocalTE_len" not in cols_to_read:
        cols_to_read.append("LocalTE_len")
        if length_col is None:
            length_col = "LocalTE_len"

    chunk_dir = None
    if args.chunk_dir:
        raw_chunk_dir = Path(args.chunk_dir)
        if not raw_chunk_dir.exists() and not raw_chunk_dir.is_absolute():
            raw_chunk_dir = resolve_output(raw_chunk_dir)
        chunk_dir = raw_chunk_dir.resolve()
        if not chunk_dir.exists():
            print(f"chunk_dir {chunk_dir} does not exist.", file=sys.stderr)
            sys.exit(1)
    else:
        # Try default split-chunk directory if present
        try_default = resolve_output("local_te_split_chunks")
        if (try_default / "metadata.json").exists():
            chunk_dir = try_default.resolve()
        else:
            # Backward-compatible fallback to legacy directory name
            legacy_default = resolve_output("local_te_chunks")
            if (legacy_default / "metadata.json").exists():
                chunk_dir = legacy_default.resolve()

    if chunk_dir is None:
        needed = ["LocalTE_bytes", "LocalTE_len", "LocalTE_dtype"]
        missing = [c for c in needed if c not in available_cols]
        if missing:
            print(f"Input parquet missing required columns: {missing}", file=sys.stderr)
            sys.exit(1)
        for col in needed:
            if col not in cols_to_read:
                cols_to_read.append(col)

    try:
        df = pd.read_parquet(args.input, columns=cols_to_read)
    except Exception as exc:
        print(f"Failed to load {args.input}: {exc}", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print("Input dataframe is empty; nothing to do.")
        return

    if length_col and length_col in df.columns:
        lengths = df[length_col].fillna(0).astype(int).to_numpy()
    elif "LocalTE_len" in df.columns:
        lengths = df["LocalTE_len"].fillna(0).astype(int).to_numpy()
    else:
        lengths = np.zeros(len(df), dtype=int)

    max_len = int(lengths.max()) if lengths.size else 0
    chunk_meta: Dict[str, object] = {}
    chunk_size = int(args.chunk_size) if args.chunk_size else 300
    if chunk_dir is not None:
        meta_path = chunk_dir / "metadata.json"
        direct_marker_path = chunk_dir / "direct_localte_done.json"
        if meta_path.exists():
            with meta_path.open() as handle:
                chunk_meta = json.load(handle)
        elif direct_marker_path.exists():
            with direct_marker_path.open() as handle:
                chunk_meta = json.load(handle)
        elif (chunk_dir.parent / "direct_localte_done.json").exists():
            with (chunk_dir.parent / "direct_localte_done.json").open() as handle:
                chunk_meta = json.load(handle)
        elif list(chunk_dir.glob("parts_w*/metadata.json")):
            with sorted(chunk_dir.glob("parts_w*/metadata.json"))[0].open() as handle:
                chunk_meta = json.load(handle)
        else:
            print(f"Expected metadata.json or direct_localte_done.json in {chunk_dir}", file=sys.stderr)
            sys.exit(1)
        if args.chunk_size is None:
            chunk_size = int(chunk_meta.get("chunk_size", chunk_size))
        max_len = max(max_len, int(chunk_meta.get("max_timepoints", 0)))
        if max_len == 0:
            num_chunks_meta = int(chunk_meta.get("num_chunks", 0))
            max_len = chunk_size * num_chunks_meta
        chunk_map_preview = _discover_chunk_files(chunk_dir)
        if max_len == 0 and chunk_map_preview:
            max_chunk_idx = max(chunk_map_preview)
            try:
                last_pf = pq.ParquetFile(chunk_map_preview[max_chunk_idx][0])
                last_len = _infer_chunk_length(last_pf)
            except Exception:
                last_len = chunk_size
            max_len = (int(max_chunk_idx) * int(chunk_size)) + int(last_len or chunk_size)
        if max_len == 0:
            max_len = chunk_size
    if args.max_timepoints is not None:
        max_len = min(max_len, int(args.max_timepoints))

    if max_len <= 0:
        print("No LocalTE entries found in the input file.")
        return

    # Optional mapping LocalIndex->OriginalIndex
    mapping = None
    if args.time_map:
        try:
            time_map_path = coerce_input_path(args.time_map)
            mdf = pd.read_parquet(time_map_path)
            if not {"LocalIndex", "OriginalIndex"}.issubset(set(mdf.columns)):
                raise ValueError("time_map missing required columns LocalIndex, OriginalIndex")
            mapping = mdf.sort_values("LocalIndex").reset_index(drop=True)
        except Exception as exc:
            print(f"Failed to load time_map {args.time_map}: {exc}", file=sys.stderr)
            sys.exit(1)

    # Build labels sequence for each LocalIndex
    if args.time_labels:
        # If mapping is present, index labels by OriginalIndex; otherwise fallback to contiguous [k:]
        if mapping is not None:
            # load full labels (for all original timepoints) and select by OriginalIndex series length max_len
            labels_path = coerce_input_path(args.time_labels)
            with open(labels_path, encoding="utf-8") as handle:
                full_labels = [line.strip() for line in handle]
            orig_idx = mapping["OriginalIndex"].to_numpy()
            try:
                labels = [str(full_labels[int(i)]) for i in orig_idx[:max_len]]
            except Exception as exc:
                print(f"Failed to map time labels via time_map: {exc}", file=sys.stderr)
                sys.exit(1)
        else:
            labels = load_time_labels(args.time_labels, args.history_length, max_len)
    else:
        # Default string labels: prefer OriginalIndex if available
        if mapping is not None:
            labels = [str(int(i)) for i in mapping["OriginalIndex"].to_numpy()[:max_len]]
        else:
            labels = [str(i) for i in range(max_len)]

    orig_idx_seq: Optional[Sequence[int]] = None
    if mapping is not None:
        arr = mapping["OriginalIndex"].to_numpy()
        if arr.size < max_len:
            print("time_map is shorter than needed for LocalTE length.", file=sys.stderr)
            sys.exit(1)
        orig_idx_seq = [int(x) for x in arr[:max_len]]

    tuples: List[Tuple[int, pd.DataFrame]] = []

    if chunk_dir is not None:
        num_chunks_eval = math.ceil(max_len / chunk_size)
        chunk_files: List[Tuple[int, Path]] = []
        temp_parent = _resolve_temp_parent(output_file_path)
        merged_parts_dir: Optional[Path] = None
        chunk_map = _discover_chunk_files(chunk_dir)
        chunk_map = {idx: files for idx, files in chunk_map.items() if idx < num_chunks_eval}
        needs_part_merge = any(len(files) > 1 for files in chunk_map.values())
        direct_part_outdegree_hist = (
            needs_part_merge
            and args.output_mode == "outdegree"
            and _env_flag("LOCAL_TE_OUTDEGREE_HIST", "on")
        )
        if needs_part_merge and direct_part_outdegree_hist:
            print(
                f"[LocalTE GRN] Detected worker-part chunks in {chunk_dir}; "
                "processing parts directly for outdegree histogram mode."
            )
            for chunk_idx in range(num_chunks_eval):
                files = chunk_map.get(chunk_idx, [])
                if files:
                    chunk_files.append((chunk_idx, list(files)))
        elif needs_part_merge:
            print(
                f"[LocalTE GRN] Detected direct worker-part chunks in {chunk_dir}; "
                "temporarily merging parts by chunk for GRN processing."
            )
            merged_parts_dir, chunk_files = _merge_chunk_parts_for_processing(chunk_map, temp_parent)
        else:
            for chunk_idx in range(num_chunks_eval):
                files = chunk_map.get(chunk_idx, [])
                if files:
                    chunk_files.append((chunk_idx, files[0]))
        def _cleanup_merged_parts() -> None:
            if merged_parts_dir is not None:
                shutil.rmtree(merged_parts_dir, ignore_errors=True)

        if not chunk_files:
            print(f"No chunk parquet files found in {chunk_dir}", file=sys.stderr)
            _cleanup_merged_parts()
            return
        chunk_total_bytes = 0
        for _, cp in chunk_files:
            paths_for_size = cp if isinstance(cp, list) else [cp]
            for part_path in paths_for_size:
                try:
                    chunk_total_bytes += Path(part_path).stat().st_size
                except OSError:
                    pass

        gene_names_list = _load_gene_names_list(str(gene_names_path))
        gene_names_arr = np.asarray(gene_names_list, dtype=object)
        batch_rows = _chunk_batch_rows()
        chunk_workers = int(getattr(args, "chunk_workers", 0) or 0)
        histogram_mode = _use_histogram_mode() and args.output_mode == "edge"
        edge_index_mode = args.output_mode == "edge_index"
        outdegree_hist_mode = args.output_mode == "outdegree" and _env_flag("LOCAL_TE_OUTDEGREE_HIST", "on")
        est_temp_per_worker = 0 if (histogram_mode or edge_index_mode or outdegree_hist_mode) else _estimate_chunk_temp_bytes(chunk_meta, chunk_size, max_len)
        n_requested = min(max(chunk_workers, 1), multiprocessing.cpu_count(), len(chunk_files))
        n_workers = _cap_workers_by_temp(n_requested, temp_parent, est_temp_per_worker)
        print(
            "[LocalTE GRN] Starting timepoint GRN generation: "
            f"chunks={_human_int(len(chunk_files))}, "
            f"chunk_size={_human_int(chunk_size)}, "
            f"timepoints={_human_int(max_len)}, "
            f"input_size={_human_bytes(chunk_total_bytes)}, "
            f"workers={_human_int(n_workers)}, "
            f"temp_dir={temp_parent}, "
            f"estimated_temp_per_worker={_human_bytes(est_temp_per_worker)}, "
            f"histogram_fdr={'on' if histogram_mode else 'off'}, "
            f"edge_index={'on' if edge_index_mode else 'off'}, "
            f"outdegree_hist={'on' if outdegree_hist_mode else 'off'}, "
            f"edge_numeric={'on' if _use_numeric_edge_output() else 'off'}, "
            f"single_pass={'on' if _use_single_pass_mode() else 'off'}, "
            f"arrow_threads={'on' if _arrow_use_threads() else 'off'}, "
            f"batch_rows={_human_int(batch_rows)}",
            flush=True,
        )
        print(
            f"[LocalTE GRN] gene_names={gene_names_path} "
            f"({ _human_int(len(gene_names_list)) } entries)",
            flush=True,
        )
        if output_file_path:
            print(f"[LocalTE GRN] output={output_file_path}", flush=True)
        _warm_hist_numba()
        temp_root_path: Optional[Path] = None
        if not histogram_mode and not edge_index_mode and not outdegree_hist_mode:
            temp_root_path = Path(tempfile.mkdtemp(prefix="grn_stream_", dir=str(temp_parent)))

        try:
            if output_file_path and (
                chunk_workers > 1 or _use_dataset_output() or args.output_mode in {"edge_index", "outdegree"}
            ):
                n = n_workers
                out_parts_dir = Path(
                    tempfile.mkdtemp(
                        prefix=f"{output_file_path.stem}.grn_parts_",
                        dir=str(output_file_path.parent),
                    )
                )
                if args.output_mode == "outdegree":
                    if _env_flag("LOCAL_TE_OUTDEGREE_HIST", "on"):
                        tasks = [
                            (
                                [str(p) for p in cp] if isinstance(cp, list) else str(cp),
                                int(cidx),
                                int(chunk_size),
                                list(gene_names_list),
                                float(args.te_cutoff),
                                float(args.fdr),
                                list(labels),
                                str(out_parts_dir),
                                int(max_len),
                                int(batch_rows),
                                float(os.environ.get("LOCALTE_TRIM_THRESHOLD", "-0.01")),
                            )
                            for cidx, cp in chunk_files
                        ]
                        task_func = _process_chunk_outdegree_hist_task
                    else:
                        tasks = [
                            (
                                str(cp),
                                int(cidx),
                                int(chunk_size),
                                list(gene_names_list),
                                float(args.te_cutoff),
                                float(args.fdr),
                                list(labels),
                                str(out_parts_dir),
                                int(max_len),
                                int(batch_rows),
                                str(temp_root_path),
                                float(os.environ.get("LOCALTE_TRIM_THRESHOLD", "-0.01")),
                            )
                            for cidx, cp in chunk_files
                        ]
                        task_func = _process_chunk_outdegree_fast_task
                    part_pattern = "outdegree_chunk_*.parquet"
                elif edge_index_mode:
                    tasks = [
                        (
                            str(cp),
                            int(cidx),
                            int(chunk_size),
                            float(args.te_cutoff),
                            float(args.fdr),
                            list(labels),
                            str(out_parts_dir),
                            int(max_len),
                            int(batch_rows),
                        )
                        for cidx, cp in chunk_files
                    ]
                    task_func = _process_chunk_edge_index_hist_task
                    part_pattern = "grn_chunk_*.parquet"
                elif histogram_mode:
                    tasks = [
                        (
                            str(cp),
                            int(cidx),
                            int(chunk_size),
                            list(gene_names_list),
                            float(args.te_cutoff),
                            float(args.fdr),
                            list(labels),
                            str(out_parts_dir),
                            int(max_len),
                            int(batch_rows),
                        )
                        for cidx, cp in chunk_files
                    ]
                    task_func = _process_chunk_histogram_file_task
                    part_pattern = "grn_chunk_*.parquet"
                else:
                    tasks = [
                        (
                            str(cp),
                            int(cidx),
                            int(chunk_size),
                            list(gene_names_list),
                            float(args.te_cutoff),
                            float(args.fdr),
                            list(labels),
                            list(orig_idx_seq) if orig_idx_seq is not None else None,
                            str(out_parts_dir),
                            int(max_len),
                            int(batch_rows),
                            str(temp_root_path),
                        )
                        for cidx, cp in chunk_files
                    ]
                    task_func = _process_chunk_file_task
                    part_pattern = "grn_chunk_*.parquet"
                print(
                    f"[LocalTE GRN] Processing chunks in parallel: workers={n}, "
                    f"temp_parts={out_parts_dir}",
                    flush=True,
                )
                with multiprocessing.Pool(processes=n) as pool:
                    iterator = pool.imap_unordered(task_func, tasks, chunksize=1)
                    written_parts = 0
                    if _tqdm is not None:
                        pbar = _tqdm(total=len(tasks), desc="LocalTE GRN chunks", unit="chunk", dynamic_ncols=True)
                        try:
                            for result in iterator:
                                if result:
                                    written_parts += 1
                                pbar.update(1)
                                pbar.set_postfix_str(f"parts={written_parts}")
                        finally:
                            pbar.close()
                    else:
                        for done, result in enumerate(iterator, start=1):
                            if result:
                                written_parts += 1
                            if done == 1 or done % 5 == 0 or done == len(tasks):
                                print(
                                    f"[LocalTE GRN] chunks {done}/{len(tasks)} completed; "
                                    f"parts={written_parts}",
                                    flush=True,
                                )
                part_parqs = sorted(Path(out_parts_dir).glob(part_pattern))
                if not part_parqs:
                    print("No GRN parts written.", file=sys.stderr)
                    shutil.rmtree(out_parts_dir, ignore_errors=True)
                    return
                out_path = output_file_path
                out_path_str = str(out_path)
                _, ext = os.path.splitext(out_path_str)
                if _use_dataset_output():
                    dataset_dir = _dataset_output_dir(output_file_path)
                    if dataset_dir.exists() and any(dataset_dir.iterdir()) and not _env_flag("LOCAL_TE_OVERWRITE_DATASET", "off"):
                        raise FileExistsError(
                            f"Dataset output directory already exists and is not empty: {dataset_dir}. "
                            "Set LOCAL_TE_OUTPUT_DATASET_DIR to a new path or LOCAL_TE_OVERWRITE_DATASET=on."
                        )
                    dataset_dir.mkdir(parents=True, exist_ok=True)
                    total_rows = 0
                    for part in _progress_iter(part_parqs, len(part_parqs), "Finalize GRN dataset", "part"):
                        pf_part = pq.ParquetFile(part)
                        total_rows += int(pf_part.metadata.num_rows)
                        final_part = dataset_dir / part.name
                        if final_part.exists() and _env_flag("LOCAL_TE_OVERWRITE_DATASET", "off"):
                            final_part.unlink()
                        shutil.move(str(part), str(final_part))
                    _write_dataset_manifest(
                        dataset_dir,
                        output_mode=args.output_mode,
                        part_pattern=part_pattern,
                        part_count=len(part_parqs),
                        row_count=total_rows,
                        fdr=float(args.fdr),
                        te_cutoff=float(args.te_cutoff),
                        chunk_size=int(chunk_size),
                        max_timepoints=int(max_len),
                    )
                    try:
                        shutil.rmtree(out_parts_dir, ignore_errors=True)
                    except Exception:
                        pass
                    print(f"Saved GRN parquet dataset to {dataset_dir}")
                    return
                if ext.lower() == ".csv":
                    with out_path.open("w", encoding="utf-8", newline="") as handle:
                        header_written = False
                        for part in _progress_iter(part_parqs, len(part_parqs), "Merge GRN parts", "part"):
                            part_df = pd.read_parquet(part)
                            if args.output_mode == "outdegree":
                                part_out = part_df.loc[:, ["Time", "Source", "Outdegree"]]
                            else:
                                part_df["Time"] = part_df["TimeLabel"]
                                part_out = part_df.loc[:, ["Source", "TE", "Target", "Time"]]
                            part_out.to_csv(handle, index=False, header=not header_written)
                            header_written = True
                else:
                    writer = None
                    try:
                        for part in _progress_iter(part_parqs, len(part_parqs), "Merge GRN parts", "part"):
                            pf_part = pq.ParquetFile(part)
                            for batch in pf_part.iter_batches():
                                table = pa.Table.from_batches([batch])
                                if writer is None:
                                    writer = pq.ParquetWriter(out_path_str, table.schema, compression="snappy")
                                writer.write_table(table)
                    finally:
                        if writer is not None:
                            writer.close()
                if args.output_mode == "outdegree":
                    summary_parts = sorted(Path(out_parts_dir).glob("summary_chunk_*.parquet"))
                    if summary_parts:
                        summary_out = output_file_path.with_suffix("").with_name(output_file_path.stem + ".summary.parquet")
                        writer = None
                        try:
                            for part in _progress_iter(summary_parts, len(summary_parts), "Merge trim summaries", "part"):
                                pf_part = pq.ParquetFile(part)
                                for batch in pf_part.iter_batches():
                                    table = pa.Table.from_batches([batch])
                                    if writer is None:
                                        writer = pq.ParquetWriter(str(summary_out), table.schema, compression="snappy")
                                    writer.write_table(table)
                        finally:
                            if writer is not None:
                                writer.close()
                elif args.output_mode == "edge_index" and _env_flag("LOCAL_TE_EDGE_INDEX_WRITE_CATALOG", "on"):
                    catalog_out = output_file_path.with_suffix("").with_name(output_file_path.stem + ".edge_catalog.parquet")
                    _write_edge_index_catalog_from_chunk(
                        chunk_files[0][1],
                        catalog_out,
                        gene_names_list,
                        batch_rows,
                    )
                    print(f"[LocalTE GRN] Edge catalog saved to {catalog_out}", flush=True)
                try:
                    shutil.rmtree(out_parts_dir, ignore_errors=True)
                except Exception:
                    pass
                print(f"Saved aggregated network to {out_path}")
                return

            if histogram_mode and output_file_path:
                out_parts_dir = Path(
                    tempfile.mkdtemp(
                        prefix=f"{output_file_path.stem}.grn_parts_",
                        dir=str(output_file_path.parent),
                    )
                )
                part_paths: List[str] = []
                for cidx, cp in _progress_iter(chunk_files, len(chunk_files), "LocalTE GRN chunks", "chunk"):
                    result = _process_chunk_histogram_file_task(
                        (
                            str(cp),
                            int(cidx),
                            int(chunk_size),
                            list(gene_names_list),
                            float(args.te_cutoff),
                            float(args.fdr),
                            list(labels),
                            str(out_parts_dir),
                            int(max_len),
                            int(batch_rows),
                        )
                    )
                    if result:
                        part_paths.append(result)
                part_parqs = sorted(Path(out_parts_dir).glob("grn_chunk_*.parquet"))
                if not part_parqs:
                    print("No GRN parts written.", file=sys.stderr)
                    shutil.rmtree(out_parts_dir, ignore_errors=True)
                    return
                writer = None
                try:
                    for part in _progress_iter(part_parqs, len(part_parqs), "Merge GRN parts", "part"):
                        pf_part = pq.ParquetFile(part)
                        for batch in pf_part.iter_batches():
                            table = pa.Table.from_batches([batch])
                            if writer is None:
                                writer = pq.ParquetWriter(str(output_file_path), table.schema, compression="snappy")
                            writer.write_table(table)
                finally:
                    if writer is not None:
                        writer.close()
                    shutil.rmtree(out_parts_dir, ignore_errors=True)
                print(f"Saved aggregated network to {output_file_path}")
                return

            for chunk_idx, chunk_path in _progress_iter(chunk_files, len(chunk_files), "LocalTE GRN chunks", "chunk"):
                frames = _process_chunk_streaming(
                    chunk_idx=chunk_idx,
                    chunk_path=chunk_path,
                    chunk_size=chunk_size,
                    te_cutoff=float(args.te_cutoff),
                    fdr=float(args.fdr),
                    labels=labels,
                    orig_idx_seq=orig_idx_seq,
                    gene_names_arr=gene_names_arr,
                    max_len=max_len,
                    batch_rows=batch_rows,
                    temp_root=temp_root_path if temp_root_path is not None else Path(tempfile.gettempdir()),
                )
                for df_time in frames:
                    if df_time is not None and not df_time.empty:
                        tuples.append((int(df_time["TimeIndex"].iat[0]), df_time))
        finally:
            if temp_root_path is not None:
                shutil.rmtree(temp_root_path, ignore_errors=True)
            _cleanup_merged_parts()
    else:
        sources = df["Source"].to_numpy()
        targets = df["Target"].to_numpy()
        matrix = decode_local_matrix(df, max_len)
        time_indices = list(range(max_len))
        workers = max(0, int(args.workers or 0))
        if workers <= 1:
            results = [
                process_timepoint(t_idx, matrix, sources, targets, labels, args.te_cutoff, args.fdr, orig_idx_seq)
                for t_idx in time_indices
            ]
        else:
            workers = min(workers, multiprocessing.cpu_count())
            with multiprocessing.Pool(
                processes=workers,
                initializer=_worker_init,
                initargs=(matrix, sources, targets, labels, orig_idx_seq, args.te_cutoff, args.fdr),
            ) as pool:
                results = list(pool.imap_unordered(_worker_task, time_indices, chunksize=1))
        for df_time in results:
            if df_time is not None and not df_time.empty:
                time_idx = int(df_time["TimeIndex"].iat[0])
                tuples.append((time_idx, df_time))
    if not tuples:
        print("No edges passed the filters; nothing was written.")
        return

    tuples.sort(key=lambda x: x[0])

    combined = pd.concat([df_time for _, df_time in tuples], ignore_index=True)
    combined["Time"] = combined["TimeLabel"]
    combined_out = combined.loc[:, ["Source", "TE", "Target", "Time"]]

    if output_file_path:
        out_path = output_file_path
    else:
        base = args.output_prefix or "localte_times"
        if base.lower().endswith((".parquet", ".csv")):
            out_path = coerce_output_path(base)
        else:
            out_path = coerce_output_path(f"{base}.parquet")
        out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.suffix.lower() == ".csv":
        combined_out.to_csv(out_path, index=False)
    else:
        combined_out.to_parquet(out_path, index=False)
    print(f"Saved aggregated network ({len(combined_out):,} rows) to {out_path}")


if __name__ == "__main__":
    main()
