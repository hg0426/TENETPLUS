#!/usr/bin/env python3
"""
Build per-timepoint GRNs from compressed LocalTE outputs (float16 bytes + len/dtype).
Replicates the legacy make_GRN_new.py logic at each time index: filter TE>cutoff, compute
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

from code.path_utils import (
    coerce_input_path,
    coerce_output_path,
    ensure_output_subdir,
    resolve_output,
)


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
    for batch in pf.iter_batches(columns=["Values"], batch_size=batch_rows, use_threads=True):
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

    for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=batch_rows, use_threads=True):
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
    for local_idx in range(chunk_len):
        handles.append((chunk_tmp / f"time_{local_idx:03d}.bin").open("wb"))

    try:
        for batch in pf.iter_batches(columns=["Source", "Target", "Values"], batch_size=batch_rows, use_threads=True):
            values = batch.column("Values")
            flat = values.values.to_numpy(zero_copy_only=False).view(np.float16)
            if flat.size == 0:
                continue
            mat_full = np.asarray(flat, dtype=np.float32).reshape(-1, raw_chunk_len)
            mat = mat_full[:, :chunk_len]
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
        triples = np.memmap(file_path, dtype=_TRIPLE_DTYPE, mode="r")
        try:
            te_vals = np.asarray(triples["TE"], dtype=np.float32)
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
    parser.add_argument("--workers", type=int, default=0, help="Worker processes (0/1 sequential, >1 multiprocessing).")
    parser.add_argument("--chunk_dir", help="Directory containing chunked LocalTE parquet files (from localte_chunk_export.py).")
    parser.add_argument("--chunk_size", type=int, help="Timepoints per chunk when using --chunk_dir (default: metadata.json value).")
    parser.add_argument("--chunk_workers", type=int, default=0, help="Parallel workers for chunk_dir processing (0/1=single-process; >1 parallel by chunk). Only used when --output_file is set.")
    args = parser.parse_args()

    output_file_path: Optional[Path] = None
    if args.output_file:
        output_file_path = coerce_output_path(args.output_file)
        output_file_path.parent.mkdir(parents=True, exist_ok=True)

    def load_gene_index_map() -> Dict[str, int]:
        gene_file = coerce_input_path("gene_names")
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
        if not meta_path.exists():
            print(f"Expected metadata.json in {chunk_dir}", file=sys.stderr)
            sys.exit(1)
        with meta_path.open() as handle:
            chunk_meta = json.load(handle)
        if args.chunk_size is None:
            chunk_size = int(chunk_meta.get("chunk_size", chunk_size))
        max_len = max(max_len, int(chunk_meta.get("max_timepoints", 0)))
        if max_len == 0:
            num_chunks_meta = int(chunk_meta.get("num_chunks", 0))
            max_len = chunk_size * num_chunks_meta
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
        for chunk_idx in range(num_chunks_eval):
            cp = chunk_dir / f"chunk_{chunk_idx:04d}.parquet"
            if cp.exists():
                chunk_files.append((chunk_idx, cp))
        if not chunk_files:
            print(f"No chunk parquet files found in {chunk_dir}", file=sys.stderr)
            return

        gene_names_list = _load_gene_names_list(str(coerce_input_path("gene_names")))
        gene_names_arr = np.asarray(gene_names_list, dtype=object)
        batch_rows = _chunk_batch_rows()
        # Use system temp so no grn_stream_* folders appear under project output
        temp_root_path = Path(tempfile.mkdtemp(prefix="grn_stream_"))

        try:
            chunk_workers = int(getattr(args, "chunk_workers", 0) or 0)
            if output_file_path and chunk_workers > 1:
                out_parts_dir = ensure_output_subdir("grn_parts")
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
                n = min(chunk_workers, multiprocessing.cpu_count(), len(tasks))
                if n <= 0:
                    n = 1
                with multiprocessing.Pool(processes=n) as pool:
                    for _ in pool.imap_unordered(_process_chunk_file_task, tasks, chunksize=1):
                        pass
                part_parqs = sorted(Path(out_parts_dir).glob("grn_chunk_*.parquet"))
                if not part_parqs:
                    print("No GRN parts written.", file=sys.stderr)
                    return
                out_path = output_file_path
                out_path_str = str(out_path)
                _, ext = os.path.splitext(out_path_str)
                if ext.lower() == ".csv":
                    with out_path.open("w", encoding="utf-8", newline="") as handle:
                        header_written = False
                        for part in part_parqs:
                            part_df = pd.read_parquet(part)
                            part_df["Time"] = part_df["TimeLabel"]
                            part_out = part_df.loc[:, ["Source", "TE", "Target", "Time"]]
                            part_out.to_csv(handle, index=False, header=not header_written)
                            header_written = True
                    # remove intermediate parts directory after successful merge
                    try:
                        shutil.rmtree(out_parts_dir, ignore_errors=True)
                    except Exception:
                        pass
                else:
                    writer = None
                    try:
                        for part in part_parqs:
                            pf_part = pq.ParquetFile(part)
                            for batch in pf_part.iter_batches():
                                table = pa.Table.from_batches([batch])
                                if writer is None:
                                    writer = pq.ParquetWriter(out_path_str, table.schema, compression="snappy")
                                writer.write_table(table)
                    finally:
                        if writer is not None:
                            writer.close()
                        # remove intermediate parts directory regardless of minor cleanup errors
                        try:
                            shutil.rmtree(out_parts_dir, ignore_errors=True)
                        except Exception:
                            pass
                print(f"Saved aggregated network to {out_path}")
                return

            for chunk_idx, chunk_path in chunk_files:
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
                    temp_root=temp_root_path,
                )
                for df_time in frames:
                    if df_time is not None and not df_time.empty:
                        tuples.append((int(df_time["TimeIndex"].iat[0]), df_time))
        finally:
            shutil.rmtree(temp_root_path, ignore_errors=True)
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
