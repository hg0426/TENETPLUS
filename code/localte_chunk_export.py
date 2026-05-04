#!/usr/bin/env python3
"""
Convert LocalTE byte payloads into time-chunked parquet files.

Each output chunk stores per-edge LocalTE segments (default window size 300)
so downstream analyses can load a limited set of timesteps at once instead of
materialising the full edge×time matrix in memory.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import queue as queue_mod
import time
import zlib
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import multiprocessing as mp

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm
import re

from code.path_utils import coerce_input_path, locate_file, resolve_output


def _localte_gene_names_file() -> Path:
    out_path = resolve_output("gene_names")
    if out_path.exists():
        return out_path
    return locate_file("gene_names")


def _set_stage(label: str) -> None:
    stage_file = os.environ.get("TENET_STAGE_FILE")
    if not stage_file:
        return
    try:
        Path(stage_file).write_text(f"{label}\n")
    except Exception:
        pass


def _as_text(value, default: str = "") -> str:
    if value is None:
        return default
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def _decode_local_te(blob: bytes, dtype_str, length: int, codec) -> np.ndarray:
    """Decode one LocalTE payload into its stored numeric array."""
    dtype = np.dtype(_as_text(dtype_str, "float16"))
    codec_str = _as_text(codec, "none").strip().lower()
    raw = zlib.decompress(blob) if codec_str == "zlib" else blob
    arr = np.frombuffer(raw, dtype=dtype, count=int(length))
    if arr.size != int(length):
        raise ValueError(
            f"Decoded LocalTE length mismatch: expected {int(length)}, "
            f"got {arr.size} (codec={codec_str or 'none'}, dtype={dtype})"
        )
    return arr


def _merge_one_task(args: Tuple[int, List[Path], Path]) -> int:
    """Top-level merge worker to combine per-worker part files for one chunk index.
    Returns the chunk index on completion (for progress accounting).
    """
    cidx_m, files_m, out_dir_m = args
    final_path_m = out_dir_m / f"chunk_{cidx_m:04d}.parquet"
    writer_m = None
    for f_m in files_m:
        pf_part = pq.ParquetFile(f_m)
        for b in pf_part.iter_batches():
            if writer_m is None:
                writer_m = pq.ParquetWriter(final_path_m, b.schema, compression="snappy")
            writer_m.write_table(pa.Table.from_batches([b]))
    if writer_m is not None:
        writer_m.close()
    return cidx_m


def _merge_one_labeled_task(args: Tuple[str, int, List[Path], Path]) -> Tuple[str, int]:
    key, cidx, files, out_dir = args
    _merge_one_task((cidx, files, out_dir))
    return key, cidx


def _sanitize_subdir_name(path: Path) -> str:
    stem = path.stem or "dataset"
    cleaned = []
    for ch in stem:
        if ch.isalnum() or ch in ("-", "_"):
            cleaned.append(ch)
        else:
            cleaned.append("_")
    return "".join(cleaned)


def _load_gene_names_map() -> Dict[str, int]:
    gene_file = _localte_gene_names_file()
    if not gene_file.exists():
        raise FileNotFoundError("gene_names file not found; required to map gene symbols to indices.")
    mapping: Dict[str, int] = {}
    with gene_file.open(encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            mapping[line.strip()] = idx
    return mapping


def _load_gene_names_list() -> List[str]:
    gene_file = _localte_gene_names_file()
    if not gene_file.exists():
        raise FileNotFoundError("gene_names file not found; required for split LocalTE export.")
    with gene_file.open(encoding="utf-8") as handle:
        return [line.strip() for line in handle]


def _coerce_index(value, mapping: Optional[Dict[str, int]]) -> int:
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, float) and value.is_integer():
        return int(value)
    string = str(value).strip()
    try:
        return int(string)
    except ValueError:
        if mapping is None:
            raise
        if string not in mapping:
            raise ValueError(f"Gene/peak name '{string}' not found in gene_names mapping.")
        return mapping[string]


def _read_selector(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path, columns=["Source", "Target"])
    if suffix in {".csv"}:
        return pd.read_csv(path, usecols=["Source", "Target"])
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", usecols=["Source", "Target"])
    raise ValueError(f"Unsupported selector format for {path}; expected parquet/csv/tsv.")


def _read_unique_source_indices(path: Path, mapping: Dict[str, int]) -> set[int]:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        con = None
        try:
            import duckdb
            con = duckdb.connect()
            path_str = str(path).replace("'", "''")
            values = con.execute(
                f"SELECT DISTINCT Source FROM read_parquet('{path_str}')"
            ).fetch_df()["Source"].tolist()
        finally:
            if con is not None:
                con.close()
    elif suffix == ".csv":
        values = pd.read_csv(path, usecols=["Source"])["Source"].drop_duplicates().tolist()
    elif suffix in {".tsv", ".txt"}:
        values = pd.read_csv(path, sep="\t", usecols=["Source"])["Source"].drop_duplicates().tolist()
    else:
        raise ValueError(f"Unsupported selector format for {path}; expected parquet/csv/tsv.")
    return {_coerce_index(v, mapping) for v in values}


def _build_allowed_map(selector_paths: Sequence[Path]) -> Dict[int, set]:
    if not selector_paths:
        return {}
    mapping: Optional[Dict[str, int]] = None
    allowed: Dict[int, set] = defaultdict(set)
    total = 0
    for sel_path in selector_paths:
        df = _read_selector(sel_path)
        if df.empty:
            continue
        for row in df.itertuples(index=False):
            try:
                src = _coerce_index(row.Source, mapping)
            except ValueError:
                if mapping is None:
                    mapping = _load_gene_names_map()
                    src = _coerce_index(row.Source, mapping)
                else:
                    raise
            try:
                tgt = _coerce_index(row.Target, mapping)
            except ValueError:
                if mapping is None:
                    mapping = _load_gene_names_map()
                    tgt = _coerce_index(row.Target, mapping)
                else:
                    raise
            allowed[src].add(tgt)
            total += 1
    if total == 0:
        print(f"[selector] Warning: combined selector contained 0 edges.")
    else:
        print(f"[selector] Loaded {total:,} edge constraints across {len(selector_paths)} file(s).")
    return allowed


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Chunk LocalTE payloads into time windows.")
    parser.add_argument("--input", required=True, nargs="+", help="One or more Parquet files containing LocalTE_bytes/len/dtype columns (or views filtered from them).")
    parser.add_argument("--output_dir", required=True, help="Directory to place chunk parquet files (when multiple inputs are provided, subdirectories are created automatically).")
    parser.add_argument("--chunk_size", type=int, default=300, help="Number of timesteps per chunk window (default 300).")
    parser.add_argument("--buffer_edges", type=int, default=2000, help="Edges per chunk to buffer before flushing (default 2000). Increase (e.g. 20000) for fewer writes and higher throughput if RAM allows.")
    parser.add_argument("--read_batch_rows", type=int, default=8192, help="Parquet batch size when scanning input (default 8192 rows). Larger is faster if RAM allows.")
    parser.add_argument("--use_threads", choices=("on","off"), default="on", help="Let Arrow use multi-threaded decode (default on).")
    parser.add_argument("--values_dtype", choices=("float16","float32"), default="float16", help="Store chunked Values as this dtype (float16 is smaller/faster; default float16).")
    parser.add_argument("--scores_output", help="Optional path to write TE scores without LocalTE payload columns.")
    parser.add_argument("--workers", type=int, default=0, help="Worker processes for parallel export (0/1 = single-process; >1 enables multiprocessing by row-group).")
    parser.add_argument("--merge_workers", type=int, default=None, help="Worker processes for merging chunk parts (default: same as --workers; 0/1 = single-process).")
    parser.add_argument(
        "--merge_parts",
        choices=("on", "off"),
        default=os.getenv("LOCAL_TE_MERGE_PARTS", "off").lower(),
        help="Merge worker part files into one chunk_XXXX.parquet file per time chunk. "
             "off is faster/lower-memory and stores an equivalent partitioned dataset (default: off).",
    )
    parser.add_argument("--compression", choices=("snappy","zstd","gzip","none"), default="snappy", help="Parquet compression codec for outputs (default snappy).")
    parser.add_argument("--zstd_level", type=int, default=None, help="ZSTD compression level when --compression zstd (optional).")
    parser.add_argument("--no_dictionary", action="store_true", help="Disable Parquet dictionary encoding (default: enabled).")
    parser.add_argument("--disable_scores", action="store_true", help="Do not write TE_result_scores.parquet (skip score output).")
    return parser


def init_chunk_buffers() -> Dict[int, Dict[str, List]]:
    return defaultdict(
        lambda: {
            "sources": [],
            "targets": [],
            "values_arrays": [],
            "offsets": [0],
        }
    )


def flush_chunk(chunk_idx: int, buffer: Dict[str, List], writers: Dict[int, pq.ParquetWriter], out_dir: Path, chunk_size: int, value_type: pa.DataType, compression: str = "snappy", zstd_level: int | None = None, use_dictionary: bool = True) -> None:
    if not buffer["sources"]:
        return
    num_rows = len(buffer["sources"])
    offsets = buffer["offsets"]
    if len(offsets) != num_rows + 1:
        raise ValueError("Offsets length mismatch when flushing chunk.")
    if buffer["values_arrays"]:
        # choose dtype according to value_type
        np_dtype = np.float16 if pa.types.is_float16(value_type) else np.float32
        flat = np.concatenate(buffer["values_arrays"]).astype(np_dtype, copy=False)
    else:
        np_dtype = np.float16 if pa.types.is_float16(value_type) else np.float32
        flat = np.empty(0, dtype=np_dtype)
    values_array = pa.array(flat, type=value_type)
    offsets_array = pa.array(offsets, type=pa.int32())
    list_array = pa.ListArray.from_arrays(offsets_array, values_array)

    table = pa.table(
        {
            "Source": pa.array(buffer["sources"]),
            "Target": pa.array(buffer["targets"]),
            "ChunkIndex": pa.array([chunk_idx] * num_rows, type=pa.int32()),
            "ChunkStart": pa.array([chunk_idx * chunk_size] * num_rows, type=pa.int32()),
            "Values": list_array,
        }
    )
    writer = writers.get(chunk_idx)
    if writer is None:
        out_path = out_dir / f"chunk_{chunk_idx:04d}.parquet"
        # Try to honor compression options if available
        try:
            writers[chunk_idx] = pq.ParquetWriter(out_path, table.schema, compression=compression, use_dictionary=use_dictionary, compression_level=zstd_level)
        except TypeError:
            writers[chunk_idx] = pq.ParquetWriter(out_path, table.schema, compression=compression, use_dictionary=use_dictionary)
        writer = writers[chunk_idx]
    writer.write_table(table)
    buffer["sources"].clear()
    buffer["targets"].clear()
    buffer["values_arrays"].clear()
    buffer["offsets"].clear()
    buffer["offsets"].append(0)


def close_writers(writers: Dict[int, pq.ParquetWriter]) -> None:
    for writer in writers.values():
        writer.close()


def _export_dataset(args, selector_paths: Sequence[Path]) -> None:
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input parquet not found: {input_path}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pf = pq.ParquetFile(input_path)
    required_cols = ["Source", "Target", "TE", "LocalTE_bytes", "LocalTE_len", "LocalTE_dtype"]
    schema_names = set(pf.schema.names)
    missing = set(required_cols) - schema_names
    if missing:
        raise ValueError(f"Input parquet missing required columns: {sorted(missing)}")
    if "LocalTE_codec" in schema_names:
        required_cols.append("LocalTE_codec")

    selector_paths = list(selector_paths)
    allowed_map: Dict[int, set] = {}
    if selector_paths:
        raw_map = _build_allowed_map(selector_paths)
        allowed_map = {src: set(targets) for src, targets in raw_map.items()}
        print(f"[selector] Filtering export to {sum(len(v) for v in allowed_map.values()):,} edges across {len(allowed_map)} sources.")

    chunk_size = max(1, int(args.chunk_size))
    buffer_edges = max(1, int(args.buffer_edges))
    read_batch = max(1, int(args.read_batch_rows))
    use_threads = (str(args.use_threads).lower() == "on")
    values_dtype = np.float16 if args.values_dtype == "float16" else np.float32
    value_type = pa.float16() if values_dtype is np.float16 else pa.float32()

    chunk_buffers = init_chunk_buffers()
    writers: Dict[int, pq.ParquetWriter] = {}

    # Stream scores to disk instead of holding all in lists (optional)
    scores_writer: pq.ParquetWriter | None = None
    max_length = 0
    min_length = math.inf
    total_rows_seen = 0

    total_rows = pf.metadata.num_rows
    # Generic fallback path. Optimized selector splits are handled by _export_known_splits_single_pass.
    # Single-process path
    progress = tqdm(total=total_rows, desc="Chunking LocalTE", unit="row")

    for batch in pf.iter_batches(columns=list(required_cols), batch_size=read_batch, use_threads=use_threads):
        batch_sources = batch.column("Source").to_numpy(zero_copy_only=False)
        batch_targets = batch.column("Target").to_numpy(zero_copy_only=False)
        batch_scores = batch.column("TE").to_numpy(zero_copy_only=False)
        batch_bytes = batch.column("LocalTE_bytes").to_pylist()
        batch_lens_np = batch.column("LocalTE_len").to_numpy(zero_copy_only=False)
        batch_dtypes = batch.column("LocalTE_dtype").to_pylist()
        batch_codecs = None
        try:
            if "LocalTE_codec" in batch.schema.names:
                batch_codecs = batch.column("LocalTE_codec").to_pylist()
        except Exception:
            batch_codecs = None

        batch_lens = [int(v) if v is not None else 0 for v in batch_lens_np]

        filt_sources: List[int] = []
        filt_targets: List[int] = []
        filt_scores: List[float] = []
        filt_lens: List[int] = []

        for idx, (s, t, score, blob, length, dtype_str) in enumerate(zip(
            batch_sources, batch_targets, batch_scores, batch_bytes, batch_lens, batch_dtypes
        )):
            progress.update(1)
            if blob is None or length is None or dtype_str is None or length == 0:
                continue
            src_i = int(s)
            tgt_i = int(t)
            if allowed_map and (src_i not in allowed_map or tgt_i not in allowed_map[src_i]):
                continue
            codec = batch_codecs[idx] if batch_codecs is not None else None
            arr = _decode_local_te(blob, dtype_str, int(length), codec)
            if arr.size == 0:
                continue
            values = arr.astype(values_dtype, copy=False) if arr.dtype != values_dtype else arr
            max_length = max(max_length, int(length))
            min_length = min(min_length, int(length))
            filt_sources.append(src_i)
            filt_targets.append(tgt_i)
            filt_scores.append(float(score))
            filt_lens.append(int(length))
            num_chunks = math.ceil(length / chunk_size)
            for chunk_idx in range(num_chunks):
                start = chunk_idx * chunk_size
                end = min(start + chunk_size, length)
                if start >= end:
                    continue
                segment = values[start:end]
                buf = chunk_buffers[chunk_idx]
                buf["sources"].append(src_i)
                buf["targets"].append(tgt_i)
                buf["values_arrays"].append(segment)
                buf["offsets"].append(buf["offsets"][-1] + segment.size)
                if len(buf["sources"]) >= buffer_edges:
                    flush_chunk(chunk_idx, buf, writers, out_dir, chunk_size, value_type)

        if not filt_sources:
            continue

        if not args.disable_scores and scores_writer is None:
            scores_schema = pa.schema(
                [
                    ("Source", pa.int32()),
                    ("Target", pa.int32()),
                    ("TE", pa.float32()),
                    ("LocalTE_total_len", pa.int32()),
                ]
            )
            scores_output = Path(args.scores_output) if args.scores_output else out_dir / "TE_result_scores.parquet"
            try:
                scores_writer = pq.ParquetWriter(scores_output, scores_schema, compression=args.compression, use_dictionary=not args.no_dictionary, compression_level=args.zstd_level)
            except TypeError:
                scores_writer = pq.ParquetWriter(scores_output, scores_schema, compression=args.compression, use_dictionary=not args.no_dictionary)

        if not args.disable_scores:
            scores_tbl = pa.table(
                {
                    "Source": pa.array(filt_sources, type=pa.int32()),
                    "Target": pa.array(filt_targets, type=pa.int32()),
                    "TE": pa.array(np.asarray(filt_scores, dtype=np.float32), type=pa.float32()),
                    "LocalTE_total_len": pa.array(filt_lens, type=pa.int32()),
                }
            )
            scores_writer.write_table(scores_tbl)
        total_rows_seen += len(filt_sources)

    progress.close()
    progress.close()

    for chunk_idx, buf in chunk_buffers.items():
        flush_chunk(chunk_idx, buf, writers, out_dir, chunk_size, value_type, args.compression, args.zstd_level, not args.no_dictionary)
    close_writers(writers)

    # Finish scores writer
    if scores_writer is not None:
        scores_writer.close()
    if total_rows_seen == 0:
        max_length = int(max_length)
        min_length = 0
    num_chunks_total = math.ceil(max_length / chunk_size) if max_length else 0

    metadata = {
        "input": str(input_path),
        "scores_output": (str(Path(args.scores_output) if args.scores_output else (out_dir / "TE_result_scores.parquet")) if not args.disable_scores else ""),
        "chunk_size": chunk_size,
        "total_edges": int(total_rows),
        "max_timepoints": max_length,
        "min_timepoints": min_length,
        "num_chunks": num_chunks_total,
        "buffer_edges": buffer_edges,
        "values_dtype": args.values_dtype,
        "use_threads": use_threads,
        "compression": args.compression,
        "zstd_level": args.zstd_level,
        "dictionary": (not args.no_dictionary),
    }
    with (out_dir / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)

    summary_scores = metadata["scores_output"] or "<disabled>"
    print(f"Wrote chunked LocalTE to {out_dir} (chunks={num_chunks_total}), scores table -> {summary_scores}")


def _new_export_state(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    return {
        "out_dir": out_dir,
        "chunk_buffers": init_chunk_buffers(),
        "writers": {},
        "scores_writer": None,
        "total": 0,
        "max_length": 0,
        "min_length": math.inf,
        "score_sources": [],
        "score_targets": [],
        "score_values": [],
        "score_lens": [],
    }


def _flush_scores_state(
    state: dict,
    scores_schema: pa.Schema,
    compression: str,
    use_dictionary: bool,
) -> None:
    if not state["score_sources"]:
        return
    if state["scores_writer"] is None:
        state["scores_writer"] = pq.ParquetWriter(
            state["out_dir"] / "TE_result_scores.parquet",
            scores_schema,
            compression=compression,
            use_dictionary=use_dictionary,
        )
    scores_tbl = pa.table(
        {
            "Source": pa.array(state["score_sources"], type=pa.int32()),
            "Target": pa.array(state["score_targets"], type=pa.int32()),
            "TE": pa.array(np.asarray(state["score_values"], dtype=np.float32), type=pa.float32()),
            "LocalTE_total_len": pa.array(state["score_lens"], type=pa.int32()),
        }
    )
    state["scores_writer"].write_table(scores_tbl)
    state["score_sources"].clear()
    state["score_targets"].clear()
    state["score_values"].clear()
    state["score_lens"].clear()


def _append_localte_to_state(
    state: dict,
    src_i: int,
    tgt_i: int,
    score: float,
    values: np.ndarray,
    length: int,
    chunk_size: int,
    buffer_edges: int,
    value_type: pa.DataType,
    compression: str,
    zstd_level: int | None,
    use_dictionary: bool,
    disable_scores: bool,
    scores_schema: pa.Schema,
) -> None:
    state["max_length"] = max(int(state["max_length"]), int(length))
    state["min_length"] = min(state["min_length"], int(length))
    state["total"] = int(state["total"]) + 1

    if not disable_scores:
        state["score_sources"].append(src_i)
        state["score_targets"].append(tgt_i)
        state["score_values"].append(float(score))
        state["score_lens"].append(int(length))
        if len(state["score_sources"]) >= buffer_edges:
            _flush_scores_state(state, scores_schema, compression, use_dictionary)

    n_chunks = math.ceil(length / chunk_size)
    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, length)
        if start >= end:
            continue
        segment = values[start:end]
        buf = state["chunk_buffers"][chunk_idx]
        buf["sources"].append(src_i)
        buf["targets"].append(tgt_i)
        buf["values_arrays"].append(segment)
        buf["offsets"].append(buf["offsets"][-1] + int(segment.size))
        if len(buf["sources"]) >= buffer_edges:
            flush_chunk(
                chunk_idx,
                buf,
                state["writers"],
                state["out_dir"],
                chunk_size,
                value_type,
                compression,
                zstd_level,
                use_dictionary,
            )


def _finalize_export_state(state: dict, args, chunk_size: int, value_type: pa.DataType) -> None:
    scores_schema = pa.schema(
        [
            ("Source", pa.int32()),
            ("Target", pa.int32()),
            ("TE", pa.float32()),
            ("LocalTE_total_len", pa.int32()),
        ]
    )
    if not args.disable_scores:
        _flush_scores_state(state, scores_schema, args.compression, not args.no_dictionary)
    for chunk_idx, buf in state["chunk_buffers"].items():
        if buf["sources"]:
            flush_chunk(
                chunk_idx,
                buf,
                state["writers"],
                state["out_dir"],
                chunk_size,
                value_type,
                args.compression,
                args.zstd_level,
                not args.no_dictionary,
            )
    close_writers(state["writers"])
    if state["scores_writer"] is not None:
        state["scores_writer"].close()
    max_length = int(state["max_length"])
    min_length = 0 if state["min_length"] == math.inf else int(state["min_length"])
    metadata = {
        "input": "",
        "scores_output": "" if args.disable_scores else str(state["out_dir"] / "TE_result_scores.parquet"),
        "chunk_size": int(chunk_size),
        "total_edges": int(state["total"]),
        "max_timepoints": max_length,
        "min_timepoints": min_length,
        "num_chunks": math.ceil(max_length / chunk_size) if max_length else 0,
        "buffer_edges": int(args.buffer_edges),
        "values_dtype": args.values_dtype,
        "use_threads": str(args.use_threads).lower() == "on",
        "compression": args.compression,
        "zstd_level": args.zstd_level,
        "dictionary": (not args.no_dictionary),
    }
    with (state["out_dir"] / "metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)


def _known_split_key(selector_paths: Sequence[Path]) -> str:
    if not selector_paths:
        return ""
    return selector_paths[0].stem


def _can_single_pass_known_splits(input_specs: Sequence[Tuple[Path, List[Path]]]) -> bool:
    if len(input_specs) <= 1:
        return False
    payloads = {str(p.resolve()) for p, _ in input_specs}
    if len(payloads) != 1:
        return False
    allowed = {"TE_TF_GN", "TE_TF_PK", "TE_PK_GN", "TE_PK_PK", "TE_GN_GN", "TE_all_features"}
    return all(paths and _known_split_key(paths) in allowed for _, paths in input_specs)


def _route_known_split_indices(src_i: int, tgt_i: int, is_peak: np.ndarray, tf_sources: set[int], keys: Sequence[str]) -> List[int]:
    src_peak = bool(src_i < is_peak.size and is_peak[src_i])
    tgt_peak = bool(tgt_i < is_peak.size and is_peak[tgt_i])
    out = []
    for out_idx, key in enumerate(keys):
        if key == "TE_TF_GN" and src_i in tf_sources and (not tgt_peak):
            out.append(out_idx)
        elif key == "TE_TF_PK" and src_i in tf_sources and tgt_peak:
            out.append(out_idx)
        elif key == "TE_PK_GN" and src_peak and (not tgt_peak):
            out.append(out_idx)
        elif key == "TE_PK_PK" and src_peak and tgt_peak:
            out.append(out_idx)
        elif key == "TE_GN_GN" and (not src_peak) and (not tgt_peak):
            out.append(out_idx)
        elif key == "TE_all_features":
            out.append(out_idx)
    return out


def _merge_scores_parts(parts: Sequence[Path], output_path: Path, compression: str, use_dictionary: bool) -> None:
    writer = None
    try:
        for part in parts:
            if not part.exists():
                continue
            part_pf = pq.ParquetFile(part)
            for batch in part_pf.iter_batches():
                table = pa.Table.from_batches([batch])
                if writer is None:
                    writer = pq.ParquetWriter(output_path, table.schema, compression=compression, use_dictionary=use_dictionary)
                writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()


def _worker_export_known_splits(args: Tuple) -> Tuple[int, List[int], List[int], List[int]]:
    (
        worker_id,
        input_path,
        base_output,
        keys,
        row_groups,
        is_peak,
        tf_sources,
        chunk_size,
        buffer_edges,
        read_batch,
        use_threads,
        values_dtype_str,
        compression,
        zstd_level,
        use_dictionary,
        disable_scores,
        progress_queue,
    ) = args
    input_path = Path(input_path)
    base_output = Path(base_output)
    keys = list(keys)
    is_peak = np.asarray(is_peak, dtype=bool)
    tf_sources = set(int(x) for x in tf_sources)
    values_dtype = np.float16 if values_dtype_str == "float16" else np.float32
    value_type = pa.float16() if values_dtype is np.float16 else pa.float32()
    scores_schema = pa.schema(
        [
            ("Source", pa.int32()),
            ("Target", pa.int32()),
            ("TE", pa.float32()),
            ("LocalTE_total_len", pa.int32()),
        ]
    )
    states = []
    for key in keys:
        states.append(_new_export_state(base_output / key / f"parts_w{worker_id:02d}"))

    pf = pq.ParquetFile(input_path)
    required_cols = ["Source", "Target", "TE", "LocalTE_bytes", "LocalTE_len", "LocalTE_dtype"]
    if "LocalTE_codec" in set(pf.schema.names):
        required_cols.append("LocalTE_codec")

    routed_counts = [0 for _ in keys]
    max_lengths = [0 for _ in keys]
    min_lengths = [math.inf for _ in keys]
    for batch in pf.iter_batches(columns=required_cols, batch_size=read_batch, use_threads=use_threads, row_groups=row_groups):
        batch_sources = batch.column("Source").to_numpy(zero_copy_only=False)
        batch_targets = batch.column("Target").to_numpy(zero_copy_only=False)
        batch_scores = batch.column("TE").to_numpy(zero_copy_only=False)
        batch_bytes = batch.column("LocalTE_bytes").to_pylist()
        batch_lens_np = batch.column("LocalTE_len").to_numpy(zero_copy_only=False)
        batch_dtypes = batch.column("LocalTE_dtype").to_pylist()
        batch_codecs = batch.column("LocalTE_codec").to_pylist() if "LocalTE_codec" in batch.schema.names else None
        for idx, (s, t, score, blob, length, dtype_str) in enumerate(zip(
            batch_sources, batch_targets, batch_scores, batch_bytes, batch_lens_np, batch_dtypes
        )):
            if blob is None or length is None or dtype_str is None or int(length) == 0:
                continue
            src_i = int(s)
            tgt_i = int(t)
            out_indices = _route_known_split_indices(src_i, tgt_i, is_peak, tf_sources, keys)
            if not out_indices:
                continue
            codec = batch_codecs[idx] if batch_codecs is not None else None
            arr = _decode_local_te(blob, dtype_str, int(length), codec)
            if arr.size == 0:
                continue
            values = arr.astype(values_dtype, copy=False) if arr.dtype != values_dtype else arr
            for out_idx in out_indices:
                routed_counts[out_idx] += 1
                max_lengths[out_idx] = max(max_lengths[out_idx], int(length))
                min_lengths[out_idx] = min(min_lengths[out_idx], int(length))
                _append_localte_to_state(
                    states[out_idx],
                    src_i,
                    tgt_i,
                    float(score),
                    values,
                    int(length),
                    chunk_size,
                    buffer_edges,
                    value_type,
                    compression,
                    zstd_level,
                    use_dictionary,
                    disable_scores,
                    scores_schema,
                )
        if progress_queue is not None:
            try:
                progress_queue.put(int(batch.num_rows))
            except Exception:
                pass

    for state in states:
        _finalize_export_state(
            state,
            argparse.Namespace(
                **{
                    "disable_scores": disable_scores,
                    "compression": compression,
                    "zstd_level": zstd_level,
                    "no_dictionary": (not use_dictionary),
                    "buffer_edges": buffer_edges,
                    "values_dtype": values_dtype_str,
                    "use_threads": "on" if use_threads else "off",
                }
            ),
            chunk_size,
            value_type,
        )
    min_lengths = [0 if x == math.inf else int(x) for x in min_lengths]
    return worker_id, routed_counts, max_lengths, min_lengths


def _export_known_splits_parallel(
    args,
    input_path: Path,
    specs: Sequence[Tuple[str, Path]],
    is_peak: np.ndarray,
    tf_sources: set[int],
) -> bool:
    n_workers = max(1, int(getattr(args, "workers", 0) or 0))
    if n_workers <= 1:
        return False
    pf = pq.ParquetFile(input_path)
    n_rg = int(pf.metadata.num_row_groups or 0)
    if n_rg <= 1:
        return False

    base_output = Path(args.output_dir)
    keys = [key for key, _ in specs]
    assignments = [[] for _ in range(n_workers)]
    for rg_idx in range(n_rg):
        assignments[rg_idx % n_workers].append(rg_idx)
    active = [(idx, row_groups) for idx, row_groups in enumerate(assignments) if row_groups]
    print(
        f"[LocalTE] Parallel single-pass export: {len(keys)} outputs, "
        f"{pf.metadata.num_rows:,} TE rows, {n_rg:,} row groups, {len(active)} workers."
    )
    _set_stage("LOCALTE_EXPORT_SPLIT")

    worker_args = [
        (
            wid,
            str(input_path),
            str(base_output),
            keys,
            row_groups,
            is_peak,
            sorted(tf_sources),
            max(1, int(args.chunk_size)),
            max(1, int(args.buffer_edges)),
            max(1, int(args.read_batch_rows)),
            str(args.use_threads).lower() == "on",
            args.values_dtype,
            args.compression,
            args.zstd_level,
            not args.no_dictionary,
            args.disable_scores,
            None,
        )
        for wid, row_groups in active
    ]

    total_counts = [0 for _ in keys]
    max_lengths = [0 for _ in keys]
    min_lengths = [math.inf for _ in keys]
    manager = mp.Manager()
    progress_queue = manager.Queue()
    worker_args = [tuple(list(wargs[:-1]) + [progress_queue]) for wargs in worker_args]
    try:
        with mp.Pool(processes=len(active)) as pool:
            async_results = [pool.apply_async(_worker_export_known_splits, (wargs,)) for wargs in worker_args]
            remaining = set(range(len(async_results)))
            with tqdm(total=int(pf.metadata.num_rows), desc="LocalTE rows", unit="row") as row_bar, tqdm(total=len(active), desc="LocalTE workers", unit="worker", leave=False) as worker_bar:
                while remaining:
                    progressed = False
                    while True:
                        try:
                            n_rows = progress_queue.get_nowait()
                        except queue_mod.Empty:
                            break
                        except Exception:
                            break
                        row_bar.update(int(n_rows))
                        progressed = True

                    done = []
                    for ridx in list(remaining):
                        if async_results[ridx].ready():
                            done.append(ridx)
                    for ridx in done:
                        _, counts, worker_max, worker_min = async_results[ridx].get()
                        for i in range(len(keys)):
                            total_counts[i] += int(counts[i])
                            max_lengths[i] = max(max_lengths[i], int(worker_max[i]))
                            if int(worker_min[i]) > 0:
                                min_lengths[i] = min(min_lengths[i], int(worker_min[i]))
                        worker_bar.update(1)
                        remaining.remove(ridx)
                        progressed = True

                    if not progressed:
                        time.sleep(0.2)

                while True:
                    try:
                        row_bar.update(int(progress_queue.get_nowait()))
                    except queue_mod.Empty:
                        break
                    except Exception:
                        break
    finally:
        manager.shutdown()

    merge_workers = int(args.merge_workers) if getattr(args, "merge_workers", None) is not None else min(len(active), max(1, mp.cpu_count() // 2))
    all_merge_tasks: List[Tuple[str, int, List[Path], Path]] = []
    for key, out_dir in specs:
        out_dir.mkdir(parents=True, exist_ok=True)
        parts_by_chunk: Dict[int, List[Path]] = {}
        pattern = re.compile(r"chunk_(\d{4})\.parquet$")
        for part_dir in sorted(out_dir.glob("parts_w*")):
            for part in part_dir.glob("chunk_*.parquet"):
                match = pattern.search(part.name)
                if match:
                    parts_by_chunk.setdefault(int(match.group(1)), []).append(part)
        all_merge_tasks.extend((key, cidx, files, out_dir) for cidx, files in sorted(parts_by_chunk.items()))

    merge_parts = str(getattr(args, "merge_parts", "off")).lower() not in ("off", "0", "false", "no")
    if all_merge_tasks and merge_parts:
        _set_stage("LOCALTE_MERGE_SPLIT")
        print(
            f"[LocalTE] Merging split chunk parts: {len(all_merge_tasks):,} chunk files "
            f"across {len(specs)} outputs with {max(1, merge_workers)} workers."
        )
        if merge_workers and merge_workers > 1 and len(all_merge_tasks) > 1:
            with mp.Pool(processes=merge_workers) as pool:
                list(
                    tqdm(
                        pool.imap_unordered(_merge_one_labeled_task, all_merge_tasks, chunksize=1),
                        total=len(all_merge_tasks),
                        desc="Merge LocalTE chunks",
                        unit="chunk",
                    )
                )
        else:
            for task in tqdm(all_merge_tasks, desc="Merge LocalTE chunks", unit="chunk"):
                _merge_one_labeled_task(task)
    elif all_merge_tasks:
        print(
            f"[LocalTE] Merge skipped: keeping {len(all_merge_tasks):,} chunk partitions "
            f"across {len(specs)} outputs (LOCAL_TE_MERGE_PARTS=off)."
        )

    for key, out_dir in specs:
        if not args.disable_scores:
            score_parts = sorted(out_dir.glob("parts_w*/TE_result_scores.parquet"))
            if score_parts and merge_parts:
                _merge_scores_parts(score_parts, out_dir / "TE_result_scores.parquet", args.compression, not args.no_dictionary)
            elif score_parts:
                print(f"[LocalTE] {key}: keeping {len(score_parts)} score part files (merge_parts=off).")

    chunk_size = max(1, int(args.chunk_size))
    for i, (key, out_dir) in enumerate(specs):
        min_len = 0 if min_lengths[i] == math.inf else int(min_lengths[i])
        metadata = {
            "input": str(input_path),
            "scores_output": "" if args.disable_scores else (
                str(out_dir / "TE_result_scores.parquet") if merge_parts else "parts_w*/TE_result_scores.parquet"
            ),
            "chunk_size": chunk_size,
            "total_edges": int(total_counts[i]),
            "max_timepoints": int(max_lengths[i]),
            "min_timepoints": min_len,
            "num_chunks": math.ceil(int(max_lengths[i]) / chunk_size) if max_lengths[i] else 0,
            "buffer_edges": int(args.buffer_edges),
            "values_dtype": args.values_dtype,
            "use_threads": str(args.use_threads).lower() == "on",
            "compression": args.compression,
            "zstd_level": args.zstd_level,
            "dictionary": (not args.no_dictionary),
            "workers": len(active),
            "layout": "merged_chunks" if merge_parts else "partitioned_parts",
            "merge_parts": merge_parts,
        }
        with (out_dir / "metadata.json").open("w") as handle:
            json.dump(metadata, handle, indent=2)
        print(f"[LocalTE] {key}: {total_counts[i]:,} edges -> {out_dir}")

    if merge_parts:
        for _, out_dir in specs:
            for part_dir in sorted(out_dir.glob("parts_w*")):
                for part in part_dir.glob("*.parquet"):
                    try:
                        part.unlink()
                    except Exception:
                        pass
                try:
                    part_dir.rmdir()
                except Exception:
                    pass
    return True


def _export_known_splits_single_pass(args, input_specs: Sequence[Tuple[Path, List[Path]]]) -> bool:
    """Export standard Matrix_generate splits in one scan of TE_result_all.parquet."""
    input_path = input_specs[0][0]
    if not input_path.exists():
        raise FileNotFoundError(f"Input parquet not found: {input_path}")
    names = _load_gene_names_list()
    name_to_idx = {name: idx for idx, name in enumerate(names, start=1)}
    is_peak = np.zeros(len(names) + 1, dtype=bool)
    for idx, name in enumerate(names, start=1):
        is_peak[idx] = str(name).startswith("chr")

    base_output = Path(args.output_dir)
    base_output.mkdir(parents=True, exist_ok=True)
    seen_names: Dict[str, int] = {}
    specs = []
    tf_sources: set[int] = set()
    for _, selector_paths in input_specs:
        key = _known_split_key(selector_paths)
        count = seen_names.get(key, 0)
        sub_name = key if count == 0 else f"{key}_{count}"
        seen_names[key] = count + 1
        out_dir = base_output / sub_name
        specs.append((key, out_dir))
        if key in {"TE_TF_GN", "TE_TF_PK"}:
            tf_sources.update(_read_unique_source_indices(selector_paths[0], name_to_idx))

    if _export_known_splits_parallel(args, input_path, specs, is_peak, tf_sources):
        return True

    states = [_new_export_state(out_dir) for _, out_dir in specs]
    pf = pq.ParquetFile(input_path)
    required_cols = ["Source", "Target", "TE", "LocalTE_bytes", "LocalTE_len", "LocalTE_dtype"]
    schema_names = set(pf.schema.names)
    missing = set(required_cols) - schema_names
    if missing:
        raise ValueError(f"Input parquet missing required columns: {sorted(missing)}")
    if "LocalTE_codec" in schema_names:
        required_cols.append("LocalTE_codec")

    chunk_size = max(1, int(args.chunk_size))
    buffer_edges = max(1, int(args.buffer_edges))
    read_batch = max(1, int(args.read_batch_rows))
    use_threads = str(args.use_threads).lower() == "on"
    values_dtype = np.float16 if args.values_dtype == "float16" else np.float32
    value_type = pa.float16() if values_dtype is np.float16 else pa.float32()
    use_dictionary = not args.no_dictionary
    scores_schema = pa.schema(
        [
            ("Source", pa.int32()),
            ("Target", pa.int32()),
            ("TE", pa.float32()),
            ("LocalTE_total_len", pa.int32()),
        ]
    )

    total_rows = pf.metadata.num_rows
    print(
        f"[LocalTE] Single-pass split export: {len(specs)} outputs, "
        f"{total_rows:,} TE rows scanned once."
    )
    progress = tqdm(total=total_rows, desc="LocalTE export", unit="row")
    routed_counts = [0 for _ in specs]
    for batch in pf.iter_batches(columns=required_cols, batch_size=read_batch, use_threads=use_threads):
        batch_sources = batch.column("Source").to_numpy(zero_copy_only=False)
        batch_targets = batch.column("Target").to_numpy(zero_copy_only=False)
        batch_scores = batch.column("TE").to_numpy(zero_copy_only=False)
        batch_bytes = batch.column("LocalTE_bytes").to_pylist()
        batch_lens_np = batch.column("LocalTE_len").to_numpy(zero_copy_only=False)
        batch_dtypes = batch.column("LocalTE_dtype").to_pylist()
        batch_codecs = batch.column("LocalTE_codec").to_pylist() if "LocalTE_codec" in batch.schema.names else None

        for idx, (s, t, score, blob, length, dtype_str) in enumerate(zip(
            batch_sources, batch_targets, batch_scores, batch_bytes, batch_lens_np, batch_dtypes
        )):
            progress.update(1)
            if blob is None or length is None or dtype_str is None or int(length) == 0:
                continue
            src_i = int(s)
            tgt_i = int(t)
            out_indices = _route_known_split_indices(src_i, tgt_i, is_peak, tf_sources, [key for key, _ in specs])
            if not out_indices:
                continue
            codec = batch_codecs[idx] if batch_codecs is not None else None
            arr = _decode_local_te(blob, dtype_str, int(length), codec)
            if arr.size == 0:
                continue
            values = arr.astype(values_dtype, copy=False) if arr.dtype != values_dtype else arr
            for out_idx in out_indices:
                routed_counts[out_idx] += 1
                _append_localte_to_state(
                    states[out_idx],
                    src_i,
                    tgt_i,
                    float(score),
                    values,
                    int(length),
                    chunk_size,
                    buffer_edges,
                    value_type,
                    args.compression,
                    args.zstd_level,
                    use_dictionary,
                    args.disable_scores,
                    scores_schema,
                )
    progress.close()

    for state in states:
        _finalize_export_state(state, args, chunk_size, value_type)
    for (key, out_dir), count in zip(specs, routed_counts):
        print(f"[LocalTE] {key}: {count:,} edges -> {out_dir}")
    return True


def main() -> None:
    parser = make_parser()
    args = parser.parse_args()

    input_specs: List[Tuple[Path, List[Path]]] = []
    for raw in args.input:
        if "=" in raw:
            payload_str, selector_str = raw.split("=", 1)
            selector_paths = [Path(sel.strip()) for sel in selector_str.split(",") if sel.strip()]
        else:
            payload_str = raw
            selector_paths = []
        input_specs.append((Path(payload_str), selector_paths))

    multi_input = len(input_specs) > 1

    if multi_input and args.scores_output and not args.disable_scores:
        raise ValueError(
            "When exporting multiple inputs, omit --scores_output or enable --disable_scores to avoid filename collisions."
        )

    base_output = Path(args.output_dir)
    if multi_input:
        base_output.mkdir(parents=True, exist_ok=True)

    if _can_single_pass_known_splits(input_specs):
        _export_known_splits_single_pass(args, input_specs)
        return

    seen_names: Dict[str, int] = {}
    total_inputs = len(input_specs)
    for idx, (input_path, selector_paths) in enumerate(input_specs, start=1):
        if not input_path.exists():
            raise FileNotFoundError(f"Input parquet not found: {input_path}")
        for sel in selector_paths:
            if not Path(sel).exists():
                raise FileNotFoundError(f"Selector file not found: {sel}")

        if multi_input:
            # Name subdir from selector file if provided; otherwise use payload name
            base_for_name = selector_paths[0] if selector_paths else input_path
            name_key = _sanitize_subdir_name(base_for_name)
            count = seen_names.get(name_key, 0)
            sub_name = name_key if count == 0 else f"{name_key}_{count}"
            seen_names[name_key] = count + 1
            out_dir = base_output / sub_name
        else:
            out_dir = base_output

        local_args = argparse.Namespace(**vars(args))
        local_args.input = str(input_path)
        local_args.output_dir = str(out_dir)
        if multi_input:
            local_args.scores_output = None

        print(f"[localte_chunk_export] ({idx}/{total_inputs}) Exporting {input_path} -> {out_dir}")
        _export_dataset(local_args, selector_paths)


if __name__ == "__main__":
    main()
