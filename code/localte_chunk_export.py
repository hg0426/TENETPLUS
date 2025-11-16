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

from code.path_utils import coerce_input_path


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
    gene_file = coerce_input_path("gene_names")
    if not gene_file.exists():
        raise FileNotFoundError("gene_names file not found; required to map gene symbols to indices.")
    mapping: Dict[str, int] = {}
    with gene_file.open(encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            mapping[line.strip()] = idx
    return mapping


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


def _worker_export_rowgroups(args: Tuple) -> Tuple[int, int, int]:
    (
        worker_id,
        input_path,
        out_dir,
        chunk_size,
        buffer_edges,
        read_batch,
        use_threads,
        values_dtype_str,
        row_groups,
        compression,
        zstd_level,
        use_dictionary,
        allowed_map,
    ) = args
    allowed_map = allowed_map or {}
    pf = pq.ParquetFile(input_path)
    # Detect optional codec column
    required_cols = ["Source", "Target", "TE", "LocalTE_bytes", "LocalTE_len", "LocalTE_dtype"]
    try:
        schema_names = set(pf.schema.names)
        if "LocalTE_codec" in schema_names:
            required_cols.append("LocalTE_codec")
    except Exception:
        pass
    value_type = pa.float16() if values_dtype_str == "float16" else pa.float32()
    writers: Dict[int, pq.ParquetWriter] = {}
    chunk_buffers = init_chunk_buffers()

    scores_parts_dir = Path(out_dir) / "scores_parts"
    scores_parts_dir.mkdir(parents=True, exist_ok=True)
    scores_writer: Optional[pq.ParquetWriter] = None
    total = 0
    max_len = 0
    min_len = math.inf

    for batch in pf.iter_batches(columns=required_cols, batch_size=read_batch, use_threads=use_threads, row_groups=row_groups):
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

        for idx, (s, t, blob, length, dtype_str, score) in enumerate(zip(batch_sources, batch_targets, batch_bytes, batch_lens, batch_dtypes, batch_scores)):
            if blob is None or length is None or dtype_str is None or length == 0:
                continue
            src_i = int(s)
            tgt_i = int(t)
            if allowed_map and (src_i not in allowed_map or tgt_i not in allowed_map[src_i]):
                continue
            if batch_codecs is not None:
                codec = batch_codecs[idx]
            else:
                codec = None
            if codec and str(codec).lower() == 'zlib':
                import zlib as _z
                try:
                    raw = _z.decompress(blob)
                except Exception:
                    raw = blob
                arr = np.frombuffer(raw, dtype=dtype_str, count=int(length))
            else:
                arr = np.frombuffer(blob, dtype=dtype_str, count=int(length))
            if arr.size == 0:
                continue
            values = arr.astype(np.float16 if values_dtype_str == "float16" else np.float32, copy=False)
            max_len = max(max_len, int(length))
            min_len = min(min_len, int(length))
            filt_sources.append(src_i)
            filt_targets.append(tgt_i)
            filt_scores.append(float(score))
            filt_lens.append(int(length))
            n_chunks = math.ceil(length / chunk_size)
            for cidx in range(n_chunks):
                start = cidx * chunk_size
                end = min(start + chunk_size, length)
                if start >= end:
                    continue
                seg = values[start:end]
                buf = chunk_buffers[cidx]
                buf["sources"].append(src_i)
                buf["targets"].append(tgt_i)
                buf["values_arrays"].append(seg)
                buf["offsets"].append(buf["offsets"][-1] + int(seg.size))
                if len(buf["sources"]) >= buffer_edges:
                    part_dir = Path(out_dir) / f"parts_w{worker_id:02d}"
                    part_dir.mkdir(parents=True, exist_ok=True)
                    writers.setdefault(cidx, None)
                    flush_chunk(cidx, buf, writers, part_dir, chunk_size, value_type, compression, zstd_level, use_dictionary)

        if not filt_sources:
            continue

        if scores_writer is None:
            scores_schema = pa.schema([
                ("Source", pa.int32()),
                ("Target", pa.int32()),
                ("TE", pa.float32()),
                ("LocalTE_total_len", pa.int32()),
            ])
            scores_writer = pq.ParquetWriter(scores_parts_dir / f"scores_w{worker_id:02d}.parquet", scores_schema, compression="snappy")
        tbl = pa.table({
            "Source": pa.array(filt_sources, type=pa.int32()),
            "Target": pa.array(filt_targets, type=pa.int32()),
            "TE": pa.array(np.asarray(filt_scores, dtype=np.float32), type=pa.float32()),
            "LocalTE_total_len": pa.array(filt_lens, type=pa.int32()),
        })
        scores_writer.write_table(tbl)
        total += len(filt_sources)

    part_dir = Path(out_dir) / f"parts_w{worker_id:02d}"
    part_dir.mkdir(parents=True, exist_ok=True)
    for cidx, buf in chunk_buffers.items():
        if buf["sources"]:
            writers.setdefault(cidx, None)
            flush_chunk(cidx, buf, writers, part_dir, chunk_size, value_type, compression, zstd_level, use_dictionary)

    for writer in writers.values():
        if writer is not None:
            writer.close()
    if scores_writer is not None:
        scores_writer.close()
    if min_len == math.inf:
        min_len = 0
    return (worker_id, total, max_len)


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
    required_cols = {"Source", "Target", "TE", "LocalTE_bytes", "LocalTE_len", "LocalTE_dtype"}
    missing = required_cols - set(pf.schema.names)
    if missing:
        raise ValueError(f"Input parquet missing required columns: {sorted(missing)}")

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
    # Multiprocessing path: split by row groups
    if int(getattr(args, "workers", 0) or 0) > 1:
        n_workers = int(args.workers)
        n_rg = int(pf.metadata.num_row_groups or 0)
        if n_rg <= 0:
            print("Parquet has 0 row groups; falling back to single-process export.")
        else:
            # Assign row groups round-robin to workers
            assignments = [[] for _ in range(n_workers)]
            for i in range(n_rg):
                assignments[i % n_workers].append(i)
            value_type_str = "float16" if values_dtype is np.float16 else "float32"
            comp = args.compression
            comp_level = args.zstd_level
            use_dict = not args.no_dictionary
            with mp.Pool(processes=n_workers) as pool:
                it = pool.imap_unordered(
                    _worker_export_rowgroups,
                    [
                        (
                            wid,
                            str(input_path),
                            str(out_dir),
                            chunk_size,
                            buffer_edges,
                            read_batch,
                            use_threads,
                            value_type_str,
                            assignments[wid],
                            comp,
                            comp_level,
                            use_dict,
                            allowed_map,
                        )
                        for wid in range(n_workers)
                    ],
                    chunksize=1,
                )
                # Track worker stats
                max_len = 0
                total_rows_seen = 0
                for wid, cnt, w_max in tqdm(it, total=n_workers, desc="Workers"):
                    total_rows_seen += cnt
                    max_len = max(max_len, int(w_max))

            # Merge per-worker scores into final
            scores_parts = sorted((out_dir / "scores_parts").glob("scores_w*.parquet"))
            if scores_parts:
                scores_output = Path(args.scores_output) if args.scores_output else out_dir / "TE_result_scores.parquet"
                writer = None
                for p in scores_parts:
                    pf_part = pq.ParquetFile(p)
                    for b in pf_part.iter_batches():
                        if writer is None:
                            writer = pq.ParquetWriter(scores_output, b.schema, compression="snappy")
                        writer.write_table(pa.Table.from_batches([b]))
                if writer is not None:
                    writer.close()

            # Merge per-worker chunk parts into final chunk files
            part_dirs = sorted(out_dir.glob("parts_w*/"))
            parts_by_chunk: Dict[int, List[Path]] = {}
            pattern = re.compile(r"chunk_(\d{4})\.parquet$")
            for d in part_dirs:
                for f in d.glob("chunk_*.parquet"):
                    m = pattern.search(f.name)
                    if not m:
                        continue
                    cidx = int(m.group(1))
                    parts_by_chunk.setdefault(cidx, []).append(f)

            merge_workers = int(args.merge_workers) if getattr(args, "merge_workers", None) is not None else n_workers
            if merge_workers and merge_workers > 1:
                with mp.Pool(processes=merge_workers) as pool:
                    list(tqdm(
                        pool.imap_unordered(
                            _merge_one_task,
                            [(c, files, out_dir) for c, files in parts_by_chunk.items()],
                            chunksize=1,
                        ),
                        total=len(parts_by_chunk), desc="Merging chunks"
                    ))
            else:
                for cidx, files in tqdm(sorted(parts_by_chunk.items()), desc="Merging chunks"):
                    _merge_one_task((cidx, files, out_dir))
            # Cleanup parts
            for d in part_dirs:
                for f in d.glob("*.parquet"):
                    try:
                        f.unlink()
                    except Exception:
                        pass
                try:
                    d.rmdir()
                except Exception:
                    pass
            # Done
            # Write metadata and return
            num_chunks_total = math.ceil(max_len / chunk_size) if max_len else 0
            metadata = {
                "input": str(input_path),
                "scores_output": str(Path(args.scores_output) if args.scores_output else (out_dir / "TE_result_scores.parquet")),
                "chunk_size": chunk_size,
                "total_edges": int(total_rows_seen),
                "max_timepoints": int(max_len),
                "min_timepoints": 0,
                "num_chunks": num_chunks_total,
                "buffer_edges": buffer_edges,
                "values_dtype": "float16" if values_dtype is np.float16 else "float32",
                "use_threads": use_threads,
                "workers": n_workers,
            }
            with (out_dir / "metadata.json").open("w") as handle:
                json.dump(metadata, handle, indent=2)
            print(f"Wrote chunked LocalTE to {out_dir} (chunks={num_chunks_total}), scores table -> {metadata['scores_output']}")
            return

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
            if batch_codecs is not None:
                codec = batch_codecs[idx]
            else:
                codec = None
            if codec and str(codec).lower() == 'zlib':
                import zlib as _z
                try:
                    raw = _z.decompress(blob)
                except Exception:
                    raw = blob
                arr = np.frombuffer(raw, dtype=dtype_str, count=int(length))
            else:
                arr = np.frombuffer(blob, dtype=dtype_str, count=int(length))
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
