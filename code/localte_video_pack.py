#!/usr/bin/env python3
"""Pack LocalTE chunk parquet files into video-like compressed floating-point tiles.

This is an experimental lossless storage layout:

  output/KEY/parts_wXX/chunk_0000.values.f16.zlib
  output/KEY/parts_wXX/chunk_0000.edges.parquet

The values file stores a dense row-major matrix of shape
edge_rows x tile_width for one chunk parquet.  The edges parquet stores row
metadata needed to interpret or reconstruct the values. Offline packing preserves
the input float16/float32/float64 dtype; streaming uses the requested value dtype.
"""

from __future__ import annotations

import argparse
import ctypes
import ctypes.util
import json
import math
import os
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm
try:
    from numba import njit
except Exception:  # pragma: no cover
    njit = None

_ZSTD_LIB = None
_ZSTD_CONTENTSIZE_ERROR = (1 << 64) - 1
_ZSTD_CONTENTSIZE_UNKNOWN = (1 << 64) - 2


def _load_zstd():
    global _ZSTD_LIB
    if _ZSTD_LIB is not None:
        return _ZSTD_LIB
    libname = ctypes.util.find_library("zstd")
    if not libname:
        raise RuntimeError("libzstd is not available")
    lib = ctypes.CDLL(libname)
    lib.ZSTD_compressBound.argtypes = [ctypes.c_size_t]
    lib.ZSTD_compressBound.restype = ctypes.c_size_t
    lib.ZSTD_compress.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.c_int,
    ]
    lib.ZSTD_compress.restype = ctypes.c_size_t
    lib.ZSTD_decompress.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.c_void_p,
        ctypes.c_size_t,
    ]
    lib.ZSTD_decompress.restype = ctypes.c_size_t
    lib.ZSTD_getFrameContentSize.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    lib.ZSTD_getFrameContentSize.restype = ctypes.c_ulonglong
    lib.ZSTD_isError.argtypes = [ctypes.c_size_t]
    lib.ZSTD_isError.restype = ctypes.c_uint
    _ZSTD_LIB = lib
    return lib


if njit is not None:
    @njit(cache=True)
    def _encode_row_codebook_u16_numba(mat):
        n, w = mat.shape
        out = np.empty(4 * (n + 1) + n * (3 + 4 * w), dtype=np.uint8)
        pos = 4 * (n + 1)
        for i in range(n + 1):
            out[4 * i] = 0
            out[4 * i + 1] = 0
            out[4 * i + 2] = 0
            out[4 * i + 3] = 0
        tmp_vals = np.empty(w, dtype=np.uint16)
        tmp_codes = np.empty(w, dtype=np.uint16)
        for i in range(n):
            u = 0
            for j in range(w):
                val = mat[i, j]
                found = -1
                for k in range(u):
                    if tmp_vals[k] == val:
                        found = k
                        break
                if found < 0:
                    tmp_vals[u] = val
                    found = u
                    u += 1
                tmp_codes[j] = found
            code_width = 1 if u <= 256 else 2
            out[pos] = code_width
            pos += 1
            out[pos] = u & 255
            out[pos + 1] = (u >> 8) & 255
            pos += 2
            for k in range(u):
                val = tmp_vals[k]
                out[pos] = val & 255
                out[pos + 1] = (val >> 8) & 255
                pos += 2
            if code_width == 1:
                for j in range(w):
                    out[pos] = tmp_codes[j] & 255
                    pos += 1
            else:
                for j in range(w):
                    c = tmp_codes[j]
                    out[pos] = c & 255
                    out[pos + 1] = (c >> 8) & 255
                    pos += 2
            off = pos - 4 * (n + 1)
            idx = i + 1
            out[4 * idx] = off & 255
            out[4 * idx + 1] = (off >> 8) & 255
            out[4 * idx + 2] = (off >> 16) & 255
            out[4 * idx + 3] = (off >> 24) & 255
        return out[:pos]
else:
    _encode_row_codebook_u16_numba = None


def warmup_video_factorized_encoder(width: int = 300) -> None:
    """Compile the numba row-codebook encoder before worker processes fork."""
    if _encode_row_codebook_u16_numba is not None:
        _encode_row_codebook_u16_numba(np.zeros((2, max(1, int(width))), dtype=np.uint16))


def _encode_row_codebook_u16(mat_u16: np.ndarray) -> bytes:
    """Encode each row as a tiny codebook plus uint8/uint16 codes."""
    if _encode_row_codebook_u16_numba is not None:
        return _encode_row_codebook_u16_numba(mat_u16).tobytes()
    chunks = []
    offsets = np.empty(mat_u16.shape[0] + 1, dtype=np.uint32)
    offsets[0] = 0
    pos = 0
    for i, row in enumerate(mat_u16):
        uniq, inv = np.unique(row, return_inverse=True)
        if uniq.size <= 256:
            part = b"\x01" + np.uint16(uniq.size).tobytes() + uniq.astype(np.uint16).tobytes() + inv.astype(np.uint8).tobytes()
        else:
            part = b"\x02" + np.uint16(uniq.size).tobytes() + uniq.astype(np.uint16).tobytes() + inv.astype(np.uint16).tobytes()
        chunks.append(part)
        pos += len(part)
        offsets[i + 1] = pos
    return offsets.tobytes() + b"".join(chunks)


def _decode_row_codebook_u16(payload: bytes, rows: int, width: int) -> np.ndarray:
    offsets = np.frombuffer(payload, dtype=np.uint32, count=rows + 1)
    base = 4 * (rows + 1)
    out = np.empty((rows, width), dtype=np.uint16)
    view = memoryview(payload)
    for i in range(rows):
        start = base + int(offsets[i])
        end = base + int(offsets[i + 1])
        rec = view[start:end]
        code_width = int(rec[0])
        n_unique = int.from_bytes(rec[1:3], "little")
        cb_start = 3
        cb_end = cb_start + 2 * n_unique
        codebook = np.frombuffer(rec[cb_start:cb_end], dtype=np.uint16)
        codes = rec[cb_end:]
        if code_width == 1:
            idx = np.frombuffer(codes, dtype=np.uint8, count=width).astype(np.int64, copy=False)
        else:
            idx = np.frombuffer(codes, dtype=np.uint16, count=width).astype(np.int64, copy=False)
        out[i, :] = codebook[idx]
    return out


def _compress_values(raw: bytes, codec: str, level: int) -> tuple[bytes, str]:
    codec = str(codec).strip().lower()
    if codec in ("none", "raw"):
        return raw, "raw"
    if codec == "zstd":
        try:
            lib = _load_zstd()
            bound = int(lib.ZSTD_compressBound(len(raw)))
            dst = ctypes.create_string_buffer(bound)
            src = ctypes.c_char_p(raw)
            n = int(lib.ZSTD_compress(dst, bound, src, len(raw), int(level)))
            if lib.ZSTD_isError(n):
                raise RuntimeError("ZSTD_compress failed")
            return bytes(dst.raw[:n]), "zstd"
        except Exception:
            # Keep the direct writer robust on systems without libzstd.
            return zlib.compress(raw, level=max(0, min(9, int(level)))), "zlib"
    return zlib.compress(raw, level=max(0, min(9, int(level)))), "zlib"


def _decompress_values(payload: bytes, codec: str, raw_size: int) -> bytes:
    codec = str(codec).strip().lower()
    if codec == "raw":
        return payload
    if codec == "zstd":
        lib = _load_zstd()
        if raw_size <= 0:
            src = ctypes.c_char_p(payload)
            frame_size = int(lib.ZSTD_getFrameContentSize(src, len(payload)))
            if frame_size in (int(_ZSTD_CONTENTSIZE_ERROR), int(_ZSTD_CONTENTSIZE_UNKNOWN)):
                raise RuntimeError("ZSTD frame content size is unavailable")
            raw_size = frame_size
        dst = ctypes.create_string_buffer(int(raw_size))
        src = ctypes.c_char_p(payload)
        n = int(lib.ZSTD_decompress(dst, int(raw_size), src, len(payload)))
        if lib.ZSTD_isError(n) or n != int(raw_size):
            raise RuntimeError("ZSTD_decompress failed")
        return bytes(dst.raw[:n])
    return zlib.decompress(payload)


def _value_ext(codec: str) -> str:
    codec = str(codec).strip().lower()
    if codec == "zstd":
        return ".values.f16.zst"
    if codec == "raw":
        return ".values.f16.raw"
    return ".values.f16.zlib"


def _value_dtype(name) -> np.dtype:
    dtype = np.dtype(name)
    if dtype.name not in ("float16", "float32", "float64"):
        raise ValueError(f"Unsupported video value dtype: {name}")
    return dtype


def _value_ext_for(codec: str, value_format: str, values_dtype: str = "float16") -> str:
    dtype = _value_dtype(values_dtype)
    if dtype != np.dtype("float16"):
        suffix = {"zstd": "zst", "raw": "raw", "none": "raw", "zlib": "zlib"}[codec]
        layout = ".rcb" if value_format == "row_codebook" else ""
        return f".values.f{dtype.itemsize * 8}{layout}.{suffix}"
    if str(value_format).lower() == "row_codebook":
        codec = str(codec).strip().lower()
        if codec == "zstd":
            return ".values.rcb.zst"
        if codec == "raw":
            return ".values.rcb.raw"
        return ".values.rcb.zlib"
    return _value_ext(codec)


def _init_video_chunk_buffers():
    from collections import defaultdict

    return defaultdict(
        lambda: {
            "sources": [],
            "targets": [],
            "lengths": [],
            "values_arrays": [],
        }
    )


def new_video_export_state(
    out_dir: Path,
    chunk_size: int,
    level: int = 1,
    values_dtype: str = "float16",
    edge_format: str = "binary",
    codec: str = "zstd",
    factorize: bool = False,
) -> dict:
    """Create state for streaming LocalTE directly into compressed video tiles."""
    values_dtype = _value_dtype(values_dtype).name
    out_dir.mkdir(parents=True, exist_ok=True)
    edge_format = str(edge_format).strip().lower()
    if edge_format not in ("binary", "parquet"):
        edge_format = "binary"
    codec = str(codec).strip().lower()
    if codec not in ("zstd", "zlib", "none", "raw"):
        codec = "zstd"
    if codec == "none":
        codec = "raw"
    return {
        "out_dir": out_dir,
        "chunk_size": int(chunk_size),
        "level": int(level),
        "values_dtype": str(values_dtype),
        "edge_format": edge_format,
        "codec": codec,
        "value_format": "row_codebook" if factorize else "dense",
        "chunk_buffers": _init_video_chunk_buffers(),
        "flush_counts": {},
        "tiles": [],
        "total": 0,
        "max_length": 0,
        "min_length": math.inf,
        "raw_value_bytes": 0,
        "compressed_value_bytes": 0,
        "edge_metadata_bytes": 0,
        "factorized_raw_bytes": 0,
    }


def _write_binary_edges(
    path: Path,
    sources: list[int],
    targets: list[int],
    chunk_idx: int,
    chunk_start: int,
    lengths: list[int],
    tile_width: int,
) -> int:
    row_count = len(sources)
    edges = np.empty((row_count, 6), dtype=np.int32)
    edges[:, 0] = np.asarray(sources, dtype=np.int32)
    edges[:, 1] = np.asarray(targets, dtype=np.int32)
    edges[:, 2] = int(chunk_idx)
    edges[:, 3] = int(chunk_start)
    edges[:, 4] = np.asarray(lengths, dtype=np.int32)
    edges[:, 5] = int(tile_width)
    path.write_bytes(edges.tobytes(order="C"))
    return int(path.stat().st_size)


def _read_edges(path: Path):
    """Read either binary or parquet video edge metadata."""
    if path.suffixes[-3:] == [".edges", ".i32", ".bin"] or path.name.endswith(".edges.i32.bin"):
        arr = np.frombuffer(path.read_bytes(), dtype=np.int32)
        if arr.size % 6 != 0:
            raise ValueError(f"Invalid binary edge metadata size: {path}")
        arr = arr.reshape((-1, 6))
        base = path.name.replace(".edges.i32.bin", "")
        matches = [path.parent / f"{base}{_value_ext_for(codec, layout, dtype)}"
                   for dtype in ("float16", "float32", "float64")
                   for layout in ("dense", "row_codebook")
                   for codec in ("zstd", "zlib", "raw")]
        matches = [candidate for candidate in matches if candidate.is_file()]
        if len(matches) != 1:
            raise ValueError(f"Expected one values tile for {path}, found {len(matches)}")
        return arr, matches[0].name
    edges = pq.read_table(path).to_pandas()
    if edges.empty:
        return np.empty((0, 6), dtype=np.int32), ""
    arr = np.column_stack(
        [
            edges["Source"].to_numpy(dtype=np.int32),
            edges["Target"].to_numpy(dtype=np.int32),
            edges["ChunkIndex"].to_numpy(dtype=np.int32),
            edges["ChunkStart"].to_numpy(dtype=np.int32),
            edges["Length"].to_numpy(dtype=np.int32),
            edges["TileWidth"].to_numpy(dtype=np.int32),
        ]
    )
    return arr.astype(np.int32, copy=False), str(edges["ValuesFile"].iloc[0])


def flush_video_chunk(state: dict, chunk_idx: int) -> None:
    """Flush one buffered LocalTE chunk while preserving its configured dtype."""
    buf = state["chunk_buffers"][chunk_idx]
    if not buf["sources"]:
        return
    out_dir = Path(state["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    flush_idx = int(state["flush_counts"].get(chunk_idx, 0))
    state["flush_counts"][chunk_idx] = flush_idx + 1
    tile_width = max(int(x) for x in buf["lengths"]) if buf["lengths"] else 0
    row_count = len(buf["sources"])
    value_dtype = _value_dtype(state["values_dtype"])
    matrix = np.zeros((row_count, tile_width), dtype=value_dtype)
    for row_idx, arr in enumerate(buf["values_arrays"]):
        vals = np.asarray(arr, dtype=value_dtype)
        matrix[row_idx, : vals.size] = vals

    value_format = str(state.get("value_format", "dense"))
    raw_dense = matrix.tobytes(order="C")
    if value_format == "row_codebook":
        raw_to_store = _encode_row_codebook_u16(matrix.view(np.uint16))
    else:
        raw_to_store = raw_dense
    compressed, codec_used = _compress_values(raw_to_store, state.get("codec", "zstd"), int(state["level"]))
    base = f"chunk_{int(chunk_idx):04d}_tile_{flush_idx:06d}"
    values_path = out_dir / f"{base}{_value_ext_for(codec_used, value_format, value_dtype.name)}"
    edge_format = str(state.get("edge_format", "binary")).lower()
    edges_path = out_dir / (f"{base}.edges.i32.bin" if edge_format == "binary" else f"{base}.edges.parquet")
    values_path.write_bytes(compressed)
    if edge_format == "binary":
        edges_bytes = _write_binary_edges(
            edges_path,
            buf["sources"],
            buf["targets"],
            int(chunk_idx),
            int(chunk_idx) * int(state["chunk_size"]),
            buf["lengths"],
            tile_width,
        )
    else:
        edges_table = pa.table(
            {
                "Source": pa.array(buf["sources"], type=pa.int32()),
                "Target": pa.array(buf["targets"], type=pa.int32()),
                "ChunkIndex": pa.array([int(chunk_idx)] * row_count, type=pa.int32()),
                "ChunkStart": pa.array([int(chunk_idx) * int(state["chunk_size"])] * row_count, type=pa.int32()),
                "Length": pa.array(buf["lengths"], type=pa.int32()),
                "TileWidth": pa.array([tile_width] * row_count, type=pa.int32()),
                "ValuesFile": pa.array([values_path.name] * row_count),
            }
        )
        pq.write_table(edges_table, edges_path, compression="snappy")
        edges_bytes = int(edges_path.stat().st_size)
    raw_bytes = len(raw_dense)
    state["factorized_raw_bytes"] += len(raw_to_store)
    compressed_bytes = len(compressed)
    state["raw_value_bytes"] += raw_bytes
    state["compressed_value_bytes"] += compressed_bytes
    state["edge_metadata_bytes"] += edges_bytes
    state["tiles"].append(
        {
            "chunk": int(chunk_idx),
            "tile": flush_idx,
            "values": str(values_path.name),
            "edges": str(edges_path.name),
            "edge_format": edge_format,
            "codec": codec_used,
            "value_format": value_format,
            "rows": row_count,
            "tile_width": int(tile_width),
            "raw_bytes": raw_bytes,
            "compressed_bytes": compressed_bytes,
            "edges_bytes": edges_bytes,
        }
    )
    buf["sources"].clear()
    buf["targets"].clear()
    buf["lengths"].clear()
    buf["values_arrays"].clear()


def append_localte_to_video_state(
    state: dict,
    src_i: int,
    tgt_i: int,
    score: float,
    values: np.ndarray,
    length: int,
    chunk_size: int,
    buffer_edges: int,
) -> None:
    """Append one edge's LocalTE vector to video-tile streaming state."""
    del score
    state["max_length"] = max(int(state["max_length"]), int(length))
    state["min_length"] = min(state["min_length"], int(length))
    state["total"] = int(state["total"]) + 1
    n_chunks = math.ceil(int(length) / int(chunk_size))
    for chunk_idx in range(n_chunks):
        start = int(chunk_idx) * int(chunk_size)
        end = min(start + int(chunk_size), int(length))
        if start >= end:
            continue
        segment = np.asarray(values[start:end], dtype=_value_dtype(state["values_dtype"]))
        buf = state["chunk_buffers"][chunk_idx]
        buf["sources"].append(int(src_i))
        buf["targets"].append(int(tgt_i))
        buf["lengths"].append(int(segment.size))
        buf["values_arrays"].append(segment)
        if len(buf["sources"]) >= int(buffer_edges):
            flush_video_chunk(state, int(chunk_idx))


def finalize_video_export_state(state: dict) -> None:
    """Flush all video tiles and write per-part metadata."""
    for chunk_idx in list(state["chunk_buffers"].keys()):
        flush_video_chunk(state, int(chunk_idx))
    min_length = 0 if state["min_length"] == math.inf else int(state["min_length"])
    raw = int(state["raw_value_bytes"])
    payload = int(state["compressed_value_bytes"]) + int(state["edge_metadata_bytes"])
    metadata = {
        "format": "localte_video_tiles_v1",
        "out_dir": str(state["out_dir"]),
        "chunk_size": int(state["chunk_size"]),
        "total_edges": int(state["total"]),
        "max_timepoints": int(state["max_length"]),
        "min_timepoints": min_length,
        "num_chunks": math.ceil(int(state["max_length"]) / int(state["chunk_size"])) if int(state["max_length"]) else 0,
        "values_dtype": state["values_dtype"],
        "compression": str(state.get("codec", "zstd")),
        "compression_level": int(state["level"]),
        "edge_format": str(state.get("edge_format", "binary")),
        "value_format": str(state.get("value_format", "dense")),
        "raw_value_bytes": raw,
        "stored_raw_value_bytes": int(state["factorized_raw_bytes"]) if state.get("value_format") == "row_codebook" else raw,
        "compressed_value_bytes": int(state["compressed_value_bytes"]),
        "edge_metadata_bytes": int(state["edge_metadata_bytes"]),
        "total_payload_bytes": payload,
        "ratio_payload_vs_raw": (float(payload) / float(raw) if raw else None),
        "tiles": state["tiles"],
    }
    (Path(state["out_dir"]) / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def _chunk_key(path: Path) -> int:
    stem = path.stem
    try:
        return int(stem.split("_")[-1])
    except Exception:
        return -1


def _relative_part_dir(path: Path, key_root: Path) -> Path:
    parent = path.parent
    try:
        return parent.relative_to(key_root)
    except Exception:
        return Path(parent.name)


def _pack_one(args_tuple):
    input_path, key_root, out_key_root, level = args_tuple
    input_path = Path(input_path)
    key_root = Path(key_root)
    out_key_root = Path(out_key_root)
    part_rel = _relative_part_dir(input_path, key_root)
    out_dir = out_key_root / part_rel
    out_dir.mkdir(parents=True, exist_ok=True)

    table = pq.read_table(input_path, columns=["Source", "Target", "ChunkIndex", "ChunkStart", "Values"])
    value_dtype = _value_dtype(table.schema.field("Values").type.value_type.to_pandas_dtype())
    data = table.to_pydict()
    values = data["Values"]
    row_count = len(values)
    lengths = [len(v) for v in values]
    max_len = max(lengths) if lengths else 0
    if row_count == 0 or max_len == 0:
        matrix = np.zeros((row_count, max_len), dtype=value_dtype)
    else:
        matrix = np.zeros((row_count, max_len), dtype=value_dtype)
        for idx, vals in enumerate(values):
            arr = np.asarray(vals, dtype=value_dtype)
            matrix[idx, : arr.size] = arr

    raw = matrix.tobytes(order="C")
    compressed = zlib.compress(raw, level=int(level))
    base = input_path.stem
    values_path = out_dir / f"{base}{_value_ext_for('zlib', 'dense', value_dtype.name)}"
    edges_path = out_dir / f"{base}.edges.parquet"
    values_path.write_bytes(compressed)

    chunk_indices = data["ChunkIndex"]
    chunk_starts = data["ChunkStart"]
    edges_table = pa.table(
        {
            "Source": pa.array(data["Source"], type=pa.int32()),
            "Target": pa.array(data["Target"], type=pa.int32()),
            "ChunkIndex": pa.array(chunk_indices, type=pa.int32()),
            "ChunkStart": pa.array(chunk_starts, type=pa.int32()),
            "Length": pa.array(lengths, type=pa.int32()),
            "TileWidth": pa.array([max_len] * row_count, type=pa.int32()),
            "ValuesFile": pa.array([values_path.name] * row_count),
        }
    )
    pq.write_table(edges_table, edges_path, compression="snappy")

    return {
        "input": str(input_path),
        "values": str(values_path),
        "edges": str(edges_path),
        "rows": row_count,
        "tile_width": int(max_len),
        "raw_bytes": len(raw),
        "compressed_bytes": len(compressed),
        "edges_bytes": int(edges_path.stat().st_size),
    }


def _discover_key_roots(root: Path, keys: list[str]) -> dict[str, Path]:
    if any((root / key).exists() for key in keys):
        return {key: root / key for key in keys if (root / key).exists()}
    return {root.name: root}


def pack(root: Path, output: Path, keys: list[str], workers: int, level: int, limit_files: int = 0) -> dict:
    key_roots = _discover_key_roots(root, keys)
    output.mkdir(parents=True, exist_ok=True)
    tasks = []
    for key, key_root in key_roots.items():
        files = sorted(key_root.rglob("chunk_*.parquet"))
        if limit_files > 0:
            files = files[:limit_files]
        out_key_root = output / key
        for path in files:
            tasks.append((str(path), str(key_root), str(out_key_root), int(level)))

    start = time.time()
    results = []
    if workers <= 1:
        iterator = (_pack_one(task) for task in tasks)
        for item in tqdm(iterator, total=len(tasks), desc="Pack LocalTE video tiles", unit="file"):
            results.append(item)
    else:
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            futures = [pool.submit(_pack_one, task) for task in tasks]
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Pack LocalTE video tiles", unit="file"):
                results.append(fut.result())

    total_raw = sum(x["raw_bytes"] for x in results)
    total_comp = sum(x["compressed_bytes"] for x in results)
    total_edges = sum(x["edges_bytes"] for x in results)
    manifest = {
        "format": "localte_video_tiles_v1",
        "source_root": str(root),
        "output_root": str(output),
        "compression": "zlib",
        "compression_level": int(level),
        "files": len(results),
        "rows": sum(x["rows"] for x in results),
        "raw_value_bytes": total_raw,
        "compressed_value_bytes": total_comp,
        "edge_metadata_bytes": total_edges,
        "total_payload_bytes": total_comp + total_edges,
        "ratio_values_vs_raw": (total_comp / total_raw if total_raw else None),
        "ratio_payload_vs_raw": ((total_comp + total_edges) / total_raw if total_raw else None),
        "elapsed_seconds": time.time() - start,
        "tiles": results,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def _parse_video_chunk_idx(path: Path) -> int:
    name = path.name
    try:
        return int(name.split("chunk_", 1)[1].split("_", 1)[0].split(".", 1)[0])
    except Exception:
        return -1


def _unpack_group(args_tuple):
    edge_paths, key_root, out_key_root = args_tuple
    edge_paths = [Path(p) for p in edge_paths]
    key_root = Path(key_root)
    out_key_root = Path(out_key_root)
    if not edge_paths:
        return {"output": None, "rows": 0, "bytes": 0}
    chunk_idx = _parse_video_chunk_idx(edge_paths[0])
    part_rel = _relative_part_dir(edge_paths[0], key_root)
    out_dir = out_key_root / part_rel
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"chunk_{chunk_idx:04d}.parquet"

    sources = []
    targets = []
    chunk_indices = []
    chunk_starts = []
    values = []
    for edges_path in sorted(edge_paths):
        edges, values_file = _read_edges(edges_path)
        if edges.size == 0:
            continue
        tile_width = int(edges[0, 5])
        value_dtype = _value_dtype("float64" if ".values.f64." in values_file else
                                   "float32" if ".values.f32." in values_file else "float16")
        dense_raw_size = int(edges.shape[0]) * int(tile_width) * value_dtype.itemsize
        codec = "zstd" if values_file.endswith(".zst") else ("raw" if values_file.endswith(".raw") else "zlib")
        is_row_codebook = ".rcb." in values_file
        raw_size = 0 if is_row_codebook and codec == "zstd" else dense_raw_size
        raw = _decompress_values((edges_path.parent / values_file).read_bytes(), codec, raw_size)
        if is_row_codebook:
            matrix_u16 = _decode_row_codebook_u16(raw, int(edges.shape[0]), int(tile_width) * value_dtype.itemsize // 2)
            matrix = matrix_u16.view(value_dtype)
        else:
            matrix = np.frombuffer(raw, dtype=value_dtype).reshape((edges.shape[0], tile_width))
        for row_idx in range(edges.shape[0]):
            length = int(edges[row_idx, 4])
            sources.append(int(edges[row_idx, 0]))
            targets.append(int(edges[row_idx, 1]))
            chunk_indices.append(int(edges[row_idx, 2]))
            chunk_starts.append(int(edges[row_idx, 3]))
            values.append(np.asarray(matrix[row_idx, :length], dtype=value_dtype))

    output_dtype = np.result_type(*(v.dtype for v in values)) if values else np.dtype("float16")
    values = [v.astype(output_dtype, copy=False) for v in values]
    table = pa.table(
        {
            "Source": pa.array(sources, type=pa.int64()),
            "Target": pa.array(targets, type=pa.int64()),
            "ChunkIndex": pa.array(chunk_indices, type=pa.int32()),
            "ChunkStart": pa.array(chunk_starts, type=pa.int32()),
            "Values": pa.array(values, type=pa.list_(pa.from_numpy_dtype(output_dtype))),
        }
    )
    pq.write_table(table, out_path, compression="snappy")
    return {"output": str(out_path), "rows": len(values), "bytes": int(out_path.stat().st_size)}


def unpack(root: Path, output: Path, keys: list[str], workers: int, limit_groups: int = 0) -> dict:
    """Reconstruct parquet LocalTE chunks from video-tile directories."""
    key_roots = _discover_key_roots(root, keys)
    output.mkdir(parents=True, exist_ok=True)
    tasks = []
    for key, key_root in key_roots.items():
        grouped = defaultdict(list)
        edge_files = list(key_root.rglob("*.edges.parquet")) + list(key_root.rglob("*.edges.i32.bin"))
        for edges_path in sorted(edge_files):
            chunk_idx = _parse_video_chunk_idx(edges_path)
            part_rel = _relative_part_dir(edges_path, key_root)
            grouped[(part_rel, chunk_idx)].append(str(edges_path))
        items = list(grouped.items())
        if limit_groups > 0:
            items = items[:limit_groups]
        out_key_root = output / key
        for (_group_key, edge_paths) in items:
            tasks.append((edge_paths, str(key_root), str(out_key_root)))

    start = time.time()
    results = []
    if workers <= 1:
        iterator = (_unpack_group(task) for task in tasks)
        for item in tqdm(iterator, total=len(tasks), desc="Unpack LocalTE video tiles", unit="chunk"):
            results.append(item)
    else:
        with ProcessPoolExecutor(max_workers=int(workers)) as pool:
            futures = [pool.submit(_unpack_group, task) for task in tasks]
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Unpack LocalTE video tiles", unit="chunk"):
                results.append(fut.result())

    manifest = {
        "format": "localte_chunk_parquet_v1",
        "source_root": str(root),
        "output_root": str(output),
        "groups": len(results),
        "rows": sum(int(x.get("rows", 0)) for x in results),
        "bytes": sum(int(x.get("bytes", 0)) for x in results),
        "elapsed_seconds": time.time() - start,
        "outputs": results,
    }
    (output / "unpack_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Pack LocalTE chunks into video-like compressed tiles.")
    parser.add_argument("root", help="local_te_split_chunks directory or one split-key directory.")
    parser.add_argument("--output", required=True, help="Output directory for video tile dataset.")
    parser.add_argument("--keys", default="TE_TF_GN,TE_TF_PK,TE_PK_GN", help="Comma-separated split keys under root.")
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--level", type=int, default=1, help="zlib compression level, usually 1 or 6.")
    parser.add_argument("--limit-files", type=int, default=0, help="Debug: only pack this many files per key.")
    parser.add_argument("--unpack", action="store_true", help="Decode video tiles back to parquet LocalTE chunks.")
    parser.add_argument("--limit-groups", type=int, default=0, help="Debug: only unpack this many part/chunk groups.")
    args = parser.parse_args()

    keys = [x.strip() for x in args.keys.split(",") if x.strip()]
    if args.unpack:
        manifest = unpack(
            Path(args.root),
            Path(args.output),
            keys,
            max(1, int(args.workers)),
            int(args.limit_groups),
        )
    else:
        manifest = pack(
            Path(args.root),
            Path(args.output),
            keys,
            max(1, int(args.workers)),
            int(args.level),
            int(args.limit_files),
        )
    print(json.dumps({k: v for k, v in manifest.items() if k != "tiles"}, indent=2))


if __name__ == "__main__":
    main()
