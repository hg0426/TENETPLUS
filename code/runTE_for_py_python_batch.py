#runTE_for_py_python_batch.py
import os
from code.feature_utils import load_peak_names
import sys
os.environ.setdefault("OMP_NUM_THREADS", "1")  # avoid oversubscription with multiprocessing
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
import pandas as pd
import argparse
import numpy as np
import multiprocessing
import time
import json
import math
import shutil
from collections import OrderedDict
from pathlib import Path
from scipy import sparse
from scipy.stats import norm
from tqdm import tqdm
try:
    from numba import njit
except Exception:
    njit = None
from .Calc_TE_python import (
    prepare_dest_context,
    compute_te_for_sources_with_context,
    _chebyshev_adjacency_packed,
)
from .path_utils import locate_file
import datetime
import logging
import duckdb
try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover
    pa = None
    pq = None
try:
    from .localte_chunk_export import (
        _append_localte_to_state,
        _decode_local_te,
        _finalize_export_state,
        _new_export_state,
        _route_known_split_indices,
    )
except Exception:  # pragma: no cover - direct LocalTE export is optional
    _append_localte_to_state = None
    _decode_local_te = None
    _finalize_export_state = None
    _new_export_state = None
    _route_known_split_indices = None
try:
    from .localte_video_pack import (
        append_localte_to_video_state,
        finalize_video_export_state,
        new_video_export_state,
        warmup_video_factorized_encoder,
    )
except Exception:  # pragma: no cover - video LocalTE export is optional
    append_localte_to_video_state = None
    finalize_video_export_state = None
    new_video_export_state = None
    warmup_video_factorized_encoder = None

# Toggle storing per-timepoint local TE
STORE_LOCAL_TE = False
LOCAL_TE_DTYPE = np.float16
LOCAL_TE_DTYPE_STR = 'float16'
LOCAL_TE_CODEC = os.getenv('TE_LOCALTE_CODEC', 'zlib').lower()
LOCAL_TE_DIRECT_CHUNKS = os.getenv("LOCAL_TE_DIRECT_CHUNKS", "off").lower()
LOCAL_TE_DIRECT_WORKERS = os.getenv("LOCAL_TE_DIRECT_WORKERS", "on").lower()
LOCAL_TE_DIRECT_SCORE_PARTS = os.getenv("LOCAL_TE_DIRECT_SCORE_PARTS", "off").lower()
LOCAL_TE_VIDEO_TILES = os.getenv("LOCAL_TE_VIDEO_TILES", "off").lower()
LOCAL_TE_VIDEO_EDGE_FORMAT = os.getenv("LOCAL_TE_VIDEO_EDGE_FORMAT", "binary").lower()
LOCAL_TE_VIDEO_CODEC = os.getenv("LOCAL_TE_VIDEO_CODEC", "zstd").lower()
LOCAL_TE_VIDEO_FACTORIZE = os.getenv("LOCAL_TE_VIDEO_FACTORIZE", "off").lower()
LOCAL_TE_VIRTUAL = os.getenv("LOCAL_TE_VIRTUAL", "off").lower()
LOCAL_TE_VIRTUAL_OUTPUT_DIR = os.getenv("LOCAL_TE_VIRTUAL_OUTPUT_DIR", "local_te_virtual_array")
LOCAL_TE_VIRTUAL_VALUE_DTYPE = os.getenv("LOCAL_TE_VIRTUAL_VALUE_DTYPE", "float32").lower()
LOCAL_TE_VIRTUAL_COMPRESSION = os.getenv("LOCAL_TE_VIRTUAL_COMPRESSION", "zstd").lower()
LOCAL_TE_VIRTUAL_BUFFER_EDGES = int(os.getenv("LOCAL_TE_VIRTUAL_BUFFER_EDGES", "5000"))
LOCAL_TE_VIRTUAL_CODE_INDEX = os.getenv("LOCAL_TE_VIRTUAL_CODE_INDEX", "off").lower()
LOCAL_TE_VIRTUAL_STATE_CODE_BITS = os.getenv("LOCAL_TE_VIRTUAL_STATE_CODE_BITS", "32").strip().lower()
LOCAL_TE_CODE_INDEX_OUTPUT_DIR = os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", "local_te_code_index")
LOCAL_TE_CODE_INDEX_COMPRESSION = os.getenv("LOCAL_TE_CODE_INDEX_COMPRESSION", "zstd").lower()
LOCAL_TE_CODE_INDEX_BATCH_FEATURES = int(os.getenv("LOCAL_TE_CODE_INDEX_BATCH_FEATURES", "512"))
LOCAL_TE_RECIPE_DIR = os.getenv("LOCAL_TE_RECIPE_DIR", "local_te_recipe")
try:
    LOCAL_TE_VIDEO_LEVEL = int(os.getenv("LOCAL_TE_VIDEO_LEVEL", "1"))
except Exception:
    LOCAL_TE_VIDEO_LEVEL = 1
DENSE_THRESHOLD = int(os.getenv("TE_DENSE_THRESHOLD", "250"))
KERNEL_BACKEND = os.getenv("TE_KERNEL_BACKEND", "auto").lower()
KERNEL_DENSE_MAX_OBS = int(os.getenv("TE_KERNEL_DENSE_MAX_OBS", "4096"))
USE_DENSE_MATRIX = os.getenv("TE_USE_DENSE_MATRIX", "auto").lower()
KERNEL_PACKED_COUNTS = os.getenv("TE_KERNEL_PACKED_COUNTS", "on").lower()
KERNEL_SOURCE_PACK_CACHE_MAX = int(os.getenv("TE_KERNEL_SOURCE_PACK_CACHE_MAX", "128"))
KERNEL_SOURCE_PACK_CACHE_BYTES = int(float(os.getenv("TE_KERNEL_SOURCE_PACK_CACHE_MB", "0")) * 1024 * 1024)
KERNEL_SOURCE_PACK_WORK_BYTES = int(float(os.getenv("TE_KERNEL_SOURCE_PACK_WORK_MB", "0")) * 1024 * 1024)
KERNEL_SOURCE_PACK_SUBBATCH = int(os.getenv("TE_KERNEL_SOURCE_PACK_SUBBATCH", "0"))
KERNEL_BATCH_WORK_BYTES = int(float(os.getenv("TE_KERNEL_BATCH_WORK_MB", "0")) * 1024 * 1024)
KERNEL_DEST_CACHE_MAX = int(os.getenv("TE_KERNEL_DEST_CACHE_MAX", "64"))
KERNEL_DEST_CACHE_BYTES = int(float(os.getenv("TE_KERNEL_DEST_CACHE_MB", "0")) * 1024 * 1024)
KERNEL_SOURCE_SCALE_CACHE_MAX = int(os.getenv("TE_KERNEL_SOURCE_SCALE_CACHE_MAX", "0"))
KERNEL_SOURCE_SCALE_CACHE_BYTES = int(float(os.getenv("TE_KERNEL_SOURCE_SCALE_CACHE_MB", "0")) * 1024 * 1024)
KERNEL_SOURCE_CODE_CACHE_MAX = int(os.getenv("TE_KERNEL_SOURCE_CODE_CACHE_MAX", "0"))
KERNEL_SOURCE_CODE_CACHE_BYTES = int(float(os.getenv("TE_KERNEL_SOURCE_CODE_CACHE_MB", "64")) * 1024 * 1024)
KERNEL_SOURCE_CODE_INTEGER = os.getenv("TE_KERNEL_SOURCE_CODE_INTEGER", "auto").lower()
KERNEL_SOURCE_CODE_INTEGER_MAX_VALUE = int(os.getenv("TE_KERNEL_SOURCE_CODE_INTEGER_MAX_VALUE", "1000000"))
KERNEL_SOURCE_CODE_INTEGER_MIN_OBS = int(os.getenv("TE_KERNEL_SOURCE_CODE_INTEGER_MIN_OBS", "1500"))
KERNEL_SOURCE_PACK_MEMMAP = os.getenv("TE_KERNEL_SOURCE_PACK_MEMMAP", "auto").lower()
KERNEL_KEEP_SOURCE_PACK_MEMMAP = os.getenv("TE_KERNEL_KEEP_SOURCE_PACK_MEMMAP", "off").lower()
KERNEL_TARGET_WORK_UNITS = os.getenv("TE_KERNEL_TARGET_WORK_UNITS", "auto").lower()
KERNEL_TARGET_MIN_SOURCES = int(os.getenv("TE_KERNEL_TARGET_MIN_SOURCES", "512"))
KERNEL_PACK_VECTOR_BLOCK = int(os.getenv("TE_KERNEL_PACK_VECTOR_BLOCK", "8"))
KERNEL_SOURCE_BLOCK_SIZE = os.getenv("TE_KERNEL_SOURCE_BLOCK_SIZE", "auto").lower()
KERNEL_SOURCE_BLOCK_TARGET_GIB = float(os.getenv("TE_KERNEL_SOURCE_BLOCK_TARGET_GIB", "4.5"))
KERNEL_GROUPED_STATE = os.getenv("TE_KERNEL_GROUPED_STATE", "on").lower()
KERNEL_GROUPED_STATE_PROBE_PAIRS = int(os.getenv("TE_KERNEL_GROUPED_STATE_PROBE_PAIRS", "64"))
KERNEL_GROUPED_STATE_MIN_COMPRESSION = float(os.getenv("TE_KERNEL_GROUPED_STATE_MIN_COMPRESSION", "4.0"))
KERNEL_GROUPED_STATE_MAX_UNIQUE = int(os.getenv("TE_KERNEL_GROUPED_STATE_MAX_UNIQUE", "8192"))
KERNEL_GROUPED_STATE_MIN_OBS = int(os.getenv("TE_KERNEL_GROUPED_STATE_MIN_OBS", "1"))
KERNEL_GROUPED_STATE_WRITE_ROWS = int(os.getenv("TE_KERNEL_GROUPED_STATE_WRITE_ROWS", "250000"))
KERNEL_GROUPED_STATE_CODED = os.getenv("TE_KERNEL_GROUPED_STATE_CODED", "on").lower()
KERNEL_GROUPED_STATE_BINCOUNT = os.getenv("TE_KERNEL_GROUPED_STATE_BINCOUNT", "auto").lower()
KERNEL_GROUPED_STATE_BINCOUNT_MAX = int(os.getenv("TE_KERNEL_GROUPED_STATE_BINCOUNT_MAX", "1000000"))
KERNEL_GROUPED_STATE_SPARSE = os.getenv("TE_KERNEL_GROUPED_STATE_SPARSE", "off").lower()
KERNEL_GROUPED_STATE_SPARSE_MIN_ZERO_FRAC = float(os.getenv("TE_KERNEL_GROUPED_STATE_SPARSE_MIN_ZERO_FRAC", "0.20"))
KERNEL_GROUPED_STATE_FACTORIZED = os.getenv("TE_KERNEL_GROUPED_STATE_FACTORIZED", "auto").lower()
KERNEL_GROUPED_STATE_FACTORIZED_MAX = int(os.getenv("TE_KERNEL_GROUPED_STATE_FACTORIZED_MAX", "1000000"))
KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS = int(os.getenv("TE_KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS", "1500"))
KERNEL_GROUPED_STATE_SOURCE_PREFIX = os.getenv("TE_KERNEL_GROUPED_STATE_SOURCE_PREFIX", "auto").lower()
KERNEL_GROUPED_STATE_PROGRESS = os.getenv("TE_KERNEL_GROUPED_STATE_PROGRESS", "on").lower()
KERNEL_GROUPED_SOURCE_CODE_MEMMAP = os.getenv("TE_KERNEL_GROUPED_SOURCE_CODE_MEMMAP", "off").lower()
KERNEL_KEEP_GROUPED_SOURCE_CODE_MEMMAP = os.getenv("TE_KEEP_GROUPED_SOURCE_CODE_MEMMAP", "off").lower()
KERNEL_GROUPED_STATE_NUMBA = os.getenv("TE_KERNEL_GROUPED_STATE_NUMBA", "off").lower()
KERNEL_GROUPED_STATE_2D_RANGE = os.getenv("TE_KERNEL_GROUPED_STATE_2D_RANGE", "on").lower()
PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK = int(os.getenv("TE_PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK", "2"))
PERM_SORT_WORK = os.getenv("TE_PERM_SORT_WORK", "on").lower()
KERNEL_PERM_TABLE_SAMPLER = os.getenv("TE_KERNEL_PERM_TABLE_SAMPLER", "off").lower()
KERNEL_PERM_TABLE_MAX_CELLS = int(os.getenv("TE_KERNEL_PERM_TABLE_MAX_CELLS", "2000000"))
POOL_CHUNKSIZE = int(os.getenv("TE_POOL_CHUNKSIZE", "4"))
def configure_logging(log_path: str = "TE_analysis.log") -> None:
    """Configure runtime logging without creating files at import time."""
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s:%(message)s")
    console_level_name = os.getenv("TE_CONSOLE_LOG_LEVEL", "WARNING").upper()
    console_level = getattr(logging, console_level_name, logging.WARNING)
    abs_log_path = os.path.abspath(log_path)

    has_file = any(
        isinstance(handler, logging.FileHandler)
        and os.path.abspath(getattr(handler, "baseFilename", "")) == abs_log_path
        for handler in root_logger.handlers
    )
    if not has_file:
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    has_console = any(
        getattr(handler, "_tenet_console_handler", False)
        for handler in root_logger.handlers
    )
    if not has_console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(console_level)
        console_handler.setFormatter(formatter)
        console_handler._tenet_console_handler = True
        root_logger.addHandler(console_handler)

GLOBAL_CELL_GENE = None
GLOBAL_CELL_GENE_DENSE = None
_SOURCE_SCALE_CACHE_KERNEL = OrderedDict()
_SOURCE_SCALE_CACHE_BYTES = 0
_SOURCE_CODE_CACHE_KERNEL = OrderedDict()
_SOURCE_CODE_CACHE_BYTES = 0
_SOURCE_SPARSE_CODE_CACHE_KERNEL = OrderedDict()
_SOURCE_SPARSE_CODE_CACHE_BYTES = 0
_GROUPED_SOURCE_CODE_MEMMAP_PATH = None
_GROUPED_SOURCE_CODE_MEMMAP_ARRAY = None
_GROUPED_SOURCE_CODE_MEMMAP_INDEX = None
_GROUPED_SOURCE_CODE_VALUES = None
_GROUPED_IMPLICIT_SOURCE_INDICES = None
_SOURCE_OK_PACK_CACHE_KERNEL = OrderedDict()
_SOURCE_OK_PACK_CACHE_BYTES = 0
_SOURCE_OK_PACK_MEMMAP_PATH = None
_SOURCE_OK_PACK_MEMMAP_SHAPE = None
_SOURCE_OK_PACK_MEMMAP_INDEX = None
_SOURCE_OK_PACK_MEMMAP_ARRAY = None
_DEST_CTX_CACHE_KERNEL = OrderedDict()
_DEST_CTX_CACHE_KERNEL_BYTES = 0
GENE_NAMES = None

# Optional global timepoint subsampling indices (applied to all series)
TIME_SUBSAMPLE_INDICES = None
VIRTUAL_LOCAL_TE_MARKER = "__VIRTUAL_LOCALTE__"


def _human_int(value) -> str:
    try:
        return f"{int(value):,}"
    except Exception:
        return str(value)


def _human_bytes(num_bytes) -> str:
    try:
        value = float(num_bytes)
    except Exception:
        return str(num_bytes)
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    unit = units[0]
    for unit in units:
        if abs(value) < 1024.0 or unit == units[-1]:
            break
        value /= 1024.0
    return f"{value:.2f} {unit}"

def _maybe_subsample(series: np.ndarray) -> np.ndarray:
    global TIME_SUBSAMPLE_INDICES
    if TIME_SUBSAMPLE_INDICES is None:
        return series
    try:
        return series[TIME_SUBSAMPLE_INDICES]
    except Exception:
        return series


def _get_series(row_idx_1based: int) -> np.ndarray:
    """Fetch one feature row as float64 without repeated sparse-to-dense copies when possible."""
    idx0 = int(row_idx_1based) - 1
    if GLOBAL_CELL_GENE_DENSE is not None:
        return np.asarray(GLOBAL_CELL_GENE_DENSE[idx0], dtype=np.float64)
    return GLOBAL_CELL_GENE[idx0].toarray().ravel().astype(np.float64, copy=False)


def _select_kernel_backend(series_len: int, history_length: int) -> str:
    """Choose an exact kernel backend; auto uses dense counts only for modest T."""
    requested = str(KERNEL_BACKEND or "auto").lower()
    if requested in ("sklearn", "ckdtree", "dense"):
        return requested
    total_obs = max(0, int(series_len) - int(history_length))
    if total_obs <= int(KERNEL_DENSE_MAX_OBS):
        return "dense"
    return "sklearn"


def _get_scaled_source_kernel(source_idx: int, historyLength: int, kernel_width: float) -> np.ndarray:
    """Cache source-only scaling because it is independent of target for kernel TE."""
    global _SOURCE_SCALE_CACHE_BYTES
    key = (int(source_idx), int(historyLength), float(kernel_width))
    cached = _SOURCE_SCALE_CACHE_KERNEL.get(key)
    if cached is not None:
        _SOURCE_SCALE_CACHE_KERNEL.move_to_end(key)
        return cached
    series = _get_series(source_idx)
    series = _maybe_subsample(series)
    values = np.asarray(series, dtype=np.float64)[int(historyLength) - 1:-1]
    epsilon = 1e-9
    kw_source = float(kernel_width) * (np.std(values, ddof=1) + epsilon)
    scaled = np.divide(
        values.astype(np.float64),
        kw_source,
        out=np.zeros_like(values, dtype=np.float64),
        where=(kw_source != 0),
    )
    if KERNEL_SOURCE_SCALE_CACHE_MAX > 0:
        _SOURCE_SCALE_CACHE_KERNEL[key] = scaled
        _SOURCE_SCALE_CACHE_KERNEL.move_to_end(key)
        _SOURCE_SCALE_CACHE_BYTES += int(scaled.nbytes)
        while (
            len(_SOURCE_SCALE_CACHE_KERNEL) > KERNEL_SOURCE_SCALE_CACHE_MAX
            or (
                KERNEL_SOURCE_SCALE_CACHE_BYTES > 0
                and _SOURCE_SCALE_CACHE_BYTES > KERNEL_SOURCE_SCALE_CACHE_BYTES
            )
        ):
            _, old = _SOURCE_SCALE_CACHE_KERNEL.popitem(last=False)
            _SOURCE_SCALE_CACHE_BYTES -= int(getattr(old, "nbytes", 0))
    return scaled


def _get_grouped_source_codes_integer_kernel(
    source_idx: int,
    historyLength: int,
    kernel_width: float,
):
    """Build exact source codes in O(N) for non-negative integer count vectors."""
    requested = str(KERNEL_SOURCE_CODE_INTEGER or "auto").lower()
    if requested in ("off", "0", "false", "no"):
        return None
    series = _get_series(source_idx)
    series = _maybe_subsample(series)
    values = np.asarray(series, dtype=np.float64)[int(historyLength) - 1:-1]
    if values.size == 0 or not np.all(np.isfinite(values)):
        return None
    if requested in ("auto", "") and int(values.size) < int(KERNEL_SOURCE_CODE_INTEGER_MIN_OBS):
        return None
    if np.any(values < 0):
        return None
    raw_int = values.astype(np.int64, copy=False)
    if not np.array_equal(values, raw_int.astype(np.float64, copy=False)):
        return None
    max_value = int(raw_int.max(initial=0))
    if max_value > int(KERNEL_SOURCE_CODE_INTEGER_MAX_VALUE):
        return None
    present = np.bincount(raw_int, minlength=max_value + 1) > 0
    unique_raw = np.flatnonzero(present).astype(np.int64, copy=False)
    lookup = np.empty(max_value + 1, dtype=np.int64)
    lookup.fill(-1)
    lookup[unique_raw] = np.arange(unique_raw.size, dtype=np.int64)
    source_codes = lookup[raw_int]

    epsilon = 1e-9
    kw_source = float(kernel_width) * (np.std(values, ddof=1) + epsilon)
    source_values = np.divide(
        unique_raw.astype(np.float64),
        kw_source,
        out=np.zeros(unique_raw.size, dtype=np.float64),
        where=(kw_source != 0),
    )
    return source_values, np.asarray(source_codes, dtype=np.int64)


def _get_grouped_source_codes_kernel(
    source_idx: int,
    historyLength: int,
    kernel_width: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return cached unique source values and inverse codes for grouped-state TE."""
    global _SOURCE_CODE_CACHE_BYTES
    key = (int(source_idx), int(historyLength), float(kernel_width))
    if _GROUPED_SOURCE_CODE_MEMMAP_INDEX is not None and _GROUPED_SOURCE_CODE_VALUES is not None:
        row_idx = _GROUPED_SOURCE_CODE_MEMMAP_INDEX.get(int(source_idx))
        if row_idx is not None:
            arr = _open_grouped_source_code_memmap_readonly()
            values = _GROUPED_SOURCE_CODE_VALUES.get(int(source_idx))
            if arr is not None and values is not None:
                return values, arr[int(row_idx)]

    cached = _SOURCE_CODE_CACHE_KERNEL.get(key)
    if cached is not None:
        _SOURCE_CODE_CACHE_KERNEL.move_to_end(key)
        return cached

    fast_codes = _get_grouped_source_codes_integer_kernel(
        int(source_idx),
        int(historyLength),
        float(kernel_width),
    )
    if fast_codes is not None:
        source_values, source_codes = fast_codes
    else:
        scaled_source = _get_scaled_source_kernel(int(source_idx), int(historyLength), float(kernel_width))
        source_values, source_codes = np.unique(scaled_source, return_inverse=True)
        source_values = np.asarray(source_values, dtype=np.float64)
        source_codes = np.asarray(source_codes, dtype=np.int64)
    result = (source_values, source_codes)

    item_bytes = int(source_values.nbytes) + int(source_codes.nbytes)
    if KERNEL_SOURCE_CODE_CACHE_MAX > 0 or KERNEL_SOURCE_CODE_CACHE_BYTES > 0:
        _SOURCE_CODE_CACHE_KERNEL[key] = result
        _SOURCE_CODE_CACHE_KERNEL.move_to_end(key)
        _SOURCE_CODE_CACHE_BYTES += item_bytes
        while _SOURCE_CODE_CACHE_KERNEL and (
            (
                KERNEL_SOURCE_CODE_CACHE_MAX > 0
                and len(_SOURCE_CODE_CACHE_KERNEL) > int(KERNEL_SOURCE_CODE_CACHE_MAX)
            )
            or (
                KERNEL_SOURCE_CODE_CACHE_BYTES > 0
                and _SOURCE_CODE_CACHE_BYTES > int(KERNEL_SOURCE_CODE_CACHE_BYTES)
            )
        ):
            _, old = _SOURCE_CODE_CACHE_KERNEL.popitem(last=False)
            _SOURCE_CODE_CACHE_BYTES -= (
                int(getattr(old[0], "nbytes", 0)) + int(getattr(old[1], "nbytes", 0))
            )
    return result


def _get_grouped_source_sparse_kernel(
    source_idx: int,
    historyLength: int,
    kernel_width: float,
):
    """Return sparse source-code representation around the exact zero/default state."""
    global _SOURCE_SPARSE_CODE_CACHE_BYTES
    sparse_mode = str(KERNEL_GROUPED_STATE_SPARSE).lower()
    if sparse_mode in ("off", "0", "false", "no"):
        return None
    key = (int(source_idx), int(historyLength), float(kernel_width))
    cached = _SOURCE_SPARSE_CODE_CACHE_KERNEL.get(key)
    if cached is not None:
        _SOURCE_SPARSE_CODE_CACHE_KERNEL.move_to_end(key)
        return cached

    scaled_source = _get_scaled_source_kernel(int(source_idx), int(historyLength), float(kernel_width))
    n_obs = int(scaled_source.size)
    if n_obs <= 0:
        return None
    source_values = np.unique(scaled_source)
    zero_hits = np.flatnonzero(source_values == 0.0)
    if zero_hits.size == 0:
        return None
    nonzero_pos = np.flatnonzero(scaled_source != 0.0)
    zero_frac = 1.0 - (float(nonzero_pos.size) / float(n_obs))
    if sparse_mode in ("auto", "") and zero_frac < float(KERNEL_GROUPED_STATE_SPARSE_MIN_ZERO_FRAC):
        return None

    nonzero_codes = np.searchsorted(source_values, scaled_source[nonzero_pos]).astype(np.int32, copy=False)
    result = (
        np.asarray(source_values, dtype=np.float64),
        int(zero_hits[0]),
        np.asarray(nonzero_pos, dtype=np.int64),
        np.asarray(nonzero_codes, dtype=np.int32),
        n_obs,
    )

    item_bytes = (
        int(result[0].nbytes)
        + int(result[2].nbytes)
        + int(result[3].nbytes)
        + 16
    )
    if KERNEL_SOURCE_CODE_CACHE_MAX > 0 or KERNEL_SOURCE_CODE_CACHE_BYTES > 0:
        _SOURCE_SPARSE_CODE_CACHE_KERNEL[key] = result
        _SOURCE_SPARSE_CODE_CACHE_KERNEL.move_to_end(key)
        _SOURCE_SPARSE_CODE_CACHE_BYTES += item_bytes
        while _SOURCE_SPARSE_CODE_CACHE_KERNEL and (
            (
                KERNEL_SOURCE_CODE_CACHE_MAX > 0
                and len(_SOURCE_SPARSE_CODE_CACHE_KERNEL) > int(KERNEL_SOURCE_CODE_CACHE_MAX)
            )
            or (
                KERNEL_SOURCE_CODE_CACHE_BYTES > 0
                and _SOURCE_SPARSE_CODE_CACHE_BYTES > int(KERNEL_SOURCE_CODE_CACHE_BYTES)
            )
        ):
            _, old = _SOURCE_SPARSE_CODE_CACHE_KERNEL.popitem(last=False)
            _SOURCE_SPARSE_CODE_CACHE_BYTES -= (
                int(getattr(old[0], "nbytes", 0))
                + int(getattr(old[2], "nbytes", 0))
                + int(getattr(old[3], "nbytes", 0))
                + 16
            )
    return result


def _grouped_source_code_memmap_enabled(total_obs: int, n_sources: int, force: bool = False) -> bool:
    requested = str(KERNEL_GROUPED_SOURCE_CODE_MEMMAP or "auto").lower()
    if force:
        return int(total_obs) > 0 and int(n_sources) > 0
    if requested in ("off", "0", "false", "no"):
        return False
    if int(total_obs) <= 0 or int(n_sources) <= 0:
        return False
    if requested in ("on", "1", "true", "yes"):
        return True
    return int(n_sources) >= 128


def _open_grouped_source_code_memmap_readonly():
    global _GROUPED_SOURCE_CODE_MEMMAP_ARRAY
    if _GROUPED_SOURCE_CODE_MEMMAP_ARRAY is None and _GROUPED_SOURCE_CODE_MEMMAP_PATH:
        _GROUPED_SOURCE_CODE_MEMMAP_ARRAY = np.lib.format.open_memmap(
            _GROUPED_SOURCE_CODE_MEMMAP_PATH,
            mode="r",
        )
    return _GROUPED_SOURCE_CODE_MEMMAP_ARRAY


def _prepare_grouped_source_code_memmap(
    source_indices,
    historyLength: int,
    kernel_width: float,
    progress_dir: str,
    force: bool = False,
):
    """Precompute source inverse codes once and share them read-only across grouped-state workers."""
    global _GROUPED_SOURCE_CODE_MEMMAP_PATH, _GROUPED_SOURCE_CODE_MEMMAP_ARRAY
    global _GROUPED_SOURCE_CODE_MEMMAP_INDEX, _GROUPED_SOURCE_CODE_VALUES
    global _SOURCE_CODE_CACHE_KERNEL, _SOURCE_CODE_CACHE_BYTES

    if GLOBAL_CELL_GENE_DENSE is not None:
        series_len = GLOBAL_CELL_GENE_DENSE.shape[1]
    elif GLOBAL_CELL_GENE is not None:
        series_len = GLOBAL_CELL_GENE.shape[1]
    else:
        return None
    total_obs = max(0, int(series_len) - int(historyLength))
    unique_sources = np.unique(np.asarray(source_indices, dtype=np.int64))
    if not _grouped_source_code_memmap_enabled(total_obs, unique_sources.size, force=force):
        return None

    os.makedirs(progress_dir, exist_ok=True)
    path = os.path.abspath(os.path.join(progress_dir, "grouped_source_codes.npy"))
    mmap = np.lib.format.open_memmap(
        path,
        mode="w+",
        dtype=np.int32,
        shape=(int(unique_sources.size), int(total_obs)),
    )
    values_by_source = {}
    index_by_source = {}
    logging.info(
        "Precomputing grouped source codes: sources=%d, total_obs=%d, size=%.2f MiB, path=%s",
        unique_sources.size,
        total_obs,
        mmap.nbytes / (1024 ** 2),
        path,
    )
    for row_idx, source_idx in enumerate(unique_sources):
        scaled_source = _get_scaled_source_kernel(int(source_idx), int(historyLength), float(kernel_width))
        source_values, source_codes = np.unique(scaled_source, return_inverse=True)
        if source_codes.size != total_obs:
            raise ValueError(
                f"Source {int(source_idx)} code length {source_codes.size} != expected {total_obs}"
            )
        mmap[row_idx, :] = np.asarray(source_codes, dtype=np.int32)
        values_by_source[int(source_idx)] = np.asarray(source_values, dtype=np.float64)
        index_by_source[int(source_idx)] = int(row_idx)
    mmap.flush()
    _GROUPED_SOURCE_CODE_MEMMAP_PATH = path
    _GROUPED_SOURCE_CODE_MEMMAP_ARRAY = None
    _GROUPED_SOURCE_CODE_MEMMAP_INDEX = index_by_source
    _GROUPED_SOURCE_CODE_VALUES = values_by_source
    _SOURCE_CODE_CACHE_KERNEL = OrderedDict()
    _SOURCE_CODE_CACHE_BYTES = 0
    logging.info("Grouped source code memmap ready: %s", path)
    return path


def _clear_grouped_source_code_memmap():
    global _GROUPED_SOURCE_CODE_MEMMAP_PATH, _GROUPED_SOURCE_CODE_MEMMAP_ARRAY
    global _GROUPED_SOURCE_CODE_MEMMAP_INDEX, _GROUPED_SOURCE_CODE_VALUES
    path = _GROUPED_SOURCE_CODE_MEMMAP_PATH
    _GROUPED_SOURCE_CODE_MEMMAP_PATH = None
    _GROUPED_SOURCE_CODE_MEMMAP_ARRAY = None
    _GROUPED_SOURCE_CODE_MEMMAP_INDEX = None
    _GROUPED_SOURCE_CODE_VALUES = None
    if path and str(KERNEL_KEEP_GROUPED_SOURCE_CODE_MEMMAP).lower() not in ("on", "1", "true", "yes"):
        try:
            os.remove(path)
            logging.info("Removed grouped source code memmap: %s", path)
        except FileNotFoundError:
            pass
        except Exception as exc:
            logging.warning("Failed to remove grouped source code memmap %s: %s", path, exc)


def _kernel_source_pack_memmap_enabled(total_obs: int, n_sources: int) -> bool:
    """Decide whether to precompute source packed masks in a shared memmap."""
    requested = str(KERNEL_SOURCE_PACK_MEMMAP or "auto").lower()
    if requested in ("off", "0", "false", "no"):
        return False
    if KERNEL_PACKED_COUNTS not in ("on", "auto", "1", "true", "yes"):
        return False
    if not hasattr(np, "bitwise_count"):
        return False
    if int(total_obs) <= 0 or int(n_sources) <= 0:
        return False
    if requested in ("on", "1", "true", "yes"):
        return True
    # Auto is aimed at large source-heavy runs where rebuilding T x T source
    # masks dominates both wall time and worker-local memory.
    return int(total_obs) >= 1500 and int(n_sources) >= 256


def _open_source_pack_memmap_readonly():
    global _SOURCE_OK_PACK_MEMMAP_ARRAY
    if _SOURCE_OK_PACK_MEMMAP_ARRAY is None and _SOURCE_OK_PACK_MEMMAP_PATH:
        _SOURCE_OK_PACK_MEMMAP_ARRAY = np.lib.format.open_memmap(
            _SOURCE_OK_PACK_MEMMAP_PATH,
            mode="r",
        )
    return _SOURCE_OK_PACK_MEMMAP_ARRAY


def _has_source_pack_memmap() -> bool:
    return bool(_SOURCE_OK_PACK_MEMMAP_PATH and _SOURCE_OK_PACK_MEMMAP_INDEX)


def _prepare_kernel_source_pack_memmap(source_indices, historyLength: int, kernel_width: float, progress_dir: str):
    """Precompute source-only packed masks once and share them read-only across workers."""
    global _SOURCE_OK_PACK_MEMMAP_PATH, _SOURCE_OK_PACK_MEMMAP_SHAPE, _SOURCE_OK_PACK_MEMMAP_INDEX
    global _SOURCE_OK_PACK_MEMMAP_ARRAY, _SOURCE_OK_PACK_CACHE_KERNEL, _SOURCE_OK_PACK_CACHE_BYTES

    if GLOBAL_CELL_GENE_DENSE is not None:
        series_len = GLOBAL_CELL_GENE_DENSE.shape[1]
    elif GLOBAL_CELL_GENE is not None:
        series_len = GLOBAL_CELL_GENE.shape[1]
    else:
        return None
    if _select_kernel_backend(series_len, historyLength) != "dense":
        return None

    total_obs = max(0, int(series_len) - int(historyLength))
    unique_sources = np.array(sorted({int(s) for s in source_indices}), dtype=np.int64)
    if not _kernel_source_pack_memmap_enabled(total_obs, unique_sources.size):
        return None

    os.makedirs(progress_dir, exist_ok=True)
    packed_cols = (((total_obs + 7) // 8) + 7) // 8
    packed_dtype = np.dtype(np.uint64)
    path = os.path.abspath(os.path.join(progress_dir, "kernel_source_ok_packed.npy"))
    shape = (int(unique_sources.size), int(total_obs), int(packed_cols))
    est_gib = (np.prod(shape, dtype=np.int64) * packed_dtype.itemsize / (1024 ** 3))
    logging.info(
        "Precomputing shared kernel source masks: sources=%d, total_obs=%d, size=%.2f GiB, path=%s",
        unique_sources.size,
        total_obs,
        est_gib,
        path,
    )
    print(f"[TE] Precomputing shared source masks ({unique_sources.size} sources, {est_gib:.2f} GiB mmap)...")

    mmap_arr = np.lib.format.open_memmap(path, mode="w+", dtype=packed_dtype, shape=shape)
    for row, source_idx in enumerate(tqdm(unique_sources, desc="Precompute source masks")):
        scaled = _get_scaled_source_kernel(int(source_idx), historyLength, kernel_width)
        packed, _ = _chebyshev_adjacency_packed(scaled)
        mmap_arr[row] = packed
    mmap_arr.flush()

    _SOURCE_OK_PACK_MEMMAP_PATH = path
    _SOURCE_OK_PACK_MEMMAP_SHAPE = shape
    _SOURCE_OK_PACK_MEMMAP_INDEX = {int(s): int(i) for i, s in enumerate(unique_sources)}
    _SOURCE_OK_PACK_MEMMAP_ARRAY = mmap_arr
    _SOURCE_OK_PACK_CACHE_KERNEL.clear()
    _SOURCE_OK_PACK_CACHE_BYTES = 0
    logging.info("Shared kernel source mask memmap ready: %s", path)
    return path


def _cleanup_kernel_source_pack_memmap():
    global _SOURCE_OK_PACK_MEMMAP_PATH, _SOURCE_OK_PACK_MEMMAP_SHAPE, _SOURCE_OK_PACK_MEMMAP_INDEX
    global _SOURCE_OK_PACK_MEMMAP_ARRAY
    path = _SOURCE_OK_PACK_MEMMAP_PATH
    _SOURCE_OK_PACK_MEMMAP_ARRAY = None
    _SOURCE_OK_PACK_MEMMAP_PATH = None
    _SOURCE_OK_PACK_MEMMAP_SHAPE = None
    _SOURCE_OK_PACK_MEMMAP_INDEX = None
    if path and str(KERNEL_KEEP_SOURCE_PACK_MEMMAP).lower() not in ("on", "1", "true", "yes"):
        try:
            os.remove(path)
            logging.info("Removed shared kernel source mask memmap: %s", path)
        except FileNotFoundError:
            pass
        except Exception as e:
            logging.warning("Could not remove shared kernel source mask memmap %s: %s", path, e)


def _get_source_ok_packed_kernel(
    source_idx: int,
    historyLength: int,
    kernel_width: float,
    scaled_source: np.ndarray | None = None,
) -> np.ndarray | None:
    """Cache packed source-neighborhood matrices for exact dense kernel counts."""
    global _SOURCE_OK_PACK_CACHE_BYTES
    if KERNEL_PACKED_COUNTS not in ("on", "auto", "1", "true", "yes"):
        return None
    if not hasattr(np, "bitwise_count"):
        return None
    if _SOURCE_OK_PACK_MEMMAP_INDEX is not None and int(source_idx) in _SOURCE_OK_PACK_MEMMAP_INDEX:
        mmap_arr = _open_source_pack_memmap_readonly()
        if mmap_arr is not None:
            return mmap_arr[_SOURCE_OK_PACK_MEMMAP_INDEX[int(source_idx)]]
    key = (int(source_idx), int(historyLength), float(kernel_width))
    cached = _SOURCE_OK_PACK_CACHE_KERNEL.get(key)
    if cached is not None:
        _SOURCE_OK_PACK_CACHE_KERNEL.move_to_end(key)
        return cached
    if scaled_source is None:
        scaled_source = _get_scaled_source_kernel(source_idx, historyLength, kernel_width)
    packed, _ = _chebyshev_adjacency_packed(scaled_source)
    if KERNEL_SOURCE_PACK_CACHE_MAX > 0:
        _SOURCE_OK_PACK_CACHE_KERNEL[key] = packed
        _SOURCE_OK_PACK_CACHE_KERNEL.move_to_end(key)
        _SOURCE_OK_PACK_CACHE_BYTES += int(packed.nbytes)
        while (
            len(_SOURCE_OK_PACK_CACHE_KERNEL) > KERNEL_SOURCE_PACK_CACHE_MAX
            or (
                KERNEL_SOURCE_PACK_CACHE_BYTES > 0
                and _SOURCE_OK_PACK_CACHE_BYTES > KERNEL_SOURCE_PACK_CACHE_BYTES
            )
        ):
            _, old = _SOURCE_OK_PACK_CACHE_KERNEL.popitem(last=False)
            _SOURCE_OK_PACK_CACHE_BYTES -= int(getattr(old, "nbytes", 0))
    return packed


def _kernel_grouped_state_enabled(mode: str, pair_mode: str, history_length: int) -> bool:
    """Use exact unique-state weighted counting for sparse/discrete k=1 kernel TE."""
    requested = str(KERNEL_GROUPED_STATE or "auto").lower()
    if requested in ("off", "0", "false", "no"):
        return False
    if mode != "kernel" or pair_mode not in ("default", "gene_only", "all_feature") or int(history_length) != 1:
        return False
    if requested in ("on", "1", "true", "yes"):
        return True
    return True


def _grouped_state_unique_count(
    scaled_source: np.ndarray,
    scaled_past: np.ndarray,
    scaled_next: np.ndarray,
) -> int:
    states = np.stack((scaled_source, scaled_past, scaled_next), axis=1)
    return int(np.unique(states, axis=0).shape[0])


def _grouped_state_te_from_scaled(
    scaled_source: np.ndarray,
    scaled_past: np.ndarray,
    scaled_next: np.ndarray,
    max_unique: int = 0,
) -> float | None:
    """Exact kernel TE for k=1 by collapsing duplicate (source,past,next) states."""
    states, weights = np.unique(
        np.stack((scaled_source, scaled_past, scaled_next), axis=1),
        axis=0,
        return_counts=True,
    )
    if states.size == 0:
        return 0.0
    if int(max_unique) > 0 and states.shape[0] > int(max_unique):
        return None
    us = states[:, 0]
    up = states[:, 1]
    uy = states[:, 2]
    w = weights.astype(np.int64, copy=False)

    past_ok = np.abs(up[:, None] - up[None, :]) <= 1.0
    next_ok = np.abs(uy[:, None] - uy[None, :]) <= 1.0
    source_ok = np.abs(us[:, None] - us[None, :]) <= 1.0

    count_past = past_ok @ w
    past_next_ok = past_ok & next_ok
    count_next_past = past_next_ok @ w
    count_past_source = (past_ok & source_ok) @ w
    count_next_past_source = (past_next_ok & source_ok) @ w

    valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
    if not np.any(valid):
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        numerator = count_next_past_source[valid] / count_past_source[valid]
        denominator = count_next_past[valid] / count_past[valid]
        log_terms = np.log(numerator / denominator)
        log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
    return float(np.sum((log_terms / np.log(2.0)) * w[valid]) / scaled_source.size)


def _target_grouped_state_codes(
    scaled_past: np.ndarray,
    scaled_next: np.ndarray,
):
    target_states, target_codes = np.unique(
        np.stack((scaled_past, scaled_next), axis=1),
        axis=0,
        return_inverse=True,
    )
    return (
        np.asarray(target_codes, dtype=np.int64),
        np.asarray(target_states[:, 0], dtype=np.float64),
        np.asarray(target_states[:, 1], dtype=np.float64),
        np.bincount(np.asarray(target_codes, dtype=np.int64), minlength=target_states.shape[0]).astype(
            np.int64,
            copy=False,
        ),
    )


def _target_grouped_state_factor_context(
    target_past_values: np.ndarray,
    target_next_values: np.ndarray,
    target_counts: np.ndarray,
):
    """Precompute target-only neighborhoods reused by every source for a target."""
    target_past_left = np.searchsorted(target_past_values, target_past_values - 1.0, side="left")
    target_past_right = np.searchsorted(target_past_values, target_past_values + 1.0, side="right")
    target_past_ok = np.abs(target_past_values[:, None] - target_past_values[None, :]) <= 1.0
    target_next_ok = np.abs(target_next_values[:, None] - target_next_values[None, :]) <= 1.0
    target_past_next_ok = (target_past_ok & target_next_ok).astype(np.int64, copy=False)
    counts = np.asarray(target_counts, dtype=np.int64)
    count_prefix = np.concatenate(
        (np.zeros(1, dtype=np.int64), np.cumsum(counts, dtype=np.int64))
    )
    target_next_order = np.argsort(target_next_values, kind="mergesort")
    target_next_sorted = target_next_values[target_next_order]
    target_next_rank = np.empty(target_next_values.size, dtype=np.int64)
    target_next_rank[target_next_order] = np.arange(target_next_values.size, dtype=np.int64)
    target_next_left = np.searchsorted(
        target_next_sorted,
        target_next_values - 1.0,
        side="left",
    ).astype(np.int64, copy=False)
    target_next_right = np.searchsorted(
        target_next_sorted,
        target_next_values + 1.0,
        side="right",
    ).astype(np.int64, copy=False)
    return (
        target_past_left,
        target_past_right,
        target_past_next_ok,
        count_prefix[target_past_right] - count_prefix[target_past_left],
        target_past_next_ok @ counts,
        target_next_rank,
        target_next_left,
        target_next_right,
    )


def _source_interval_weighted_counts(
    source_values: np.ndarray,
    joint_weights: np.ndarray,
) -> np.ndarray:
    """Exact replacement for `(abs(source_i-source_j)<=1) @ joint_weights`.

    `source_values` is sorted by construction (`np.unique` or sorted integer
    values).  Each source neighborhood is therefore a contiguous interval, so a
    row-wise prefix sum gives the same integer counts without materialising the
    source-by-source boolean matrix.
    """
    n_source_states = int(source_values.size)
    if n_source_states <= 0:
        return np.zeros_like(joint_weights)
    left = np.searchsorted(source_values, source_values - 1.0, side="left")
    right = np.searchsorted(source_values, source_values + 1.0, side="right")
    prefix = np.vstack(
        (
            np.zeros((1, joint_weights.shape[1]), dtype=np.int64),
            np.cumsum(np.asarray(joint_weights, dtype=np.int64), axis=0),
        )
    )
    return prefix[right] - prefix[left]


def _target_past_interval_weighted_counts(
    target_past_left: np.ndarray,
    target_past_right: np.ndarray,
    joint_weights: np.ndarray,
) -> np.ndarray:
    """Exact replacement for `joint_weights @ target_past_ok.T`."""
    if joint_weights.size == 0:
        return np.zeros_like(joint_weights)
    prefix = np.concatenate(
        (
            np.zeros((joint_weights.shape[0], 1), dtype=np.int64),
            np.cumsum(np.asarray(joint_weights, dtype=np.int64), axis=1),
        ),
        axis=1,
    )
    return prefix[:, target_past_right] - prefix[:, target_past_left]


if njit is not None:
    @njit(cache=True)
    def _target_past_next_interval_weighted_counts_numba(
        joint_weights,
        target_past_left,
        target_past_right,
        target_next_rank,
        target_next_left,
        target_next_right,
    ):
        n_rows = joint_weights.shape[0]
        n_targets = joint_weights.shape[1]
        out = np.zeros((n_rows, n_targets), dtype=np.int64)
        for row in range(n_rows):
            bit = np.zeros(n_targets + 1, dtype=np.int64)
            left_cursor = 0
            right_cursor = 0
            for target_idx in range(n_targets):
                while right_cursor < target_past_right[target_idx]:
                    val = joint_weights[row, right_cursor]
                    if val != 0:
                        bit_idx = target_next_rank[right_cursor] + 1
                        while bit_idx <= n_targets:
                            bit[bit_idx] += val
                            bit_idx += bit_idx & -bit_idx
                    right_cursor += 1
                while left_cursor < target_past_left[target_idx]:
                    val = joint_weights[row, left_cursor]
                    if val != 0:
                        bit_idx = target_next_rank[left_cursor] + 1
                        while bit_idx <= n_targets:
                            bit[bit_idx] -= val
                            bit_idx += bit_idx & -bit_idx
                    left_cursor += 1

                total_right = 0
                bit_idx = target_next_right[target_idx]
                while bit_idx > 0:
                    total_right += bit[bit_idx]
                    bit_idx -= bit_idx & -bit_idx
                total_left = 0
                bit_idx = target_next_left[target_idx]
                while bit_idx > 0:
                    total_left += bit[bit_idx]
                    bit_idx -= bit_idx & -bit_idx
                out[row, target_idx] = total_right - total_left
        return out
else:
    _target_past_next_interval_weighted_counts_numba = None


def _target_past_next_interval_weighted_counts(
    target_factor_context,
    joint_weights: np.ndarray,
) -> np.ndarray:
    """Exact replacement for `joint_weights @ target_past_next_ok.T`.

    The target states are sorted by target-past value, so a sliding past-window
    plus a Fenwick tree over target-next ranks computes the same 2D rectangle
    counts without materialising a large target-by-target matrix product for
    every permutation sample.
    """
    if (
        _target_past_next_interval_weighted_counts_numba is not None
        and str(KERNEL_GROUPED_STATE_2D_RANGE).lower() not in ("off", "0", "false", "no")
        and len(target_factor_context) >= 8
    ):
        return _target_past_next_interval_weighted_counts_numba(
            np.asarray(joint_weights, dtype=np.int64),
            np.asarray(target_factor_context[0], dtype=np.int64),
            np.asarray(target_factor_context[1], dtype=np.int64),
            np.asarray(target_factor_context[5], dtype=np.int64),
            np.asarray(target_factor_context[6], dtype=np.int64),
            np.asarray(target_factor_context[7], dtype=np.int64),
        )
    return np.asarray(joint_weights, dtype=np.int64) @ target_factor_context[2].T


if njit is not None:
    @njit(cache=True)
    def _grouped_state_te_numba(us, up, uy, w, n_obs):
        total = 0.0
        log2 = np.log(2.0)
        m = us.size
        for i in range(m):
            count_past = 0
            count_next_past = 0
            count_past_source = 0
            count_next_past_source = 0
            us_i = us[i]
            up_i = up[i]
            uy_i = uy[i]
            for j in range(m):
                past_ok = abs(up[j] - up_i) <= 1.0
                if past_ok:
                    wj = w[j]
                    count_past += wj
                    next_ok = abs(uy[j] - uy_i) <= 1.0
                    if next_ok:
                        count_next_past += wj
                    if abs(us[j] - us_i) <= 1.0:
                        count_past_source += wj
                        if next_ok:
                            count_next_past_source += wj
            if count_past_source > 0 and count_past > 0 and count_next_past > 0:
                numerator = count_next_past_source / count_past_source
                denominator = count_next_past / count_past
                ratio = numerator / denominator
                if ratio > 0.0:
                    total += (np.log(ratio) / log2) * w[i]
        return total / n_obs
else:
    _grouped_state_te_numba = None


def _grouped_state_te_from_coded(
    source_values: np.ndarray,
    source_codes: np.ndarray,
    target_codes: np.ndarray,
    target_past_values: np.ndarray,
    target_next_values: np.ndarray,
    max_unique: int = 0,
) -> float | None:
    """Exact grouped-state TE using integer-coded unique states.

    This is algebraically the same as unique rows of
    (source, target_past, target_next), but avoids repeatedly materialising and
    sorting a 3-column float matrix for every source-target pair.
    """
    if source_values.size == 0:
        return 0.0
    n_target_states = int(target_past_values.size)
    if n_target_states <= 0:
        return 0.0
    n_source_states = int(source_values.size)
    code_space = n_source_states * n_target_states
    combined_codes = np.asarray(source_codes, dtype=np.int64) * n_target_states + target_codes
    use_bincount = str(KERNEL_GROUPED_STATE_BINCOUNT).lower() not in ("off", "0", "false", "no")
    if use_bincount and int(KERNEL_GROUPED_STATE_BINCOUNT_MAX) > 0:
        use_bincount = code_space <= int(KERNEL_GROUPED_STATE_BINCOUNT_MAX)
    if use_bincount:
        counts = np.bincount(combined_codes, minlength=code_space)
        combined_unique = np.flatnonzero(counts)
        weights = counts[combined_unique]
    else:
        combined_unique, weights = np.unique(combined_codes, return_counts=True)
    if combined_unique.size == 0:
        return 0.0
    if int(max_unique) > 0 and combined_unique.size > int(max_unique):
        return None

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    us = source_values[src_code]
    up = target_past_values[tgt_code]
    uy = target_next_values[tgt_code]
    w = weights.astype(np.int64, copy=False)

    use_numba = (
        _grouped_state_te_numba is not None
        and str(KERNEL_GROUPED_STATE_NUMBA).lower() not in ("off", "0", "false", "no")
    )
    if use_numba:
        return float(_grouped_state_te_numba(us, up, uy, w, int(source_codes.size)))

    past_ok = np.abs(up[:, None] - up[None, :]) <= 1.0
    next_ok = np.abs(uy[:, None] - uy[None, :]) <= 1.0
    source_ok = np.abs(us[:, None] - us[None, :]) <= 1.0

    count_past = past_ok @ w
    past_next_ok = past_ok & next_ok
    count_next_past = past_next_ok @ w
    count_past_source = (past_ok & source_ok) @ w
    count_next_past_source = (past_next_ok & source_ok) @ w

    valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
    if not np.any(valid):
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        numerator = count_next_past_source[valid] / count_past_source[valid]
        denominator = count_next_past[valid] / count_past[valid]
        log_terms = np.log(numerator / denominator)
        log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
    return float(np.sum((log_terms / np.log(2.0)) * w[valid]) / source_codes.size)


def _local_terms_from_combined_codes(
    combined_codes: np.ndarray,
    combined_unique: np.ndarray,
    local_by_unique: np.ndarray,
    code_space: int,
) -> np.ndarray:
    """Expand grouped local TE values back to the original observation order."""
    local_by_unique = np.asarray(local_by_unique)
    if int(code_space) > 0 and int(code_space) <= int(KERNEL_GROUPED_STATE_BINCOUNT_MAX):
        local_by_code = np.zeros(int(code_space), dtype=local_by_unique.dtype)
        local_by_code[np.asarray(combined_unique, dtype=np.int64)] = local_by_unique
        return local_by_code[np.asarray(combined_codes, dtype=np.int64)]
    positions = np.searchsorted(combined_unique, combined_codes)
    return local_by_unique[positions]


def _grouped_state_te_local_from_coded(
    source_values: np.ndarray,
    source_codes: np.ndarray,
    target_codes: np.ndarray,
    target_past_values: np.ndarray,
    target_next_values: np.ndarray,
    target_factor_context=None,
    max_unique: int = 0,
) -> tuple[float, np.ndarray] | None:
    """Exact grouped-state TE plus per-observation LocalTE for k=1 kernel mode."""
    if source_values.size == 0:
        return 0.0, np.zeros(int(source_codes.size), dtype=np.float32)
    n_target_states = int(target_past_values.size)
    if n_target_states <= 0:
        return 0.0, np.zeros(int(source_codes.size), dtype=np.float32)
    n_source_states = int(source_values.size)
    code_space = n_source_states * n_target_states
    combined_codes = np.asarray(source_codes, dtype=np.int64) * n_target_states + np.asarray(
        target_codes,
        dtype=np.int64,
    )
    use_bincount = str(KERNEL_GROUPED_STATE_BINCOUNT).lower() not in ("off", "0", "false", "no")
    if use_bincount and int(KERNEL_GROUPED_STATE_BINCOUNT_MAX) > 0:
        use_bincount = code_space <= int(KERNEL_GROUPED_STATE_BINCOUNT_MAX)
    if use_bincount:
        counts = np.bincount(combined_codes, minlength=code_space)
        combined_unique = np.flatnonzero(counts)
        weights = counts[combined_unique]
    else:
        combined_unique, weights = np.unique(combined_codes, return_counts=True)
    if combined_unique.size == 0:
        return 0.0, np.zeros(int(source_codes.size), dtype=np.float32)
    if int(max_unique) > 0 and combined_unique.size > int(max_unique):
        return None

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    w = weights.astype(np.int64, copy=False)

    if target_factor_context is not None:
        requested = str(KERNEL_GROUPED_STATE_FACTORIZED).lower()
        if (
            requested in ("off", "0", "false", "no")
            or (requested in ("auto", "") and int(source_codes.size) < int(KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS))
            or (int(KERNEL_GROUPED_STATE_FACTORIZED_MAX) > 0 and code_space > int(KERNEL_GROUPED_STATE_FACTORIZED_MAX))
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
        use_source_prefix = str(KERNEL_GROUPED_STATE_SOURCE_PREFIX).lower() not in ("off", "0", "false", "no")
        if use_source_prefix:
            source_weighted = _source_interval_weighted_counts(source_values, W)
        else:
            source_ok = (np.abs(source_values[:, None] - source_values[None, :]) <= 1.0).astype(
                np.int64,
                copy=False,
            )
            source_weighted = source_ok @ W
        count_past_source_matrix = _target_past_interval_weighted_counts(
            target_past_left,
            target_past_right,
            source_weighted,
        )
        count_next_past_source_matrix = _target_past_next_interval_weighted_counts(
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
    local_values64 = _local_terms_from_combined_codes(
        combined_codes,
        combined_unique,
        local_by_unique,
        code_space,
    )
    # Match the scalar grouped-state reduction before storing local values.
    te = float(np.sum(local_by_unique[valid] * w[valid]) / source_codes.size)
    return te, local_values64.astype(np.float32, copy=False)


def _virtual_localte_value_type(dtype: str):
    dtype = str(dtype or "float32").lower()
    if dtype in ("float16", "fp16", "half"):
        return pa.float16(), np.float16, "float16"
    return pa.float32(), np.float32, "float32"


def _virtual_state_code_type():
    bits = str(os.getenv("LOCAL_TE_VIRTUAL_STATE_CODE_BITS", LOCAL_TE_VIRTUAL_STATE_CODE_BITS)).strip().lower()
    if bits in ("32", "int32", "i32"):
        return pa.int32(), np.int32, "int32"
    return pa.int64(), np.int64, "int64"


def _virtual_localte_schema(value_type) -> pa.Schema:
    code_type, _code_np_dtype, _code_dtype_str = _virtual_state_code_type()
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


def _grouped_state_te_virtual_from_coded(
    source_values: np.ndarray,
    source_codes: np.ndarray,
    target_codes: np.ndarray,
    target_past_values: np.ndarray,
    target_next_values: np.ndarray,
    target_factor_context=None,
    max_unique: int = 0,
    value_dtype=np.float32,
) -> tuple[float, int, int, np.ndarray, np.ndarray] | None:
    """Exact TE plus compact LocalTE state table for virtual LocalTE storage."""
    n_obs = int(source_codes.size)
    if source_values.size == 0:
        return 0.0, 0, n_obs, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_dtype)
    n_target_states = int(target_past_values.size)
    if n_target_states <= 0:
        return 0.0, n_target_states, n_obs, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_dtype)
    n_source_states = int(source_values.size)
    code_space = n_source_states * n_target_states
    combined_codes = np.asarray(source_codes, dtype=np.int64) * n_target_states + np.asarray(
        target_codes,
        dtype=np.int64,
    )
    use_bincount = str(KERNEL_GROUPED_STATE_BINCOUNT).lower() not in ("off", "0", "false", "no")
    if use_bincount and int(KERNEL_GROUPED_STATE_BINCOUNT_MAX) > 0:
        use_bincount = code_space <= int(KERNEL_GROUPED_STATE_BINCOUNT_MAX)
    if use_bincount:
        counts = np.bincount(combined_codes, minlength=code_space)
        combined_unique = np.flatnonzero(counts)
        weights = counts[combined_unique]
    else:
        combined_unique, weights = np.unique(combined_codes, return_counts=True)
    if combined_unique.size == 0:
        return 0.0, n_target_states, n_obs, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=value_dtype)
    if int(max_unique) > 0 and combined_unique.size > int(max_unique):
        return None

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    w = weights.astype(np.int64, copy=False)

    if target_factor_context is not None:
        requested = str(KERNEL_GROUPED_STATE_FACTORIZED).lower()
        if (
            requested in ("off", "0", "false", "no")
            or (requested in ("auto", "") and int(source_codes.size) < int(KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS))
            or (int(KERNEL_GROUPED_STATE_FACTORIZED_MAX) > 0 and code_space > int(KERNEL_GROUPED_STATE_FACTORIZED_MAX))
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
        use_source_prefix = str(KERNEL_GROUPED_STATE_SOURCE_PREFIX).lower() not in ("off", "0", "false", "no")
        if use_source_prefix:
            source_weighted = _source_interval_weighted_counts(source_values, W)
        else:
            source_ok = (np.abs(source_values[:, None] - source_values[None, :]) <= 1.0).astype(
                np.int64,
                copy=False,
            )
            source_weighted = source_ok @ W
        count_past_source_matrix = _target_past_interval_weighted_counts(
            target_past_left,
            target_past_right,
            source_weighted,
        )
        count_next_past_source_matrix = _target_past_next_interval_weighted_counts(
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
    te = float(np.sum(local_by_unique[valid] * w[valid]) / n_obs)
    return (
        te,
        n_target_states,
        n_obs,
        np.asarray(combined_unique, dtype=np.int64),
        np.asarray(local_by_unique, dtype=value_dtype),
    )


def _grouped_state_te_from_factorized(
    source_values: np.ndarray,
    source_codes: np.ndarray,
    target_codes: np.ndarray,
    target_counts: np.ndarray,
    target_past_values: np.ndarray,
    target_next_values: np.ndarray,
    target_factor_context,
    max_unique: int = 0,
) -> float | None:
    """Exact grouped-state TE using source x target count matrices.

    This computes the same joint-state counts as `_grouped_state_te_from_coded`,
    but factorizes neighborhoods into source-only and target-only matrices.
    """
    if source_values.size == 0:
        return 0.0
    n_source_states = int(source_values.size)
    n_target_states = int(target_past_values.size)
    if n_target_states <= 0:
        return 0.0
    code_space = n_source_states * n_target_states
    requested = str(KERNEL_GROUPED_STATE_FACTORIZED).lower()
    if requested in ("off", "0", "false", "no"):
        return None
    if requested in ("auto", "") and int(source_codes.size) < int(KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS):
        return None
    if int(KERNEL_GROUPED_STATE_FACTORIZED_MAX) > 0 and code_space > int(KERNEL_GROUPED_STATE_FACTORIZED_MAX):
        return None

    combined_codes = np.asarray(source_codes, dtype=np.int64) * n_target_states + target_codes
    counts = np.bincount(combined_codes, minlength=code_space)
    combined_unique = np.flatnonzero(counts)
    if combined_unique.size == 0:
        return 0.0
    if int(max_unique) > 0 and combined_unique.size > int(max_unique):
        return None
    W = counts.reshape(n_source_states, n_target_states)

    target_past_left = target_factor_context[0]
    target_past_right = target_factor_context[1]
    count_past_by_target = target_factor_context[3]
    count_next_past_by_target = target_factor_context[4]
    use_source_prefix = str(KERNEL_GROUPED_STATE_SOURCE_PREFIX).lower() not in ("off", "0", "false", "no")
    if use_source_prefix:
        source_weighted = _source_interval_weighted_counts(source_values, W)
    else:
        source_ok = (np.abs(source_values[:, None] - source_values[None, :]) <= 1.0).astype(np.int64, copy=False)
        source_weighted = source_ok @ W
    count_past_source_matrix = _target_past_interval_weighted_counts(
        target_past_left,
        target_past_right,
        source_weighted,
    )

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    w = counts[combined_unique].astype(np.int64, copy=False)

    count_next_past_source_matrix = _target_past_next_interval_weighted_counts(
        target_factor_context,
        source_weighted,
    )
    count_past_source = count_past_source_matrix[src_code, tgt_code]
    count_next_past_source = count_next_past_source_matrix[src_code, tgt_code]
    count_past = count_past_by_target[tgt_code]
    count_next_past = count_next_past_by_target[tgt_code]

    valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
    if not np.any(valid):
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        numerator = count_next_past_source[valid] / count_past_source[valid]
        denominator = count_next_past[valid] / count_past[valid]
        log_terms = np.log(numerator / denominator)
        log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
    return float(np.sum((log_terms / np.log(2.0)) * w[valid]) / source_codes.size)


def _grouped_state_te_from_weight_matrix(
    source_values: np.ndarray,
    joint_weights: np.ndarray,
    target_context: dict,
    max_unique: int = 0,
) -> float | None:
    """Exact grouped-state TE from a source-state x target-state count table.

    Permuting source labels only changes the joint count table between source
    states and target states.  Working directly on that table avoids rebuilding
    a length-N permuted code vector and then counting it again for every null
    sample.
    """
    W = np.asarray(joint_weights, dtype=np.int64)
    if W.size == 0 or W.sum() <= 0:
        return 0.0
    n_source_states, n_target_states = W.shape
    if int(source_values.size) != int(n_source_states):
        return None
    target_past_values = target_context["target_past_values"]
    target_next_values = target_context["target_next_values"]
    if int(target_past_values.size) != int(n_target_states):
        return None

    flat_weights = W.ravel()
    combined_unique = np.flatnonzero(flat_weights)
    if combined_unique.size == 0:
        return 0.0
    if int(max_unique) > 0 and combined_unique.size > int(max_unique):
        return None

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    w = flat_weights[combined_unique].astype(np.int64, copy=False)
    target_factor_context = target_context.get("target_factor_context")

    if target_factor_context is not None:
        target_past_left = target_factor_context[0]
        target_past_right = target_factor_context[1]
        count_past_by_target = target_factor_context[3]
        count_next_past_by_target = target_factor_context[4]
        use_source_prefix = str(KERNEL_GROUPED_STATE_SOURCE_PREFIX).lower() not in (
            "off",
            "0",
            "false",
            "no",
        )
        if use_source_prefix:
            source_weighted = _source_interval_weighted_counts(source_values, W)
        else:
            source_ok = (
                np.abs(source_values[:, None] - source_values[None, :]) <= 1.0
            ).astype(np.int64, copy=False)
            source_weighted = source_ok @ W
        count_past_source_matrix = _target_past_interval_weighted_counts(
            target_past_left,
            target_past_right,
            source_weighted,
        )
        count_next_past_source_matrix = _target_past_next_interval_weighted_counts(
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

    valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
    if not np.any(valid):
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        numerator = count_next_past_source[valid] / count_past_source[valid]
        denominator = count_next_past[valid] / count_past[valid]
        log_terms = np.log(numerator / denominator)
        log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
    return float(np.sum((log_terms / np.log(2.0)) * w[valid]) / W.sum())


def _rng_multivariate_hypergeometric(rng, colors: np.ndarray, nsample: int) -> np.ndarray:
    """Draw one multivariate hypergeometric sample with a NumPy-version fallback."""
    colors = np.asarray(colors, dtype=np.int64)
    nsample = int(nsample)
    if nsample <= 0:
        return np.zeros_like(colors, dtype=np.int64)
    if hasattr(rng, "multivariate_hypergeometric"):
        return np.asarray(rng.multivariate_hypergeometric(colors, nsample), dtype=np.int64)

    draws = np.zeros_like(colors, dtype=np.int64)
    remaining_total = int(colors.sum())
    remaining_sample = nsample
    for i in range(max(0, colors.size - 1)):
        if remaining_sample <= 0 or remaining_total <= 0:
            break
        good = int(colors[i])
        bad = int(remaining_total - good)
        draw = int(rng.hypergeometric(good, bad, remaining_sample))
        draws[i] = draw
        remaining_sample -= draw
        remaining_total -= good
    if colors.size:
        draws[-1] = remaining_sample
    return draws


def _sample_permuted_joint_weights(
    source_counts: np.ndarray,
    target_counts: np.ndarray,
    rng,
) -> np.ndarray:
    """Sample the joint source-target state table induced by a random shuffle."""
    source_counts = np.asarray(source_counts, dtype=np.int64)
    target_counts = np.asarray(target_counts, dtype=np.int64)
    W = np.zeros((source_counts.size, target_counts.size), dtype=np.int64)
    active_rows = np.flatnonzero(source_counts > 0)
    if active_rows.size == 0:
        return W
    remaining_targets = target_counts.copy()
    for row in active_rows[:-1]:
        draw = _rng_multivariate_hypergeometric(
            rng,
            remaining_targets,
            int(source_counts[row]),
        )
        W[row, :] = draw
        remaining_targets -= draw
    W[int(active_rows[-1]), :] = remaining_targets
    return W


def _grouped_state_te_from_sparse_coded(
    source_values: np.ndarray,
    zero_code: int,
    nonzero_pos: np.ndarray,
    nonzero_codes: np.ndarray,
    n_obs: int,
    target_codes: np.ndarray,
    target_counts: np.ndarray,
    target_past_values: np.ndarray,
    target_next_values: np.ndarray,
    max_unique: int = 0,
) -> float | None:
    """Exact grouped-state TE using zero/default source state plus sparse corrections."""
    n_source_states = int(source_values.size)
    n_target_states = int(target_past_values.size)
    if n_source_states <= 0 or n_target_states <= 0 or int(n_obs) <= 0:
        return 0.0

    nonzero_target_codes = target_codes[nonzero_pos]
    nonzero_target_counts = np.bincount(nonzero_target_codes, minlength=n_target_states)
    default_weights = np.asarray(target_counts, dtype=np.int64) - nonzero_target_counts.astype(np.int64, copy=False)
    default_codes = int(zero_code) * n_target_states + np.flatnonzero(default_weights > 0)
    default_w = default_weights[default_weights > 0]

    nonzero_combined = np.asarray(nonzero_codes, dtype=np.int64) * n_target_states + nonzero_target_codes
    if nonzero_combined.size:
        code_space = n_source_states * n_target_states
        use_bincount = str(KERNEL_GROUPED_STATE_BINCOUNT).lower() not in ("off", "0", "false", "no")
        if use_bincount and int(KERNEL_GROUPED_STATE_BINCOUNT_MAX) > 0:
            use_bincount = code_space <= int(KERNEL_GROUPED_STATE_BINCOUNT_MAX)
        if use_bincount:
            counts = np.bincount(nonzero_combined, minlength=code_space)
            nonzero_unique = np.flatnonzero(counts)
            nonzero_w = counts[nonzero_unique]
        else:
            nonzero_unique, nonzero_w = np.unique(nonzero_combined, return_counts=True)
    else:
        nonzero_unique = np.empty(0, dtype=np.int64)
        nonzero_w = np.empty(0, dtype=np.int64)

    if default_codes.size and nonzero_unique.size:
        combined_unique = np.concatenate((default_codes.astype(np.int64, copy=False), nonzero_unique))
        weights = np.concatenate((default_w.astype(np.int64, copy=False), nonzero_w.astype(np.int64, copy=False)))
        order = np.argsort(combined_unique)
        combined_unique = combined_unique[order]
        weights = weights[order]
    elif default_codes.size:
        combined_unique = default_codes.astype(np.int64, copy=False)
        weights = default_w.astype(np.int64, copy=False)
    else:
        combined_unique = nonzero_unique.astype(np.int64, copy=False)
        weights = nonzero_w.astype(np.int64, copy=False)

    if combined_unique.size == 0:
        return 0.0
    if int(max_unique) > 0 and combined_unique.size > int(max_unique):
        return None

    src_code = combined_unique // n_target_states
    tgt_code = combined_unique - (src_code * n_target_states)
    us = source_values[src_code]
    up = target_past_values[tgt_code]
    uy = target_next_values[tgt_code]
    w = weights.astype(np.int64, copy=False)

    past_ok = np.abs(up[:, None] - up[None, :]) <= 1.0
    next_ok = np.abs(uy[:, None] - uy[None, :]) <= 1.0
    source_ok = np.abs(us[:, None] - us[None, :]) <= 1.0

    count_past = past_ok @ w
    past_next_ok = past_ok & next_ok
    count_next_past = past_next_ok @ w
    count_past_source = (past_ok & source_ok) @ w
    count_next_past_source = (past_next_ok & source_ok) @ w

    valid = (count_past_source > 0) & (count_past > 0) & (count_next_past > 0)
    if not np.any(valid):
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        numerator = count_next_past_source[valid] / count_past_source[valid]
        denominator = count_next_past[valid] / count_past[valid]
        log_terms = np.log(numerator / denominator)
        log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
    return float(np.sum((log_terms / np.log(2.0)) * w[valid]) / int(n_obs))


_PREFIX_DATA = None

class _PrefixColumns:
    """Preserve integer IDs and float64 scores without per-pair Python objects."""
    def __init__(self, sources, target, scores):
        self.sources = np.asarray(sources, dtype=np.int64)
        self.target = int(target)
        self.scores = scores

    def __len__(self):
        return len(self.sources)

def _process_prefix_target(args):
    from .exact_prefix import prepare_target, score_target
    target_idx, source_indices, history_length, kernel_width = args
    packed, values, codes = _PREFIX_DATA
    target, context = prepare_target(_get_series(int(target_idx)), sys.modules[__name__])
    source_indices_array = np.asarray(source_indices, dtype=np.int64) - 1
    scores, _ = score_target(source_indices_array, packed, target, context, values, codes,
                             sys.modules[__name__])
    return _PrefixColumns(source_indices, target_idx, scores)

def _process_grouped_state_target(args):
    target_idx, source_indices, historyLength, kernel_width = args
    try:
        dest = _maybe_subsample(_get_series(int(target_idx)))
        if dest.size <= int(historyLength):
            return [(int(s), int(target_idx), 0.0) for s in source_indices]
        dest_past = dest[:-1].astype(np.float64, copy=False)
        dest_next = dest[1:].astype(np.float64, copy=False)
        epsilon = 1e-9
        kw_past = float(kernel_width) * (np.std(dest_past, ddof=1) + epsilon)
        kw_next = float(kernel_width) * (np.std(dest_next, ddof=1) + epsilon)
        scaled_past = np.divide(
            dest_past,
            kw_past,
            out=np.zeros_like(dest_past, dtype=np.float64),
            where=(kw_past != 0),
        )
        scaled_next = np.divide(
            dest_next,
            kw_next,
            out=np.zeros_like(dest_next, dtype=np.float64),
            where=(kw_next != 0),
        )
        rows = []
        fallback_ctx = None
        max_unique = int(KERNEL_GROUPED_STATE_MAX_UNIQUE)
        virtual_localte = _is_on(LOCAL_TE_VIRTUAL)
        virtual_np_dtype = _virtual_localte_value_type(LOCAL_TE_VIRTUAL_VALUE_DTYPE)[1] if virtual_localte else None
        use_coded = str(KERNEL_GROUPED_STATE_CODED).lower() not in ("off", "0", "false", "no")
        if STORE_LOCAL_TE or virtual_localte:
            use_coded = True
        if use_coded:
            target_codes, target_past_values, target_next_values, target_counts = _target_grouped_state_codes(
                scaled_past,
                scaled_next,
            )
            target_factor_context = None
            use_factorized_target = str(KERNEL_GROUPED_STATE_FACTORIZED).lower() not in ("off", "0", "false", "no")
            if (
                use_factorized_target
                and not (
                    str(KERNEL_GROUPED_STATE_FACTORIZED).lower() in ("auto", "")
                    and int(scaled_past.shape[0]) < int(KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS)
                )
            ):
                target_factor_context = _target_grouped_state_factor_context(
                    target_past_values,
                    target_next_values,
                    target_counts,
                )
            use_sparse_source = (
                str(KERNEL_GROUPED_STATE_SPARSE).lower() not in ("off", "0", "false", "no")
                and not STORE_LOCAL_TE
                and not virtual_localte
            )
        for source_idx in source_indices:
            if use_coded:
                sparse_rep = None
                if use_sparse_source:
                    sparse_rep = _get_grouped_source_sparse_kernel(
                        int(source_idx),
                        int(historyLength),
                        float(kernel_width),
                    )
                if sparse_rep is not None and int(sparse_rep[4]) == int(scaled_past.shape[0]):
                    te = _grouped_state_te_from_sparse_coded(
                        sparse_rep[0],
                        sparse_rep[1],
                        sparse_rep[2],
                        sparse_rep[3],
                        sparse_rep[4],
                        target_codes,
                        target_counts,
                        target_past_values,
                        target_next_values,
                        max_unique=max_unique,
                    )
                else:
                    source_values, source_codes = _get_grouped_source_codes_kernel(
                        int(source_idx),
                        int(historyLength),
                        float(kernel_width),
                    )
                    if source_codes.shape[0] != scaled_past.shape[0]:
                        if STORE_LOCAL_TE:
                            empty_local = np.zeros(int(scaled_past.shape[0]), dtype=np.float32)
                            local_bytes, local_len, local_dtype, local_codec = _encode_local_te(empty_local)
                            rows.append((
                                int(source_idx),
                                int(target_idx),
                                0.0,
                                local_bytes,
                                local_len,
                                local_dtype,
                                local_codec,
                            ))
                        elif virtual_localte:
                            rows.append((
                                int(source_idx),
                                int(target_idx),
                                0.0,
                                VIRTUAL_LOCAL_TE_MARKER,
                                0,
                                int(scaled_past.shape[0]),
                                np.zeros(0, dtype=np.int64),
                                np.zeros(0, dtype=virtual_np_dtype),
                            ))
                        else:
                            rows.append((int(source_idx), int(target_idx), 0.0))
                        continue
                    te = None
                    if virtual_localte:
                        virtual_result = _grouped_state_te_virtual_from_coded(
                            source_values,
                            source_codes,
                            target_codes,
                            target_past_values,
                            target_next_values,
                            target_factor_context=target_factor_context,
                            max_unique=max_unique,
                            value_dtype=virtual_np_dtype,
                        )
                        if virtual_result is not None:
                            te, n_target_states, n_obs, state_codes, local_values = virtual_result
                            rows.append((
                                int(source_idx),
                                int(target_idx),
                                te,
                                VIRTUAL_LOCAL_TE_MARKER,
                                int(n_target_states),
                                int(n_obs),
                                state_codes,
                                local_values,
                            ))
                            continue
                    elif STORE_LOCAL_TE:
                        local_result = _grouped_state_te_local_from_coded(
                            source_values,
                            source_codes,
                            target_codes,
                            target_past_values,
                            target_next_values,
                            target_factor_context=target_factor_context,
                            max_unique=max_unique,
                        )
                        if local_result is not None:
                            te, local_values = local_result
                            local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local_values)
                            rows.append((
                                int(source_idx),
                                int(target_idx),
                                te,
                                local_bytes,
                                local_len,
                                local_dtype,
                                local_codec,
                            ))
                            continue
                    elif target_factor_context is not None:
                        te = _grouped_state_te_from_factorized(
                            source_values,
                            source_codes,
                            target_codes,
                            target_counts,
                            target_past_values,
                            target_next_values,
                            target_factor_context,
                            max_unique=max_unique,
                        )
                    if te is None:
                        te = _grouped_state_te_from_coded(
                            source_values,
                            source_codes,
                            target_codes,
                            target_past_values,
                            target_next_values,
                            max_unique=max_unique,
                        )
            else:
                scaled_source = _get_scaled_source_kernel(int(source_idx), int(historyLength), float(kernel_width))
                if scaled_source.shape[0] != scaled_past.shape[0]:
                    if STORE_LOCAL_TE:
                        empty_local = np.zeros(int(scaled_past.shape[0]), dtype=np.float32)
                        local_bytes, local_len, local_dtype, local_codec = _encode_local_te(empty_local)
                        rows.append((
                            int(source_idx),
                            int(target_idx),
                            0.0,
                            local_bytes,
                            local_len,
                            local_dtype,
                            local_codec,
                        ))
                    elif virtual_localte:
                        rows.append((
                            int(source_idx),
                            int(target_idx),
                            0.0,
                            VIRTUAL_LOCAL_TE_MARKER,
                            0,
                            int(scaled_past.shape[0]),
                            np.zeros(0, dtype=np.int64),
                            np.zeros(0, dtype=virtual_np_dtype),
                        ))
                    else:
                        rows.append((int(source_idx), int(target_idx), 0.0))
                    continue
                te = _grouped_state_te_from_scaled(scaled_source, scaled_past, scaled_next, max_unique=max_unique)
            if te is None:
                # Rare safety path for near-continuous data: keep exactness by
                # falling back to the existing per-pair kernel estimator instead
                # of materialising an oversized unique-state matrix.
                if fallback_ctx is None:
                    fallback_ctx = prepare_dest_context(
                        dest,
                        int(historyLength),
                        kernel_width=float(kernel_width),
                        normalise=True,
                        backend=_select_kernel_backend(dest.size, int(historyLength)),
                        workers=1,
                    )
                scaled_source = _get_scaled_source_kernel(int(source_idx), int(historyLength), float(kernel_width))
                source_series = _maybe_subsample(_get_series(int(source_idx)))
                if virtual_localte:
                    te_vals, local_vals = compute_te_for_sources_with_context(
                        fallback_ctx,
                        [source_series],
                        reuse_dest_neighbors=False,
                        dense_threshold=DENSE_THRESHOLD,
                        return_local=True,
                        scaled_sources=[scaled_source],
                        source_ok_packed=[None],
                        packed_vector_block=1,
                    )
                    dense_values = np.asarray(local_vals[0], dtype=virtual_np_dtype)
                    rows.append((
                        int(source_idx),
                        int(target_idx),
                        float(te_vals[0]),
                        VIRTUAL_LOCAL_TE_MARKER,
                        -1,
                        int(dense_values.size),
                        np.zeros(0, dtype=np.int64),
                        dense_values,
                    ))
                    continue
                elif STORE_LOCAL_TE:
                    te_vals, local_vals = compute_te_for_sources_with_context(
                        fallback_ctx,
                        [source_series],
                        reuse_dest_neighbors=False,
                        dense_threshold=DENSE_THRESHOLD,
                        return_local=True,
                        scaled_sources=[scaled_source],
                        source_ok_packed=[None],
                        packed_vector_block=1,
                    )
                    local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local_vals[0])
                    rows.append((
                        int(source_idx),
                        int(target_idx),
                        float(te_vals[0]),
                        local_bytes,
                        local_len,
                        local_dtype,
                        local_codec,
                    ))
                    continue
                else:
                    te = compute_te_for_sources_with_context(
                        fallback_ctx,
                        [source_series],
                        reuse_dest_neighbors=False,
                        dense_threshold=DENSE_THRESHOLD,
                        return_local=False,
                        scaled_sources=[scaled_source],
                        source_ok_packed=[None],
                        packed_vector_block=1,
                    )[0]
            rows.append((int(source_idx), int(target_idx), te))
        return rows
    except Exception as e:
        logging.warning(
            "Grouped-state kernel failed for target %s; falling back to exact per-target kernel path: %s",
            target_idx,
            e,
        )
        try:
            return _process_chunk((target_idx, source_indices, historyLength, kernel_width, "kernel"))
        except Exception as fallback_exc:
            logging.error(
                "Fallback exact kernel also failed for target %s: %s",
                target_idx,
                fallback_exc,
            )
            if STORE_LOCAL_TE:
                return [
                    (
                        int(s),
                        int(target_idx),
                        0.0,
                        *_encode_local_te(np.zeros(0, dtype=np.float32)),
                    )
                    for s in source_indices
                ]
            if _is_on(LOCAL_TE_VIRTUAL):
                return [
                    (
                        int(s),
                        int(target_idx),
                        0.0,
                        VIRTUAL_LOCAL_TE_MARKER,
                        0,
                        0,
                        np.zeros(0, dtype=np.int64),
                        np.zeros(0, dtype=_virtual_localte_value_type(LOCAL_TE_VIRTUAL_VALUE_DTYPE)[1]),
                    )
                    for s in source_indices
                ]
            return [(int(s), int(target_idx), 0.0) for s in source_indices]


def _init_grouped_state_implicit(source_indices):
    """Share implicit all-pair source indices once per worker, not per task."""
    global _GROUPED_IMPLICIT_SOURCE_INDICES
    _GROUPED_IMPLICIT_SOURCE_INDICES = np.asarray(source_indices, dtype=np.int64)


def _process_grouped_state_target_implicit(args):
    target_idx, historyLength, kernel_width = args
    indices = _GROUPED_IMPLICIT_SOURCE_INDICES
    if indices is None:
        raise RuntimeError("Implicit grouped-state source indices were not initialised.")
    target_idx = int(target_idx)
    source_indices = indices[indices != target_idx]
    return _process_grouped_state_target((target_idx, source_indices, historyLength, kernel_width))


def _grouped_state_auto_probe(list_pairs, historyLength: int, kernel_width: float) -> bool:
    """Enable grouped-state automatically only when repeated states make it worthwhile."""
    requested = str(KERNEL_GROUPED_STATE or "auto").lower()
    if requested in ("on", "1", "true", "yes"):
        return True
    if requested not in ("auto", ""):
        return False
    n_probe = min(int(KERNEL_GROUPED_STATE_PROBE_PAIRS), int(len(list_pairs)))
    if n_probe <= 0:
        return False

    if GLOBAL_CELL_GENE_DENSE is not None:
        total_obs = max(0, int(GLOBAL_CELL_GENE_DENSE.shape[1]) - int(historyLength))
    elif GLOBAL_CELL_GENE is not None:
        total_obs = max(0, int(GLOBAL_CELL_GENE.shape[1]) - int(historyLength))
    else:
        total_obs = 0
    if int(KERNEL_GROUPED_STATE_MIN_OBS) > 0 and total_obs < int(KERNEL_GROUPED_STATE_MIN_OBS):
        logging.info(
            "Grouped-state auto probe: disabled because total_obs=%d < min_obs=%d.",
            total_obs,
            int(KERNEL_GROUPED_STATE_MIN_OBS),
        )
        print(
            "[TE] Grouped-state auto probe: "
            f"enabled=False, total_obs={total_obs} < min_obs={int(KERNEL_GROUPED_STATE_MIN_OBS)}"
        )
        return False

    # Spread the probe over the pair list so one unusually sparse target/source
    # block does not decide the backend for the whole run.
    probe_indices = np.linspace(0, len(list_pairs) - 1, num=n_probe, dtype=np.int64)
    target_cache = {}
    compression = []
    unique_counts = []
    for pair_i in probe_indices:
        source_idx, target_idx = list_pairs[int(pair_i)]
        target_idx = int(target_idx)
        if target_idx not in target_cache:
            dest = _maybe_subsample(_get_series(target_idx))
            if dest.size <= int(historyLength):
                continue
            dest_past = dest[:-1].astype(np.float64, copy=False)
            dest_next = dest[1:].astype(np.float64, copy=False)
            epsilon = 1e-9
            kw_past = float(kernel_width) * (np.std(dest_past, ddof=1) + epsilon)
            kw_next = float(kernel_width) * (np.std(dest_next, ddof=1) + epsilon)
            scaled_past = np.divide(
                dest_past,
                kw_past,
                out=np.zeros_like(dest_past, dtype=np.float64),
                where=(kw_past != 0),
            )
            scaled_next = np.divide(
                dest_next,
                kw_next,
                out=np.zeros_like(dest_next, dtype=np.float64),
                where=(kw_next != 0),
            )
            target_cache[target_idx] = (scaled_past, scaled_next)
        scaled_past, scaled_next = target_cache[target_idx]
        scaled_source = _get_scaled_source_kernel(int(source_idx), int(historyLength), float(kernel_width))
        if scaled_source.shape[0] != scaled_past.shape[0]:
            continue
        n_unique = _grouped_state_unique_count(scaled_source, scaled_past, scaled_next)
        if n_unique <= 0:
            continue
        unique_counts.append(n_unique)
        compression.append(float(scaled_source.size) / float(n_unique))

    if not compression:
        return False
    median_compression = float(np.median(compression))
    max_unique = int(np.max(unique_counts))
    enabled = (
        median_compression >= float(KERNEL_GROUPED_STATE_MIN_COMPRESSION)
        and (int(KERNEL_GROUPED_STATE_MAX_UNIQUE) <= 0 or max_unique <= int(KERNEL_GROUPED_STATE_MAX_UNIQUE))
    )
    logging.info(
        "Grouped-state auto probe: enabled=%s, pairs=%d, median_compression=%.2f, max_unique=%d",
        enabled,
        len(compression),
        median_compression,
        max_unique,
    )
    print(
        "[TE] Grouped-state auto probe: "
        f"enabled={enabled}, median_compression={median_compression:.2f}x, max_unique={max_unique}"
    )
    return bool(enabled)


def _grouped_state_auto_probe_implicit(indices, historyLength: int, kernel_width: float) -> bool:
    """Probe implicit all-pair mode without materialising the full pair set."""
    requested = str(KERNEL_GROUPED_STATE or "auto").lower()
    if requested in ("on", "1", "true", "yes"):
        return True
    if requested not in ("auto", ""):
        return False

    indices = np.asarray(indices, dtype=np.int64)
    if indices.size <= 1:
        return False
    n_probe = min(int(KERNEL_GROUPED_STATE_PROBE_PAIRS), int(indices.size) * int(indices.size - 1))
    if n_probe <= 0:
        return False

    probe_pairs = []
    for tgt in indices:
        for src in indices:
            if int(src) == int(tgt):
                continue
            probe_pairs.append((int(src), int(tgt)))
            if len(probe_pairs) >= n_probe:
                return _grouped_state_auto_probe(
                    np.asarray(probe_pairs, dtype=np.int64),
                    int(historyLength),
                    float(kernel_width),
                )
    return _grouped_state_auto_probe(
        np.asarray(probe_pairs, dtype=np.int64),
        int(historyLength),
        float(kernel_width),
    )


def _bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values mapped back to input order."""
    pvals = np.asarray(pvals, dtype=np.float64)
    m = int(pvals.size)
    if m == 0:
        return pvals.copy()
    order = np.argsort(pvals)
    ranked = pvals[order] * float(m) / np.arange(1, m + 1, dtype=np.float64)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty_like(ranked)
    qvals[order] = np.clip(ranked, 0.0, 1.0)
    return qvals


def _select_grn_fdr_candidates(
    df_fast: pd.DataFrame,
    alpha: float,
    te_cutoff: float,
) -> pd.DataFrame:
    """Select GRN global z-score + BH-FDR candidates."""
    cutoff = float(te_cutoff)
    subset = df_fast.loc[df_fast["TE"].astype(float) > cutoff, ["Source", "Target", "TE"]].copy()
    if subset.shape[0] < 2:
        return subset.iloc[0:0]
    te_arr = subset["TE"].to_numpy(dtype=np.float64, copy=False)
    std = float(np.std(te_arr, ddof=0))
    if np.allclose(std, 0.0):
        return subset.iloc[0:0]
    mean = float(np.mean(te_arr))
    pvals = norm.sf((te_arr - mean) / std)
    qvals = _bh_qvalues(pvals)
    return subset.loc[qvals < float(alpha), ["Source", "Target", "TE"]].copy()


def _kernel_source_pack_subbatch_size(dest_ctx: dict, requested_count: int) -> int:
    """Limit simultaneously-live packed source masks by bytes, not just pair batch size."""
    requested_count = int(requested_count)
    if requested_count <= 1:
        return max(1, requested_count)
    if dest_ctx.get("backend") != "dense":
        return requested_count
    if KERNEL_PACKED_COUNTS not in ("on", "auto", "1", "true", "yes"):
        return requested_count
    if not hasattr(np, "bitwise_count"):
        return requested_count

    total_obs = int(dest_ctx.get("total_obs", 0))
    if total_obs <= 0:
        return requested_count
    packed_bytes_per_source = total_obs * ((((total_obs + 7) // 8) + 7) // 8) * np.dtype(np.uint64).itemsize

    limit = requested_count
    if KERNEL_SOURCE_PACK_SUBBATCH > 0:
        limit = min(limit, int(KERNEL_SOURCE_PACK_SUBBATCH))
    if KERNEL_SOURCE_PACK_WORK_BYTES > 0:
        by_bytes = max(1, int(KERNEL_SOURCE_PACK_WORK_BYTES) // max(1, packed_bytes_per_source))
        limit = min(limit, by_bytes)
    return max(1, min(requested_count, limit))


def _maybe_cap_kernel_batch_size(batch_size: int, historyLength: int, mode: str) -> int:
    """Optionally cap kernel source batch size by per-worker packed-mask memory."""
    batch_size = max(1, int(batch_size))
    if mode != "kernel" or KERNEL_BATCH_WORK_BYTES <= 0:
        return batch_size
    if KERNEL_PACKED_COUNTS not in ("on", "auto", "1", "true", "yes"):
        return batch_size
    if not hasattr(np, "bitwise_count"):
        return batch_size
    if GLOBAL_CELL_GENE_DENSE is not None:
        series_len = GLOBAL_CELL_GENE_DENSE.shape[1]
    elif GLOBAL_CELL_GENE is not None:
        series_len = GLOBAL_CELL_GENE.shape[1]
    else:
        return batch_size
    if _select_kernel_backend(series_len, historyLength) != "dense":
        return batch_size

    total_obs = max(0, int(series_len) - int(historyLength))
    if total_obs <= 0:
        return batch_size
    packed_bytes_per_source = total_obs * ((((total_obs + 7) // 8) + 7) // 8) * np.dtype(np.uint64).itemsize
    capped = max(1, int(KERNEL_BATCH_WORK_BYTES) // max(1, packed_bytes_per_source))
    if capped < batch_size:
        logging.info(
            "Capping kernel batch_size from %d to %d by TE_KERNEL_BATCH_WORK_MB "
            "(%.2f MiB per worker; %.2f MiB per packed source mask).",
            batch_size,
            capped,
            KERNEL_BATCH_WORK_BYTES / (1024 * 1024),
            packed_bytes_per_source / (1024 * 1024),
        )
        print(
            f"[TE] Capping kernel batch_size {batch_size} -> {capped} "
            f"by TE_KERNEL_BATCH_WORK_MB={KERNEL_BATCH_WORK_BYTES / (1024 * 1024):.2f}"
        )
        return capped
    return batch_size


def _resolve_kernel_source_block_size(list_pairs, historyLength: int) -> int:
    """Resolve explicit/auto source-block size for exact packed kernel mode."""
    requested = str(KERNEL_SOURCE_BLOCK_SIZE or "0").strip().lower()
    if requested in ("0", "off", "false", "no", "none"):
        return 0
    if requested not in ("auto", "on", "true", "yes"):
        try:
            return max(0, int(float(requested)))
        except Exception:
            logging.warning("Invalid TE_KERNEL_SOURCE_BLOCK_SIZE=%r; falling back to auto.", requested)

    if list_pairs is None or len(list_pairs) == 0:
        return 0
    if GLOBAL_CELL_GENE_DENSE is not None:
        series_len = GLOBAL_CELL_GENE_DENSE.shape[1]
    elif GLOBAL_CELL_GENE is not None:
        series_len = GLOBAL_CELL_GENE.shape[1]
    else:
        return 0
    if _select_kernel_backend(series_len, historyLength) != "dense":
        return 0

    total_obs = max(0, int(series_len) - int(historyLength))
    if total_obs <= 0:
        return 0
    unique_sources = np.unique(list_pairs[:, 0].astype(np.int64))
    if unique_sources.size <= 0:
        return 0
    packed_cols = (((total_obs + 7) // 8) + 7) // 8
    bytes_per_source = total_obs * packed_cols * np.dtype(np.uint64).itemsize
    target_bytes = max(1.0, float(KERNEL_SOURCE_BLOCK_TARGET_GIB)) * (1024 ** 3)
    block_size = max(1, int(target_bytes // max(1, bytes_per_source)))
    if unique_sources.size <= block_size:
        return 0
    logging.info(
        "Auto source-block size resolved to %d sources (target %.2f GiB, %.2f MiB/source, total sources=%d).",
        block_size,
        float(KERNEL_SOURCE_BLOCK_TARGET_GIB),
        bytes_per_source / (1024 ** 2),
        unique_sources.size,
    )
    return block_size


def _kernel_ctx_nbytes(ctx: dict) -> int:
    total = 0
    for key in (
        "scaled_past",
        "scaled_next",
        "count_past",
        "count_next_past",
        "dense_past",
        "dense_next_past",
        "dense_past_packed",
        "dense_next_past_packed",
    ):
        arr = ctx.get(key)
        if arr is not None:
            total += int(getattr(arr, "nbytes", 0))
    return total


def ensure_gene_names(path='gene_names'):
    """
    Load gene names from a text file (one per line) if available.
    Caches the result globally to avoid repeated disk reads.
    """
    global GENE_NAMES
    if GENE_NAMES is not None:
        return GENE_NAMES
    try:
        with open(path, 'r') as fh:
            GENE_NAMES = [line.rstrip('\r\n') for line in fh]
        if GENE_NAMES:
            logging.info(f"Loaded {len(GENE_NAMES)} gene names from {path}.")
        else:
            logging.warning(f"Gene names file {path} is empty; outputs will use indices.")
    except FileNotFoundError:
        logging.warning(f"Gene names file {path} not found; outputs will use indices.")
        GENE_NAMES = []
    except Exception as e:
        logging.error(f"Error loading gene names from {path}: {e}")
        GENE_NAMES = []
    return GENE_NAMES


def gene_name_from_index(idx):
    """
    Translate a 1-based gene index into its name using the cached list.
    Falls back to the index string if lookup is unavailable.
    """
    names = ensure_gene_names()
    if not names:
        return str(int(idx)) if idx is not None else None
    try:
        return names[int(idx) - 1]
    except (IndexError, ValueError, TypeError):
        return str(int(idx)) if idx is not None else None


def _is_on(value) -> bool:
    return str(value).strip().lower() in ("on", "1", "true", "yes", "y")


def _direct_localte_keys(pair_mode: str, split_type: str) -> list[str]:
    """Return LocalTE split outputs that match Matrix_generate outputs."""
    pair_mode = str(pair_mode or "default").strip().lower()
    split_type = str(split_type or "1").strip()
    if pair_mode in ("gene_only", "all_gene"):
        return ["TE_GN_GN"]
    if pair_mode in ("all_feature", "all", "full", "gene_peak"):
        return ["TE_all_features"]
    if split_type == "1":
        return ["TE_TF_GN", "TE_TF_PK", "TE_PK_GN"]
    if split_type == "2":
        return ["TE_TF_GN"]
    if split_type == "3":
        return ["TE_TF_PK"]
    if split_type == "4":
        return ["TE_TF_GN", "TE_TF_PK"]
    if split_type == "5":
        return ["TE_PK_GN"]
    if split_type == "6":
        return ["TE_PK_PK"]
    return ["TE_TF_GN"]


def _read_name_set(path) -> set[str]:
    try:
        resolved = locate_file(str(path))
    except Exception:
        resolved = path
    try:
        with open(resolved, "r", encoding="utf-8") as handle:
            return {line.strip() for line in handle if line.strip()}
    except Exception:
        return set()


def _new_virtual_localte_state(out_dir: Path, value_type, compression: str, buffer_edges: int, value_dtype: str) -> dict:
    part_dir = Path(out_dir)
    part_dir.mkdir(parents=True, exist_ok=True)
    return {
        "part_dir": part_dir,
        "rows": [],
        "part_idx": 0,
        "edges": 0,
        "states": 0,
        "value_type": value_type,
        "compression": None if str(compression).lower() == "none" else str(compression),
        "buffer_edges": max(1, int(buffer_edges)),
        "value_dtype": str(value_dtype),
    }


def _write_virtual_localte_part(state: dict) -> None:
    rows = state.get("rows") or []
    if not rows:
        return
    value_type = state["value_type"]
    table = pa.table(
        {
            "Source": pa.array([r["Source"] for r in rows], type=pa.int64()),
            "Target": pa.array([r["Target"] for r in rows], type=pa.int64()),
            "TE": pa.array([r["TE"] for r in rows], type=pa.float64()),
            "NTargetStates": pa.array([r["NTargetStates"] for r in rows], type=pa.int32()),
            "NObs": pa.array([r["NObs"] for r in rows], type=pa.int32()),
            "NLocalStates": pa.array([r["NLocalStates"] for r in rows], type=pa.int32()),
            "StateCodes": pa.array([r["StateCodes"] for r in rows], type=pa.list_(pa.int64())),
            "LocalValues": pa.array([r["LocalValues"] for r in rows], type=pa.list_(value_type)),
        },
        schema=_virtual_localte_schema(value_type),
    )
    part_path = Path(state["part_dir"]) / f"part_{int(state['part_idx']):06d}.parquet"
    pq.write_table(table, part_path, compression=state["compression"])
    state["part_idx"] = int(state["part_idx"]) + 1
    state["rows"] = []


def _append_virtual_localte_row(
    state: dict,
    source: int,
    target: int,
    score: float,
    n_target_states: int,
    n_obs: int,
    state_codes,
    local_values,
) -> None:
    _code_type, code_np_dtype, _code_dtype_str = _virtual_state_code_type()
    if code_np_dtype is np.int32 and np.asarray(state_codes).size:
        max_code = int(np.max(state_codes))
        min_code = int(np.min(state_codes))
        if min_code < np.iinfo(np.int32).min or max_code > np.iinfo(np.int32).max:
            raise OverflowError(
                "StateCodes exceed int32 range; rerun with LOCAL_TE_VIRTUAL_STATE_CODE_BITS=64."
            )
    state_codes = np.asarray(state_codes, dtype=code_np_dtype)
    local_values = np.asarray(local_values)
    state["rows"].append(
        {
            "Source": int(source),
            "Target": int(target),
            "TE": float(score),
            "NTargetStates": int(n_target_states),
            "NObs": int(n_obs),
            "NLocalStates": int(local_values.size if int(n_target_states) < 0 else state_codes.size),
            "StateCodes": state_codes,
            "LocalValues": local_values,
        }
    )
    state["edges"] = int(state.get("edges", 0)) + 1
    state["states"] = int(state.get("states", 0)) + int(local_values.size if int(n_target_states) < 0 else state_codes.size)
    if len(state["rows"]) >= int(state["buffer_edges"]):
        _write_virtual_localte_part(state)


def _finalize_virtual_localte_state(state: dict) -> dict:
    _write_virtual_localte_part(state)
    return {
        "edges": int(state.get("edges", 0)),
        "states": int(state.get("states", 0)),
        "parts": int(state.get("part_idx", 0)),
    }


def _virtual_code_index_enabled() -> bool:
    return _is_on(LOCAL_TE_VIRTUAL) and _is_on(os.getenv("LOCAL_TE_VIRTUAL_CODE_INDEX", LOCAL_TE_VIRTUAL_CODE_INDEX))


def _virtual_code_index_feature_sets(pair_mode: str, list_pairs, pair_index_info) -> tuple[np.ndarray, np.ndarray]:
    if pair_mode == "default" and list_pairs is not None and len(list_pairs) > 0:
        arr = np.asarray(list_pairs, dtype=np.int64)
        return np.unique(arr[:, 0]), np.unique(arr[:, 1])
    if pair_index_info is not None and "indices" in pair_index_info:
        indices = np.asarray(pair_index_info["indices"], dtype=np.int64)
        indices = np.unique(indices)
        return indices, indices
    return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)


def _target_codes_for_code_index(target_idx: int, history_length: int, kernel_width: float) -> tuple[np.ndarray, int, int]:
    dest = _maybe_subsample(_get_series(int(target_idx)))
    if dest.size <= int(history_length):
        return np.zeros(0, dtype=np.int32), 0, 0
    dest_past = dest[:-1].astype(np.float64, copy=False)
    dest_next = dest[1:].astype(np.float64, copy=False)
    epsilon = 1e-9
    kw_past = float(kernel_width) * (np.std(dest_past, ddof=1) + epsilon)
    kw_next = float(kernel_width) * (np.std(dest_next, ddof=1) + epsilon)
    scaled_past = np.divide(dest_past, kw_past, out=np.zeros_like(dest_past), where=(kw_past != 0))
    scaled_next = np.divide(dest_next, kw_next, out=np.zeros_like(dest_next), where=(kw_next != 0))
    target_codes, target_past_values, _target_next_values, _target_counts = _target_grouped_state_codes(
        scaled_past,
        scaled_next,
    )
    return (
        np.asarray(target_codes, dtype=np.int32),
        int(target_past_values.size),
        int(target_codes.size),
    )


def _write_code_index_rows(path: Path, rows: list[dict], schema: pa.Schema, compression: str, writer) -> object:
    if not rows:
        return writer
    names = set(schema.names)
    payload = {
        "Feature": pa.array([r["Feature"] for r in rows], type=pa.int64()),
        "NObs": pa.array([r["NObs"] for r in rows], type=pa.int32()),
    }
    if "NSourceStates" in names:
        payload["NSourceStates"] = pa.array([r["NSourceStates"] for r in rows], type=pa.int32())
        payload["SourceCodes"] = pa.array([r["SourceCodes"] for r in rows], type=pa.list_(pa.int32()))
    else:
        payload["NTargetStates"] = pa.array([r["NTargetStates"] for r in rows], type=pa.int32())
        payload["TargetCodes"] = pa.array([r["TargetCodes"] for r in rows], type=pa.list_(pa.int32()))
    table = pa.table(payload, schema=schema)
    if writer is None:
        writer = pq.ParquetWriter(
            path,
            schema,
            compression=None if str(compression).lower() == "none" else str(compression),
        )
    writer.write_table(table)
    rows.clear()
    return writer


def _write_virtual_code_index(
    source_indices,
    target_indices,
    history_length: int,
    kernel_width: float = 0.5,
    output_dir: str | None = None,
) -> dict | None:
    """Write reusable source/target code sequences for virtual LocalTE export.

    The virtual LocalTE table stores edge-specific state values.  Exporting full
    timepoint arrays also needs source and target code sequences; storing them
    once per feature avoids recomputing the same sequences for every edge.
    """
    if not _virtual_code_index_enabled():
        return None
    if pa is None or pq is None:
        logging.warning("Virtual LocalTE code index skipped because pyarrow is unavailable.")
        return None

    sources = np.asarray(source_indices, dtype=np.int64)
    targets = np.asarray(target_indices, dtype=np.int64)
    sources = np.unique(sources[sources > 0])
    targets = np.unique(targets[targets > 0])
    if sources.size == 0 and targets.size == 0:
        return None

    out_dir = Path(output_dir or os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", LOCAL_TE_CODE_INDEX_OUTPUT_DIR))
    compression = os.getenv("LOCAL_TE_CODE_INDEX_COMPRESSION", LOCAL_TE_CODE_INDEX_COMPRESSION).strip().lower()
    batch_features = max(1, int(os.getenv("LOCAL_TE_CODE_INDEX_BATCH_FEATURES", str(LOCAL_TE_CODE_INDEX_BATCH_FEATURES))))
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in ("source_codes.parquet", "target_codes.parquet", "metadata.json"):
        try:
            (out_dir / stale).unlink()
        except FileNotFoundError:
            pass

    source_schema = pa.schema(
        [
            ("Feature", pa.int64()),
            ("NObs", pa.int32()),
            ("NSourceStates", pa.int32()),
            ("SourceCodes", pa.list_(pa.int32())),
        ]
    )
    target_schema = pa.schema(
        [
            ("Feature", pa.int64()),
            ("NObs", pa.int32()),
            ("NTargetStates", pa.int32()),
            ("TargetCodes", pa.list_(pa.int32())),
        ]
    )
    start = time.time()
    print(
        f"[LocalTE] Building virtual code index: sources={_human_int(sources.size)}, "
        f"targets={_human_int(targets.size)}, out={out_dir}"
    )

    source_writer = None
    source_rows: list[dict] = []
    try:
        for feature in sources:
            source_values, source_codes = _get_grouped_source_codes_kernel(
                int(feature),
                int(history_length),
                float(kernel_width),
            )
            source_rows.append(
                {
                    "Feature": int(feature),
                    "NObs": int(np.asarray(source_codes).size),
                    "NSourceStates": int(np.asarray(source_values).size),
                    "SourceCodes": np.asarray(source_codes, dtype=np.int32),
                }
            )
            if len(source_rows) >= batch_features:
                source_writer = _write_code_index_rows(
                    out_dir / "source_codes.parquet",
                    source_rows,
                    source_schema,
                    compression,
                    source_writer,
                )
        source_writer = _write_code_index_rows(
            out_dir / "source_codes.parquet",
            source_rows,
            source_schema,
            compression,
            source_writer,
        )
    finally:
        if source_writer is not None:
            source_writer.close()

    target_writer = None
    target_rows: list[dict] = []
    try:
        for feature in targets:
            target_codes, n_target_states, n_obs = _target_codes_for_code_index(
                int(feature),
                int(history_length),
                float(kernel_width),
            )
            target_rows.append(
                {
                    "Feature": int(feature),
                    "NObs": int(n_obs),
                    "NTargetStates": int(n_target_states),
                    "TargetCodes": np.asarray(target_codes, dtype=np.int32),
                }
            )
            if len(target_rows) >= batch_features:
                target_writer = _write_code_index_rows(
                    out_dir / "target_codes.parquet",
                    target_rows,
                    target_schema,
                    compression,
                    target_writer,
                )
        target_writer = _write_code_index_rows(
            out_dir / "target_codes.parquet",
            target_rows,
            target_schema,
            compression,
            target_writer,
        )
    finally:
        if target_writer is not None:
            target_writer.close()

    elapsed = time.time() - start
    payload = {
        "format": "tenet_localte_code_index_v1",
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "history_length": int(history_length),
        "kernel_width": float(kernel_width),
        "sources": int(sources.size),
        "targets": int(targets.size),
        "compression": str(compression),
        "batch_features": int(batch_features),
        "elapsed_seconds": float(elapsed),
    }
    (out_dir / "metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[LocalTE] Code index ready in {elapsed:.2f} seconds: {out_dir}")
    return payload


def _build_direct_localte_exporter(args, pair_mode: str, total_pairs: int):
    """Prepare state for direct LocalTE time-chunk writing."""
    if _is_on(LOCAL_TE_VIRTUAL):
        if getattr(args, "permute", False):
            return None
        if pa is None or pq is None:
            logging.warning("Virtual LocalTE disabled because pyarrow is unavailable.")
            return None
        out_root = os.getenv("LOCAL_TE_VIRTUAL_OUTPUT_DIR", LOCAL_TE_VIRTUAL_OUTPUT_DIR)
        value_type, value_np_dtype, value_dtype_str = _virtual_localte_value_type(
            os.getenv("LOCAL_TE_VIRTUAL_VALUE_DTYPE", LOCAL_TE_VIRTUAL_VALUE_DTYPE)
        )
        buffer_edges = max(1, int(os.getenv("LOCAL_TE_VIRTUAL_BUFFER_EDGES", str(LOCAL_TE_VIRTUAL_BUFFER_EDGES))))
        compression = os.getenv("LOCAL_TE_VIRTUAL_COMPRESSION", LOCAL_TE_VIRTUAL_COMPRESSION).strip().lower()
        out_path = Path(out_root)
        if (out_path / "parts").exists():
            shutil.rmtree(out_path / "parts")
        for stale in ("metadata.json", "direct_virtual_localte_done.json"):
            try:
                (out_path / stale).unlink()
            except FileNotFoundError:
                pass
        out_path.mkdir(parents=True, exist_ok=True)
        print(
            f"[LocalTE] Virtual exact array: out={out_root}, dtype={value_dtype_str}, "
            f"buffer_edges={buffer_edges}, compression={compression}"
        )
        return {
            "virtual_array": True,
            "out_root": out_root,
            "keys": ["ALL"],
            "specs": [("ALL", out_root)],
            "states": [],
            "value_type": value_type,
            "value_np_dtype": value_np_dtype,
            "value_dtype": value_dtype_str,
            "compression": compression,
            "buffer_edges": buffer_edges,
            "routed_counts": [0],
            "total_pairs": int(total_pairs),
            "args": args,
        }
    if not STORE_LOCAL_TE:
        return None
    if not _is_on(LOCAL_TE_DIRECT_CHUNKS):
        return None
    if getattr(args, "permute", False):
        return None
    if not _is_on(os.getenv("LOCAL_TE_SPLIT_EXPORT", "on")):
        return None
    use_video_tiles = _is_on(LOCAL_TE_VIDEO_TILES)
    if any(
        helper is None
        for helper in (
            _append_localte_to_state,
            _decode_local_te,
            _finalize_export_state,
            _new_export_state,
            _route_known_split_indices,
        )
    ):
        logging.warning("Direct LocalTE chunks disabled because chunk-export helpers are unavailable.")
        return None
    if use_video_tiles and any(
        helper is None
        for helper in (
            append_localte_to_video_state,
            finalize_video_export_state,
            new_video_export_state,
            warmup_video_factorized_encoder,
        )
    ):
        logging.warning("Direct LocalTE video tiles disabled because video helpers are unavailable.")
        return None
    if pa is None or pq is None:
        logging.warning("Direct LocalTE chunks disabled because pyarrow is unavailable.")
        return None

    names = ensure_gene_names()
    if not names:
        logging.warning("Direct LocalTE chunks disabled because gene_names is unavailable.")
        return None

    split_type = str(getattr(args, "split_type", "1") or "1")
    species = str(getattr(args, "species", "") or "").strip()
    keys = _direct_localte_keys(pair_mode, split_type)
    if not keys:
        return None

    name_to_idx = {name: idx for idx, name in enumerate(names, start=1)}
    is_peak = np.zeros(len(names) + 1, dtype=bool)
    peak_names = load_peak_names(names, locate_file("TE_peak_list.txt"))
    for idx, name in enumerate(names, start=1):
        if name in peak_names:
            is_peak[idx] = True

    tf_sources: set[int] = set()
    if any(key.startswith("TE_TF_") for key in keys):
        from code.species_support import load_tf_names
        tf_names = set(load_tf_names(species))
        tf_sources = {name_to_idx[name] for name in tf_names if name in name_to_idx}
        if not tf_sources:
            logging.warning("Direct LocalTE chunks disabled because TF source list is empty.")
            return None

    if use_video_tiles:
        out_root = os.getenv("LOCAL_TE_VIDEO_OUTPUT_DIR", "local_te_video_tiles")
    else:
        out_root = os.getenv("LOCAL_TE_SPLIT_OUTPUT_DIR", "local_te_split_chunks")
    chunk_size = max(1, int(os.getenv("LOCAL_TE_CHUNK_SIZE", "300")))
    direct_buffer_default = os.getenv("LOCAL_TE_DIRECT_BUFFER_EDGES", "").strip()
    if not direct_buffer_default or direct_buffer_default.lower() == "auto":
        direct_buffer_default = (
            "5000"
            if use_video_tiles and _is_on(LOCAL_TE_DIRECT_WORKERS)
            else ("500" if _is_on(LOCAL_TE_DIRECT_WORKERS) else os.getenv("LOCAL_TE_BUFFER_EDGES", "10000"))
        )
    buffer_edges = max(1, int(direct_buffer_default))
    values_dtype_str = os.getenv("LOCAL_TE_VALUES_DTYPE", "float16").strip().lower()
    if values_dtype_str not in ("float16", "float32"):
        values_dtype_str = "float16"
    values_dtype = np.float16 if values_dtype_str == "float16" else np.float32
    value_type = pa.float16() if values_dtype is np.float16 else pa.float32()
    compression = os.getenv("LOCAL_TE_CHUNK_COMPRESSION", "snappy").strip().lower()
    if compression == "none":
        compression = None
    zstd_level_env = os.getenv("LOCAL_TE_ZSTD_LEVEL", "")
    try:
        zstd_level = int(zstd_level_env) if zstd_level_env else None
    except Exception:
        zstd_level = None
    use_dictionary = not _is_on(os.getenv("LOCAL_TE_NO_DICTIONARY", "off"))
    disable_chunk_scores = _is_on(os.getenv("LOCAL_TE_DISABLE_CHUNK_SCORES", "off"))

    export_args = argparse.Namespace(
        disable_scores=disable_chunk_scores,
        compression=compression or "none",
        zstd_level=zstd_level,
        no_dictionary=(not use_dictionary),
        buffer_edges=buffer_edges,
        values_dtype=values_dtype_str,
        use_threads=os.getenv("LOCAL_TE_USE_THREADS", "on"),
    )
    states = []
    specs = []
    for key in keys:
        out_dir = os.path.join(out_root, key)
        if use_video_tiles:
            states.append(
                new_video_export_state(
                    Path(out_dir),
                    chunk_size,
                    LOCAL_TE_VIDEO_LEVEL,
                    values_dtype_str,
                    LOCAL_TE_VIDEO_EDGE_FORMAT,
                    LOCAL_TE_VIDEO_CODEC,
                    _is_on(LOCAL_TE_VIDEO_FACTORIZE),
                )
            )
        else:
            states.append(_new_export_state(Path(out_dir)))
        specs.append((key, out_dir))

    scores_schema = pa.schema(
        [
            ("Source", pa.int32()),
            ("Target", pa.int32()),
            ("TE", pa.float32()),
            ("LocalTE_total_len", pa.int32()),
        ]
    )
    Path(out_root).mkdir(parents=True, exist_ok=True)
    if use_video_tiles and _is_on(LOCAL_TE_VIDEO_FACTORIZE):
        warmup_video_factorized_encoder(chunk_size)
    if use_video_tiles:
        print(
            f"[LocalTE] Direct video-tile streaming: {len(keys)} outputs -> {out_root} "
            f"(chunk={chunk_size}, dtype={values_dtype_str}, buffer_edges={buffer_edges}, "
            f"codec={LOCAL_TE_VIDEO_CODEC}-{LOCAL_TE_VIDEO_LEVEL}, edges={LOCAL_TE_VIDEO_EDGE_FORMAT}, "
            f"factorize={LOCAL_TE_VIDEO_FACTORIZE})"
        )
    else:
        print(
            f"[LocalTE] Direct chunk streaming: {len(keys)} outputs -> {out_root} "
            f"(chunk={chunk_size}, dtype={values_dtype_str}, buffer_edges={buffer_edges})"
        )
    logging.info(
        "Direct LocalTE %s enabled: keys=%s, total_pairs=%d, out=%s.",
        "video tiles" if use_video_tiles else "chunks",
        ",".join(keys),
        int(total_pairs),
        out_root,
    )
    return {
        "out_root": out_root,
        "keys": keys,
        "specs": specs,
        "states": states,
        "is_peak": is_peak,
        "tf_sources": tf_sources,
        "chunk_size": chunk_size,
        "buffer_edges": buffer_edges,
        "values_dtype": values_dtype,
        "value_type": value_type,
        "compression": compression or "none",
        "zstd_level": zstd_level,
        "use_dictionary": use_dictionary,
        "video_tiles": use_video_tiles,
        "video_level": int(LOCAL_TE_VIDEO_LEVEL),
        "video_edge_format": str(LOCAL_TE_VIDEO_EDGE_FORMAT),
        "video_codec": str(LOCAL_TE_VIDEO_CODEC),
        "video_factorize": bool(_is_on(LOCAL_TE_VIDEO_FACTORIZE)),
        "scores_schema": scores_schema,
        "args": export_args,
        "routed_counts": [0 for _ in keys],
    }


def _direct_localte_append_rows(exporter: dict | None, rows) -> list[tuple[int, int, float]]:
    """Route LocalTE payloads to direct chunk writers and return score-only rows."""
    if exporter is None:
        return rows
    score_rows = []
    keys = exporter["keys"]
    for row in rows:
        if exporter.get("virtual_array") and len(row) >= 8 and row[3] == VIRTUAL_LOCAL_TE_MARKER:
            src_i = int(row[0])
            tgt_i = int(row[1])
            score = float(row[2])
            score_rows.append((src_i, tgt_i, score))
            _append_virtual_localte_row(
                exporter["states"][0],
                src_i,
                tgt_i,
                score,
                int(row[4]),
                int(row[5]),
                row[6],
                row[7],
            )
            exporter["routed_counts"][0] += 1
            continue
        if len(row) < 7:
            score_rows.append((int(row[0]), int(row[1]), float(row[2])))
            continue
        src_i = int(row[0])
        tgt_i = int(row[1])
        score = float(row[2])
        score_rows.append((src_i, tgt_i, score))
        blob, length, dtype_str, codec = row[3], row[4], row[5], row[6]
        if blob is None or length is None or int(length) == 0:
            continue
        out_indices = _route_known_split_indices(
            src_i,
            tgt_i,
            exporter["is_peak"],
            exporter["tf_sources"],
            keys,
        )
        if not out_indices:
            continue
        values = _decode_local_te(blob, dtype_str, int(length), codec)
        if values.size == 0:
            continue
        if values.dtype != exporter["values_dtype"]:
            values = values.astype(exporter["values_dtype"], copy=False)
        for out_idx in out_indices:
            exporter["routed_counts"][out_idx] += 1
            if exporter.get("video_tiles"):
                append_localte_to_video_state(
                    exporter["states"][out_idx],
                    src_i,
                    tgt_i,
                    score,
                    values,
                    int(length),
                    exporter["chunk_size"],
                    exporter["buffer_edges"],
                )
            else:
                _append_localte_to_state(
                    exporter["states"][out_idx],
                    src_i,
                    tgt_i,
                    score,
                    values,
                    int(length),
                    exporter["chunk_size"],
                    exporter["buffer_edges"],
                    exporter["value_type"],
                    exporter["compression"],
                    exporter["zstd_level"],
                    exporter["use_dictionary"],
                    bool(exporter["args"].disable_scores),
                    exporter["scores_schema"],
                )
    return score_rows


def _direct_localte_finalize(exporter: dict | None) -> None:
    if exporter is None:
        return
    if exporter.get("virtual_array"):
        summaries = [_finalize_virtual_localte_state(state) for state in exporter["states"]]
        total_edges = sum(int(s["edges"]) for s in summaries)
        total_states = sum(int(s["states"]) for s in summaries)
        out_root = Path(exporter["out_root"])
        payload = {
            "format": "tenet_localte_virtual_array_v1",
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "recipe": os.path.abspath(os.path.join(os.getenv("LOCAL_TE_RECIPE_DIR", LOCAL_TE_RECIPE_DIR), "recipe.json")),
            "edges": os.path.abspath("TE_result_all.parquet"),
            "matrix_parquet": os.path.abspath("cell_gene_trsps.parquet"),
            "gene_names": os.path.abspath("gene_names"),
            "time_index_map": os.path.abspath("time_index_map.parquet") if os.path.exists("time_index_map.parquet") else None,
            "all_pairs": os.path.abspath(getattr(exporter.get("args"), "input_file", "all_pairs.csv")),
            "history_length": int(getattr(exporter.get("args"), "history_length", 1)),
            "kernel_width": 0.5,
            "code_index": os.path.abspath(os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", LOCAL_TE_CODE_INDEX_OUTPUT_DIR))
            if os.path.exists(os.path.join(os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", LOCAL_TE_CODE_INDEX_OUTPUT_DIR), "metadata.json"))
            else None,
            "parts_dir": "parts",
            "edges_stored": int(total_edges),
            "state_values_stored": int(total_states),
            "mean_states_per_edge": float(total_states / total_edges) if total_edges else 0.0,
            "value_dtype": str(exporter["value_dtype"]),
            "compression": str(exporter["compression"]),
            "buffer_edges": int(exporter["buffer_edges"]),
            "worker_results": summaries,
            "note": "Virtual LocalTE written during kernel TE. Dense fallback edges use NTargetStates=-1.",
        }
        (out_root / "metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        (out_root / "direct_virtual_localte_done.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(
            f"[LocalTE] Virtual array completed: {total_edges:,} edges, "
            f"{total_states:,} stored values -> {out_root}"
        )
        return
    for state in exporter["states"]:
        if exporter.get("video_tiles"):
            finalize_video_export_state(state)
        else:
            _finalize_export_state(
                state,
                exporter["args"],
                exporter["chunk_size"],
                exporter["value_type"],
            )
    marker = Path(exporter["out_root"]) / "direct_localte_done.json"
    payload = {
        "mode": "direct_video_tile_stream" if exporter.get("video_tiles") else "direct_chunk_stream",
        "format": "localte_video_tiles_v1" if exporter.get("video_tiles") else "localte_chunk_parquet_v1",
        "keys": exporter["keys"],
        "outputs": [
            {"key": key, "out_dir": out_dir, "edges": int(count)}
            for (key, out_dir), count in zip(exporter["specs"], exporter["routed_counts"])
        ],
        "chunk_size": int(exporter["chunk_size"]),
        "values_dtype": str(exporter["args"].values_dtype),
        "compression": str(exporter.get("video_codec", "zstd")) if exporter.get("video_tiles") else str(exporter["compression"]),
        "video_tiles": bool(exporter.get("video_tiles")),
        "video_level": int(exporter.get("video_level", 0)),
        "video_edge_format": str(exporter.get("video_edge_format", "")),
        "video_codec": str(exporter.get("video_codec", "")),
        "video_factorize": bool(exporter.get("video_factorize", False)),
        "chunk_scores": False if exporter.get("video_tiles") else (not bool(exporter["args"].disable_scores)),
    }
    marker.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    for item in payload["outputs"]:
        print(f"[LocalTE] Direct {item['key']}: {item['edges']:,} edges -> {item['out_dir']}")
    label = "video-tile" if exporter.get("video_tiles") else "chunk"
    print(f"[LocalTE] Direct {label} streaming completed: {marker}")


def _direct_localte_make_worker_exporter(exporter: dict, worker_id: int) -> dict:
    """Create per-worker LocalTE chunk states under parts_wXX directories."""
    worker = dict(exporter)
    worker_specs = []
    worker_states = []
    for key, out_dir in exporter["specs"]:
        if exporter.get("virtual_array"):
            part_dir = Path(out_dir) / "parts" / f"worker_{int(worker_id):03d}"
            worker_specs.append((key, str(part_dir)))
            worker_states.append(
                _new_virtual_localte_state(
                    part_dir,
                    exporter["value_type"],
                    exporter["compression"],
                    exporter["buffer_edges"],
                    exporter["value_dtype"],
                )
            )
            continue
        part_dir = Path(out_dir) / f"parts_w{int(worker_id):02d}"
        worker_specs.append((key, str(part_dir)))
        if exporter.get("video_tiles"):
            worker_states.append(
                new_video_export_state(
                    part_dir,
                    exporter["chunk_size"],
                    exporter.get("video_level", 1),
                    str(exporter["args"].values_dtype),
                    str(exporter.get("video_edge_format", "binary")),
                    str(exporter.get("video_codec", "zstd")),
                    bool(exporter.get("video_factorize", False)),
                )
            )
        else:
            worker_states.append(_new_export_state(part_dir))
    worker["specs"] = worker_specs
    worker["states"] = worker_states
    worker["routed_counts"] = [0 for _ in exporter["keys"]]
    return worker


def _direct_localte_finalize_worker(exporter: dict | None) -> list[int]:
    """Finalize one worker's chunk writers and return routed counts."""
    if exporter is None:
        return []
    for state in exporter["states"]:
        if exporter.get("virtual_array"):
            _finalize_virtual_localte_state(state)
            continue
        if exporter.get("video_tiles"):
            finalize_video_export_state(state)
        else:
            _finalize_export_state(
                state,
                exporter["args"],
                exporter["chunk_size"],
                exporter["value_type"],
            )
    return [int(x) for x in exporter.get("routed_counts", [])]


def _direct_localte_write_marker_from_worker_counts(exporter: dict | None, worker_counts: list[list[int]]) -> None:
    """Write the direct LocalTE completion marker for worker-owned part outputs."""
    if exporter is None:
        return
    totals = [0 for _ in exporter["keys"]]
    for counts in worker_counts:
        for i, value in enumerate(counts):
            if i < len(totals):
                totals[i] += int(value)
    marker = Path(exporter["out_root"]) / "direct_localte_done.json"
    if exporter.get("virtual_array"):
        out_root = Path(exporter["out_root"])
        parts = sorted((out_root / "parts").rglob("part_*.parquet"))
        total_states = 0
        for part in parts:
            try:
                table = pq.read_table(str(part), columns=["NLocalStates"])
                total_states += int(table.column("NLocalStates").combine_chunks().to_numpy().sum())
            except Exception:
                pass
        total_edges = int(totals[0]) if totals else 0
        payload = {
            "format": "tenet_localte_virtual_array_v1",
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "recipe": os.path.abspath(os.path.join(os.getenv("LOCAL_TE_RECIPE_DIR", LOCAL_TE_RECIPE_DIR), "recipe.json")),
            "edges": os.path.abspath("TE_result_all.parquet"),
            "matrix_parquet": os.path.abspath("cell_gene_trsps.parquet"),
            "gene_names": os.path.abspath("gene_names"),
            "time_index_map": os.path.abspath("time_index_map.parquet") if os.path.exists("time_index_map.parquet") else None,
            "all_pairs": os.path.abspath(getattr(exporter.get("args"), "input_file", "all_pairs.csv")),
            "history_length": int(getattr(exporter.get("args"), "history_length", 1)),
            "kernel_width": 0.5,
            "code_index": os.path.abspath(os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", LOCAL_TE_CODE_INDEX_OUTPUT_DIR))
            if os.path.exists(os.path.join(os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", LOCAL_TE_CODE_INDEX_OUTPUT_DIR), "metadata.json"))
            else None,
            "parts_dir": "parts",
            "edges_stored": total_edges,
            "state_values_stored": int(total_states),
            "mean_states_per_edge": float(total_states / total_edges) if total_edges else 0.0,
            "value_dtype": str(exporter["value_dtype"]),
            "compression": str(exporter["compression"]),
            "buffer_edges": int(exporter["buffer_edges"]),
            "workers": len(worker_counts),
            "worker_counts": worker_counts,
            "note": "Virtual LocalTE written during kernel TE. Dense fallback edges use NTargetStates=-1.",
        }
        (out_root / "metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        (out_root / "direct_virtual_localte_done.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(
            f"[LocalTE] Virtual array completed: {total_edges:,} edges, "
            f"{int(total_states):,} stored values -> {out_root}"
        )
        return
    payload = {
        "mode": (
            "direct_video_tile_stream_worker_parts"
            if exporter.get("video_tiles")
            else "direct_chunk_stream_worker_parts"
        ),
        "format": "localte_video_tiles_v1" if exporter.get("video_tiles") else "localte_chunk_parquet_v1",
        "keys": exporter["keys"],
        "outputs": [
            {"key": key, "out_dir": out_dir, "edges": int(count)}
            for (key, out_dir), count in zip(exporter["specs"], totals)
        ],
        "chunk_size": int(exporter["chunk_size"]),
        "values_dtype": str(exporter["args"].values_dtype),
        "compression": str(exporter.get("video_codec", "zstd")) if exporter.get("video_tiles") else str(exporter["compression"]),
        "video_tiles": bool(exporter.get("video_tiles")),
        "video_level": int(exporter.get("video_level", 0)),
        "video_edge_format": str(exporter.get("video_edge_format", "")),
        "video_codec": str(exporter.get("video_codec", "")),
        "video_factorize": bool(exporter.get("video_factorize", False)),
        "merge_parts": False,
        "chunk_scores": False if exporter.get("video_tiles") else (not bool(exporter["args"].disable_scores)),
    }
    marker.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    for item in payload["outputs"]:
        print(f"[LocalTE] Direct {item['key']}: {item['edges']:,} edges -> {item['out_dir']}")
    label = "video-tile" if exporter.get("video_tiles") else "chunk"
    print(f"[LocalTE] Direct worker {label} streaming completed: {marker}")


def _direct_localte_worker_loop(
    worker_id: int,
    task_queue,
    result_queue,
    exporter: dict,
    score_parts_dir: str | None = None,
    score_buffer_rows: int = 5000,
) -> None:
    """Compute grouped-state target tasks and write LocalTE chunks inside one worker."""
    worker_exporter = _direct_localte_make_worker_exporter(exporter, worker_id)
    score_writer_state = {"writer": None, "rows": 0}
    score_buffer = []
    score_columns = ["Source", "Target", "TE"]
    score_output = None
    if score_parts_dir:
        score_output = str(Path(score_parts_dir) / f"part_w{int(worker_id):02d}.parquet")

    def flush_scores():
        nonlocal score_buffer
        if not score_output or not score_buffer:
            return
        _append_grouped_rows_to_parquet_writer(
            score_buffer,
            score_columns,
            score_writer_state,
            score_output,
        )
        score_buffer = []

    try:
        while True:
            task = task_queue.get()
            if task is None:
                break
            target_rows = _process_grouped_state_target(task)
            score_rows = _direct_localte_append_rows(worker_exporter, target_rows)
            if score_output:
                score_buffer.extend(score_rows)
                if len(score_buffer) >= max(1, int(score_buffer_rows)):
                    flush_scores()
                result_queue.put(("progress", int(len(score_rows))))
            else:
                result_queue.put(("rows", int(len(score_rows)), score_rows))
    except BaseException as exc:
        result_queue.put(("error", int(worker_id), repr(exc)))
    finally:
        try:
            flush_scores()
        except BaseException as exc:
            result_queue.put(("error", int(worker_id), f"score finalize failed: {exc!r}"))
        finally:
            if score_writer_state.get("writer") is not None:
                score_writer_state["writer"].close()
                score_writer_state["writer"] = None
        try:
            counts = _direct_localte_finalize_worker(worker_exporter)
        except BaseException as exc:
            result_queue.put(("error", int(worker_id), f"finalize failed: {exc!r}"))
            counts = [0 for _ in exporter.get("keys", [])]
        result_queue.put(("done", int(worker_id), counts))


def _run_grouped_state_direct_worker_writers(
    work_units,
    total_pairs: int,
    num_cpus: int,
    output_file: str,
    direct_localte_exporter: dict,
    result_buffer_rows: int,
) -> bool:
    """Run grouped-state TE with worker-owned LocalTE chunk writers.

    This avoids sending large LocalTE byte payloads back to the parent process.
    Workers return only score rows; each worker writes its own chunk partitions.
    """
    if direct_localte_exporter is None:
        return False
    n_workers = max(1, min(int(num_cpus), int(len(work_units))))
    use_score_parts = _is_on(LOCAL_TE_DIRECT_SCORE_PARTS)
    score_parts_dir = None
    if use_score_parts:
        score_parts_dir = str(Path(output_file).with_suffix("").parent / "TE_score_parts")
        Path(score_parts_dir).mkdir(parents=True, exist_ok=True)
        for stale_part in Path(score_parts_dir).glob("part_w*.parquet"):
            try:
                stale_part.unlink()
            except Exception:
                pass
    if direct_localte_exporter.get("virtual_array"):
        localte_format = "virtual array parts"
    else:
        localte_format = "video-tile parts" if direct_localte_exporter.get("video_tiles") else "chunk parts"
    print(
        f"[TE] LocalTE direct worker writers: {n_workers} workers write {localte_format}; "
        + (
            "workers also write score parts."
            if use_score_parts
            else "main process writes score-only TE_result_all.parquet."
        )
    )
    ctx = multiprocessing.get_context("fork")
    task_queue = ctx.Queue()
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=_direct_localte_worker_loop,
            args=(
                worker_id,
                task_queue,
                result_queue,
                direct_localte_exporter,
                score_parts_dir,
                int(result_buffer_rows or 5000),
            ),
        )
        for worker_id in range(n_workers)
    ]

    for proc in processes:
        proc.start()
    for task in work_units:
        task_queue.put(task)
    for _ in processes:
        task_queue.put(None)

    direct_writer_state = {"writer": None, "rows": 0}
    columns = ["Source", "Target", "TE"]
    row_buffer = []
    buffered_count = 0
    worker_counts: list[list[int]] = []
    done_workers = 0
    processed_pairs = 0
    start_time = time.time()
    show_progress = str(KERNEL_GROUPED_STATE_PROGRESS).lower() not in ("off", "0", "false", "no")
    pair_pbar = tqdm(
        total=int(total_pairs),
        desc="TE pairs (compute+write)",
        unit="pair",
        dynamic_ncols=True,
        mininterval=1.0,
        leave=True,
    ) if show_progress else None

    def flush_score_rows(buffered_rows):
        if not buffered_rows:
            return []
        _append_grouped_rows_to_parquet_writer(
            buffered_rows,
            columns,
            direct_writer_state,
            output_file,
        )
        return []

    error_messages = []
    try:
        while done_workers < n_workers:
            msg = result_queue.get()
            tag = msg[0]
            if tag == "rows":
                n_rows = int(msg[1])
                score_rows = msg[2]
                processed_pairs += n_rows
                if pair_pbar is not None and n_rows:
                    pair_pbar.update(n_rows)
                if score_rows:
                    row_buffer.extend(score_rows)
                    buffered_count += n_rows
                    if buffered_count >= max(1, int(result_buffer_rows or 5000)):
                        row_buffer = flush_score_rows(row_buffer)
                        buffered_count = 0
            elif tag == "progress":
                n_rows = int(msg[1])
                processed_pairs += n_rows
                if pair_pbar is not None and n_rows:
                    pair_pbar.update(n_rows)
            elif tag == "done":
                done_workers += 1
                worker_counts.append([int(x) for x in msg[2]])
            elif tag == "error":
                error_messages.append(f"worker {msg[1]}: {msg[2]}")
        if row_buffer:
            flush_score_rows(row_buffer)
    finally:
        if pair_pbar is not None:
            pair_pbar.close()
        for proc in processes:
            proc.join(timeout=2)
        for proc in processes:
            if proc.is_alive():
                proc.terminate()
                proc.join(timeout=2)
        if direct_writer_state.get("writer") is not None:
            direct_writer_state["writer"].close()
            direct_writer_state["writer"] = None

    if error_messages:
        logging.error("Direct LocalTE worker writer failed: %s", "; ".join(error_messages[:5]))
        return False
    if use_score_parts:
        score_parts = sorted(str(p) for p in Path(score_parts_dir).glob("part_w*.parquet"))
        if not score_parts and int(total_pairs) > 0:
            return False
        if not concat_parquet_files_duckdb(score_parts, output_file):
            return False
    elif int(direct_writer_state.get("rows", 0)) <= 0 and int(total_pairs) > 0:
        return False
    _direct_localte_write_marker_from_worker_counts(direct_localte_exporter, worker_counts)
    elapsed = time.time() - start_time
    logging.info(
        "Grouped-state direct-worker LocalTE results saved to %s in %.2f seconds (%d pairs).",
        output_file,
        elapsed,
        processed_pairs,
    )
    return True


def _get_dest_ctx_kernel(target_idx: int, historyLength: int, kernel_width: float):
    global _DEST_CTX_CACHE_KERNEL_BYTES
    series_len = GLOBAL_CELL_GENE_DENSE.shape[1] if GLOBAL_CELL_GENE_DENSE is not None else GLOBAL_CELL_GENE.shape[1]
    backend = _select_kernel_backend(series_len, historyLength)
    key = (target_idx, historyLength, float(kernel_width), backend)
    ctx = _DEST_CTX_CACHE_KERNEL.get(key)
    if ctx is not None:
        _DEST_CTX_CACHE_KERNEL.move_to_end(key)
        return ctx
    dest = _get_series(target_idx)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context(
        dest,
        k=historyLength,
        kernel_width=kernel_width,
        normalise=True,
        precompute_neighbors=False,
        backend=backend,
        workers=1,
    )
    _DEST_CTX_CACHE_KERNEL[key] = ctx
    _DEST_CTX_CACHE_KERNEL.move_to_end(key)
    _DEST_CTX_CACHE_KERNEL_BYTES += _kernel_ctx_nbytes(ctx)
    while (
        len(_DEST_CTX_CACHE_KERNEL) > max(1, int(KERNEL_DEST_CACHE_MAX))
        or (
            KERNEL_DEST_CACHE_BYTES > 0
            and _DEST_CTX_CACHE_KERNEL_BYTES > KERNEL_DEST_CACHE_BYTES
        )
    ):
        _, old = _DEST_CTX_CACHE_KERNEL.popitem(last=False)
        _DEST_CTX_CACHE_KERNEL_BYTES -= _kernel_ctx_nbytes(old)
    return ctx


def _encode_local_te(local_array) -> tuple[bytes, int, str, str]:
    """Encode LocalTE array to bytes with optional compression.
    Returns (bytes, length, dtype_str, codec).
    """
    arr = np.asarray(local_array, dtype=LOCAL_TE_DTYPE, order='C')
    raw = arr.tobytes()
    codec = LOCAL_TE_CODEC if LOCAL_TE_CODEC in ('zlib', 'none') else 'zlib'
    if codec == 'zlib':
        try:
            import zlib as _z
            level = int(os.getenv('TE_LOCALTE_ZLIB_LEVEL', '3'))
            comp = _z.compress(raw, level=level)
            return comp, int(arr.size), LOCAL_TE_DTYPE_STR, 'zlib'
        except Exception:
            # Fallback to raw bytes on any compression error
            return raw, int(arr.size), LOCAL_TE_DTYPE_STR, 'none'
    return raw, int(arr.size), LOCAL_TE_DTYPE_STR, 'none'

def _process_chunk(args):
    """Worker: compute exact kernel TE for a single target and a chunk of sources."""
    target_idx, source_indices, historyLength, kernel_width, mode = args
    if mode != 'kernel':
        raise ValueError("This TENETPLUS runner only supports kernel mode.")
    results = []
    try:
        dest_ctx = _get_dest_ctx_kernel(target_idx, historyLength, kernel_width)
    except Exception as e:
        logging.error(f"Failed to build kernel context for target {target_idx}: {e}")
        return results
    if dest_ctx.get('total_obs', 0) <= 0:
        return results

    # Packed source-neighborhood masks scale with cells^2.  Process sources in
    # memory-bounded sub-batches while reusing the destination context.
    pack_subbatch = _kernel_source_pack_subbatch_size(dest_ctx, len(source_indices))
    for start in range(0, len(source_indices), pack_subbatch):
        source_chunk = source_indices[start : start + pack_subbatch]
        src_arrays = []
        scaled_src_arrays = []
        source_ok_packed_arrays = [] if dest_ctx.get('backend') == 'dense' else None
        valid_sources = []
        for sidx in source_chunk:
            try:
                scaled = _get_scaled_source_kernel(sidx, historyLength, kernel_width)
                scaled_src_arrays.append(scaled)
                if source_ok_packed_arrays is not None:
                    source_ok_packed_arrays.append(
                        _get_source_ok_packed_kernel(sidx, historyLength, kernel_width, scaled)
                    )
                # The optimized kernel path consumes scaled_sources directly.
                src_arrays.append(None)
                valid_sources.append(sidx)
            except Exception as e:
                logging.error(f"Failed to load source {sidx}: {e}")
        if not src_arrays:
            continue

        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_te_for_sources_with_context(
                dest_ctx,
                src_arrays,
                reuse_dest_neighbors=False,
                dense_threshold=DENSE_THRESHOLD,
                return_local=True,
                scaled_sources=scaled_src_arrays,
                source_ok_packed=source_ok_packed_arrays,
                packed_vector_block=KERNEL_PACK_VECTOR_BLOCK,
            )
            for sidx, te, local in zip(valid_sources, te_vals, local_vals):
                local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local)
                results.append((sidx, target_idx, te, local_bytes, local_len, local_dtype, local_codec))
        else:
            te_vals = compute_te_for_sources_with_context(
                dest_ctx,
                src_arrays,
                reuse_dest_neighbors=False,
                dense_threshold=DENSE_THRESHOLD,
                scaled_sources=scaled_src_arrays,
                source_ok_packed=source_ok_packed_arrays,
                packed_vector_block=KERNEL_PACK_VECTOR_BLOCK,
            )
            for sidx, te in zip(valid_sources, te_vals):
                results.append((sidx, target_idx, te))
    return results


def _kernel_te_perm_pvalue(
    dest_ctx: dict,
    series: np.ndarray,
    k_hist: int,
    num_perms: int = 100,
    seed: int = 0,
    store_local: bool = False,
) -> tuple:
    """Compute kernel TE for one series and a permutation p-value by shuffling source_values = series[k-1:-1]."""
    if store_local:
        te_vals, local_vals = compute_te_for_sources_with_context(
            dest_ctx,
            [series],
            reuse_dest_neighbors=False,
            dense_threshold=DENSE_THRESHOLD,
            return_local=True,
        )
    else:
        te_vals = compute_te_for_sources_with_context(
            dest_ctx,
            [series],
            reuse_dest_neighbors=False,
            dense_threshold=DENSE_THRESHOLD,
        )
        local_vals = None
    if te_vals:
        te_orig = float(te_vals[0])
        if store_local:
            local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local_vals[0])
        else:
            local_bytes = local_len = local_dtype = local_codec = None
    else:
        te_orig = 0.0
        if store_local:
            local_bytes, local_len, local_dtype, local_codec = b'', 0, LOCAL_TE_DTYPE_STR, 'none'
        else:
            local_bytes = local_len = local_dtype = local_codec = None
    if num_perms <= 0:
        if store_local:
            return te_orig, 1.0, local_bytes, local_len, local_dtype, local_codec
        return te_orig, 1.0
    rng = np.random.default_rng(int(seed))
    s = np.asarray(series, dtype=np.float64).copy()
    middle = s[k_hist - 1: -1]
    n = middle.shape[0]
    count = 0
    for _ in range(int(num_perms)):
        order = rng.permutation(n)
        s_perm = s.copy()
        s_perm[k_hist - 1: -1] = middle[order]
        te_p = compute_te_for_sources_with_context(
            dest_ctx,
            [s_perm],
            reuse_dest_neighbors=False,
            dense_threshold=DENSE_THRESHOLD,
        )
        te_perm = float(te_p[0]) if te_p else 0.0
        if te_perm >= te_orig:
            count += 1
    pval = (count + 1.0) / (num_perms + 1.0)
    if store_local:
        return te_orig, float(pval), local_bytes, local_len, local_dtype, local_codec
    return te_orig, float(pval)


def _kernel_grouped_target_context_for_permutation(
    target_idx: int,
    historyLength: int,
    kernel_width: float,
):
    """Build the same target grouped-state context used by fast kernel TE."""
    dest = _maybe_subsample(_get_series(int(target_idx)))
    if dest.size <= int(historyLength):
        return None
    dest_past = dest[:-1].astype(np.float64, copy=False)
    dest_next = dest[1:].astype(np.float64, copy=False)
    epsilon = 1e-9
    kw_past = float(kernel_width) * (np.std(dest_past, ddof=1) + epsilon)
    kw_next = float(kernel_width) * (np.std(dest_next, ddof=1) + epsilon)
    scaled_past = np.divide(
        dest_past,
        kw_past,
        out=np.zeros_like(dest_past, dtype=np.float64),
        where=(kw_past != 0),
    )
    scaled_next = np.divide(
        dest_next,
        kw_next,
        out=np.zeros_like(dest_next, dtype=np.float64),
        where=(kw_next != 0),
    )
    target_codes, target_past_values, target_next_values, target_counts = _target_grouped_state_codes(
        scaled_past,
        scaled_next,
    )
    target_factor_context = None
    use_factorized_target = str(KERNEL_GROUPED_STATE_FACTORIZED).lower() not in ("off", "0", "false", "no")
    if (
        use_factorized_target
        and not (
            str(KERNEL_GROUPED_STATE_FACTORIZED).lower() in ("auto", "")
            and int(scaled_past.shape[0]) < int(KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS)
        )
    ):
        target_factor_context = _target_grouped_state_factor_context(
            target_past_values,
            target_next_values,
            target_counts,
        )
    return {
        "target_codes": target_codes,
        "target_past_values": target_past_values,
        "target_next_values": target_next_values,
        "target_counts": target_counts,
        "target_factor_context": target_factor_context,
        "n_obs": int(scaled_past.shape[0]),
    }


def _kernel_grouped_target_state_count(
    target_idx: int,
    historyLength: int,
    kernel_width: float,
) -> int:
    """Cheap target complexity proxy used only for permutation work scheduling."""
    dest = _maybe_subsample(_get_series(int(target_idx)))
    if dest.size <= int(historyLength):
        return 0
    dest_past = dest[:-1].astype(np.float64, copy=False)
    dest_next = dest[1:].astype(np.float64, copy=False)
    epsilon = 1e-9
    kw_past = float(kernel_width) * (np.std(dest_past, ddof=1) + epsilon)
    kw_next = float(kernel_width) * (np.std(dest_next, ddof=1) + epsilon)
    scaled_past = np.divide(
        dest_past,
        kw_past,
        out=np.zeros_like(dest_past, dtype=np.float64),
        where=(kw_past != 0),
    )
    scaled_next = np.divide(
        dest_next,
        kw_next,
        out=np.zeros_like(dest_next, dtype=np.float64),
        where=(kw_next != 0),
    )
    return int(np.unique(np.stack((scaled_past, scaled_next), axis=1), axis=0).shape[0])


def _kernel_grouped_te_for_source_codes(
    source_values: np.ndarray,
    source_codes: np.ndarray,
    target_context: dict,
) -> float | None:
    """Compute exact grouped-state kernel TE from already-coded source values."""
    if int(source_codes.shape[0]) != int(target_context["n_obs"]):
        return 0.0
    max_unique = int(KERNEL_GROUPED_STATE_MAX_UNIQUE)
    te = None
    target_factor_context = target_context.get("target_factor_context")
    if target_factor_context is not None:
        te = _grouped_state_te_from_factorized(
            source_values,
            source_codes,
            target_context["target_codes"],
            target_context["target_counts"],
            target_context["target_past_values"],
            target_context["target_next_values"],
            target_factor_context,
            max_unique=max_unique,
        )
    if te is None:
        te = _grouped_state_te_from_coded(
            source_values,
            source_codes,
            target_context["target_codes"],
            target_context["target_past_values"],
            target_context["target_next_values"],
            max_unique=max_unique,
        )
    return te


def _kernel_te_perm_pvalue_grouped(
    target_context: dict,
    source_idx: int,
    k_hist: int,
    kernel_width: float,
    num_perms: int = 100,
    seed: int = 0,
) -> tuple[float, float]:
    """Permutation p-value using the fast exact grouped-state kernel path.

    This is equivalent to shuffling source_values = series[k-1:-1] because the
    source scaling and unique-value coding are source-only; a permutation only
    reorders the already-computed source code vector.
    """
    source_values, source_codes = _get_grouped_source_codes_kernel(
        int(source_idx),
        int(k_hist),
        float(kernel_width),
    )
    if source_codes.shape[0] != int(target_context["n_obs"]):
        return 0.0, 1.0
    te_orig = _kernel_grouped_te_for_source_codes(source_values, source_codes, target_context)
    if te_orig is None:
        raise RuntimeError(
            "Grouped-state kernel permutation exceeded configured unique-state limits; "
            "raise TE_KERNEL_GROUPED_STATE_MAX_UNIQUE or use fewer/less-continuous cells."
        )
    te_orig = float(te_orig)
    if int(num_perms) <= 0:
        return te_orig, 1.0

    rng = np.random.default_rng(int(seed))
    n_obs = int(source_codes.shape[0])
    use_table_sampler = str(KERNEL_PERM_TABLE_SAMPLER).lower() not in (
        "off",
        "0",
        "false",
        "no",
    )
    n_source_states = int(source_values.size)
    n_target_states = int(target_context["target_counts"].size)
    if (
        use_table_sampler
        and n_source_states > 0
        and n_target_states > 0
        and n_source_states * n_target_states <= int(KERNEL_PERM_TABLE_MAX_CELLS)
    ):
        source_counts = np.bincount(
            np.asarray(source_codes, dtype=np.int64),
            minlength=n_source_states,
        ).astype(np.int64, copy=False)
        target_counts = np.asarray(target_context["target_counts"], dtype=np.int64)
        if int(source_counts.sum()) == n_obs and int(target_counts.sum()) == n_obs:
            count = 0
            for _ in range(int(num_perms)):
                W_perm = _sample_permuted_joint_weights(source_counts, target_counts, rng)
                te_perm = _grouped_state_te_from_weight_matrix(
                    source_values,
                    W_perm,
                    target_context,
                    max_unique=0,
                )
                if te_perm is None:
                    break
                if float(te_perm) >= te_orig:
                    count += 1
            else:
                pval = (count + 1.0) / (int(num_perms) + 1.0)
                return te_orig, float(pval)

    count = 0
    for _ in range(int(num_perms)):
        perm_codes = source_codes[rng.permutation(n_obs)]
        te_perm = _kernel_grouped_te_for_source_codes(source_values, perm_codes, target_context)
        if te_perm is None:
            raise RuntimeError(
                "Grouped-state kernel permutation exceeded configured unique-state limits "
                "for a permuted source."
            )
        if float(te_perm) >= te_orig:
            count += 1
    pval = (count + 1.0) / (int(num_perms) + 1.0)
    return te_orig, float(pval)



def _perm_worker_kernel(args_w):
    """Permutation worker for kernel mode.
    Args tuple: (target_idx, source_indices, historyLength, perm_n, base_seed, kernel_width)
    Returns list of (Source, Target, TE, p_value)
    """
    tgt, srcs, k_hist, perm_n, base_seed, kw = args_w
    res = []
    # Payload storage must not choose the observed/null scoring engine.
    use_grouped_perm = _kernel_grouped_state_enabled("kernel", "default", int(k_hist))
    if use_grouped_perm:
        try:
            grouped_ctx = _kernel_grouped_target_context_for_permutation(tgt, k_hist, kw)
        except Exception as e:
            logging.error(f"Failed to build grouped kernel permutation context for target {tgt}: {e}")
            return res
        if grouped_ctx is None or int(grouped_ctx.get("n_obs", 0)) <= 0:
            return res
    else:
        try:
            ctx = _get_dest_ctx_kernel(tgt, k_hist, kw)
        except Exception as e:
            logging.error(f"Failed to build kernel context for target {tgt}: {e}")
            return res
        if ctx is None or ctx.get('total_obs', 0) <= 0:
            return res
    for s in srcs:
        try:
            if use_grouped_perm:
                te, p = _kernel_te_perm_pvalue_grouped(
                    grouped_ctx,
                    int(s),
                    int(k_hist),
                    float(kw),
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                )
                if STORE_LOCAL_TE:
                    values, codes = _get_grouped_source_codes_kernel(int(s), int(k_hist), float(kw))
                    local_result = _grouped_state_te_local_from_coded(
                        values, codes, grouped_ctx["target_codes"],
                        grouped_ctx["target_past_values"], grouped_ctx["target_next_values"],
                        max_unique=int(KERNEL_GROUPED_STATE_MAX_UNIQUE),
                        target_factor_context=grouped_ctx.get("target_factor_context"),
                    )
                    if local_result is None:
                        raise RuntimeError("LocalTE payload exceeds the grouped-state limit")
                    local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local_result[1])
                    res.append((s, tgt, te, p, local_bytes, local_len, local_dtype, local_codec))
                else:
                    res.append((s, tgt, te, p))
            elif STORE_LOCAL_TE:
                series = _get_series(s)
                series = _maybe_subsample(series)
                te, p, local_bytes, local_len, local_dtype, local_codec = _kernel_te_perm_pvalue(
                    ctx,
                    series,
                    k_hist,
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                    store_local=True,
                )
                res.append((s, tgt, te, p, local_bytes, local_len, local_dtype, local_codec))
            else:
                series = _get_series(s)
                series = _maybe_subsample(series)
                te, p = _kernel_te_perm_pvalue(
                    ctx,
                    series,
                    k_hist,
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                    store_local=False,
                )
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (kernel) failed for pair ({s},{tgt}): {e}")
    return res



def merge_parquet_files(progress_dir, merged_filename_prefix='merged_TE_progress'):
    """
    Merges all batch Parquet files in the progress directory into a single merged Parquet file using DuckDB.
    Handles incomplete/corrupted files by excluding them from the merge and deleting them.

    Args:
        progress_dir (str): Directory containing Parquet files to merge.
        merged_filename_prefix (str): Prefix for the merged Parquet file.
    """
    start_time = time.time()
    con = duckdb.connect()

    # Use high-resolution + random suffix to avoid filename collisions when merges
    # happen multiple times within the same second (which can lead to overwriting
    # and loss of previously merged content).
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    try:
        import uuid as _uuid
        rand = _uuid.uuid4().hex[:8]
    except Exception:
        rand = f"{int(time.time()*1e6)%1000000:06d}"
    merged_filename = f"{merged_filename_prefix}_{timestamp}_{rand}.parquet"
    merged_filepath = os.path.join(progress_dir, merged_filename)

    valid_files = []
    invalid_files = []  # Keep track of invalid files for deletion
    for fname in os.listdir(progress_dir):
        if fname.endswith('.parquet') and fname.startswith('batch_'):
            filepath = os.path.join(progress_dir, fname)
            try:
                # Attempt to read the Parquet file to check for validity
                con.execute(f"SELECT * FROM read_parquet('{filepath}') LIMIT 1")
                valid_files.append(filepath)
                logging.info(f"File {fname} is valid.")
            except Exception as e:
                logging.warning(f"Skipping invalid/incomplete file {fname}: {e}")
                invalid_files.append(filepath)  # Add invalid file to the list

    if not valid_files:
        logging.warning("No valid batch Parquet files found to merge.")
        con.close()
        return

    try:
        # Include only valid batch Parquet files in the merge
        file_list_str = ', '.join([f"'{file}'" for file in valid_files])
        query = f"""
        COPY (
            SELECT * FROM read_parquet([{file_list_str}])
        ) TO '{merged_filepath}' (FORMAT PARQUET);
        """
        con.execute(query)
        logging.info(f"Merged {len(valid_files)} batch Parquet files into {merged_filename} using DuckDB.")
    except Exception as e:
        logging.error(f"Error merging Parquet files: {e}")
    finally:
        con.close()
        end_time = time.time()
        logging.info(f"Parquet merging took {end_time - start_time:.2f} seconds.")

    # Delete invalid files
    for filepath in invalid_files:
        try:
            os.remove(filepath)
            logging.info(f"Deleted invalid file: {filepath}")
        except Exception as e:
            logging.error(f"Error deleting invalid file {filepath}: {e}")

    # Delete valid original batch files after merging
    for filepath in valid_files:
        try:
            os.remove(filepath)
            logging.info(f"Deleted batch file after merging: {filepath}")
        except Exception as e:
            logging.error(f"Error deleting batch file {filepath}: {e}")


def consolidate_merged_results(progress_dir, output_path, delete_after=False):
    """
    Combine all merged progress Parquet files into a single Parquet file.

    Args:
        progress_dir (str): Directory containing merged progress files.
        output_path (str): Destination Parquet filename.
        delete_after (bool): When True, delete merged inputs after consolidation.

    Returns:
        bool: True if consolidation wrote an output file, False otherwise.
    """
    if not progress_dir:
        logging.warning("No progress directory provided for consolidation.")
        return False

    merged_files = [
        os.path.join(progress_dir, fname)
        for fname in sorted(os.listdir(progress_dir))
        if fname.endswith('.parquet') and fname.startswith('merged_')
    ]
    if not merged_files:
        logging.warning("No merged Parquet files found in progress directory for consolidation.")
        return False

    output_path_abs = os.path.abspath(output_path)

    if pq is None or pa is None:
        logging.warning("pyarrow not available; falling back to DuckDB consolidation (higher memory usage).")
        con = duckdb.connect()
        try:
            file_list_str = ', '.join(f"'{path}'" for path in merged_files)
            query = f"""
            COPY (
                SELECT * FROM read_parquet([{file_list_str}])
            ) TO '{output_path_abs}' (FORMAT PARQUET);
            """
            con.execute(query)
            logging.info(f"Wrote consolidated results to {output_path_abs} from {len(merged_files)} merged files.")
        except Exception as e:
            logging.error(f"Error consolidating merged Parquet files into {output_path_abs}: {e}")
            raise
        finally:
            con.close()
    else:
        writer = None
        total_rows = 0
        start_time = time.time()
        try:
            for filepath in merged_files:
                try:
                    pf = pq.ParquetFile(filepath)
                except Exception as e:
                    logging.error(f"Failed to open merged parquet {filepath}: {e}")
                    raise
                for batch in pf.iter_batches(use_threads=False):
                    if batch.num_rows == 0:
                        continue
                    if writer is None:
                        writer = pq.ParquetWriter(output_path_abs, batch.schema, compression="snappy")
                    writer.write_table(pa.Table.from_batches([batch]))
                    total_rows += batch.num_rows
            if writer is None:
                logging.warning("Merged parquet files were empty; no output written.")
                return False
        except Exception as e:
            logging.error(f"Streaming consolidation failed: {e}")
            if writer is not None:
                writer.close()
            raise
        else:
            writer.close()
            elapsed = time.time() - start_time
            logging.info(f"Wrote consolidated results to {output_path_abs} ({total_rows:,} rows) from {len(merged_files)} merged files in {elapsed:.2f}s.")

    if delete_after:
        for filepath in merged_files:
            try:
                os.remove(filepath)
                logging.info(f"Deleted merged file {filepath}.")
            except Exception as e:
                logging.error(f"Error deleting merged file {filepath}: {e}")
    return True



def concat_parquet_files_duckdb(input_files: list[str], output_file: str) -> bool:
    """Concatenate Parquet files with DuckDB without materialising them in pandas."""
    files = [os.path.abspath(p) for p in input_files if p and os.path.exists(p)]
    if not files:
        return False
    out_path = os.path.abspath(output_file)
    file_list = "[" + ", ".join("'" + p.replace("'", "''") + "'" for p in files) + "]"
    con = duckdb.connect()
    try:
        con.execute(f"COPY (SELECT * FROM read_parquet({file_list})) TO '{out_path}' (FORMAT PARQUET)")
    finally:
        con.close()
    return True


def _is_auto_localte_stream_dir(progress_dir: str | None) -> bool:
    """True for the internal LocalTE temp dir used when user resume is off."""
    if not STORE_LOCAL_TE or not progress_dir:
        return False
    return os.path.basename(os.path.abspath(progress_dir)).startswith("TE_progress_localte_stream_")


def _grouped_result_columns() -> list[str]:
    if STORE_LOCAL_TE:
        return [
            "Source",
            "Target",
            "TE",
            "LocalTE_bytes",
            "LocalTE_len",
            "LocalTE_dtype",
            "LocalTE_codec",
        ]
    return ["Source", "Target", "TE"]


def _append_grouped_rows_to_parquet_writer(
    buffered_rows,
    columns: list[str],
    writer_state: dict,
    output_file: str,
) -> None:
    """Append one grouped-state row buffer directly as a Parquet row group."""
    if pa is None or pq is None:
        raise RuntimeError("pyarrow is required for direct LocalTE streaming output.")
    batch_df = pd.DataFrame.from_records(buffered_rows, columns=columns)
    table = pa.Table.from_pandas(batch_df, preserve_index=False)
    if writer_state.get("writer") is None:
        writer_state["writer"] = pq.ParquetWriter(output_file, table.schema, compression="snappy")
    writer_state["writer"].write_table(table)
    writer_state["rows"] = int(writer_state.get("rows", 0)) + int(len(batch_df))


def run_kernel_grouped_state(
    list_pairs,
    historyLength,
    kernel_width,
    num_cpus,
    progress_dir,
    merge_threshold,
    enable_intermediate_save,
    output_file,
    pool_chunksize,
    result_buffer_rows,
    direct_localte_exporter=None,
):
    """Run exact k=1 kernel TE by weighted unique-state counting."""
    if list_pairs is None or len(list_pairs) == 0:
        return False
    if int(historyLength) != 1:
        return False
    if not _grouped_state_auto_probe(list_pairs, int(historyLength), float(kernel_width)):
        return False
    if enable_intermediate_save:
        os.makedirs(progress_dir, exist_ok=True)

    # Stable grouping preserves both first-seen target order and within-target
    # source order, without 20 million Python pair unpack/dictionary operations.
    group_order = np.argsort(list_pairs[:, 1], kind="stable")
    sorted_targets = list_pairs[group_order, 1]
    sorted_sources = list_pairs[group_order, 0]
    group_starts = np.r_[0, np.flatnonzero(sorted_targets[1:] != sorted_targets[:-1]) + 1]
    group_ends = np.r_[group_starts[1:], len(group_order)]
    first_seen_order = np.argsort(group_order[group_starts], kind="stable")
    work_units = [
        (int(sorted_targets[group_starts[g]]), sorted_sources[group_starts[g]:group_ends[g]],
         int(historyLength), float(kernel_width))
        for g in first_seen_order
    ]
    del group_order, sorted_targets
    logging.info(
        "Using grouped-state exact kernel mode: targets=%d, pairs=%d.",
        len(work_units),
        len(list_pairs),
    )
    global _PREFIX_DATA
    prefix_enabled = (
        os.getenv("TE_EXACT_PREFIX", "on").lower() == "on"
        and not STORE_LOCAL_TE and not _is_on(LOCAL_TE_VIRTUAL)
        and direct_localte_exporter is None
        and not _is_auto_localte_stream_dir(progress_dir)
        and GLOBAL_CELL_GENE_DENSE is not None
        and TIME_SUBSAMPLE_INDICES is None
        and GLOBAL_CELL_GENE_DENSE.shape[1] - 1 >= KERNEL_GROUPED_STATE_FACTORIZED_MIN_OBS
        and str(KERNEL_GROUPED_STATE_FACTORIZED).lower() not in ("off", "0", "false", "no")
        and str(KERNEL_GROUPED_STATE_SOURCE_PREFIX).lower() not in ("off", "0", "false", "no")
        and str(KERNEL_GROUPED_STATE_2D_RANGE).lower() not in ("off", "0", "false", "no")
        and str(KERNEL_GROUPED_STATE_CODED).lower() not in ("off", "0", "false", "no")
        and str(KERNEL_GROUPED_STATE_SPARSE).lower() in ("off", "0", "false", "no")
        and GLOBAL_CELL_GENE_DENSE.shape[1] - 1 <= min(8192, KERNEL_GROUPED_STATE_MAX_UNIQUE)
        and KERNEL_GROUPED_STATE_FACTORIZED_MAX >= 250000
        and np.all(np.isfinite(GLOBAL_CELL_GENE_DENSE))
    )
    if prefix_enabled:
        from .exact_prefix import prepare_sources, prepare_target, count_batch
        prefix_start = time.time()
        _PREFIX_DATA = prepare_sources(GLOBAL_CELL_GENE_DENSE, sys.modules[__name__])
        warm_target, _ = prepare_target(GLOBAL_CELL_GENE_DENSE[0], sys.modules[__name__])
        count_batch(np.asarray([0], dtype=np.int64), _PREFIX_DATA[0], warm_target)
        logging.info("Exact sparse-prefix preparation took %.3f seconds.", time.time() - prefix_start)
        print("[TE] Exact sparse-prefix counts enabled; original NumPy scoring retained.")
    local_label = "LocalTE on" if STORE_LOCAL_TE else "LocalTE off"
    print(
        f"[TE] Engine: grouped-state exact kernel ({local_label}); "
        f"targets={_human_int(len(work_units))}, pairs={_human_int(len(list_pairs))}"
    )

    grouped_write_rows = max(0, int(KERNEL_GROUPED_STATE_WRITE_ROWS))
    if STORE_LOCAL_TE and result_buffer_rows is not None:
        grouped_write_rows = min(grouped_write_rows, max(1, int(result_buffer_rows)))
    if grouped_write_rows > 0:
        logging.info("Grouped-state result write buffer set to %d rows.", grouped_write_rows)
        if STORE_LOCAL_TE:
            print(f"[TE] LocalTE streaming buffer: {_human_int(grouped_write_rows)} rows per Parquet row group")
    else:
        logging.info("Grouped-state result write buffer disabled; writing one target batch at a time.")
    direct_localte_output = _is_auto_localte_stream_dir(progress_dir)
    if direct_localte_output:
        print("[TE] LocalTE output mode: direct streaming to TE_result_all.parquet; final merge skipped.")
    if (
        _grouped_state_te_numba is not None
        and str(KERNEL_GROUPED_STATE_NUMBA).lower() not in ("off", "0", "false", "no")
    ):
        try:
            _grouped_state_te_numba(
                np.asarray([0.0], dtype=np.float64),
                np.asarray([0.0], dtype=np.float64),
                np.asarray([0.0], dtype=np.float64),
                np.asarray([1], dtype=np.int64),
                1,
            )
            logging.info("Grouped-state numba TE loop is enabled.")
        except Exception as exc:
            logging.warning("Grouped-state numba warmup failed; continuing with NumPy path: %s", exc)
    grouped_source_code_memmap_path = None
    try:
        grouped_source_code_memmap_path = _prepare_grouped_source_code_memmap(
            list_pairs[:, 0],
            int(historyLength),
            float(kernel_width),
            progress_dir,
        )
    except Exception as exc:
        logging.warning("Grouped source code memmap disabled after preparation failure: %s", exc)
        _clear_grouped_source_code_memmap()

    if direct_localte_exporter is not None and _is_on(LOCAL_TE_DIRECT_WORKERS):
        try:
            return _run_grouped_state_direct_worker_writers(
                work_units,
                int(len(list_pairs)),
                int(num_cpus),
                output_file,
                direct_localte_exporter,
                int(result_buffer_rows or grouped_write_rows or 5000),
            )
        finally:
            if grouped_source_code_memmap_path:
                _clear_grouped_source_code_memmap()

    direct_writer_state = {"writer": None, "rows": 0}
    columns = ["Source", "Target", "TE"] if direct_localte_exporter is not None else _grouped_result_columns()

    def flush_grouped_rows(buffered_rows):
        if not buffered_rows:
            return []
        if prefix_enabled:
            batch_df = pd.DataFrame({
                "Source": np.concatenate([b.sources for b in buffered_rows]),
                "Target": np.concatenate([np.full(len(b), b.target, dtype=np.int64) for b in buffered_rows]),
                "TE": np.concatenate([b.scores for b in buffered_rows]),
            })
            rows_to_write = None
        else:
            rows_to_write = _direct_localte_append_rows(direct_localte_exporter, buffered_rows)
        if direct_localte_output:
            _append_grouped_rows_to_parquet_writer(
                rows_to_write,
                columns,
                direct_writer_state,
                output_file,
            )
            return []
        if not prefix_enabled:
            batch_df = pd.DataFrame.from_records(rows_to_write, columns=columns)
        if enable_intermediate_save:
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            try:
                import uuid as _uuid
                suf = _uuid.uuid4().hex[:8]
            except Exception:
                suf = f"{int(time.time()*1e6)%1000000:06d}"
            batch_filename = os.path.join(progress_dir, f"batch_grouped_{ts}_{suf}.parquet")
            batch_df.to_parquet(batch_filename, index=False)
            logging.info("Written %s (%d rows)", batch_filename, len(batch_df))
            current_batch_file_count = len([
                fname for fname in os.listdir(progress_dir)
                if fname.endswith(".parquet") and fname.startswith("batch_")
            ])
            if current_batch_file_count >= merge_threshold:
                merge_parquet_files(progress_dir)
        else:
            all_results.append(batch_df)
        return []

    all_results = []
    row_buffer = []
    buffered_count = 0
    start_time = time.time()
    grouped_chunksize = 1
    show_progress = str(KERNEL_GROUPED_STATE_PROGRESS).lower() not in ("off", "0", "false", "no")
    processed_targets = 0
    processed_pairs = 0
    try:
        with multiprocessing.Pool(processes=num_cpus) as pool:
            iterator = pool.imap_unordered(
                _process_prefix_target if prefix_enabled else _process_grouped_state_target,
                work_units,
                chunksize=grouped_chunksize,
            )
            target_pbar = None
            pair_pbar = None
            if show_progress:
                target_pbar = tqdm(
                    total=len(work_units),
                    desc="TE targets",
                    unit="target",
                    dynamic_ncols=True,
                    mininterval=1.0,
                    leave=True,
                    bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}",
                )
                pair_pbar = tqdm(
                    total=len(list_pairs),
                    desc="TE pairs (compute)",
                    unit="pair",
                    dynamic_ncols=True,
                    mininterval=1.0,
                    leave=True,
                )
            try:
                for target_rows in iterator:
                    n_rows = len(target_rows) if target_rows else 0
                    processed_targets += 1
                    processed_pairs += n_rows
                    if target_pbar is not None:
                        target_pbar.update(1)
                    if pair_pbar is not None and n_rows:
                        pair_pbar.update(n_rows)
                    if not target_rows:
                        continue
                    if grouped_write_rows > 0:
                        if prefix_enabled:
                            row_buffer.append(target_rows)
                        else:
                            row_buffer.extend(target_rows)
                        buffered_count += n_rows
                        if buffered_count >= grouped_write_rows:
                            row_buffer = flush_grouped_rows(row_buffer)
                            buffered_count = 0
                    else:
                        flush_grouped_rows([target_rows] if prefix_enabled else target_rows)
            finally:
                if target_pbar is not None:
                    target_pbar.close()
                if pair_pbar is not None:
                    pair_pbar.close()
    finally:
        if grouped_source_code_memmap_path:
            _clear_grouped_source_code_memmap()

    if row_buffer:
        row_buffer = flush_grouped_rows(row_buffer)

    if direct_writer_state.get("writer") is not None:
        direct_writer_state["writer"].close()
        direct_writer_state["writer"] = None
    _direct_localte_finalize(direct_localte_exporter)

    if direct_localte_output:
        if int(direct_writer_state.get("rows", 0)) <= 0:
            return False
    elif enable_intermediate_save:
        current_batch_file_count = len([
            fname for fname in os.listdir(progress_dir)
            if fname.endswith(".parquet") and fname.startswith("batch_")
        ])
        if current_batch_file_count > 0:
            merge_parquet_files(progress_dir)
        consolidated = consolidate_merged_results(progress_dir, output_file, delete_after=True)
        if not consolidated:
            return False
    elif all_results:
        pd.concat(all_results, ignore_index=True).to_parquet(output_file, index=False)
    else:
        return False

    elapsed = time.time() - start_time
    logging.info("Grouped-state exact kernel results saved to %s in %.2f seconds.", output_file, elapsed)
    return True


def run_kernel_grouped_state_implicit(
    indices,
    pair_mode,
    historyLength,
    kernel_width,
    num_cpus,
    progress_dir,
    merge_threshold,
    enable_intermediate_save,
    output_file,
    pool_chunksize,
    result_buffer_rows,
    direct_localte_exporter=None,
):
    """Run exact grouped-state kernel TE for implicit all-pair modes.

    This keeps the estimator identical to `run_kernel_grouped_state`, but avoids
    both materialising all pairs and sending a full source list with every task.
    """
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size <= 1 or int(historyLength) != 1:
        return False
    if not _grouped_state_auto_probe_implicit(indices, int(historyLength), float(kernel_width)):
        return False

    total_pairs = int(indices.size) * int(indices.size - 1)
    effective_progress_dir = progress_dir
    if enable_intermediate_save:
        # Keep implicit all-pair batches isolated from any interrupted fallback
        # run that may have left batch_*.parquet in the default progress dir.
        effective_progress_dir = os.path.join(progress_dir, f"grouped_implicit_{pair_mode}")
        os.makedirs(effective_progress_dir, exist_ok=True)

    logging.info(
        "Using implicit grouped-state exact kernel mode: pair_mode=%s, targets=%d, pairs=%d.",
        pair_mode,
        indices.size,
        total_pairs,
    )
    print(
        f"[TE] Engine: implicit grouped-state exact kernel ({pair_mode}, "
        f"{'LocalTE on' if STORE_LOCAL_TE else 'LocalTE off'}); "
        f"targets={_human_int(indices.size)}, pairs={_human_int(total_pairs)}"
    )

    grouped_write_rows = max(0, int(KERNEL_GROUPED_STATE_WRITE_ROWS))
    if STORE_LOCAL_TE and result_buffer_rows is not None:
        grouped_write_rows = min(grouped_write_rows, max(1, int(result_buffer_rows)))
    if grouped_write_rows > 0:
        logging.info("Grouped-state result write buffer set to %d rows.", grouped_write_rows)
        if STORE_LOCAL_TE:
            print(f"[TE] LocalTE streaming buffer: {_human_int(grouped_write_rows)} rows per Parquet row group")
    else:
        logging.info("Grouped-state result write buffer disabled; writing one target batch at a time.")
    direct_localte_output = _is_auto_localte_stream_dir(effective_progress_dir)
    if direct_localte_output:
        print("[TE] LocalTE output mode: direct streaming to TE_result_all.parquet; final merge skipped.")

    if (
        _grouped_state_te_numba is not None
        and str(KERNEL_GROUPED_STATE_NUMBA).lower() not in ("off", "0", "false", "no")
    ):
        try:
            _grouped_state_te_numba(
                np.asarray([0.0], dtype=np.float64),
                np.asarray([0.0], dtype=np.float64),
                np.asarray([0.0], dtype=np.float64),
                np.asarray([1], dtype=np.int64),
                1,
            )
            logging.info("Grouped-state numba TE loop is enabled.")
        except Exception as exc:
            logging.warning("Grouped-state numba warmup failed; continuing with NumPy path: %s", exc)

    grouped_source_code_memmap_path = None
    try:
        grouped_source_code_memmap_path = _prepare_grouped_source_code_memmap(
            indices,
            int(historyLength),
            float(kernel_width),
            effective_progress_dir or ".",
            force=True,
        )
    except Exception as exc:
        logging.warning("Grouped source code memmap disabled after preparation failure: %s", exc)
        _clear_grouped_source_code_memmap()

    direct_writer_state = {"writer": None, "rows": 0}
    columns = ["Source", "Target", "TE"] if direct_localte_exporter is not None else _grouped_result_columns()

    def flush_grouped_rows(buffered_rows):
        if not buffered_rows:
            return []
        rows_to_write = _direct_localte_append_rows(direct_localte_exporter, buffered_rows)
        if direct_localte_output:
            _append_grouped_rows_to_parquet_writer(
                rows_to_write,
                columns,
                direct_writer_state,
                output_file,
            )
            return []
        batch_df = pd.DataFrame.from_records(rows_to_write, columns=columns)
        if enable_intermediate_save:
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            try:
                import uuid as _uuid
                suf = _uuid.uuid4().hex[:8]
            except Exception:
                suf = f"{int(time.time()*1e6)%1000000:06d}"
            batch_filename = os.path.join(effective_progress_dir, f"batch_grouped_{ts}_{suf}.parquet")
            batch_df.to_parquet(batch_filename, index=False)
            logging.info("Written %s (%d rows)", batch_filename, len(batch_df))
            current_batch_file_count = len([
                fname for fname in os.listdir(effective_progress_dir)
                if fname.endswith(".parquet") and fname.startswith("batch_")
            ])
            if current_batch_file_count >= merge_threshold:
                merge_parquet_files(effective_progress_dir)
        else:
            all_results.append(batch_df)
        return []

    all_results = []
    row_buffer = []
    buffered_count = 0
    start_time = time.time()
    show_progress = str(KERNEL_GROUPED_STATE_PROGRESS).lower() not in ("off", "0", "false", "no")
    processed_targets = 0
    processed_pairs = 0
    work_units = [(int(tgt), int(historyLength), float(kernel_width)) for tgt in indices]

    try:
        with multiprocessing.Pool(
            processes=num_cpus,
            initializer=_init_grouped_state_implicit,
            initargs=(indices,),
        ) as pool:
            iterator = pool.imap_unordered(
                _process_grouped_state_target_implicit,
                work_units,
                chunksize=max(1, int(pool_chunksize)),
            )
            target_pbar = None
            pair_pbar = None
            if show_progress:
                target_pbar = tqdm(
                    total=len(work_units),
                    desc="TE targets",
                    unit="target",
                    dynamic_ncols=True,
                    mininterval=1.0,
                    leave=True,
                    bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}",
                )
                pair_pbar = tqdm(
                    total=total_pairs,
                    desc="TE pairs (compute)",
                    unit="pair",
                    dynamic_ncols=True,
                    mininterval=1.0,
                    leave=True,
                )
            try:
                for target_rows in iterator:
                    n_rows = len(target_rows) if target_rows else 0
                    processed_targets += 1
                    processed_pairs += n_rows
                    if target_pbar is not None:
                        target_pbar.update(1)
                    if pair_pbar is not None and n_rows:
                        pair_pbar.update(n_rows)
                    if not target_rows:
                        continue
                    if grouped_write_rows > 0:
                        row_buffer.extend(target_rows)
                        buffered_count += n_rows
                        if buffered_count >= grouped_write_rows:
                            row_buffer = flush_grouped_rows(row_buffer)
                            buffered_count = 0
                    else:
                        flush_grouped_rows(target_rows)
            finally:
                if target_pbar is not None:
                    target_pbar.close()
                if pair_pbar is not None:
                    pair_pbar.close()
    finally:
        if grouped_source_code_memmap_path:
            _clear_grouped_source_code_memmap()

    if row_buffer:
        row_buffer = flush_grouped_rows(row_buffer)

    if direct_writer_state.get("writer") is not None:
        direct_writer_state["writer"].close()
        direct_writer_state["writer"] = None
    _direct_localte_finalize(direct_localte_exporter)

    if direct_localte_output:
        if int(direct_writer_state.get("rows", 0)) <= 0:
            return False
    elif enable_intermediate_save:
        current_batch_file_count = len([
            fname for fname in os.listdir(effective_progress_dir)
            if fname.endswith(".parquet") and fname.startswith("batch_")
        ])
        if current_batch_file_count > 0:
            merge_parquet_files(effective_progress_dir)
        consolidated = consolidate_merged_results(effective_progress_dir, output_file, delete_after=True)
        if not consolidated:
            return False
    elif all_results:
        pd.concat(all_results, ignore_index=True).to_parquet(output_file, index=False)
    else:
        return False

    elapsed = time.time() - start_time
    logging.info(
        "Implicit grouped-state exact kernel results saved to %s in %.2f seconds "
        "(targets=%d, pairs=%d).",
        output_file,
        elapsed,
        processed_targets,
        processed_pairs,
    )
    return True


def run_kernel_source_blocked(
    list_pairs,
    cell_gene_all,
    historyLength,
    kernel_width,
    num_cpus,
    batch_size,
    progress_dir,
    buffer_size,
    merge_threshold,
    enable_intermediate_save,
    output_file,
    pool_chunksize,
    source_block_size,
):
    """Run exact kernel TE in bounded source blocks to share masks without a huge mmap."""
    if list_pairs is None or len(list_pairs) == 0:
        return False
    source_block_size = int(source_block_size)
    if source_block_size <= 0:
        return False

    unique_sources, source_counts = np.unique(list_pairs[:, 0].astype(np.int64), return_counts=True)
    if unique_sources.size <= source_block_size:
        return False

    effective_progress_dir = progress_dir or "."
    block_root = os.path.join(effective_progress_dir, "kernel_source_blocks")
    os.makedirs(block_root, exist_ok=True)
    block_outputs: list[str] = []
    total_blocks = int((unique_sources.size + source_block_size - 1) // source_block_size)
    # Source-heavy pair lists can be highly skewed by source index. Balance
    # blocks by pair count so one huge source block does not dominate the run.
    blocks: list[list[int]] = [[] for _ in range(total_blocks)]
    block_pair_counts = [0 for _ in range(total_blocks)]
    source_order = np.argsort(-source_counts)
    for source_pos in source_order:
        available = [
            i for i in range(total_blocks)
            if len(blocks[i]) < source_block_size
        ]
        block_idx = min(available, key=lambda i: block_pair_counts[i])
        blocks[block_idx].append(int(unique_sources[source_pos]))
        block_pair_counts[block_idx] += int(source_counts[source_pos])
    logging.info(
        "Using source-blocked kernel mode: sources=%d, block_size=%d, blocks=%d.",
        unique_sources.size,
        source_block_size,
        total_blocks,
    )
    print(
        f"[TE] Source-blocked kernel mode: {unique_sources.size} sources, "
        f"block_size={source_block_size}, blocks={total_blocks}"
    )

    for block_idx, block_source_list in enumerate(blocks, start=1):
        block_sources = np.array(sorted(block_source_list), dtype=np.int64)
        source_mask = np.isin(list_pairs[:, 0], block_sources)
        block_pairs = list_pairs[source_mask]
        if len(block_pairs) == 0:
            continue

        block_progress_dir = os.path.join(block_root, f"block_{block_idx:04d}")
        os.makedirs(block_progress_dir, exist_ok=True)
        block_output = os.path.abspath(os.path.join(block_root, f"TE_source_block_{block_idx:04d}.parquet"))
        logging.info(
            "Running source block %d/%d: sources=%d, pairs=%d.",
            block_idx,
            total_blocks,
            block_sources.size,
            len(block_pairs),
        )
        print(
            f"[TE] Source block {block_idx}/{total_blocks}: "
            f"{block_sources.size} sources, {len(block_pairs)} pairs"
        )

        try:
            _prepare_kernel_source_pack_memmap(
                block_sources,
                int(historyLength),
                float(kernel_width),
                block_progress_dir,
            )
            run_parallel_batches(
                list_pairs=block_pairs,
                cell_gene_all=cell_gene_all,
                historyLength=historyLength,
                kernel_width=kernel_width,
                num_cpus=num_cpus,
                batch_size=batch_size,
                progress_dir=block_progress_dir,
                buffer_size=buffer_size,
                merge_threshold=merge_threshold,
                enable_intermediate_save=enable_intermediate_save,
                mode="kernel",
                output_file=block_output,
                pair_mode="default",
                pair_index_info=None,
                pool_chunksize=pool_chunksize,
            )
            if enable_intermediate_save and not os.path.exists(block_output):
                consolidated = consolidate_merged_results(block_progress_dir, block_output, delete_after=False)
                if not consolidated:
                    raise RuntimeError(f"No source-block progress files available for block {block_idx}.")
            block_outputs.append(block_output)
        finally:
            _cleanup_kernel_source_pack_memmap()

    if not block_outputs:
        return False

    if not concat_parquet_files_duckdb(block_outputs, output_file):
        raise RuntimeError("Failed to concatenate source-block kernel outputs.")
    print(f"Fast-mode results saved to {output_file}.")
    logging.info("Source-blocked kernel results saved to %s from %d block files.", output_file, len(block_outputs))
    return True

def run_parallel_batches(
    list_pairs,
    cell_gene_all,
    historyLength,
    kernel_width,
    num_cpus,
    batch_size,
    progress_dir,
    buffer_size=100000,
    merge_threshold=50,
    enable_intermediate_save=True,
    mode='kernel',
    output_file='TE_result_all.parquet',
    pair_mode: str = "default",
    pair_index_info: dict | None = None,
    pool_chunksize: int = 1,
):
    """
    Processes the list of pairs in batches using multiprocessing Pool.
    Writes results to Parquet files incrementally in the specified progress directory if enable_intermediate_save is True.
    Merges only the new batch Parquet files when the number of batch files exceeds the merge_threshold.

    Args:
        list_pairs (np.ndarray): Array of pairs to process.
        cell_gene_all (scipy.sparse.csr_matrix): Sparse matrix of gene expressions.
        historyLength (int): History length (k) for the analysis.
        num_cpus (int): Number of parallel CPUs to use.
        batch_size (int): Number of pairs per batch.
        progress_dir (str): Directory to store intermediate Parquet files.
        buffer_size (int): Number of batches to accumulate before writing.
        merge_threshold (int): Number of batch Parquet files to trigger a merge.
        enable_intermediate_save (bool): Flag to enable/disable intermediate saving.
    """

    if enable_intermediate_save:
        # Ensure the progress directory exists
        os.makedirs(progress_dir, exist_ok=True)

    # Build work units. For the default mode, we group the explicit list_pairs by
    # target. For all-pair streaming modes, we construct (target, source_chunk)
    # implicitly from the feature indices to avoid materialising all pairs.
    work_units = None

    if pair_mode == "default":
        # Group by target to reuse destination computations
        target_to_sources = {}
        for src, tgt in list_pairs:
            target_to_sources.setdefault(int(tgt), []).append(int(src))

        avg_sources_per_target = (
            float(sum(len(srcs) for srcs in target_to_sources.values())) / max(1, len(target_to_sources))
        )
        target_work_setting = str(KERNEL_TARGET_WORK_UNITS or "auto").lower()
        target_work_units = (
            mode == "kernel"
            and target_work_setting not in ("off", "0", "false", "no")
            and _has_source_pack_memmap()
            and (
                target_work_setting in ("on", "1", "true", "yes")
                or avg_sources_per_target >= max(1, int(KERNEL_TARGET_MIN_SOURCES))
            )
        )
        if target_work_units:
            work_units = [
                (tgt, srcs, historyLength, kernel_width, mode)
                for tgt, srcs in target_to_sources.items()
            ]
            logging.info(
                "Using target-level kernel work units with shared source masks: targets=%d.",
                len(work_units),
            )
        else:
            if mode == "kernel" and _has_source_pack_memmap():
                logging.info(
                    "Using chunk-level kernel work units with shared source masks: targets=%d, avg_sources_per_target=%.2f.",
                    len(target_to_sources),
                    avg_sources_per_target,
                )
            work_units = [
                (tgt, srcs[i : i + batch_size], historyLength, kernel_width, mode)
                for tgt, srcs in target_to_sources.items()
                for i in range(0, len(srcs), batch_size)
            ]
        total_batches = len(work_units)

        def work_iter():
            for w in work_units:
                yield w

    else:
        if pair_index_info is None or "indices" not in pair_index_info:
            raise ValueError("pair_index_info with 'indices' is required for all-pair modes.")
        indices = list(pair_index_info["indices"])
        n = len(indices)
        if n <= 1:
            return
        target_work_setting = str(KERNEL_TARGET_WORK_UNITS or "auto").lower()
        # Implicit all-pair modes have n-1 sources per target, so the target
        # owner strategy is only worthwhile once that count is large enough.
        target_work_units = (
            mode == "kernel"
            and target_work_setting not in ("off", "0", "false", "no")
            and _has_source_pack_memmap()
            and (
                target_work_setting in ("on", "1", "true", "yes")
                or (n - 1) >= max(1, int(KERNEL_TARGET_MIN_SOURCES))
            )
        )
        # Number of chunks per target: ceil((n-1)/batch_size), unless a
        # shared source-mask memmap lets one worker own the full target.
        per_target_chunks = 1 if target_work_units else (n - 1 + batch_size - 1) // batch_size
        total_batches = n * per_target_chunks

        def work_iter():
            for tgt in indices:
                chunk = []
                for src in indices:
                    if src == tgt:
                        continue
                    chunk.append(src)
                    if (not target_work_units) and len(chunk) >= batch_size:
                        yield (tgt, list(chunk), historyLength, kernel_width, mode)
                        chunk.clear()
                if chunk:
                    yield (tgt, list(chunk), historyLength, kernel_width, mode)

    # Initialize buffer for accumulating batches
    buffered_batches = []
    buffered_rows = 0
    all_results = []

    # Use multiprocessing Pool for parallel processing
    with multiprocessing.Pool(processes=num_cpus) as pool:
        with tqdm(total=total_batches, desc="Processing Chunks") as pbar:
            for batch_result in pool.imap_unordered(_process_chunk, work_iter(), chunksize=max(1, int(pool_chunksize))):
                if not batch_result:
                    pbar.update(1)
                    continue

                # Build columns directly; per-row dict construction is costly for LocalTE payloads.
                if STORE_LOCAL_TE:
                    sources, targets, te_values, local_bytes, local_lens, local_dtypes, local_codecs = zip(*batch_result)
                    batch_df = pd.DataFrame(
                        {
                            'Source': sources,
                            'Target': targets,
                            'TE': te_values,
                            'LocalTE_bytes': local_bytes,
                            'LocalTE_len': local_lens,
                            'LocalTE_dtype': local_dtypes,
                            'LocalTE_codec': local_codecs,
                        }
                    )
                else:
                    sources, targets, te_values = zip(*batch_result)
                    batch_df = pd.DataFrame({'Source': sources, 'Target': targets, 'TE': te_values})

                if enable_intermediate_save:
                    buffered_batches.append(batch_df)
                    buffered_rows += len(batch_df)
                    # If buffer is full by rows, write to a batch Parquet file
                    if buffered_rows >= buffer_size:
                        combined_df = pd.concat(buffered_batches, ignore_index=True)
                        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                        try:
                            import uuid as _uuid
                            suf = _uuid.uuid4().hex[:8]
                        except Exception:
                            suf = f"{int(time.time()*1e6)%1000000:06d}"
                        batch_filename = os.path.join(progress_dir, f'batch_{ts}_{suf}.parquet')
                        try:
                            combined_df.to_parquet(batch_filename, index=False)
                            logging.info(f"Written {batch_filename}")
                        except Exception as e:
                            logging.error(f"Error writing {batch_filename}: {e}")
                            raise e
                        buffered_batches = []  # Reset buffer
                        buffered_rows = 0

                        # Check if merging is needed
                        current_batch_file_count = len([
                            fname for fname in os.listdir(progress_dir)
                            if fname.endswith('.parquet') and fname.startswith('batch_')
                        ])
                        if current_batch_file_count >= merge_threshold:
                            logging.info(f"Batch file count {current_batch_file_count} reached merge threshold {merge_threshold}. Merging batch files.")
                            merge_parquet_files(progress_dir)
                else:
                    all_results.append(batch_df)

                pbar.update(1)

    if enable_intermediate_save:
        # Write any remaining batches in the buffer
        if buffered_batches:
            combined_df = pd.concat(buffered_batches, ignore_index=True)
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            try:
                import uuid as _uuid
                suf = _uuid.uuid4().hex[:8]
            except Exception:
                suf = f"{int(time.time()*1e6)%1000000:06d}"
            batch_filename = os.path.join(progress_dir, f'batch_{ts}_{suf}.parquet')
            try:
                combined_df.to_parquet(batch_filename, index=False)
                logging.info(f"Written {batch_filename}")
            except Exception as e:
                logging.error(f"Error writing {batch_filename}: {e}")
                raise e

        # Final check for merging remaining batch files
        current_batch_file_count = len([
            fname for fname in os.listdir(progress_dir)
            if fname.endswith('.parquet') and fname.startswith('batch_')
        ])
        if current_batch_file_count > 0:
            logging.info(f"Final batch file count {current_batch_file_count}. Merging remaining batch files.")
            merge_parquet_files(progress_dir)
    else:
        # Combine all results if intermediate saving is disabled
        if all_results:
            final_results_df = pd.concat(all_results, ignore_index=True)
            final_results_df.to_parquet(output_file, index=False)
            logging.info(f"Final results saved to {output_file}.")

def load_progress(input_csv, progress_dir):
    """
    Loads progress from the progress directory if it exists.
    Returns the remaining pairs to process.
    Handles incomplete/corrupted files by skipping them.

    Args:
        input_csv (str): Path to the input CSV file containing all pairs.
        progress_dir (str): Directory where intermediate Parquet files are stored.

    Returns:
        np.ndarray: Array of remaining pairs to process.
    """
    load_start_time = time.time()

    # Ensure progress directory exists
    if not os.path.isdir(progress_dir):
        os.makedirs(progress_dir, exist_ok=True)

    # Check for existing merged Parquet files and load them
    merged_files = [
        fname for fname in os.listdir(progress_dir)
        if fname.endswith('.parquet') and fname.startswith('merged_')
    ]

    if len(merged_files) > 1:
        logging.info(f"Multiple merged Parquet files detected at startup. Merging them.")
        merge_parquet_files(progress_dir)
        # After merging, there should be only one merged file
        merged_files = [
            fname for fname in os.listdir(progress_dir)
            if fname.endswith('.parquet') and fname.startswith('merged_')
        ]
    elif len(merged_files) == 1:
        logging.info(f"Single merged Parquet file detected: {merged_files[0]}.")
    else:
        logging.info("No merged Parquet files detected in progress directory.")

    # Load all pairs from input CSV
    try:
        list_pairs = pd.read_csv(input_csv, delimiter=',', header=None).to_numpy().astype(int)
        logging.info(f"Loaded {len(list_pairs)} pairs from {input_csv}.")
    except Exception as e:
        logging.error(f"Error loading input CSV {input_csv}: {e}")
        raise

    pairs_df = pd.DataFrame(list_pairs, columns=['Source', 'Target'])

    # Load processed pairs from merged Parquet files, handling potential errors
    processed_pairs = []
    for fname in merged_files:
        filepath = os.path.join(progress_dir, fname)
        try:
            df = pd.read_parquet(filepath)
            processed_pairs.append(df[['Source', 'Target']])
            logging.info(f"Loaded processed pairs from {fname}.")
        except Exception as e:
            logging.error(f"Error reading progress file {fname}: {e}")
            raise

    if processed_pairs:
        processed_pairs_df = pd.concat(processed_pairs, ignore_index=True).drop_duplicates()
        logging.info(f"Total processed pairs: {len(processed_pairs_df)}.")

        # Identify remaining pairs by excluding processed ones
        remaining_pairs_df = pairs_df.merge(
            processed_pairs_df,
            on=['Source', 'Target'],
            how='left',
            indicator=True
        )
        remaining_pairs = remaining_pairs_df[remaining_pairs_df['_merge'] == 'left_only'][['Source', 'Target']].to_numpy()
        logging.info(f"Remaining pairs to process: {len(remaining_pairs)}.")

    else:
        logging.info("No processed pairs found in progress directory. Starting fresh.")
        remaining_pairs = list_pairs

    load_end_time = time.time()
    logging.info(f"Loading progress took {load_end_time - load_start_time:.2f} seconds.")
    return remaining_pairs

def main(args):
    """
    Main function to execute the TE analysis workflow.

    Args:
        args (Namespace): Parsed command-line arguments.
    """
    configure_logging()
    start_time = time.time()
    print("[TE] Starting kernel TE calculation.")
    logging.info('Starting Calculate_TE')
    if getattr(args, 'mode', 'kernel') != 'kernel':
        raise SystemExit(
            "This default runner is kernel-only. "
            "This standalone package only supports kernel TE."
        )

    # Load expression data
    print("[TE] Loading expression matrix...")
    logging.info('Loading expression data...')
    load_data_start_time = time.time()
    try:
        cell_gene_dense = pd.read_parquet('cell_gene_trsps.parquet').to_numpy(dtype=np.float64)
        dense_matrix_enabled = str(USE_DENSE_MATRIX).lower() in ("on", "auto")
        if dense_matrix_enabled:
            # Kernel fast paths read rows from GLOBAL_CELL_GENE_DENSE directly.
            # Keep only a shape-compatible sparse placeholder to avoid duplicating
            # the whole expression matrix in CSR form.
            cell_gene_all = sparse.csr_matrix(cell_gene_dense.shape, dtype=np.float64)
        else:
            cell_gene_all = sparse.csr_matrix(cell_gene_dense)
        logging.info('Expression data loaded successfully.')
    except Exception as e:
        print(f"Error loading expression data: {e}")
        logging.error(f"Error loading expression data: {e}")
        raise
    load_data_end_time = time.time()
    logging.info(f"Loading expression data took {load_data_end_time - load_data_start_time:.2f} seconds.")

    # Expose matrix as global for worker processes (avoids pickling per task)
    global GLOBAL_CELL_GENE, GLOBAL_CELL_GENE_DENSE
    GLOBAL_CELL_GENE = cell_gene_all
    GLOBAL_CELL_GENE_DENSE = cell_gene_dense if dense_matrix_enabled else None
    if GLOBAL_CELL_GENE_DENSE is not None:
        logging.info(
            "Dense row matrix enabled for %s mode: shape=%s, nbytes=%.2f MiB",
            getattr(args, 'mode', 'unknown'),
            tuple(GLOBAL_CELL_GENE_DENSE.shape),
            GLOBAL_CELL_GENE_DENSE.nbytes / (1024 * 1024),
        )
        print(
            f"[TE] Matrix loaded: features={_human_int(GLOBAL_CELL_GENE_DENSE.shape[0])}, "
            f"cells={_human_int(GLOBAL_CELL_GENE_DENSE.shape[1])}, "
            f"dense={_human_bytes(GLOBAL_CELL_GENE_DENSE.nbytes)}"
        )

    global STORE_LOCAL_TE
    STORE_LOCAL_TE = bool(getattr(args, 'store_local_te', False))
    if _is_on(LOCAL_TE_VIRTUAL):
        print("[TE] LocalTE: VIRTUAL exact array. TE scores stay score-only; timepoint values are reconstructed from compact state tables.")
        logging.info("Virtual LocalTE storage enabled; outputs will include compact state tables.")
    elif STORE_LOCAL_TE:
        logging.info("Local TE storage enabled; outputs will include LocalTE arrays.")
        print("[TE] LocalTE: ON. Per-pair timepoint arrays will be written to the TE parquet.")
    else:
        print("[TE] LocalTE: OFF.")
    # Configure LocalTE codec
    global LOCAL_TE_CODEC
    try:
        codec_arg = getattr(args, 'localte_codec', None)
        if codec_arg:
            LOCAL_TE_CODEC = str(codec_arg).lower()
        else:
            LOCAL_TE_CODEC = os.getenv('TE_LOCALTE_CODEC', LOCAL_TE_CODEC or 'zlib').lower()
    except Exception:
        LOCAL_TE_CODEC = os.getenv('TE_LOCALTE_CODEC', 'zlib').lower()

    global TIME_SUBSAMPLE_INDICES
    TIME_SUBSAMPLE_INDICES = None
    try:
        probe = _get_series(1)
        N_series = int(probe.size)
    except Exception:
        N_series = None

    # Write time index map for LocalTE consumers (maps LocalIndex -> OriginalIndex)
    try:
        indices = TIME_SUBSAMPLE_INDICES if TIME_SUBSAMPLE_INDICES is not None else np.arange(N_series or 0, dtype=int)
        k_hist = int(args.history_length)
        if indices.size > k_hist:
            local_len = int(indices.size - k_hist)
            df_idx = pd.DataFrame({
                'LocalIndex': np.arange(local_len, dtype=int),
                'OriginalIndex': indices[k_hist:].astype(int),
                'HistoryLength': np.full(local_len, k_hist, dtype=int),
            })
            df_idx.to_parquet('time_index_map.parquet', index=False)
            logging.info(f"Wrote time_index_map.parquet with {local_len} rows (k={k_hist}).")
    except Exception as e:
        logging.warning(f"Failed to write time_index_map.parquet: {e}")


    gene_names_available = bool(ensure_gene_names())

    # Define progress directory. LocalTE payloads are too large to keep in RAM
    # until the final write, so LocalTE runs always stream temporary batches.
    stream_localte_batches = bool(STORE_LOCAL_TE)
    if args.enable_intermediate_save:
        progress_dir = 'TE_progress_parquet'
    elif stream_localte_batches:
        try:
            import uuid as _uuid
            suffix = _uuid.uuid4().hex[:8]
        except Exception:
            suffix = str(int(time.time() * 1000))
        progress_dir = f"TE_progress_localte_stream_{suffix}"
    else:
        progress_dir = None

    # Create progress dir if intermediate saving is enabled (all modes)
    if args.enable_intermediate_save or stream_localte_batches:
        os.makedirs(progress_dir, exist_ok=True)
    if stream_localte_batches and not args.enable_intermediate_save:
        print(f"[TE] LocalTE streaming temp dir: {progress_dir} (batch files cleaned after final parquet)")

    pair_mode = getattr(args, "pair_mode", "default")
    pair_index_info = None

    if pair_mode == "default":
        # Load or initialize list of pairs from CSV/progress directory
        if args.enable_intermediate_save:
            list_pairs = load_progress(args.input_csv, progress_dir)
        else:
            try:
                list_pairs = pd.read_csv(args.input_csv, delimiter=',', header=None).to_numpy().astype(int)
                logging.info(f"Loaded {len(list_pairs)} pairs from {args.input_csv}.")
            except Exception as e:
                logging.error(f"Error loading input CSV {args.input_csv}: {e}")
                raise

        total_pairs = len(list_pairs)
        print(f"[TE] Pairs: {_human_int(total_pairs)} from all_pairs.csv")
        logging.info(f"Total pairs to process: {total_pairs}")

        if total_pairs == 0 and args.enable_intermediate_save:
            print("All pairs have already been processed.")
            logging.info("All pairs have already been processed.")

            # If intermediate files exist, consolidate them into the final result
            all_batches = []
            for fname in sorted(os.listdir(progress_dir)):
                if fname.endswith('.parquet') and fname.startswith('merged_'):
                    try:
                        batch_df = pd.read_parquet(os.path.join(progress_dir, fname))
                        all_batches.append(batch_df)
                        logging.info(f"Loaded batch {fname} for consolidation.")
                    except Exception as e:
                        print(f"Error reading file {fname}: {e}")
                        logging.error(f"Error reading file {fname}: {e}")
                        raise

            if all_batches:
                try:
                    final_results_df = pd.concat(all_batches, ignore_index=True)
                    # output_file not accessible here; keep default filename for intermediate mode
                    final_results_df.to_parquet('TE_result_all.parquet', index=False)
                    print("Final results saved to TE_result_all.parquet.")
                    logging.info("Final results saved to TE_result_all.parquet.")
                except Exception as e:
                    print(f"Error saving final results: {e}")
                    logging.error(f"Error saving final results: {e}")
                    raise

                # Optionally, delete intermediate merged files
                for fname in os.listdir(progress_dir):
                    if fname.endswith('.parquet') and fname.startswith('merged_'):
                        try:
                            os.remove(os.path.join(progress_dir, fname))
                            logging.info(f"Deleted merged file {fname}.")
                        except Exception as e:
                            print(f"Error deleting file {fname}: {e}")
                            logging.error(f"Error deleting file {fname}: {e}")
                print("Intermediate merged progress files have been deleted.")
                logging.info("Intermediate merged progress files have been deleted.")
            else:
                print("No progress files found to combine.")
                logging.info("No progress files found to combine.")
            return
    else:
        # Implicit all-pair modes: derive indices from gene_names and/or matrix rows
        names = ensure_gene_names()
        n_features = GLOBAL_CELL_GENE.shape[0]
        if names and len(names) == n_features:
            all_idx = np.arange(1, n_features + 1, dtype=int)
            if pair_mode == "gene_only":
                peak_names = load_peak_names(names, locate_file("TE_peak_list.txt"))
                mask = np.array([n not in peak_names for n in names], dtype=bool)
                indices = all_idx[mask]
            elif pair_mode == "all_feature":
                indices = all_idx
            else:
                indices = all_idx
        else:
            all_idx = np.arange(1, n_features + 1, dtype=int)
            indices = all_idx
        if indices.size <= 1:
            print("[TE] Not enough features to build all-pair TE.")
            logging.warning("Not enough features to build all-pair TE.")
            return
        pair_index_info = {"indices": indices}
        list_pairs = None
        total_pairs = int(indices.size) * int(indices.size - 1)
        print(f"[TE] Pairs: {_human_int(total_pairs)} implicit {pair_mode}")
        logging.info(f"Total pairs to process (implicit {pair_mode}): {total_pairs}")

    # Define source chunk size per target (smaller for faster first results)
    try:
        batch_size = int(getattr(args, "batch_size", 100) or 100)
        if batch_size <= 0:
            batch_size = 100
    except Exception:
        batch_size = 100  # Adjust based on memory and performance considerations
    batch_size = _maybe_cap_kernel_batch_size(batch_size, int(args.history_length), "kernel")

    # Define buffer size (number of result rows to accumulate before writing)
    if args.results_buffer_rows is None:
        buffer_size = 5000 if STORE_LOCAL_TE else 200000
    else:
        try:
            buffer_size = max(500, int(args.results_buffer_rows))
        except Exception:
            buffer_size = 5000 if STORE_LOCAL_TE else 200000
            logging.warning(f"Invalid --results_buffer_rows value; falling back to {buffer_size}.")
    logging.info(f"Result buffering threshold set to {buffer_size} rows (local TE {'enabled' if STORE_LOCAL_TE else 'disabled'}).")
    print(f"[TE] Output buffer: {_human_int(buffer_size)} rows")
    if STORE_LOCAL_TE:
        local_len = 0
        try:
            local_len = max(0, int(GLOBAL_CELL_GENE.shape[1]) - int(args.history_length))
            if TIME_SUBSAMPLE_INDICES is not None:
                local_len = max(0, int(len(TIME_SUBSAMPLE_INDICES)) - int(args.history_length))
        except Exception:
            local_len = 0
        dtype_bytes = int(np.dtype(LOCAL_TE_DTYPE).itemsize)
        raw_bytes = int(total_pairs) * int(local_len) * dtype_bytes
        print(
            f"[TE] LocalTE payload estimate: {_human_int(total_pairs)} pairs x "
            f"{_human_int(local_len)} timepoints x {LOCAL_TE_DTYPE_STR} "
            f"≈ {_human_bytes(raw_bytes)} before Parquet/zlib compression."
        )
        print("[TE] Note: progress bars track computation; LocalTE compression and disk writes add extra wall time.")

    # Define merge threshold (number of batch Parquet files to trigger a merge)
    merge_threshold = 20  # Adjust based on desired maximum number of files

    # Run exact kernel TE.
    stage1_enable_save = args.enable_intermediate_save or stream_localte_batches
    fast_output = "TE_result_all.parquet"
    direct_localte_exporter = _build_direct_localte_exporter(args, pair_mode, total_pairs)
    if direct_localte_exporter is not None:
        if direct_localte_exporter.get("virtual_array"):
            print("[TE] LocalTE output mode: direct virtual array; TE_result_all.parquet will remain score-only.")
        else:
            print("[TE] LocalTE output mode: direct time-chunk streaming; TE_result_all.parquet will remain score-only.")

    source_blocked_done = False
    kernel_source_block_size = 0
    grouped_state_done = False
    if (
        pair_mode == "default"
        and list_pairs is not None
        and _kernel_grouped_state_enabled("kernel", pair_mode, int(args.history_length))
    ):
        grouped_state_done = run_kernel_grouped_state(
            list_pairs=list_pairs,
            historyLength=int(args.history_length),
            kernel_width=0.5,
            num_cpus=args.num_jobs,
            progress_dir=progress_dir,
            merge_threshold=merge_threshold,
            enable_intermediate_save=stage1_enable_save,
            output_file=fast_output,
            pool_chunksize=max(1, int(POOL_CHUNKSIZE)),
            result_buffer_rows=buffer_size,
            direct_localte_exporter=direct_localte_exporter,
        )

    if (
        not grouped_state_done
        and pair_mode in ("gene_only", "all_feature")
        and pair_index_info is not None
        and "indices" in pair_index_info
        and _kernel_grouped_state_enabled("kernel", pair_mode, int(args.history_length))
    ):
        grouped_state_done = run_kernel_grouped_state_implicit(
            indices=pair_index_info["indices"],
            pair_mode=pair_mode,
            historyLength=int(args.history_length),
            kernel_width=0.5,
            num_cpus=args.num_jobs,
            progress_dir=progress_dir,
            merge_threshold=merge_threshold,
            enable_intermediate_save=stage1_enable_save,
            output_file=fast_output,
            pool_chunksize=max(1, int(POOL_CHUNKSIZE)),
            result_buffer_rows=buffer_size,
            direct_localte_exporter=direct_localte_exporter,
        )

    if (not grouped_state_done) and pair_mode == "default" and list_pairs is not None:
        kernel_source_block_size = _resolve_kernel_source_block_size(list_pairs, int(args.history_length))
    if (
        not grouped_state_done
        and
        pair_mode == "default"
        and list_pairs is not None
        and kernel_source_block_size > 0
    ):
        source_blocked_done = run_kernel_source_blocked(
            list_pairs=list_pairs,
            cell_gene_all=cell_gene_all,
            historyLength=int(args.history_length),
            kernel_width=0.5,
            num_cpus=args.num_jobs,
            batch_size=batch_size,
            progress_dir=progress_dir,
            buffer_size=buffer_size,
            merge_threshold=merge_threshold,
            enable_intermediate_save=stage1_enable_save,
            output_file=fast_output,
            pool_chunksize=max(1, int(POOL_CHUNKSIZE)),
            source_block_size=kernel_source_block_size,
        )

    if (not grouped_state_done) and (not source_blocked_done):
        if pair_mode == "default" and list_pairs is not None and len(list_pairs) > 0:
            source_indices_for_pack = list_pairs[:, 0]
        elif pair_index_info is not None and "indices" in pair_index_info:
            source_indices_for_pack = pair_index_info["indices"]
        else:
            source_indices_for_pack = []
        _prepare_kernel_source_pack_memmap(
            source_indices_for_pack,
            int(args.history_length),
            0.5,
            progress_dir or ".",
        )

    if (not grouped_state_done) and (not source_blocked_done):
        try:
            run_parallel_batches(
                list_pairs=list_pairs,
                cell_gene_all=cell_gene_all,
                historyLength=int(args.history_length),
                kernel_width=0.5,
                num_cpus=args.num_jobs,
                batch_size=batch_size,
                progress_dir=progress_dir,
                buffer_size=buffer_size,
                merge_threshold=merge_threshold,
                enable_intermediate_save=stage1_enable_save,
                mode="kernel",
                output_file=fast_output,
                pair_mode=pair_mode,
                pair_index_info=pair_index_info,
                pool_chunksize=max(1, int(POOL_CHUNKSIZE)),
            )
        finally:
            _cleanup_kernel_source_pack_memmap()

    if stage1_enable_save and not os.path.exists(fast_output):
        try:
            consolidated = consolidate_merged_results(progress_dir, fast_output, delete_after=False)
        except Exception:
            logging.error(f"Failed to consolidate progress files into {fast_output}.")
            raise SystemExit(1)
        if not consolidated:
            logging.error(f"No merged progress files available to build {fast_output}.")
            raise SystemExit(1)
        print(f"Fast-mode results saved to {fast_output}.")
        logging.info(f"Fast-mode results saved to {fast_output}.")

    if _virtual_code_index_enabled() and os.path.exists(fast_output):
        try:
            code_sources, code_targets = _virtual_code_index_feature_sets(pair_mode, list_pairs, pair_index_info)
            _write_virtual_code_index(
                code_sources,
                code_targets,
                int(args.history_length),
                0.5,
                os.getenv("LOCAL_TE_CODE_INDEX_OUTPUT_DIR", LOCAL_TE_CODE_INDEX_OUTPUT_DIR),
            )
        except Exception as exc:
            logging.warning("Virtual LocalTE code index generation failed: %s", exc)
            print(f"[LocalTE] Warning: code index generation failed ({exc}). Virtual LocalTE data is still valid.")

    # Permutation test to build network
    if getattr(args, 'permute', False):
        try:
            df_fast = pd.read_parquet(fast_output)
        except Exception as e:
            logging.error(f"Failed to load fast output {fast_output} for permutation testing: {e}")
            raise SystemExit(1)

        # Select pairs to test
        if args.perm_candidate_grn_fdr and args.perm_candidate_grn_fdr > 0:
            candidate_subset = _select_grn_fdr_candidates(
                df_fast,
                alpha=float(args.perm_candidate_grn_fdr),
                te_cutoff=float(args.perm_candidate_te_cutoff),
            )
            perm_pairs = candidate_subset[['Source', 'Target']].to_numpy(dtype=int)
            logging.info(
                "Permutation candidate selection: grn_fdr=%.4g, te_cutoff=%.4g, pairs=%d.",
                float(args.perm_candidate_grn_fdr),
                float(args.perm_candidate_te_cutoff),
                len(perm_pairs),
            )
            print(
                "[TE] Permutation candidates from GRN FDR: "
                f"alpha={float(args.perm_candidate_grn_fdr):g}, "
                f"te_cutoff={float(args.perm_candidate_te_cutoff):g}, "
                f"pairs={len(perm_pairs)}"
            )
        else:
            perm_pairs = df_fast[['Source', 'Target']].to_numpy(dtype=int)

        logging.info(f"Permutation test (kernel): {len(perm_pairs)} pairs selected.")

        kernel_grouped_perm_requested = (
            _kernel_grouped_state_enabled("kernel", "default", int(args.history_length))
            and len(perm_pairs) > 0
        )
        perm_grouped_source_memmap_path = None
        if (
            kernel_grouped_perm_requested
        ):
            try:
                perm_grouped_dir = os.path.join(progress_dir or ".", "perm_grouped_kernel")
                os.makedirs(perm_grouped_dir, exist_ok=True)
                perm_grouped_source_memmap_path = _prepare_grouped_source_code_memmap(
                    perm_pairs[:, 0],
                    int(args.history_length),
                    0.5,
                    perm_grouped_dir,
                    force=True,
                )
                logging.info(
                    "Kernel permutation using grouped-state exact path%s.",
                    " with source-code memmap" if perm_grouped_source_memmap_path else "",
                )
                print(
                    "[TE] Kernel permutation grouped-state exact path"
                    + (" with source-code memmap" if perm_grouped_source_memmap_path else "")
                )
                if str(KERNEL_PERM_TABLE_SAMPLER).lower() not in ("off", "0", "false", "no"):
                    logging.info(
                        "Kernel permutation table sampler enabled "
                        "(max_cells=%d).",
                        int(KERNEL_PERM_TABLE_MAX_CELLS),
                    )
                    print(
                        "[TE] Kernel permutation table sampler enabled "
                        f"(max_cells={int(KERNEL_PERM_TABLE_MAX_CELLS)})"
                    )
            except Exception as exc:
                logging.warning("Grouped source code memmap disabled for permutation: %s", exc)
                _clear_grouped_source_code_memmap()

        # Group by target and process in chunks
        target_to_sources = {}
        for src, tgt in perm_pairs:
            target_to_sources.setdefault(int(tgt), []).append(int(src))


        work = []
        # Allow user to tune granularity for smoother progress (default as before)
        if getattr(args, 'perm_srcs_per_chunk', None):
            per_chunk = max(1, int(args.perm_srcs_per_chunk))
        elif kernel_grouped_perm_requested and int(PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK) > 0:
            per_chunk = max(1, int(PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK))
        else:
            per_chunk = max(50, 1000 // max(1, args.num_jobs))
        for tgt, srcs in target_to_sources.items():
            for i in range(0, len(srcs), per_chunk):
                work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed), 0.5))

        if (
            kernel_grouped_perm_requested
            and str(PERM_SORT_WORK).lower() not in ("off", "0", "false", "no")
            and work
        ):
            target_cost_cache: dict[int, int] = {}
            source_cost_cache: dict[int, int] = {}

            def kernel_perm_work_cost(item) -> int:
                tgt = int(item[0])
                srcs = item[1]
                if tgt not in target_cost_cache:
                    target_cost_cache[tgt] = max(
                        1,
                        _kernel_grouped_target_state_count(
                            tgt,
                            int(args.history_length),
                            0.5,
                        ),
                    )
                source_cost_sum = 0
                for src in srcs:
                    src = int(src)
                    if src not in source_cost_cache:
                        try:
                            source_cost_cache[src] = max(
                                1,
                                int(
                                    _get_grouped_source_codes_kernel(
                                        src,
                                        int(args.history_length),
                                        0.5,
                                    )[0].size
                                ),
                            )
                        except Exception:
                            source_cost_cache[src] = 1
                    source_cost_sum += source_cost_cache[src]
                return int(target_cost_cache[tgt]) * max(1, int(source_cost_sum))

            ranked_work = sorted(work, key=kernel_perm_work_cost)
            # Start with cheap chunks, then alternate expensive/cheap chunks.  If
            # the first scheduled wave is all hard targets, tqdm can sit at 0%
            # for a long time even though all workers are busy.
            balanced_work = []
            left = 0
            right = len(ranked_work) - 1
            while left <= right:
                balanced_work.append(ranked_work[left])
                left += 1
                if left <= right:
                    balanced_work.append(ranked_work[right])
                    right -= 1
            work = balanced_work
            logging.info(
                "Kernel permutation work balanced by grouped-state complexity "
                "(chunks=%d, sources_per_chunk=%d).",
                len(work),
                per_chunk,
            )
            print(
                "[TE] Kernel permutation work balanced by grouped-state complexity "
                f"(chunks={len(work)}, sources_per_chunk={per_chunk})"
            )

        import multiprocessing
        out_rows = []
        total_chunks = len(work)
        try:
            total_pairs = sum(len(w[1]) for w in work)
        except Exception:
            total_pairs = None
        logging.info(f"Permutation workload: {total_chunks} chunks" + (f", {total_pairs} pairs" if total_pairs is not None else ""))
        try:
            with multiprocessing.Pool(processes=args.num_jobs) as pool:
                # Use small imap chunksize for smoother progress updates (override-able via --perm_imap_chunksize)
                try:
                    cs = int(getattr(args, 'perm_imap_chunksize', 1) or 1)
                except Exception:
                    cs = 1
                if total_pairs is not None:
                    from tqdm import tqdm
                    with tqdm(total=total_chunks, desc="Perm chunks") as pbar_c, tqdm(total=total_pairs, desc="Perm pairs") as pbar_p:
                        worker = _perm_worker_kernel
                        for rows in pool.imap_unordered(worker, work, chunksize=cs):
                            out_rows.extend(rows)
                            pbar_c.update(1)
                            pbar_p.update(len(rows))
                else:
                    from tqdm import tqdm
                    with tqdm(total=total_chunks, desc="Perm chunks") as pbar_c:
                        worker = _perm_worker_kernel
                        for rows in pool.imap_unordered(worker, work, chunksize=cs):
                            out_rows.extend(rows)
                            pbar_c.update(1)
        finally:
            if perm_grouped_source_memmap_path:
                _clear_grouped_source_code_memmap()

        if not out_rows:
            logging.warning("No permutation results produced.")
            return

        if STORE_LOCAL_TE:
            # Include codec column when LocalTE is present
            has_codec = any(len(r) == 8 for r in out_rows)
            if has_codec:
                df_perm = pd.DataFrame(out_rows, columns=['Source', 'Target', 'TE', 'p_value', 'LocalTE_bytes', 'LocalTE_len', 'LocalTE_dtype', 'LocalTE_codec'])
            else:
                df_perm = pd.DataFrame(out_rows, columns=['Source', 'Target', 'TE', 'p_value', 'LocalTE_bytes', 'LocalTE_len', 'LocalTE_dtype'])
        else:
            df_perm = pd.DataFrame(out_rows, columns=['Source', 'Target', 'TE', 'p_value'])
        # Benjamini-Hochberg FDR (monotone, mapped back to original order)
        if args.use_fdr:
            m = len(df_perm)
            pvals = df_perm['p_value'].to_numpy()
            order = np.argsort(pvals)
            ranks = np.empty(m, dtype=int); ranks[order] = np.arange(1, m+1)
            q = pvals * m / ranks
            # enforce monotonicity on the sorted sequence, then map back by 'order'
            q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
            q_final = np.empty_like(q)
            q_final[order] = q_sorted
            df_perm['q_value'] = np.clip(q_final, 0, 1)
        if gene_names_available:
            df_perm.insert(0, 'Source_Name', df_perm['Source'].map(gene_name_from_index))
            df_perm.insert(1, 'Target_Name', df_perm['Target'].map(gene_name_from_index))
        try:
            df_perm.to_parquet(args.perm_output, index=False)
            logging.info(f"Permutation results saved to {args.perm_output}")
        except Exception as e:
            logging.error(f"Failed to write permutation output {args.perm_output}: {e}")
            raise

        # Build network by p-value threshold
        try:
            base_cols = ['Source', 'Target', 'TE', 'p_value']
            if STORE_LOCAL_TE:
                for col in ('LocalTE_bytes', 'LocalTE_len', 'LocalTE_dtype', 'LocalTE_codec'):
                    if col in df_perm.columns:
                        base_cols.append(col)
            if args.use_fdr and 'q_value' in df_perm.columns:
                cols = base_cols + ['q_value']
                edges = df_perm[df_perm['q_value'] < float(args.perm_q_alpha)][cols]
                edge_msg = f"Network edges saved to {args.network_output} (q_alpha={args.perm_q_alpha})."
            else:
                edges = df_perm[df_perm['p_value'] < float(args.perm_alpha)][base_cols]
                edge_msg = f"Network edges saved to {args.network_output} (alpha={args.perm_alpha})."
            if gene_names_available:
                edges.insert(0, 'Source_Name', edges['Source'].map(gene_name_from_index))
                edges.insert(1, 'Target_Name', edges['Target'].map(gene_name_from_index))
            edges.to_parquet(args.network_output, index=False)
            print(edge_msg)
            logging.info(edge_msg)
        except Exception as e:
            logging.error(f"Failed to write network edges: {e}")
            raise
    print(f"[TE] Kernel TE completed in {time.time() - start_time:.2f} seconds.")
    logging.info(f"Calculate TE execution time: {time.time() - start_time:.2f} seconds")

    if stage1_enable_save and (not grouped_state_done) and (not source_blocked_done):
        combine_start_time = time.time()
        try:
            combined = consolidate_merged_results(progress_dir, 'TE_result_all.parquet', delete_after=True)
        except Exception:
            logging.error("Failed to consolidate progress files into final output.")
            raise SystemExit(1)

        if combined:
            print("Final results saved to TE_result_all.parquet.")
            logging.info("Final results saved to TE_result_all.parquet.")
            print("Intermediate merged progress files have been deleted.")
            logging.info("Intermediate merged progress files have been deleted.")
        else:
            print("No progress files found to combine.")
            logging.info("No progress files found to combine.")

        combine_end_time = time.time()
        logging.info(f"Final consolidation took {combine_end_time - combine_start_time:.2f} seconds.")

    if stream_localte_batches and not args.enable_intermediate_save and progress_dir:
        try:
            os.rmdir(progress_dir)
        except OSError:
            # The directory is intentionally left in place if a backend kept
            # non-batch diagnostic files there.
            pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TE analysis with parallel processing and batch management.")
    parser.add_argument('input_csv', type=str, help="Input CSV file containing all pairs.")
    parser.add_argument('num_jobs', type=int, help="Number of parallel jobs.")
    parser.add_argument('history_length', type=str, help="History length (k) for the analysis.")
    parser.add_argument('--enable_intermediate_save', action='store_true', help="Enable intermediate saving of results.")
    parser.add_argument('--mode', choices=['kernel'], default='kernel', help="TE estimator. The default runner is kernel-only.")
    parser.add_argument('--permute', action='store_true', help="Enable permutation testing for kernel mode.")
    parser.add_argument('--perm_imap_chunksize', type=int, default=1, help="imap_unordered chunksize for permutation progress smoothness (default=1 for smooth updates).")
    parser.add_argument('--perm_srcs_per_chunk', type=int, default=None, help="Sources per permutation work unit (smaller => more updates, more overhead).")
    parser.add_argument('--perm_n', type=int, default=100, help="Number of permutations per pair.")
    parser.add_argument('--perm_seed', type=int, default=42, help="Base RNG seed for permutations.")
    parser.add_argument('--perm_candidate_grn_fdr', dest='perm_candidate_grn_fdr', type=float, default=0.0, help="If >0, select permutation candidates by GRN global z-score + BH-FDR at this alpha.")
    parser.add_argument('--perm_candidate_te_cutoff', type=float, default=0.0, help="TE cutoff applied before GRN-FDR permutation candidate selection.")
    parser.add_argument('--perm_alpha', type=float, default=0.01, help="Significance threshold for edges (p-value).")
    parser.add_argument('--perm_output', type=str, default='TE_kernel_perm.parquet', help="Output file for permutation results.")
    parser.add_argument('--network_output', type=str, default='network_edges.parquet', help="Output file for significant network edges.")
    parser.add_argument('--use_fdr', action='store_true', help="Use BH-FDR to threshold permutation results.")
    parser.add_argument('--perm_q_alpha', type=float, default=0.05, help="FDR q-value threshold when --use_fdr is set.")
    parser.add_argument('--store_local_te', action='store_true', help="Store per-timepoint local TE arrays in outputs.")
    parser.add_argument('--localte_codec', choices=['none','zlib'], default=None, help="Compression for LocalTE_bytes column (default: TE_LOCALTE_CODEC or zlib).")
    parser.add_argument('--pool_maxtasks', type=int, default=None, help="Recycle worker processes after N tasks to reduce RSS growth.")
    parser.add_argument('--results_buffer_rows', type=int, default=None, help="Maximum TE rows to buffer in memory before flushing (default: 200000, or 5000 when storing local TE).")
    parser.add_argument('--batch_size', type=int, default=100, help="Sources per target chunk for TE computation (default: 100).")
    parser.add_argument('--species', type=str, default="", help="Species used to resolve TF list for direct LocalTE chunk export.")
    parser.add_argument('--split_type', type=str, default="1", help="Network mode/split type used by Matrix_generate.")
    parser.add_argument(
        '--pair_mode',
        choices=['default', 'gene_only', 'all_feature'],
        default='default',
        help="Pair generation mode: 'default' uses pairs from input_csv; "
             "'gene_only' and 'all_feature' enumerate all pairs implicitly "
             "in a streaming fashion (avoids materialising all_pairs.csv).",
    )
    args = parser.parse_args()
    main(args)
    
