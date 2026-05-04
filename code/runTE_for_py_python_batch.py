#runTE_for_py_python_batch.py
import os
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
from collections import OrderedDict
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
import datetime
import logging
import duckdb
try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover
    pa = None
    pq = None

# Toggle storing per-timepoint local TE
STORE_LOCAL_TE = False
LOCAL_TE_DTYPE = np.float16
LOCAL_TE_DTYPE_STR = 'float16'
LOCAL_TE_CODEC = os.getenv('TE_LOCALTE_CODEC', 'zlib').lower()
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
    # Auto is aimed at the CD4-like regime where rebuilding T x T source masks
    # dominates both wall time and worker-local memory.
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
    te = float(np.sum(local_values64) / source_codes.size)
    return te, local_values64.astype(np.float32, copy=False)


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
        use_coded = str(KERNEL_GROUPED_STATE_CODED).lower() not in ("off", "0", "false", "no")
        if STORE_LOCAL_TE:
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
                        else:
                            rows.append((int(source_idx), int(target_idx), 0.0))
                        continue
                    te = None
                    if STORE_LOCAL_TE:
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
                if STORE_LOCAL_TE:
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
            GENE_NAMES = [line.strip() for line in fh if line.strip()]
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
        raise ValueError("TENETPLUS_KERNEL runner only supports kernel mode.")
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
    use_grouped_perm = _kernel_grouped_state_enabled("kernel", "default", int(k_hist)) and not STORE_LOCAL_TE
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

    target_to_sources: dict[int, list[int]] = {}
    for src, tgt in list_pairs:
        target_to_sources.setdefault(int(tgt), []).append(int(src))
    work_units = [
        (tgt, srcs, int(historyLength), float(kernel_width))
        for tgt, srcs in target_to_sources.items()
    ]
    logging.info(
        "Using grouped-state exact kernel mode: targets=%d, pairs=%d.",
        len(work_units),
        len(list_pairs),
    )
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

    direct_writer_state = {"writer": None, "rows": 0}
    columns = _grouped_result_columns()

    def flush_grouped_rows(buffered_rows):
        if not buffered_rows:
            return []
        if direct_localte_output:
            _append_grouped_rows_to_parquet_writer(
                buffered_rows,
                columns,
                direct_writer_state,
                output_file,
            )
            return []
        batch_df = pd.DataFrame.from_records(buffered_rows, columns=columns)
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
                _process_grouped_state_target,
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
    columns = _grouped_result_columns()

    def flush_grouped_rows(buffered_rows):
        if not buffered_rows:
            return []
        if direct_localte_output:
            _append_grouped_rows_to_parquet_writer(
                buffered_rows,
                columns,
                direct_writer_state,
                output_file,
            )
            return []
        batch_df = pd.DataFrame.from_records(buffered_rows, columns=columns)
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
    # CD4-like pair lists can be highly skewed by source index. Balance blocks
    # by pair count so one huge source block does not dominate the full run.
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
        return np.array([])

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
            logging.error(f"Error reading file {fname}: {e}. Skipping this file.")

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
        return
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
    if STORE_LOCAL_TE:
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

    # Optional: build global timepoint subsample indices
    global TIME_SUBSAMPLE_INDICES
    TIME_SUBSAMPLE_INDICES = None
    try:
        stride = int(getattr(args, 'time_stride', 1) or 1)
        pct = float(getattr(args, 'time_pct', 100.0) or 100.0)
        seed = int(getattr(args, 'time_seed', 42) or 42)
    except Exception:
        stride, pct, seed = 1, 100.0, 42
    # Determine series length from the first row
    try:
        probe = _get_series(1)
        N_series = int(probe.size)
    except Exception:
        N_series = None
    if N_series and N_series > 0:
        if stride > 1:
            idx = np.arange(0, N_series, stride, dtype=int)
            if idx.size > 0:
                TIME_SUBSAMPLE_INDICES = idx
        elif 0.0 < pct < 100.0:
            m = max(0, min(N_series, int(np.ceil(N_series * (pct / 100.0)))))
            if m > 0:
                rng = np.random.default_rng(seed)
                idx = np.sort(rng.choice(N_series, size=m, replace=False)).astype(int)
                TIME_SUBSAMPLE_INDICES = idx
        if TIME_SUBSAMPLE_INDICES is not None:
            logging.info(f"Time subsampling enabled: using {TIME_SUBSAMPLE_INDICES.size}/{N_series} timepoints (stride={stride}, pct={pct}).")
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
                return np.array([])

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
                    return

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
                mask = np.array([not str(n).startswith("chr") for n in names], dtype=bool)
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
        elif args.permute_topk_per_target and args.permute_topk_per_target > 0:
            k = int(args.permute_topk_per_target)
            df_fast_sorted = df_fast.sort_values(['Target', 'TE'], ascending=[True, False])
            candidate_subset = df_fast_sorted.groupby('Target', as_index=False).head(k)
            perm_pairs = candidate_subset[['Source', 'Target']].to_numpy(dtype=int)
        elif args.permute_top_pct and args.permute_top_pct > 0.0:
            pct = float(args.permute_top_pct)
            thresh = df_fast['TE'].quantile(max(min(pct, 100.0), 0.0) / 100.0)
            candidate_subset = df_fast[df_fast['TE'] >= thresh]
            perm_pairs = candidate_subset[['Source', 'Target']].to_numpy(dtype=int)
        else:
            # all pairs
            perm_pairs = df_fast[['Source', 'Target']].to_numpy(dtype=int)

        logging.info(f"Permutation test (kernel): {len(perm_pairs)} pairs selected.")

        kernel_grouped_perm_requested = (
            _kernel_grouped_state_enabled("kernel", "default", int(args.history_length))
            and not STORE_LOCAL_TE
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
    parser.add_argument('--permute_topk_per_target', type=int, default=0, help="If >0, run permutations for top-K pairs per target by fast TE.")
    parser.add_argument('--permute_top_pct', type=float, default=0.0, help="If >0, run permutations for global top percentile [0-100].")
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
    parser.add_argument(
        '--pair_mode',
        choices=['default', 'gene_only', 'all_feature'],
        default='default',
        help="Pair generation mode: 'default' uses pairs from input_csv; "
             "'gene_only' and 'all_feature' enumerate all pairs implicitly "
             "in a streaming fashion (avoids materialising all_pairs.csv).",
    )
    # Time subsampling options
    parser.add_argument('--time_stride', type=int, default=1, help="Use every N-th timepoint (applied before windowing). Overrides time_pct when >1.")
    parser.add_argument('--time_pct', type=float, default=100.0, help="Randomly sample this percent of timepoints per series (0-100]. Used only when stride==1.")
    parser.add_argument('--time_seed', type=int, default=42, help="RNG seed for time_pct sampling.")

    args = parser.parse_args()
    main(args)
    
