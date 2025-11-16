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
from scipy import sparse
from tqdm import tqdm
from .Calc_TE_python import (
    TransferEntropyCalculatorKernel,
    KernelEstimatorTransferEntropy,
    compute_te_for_target,
    prepare_dest_context,
    compute_te_for_sources_with_context,
)
from .te_linear import (
    prepare_dest_context_linear,
    compute_linear_te_for_sources,
    linear_te_permutation_pvalue,
)
from .te_ksg import (
    prepare_dest_context_ksg,
    compute_ksg_te_for_sources,
)
from .te_poly import (
    prepare_dest_context_poly,
    compute_poly_te_for_sources,
)
from .te_gcmi import (
    prepare_dest_context_gcmi,
    compute_gcmi_te_for_sources,
)
from .te_discrete import (
    prepare_dest_context_discrete_quantile,
    prepare_dest_context_ordinal,
    compute_discrete_te_for_sources,
)
from .te_kernel_grid import (
    prepare_dest_context_kernel_grid,
    compute_kernel_grid_te_for_sources,
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

# cKDTree workers per process (tune via env)
CKD_WORKERS = int(os.getenv("CKD_WORKERS", "1"))
# Toggle storing per-timepoint local TE
STORE_LOCAL_TE = False
LOCAL_TE_DTYPE = np.float16
LOCAL_TE_DTYPE_STR = 'float16'
LOCAL_TE_CODEC = os.getenv('TE_LOCALTE_CODEC', 'zlib').lower()
DENSE_THRESHOLD = int(os.getenv("TE_DENSE_THRESHOLD", "250"))
DISC_BINS = 6
DISC_BIAS = 'miller'
ORDINAL_KX = None

# Configure logging to track the progress and debug if needed
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create file handler which logs even debug messages
fh = logging.FileHandler('TE_analysis.log')
fh.setLevel(logging.INFO)

# Create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# Create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s %(levelname)s:%(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

GLOBAL_CELL_GENE = None
_DEST_CTX_CACHE_KERNEL = {}
_DEST_CTX_CACHE_LINEAR = {}
_DEST_CTX_CACHE_POLY = {}

_DEST_CTX_CACHE_KSG = {}
KSG_K = 4

_DEST_CTX_CACHE_GCMI = {}
_DEST_CTX_CACHE_DISC = {}
_DEST_CTX_CACHE_ORD = {}
_DEST_CTX_CACHE_KGRID = {}

GENE_NAMES = None

# Optional global timepoint subsampling indices (applied to all series)
TIME_SUBSAMPLE_INDICES = None

def _maybe_subsample(series: np.ndarray) -> np.ndarray:
    global TIME_SUBSAMPLE_INDICES
    if TIME_SUBSAMPLE_INDICES is None:
        return series
    try:
        return series[TIME_SUBSAMPLE_INDICES]
    except Exception:
        return series


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
    key = (target_idx, historyLength, float(kernel_width))
    ctx = _DEST_CTX_CACHE_KERNEL.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context(
        dest,
        k=historyLength,
        kernel_width=kernel_width,
        normalise=True,
        precompute_neighbors=False,
        backend="sklearn",
        workers=1,
    )
    if len(_DEST_CTX_CACHE_KERNEL) > 64:
        _DEST_CTX_CACHE_KERNEL.clear()
    _DEST_CTX_CACHE_KERNEL[key] = ctx
    return ctx

_LINEAR_CACHE_MAX = int(os.getenv("DEST_CTX_CACHE_LINEAR_MAX", "64"))

def _get_dest_ctx_linear(target_idx: int, historyLength: int):
    key = (target_idx, historyLength)
    ctx = _DEST_CTX_CACHE_LINEAR.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context_linear(dest, k=historyLength)
    if len(_DEST_CTX_CACHE_LINEAR) > _LINEAR_CACHE_MAX:
        _DEST_CTX_CACHE_LINEAR.clear()
    _DEST_CTX_CACHE_LINEAR[key] = ctx
    return ctx

def _get_dest_ctx_poly(target_idx: int, historyLength: int):
    key = (target_idx, historyLength)
    ctx = _DEST_CTX_CACHE_POLY.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context_poly(dest, k=historyLength)
    if len(_DEST_CTX_CACHE_POLY) > 256:
        _DEST_CTX_CACHE_POLY.clear()
    _DEST_CTX_CACHE_POLY[key] = ctx
    return ctx


def _get_dest_ctx_ksg(target_idx: int, historyLength: int):
    key = (target_idx, historyLength)
    ctx = _DEST_CTX_CACHE_KSG.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context_ksg(dest, k_hist=historyLength)
    if len(_DEST_CTX_CACHE_KSG) > 256:
        _DEST_CTX_CACHE_KSG.clear()
    _DEST_CTX_CACHE_KSG[key] = ctx
    return ctx

def _get_dest_ctx_gcmi(target_idx: int, historyLength: int):
    key = (target_idx, historyLength)
    ctx = _DEST_CTX_CACHE_GCMI.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context_gcmi(dest, k_hist=historyLength)
    if len(_DEST_CTX_CACHE_GCMI) > 256:
        _DEST_CTX_CACHE_GCMI.clear()
    _DEST_CTX_CACHE_GCMI[key] = ctx
    return ctx

def _get_dest_ctx_disc(target_idx: int, historyLength: int, nbins: int):
    key = (target_idx, historyLength, int(nbins))
    ctx = _DEST_CTX_CACHE_DISC.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context_discrete_quantile(dest, k_hist=historyLength, nbins=int(nbins))
    if len(_DEST_CTX_CACHE_DISC) > 256:
        _DEST_CTX_CACHE_DISC.clear()
    _DEST_CTX_CACHE_DISC[key] = ctx
    return ctx

def _get_dest_ctx_ordinal(target_idx: int, historyLength: int, nbins_y: int):
    key = (target_idx, historyLength, int(nbins_y))
    ctx = _DEST_CTX_CACHE_ORD.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    dest = _maybe_subsample(dest)
    ctx = prepare_dest_context_ordinal(dest, k_hist=historyLength, nbins_y=int(nbins_y))
    if len(_DEST_CTX_CACHE_ORD) > 256:
        _DEST_CTX_CACHE_ORD.clear()
    _DEST_CTX_CACHE_ORD[key] = ctx
    return ctx

def _get_dest_ctx_kgrid(target_idx: int, historyLength: int, kernel_width: float):
    key = (target_idx, historyLength, float(kernel_width))
    ctx = _DEST_CTX_CACHE_KGRID.get(key)
    if ctx is not None:
        return ctx
    dest = GLOBAL_CELL_GENE[target_idx - 1].toarray().ravel().astype(np.float64, copy=False)
    ctx = prepare_dest_context_kernel_grid(dest, k=historyLength, kernel_width=kernel_width, normalise=True)
    if len(_DEST_CTX_CACHE_KGRID) > 64:
        _DEST_CTX_CACHE_KGRID.clear()
    _DEST_CTX_CACHE_KGRID[key] = ctx
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
    """Worker: compute TE for a single target and a chunk of its sources."""
    target_idx, source_indices, historyLength, kernel_width, mode = args
    global GLOBAL_CELL_GENE
    global DISC_BINS, DISC_BIAS, ORDINAL_KX
    results = []
    try:
        if mode == 'linear':
            dest_ctx = _get_dest_ctx_linear(target_idx, historyLength)
        elif mode == 'poly':
            dest_ctx = _get_dest_ctx_poly(target_idx, historyLength)
        elif mode == 'ksg':
            dest_ctx = _get_dest_ctx_ksg(target_idx, historyLength)
        elif mode == 'gcmi':
            dest_ctx = _get_dest_ctx_gcmi(target_idx, historyLength)
        elif mode == 'disc':
            dest_ctx = _get_dest_ctx_disc(target_idx, historyLength, int(DISC_BINS))
        elif mode == 'ordinal':
            dest_ctx = _get_dest_ctx_ordinal(target_idx, historyLength, int(DISC_BINS))
        elif mode == 'kernel_grid':
            dest_ctx = _get_dest_ctx_kgrid(target_idx, historyLength, kernel_width)
        else:
            dest_ctx = _get_dest_ctx_kernel(target_idx, historyLength, kernel_width)
    except Exception as e:
        logging.error(f"Failed to build context for target {target_idx}: {e}")
        return results
    # Validate context by mode
    if mode in ('linear', 'poly'):
        if (dest_ctx is None) or getattr(dest_ctx, 'nobs', 0) <= 0:
            return results
    elif mode == 'ksg':
        if not isinstance(dest_ctx, dict) or dest_ctx.get('nobs', 0) <= 0:
            return results
    elif mode in ('gcmi', 'disc', 'ordinal', 'kernel_grid'):
        if not isinstance(dest_ctx, dict):
            return results
        nobs_like = dest_ctx.get('nobs', dest_ctx.get('total_obs', 0))
        if int(nobs_like) <= 0:
            return results
    else:
        if dest_ctx.get('total_obs', 0) <= 0:
            return results

    src_arrays = []
    valid_sources = []
    for sidx in source_indices:
        try:
            s = GLOBAL_CELL_GENE[sidx - 1].toarray().ravel().astype(np.float64, copy=False)
            s = _maybe_subsample(s)
            src_arrays.append(s)
            valid_sources.append(sidx)
        except Exception as e:
            logging.error(f"Failed to load source {sidx}: {e}")
    if not src_arrays:
        return results

    if mode == 'linear':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_linear_te_for_sources(dest_ctx, src_arrays, return_local=True)
        else:
            te_vals = compute_linear_te_for_sources(dest_ctx, src_arrays)
            local_vals = None
    elif mode == 'poly':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_poly_te_for_sources(dest_ctx, src_arrays, return_local=True)
        else:
            te_vals = compute_poly_te_for_sources(dest_ctx, src_arrays)
            local_vals = None
    elif mode == 'ksg':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_ksg_te_for_sources(dest_ctx, src_arrays, historyLength, k_nn=KSG_K, return_local=True)
        else:
            te_vals = compute_ksg_te_for_sources(dest_ctx, src_arrays, historyLength, k_nn=KSG_K)
            local_vals = None
    elif mode == 'gcmi':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_gcmi_te_for_sources(dest_ctx, src_arrays, return_local=True)
        else:
            te_vals = compute_gcmi_te_for_sources(dest_ctx, src_arrays, return_local=False)
            local_vals = None
    elif mode == 'disc':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_discrete_te_for_sources(
                dest_ctx,
                src_arrays,
                method='disc',
                nbins=int(DISC_BINS),
                bias=str(DISC_BIAS),
                return_local=True,
            )
        else:
            te_vals = compute_discrete_te_for_sources(
                dest_ctx,
                src_arrays,
                method='disc',
                nbins=int(DISC_BINS),
                bias=str(DISC_BIAS),
                return_local=False,
            )
            local_vals = None
    elif mode == 'ordinal':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_discrete_te_for_sources(
                dest_ctx,
                src_arrays,
                method='ordinal',
                ord_kx=int(ORDINAL_KX) if ORDINAL_KX is not None else int(historyLength),
                bias=str(DISC_BIAS),
                return_local=True,
            )
        else:
            te_vals = compute_discrete_te_for_sources(
                dest_ctx,
                src_arrays,
                method='ordinal',
                ord_kx=int(ORDINAL_KX) if ORDINAL_KX is not None else int(historyLength),
                bias=str(DISC_BIAS),
                return_local=False,
            )
            local_vals = None
    elif mode == 'kernel_grid':
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_kernel_grid_te_for_sources(dest_ctx, src_arrays, return_local=True)
        else:
            te_vals = compute_kernel_grid_te_for_sources(dest_ctx, src_arrays, return_local=False)
            local_vals = None
    else:
        if STORE_LOCAL_TE:
            te_vals, local_vals = compute_te_for_sources_with_context(
                dest_ctx,
                src_arrays,
                reuse_dest_neighbors=False,
                dense_threshold=DENSE_THRESHOLD,
                return_local=True,
            )
        else:
            te_vals = compute_te_for_sources_with_context(
                dest_ctx,
                src_arrays,
                reuse_dest_neighbors=False,
                dense_threshold=DENSE_THRESHOLD,
            )
            local_vals = None
    if STORE_LOCAL_TE:
        for sidx, te, local in zip(valid_sources, te_vals, local_vals):
            local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local)
            results.append((sidx, target_idx, te, local_bytes, local_len, local_dtype, local_codec))
    else:
        for sidx, te in zip(valid_sources, te_vals):
            results.append((sidx, target_idx, te))

    return results


def _perm_worker_linear(args_w):
    """Top-level worker for linear-mode permutation testing (picklable).
    Args tuple: (target_idx, source_indices, historyLength, perm_n, base_seed)
    Returns list of (Source, Target, TE, p_value)
    """
    tgt, srcs, k_hist, perm_n, base_seed = args_w
    global GLOBAL_CELL_GENE
    res = []
    try:
        ctx = _get_dest_ctx_linear(tgt, k_hist)
    except Exception as e:
        logging.error(f"Failed to build linear context for target {tgt}: {e}")
        return res
    if ctx is None or getattr(ctx, 'nobs', 0) <= 0:
        return res
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            series = _maybe_subsample(series)
            if STORE_LOCAL_TE:
                te, p, local = linear_te_permutation_pvalue(
                    ctx,
                    series,
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                    return_local=True,
                )
                local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local)
                res.append((s, tgt, te, p, local_bytes, local_len, local_dtype, local_codec))
            else:
                te, p = linear_te_permutation_pvalue(
                    ctx,
                    series,
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                )
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation failed for pair ({s},{tgt}): {e}")
    return res


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


def _perm_worker_poly(args_w):
    """Permutation worker for poly mode.
    Args tuple: (target_idx, source_indices, historyLength, perm_n, base_seed)
    Returns list of (Source, Target, TE, p_value)
    """
    from .te_poly import poly_te_permutation_pvalue
    tgt, srcs, k_hist, perm_n, base_seed = args_w
    res = []
    try:
        ctx = _get_dest_ctx_poly(tgt, k_hist)
    except Exception as e:
        logging.error(f"Failed to build poly context for target {tgt}: {e}")
        return res
    if ctx is None or getattr(ctx, 'nobs', 0) <= 0:
        return res
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            series = _maybe_subsample(series)
            if STORE_LOCAL_TE:
                te, p, local = poly_te_permutation_pvalue(
                    ctx,
                    series,
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                    return_local=True,
                )
                local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local)
                res.append((s, tgt, te, p, local_bytes, local_len, local_dtype, local_codec))
            else:
                te, p = poly_te_permutation_pvalue(
                    ctx,
                    series,
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                )
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (poly) failed for pair ({s},{tgt}): {e}")
    return res


def _perm_worker_kernel(args_w):
    """Permutation worker for kernel mode.
    Args tuple: (target_idx, source_indices, historyLength, perm_n, base_seed, kernel_width)
    Returns list of (Source, Target, TE, p_value)
    """
    tgt, srcs, k_hist, perm_n, base_seed, kw = args_w
    res = []
    try:
        ctx = _get_dest_ctx_kernel(tgt, k_hist, kw)
    except Exception as e:
        logging.error(f"Failed to build kernel context for target {tgt}: {e}")
        return res
    if ctx is None or ctx.get('total_obs', 0) <= 0:
        return res
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            series = _maybe_subsample(series)
            if STORE_LOCAL_TE:
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


def _perm_worker_ksg(args_w):
    """Permutation worker for ksg mode.
    Args tuple: (target_idx, source_indices, historyLength, perm_n, base_seed, k_nn)
    Returns list of (Source, Target, TE, p_value)
    """
    from .te_ksg import ksg_te_permutation_pvalue
    tgt, srcs, k_hist, perm_n, base_seed, k_nn = args_w
    res = []
    try:
        ctx = _get_dest_ctx_ksg(tgt, k_hist)
    except Exception as e:
        logging.error(f"Failed to build KSG context for target {tgt}: {e}")
        return res
    if (not isinstance(ctx, dict)) or ctx.get('nobs', 0) <= 0:
        return res
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            series = _maybe_subsample(series)
            if STORE_LOCAL_TE:
                te, p, local = ksg_te_permutation_pvalue(
                    ctx,
                    series,
                    int(k_hist),
                    k_nn=int(k_nn),
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                    return_local=True,
                )
                local_bytes, local_len, local_dtype, local_codec = _encode_local_te(local)
                res.append((s, tgt, te, p, local_bytes, local_len, local_dtype, local_codec))
            else:
                te, p = ksg_te_permutation_pvalue(
                    ctx,
                    series,
                    int(k_hist),
                    k_nn=int(k_nn),
                    num_perms=int(perm_n),
                    seed=int(base_seed) + int(s),
                )
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (ksg) failed for pair ({s},{tgt}): {e}")
    return res


def _kgrid_te_perm_pvalue(dest_ctx: dict, series: np.ndarray, k_hist: int, num_perms: int = 100, seed: int = 0, store_local: bool = False):
    from .te_kernel_grid import compute_kernel_grid_te_for_sources
    x = series
    if store_local:
        te_vals, local_vals = compute_kernel_grid_te_for_sources(dest_ctx, [x], return_local=True)
        te = float(te_vals[0])
    else:
        te_vals = compute_kernel_grid_te_for_sources(dest_ctx, [x], return_local=False)
        te = float(te_vals[0])
        local_vals = None
    rng = np.random.RandomState(int(seed))
    ge = 0
    perms = 0
    x_vec = x.copy()
    for _ in range(int(num_perms)):
        rng.shuffle(x_vec)
        te_p = float(compute_kernel_grid_te_for_sources(dest_ctx, [x_vec], return_local=False)[0])
        if te_p >= te:
            ge += 1
        perms += 1
    pval = (ge + 1.0) / (perms + 1.0)
    if store_local:
        b, n, dt, cc = _encode_local_te(local_vals[0])
        return te, pval, b, n, dt, cc
    return te, pval


def _perm_worker_kgrid(args_w):
    tgt, srcs, k_hist, perm_n, base_seed, kw = args_w
    global GLOBAL_CELL_GENE
    res = []
    try:
        ctx = _get_dest_ctx_kgrid(tgt, k_hist, kw)
    except Exception as e:
        logging.error(f"Failed to build kernel_grid context for target {tgt}: {e}")
        return res
    if ctx is None or ctx.get('total_obs', 0) <= 0:
        return res
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            series = _maybe_subsample(series)
            if STORE_LOCAL_TE:
                te, p, b, n, dt, cc = _kgrid_te_perm_pvalue(ctx, series, k_hist, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=True)
                res.append((s, tgt, te, p, b, n, dt, cc))
            else:
                te, p = _kgrid_te_perm_pvalue(ctx, series, k_hist, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=False)
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (kernel_grid) failed for pair ({s},{tgt}): {e}")
    return res


def _gcmi_te_perm_pvalue(dest_ctx: dict, series: np.ndarray, k_hist: int, num_perms: int = 100, seed: int = 0, store_local: bool = False):
    """Compute GCMI TE and permutation p-value by shuffling source values x = series[k-1:-1]."""
    from .te_gcmi import compute_gcmi_te_for_sources
    x = series
    if store_local:
        te_vals, local_vals = compute_gcmi_te_for_sources(dest_ctx, [x], return_local=True)
        te = float(te_vals[0])
    else:
        te_vals = compute_gcmi_te_for_sources(dest_ctx, [x], return_local=False)
        te = float(te_vals[0])
        local_vals = None
    rng = np.random.RandomState(int(seed))
    perms = 0
    ge = 0
    x_vec = x.copy()
    for _ in range(int(num_perms)):
        rng.shuffle(x_vec)
        te_p = float(compute_gcmi_te_for_sources(dest_ctx, [x_vec], return_local=False)[0])
        if te_p >= te:
            ge += 1
        perms += 1
    pval = (ge + 1.0) / (perms + 1.0)
    if store_local:
        loc = local_vals[0]
        b, n, dt, cc = _encode_local_te(loc)
        return te, pval, b, n, dt, cc
    return te, pval


def _perm_worker_gcmi(args_w):
    tgt, srcs, k_hist, perm_n, base_seed = args_w
    global GLOBAL_CELL_GENE
    res = []
    try:
        ctx = _get_dest_ctx_gcmi(tgt, k_hist)
    except Exception as e:
        logging.error(f"Failed to build GCMI context for target {tgt}: {e}")
        return res
    if ctx is None or ctx.get('nobs', 0) <= 0:
        return res
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            series = _maybe_subsample(series)
            if STORE_LOCAL_TE:
                te, p, b, n, dt, cc = _gcmi_te_perm_pvalue(ctx, series, k_hist, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=True)
                res.append((s, tgt, te, p, b, n, dt, cc))
            else:
                te, p = _gcmi_te_perm_pvalue(ctx, series, k_hist, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=False)
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (gcmi) failed for pair ({s},{tgt}): {e}")
    return res


def _disc_te_perm_pvalue(ctx: dict, x_codes: np.ndarray, nstates_z: int, nstates_x: int, nstates_y: int, num_perms: int = 100, seed: int = 0, store_local: bool = False):
    from .te_discrete import discrete_te, discrete_local_te
    z_codes = ctx['z_codes']
    y_codes = ctx['y_codes']
    te = float(discrete_te(z_codes, x_codes, y_codes, nstates_z, nstates_x, nstates_y))
    if store_local:
        loc = discrete_local_te(z_codes, x_codes, y_codes, nstates_z, nstates_x, nstates_y, alpha=1.0)
    rng = np.random.RandomState(int(seed))
    ge = 0
    perms = 0
    x_perm = x_codes.copy()
    for _ in range(int(num_perms)):
        rng.shuffle(x_perm)
        te_p = float(discrete_te(z_codes, x_perm, y_codes, nstates_z, nstates_x, nstates_y))
        if te_p >= te:
            ge += 1
        perms += 1
    pval = (ge + 1.0) / (perms + 1.0)
    if store_local:
        b, n, dt, cc = _encode_local_te(loc)
        return te, pval, b, n, dt, cc
    return te, pval


def _perm_worker_disc(args_w):
    tgt, srcs, k_hist, perm_n, base_seed, nbins = args_w
    global GLOBAL_CELL_GENE
    res = []
    try:
        ctx = _get_dest_ctx_disc(tgt, k_hist, int(nbins))
    except Exception as e:
        logging.error(f"Failed to build discrete context for target {tgt}: {e}")
        return res
    if ctx is None or ctx.get('nobs', 0) <= 0:
        return res
    nstates_z = int(ctx['nstates_z'])
    nstates_y = int(ctx['nstates_y'])
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            x = series[k_hist - 1:-1]
            from .te_discrete import discretize_quantile
            x_codes = discretize_quantile(x, int(nbins)).astype(np.int32)
            if STORE_LOCAL_TE:
                te, p, b, n, dt, cc = _disc_te_perm_pvalue(ctx, x_codes, nstates_z, int(nbins), nstates_y, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=True)
                res.append((s, tgt, te, p, b, n, dt, cc))
            else:
                te, p = _disc_te_perm_pvalue(ctx, x_codes, nstates_z, int(nbins), nstates_y, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=False)
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (disc) failed for pair ({s},{tgt}): {e}")
    return res


def _ordinal_te_perm_pvalue(ctx: dict, x_codes: np.ndarray, nstates_z: int, nstates_x: int, nstates_y: int, num_perms: int = 100, seed: int = 0, store_local: bool = False):
    from .te_discrete import discrete_te, discrete_local_te
    z_codes = ctx['z_codes']
    y_codes = ctx['y_codes']
    te = float(discrete_te(z_codes, x_codes, y_codes, nstates_z, nstates_x, nstates_y))
    if store_local:
        loc = discrete_local_te(z_codes, x_codes, y_codes, nstates_z, nstates_x, nstates_y, alpha=1.0)
    rng = np.random.RandomState(int(seed))
    ge = 0
    perms = 0
    x_perm = x_codes.copy()
    for _ in range(int(num_perms)):
        rng.shuffle(x_perm)
        te_p = float(discrete_te(z_codes, x_perm, y_codes, nstates_z, nstates_x, nstates_y))
        if te_p >= te:
            ge += 1
        perms += 1
    pval = (ge + 1.0) / (perms + 1.0)
    if store_local:
        b, n, dt, cc = _encode_local_te(loc)
        return te, pval, b, n, dt, cc
    return te, pval


def _perm_worker_ordinal(args_w):
    tgt, srcs, k_hist, perm_n, base_seed, ord_kx, nbins_y = args_w
    global GLOBAL_CELL_GENE
    res = []
    try:
        ctx = _get_dest_ctx_ordinal(tgt, k_hist, int(nbins_y))
    except Exception as e:
        logging.error(f"Failed to build ordinal context for target {tgt}: {e}")
        return res
    if ctx is None or ctx.get('nobs', 0) <= 0:
        return res
    nstates_z = int(ctx['nstates_z'])
    nstates_y = int(ctx['nstates_y'])
    from .te_discrete import ordinal_lehmer_codes
    for s in srcs:
        try:
            series = GLOBAL_CELL_GENE[s - 1].toarray().ravel().astype(np.float64, copy=False)
            # build source ordinal codes of length nobs
            try:
                from numpy.lib.stride_tricks import sliding_window_view
                win = sliding_window_view(series, int(ord_kx))[:-1]
            except Exception:
                n = series.size
                win = np.empty((n - int(ord_kx), int(ord_kx)), dtype=np.float64)
                for r in range(n - int(ord_kx)):
                    win[r, :] = series[r:r + int(ord_kx)]
                win = win[:-1]
            x_codes = ordinal_lehmer_codes(win)
            nstates_x = int(np.math.factorial(int(ord_kx)))
            if STORE_LOCAL_TE:
                te, p, b, n, dt, cc = _ordinal_te_perm_pvalue(ctx, x_codes, nstates_z, nstates_x, nstates_y, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=True)
                res.append((s, tgt, te, p, b, n, dt, cc))
            else:
                te, p = _ordinal_te_perm_pvalue(ctx, x_codes, nstates_z, nstates_x, nstates_y, num_perms=int(perm_n), seed=int(base_seed) + int(s), store_local=False)
                res.append((s, tgt, te, p))
        except Exception as e:
            logging.error(f"Permutation (ordinal) failed for pair ({s},{tgt}): {e}")
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


def select_refine_pairs_duckdb(fast_path: str, topk: int = 0, top_pct: float = 0.0):
    """Select refine pairs from a parquet file using DuckDB without loading the full table into RAM.
    Returns numpy ndarray of shape (N,2) with columns [Source, Target]."""
    con = duckdb.connect()
    try:
        if topk and topk > 0:
            q = f"""
                SELECT Source, Target FROM (
                  SELECT Source, Target, TE,
                         ROW_NUMBER() OVER (PARTITION BY Target ORDER BY TE DESC) AS rn
                  FROM read_parquet('{fast_path}')
                ) WHERE rn <= {int(topk)}
            """
            df = con.execute(q).fetchdf()
            return df[['Source', 'Target']].to_numpy(dtype=int)
        elif top_pct and top_pct > 0.0:
            p = max(min(float(top_pct), 100.0), 0.0) / 100.0
            try:
                q_th = con.execute(f"SELECT quantile_cont(TE, {p}) AS q FROM read_parquet('{fast_path}')").fetchone()[0]
            except Exception:
                q_th = con.execute(f"SELECT approx_quantile(TE, {p}) AS q FROM read_parquet('{fast_path}')").fetchone()[0]
            df = con.execute(f"SELECT Source, Target FROM read_parquet('{fast_path}') WHERE TE >= {q_th}").fetchdf()
            return df[['Source', 'Target']].to_numpy(dtype=int)
        else:
            df = con.execute(f"SELECT Source, Target FROM read_parquet('{fast_path}')").fetchdf()
            return df[['Source', 'Target']].to_numpy(dtype=int)
    finally:
        con.close()


def copy_parquet_duckdb(in_path: str, out_path: str):
    con = duckdb.connect()
    try:
        con.execute(f"COPY (SELECT * FROM read_parquet('{in_path}')) TO '{out_path}' (FORMAT PARQUET)")
    finally:
        con.close()


def _duckdb_cast_for_col(col_name: str) -> str:
    name = col_name.lower()
    if 'bytes' in name:
        return 'BLOB'
    if 'len' in name:
        return 'INTEGER'
    return 'VARCHAR'


def merge_fast_and_refined_duckdb(fast_path: str, refined_path: str, out_path: str):
    """Merge refined scores into fast results on disk with DuckDB.
    Coalesces TE and LocalTE_* columns without loading full tables into pandas.
    """
    con = duckdb.connect()
    try:
        fast_cols = list(con.execute(f"SELECT * FROM read_parquet('{fast_path}') LIMIT 0").fetchdf().columns)
        ref_cols = list(con.execute(f"SELECT * FROM read_parquet('{refined_path}') LIMIT 0").fetchdf().columns)
        local_cols_all = [c for c in ('LocalTE_bytes','LocalTE_len','LocalTE_dtype','LocalTE_codec') if c in set(fast_cols) | set(ref_cols)]

        def select_with_alias(path: str, cols: list[str], prefix: str | None = None):
            sels = ["Source", "Target", ("TE AS TE_r" if prefix == 'r' else "TE")]
            for c in local_cols_all:
                if c in cols:
                    if prefix == 'r':
                        sels.append(f"{c} AS {c}_r")
                    else:
                        sels.append(c)
                else:
                    cast_t = _duckdb_cast_for_col(c)
                    if prefix == 'r':
                        sels.append(f"CAST(NULL AS {cast_t}) AS {c}_r")
                    else:
                        sels.append(f"CAST(NULL AS {cast_t}) AS {c}")
            return f"SELECT {', '.join(sels)} FROM read_parquet('{path}')"

        fast_sub = select_with_alias(fast_path, fast_cols, prefix=None)
        ref_sub = select_with_alias(refined_path, ref_cols, prefix='r')

        assigns = ["COALESCE(r.TE_r, f.TE) AS TE"]
        for c in local_cols_all:
            assigns.append(f"COALESCE(r.{c}_r, f.{c}) AS {c}")
        sql = f"""
            WITH f AS ({fast_sub}), r AS ({ref_sub})
            SELECT f.Source, f.Target, {', '.join(assigns)}
            FROM f LEFT JOIN r ON (f.Source=r.Source AND f.Target=r.Target)
        """
        con.execute(f"COPY ({sql}) TO '{out_path}' (FORMAT PARQUET)")
    finally:
        con.close()

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
    mode='linear',
    output_file='TE_result_all.parquet',
    pair_mode: str = "default",
    pair_index_info: dict | None = None,
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
        # Number of chunks per target: ceil((n-1)/batch_size)
        per_target_chunks = (n - 1 + batch_size - 1) // batch_size
        total_batches = n * per_target_chunks

        def work_iter():
            for tgt in indices:
                chunk = []
                for src in indices:
                    if src == tgt:
                        continue
                    chunk.append(src)
                    if len(chunk) >= batch_size:
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
            for batch_result in pool.imap_unordered(_process_chunk, work_iter(), chunksize=1):
                # Prepare DataFrame for the batch
                batch_records = []
                for row in batch_result:
                    if STORE_LOCAL_TE:
                        source, target, te_value, local_bytes, local_len, local_dtype, local_codec = row
                        record = {
                            'Source': source,
                            'Target': target,
                            'TE': te_value,
                            'LocalTE_bytes': local_bytes,
                            'LocalTE_len': local_len,
                            'LocalTE_dtype': local_dtype,
                            'LocalTE_codec': local_codec,
                        }
                    else:
                        source, target, te_value = row
                        record = {
                            'Source': source,
                            'Target': target,
                            'TE': te_value,
                        }
                    batch_records.append(record)

                batch_df = pd.DataFrame(batch_records)

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
    start_time = time.time()
    print('Starting Calculate_TE')
    logging.info('Starting Calculate_TE')

    # Load expression data
    print('Loading expression data...')
    logging.info('Loading expression data...')
    load_data_start_time = time.time()
    try:
        cell_gene_all = sparse.csr_matrix(pd.read_parquet('cell_gene_trsps.parquet').to_numpy(dtype=np.float64))
        print('Expression data loaded successfully.')
        logging.info('Expression data loaded successfully.')
    except Exception as e:
        print(f"Error loading expression data: {e}")
        logging.error(f"Error loading expression data: {e}")
        return
    load_data_end_time = time.time()
    logging.info(f"Loading expression data took {load_data_end_time - load_data_start_time:.2f} seconds.")

    # Expose matrix as global for worker processes (avoids pickling per task)
    global GLOBAL_CELL_GENE
    GLOBAL_CELL_GENE = cell_gene_all

    global STORE_LOCAL_TE
    STORE_LOCAL_TE = bool(getattr(args, 'store_local_te', False))
    if STORE_LOCAL_TE:
        logging.info("Local TE storage enabled; outputs will include LocalTE arrays.")
    # Configure LocalTE codec
    global LOCAL_TE_CODEC
    try:
        LOCAL_TE_CODEC = str(getattr(args, 'localte_codec', LOCAL_TE_CODEC or 'zlib')).lower()
    except Exception:
        LOCAL_TE_CODEC = os.getenv('TE_LOCALTE_CODEC', 'zlib').lower()
    # Optional profiling for linear TE memory
    try:
        prof_lin = bool(getattr(args, 'profile_linear_mem', False))
    except Exception:
        prof_lin = False
    if prof_lin and getattr(args, 'mode', 'linear') == 'linear':
        os.environ['TE_LINEAR_PROFILE'] = '1'
        logging.info("TE_LINEAR_PROFILE enabled (env set) for detailed memory tracing.")

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
        probe = GLOBAL_CELL_GENE[0].toarray().ravel()
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

    # Set globals for discrete/ordinal configs
    global DISC_BINS, DISC_BIAS, ORDINAL_KX
    DISC_BINS = int(getattr(args, 'disc_bins', 6))
    DISC_BIAS = str(getattr(args, 'disc_bias', 'miller'))
    try:
        ORDINAL_KX = int(args.ordinal_kx) if args.ordinal_kx is not None else None
    except Exception:
        ORDINAL_KX = None

    gene_names_available = bool(ensure_gene_names())

    # Define progress directory
    progress_dir = 'TE_progress_parquet' if args.enable_intermediate_save else None

    # Create progress dir if intermediate saving is enabled (all modes)
    if args.enable_intermediate_save:
        os.makedirs(progress_dir, exist_ok=True)

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
        print(f"Total pairs to process: {total_pairs}")
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
        print(f"Total pairs to process (implicit {pair_mode}): {total_pairs}")
        logging.info(f"Total pairs to process (implicit {pair_mode}): {total_pairs}")

    # Define source chunk size per target (smaller for faster first results)
    try:
        batch_size = int(getattr(args, "batch_size", 100) or 100)
        if batch_size <= 0:
            batch_size = 100
    except Exception:
        batch_size = 100  # Adjust based on memory and performance considerations

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

    # Define merge threshold (number of batch Parquet files to trigger a merge)
    merge_threshold = 20  # Adjust based on desired maximum number of files

    # Run parallel batch processing (fast mode)
    refine_requested = (
        args.mode in ("linear", "poly", "gcmi", "disc", "ordinal", "kernel_grid")
        and (args.hybrid_refine_topk_per_target > 0 or args.hybrid_refine_top_pct > 0.0)
    )
    stage1_enable_save = False if refine_requested else args.enable_intermediate_save
    fast_output = "TE_fast.parquet" if (
        args.mode in ("linear", "poly", "gcmi", "disc", "ordinal", "kernel_grid")
    ) else "TE_result_all.parquet"
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
        mode=args.mode,
        output_file=fast_output,
        pair_mode=pair_mode,
        pair_index_info=pair_index_info,
    )

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

    # Optional hybrid refinement: re-score a subset with kernel/KSG and replace
    if refine_requested:
        # Determine subset to refine using DuckDB to avoid full in-RAM DataFrame
        try:
            refine_pairs = select_refine_pairs_duckdb(
                fast_output,
                int(args.hybrid_refine_topk_per_target) if args.hybrid_refine_topk_per_target else 0,
                float(args.hybrid_refine_top_pct) if args.hybrid_refine_top_pct else 0.0,
            )
        except Exception as e:
            logging.error(f"Failed to select refine pairs via DuckDB: {e}")
            return

        if refine_pairs is None or len(refine_pairs) == 0:
            logging.info("No pairs selected for kernel refinement. Writing fast results as final output.")
            try:
                copy_parquet_duckdb(fast_output, args.hybrid_output)
            except Exception as e:
                logging.error(f"Failed to write final output {args.hybrid_output}: {e}")
            return

        method = getattr(args, "hybrid_refine_method", "kernel")
        logging.info(f"Selected {len(refine_pairs)} pairs for {method} refinement.")

        refine_progress_dir = f"{progress_dir}_{method}"
        os.makedirs(refine_progress_dir, exist_ok=True)
        out_name = f"TE_refined_{method}.parquet"
        run_parallel_batches(
            list_pairs=refine_pairs,
            cell_gene_all=cell_gene_all,
            historyLength=int(args.history_length),
            kernel_width=0.5,
            num_cpus=args.num_jobs,
            batch_size=batch_size,
            progress_dir=refine_progress_dir,
            buffer_size=buffer_size,
            merge_threshold=merge_threshold,
            enable_intermediate_save=True,
            mode=("ksg" if method == "ksg" else "kernel"),
            output_file=out_name,
        )

        # If the refined output file does not exist yet, consolidate
        # any progress files into it (mirrors the fast-mode behaviour).
        if not os.path.exists(out_name):
            try:
                combined_refine = consolidate_merged_results(
                    refine_progress_dir,
                    out_name,
                    delete_after=True,
                )
            except Exception as e:
                logging.error(
                    f"Failed to consolidate refined progress files into {out_name}: {e}"
                )
                combined_refine = False
            if not combined_refine:
                logging.warning(
                    "No refined progress files found; falling back to fast "
                    f"results as final output {args.hybrid_output}."
                )
                try:
                    copy_parquet_duckdb(fast_output, args.hybrid_output)
                    print(
                        f"Final results saved to {args.hybrid_output} "
                        "(fast-only fallback)."
                    )
                    logging.info(
                        f"Final results saved to {args.hybrid_output} "
                        "(fast-only fallback)."
                    )
                except Exception as e2:
                    logging.error(
                        f"Failed to copy fast output {fast_output} "
                        f"to {args.hybrid_output}: {e2}"
                    )
                    raise SystemExit(1)
                return

        # Merge fast + refined using DuckDB to keep memory flat
        try:
            merge_fast_and_refined_duckdb(fast_output, out_name, args.hybrid_output)
            print(f"Final results saved to {args.hybrid_output}.")
            logging.info(f"Final results saved to {args.hybrid_output}.")
        except Exception as e:
            logging.error(
                f"Failed to merge fast and refined outputs via DuckDB: {e}"
            )
            # On merge failure, keep the fast results as the final TE output.
            try:
                copy_parquet_duckdb(fast_output, args.hybrid_output)
                print(
                    f"Final results saved to {args.hybrid_output} "
                    "(fast-only fallback after merge error)."
                )
                logging.info(
                    f"Final results saved to {args.hybrid_output} "
                    "(fast-only fallback after merge error)."
                )
            except Exception as e2:
                logging.error(
                    f"Failed to copy fast output {fast_output} "
                    f"to {args.hybrid_output} after merge error: {e2}"
                )
                raise SystemExit(1)
        return



    # Permutation test to build network
    if (args.mode == 'linear' and args.permute_linear) or getattr(args, 'permute', False):
        try:
            df_fast = pd.read_parquet(fast_output)
        except Exception as e:
            logging.error(f"Failed to load fast output {fast_output} for permutation testing: {e}")
            raise SystemExit(1)

        # Select pairs to test
        if args.permute_topk_per_target and args.permute_topk_per_target > 0:
            k = int(args.permute_topk_per_target)
            df_fast_sorted = df_fast.sort_values(['Target', 'TE'], ascending=[True, False])
            refine_subset = df_fast_sorted.groupby('Target', as_index=False).head(k)
            refine_pairs = refine_subset[['Source', 'Target']].to_numpy(dtype=int)
        elif args.permute_top_pct and args.permute_top_pct > 0.0:
            pct = float(args.permute_top_pct)
            thresh = df_fast['TE'].quantile(max(min(pct, 100.0), 0.0) / 100.0)
            refine_subset = df_fast[df_fast['TE'] >= thresh]
            refine_pairs = refine_subset[['Source', 'Target']].to_numpy(dtype=int)
        else:
            # all pairs
            refine_pairs = df_fast[['Source', 'Target']].to_numpy(dtype=int)

        logging.info(f"Permutation test ({args.mode}): {len(refine_pairs)} pairs selected.")

        # Group by target and process in chunks
        target_to_sources = {}
        for src, tgt in refine_pairs:
            target_to_sources.setdefault(int(tgt), []).append(int(src))


        work = []
        # Allow user to tune granularity for smoother progress (default as before)
        if getattr(args, 'perm_srcs_per_chunk', None):
            per_chunk = max(1, int(args.perm_srcs_per_chunk))
        else:
            per_chunk = max(50, 1000 // max(1, args.num_jobs))
        for tgt, srcs in target_to_sources.items():
            for i in range(0, len(srcs), per_chunk):
                if args.mode == 'linear':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed)))
                elif args.mode == 'poly':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed)))
                elif args.mode == 'kernel':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed), 0.5))
                elif args.mode == 'ksg':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed), int(getattr(args, 'ksg_k', 4))))
                elif args.mode == 'gcmi':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed)))
                elif args.mode == 'disc':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed), int(DISC_BINS)))
                elif args.mode == 'kernel_grid':
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed), 0.5))
                else:  # ordinal
                    kx = int(ORDINAL_KX) if ORDINAL_KX is not None else int(args.history_length)
                    work.append((tgt, srcs[i:i + per_chunk], int(args.history_length), int(args.perm_n), int(args.perm_seed), kx, int(DISC_BINS)))

        import multiprocessing
        out_rows = []
        total_chunks = len(work)
        try:
            total_pairs = sum(len(w[1]) for w in work)
        except Exception:
            total_pairs = None
        logging.info(f"Permutation workload: {total_chunks} chunks" + (f", {total_pairs} pairs" if total_pairs is not None else ""))
        with multiprocessing.Pool(processes=args.num_jobs) as pool:
            # Use small imap chunksize for smoother progress updates (override-able via --perm_imap_chunksize)
            try:
                cs = int(getattr(args, 'perm_imap_chunksize', 1) or 1)
            except Exception:
                cs = 1
            if total_pairs is not None:
                from tqdm import tqdm
                with tqdm(total=total_chunks, desc="Perm chunks") as pbar_c, tqdm(total=total_pairs, desc="Perm pairs") as pbar_p:
                    if args.mode == 'linear':
                        worker = _perm_worker_linear
                    elif args.mode == 'poly':
                        worker = _perm_worker_poly
                    elif args.mode == 'kernel':
                        worker = _perm_worker_kernel
                    elif args.mode == 'ksg':
                        worker = _perm_worker_ksg
                    elif args.mode == 'gcmi':
                        worker = _perm_worker_gcmi
                    elif args.mode == 'disc':
                        worker = _perm_worker_disc
                    elif args.mode == 'kernel_grid':
                        worker = _perm_worker_kgrid
                    else:
                        worker = _perm_worker_ordinal
                    for rows in pool.imap_unordered(worker, work, chunksize=cs):
                        out_rows.extend(rows)
                        pbar_c.update(1)
                        pbar_p.update(len(rows))
            else:
                from tqdm import tqdm
                with tqdm(total=total_chunks, desc="Perm chunks") as pbar_c:
                    if args.mode == 'linear':
                        worker = _perm_worker_linear
                    elif args.mode == 'poly':
                        worker = _perm_worker_poly
                    elif args.mode == 'kernel':
                        worker = _perm_worker_kernel
                    elif args.mode == 'ksg':
                        worker = _perm_worker_ksg
                    elif args.mode == 'gcmi':
                        worker = _perm_worker_gcmi
                    elif args.mode == 'disc':
                        worker = _perm_worker_disc
                    elif args.mode == 'kernel_grid':
                        worker = _perm_worker_kgrid
                    else:
                        worker = _perm_worker_ordinal
                    for rows in pool.imap_unordered(worker, work, chunksize=cs):
                        out_rows.extend(rows)
                        pbar_c.update(1)

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
    print("--- Calculate TE execution time: %s seconds ---" % (time.time() - start_time))
    logging.info(f"Calculate TE execution time: {time.time() - start_time:.2f} seconds")

    if args.enable_intermediate_save:
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TE analysis with parallel processing and batch management.")
    parser.add_argument('input_csv', type=str, help="Input CSV file containing all pairs.")
    parser.add_argument('num_jobs', type=int, help="Number of parallel jobs.")
    parser.add_argument('history_length', type=str, help="History length (k) for the analysis.")
    parser.add_argument('--enable_intermediate_save', action='store_true', help="Enable intermediate saving of results.")
    parser.add_argument('--mode', choices=['linear', 'poly', 'kernel', 'ksg', 'gcmi', 'disc', 'ordinal', 'kernel_grid'], default='linear', help="TE estimator: linear, polynomial, kernel, KSG, GCMI (gcmi), discrete (disc), ordinal, or kernel_grid (grid-hash approximation to kernel).")
    parser.add_argument('--ksg_k', type=int, default=4, help="K for KSG kNN-TE (conditional MI).")
    # Discrete/ordinal options
    parser.add_argument('--disc_bins', type=int, default=6, help="Bins for discrete quantile TE and y discretization in ordinal mode.")
    parser.add_argument('--disc_bias', choices=['none','miller'], default='miller', help="Bias correction for discrete entropies.")
    parser.add_argument('--ordinal_kx', type=int, default=None, help="Embedding dimension for source in ordinal TE (default=history_length).")
    parser.add_argument('--hybrid_refine_topk_per_target', type=int, default=0, help="If >0 and mode in {linear,poly,gcmi,disc,ordinal}, refine top-K per target with selected method and replace scores.")
    parser.add_argument('--hybrid_refine_top_pct', type=float, default=0.0, help="If >0 refine global top percentile (0-100] with selected method and replace scores (for fast modes).")
    parser.add_argument('--hybrid_refine_method', choices=['kernel','ksg'], default='kernel', help="Method used for refinement of top pairs.")
    parser.add_argument('--hybrid_output', type=str, default='TE_result_all.parquet', help="Final output filename after optional refinement.")
    # Linear-mode permutation testing options
    parser.add_argument('--permute_linear', action='store_true', help="Enable permutation testing (linear mode only).")
    parser.add_argument('--permute', action='store_true', help="Enable permutation testing for current mode (linear, poly, kernel, ksg, gcmi, disc, ordinal).")
    parser.add_argument('--perm_imap_chunksize', type=int, default=1, help="imap_unordered chunksize for permutation progress smoothness (default=1 for smooth updates).")
    parser.add_argument('--perm_srcs_per_chunk', type=int, default=None, help="Sources per permutation work unit (smaller => more updates, more overhead).")
    parser.add_argument('--perm_n', type=int, default=100, help="Number of permutations per pair (linear mode).")
    parser.add_argument('--perm_seed', type=int, default=42, help="Base RNG seed for permutations.")
    parser.add_argument('--permute_topk_per_target', type=int, default=0, help="If >0, run permutations for top-K pairs per target by fast TE.")
    parser.add_argument('--permute_top_pct', type=float, default=0.0, help="If >0, run permutations for global top percentile [0-100].")
    parser.add_argument('--perm_alpha', type=float, default=0.01, help="Significance threshold for edges (p-value).")
    parser.add_argument('--perm_output', type=str, default='TE_linear_perm.parquet', help="Output file for permutation results.")
    parser.add_argument('--network_output', type=str, default='network_edges.parquet', help="Output file for significant network edges.")
    parser.add_argument('--use_fdr', action='store_true', help="Use BH-FDR to threshold permutation results.")
    parser.add_argument('--perm_q_alpha', type=float, default=0.05, help="FDR q-value threshold when --use_fdr is set.")
    parser.add_argument('--store_local_te', action='store_true', help="Store per-timepoint local TE arrays in outputs.")
    parser.add_argument('--localte_codec', choices=['none','zlib'], default='zlib', help="Compression for LocalTE_bytes column (default: zlib).")
    parser.add_argument('--profile_linear_mem', action='store_true', help="Enable detailed memory tracing for linear TE (prints to stderr).")
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
    KSG_K = int(getattr(args, "ksg_k", 4))
    main(args)
    
