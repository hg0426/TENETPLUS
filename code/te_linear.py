import os
import sys
import time
import numpy as np
from typing import Callable

try:
    from numpy.lib.stride_tricks import sliding_window_view as _sliding_window_view
except Exception:
    _sliding_window_view = None


class LinearTEContext:
    def __init__(self, k: int, Q: np.ndarray, r_b: np.ndarray, var_bX: float, nobs: int):
        self.k = int(k)
        self.Q = Q  # (nobs, k) orthonormal columns
        self.r_b = r_b  # (nobs,)
        self.var_bX = float(var_bX)
        self.nobs = int(nobs)


def prepare_dest_context_linear(dest_series: np.ndarray, k: int) -> LinearTEContext | None:
    dest_series = np.asarray(dest_series, dtype=np.float64)
    N = dest_series.shape[0]
    if N <= k:
        return None
    nobs = N - k
    if _sliding_window_view is not None:
        Xp = _sliding_window_view(dest_series, window_shape=k)[:-1]
    else:
        Xp = np.empty((nobs, k), dtype=np.float64)
        for i in range(nobs):
            Xp[i] = dest_series[i:i + k]
    y = dest_series[k:]

    # Orthonormal basis of X_past via thin QR
    Q, _ = np.linalg.qr(Xp, mode='reduced')  # (nobs, k)
    # Residual of y given X_past
    proj_y = Q @ (Q.T @ y)
    r_b = y - proj_y
    var_bX = float(np.mean(r_b * r_b))

    return LinearTEContext(k=k, Q=Q, r_b=r_b, var_bX=var_bX, nobs=nobs)


def _source_past(series: np.ndarray, k: int) -> np.ndarray:
    series = np.asarray(series, dtype=np.float64)
    N = series.shape[0]
    nobs = N - k
    if _sliding_window_view is not None:
        return _sliding_window_view(series, window_shape=k)[:-1]
    X = np.empty((nobs, k), dtype=np.float64)
    for i in range(nobs):
        X[i] = series[i:i + k]
    return X


def _rss_mb() -> float:
    """Return current resident set size in MB (best-effort, Linux-friendly)."""
    # Try /proc/self/statm -> RSS pages * page_size
    try:
        with open('/proc/self/statm', 'r') as f:
            parts = f.read().split()
            if len(parts) >= 2:
                rss_pages = int(parts[1])
                page_size = os.sysconf('SC_PAGE_SIZE') if hasattr(os, 'sysconf') else 4096
                return (rss_pages * page_size) / (1024.0 * 1024.0)
    except Exception:
        pass
    # Fallback to parsing /proc/self/status VmRSS
    try:
        with open('/proc/self/status', 'r') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    parts = line.split()
                    if len(parts) >= 2:
                        kb = float(parts[1])
                        return kb / 1024.0
    except Exception:
        pass
    return 0.0


def compute_linear_te_for_sources(
    ctx: LinearTEContext,
    sources: list[np.ndarray],
    ridge: float = 1e-9,
    batch_size: int = 256,
    return_local: bool = False,
    # Memory control knobs (soft caps; safe defaults)
    max_stack_mb: float = 128.0,
    proj_block_mb: float = 64.0,
    # Optional profiling callback; when None, honors env TE_LINEAR_PROFILE=1 to print progress
    progress_hook: Callable | None = None,
) -> list[float] | tuple[list[float], list[np.ndarray]]:
    """
    Batched linear TE computation using GEMM + einsum to reduce Python loop overhead.
    """
    if ctx is None or ctx.nobs <= 0:
        if return_local:
            S = len(sources)
            return [0.0 for _ in sources], [np.zeros(0, dtype=np.float32) for _ in range(S)]
        return [0.0 for _ in sources]
    Q = ctx.Q
    QT = Q.T  # (k, nobs)
    r_b = ctx.r_b
    var_bX = ctx.var_bX
    nobs = float(ctx.nobs)
    k = ctx.k

    if var_bX <= 0.0:
        if return_local:
            S = len(sources)
            return [0.0 for _ in sources], [np.zeros(ctx.nobs, dtype=np.float32) for _ in range(S)]
        return [0.0 for _ in sources]

    S = len(sources)
    te_values = [0.0] * S
    local_values: list[np.ndarray] | None = [] if return_local else None
    LOG2 = np.log(2.0)

    # Optional profiling set up (no overhead unless enabled)
    prof_env = os.getenv('TE_LINEAR_PROFILE')
    if progress_hook is None and prof_env:
        def _default_hook(ev: str, **info):
            try:
                rss = _rss_mb()
                ts = time.strftime('%H:%M:%S')
                msg = {
                    'event': ev,
                    'rss_mb': round(rss, 2),
                    **{k: (int(v) if isinstance(v, (int, np.integer)) else v) for k, v in info.items()}
                }
                print(f"[TE_LINEAR_PROFILE {ts}] {msg}", file=sys.stderr, flush=True)
            except Exception:
                pass
        progress_hook = _default_hook

    base = 0
    while base < S:
        chunk = sources[base: base + batch_size]
        bs = len(chunk)
        # Build stacked past windows for the chunk: (nobs, k*bs)
        stack = np.empty((ctx.nobs, k * bs), dtype=np.float64)
        for j, src in enumerate(chunk):
            Xs = _source_past(src, k)  # (nobs, k)
            stack[:, j * k:(j + 1) * k] = Xs
        if progress_hook is not None:
            progress_hook('chunk_start', base=base, bs=bs, nobs=ctx.nobs, k=k, stack_cols=stack.shape[1])

        # Project and residualize in one shot for speed
        proj = Q @ (QT @ stack)
        Xr_stack = stack - proj
        Xr = Xr_stack.reshape(ctx.nobs, k, bs)  # (nobs, k, bs)

        # Normal equations for all sources in chunk
        G = np.einsum('nib,njb->ijb', Xr, Xr, optimize=True)  # (k,k,bs)
        yv = np.einsum('nib,n->ib', Xr, r_b, optimize=True)   # (k,bs)

        for j in range(bs):
            Gj = G[:, :, j]
            yvj = yv[:, j]
            if ridge > 0:
                Gj = Gj.copy()
                Gj.flat[::Gj.shape[0] + 1] += ridge
            try:
                beta = np.linalg.solve(Gj, yvj)
            except np.linalg.LinAlgError:
                beta = np.linalg.lstsq(Gj, yvj, rcond=None)[0]
            # Variance reduction and TE
            delta = float(yvj.T @ beta) / nobs
            var_bXY = max(var_bX - delta, 1e-15)
            te = float(0.5 * np.log(var_bX / var_bXY) / LOG2)
            te_values[base + j] = te
            if return_local and local_values is not None:
                # Compute pointwise (local) TE via Gaussian log-density ratio
                # r_b: residuals of y|X; r_f: residuals of y|X,S
                Xrj = Xr[:, :, j]
                r_f = r_b - (Xrj @ beta)
                var_f = max(float(np.mean(r_f * r_f)), 1e-15)
                # local_t = 0.5 * [ ln(σ_b^2/σ_f^2) + (r_f^2/σ_f^2) - (r_b^2/σ_b^2) ] / ln 2
                const = 0.5 * (np.log(var_bX / var_f)) / LOG2
                term = 0.5 * ((r_f * r_f) / var_f - (r_b * r_b) / var_bX) / LOG2
                local = (const + term).astype(np.float32)
                local_values.append(local)
        if progress_hook is not None:
            progress_hook('chunk_done', base=base, consumed=bs)

        base += bs

    if return_local and local_values is not None:
        return te_values, local_values
    return te_values


def _residualize_against_Q(Q: np.ndarray, QT: np.ndarray, X: np.ndarray) -> np.ndarray:
    return X - Q @ (QT @ X)


def _te_from_Xs(ctx: LinearTEContext, Xs: np.ndarray, ridge: float = 1e-9) -> float:
    Q = ctx.Q
    QT = Q.T
    r_b = ctx.r_b
    var_bX = ctx.var_bX
    nobs = float(ctx.nobs)
    LOG2 = np.log(2.0)
    eps = 1e-15
    # If destination residual variance given its own past is ~0, TE is 0
    if var_bX <= eps:
        return 0.0
    Xr = _residualize_against_Q(Q, QT, Xs)
    G = Xr.T @ Xr
    if ridge > 0:
        G.flat[:: G.shape[0] + 1] += ridge
    yv = Xr.T @ r_b
    try:
        beta = np.linalg.solve(G, yv)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(G, yv, rcond=None)[0]
    delta = float(yv.T @ beta) / nobs
    var_bXY = max(var_bX - delta, 1e-15)
    ratio = max(var_bX, 1e-15) / max(var_bXY, 1e-15)
    return float(0.5 * np.log(ratio) / LOG2)


def linear_te_permutation_pvalue(
    ctx: LinearTEContext,
    source_series: np.ndarray,
    num_perms: int = 100,
    seed: int = 0,
    ridge: float = 1e-9,
    return_local: bool = False,
) -> tuple[float, float] | tuple[float, float, np.ndarray]:
    """
    Compute linear TE and a permutation p-value by permuting source past windows (row-shuffle).
    Returns (te_original, p_value).
    """
    if ctx is None or ctx.nobs <= 0:
        if return_local:
            return 0.0, 1.0, np.zeros(0, dtype=np.float32)
        return 0.0, 1.0
    Xs = _source_past(np.asarray(source_series, dtype=np.float64), ctx.k)
    # Compute TE and local contributions for the original ordering
    Q = ctx.Q
    QT = Q.T
    r_b = ctx.r_b
    var_bX = ctx.var_bX
    LOG2 = np.log(2.0)
    Xr = _residualize_against_Q(Q, QT, Xs)
    G = Xr.T @ Xr
    if ridge > 0:
        G.flat[:: G.shape[0] + 1] += ridge
    yv = Xr.T @ r_b
    try:
        beta = np.linalg.solve(G, yv)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(G, yv, rcond=None)[0]
    delta = float(yv.T @ beta) / float(ctx.nobs)
    var_bXY = max(var_bX - delta, 1e-15)
    te_orig = float(0.5 * np.log(max(var_bX, 1e-15) / max(var_bXY, 1e-15)) / LOG2)
    if return_local:
        r_f = r_b - (Xr @ beta)
        var_f = max(float(np.mean(r_f * r_f)), 1e-15)
        const = 0.5 * (np.log(var_bX / var_f)) / LOG2
        term = 0.5 * ((r_f * r_f) / var_f - (r_b * r_b) / var_bX) / LOG2
        local_vec = (const + term).astype(np.float32)
    if num_perms <= 0:
        return (te_orig, 1.0, local_vec) if return_local else (te_orig, 1.0)
    rng = np.random.default_rng(seed)
    n = Xs.shape[0]
    count = 0
    for i in range(num_perms):
        order = rng.permutation(n)
        te_perm = _te_from_Xs(ctx, Xs[order], ridge=ridge)
        if te_perm >= te_orig:
            count += 1
    pval = (count + 1.0) / (num_perms + 1.0)
    if return_local:
        return te_orig, float(pval), local_vec
    return te_orig, float(pval)
