import numpy as np


def _rank_gauss(x: np.ndarray) -> np.ndarray:
    """Rank-normalize to approx. standard normal via inverse CDF.
    Handles ties by average rank; adds small jitter to avoid infs at 0/1.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if n == 0:
        return x
    ranks = np.argsort(np.argsort(x, kind="mergesort"), kind="mergesort").astype(np.float64)
    # average rank for ties via stable sort trick: convert to 1..n
    ranks += 1.0
    u = ranks / (n + 1.0)
    # protect against 0/1
    eps = 1e-6
    u = np.clip(u, eps, 1 - eps)
    from scipy.stats import norm
    z = norm.ppf(u)
    # standardize to unit variance as a precaution
    z -= z.mean()
    std = z.std(ddof=1)
    if std > 0:
        z /= std
    return z


def _sliding_windows_1d(x: np.ndarray, k: int) -> np.ndarray:
    """Return view of shape (n-k, k) with consecutive windows."""
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if k <= 0 or n <= k:
        return np.empty((0, k), dtype=np.float64)
    try:
        from numpy.lib.stride_tricks import sliding_window_view
        return sliding_window_view(x, k)[:-1]
    except Exception:
        # fallback copy
        out = np.empty((n - k, k), dtype=np.float64)
        for i in range(n - k):
            out[i, :] = x[i:i + k]
        return out


def prepare_dest_context_gcmi(dest: np.ndarray, k_hist: int) -> dict:
    """Prepare destination context for GCMI TE.

    - Rank-Gaussianize destination series
    - Build Z = windows of length k_hist for Y past
    - y = destination at t aligned with Z rows
    - Precompute QR factor Q for fast residualization
    """
    dest = np.asarray(dest, dtype=np.float64)
    if k_hist <= 0 or dest.size <= k_hist:
        return {"nobs": 0}
    dest_z = _rank_gauss(dest)
    Z = _sliding_windows_1d(dest_z, k_hist)
    y = dest_z[k_hist:]
    # add intercept column
    Z_aug = np.hstack([Z, np.ones((Z.shape[0], 1), dtype=np.float64)])
    # reduced QR
    Q, R = np.linalg.qr(Z_aug, mode="reduced")
    # precompute y residual
    y_proj = Q @ (Q.T @ y)
    y_resid = y - y_proj
    return {
        "Q": Q,
        "y_resid": y_resid,
        "nobs": int(y_resid.size),
        "k": int(k_hist),
    }


from typing import List, Tuple


def compute_gcmi_te_for_sources(ctx: dict, src_arrays: 'List[np.ndarray]', return_local: bool = False) -> np.ndarray | Tuple[np.ndarray, List[np.ndarray]]:
    """Compute TE (bits) for each source using Gaussian-copula partial correlation.
    TE = -0.5 * log(1 - r^2) / ln(2), where r = corr(resid_x, resid_y)

    If return_local is True, returns (te_vals, locals) where each local array
    is constant per observation such that mean(local) == TE.
    """
    Q = ctx.get("Q")
    y_resid = ctx.get("y_resid")
    k = int(ctx.get("k", 1))
    nobs = int(ctx.get("nobs", 0))
    if nobs <= 0:
        return np.zeros(len(src_arrays), dtype=np.float64)
    te_vals = np.zeros(len(src_arrays), dtype=np.float64)
    locals_list: List[np.ndarray] | None = [] if return_local else None
    for i, s in enumerate(src_arrays):
        s = np.asarray(s, dtype=np.float64)
        if s.size <= k:
            te_vals[i] = 0.0
            continue
        s_z = _rank_gauss(s)
        x = s_z[k - 1:-1]
        # residualize x by Z via Q
        x_proj = Q @ (Q.T @ x)
        x_resid = x - x_proj
        # correlation of residuals
        xr = x_resid - x_resid.mean()
        yr = y_resid - y_resid.mean()
        denom = np.sqrt(np.sum(xr * xr) * np.sum(yr * yr))
        if denom <= 0:
            te_vals[i] = 0.0
        else:
            r = float(np.sum(xr * yr) / denom)
            r2 = max(0.0, min(1.0, r * r))
            # bits
            te_vals[i] = -0.5 * np.log(max(1e-15, 1.0 - r2)) / np.log(2.0)
        if return_local:
            # constant per-sample local contribution whose mean equals TE
            locals_list.append(np.full(nobs, te_vals[i], dtype=np.float64))
    if return_local:
        return te_vals, locals_list  # type: ignore[return-value]
    return te_vals
