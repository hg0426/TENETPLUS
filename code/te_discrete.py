import numpy as np
from math import factorial
from typing import List, Optional


def discretize_quantile(x: np.ndarray, nbins: int) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if nbins <= 1 or x.size == 0:
        return np.zeros_like(x, dtype=np.int32)
    q = np.linspace(0.0, 1.0, nbins + 1)
    q[0] = 0.0
    q[-1] = 1.0
    edges = np.quantile(x, q)
    # ensure monotonic edges (handle ties)
    edges = np.unique(edges)
    if edges.size <= 1:
        return np.zeros_like(x, dtype=np.int32)
    bins = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, edges.size - 2)
    return bins.astype(np.int32)


def _sliding_windows_1d(x: np.ndarray, k: int) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if k <= 0 or n <= k:
        return np.empty((0, k), dtype=np.float64)
    try:
        from numpy.lib.stride_tricks import sliding_window_view
        return sliding_window_view(x, k)[:-1]
    except Exception:
        out = np.empty((n - k, k), dtype=np.float64)
        for i in range(n - k):
            out[i, :] = x[i:i + k]
        return out


def ordinal_lehmer_codes(windows: np.ndarray) -> np.ndarray:
    """Compute ordinal pattern Lehmer codes for each row window.
    Ties are broken by stable ranking (order of appearance).
    """
    if windows.size == 0:
        return np.zeros((0,), dtype=np.int32)
    n, k = windows.shape
    # stable argsort twice to get ranks with ties respecting order
    order = np.argsort(windows, axis=1, kind="mergesort")
    ranks = np.argsort(order, axis=1, kind="mergesort")
    # Compute Lehmer code per row
    codes = np.zeros(n, dtype=np.int64)
    # precompute factorial weights
    weights = np.array([factorial(k - 1 - j) for j in range(k)], dtype=np.int64)
    # convert to Lehmer digits by counting smaller elements to the right
    for i in range(n):
        row = ranks[i]
        # O(k^2) but k is small
        csum = 0
        for j in range(k):
            val = row[j]
            # count how many to the right are < val
            cnt = np.count_nonzero(row[j + 1:] < val)
            csum += cnt * int(weights[j])
        codes[i] = csum
    return codes.astype(np.int32)


def _combine_codes(z: np.ndarray, x: np.ndarray, y: np.ndarray, bases: tuple[int, int, int]) -> np.ndarray:
    bz, bx, by = bases
    return (z.astype(np.int64) * (bx * by) + x.astype(np.int64) * by + y.astype(np.int64)).astype(np.int64)


def discrete_te(z_codes: np.ndarray, x_codes: np.ndarray, y_codes: np.ndarray,
                nstates_z: int, nstates_x: int, nstates_y: int,
                bias: str = "miller") -> float:
    """Compute TE = I(X;Y|Z) (bits) from discrete codes.
    Uses unique/counts; optional Miller–Madow small-sample bias for entropies.
    """
    z = z_codes.astype(np.int64)
    x = x_codes.astype(np.int64)
    y = y_codes.astype(np.int64)
    n = z.size
    if n == 0:
        return 0.0

    # Helper to entropy from counts
    def H_from_counts(counts: np.ndarray) -> float:
        p = counts.astype(np.float64)
        p = p[p > 0]
        p /= n
        H = -np.sum(p * (np.log(p) / np.log(2.0)))
        if bias == "miller":
            # k = number of non-empty bins
            k = p.size
            H += (k - 1) / (2.0 * n * np.log(2.0))
        return float(H)

    # Entropies: H(Y,Z), H(Z), H(X,Y,Z), H(X,Z)
    # Use combined codes to count efficiently
    yz = (z * nstates_y + y)
    xz = (z * nstates_x + x)
    xyz = (z * (nstates_x * nstates_y) + x * nstates_y + y)

    # counts via unique
    def counts_of(arr: np.ndarray) -> np.ndarray:
        _, c = np.unique(arr, return_counts=True)
        return c

    H_yz = H_from_counts(counts_of(yz))
    H_z = H_from_counts(counts_of(z))
    H_xyz = H_from_counts(counts_of(xyz))
    H_xz = H_from_counts(counts_of(xz))

    # TE = H(Y|Z) - H(Y|X,Z) = H(Y,Z) - H(Z) - H(X,Y,Z) + H(X,Z)
    te = H_yz - H_z - H_xyz + H_xz
    return max(0.0, te)


def discrete_local_te(z_codes: np.ndarray, x_codes: np.ndarray, y_codes: np.ndarray,
                      nstates_z: int, nstates_x: int, nstates_y: int,
                      alpha: float = 1.0) -> np.ndarray:
    """Per-sample local TE = log2 p(y|x,z) - log2 p(y|z) with Laplace smoothing alpha.
    Returns array length nobs.
    """
    z = z_codes.astype(np.int64)
    x = x_codes.astype(np.int64)
    y = y_codes.astype(np.int64)
    n = z.size
    if n == 0:
        return np.zeros((0,), dtype=np.float64)

    yz = (z * nstates_y + y)
    xz = (z * nstates_x + x)
    xyz = (z * (nstates_x * nstates_y) + x * nstates_y + y)

    # counts
    def counts_indexed(max_state: int, codes: np.ndarray) -> np.ndarray:
        cnt = np.bincount(codes, minlength=max_state)
        return cnt.astype(np.float64)

    counts_z = counts_indexed(nstates_z, z)
    counts_yz = counts_indexed(nstates_z * nstates_y, yz)
    counts_xz = counts_indexed(nstates_z * nstates_x, xz)
    counts_xyz = counts_indexed(nstates_z * nstates_x * nstates_y, xyz)

    # local probabilities with Laplace smoothing across Y categories
    # p(y|z) ~ (c_yz + alpha)/(c_z + alpha*nstates_y)
    p_y_given_z = (counts_yz[yz] + alpha) / (counts_z[z] + alpha * nstates_y)
    p_y_given_xz = (counts_xyz[xyz] + alpha) / (counts_xz[xz] + alpha * nstates_y)
    with np.errstate(divide='ignore', invalid='ignore'):
        loc = np.log2(np.maximum(1e-300, p_y_given_xz)) - np.log2(np.maximum(1e-300, p_y_given_z))
    return loc.astype(np.float64)


def prepare_dest_context_discrete_quantile(dest: np.ndarray, k_hist: int, nbins: int) -> dict:
    dest = np.asarray(dest, dtype=np.float64)
    if k_hist <= 0 or dest.size <= k_hist:
        return {"nobs": 0}
    Z = _sliding_windows_1d(dest, k_hist)
    y = dest[k_hist:]
    # discretize per column then encode base-nbins
    z_codes_cols = [discretize_quantile(Z[:, j], nbins) for j in range(Z.shape[1])]
    # base encoding
    z_codes = np.zeros(Z.shape[0], dtype=np.int64)
    base = 1
    for j in range(Z.shape[1] - 1, -1, -1):
        z_codes += z_codes_cols[j].astype(np.int64) * base
        base *= nbins
    y_codes = discretize_quantile(y, nbins).astype(np.int32)
    return {
        "z_codes": z_codes.astype(np.int32),
        "y_codes": y_codes,
        "nstates_z": int(base),
        "nstates_x": int(nbins),
        "nstates_y": int(nbins),
        "nobs": int(y_codes.size),
        "k": int(k_hist),
        "nbins": int(nbins),
    }


def prepare_dest_context_ordinal(dest: np.ndarray, k_hist: int, nbins_y: int) -> dict:
    dest = np.asarray(dest, dtype=np.float64)
    if k_hist <= 0 or dest.size <= k_hist:
        return {"nobs": 0}
    Z = _sliding_windows_1d(dest, k_hist)
    y = dest[k_hist:]
    z_codes = ordinal_lehmer_codes(Z)
    y_codes = discretize_quantile(y, nbins_y).astype(np.int32)
    # number of permutations k!
    nstates_z = int(factorial(k_hist))
    return {
        "z_codes": z_codes.astype(np.int32),
        "y_codes": y_codes,
        "nstates_z": nstates_z,
        "nstates_x": nstates_z,  # default if source also uses same k to ordinal-encode
        "nstates_y": int(nbins_y),
        "nobs": int(y_codes.size),
        "k": int(k_hist),
        "nbins_y": int(nbins_y),
    }


def compute_discrete_te_for_sources(ctx: dict, src_arrays: 'List[np.ndarray]', *, method: str,
                                    nbins: int = 6, ord_kx: Optional[int] = None,
                                    bias: str = "miller", return_local: bool = False) -> np.ndarray | tuple[np.ndarray, List[np.ndarray]]:
    z_codes = ctx.get("z_codes")
    y_codes = ctx.get("y_codes")
    nobs = int(ctx.get("nobs", 0))
    k = int(ctx.get("k", 1))
    if nobs <= 0:
        return np.zeros(len(src_arrays), dtype=np.float64)

    te_vals = np.zeros(len(src_arrays), dtype=np.float64)
    locals_list: List[np.ndarray] | None = [] if return_local else None
    if method == "disc":
        nstates_z = int(ctx.get("nstates_z"))
        nstates_y = int(ctx.get("nstates_y"))
        for i, s in enumerate(src_arrays):
            s = np.asarray(s, dtype=np.float64)
            if s.size <= k:
                te_vals[i] = 0.0
                continue
            x = s[k - 1:-1]
            x_codes = discretize_quantile(x, nbins)
            te_vals[i] = discrete_te(z_codes, x_codes, y_codes, nstates_z, nbins, nstates_y, bias=bias)
            if return_local:
                locals_list.append(discrete_local_te(z_codes, x_codes, y_codes, nstates_z, nbins, nstates_y, alpha=1.0))
    else:  # ordinal
        if ord_kx is None:
            ord_kx = k
        nstates_z = int(ctx.get("nstates_z"))  # = k!
        nstates_y = int(ctx.get("nstates_y"))
        for i, s in enumerate(src_arrays):
            s = np.asarray(s, dtype=np.float64)
            if s.size <= ord_kx:
                te_vals[i] = 0.0
                continue
            # source windows end at t-1 -> use same slicing as dest windows alignment
            # Build windows of length ord_kx and drop the last row to align with y
            try:
                from numpy.lib.stride_tricks import sliding_window_view
                win = sliding_window_view(s, ord_kx)[:-1]
            except Exception:
                n = s.size
                win = np.empty((n - ord_kx, ord_kx), dtype=np.float64)
                for r in range(n - ord_kx):
                    win[r, :] = s[r:r + ord_kx]
                win = win[:-1]
            x_codes = ordinal_lehmer_codes(win)
            # If ord_kx != k, nstates_x = ord_kx!
            nstates_x = int(factorial(ord_kx))
            te_vals[i] = discrete_te(z_codes, x_codes, y_codes, nstates_z, nstates_x, nstates_y, bias=bias)
            if return_local:
                locals_list.append(discrete_local_te(z_codes, x_codes, y_codes, nstates_z, nstates_x, nstates_y, alpha=1.0))
    if return_local:
        return te_vals, locals_list  # type: ignore[return-value]
    return te_vals
