import numpy as np
from typing import Dict, Tuple, List


def _rint_grid(arr: np.ndarray) -> np.ndarray:
    return np.rint(arr).astype(np.int64, copy=False)


def _adjacent_offsets(dim: int) -> np.ndarray:
    """Return all offsets in {-1,0,1}^dim as (3^dim, dim) array.
    For dim==0, return shape (1,0).
    """
    if dim <= 0:
        return np.zeros((1, 0), dtype=np.int8)
    base = np.array([-1, 0, 1], dtype=np.int8)
    grids = np.meshgrid(*([base] * dim), indexing='ij')
    return np.stack(grids, axis=-1).reshape(-1, dim)


def _count_codes(coords: np.ndarray) -> Dict[Tuple[int, ...], int]:
    counts: Dict[Tuple[int, ...], int] = {}
    for row in coords:
        key = tuple(int(v) for v in row)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _neighbor_count_map(unique_coords: np.ndarray, counts_map: Dict[Tuple[int, ...], int], offsets: np.ndarray) -> Dict[Tuple[int, ...], int]:
    out: Dict[Tuple[int, ...], int] = {}
    for row in unique_coords:
        base = tuple(int(v) for v in row)
        s = 0
        for off in offsets:
            key = tuple(int(b + int(o)) for b, o in zip(base, off))
            s += counts_map.get(key, 0)
        out[base] = s
    return out


def prepare_dest_context_kernel_grid(dest: np.ndarray, k: int, kernel_width: float = 0.5, normalise: bool = True) -> dict:
    dest = np.asarray(dest, dtype=np.float64)
    N = dest.size
    if N <= k:
        return {"total_obs": 0}
    total_obs = N - k
    # Scale as in kernel implementation
    if normalise:
        eps = 1e-9
        from numpy.lib.stride_tricks import sliding_window_view
        past = sliding_window_view(dest, k)[:-1]
        nxt = dest[k:]
        std_p = np.std(past, axis=0, ddof=1)
        std_n = np.std(nxt, ddof=1)
        spast = np.divide(past, kernel_width * (std_p + eps), out=np.zeros_like(past), where=(std_p + eps) != 0)
        snext = np.divide(nxt.reshape(-1, 1), kernel_width * (std_n + eps), out=np.zeros((total_obs, 1)), where=(std_n + eps) != 0)
    else:
        from numpy.lib.stride_tricks import sliding_window_view
        past = sliding_window_view(dest, k)[:-1]
        nxt = dest[k:]
        spast = past / kernel_width
        snext = nxt.reshape(-1, 1) / kernel_width

    gp = _rint_grid(spast)  # (nobs, k)
    gpn = _rint_grid(np.hstack([spast, snext]))  # (nobs, k+1)

    # Build counts and neighbor-count maps for destination-only contexts
    counts_p = _count_codes(gp)
    counts_pn = _count_codes(gpn)
    off_k = _adjacent_offsets(k)
    off_kn = _adjacent_offsets(k + 1)
    unique_gp = np.unique(gp, axis=0)
    unique_gpn = np.unique(gpn, axis=0)
    nbrmap_p = _neighbor_count_map(unique_gp, counts_p, off_k)
    nbrmap_pn = _neighbor_count_map(unique_gpn, counts_pn, off_kn)

    # Resolve per-row neighbor counts
    countPast = np.fromiter((nbrmap_p[tuple(row)] for row in gp), dtype=np.int32, count=gp.shape[0])
    countNextPast = np.fromiter((nbrmap_pn[tuple(row)] for row in gpn), dtype=np.int32, count=gpn.shape[0])

    return {
        "N": N,
        "k": k,
        "kernel_width": float(kernel_width),
        "normalise": bool(normalise),
        "total_obs": int(total_obs),
        "spast": spast,
        "snext": snext,
        "gp": gp,
        "gpn": gpn,
        "countPast": countPast,
        "countNextPast": countNextPast,
        "off_k": off_k,
        "off_kn": off_kn,
    }


def compute_kernel_grid_te_for_sources(ctx: dict, src_arrays: List[np.ndarray], return_local: bool = False):
    if ctx.get("total_obs", 0) <= 0:
        return ([0.0 for _ in src_arrays], [np.zeros((0,), dtype=np.float32) for _ in src_arrays]) if return_local else [0.0 for _ in src_arrays]

    k = int(ctx["k"])
    total_obs = int(ctx["total_obs"])
    spast = ctx["spast"]
    snext = ctx["snext"]
    countPast = ctx["countPast"]
    countNextPast = ctx["countNextPast"]
    gp = ctx["gp"]
    off_k = ctx["off_k"]
    off_kn = ctx["off_kn"]
    LOG2 = np.log(2.0)

    te_vals: List[float] = []
    local_list: List[np.ndarray] = []

    # Precompute integer next grid for concat
    gn = _rint_grid(snext)

    for s in src_arrays:
        s = np.asarray(s, dtype=np.float64)
        if s.size <= k:
            te_vals.append(0.0)
            if return_local:
                local_list.append(np.zeros((0,), dtype=np.float32))
            continue
        x = s[k - 1:-1]
        if ctx["normalise"]:
            eps = 1e-9
            std_x = np.std(x, ddof=1)
            sx = np.divide(x, ctx["kernel_width"] * (std_x + eps), out=np.zeros_like(x), where=(std_x + eps) != 0)
        else:
            sx = x / ctx["kernel_width"]
        gs = _rint_grid(sx.reshape(-1, 1))  # (nobs,1)

        # Build combined grid arrays
        gps = np.hstack([gp, gs])            # (nobs, k+1)
        gpsn = np.hstack([gp, gs, gn])       # (nobs, k+2)

        # Count maps for combined
        cnt_ps = _count_codes(gps)
        cnt_psn = _count_codes(gpsn)
        # neighbor count maps
        uniq_ps = np.unique(gps, axis=0)
        uniq_psn = np.unique(gpsn, axis=0)
        nbr_ps = _neighbor_count_map(uniq_ps, cnt_ps, _adjacent_offsets(k + 1))
        nbr_psn = _neighbor_count_map(uniq_psn, cnt_psn, _adjacent_offsets(k + 2))
        # per-row counts
        c_ps = np.fromiter((nbr_ps[tuple(row)] for row in gps), dtype=np.int32, count=gps.shape[0])
        c_psn = np.fromiter((nbr_psn[tuple(row)] for row in gpsn), dtype=np.int32, count=gpsn.shape[0])

        valid = (c_ps > 0) & (countPast > 0) & (countNextPast > 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            num = c_psn[valid] / c_ps[valid]
            den = countNextPast[valid] / countPast[valid]
            log_terms = np.log(num / den)
            log_terms = np.where(np.isfinite(log_terms), log_terms, 0.0)
        te = float(np.sum(log_terms) / total_obs / LOG2)
        te_vals.append(te)
        if return_local:
            local_bits = np.zeros(total_obs, dtype=np.float32)
            if np.any(valid):
                local_bits[valid] = (log_terms / LOG2).astype(np.float32)
            local_list.append(local_bits)

    if return_local:
        return te_vals, local_list
    return te_vals
