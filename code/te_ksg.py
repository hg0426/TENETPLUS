import numpy as np
from sklearn.neighbors import KDTree
from scipy.special import digamma


def _windows(series: np.ndarray, k: int) -> np.ndarray:
    series = np.asarray(series, dtype=float)
    N = len(series)
    nobs = N - k
    X = np.empty((nobs, k), dtype=float)
    for i in range(nobs):
        X[i] = series[i:i + k]
    return X


def prepare_dest_context_ksg(dest_series: np.ndarray, k_hist: int):
    """Prepare destination-only arrays/trees reused across sources for KSG TE.
    TE = I(X_past; Y_next | Y_past)
    We reuse Z = Y_past and YZ = [Y_next, Y_past].
    """
    y = np.asarray(dest_series, dtype=float)
    N = len(y)
    if N <= k_hist:
        return {"nobs": 0}
    nobs = N - k_hist
    Z = _windows(y, k_hist)               # (nobs, k)
    Ynext = y[k_hist:].reshape(-1, 1)     # (nobs, 1)
    YZ = np.hstack((Ynext, Z))            # (nobs, 1+k)
    return {
        "nobs": nobs,
        "Z": Z,
        "YZ": YZ,
        # Build KDTree lazily in compute to avoid pickling big trees in parent
    }

def _count_ball(tree, pts, radii):
    """Return counts of neighbors within radii for each point in pts (exclude self).
    Uses vector query when available; falls back to per-point loop otherwise.
    """
    try:
        idx_list = tree.query_ball_point(pts, r=radii, p=np.inf)
        return np.fromiter((max(len(a) - 1, 0) for a in idx_list), dtype=int, count=len(pts))
    except TypeError:
        out = np.empty(len(pts), dtype=int)
        for i in range(len(pts)):
            idx = tree.query_ball_point(pts[i], r=float(radii[i]), p=np.inf)
            out[i] = max(len(idx) - 1, 0)
        return out



def compute_ksg_te_for_sources(
    ctx: dict,
    sources: list[np.ndarray],
    k_hist: int,
    k_nn: int = 4,
    return_local: bool = False,
) -> list[float] | tuple[list[float], list[np.ndarray]]:
    """KSG estimator for TE via CMI: I(X;Y|Z). Uses Chebyshev (max) norm.
    For each source, builds trees on XZ and XYZ; reuses Z and YZ from context.
    """
    if ctx.get("nobs", 0) <= 0:
        return [0.0 for _ in sources]

    Z = ctx["Z"]               # (nobs, k)
    YZ = ctx["YZ"]             # (nobs, 1+k)
    nobs = ctx["nobs"]
    k = int(k_hist)

    # Reuse destination-only trees
    tree_Z = KDTree(Z)
    tree_YZ = KDTree(YZ)

    te_vals: list[float] = []
    local_vals: list[np.ndarray] | None = [] if return_local else None
    LOG2 = np.log(2.0)
    digamma_k = digamma(k_nn).item()
    for src in sources:
        X = np.asarray(src, dtype=float)
        if len(X) != nobs + k:
            te_vals.append(0.0)
            continue
        Xp = _windows(X, k)            # (nobs, k)
        XYZ = np.hstack((Xp, YZ))      # (nobs, k + 1 + k) = (nobs, 2k+1)
        XZ = np.hstack((Xp, Z))        # (nobs, 2k)

        tree_XYZ = KDTree(XYZ)
        tree_XZ = KDTree(XZ)

        # Distance to k-th neighbor in joint space (exclude self => ask for k+1, take index k)
        dist, _ = tree_XYZ.query(XYZ, k=k_nn + 1)
        # Use non-negative radii; if k-th distance is 0, open ball radius becomes 0
        # (only self included in marginals; we exclude self below → count 0)
        eps = np.maximum(dist[:, k_nn] - 1e-12, 0.0)

        # Counts within eps in marginals (exclude self) — use helper for per-point radii
        n_z = np.maximum(tree_Z.query_radius(Z, r=eps, count_only=True) - 1, 0)
        n_yz = np.maximum(tree_YZ.query_radius(YZ, r=eps, count_only=True) - 1, 0)
        n_xz = np.maximum(tree_XZ.query_radius(XZ, r=eps, count_only=True) - 1, 0)

        # KSG CMI estimator (Frenzel & Pompe style, +1 correction)
        # Safe digamma on integers >= 1
        nz1 = np.maximum(n_z, 0) + 1
        nxz1 = np.maximum(n_xz, 0) + 1
        nyz1 = np.maximum(n_yz, 0) + 1
        local_nat = digamma_k + (digamma(nz1) - digamma(nxz1) - digamma(nyz1))
        local_bits = local_nat / LOG2
        te_vals.append(float(np.mean(local_bits)))
        if return_local and local_vals is not None:
            local_vals.append(local_bits.astype(np.float32))

    if return_local and local_vals is not None:
        return te_vals, local_vals
    return te_vals


def ksg_te_permutation_pvalue(
    ctx: dict,
    source_series: np.ndarray,
    k_hist: int,
    k_nn: int = 4,
    num_perms: int = 100,
    seed: int = 0,
    return_local: bool = False,
) -> tuple[float, float] | tuple[float, float, np.ndarray]:
    """Compute KSG-TE and a permutation p-value by shuffling source past windows.
    Returns (te_original, p_value).
    """
    import numpy as _np
    from sklearn.neighbors import KDTree as _KDTree
    from scipy.special import digamma as _digamma
    if ctx.get("nobs", 0) <= 0:
        return 0.0, 1.0
    Z = ctx["Z"]
    YZ = ctx["YZ"]
    nobs = int(ctx["nobs"])
    k = int(k_hist)

    tree_Z = _KDTree(Z)
    tree_YZ = _KDTree(YZ)

    X = _np.asarray(source_series, dtype=float)
    if len(X) != nobs + k:
        return 0.0, 1.0
    Xp = _windows(X, k)

    LOG2 = _np.log(2.0)
    digamma_k = _digamma(k_nn).item()

    def _te_from_Xp(Xp_local: _np.ndarray, with_local: bool = False):
        XYZ = _np.hstack((Xp_local, YZ))
        XZ = _np.hstack((Xp_local, Z))
        tree_XYZ = _KDTree(XYZ)
        tree_XZ = _KDTree(XZ)
        dist, _ = tree_XYZ.query(XYZ, k=k_nn + 1)
        eps = _np.maximum(dist[:, k_nn] - 1e-12, 0.0)
        n_z = _np.maximum(tree_Z.query_radius(Z, r=eps, count_only=True) - 1, 0)
        n_yz = _np.maximum(tree_YZ.query_radius(YZ, r=eps, count_only=True) - 1, 0)
        n_xz = _np.maximum(tree_XZ.query_radius(XZ, r=eps, count_only=True) - 1, 0)
        nz1 = _np.maximum(n_z, 0) + 1
        nxz1 = _np.maximum(n_xz, 0) + 1
        nyz1 = _np.maximum(n_yz, 0) + 1
        local_nat = digamma_k + (_digamma(nz1) - _digamma(nxz1) - _digamma(nyz1))
        local_bits = local_nat / LOG2
        te_val = float(_np.mean(local_bits))
        if with_local:
            return te_val, local_bits.astype(_np.float32)
        return te_val

    if return_local:
        te_orig, local_bits = _te_from_Xp(Xp, with_local=True)
    else:
        te_orig = _te_from_Xp(Xp, with_local=False)
    if num_perms <= 0:
        if return_local:
            return te_orig, 1.0, local_bits
        return te_orig, 1.0
    rng = _np.random.default_rng(int(seed))
    n = Xp.shape[0]
    count = 0
    for _ in range(int(num_perms)):
        order = rng.permutation(n)
        te_perm = _te_from_Xp(Xp[order], with_local=False)
        if te_perm >= te_orig:
            count += 1
    pval = (count + 1.0) / (num_perms + 1.0)
    if return_local:
        return te_orig, float(pval), local_bits
    return te_orig, float(pval)
