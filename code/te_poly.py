import numpy as np

try:
    from numpy.lib.stride_tricks import sliding_window_view as _sliding_window_view
except Exception:
    _sliding_window_view = None


class PolyTEContext:
    def __init__(self, k: int, Q: np.ndarray, QT: np.ndarray, r_b: np.ndarray,
                 var_bX: float, zXp: np.ndarray, mu_past: np.ndarray, std_past: np.ndarray, nobs: int):
        self.k = int(k)
        self.Q = Q            # (nobs, k)
        self.QT = QT          # (k, nobs)
        self.r_b = r_b        # (nobs,)
        self.var_bX = float(var_bX)
        self.zXp = zXp        # (nobs, k) standardized dest past windows
        self.mu_past = mu_past  # (k,)
        self.std_past = std_past  # (k,)
        self.nobs = int(nobs)


def _windows(series: np.ndarray, k: int) -> np.ndarray:
    series = np.asarray(series, dtype=np.float64)
    if _sliding_window_view is not None:
        return _sliding_window_view(series, window_shape=k)[:-1]
    N = series.shape[0]
    nobs = N - k
    X = np.empty((nobs, k), dtype=np.float64)
    for i in range(nobs):
        X[i] = series[i:i + k]
    return X


def prepare_dest_context_poly(dest_series: np.ndarray, k: int, epsilon: float = 1e-9) -> PolyTEContext | None:
    dest_series = np.asarray(dest_series, dtype=np.float64)
    N = dest_series.shape[0]
    if N <= k:
        return None
    nobs = N - k
    Xp = _windows(dest_series, k)             # (nobs, k)
    y = dest_series[k:]                       # (nobs,)

    # Standardize dest past per lag to mimic kernel normalise
    mu_past = Xp.mean(axis=0)
    std_past = Xp.std(axis=0, ddof=1) + epsilon
    zXp = (Xp - mu_past) / std_past

    # Orthonormal basis of standardized X_past via thin QR
    Q, _ = np.linalg.qr(zXp, mode='reduced')  # (nobs, k)
    QT = Q.T

    # Center y; residual of y given Q
    y_c = y - y.mean()
    proj_y = Q @ (QT @ y_c)
    r_b = y_c - proj_y
    var_bX = float(np.mean(r_b * r_b))

    return PolyTEContext(k=k, Q=Q, QT=QT, r_b=r_b, var_bX=var_bX,
                         zXp=zXp, mu_past=mu_past, std_past=std_past, nobs=nobs)


def compute_poly_te_for_sources(
    ctx: PolyTEContext,
    sources: list[np.ndarray],
    use_square: bool = True,
    use_same_lag_interact: bool = True,
    ridge: float = 1e-8,
    epsilon: float = 1e-9,
    return_local: bool = False,
) -> list[float] | tuple[list[float], list[np.ndarray]]:
    if ctx is None or ctx.nobs <= 0:
        if return_local:
            S = len(sources)
            return [0.0 for _ in sources], [np.zeros(0, dtype=np.float32) for _ in range(S)]
        return [0.0 for _ in sources]

    k = ctx.k
    Q = ctx.Q
    QT = ctx.QT
    r_b = ctx.r_b
    var_bX = ctx.var_bX
    zXp = ctx.zXp
    nobs = float(ctx.nobs)

    if var_bX <= 0.0:
        if return_local:
            S = len(sources)
            return [0.0 for _ in sources], [np.zeros(ctx.nobs, dtype=np.float32) for _ in range(S)]
        return [0.0 for _ in sources]

    te_values: list[float] = []
    local_values: list[np.ndarray] | None = [] if return_local else None
    LOG2 = np.log(2.0)

    for src in sources:
        Xs = _windows(np.asarray(src, dtype=np.float64), k)  # (nobs, k)
        # Standardize source past per lag to match kernel normalise
        mu_s = Xs.mean(axis=0)
        std_s = Xs.std(axis=0, ddof=1) + epsilon
        zXs = (Xs - mu_s) / std_s

        feats = [zXs]
        if use_square:
            feats.append(zXs * zXs)
        if use_same_lag_interact:
            # same-lag interactions between dest and source past
            feats.append(zXp * zXs)
        F = np.hstack(feats)  # (nobs, d) where d = k * (#enabled blocks)

        # Residualize features against Q (remove dest past effect)
        F_proj = Q @ (QT @ F)
        Fr = F - F_proj

        # Normal equations with ridge (d is small: <= 3k)
        G = Fr.T @ Fr
        if ridge > 0:
            # add ridge to diagonal
            diag_idx = np.arange(G.shape[0])
            G[diag_idx, diag_idx] += ridge
        yv = Fr.T @ r_b
        try:
            beta = np.linalg.solve(G, yv)
        except np.linalg.LinAlgError:
            beta = np.linalg.lstsq(G, yv, rcond=None)[0]

        # Variance reduction provided by source features
        delta = float(yv.T @ beta) / nobs
        var_bXY = max(var_bX - delta, 1e-15)
        te = 0.5 * np.log(var_bX / var_bXY) / LOG2
        te_values.append(float(te))
        if return_local and local_values is not None:
            # Local TE via Gaussian log-density ratio using residuals
            r_f = r_b - (Fr @ beta)
            var_f = max(float(np.mean(r_f * r_f)), 1e-15)
            const = 0.5 * (np.log(var_bX / var_f)) / LOG2
            term = 0.5 * ((r_f * r_f) / var_f - (r_b * r_b) / var_bX) / LOG2
            local = (const + term).astype(np.float32)
            local_values.append(local)

    if return_local and local_values is not None:
        return te_values, local_values
    return te_values


def poly_te_permutation_pvalue(
    ctx: PolyTEContext,
    source_series: np.ndarray,
    num_perms: int = 100,
    seed: int = 0,
    use_square: bool = True,
    use_same_lag_interact: bool = True,
    ridge: float = 1e-8,
    epsilon: float = 1e-9,
    return_local: bool = False,
) -> tuple[float, float] | tuple[float, float, np.ndarray]:
    """Compute polynomial-approx TE and a permutation p-value by shuffling source past windows.
    Returns (te_original, p_value).
    """
    import numpy as _np
    if ctx is None or getattr(ctx, 'nobs', 0) <= 0:
        if return_local:
            return 0.0, 1.0, np.zeros(0, dtype=np.float32)
        return 0.0, 1.0
    k = ctx.k
    Q = ctx.Q
    QT = ctx.QT
    r_b = ctx.r_b
    var_bX = ctx.var_bX
    zXp = ctx.zXp
    nobs = int(ctx.nobs)
    LOG2 = _np.log(2.0)
    eps = max(float(epsilon), 1e-15)

    # If baseline residual variance is zero or invalid, TE is zero; avoid divide-by-zero in locals
    if not _np.isfinite(var_bX) or var_bX <= 0.0:
        if return_local:
            return 0.0, 1.0, _np.zeros(nobs, dtype=_np.float32)
        return 0.0, 1.0

    def _te_from_Xs(Xs: _np.ndarray) -> tuple[float, _np.ndarray | None]:
        mu_s = Xs.mean(axis=0)
        std_s = Xs.std(axis=0, ddof=1) + epsilon
        zXs = (Xs - mu_s) / std_s
        feats = [zXs]
        if use_square:
            feats.append(zXs * zXs)
        if use_same_lag_interact:
            feats.append(zXp * zXs)
        F = _np.hstack(feats)
        F_proj = Q @ (QT @ F)
        Fr = F - F_proj
        G = Fr.T @ Fr
        if ridge > 0:
            idx = _np.arange(G.shape[0])
            G[idx, idx] += ridge
        yv = Fr.T @ r_b
        try:
            beta = _np.linalg.solve(G, yv)
        except _np.linalg.LinAlgError:
            beta = _np.linalg.lstsq(G, yv, rcond=None)[0]
        delta = float(yv.T @ beta) / float(nobs)
        var_bXY = max(var_bX - delta, eps)
        te = float(0.5 * _np.log(max(var_bX, eps) / max(var_bXY, eps)) / LOG2)
        if return_local:
            r_f = r_b - (Fr @ beta)
            var_f = max(float(_np.mean(r_f * r_f)), eps)
            const = 0.5 * (_np.log(max(var_bX, eps) / max(var_f, eps))) / LOG2
            term = 0.5 * ((r_f * r_f) / max(var_f, eps) - (r_b * r_b) / max(var_bX, eps)) / LOG2
            local = (const + term).astype(_np.float32)
            return te, local
        return te, None

    Xs = _windows(_np.asarray(source_series, dtype=_np.float64), k)
    te_orig, local_vec = _te_from_Xs(Xs)
    if num_perms <= 0:
        return (te_orig, 1.0, local_vec) if return_local else (te_orig, 1.0)
    rng = _np.random.default_rng(int(seed))
    n = Xs.shape[0]
    count = 0
    for _ in range(int(num_perms)):
        order = rng.permutation(n)
        te_perm, _ = _te_from_Xs(Xs[order])
        if te_perm >= te_orig:
            count += 1
    pval = (count + 1.0) / (num_perms + 1.0)
    if return_local:
        return te_orig, float(pval), local_vec
    return te_orig, float(pval)
