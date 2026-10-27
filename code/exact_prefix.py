"""Exact k=1 box-kernel counts using sparse corrections and summed volumes.

No quantization: axes contain the original distinct scaled values. Integer
counts are exact. Float scoring deliberately uses the original NumPy expression
and per-pair reduction order. Oversized grids use the original implementation.
"""
import numpy as np
from numba import njit


def bounds(values):
    return (np.searchsorted(values, values - 1., side='left').astype(np.int64),
            np.searchsorted(values, values + 1., side='right').astype(np.int64))


def prepare_sources(matrix, reference):
    reference.GLOBAL_CELL_GENE_DENSE = matrix
    values_list, codes_list = [], []
    positions, exceptions, lefts, rights = [], [], [], []
    offsets, sparse_offsets, defaults = [0], [0], []
    for i in range(matrix.shape[0]):
        values, codes = reference._get_grouped_source_codes_kernel(i + 1, 1, .5)
        values_list.append(values)
        codes_list.append(codes.astype(np.int32))
        default = int(np.argmax(np.bincount(codes)))
        pos = np.flatnonzero(codes != default)
        lo, hi = bounds(values)
        positions.append(pos.astype(np.int32))
        exceptions.append(codes[pos].astype(np.int32))
        lefts.append(lo); rights.append(hi)
        defaults.append(default)
        offsets.append(offsets[-1] + len(values))
        sparse_offsets.append(sparse_offsets[-1] + len(pos))
    packed = (np.array(offsets, dtype=np.int64), np.array(sparse_offsets, dtype=np.int64),
              np.concatenate(positions), np.concatenate(exceptions),
              np.array(defaults, dtype=np.int64), np.concatenate(lefts), np.concatenate(rights))
    return packed, values_list, codes_list


def prepare_target(series, reference):
    past, future = series[:-1], series[1:]
    sp = past / (.5 * (np.std(past, ddof=1) + 1e-9))
    sf = future / (.5 * (np.std(future, ddof=1) + 1e-9))
    codes, pv, fv, counts = reference._target_grouped_state_codes(sp, sf)
    context = reference._target_grouped_state_factor_context(pv, fv, counts)
    pvalues, pindex = np.unique(pv, return_inverse=True)
    fvalues, findex = np.unique(fv, return_inverse=True)
    pl, pr = bounds(pvalues)
    fl, fr = bounds(fvalues)
    target = (codes, pindex, findex, counts, pl, pr, fl, fr, context[3], context[4])
    return target, (codes, pv, fv, counts, context)


@njit(cache=True, inline='always')
def rectangle(v, sl, sr, pl, pr, fl, fr):
    return (v[sr, pr, fr] - v[sl, pr, fr] - v[sr, pl, fr] - v[sr, pr, fl]
            + v[sl, pl, fr] + v[sl, pr, fl] + v[sr, pl, fl] - v[sl, pl, fl])


def count_batch(sources, packed, target, grid_limit=250000, max_unique=8192):
    # Pass concrete defaults to avoid repeated Numba omitted-argument typing.
    return _count_batch(sources, packed, target, grid_limit, max_unique)


@njit(cache=True)
def _count_batch(sources, packed, target, grid_limit, max_unique):
    offsets, sparse_offsets, positions, exceptions, defaults, lefts, rights = packed
    tcodes, pindex, findex, tc, pl, pr, fl, fr, cp, cpn = target
    nt, npast, nfuture = len(tc), len(pl), len(fl)
    capacity = 0
    for src in sources:
        ns = offsets[src + 1] - offsets[src]
        if (ns + 1) * (npast + 1) * (nfuture + 1) <= grid_limit:
            capacity += min(ns * nt, len(tcodes))
    # Separate contiguous columns retain NumPy's original vector arithmetic.
    out = np.empty((5, capacity), dtype=np.int64)
    ends = np.zeros(len(sources) + 1, dtype=np.int64)
    fallback = np.zeros(len(sources), dtype=np.bool_)
    cursor = 0
    for ii in range(len(sources)):
        src = sources[ii]
        begin = offsets[src]
        ns = offsets[src + 1] - begin
        if (ns + 1) * (npast + 1) * (nfuture + 1) > grid_limit:
            fallback[ii] = True
            ends[ii + 1] = cursor
            continue
        # Assume every observation has the modal source state, then correct
        # only exceptions. This includes zeros exactly without visiting them.
        w = np.zeros((ns, nt), dtype=np.int64)
        default = defaults[src]
        w[default, :] = tc
        for j in range(sparse_offsets[src], sparse_offsets[src + 1]):
            target_code = tcodes[positions[j]]
            w[default, target_code] -= 1
            w[exceptions[j], target_code] += 1
        occupied = 0
        for s in range(ns):
            for t in range(nt):
                if w[s, t] > 0:
                    occupied += 1
        if occupied > max_unique:
            fallback[ii] = True
            ends[ii + 1] = cursor
            continue
        volume = np.zeros((ns + 1, npast + 1, nfuture + 1), dtype=np.int64)
        for s in range(ns):
            for t in range(nt):
                volume[s + 1, pindex[t] + 1, findex[t] + 1] = w[s, t]
        for s in range(1, ns + 1):
            for p in range(1, npast + 1):
                for f in range(1, nfuture + 1):
                    volume[s, p, f] += volume[s - 1, p, f]
        for s in range(1, ns + 1):
            for p in range(1, npast + 1):
                for f in range(1, nfuture + 1):
                    volume[s, p, f] += volume[s, p - 1, f]
        for s in range(1, ns + 1):
            for p in range(1, npast + 1):
                for f in range(1, nfuture + 1):
                    volume[s, p, f] += volume[s, p, f - 1]
        # Match baseline np.flatnonzero(W.ravel()) order exactly.
        for s in range(ns):
            sl, sr = lefts[begin + s], rights[begin + s]
            for t in range(nt):
                if w[s, t] == 0:
                    continue
                p, f = pindex[t], findex[t]
                out[0, cursor] = w[s, t]
                out[1, cursor] = rectangle(volume, sl, sr, pl[p], pr[p], 0, nfuture)
                out[2, cursor] = rectangle(volume, sl, sr, pl[p], pr[p], fl[f], fr[f])
                out[3, cursor] = cp[t]
                out[4, cursor] = cpn[t]
                cursor += 1
        ends[ii + 1] = cursor
    return out[:, :cursor], ends, fallback


def scores_from_counts(counts, ends, nobs):
    w, cps, cpns, cp, cpn = counts
    # Each occupied state contains itself, hence the baseline validity mask is
    # all true. Check this rather than silently changing the estimator.
    if np.any(cps <= 0) or np.any(cp <= 0) or np.any(cpn <= 0):
        raise ValueError('Invalid neighborhood count')
    with np.errstate(divide='ignore', invalid='ignore'):
        numerator = cpns / cps
        denominator = cpn / cp
        terms = np.log(numerator / denominator)
        terms = np.where(np.isfinite(terms), terms, 0.)
    weighted = (terms / np.log(2.)) * w
    return np.asarray([float(np.sum(weighted[ends[i]:ends[i+1]]) / nobs)
                       for i in range(len(ends) - 1)])


def score_target(sources, packed, target, reference_context, values, codes, reference,
                 batch_size=128):
    results = []
    n_fallback = 0
    tcodes, pv, fv, tc, context = reference_context
    for start in range(0, len(sources), batch_size):
        batch = np.asarray(sources[start:start + batch_size], dtype=np.int64)
        counts, ends, fallback = count_batch(batch, packed, target)
        scores = scores_from_counts(counts, ends, len(tcodes))
        for i in np.flatnonzero(fallback):
            src = batch[i]
            score = reference._grouped_state_te_from_factorized(
                values[src], codes[src], tcodes, tc, pv, fv, context, max_unique=8192)
            if score is None:
                score = reference._grouped_state_te_from_coded(
                    values[src], codes[src], tcodes, pv, fv, max_unique=0)
            scores[i] = score
        n_fallback += int(fallback.sum())
        results.append(scores)
    return np.concatenate(results), n_fallback
