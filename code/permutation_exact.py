"""Same-seed source permutations reused across targets; exact batched counts.

Reduction follows NumPy float64 pairwise order (block 128, eight lanes).
Algorithm reference: numpy/core/src/umath/loops_utils.h.src, NumPy 1.26.4.
No fastmath; logarithms and elementwise arithmetic remain NumPy operations.
"""
import numpy as np
from numba import njit
from . import perm_counts as counts_engine


@njit(cache=True)
def pairwise(a, lo, n):
    if n < 8:
        total = -0.0
        for j in range(n):
            total += a[lo+j]
        return total
    if n <= 128:
        r0,r1,r2,r3=a[lo],a[lo+1],a[lo+2],a[lo+3]
        r4,r5,r6,r7=a[lo+4],a[lo+5],a[lo+6],a[lo+7]
        stop=n-n%8
        for j in range(8,stop,8):
            r0+=a[lo+j];r1+=a[lo+j+1];r2+=a[lo+j+2];r3+=a[lo+j+3]
            r4+=a[lo+j+4];r5+=a[lo+j+5];r6+=a[lo+j+6];r7+=a[lo+j+7]
        total=((r0+r1)+(r2+r3))+((r4+r5)+(r6+r7))
        for j in range(stop,n):
            total+=a[lo+j]
        return total
    half=n//2
    half-=half%8
    return pairwise(a,lo,half)+pairwise(a,lo+half,n-half)


@njit(cache=True)
def reduce_segments(weighted, ends, nobs):
    result=np.empty(len(ends)-1,dtype=np.float64)
    for i in range(len(result)):
        # NumPy reduce starts with the additive identity +0.
        result[i]=(0.0+pairwise(weighted,ends[i],ends[i+1]-ends[i]))/nobs
    return result


def scores_from_counts(counts, ends, nobs):
    w,cps,cpns,cp,cpn=counts
    if np.any(cps<=0) or np.any(cp<=0) or np.any(cpn<=0):
        raise ValueError('Invalid neighborhood count')
    with np.errstate(divide='ignore',invalid='ignore'):
        numerator=cpns/cps
        denominator=cpn/cp
        terms=np.log(numerator/denominator)
        terms=np.where(np.isfinite(terms),terms,0.)
    weighted=(terms/np.log(2.))*w
    return reduce_segments(weighted,ends,nobs)


def prepare_permutations(values, codes, source_id, nperm, base_seed=42):
    """Generate EXACTLY rng(seed+source).permutation(n) in original order."""
    rng=np.random.default_rng(int(base_seed)+int(source_id))
    n=len(codes);ns=len(values)
    all_codes=np.empty((nperm,n),dtype=np.int32)
    hist=np.bincount(codes,minlength=ns).astype(np.int64)
    default=int(np.argmax(hist))
    nonmodal=n-int(hist[default])
    positions=np.empty(nperm*nonmodal,dtype=np.int32)
    exceptions=np.empty_like(positions)
    for b in range(nperm):
        all_codes[b]=codes[rng.permutation(n)]
        pos=np.flatnonzero(all_codes[b]!=default)
        positions[b*nonmodal:(b+1)*nonmodal]=pos
        exceptions[b*nonmodal:(b+1)*nonmodal]=all_codes[b,pos]
    left,right=counts_engine.bounds(values)
    return (np.arange(nperm+1,dtype=np.int64)*ns,
            np.arange(nperm+1,dtype=np.int64)*nonmodal,positions,exceptions,
            np.full(nperm,default,dtype=np.int64),np.tile(left,nperm),np.tile(right,nperm),
            all_codes,np.tile(hist,nperm),
            np.full(nperm,counts_engine.isolated(values,left,right),dtype=np.bool_))


def score_permutations(packed, values, target, context, reference, batch_size=256):
    nperm=len(packed[4]);result=np.empty(nperm,dtype=np.float64)
    tc,pv,fv,tw,factor=context
    for begin in range(0,nperm,batch_size):
        ids=np.arange(begin,min(begin+batch_size,nperm),dtype=np.int64)
        counts,ends,fallback=counts_engine.count_batch(ids,packed,target)
        scores=scores_from_counts(counts,ends,len(tc))
        for j in np.flatnonzero(fallback):
            score=reference._grouped_state_te_from_factorized(values,packed[7][ids[j]],tc,tw,pv,fv,factor,max_unique=8192)
            if score is None:
                score=reference._grouped_state_te_from_coded(values,packed[7][ids[j]],tc,pv,fv,max_unique=8192)
            if score is None:
                raise RuntimeError('Original grouped permutation state limit exceeded')
            scores[j]=score
        result[begin:begin+len(ids)]=scores
    return result
