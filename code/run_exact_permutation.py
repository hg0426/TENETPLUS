"""Resumable exact permutation runner for prepared TENET inputs, <=32 workers.

Baseline calls the original same-seed permutation kernel, with shared prepared
contexts. Optimized reuses each source's exact shuffles across its targets.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import sys,time,json,argparse,hashlib,multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd
from . import runTE_for_py_python_batch as ref
from . import permutation_exact as fast
from . import perm_counts as engine
from .reference_bundle import sha256
from .feature_utils import read_names, load_manifest

def save(path,obj):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2));temp.replace(path)

def worker(task):
    index,source,start,end=task
    t=time.process_time()
    values,codes=VALUES[source-1],CODES[source-1]
    result=np.empty(end-start,dtype=np.int32)
    if MODE=='optimized':pack=fast.prepare_permutations(values,codes,source,NPERM,SEED)
    prep=time.process_time()-t
    for j,row in enumerate(range(start,end)):
        target_id=TARGETS[row]
        if MODE=='optimized':
            target,context=CONTEXTS[target_id-1]
            tc,pv,fv,tw,factor=context
            ctx={'target_codes':tc,'target_past_values':pv,'target_next_values':fv,'target_counts':tw,'target_factor_context':factor,'n_obs':len(tc)}
            original=ref._kernel_grouped_te_for_source_codes(values,codes,ctx)
            if np.float64(original).view(np.uint64)!=OBSERVED[row].view(np.uint64):
                raise ValueError(('Observed TE differs',source,int(target_id)))
            scores=fast.score_permutations(pack,values,target,context,ref)
            result[j]=np.count_nonzero(scores>=OBSERVED[row])
        else:
            target,context=CONTEXTS[target_id-1]
            tc,pv,fv,tw,factor=context
            ctx={'target_codes':tc,'target_past_values':pv,'target_next_values':fv,
                 'target_counts':tw,'target_factor_context':factor,'n_obs':len(tc)}
            original,p=ref._kernel_te_perm_pvalue_grouped(ctx,source,1,.5,NPERM,SEED+source)
            if np.float64(original).view(np.uint64)!=OBSERVED[row].view(np.uint64):
                raise ValueError(('Observed TE differs',source,int(target_id)))
            result[j]=round(p*(NPERM+1)-1)
    return index,start,end,result,time.process_time()-t,prep

def main():
    global VALUES,CODES,CONTEXTS,TARGETS,OBSERVED,MODE,NPERM,SEED
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['baseline','optimized'],default='optimized')
    parser.add_argument('--perms',type=int,required=True)
    parser.add_argument('--jobs',type=int,default=32,choices=range(1,33))
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--matrix', required=True, help='Prepared feature-by-cell Parquet, k=1 and kernel width=0.5')
    parser.add_argument('--observed', required=True, help='Integer Source/Target and observed TE Parquet')
    parser.add_argument('--gene-names', required=True)
    parser.add_argument('--output', required=True)
    args=parser.parse_args()
    if args.perms<1:parser.error('perms must be positive')
    if args.seed < 0:parser.error('seed must be nonnegative')
    if str(ref.KERNEL_PERM_TABLE_SAMPLER).lower() not in ('off','0','false','no'):
        parser.error('Exact shuffle runner requires TE_KERNEL_PERM_TABLE_SAMPLER=off')
    MODE,NPERM,SEED=args.mode,args.perms,args.seed
    cpus=sorted(os.sched_getaffinity(0))[-args.jobs:];os.sched_setaffinity(0,cpus)
    output=Path(args.output).resolve();output.mkdir(parents=True,exist_ok=True)
    matrix_path=Path(args.matrix).resolve();observed_path=Path(args.observed).resolve()
    names_path=Path(args.gene_names).resolve()
    os.environ['TENET_OUTPUT_DIR']=str(names_path.parent)
    names=read_names(names_path);manifest=load_manifest(names)
    start=time.perf_counter()
    save(output/'STATUS.json',{'state':'PREPARING','mode':MODE,'permutations':NPERM,'cpus':cpus})
    matrix=pd.read_parquet(matrix_path).to_numpy(np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != len(names) or not 3 <= matrix.shape[1] <= 8193 or not np.isfinite(matrix).all():
        raise ValueError('Exact permutation requires finite feature-by-cell data, 3..8193 cells, matching gene_names')
    ref.GLOBAL_CELL_GENE_DENSE=matrix
    # Each worker verifies its observed TE against the original kernel.
    frame=pd.read_parquet(observed_path,columns=['Source','Target','TE'])
    if frame.empty or frame.duplicated(['Source','Target']).any() or not np.isfinite(frame.TE).all():
        raise ValueError('Observed TE must be nonempty, finite, and have unique pairs')
    for col in ['Source','Target']:
        if frame[col].dtype.kind not in 'iu' or not frame[col].between(1,len(names)).all():
            raise ValueError('Observed endpoints must be valid 1-based integer feature IDs')
    frame=frame.sort_values(['Source','Target'],kind='stable').reset_index(drop=True)
    SOURCES=frame.Source.to_numpy(np.int64);TARGETS=frame.Target.to_numpy(np.int64);OBSERVED=frame.TE.to_numpy(np.float64)
    starts=np.r_[0,np.flatnonzero(np.diff(SOURCES))+1];ends=np.r_[starts[1:],len(frame)]
    tasks=[(i,int(SOURCES[a]),int(a),int(b)) for i,(a,b) in enumerate(zip(starts,ends))]
    VALUES=[];CODES=[]
    for i in range(len(matrix)):
        values,codes=ref._get_grouped_source_codes_kernel(i+1,1,.5)
        VALUES.append(values);CODES.append(codes)
    CONTEXTS=[engine.prepare_target(series,ref) for series in matrix]
    warm=fast.prepare_permutations(VALUES[0],CODES[0],1,2,SEED)
    fast.score_permutations(warm,VALUES[0],*CONTEXTS[0],ref)
    config={'mode':MODE,'permutations':NPERM,'seed':SEED,'pairs':len(frame),'sources':len(tasks),
            'pair_order_sha256':hashlib.sha256(np.column_stack((SOURCES,TARGETS)).tobytes()).hexdigest(),
            'code_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}}
    config['input_hashes']={str(p):sha256(p) for p in [matrix_path,observed_path,names_path,names_path.parent/'feature_manifest.json']}
    config['numpy_version']=np.__version__
    config['kernel_environment']={k:v for k,v in os.environ.items() if k.startswith('TE_KERNEL_')}
    if (output/'CONFIG.json').exists():
        if json.loads((output/'CONFIG.json').read_text()) != config:
            raise ValueError('Resume configuration differs')
    else:save(output/'CONFIG.json',config)
    if (output/'FINAL.json').exists():
        final=json.loads((output/'FINAL.json').read_text())
        if sha256(output/'TE_kernel_perm.parquet') != final['result_sha256']:
            raise ValueError('Completed permutation result checksum mismatch')
        save(output/'STATUS.json',dict(state='COMPLETE',**final))
        print(json.dumps(final));return
    countpath=output/'exceedances.npy';donepath=output/'done.npy'
    if countpath.exists():
        counts=np.load(countpath,mmap_mode='r+');done=np.load(donepath)
    else:
        counts=np.lib.format.open_memmap(countpath,mode='w+',dtype=np.int32,shape=(len(frame),));counts[:]=-1;counts.flush()
        done=np.zeros(len(tasks),bool);np.save(donepath,done)
    pending=[x for x in tasks if not done[x[0]]]
    pending.sort(key=lambda x:-(x[3]-x[2])*max(1,len(VALUES[x[1]-1])))
    prepwall=time.perf_counter()-start;compute_start=time.perf_counter();cpu_sum=0.;sourceprep=0.;nreturned=0
    completed_pairs=int(sum(b-a for i,s,a,b in tasks if done[i]));initial_pairs=completed_pairs
    def checkpoint():
        counts.flush();np.save(output/'done.next.npy',done);(output/'done.next.npy').replace(donepath)
        elapsed=time.perf_counter()-compute_start;newpairs=completed_pairs-initial_pairs
        save(output/'STATUS.json',{'state':'RUNNING','mode':MODE,'permutations':NPERM,'jobs':args.jobs,'cpus':cpus,'completed_pairs':completed_pairs,'total_pairs':len(frame),
             'completed_sources':int(done.sum()),'total_sources':len(tasks),'preparation_seconds':prepwall,'compute_seconds_this_session':elapsed,
             'worker_cpu_seconds_this_session':cpu_sum,'estimated_remaining_seconds':elapsed*(len(frame)-completed_pairs)/newpairs if newpairs else None})
    checkpoint()
    with mp.get_context('fork').Pool(args.jobs) as pool:
        for i,a,b,values,cpu,scpu in pool.imap_unordered(worker,pending,chunksize=1):
            counts[a:b]=values;done[i]=True;cpu_sum+=cpu;sourceprep+=scpu;completed_pairs+=b-a;nreturned+=1
            if nreturned%25==0:checkpoint()
    checkpoint();compute_seconds=time.perf_counter()-compute_start
    if not (done.all() and np.all(counts>=0) and np.all(counts<=NPERM)):
        raise ValueError('Incomplete or invalid permutation exceedance counts')
    frame['exceedances']=np.asarray(counts)
    frame['p_value']=(np.asarray(counts,dtype=np.float64)+1.)/(NPERM+1.)
    p=frame.p_value.to_numpy();order=np.argsort(p)
    ranks=np.empty(len(p),np.int64);ranks[order]=np.arange(1,len(p)+1)
    q=p*len(p)/ranks
    q_sorted=np.minimum.accumulate(q[order][::-1])[::-1]
    q_final=np.empty_like(q);q_final[order]=q_sorted
    frame['q_value']=np.clip(q_final,0,1)
    frame.to_parquet(output/'TE_kernel_perm.parquet',index=False)
    network_counts={}
    for alpha in [.01,.05]:
        edges=frame[frame.q_value<alpha]
        edges.to_parquet(output/f'network_q{alpha}.parquet',index=False)
        network_counts[str(alpha)]=len(edges)
    result=dict(result_sha256=sha256(output/'TE_kernel_perm.parquet'),status='PASS',mode=MODE,permutations=NPERM,pairs=len(frame),jobs=args.jobs,cpus=cpus,
                preparation_seconds=prepwall,compute_seconds=compute_seconds,worker_cpu_seconds=cpu_sum,source_prepare_cpu_seconds=sourceprep,
                total_seconds=time.perf_counter()-start,networks=network_counts,resumed_initial_pairs=initial_pairs,
                note='Exact source shuffles for all supplied observed pairs; no early stopping or table sampling.')
    save(output/'FINAL.json',result);save(output/'STATUS.json',dict(state='COMPLETE',**result));print(json.dumps(result),flush=True)

if __name__=='__main__':main()
