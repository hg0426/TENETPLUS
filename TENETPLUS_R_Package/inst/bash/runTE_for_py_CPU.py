from jpype import startJVM, getDefaultJVMPath, JPackage, JArray, JDouble
import os
import pandas as pd
import argparse
import numpy as np
import multiprocessing
import time
from scipy import sparse
from multiprocessing import Process, Array
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
    
def calculate_TE(params, progress, idx_offset):
    jarLocation = "infodynamics.jar"
    startJVM(getDefaultJVMPath(), "-ea", f"-Djava.class.path={jarLocation}", "-Xmx16G")
    idx, list_pairs, cell_gene_all, historyLength = params
    teCalcClass = JPackage("infodynamics.measures.continuous.kernel").TransferEntropyCalculatorKernel
    teCalc = teCalcClass()
    teCalc.setProperty("NORMALISE", "true")
    teCalc.initialise(historyLength, 0.5)
    
    results = []
    for pair in list_pairs:
        expression_data = [list(cell_gene_all[pair[0] - 1].toarray().ravel()),
                               list(cell_gene_all[pair[1] - 1].toarray().ravel())]
        #expression_data = [list(cell_gene_all[pair[0] - 1]),
        #                       list(cell_gene_all[pair[1] - 1])]
        teCalc.setObservations(JArray(JDouble, 1)(expression_data[0]), JArray(JDouble, 1)(expression_data[1]))
        result = teCalc.computeAverageLocalOfObservations()
        results.append(result)
        with open('TE_progress.csv', 'a') as file:
            file.write(f"{idx},{','.join(map(str, pair))},{result}\n")
        progress[idx_offset] += 1  # Update progress


def run_parallel_process(params_list, num_cpus):
    manager = multiprocessing.Manager()
    progress = manager.Array('i', [0] * num_cpus)  # Shared array to track progress
    processes = []

    # Create and start a process for each CPU
    for i in range(num_cpus):
        ctx = multiprocessing.get_context("spawn")
        process = ctx.Process(target=calculate_TE, args=(params_list[i], progress, i))
        processes.append(process)
        process.start()

    # Monitor progress
    with tqdm(total=sum(len(params[1]) for params in params_list)) as pbar:
        current_progress = 0
        while any(p.is_alive() for p in processes):
            total_progress = sum(progress)
            pbar.update(total_progress - current_progress)
            current_progress = total_progress
            time.sleep(0.5)  # Update every half second

    # Ensure all processes finish
    for p in processes:
        p.join()


def load_progress():
    if os.path.exists('TE_progress.csv') and os.stat('TE_progress.csv').st_size > 0:
        list_pairs = pd.read_csv(args.input_csv, delimiter=',', header=None).to_numpy().astype(int)
        print('TE_progress.csv was found and continues with existing results. If you want to new proceed, please delete TE_progress.csv.')
        progress_df = pd.read_csv('TE_progress.csv', header=None,dtype={3:str})
        pairs_df = pd.DataFrame(list_pairs, columns=['Source', 'Target'])
        progress_df.columns = ['Index', 'Source', 'Target', 'Result']
        merged_df = pd.merge(pairs_df, progress_df, on=['Source', 'Target'], how='left', indicator=True)
        filtered_df = merged_df[merged_df['_merge'] == 'left_only']
        remaining_pairs = filtered_df[['Source', 'Target']].to_numpy()
        return remaining_pairs
    else:
        with open('TE_progress.csv', 'w') as file:
            file.write('')
        list_pairs = pd.read_csv(args.input_csv, delimiter=',', header=None).to_numpy().astype(int)
        return list_pairs

        

def main(args):
    start_time = time.time()
    print('Starting Calculate_TE')
    start_time2 = time.time()
    cell_gene_all = sparse.csr_matrix(pd.read_csv('cell_gene_trsps.csv', delimiter=',', header=None).to_numpy())
    
    #print('pairs')
    
    list_pairs = load_progress()
    
    #print('repro')
    # Prepare parameters and split work among CPUs
    num_cpus = args.num_jobs
    pairs_per_cpu = np.array_split(list_pairs, num_cpus)
    params_list = [(i, pairs, cell_gene_all, int(args.extra_param)) for i, pairs in enumerate(pairs_per_cpu)]
    
    print("---prepare : %s seconds ---" % (time.time() - start_time2))
    
    # Run parallel processing
    run_parallel_process(params_list, num_cpus)

    print("---Calculate TE time : %s seconds ---" % (time.time() - start_time))
    
    if os.path.exists('TE_progress.csv'):
        final_results_df = pd.read_csv('TE_progress.csv', header=None,dtype={3:str})
        final_results_df.to_csv('TE_result_all.csv', index=False, header=False,columns=[1,2,3])
        os.rename('TE_progress.csv','TE_progress_remain.csv')

if __name__ == "__main__":
    #multiprocessing.set_start_method('spawn')
    parser = argparse.ArgumentParser(description="Run parallel processing on TE analysis.")
    parser.add_argument('input_csv', type=str, help="The input CSV file containing all pairs.")
    parser.add_argument('num_jobs', type=int, help="The number of parallel jobs.")
    parser.add_argument('extra_param', type=str, help="Additional parameter for the analysis.")
    
    args = parser.parse_args()
    main(args)