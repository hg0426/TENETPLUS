# main_script.py

import os
import pandas as pd
import argparse
import numpy as np
import multiprocessing
import time
from scipy import sparse
from tqdm import tqdm
from Calc_TE_python_localTE import TransferEntropyCalculatorKernel, EmpiricalMeasurementDistribution
import datetime
import logging
import duckdb

# Configure logging to track the progress and debug if needed
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create file handler which logs even debug messages
fh = logging.FileHandler('TE_analysis.log')
fh.setLevel(logging.INFO)

# Create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# Create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s %(levelname)s:%(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

def calculate_TE_batch(batch_params):
    """
    Worker function to calculate local Transfer Entropy (TE) for a batch of pairs.

    Args:
        batch_params (tuple): Contains list_pairs, cell_gene_all, historyLength, kernel_width.

    Returns:
        list: A list of tuples in the form (source, target, [t1, t2, ..., tN]).
    """
    list_pairs, cell_gene_all, historyLength, kernel_width = batch_params

    results = []
    for pair in list_pairs:
        source_idx, target_idx = pair
        te_calculator = TransferEntropyCalculatorKernel()
        te_calculator.set_debug(False)
        te_calculator.set_property("KERNEL_WIDTH", str(kernel_width))
        te_calculator.set_property("NORMALISE", "True")
        te_calculator.initialise(k=historyLength)

        # Adjust indices for 0-based indexing
        try:
            expression_data_source = cell_gene_all[source_idx - 1].toarray().ravel()
            expression_data_target = cell_gene_all[target_idx - 1].toarray().ravel()
        except IndexError as e:
            logging.error(f"IndexError for pair ({source_idx}, {target_idx}): {e}")
            continue  # Skip this pair if indices are out of bounds

        te_calculator.add_observations(source=expression_data_source, destination=expression_data_target)
        te_calculator.finalise_add_observations()
        local_te = te_calculator.compute_local_TE()

        results.append((source_idx, target_idx, local_te.tolist()))

    return results

def merge_parquet_files(progress_dir, merged_filename_prefix='merged_TE_progress'):
    """
    Merges all batch Parquet files in the progress directory into a single merged Parquet file using DuckDB.
    Handles incomplete/corrupted files by excluding them from the merge and deleting them.

    Args:
        progress_dir (str): Directory containing Parquet files to merge.
        merged_filename_prefix (str): Prefix for the merged Parquet file.
    """
    start_time = time.time()
    con = duckdb.connect()

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    merged_filename = f"{merged_filename_prefix}_{timestamp}.parquet"
    merged_filepath = os.path.join(progress_dir, merged_filename)

    valid_files = []
    invalid_files = []  # Keep track of invalid files for deletion
    for fname in os.listdir(progress_dir):
        if fname.endswith('.parquet') and fname.startswith('batch_'):
            filepath = os.path.join(progress_dir, fname)
            try:
                # Attempt to read the Parquet file to check for validity
                con.execute(f"SELECT * FROM read_parquet('{filepath}') LIMIT 1")
                valid_files.append(filepath)
                logging.info(f"File {fname} is valid.")
            except Exception as e:
                logging.warning(f"Skipping invalid/incomplete file {fname}: {e}")
                invalid_files.append(filepath)  # Add invalid file to the list

    if not valid_files:
        logging.warning("No valid batch Parquet files found to merge.")
        con.close()
        return

    try:
        # Include only valid batch Parquet files in the merge
        file_list_str = ', '.join([f"'{file}'" for file in valid_files])
        query = f"""
        COPY (
            SELECT * FROM read_parquet([{file_list_str}])
        ) TO '{merged_filepath}' (FORMAT PARQUET);
        """
        con.execute(query)
        logging.info(f"Merged {len(valid_files)} batch Parquet files into {merged_filename} using DuckDB.")
    except Exception as e:
        logging.error(f"Error merging Parquet files: {e}")
    finally:
        con.close()
        end_time = time.time()
        logging.info(f"Parquet merging took {end_time - start_time:.2f} seconds.")

    # Delete invalid files
    for filepath in invalid_files:
        try:
            os.remove(filepath)
            logging.info(f"Deleted invalid file: {filepath}")
        except Exception as e:
            logging.error(f"Error deleting invalid file {filepath}: {e}")

    # Delete valid original batch files after merging
    for filepath in valid_files:
        try:
            os.remove(filepath)
            logging.info(f"Deleted batch file after merging: {filepath}")
        except Exception as e:
            logging.error(f"Error deleting batch file {filepath}: {e}")

def run_parallel_batches(list_pairs, cell_gene_all, historyLength, kernel_width, num_cpus, batch_size, progress_dir, buffer_size=10, merge_threshold=50, enable_intermediate_save=True):
    """
    Processes the list of pairs in batches using multiprocessing Pool.
    Writes results to Parquet files incrementally in the specified progress directory if enable_intermediate_save is True.
    Merges only the new batch Parquet files when the number of batch files exceeds the merge_threshold.

    Args:
        list_pairs (np.ndarray): Array of pairs to process.
        cell_gene_all (scipy.sparse.csr_matrix): Sparse matrix of gene expressions.
        historyLength (int): History length (k) for the analysis.
        kernel_width (float): Kernel width parameter.
        num_cpus (int): Number of parallel CPUs to use.
        batch_size (int): Number of pairs per batch.
        progress_dir (str): Directory to store intermediate Parquet files.
        buffer_size (int): Number of batches to accumulate before writing.
        merge_threshold (int): Number of batch Parquet files to trigger a merge.
        enable_intermediate_save (bool): Flag to enable/disable intermediate saving.
    """
    if enable_intermediate_save:
        # Ensure the progress directory exists
        os.makedirs(progress_dir, exist_ok=True)

    # Split list_pairs into batches
    batches = [list_pairs[i:i + batch_size] for i in range(0, len(list_pairs), batch_size)]
    total_batches = len(batches)

    # Prepare parameters for each batch
    batch_params = [(batch, cell_gene_all, historyLength, kernel_width) for batch in batches]

    # Determine number of time points based on the first pair
    if len(list_pairs) > 0:
        first_pair = list_pairs[0]
        _, target_idx = first_pair
        try:
            expression_data_target = cell_gene_all[target_idx - 1].toarray().ravel()
            num_time_points = len(expression_data_target) - historyLength
            columns = ['Source', 'Target'] + [f't{t+1}' for t in range(num_time_points)]
        except IndexError as e:
            logging.error(f"IndexError when determining number of time points: {e}")
            columns = ['Source', 'Target']  # Fallback in case of error
    else:
        columns = ['Source', 'Target']  # Handle empty case

    # Initialize buffer for accumulating batches
    buffered_batches = []
    all_results = []

    # Use multiprocessing Pool for parallel processing
    with multiprocessing.Pool(processes=num_cpus) as pool:
        with tqdm(total=total_batches, desc="Processing Batches") as pbar:
            for batch_result in pool.imap_unordered(calculate_TE_batch, batch_params):
                # Prepare DataFrame for the batch
                batch_records = []
                for source, target, te_values in batch_result:
                    record = {'Source': source, 'Target': target}
                    for t, te in enumerate(te_values):
                        record[f't{t+1}'] = te
                    batch_records.append(record)

                batch_df = pd.DataFrame(batch_records, columns=columns)
                
                if enable_intermediate_save:
                    buffered_batches.append(batch_df)
                    # If buffer is full, write to a batch Parquet file
                    if len(buffered_batches) >= buffer_size:
                        combined_df = pd.concat(buffered_batches, ignore_index=True)
                        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                        batch_filename = os.path.join(progress_dir, f'batch_{timestamp}.parquet')
                        try:
                            combined_df.to_parquet(batch_filename, index=False)
                            logging.info(f"Written {batch_filename}")
                        except Exception as e:
                            logging.error(f"Error writing {batch_filename}: {e}")
                            raise e
                        buffered_batches = []  # Reset buffer

                        # Check if merging is needed
                        current_batch_file_count = len([
                            fname for fname in os.listdir(progress_dir)
                            if fname.endswith('.parquet') and fname.startswith('batch_')
                        ])
                        if current_batch_file_count >= merge_threshold:
                            logging.info(f"Batch file count {current_batch_file_count} reached merge threshold {merge_threshold}. Merging batch files.")
                            merge_parquet_files(progress_dir)
                else:
                    all_results.append(batch_df)

                pbar.update(1)
    if enable_intermediate_save:
        # Write any remaining batches in the buffer
        if buffered_batches:
            combined_df = pd.concat(buffered_batches, ignore_index=True)
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            batch_filename = os.path.join(progress_dir, f'batch_{timestamp}.parquet')
            try:
                combined_df.to_parquet(batch_filename, index=False)
                logging.info(f"Written {batch_filename}")
            except Exception as e:
                logging.error(f"Error writing {batch_filename}: {e}")
                raise e

        # Final check for merging remaining batch files
        current_batch_file_count = len([
            fname for fname in os.listdir(progress_dir)
            if fname.endswith('.parquet') and fname.startswith('batch_')
        ])
        if current_batch_file_count > 0:
            logging.info(f"Final batch file count {current_batch_file_count}. Merging remaining batch files.")
            merge_parquet_files(progress_dir)
    else:
        # Combine all results if intermediate saving is disabled
        if all_results:
            final_results_df = pd.concat(all_results, ignore_index=True)
            final_results_df.to_parquet('TE_result_all.parquet', index=False)
            logging.info("Final results saved to TE_result_all.parquet.")

def load_progress(input_csv, progress_dir):
    """
    Loads progress from the progress directory if it exists.
    Returns the remaining pairs to process.
    Handles incomplete/corrupted files by skipping them.

    Args:
        input_csv (str): Path to the input CSV file containing all pairs.
        progress_dir (str): Directory where intermediate Parquet files are stored.

    Returns:
        np.ndarray: Array of remaining pairs to process.
    """
    load_start_time = time.time()

    # Check for existing merged Parquet files and load them
    merged_files = [
        fname for fname in os.listdir(progress_dir)
        if fname.endswith('.parquet') and fname.startswith('merged_')
    ]

    if len(merged_files) > 1:
        logging.info(f"Multiple merged Parquet files detected at startup. Merging them.")
        merge_parquet_files(progress_dir)
        # After merging, there should be only one merged file
        merged_files = [
            fname for fname in os.listdir(progress_dir)
            if fname.endswith('.parquet') and fname.startswith('merged_')
        ]
    elif len(merged_files) == 1:
        logging.info(f"Single merged Parquet file detected: {merged_files[0]}.")
    else:
        logging.info("No merged Parquet files detected in progress directory.")

    # Load all pairs from input CSV
    try:
        list_pairs = pd.read_csv(input_csv, delimiter=',', header=None).to_numpy().astype(int)
        logging.info(f"Loaded {len(list_pairs)} pairs from {input_csv}.")
    except Exception as e:
        logging.error(f"Error loading input CSV {input_csv}: {e}")
        return np.array([])

    pairs_df = pd.DataFrame(list_pairs, columns=['Source', 'Target'])

    # Load processed pairs from merged Parquet files, handling potential errors
    processed_pairs = []
    for fname in merged_files:
        filepath = os.path.join(progress_dir, fname)
        try:
            df = pd.read_parquet(filepath)
            processed_pairs.append(df[['Source', 'Target']])
            logging.info(f"Loaded processed pairs from {fname}.")
        except Exception as e:
            logging.error(f"Error reading file {fname}: {e}. Skipping this file.")

    if processed_pairs:
        processed_pairs_df = pd.concat(processed_pairs, ignore_index=True).drop_duplicates()
        logging.info(f"Total processed pairs: {len(processed_pairs_df)}.")

        # Identify remaining pairs by excluding processed ones
        remaining_pairs_df = pairs_df.merge(
            processed_pairs_df,
            on=['Source', 'Target'],
            how='left',
            indicator=True
        )
        remaining_pairs = remaining_pairs_df[remaining_pairs_df['_merge'] == 'left_only'][['Source', 'Target']].to_numpy()
        logging.info(f"Remaining pairs to process: {len(remaining_pairs)}.")

    else:
        logging.info("No processed pairs found in progress directory. Starting fresh.")
        remaining_pairs = list_pairs

    load_end_time = time.time()
    logging.info(f"Loading progress took {load_end_time - load_start_time:.2f} seconds.")
    return remaining_pairs

def main(args):
    """
    Main function to execute the TE analysis workflow.

    Args:
        args (Namespace): Parsed command-line arguments.
    """
    start_time = time.time()
    print('Starting Calculate_TE')
    logging.info('Starting Calculate_TE')

    # Load expression data
    print('Loading expression data...')
    logging.info('Loading expression data...')
    load_data_start_time = time.time()
    try:
        cell_gene_all = sparse.csr_matrix(pd.read_parquet('cell_gene_trsps.parquet').to_numpy())
        print('Expression data loaded successfully.')
        logging.info('Expression data loaded successfully.')
    except Exception as e:
        print(f"Error loading expression data: {e}")
        logging.error(f"Error loading expression data: {e}")
        return
    load_data_end_time = time.time()
    logging.info(f"Loading expression data took {load_data_end_time - load_data_start_time:.2f} seconds.")

    # Define progress directory
    progress_dir = 'TE_progress_parquet' if args.enable_intermediate_save else None
    
    # Create progress dir if intermediate saving is enabled
    if args.enable_intermediate_save:
        os.makedirs(progress_dir, exist_ok=True)
        
    # Load or initialize list of pairs
    if args.enable_intermediate_save:
        list_pairs = load_progress(args.input_csv, progress_dir)
    else:
        try:
            list_pairs = pd.read_csv(args.input_csv, delimiter=',', header=None).to_numpy().astype(int)
            logging.info(f"Loaded {len(list_pairs)} pairs from {args.input_csv}.")
        except Exception as e:
            logging.error(f"Error loading input CSV {args.input_csv}: {e}")
            return np.array([])

    total_pairs = len(list_pairs)
    print(f"Total pairs to process: {total_pairs}")
    logging.info(f"Total pairs to process: {total_pairs}")

    if total_pairs == 0 and args.enable_intermediate_save:
        print("All pairs have already been processed.")
        logging.info("All pairs have already been processed.")

        # If intermediate files exist, consolidate them into the final result
        all_batches = []
        for fname in sorted(os.listdir(progress_dir)):
            if fname.endswith('.parquet') and fname.startswith('merged_'):
                try:
                    batch_df = pd.read_parquet(os.path.join(progress_dir, fname))
                    all_batches.append(batch_df)
                    logging.info(f"Loaded batch {fname} for consolidation.")
                except Exception as e:
                    print(f"Error reading file {fname}: {e}")
                    logging.error(f"Error reading file {fname}: {e}")

        if all_batches:
            try:
                final_results_df = pd.concat(all_batches, ignore_index=True)
                final_results_df.to_parquet('TE_result_all.parquet', index=False)
                print("Final results saved to TE_result_all.parquet.")
                logging.info("Final results saved to TE_result_all.parquet.")
            except Exception as e:
                print(f"Error saving final results: {e}")
                logging.error(f"Error saving final results: {e}")
                return

            # Optionally, delete intermediate merged files
            for fname in os.listdir(progress_dir):
                if fname.endswith('.parquet') and fname.startswith('merged_'):
                    try:
                        os.remove(os.path.join(progress_dir, fname))
                        logging.info(f"Deleted merged file {fname}.")
                    except Exception as e:
                        print(f"Error deleting file {fname}: {e}")
                        logging.error(f"Error deleting file {fname}: {e}")
            print("Intermediate merged progress files have been deleted.")
            logging.info("Intermediate merged progress files have been deleted.")
        else:
            print("No progress files found to combine.")
            logging.info("No progress files found to combine.")
        return

    # Define batch size
    batch_size = 500  # Adjust based on memory and performance considerations

    # Define buffer size (number of batches to accumulate before writing)
    buffer_size = 5000  # For example: accumulate 100 batches before writing to Parquet

    # Define merge threshold (number of batch Parquet files to trigger a merge)
    merge_threshold = 30  # Adjust based on desired maximum number of files

    # Run parallel batch processing
    run_parallel_batches(
        list_pairs=list_pairs,
        cell_gene_all=cell_gene_all,
        historyLength=int(args.history_length),
        kernel_width=0.5,  # Adjust as needed
        num_cpus=args.num_jobs,
        batch_size=batch_size,
        progress_dir=progress_dir,
        buffer_size=buffer_size,
        merge_threshold=merge_threshold,
        enable_intermediate_save=args.enable_intermediate_save
    )

    print("--- Calculate TE execution time: %s seconds ---" % (time.time() - start_time))
    logging.info(f"Calculate TE execution time: {time.time() - start_time:.2f} seconds")

    if args.enable_intermediate_save:
        # Combine all merged Parquet files into the final result
        all_batches = []
        combine_start_time = time.time()
        for fname in sorted(os.listdir(progress_dir)):
            if fname.endswith('.parquet') and fname.startswith('merged_'):
                try:
                    batch_df = pd.read_parquet(os.path.join(progress_dir, fname))
                    all_batches.append(batch_df)
                    logging.info(f"Loaded merged batch {fname} for final consolidation.")
                except Exception as e:
                    print(f"Error combining file {fname}: {e}")
                    logging.error(f"Error combining file {fname}: {e}")

        if all_batches:
            try:
                final_results_df = pd.concat(all_batches, ignore_index=True)
                final_results_df.to_parquet('TE_result_all.parquet', index=False)
                print("Final results saved to TE_result_all.parquet.")
                logging.info("Final results saved to TE_result_all.parquet.")
            except Exception as e:
                print(f"Error saving final results: {e}")
                logging.error(f"Error saving final results: {e}")
                return

            # Optionally, delete intermediate merged files
            for fname in os.listdir(progress_dir):
                if fname.endswith('.parquet') and fname.startswith('merged_'):
                    try:
                        os.remove(os.path.join(progress_dir, fname))
                        logging.info(f"Deleted merged file {fname}.")
                    except Exception as e:
                        print(f"Error deleting file {fname}: {e}")
                        logging.error(f"Error deleting file {fname}: {e}")
            print("Intermediate merged progress files have been deleted.")
            logging.info("Intermediate merged progress files have been deleted.")
        else:
            print("No progress files found to combine.")
            logging.info("No progress files found to combine.")

        combine_end_time = time.time()
        logging.info(f"Final consolidation took {combine_end_time - combine_start_time:.2f} seconds.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TE analysis with parallel processing and batch management.")
    parser.add_argument('input_csv', type=str, help="Input CSV file containing all pairs.")
    parser.add_argument('num_jobs', type=int, help="Number of parallel jobs.")
    parser.add_argument('history_length', type=str, help="History length (k) for the analysis.")
    parser.add_argument('--enable_intermediate_save', action='store_true', help="Enable intermediate saving of results.")

    args = parser.parse_args()
    main(args)
