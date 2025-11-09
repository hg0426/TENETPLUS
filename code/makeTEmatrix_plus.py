import sys
import pandas as pd
import numpy as np

from code.path_utils import locate_file, resolve_output

def load_gene_names(file_path):
    """Reads gene names from a file and returns them as a list."""
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            gene_names = [line.strip() for line in file]
        return gene_names
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")

def initialize_te_matrix(gene_count):
    """Creates a TE matrix initialized with zeros based on the number of genes."""
    return pd.DataFrame(
        np.zeros((gene_count, gene_count)),
        index=gene_names,
        columns=gene_names
    )

def populate_te_matrix(te_matrix, csv_path):
    """Fills the TE matrix with data from a CSV file."""
    try:
        with open(csv_path, "r", encoding="utf-8") as file:
            for line in file:
                parts = line.strip().split(',')
                if len(parts) < 3:
                    continue  # Skip if data is insufficient
                row_idx = int(parts[0]) - 1
                col_idx = int(parts[1]) - 1
                te_matrix.iat[row_idx, col_idx] = float(parts[2])
                if len(parts) > 3:
                    te_matrix.iat[col_idx, row_idx] = float(parts[3])
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {csv_path}")
    except ValueError as e:
        raise ValueError(f"Data format error: {e}")

def categorize_genes(gene_names, species, split_type):
    """
    Returns a list of genes categorized based on the split_type.
    """
    if split_type == 'rowPeak':
        peak_list_file = locate_file("TE_peak_list.txt")
        try:
            with open(peak_list_file, "r", encoding="utf-8") as file:
                peaks = set(line.strip() for line in file)
            selected_genes = [gene for gene in gene_names if gene in peaks]
            return selected_genes
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {peak_list_file}")
    elif split_type == 'rowTF':
        tf_list_file = locate_file(
            f"GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
        )
        try:
            with open(tf_list_file, "r", encoding="utf-8") as file:
                tfs = set(line.strip() for line in file)
            selected_genes = [gene for gene in gene_names if gene in tfs]
            return selected_genes
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {tf_list_file}")
    else:
        return gene_names

def split_columns(gene_names, te_matrix, split_type):
    """
    Categorizes columns and returns the filtered matrix.
    """
    if split_type == 'colGN':
        # Genomic Names (e.g., "gene1", "gene2", ...)
        col_genes = [gene for gene in gene_names if not gene.startswith('chr')]
        return te_matrix[col_genes]
    elif split_type == 'colPK':
        # Peak Columns (e.g., "chr1-12345-67890", ...)
        col_peaks = [gene for gene in gene_names if gene.startswith('chr')]
        return te_matrix[col_peaks]
    else:
        return te_matrix

def save_matrix(df, output_filename):
    """Saves the DataFrame to a file."""
    try:
        output_path = resolve_output(output_filename)
        df.to_csv(output_path, sep='\t', index=True, header=True)
    except IOError as e:
        raise IOError(f"Error saving file: {e}")

input_csv = locate_file(sys.argv[1])
species = sys.argv[2]
split_type = sys.argv[3] 

# Load gene names
gene_names_file = locate_file("gene_names")
gene_names = load_gene_names(gene_names_file)

# Initialize TE matrix
te_matrix = initialize_te_matrix(len(gene_names))

# Populate TE matrix with data
populate_te_matrix(te_matrix, input_csv)

te_matrix.index.name = 'TE'

# Split the matrix based on conditions
if split_type in ['1', '5']:
    print('Splitting Matrix for TENET_Plus(rowPeak_colGN)')
    # Split by rowPeak
    row_peak_genes = categorize_genes(gene_names, species, 'rowPeak')
    te_matrix_row_peak = te_matrix.loc[row_peak_genes]

    # Split by colGN
    te_matrix_row_peak_col_gn = split_columns(gene_names, te_matrix_row_peak, 'colGN')
    save_matrix(te_matrix_row_peak_col_gn, "TE_result_matrix_rowPeak_colGN.txt")

if split_type != '5':
    print('Splitting Matrix for rowTF')
    # Split by rowTF
    row_tf_genes = categorize_genes(gene_names, species, 'rowTF')
    te_matrix_row_tf = te_matrix.loc[row_tf_genes]

    if split_type not in ['3']:
        print('Splitting Matrix by colGN')
        te_matrix_row_tf_col_gn = split_columns(gene_names, te_matrix_row_tf, 'colGN')
        save_matrix(te_matrix_row_tf_col_gn, "TE_result_matrix_rowTF_colGN.txt")

    if split_type in ['1', '3', '4']:
        print('Splitting Matrix by colPK')
        te_matrix_row_tf_col_pk = split_columns(gene_names, te_matrix_row_tf, 'colPK')
        save_matrix(te_matrix_row_tf_col_pk, "TE_result_matrix_rowTF_colPK.txt")

print("Matrix splitting completed.")
