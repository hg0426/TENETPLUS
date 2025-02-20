import csv
import numpy as np
import time
import sys

start_time = time.time()

def load_file_to_list(filename, delimiter="\n"):
    with open(filename, "r") as f:
        return [line.replace("\n", "") for line in f]

def load_and_transpose_tsv_to_csv(input_tsv, output_csv):
    with open(input_tsv, "r") as f_in:
        data = np.array(list(csv.reader(f_in, delimiter=" "))).astype("float").T
    with open(output_csv, "w") as f_out:
        writer = csv.writer(f_out)
        writer.writerows(data)
    return data

def create_gene_pairs(gene_names_dict, tf_list):
    gene_pairs = []
    
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend([(tf_idx + 1, gene_names_dict[gene] + 1) for gene in gene_names_dict if gene_names_dict[gene] != tf_idx])
    
    return gene_pairs

# Load and transpose data
expression_data = load_and_transpose_tsv_to_csv("cell_gene.tsv", "cell_gene_trsps.csv")

# Load gene names, and TF list
gene_names = load_file_to_list("gene_names")
tf_list = load_file_to_list("GO_symbol_"+sys.argv[1]+"_regulation_of_transcription+sequence-specific_DNA_binding_list.txt")

gene_names_dict = {name: idx for idx, name in enumerate(gene_names)}


# Create gene pairs and save to all_pairs.csv
gene_pairs_indices = create_gene_pairs(gene_names_dict, tf_list)
with open("all_pairs.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(gene_pairs_indices)
    
print("---Preprocess time : %s seconds ---" % (time.time() - start_time))
