import csv
import os
import sys
import numpy as np
from multiprocessing import Pool, cpu_count
import time
start_time = time.time()

abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
os.chdir(dname)


def load_gene_chr_mapping(filename):
    with open(filename, "r") as f:
        return {gene: chromosome for chromosome, gene in (line.strip().split() for line in f) if not gene.startswith('chr')}

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

def create_pairs_for_peak(peak):
    gene_pairs = []
    peak_chr = peak.split('-')[0]
    for gene in genes_per_chromosome.get(peak_chr, []):
        if gene != peak and gene in gene_names_dict and peak in gene_names_dict:
            gene_pairs.append((peak, gene))
    return [[gene_names_dict[pair[0]] + 1, gene_names_dict[pair[1]] + 1] for pair in gene_pairs]

def create_gene_pairs_partial(tf):
    gene_pairs = []
    for gene in gene_names:
        if gene != tf and tf in gene_names_dict and gene in gene_names_dict:
            gene_pairs.append((tf, gene))
    return [[gene_names_dict[pair[0]] + 1, gene_names_dict[pair[1]] + 1] for pair in gene_pairs]

def create_gene_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list):
    gene_pairs = []
    
    
    with Pool(num_processes) as pool: 
        gene_pairs += [pair for sublist in pool.map(create_pairs_for_peak, Peaks) for pair in sublist]
                
    
    with Pool(num_processes) as pool:  
        gene_pairs += [pair for sublist in pool.map(create_gene_pairs_partial, tf_list) for pair in sublist]
                
    return gene_pairs

# Load gene-chromosome mapping
gene_chr_mapping = load_gene_chr_mapping("gene_chr.txt")

# Load and transpose data
expression_data = load_and_transpose_tsv_to_csv("cell_gene.tsv", "cell_gene_trsps.csv")

# Load gene names, Peaks, and TF list
gene_names = load_file_to_list("gene_names")
Peaks = load_file_to_list(f"TE_peak_list.txt")
tf_list = load_file_to_list("GO_symbol_human_regulation_of_transcription+sequence-specific_DNA_binding_list.txt")


gene_names_dict = {name: idx for idx, name in enumerate(gene_names)}

# Create a dictionary mapping chromosome to genes
genes_per_chromosome = {}
for gene, chromosome in gene_chr_mapping.items():
    if gene in gene_names_dict:
        genes_per_chromosome.setdefault(chromosome, []).append(gene)

num_processes = int(sys.argv[1])

# Create gene pairs and save to all_pairs.csv
gene_pairs_indices = create_gene_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list)
with open("all_pairs.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(gene_pairs_indices)
    
print("---Preprocess time : %s seconds ---" % (time.time() - start_time))
