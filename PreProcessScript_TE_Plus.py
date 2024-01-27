import csv
import numpy as np
import time
import sys

start_time = time.time()

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

def create_TENET_Plus_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list):
    gene_pairs = []
    
    for peak in Peaks:
        peak_chr = peak.split('-')[0]
        peak_idx = gene_names_dict.get(peak, None)
        if peak_idx is not None:
            gene_pairs.extend([(peak_idx + 1, gene_names_dict[gene] + 1) for gene in genes_per_chromosome.get(peak_chr, []) if gene in gene_names_dict])
    
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend([(tf_idx + 1, gene_names_dict[gene] + 1) for gene in gene_names_dict if gene_names_dict[gene] != tf_idx])
    
    return gene_pairs
    
def create_tf_peak_pairs(gene_names_dict,only_peak_names_dict, tf_list):
    peak_pairs = []
    
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        peak_pairs.extend([(tf_idx + 1, only_peak_names_dict[peak] + 1) for peak in only_peak_names_dict if only_peak_names_dict[peak] != tf_idx])
    
    return peak_pairs
    
def create_tf_gene_pairs(only_gene_names_dict, tf_list):
    gene_pairs = []
    
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend([(tf_idx + 1, only_gene_names_dict[gene] + 1) for gene in only_gene_names_dict if only_gene_names_dict[gene] != tf_idx])
    
    return gene_pairs

def create_cis_peaksource_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list):
    gene_pairs = []
    
    for peak in Peaks:
        peak_chr = peak.split('-')[0]
        peak_idx = gene_names_dict.get(peak, None)
        if peak_idx is not None:
            gene_pairs.extend([(peak_idx + 1, gene_names_dict[gene] + 1) for gene in genes_per_chromosome.get(peak_chr, []) if gene in gene_names_dict])
            
    return gene_pairs

# Load gene-chromosome mapping
gene_chr_mapping = load_gene_chr_mapping("gene_chr.txt")

# Load and transpose data
expression_data = load_and_transpose_tsv_to_csv("cell_gene.tsv", "cell_gene_trsps.csv")

# Load gene names, Peaks, and TF list
gene_names = load_file_to_list("gene_names")
Peaks = load_file_to_list(f"TE_peak_list.txt")
tf_list = load_file_to_list("GO_symbol_"+sys.argv[1]+"_regulation_of_transcription+sequence-specific_DNA_binding_list.txt")

gene_names_dict = {name: idx for idx, name in enumerate(gene_names)}
only_peak_names_dict = {gene: idx for gene, idx in gene_names_dict.items() if gene in Peaks}
only_gene_names_dict = {gene: idx for gene, idx in gene_names_dict.items() if gene not in Peaks}


# Create a dictionary mapping chromosome to genes
genes_per_chromosome = {}
for gene, chromosome in gene_chr_mapping.items():
    if gene in gene_names_dict:
        genes_per_chromosome.setdefault(chromosome, []).append(gene)


# Create gene pairs and save to all_pairs.csv

if sys.argv[2] == "1":
    gene_pairs_indices = create_TENET_Plus_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list)
elif sys.argv[2] == "2":
    gene_pairs_indices = create_tf_gene_pairs(only_gene_names_dict, tf_list)
elif sys.argv[2] == "3":
    gene_pairs_indices = create_tf_peak_pairs(gene_names_dict,only_peak_names_dict, tf_list)
elif sys.argv[2] == "4":
    gene_pairs_indices = create_cis_peaksource_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list)




with open("all_pairs.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(gene_pairs_indices)
    
print("---Preprocess time : %s seconds ---" % (time.time() - start_time))
