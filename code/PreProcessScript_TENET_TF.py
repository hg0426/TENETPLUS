import csv
import sys
import time

from code.path_utils import coerce_input_path, coerce_output_path

start_time = time.time()


def load_file_to_list(filename, delimiter="\n"):
    path = coerce_input_path(filename)
    with open(path, "r", encoding="utf-8") as f:
        return [line.replace("\n", "") for line in f]


def create_gene_pairs(gene_names_dict, tf_list):
    gene_pairs = []

    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend(
            [
                (tf_idx + 1, gene_names_dict[gene] + 1)
                for gene in gene_names_dict
                if gene_names_dict[gene] != tf_idx
            ]
        )

    return gene_pairs


# Load gene names, and TF list
gene_names = load_file_to_list("gene_names")
tf_list = load_file_to_list(
    f"GO_symbol_{sys.argv[1]}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
)

gene_names_dict = {name: idx for idx, name in enumerate(gene_names)}

# Create gene pairs and save to all_pairs.csv
gene_pairs_indices = create_gene_pairs(gene_names_dict, tf_list)
output_path = coerce_output_path("all_pairs.csv")
with open(output_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(gene_pairs_indices)

print("---Preprocess time : %s seconds ---" % (time.time() - start_time))
