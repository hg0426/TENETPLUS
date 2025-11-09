import csv
import sys

import numpy

from code.path_utils import coerce_input_path, coerce_output_path


def load_matrix(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=" ")
        return numpy.array(list(reader)).astype(float)


def load_list(path):
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def main() -> None:
    matrix_path = coerce_input_path("cell_gene.tsv")
    expression_data = load_matrix(matrix_path).T
    transposed_path = coerce_output_path("cell_gene_trsps.csv")
    with open(transposed_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(expression_data)

    gene_names = load_list(coerce_input_path("gene_names"))
    tfs = load_list(
        coerce_input_path(
            f"GO_symbol_{sys.argv[1]}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
        )
    )
    tf_index = {tf: idx + 1 for idx, tf in enumerate(gene_names) if tf in tfs}

    indices = []
    for tf_idx in tf_index.values():
        for gene_idx in range(1, len(gene_names) + 1):
            if tf_idx != gene_idx:
                indices.append([tf_idx, gene_idx])

    all_pairs_path = coerce_output_path("all_pairs.csv")
    with open(all_pairs_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(indices)


if __name__ == "__main__":
    main()
