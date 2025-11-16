import csv
import sys
import time

from code.path_utils import coerce_input_path, coerce_output_path


def load_file_to_list(filename, delimiter: str = "\n"):
    path = coerce_input_path(filename)
    with open(path, "r", encoding="utf-8") as f:
        return [line.replace(delimiter, "") for line in f]


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


def run_preprocess(species: str) -> None:
    """Build all_pairs.csv for the original TENET_TF mode."""
    start_time = time.time()

    gene_names = load_file_to_list("gene_names")
    tf_list = load_file_to_list(
        f"GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
    )

    gene_names_dict = {name: idx for idx, name in enumerate(gene_names)}

    gene_pairs_indices = create_gene_pairs(gene_names_dict, tf_list)
    output_path = coerce_output_path("all_pairs.csv")
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(gene_pairs_indices)

    print(
        "---Preprocess time : %s seconds ---"
        % (time.time() - start_time)
    )


def main(argv: list[str] | None = None) -> None:
    """Entry point for `python -m code.PreProcessScript_TENET_TF`."""
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print("Usage: python -m code.PreProcessScript_TENET_TF <species>", file=sys.stderr)
        sys.exit(1)
    species = argv[0]
    run_preprocess(species)


if __name__ == "__main__":
    main()
