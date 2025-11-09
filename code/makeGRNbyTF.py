import sys
from pathlib import Path

import numpy
import scipy.stats
import statsmodels.sandbox.stats.multicomp

from code.path_utils import coerce_input_path, coerce_output_path


def load_tf_list(species: str) -> set[str]:
    tf_path = coerce_input_path(
        f"GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
    )
    with open(tf_path, encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def load_matrix(matrix_path: Path, tf_list: set[str]):
    sources = []
    targets = []
    te_vals = []
    with matrix_path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        gene_names = header[1:]
        for source_idx, line in enumerate(handle):
            if source_idx >= len(gene_names):
                break
            if gene_names[source_idx] not in tf_list:
                continue
            parts = line.split()
            for target_idx in range(len(parts) - 1):
                value = float(parts[target_idx + 1])
                if value > 0:
                    sources.append(gene_names[source_idx])
                    targets.append(gene_names[target_idx])
                    te_vals.append(value)
    return sources, targets, te_vals


def main() -> None:
    species = sys.argv[1]
    tf_list = load_tf_list(species)

    matrix_path = coerce_output_path("TE_result_matrix.txt")
    if not matrix_path.exists():
        raise FileNotFoundError(f"TE matrix not found at {matrix_path}")

    sources, targets, te_vals = load_matrix(matrix_path, tf_list)
    te_array = numpy.asarray(te_vals, dtype=float)
    if te_array.size == 0:
        print("No TE entries for selected TFs; nothing to write.")
        return

    zscores = (te_array - te_array.mean()) / te_array.std(ddof=0)
    pvals = 1 - scipy.stats.norm.cdf(zscores)
    _, qvals, _, _ = statsmodels.sandbox.stats.multicomp.multipletests(
        pvals,
        alpha=0.05,
        method="fdr_bh",
    )

    fdr_cutoff = float(sys.argv[2])
    out_name = matrix_path.name.replace(".txt", ".byGRN.fdr") + f"{fdr_cutoff}.sif"
    out_path = matrix_path.with_name(out_name)
    with out_path.open("w", encoding="utf-8") as handle:
        for src, tgt, te_val, q in zip(sources, targets, te_array, qvals):
            if q < fdr_cutoff:
                handle.write(f"{src}\t{te_val}\t{tgt}\n")


if __name__ == "__main__":
    main()
