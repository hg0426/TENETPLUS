import sys
from pathlib import Path

import numpy
import scipy.stats
import statsmodels.sandbox.stats.multicomp

from code.path_utils import coerce_output_path


def load_matrix(matrix_path: Path):
    with matrix_path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        gene_name = header[1:]
        cut_off = 0.0
        source = []
        te_vals = []
        targets = []
        for source_index, line in enumerate(handle):
            temp = line.split()
            for target_index in range(len(temp) - 1):
                value = float(temp[target_index + 1])
                if value > cut_off:
                    source.append(gene_name[source_index])
                    te_vals.append(value)
                    targets.append(gene_name[target_index])
    return gene_name, source, targets, te_vals


def main() -> None:
    matrix_path = coerce_output_path("TE_result_matrix.txt")
    if not matrix_path.exists():
        raise FileNotFoundError(f"TE matrix not found at {matrix_path}")
    _, source, target, te_vals = load_matrix(matrix_path)

    te_array = numpy.asarray(te_vals, dtype=float)
    if te_array.size == 0:
        print("No TE entries greater than zero; nothing to write.")
        return

    zscores = (te_array - te_array.mean()) / te_array.std(ddof=0)
    pvals = 1 - scipy.stats.norm.cdf(zscores)
    _, qvals, _, _ = statsmodels.sandbox.stats.multicomp.multipletests(
        pvals,
        alpha=0.05,
        method="fdr_bh",
    )

    fdr_cutoff = float(sys.argv[1])
    out_name = matrix_path.name.replace(".txt", ".fdr") + f"{fdr_cutoff}.sif"
    out_path = matrix_path.with_name(out_name)
    with out_path.open("w", encoding="utf-8") as handle:
        for src, tgt, te_val, q in zip(source, target, te_array, qvals):
            if q < fdr_cutoff:
                handle.write(f"{src}\t{te_val}\t{tgt}\n")


if __name__ == "__main__":
    main()
