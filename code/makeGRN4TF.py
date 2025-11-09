import sys
from pathlib import Path

import numpy
import scipy.stats
import statsmodels.sandbox.stats.multicomp

from code.path_utils import coerce_output_path, locate_file


def load_gene_names() -> list[str]:
    path = locate_file("gene_names")
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def load_te_edges(path: Path):
    cut_off = 0.0
    sources, targets, tes = [], [], []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            temp = line.strip().split(",")
            if len(temp) < 3:
                continue
            dst_val = float(temp[1])
            if dst_val > cut_off:
                sources.append(int(temp[0]))
                tes.append(float(temp[2]))
                targets.append(int(temp[1]))
    return sources, targets, tes


def main() -> None:
    gene_names = load_gene_names()
    te_path = locate_file("TE_result_all.csv")
    sources, targets, tes = load_te_edges(te_path)

    if not tes:
        print("No TE entries above cutoff; nothing to write.")
        return

    te_arr = numpy.asarray(tes, dtype=float)
    zscores = (te_arr - te_arr.mean()) / te_arr.std(ddof=0)
    pvals = 1 - scipy.stats.norm.cdf(zscores)
    _, qvals, _, _ = statsmodels.sandbox.stats.multicomp.multipletests(
        pvals,
        alpha=0.05,
        method='fdr_bh',
    )

    fdr_cutoff = float(sys.argv[1])
    out_name = te_path.name.replace(".csv", ".fdr") + f"{fdr_cutoff}.sif"
    out_path = coerce_output_path(out_name)
    with out_path.open("w", encoding="utf-8") as handle:
        for src, tgt, te_value, q in zip(sources, targets, te_arr, qvals):
            if q < fdr_cutoff:
                handle.write(f"{gene_names[src-1]}\t{te_value}\t{gene_names[tgt-1]}\n")


if __name__ == "__main__":
    main()
