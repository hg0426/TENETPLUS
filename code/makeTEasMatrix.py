from __future__ import annotations

from code.path_utils import locate_file, resolve_output


def load_gene_names(path) -> list[str]:
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def main() -> None:
    gene_names_path = locate_file("gene_names")
    gene_name = load_gene_names(gene_names_path)

    te_matrix = [[0 for _ in gene_name] for _ in gene_name]

    te_input_path = locate_file("TE_result_all.csv")
    with open(te_input_path, encoding="utf-8") as handle:
        for line in handle:
            temp = line.strip().split(",")
            if len(temp) < 3:
                continue
            src = int(temp[0]) - 1
            dst = int(temp[1]) - 1
            value = float(temp[2])
            te_matrix[src][dst] = value
            if len(temp) > 3:
                te_matrix[dst][src] = float(temp[3])

    out_path = resolve_output("TE_result_matrix.txt")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("TE")
        for name in gene_name:
            handle.write(f"\t{name}")
        for i, row in enumerate(te_matrix):
            handle.write(f"\n{gene_name[i]}")
            for value in row:
                handle.write(f"\t{value}")


if __name__ == "__main__":
    main()
