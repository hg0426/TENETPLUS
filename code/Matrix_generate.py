import sys
import pandas as pd
import numpy as np
import duckdb

from code.path_utils import locate_file, resolve_output


def load_gene_names(file_path):
    """Reads gene names from a file and returns them as a list."""
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            gene_names = [line.strip() for line in file]
        return gene_names
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")


def split_columns(gene_names, te_matrix, split_type):
    """
    Categorizes columns and returns the filtered matrix.
    """
    if split_type == "colGN":
        col_genes = [gene for gene in gene_names if not gene.startswith("chr")]
        return te_matrix[te_matrix["Target"].isin(col_genes)]
    if split_type == "colPK":
        col_peaks = [gene for gene in gene_names if gene.startswith("chr")]
        return te_matrix[te_matrix["Target"].isin(col_peaks)]
    return te_matrix


def categorize_genes(gene_names, species, split_type):
    """
    Returns a list of genes categorized based on the split_type.
    """
    if split_type == "rowPeak":
        peak_list_file = locate_file("TE_peak_list.txt")
        try:
            with open(peak_list_file, "r", encoding="utf-8") as file:
                peaks = set(line.strip() for line in file)
            selected_genes = [gene for gene in gene_names if gene in peaks]
            return selected_genes
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {peak_list_file}")
    if split_type == "rowTF":
        tf_list_file = locate_file(
            f"GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
        )
        try:
            with open(tf_list_file, "r", encoding="utf-8") as file:
                tfs = set(line.strip() for line in file)
            selected_genes = [gene for gene in gene_names if gene in tfs]
            return selected_genes
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {tf_list_file}")
    return gene_names


def save_matrix(df, output_filename):
    """Saves the DataFrame to a Parquet file in the output directory."""
    output_path = resolve_output(output_filename)
    df.to_parquet(output_path)


def generate_matrices(input_parquet: str, species: str, split_type: str) -> None:
    """Split TE_result_all.parquet into mode-specific matrices."""
    input_csv_path = locate_file(input_parquet)

    gene_names_file = locate_file("gene_names")
    gene_names = load_gene_names(gene_names_file)

    gene_names_df = pd.read_csv(gene_names_file, header=None)
    gene_names_df.reset_index(drop=True, inplace=True)
    gene_names_df.index += 1
    gene_names_df.index.name = "id"
    gene_names_df.columns = ["list"]

    con = duckdb.connect()
    con.register("gene_map", gene_names_df.reset_index())
    te_df = con.execute(
        """
        WITH base AS (
            SELECT Source, Target, TE
            FROM read_parquet(?)
        )
        SELECT
            COALESCE(g_src.list, CAST(base.Source AS VARCHAR)) AS Source,
            COALESCE(g_tgt.list, CAST(base.Target AS VARCHAR)) AS Target,
            base.TE
        FROM base
        LEFT JOIN gene_map g_src ON try_cast(base.Source AS BIGINT) = g_src.id
        LEFT JOIN gene_map g_tgt ON try_cast(base.Target AS BIGINT) = g_tgt.id
        """,
        [str(input_csv_path)],
    ).fetchdf()
    con.close()

    if split_type in ["1", "5"]:
        print("Splitting Matrix for TENET_Plus(rowPeak_colGN)")
        row_peak_genes = categorize_genes(gene_names, species, "rowPeak")
        te_matrix_row_peak = te_df[te_df["Source"].isin(row_peak_genes)]
        te_matrix_row_peak_col_gn = split_columns(
            gene_names, te_matrix_row_peak, "colGN"
        )
        save_matrix(
            te_matrix_row_peak_col_gn,
            "Local_TE_result_matrix_rowPeak_colGN.parquet",
        )

    if split_type in ["6"]:
        print("Splitting Matrix for TENET_Plus(rowPeak_colGN)")
        row_peak_genes = categorize_genes(gene_names, species, "rowPeak")
        te_matrix_row_peak = te_df[te_df["Source"].isin(row_peak_genes)]
        te_matrix_row_peak_col_pk = split_columns(
            gene_names, te_matrix_row_peak, "colPK"
        )
        save_matrix(
            te_matrix_row_peak_col_pk,
            "Local_TE_result_matrix_rowPeak_colPK.parquet",
        )

    if split_type != "5":
        print("Splitting Matrix for rowTF")
        row_tf_genes = categorize_genes(gene_names, species, "rowTF")
        te_matrix_row_tf = te_df[te_df["Source"].isin(row_tf_genes)]

        if split_type not in ["3"]:
            print("Splitting Matrix by colGN")
            te_matrix_row_tf_col_gn = split_columns(
                gene_names, te_matrix_row_tf, "colGN"
            )
            save_matrix(
                te_matrix_row_tf_col_gn,
                "Local_TE_result_matrix_rowTF_colGN.parquet",
            )

        if split_type in ["1", "3", "4"]:
            print("Splitting Matrix by colPK")
            te_matrix_row_tf_col_pk = split_columns(
                gene_names, te_matrix_row_tf, "colPK"
            )
            save_matrix(
                te_matrix_row_tf_col_pk,
                "Local_TE_result_matrix_rowTF_colPK.parquet",
            )

    print("Matrix splitting completed.")


def main(argv: list[str] | None = None) -> None:
    """Entry point for `python -m code.Matrix_generate`."""
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) != 3:
        print(
            "Usage: python -m code.Matrix_generate "
            "<input_parquet> <species> <split_type>",
            file=sys.stderr,
        )
        sys.exit(1)
    input_parquet, species, split_type = argv
    generate_matrices(input_parquet, species, split_type)


if __name__ == "__main__":
    main()
