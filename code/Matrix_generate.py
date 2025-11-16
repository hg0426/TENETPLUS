import os
import sys

import duckdb
import pandas as pd

from code.path_utils import locate_file, resolve_output


def generate_matrices(input_parquet: str, species: str, split_type: str) -> None:
    """Split TE_result_all.parquet into mode- or pair-mode-specific matrices."""
    # Resolve input; if the requested file is missing, fall back to TE_fast.parquet.
    input_parquet_path = locate_file(input_parquet)
    if not input_parquet_path.exists():
        fallback = locate_file("TE_fast.parquet")
        if fallback.exists():
            print(
                f"[Matrix_generate] {input_parquet} not found; falling back to {fallback.name}."
            )
            input_parquet_path = fallback
        else:
            raise FileNotFoundError(
                f"No TE result file found: tried {input_parquet} and TE_fast.parquet"
            )
    gene_names_file = locate_file("gene_names")

    # Build a 1-based gene index mapping (id -> name) for joining.
    gene_names_df = pd.read_csv(gene_names_file, header=None)
    gene_names_df.reset_index(drop=True, inplace=True)
    gene_names_df.index += 1
    gene_names_df.index.name = "id"
    gene_names_df.columns = ["list"]

    pair_mode = os.getenv("TENET_PAIR_MODE", "default").strip().lower()

    # Use DuckDB to stream from Parquet and write filtered Parquet files
    # without materialising the full TE_result_all in pandas.
    con = duckdb.connect()
    try:
        con.register("gene_map", gene_names_df.reset_index())

        input_str = str(input_parquet_path).replace("'", "''")

        def copy_full_matrix(output_filename: str) -> None:
            """Write full Source×Target matrix, mapping indices to names."""
            out_path = resolve_output(output_filename)
            out_str = str(out_path).replace("'", "''")
            con.execute(
                f"""
                COPY (
                    SELECT
                        COALESCE(g_src.list, CAST(base.Source AS VARCHAR)) AS Source,
                        COALESCE(g_tgt.list, CAST(base.Target AS VARCHAR)) AS Target,
                        base.TE
                    FROM read_parquet('{input_str}') AS base
                    LEFT JOIN gene_map AS g_src
                        ON try_cast(base.Source AS BIGINT) = g_src.id
                    LEFT JOIN gene_map AS g_tgt
                        ON try_cast(base.Target AS BIGINT) = g_tgt.id
                )
                TO '{out_str}'
                (FORMAT 'parquet')
                """
            )

        # When using full pair modes, avoid TF/peak-specific rowTF/rowPeak naming
        # and instead write concise, mode-aware outputs.
        if pair_mode in ("all_gene", "gene_only"):
            print(
                "[Matrix_generate] pair_mode=all_gene: writing full gene×gene matrix."
            )
            copy_full_matrix("TE_GN_GN.parquet")
            print("Matrix splitting completed (GN_GN).")
            return
        if pair_mode in ("all_feature", "all", "full", "gene_peak"):
            print(
                "[Matrix_generate] pair_mode=all_feature: writing full feature×feature matrix."
            )
            copy_full_matrix("TE_all_features.parquet")
            print("Matrix splitting completed (all_features).")
            return

        # Default pair_mode: materialise a compact in-memory table once and
        # slice it into mode-specific matrices using SQL.
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE te_named AS
            SELECT
                COALESCE(g_src.list, CAST(base.Source AS VARCHAR)) AS Source,
                COALESCE(g_tgt.list, CAST(base.Target AS VARCHAR)) AS Target,
                base.TE
            FROM read_parquet('{input_str}') AS base
            LEFT JOIN gene_map AS g_src
                ON try_cast(base.Source AS BIGINT) = g_src.id
            LEFT JOIN gene_map AS g_tgt
                ON try_cast(base.Target AS BIGINT) = g_tgt.id
            """
        )

        def copy_query(query_sql: str, output_filename: str) -> None:
            out_path = resolve_output(output_filename)
            out_str = str(out_path).replace("'", "''")
            con.execute(
                f"""
                COPY (
                    {query_sql}
                )
                TO '{out_str}'
                (FORMAT 'parquet')
                """
            )

        peak_df = None
        tf_df = None

        # Peak list is required for rowPeak / peak->peak modes.
        if split_type in {"1", "5", "6"}:
            peak_list_file = locate_file("TE_peak_list.txt")
            try:
                with open(peak_list_file, "r", encoding="utf-8") as file:
                    peaks = [line.strip() for line in file if line.strip()]
            except FileNotFoundError:
                raise FileNotFoundError(f"File not found: {peak_list_file}")
            peak_df = pd.DataFrame({"name": peaks})
            con.register("peak_list", peak_df)

        # TF list is required for all modes except pure rowPeak (mode 5).
        if split_type != "5":
            tf_list_file = locate_file(
                f"GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
            )
            try:
                with open(tf_list_file, "r", encoding="utf-8") as file:
                    tfs = [line.strip() for line in file if line.strip()]
            except FileNotFoundError:
                raise FileNotFoundError(f"File not found: {tf_list_file}")
            tf_df = pd.DataFrame({"name": tfs})
            con.register("tf_list", tf_df)

        # rowPeak splits (peak as source)
        if split_type in {"1", "5"} and peak_df is not None:
            print("Splitting Matrix for TENET_Plus(rowPeak_colGN)")
            copy_query(
                """
                SELECT te.Source, te.Target, te.TE
                FROM te_named AS te
                JOIN peak_list AS pl ON te.Source = pl.name
                WHERE te.Target NOT LIKE 'chr%'
                """,
                "TE_PK_GN.parquet",
            )

        if split_type in {"6"} and peak_df is not None:
            print("Splitting Matrix for TENET_Plus(rowPeak_colPK)")
            copy_query(
                """
                SELECT te.Source, te.Target, te.TE
                FROM te_named AS te
                JOIN peak_list AS pl ON te.Source = pl.name
                WHERE te.Target LIKE 'chr%'
                """,
                "TE_PK_PK.parquet",
            )

        # rowTF splits (TF as source)
        if split_type != "5" and tf_df is not None:
            print("Splitting Matrix for rowTF")
            if split_type not in {"3"}:
                print("Splitting Matrix by colGN")
                copy_query(
                    """
                    SELECT te.Source, te.Target, te.TE
                    FROM te_named AS te
                    JOIN tf_list AS tf ON te.Source = tf.name
                    WHERE te.Target NOT LIKE 'chr%'
                    """,
                    "TE_TF_GN.parquet",
                )

            if split_type in {"1", "3", "4"}:
                print("Splitting Matrix by colPK")
                copy_query(
                    """
                    SELECT te.Source, te.Target, te.TE
                    FROM te_named AS te
                    JOIN tf_list AS tf ON te.Source = tf.name
                    WHERE te.Target LIKE 'chr%'
                    """,
                    "TE_TF_PK.parquet",
                )

        print("Matrix splitting completed.")
    finally:
        con.close()


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
