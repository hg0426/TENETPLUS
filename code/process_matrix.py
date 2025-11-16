import os
import sys
import time

import pandas as pd

from code.path_utils import coerce_input_path, coerce_output_path


def print_error(message: str) -> None:
    """Print an error message and exit."""
    print(message, file=sys.stderr)
    sys.exit(1)


def detect_file_format(file_path):
    """Return a simple format label from a file path."""
    _, ext = os.path.splitext(str(file_path).lower())
    if ext == ".csv":
        return "csv"
    if ext == ".tsv":
        return "tsv"
    if ext == ".txt":
        return "txt"
    if ext == ".parquet":
        return "parquet"
    print_error(
        f"Unsupported file format: '{ext}'. "
        "Supported formats are CSV, TSV, TXT, and Parquet."
    )


def read_matrix(file_path, file_format):
    """Load the input matrix according to its format."""
    try:
        if file_format in ["txt", "tsv"]:
            return pd.read_table(file_path, index_col=0)
        if file_format == "csv":
            return pd.read_csv(file_path, index_col=0)
        if file_format == "parquet":
            return pd.read_parquet(file_path)
    except FileNotFoundError:
        print_error(f"Input file '{file_path}' not found. Please check the file path.")
    except Exception as e:  # pragma: no cover - defensive
        print_error(f"Error reading input file: {e}")
    return None


def process_matrix(matrix_file: str, trajectory_file: str, cell_select_file: str) -> None:
    """Main processing routine used by TENET_Plus_for_py.sh.

    - Reads the input matrix, trajectory, and cell_select.
    - Filters cells where cell_select == 1.
    - Sorts by trajectory.
    - Writes filtered matrix/trajectory/cell_select and gene_names.
    """
    start_time = time.time()

    matrix_file_path = coerce_input_path(matrix_file)
    trajectory_file_path = coerce_input_path(trajectory_file)
    cell_select_file_path = coerce_input_path(cell_select_file)

    file_format = detect_file_format(matrix_file_path)

    df = read_matrix(matrix_file_path, file_format)
    if df is None:
        return

    trajectory = pd.read_csv(trajectory_file_path, header=None).squeeze("columns")
    cell_select = pd.read_csv(cell_select_file_path, header=None).squeeze("columns")

    print(len(df))
    print(len(cell_select))

    if len(df) != len(cell_select):
        raise ValueError(
            "The length of CSV file and cell_select file do not match."
        )

    columns_list = df.columns.tolist()

    gene_names_input_path = coerce_input_path("gene_names")
    gene_names_output_path = coerce_output_path("gene_names")
    pd.DataFrame(columns_list).to_csv(
        gene_names_input_path, index=False, header=False
    )
    pd.DataFrame(columns_list).to_csv(
        gene_names_output_path, index=False, header=False
    )

    cell_select.index = df.index
    trajectory.index = df.index

    filtered_df = df[cell_select == 1]
    filtered_trajectory = trajectory[cell_select == 1]
    trajectory_sorted = filtered_trajectory.sort_values(ascending=True).copy()
    df_sorted = filtered_df.reindex(trajectory_sorted.index)

    filtered_matrix_path = coerce_output_path("filtered_matrix.parquet")
    df_sorted.to_parquet(filtered_matrix_path)
    filtered_trajectory_path = coerce_output_path("filtered_trajectory.txt")
    trajectory_sorted.to_csv(filtered_trajectory_path, index=False, header=False)

    filtered_cellselect_path = coerce_output_path("filtered_cellselect.txt")
    with open(filtered_cellselect_path, "w", encoding="utf-8") as file:
        for _ in range(len(filtered_df)):
            file.write("1\n")

    print(f"---matrix process time : {time.time() - start_time} seconds ---")


def main(argv: list[str] | None = None) -> None:
    """Entry point for `python -m code.process_matrix`."""
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) != 3:
        print_error(
            "Usage: python -m code.process_matrix "
            "<matrix_file> <trajectory_file> <cell_select_file>"
        )
    matrix_file, trajectory_file, cell_select_file = argv
    process_matrix(matrix_file, trajectory_file, cell_select_file)


if __name__ == "__main__":
    main()
