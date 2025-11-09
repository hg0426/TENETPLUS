import os
import sys
import time

import pandas as pd

from code.path_utils import coerce_input_path, coerce_output_path

start_time = time.time()


def print_error(message: str) -> None:
    """Print an error message and exit."""
    print(message, file=sys.stderr)
    sys.exit(1)

def detect_file_format(file_path):
    _, ext = os.path.splitext(str(file_path).lower())
    if ext == '.csv':
        return 'csv'
    elif ext == '.tsv':
        return 'tsv'
    elif ext == '.txt':
        return 'txt'
    elif ext == '.parquet':
        return 'parquet'
    else:
        print_error(f"Unsupported file format: '{ext}'. Supported formats are CSV, TSV, TXT, and Parquet.")


def read_matrix(file_path, file_format):
    try:
        if file_format in ['txt', 'tsv']:
            df = pd.read_table(file_path, index_col=0)
        elif file_format == 'csv':
            df = pd.read_csv(file_path, index_col=0)
        elif file_format == 'parquet':
            df = pd.read_parquet(file_path)
        return df
    except FileNotFoundError:
        print_error(f"Input file '{file_path}' not found. Please check the file path.")
    except Exception as e:
        print_error(f"Error reading input file: {e}")
        

# input file
matrix_file_path = coerce_input_path(sys.argv[1])

trajectory_file_path = coerce_input_path(sys.argv[2])

cell_select_file_path = coerce_input_path(sys.argv[3])

file_format = detect_file_format(matrix_file_path)

# read matrix
df = read_matrix(matrix_file_path, file_format)
#df = pd.read_csv(matrix_file_path, index_col=0)

# read trajectory
trajectory = pd.read_csv(trajectory_file_path, header=None).squeeze('columns')

# read cell_select
cell_select = pd.read_csv(cell_select_file_path, header=None).squeeze('columns')

print(len(df))
print(len(cell_select))

# validate
if len(df) != len(cell_select):
    raise ValueError("The length of CSV file and cell_select file do not match.")

columns_list = df.columns.tolist()

# Write gene_names into both input and output directories for robustness
gene_names_input_path = coerce_input_path("gene_names")
gene_names_output_path = coerce_output_path("gene_names")
pd.DataFrame(columns_list).to_csv(gene_names_input_path, index=False, header=False)
pd.DataFrame(columns_list).to_csv(gene_names_output_path, index=False, header=False)

# matching index
cell_select.index = df.index
trajectory.index = df.index

# filtering with cell_select
filtered_df = df[cell_select == 1]
filtered_trajectory =trajectory[cell_select == 1]
trajectory_sorted = filtered_trajectory.sort_values(ascending=True).copy()
df_sorted = filtered_df.reindex(trajectory_sorted.index)





# save result
#df_sorted.to_csv('filtered_matrix.csv', index=True, index_label=None)
filtered_matrix_path = coerce_output_path('filtered_matrix.parquet')
df_sorted.to_parquet(filtered_matrix_path)
filtered_trajectory_path = coerce_output_path('filtered_trajectory.txt')
trajectory_sorted.to_csv(filtered_trajectory_path, index=False, header=False)

# make new cell_select(all values 1 )
filtered_cellselect_path = coerce_output_path('filtered_cellselect.txt')
with open(filtered_cellselect_path, 'w', encoding='utf-8') as file:
    for value in [1] * len(filtered_df):
        file.write(str(value) + '\n')





# num_bins = 200
# # num_bins = len(df_sorted)
# trajectory_binned, bin_edges = pd.qcut(trajectory_sorted, q=num_bins, retbins=True, labels=False, duplicates='drop')
# trajectory_binned += 1  # Optional: start bin numbering at 1
# actual_bins = len(bin_edges) - 1
# print(f"Number of bins after qcut: {actual_bins}")
# df_sorted['Trajectory_Bin'] = trajectory_binned

# # 예시: 각 bin별로 평균 발현량 계산
# final_df = df_sorted.groupby('Trajectory_Bin').mean()

# final_df.to_csv('filtered_matrix.csv', index=True, index_label=None)
# pd.DataFrame(final_df.index).to_csv('filtered_trajectory.txt',index=False,header=False)

# # make new cell_select(all values 1 )
# with open('filtered_cellselect.txt', 'w') as file:
#     for value in [1] * len(final_df):
#         file.write(str(value) + '\n')




print("---matrix process time : %s seconds ---" % (time.time() - start_time))
