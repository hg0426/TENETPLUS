import pandas as pd
import time
import sys

start_time = time.time()

# input file
matrix_file_path = sys.argv[1]

trajectory_file_path = sys.argv[2]

cell_select_file_path = sys.argv[3]

# read matrix
df = pd.read_csv(matrix_file_path, index_col=0)

# read trajectory
trajectory = pd.read_csv(trajectory_file_path, header=None).squeeze('columns')

# read cell_select
cell_select = pd.read_csv(cell_select_file_path, header=None).squeeze('columns')

print(len(df))
print(len(cell_select))

# validate
if len(df) != len(cell_select):
    raise ValueError("The length of CSV file and cell_select file do not match.")

# matching index
cell_select.index = df.index
trajectory.index = df.index

# filtering with cell_select
filtered_df = df[cell_select == 1]
filtered_trajectory =trajectory[cell_select == 1]

# save result
filtered_df.to_csv('filtered_matrix.csv', index=True, index_label=None)
filtered_trajectory.to_csv('filtered_trajectory.txt',index=False,header=False)

# make new cell_select(all values 1 )
with open('new_cell_select.txt', 'w') as file:
    for value in [1] * len(filtered_df):
        file.write(str(value) + '\n')

print("---matrix process time : %s seconds ---" % (time.time() - start_time))
