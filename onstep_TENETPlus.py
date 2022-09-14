import os, glob

tenet = input("Enter as example :  ./TENET [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] \n                  ")
cmd = tenet
print(cmd) 
os.system(cmd)

result_dic = os.getcwd() + "/TE_result_matrixl.txt"

cmd = "sh getMatrix_rowTF_colPK_C-matrix.sh " + result_dic
print(cmd)
os.system(cmd)

cmd = "sh getMatrix_rowTF_colGN_AB-matrix.sh " result_dic
print(cmd)
os.system(cmd)

cmd = "python make_GRN_new.py " + "TE_result_matrix_rowTF_colGN.txt"
print(cmd)
os.system(cmd)

cmd = "python make_GRN_new.py " + "TE_result_matrix_rowTF_colPK.txt"
print(cmd)
os.system(cmd)

cmd = "python countOutdegree.py TE_result_matrix_rowTF_colPK.sif"
print(cmd)
os.system(cmd)

cmd = "python countOutdegree.py TE_result_matrix_rowTF_colGN.sif"
print(cmd)
os.system(cmd)
