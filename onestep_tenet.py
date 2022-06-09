import os, glob

cmd = "./TENET expression_data.csv 10 trajectory.txt cell_select.txt 1"
print(cmd) 
os.system(cmd)

cmd = "python makeGRN.py 0.01"
print(cmd)
os.system(cmd)

#cmd = "python trim_indirect.py TE_result_matrix.fdr0.01.sif 0"
#print(cmd)
#os.system(cmd)

cmd = "python countOutdegree.py TE_result_matrix.fdr0.01.trimIndirect0.0.sif"
print(cmd)
os.system(cmd)

