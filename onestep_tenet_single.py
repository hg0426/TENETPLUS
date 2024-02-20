
import os, glob

cmd = "./TENET_Plus_single ExprMatMerged_Counts_DEGDAR.csv 1 Tenet_pseudo.txt Tenet_cellselect.txt 1 human 1"
#cmd = "./TENET_Plus_single merged_expression_data.csv 1 trajectory.txt cell_select.txt 1 human 1"
print(cmd) 
os.system(cmd)

cmd = "python make_GRN_new.py TE_result_matrix_rowTF_colGN.txt 0.01"
print(cmd)
os.system(cmd)

cmd = "python trim_indirect.py TE_result_matrix_rowTF_colGN.fdr0.01.sif -0.01"
print(cmd)
os.system(cmd)

cmd = "python countOutdegree.py TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect-0.01.sif"
print(cmd)
os.system(cmd)

#===========================================================================================

cmd = "python make_GRN_new.py TE_result_matrix_rowPeak_colGN.txt 0.01"
print(cmd)
os.system(cmd)
