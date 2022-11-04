# TENETPLUS

A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data 

<div>
<img src="https://user-images.githubusercontent.com/61915842/190166344-069614e4-4d40-4591-b212-178cd6becb02.png" width="90%"></img>
</div>

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014
https://github.com/neocaleb/TENET


## 1. Run TENET using expression data in a csv file and pseudotime result in a text file
#### Usage

	./TENET [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length]
	
#### example

	./TENET expression_data.csv 10 trajectory.txt cell_select.txt 1

## 1-1. Run TENET from TF to target using expression data in a csv file and pseudotime result in a text file
#### Usage

 	./TENET_TF [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [species]
	
#### example

	 ./TENET_TF expression_data.csv 10 trajectory.txt cell_select.txt 1 human
	 
## 1-2. Run TENET from selected gene peak to target using expression data in a csv file and pseudotime result in a text file
#### Usage

 	./TENET_select [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [selcet_list]
	
#### example

	 ./TENET_select expression_data.csv 10 trajectory.txt cell_select.txt 1 select_list.txt
	
## 2 split result matrix (detect peak standard : xxx - xxxx - xxxx)

#### Usage

 	sh getMatrix_rowTF_colGN_AB-matrix.sh [TENET result matrix]  
 	sh getMatrix_rowTF_colPK_C-matrix.sh [TENET result matrix]
	
#### example

	 sh getMatrix_rowTF_colGN_AB-matrix.sh TE_result_matrix.txt  
	 sh getMatrix_rowTF_colPK_C-matrix.sh TE_result_matrix.txt
	
## 3 Reconstructing GRN

#### Usage

 	python make_GRN_new.py [AB matrix]
	
#### example

 	python make_GRN_new.py TE_result_matrix_rowTF_colGN.txt  
 	python make_GRN_new.py TE_result_matrix_rowTF_colPK.txt
	
## 4 Counting out-degree of a given GRN

#### Usage
 	python countOutdegree.py [name of GRN]
#### Example
 	python countOutdegree.py TE_result_matrix_rowTF_colGN.sif  
 	python countOutdegree.py TE_result_matrix_rowTF_colPK.sif

###### Output file
	 TE_result_matrix_rowTF_colPK.sif.outdegree.txt

	
	
	
	
	
	
	
	
	
	
	
