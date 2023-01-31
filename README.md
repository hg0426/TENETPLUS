# TENETPLUS

A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data 

## multiome version (non-multiome version is under construction)
## TENET Time Predict in SingleCore

-0.0006555 x Gene^2(or TF x Gene)+0.000006553 x Gene^2 (or TF x Gene) x selected cell count

<div>
	
![TENET_PLUSE_workflow](https://user-images.githubusercontent.com/61915842/200551135-54f726db-6863-4d52-8322-98246300a13d.PNG)


</div>


## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014

https://github.com/neocaleb/TENET


## 1. Run TENET using expression data in a csv file and pseudotime result in a text file
#### Usage

	./TENET [merged_expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length]
	
#### example

	./TENET merged_expression_data.csv 10 trajectory.txt cell_select.txt 1

## 1-1. Run TENET from TF to target using expression data in a csv file and pseudotime result in a text file
#### Usage

 	./TENET_TF [merged_expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [species]
	
#### example

	 ./TENET_TF merged_expression_data.csv 10 trajectory.txt cell_select.txt 1 human
	 
## 1-2. Run TENET from selected gene peak to target using expression data in a csv file and pseudotime result in a text file
#### Usage

 	./TENET_select [merged_expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [selcet_list]
	
#### example

	 ./TENET_select merged_expression_data.csv 10 trajectory.txt cell_select.txt 1 select_list.txt
	
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

	
	
	
	
	
	
	
	
	
	
	
