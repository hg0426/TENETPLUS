# TENETPLUS

A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data 

## multiome version (non-multiome version is under construction)

## Method
<div>


![tenetplus_workflow](https://github.com/hg0426/TENETPLUS/assets/61915842/705322aa-1e17-4deb-baef-b6fb4395f0c4)



</div>


## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014

https://github.com/neocaleb/TENET


## 1. Run TENETPLUS from TF to target and peak source once using expression data in a csv file and pseudotime result in a text file after make result matrix split matrix
#### Usage

 	./TENET_Plus [merged_expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [species]
	
#### example

	./TENET_Plus merged_expression_data.csv 10 trajectory.txt cell_select.txt 1 human
	
#### output

	TE_result_matrix_rowTF_colPK.txt
	TE_result_matrix_rowTF_colGN.txt
	TE_result_matrix_rowPeak_colGN.txt
	TE_result_matrix.txt
 
## 2. Reconstructing GRN

#### Usage

 	python make_GRN_new.py [AB matrix] [cut-off]
	
#### example

 	python make_GRN_new.py TE_result_matrix_rowTF_colGN.txt 0.01 
 	python make_GRN_new.py TE_result_matrix_rowTF_colPK.txt 0.01

###### Output file
	TE_result_matrix_rowTF_colPK.sif
	TE_result_matrix_rowTF_colGN.sif
 	TE_result_matrix_rowPeak_colGN.sif
	
## 3. Counting out-degree of a given GRN

#### Usage
 	python countOutdegree.py [name of GRN]
#### Example
 	python countOutdegree.py TE_result_matrix_rowTF_colGN.sif  
 	python countOutdegree.py TE_result_matrix_rowTF_colPK.sif

###### Output file
	 TE_result_matrix_rowTF_colPK.sif.outdegree.txt

	
	
	
	
	
	
	
	
	
	
	
