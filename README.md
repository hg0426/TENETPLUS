# TENETPLUS

A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic and epigenetic data 

## multiome version (non-multiome version is under construction)

## Method
<div>


![tenetplus_workflow](https://github.com/user-attachments/assets/36e6c47a-4f12-4a6a-a74e-1da6eea6ab32)


</div>

## Dependency

	python3
	openmpi (>4.0)
	JPype

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014

https://github.com/neocaleb/TENET


## 1. Run TENETPLUS from TF to target and peak source once using expression data in a csv file and pseudotime result in a text file after make result matrix split matrix
#### Usage

 	./TENET_Plus [merged_expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [species] [options]
##### options : 0: TENET_TF(Only RNA), 1: TENET_Plus(RNA + ATAC) all matrix, 2: TENET_Plus only rowTF colGN, 3: TENET_Plus only rowTF colPK, 4: TENET_Plus only rowPeak(cis-peaksource)
 
#### example

	./TENET_Plus merged_expression_data.csv 10 trajectory.txt cell_select.txt 1 human 1
	
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
	TE_result_matrix_rowTF_colPK.fdr0.01.sif
	TE_result_matrix_rowTF_colGN.fdr0.01.sif
 	TE_result_matrix_rowPeak_colGN.fdr0.01.sif
	

## 3. Trimming indirect edges for TF_Gene result

##### Trimming only AB matrix (original_TENET)
#### Usage
 	python trim_indirect.py [AB matrix] [cutoff]
	
#### example
 	python trim_indirect.py TE_result_matrix_rowTF_colGN.fdr0.01.sif -0.01

###### Output file
	TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect-0.01.sif
###### Parameter
	[cutoff] - A cutoff value for trimming indirect edges. Recommended range is -0.1 to 0.1


## 4. Counting out-degree of a given GRN

#### Usage
 	python countOutdegree.py [name of GRN]
#### Example
 	python countOutdegree.py TE_result_matrix_rowTF_colPK.fdr0.01.sif 
 	python countOutdegree.py TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect-0.01.sif

###### Output file
	TE_result_matrix_rowTF_colPK.fdr0.01.sif.outdegree.txt
  	TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect-0.01.sif.outdegree.txt
  

	
	
	
	
	
	
	
	
	
	
	
