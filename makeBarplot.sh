#!/bin/bash

### USAGE: sh makeBarplot.sh
### Key commandline: Rscript makeBarplot.R [fullPath_degree.txt]

Path=/home/Data/tenet_result_kyu
cd $Path
 
#for DegreeTxt in `ls *normal/TE_result_matrix_colPK_rowTF.fdr0.01.trimIndirect0.0.sif.*degree.txt`
for DegreeTxt in `ls *normal/TE_result_matrix*.fdr0.01.trimIndirect0.0.sif.*degree.txt`
do
	Rscript ./script/makeBarplot.R $DegreeTxt
done

