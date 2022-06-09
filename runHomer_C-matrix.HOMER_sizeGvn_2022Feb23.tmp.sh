#!/bin/bash

### Please activate conda environment before running this script!!
### USAGE: sh XXX.sh

### scripts to create a conda env and install tools
#conda create -n bedtools_meme_homer
#conda activate bedtools_meme_homer
#conda install -c bioconda bedtools meme homer
#conda activate bedtools_meme_homer

Path=/home/Data/tenet_result_kyu
Job=HOMER_sizeGvn_2022Feb23_tmp
#Job=HOMER_size200_2022Feb23

cd $Path 

#for Celltype in `ls -d *normal/ | cut -f1 -d'/'`
for Celltype in `ls -d Bcell_merge_normal/ | cut -f1 -d'/'`
do
	cd $Celltype
	mkdir -p $Job

	OutD=TE_result_matrix_rowTF_colPK.fdr0.01.trimIndirect0.0.sif.outdegree.txt
	Smpl=`echo $OutD| sed 's/\//./g'| sed 's/.fdr0.01.trimIndirect0.0.sif.outdegree.txt//g' | sed 's/TE_result_matrix_//g'`
	SIF=${OutD%.outdegree.txt}

	echo -e "=================================== START \c"; date "+%T"
	echo -e $Celltype"\t:\t"$Smpl"\t:\t"$SIF

	#cut -f 1: C-matrix (rowTF-colPK); cut -f 3: GH-matrix (rowPK-colGN);
	TFs=`head -n 5 $OutD| cut -f 1`
	#TFs=`awk '{if($2>400) print $1}' $OutD`  
	echo $TFs

	for TF in $TFs
	do
		echo -e "---------------- \c"; date "+%T"
		echo $TF

		BED=$Path/$Celltype/$Job/$Smpl"-"$TF.bed
		awk -v var="$TF" '{if($1==var) print $3}' $SIF| awk -F"-" '{print $1"\t"$2"\t"$3}' > $BED

		DIRN=$Path/$Celltype/$Job/"HOMER-"$Smpl"-"$TF
		echo $DIRN
		mkdir -p $DIRN
		
		findMotifsGenome.pl $BED hg38 $DIRN -size given -mset all 2>$Path/$Celltype/$Job/"HOMER-"$Smpl"-"$TF.err &
		#findMotifsGenome.pl $BED hg38 $DIRN -mset vert
		
	done
	cd $Path
	echo -e "===================================== END \c"; date "+%T"
	#sleep 4h
done

