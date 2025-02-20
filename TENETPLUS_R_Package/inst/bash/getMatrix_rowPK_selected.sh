
##USAGE: sh splitMatrix_rowTF-colPK_C-matrix.sh

#cd /home/Data/tenet_result_kyu/

for TE_mtrx in `ls $1`; do
	echo "-------------------"
	echo $TE_mtrx
	selected="$2"
  echo $selected

        bash splitMatrix_rowPK_select ${TE_mtrx} ${selected}
        #===> ${TE_mtrx%.txt}_rowTF.txt

	#bash splitMatrix_colPK ${TE_mtrx%.txt}_rowTF.txt
	#===> ${TE_mtrx%.txt}_rowTF_colPK.txt

	#rm ${TE_mtrx%.txt}_rowTF.txt
done
