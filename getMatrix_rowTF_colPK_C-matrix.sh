
##USAGE: sh splitMatrix_rowTF-colPK_C-matrix.sh

#cd /home/Data/tenet_result_kyu/

for TE_mtrx in `ls /home/Data_Drive_8TB/Data/tenet_result_kyu/cell*_normal/TE_result_matrix.txt`; do
	echo "-------------------"
	echo $TE_mtrx

        bash splitMatrix_rowTF ${TE_mtrx}
        #===> ${TE_mtrx%.txt}_rowTF.txt

	bash splitMatrix_colPK ${TE_mtrx%.txt}_rowTF.txt
	#===> ${TE_mtrx%.txt}_rowTF_colPK.txt

	rm ${TE_mtrx%.txt}_rowTF.txt
done
