##USAGE: sh getMatrix_rowPK_colGN_GH-matrix.sh

#cd /home/Data/tenet_result_kyu/

for TE_mtrx in `ls /home/Data_Drive_8TB/Data/tenet_result_kyu/cell*_normal/TE_result_matrix.txt`; do
	echo "-------------------"
	echo $TE_mtrx

	bash splitMatrix_rowTF $TE_mtrx
	#===> ${TE_mtrx%.txt}_rowPK.txt

	bash splitMatrix_colGN ${TE_mtrx%.txt}_rowTF.txt
	#===> ${TE_mtrx%.txt}_rowPK_colGN.txt

	rm ${TE_mtrx%.txt}_rowTF.txt
done
