#!/usr/bin/sh

### USAGE: sh XXX.sh

Path=/home/Data_Drive_8TB/Data/tenet_result_kyu
cd $Path 

#homerHTML=/home/Data/tenet_result_kyu/Bcell_merge_normal/HOMER_sizeGvn_2022Feb23/HOMER-colPK_rowTF-ARID5B/homerResults.html
i=0

for homerHTML in `ls $Path/*/*/HOMER-*/homerResults.html`
#for homerHTML in `ls $Path/Bcell_merge_normal/HOMER_sizeGvn_2022Feb23_tmp/HOMER-*/homerResults.html`
do
	indvMotif=${homerHTML%.html}/motif1.info.html
	if [ -f "$indvMotif" ]
	then
		i=$((i+1))
		bash $Path/script/convertHomerResultsHtml2Txt $homerHTML	
		#echo "$i: $homerHTML"
		echo "$i: ${homerHTML%.html}.txt"
	fi
done

