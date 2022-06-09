#!/usr/bin/bash

### USAGE: sh XXX.sh

Path=/home/Data_Drive_8TB/Data/tenet_result_kyu
cd $Path

for kResult in `ls */HOMER_sizeGvn_2022Mar04*/HOMER-*/knownResults.txt`
do
	#kResult=Bcell_merge_TF/HOMER_sizeGvn_2022Feb23_tmp/HOMER-rowTF_colPK-ARID5B/knownResults.txt
	
	TF=`echo $kResult| sed 's/\/knownResults.txt//g' | grep -o '[^-]*$' `
	#grep -i "^$TF" $kResult #results starting with $TF
	tf=`echo $TF | sed 's/[0-9]$//'`; #tf=${TF%%+([[:digit:]])} #removed the last digits to search for TF families ex) IRF3 --> IRF* 

	echo "---------- $TF || $tf "'\t'$kResult
	grep --color -i $tf $kResult
done
