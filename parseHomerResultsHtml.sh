#homerHTML=/home/Data/tenet_result_kyu/Bcell_merge_normal/HOMER_sizeGvn_2022Feb23/HOMER-colPK_rowTF-ARID5B/homerResults.html 

#homerHTML=/Users/jylee43/SSU/TENET+/4_MotifSearch_sizegiven_2022Feb07_HyunkKim/HOMER_Bcell_peak.fdr0.01_BCL11A_sizegiven/homerResults.html

#awk -F"\t" '{if ($1 ~ /^<\/TD><TD>/ && $1 !~ /^<\/TD><TD><svg/) print $0; else if ($1 ~ /\*<\/FONT><\/TD><TD>/) exit}' $homerHTML |sed 's/<\/TD><TD>/\t/g' | cut -d"<" -f 1 | awk '{print NR$0}'


awk -F"\t" '{if ($1 ~ /^<\/TD><TD>/ && $1 !~ /^<\/TD><TD><svg/) print $0; else if ($1 ~ /\*<\/FONT><\/TD><TD>/) exit}' $homerHTML |sed 's/<\/TD><TD>/\t/g' | cut -d"<" -f 1 | awk '{print NR$0}' > ${homerHTML%.html}.txt
#filtering p-value cutoff

for num in `cut -f 1 ${homerHTML%.html}.txt`
do
	echo -e "$num\c"
	indvMotif=${homerHTML%.html}"/motif"$num".info.html"
	echo $indvMotif

