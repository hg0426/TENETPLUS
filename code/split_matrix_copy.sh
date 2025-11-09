#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="$REPO_ROOT/input"
OUTPUT_DIR="$REPO_ROOT/output"

cd "$OUTPUT_DIR"

echo 'Split Matrix for TENET_Plus(rowPeak_colGN)'
#split rowPeaks in list
TE_mtrx="TE_result_matrix.txt"

cut -f 1 $TE_mtrx | sort > TE_Matrix.sort.tmp

sort "$INPUT_DIR/TE_peak_list.txt" > TE_peak_list.sort.tmp
comm -12 TE_Matrix.sort.tmp TE_peak_list.sort.tmp > TE_peak_list.comm.tmp
wc -l *tmp
cut -d$'\t' -f 1 $TE_mtrx > rowNames.tmp

for i in `cat  TE_peak_list.comm.tmp`; do grep -n -w ^$i$ rowNames.tmp; done | cut -d":" -f 1 > LN.tmp
TE_peak_list=`perl -pe 's/\n/$1,/' LN.tmp |sed 's/.$//'`
TE_peak_list_Head=`echo "1,"$TE_peak_list`

#https://stackoverflow.com/questions/40842008/subset-a-file-by-row-and-column-numbers
awk -v rows="$TE_peak_list_Head" 'BEGIN{split(rows, a, ","); for (i in a) b[a[i]]} NR in b' $TE_mtrx > ${TE_mtrx%.txt}_rowPeak.txt

rm *tmp

#split col GN
TE_mtrx="TE_result_matrix_rowPeak.txt"

head -n 1 $TE_mtrx |sed 's/\t/\n/g' |grep -vE "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |wc -l

head -n 1 $TE_mtrx |sed 's/\t/\n/g' |grep -vE "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |awk '{ printf "%s,", $0 }'|sed 's/.$//' > TeGeneCols.txt

numGenes=`head -n 1 $TE_mtrx |sed 's/\t/\n/g' | grep -vE "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |wc -l`
divisor=2500
quotient=`expr $numGenes / $divisor`

echo $numGenes $divisor $quotient

for ((i=0;i<=$quotient;i++))
do
  echo $i
  cut -d "," -f $((1+$divisor*$i))-$(($divisor*$((i+1)))) TeGeneCols.txt > tmp$i
  awk -v cols=`cat tmp$i` 'BEGIN{FS=OFS="\t"; nc=split(cols, a, ",")} NR==1{for (i=1; i<=NF; i++) hdr[$i]=i} {for (i=1; i<=nc; i++) if (a[i] in hdr) printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)}' $TE_mtrx >TMP$i
done

paste TMP* > ${TE_mtrx%.txt}_colGN.txt
awk '{print NF; exit}' ${TE_mtrx%.txt}_colGN.txt

rm tmp* TMP* TeGeneCols.txt
rm TE_result_matrix_rowPeak.txt

echo 'Split Matrix rowTF'
#split row TF
TE_mtrx="TE_result_matrix.txt"

TF_list="$INPUT_DIR/GO_symbol_human_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"

cut -f 1 $TE_mtrx | sort > TE_Matrix.sort.tmp
sort $TF_list > TFs.sort.tmp
comm -12 TE_Matrix.sort.tmp TFs.sort.tmp > TFs.comm.tmp
wc -l *tmp
cut -d$'\t' -f 1 $TE_mtrx > rowNames.tmp

for i in `cat  TFs.comm.tmp`; do grep -n -w ^$i$ rowNames.tmp; done | cut -d":" -f 1 > LN.tmp
TFs=`perl -pe 's/\n/$1,/' LN.tmp |sed 's/.$//'`
HeadTFs=`echo "1,"$TFs`

awk -v rows="$HeadTFs" 'BEGIN{split(rows, a, ","); for (i in a) b[a[i]]} NR in b' $TE_mtrx > ${TE_mtrx%.txt}_rowTF.txt

rm *tmp

echo 'Split Matrix colGN'
#split col GN
TE_mtrx="TE_result_matrix_rowTF.txt"

head -n 1 $TE_mtrx |sed 's/\t/\n/g' |grep -vE "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |wc -l

head -n 1 $TE_mtrx |sed 's/\t/\n/g' |grep -vE "^chr[0-9XY]+-([0-9]+)-([0-9]+)$"|awk '{ printf "%s,", $0 }'|sed 's/.$//' > TeGeneCols.txt

numGenes=`head -n 1 $TE_mtrx |sed 's/\t/\n/g' | grep -vE "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |wc -l`
divisor=2500
quotient=`expr $numGenes / $divisor`

echo $numGenes $divisor $quotient

for ((i=0;i<=$quotient;i++))
do
  echo $i
  cut -d "," -f $((1+$divisor*$i))-$(($divisor*$((i+1)))) TeGeneCols.txt > tmp$i
  awk -v cols=`cat tmp$i` 'BEGIN{FS=OFS="\t"; nc=split(cols, a, ",")} NR==1{for (i=1; i<=NF; i++) hdr[$i]=i} {for (i=1; i<=nc; i++) if (a[i] in hdr) printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)}' $TE_mtrx >TMP$i
done

paste TMP* > ${TE_mtrx%.txt}_colGN.txt
awk '{print NF; exit}' ${TE_mtrx%.txt}_colGN.txt

rm tmp* TMP* TeGeneCols.txt

echo 'Split Matrix colPK'
#split col PK
TE_mtrx="TE_result_matrix_rowTF.txt"
head -n 1 $TE_mtrx |sed 's/\t/\n/g' |grep -P "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |wc -l

PeakCols=`head -n 1 $TE_mtrx |sed 's/\t/\n/g' |grep -P "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |awk '{ printf "%s,", $0 }'|sed 's/.$//'`
TE='TE'
TePeakCols="$TE,$PeakCols"
echo $TePeakCols > TePeakCols.txt

numPeaks=`head -n 1 $TE_mtrx |sed 's/\t/\n/g' | grep -P "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" |wc -l`
divisor=2500
quotient=`expr $numPeaks / $divisor`

echo $numPeaks $divisor $quotient

for ((i=0;i<=$quotient;i++))
do
  echo $i
  cut -d "," -f $((1+$divisor*$i))-$(($divisor*$((i+1)))) TePeakCols.txt > tmp$i
  awk -v cols=`cat tmp$i` 'BEGIN{FS=OFS="\t"; nc=split(cols, a, ",")} NR==1{for (i=1; i<=NF; i++) hdr[$i]=i} {for (i=1; i<=nc; i++) if (a[i] in hdr) printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)}' $TE_mtrx >TMP$i
done

paste TMP* > ${TE_mtrx%.txt}_colPK.txt
awk '{print NF; exit}' ${TE_mtrx%.txt}_colPK.txt

rm tmp* TMP* TePeakCols.txt

rm TE_result_matrix_rowTF.txt
