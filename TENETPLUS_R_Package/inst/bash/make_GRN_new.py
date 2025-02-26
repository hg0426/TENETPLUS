import numpy
import statsmodels.sandbox.stats.multicomp
import scipy.stats
import sys

file_name=sys.argv[1]

ifile = open(file_name)
line = ifile.readline()
lines_target = ifile.read().split("\n")
temp = line.split()
target_name=[]
source_name=[]
#print(len(lines_target))
#print(len(temp))
for i in range(len(temp)-1):
    target_name.append(temp[i+1])

for i in range(len(lines_target)):
    lines2 = lines_target[i].split("\t")
    source_name.append(lines2[0])

cutOff=0
sourceIndex=0
TEnetwork=[]
source=[]
TE=[]
target=[]
#print(len(source_name))
#print(len(target_name))

file_name=sys.argv[1]

ifile = open(file_name)
line = ifile.readline()

for line in ifile:
    temp = line.split()
    for targetIndex in range(len(temp)-1):
        if float(temp[targetIndex+1])>cutOff:            
            source.append(source_name[sourceIndex])
            TE.append(float(temp[targetIndex+1]))
            target.append(target_name[targetIndex])
    sourceIndex=sourceIndex+1
ifile.close()

TEzscore=(TE-numpy.mean(TE))/numpy.std(TE)
TEpvalue=1-scipy.stats.norm.cdf(TEzscore)
TEfdr=statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue,alpha=0.05,method='fdr_bh')

fdrCutoff=float(sys.argv[2])
ofile = open(file_name.replace(".txt",".fdr")+str(fdrCutoff)+".sif","w")
for i in range(len(source)):
    if TEfdr[1][i]<fdrCutoff:
        ofile.write(source[i]+"\t"+str(TE[i])+"\t"+target[i]+"\n")
ofile.close()
