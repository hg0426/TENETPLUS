library(stringr)
library(dplyr)

source("~/function_PAG_trimming.r")

total_peak_gene <- read.delim("~/total_peak_gene_hg38.txt")
target <- c("CD4","CD8","mono","Bcell")

for (i in target){
  print(i)
  peak_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/0421/",i,"/TE_result_matrix_rowTF_colPK.fdr0.01.sif"), header=FALSE)
  gene_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/0421/",i,"/TE_result_matrix_rowTF_colGN.fdr0.01.sif"), header=FALSE)
  PAG_list <- trimming_PAG(total_peak_gene = total_peak_gene,peak_list = peak_list,gene_list = gene_list)
}


att <- as.data.frame(table(gene_list$V1))[order(as.data.frame(table(gene_list$V1),)[,2], decreasing = T),]
