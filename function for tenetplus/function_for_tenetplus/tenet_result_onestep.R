library(stringr)
library(dplyr)
library(ggplot2)
library(gplots)
library(VennDiagram)
library(RColorBrewer)
library(pheatmap)
source("C:/Users/LSCS/Desktop/r/function_for_tenet.R")

total_peak_gene <- read.delim("~/total_peak_gene_hg38.txt")
target <- c("CD4","CD8","mono","Bcell")

for (i in target){
  print(i)
  peak_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/0421/",i,"/TE_result_matrix_rowTF_colPK.fdr0.01.sif"), header=FALSE)
  gene_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/0421/",i,"/TE_result_matrix_rowTF_colGN.fdr0.01.sif"), header=FALSE)
  PAG_list <- trimming_PAG(total_peak_gene = total_peak_gene,peak_list = peak_list,gene_list = gene_list)
}

TE_barplot(gene_list = gene_list, peak_list = peak_list, highlight_gene = c("LRRFIP1","KLF2","ELF1","Adasd"),fill_col = "red", highlight_col = "blue")

sample_matrix <- read.table("~/sample_matix.txt", quote="\"", comment.char="")
sample_pseudotime <- read.table("~/pseudotime.txt", quote="\"", comment.char="")
sample_gene_list <- read.table("~/gene_list.txt", quote="\"", comment.char="")
sample_gene_list <- as.character(unlist(sample_gene_list))

aa <- count_degree(gene_list,decreasing = F,degree = "out")


sample_gene_list <- c("TCF7","ELF1")
pseudotime_heatmap(matrix = sample_matrix,selected_gene = sample_gene_list,pseudotime = sample_pseudotime,
                   span = 0.6,use_pseudotime_origin = F,use_z_score = T,
                   order_pseudo_score = T,max_min = F, p_max = 1, p_min = -0.8, out_result = F, filename = "test")


