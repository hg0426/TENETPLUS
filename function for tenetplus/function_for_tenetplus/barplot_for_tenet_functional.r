library(ggplot2)
library(gplots)
library(VennDiagram)
library(RColorBrewer)
source("~/function_tenet_barplot.r")
target <- c("CD4","CD8","Bcell","mono")
i <- "CD4"
peak_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/0421/",i,"/TE_result_matrix_rowTF_colPK.fdr0.01.sif"), header=FALSE)
gene_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/0421/",i,"/TE_result_matrix_rowTF_colGN.fdr0.01.sif"), header=FALSE)

TE_barplot(gene_list = gene_list, peak_list = peak_list, highlight_gene = c("LRRFIP1","KLF2","ELF1","Adasd"),fill_col = "red", highlight_col = "blue")

