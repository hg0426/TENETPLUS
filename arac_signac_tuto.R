library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

pbmc <- readRDS("C:/Users/LSCS/Desktop/kyu/R_proj/ATAC/RDS/pbmc_atac.rds")
annotations <- readRDS("C:/Users/LSCS/Desktop/kyu/R_proj/ATAC/RDS/annotations.rds")

pbmc
pbmc[['ATAC']]
granges(pbmc)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(pbmc) <- annotations

