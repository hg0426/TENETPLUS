library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(SeuratData)


pbmc_atac2 <- readRDS("C:/Users/LSCS/Desktop/kyu/R_proj/ATAC/RDS/pbmc_atac.rds")
pbmc_atac2 <- RunTFIDF(pbmc_atac2)
pbmc_atac2 <- FindTopFeatures(pbmc_atac2, min.cutoff = "q0")
pbmc_atac2 <- RunSVD(pbmc_atac2)
pbmc_atac2 <- RunUMAP(pbmc_atac2, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc_atac2 <- FindNeighbors(object = pbmc_atac2, reduction = 'lsi', dims = 2:30)
pbmc_atac2 <- FindClusters(object = pbmc_atac2, verbose = FALSE, algorithm = 3)

DimPlot(pbmc_atac2)
DimPlot(pbmc_atac2, group.by = "seurat_annotations")


pbmc.markers <- FindAllMarkers(pbmc_atac2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers2 <- FindAllMarkers(pbmc_atac2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
pbmc.markers3 <- FindAllMarkers(pbmc_atac2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)


matrix_file <- pbmc_atac2@assays[["ATAC"]]@counts
dim(matrix_file)
matrix_file<-matrix_file[unique(pbmc.markers2$gene),]
dim(matrix_file)

matrix_file <- as.matrix(t(matrix_file))

dim(matrix_file)

length(colnames(matrix_file))

write.table(matrix_file,"matrix_file_markergene.csv",quote = F, sep = ",",  row.names = F)

matrix_gene <- read_csv("C:/Users/LSCS/Desktop/tenet/matrix.csv")

matrix_gene <- as.matrix(matrix_gene)

dim(matrix_gene)
dim(matrix_file)

merge_matrix <- cbind(matrix_gene,matrix_file)
colnames(merge_matrix)[1] <- "x"
dim(merge_matrix)

write.table(merge_matrix,"merged_matrix.csv",quote = F, sep = ",",  row.names = F)
