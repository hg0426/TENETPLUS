library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(SeuratData)

coembed <- readRDS("D:/shared/coembed.rds")

pbmc.rna <- subset(coembed, orig.ident == "RNA")

pbmc.rna <- FindNeighbors(pbmc.rna, dims = 1:50) # do it
pbmc.rna <- FindClusters(pbmc.rna, resolution = 0.8)
DimPlot(pbmc.rna, group.by = c("RNA_snn_res.0.8"))

pbmc.rna_2 <- subset(pbmc.rna, RNA_snn_res.0.8 == "2")
pbmc.rna_2 <- subset(pbmc.rna_2, seurat_annotations == "CD8 Naive")

p1 <- DimPlot(coembed, group.by = c("seurat_annotations"), label = T)
p2 <- DimPlot(pbmc.rna_2, group.by = c("seurat_annotations"))

p1 | p2

cluster.cells <- WhichCells(object = pbmc.rna_2)


#================================================ATAC====================

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))


atac <- subset(coembed, orig.ident == "ATAC")

cds <- as.cell_data_set(atac)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
#cds <- order_cells(cds)
cds <- order_cells(cds, root_cells = cluster.cells)


plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 

colnames(pbmc.pseudo)
dim(pbmc.pseudo)

file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)



write.table(file,"trajecory_cd8_atac.txt",row.names = F,col.names = F,quote = F)

trajecory_cd8_atac <- file

#================================================RNA====================

rna <- subset(coembed, orig.ident == "RNA")

cds <- as.cell_data_set(rna)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds)
cds <- order_cells(cds, root_cells = cluster.cells)

plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 

colnames(pbmc.pseudo)
dim(pbmc.pseudo)

file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)

file$cell_select <- ifelse(file!=Inf,1,0)

write.table(file$pbmc.pseudo.meta.data.monocle3_pseudotime,"trajecory_mono.txt",row.names = F,col.names = F,quote = F)
write.table(file$cell_select,"cell_select_mono.txt",row.names = F,col.names = F,quote = F)

trajecory_cd8_rna <- file

#================================================RNA + ATAC====================

cds <- as.cell_data_set(coembed)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
#cds <- order_cells(cds)
cds <- order_cells(cds, root_cells = cluster.cells)


plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 

colnames(pbmc.pseudo)
dim(pbmc.pseudo)

file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)

write.table(file,"trajecory_cd8_atac.txt",row.names = F,col.names = F,quote = F)
trajecory_cd8_rna_atac <- file

#==========================================

library(pheatmap)

a <- c("RNA","ATAC","m_RNA","m_ATAC")
b <- substr(colnames(pbmc.pseudo),1,16)

merged <- cbind(trajecory_cd8_atac,trajecory_cd8_rna, trajecory_cd8_rna_atac[1:10412,],trajecory_cd8_rna_atac[10413:20824,])
colnames(merged) <- a
rownames(merged) <- b[1:10412]

write.csv(merged,"merged.csv")

x <- data.frame(merged)
x1 <-as.numeric(gsub("Inf",0,x$RNA))
x2 <-as.numeric(gsub("Inf",0,x$ATAC))
x3 <-as.numeric(gsub("Inf",0,x$m_RNA))
x4 <-as.numeric(gsub("Inf",0,x$m_ATAC))

x__ <- cbind(x1,x2,x3,x4)

class(x1)

colnames(x__) <- a
rownames(x__) <- b[1:10412]
x_done <- as.matrix(x__)

write.csv(x_done,"aa.csv")

p1 <- heatmap(x_done, labels_col = F)

pheatmap(x_done, cluster_rows = F, scale = "none" )
pheatmap(x_done, scale = "none" , show_rownames=F)


x <- data.frame(merged)
x1 <-as.numeric(gsub("Inf",0,x$RNA))
x2 <-as.numeric(gsub("Inf",0,x$ATAC))
x3 <-as.numeric(gsub("Inf",0,x$RNA_mer))
x4 <-as.numeric(gsub("Inf",0,x$ATAC_mer))


b <- as.data.frame(x_done)
x_filter <- b[b$RNA > 0 & b$ATAC >0 & b$m_RNA >0 & b$m_ATAC >0, ]

x_filter <- as.matrix(x_filter)

p3 <- pheatmap(x_filter, scale = "none",show_rownames=F)

b <- as.data.frame(x_filter)

b1 <- order(b$RNA)

b_order <- b[b1,]

x_order <- as.matrix(b_order)

p3 <- pheatmap(x_order)
p3 <- pheatmap(x_order, scale = "none", cluster_rows = F,show_rownames=F)


#=====================================================================

rna_ <-as.numeric(gsub("Inf",0,x$RNA))
atac_ <-as.numeric(gsub("Inf",0,x$ATAC))
rna_mer <-as.numeric(gsub("Inf",0,x$RNA_mer))
atac_mer <-as.numeric(gsub("Inf",0,x$ATAC_mer))

p1
par(mfrow=c(2,2))

p2 <- hist(atac_mer)

