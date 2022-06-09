library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(SeuratData)

coembed <- readRDS("D:/shared/coembed.rds")

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))


#================================================ATAC====================
atac <- subset(coembed, orig.ident == "ATAC")

cds <- as.cell_data_set(atac)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)

#===================================================CD8================
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 


file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)



write.table(file,"trajecory_cd8_atac.txt",row.names = F,col.names = F,quote = F)
#===================================================CD4================
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 


file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)



write.table(file,"trajecory_cd4_atac.txt",row.names = F,col.names = F,quote = F)
#===================================================Bcell================
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 


file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)



write.table(file,"trajecory_Bcell_atac.txt",row.names = F,col.names = F,quote = F)
#===================================================mono================
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
pbmc.pseudo <- as.Seurat(cds, assay = NULL)

pbmc.pseudo@meta.data$monocle3_pseudotime 


file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)



write.table(file,"trajecory_mono_atac.txt",row.names = F,col.names = F,quote = F)




#================================================RNA====================

rna <- subset(coembed, cells = cluster.cells)

cds <- as.cell_data_set(rna)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
pbmc.pseudo <- as.Seurat(cds, assay = NULL)
pbmc.pseudo@meta.data$monocle3_pseudotime 

colnames(pbmc.pseudo)
dim(pbmc.pseudo)

file <- data.frame(pbmc.pseudo@meta.data$monocle3_pseudotime)


write.table(file,"trajecory_cd8_rna.txt",row.names = F,col.names = F,quote = F)

trajecory_cd8_rna <- file
