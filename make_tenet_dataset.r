library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
#####
library(SeuratData)

coembed <- readRDS("D:/shared/coembed.rds")

rna <- subset(coembed, orig.ident == "RNA")


###
DimPlot(rna, label = TRUE,label.size = 5)
DimPlot(rna, group.by = "seurat_annotations",label = TRUE)

cds <- as.cell_data_set(rna)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds2)
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

rna_pseudo <- as.Seurat(cds, assay = NULL)

file <- data.frame(barcode=colnames(rna_pseudo), # barcodes
                   pseudotime=rna_pseudo@meta.data$monocle3_pseudotime) # pseudotime

file$cell_select <- ifelse(file$pseudotime!=Inf,1,0)

write.table(file$pseudotime,"trajecory_mono.txt",row.names = F,col.names = F,quote = F)
write.table(file$cell_select,"cell_select_mono.txt",row.names = F,col.names = F,quote = F)
 


#check
10412 - length(which(trajecory_mono$V1=='Inf'))
10412 - length(which(cell_select_mono$V1==0))

10412 - length(which(trajecory_bcell$V1=='Inf'))
10412 - length(which(`cell_select_bcell.(2)`$V1==0))

10412 - length(which(trajecory_cd4$V1=='Inf'))
10412 - length(which(cell_select_cd4$V1==0))

10412 - length(which(trajecory_cd8$V1=='Inf'))
10412 - length(which(cell_select_cd8$V1==0))

write.table(,"twaat.txt",row.names = F,col.names = F,quote = F)

length(which(trajecory_bcell$V1!='Inf'))
length(which(cell_select_bcell$V1==1))
length(which(trajecory_bcell$V1!='Inf' & cell_select_bcell$V1==1))

length(which(trajecory_cd8$V1!='Inf'))
length(which(cell_select_cd8$V1==1))
length(which(trajecory_cd8$V1!='Inf' & cell_select_cd8$V1==1))

length(which(trajecory_cd4$V1!='Inf'))
length(which(cell_select_cd4$V1==1))
length(which(trajecory_cd4$V1!='Inf' & cell_select_cd4$V1==1))

length(which(trajecory_mono$V1!='Inf'))
length(which(cell_select_mono$V1==1))
length(which(trajecory_mono$V1!='Inf' & cell_select_mono$V1==1))

