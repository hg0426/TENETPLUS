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

cds <- as.cell_data_set(coembed)
cds <- cluster_cells(cds,reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds)
rna_pseudo <- as.Seurat(cds, assay = NULL)


merge_trajectory <- read.table("C:/Users/LSCS/Desktop/tenet/merge_trajectory2.txt", quote="\"", comment.char="")

merge_trajectory <- as.numeric(merge_trajectory$V1)

cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- merge_trajectory

plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


rna_pseudo@meta.data$merge_trajectory <- trajecory_merge
FeaturePlot(rna_pseudo, "merge_trajectory", pt.size = 0.1)& scale_color_viridis_c()

rna_pseudo@meta.data$bcell_trajectory <- trajecory_bcell
FeaturePlot(rna_pseudo, "bcell_trajectory", pt.size = 0.1)& scale_color_viridis_c()

rna_pseudo@meta.data$mono_trajectory <- trajecory_mono
FeaturePlot(rna_pseudo, "mono_trajectory", pt.size = 0.1)& scale_color_viridis_c()

rna_pseudo@meta.data$cd4_trajectory <- trajecory_cd4
FeaturePlot(rna_pseudo, "cd4_trajectory", pt.size = 0.1)& scale_color_viridis_c()    

rna_pseudo@meta.data$cd8_trajectory <- trajecory_cd8
FeaturePlot(rna_pseudo, "cd8_trajectory", pt.size = 0.1)& scale_color_viridis_c()

library(RColorBrewer)
rna@meta.data$merge_trajectory <- merge_trajectory

FeaturePlot(rna, "merge_trajectory", pt.size = 0.1)

FeaturePlot(rna, "merge_trajectory") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 15, name = "Rd")))

FeaturePlot(rna_pseudo, "merge_trajectory") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 15, name = "RdBu")))


rna <- subset(coembed, orig.ident == "RNA")
cd8 <- subset(coembed, seurat_annotations == "CD8 Naive" | seurat_annotations == "CD8 TEM_1" | seurat_annotations == "CD8 TEM_2"|seurat_annotations=="MAIT"|seurat_annotations=="gdt"|seurat_annotations=="NK")
DimPlot(cd8, group.by = "seurat_annotations",label = TRUE)
cd8_barcode <- WhichCells(cd8)

file <- data.frame(barcode=colnames(rna_pseudo), # barcodes
                   pseudotime=rna_pseudo@meta.data$monocle3_pseudotime)

pseudotime_init <- subset(file$barcode, file$pseudotime!='Inf')

merge_barcode <- append(cd8_barcode,pseudotime_init)

sort_barcoed <- list()

length(merge_barcode)

for (i in 1:length(merge_barcode)){
  if (duplicated(merge_barcode)[i] == TRUE){
    sort_barcoed[i] <- merge_barcode[i]
    print(i)
  }
}
library(purrr)

compact(sort_barcoed)
sort_barcoed <- as.character(sort_barcoed)
sort_barcoed <- as.data.frame(sort_barcoed)
sort_barcoed <- drop_na(sort_barcoed)
write.table(sort_barcoed,"barcode_sort.txt",row.names = F,col.names = F,quote = F) 
which(duplicated(merge_barcode),TRUE)
sort_barcoed[5379]
sort_barcoed <- as.data.frame(sorted_barcode)
sort_barcoed <- as.character(sorted_barcode$V1)

sort_barcoed2 <- as.character(sorted_barcode2$V1)
barcodes_both <- as.character(barcodes_both$V1)
barcodes_only_annotation <- as.character(barcodes_only_annotation$V1)
barcodes_only_pseudo <- as.character(barcodes_only_pseudo$V1)
barcodes_annotation <- as.character(barcodes_annotation$V1)
barcodes_pseudo <- as.character(barcodes_pseudo$V1)


DimPlot(coembed, cells.highlight = list(barcodes_annotation), cols.highlight = c("royalblue4"), cols = "gray", order = TRUE)
DimPlot(coembed, cells.highlight = list(barcodes_pseudo), cols.highlight = c("salmon2"), cols = "gray", order = TRUE)

DimPlot(coembed, cells.highlight = list(barcodes_both,barcodes_only_annotation,barcodes_only_pseudo), cols.highlight = c("salmon2","royalblue4","cyan4"), cols = "gray", order = TRUE)



only_pseudo <- list()

merge_psuedo_cd8 <- append(cd8_barcode,pseudotime_init)

for (i in 1:length(merge_barcode)){
  if (duplicated(merge_barcode)[i] == FALSE){
    only_pseudo[i] <- merge_barcode[i]
    print(i)
  }
}
only_pseudo <- as.character(only_pseudo)
only_pseudo <- as.data.frame(only_pseudo)

write.table(pseudotime_init,"barcode_sort2.txt",row.names = F,col.names = F,quote = F) 


all_barcoded <- WhichCells(coembed)








