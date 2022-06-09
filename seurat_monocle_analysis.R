###### Library
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

##### Load Data
## pbmc
setwd("C:/Users/LSCS/Desktop/tenet/files_0504")
pbmc.rna <- readRDS('final_pbmc.rna')
pbmc.markers <- readRDS('pbmc_marker.rds')
pbmc.atac <- readRDS('final_pbmc.atac')
coembed <- readRDS('coembed')
pbmc.rna <- coembed


###### Processing # for atlas data
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna,dims=1:10)
DimPlot(pbmc.rna)

pbmc.rna <- FindNeighbors(pbmc.rna, dims = 1:10) # do it pbmc data
pbmc.rna <- FindClusters(pbmc.rna,resolution = 0.8) # do it pbmc data 

###### Plotting
DimPlot(pbmc.rna, group.by = "ident",label = TRUE,label.size = 5)
DimPlot(pbmc.rna, group.by = "seurat_annotations",label = TRUE)


####### Markersu
# RNA.markers
pbmc.markers <- FindAllMarkers(pbmc.rna, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

dim(pbmc.markers) # 15480     7

length(unique(pbmc.markers$gene)) # 4991ê°? 
# loaddata and readRDSë¡œë?€?„° ?‹œ?ž‘?•¨, findcluster(resolution=0.5) -> 5238ê°?
# 5150ê°?
dim(pbmc.markers)
unique(pbmc.markers$cluster)

# ATAC.markers (peak/gene matrix) ???????????
DefaultAssay(coembed)<-"ATAC"
DefaultAssay(coembed)
da_peaks <- FindAllMarkers(object =  coembed, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(da_peaks,"E:/labs_data/tenet_data/da_peaks.rds")
da_peaks <- readRDS('E:/labs_data/tenet_data/da_peaks.rds')
dim(da_peaks) # 9560    7
length(unique(da_peaks$gene)) # 9560ê°?

###### subset RNA or ATAC
pbmc.rna<- subset(pbmc.rna, orig.ident == "RNA")
View(pbmc.rna)
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)

######## Monocle
cds <- as.cell_data_set(pbmc.rna)
cds <- cluster_cells(cds,reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "orig.ident", show_trajectory_graph = FALSE,group_label_size = 5)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size = 5)
plot_cells(cds, color_cells_by = "seurat_annotations", show_trajectory_graph = FALSE,group_label_size = 5)
plot_cells(cds, color_cells_by = "seurat_annotations", show_trajectory_graph = FALSE,group_label_size = 5,labels_per_group = "seurat_annotations")

###### learn_graph
cds <- learn_graph(cds)

####### Ordering Cells
## 1) root <-- click
cds <- order_cells(cds)

## 2) root <-- idents
cluster.cells <- WhichCells(object = pbmc.rna, idents = 3) # Seurat object
DimPlot(pbmc.rna, group.by = "ident",label = TRUE,label.size = 5)
DimPlot(pbmc.rna,cells.highlight = cluster.cells)
cds <- order_cells(cds, root_cells = cluster.cells)

## 3) root <-- annotations
index=which(unlist(FetchData(pbmc.rna,"seurat_annotations"))=="CD8 Naive")
length(index)
index=colnames(pbmc.rna)[index]
index
length(index)
DimPlot(pbmc.rna,cells.highlight = index)
cluster.cells=index
cds <- order_cells(cds, root_cells = cluster.cells)

## 4) root <-- gene
max.IL7R <- which.max(unlist(FetchData(integrated.sub, "IL7R")))
max.IL7R <- colnames(integrated.sub)[max.IL7R]
cds <- order_cells(cds, root_cells = max.IL7R)

######## Plotting Pseudotime
plot_cells(cds, color_cells_by = "pseudotime",cell_size = 1,trajectory_graph_color = "blue",label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds,color_cells_by = "seurat_annotations",group_label_size = 5,labels_per_group = "seurat_annotations",cell_size = 1,trajectory_graph_color = "BLUE",label_branch_points = FALSE,label_principal_points = FALSE)
plot_cells(cds,color_cells_by = "seurat_annotations",group_label_size = 5,cell_size = 1,trajectory_graph_color = "BLUE",label_branch_points = TRUE,label_principal_points = TRUE)

######## CDS --> Seurat object
pbmc.rna_pseudo <- as.Seurat(cds, assay = NULL) # pseudotime?´ ?žˆ?Š” cds
head(pbmc.rna_pseudo@meta.data$monocle3_pseudotime)
head(pseudotime(cds))
# pbmc.rna_pseudo@meta.data$monocle3_pseudotime <- pseudotime(cds) <if pseoduotime not in object>

## plotting using Seurat <unnecessary>
pbmc.rna_pseudo@meta.data$monocle3_pseudotime <- ifelse(pbmc.rna_pseudo@meta.data$monocle3_pseudotime==Inf,NA,pbmc.rna_pseudo@meta.data$monocle3_pseudotime)
FeaturePlot(pbmc.rna_pseudo, "monocle3_pseudotime", pt.size = 0.1)& scale_color_viridis_c()
pbmc.rna_pseudo <- as.Seurat(cds, assay = NULL)

####### Making Matrix file / markergene
### counts=raw count
## for coembed(AAACAGCCAAGGAATC-1_1)
matrix_file <- pbmc.rna_pseudo@assays[["RNA"]]@counts
# head(matrix_file)
colnames(matrix_file)<- substr(colnames(matrix_file),1,18)
# head(t(matrix_file))
# dim(t(matrix_file)) # RNA: 10412 36601/ coembed: 20824 36601
matrix_file <- as.matrix(t(matrix_file))
# head(matrix_file)
# matrix_file <- as.data.frame(matrix_file)
# dim(matrix_file)# RNA: 10412/ coembed: 20824  4991
# head(matrix_file)
# saveRDS(matrix_file,'E:/labs_data/monocle_result/seuratwrapper/matrix_file.rds')

## Using Marker genes
pbmc.markers <- readRDS('E://labs_data/tenet_data/pbmc_marker.rds')
# length(pbmc.markers$gene) 
# table(pbmc.markers$gene)
# length(unique(pbmc.markers$gene))# 4991ê°? 

matrix_file<-matrix_file[,unique(pbmc.markers$gene)]
# dim(matrix_file) # RNA: 10412  4991

write.csv(matrix_file,"E:/labs_data/monocle_result/seuratwrapper/matrix_file_markergene.csv",quote = F)

####### Making tenet file (Use Method (1))
pbmc.rna_pseudo@meta.data$monocle3_pseudotime 
colnames(pbmc.rna_pseudo)
dim(pbmc.rna_pseudo)
length(substr(colnames(pbmc.rna_pseudo),1,18))

# for pbmc.rna
file <- data.frame(barcode=colnames(pbmc.rna_pseudo), # barcodes
                   pseudotime=pbmc.rna_pseudo@meta.data$monocle3_pseudotime) # pseudotime
# for coembed
file <- data.frame(barcode=substr(colnames(pbmc.rna_pseudo),1,18), # barcodes
                   pseudotime=pbmc.rna_pseudo@meta.data$monocle3_pseudotime) # pseudotime

## 1) from pseudotime 
file$cell_select <- ifelse(file$pseudotime!=Inf,1,0) #cell_type
table(file$cell_select)
head(file)

## 2) from annotations
cd4 <- subset(pbmc.rna_pseudo, seurat_annotations == "CD4 Naive" | seurat_annotations == "CD4 TEM" | seurat_annotations == "CD4 TCM"|seurat_annotations=="Treg"|seurat_annotations=="HSPC")
cd8 <- subset(pbmc.rna, seurat_annotations == "CD8 Naive" | seurat_annotations == "CD8 TEM_1" | seurat_annotations == "CD8 TEM_2"|seurat_annotations=="MAIT"|seurat_annotations=="gdt"|seurat_annotations=="NK")
bcell <- subset(pbmc.rna, seurat_annotations == "Intermediate B" | seurat_annotations == "Naive B" | seurat_annotations == "Memory B")
#
length(colnames(cd4)) # 3054 # 2698 # 1732
colnames(cd4)
sub <- substr(colnames(cd4),1,18)
length(sub)
file$cell_select <- ifelse(file$barcode %in% sub,1,0)
table(file$cell_select)
head(file)

## Store tenet files
nrow(file[which(file$cell_select==0),])
nrow(file[which(file$cell_select==1),])
write.table(file$barcode,"./tenet+/barcode_cd4.txt",row.names = F,col.names = F,quote = F) 
write.table(file$pseudotime,"./tenet+/trajecory_cd4.txt",row.names = F,col.names = F,quote = F)
write.table(file$cell_select,"./tenet+/cell_select_cd4.txt",row.names = F,col.names = F,quote = F)


