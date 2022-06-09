pseudo_heatmap <- function(x,y){
  cell <- x
  cell_data <- log2(cell[,2:dim(cell)[2]]+1)
  cell_data <- cbind(cell[,1],cell_data)
  temp=array()
  for (i in 2:dim(cell_data)[2]){
    lo <- loess(cell_data[,i]~cell_data[,1], span = y)
    kk <- lo$fitted
    temp <- cbind(temp,kk)
    colnames(temp)[i] <- colnames(cell)[i]
  }
  temp[,1] <- x[,1]
  temp <- t(temp)
  colnames(temp) <- round(temp[1,],2)
  return(temp)
}
order_peusodata2 <- function(x,y,z){
  
  data_sort <- x[,y]
  data_sort <- cbind(z,data_sort)
  data_sort <- subset(data_sort,data_sort$V1!="Inf")
  
  b <- as.data.frame(data_sort)
  
  b1 <- order(b$V1)
  
  b_order <- b[b1,]
  
  Data_order <- as.matrix(b_order)
  
  print(dim(Data_order))
  return(Data_order)
}
order_peusodata_peak <- function(x,y,z){
  
  data_sort <- x[,y]
  colnames(data_sort) <- c(select_gene,paste0(intersection$chr," (",intersection$Gene.Name,")"))
  data_sort <- cbind(z,data_sort)
  data_sort <- subset(data_sort,data_sort$V1!="Inf")
  
  b <- as.data.frame(data_sort)
  
  b1 <- order(b$V1)
  
  b_order <- b[b1,]
  
  Data_order <- as.matrix(b_order)
  
  print(dim(Data_order))
  return(Data_order)
}
library(pheatmap)
merge_matrix <- read.csv("C:/Users/LSCS/Desktop/tenet/merged_matrix.csv")
total_peak_gene <- read.delim("~/total_peak_gene.txt")
total_gene_list <- colnames(merge_matrix)[2:4992]
total_peak_list <- colnames(merge_matrix)[4993:11069]
total_peak_list <- as.data.frame(total_peak_list)
colnames(total_peak_list) <- "test2"
total_peak_gene$test2 <- gsub('-',".",total_peak_gene$chr)
total_peak_list2 <- merge(total_peak_list,total_peak_gene, by="test2")
total_gene_list <- unique(total_gene_list)
total_peak_list <- unique(total_peak_list2$Gene.Name)

length(celltype_list)

celltype_list <- c("CD4","CD8","Bcell")
for (i in 1:length(celltype_list)){
  select_gene <- "KLF2"
  celltype <- celltype_list[i]
  trajectory<- read.table(paste0("C:/Users/LSCS/Desktop/tenet/tenet_files/trajecory_",celltype,".txt"), quote="\"", comment.char="")
  peak_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/new/",celltype,"_merge_TF/TE_result_matrix_rowTF_colPK.fdr0.01.sif"), header=FALSE)
  gene_list <- read.delim(paste0("C:/Users/LSCS/Desktop/tenet/new/",celltype,"_merge_TF/TE_result_matrix_rowTF_colGN.fdr0.01.sif"), header=FALSE)
  
  
  
  peak_target <- subset(peak_list,peak_list$V1 == select_gene)
  gene_target <- subset(gene_list,gene_list$V1 == select_gene)
  total_peak_gene$V3 <- total_peak_gene$chr 
  a <- merge(peak_target,total_peak_gene, by = "V3")
  peak_target_ann <- a[,c(2:5)]
  peak_target_ann$V3 <- peak_target_ann$Gene.Name
  intersection <- merge(gene_target,peak_target_ann, by = "V3")
  
  #==========================gene==============
  
  name <- paste0(celltype,"_",select_gene,"_related_gene_by_TENET pseudotime_expression")
  sif_gene <- subset(gene_list,gene_list$V1==select_gene)
  order_sif <- order(sif_gene$V2,decreasing = T)
  sort_sif <- sif_gene[order_sif,]
  head(sort_sif)
  top_20_sif <- unique(intersection$Gene.Name)
  top_20_sif <- c(select_gene,top_20_sif)
  
  temp2 <- order_peusodata2(merge_matrix,top_20_sif,trajectory)
  temp2 <- pseudo_heatmap(temp2,0.7)
  dim(temp2)
  pheatmap(temp2[2:dim(temp2)[1],], show_colnames=F, cluster_rows = F, cluster_cols = F,main=name)
  
  #==========================peak==============
  name <- paste0(celltype,"_",select_gene,"_related_peak_by_TENET pseudotime_expression")
  sif_gene <- subset(peak_list,peak_list$V1==select_gene)
  order_sif <- order(sif_gene$V2,decreasing = T)
  sort_sif <- sif_gene[order_sif,]
  head(sort_sif)
  top_20_sif <- unique(gsub("-",".",intersection$chr))
  top_20_sif <- c(select_gene,top_20_sif)
  
  temp2 <- order_peusodata_peak(merge_matrix,top_20_sif,trajectory)
  temp2 <- pseudo_heatmap(temp2,0.8)
  dim(temp2)
  pheatmap(temp2[2:dim(temp2)[1],], show_colnames=F, cluster_rows = F, cluster_cols = F,main=name)
}
