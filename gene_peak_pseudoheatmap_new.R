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
library(pheatmap)
TF_gene <- read.delim("C:/Users/LSCS/Desktop/tenet/CD8_merge_TF/TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect0.0.sif", header=FALSE)
TF_peak <- read.delim("C:/Users/LSCS/Desktop/tenet/CD8_merge_TF/TE_result_matrix_rowTF_colPK.fdr0.01.trimIndirect0.0.sif", header=FALSE)
TF_peak$V3 <- gsub("-",".",TF_peak$V3)
merge_matrix <- read.csv("C:/Users/LSCS/Desktop/tenet/merged_matrix.csv")
trajectory<- read.table("C:/Users/LSCS/Desktop/tenet/tenet_files/trajecory_cd8.txt", quote="\"", comment.char="")
#==========================gene==============
select_gene <- "LEF1"
name <- paste0(select_gene,"_related_top20_gene_by_TENET pseudotime_expression")
sif_gene <- subset(TF_gene,TF_gene$V1==select_gene)
order_sif <- order(sif_gene$V2,decreasing = T)
sort_sif <- sif_gene[order_sif,]
head(sort_sif)
top_20_sif <- sort_sif[1:20,3]
top_20_sif <- c(select_gene,top_20_sif)
#===============heatmap===============================

temp2 <- order_peusodata2(merge_matrix,top_20_sif,trajectory)
temp2 <- pseudo_heatmap(temp2,0.7)
dim(temp2)
pheatmap(temp2[2:22,], show_colnames=F, cluster_rows = F, cluster_cols = F,main=name)

#==========================peak==============
select_gene <- "KLF6"
name <- paste0(select_gene,"_related_top20_peak_by_TENET pseudotime_expression")
sif_gene <- subset(TF_peak,TF_peak$V1==select_gene)
order_sif <- order(sif_gene$V2,decreasing = T)
sort_sif <- sif_gene[order_sif,]
head(sort_sif)
top_20_sif <- sort_sif[1:20,3]
top_20_sif
top_20_sif <- c(select_gene,top_20_sif)

#===============heatmap===============================

temp2 <- order_peusodata2(merge_matrix,top_20_sif,trajectory)
temp2 <- pseudo_heatmap(temp2,0.7)
dim(temp2)
pheatmap(temp2[2:22,], show_colnames=F, cluster_rows = F, cluster_cols = F,main=name)

