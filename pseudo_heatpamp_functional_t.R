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

pseudo_heatmap_log_none <- function(x,y){
  cell_data <- x
  temp=array()
  for (i in 2:dim(cell_data)[2]){
    lo <- loess(cell_data[,i]~cell_data[,1], span = y)
    kk <- lo$fitted
    temp <- cbind(temp,kk)
    colnames(temp)[i] <- colnames(cell)[i]
  }
  return(temp)
}

log_peuso_data <- function(x){
  cell <- x
  cell_data <- log2(cell[,2:dim(cell)[2]]+1)
  cell_data <- cbind(cell[,1],cell_data)
  return(cell_data)
}

order_peusodata <- function(x,y,z){
  
  data_sort <- x[,y[,1]]
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
library(readr)


matrix <- read_csv("C:/Users/LSCS/Desktop/tenet/matrix.csv")
trajecory_cd8 <- read.table("C:/Users/LSCS/Desktop/tenet/trajecory_cd8.txt", quote="\"", comment.char="")
CD8_tenet <- read.delim("C:/Users/LSCS/Desktop/tenet/tenet_result/CD8/TE_result_matrix.byGRN.fdr0.01.trimIndirect0.0.sif.outdegree.txt", header=FALSE)


#=====================================Data preparation=================================

#Data_order <- pseudo_order(total_matrix,TENET_result,trajecory_data)

CD8_order = order_peusodata(matrix,CD8_tenet,trajecory_cd8)


#====================================Function Version===================================

#object_name = pseudo_heatmap(orderdata, span( 0 > The larger the more accurate))


temp = pseudo_heatmap(CD8_order,0.7)

pheatmap(temp[2:10,], show_colnames=F, cluster_rows = F, cluster_cols = F, main="Pseudotime -->", xlabel  = " a")
pheatmap(temp[1050:1100,], show_colnames=F, cluster_rows = F, cluster_cols = F, main="")

#=====================================verification======================================

cell_data = log_peuso_data(CD8_order)
cell <- CD8_order

par(mfrow=c(2,1))
par(mar=c(1,1,1,1))

p = 7

p1 <- plot(cell[,1], cell[,p], type="l",xlab="",ylab="a")
lines(cell_data[,1],temp[,p],col="red",lwd=2)
p2 <- plot(cell_data[,1], cell_data[,p], type="l",,xlab="",ylab="")
lines(cell_data[,1],temp[,p],col="red",lwd=2)




