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
  return(temp)
}

pseudo_order <- function(x,y,z){

  data_sort <- x[,y[,2]]
  data_sort <- cbind(z,data_sort)
  data_sort <- subset(data_sort,data_sort$V1!="Inf")
  
  b <- as.data.frame(data_sort)
  
  b1 <- order(b$V1)
  
  b_order <- b[b1,]
  
  Data_order <- as.matrix(b_order)
  
  print(dim(Data_order))
  return(Data_order)
}

#=====================================Data preparation=================================

#Data_sort <- pseudo_order(total_matrix,TENET_result,trajecory_data)

CD4_order = pseudo_order(matrix,CD4_tenet,trajecory_cd4)


#====================================Function Version===================================

#object_name = pseudo_heatmap(orderdata, span( 0 > The larger the more accurate))

temp2 = pseudo_heatmap(CD4_order,1)


pheatmap(temp2[,2:10], show_rownames=F, cluster_rows = F, cluster_cols = F, main="")


#===================================================================================

cell <- CD4_order
cell_data <- log2(cell[,2:dim(cell)[2]]+1)
cell_data <- cbind(cell[,1],cell_data)

temp=array()

for (i in 2:dim(cell_data)[2]){
  lo <- loess(cell_data[,i]~cell_data[,1],span = 1)
  kk <- lo$fitted
  temp <- cbind(temp,kk)
  colnames(temp)[i] <- colnames(cell)[i]
}



dim(temp)

temp[,1] <- cell[,1]

pheatmap(temp[,2:dim(cell_data)[2]], show_rownames=F, cluster_rows = F, cluster_cols = F, main="")
pheatmap(temp[,2:30], show_rownames=F, cluster_rows = F, cluster_cols = F, main="")
pheatmap(temp[,80:110], show_rownames=F, cluster_rows = F, cluster_cols = F, main="")

par(mfrow=c(2,1))
par(mar=c(1,1,1,1))

p = 7

p1 <- plot(cell[,1], cell[,p], type="l",xlab="",ylab="a")
lines(cell_data[,1],temp[,p],col="red",lwd=2)
p2 <- plot(cell_data[,1], cell_data[,p], type="l",,xlab="",ylab="")
lines(cell_data[,1],temp[,p],col="red",lwd=2)

