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

temp1 = sudo(CD4_order,0.2)

pheatmap(temp1[,2:10], show_rownames=F, cluster_rows = F, cluster_cols = F, main="")
