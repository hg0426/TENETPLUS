pseudotime_heatmap<- function(matrix,selected_gene,gene_list,peak_list,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                              order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",
                              max_min = T,out_result =F,target_average = F, target){
  #====================subset & order by pseudotime =================
  
  sub_matrix <- matrix[,selected_gene]
  sub_matrix <- cbind(pseudotime,sub_matrix)
  sub_matrix <- subset(sub_matrix,sub_matrix[,1]!="Inf" & sub_matrix[,1]!="NA")
  
  sub_matrix <- as.data.frame(sub_matrix)
  
  sub_matrix_order <- order(sub_matrix[,1])
  
  sub_matrix_sorted <- sub_matrix[sub_matrix_order,]
  
  sub_matrix_sorted <- as.matrix(sub_matrix_sorted)
  
  print(dim(sub_matrix_sorted))
  
  #====================z-score =================
  
  if(use_z_score == T){
    for (i in 2:dim(sub_matrix_sorted)[2]){
      #ave_merge_z[,i] <- log2(ave_merge_z[,i]+1)
      sub_matrix_sorted[,i] <-(sub_matrix_sorted[,i]-mean(sub_matrix_sorted[,i]))/sd(sub_matrix_sorted[,i])
    }
  }
  #====================use origin pseudotime or use order pseudotime =================
  if(use_pseudotime_origin == T){
  }else{
    expanded_pseudotime <- data.frame(1:dim(sub_matrix_sorted)[1])
    sub_matrix_sorted[,1] <- expanded_pseudotime[,1]
  }
  
  if (target_average == T){
    average_matrix <- sub_matrix_sorted[,1]
    for (j in 1:length(selected_gene)){
      peak_target <- subset(peak_list,peak_list[,1] == selected_gene[j])
      colnames(peak_target) <- c("source","score","chr")
      gene_target <- subset(gene_list,gene_list[,1] == selected_gene[j])
      colnames(gene_target) <- c("source","score","Gene.Name")
      PAG_target <- merge(peak_target,total_peak_gene, by = "chr")
      intersection <- merge(gene_target,PAG_target, by = "Gene.Name")
      if(target == "peak"){
        target_list <- gsub("-",".",unique(intersection$chr))
      }
      if(target == "gene"){
        target_list <- gsub("-",".",unique(intersection$Gene.Name))
      }
      
      temp_matrix <- matrix[,target_list]
      temp_matrix <- cbind(pseudotime,temp_matrix)
      temp_matrix <- subset(temp_matrix,temp_matrix[,1]!="Inf" & temp_matrix[,1]!="NA")
      
      temp_matrix <- as.data.frame(temp_matrix)
      
      temp_matrix_order <- order(temp_matrix[,1])
      
      temp_matrix_sorted <- temp_matrix[temp_matrix_order,]
      
      temp_matrix_sorted <- as.matrix(temp_matrix_sorted)
      
      temp_average <- apply(temp_matrix_sorted[,-1],1,mean)
      temp_average <- as.data.frame(temp_average)
      colnames(temp_average) <- paste0(selected_gene[j]," (",dim(temp_matrix_sorted)[2]-1,")")
      average_matrix <- cbind(average_matrix,temp_average)
    }
    sub_matrix_sorted <- average_matrix
  }
  #====================normalization =================
  temp=array()
  for (i in 2:dim(sub_matrix_sorted)[2]){
    j <- sub_matrix_sorted[,i]
    k <- sub_matrix_sorted[,1]
    lo <- loess(j~k, span = span)
    xl <- seq(min(k),max(k), (max(k) - min(k))/(length(k)-1))
    out = predict(lo,xl)
    temp <- cbind(temp,out)
    colnames(temp)[i] <- colnames(sub_matrix_sorted)[i]
  }
  temp[,1] <- sub_matrix_sorted[,1]
  colnames(temp)[1] <- "pseudotime"
  temp <- t(temp)
  
  normalized_sub_matrix_sorted <- temp
  #====================pseudo_score =================
  if(order_pseudo_score==T){
    pseudo_score <- data.frame()
    
    for (i in 2:dim(normalized_sub_matrix_sorted)[1]){
      pseudo_score[i,1] <- sum(normalized_sub_matrix_sorted[1,]*normalized_sub_matrix_sorted[i,])
    }
    
    normalized_sub_matrix_sorted_pseudo_ordered <- cbind(pseudo_score,normalized_sub_matrix_sorted)
    rownames(normalized_sub_matrix_sorted_pseudo_ordered) <- rownames(normalized_sub_matrix_sorted)
    temp <- normalized_sub_matrix_sorted_pseudo_ordered[-1,]
    temp <- temp[order(temp[,1]),]
    normalized_sub_matrix_sorted_pseudo_ordered <- rbind(normalized_sub_matrix_sorted_pseudo_ordered[1,],temp)
    final_matrix <- normalized_sub_matrix_sorted_pseudo_ordered
    pheat_start <-2
  }else{
    final_matrix <- normalized_sub_matrix_sorted
    pheat_start <-1
  }
  
  #====================make pheatmap =================
  if(max_min==T){
    
    s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1],pheat_start:dim(final_matrix)[2]],
                   breaks = seq(p_min,p_max,length.out=100),show_colnames=F, cluster_rows = F, cluster_cols = F,legend = p_legend)
  }else{
    s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1],pheat_start:dim(final_matrix)[2]],
                   show_colnames=F, cluster_rows = F, cluster_cols = F,legend = p_legend)
    
  }
  
  png(paste0(filename,".png"),width=6000,height=4000,res=500)
  print(s1)
  dev.off()
  if (out_result==T){
    return(final_matrix)
  }
}



sample_gene_list <- c("ELF1","TCF7","LRRFIP1")
  
pseudotime_heatmap(matrix = matrix,selected_gene = sample_gene_list,pseudotime = sample_pseudotime,
                   span = 0.6,use_pseudotime_origin = F,use_z_score = T,
                   order_pseudo_score = T,max_min = F, p_max = 1, p_min = -0.8, out_result = F, filename = "test",
                   target_average = T,peak_list = peak_list,gene_list = gene_list, target = "gene")

confirm <- pseudotime_heatmap_for_confirm(matrix = matrix,selected_gene = sample_gene_list,
                                          pseudotime = sample_pseudotime,span = 0.6,use_pseudotime_origin = F,
                                          use_z_score = T,order_pseudo_score = F,target_average = T,peak_list = peak_list,gene_list = gene_list, target = "gene")

par(mfrow=c(3,1))

# ================for confirm Check by changing the numbers in the loop===============

for (i in 2:4){
  x <- as.numeric(confirm[1,])
  y <- as.numeric(confirm[i,])
  lo <- loess(y~x, span= 0.8)
  xl <- seq(min(x),max(x), (max(x) - min(x))/(length(y)-1))
  out = predict(lo,xl)
  plot(x,y,type="p", main = rownames(confirm)[i],ylim = c(-0.8,0.8),xlab = "",ylab = "")
  lines(xl, out, col='red', lwd=2)
}
