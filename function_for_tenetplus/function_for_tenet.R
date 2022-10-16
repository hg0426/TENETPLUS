print("functions_for_TENET_2022.10.17")

TE_barplot <- function(gene_list=as.data.frame("non","non"),peak_list=as.data.frame("non","non"),
                       PAG_list=as.data.frame("non","non"),top_gene_num=20,highlight_gene="",ven = F,fill_col = "grey19",
                       highlight_col = "brown4",save_files = F,
                       name_RNA="TE_result_RNA_count_outdegree",
                       name_ATAC="TE_result_ATAC_count_outdegree",
                       name_trim="TE_result_trim_count_outdegree") {
  #======================RNA_Gene=======================================
  if (gene_list[1,1] != "non") {
    att <- as.data.frame(table(gene_list[,1]))[order(as.data.frame(table(gene_list[,1]),)[,2], decreasing = T),]
    top_genes <- att[1:top_gene_num,]
    top_RNA <- top_genes
    cols <- rep(fill_col,times=top_gene_num)
    if(length(highlight_gene)>0){
      for (z in 1:length(top_genes[,1])){
        if (top_genes[z,1] %in% highlight_gene){
          cols[z] <- highlight_col
        }
      }
      
    }
    p <- ggplot(data=top_genes, aes(x=reorder(Var1, Freq), y=Freq)) + geom_bar(stat="identity",fill=cols) +
      theme(plot.margin=unit(c(0.5,5,0.5,5), 'cm'))+
      xlab("") + ylab("") + ggtitle(name_RNA)  +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    p + coord_flip()
    if (save_files == T) {
      ggsave(paste0(name_RNA,".png"))
    }
  
    return(p + coord_flip())  
  }
  
  #===================ATAC_Peak==================================
  if (peak_list[1,1] != "non") {
    att <- as.data.frame(table(peak_list[,1]))[order(as.data.frame(table(peak_list[,1]),)[,2], decreasing = T),]
    top_genes <- att[1:top_gene_num,]
    top_ATAC <- top_genes
    cols <- rep(fill_col,times=top_gene_num)
    if(length(highlight_gene)>0){
      for (z in 1:length(top_genes[,1])){
        if (top_genes[z,1] %in% highlight_gene){
          cols[z] <- highlight_col
        }
      }
      
    }
    p <- ggplot(data=top_genes, aes(x=reorder(Var1, Freq), y=Freq)) + geom_bar(stat="identity",fill=cols) +
      theme(plot.margin=unit(c(0.5,5,0.5,5), 'cm'))+
      xlab("") + ylab("") + ggtitle(name_ATAC)  +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    p + coord_flip()
    if (save_files == T) {
      ggsave(paste0(name_ATAC,".png"))
    }
    return(p + coord_flip())
  }
  #=================trimmed_by_PAG======================
  if (PAG_list[1,1] != "non") {
    att <- as.data.frame(table(PAG_list[,1]))[order(as.data.frame(table(PAG_list[,1]),)[,2], decreasing = T),]
    top_genes <- att[1:top_gene_num,]
    top_Trim <- top_genes
    cols <- rep(fill_col,times=top_gene_num)
    if(length(highlight_gene)>0){
      for (z in 1:length(top_genes[,1])){
        if (top_genes[z,1] %in% highlight_gene){
          cols[z] <- highlight_col
        }
      }
      
    }
    p <- ggplot(data=top_genes, aes(x=reorder(Var1, Freq), y=Freq)) + geom_bar(stat="identity",fill=cols) +
      theme(plot.margin=unit(c(0.5,5,0.5,5), 'cm'))+
      xlab("") + ylab("") + ggtitle(name_trim)  +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    p + coord_flip()
    if (save_files == T) {
      ggsave(paste0(name_trim,".png"))
    }
    
    return(p + coord_flip())
  }
  
  if (ven == T) {
    myCol <- brewer.pal(3, "Pastel2")
    venn.diagram(
      x = list(top_RNA[1:top_gene_num,][1], top_ATAC[1:top_gene_num,][1], top_Trim[1:top_gene_num,][1]),
      category.names = c("ATAC" , "PK " , "GN"),
      filename = paste0(i,"_ven.png"),
      output=TRUE,
      fill = myCol,
      cex = 2,
      fontface = "bold",
      fontfamily = "sans",
      cat.cex = 2,
      cat.fontface = "bold",
      cat.default.pos = "outer",
    )
    
  }
}

trimming_PAG <- function(total_peak_gene,peak_list,gene_list,name = "trimmed_by_PAG_TE_result_matrix.sif"){
  total_peak_gene$target <- total_peak_gene$chr
  colnames(peak_list) <- c("source","score","target")
  colnames(gene_list) <- c("source","score","target")
  PAG_list <- merge(peak_list,total_peak_gene, by = "target")
  gene_list$Gene.Name <- gene_list$target
  trimmed <- merge(PAG_list,gene_list, by = c("source","Gene.Name"))
  trimmed_sif <- cbind(trimmed$source,trimmed$score.x,trimmed$Gene.Name,trimmed$chr)
  write.table(trimmed_sif,name,row.names = F,col.names = F,quote = F,sep = "\t")
  output <- as.data.frame(trimmed_sif)
  colnames(output) <- c("source","score","target_gene","target_peak")
  return(output)
}

pseudotime_heatmap<- function(matrix,selected_gene,gene_list,peak_list=F,total_peak_gene,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                              order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",
                              max_min = T,out_result =F,target_average = F, target,pseudo_decrease=T){
  #====================subset & order by pseudotime =================
  
  sub_matrix <- matrix[,selected_gene]
  sub_matrix <- cbind(pseudotime,sub_matrix)
  sub_matrix <- subset(sub_matrix,sub_matrix[,1]!="Inf" & sub_matrix[,1]!="NA")
  
  sub_matrix <- as.data.frame(sub_matrix)
  
  sub_matrix_order <- order(sub_matrix[,1])
  
  sub_matrix_sorted <- sub_matrix[sub_matrix_order,]
  
  sub_matrix_sorted <- as.matrix(sub_matrix_sorted)
  
  print(dim(sub_matrix_sorted))
  
  #====================use origin pseudotime or use order pseudotime =================
  if(use_pseudotime_origin == T){
  }else{
    expanded_pseudotime <- data.frame(1:dim(sub_matrix_sorted)[1])
    sub_matrix_sorted[,1] <- expanded_pseudotime[,1]
  }
  #=================== target average heatmap =================
  if (target_average == T){
    average_matrix <- sub_matrix_sorted[,1]
    for (j in 1:length(selected_gene)){
      if (peak_list != F) {
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
        
      }
      target_list <- subset(gene_list,gene_list[,1] == selected_gene[j])[,3]
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
  #====================z-score =================
  
  if(use_z_score == T){
    for (i in 2:dim(sub_matrix_sorted)[2]){
      #ave_merge_z[,i] <- log2(ave_merge_z[,i]+1)
      sub_matrix_sorted[,i] <-(sub_matrix_sorted[,i]-mean(sub_matrix_sorted[,i]))/sd(sub_matrix_sorted[,i])
    }
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
    temp <- temp[order(temp[,1],decreasing = pseudo_decrease),]
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

pseudotime_heatmap_for_confirm<- function(matrix,selected_gene,gene_list,total_peak_gene,peak_list,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                                          order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",
                                          max_min = T,target_average = F, target,pseudo_decrease=T){
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
  
  #====================target average heatmap =================
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
  normalized_sub_matrix_sorted <- t(sub_matrix_sorted)
  #====================pseudo_score =================
  if(order_pseudo_score==T){
    pseudo_score <- data.frame()
    
    for (i in 2:dim(normalized_sub_matrix_sorted)[1]){
      pseudo_score[i,1] <- sum(normalized_sub_matrix_sorted[1,]*normalized_sub_matrix_sorted[i,])
    }
    
    normalized_sub_matrix_sorted_pseudo_ordered <- cbind(pseudo_score,normalized_sub_matrix_sorted)
    rownames(normalized_sub_matrix_sorted_pseudo_ordered) <- rownames(normalized_sub_matrix_sorted)
    temp <- normalized_sub_matrix_sorted_pseudo_ordered[-1,]
    temp <- temp[order(temp[,1],decreasing = pseudo_decrease),]
    normalized_sub_matrix_sorted_pseudo_ordered <- rbind(normalized_sub_matrix_sorted_pseudo_ordered[1,],temp)
    final_matrix <- normalized_sub_matrix_sorted_pseudo_ordered
  }else{
    final_matrix <- normalized_sub_matrix_sorted
  }
  return(final_matrix)
}


count_degree <- function(data,decreasing =T,degree = "out"){
  if (degree == "out") {
    att <- as.data.frame(table(data[,1]))[order(as.data.frame(table(data[,1]),)[,2], decreasing = decreasing),]
  }
  if (degree == "in") {
    att <- as.data.frame(table(data[,3]))[order(as.data.frame(table(data[,3]),)[,2], decreasing = decreasing),]
  }
  
  return(att)
}

define_source_target <- function(data,source ="TF",target="gene"){
  source_list <- unique(data[,1])
  target_list <- unique(data[,3])
  target_sorted <- target_list[!(target_list %in% source_list)]
  source_data <- cbind(source_list,source)
  target_data <- cbind(target_sorted,target)
  data <- rbind(source_data,target_data)
  colnames(data) <- c("unique","define")
  return(data)
}
