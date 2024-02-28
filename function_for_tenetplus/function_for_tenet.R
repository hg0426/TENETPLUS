print("functions_for_TENET_attached_ver.2024.02.21")
library(ggplot2)
library(pheatmap)
library(stringr)
library(igraph)

TE_barplot <- function(gene_list=as.data.frame("non","non"),peak_list=as.data.frame("non","non"),
                       PAG_list=as.data.frame("non","non"),top_gene_num=20,highlight_gene="",ven = F,fill_col = "grey19",
                       highlight_col = "brown4",save_files = F,font_size=10,
                       name_RNA="TE_result_RNA_count_outdegree",
                       name_ATAC="TE_result_ATAC_count_outdegree",
                       name_trim="TE_result_trim_count_outdegree",width = 18, height = 10) {
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
      theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(face="bold"),text = element_text(size = font_size))
    
    
    p + coord_flip()
    if (save_files == T) {
      ggsave(paste0(name_RNA,".png"),width = width, height = width, units= "cm", dpi = 300)
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
      theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(face="bold"),text = element_text(size = font_size))
    
    
    p + coord_flip()
    if (save_files == T) {
      ggsave(paste0(name_ATAC,".png"),width = width, height = width, units= "cm", dpi = 300)
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
      theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(face="bold"),text = element_text(size = font_size))
    
    
    p + coord_flip()
    if (save_files == T) {
      ggsave(paste0(name_trim,".png"),width = width, height = width, units= "cm", dpi = 300)
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

trimming_PAG <- function(total_peak_gene,peak_list,gene_list,name = "trimmed_by_PAG_TE_result_matrix.sif",anno=F,distance=F){
  total_peak_gene$target <- total_peak_gene$chr
  colnames(peak_list) <- c("source","score","target")
  colnames(gene_list) <- c("source","score","target")
  PAG_list <- merge(peak_list,total_peak_gene, by = "target")
  gene_list$Gene.Name <- gene_list$target
  trimmed <- merge(PAG_list,gene_list, by = c("source","Gene.Name"))
  trimmed_sif <- cbind(trimmed$source,trimmed$score.x,trimmed$Gene.Name,trimmed$chr)
  if (anno==T) {
    trimmed_sif <- cbind(trimmed_sif,trimmed$anno)
  }
  if (distance==T) {
    trimmed_sif <- cbind(trimmed_sif,trimmed$Distance.to.TSS)
  }
  write.table(trimmed_sif,name,row.names = F,col.names = F,quote = F,sep = "\t")
  output <- as.data.frame(trimmed_sif)
  colnames_list <- c("source","score","target_gene","target_peak")
  if (anno==T) {
    colnames_list <- c(colnames_list,"anno")
  }
  if (distance==T) {
    colnames_list <- c(colnames_list,"distance")
  }
  colnames(output) <- colnames_list
  return(output)
}

pseudotime_heatmap<- function(matrix,selected_gene,gene_list,peak_list=F,ues_cell_select=F,cell_select=F,total_peak_gene,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                              order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",
                              max_min = T,out_result =F,target_average = F, target,pseudo_decrease=T,save_file=T){
  #====================subset & order by pseudotime =================
  
  sub_matrix <- matrix[,selected_gene]
  sub_matrix <- cbind(pseudotime,sub_matrix)
  sub_matrix <- subset(sub_matrix,sub_matrix[,1]!="Inf" & sub_matrix[,1]!="NA")
  if(ues_cell_select!=F){
    cell_select <- as.data.frame(cell_select)
    sub_matrix <- sub_matrix[which(sample_cell_select==1),]
  }
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
      if (peak_list[1,1] != F) {
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
        
      } else{
      target_list <- subset(gene_list,gene_list[,1] == gsub("[.]","-",selected_gene[j]))[,3]
      }
      target_list <- gsub("-",".",target_list)
      temp_matrix <- as.data.frame(matrix[,target_list])
      if (length(target_list)==1) {
        colnames(temp_matrix) <- target_list
        
      }
      temp_matrix <- cbind(pseudotime,temp_matrix)
      temp_matrix <- subset(temp_matrix,temp_matrix[,1]!="Inf" & temp_matrix[,1]!="NA")
      
      temp_matrix <- as.data.frame(temp_matrix)
      
      temp_matrix_order <- order(temp_matrix[,1])
      
      temp_matrix_sorted <- temp_matrix[temp_matrix_order,]
      
      temp_matrix_sorted <- as.matrix(temp_matrix_sorted)
      
      if (length(target_list)!=1) {
        temp_average <- apply(temp_matrix_sorted[,-1],1,mean)
      }else{
        temp_average <- temp_matrix_sorted[,2]
      }
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
  temp <- cbind(sub_matrix_sorted[,1], apply(sub_matrix_sorted[,-1], 2, function(j) {
    k <- sub_matrix_sorted[,1]
    lo <- loess(j~k, span = span)
    xl <- seq(min(k),max(k), (max(k) - min(k))/(length(k)-1))
    predict(lo,xl)
  }))
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
  if (save_file==T) {
    png(paste0(filename,".png"),width=6000,height=4000,res=500)
    print(s1)
    dev.off()
  }

  if (out_result==T){
    return(final_matrix)
  }
}


pseudotime_heatmap_for_confirm<- function(matrix,selected_gene,gene_list,total_peak_gene,peak_list,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                                          order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",
                                          max_min = T,target_average = F, target, pseudo_decrease=T){
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

data_count_for_cytoscape <- function(GRN_list,celltype_group=GRN_list){
  print(celltype_group)
  merge_uniq_geme <- c()
  for (i in 1:length(GRN_list)) {
    gene_list <- read.delim(GRN_list[i], header=FALSE)
    uniq_gene <- unique(gene_list[,1])
    merge_uniq_geme <- unique(c(uniq_gene,merge_uniq_geme))
  }
  merge_uniq_geme <- as.data.frame(merge_uniq_geme)
  dim(merge_uniq_geme)
  
  for (q in 1:length(GRN_list)) {
    gene_list <- read.delim(GRN_list[q], header=FALSE)
    for (j in 1:dim(merge_uniq_geme)[1]) {
      celltype <- celltype_list[q]
      merge_uniq_geme[j,q+1] <- length(which(merge_uniq_geme[j,1]==gene_list$V1))
    }
    colnames(merge_uniq_geme)[q+1] <- celltype_group[q]
  }
  return(merge_uniq_geme)
}


split_mtx_AB <- function(TE_matrix,species){
  system(paste0("sh /home/Data_Drive_8TB_2/TENET_workshop/TENETPLUS/function_for_tenetplus/script/getMatrix_rowTF_colGN_AB-matrix.sh ",TE_matrix," ",species))
}

split_mtx_C <- function(TE_matrix,species){
  system(paste0("sh /home/Data_Drive_8TB_2/TENET_workshop/TENETPLUS/function_for_tenetplus/script/getMatrix_rowTF_colPK_C-matrix.sh ",TE_matrix," ",species))
}

convertNumberToUnit <- function(x) {
  if(x >= 1000 & x < 1000000) {
    return(paste(x / 1000, "K"))
  } else if(x >= 1000000) {
    return(paste(x / 1000000, "M"))
  } else {
    return(as.character(x))
  }
}

function_PeakSource_Distance <- function(Tenet_result_dir=Tenet_result_dir,
                                         gene_anno=gene_anno,
                                         trim_indirect=T,
                                         trim_distance=1000000,
                                         save=F
){
  
  # Input files -------------------------------------------------------------
  print(paste0("##### Input files"))
  print(paste0("- Tenet_result_dir: ",Tenet_result_dir))
  rowPeak_colGN_file=paste0(Tenet_result_dir,"/TE_result_matrix_rowPeak_colGN.fdr0.01.sif")
  rowTF_colPK_file=paste0(Tenet_result_dir,"/TE_result_matrix_rowTF_colPK.fdr0.01.sif")
  if (trim_indirect==T){
    # Trim indirect 여러개 일수도
    rowTF_colGN_file=paste0(Tenet_result_dir,"/TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect-0.01.sif")
    
  } else if(trim_indirect==F){
    rowTF_colGN_file=paste0(Tenet_result_dir,"/TE_result_matrix_rowTF_colGN.fdr0.01.sif")
  }
  if (file.exists(rowPeak_colGN_file)){
    print(paste0("- rowPeak_colGN_file: ",basename(rowPeak_colGN_file)))
  } else{
    print(paste0("- rowPeak_colGN_file not in directory"))
  }
  if (file.exists(rowTF_colPK_file)){
    print(paste0("- rowTF_colPK_file: ",basename(rowTF_colPK_file)))
  } else{
    print(paste0("- rowTF_colPK_file not in directory"))
  }
  if (file.exists(rowTF_colGN_file)){
    print(paste0("- rowTF_colGN_file: ",basename(rowTF_colGN_file)))
  } else{
    print(paste0("- rowTF_colGN_file not in directory"))
  }
  cat('\n')
  
  # Calculate PeakSource Distance -------------------------------------------------------------
  print(paste0("##### Calculate PeakSource Distance"))
  rowPeak_colGN <- read.delim(paste0(rowPeak_colGN_file),header = F) ; colnames(rowPeak_colGN) <- c("source","score","target")
  dim(rowPeak_colGN) ; head(rowPeak_colGN)
  
  rowPeak_colGN$peak_chr <- str_split(rowPeak_colGN$source,'-',simplify = T)[,1]
  rowPeak_colGN$peak_start <- as.numeric(str_split(rowPeak_colGN$source,'-',simplify = T)[,2])
  rowPeak_colGN$peak_end <- as.numeric(str_split(rowPeak_colGN$source,'-',simplify = T)[,3])
  dim(rowPeak_colGN) ; head(rowPeak_colGN)
  
  ### Process gene_anno 
  gene_anno$gene_start <- as.numeric(gene_anno$gene_start);gene_anno$gene_end <- as.numeric(gene_anno$gene_end)
  gene_dupli <- unique(gene_anno$gene[duplicated(gene_anno$gene)]) ; length(gene_dupli)
  dim(gene_anno);head(gene_anno)
  gene_anno <- gene_anno[!gene_anno$gene %in% gene_dupli,]
  dim(gene_anno);head(gene_anno);length(unique(gene_anno$gene))
  table(gene_anno$gene_chr,useNA = "always")
  
  length(unique(rowPeak_colGN$target)) 
  length(unique(gene_anno$gene))
  length(intersect(rowPeak_colGN$target,gene_anno$gene)) 
  intersect_gene <- length(intersect(unique(rowPeak_colGN$target),gene_anno$gene))
  rowPeak_colGN_gene<- length(unique(rowPeak_colGN$target))
  print(paste0("- annotated gene/total gene: ",intersect_gene,"/",rowPeak_colGN_gene))
  if (rowPeak_colGN_gene>intersect_gene){
    print(paste0("Genes not in Bed File: ",unique(rowPeak_colGN$target)[!unique(rowPeak_colGN$target) %in% gene_anno$gene]))
  }
  rowPeak_colGN <- merge(rowPeak_colGN,gene_anno,all.x=T,by.x="target",by.y="gene")
  rowPeak_colGN$chr_info <- ifelse(rowPeak_colGN$peak_chr == rowPeak_colGN$gene_chr , "same", "diff")
  dim(rowPeak_colGN) ; head(rowPeak_colGN)
  table(rowPeak_colGN$chr_info,useNA = "always")
  table(rowPeak_colGN$peak_chr,useNA = "always")
  table(rowPeak_colGN$gene_chr,useNA = "always")
  
  ### Calculate Distance
  calc_dist <- function(a1,a2,b1,b2){
    if (a2 < b1) {
      distance <- b1 - a2
    } else if (b2 < a1) {
      distance <- a1 - b2
    } else {
      distance <- 0  
    }
    return(distance)
  }
  rowPeak_colGN$distance <- apply(rowPeak_colGN,MARGIN=1,function(x) calc_dist(as.numeric(x["peak_start"]),as.numeric(x["peak_end"]),as.numeric(x["gene_start"]),as.numeric(x["gene_end"])))
  head(rowPeak_colGN);dim(rowPeak_colGN)
  
  rowPeak_colGN_dist_file <- paste0(Tenet_result_dir,"/TE_result_matrix_rowPeak_colGN.fdr0.01_AddDist.sif")
  write.table(rowPeak_colGN,rowPeak_colGN_dist_file,row.names = F,quote = F,col.names = T)
  cat('\n')
  
  # Trimming -------------------------------------------------------------
  print(paste0("##### Trimming"))
  gene_list <- read.delim(rowTF_colGN_file, header=FALSE);colnames(gene_list) <- c("TF","TF_GN_TE","gene")
  peak_list <- read.delim(rowTF_colPK_file, header=FALSE);colnames(peak_list) <- c("TF","TF_PK_TE","peak")
  PK_gene_list <- read.table(rowPeak_colGN_dist_file,header = TRUE)
  PK_gene_list=dplyr::rename(PK_gene_list,"PK_GN_TE"="score")
  PK_gene_list=dplyr::rename(PK_gene_list,"peak"="source")
  PK_gene_list=dplyr::rename(PK_gene_list,"gene"="target")
  PK_gene_list$distance_category <- PK_gene_list$chr_info
  PK_gene_list$distance_category <- ifelse(PK_gene_list$distance_category=="same",
                                           ifelse(PK_gene_list$distance<=1000,"promoter",
                                                  ifelse(PK_gene_list$distance<=1000000,"enhancer","same")),PK_gene_list$distance_category)
  table(PK_gene_list$distance_category,useNA = "always")
  
  ## Step 1: PK -> GN
  print(paste0("- distance: ",trim_distance))
  if (trim_distance != 0){
    PK_gene_list_enhan=subset(PK_gene_list,distance < trim_distance)
  } else if (trim_distance == 0) {
    PK_gene_list_enhan=PK_gene_list
  }
  dim(PK_gene_list_enhan) ;head(PK_gene_list_enhan,3)
  
  ## Step 2: TF -> PK
  PK_gene_list_enhan_peak <- merge(PK_gene_list_enhan,peak_list,by="peak")
  dim(PK_gene_list_enhan_peak) ; head(PK_gene_list_enhan_peak)
  
  ## Step 3: TF -> GN
  trim_result <- merge(PK_gene_list_enhan_peak,gene_list,by=c("TF","gene"))
  trim_result <- trim_result[,c("TF","gene","peak","TF_GN_TE","TF_PK_TE","PK_GN_TE","distance")]
  dim(trim_result);head(trim_result)
  length(paste0(trim_result$TF,"_",trim_result$gene,"_",trim_result$peak))
  length(unique(paste0(trim_result$TF,"_",trim_result$gene,"_",trim_result$peak)))
  print(paste("- trim_result's nrow: ",nrow(trim_result)))
  
  outdegree_trim_new <- as.data.frame(table(trim_result$TF)) ; outdegree_trim_new <- outdegree_trim_new[order(outdegree_trim_new$Freq,decreasing = T),] ; names(outdegree_trim_new) <- c("TF","outdegree")
  outdegree_plot <- ggplot(head(outdegree_trim_new,20))+
    geom_col(mapping=aes(x=reorder(TF,outdegree),y=outdegree))+
    xlab("") + ylab("") +
    theme_classic()+
    ggtitle(paste0("Trimm Outdegree"))  +
    theme(plot.title = element_text(hjust = 0.5,size=20),axis.text = element_text(face="bold"))+
    coord_flip()
  if (save == T){
    trim_distance_temp <- convertNumberToUnit(trim_distance)
    if (trim_indirect == T) {
      result_name=paste0(Tenet_result_dir,"/trimmed_by_PeaksourceDistance",trim_distance_temp,"_Indirect-0.01_TE_result_matrix.sif")
    } else if (trim_indirect == F) {
       result_name=paste0(Tenet_result_dir,"/trimmed_by_PeaksourceDistance",trim_distance_temp,"_TE_result_matrix.sif")
      }
    write.table(trim_result,result_name,quote = F,row.names = F,col.names = T)
    plot_name=gsub("_TE_result_matrix.sif",".pdf",basename(result_name))
    plot_name=gsub("trimmed_by_","Outdegree_",plot_name)
    pdf(paste0(Tenet_result_dir,"/",plot_name))
    print(outdegree_plot)
    dev.off()
    }
  cat('\n')
  return(trim_result)
                                  
}


pseudotime_heatmap2 <- function(matrix,selected_gene,gene_list,peak_list=F,ues_cell_select=F,cell_select=F,total_peak_gene,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                                order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",out_dataset=NULL,
                                max_min = T,out_result =F,target_average = F, target,pseudo_decrease=T,save_file=T,out_plot=T,fontsize=15){
  #====================subset & order by pseudotime =================
  sub_matrix <- matrix[,selected_gene]
  sub_matrix <- cbind(pseudotime,sub_matrix)
  sub_matrix <- subset(sub_matrix,sub_matrix[,1]!="Inf" & sub_matrix[,1]!="NA")
  if(ues_cell_select!=F){
    cell_select <- as.data.frame(cell_select)
    sub_matrix <- sub_matrix[which(sample_cell_select==1),]
  }
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
      if(target == "peak"){
        temp_peak_list <- subset(peak_list,peak_list[,1]==selected_gene[j])
        target_list <- gsub("-",".",unique(temp_peak_list[,3]))
      }
      if(target == "gene"){
        temp_gene_list <- subset(gene_list,gene_list[,1]==selected_gene[j])
        target_list <- gsub("-",".",unique(temp_gene_list[,3]))
      }
      target_list <- gsub("-",".",target_list)
      temp_matrix <- as.data.frame(matrix[,target_list])
      if (length(target_list)==1) {
        colnames(temp_matrix) <- target_list
        
      }
      temp_matrix <- cbind(pseudotime,temp_matrix)
      temp_matrix <- subset(temp_matrix,temp_matrix[,1]!="Inf" & temp_matrix[,1]!="NA")
      
      temp_matrix <- as.data.frame(temp_matrix)
      
      temp_matrix_order <- order(temp_matrix[,1])
      
      temp_matrix_sorted <- temp_matrix[temp_matrix_order,]
      
      temp_matrix_sorted <- as.matrix(temp_matrix_sorted)
      
      if (length(target_list)!=1) {
        temp_average <- apply(temp_matrix_sorted[,-1],1,mean)
      }else{
        temp_average <- temp_matrix_sorted[,2]
      }
      temp_average <- as.data.frame(temp_average)
      colnames(temp_average) <- paste0(selected_gene[j])
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
  temp <- cbind(sub_matrix_sorted[,1], apply(sub_matrix_sorted[,-1], 2, function(j) {
    k <- sub_matrix_sorted[,1]
    lo <- loess(j~k, span = span)
    xl <- seq(min(k),max(k), (max(k) - min(k))/(length(k)-1))
    predict(lo,xl)
  }))
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
                   breaks = seq(p_min,p_max,length.out=100),show_colnames=F, cluster_rows = F, cluster_cols = F,legend = p_legend,main=filename,fontsize=fontsize)
  }else{
    s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1],pheat_start:dim(final_matrix)[2]],
                   show_colnames=F, cluster_rows = F, cluster_cols = F,legend = p_legend,main=filename,fontsize=fontsize)
    
  }
  if (save_file==T) {
    png(paste0(filename,".png"),width=6000,height=4000,res=500)
    print(s1)
    dev.off()
  }
  if (is.null(out_dataset)!=T){
    assign(paste0(out_dataset),final_matrix, envir = .GlobalEnv)
    #print("out_dataset")
  }
  
  if (out_plot==T){
    return(s1)
  }

}

mytriangle <- function(coords, v = 
                         NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(
    x = coords[, 1], y = coords[, 2], bg = vertex.color,
    stars = cbind(vertex.size, vertex.size, vertex.size),
    add = TRUE, inches = FALSE
  )
}
add_shape("triangle",
          clip = shapes("circle")$clip,
          plot = mytriangle
)

make_GRN_graph <- function(input_list, nodes="nodes",edges="edges",node_sizes="node_sizes",
                           node_label_sizes = "node_label_sizes",node_colors = "node_colors",
                           node_shapes = "node_shapes", min_node_size = 5, min_size_cutoff = 5,
                           min_size = 5,min_node_label_size = 1,node_size_adjust=1){
  temp_nodes <- unique(c(input_list$TF,input_list$gene,input_list$peak))
  assign(paste0(nodes),temp_nodes, envir = .GlobalEnv)
  
  node_shapes_temp <- rep("circle",length(temp_nodes))
  for(i in 1:length(temp_nodes)){
    if (temp_nodes[i] %in% input_list$gene){
      node_shapes_temp[i] <- "triangle"
      if(temp_nodes[i] %in% input_list$TF){
        node_shapes_temp[i] <- "circle"
      }
    }else if (temp_nodes[i] %in% input_list$TF){
      node_shapes_temp[i] <- "circle"
    }else if (temp_nodes[i] %in% input_list$peak){
      node_shapes_temp[i] <- "square"
    }
    #print(paste0(nodes[i],"_",node_shapes_temp[i]))
  }

  assign(paste0(node_shapes),node_shapes_temp, envir = .GlobalEnv)  
  
  temp_edges <- cbind(c(input_list$TF),c(input_list$gene))
  temp_edges2 <- cbind(c(input_list$peak),c(input_list$gene))
  temp_edges3 <- cbind(c(input_list$TF),c(input_list$peak))
  total_edges <- rbind(temp_edges,temp_edges2,temp_edges3)
  total_edges <- as.matrix(total_edges)
  assign(edges,total_edges, envir = .GlobalEnv)  
  
  
  input_count <- count_degree(input_list)
  size_dict <- as.list(input_count$Freq)
  names(size_dict) <- input_count$Var1
  
  # 각 노드의 크기 설정
  node_sizes_temp <- sapply(temp_nodes, function(temp_nodes) {
    if (is.null(size_dict[[temp_nodes]])) {
      # NULL 값인 경우 기본값으로 설정
      return(min_node_size)  # 또는 다른 기본값 지정 가능
    }else if(size_dict[[temp_nodes]] <= min_size_cutoff){
      return(min_node_size)
    }else {
      # NULL 값이 아닌 경우 해당 크기 값 반환
      return(size_dict[[temp_nodes]]*node_size_adjust)
    }
  })
  
  assign(paste0(node_sizes),node_sizes_temp, envir = .GlobalEnv)
  
  node_label_sizes_temp <- sapply(temp_nodes, function(temp_nodes) {
    if (is.null(size_dict[[temp_nodes]])) {
      # NULL 값인 경우 기본값으로 설정
      return(min_node_label_size)  # 또는 다른 기본값 지정 가능
    }else if(size_dict[[temp_nodes]] <= min_size_cutoff){
      return(min_node_label_size)
    }else {
      # NULL 값이 아닌 경우 해당 크기 값 반환
      #return(2.5)
      return(size_dict[[temp_nodes]]/5)
    }
  })
  assign(paste0(node_label_sizes),node_label_sizes_temp, envir = .GlobalEnv)
  
  node_colors_temp <- rep("orange",length(temp_nodes))
  for(i in 1:length(temp_nodes)){
    if (temp_nodes[i] %in% input_list$gene){
      node_colors_temp[i] <- "green"
      if(temp_nodes[i] %in% input_list$TF){
        node_colors_temp[i] <- "orange"
      }
    }else if (temp_nodes[i] %in% input_list$TF){
      node_colors_temp[i] <- "orange"
    }else if (temp_nodes[i] %in% input_list$peak){
      node_colors_temp[i] <- "skyblue"
    }
  }
  
  assign(paste0(node_colors),node_colors_temp, envir = .GlobalEnv)
  
  graph <- graph_from_data_frame(total_edges, directed = TRUE, vertices = as.character(temp_nodes))
  return(graph)
}
