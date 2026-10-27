# TENETPLUS downstream R utilities.
# Gene coordinates must match the species, assembly, and identifiers of the input.
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


# Compatibility entry points for external legacy matrix-splitting scripts.
# Current Parquet splitting is performed by TENET_Plus_for_py.sh.
.run_legacy_split <- function(TE_matrix, species, script) {
  if (!nzchar(script) || !file.exists(script))
    stop("Supply the legacy split script path; for current Parquet outputs use TENET_Plus_for_py.sh.")
  status <- system2("sh", c(shQuote(normalizePath(script)), shQuote(TE_matrix), shQuote(species)))
  if (status != 0L) stop("Legacy matrix splitting failed with exit status ", status)
  invisible(status)
}
split_mtx_AB <- function(TE_matrix, species, script = Sys.getenv("TENET_LEGACY_SPLIT_AB")) {
  .run_legacy_split(TE_matrix, species, script)
}
split_mtx_C <- function(TE_matrix, species, script = Sys.getenv("TENET_LEGACY_SPLIT_C")) {
  .run_legacy_split(TE_matrix, species, script)
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

function_PeakSource_Distance <- function(Tenet_result_dir = NULL,
                                        gene_anno = NULL,
                                        trim_indirect = TRUE,
                                        trim_distance = 1000000,
                                        save = FALSE,
                                        rowPeak_colGN = NA,
                                        gene_list = NA,
                                        peak_list = NA,
                                        fdr = 0.01,
                                        contig_aliases = NULL) {
  if (is.null(gene_anno)) stop("Supply gene_anno with matching gene coordinates.")
  if (length(trim_distance) != 1L || !is.finite(trim_distance) || trim_distance < 0)
    stop("trim_distance must be a finite non-negative number.")
  if (length(fdr) != 1L || !is.finite(fdr) || fdr <= 0 || fdr > 1)
    stop("fdr must be in (0, 1].")
  alpha <- format(fdr, scientific = FALSE, trim = TRUE)
  root <- if (is.null(Tenet_result_dir)) NULL else normalizePath(Tenet_result_dir, mustWork = TRUE)
  manifest_path <- if (!is.null(root)) file.path(root, "grn_network_outputs.tsv") else ""
  manifest <- if (nzchar(manifest_path) && file.exists(manifest_path))
    read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE,
               colClasses = "character", quote = "", comment.char = "") else NULL
  canonical <- function(x) {
    x <- as.character(x)
    if (!is.null(contig_aliases)) {
      if (is.null(names(contig_aliases)) || anyDuplicated(names(contig_aliases)) ||
          anyNA(contig_aliases) || anyNA(names(contig_aliases)) ||
          any(!nzchar(names(contig_aliases))) || any(!nzchar(contig_aliases)))
        stop("contig_aliases must be a named vector of explicit alias-to-canonical mappings.")
      i <- match(x, names(contig_aliases)); use <- !is.na(i)
      x[use] <- as.character(contig_aliases[i[use]])
    }
    x
  }
  matches_threshold <- function(paths, prefix, suffix) {
    labels <- basename(as.character(paths))
    valid <- !is.na(labels) & startsWith(labels, prefix) & endsWith(labels, suffix)
    result <- rep(FALSE, length(labels))
    idx <- which(valid)
    values <- suppressWarnings(as.numeric(substring(labels[idx], nchar(prefix) + 1,
                                                     nchar(labels[idx]) - nchar(suffix))))
    # The Python Parquet exporter uses six significant digits in filenames.
    result[idx] <- !is.na(values) & values %in% c(fdr, as.numeric(sprintf("%.6g", fdr)))
    result
  }
  network_path <- function(key, legacy_key, trimmed = FALSE) {
    if (is.null(root)) stop("Supply Tenet_result_dir or all three network data frames.")
    suffix <- if (trimmed) ".trimIndirect-0.01.sif" else ".sif"
    if (!is.null(manifest)) {
      required <- c("te_parquet", "sif_before_trim", "trimmed_sif")
      if (!all(required %in% names(manifest))) stop("Invalid grn_network_outputs.tsv schema.")
      hit <- which(basename(manifest$te_parquet) == paste0(key, ".parquet") &
                   matches_threshold(manifest$sif_before_trim, paste0(key, "_fdr"), ".sif"))
      if (length(hit) != 1L) stop("Expected one manifest entry for ", key, " at FDR ", alpha)
      value <- manifest[[if (trimmed) "trimmed_sif" else "sif_before_trim"]][hit]
      if (is.na(value) || !nzchar(value)) stop("Requested network was not generated: ", key, suffix)
      path <- if (grepl("^/", value)) value else file.path(root, value)
    } else {
      entries <- list.files(root, full.names = TRUE)
      dirs <- entries[matches_threshold(entries, paste0("network_", key, "_fdr"), "") & dir.exists(entries)]
      paths <- unlist(lapply(dirs, list.files, full.names = TRUE), use.names = FALSE)
      found <- paths[matches_threshold(paths, paste0(key, "_fdr"), suffix)]
      if (!length(found))
        found <- entries[matches_threshold(entries, paste0("TE_result_matrix_", legacy_key, ".fdr"), suffix)]
      if (length(found) != 1L)
        stop("Expected one network file for ", key, " at FDR ", alpha, "; found ", length(found))
      path <- found[1]
    }
    if (!file.exists(path)) stop("Network file not found: ", path)
    path
  }
  network <- function(value, key, legacy_key, cols, trimmed = FALSE) {
    if (!is.data.frame(value)) {
      path <- if (is.character(value) && length(value) == 1L && !is.na(value)) value else
        network_path(key, legacy_key, trimmed)
      if (!file.exists(path)) stop("Network file not found: ", path)
      if (file.info(path)$size == 0) {
        value <- data.frame(source = character(), score = numeric(), target = character())
      } else {
        value <- read.delim(path, header = FALSE, quote = "", comment.char = "",
                            colClasses = "character", check.names = FALSE)
      }
    }
    if (all(c("Source", "TE", "Target") %in% names(value)))
      value <- value[, c("Source", "TE", "Target"), drop = FALSE]
    if (ncol(value) != 3L) stop("Network must have source, score, target columns: ", key)
    names(value) <- cols
    value[[1]] <- as.character(value[[1]]); value[[3]] <- as.character(value[[3]])
    value[[2]] <- suppressWarnings(as.numeric(as.character(value[[2]])))
    if (anyNA(value) || any(!nzchar(value[[1]])) || any(!nzchar(value[[3]])) ||
        any(!is.finite(value[[2]]))) stop("Invalid network values: ", key)
    value
  }
  pkgn <- network(rowPeak_colGN, "TE_PK_GN", "rowPeak_colGN", c("peak", "PK_GN_TE", "gene"))
  tfgn <- network(gene_list, "TE_TF_GN", "rowTF_colGN", c("TF", "TF_GN_TE", "gene"), trim_indirect)
  tfpk <- network(peak_list, "TE_TF_PK", "rowTF_colPK", c("TF", "TF_PK_TE", "peak"))
  if (is.character(gene_anno) && length(gene_anno) == 1L)
    gene_anno <- read.table(gene_anno, header = TRUE, quote = "", comment.char = "",
                            colClasses = "character", check.names = FALSE)
  fields <- c("gene_chr", "gene_start", "gene_end", "gene")
  if (!is.data.frame(gene_anno) || !all(fields %in% names(gene_anno)))
    stop("gene_anno must contain gene_chr, gene_start, gene_end, gene.")
  anno <- gene_anno[, fields, drop = FALSE]
  anno$gene <- as.character(anno$gene); anno$gene_chr <- canonical(anno$gene_chr)
  for (col in c("gene_start", "gene_end")) anno[[col]] <- suppressWarnings(as.numeric(as.character(anno[[col]])))
  if (anyNA(anno) || any(!nzchar(anno$gene)) || any(!nzchar(anno$gene_chr)) ||
      any(!is.finite(anno$gene_start)) || any(!is.finite(anno$gene_end)) ||
      any(anno$gene_start < 0 | anno$gene_end <= anno$gene_start) ||
      any(anno$gene_start != floor(anno$gene_start) | anno$gene_end != floor(anno$gene_end)))
    stop("Invalid gene intervals; use 0-based, half-open coordinates.")
  duplicates <- unique(anno$gene[duplicated(anno$gene)])
  # Preserve the original R policy: exclude every row for duplicate gene names.
  anno <- anno[!anno$gene %in% duplicates, , drop = FALSE]
  peaks <- unique(c(pkgn$peak, tfpk$peak))
  peak_info <- data.frame(peak = peaks, peak_chr = character(length(peaks)),
                          peak_start = numeric(length(peaks)), peak_end = numeric(length(peaks)))
  feature_path <- if (!is.null(root)) file.path(root, "feature_manifest.json") else ""
  if (nzchar(feature_path) && file.exists(feature_path)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite to read feature_manifest.json.")
    features <- jsonlite::fromJSON(feature_path, simplifyVector = FALSE)$features
    ids <- vapply(features, function(x) x$feature_id, character(1))
    if (anyDuplicated(ids)) stop("Duplicate IDs in feature_manifest.json.")
    positions <- match(peaks, ids)
    if (anyNA(positions)) stop("Peak missing from feature_manifest.json.")
    for (i in seq_along(peaks)) {
      f <- features[[positions[i]]]
      if (!identical(f$feature_type, "peak")) stop("Non-peak endpoint in peak network: ", peaks[i])
      peak_info$peak_chr[i] <- f$contig; peak_info$peak_start[i] <- f$start; peak_info$peak_end[i] <- f$end
    }
    rna_ids <- unique(c(pkgn$gene, tfgn$TF, tfgn$gene, tfpk$TF))
    pos <- match(rna_ids, ids)
    if (anyNA(pos) || any(vapply(features[pos], function(x) !identical(x$feature_type, "RNA"), logical(1))))
      stop("Unknown or non-RNA endpoint in RNA network.")
  } else {
    parsed <- regmatches(peaks, regexec("^(.+)[:-]([0-9]+)-([0-9]+)$", peaks))
    if (any(lengths(parsed) != 4L)) stop("Peak ID must be contig-start-end or contig:start-end.")
    for (i in seq_along(peaks)) {
      peak_info$peak_chr[i] <- parsed[[i]][2]
      peak_info$peak_start[i] <- as.numeric(parsed[[i]][3]); peak_info$peak_end[i] <- as.numeric(parsed[[i]][4])
    }
  }
  peak_info$peak_chr <- canonical(peak_info$peak_chr)
  if (anyNA(peak_info) || any(!nzchar(peak_info$peak_chr)) ||
      any(!is.finite(peak_info$peak_start)) || any(!is.finite(peak_info$peak_end)) ||
      any(peak_info$peak_start < 0 | peak_info$peak_start >= peak_info$peak_end) ||
      any(peak_info$peak_start != floor(peak_info$peak_start) | peak_info$peak_end != floor(peak_info$peak_end)))
    stop("Invalid peak coordinates.")
  annotated <- merge(pkgn, peak_info, by = "peak", all.x = TRUE, sort = FALSE)
  annotated <- merge(annotated, anno, by = "gene", all.x = TRUE, sort = FALSE)
  annotated$status <- ifelse(annotated$gene %in% duplicates, "duplicate_gene",
                            ifelse(is.na(annotated$gene_chr), "unmapped_gene",
                                   ifelse(annotated$peak_chr != annotated$gene_chr, "different_contig", "same_contig")))
  annotated$distance <- rep(NA_real_, nrow(annotated))
  same <- annotated$status == "same_contig"
  annotated$distance[same] <- pmax(annotated$gene_start[same] - annotated$peak_end[same],
                                    annotated$peak_start[same] - annotated$gene_end[same], 0)
  for (reason in c("duplicate_gene", "unmapped_gene", "different_contig")) {
    n <- sum(annotated$status == reason)
    if (n) warning(n, " peak-to-gene edges excluded: ", reason, call. = FALSE)
  }
  keep <- same & (trim_distance == 0 | (!is.na(annotated$distance) & annotated$distance < trim_distance))
  annotated$status[same & !keep] <- "outside_distance"
  trim_result <- merge(merge(annotated[keep, , drop = FALSE], tfpk, by = "peak"), tfgn, by = c("TF", "gene"))
  trim_result <- trim_result[, c("TF", "gene", "peak", "TF_GN_TE", "TF_PK_TE", "PK_GN_TE", "distance"), drop = FALSE]
  trim_result <- trim_result[order(trim_result$TF, trim_result$gene, trim_result$peak), , drop = FALSE]
  rownames(trim_result) <- NULL
  if (isTRUE(save)) {
    if (is.null(root)) stop("Tenet_result_dir is required when save=TRUE.")
    # Retain the original default result filename for downstream scripts.
    stem <- paste0("trimmed_by_PeaksourceDistance", convertNumberToUnit(trim_distance),
                   if (trim_indirect) "_Indirect-0.01" else "", if (fdr != 0.01) paste0("_fdr", alpha) else "")
    write.table(trim_result, file.path(root, paste0(stem, "_TE_result_matrix.sif")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(annotated, file.path(root, paste0(stem, "_peak_gene_annotations.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    degrees <- if (nrow(trim_result)) {
      counts <- table(trim_result$TF)
      data.frame(TF = names(counts), outdegree = as.integer(counts), stringsAsFactors = FALSE)
    } else data.frame(TF = character(), outdegree = integer())
    degrees <- head(degrees[order(-degrees$outdegree), , drop = FALSE], 20)
    plot <- ggplot(degrees, aes(x = reorder(TF, outdegree), y = outdegree)) +
      geom_col() + coord_flip() + theme_classic() + labs(x = NULL, y = "Triplet count")
    ggsave(file.path(root, paste0(sub("^trimmed_by_", "Outdegree_", stem), ".pdf")),
           plot = plot, width = 7, height = 5)
  }
  attr(trim_result, "peak_gene_annotations") <- annotated
  trim_result
}


pseudotime_heatmap2 <- function(matrix,selected_gene,gene_list,peak_list=F,ues_cell_select=F,cell_select=F,total_peak_gene,pseudotime, span=0.5,use_pseudotime_origin = T, use_z_score = T,
                                order_pseudo_score=F,p_min=-0.7,p_max=0.7,p_legend=T,filename = "pseudotime_heatmap",out_dataset=NULL,cell_aanotaion=NULL,
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
  if(is.null(cell_aanotaion)==F){
    cell_aanotaion <- cbind(pseudotime,cell_aanotaion)
    cell_aanotaion <- subset(cell_aanotaion,cell_aanotaion[,1]!="Inf" & cell_aanotaion[,1]!="NA")
    cell_pseudo_order <- order(cell_aanotaion[,1])
    cell_aanotaion <- cell_aanotaion[cell_pseudo_order,]
    cell_aanotaion_list <- as.data.frame(cell_aanotaion[,2])
    rownames(cell_aanotaion_list) <- rownames(cell_aanotaion)
    cell_type_list <- unique(cell_aanotaion[,2])
    cell_type_list
    color_list <- c('#F68282','#31C53F','#1FA195','#B95FBB','#D4D915','#28CECA','#ff9a36','#2FF18B','#aeadb3','#faf4cf',
                    '#CCB1F1','#25aff5','#A4DFF2','#4B4BF7','#AC8F14','#E6C122')
    aano_colors <- c()
    kk <- 1
    for (z in cell_type_list) {
      aano_colors_temp <- c(paste0('"', z, '" = "', color_list[kk], '"'))
      aano_colors <- c(aano_colors,aano_colors_temp)
      kk <- kk+1
    }
    aano_colors
    aano_colors <- as.list(aano_colors)
    aano_colors
    if(max_min==T){
      s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1],1:dim(final_matrix)[2]],
                     show_colnames=F, cluster_rows = F, cluster_cols = F,legend = p_legend,annotation = cell_aanotaion_list,
                     annotation_colors = anno_colors,breaks = seq(p_min,p_max,length.out=100),main=filename,fontsize=fontsize)
    }else{
      s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1],1:dim(final_matrix)[2]],
                     show_colnames=F, cluster_rows = F, cluster_cols = F,legend = p_legend,annotation = cell_aanotaion_list,
                     annotation_colors = anno_colors,breaks = seq(-0.4,0.4,length.out=100),main=filename,fontsize=fontsize) 
    }
    
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
