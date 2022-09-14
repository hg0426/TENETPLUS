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