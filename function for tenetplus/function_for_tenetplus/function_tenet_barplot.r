TE_barplot <- function(gene_list=as.data.frame("non","non"),peak_list=as.data.frame("non","non"),
                       PAG_list=as.data.frame("non","non"),top_gene_num=20,highlight_gene="",ven = F,fill_col = "grey",
                       highlight_col = "brown4",save_files = F,
                       name_RNA="TE_result_RNA_count_outdegree",
                       name_ATAC="TE_result_ATAC_count_outdegree",
                       name_trim="TE_result_trim_count_outdegree") {
  #======================RNA_Gene=======================================
  if (gene_list[1,1] != "non") {
    att <- as.data.frame(table(gene_list[,1]))[order(as.data.frame(table(gene_list[,1]),)[,2], decreasing = T),]
    top_genes <- att[1:top_gene_num,]
    top_RNA <- top_genes
    cols <- rep(fill_col,times=20)
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
    ggsave(paste0(name_RNA,".png"))
    
  }

  #===================ATAC_Peak==================================
  if (peak_list[1,1] != "non") {
    att <- as.data.frame(table(peak_list[,1]))[order(as.data.frame(table(peak_list[,1]),)[,2], decreasing = T),]
    top_genes <- att[1:top_gene_num,]
    top_ATAC <- top_genes
    cols <- rep(fill_col,times=20)
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
    ggsave(paste0(name_ATAC,".png"))
  }
  #=================trimmed_by_PAG======================
  if (PAG_list[1,1] != "non") {
    att <- as.data.frame(table(PAG_list[,1]))[order(as.data.frame(table(PAG_list[,1]),)[,2], decreasing = T),]
    top_genes <- att[1:top_gene_num,]
    top_Trim <- top_genes
    cols <- rep(fill_col,times=20)
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
    ggsave(paste0(name_trim,".png"))
  
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

      