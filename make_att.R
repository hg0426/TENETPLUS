library(dplyr)
library(tidyverse)

TF_gene <- read.delim("C:/Users/LSCS/Desktop/tenet/new/CD8_merge_TF/TF_gene.sif", header=FALSE)
TF_peak <- read.delim("C:/Users/LSCS/Desktop/tenet/new/CD8_merge_TF/TF_peak.sif", header=FALSE)

#==================unique_TF=====================
TF <- unique(TF_gene$V1)
TF2 <- unique(TF_peak$V1)
mer_TF <- c(TF,TF2)
unique_TF <- unique(mer_TF)
unique_TF <- as.data.frame(unique_TF)
unique_TF <- cbind(unique_TF,"TF")
colnames(unique_TF) <- c("gene,paek","stat")
#==================unique_gene=====================
gene_mer <- unique(TF_gene$V3)
gene <- list()
for (i in 1:length(gene_mer)) {
  if (any(gene_mer[i] %in% mer_TF)) {
    
  }else{
    gene <- append(gene,gene_mer[i])
  }
}
unique_gene <- as.character(gene)
unique_gene <- as.data.frame(unique_gene)
unique_gene <- cbind(unique_gene,"gene")
colnames(unique_gene) <- c("gene,paek","stat")
#==================unique_peak=====================
unique_peak <- unique(TF_peak$V3)
unique_peak <- as.data.frame(unique_peak)
unique_peak <- cbind(unique_peak,"peak")
colnames(unique_peak) <- c("gene,paek","stat")
#==================merge data====================
atribute <- rbind(unique_TF,unique_gene,unique_peak)

write.table(atribute,"gene_paek_stat.txt", quote = F, sep = "\t",  row.names = FALSE, col.names = FALSE)
write.table(atribute,"gene_paek_stat.att", quote = F, sep = "\t",  row.names = FALSE, col.names = FALSE)

