library(dplyr)
trim_sif <- read.delim("C:/Users/LSCS/Desktop/tenet/new/Bcell_merge_TF/TE_result_matrix.fdr0.01.trimIndirect0.0.sif", header=FALSE)

TF_peak <- dplyr::filter(trim_sif, grepl("chr",trim_sif$V3))
TF_gene <- dplyr::filter(trim_sif, !grepl("chr",trim_sif$V3))

write.table(TF_peak,"TF_peak.sif", quote = F, sep = "\t",  row.names = FALSE, col.names = FALSE)
write.table(TF_gene,"TF_gene.sif", quote = F, sep = "\t",  row.names = FALSE, col.names = FALSE)



