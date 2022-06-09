library(ggplot2)
data <- read.delim("C:/Users/hg042/Downloads/TENET_TF_pbmc/fdr0.01/mono_merge_TF/TF_gene.sif.outdegree.txt", header=FALSE)
data2 <- read.delim("C:/Users/hg042/Downloads/TENET_TF_pbmc/fdr0.01/mono_merge_TF/TF_peak.sif.outdegree.txt", header=FALSE)

top20 <- data[1:20,]
top20_2 <- data2[1:8,]

name <- "TENET_TF\nmono_gene_top20"
p <- ggplot(data=top20, aes(x=reorder(V1, V2), y=V2)) + geom_bar(stat="identity") +
  theme(plot.margin=unit(c(0.5,5,0.5,5), 'cm'))+
  xlab("") + ylab("") + ggtitle(name)  +
  theme(plot.title = element_text(hjust = 0.5))


p + coord_flip()

ggsave("image/mono_gene_top20.png")

name <- "TENET_TF\nmono_peak_top20"
p <- ggplot(data=top20_2, aes(x=reorder(V1, V2), y=V2)) + geom_bar(stat="identity") +
  theme(plot.margin=unit(c(0.5,5,0.5,5), 'cm'))+
  xlab("") + ylab("") + ggtitle(name)  +
  theme(plot.title = element_text(hjust = 0.5))


p + coord_flip()

ggsave("image/mono_peak_top20.png")
