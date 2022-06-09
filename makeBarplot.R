#!/usr/bin/env Rscript

### USAGE: Rscript makeBarplot.R [fullPath_out|indegree.txt]
### ex) Rscript makeBarplot.R /home/Data/tenet_result_kyu/Bcell_merge_normal/TE_result_matrix_colPK_rowTF.fdr0.01.trimIndirect0.0.sif.outdegree.txt

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)

inFnm=args[1] #TE_result_matrix_colPK_rowTF.fdr0.01.trimIndirect0.0.sif.outdegree.txt
outFnm=gsub(".txt", "_top20Ftr.png", inFnm)
print(outFnm)

dt=read.table(inFnm)
top20Ftr=dt[1:20,] #top20Ftr=na.omit(dt[1:20,])

p <- ggplot(data=top20Ftr, aes(x=reorder(V1, V2), y=V2)) + geom_bar(stat="identity") +
  theme(plot.margin=unit(c(0.5,0.5,0,0), 'cm'))+
  xlab("") + ylab("") + ggtitle(inFnm)  +
  theme(plot.title = element_text(hjust = 0.5, size = 7))
p + coord_flip()

ggsave(filename= outFnm, width=5, height=7)

