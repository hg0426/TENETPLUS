source("~/function_pseudo_heatmap.r")

# matrix in put (row = barcode, col = Gene)
# psuedotime (row = barcodes, col = pseudotime [one column], * make sure if Pseudotime is not defined should be set to NA or Inf)
# gene_list is just character

# span = weight / default = 0.6
# use_pseudotime_origin = T -> use original pseudotime (0.1,0.2,0.23....) /F use pseudotime order (1,2,3....) /  default = T
# use_z_score (when use F option max_min must be "F" and also recommended to lower the span value) / default = T
# order_pseudo_score Sort the heatmap using the concept of pseudoscore / default = F
# max_min pheatmap Specifies the maximum and minimum values in the pheatmap. / default = T
# p_min,p_max Only works when max_min = T / default = +- 0.7
# p_legend pheatmap legend / default = T
# out_result for make object *how to use (object name <- pseudotime_heatmap(~~~~, out_result = T)) / default = F
# filename save file name "xxx" -> xxx.png / defualt "pseudotime_heatmap"

library(pheatmap)

sample_matrix <- read.table("~/sample_matix.txt", quote="\"", comment.char="")
sample_pseudotime <- read.table("~/pseudotime.txt", quote="\"", comment.char="")
sample_gene_list <- read.table("~/gene_list.txt", quote="\"", comment.char="")
sample_gene_list <- as.character(unlist(sample_gene_list))

pseudotime_heatmap(matrix = sample_matrix,selected_gene = sample_gene_list,pseudotime = sample_pseudotime,
                   span = 0.6,use_pseudotime_origin = F,use_z_score = T,
                   order_pseudo_score = T,max_min = F, p_max = 1, p_min = -0.8, out_result = F, filename = "test",pseudo_decrease = F
                   ,target_average = T,gene_list = gene_list, peak_list = peak_list, total_peak_gene = total_peak_gene, target = "gene")






# ==============for confirm ==============

confirm <- pseudotime_heatmap_for_confirm(matrix = sample_matrix,selected_gene = sample_gene_list,gene_list = gene_list,
                                          peak_list = peak_list,total_peak_gene = total_peak_gene, pseudotime = sample_pseudotime,span = 0.6,use_pseudotime_origin = F,
                                          use_z_score = T,order_pseudo_score = F,target_average = T, target = "gene")

# ======== how many gene shows in one slide========
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
# ================for confirm when use order_pseudo_score = T ===============
for (i in 2:4){
  x <- as.numeric(confirm[1,2:dim(confirm)[2]])
  y <- as.numeric(confirm[i,2:dim(confirm)[2]])
  lo <- loess(y~x, span= 0.8)
  xl <- seq(min(x),max(x), (max(x) - min(x))/(length(y)-1))
  out = predict(lo,xl)
  plot(x,y,type="p", main = rownames(confirm)[i],ylim = c(-0.8,0.8),xlab = "",ylab = "")
  lines(xl, out, col='red', lwd=2)
}
