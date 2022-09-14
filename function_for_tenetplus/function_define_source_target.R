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
