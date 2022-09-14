count_degree <- function(data,decreasing =T,degree = "out"){
  if (degree == "out") {
    att <- as.data.frame(table(data[,1]))[order(as.data.frame(table(data[,1]),)[,2], decreasing = decreasing),]
  }
  if (degree == "in") {
    att <- as.data.frame(table(data[,3]))[order(as.data.frame(table(data[,3]),)[,2], decreasing = decreasing),]
  }
  
  return(att)
}
