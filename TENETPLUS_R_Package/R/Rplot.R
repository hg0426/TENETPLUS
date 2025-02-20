#' Plot Degree and Expression Correlation
#'
#' This function calculates the degree distribution from a specified layer in a regulatory network object and correlates it with gene expression levels.
#' It extracts the nodes from the layer, computes their degree (either outdegree or indegree), retrieves the corresponding gene expression data,
#' and then visualizes the relationship between degree and expression using a scatter plot. The top \code{top_n} nodes by expression and degree are highlighted.
#'
#' @param object An object containing network data and gene expression data. It should include:
#'   \itemize{
#'     \item \code{object@output$result}: A list of data frames representing different network layers. The specified \code{layer} should be present here.
#'     \item \code{object@input$count$rna_counts}: A matrix or data frame of gene expression counts with gene names as rownames.
#'   }
#'
#' @param layer Character. The name of the layer in \code{object@output$result} to be used. Default is \code{"Trimm"}.
#'
#' @param degree_type Character. Specifies whether to compute \code{"outdegree"} or \code{"indegree"}.
#'   For the \code{"Trimm"} layer, \code{"outdegree"} is computed using the \code{TF} column and \code{"indegree"} using the \code{TG} column.
#'   For other layers, \code{"outdegree"} uses the \code{source} column and \code{"indegree"} uses the \code{target} column.
#'   Default is \code{"outdegree"}.
#'
#' @param top_n Numeric. The number of top nodes (by expression and degree) to highlight in the plot. Default is 5.
#'
#' @param return_df Logical. If \code{TRUE}, the function returns a data frame containing the node names, degree, and expression values.
#'   Otherwise, the function only prints the plot. Default is \code{FALSE}.
#'
#' @return A scatter plot is printed displaying the relationship between node degree and gene expression.
#'   If \code{return_df} is \code{TRUE}, a data frame with columns \code{node}, \code{degree}, and \code{expression} is returned.
#'
#' @details
#' The function operates as follows:
#' \enumerate{
#'   \item Validates that the specified \code{layer} exists in \code{object@output$result}.
#'   \item Extracts the relevant degree information based on the layer and \code{degree_type}.
#'   \item Retrieves gene expression data from \code{object@input$count$rna_counts} for the unique nodes found in the degree calculation.
#'   \item Aggregates and normalizes the expression data using \code{NormalizeData}.
#'   \item Highlights the top \code{top_n} nodes based on both expression and degree.
#'   \item Generates a scatter plot with \code{ggplot2} and uses \code{geom_text_repel} to label the highlighted nodes.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'net_obj' is a properly formatted object:
#' Rplot(net_obj, layer = "Trimm", degree_type = "outdegree", top_n = 5)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import stringr
#' @import ggrepel
#' @import Matrix
#' @import Seurat
#'
#' @export
Rplot <- function(object,
                  layer = "Trimm",
                  degree_type = "outdegree",
                  top_n = 5,
                  return_df = FALSE) {
  # 예외 처리: 해당 layer가 존재하지 않을 경우 오류 메시지 출력
  if (!(layer %in% names(object@output$result))) {
    stop(paste("Error: Layer", layer, "does not exist in object@output$result"))
  }
  
  # 결과 데이터 가져오기
  result <- object@output$result[[layer]]
  
  # Degree 계산
  if (layer == "Trimm") {
    if (degree_type == "outdegree") {
      degree <- as.data.frame(sort(table(result$TF), decreasing = TRUE))
      var_id <- "TF"
    } else if (degree_type == "indegree") {
      degree <- as.data.frame(sort(table(result$TG), decreasing = TRUE))
      var_id <- "TG"
    }
  } else {
    if (degree_type == "outdegree") {
      degree <- as.data.frame(sort(table(result$source), decreasing = TRUE))
      var_id <- strsplit(layer, "_")[[1]][1]
    } else if (degree_type == "indegree") {
      degree <- as.data.frame(sort(table(result$target), decreasing = TRUE))
      var_id <- strsplit(layer, "_")[[1]][2]
    }
  }
  
  # 컬럼명 설정
  names(degree) <- c("Var", "Freq")
  unique_feature <- as.character(unique(degree$Var))
  
  # Expression 데이터 가져오기
  table(unique_feature %in% rownames(object@input$count$rna_counts))
  if (all(unique_feature %in% rownames(object@input$count$rna_counts)) == FALSE){
    stop("Error: Some genes are not found in object@input$count$rna_counts.")
  }
  exp_data <- as.matrix(object@input$count$rna_counts[unique_feature, ])
  aggregated_exp_data <- Matrix::rowSums(exp_data)
  
  # 정규화
  normalized_data <- NormalizeData(as.data.frame(aggregated_exp_data))
  normalized_data <- normalized_data[unique_feature, ]
  degree_exp_data <- data.frame(Node = unique_feature, Degree = degree$Freq, Expression = normalized_data)
  
  # Expression 값이 가장 높은 5개
  top_exp <- degree_exp_data[order(-degree_exp_data$Expression), ][1:top_n, ]
  
  # Outdegree 값이 가장 높은 5개
  top_degree <- degree_exp_data[order(-degree_exp_data$Degree), ][1:top_n, ]
  
  # 두 개 합치기 (중복 제거)
  highlight_nodes <- unique(rbind(top_exp, top_degree))
  
  # ggplot 시각화
  p <- ggplot(degree_exp_data, aes(x = Degree, y = Expression)) +
    geom_point(alpha = 0.6,size=2) +  
    geom_text_repel(data = highlight_nodes, aes(label = Node), size = 5, color = "red") +
    theme_bw() +
    labs(
      title = "R Plot",
      x = degree_type,
      y = "expression"
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 13, color = 'black'),
      axis.text.y = element_text(size = 13, color = 'black')
    )
  
  # 그래프 출력
  print(p)
  
  if (return_df) {
    names(degree_exp_data) <- c("node","degree","expression")
    return(degree_exp_data)
  }
}
