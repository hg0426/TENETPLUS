#' Plot Degree Distribution of Regulatory Network Layers
#'
#' This function creates a bar plot displaying the degree distribution (either outdegree or indegree)
#' for a specified layer within a regulatory network stored in an object. The function extracts the
#' relevant data from the chosen layer, computes the frequency of occurrence for each node (e.g., transcription
#' factors or target genes), and visualizes the top \code{top_n} nodes using \code{ggplot2}.
#'
#' @param object An object containing the regulatory network results in its \code{output$result} slot.
#'   Each layer within this slot should be a data frame with the appropriate columns.
#'
#' @param layer Character. The name of the layer within \code{object@output$result} to be used for plotting.
#'   For example, \code{"Trimm"} is commonly used to represent the trimmed network layer.
#'
#' @param degree_type Character. Specifies whether to compute the \code{"outdegree"} or \code{"indegree"}.
#'   For the \code{"Trimm"} layer, \code{"outdegree"} calculates the frequency of the \code{TF} column,
#'   while \code{"indegree"} uses the \code{TG} column. For other layers, \code{"outdegree"} uses the \code{source}
#'   column and \code{"indegree"} uses the \code{target} column.
#'
#' @param top_n Numeric. The number of top nodes (by frequency) to display in the plot. The default is 20.
#'
#' @param return_df Logical. If \code{TRUE}, the function returns the computed degree data frame in addition to
#'   printing the plot. The data frame contains two columns: \code{"Var"} (the node) and \code{"Freq"} (its frequency).
#'
#' @return If \code{return_df} is \code{TRUE}, a data frame containing the degree distribution is returned.
#'   Otherwise, the function prints the plot without returning any value.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that the specified \code{layer} exists in \code{object@output$result}. If not, an error is thrown.
#'   \item Extracts the data for the selected layer.
#'   \item Computes the degree distribution:
#'     \itemize{
#'       \item For the \code{"Trimm"} layer:
#'         \itemize{
#'           \item \code{"outdegree"} is computed from the \code{TF} column.
#'           \item \code{"indegree"} is computed from the \code{TG} column.
#'         }
#'       \item For other layers:
#'         \itemize{
#'           \item \code{"outdegree"} is computed from the \code{source} column.
#'           \item \code{"indegree"} is computed from the \code{target} column.
#'         }
#'       \item The variable identifier (\code{var_id}) is set based on the chosen column.
#'     }
#'   \item Generates a bar plot of the top \code{top_n} nodes using \code{ggplot2}, with a flipped coordinate system
#'         for better readability.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'net_obj' is an object with appropriate regulatory network data:
#' DegreePlot(net_obj, layer = "Trimm", degree_type = "outdegree", top_n = 15)
#'
#' # To return the degree data frame instead of just plotting:
#' degree_df <- DegreePlot(net_obj, layer = "Trimm", degree_type = "indegree", return_df = TRUE)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import stringr
#'
#' @export 
DegreePlot <- function(object,
                       layer = "Trimm",
                       degree_type = "outdegree",
                       top_n = 20,
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
      var_id <- str_split(layer, "_", simplify = TRUE)[, 1]
    } else if (degree_type == "indegree") {
      degree <- as.data.frame(sort(table(result$target), decreasing = TRUE))
      var_id <- str_split(layer, "_", simplify = TRUE)[, 2]
    }
  }
  # 컬럼명 설정
  names(degree) <- c("Var", "Freq")
  degree_type2 <- paste0(toupper(substr(degree_type, 1, 1)), substr(degree_type, 2, nchar(degree_type)))
  
  # ggplot 그래프 생성
  p <- head(degree, top_n) %>% 
    ggplot() +
    geom_col(aes(x = reorder(Var, Freq), y = Freq), fill = "steelblue") +  
    coord_flip() +
    theme_bw() +
    labs(
      title = paste0("Top ", degree_type2, " (", layer, ")"),  # 그래프 제목
      x = var_id,  # x축 라벨
      y = "Degree"  # y축 라벨
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # 제목 크기, 굵기, 중앙 정렬
      axis.title.x = element_text(size = 15, face = "bold"),  # x축 제목 크기 및 굵기
      axis.title.y = element_text(size = 15, face = "bold"),  # y축 제목 크기 및 굵기
      axis.text.x = element_text(size = 13, color = 'black'),  # x축 눈금 크기
      axis.text.y = element_text(size = 13, color = 'black')   # y축 눈금 크기
    )
  
  # 그래프 출력
  print(p)
  
  if (return_df == TRUE){
    return(degree) 
  }
}
