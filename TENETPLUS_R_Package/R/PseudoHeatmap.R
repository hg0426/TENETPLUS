#' Generate Pseudotime Heatmaps for Transcription Factors (TF), Target Genes (TG), and Target Accessible Regions (TAR)
#'
#' This function generates a pseudotime-ordered heatmap for selected transcription factors (TFs) along with their 
#' target genes (TG) and target accessible regions (TAR) based on gene expression (RNA-seq) and chromatin accessibility (ATAC-seq) data.
#'
#' @param object A list or S4 object containing input and output data. Expected structure:
#'   \itemize{
#'     \item \code{object@input$count$rna_counts}: A matrix of RNA expression counts (rows: genes, columns: cells).
#'     \item \code{object@input$count$atac_counts}: A matrix of ATAC accessibility counts (rows: peaks, columns: cells).
#'     \item \code{object@input$metadata}: A data frame containing metadata with \code{pseudotime} and \code{select} columns.
#'     \item \code{object@output$result$Trimm}: A data frame with TF-TG-TAR relationships.
#'   }
#'
#' @param features Character vector of transcription factors (TFs) to include in the heatmap.
#' @param group Character (optional). Specifies the grouping variable in \code{object@input$metadata} for heatmap annotation.
#' @param use_result Logical. If \code{TRUE}, includes target genes (TG) and target accessible regions (TAR) heatmaps. Default is \code{TRUE}.
#' @param span Numeric. The span parameter for LOESS smoothing. Default is \code{0.7}.
#'
#' @return Generates heatmaps showing pseudotime-ordered expression patterns of TFs, TGs, and TARs.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts metadata and pseudotime information, filtering cells where \code{select == 1}.
#'   \item Retrieves RNA expression data for selected TFs, standardizes it (Z-score transformation), and applies LOESS smoothing.
#'   \item Generates a heatmap of TF expression over pseudotime.
#'   \item If \code{use_result} is \code{TRUE}:
#'     \itemize{
#'       \item Extracts TG expression data, computes average expression across targets, and applies LOESS smoothing.
#'       \item Extracts TAR accessibility data from ATAC-seq, processes it similarly to TG expression.
#'       \item Generates heatmaps for TG and TAR.
#'       \item Merges the TF, TG, and TAR heatmaps vertically (\code{\%v\%} operator).
#'     }
#'   \item If \code{group} is provided, adds cell-type annotation to the heatmap.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' PseudoHeatmap(object, features = c("FOXO1", "TP53"), group = "cell_type", use_result = TRUE, span = 0.7)
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import stringr
#' @import circlize
#' @import grid
#'
#' @export
PseudoHeatmap <- function(object,
                          features,
                          group = NULL,
                          use_result = TRUE,
                          span=0.7
                          ) {
  # 1. 예외 처리: features가 object의 rna_counts에 모두 포함되어 있는지 확인
  if (!all(features %in% rownames(object@input$count$rna_counts))) {
    stop("Error: Some features are not present in object@input$count$rna_counts")
  }
  
  # 3. select가 1인 것 subset, pseudotime ordering
  select_idx <- rownames(object@input$metadata)[object@input$metadata$select == 1]
  meta_sub <- object@input$metadata[select_idx,]
  meta_sub <- meta_sub[order(meta_sub$pseudotime),]
  
  # 4. count matrix (cell subset)
  count_matrix <- t(as.matrix(object@input$count$rna_counts[,rownames(meta_sub) , drop = FALSE]))
  count_matrix[1:3,1:3];dim(count_matrix)
  
  # 5. pseudotime을 1부터 순차적으로 재정의
  pseudo_idx <- 1:nrow(count_matrix)
  count_matrix <- cbind(pseudo_idx = pseudo_idx, count_matrix)
  count_matrix[1:3,1:3];dim(count_matrix)
  
  # 6. TF 자체 expression
  count_matrix_tf <- count_matrix[,c("pseudo_idx",features)]
  head(count_matrix_tf)
  count_matrix_tf[, -1] <- apply(count_matrix_tf[, -1, drop = FALSE], 2, function(x) as.numeric(scale(x)))
  count_matrix_tf[1:3,1:3];dim(count_matrix_tf)
  
  # 9. LOESS smoothing: 각 feature에 대해 pseudotime을 이용하여 스무딩 적용
  count_matrix_tf[, -1] <- apply(count_matrix_tf[, -1, drop = FALSE], 2, function(y) {
    lo <- loess(y ~ count_matrix_tf[, 1], span = span)
    predict(lo, newdata = count_matrix_tf[, 1])
  })
  count_matrix_tf[1:3,1:3];dim(count_matrix_tf)
  
  if (use_result){
    col_title = NULL
    row_title = "TF"
  } else {
    col_title = "pseudotime"
    row_title = NULL
  }
  # 9. Pseudo-heatmap
  if (is.null(group)){
    ht <- ComplexHeatmap::Heatmap(t(count_matrix_tf[,-1]), use_raster = FALSE,
                            cluster_columns = F,
                            cluster_rows = F,
                            show_column_names = F,
                            show_row_dend = F,
                            column_title = col_title, 
                            column_title_side = "bottom",
                            row_title = row_title,
                            name="value",
                            column_title_gp = gpar(fontsize = 10),
                            row_names_gp = gpar(fontsize = 10),
                            show_heatmap_legend = T)
  } else {
    orders <- meta_sub %>% 
      group_by(get(group)) %>% 
      summarise(median=median(pseudotime)) %>% 
      arrange(median)
    pseudo_type <- meta_sub[[group]]
    colors <- brewer.pal(length(unique(pseudo_type)), "Set1")
    names(colors) <- unique(pseudo_type)
    ha = HeatmapAnnotation("celltype" = pseudo_type,
                           col = list("celltype" = colors),
                           annotation_legend_param = list(
                             "celltype" = list(
                               title = "celltype", 
                               title_position = "topcenter",
                               at = orders[[1]])
                           ),
                           simple_anno_size=unit(0.3, "cm"),
                           annotation_name_gp= gpar(fontsize = 10)
    )
    ht <- ComplexHeatmap::Heatmap(t(count_matrix_tf[,-1]), use_raster = FALSE,
                            cluster_columns = F,
                            cluster_rows = F,
                            show_column_names = F,
                            show_row_dend = F,
                            column_title = col_title,
                            column_title_side = "bottom",
                            row_title = row_title,
                            name="value",
                            column_title_gp = gpar(fontsize = 10),
                            row_names_gp = gpar(fontsize = 10),
                            show_heatmap_legend = T,
                            ,top_annotation = ha)
  }

  # Use TENET+ result (TG) -------------------------------------------------------
  if (use_result){
    # 2. features가 존재하는 Trimm TENET 결과 가져오기
    result <- object@output$result$Trimm
    result <- result[result$TF %in% features, ]
    head(result);dim(result)
    
    # if (!all(unique(result$TG) %in% rownames(object@input$count$rna_counts))) {
    #   stop("Error: Some TGs are not present in object@input$count$rna_counts")
    # }
    
    count_matrix_avg <- data.frame(pseudo_idx = count_matrix[, "pseudo_idx"])
    for (feature in features) {
      #print(feature)
      result_sub <- result[result$TF == feature,]
      head(result_sub);dim(result_sub)
      targets <- unique(result_sub$TG)
      length(targets)
      
      ### Tmp code
      # targets <- targets[targets %in% rownames(object@input$count$rna_counts)]
      # targets <- "FOXO1"
      
      count_matrix_sub <- count_matrix[,c("pseudo_idx",targets)]
      
      avg_expr <- if (ncol(count_matrix_sub) > 1) rowMeans(count_matrix_sub) else as.vector(count_matrix_sub)
      count_matrix_avg[[feature]] <- avg_expr
    }
    count_matrix_avg[1:3,1:3];dim(count_matrix_avg)
    
    # 6. TF 자체 expression
    count_matrix_avg[, -1] <- apply(count_matrix_avg[, -1, drop = FALSE], 2, function(x) as.numeric(scale(x)))
    count_matrix_avg[1:3,1:3];dim(count_matrix_avg)
    
    # 9. LOESS smoothing: 각 feature에 대해 pseudotime을 이용하여 스무딩 적용
    count_matrix_avg[, -1] <- apply(count_matrix_avg[, -1, drop = FALSE], 2, function(y) {
      lo <- loess(y ~ count_matrix_avg[, 1], span = span)
      predict(lo, newdata = count_matrix_avg[, 1])
    })
    count_matrix_avg[1:3,1:3];dim(count_matrix_avg)
    
    ht2 <- ComplexHeatmap::Heatmap(t(count_matrix_avg[,-1]), use_raster = FALSE,
                                   cluster_columns = F,
                                   cluster_rows = F,
                                   show_column_names = F,
                                   show_row_dend = F,
                                   column_title = "pseudotime", 
                                   #column_title = NULL, 
                                   column_title_side = "bottom",
                                   name="value",
                                   column_title_gp = gpar(fontsize = 10),
                                   row_names_gp = gpar(fontsize = 10),
                                   show_heatmap_legend = T,
                                   row_title = "TG")
  } 
  
  # Use TENET+ result (TAR) -------------------------------------------------------
  if (use_result){
    # 2. features가 존재하는 Trimm TENET 결과 가져오기
    result <- object@output$result$Trimm
    result <- result[result$TF %in% features, ]
    head(result)
    
    if (!all(unique(result$TAR) %in% rownames(object@input$count$atac_counts))) {
      stop("Error: Some TARs are not present in object@input$count$atac_counts")
    }
    
    # 4. count matrix (cell subset)
    count_matrix_atac <- t(as.matrix(object@input$count$atac_counts[,rownames(meta_sub) , drop = FALSE]))
    count_matrix_atac[1:3,1:3];dim(count_matrix_atac)
    
    # 5. pseudotime을 1부터 순차적으로 재정의
    pseudo_idx <- 1:nrow(count_matrix_atac)
    count_matrix_atac <- cbind(pseudo_idx = pseudo_idx, count_matrix_atac)
    count_matrix_atac[1:3,1:3];dim(count_matrix_atac)
    
    count_matrix_avg_atac <- data.frame(pseudo_idx = count_matrix_atac[, "pseudo_idx"])
    for (feature in features) {
      #print(feature)
      result_sub <- result[result$TF == feature,]
      head(result_sub);dim(result_sub)
      targets <- unique(result_sub$TAR)
      length(targets)
      
      ### Tmp code
      #targets <- targets[targets %in% rownames(object@input$count$atac_counts)]
      #targets <- "FOXO1"
      
      count_matrix_sub <- count_matrix_atac[,c("pseudo_idx",targets)]
      
      avg_expr <- if (ncol(count_matrix_sub) > 1) rowMeans(count_matrix_sub) else as.vector(count_matrix_sub)
      count_matrix_avg_atac[[feature]] <- avg_expr
    }
    count_matrix_avg_atac[1:3,1:3];dim(count_matrix_avg_atac)
    
    # 6. TF 자체 expression
    count_matrix_avg_atac[, -1] <- apply(count_matrix_avg_atac[, -1, drop = FALSE], 2, function(x) as.numeric(scale(x)))
    count_matrix_avg_atac[1:3,1:3];dim(count_matrix_avg_atac)
    
    # 9. LOESS smoothing: 각 feature에 대해 pseudotime을 이용하여 스무딩 적용
    count_matrix_avg_atac[, -1] <- apply(count_matrix_avg_atac[, -1, drop = FALSE], 2, function(y) {
      lo <- loess(y ~ count_matrix_avg_atac[, 1], span = span)
      predict(lo, newdata = count_matrix_avg_atac[, 1])
    })
    count_matrix_avg_atac[1:3,1:3];dim(count_matrix_avg_atac)
    
    ht3 <- ComplexHeatmap::Heatmap(t(count_matrix_avg_atac[,-1]), use_raster = FALSE,
                                   cluster_columns = F,
                                   cluster_rows = F,
                                   show_column_names = F,
                                   show_row_dend = F,
                                   column_title = "pseudotime", 
                                   #column_title = NULL, 
                                   column_title_side = "bottom",
                                   name="value",
                                   column_title_gp = gpar(fontsize = 10),
                                   row_names_gp = gpar(fontsize = 10),
                                   show_heatmap_legend = T,
                                   row_title = "TAR")
    draw(ht %v% ht2 %v% ht3,merge_legend = TRUE)
  } else{
    draw(ht,merge_legend = TRUE)
  }
}
