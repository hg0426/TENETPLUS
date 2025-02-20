#' Triplet Trimming for Regulatory Network Analysis
#'
#' This function performs triplet trimming on regulatory network data stored in a given object.
#' It calculates the genomic distance between peak regions and gene annotations, and then trims
#' the triplet network based on a specified distance threshold and the choice between direct or
#' indirect transcription factor (TF) to target gene (TG) interactions.
#'
#' @param object An object (often an S4 or list) containing the necessary input and output data.
#'   It is expected to have:
#'   \itemize{
#'     \item \code{object@input$feature_list$gtf}: A data.frame with gene annotations, including columns
#'           \code{gene}, \code{gene_chr}, \code{gene_start}, and \code{gene_end}.
#'     \item \code{object@output$result$TAR_TG}: A data.frame with peak data. The \code{source} column should
#'           contain strings formatted as "chr-start-end", along with columns such as \code{score} and \code{target}.
#'     \item \code{object@output$result$TF_TG_indirect}: A data.frame of indirect TF-target interactions.
#'     \item \code{object@output$result$TF_TG}: A data.frame of direct TF-target interactions.
#'     \item \code{object@output$result$TF_TAR}: A data.frame of TF-peak interactions, with a \code{TAR} column.
#'   }
#'
#' @param trim_indirect Logical; if \code{TRUE}, indirect TF-target interactions (\code{TF_TG_indirect}) will
#'   be used for trimming. If \code{FALSE}, direct interactions (\code{TF_TG}) are used. Default is \code{TRUE}.
#'
#' @param trim_distance Numeric or character; specifies the maximum allowed distance between peaks and gene
#'   annotations for a valid interaction. When a numeric value is provided, only peak-gene pairs with a
#'   calculated distance less than \code{trim_distance} are retained. If set to \code{"None"}, no distance-based
#'   trimming is applied. Default is \code{1000000}.
#'
#' @return The input \code{object} with updated results in its \code{@output$result} slot:
#'   \itemize{
#'     \item \code{$TAR_TG}: A data.frame with added columns for peak chromosome (\code{peak_chr}),
#'           start (\code{peak_start}), end (\code{peak_end}), and the calculated \code{distance}
#'           between the peak and gene annotation.
#'     \item \code{$Trimm}: A data.frame representing the trimmed triplet regulatory network. It includes
#'           columns \code{TF}, \code{TG}, \code{TAR}, \code{TF_TG_TE}, \code{TF_TAR_TE}, \code{TAR_TG_TE},
#'           and \code{distance}.
#'   }
#'
#' @details
#' The function executes the following steps:
#' \enumerate{
#'   \item Extracts and parses the peak data from the \code{source} column to retrieve chromosome,
#'         start, and end positions.
#'   \item Processes the gene annotation data by ensuring numeric types for start and end positions,
#'         and removes duplicate gene entries.
#'   \item Merges the peak data with gene annotation data and verifies that the chromosome information
#'         matches. An error is raised if mismatches are found.
#'   \item Calculates the distance between each peak and the corresponding gene using a custom
#'         distance function.
#'   \item Depending on the \code{trim_indirect} flag, selects either indirect or direct TF-target
#'         interaction data for further processing.
#'   \item Merges the filtered peak-gene list with TF-peak and TF-target interaction data to form a
#'         trimmed triplet network based on the specified \code{trim_distance}.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_object' is a properly formatted object:
#' result_object <- TripletTrimming(my_object, trim_indirect = TRUE, trim_distance = 500000)
#' }
#'
#' @import stringr
#' @import dplyr
#' 
#' @export
TripletTrimming <- function(object,
                            trim_indirect=TRUE,
                            trim_distance=1000000
){
  # Calculate PeakSource Distance -------------------------------------------------------------
  print(paste0("##### Calculate TAR_TG Distance"))
  rowPeak_colGN <- object@output$result$TAR_TG
  rowPeak_colGN$peak_chr <- str_split(rowPeak_colGN$source,'-',simplify = T)[,1]
  rowPeak_colGN$peak_start <- as.numeric(str_split(rowPeak_colGN$source,'-',simplify = T)[,2])
  rowPeak_colGN$peak_end <- as.numeric(str_split(rowPeak_colGN$source,'-',simplify = T)[,3])
  
  ### Process gene_anno 
  gene_anno <- object@input$feature_list$gtf
  gene_anno$gene_start <- as.numeric(gene_anno$gene_start);gene_anno$gene_end <- as.numeric(gene_anno$gene_end)
  gene_anno <- gene_anno[!duplicated(gene_anno$gene),] 
  
  intersect_gene <- length(intersect(unique(rowPeak_colGN$target),gene_anno$gene))
  rowPeak_colGN_gene<- length(unique(rowPeak_colGN$target))
  print(paste0("- annotated gene/total gene: ",intersect_gene,"/",rowPeak_colGN_gene))
  if (rowPeak_colGN_gene>intersect_gene){
    tmp_gene <- unique(rowPeak_colGN$target)[!unique(rowPeak_colGN$target) %in% gene_anno$gene]
    print(paste0("Genes are not in gtf file: ",paste0(paste(head(tmp_gene), collapse = ", "),"...")))
  }
  rowPeak_colGN <- merge(rowPeak_colGN,gene_anno,all.x=T,by.x="target",by.y="gene",sort=F)
  if (identical(rowPeak_colGN$peak_chr, rowPeak_colGN$gene_chr) == FALSE){
    object@output$result$TAR_TG <- rowPeak_colGN
    return(object)
    stop(paste("Error: Chromosomes do not match between peak and gene in TAR_TG. Check 'peak_chr' and 'gene_chr' columns in TAR_TG"))
  }
  
  ### Calculate Distance
  calc_dist <- function(a1,a2,b1,b2){
    if (a2 < b1) {
      distance <- b1 - a2
    } else if (b2 < a1) {
      distance <- a1 - b2
    } else {
      distance <- 0  
    }
    return(distance)
  }
  rowPeak_colGN$distance <- apply(rowPeak_colGN,MARGIN=1,function(x) calc_dist(as.numeric(x["peak_start"]),as.numeric(x["peak_end"]),as.numeric(x["gene_start"]),as.numeric(x["gene_end"])))
  rowPeak_colGN %>% 
    dplyr::select(c("source","score","target","distance")) -> rowPeak_colGN
  object@output$result$TAR_TG <- rowPeak_colGN
  cat('\n')
  
  # TripletTrimming -------------------------------------------------------------
  print(paste0("##### TripletTrimming"))
  if (trim_indirect == TRUE){
    gene_list <- object@output$result$TF_TG_indirect
  } else {
    gene_list <- object@output$result$TF_TG
  }
  peak_list <- object@output$result$TF_TAR
  PK_gene_list <- object@output$result$TAR_TG
  colnames(gene_list) <- c("TF","TF_TG_TE","TG")
  colnames(peak_list) <- c("TF","TF_TAR_TE","TAR")
  colnames(PK_gene_list) <- c("TAR","TAR_TG_TE","TG","distance")
  
  ## Step 1: PK -> GN
  print(paste0("- distance: ", format(trim_distance, scientific = FALSE)))
  if (trim_distance != "None"){
    PK_gene_list_enhan=subset(PK_gene_list,distance < trim_distance)
  } else if (trim_distance == "None") {
    PK_gene_list_enhan=PK_gene_list
  }
  
  ## Step 2: TF -> PK
  PK_gene_list_enhan_peak <- merge(PK_gene_list_enhan,peak_list,by="TAR")
  
  ## Step 3: TF -> GN
  trim_result <- merge(PK_gene_list_enhan_peak,gene_list,by=c("TF","TG"))
  trim_result <- trim_result[,c("TF","TG","TAR","TF_TG_TE","TF_TAR_TE","TAR_TG_TE","distance")]
  length(paste0(trim_result$TF,"_",trim_result$gene,"_",trim_result$peak))
  length(unique(paste0(trim_result$TF,"_",trim_result$gene,"_",trim_result$peak)))
  print(paste("- Trimm's nrow: ",nrow(trim_result)))
  print(paste("- Trimm's unique TF: ",length(unique(trim_result$TF))))
  object@output$result$Trimm = trim_result
  cat('\n')
  return(object)
}



