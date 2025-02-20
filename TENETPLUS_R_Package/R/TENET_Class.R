#' TENET Class Definition
#'
#' Defines the S4 class TENET with slots for input and output.
#'
#' @name TENET
#' @slot input A list containing RNA/ATAC counts, metadata, and feature list.
#' @slot output A list to store processed data/results.
#' @export
setClass(
  "TENET",
  slots = list(
    input = "list",
    output = "list"
  )
)

#' Validity Check for TENET Class
#'
#' Ensures that the TENET object contains required elements such as count matrices,
#' metadata, and feature lists, and that they have the correct format.
#'
#' @param object A TENET object.
#' @return Returns TRUE if the object is valid; otherwise, returns an error message.
#' @export
setValidity("TENET", function(object) {
  required_input_keys <- c("count", "metadata", "feature_list")
  if (!all(required_input_keys %in% names(object@input))) {
    return(paste("input list must contain:", paste(required_input_keys, collapse = ", ")))
  }

  # Check RNA and ATAC count matrices
  count_list <- object@input$count
  if (!all(c("rna_counts", "atac_counts") %in% names(count_list))) {
    return("input$count must contain 'rna_counts' and 'atac_counts'.")
  }

  rna <- count_list$rna_counts
  atac <- count_list$atac_counts
  if (!(is.matrix(rna) || inherits(rna, "dgCMatrix"))) {
    return("rna_counts must be a matrix or a dgCMatrix.")
  }
  if (!(is.matrix(atac) || inherits(atac, "dgCMatrix"))) {
    return("atac_counts must be a matrix or a dgCMatrix.")
  }
  if (ncol(rna) != ncol(atac)) {
    return("rna_counts and atac_counts must have the same number of columns (cells).")
  }

  # Check metadata
  metadata_df <- object@input$metadata
  if (!is.data.frame(metadata_df)) {
    return("metadata must be a data.frame.")
  }
  if (is.null(rownames(metadata_df)) || nchar(rownames(metadata_df)[1]) == 0) {
    return("metadata must have row names corresponding to cell names.")
  }

  return(TRUE)
})

