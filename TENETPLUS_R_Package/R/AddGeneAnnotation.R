#' Add Gene Annotation to TENET Object
#'
#' This function reads a gene annotation file from the package's `inst/bash/`
#' directory and adds it to the `feature_list$gtf` field of the TENET object.
#'
#' @param obj A TENET object.
#' @param species The species for which the annotation is required. Default is `"human"`.
#' @param version The genome version (e.g., `"hg38"`, `"hg19"`, `"mm10"`). Default is `"hg38"`.
#' @return The TENET object with updated `feature_list$gtf`.
#' @export
AddGeneAnnotation <- function(obj, species = "human", version = "hg38") {
  
  # Validate input
  if (!inherits(obj, "TENET")) {
    stop("obj must be a TENET object.")
  }

  # Construct the file name dynamically
  annotation_filename <- paste0(species, "_", version, "_gtf.bed")
  annotation_file <- system.file("bash", annotation_filename, package = "TENETPLUS")

  # Check if the file exists
  if (!file.exists(annotation_file)) {
    stop("Annotation file not found: ", annotation_file)
  }

  # Read the gene annotation file
  gene_anno <- tryCatch(
    {
      read.table(annotation_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    },
    error = function(e) {
      stop("Failed to read annotation file: ", annotation_file)
    }
  )

  # Ensure the file contains at least 4 columns
  if (ncol(gene_anno) < 4) {
    stop("Invalid annotation file format. Expected at least 4 columns.")
  }

  # Keep only relevant columns and rename them
  gene_anno <- gene_anno[, 1:4]
  colnames(gene_anno) <- c("gene_chr", "gene_start", "gene_end", "gene")

  # Store in TENET object
  obj@input$feature_list$gtf <- gene_anno

  message("Successfully added gene annotation for ", species, " (", version, ").")

  return(obj)
}

