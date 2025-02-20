#' Convert Seurat Object to TENET Object
#'
#' Converts a Seurat object (with RNA and optional ATAC assays) into a TENET object.
#'
#' @param seurat_obj A Seurat object.
#' @param Species A character string indicating the species.
#' @param DEG A list or data frame for DEG information.
#' @param DAR A list or data frame for DAR information.
#' @param TF_list A list or data frame for TF_list information.
#' @param pseudotime Column name or data frame with pseudotime values.
#' @param cell_select Column name or data frame indicating cell selection.
#'
#' @return A TENET object.
#' @export
SeuratToTENET <- function(seurat_obj,
                          Species = "human",
                          DEG = list(), DAR = list(), TF_list = list(),
                          pseudotime = "pseudotime",
                          cell_select = "cell_select") {
  if (!("RNA" %in% names(seurat_obj@assays))) {
    stop("Seurat object does not have an RNA assay.")
  }
  rna_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  if ("ATAC" %in% names(seurat_obj@assays)) {
    atac_counts <- GetAssayData(seurat_obj, assay = "ATAC", layer = "counts")
  } else {
    warning("Seurat object does not have an ATAC assay; using an empty matrix.")
    atac_counts <- matrix(0, nrow = 0, ncol = ncol(rna_counts))
    colnames(atac_counts) <- colnames(rna_counts)
  }
  meta <- seurat_obj@meta.data
  cells <- rownames(meta)
  
  if (is.data.frame(pseudotime)) {
    ext_pseudo <- pseudotime
  } else {
    if (pseudotime %in% colnames(meta)) {
      ext_pseudo <- data.frame(pseudotime = meta[[pseudotime]], row.names = cells)
    } else {
      warning(paste("Pseudotime column", pseudotime, "not found in meta.data; using NA values for pseudotime."))
      ext_pseudo <- data.frame(pseudotime = rep(NA, length(cells)), row.names = cells)
    }
  }
  ext_pseudo$pseudotime <- as.numeric(ext_pseudo$pseudotime)
  
  if (is.data.frame(cell_select)) {
    if (ncol(cell_select) > 1) {
      warning("cell_select data frame has multiple columns. Only the first column will be used and renamed to 'select'.")
      cell_select <- cell_select[, 1, drop = FALSE]
    }
    colnames(cell_select) <- "select"
    ext_cell_select <- cell_select
  } else {
    if (cell_select %in% colnames(meta)) {
      ext_cell_select <- data.frame(cell_select = meta[[cell_select]], row.names = cells)
      colnames(ext_cell_select) <- "select"
    } else {
      warning(paste("cell_select column", cell_select, "not found in meta.data; using default value 1 for all cells."))
      ext_cell_select <- data.frame(select = rep(1, length(cells)), row.names = cells)
    }
  }
  ext_cell_select$select <- as.numeric(ext_cell_select$select)
  
  new_meta <- meta
  new_meta$pseudotime <- ext_pseudo$pseudotime
  new_meta$cell_select <- ext_cell_select$select
  new_meta$Species <- Species
  
  new_meta <- new_meta[order(new_meta$pseudotime, na.last = TRUE), , drop = FALSE]
  order_cells <- rownames(new_meta)
  
  rna_counts <- rna_counts[, order_cells, drop = FALSE]
  atac_counts <- atac_counts[, order_cells, drop = FALSE]
  
  metadata_df <- new_meta
  
  tenet_obj <- CreateTENET(
    rna_counts = rna_counts,
    atac_counts = atac_counts,
    pseudo = ext_pseudo,
    cell_select = ext_cell_select,
    processedData = NULL,
    DEG = DEG,
    DAR = DAR,
    TF_list = TF_list,
    Species = Species,
    metadata = metadata_df
  )
  
  return(tenet_obj)
}

