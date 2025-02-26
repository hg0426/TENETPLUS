CreateTENET <- function(rna_counts, atac_counts = NULL, pseudo, cell_select,
                        processedData = NULL,
                        DEG = list(), DAR = list(), TF_list = list(),
                        feature_list = NULL,
                        Species,
                        metadata = NULL,
                        ident = "orig.ident") {
  if (is.null(atac_counts)) {
    atac_counts <- matrix(0, nrow = 0, ncol = ncol(rna_counts))
    if (!is.null(colnames(rna_counts))) {
      colnames(atac_counts) <- colnames(rna_counts)
    }
  }
  
  if (is.null(processedData)) {
    output_list <- list()
  } else {
    output_list <- list(processedData = processedData)
  }
  
  if (!is.null(metadata)) {
    metadata_df <- metadata
  } else {
    if (!is.data.frame(pseudo)) {
      stop("pseudo must be a data.frame.")
    }
    if (ncol(pseudo) < 1) {
      stop("pseudo must have at least one column.")
    }
    if (!("pseudotime" %in% colnames(pseudo))) {
      if (ncol(pseudo) > 1) {
        warning("pseudo data.frame has multiple columns but no 'pseudotime' column. Using the first column as pseudotime.")
      }
      pseudo <- pseudo[, 1, drop = FALSE]
      colnames(pseudo)[1] <- "pseudotime"
    }
    if (is.null(rownames(pseudo))) {
      stop("pseudo must have row names corresponding to cell names.")
    }
    pseudo <- as.data.frame(pseudo)
    rownames(pseudo) <- rownames(metadata_df)
    pseudo <- pseudo[colnames(rna_counts), , drop = FALSE]
    
    pseudo <- pseudo[order(pseudo$pseudotime), , drop = FALSE]
    order_cells <- rownames(pseudo)
    
    if (!is.data.frame(cell_select)) {
      stop("cell_select must be a data.frame.")
    }
    if (ncol(cell_select) > 1) {
      warning("cell_select data frame has multiple columns. Only the first column will be used and renamed to 'select'.")
      cell_select <- cell_select[, 1, drop = FALSE]
    }
    colnames(cell_select) <- "select"
    cell_select <- cell_select[order_cells, , drop = FALSE]
    
    metadata_df <- cbind(pseudo, cell_select)
    metadata_df$Species <- Species
  }
  
  meta_rna <- metadata_df[metadata_df$ident == "RNA", ]
  meta_atac <- metadata_df[metadata_df$ident == "ATAC", ]
  rna_order_cells <- rownames(meta_rna)
  atac_order_cells <- rownames(meta_atac)

  rna_counts <- rna_counts[, rna_order_cells, drop = FALSE]
  atac_counts <- atac_counts[, atac_order_cells, drop = FALSE]
  
  if (is.data.frame(DEG)) {
    if (ncol(DEG) > 1) {
      warning("DEG data frame has multiple columns. Only the first column will be used and renamed to 'DEG'.")
      DEG <- DEG[, 1, drop = FALSE]
    }
    colnames(DEG) <- "DEG"
  }
  if (is.data.frame(DAR)) {
    if (ncol(DAR) > 1) {
      warning("DAR data frame has multiple columns. Only the first column will be used and renamed to 'DAR'.")
      DAR <- DAR[, 1, drop = FALSE]
    }
    colnames(DAR) <- "DAR"
  }
  if (is.data.frame(TF_list)) {
    if (ncol(TF_list) > 1) {
      warning("TF_list data frame has multiple columns. Only the first column will be used and renamed to 'TF_list'.")
      TF_list <- TF_list[, 1, drop = FALSE]
    }
    colnames(TF_list) <- "TF_list"
  }
  
  if (is.null(feature_list)) {
    feature_list <- list(DEG = DEG, DAR = DAR, TF_list = TF_list)
  } else {
    if (!("DEG" %in% names(feature_list))) {
      feature_list$DEG <- DEG
    }
    if (!("DAR" %in% names(feature_list))) {
      feature_list$DAR <- DAR
    }
    if (!("TF_list" %in% names(feature_list))) {
      feature_list$TF_list <- TF_list
    }
  }
  
  input_list <- list(
    count = list(
      rna_counts = rna_counts,
      atac_counts = atac_counts
    ),
    metadata = metadata_df,
    feature_list = feature_list
  )
  
  new("TENET", input = input_list, output = output_list)
}

