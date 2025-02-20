#' Show Method for TENET Objects
#'
#' Custom display of a TENET object in the R console.
#'
#' @param object A TENET object.
#' @export
setMethod("show", "TENET", function(object) {
 # cat("   _________   ______    ___    ___   ______   _________              \n")
 # cat("  /___  ___/  / ____/   /   |  /  /  / ____/  /___  ___/   _          \n")
 # cat("     / /     / __      /  | | /  /  / __        / /      _| |_        \n")
 # cat("    / /     / __/     /   | |/  /  /  _/       / /      |_   _|       \n")
 # cat("   / /     / ___     /  / |    /  / ___       / /         |_|         \n")
 # cat("  /_/     /_____/   /__/  |_._/  /_____/     /_/                      \n")
 # cat("\n")
 # cat("----- TENET Object Summary -----\n")
  
  rna <- object@input$count$rna_counts
  atac <- object@input$count$atac_counts
  meta <- object@input$metadata
  species <- meta$Species[1]
  
  cell_count <- ncol(rna)
  gene_count <- nrow(rna)
  peak_count <- nrow(atac)
  
  cat("Cells:", cell_count, "\n")
  cat("Genes (RNA):", gene_count, "\n")
  cat("Peaks (ATAC):", peak_count, "\n")
  cat("Pseudo entries:", nrow(meta), "\n")
  selected_count <- sum(meta$select == 1, na.rm = TRUE)
  cat("Pseudo entries (selected cells):", selected_count, "\n")
  other_cols <- setdiff(colnames(meta), c("pseudotime", "Species"))
  if (length(other_cols) > 0) {
    cat("Additional metadata columns:", paste(other_cols, collapse = ", "), "\n")
  } else {
    cat("No additional metadata columns.\n")
  }
  
  cat("Species:", species, "\n")
  
  deg <- object@input$feature_list$DEG
  if (length(deg) > 0) {
    if (is.list(deg)) {
      deg_keys <- names(deg)
      if (is.null(deg_keys) || all(deg_keys == "")) {
        deg_keys <- rep("", length(deg))
      }
      deg_counts <- sapply(deg, length)
      cat("Number of DEG entries: [", paste(deg_counts, collapse = ", "), "]\n")
    } else if (is.vector(deg)) {
      cat("Number of DEG entries:", length(deg), "\n")
      cat(" DEG: DEG\n")
    }
  } else {
    cat("No DEG entries found.\n")
  }
  
  dar <- object@input$feature_list$DAR
  if (length(dar) > 0) {
    if (is.list(dar)) {
      dar_keys <- names(dar)
      if (is.null(dar_keys) || all(dar_keys == "")) {
        dar_keys <- rep("", length(dar))
      }
      dar_counts <- sapply(dar, length)
      cat("Number of DAR entries: [", paste(dar_counts, collapse = ", "), "]\n")
    } else if (is.vector(dar)) {
      cat("Number of DAR entries:", length(dar), "\n")
      cat(" DAR: DAR\n")
    }
  } else {
    cat("No DAR entries found.\n")
  }
  
  tf <- object@input$feature_list$TF_list
  if (!is.null(tf) && length(tf) > 0) {
    if (is.list(tf)) {
      tf_keys <- names(tf)
      if (is.null(tf_keys) || all(tf_keys == "")) {
        tf_keys <- rep("", length(tf))
      }
      tf_counts <- sapply(tf, length)
      cat("Number of TF entries: [", paste(tf_counts, collapse = ", "), "]\n")
    } else if (is.vector(tf)) {
      cat("Number of TF entries:", length(tf), "\n")
      cat(" TF:", paste(tf, collapse = ", "), "\n")
    }
  }
  
  if ("processedData" %in% names(object@output)) {
    pd <- object@output$processedData
    cat("Processed Data:", nrow(pd), "rows,", ncol(pd), "columns\n")
  } else if ("result" %in% names(object@output)) {
    res <- object@output$result
    if (inherits(res, "dgCMatrix")) {
      cat("Result: dgCMatrix with", nrow(res), "rows,", ncol(res), "columns\n")
    } else if (is.data.frame(res) || is.matrix(res)) {
      cat("Result:", nrow(res), "rows,", ncol(res), "columns\n")
    } else if (is.list(res)) {
      cat("Result contains multiple elements:\n")
      for (name in names(res)) {
        elem <- res[[name]]
        if (inherits(elem, "dgCMatrix")) {
          cat(" -", name, ": dgCMatrix with", nrow(elem), "rows,", ncol(elem), "columns\n")
        } else if (is.data.frame(elem) || is.matrix(elem)) {
          cat(" -", name, ":", nrow(elem), "rows,", ncol(elem), "columns\n")
        } else {
          cat(" -", name, ":", elem, "\n")
        }
      }
    } else {
      cat("Result:", res, "\n")
    }
  } else if (length(object@output) > 0) {
    cat("Output:\n")
    for (name in names(object@output)) {
      out <- object@output[[name]]
      if (inherits(out, "dgCMatrix")) {
        cat(name, ": dgCMatrix with", nrow(out), "rows,", ncol(out), "columns\n")
      } else if (is.data.frame(out) || is.matrix(out)) {
        cat(name, ":", nrow(out), "rows,", ncol(out), "columns\n")
      } else {
        cat(name, ":", out, "\n")
      }
    }
  }
  cat("----------------------------------------------\n")
})

