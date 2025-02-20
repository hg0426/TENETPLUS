setwd('/home/Data_Drive_8TB/chs1151/Vesalius/Renv3/')

getwd()

source('renv/activate.R')

.libPaths()
.libPaths('/home/Data_Drive_8TB/chs1151/Vesalius/Renv3/renv/library/R-4.3/x86_64-pc-linux-gnu')

# Load required library
library(Seurat)

rna <- read.table('/home/NAS/TENETPLUS/TEP_obj/test_matrix_rna.csv',sep = ',',header = T,row.names = 1,check.names = F)
atac <- read.table('/home/NAS/TENETPLUS/TEP_obj/test_matrix_atac.csv',sep = ',',header = T,row.names = 1,check.names = F)

cell_select <- read.table('/home/NAS/TENETPLUS/TEP_obj//test_cellselect.txt',sep='\t')
cell_select <- as.data.frame(cell_select,row.names = rownames(rna))
colnames(cell_select) <- 'select'

pseudotime <- read.table('/home/NAS/TENETPLUS/TEP_obj//test_trajectory.txt',sep='\t')
pseudotime <- as.data.frame(pseudotime,row.names = rownames(rna))
# colnames(pseudotime) <- 'pseudotime'

# --------------------------------------------------------------------------------
#  TENET Object Complete Code (Latest Version, Preserving Seurat Metadata)
# --------------------------------------------------------------------------------

# Load Required Packages
library(Matrix)
library(Seurat)

# 1. TENET Class Definition
setClass(
  "TENET",
  slots = list(
    input = "list",
    output = "list"
  )
)

# 2. Validity Check Function
setValidity("TENET", function(object) {
  # Check if input list contains count, metadata, feature_list
  required_input_keys <- c("count", "metadata", "feature_list")
  if (!all(required_input_keys %in% names(object@input))) {
    return(paste("input list must contain:", paste(required_input_keys, collapse = ", ")))
  }

  # Check if count list contains rna_counts and atac_counts
  count_list <- object@input$count
  if (!all(c("rna_counts", "atac_counts") %in% names(count_list))) {
    return("input$count must contain 'rna_counts' and 'atac_counts'.")
  }
  rna <- count_list$rna_counts
  atac <- count_list$atac_counts

  # Check types of rna_counts and atac_counts
  if (!(is.matrix(rna) || inherits(rna, "dgCMatrix"))) {
    return("rna_counts must be a matrix or a dgCMatrix.")
  }
  if (!(is.matrix(atac) || inherits(atac, "dgCMatrix"))) {
    return("atac_counts must be a matrix or a dgCMatrix.")
  }
  # Ensure both have same number of columns
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
  if (!("pseudotime" %in% colnames(metadata_df))) {
    return("metadata must contain a 'pseudotime' column.")
  }
  if (!("Species" %in% colnames(metadata_df))) {
    return("metadata must contain a 'Species' column.")
  }

  # Check feature_list contents
  feature_list <- object@input$feature_list
  required_feature_keys <- c("DEG", "DAR", "TF_list")
  missing_feat <- setdiff(required_feature_keys, names(feature_list))
  if (length(missing_feat) > 0) {
    for (key in missing_feat) {
      feature_list[[key]] <- list()
    }
    object@input$feature_list <- feature_list
  }

  # If processedData exists in the output, validate it
  if ("processedData" %in% names(object@output)) {
    pd <- object@output$processedData
    if (!is.data.frame(pd)) {
      return("processedData must be a data.frame.")
    }
    required_pd_cols <- c("source", "TE score", "target")
    if (!all(required_pd_cols %in% colnames(pd))) {
      return(paste("processedData must contain columns:", paste(required_pd_cols, collapse = ", ")))
    }
  }

  return(TRUE)
})

# 3. CreateTENET Function (Constructor)
CreateTENET <- function(rna_counts, atac_counts, pseudo, cell_select,
                        processedData = NULL,
                        DEG = list(), DAR = list(), TF_list = list(),
                        feature_list = NULL,
                        Species,
                        metadata = NULL) {

  # If processedData is NULL, set the output slot to an empty list
  if (is.null(processedData)) {
    output_list <- list()
  } else {
    output_list <- list(processedData = processedData)
  }

  # If metadata is not provided, build it from pseudo and cell_select
  if (!is.null(metadata)) {
    metadata_df <- metadata
  } else {
    if (!is.data.frame(pseudo)) {
      stop("pseudo must be a data.frame.")
    }
    if (ncol(pseudo) < 1) {
      stop("pseudo must have at least one column.")
    }
    # If 'pseudotime' column does not exist, rename the first column
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
    pseudo <- pseudo[order(pseudo$pseudotime), , drop = FALSE]
    order_cells <- rownames(pseudo)

    # cell_select data frame must have only one column, rename it to "select"
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

  # Handle DEG, DAR, TF_list
  # If they are data frames with multiple columns, only use the first column
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

  # If feature_list is NULL, build one; otherwise add missing keys
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

# 4. Define the show Method for TENET Objects
setMethod("show", "TENET", function(object) {
  cat("   _________   ______    ___    ___   ______   _________              \n")
  cat("  /___  ___/  / ____/   /   |  /  /  / ____/  /___  ___/   _          \n")
  cat("     / /     / __      /  | | /  /  / __        / /      _| |_        \n")
  cat("    / /     / __/     /   | |/  /  /  _/       / /      |_   _|       \n")
  cat("   / /     / ___     /  / |    /  / ___       / /         |_|         \n")
  cat("  /_/     /_____/   /__/  |_._/  /_____/     /_/                      \n")
  cat("\n")
  cat("----- TENET Object Summary -----\n")

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

  other_cols <- setdiff(colnames(meta), c("pseudotime", "Species"))
  if (length(other_cols) > 0) {
    cat("Additional metadata columns:", paste(other_cols, collapse = ", "), "\n")
  } else {
    cat("No additional metadata columns.\n")
  }

  cat("Species:", species, "\n")

  # Display DEG entries
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

  # Display DAR entries
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

  # Display TF_list entries
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

  # Display output information
  if ("processedData" %in% names(object@output)) {
    pd <- object@output$processedData
    cat("Processed Data:", nrow(pd), "rows,", ncol(pd), "columns\n")

  } else if ("result" %in% names(object@output)) {
    res <- object@output$result

    # Check for dgCMatrix
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

      # Handle dgCMatrix inside the output as well
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

# 5. Convert a Seurat Object to a TENET Object
SeuratToTENET <- function(seurat_obj,
                          Species = "human",
                          DEG = list(), DAR = list(), TF_list = list(),
                          pseudotime = "pseudotime",
                          cell_select = "cell_select") {

  # Check if the RNA assay is present
  if (!("RNA" %in% names(seurat_obj@assays))) {
    stop("Seurat object does not have an RNA assay.")
  }
  rna_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

  # Handle ATAC assay (if missing, create an empty matrix)
  if ("ATAC" %in% names(seurat_obj@assays)) {
    atac_counts <- GetAssayData(seurat_obj, assay = "ATAC", layer = "counts")
  } else {
    warning("Seurat object does not have an ATAC assay; using an empty matrix.")
    atac_counts <- matrix(0, nrow = 0, ncol = ncol(rna_counts))
    colnames(atac_counts) <- colnames(rna_counts)
  }

  meta <- seurat_obj@meta.data
  cells <- rownames(meta)

  # Handle pseudotime
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

  # Handle cell_select
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

  # Build new metadata: preserve all columns of meta, override pseudotime, cell_select, Species
  new_meta <- meta
  new_meta$pseudotime <- ext_pseudo$pseudotime
  new_meta$cell_select <- ext_cell_select$select
  new_meta$Species <- Species

  # Reorder cells based on pseudotime
  new_meta <- new_meta[order(new_meta$pseudotime, na.last = TRUE), , drop = FALSE]
  order_cells <- rownames(new_meta)

  # Reorder RNA and ATAC counts accordingly
  rna_counts <- rna_counts[, order_cells, drop = FALSE]
  atac_counts <- atac_counts[, order_cells, drop = FALSE]

  # Use the fully updated metadata
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


tenetObj2 <- CreateTENET(
  rna_counts   = t(rna),
  atac_counts  = t(atac),
  pseudo       = pseudotime,
  cell_select  = cell_select,
  # processedData = proc_df,  # 선택적: 생략 가능
  # DEG          = deg_obj,    # 값이 없으면 list()로 넘기면 됨
  # DAR          = dar_obj,    # 값이 없으면 list()로 넘기면 됨
  # TF_list      = tf_obj,     # 값이 없으면 list()로 넘기면 됨
  Species      = "human"
)

tenetObj2@input$count$rna_counts[1:3,1:3]

tenetObj
tenetObj@input$metadata

########################################
#  runTENETPlus : TENET+ 파이프라인 실행 (개선 버전)
########################################
#' runTENETPlus
#'
#' This function takes a TENET object (with RNA and ATAC data) and runs the
#' TENET+ pipeline by writing out the necessary input files, calling the
#' bash+python scripts, and reading all generated results back into the
#' \code{TENET} object.
#'
#' \describe{
#'   \item{\code{output$TE_mat}}{Stores \code{TE_result_matrix.txt} as a sparse matrix}
#'   \item{\code{output$TF_TG}}{From \code{TE_result_matrix_rowTF_colGN.txt}, pivoted into (source, score, target)}
#'   \item{\code{output$TF_TAR}}{From \code{TE_result_matrix_rowTF_colPK.txt}, pivoted into (source, score, target)}
#'   \item{\code{output$TAR_TG}}{From \code{TE_result_matrix_rowPeak_colGN.txt}, pivoted into (source, score, target)}
#' }
#'
#' @param obj A TENET object
#' @param outdir A directory to store intermediate / result files
#' @param ncores Number of parallel jobs
#' @param historyLength (k) param
#' @param species Force using this species (if empty, use obj@input$metadata$Species)
#' @param mode 0,1,2,... (TENET_Plus_for_py.sh에서 정의하는 모드)
#' @param bashFile If NULL, use package-included script; otherwise use user-provided path
#'
#' @return The TENET object, with several new fields in \code{obj@output}.
#'
#' @importFrom utils write.csv write.table read.table read.csv
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix
runTENETPlus <- function(obj,
                         outdir = ".",
                         ncores = 2,
                         historyLength = 1,
                         species = NULL,
                         mode = 1,
                         bashFile = NULL) {
  #################################
  # 1) Basic checks
  #################################
  if (!inherits(obj, "TENET")) {
    stop("obj must be a TENET object.")
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(bashFile)) {
    stop("Please provide 'bashFile' that points to TENET_Plus_for_py.sh.")
  }
  if (!file.exists(bashFile)) {
    stop("The specified bashFile does not exist: ", bashFile)
  }

  #################################
  # 2) Prepare data
  #################################
  rna_mat <- obj@input$count$rna_counts
  atac_mat <- obj@input$count$atac_counts
  meta <- obj@input$metadata
  feat <- obj@input$feature_list

  # Determine species
  if (is.null(species) || !nzchar(species)) {
    sp_all <- unique(meta$Species)
    if (length(sp_all) > 1) {
      message("Multiple species found. Using the first: ", sp_all[1])
      sp_all <- sp_all[1]
    }
    sp <- sp_all
  } else {
    sp <- species
  }

  # 3) Optional filter by DEG/DAR
  deg_vec <- feat$DEG
  dar_vec <- feat$DAR
  gene_names <- rownames(rna_mat)
  peak_names <- rownames(atac_mat)

  use_genes <- gene_names
  use_peaks <- peak_names
  union_deg_dar <- unique(c(deg_vec, dar_vec))
  if (length(union_deg_dar) > 0) {
    use_genes <- intersect(gene_names, union_deg_dar)
    use_peaks <- intersect(peak_names, union_deg_dar)
    if (length(use_genes) == 0 && length(use_peaks) == 0) {
      warning("Filtering leaves no features. Using all features instead.")
      use_genes <- gene_names
      use_peaks <- peak_names
    }
  }

  rna_mat_sub <- rna_mat[use_genes, , drop = FALSE]
  atac_mat_sub <- atac_mat[use_peaks, , drop = FALSE]

  #################################
  # 4) Combine RNA + ATAC into one (#cells x #features) matrix
  #################################
  mat_rna_t  <- t(rna_mat_sub)   # cells x genes
  mat_atac_t <- t(atac_mat_sub)  # cells x peaks
  combined   <- cbind(mat_rna_t, mat_atac_t)

  # Rename peak columns to have "chr" prefix if missing
  gene_names_sub <- rownames(rna_mat_sub)
  peak_names_sub <- rownames(atac_mat_sub)
  chr_peaks <- vapply(seq_along(peak_names_sub), function(i) {
    pk <- peak_names_sub[i]
    if (!grepl("^chr", pk)) {
      paste0("chrX-", i*100, "-", i*100+50)
    } else {
      pk
    }
  }, FUN.VALUE = character(1))
  colnames(combined) <- c(gene_names_sub, chr_peaks)

  # Sanity check
  if (!identical(rownames(combined), rownames(meta))) {
    stop("Row names of combined data do not match rownames of metadata.")
  }

  #################################
  # 5) Write the input files to outdir
  #################################
  matrix_csv_file <- file.path(outdir, "tenet_plus_input.csv")
  df_for_write <- as.data.frame(combined, check.names = FALSE)
  write.csv(df_for_write, matrix_csv_file, row.names = TRUE, quote = FALSE)

  pseudotime_file <- file.path(outdir, "pseudotime.txt")
  write.table(meta$pseudotime, pseudotime_file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  cell_select_file <- file.path(outdir, "cell_select.txt")
  if ("cell_select" %in% colnames(meta)) {
    write.table(meta$cell_select, cell_select_file,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else {
    write.table(rep(1, nrow(meta)), cell_select_file,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  #################################
  # 6) Construct the command
  #################################
  scriptDir <- dirname(bashFile)
  scriptName <- basename(bashFile)
  cmd <- sprintf(
    "cd %s && ./%s %s %d %s %s %d %s %d",
    shQuote(scriptDir),
    shQuote(scriptName),
    shQuote(normalizePath(matrix_csv_file)),
    ncores,
    shQuote(normalizePath(pseudotime_file)),
    shQuote(normalizePath(cell_select_file)),
    historyLength,
    shQuote(sp),
    mode
  )
  message("Running TENET+ command:\n", cmd)

  #################################
  # 7) Run system command
  #################################
  exit_code <- system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  message("TENET+ exit code: ", exit_code)

  #################################
  # Helper functions
  #################################

  # A) Read TE_result_matrix.* as a sparse matrix
  readAsSparseMatrix <- function(filePath) {
    # read.table(..., row.names=1) => interpret first col as rownames
    df <- tryCatch(
      {
        read.table(filePath, header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      },
      error = function(e) {
        message("Failed to read file as matrix: ", filePath)
        return(NULL)
      }
    )
    if (is.null(df) || ncol(df) == 0) {
      return(NULL)
    }
    mat <- as.matrix(df)
    # Convert to sparse
    mat_sp <- methods::as(mat, "sparseMatrix")
    return(mat_sp)
  }

  # B) Pivot splitted matrix files into (source, score, target)
  pivotMatrixFileToDF <- function(filePath) {
    df <- tryCatch(
      {
        read.table(filePath, header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      },
      error = function(e) {
        message("Failed to read file as splitted matrix: ", filePath)
        return(NULL)
      }
    )
    if (is.null(df) || ncol(df) == 0) {
      return(NULL)
    }
    mat <- as.matrix(df)
    row_ids <- rownames(mat)
    col_ids <- colnames(mat)
    # pivot
    out_list <- vector("list", length = nrow(mat)*ncol(mat))
    idx <- 1
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        # (source, score, target)
        out_list[[idx]] <- c(row_ids[i], as.character(mat[i, j]), col_ids[j])
        idx <- idx + 1
      }
    }
    outdf <- data.frame(do.call(rbind, out_list), stringsAsFactors = FALSE)
    colnames(outdf) <- c("source", "score", "target")
    # Convert score to numeric if possible
    outdf$score <- suppressWarnings(as.numeric(outdf$score))
    return(outdf)
  }

  #################################
  # 8) Gather outputs
  #################################
  # We'll store them in obj@output

  obj@output$TE_mat  <- NULL
  obj@output$TF_TG   <- NULL
  obj@output$TF_TAR  <- NULL
  obj@output$TAR_TG  <- NULL

  # Move/copy relevant output files from scriptDir to outdir
  possible_files_dir <- unique(c(scriptDir, outdir))
  for (d in possible_files_dir) {
    all_files <- list.files(d, pattern = "^TE_result_matrix", full.names = TRUE)
    for (f in all_files) {
      fname <- basename(f)
      dest_f <- file.path(outdir, fname)
      if (!file.exists(dest_f)) {
        file.copy(f, dest_f, overwrite = TRUE)
      }
    }
  }

  # (1) TE_result_matrix.txt => obj@output$TE_mat (sparse matrix)
  mainMatFile <- file.path(outdir, "TE_result_matrix.txt")
  if (file.exists(mainMatFile)) {
    sparse_mat <- readAsSparseMatrix(mainMatFile)
    if (!is.null(sparse_mat)) {
      obj@output$TE_mat <- sparse_mat
    }
  }

  # (2) TE_result_matrix_rowTF_colGN.txt => obj@output$TF_TG (df: source, score, target)
  tfTgFile <- file.path(outdir, "TE_result_matrix_rowTF_colGN.txt")
  if (file.exists(tfTgFile)) {
    tfTgDF <- pivotMatrixFileToDF(tfTgFile)
    if (!is.null(tfTgDF)) {
      obj@output$TF_TG <- tfTgDF
    }
  }

  # (3) TE_result_matrix_rowTF_colPK.txt => obj@output$TF_TAR (df: source, score, target)
  tfTarFile <- file.path(outdir, "TE_result_matrix_rowTF_colPK.txt")
  if (file.exists(tfTarFile)) {
    tfTarDF <- pivotMatrixFileToDF(tfTarFile)
    if (!is.null(tfTarDF)) {
      obj@output$TF_TAR <- tfTarDF
    }
  }

  # (4) TE_result_matrix_rowPeak_colGN.txt => obj@output$TAR_TG (df: source, score, target)
  tarTgFile <- file.path(outdir, "TE_result_matrix_rowPeak_colGN.txt")
  if (file.exists(tarTgFile)) {
    tarTgDF <- pivotMatrixFileToDF(tarTgFile)
    if (!is.null(tarTgDF)) {
      obj@output$TAR_TG <- tarTgDF
    }
  }

  invisible(obj)
}



#setwd('/home/Data_Drive_8TB/chs1151/TEP_obj/')
library(reticulate)
#options(reticulate.conda_binary = "/usr/local/anaconda3")
use_python('/home/Data_Drive_8TB/hg0426/.conda/envs/TENETPLUS/bin/python',required=TRUE)
py_config()


# (2) TENET+ 실행
tenetObj <- runTENETPlus(tenetObj,
                         outdir = "myTENETout",
                         ncores = 32,
                         historyLength = 1,
                         species = 'human',   # or "human"
                         mode = 1,
                         bashFile = '/home/Data_Drive_8TB_2/TENETPLUS/TENETPLUSpack/TENETPLUS/TENET_Plus_for_py')  # or "path/to/custom.sh"

str(tenetObj@output$TE_mat)
#tenetObj@output$TE_mat <- NULL
tenetObj

str(tenetObj@output)

tenetObj <- CreateTENET(
  rna_counts   = t(rna),
  #atac_counts  = t(atac),
  pseudo       = pseudotime,
  cell_select  = cell_select,
  # processedData = proc_df,  # 선택적: 생략 가능
  # DEG          = deg_obj,    # 값이 없으면 list()로 넘기면 됨
  # DAR          = dar_obj,    # 값이 없으면 list()로 넘기면 됨
  # TF_list      = tf_obj,     # 값이 없으면 list()로 넘기면 됨
  Species      = "human"
)
tenetObj

# (2) TENET+ 실행
tenetObj <- runTENETPlus(tenetObj,
                         outdir = "myTENETout",
                         ncores = 32,
                         historyLength = 1,
                         species = 'human',   # or "human"
                         mode = 1,
                         bashFile = '/home/Data_Drive_8TB_2/TENETPLUS/TENETPLUSpack/TENETPLUS/TENET_Plus_for_py')






saveRDS(tenetObj,'TE_pack_result.RDS')

