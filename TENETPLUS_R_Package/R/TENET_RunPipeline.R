#' Run TENET+ Pipeline
#'
#' Executes the TENET+ pipeline by writing out the necessary input files,
#' invoking an external bash script (and Python scripts), and reading the results
#' back into the TENET object.
#'
#' @param obj A TENET object.
#' @param outdir Directory to store intermediate and result files.
#' @param ncores Number of parallel jobs.
#' @param historyLength Parameter for history length.
#' @param species Species to use (if NULL, uses obj@input$metadata$Species).
#' @param mode Mode to run (as defined in TENET_Plus_for_py.sh).
#' @param bashFile Path to the bash script. If using an internal script, retrieve its path via system.file().
#'
#' @return Updated TENET object with output stored in obj@output.
#' @export

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
    bashFile <- system.file("bash", "TENET_Plus_for_py", package = "TENETPLUS")
    if (file.exists(bashFile)) {
      # Set executable permission
      system(paste("chmod +x", shQuote(bashFile)))
    } else {
      warning("TENET_Plus_for_py script not found in package.")
    }
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
  if ("select" %in% colnames(meta)) {
    write.table(meta$select, cell_select_file,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else {
    write.table(rep(1, nrow(meta)), cell_select_file,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  #################################
  # 6) Copy bashFile, dependent python scripts, and reference files to outdir
  #################################
  # Directory containing TENET_Plus_for_py.sh (copy necessary files from here)
  source_dir <- dirname(bashFile)

  # Copy bashFile
  bashFile_new <- file.path(outdir, basename(bashFile))
  if (!file.exists(bashFile_new)) {
    file.copy(bashFile, bashFile_new, overwrite = TRUE)
  }

  # Copy memory_check.sh if it exists
  mem_script <- file.path(source_dir, "memory_check.sh")
  if (file.exists(mem_script)) {
    mem_script_new <- file.path(outdir, "memory_check.sh")
    if (!file.exists(mem_script_new)) {
      file.copy(mem_script, mem_script_new, overwrite = TRUE)
    }
  }

  # List of necessary python scripts (called from the bash script using relative paths)
  python_scripts <- c("process_matrix.py",
                      "PreProcessScript_TE_Plus.py",
                      "PreProcessScript_TENET_TF.py",
                      "runTE_for_py_CPU.py",
                      "makeTEasMatrix.py",
                      "trim_indirect.py",
                      "make_GRN_new.py")
  for (ps in python_scripts) {
    src_ps <- file.path(source_dir, ps)
    dest_ps <- file.path(outdir, ps)
    if (file.exists(src_ps)) {
      file.copy(src_ps, dest_ps, overwrite = TRUE)
    } else {
      message("Warning: ", src_ps, " not found. Please ensure this file exists.")
    }
  }

  # List of necessary reference files (example)
  ref_files <- c("gene_chr.txt",
                 "GO_symbol_human_regulation_of_transcription+sequence-specific_DNA_binding_list.txt",
                 "infodynamics.jar",
                 "GO_symbol_mouse_regulation_of_transcription+sequence-specific_DNA_binding_list.txt")
  for (rf in ref_files) {
    src_rf <- file.path(source_dir, rf)
    dest_rf <- file.path(outdir, rf)
    if (file.exists(src_rf)) {
      file.copy(src_rf, dest_rf, overwrite = TRUE)
    } else {
      message("Warning: Reference file ", src_rf, " not found. Please ensure this file exists.")
    }
  }

  #################################
  # 7) Construct the command with working directory set to outdir
  #################################
  # For Linux, add 'stdbuf -oL' to force line buffering so that output is displayed in real-time.
  cmd <- sprintf(
    "cd %s && stdbuf -oL ./%s %s %d %s %s %d %s %d 2>&1",
    shQuote(normalizePath(outdir)),
    shQuote(basename(bashFile_new)),
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
  # 8) Run system command and display output in real-time
  #################################
  # Using intern = FALSE streams the output directly to the console in real-time.
  exit_code <- system(cmd, intern = FALSE)
  message("TENET+ exit code: ", exit_code)

  #################################
  # Helper functions
  #################################

  # A) Read TE_result_matrix.* as a sparse matrix
  readAsSparseMatrix <- function(filePath) {
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
    mat_sp <- methods::as(mat, "sparseMatrix")
    return(mat_sp)
  }

  # B) New function to read SIF files as data frames with source, score, target columns
  readSifToDF <- function(filePath) {
    df <- tryCatch(
      {
        # Read the SIF file assuming tab-separated values without header
        read.table(filePath, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      },
      error = function(e) {
        message("Failed to read SIF file: ", filePath)
        return(NULL)
      }
    )
    # Check if the file has at least 3 columns
    if (!is.null(df) && ncol(df) >= 3) {
      colnames(df)[1:3] <- c("source", "score", "target")
      df$score <- suppressWarnings(as.numeric(df$score))
    } else {
      message("File does not have the expected format: ", filePath)
      return(NULL)
    }
    return(df)
  }

  #################################
  # 9) Gather outputs -> obj@output$result
  #################################
  all_results <- list(
    TE_mat  = NULL,
    TF_TG   = NULL,
    TF_TAR  = NULL,
    TAR_TG  = NULL,
    TF_TG_indirect  = NULL
  )

  # Copy result files from source_dir and outdir if necessary
  possible_files_dir <- unique(c(source_dir, outdir))
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

  # (1) TE_result_matrix.txt => TE_mat (sparse matrix)
  mainMatFile <- file.path(outdir, "TE_result_matrix.txt")
  if (file.exists(mainMatFile)) {
    sparse_mat <- readAsSparseMatrix(mainMatFile)
    if (!is.null(sparse_mat)) {
      all_results$TE_mat <- sparse_mat
    }
  }

  # (2) TE_result_matrix_rowTF_colGN.sif => TF_TG
  tfTgFile <- file.path(outdir, "TE_result_matrix_rowTF_colGN.fdr0.01.sif")
  if (file.exists(tfTgFile)) {
    tfTgDF <- readSifToDF(tfTgFile)
    if (!is.null(tfTgDF)) {
      all_results$TF_TG <- tfTgDF
    }
  }

  # (3) TE_result_matrix_rowTF_colPK.sif => TF_TAR
  tfTarFile <- file.path(outdir, "TE_result_matrix_rowTF_colPK.fdr0.01.sif")
  if (file.exists(tfTarFile)) {
    tfTarDF <- readSifToDF(tfTarFile)
    if (!is.null(tfTarDF)) {
      all_results$TF_TAR <- tfTarDF
    }
  }

  # (4) TE_result_matrix_rowPeak_colGN.sif => TAR_TG
  tarTgFile <- file.path(outdir, "TE_result_matrix_rowPeak_colGN.fdr0.01.sif")
  if (file.exists(tarTgFile)) {
    tarTgDF <- readSifToDF(tarTgFile)
    if (!is.null(tarTgDF)) {
      all_results$TAR_TG <- tarTgDF
    }
  }

  # (5) TE_result_matrix_rowTF_colGN.sif => TF_TG_indirect
  tfTgFile <- file.path(outdir, "TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect0.0.sif")
  if (file.exists(tfTgFile)) {
    tfTgDF <- readSifToDF(tfTgFile)
    if (!is.null(tfTgDF)) {
      all_results$TF_TG_indirect <- tfTgDF
    }
  }

  obj@output$result <- all_results

  invisible(obj)
}
