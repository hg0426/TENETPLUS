# Load necessary packages
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Function: Load Seurat object from file
load_seurat_object <- function(file_path) {
  seurat_obj <- readRDS(file_path)
  return(seurat_obj)
}

# Function: Aggregate expression data (RNA or ATAC)
aggregate_expression_data <- function(seurat_obj, group_by = "seurat_annotations_large", assay = "RNA") {
  # Check that the group_by column exists in meta.data
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("The specified group_by column does not exist in the Seurat object's meta.data.")
  }
  
  # Get unique groups and the counts matrix
  unique_groups <- unique(seurat_obj@meta.data[[group_by]])
  expression_data <- seurat_obj@assays[[assay]]$counts
  
  
  # Map each cell to its group
  group_indices <- split(seq_len(ncol(expression_data)), seurat_obj@meta.data[[group_by]])
  
  # Efficient rowSums computation for each group
  aggregated_expression <- sapply(group_indices, function(cell_indices) {
    Matrix::rowSums(expression_data[, cell_indices, drop = FALSE])
  })
  
  # Convert to data frame and set column names
  aggregated_expression_df <- as.data.frame(as.matrix(aggregated_expression))
  colnames(aggregated_expression_df) <- unique_groups
  
  # Normalize or apply TF-IDF based on assay type and unique_groups length
  if (assay == "RNA") {
    normalized_data <- as.data.frame(NormalizeData(aggregated_expression_df))
  } else if (assay == "ATAC") {
    if (length(unique_groups) == 1) {
      normalized_data <- as.data.frame(NormalizeData(aggregated_expression_df))
    } else {
      normalized_data <- as.data.frame(RunTFIDF(aggregated_expression_df))
    }
  } else {
    stop("Unknown assay type. Assay must be either 'RNA' or 'ATAC'.")
  }
  
  # Add feature names as a column
  normalized_data$Feature <- rownames(normalized_data)
  
  return(normalized_data)
}
process_triplet <- function(type, triplet_list) {
  # Read triplet data file
  # triplet_file <- paste0(base_path, type, "/trimmed_by_PeaksourceDistance1M_Indirect-0.01_", type, "_TE_result_matrix.sif")
  # triplet_list <- read.delim(triplet_file, header = TRUE, sep = " ", stringsAsFactors = FALSE)
  triplet_list <- triplet_list
  # Collect lists of nodes
  TFs <- unique(triplet_list$TF)
  peaks <- unique(triplet_list$peak)
  genes <- unique(triplet_list$gene)
  
  # Construct edge lists
  edge_list_TF_peak <- data.frame(Source = triplet_list$TF, Target = triplet_list$peak, stringsAsFactors = FALSE)
  edge_list_TF_gene <- data.frame(Source = triplet_list$TF, Target = triplet_list$gene, stringsAsFactors = FALSE)
  edge_list_peak_gene <- data.frame(Source = triplet_list$peak, Target = triplet_list$gene, stringsAsFactors = FALSE)
  
  # Compute TF OutDegree (TF-peak and TF-gene)
  tf_peak_counts <- aggregate(Target ~ Source, data = edge_list_TF_peak, FUN = function(x) length(unique(x)))
  colnames(tf_peak_counts) <- c("Node", "TF_Peak_Count")
  
  tf_gene_counts <- aggregate(Target ~ Source, data = edge_list_TF_gene, FUN = function(x) length(unique(x)))
  colnames(tf_gene_counts) <- c("Node", "TF_Gene_Count")
  
  tf_degree <- full_join(tf_peak_counts, tf_gene_counts, by = "Node")
  tf_degree[is.na(tf_degree)] <- 0
  tf_degree$TF_OutDegree <- tf_degree$TF_Peak_Count + tf_degree$TF_Gene_Count
  
  # Compute gene InDegree (TF-gene and peak-gene)
  gene_tf_counts <- aggregate(Source ~ Target, data = edge_list_TF_gene, FUN = function(x) length(unique(x)))
  colnames(gene_tf_counts) <- c("Node", "Gene_TF_Count")
  
  gene_peak_counts <- aggregate(Source ~ Target, data = edge_list_peak_gene, FUN = function(x) length(unique(x)))
  colnames(gene_peak_counts) <- c("Node", "Gene_Peak_Count")
  
  gene_degree <- full_join(gene_tf_counts, gene_peak_counts, by = "Node")
  gene_degree[is.na(gene_degree)] <- 0
  gene_degree$Gene_InDegree <- gene_degree$Gene_TF_Count + gene_degree$Gene_Peak_Count
  
  # Compute peak OutDegree (peak-gene) and InDegree (TF-peak)
  peak_out_degree <- aggregate(Target ~ Source, data = edge_list_peak_gene, FUN = function(x) length(unique(x)))
  colnames(peak_out_degree) <- c("Node", "Peak_OutDegree")
  
  peak_in_degree <- aggregate(Source ~ Target, data = edge_list_TF_peak, FUN = function(x) length(unique(x)))
  colnames(peak_in_degree) <- c("Node", "Peak_InDegree")
  
  # Merge all degrees
  degree_counts <- full_join(tf_degree[, c("Node", "TF_OutDegree")], gene_degree[, c("Node", "Gene_InDegree")], by = "Node")
  degree_counts <- full_join(degree_counts, peak_out_degree, by = "Node")
  degree_counts <- full_join(degree_counts, peak_in_degree, by = "Node")
  
  # Replace NA with 0
  degree_counts[is.na(degree_counts)] <- 0
  
  # Determine node types
  degree_counts$NodeType <- ifelse(degree_counts$Node %in% TFs, "TF",
                                   ifelse(degree_counts$Node %in% genes, "gene",
                                          ifelse(degree_counts$Node %in% peaks, "peak", "unknown")))
  
  # Add cell type column
  degree_counts$CellType <- type
  
  return(degree_counts)
}


# Function: Create degree list for multiple cell types
get_degree_list <- function(cell_types, triplet_list) {
  degree_list <- list()
  node_set <- character()
  
  for (type in cell_types) {
    degree_counts <- process_triplet(type, triplet_list)
    degree_list[[type]] <- degree_counts
    node_set <- c(node_set, as.character(degree_counts$Node))
  }
  
  Node_unique <- unique(node_set)
  Node_df <- data.frame(Node = Node_unique, stringsAsFactors = FALSE)
  
  return(list(Node_df = Node_df, degree_list = degree_list))
}

# Function: Merge degrees with expression data
merge_degree_with_expression <- function(degree_counts_list, exp_rna = NULL, exp_atac = NULL, cell_types) {
  merged_data <- data.frame()
  
  for (type in cell_types) {
    degree_counts <- degree_counts_list[[type]]
    
    # Split degree_counts into TFs, genes, peaks
    TFs <- degree_counts[degree_counts$NodeType == "TF", ]
    genes <- degree_counts[degree_counts$NodeType == "gene", ]
    peaks <- degree_counts[degree_counts$NodeType == "peak", ]
    
    # Merge with expression data where appropriate
    if (!is.null(exp_rna)) {
      # Get expression data for the cell type
      exp_rna_subset <- exp_rna[, c("Feature", type)]
      colnames(exp_rna_subset) <- c("Feature", "Expression")
      
      # Merge TFs and genes with RNA expression
      TFs_exp <- merge(TFs, exp_rna_subset, by.x = "Node", by.y = "Feature", all.x = TRUE)
      genes_exp <- merge(genes, exp_rna_subset, by.x = "Node", by.y = "Feature", all.x = TRUE)
      
      TFs_exp$FeatureType <- "TF"
      genes_exp$FeatureType <- "gene"
      
      # Add missing columns to TFs and genes
      TFs_exp$Accessibility <- NA
      genes_exp$Accessibility <- NA
      
      # Combine TFs and genes
      merged_data <- bind_rows(merged_data, TFs_exp, genes_exp)
    }
    
    if (!is.null(exp_atac)) {
      # Get accessibility data for the cell type
      exp_atac_subset <- exp_atac[, c("Feature", type)]
      colnames(exp_atac_subset) <- c("Feature", "Accessibility")
      
      # Merge peaks with ATAC data
      peaks_exp <- merge(peaks, exp_atac_subset, by.x = "Node", by.y = "Feature", all.x = TRUE)
      peaks_exp$FeatureType <- "peak"
      
      # Add missing columns to peaks
      peaks_exp$Expression <- NA
      
      # Combine peaks
      merged_data <- bind_rows(merged_data, peaks_exp)
    }
  }
  
  return(merged_data)
}

# Function: Visualization with separated degrees for peaks
plot_degree_expression <- function(data, cell_type, node_type = "TF", degree_type = "OutDegree") {
  # Filter data for the specific cell type and node type
  data_filtered <- data[data$CellType == cell_type & data$NodeType == node_type, ]
  
  # Determine y-axis variable and label
  if (node_type == "peak") {
    y_var <- "Accessibility"
    y_label <- "Peak Accessibility"
    
    # Use specific degree columns for peaks
    if (degree_type == "InDegree") {
      x_var <- "Peak_InDegree"
      x_label <- "TF InDegree"
    } else if (degree_type == "OutDegree") {
      x_var <- "Peak_OutDegree"
      x_label <- "Gene OutDegree"
    } else {
      stop("degree_type must be 'InDegree' or 'OutDegree'")
    }
  } else {
    y_var <- "Expression"
    y_label <- "Normalized Expression"
    
    # Use general degree columns for TFs and genes
    if (degree_type == "InDegree") {
      x_var <- "Gene_InDegree"
      x_label <- "InDegree"
    } else if (degree_type == "OutDegree") {
      x_var <- "TF_OutDegree"
      x_label <- "OutDegree"
    } else {
      stop("degree_type must be 'InDegree' or 'OutDegree'")
    }
  }
  
  # Check if degree and expression/accessibility data are available
  if (!(x_var %in% colnames(data_filtered)) || !(y_var %in% colnames(data_filtered))) {
    stop(paste("Required data columns not found for node type", node_type))
  }
  
  # Plot
  p <- ggplot(data_filtered, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.6) +  # Adjust point transparency
    geom_text_repel(aes(label = Node), size = 3) +
    theme_minimal() +
    labs(
      x = x_label,
      y = y_label,
      title = paste(cell_type, node_type, degree_type, "-", y_label, "Plot")
    )
  
  return(p)
}
