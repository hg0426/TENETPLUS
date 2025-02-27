# 1. Register custom vertex shape "triangle" if not already registered
if (!"triangle" %in% names(igraph::vertex.shapes())) {
  # Custom plotting function for triangle shape
  mytriangle <- function(coords, v=NULL, params) {
    # Get vertex color and size from parameters
    vertex.color <- params("vertex", "color")
    vertex.size <- 1/200 * params("vertex", "size")
    
    # Draw a triangle for each vertex
    for (i in 1:nrow(coords)) {
      x <- coords[i, 1]
      y <- coords[i, 2]
      r <- vertex.size[i]
      # Define three angles (in radians) for the triangle vertices.
      angles <- (90 + c(0, 120, 240)) * pi / 180
      xs <- x + r * cos(angles)
      ys <- y + r * sin(angles)
      polygon(xs, ys, col = vertex.color[i], border = "black")
    }
  }
  # Register the new vertex shape "triangle"
  igraph::add.vertex.shape(
    name = "triangle",
    clip = igraph::shape.noclip,
    plot = mytriangle
  )
}

# 2. Function to create a GRN graph from TENETobj@output$result$Trimm
make_GRN_graph <- function(trimm_df, nodes = "nodes", edges = "edges", node_sizes = "node_sizes",
                           node_label_sizes = "node_label_sizes", node_colors = "node_colors",
                           node_shapes = "node_shapes", min_node_size = 5, min_size_cutoff = 5,
                           min_node_label_size = 1, node_size_adjust = 1) {
  # Create unique node list from TF, TG, and TAR columns
  temp_nodes <- unique(c(trimm_df$TF, trimm_df$TG, trimm_df$TAR))
  assign(nodes, temp_nodes, envir = .GlobalEnv)
  
  # Define node shapes:
  # - 기본값 "circle"
  # - TG (여기서는 gene) -> "triangle" (단, 동시에 TF에 속하면 "circle")
  # - TAR -> "square"
  node_shapes_temp <- rep("circle", length(temp_nodes))
  for (i in seq_along(temp_nodes)) {
    if (temp_nodes[i] %in% trimm_df$TG) {
      node_shapes_temp[i] <- "triangle"
      if (temp_nodes[i] %in% trimm_df$TF) {
        node_shapes_temp[i] <- "circle"
      }
    } else if (temp_nodes[i] %in% trimm_df$TF) {
      node_shapes_temp[i] <- "circle"
    } else if (temp_nodes[i] %in% trimm_df$TAR) {
      node_shapes_temp[i] <- "square"
    }
  }
  assign(node_shapes, node_shapes_temp, envir = .GlobalEnv)
  
  # Create edge lists as 2-column matrices:
  # (a) TF -> TG
  temp_edges  <- matrix(c(as.character(trimm_df$TF), as.character(trimm_df$TG)),
                        ncol = 2)
  # (b) TF -> TAR
  temp_edges2 <- matrix(c(as.character(trimm_df$TF), as.character(trimm_df$TAR)),
                        ncol = 2)
  # (c) TAR -> TG
  temp_edges3 <- matrix(c(as.character(trimm_df$TAR), as.character(trimm_df$TG)),
                        ncol = 2)
  total_edges <- rbind(temp_edges, temp_edges2, temp_edges3)
  assign(edges, total_edges, envir = .GlobalEnv)
  
  # Count node occurrences (degree) across TF, TG, TAR
  count_degree <- function(df) {
    all_nodes <- c(as.character(df$TF), as.character(df$TG), as.character(df$TAR))
    tbl <- table(all_nodes)
    return(data.frame(Var1 = names(tbl), Freq = as.numeric(tbl)))
  }
  input_count <- count_degree(trimm_df)
  size_dict <- as.list(input_count$Freq)
  names(size_dict) <- input_count$Var1
  
  # Set node sizes based on frequency
  node_sizes_temp <- sapply(temp_nodes, function(x) {
    if (is.null(size_dict[[x]])) {
      return(min_node_size)
    } else if (size_dict[[x]] <= min_size_cutoff) {
      return(min_node_size)
    } else {
      return(size_dict[[x]] * node_size_adjust)
    }
  })
  assign(node_sizes, node_sizes_temp, envir = .GlobalEnv)
  
  # Set node label sizes based on frequency
  node_label_sizes_temp <- sapply(temp_nodes, function(x) {
    if (is.null(size_dict[[x]])) {
      return(min_node_label_size)
    } else if (size_dict[[x]] <= min_size_cutoff) {
      return(min_node_label_size)
    } else {
      return(size_dict[[x]] / 5)
    }
  })
  assign(node_label_sizes, node_label_sizes_temp, envir = .GlobalEnv)
  
  # Define node colors:
  # - 기본값 "orange"
  # - TG -> "green" (단, 동시에 TF에 속하면 "orange")
  # - TAR -> "skyblue"
  node_colors_temp <- rep("orange", length(temp_nodes))
  for (i in seq_along(temp_nodes)) {
    if (temp_nodes[i] %in% trimm_df$TG) {
      node_colors_temp[i] <- "green"
      if (temp_nodes[i] %in% trimm_df$TF) {
        node_colors_temp[i] <- "orange"
      }
    } else if (temp_nodes[i] %in% trimm_df$TF) {
      node_colors_temp[i] <- "orange"
    } else if (temp_nodes[i] %in% trimm_df$TAR) {
      node_colors_temp[i] <- "skyblue"
    }
  }
  assign(node_colors, node_colors_temp, envir = .GlobalEnv)
  
  # Create a directed graph from total_edges with vertices from temp_nodes using igraph::
  graph <- igraph::graph_from_data_frame(d = total_edges, directed = TRUE, vertices = as.character(temp_nodes))
  return(graph)
}

# # ---------------------------------------------------------------------------
# # Example usage:
# # Assuming TENETobj@output$result$Trimm is your data frame input.
# trimm_df <- TENETobj@output$result$Trimm

# # Create the network graph using the function
# grn_graph <- make_GRN_graph(trimm_df)

# # Calculate layout using Fruchterman-Reingold algorithm
# layout <- igraph::layout_with_fr(grn_graph)

# # Set reproducibility and plot options
# set.seed(123)
# options(repr.plot.width = 10, repr.plot.height = 10)

# # Plot the graph using global variables for node attributes
# plot(grn_graph,
#      vertex.label = get("nodes", envir = .GlobalEnv),
#      vertex.size = get("node_sizes", envir = .GlobalEnv),
#      edge.arrow.size = 0.2,
#      layout = layout,
#      vertex.color = get("node_colors", envir = .GlobalEnv),
#      vertex.shape = get("node_shapes", envir = .GlobalEnv),
#      vertex.label.cex = get("node_label_sizes", envir = .GlobalEnv),
#      vertex.label.color = "black")
