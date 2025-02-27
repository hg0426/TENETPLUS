make_GRN_graph <- function(trimm_df, nodes = "nodes", edges = "edges", node_sizes = "node_sizes",
                           node_label_sizes = "node_label_sizes", node_colors = "node_colors",
                           node_shapes = "node_shapes", min_node_size = 5, min_size_cutoff = 5,
                           min_node_label_size = 1, node_size_adjust = 1) {
  # 1. TF, TG, TAR 열을 합쳐 고유 노드 목록 생성
  temp_nodes <- unique(c(trimm_df$TF, trimm_df$TG, trimm_df$TAR))
  assign(nodes, temp_nodes, envir = .GlobalEnv)
  
  # 2. 노드 모양 지정:
  # 기본값은 "circle"
  # - TG 노드는 "csquare"로 지정 (단, 동시에 TF에 속하면 "circle")
  # - TAR 노드는 "square"
  node_shapes_temp <- rep("circle", length(temp_nodes))
  for (i in seq_along(temp_nodes)) {
    if (temp_nodes[i] %in% trimm_df$TG) {
      node_shapes_temp[i] <- "csquare"
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
  
  # 3. 엣지 생성: (a) TF -> TG, (b) TF -> TAR, (c) TAR -> TG
  temp_edges  <- matrix(c(as.character(trimm_df$TF), as.character(trimm_df$TG)), ncol = 2)
  temp_edges2 <- matrix(c(as.character(trimm_df$TF), as.character(trimm_df$TAR)), ncol = 2)
  temp_edges3 <- matrix(c(as.character(trimm_df$TAR), as.character(trimm_df$TG)), ncol = 2)
  total_edges <- rbind(temp_edges, temp_edges2, temp_edges3)
  assign(edges, total_edges, envir = .GlobalEnv)
  
  # 4. 각 노드의 출현 빈도(degree)를 계산 (TF, TG, TAR 모두 고려)
  count_degree <- function(df) {
    all_nodes <- c(as.character(df$TF), as.character(df$TG), as.character(df$TAR))
    tbl <- table(all_nodes)
    data.frame(Var1 = names(tbl), Freq = as.numeric(tbl))
  }
  input_count <- count_degree(trimm_df)
  size_dict <- as.list(input_count$Freq)
  names(size_dict) <- input_count$Var1
  
  # 5. 노드 크기 결정: 출현 빈도 기반, 최소값 적용
  node_sizes_temp <- sapply(temp_nodes, function(x) {
    if (is.null(size_dict[[x]])) {
      min_node_size
    } else if (size_dict[[x]] <= min_size_cutoff) {
      min_node_size
    } else {
      size_dict[[x]] * node_size_adjust
    }
  })
  assign(node_sizes, node_sizes_temp, envir = .GlobalEnv)
  
  # 6. 노드 레이블 크기 결정: 출현 빈도에 따라 조정
  node_label_sizes_temp <- sapply(temp_nodes, function(x) {
    if (is.null(size_dict[[x]])) {
      min_node_label_size
    } else if (size_dict[[x]] <= min_size_cutoff) {
      min_node_label_size
    } else {
      size_dict[[x]] / 5
    }
  })
  assign(node_label_sizes, node_label_sizes_temp, envir = .GlobalEnv)
  
  # 7. 노드 색상 지정:
  # 기본은 "orange"
  # - TG 노드는 "green" (단, TF에도 속하면 "orange")
  # - TAR 노드는 "skyblue"
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
  
  # 8. igraph 객체 생성: 노드 목록과 엣지를 사용하여 directed graph 생성
  graph <- igraph::graph_from_data_frame(d = total_edges, directed = TRUE, vertices = as.character(temp_nodes))
  return(graph)
}
