#' Read network and node group tables
#'
#' @param edge_file Path to edges table with columns: edge_id, weight, node1, node2
#' @param node_group_file Path to node groups table with columns including id, clade, subclade
#' @return A list with elements: network_data (data.frame), node_groups (data.frame)
#' @export
read_network_tables <- function(edge_file, node_group_file) {
  # 读取边文件，自动检测列数与格式
  network_data <- utils::read.table(edge_file, header = FALSE, sep = "\t",
                                    stringsAsFactors = FALSE, fill = TRUE, quote = "\"", comment.char = "")
  n_cols <- ncol(network_data)
  if (n_cols >= 4) {
    # 常见格式：edge_id, weight, node1, node2
    colnames(network_data)[1:4] <- c("edge_id", "weight", "node1", "node2")
  } else if (n_cols == 3) {
    # 退化格式：node1, node2, weight
    colnames(network_data)[1:3] <- c("node1", "node2", "weight")
  } else if (n_cols == 2) {
    # 最简格式：node1, node2（补权重=1）
    colnames(network_data)[1:2] <- c("node1", "node2")
    network_data$weight <- 1
  } else {
    stop("Edge file must have at least 2 columns (node1, node2)")
  }

  # 标准化节点与权重类型
  if ("node1" %in% names(network_data)) network_data$node1 <- as.character(network_data$node1)
  if ("node2" %in% names(network_data)) network_data$node2 <- as.character(network_data$node2)
  if ("weight" %in% names(network_data)) {
    network_data$weight <- suppressWarnings(as.numeric(network_data$weight))
    network_data$weight[is.na(network_data$weight)] <- 1
  } else {
    network_data$weight <- 1
  }
  # 去除节点名首尾空白
  network_data$node1 <- trimws(network_data$node1)
  network_data$node2 <- trimws(network_data$node2)
  
  # 读取节点分组文件
  node_groups <- utils::read.table(node_group_file, header = TRUE, sep = "\t",
                                   stringsAsFactors = FALSE, fill = TRUE, quote = "\"", comment.char = "")
  
  # 确保所有字符串列都是字符类型
  if ("id" %in% colnames(node_groups)) {
    node_groups$id <- as.character(node_groups$id)
    node_groups$id <- trimws(node_groups$id)
  }
  if ("clade" %in% colnames(node_groups)) {
    node_groups$clade <- as.character(node_groups$clade)
  }
  if ("subclade" %in% colnames(node_groups)) {
    node_groups$subclade <- as.character(node_groups$subclade)
  }
  
  # 确保边的节点名称也是字符类型
  network_data$node1 <- as.character(network_data$node1)
  network_data$node2 <- as.character(network_data$node2)
  network_data$weight <- as.numeric(network_data$weight)
  network_data$weight[is.na(network_data$weight)] <- 1
  
  list(network_data = network_data, node_groups = node_groups)
}

#' Build igraph from tables and attach attributes
#'
#' @param network_data Data frame with columns node1, node2, weight
#' @param node_groups Data frame with columns id, clade, subclade
#' @return igraph object with vertex attributes clade, subclade
#' @export
build_fab_graph <- function(network_data, node_groups) {
  # 确保数据格式正确
  if (!"node1" %in% colnames(network_data) || !"node2" %in% colnames(network_data)) {
    stop("network_data must have columns 'node1' and 'node2'")
  }
  
  # 确保节点名称是字符类型
  network_data$node1 <- as.character(network_data$node1)
  network_data$node2 <- as.character(network_data$node2)
  
  # 获取所有唯一的节点名称
  all_nodes <- unique(c(network_data$node1, network_data$node2))
  all_nodes <- as.character(all_nodes)
  n_nodes <- length(all_nodes)
  
  # 初始化属性向量
  node_attrs_clade <- character(n_nodes)
  node_attrs_subclade <- character(n_nodes)
  
  # 检查 node_groups 是否有必要的列，如果没有则创建默认值
  if (!"id" %in% colnames(node_groups)) {
    warning("node_groups missing 'id' column, creating default attributes")
    node_attrs_clade[] <- "Unknown"
    node_attrs_subclade[] <- "Unknown"
  } else {
    # 确保 id 列是字符类型
    node_groups$id <- as.character(node_groups$id)
    
    # 检查是否有 clade 和 subclade 列
    has_clade <- "clade" %in% colnames(node_groups)
    has_subclade <- "subclade" %in% colnames(node_groups)
    
    # 创建匹配索引
    matched_idx <- match(all_nodes, node_groups$id)
    
    # 逐个节点设置属性
    for (i in seq_len(n_nodes)) {
      idx <- matched_idx[i]
      
      if (!is.na(idx) && idx > 0 && idx <= nrow(node_groups)) {
        # 匹配成功
        if (has_clade) {
          clade_val <- node_groups$clade[idx]
          if (is.na(clade_val) || clade_val == "" || length(clade_val) == 0) {
            node_attrs_clade[i] <- "Unknown"
          } else {
            node_attrs_clade[i] <- as.character(clade_val)[1]
          }
        } else {
          node_attrs_clade[i] <- "Unknown"
        }
        
        if (has_subclade) {
          subclade_val <- node_groups$subclade[idx]
          if (is.na(subclade_val) || subclade_val == "" || length(subclade_val) == 0) {
            node_attrs_subclade[i] <- node_attrs_clade[i]
          } else {
            node_attrs_subclade[i] <- as.character(subclade_val)[1]
          }
        } else {
          node_attrs_subclade[i] <- node_attrs_clade[i]
        }
      } else {
        # 未匹配到的节点
        node_attrs_clade[i] <- "Unknown"
        node_attrs_subclade[i] <- "Unknown"
      }
    }
  }
  
  # 确保所有属性都是有效的字符向量，不是因子
  node_attrs_clade <- as.character(node_attrs_clade)
  node_attrs_subclade <- as.character(node_attrs_subclade)
  node_attrs_clade[is.na(node_attrs_clade) | node_attrs_clade == ""] <- "Unknown"
  node_attrs_subclade[is.na(node_attrs_subclade) | node_attrs_subclade == ""] <- "Unknown"
  
  # 确保所有节点的名称都是字符类型，不是因子
  all_nodes <- as.character(all_nodes)
  
  # 验证长度完全匹配
  if (length(all_nodes) != length(node_attrs_clade) || 
      length(all_nodes) != length(node_attrs_subclade)) {
    stop(sprintf("Attribute length mismatch: nodes=%d, clade=%d, subclade=%d", 
                 length(all_nodes), length(node_attrs_clade), length(node_attrs_subclade)))
  }
  
  # 创建顶点数据框，在构建图时直接添加属性（这是最安全的方法）
  # 确保所有列都是简单类型，不是因子
  vertices_df <- data.frame(
    name = all_nodes,
    clade = node_attrs_clade,
    subclade = node_attrs_subclade,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # 再次验证，确保没有因子类型
  if (is.factor(vertices_df$name) || is.factor(vertices_df$clade) || is.factor(vertices_df$subclade)) {
    vertices_df$name <- as.character(vertices_df$name)
    vertices_df$clade <- as.character(vertices_df$clade)
    vertices_df$subclade <- as.character(vertices_df$subclade)
  }
  
  # 确保边数据框的节点名称也是字符类型
  network_data$node1 <- as.character(network_data$node1)
  network_data$node2 <- as.character(network_data$node2)
  
  # 构建图，直接在构建时添加顶点属性
  # 这样可以完全避免后续设置属性时的索引问题
  g <- igraph::graph_from_data_frame(
    d = network_data[, c("node1", "node2", "weight")],
    directed = FALSE,
    vertices = vertices_df
  )
  
  # 验证图的顶点数量
  if (igraph::vcount(g) != length(all_nodes)) {
    stop(sprintf("Graph construction failed: expected %d vertices, got %d vertices.", 
                 length(all_nodes), igraph::vcount(g)))
  }
  
  g
}


