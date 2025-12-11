#' Cluster graph using multiple algorithms
#'
#' @param g igraph object
#' @param method One of c("infomap","louvain","walktrap","leiden")
#' @param weights Edge weights to use (default: E(g)$weight)
#' @param ... Extra parameters passed to the underlying clustering function
#' @return The graph object with community attribute added
#' @export
cluster_graph <- function(g, method = c("infomap","louvain","walktrap","leiden"),
                          weights = igraph::E(g)$weight, ...) {
  method <- match.arg(method)
  
  # 在进行聚类前，先获取图的完整数据，确保顶点顺序固定
  g_list <- igraph::as_data_frame(g, what = "both")
  
  # 确保顶点数据框存在，如果不存在则创建
  if (is.null(g_list$vertices)) {
    # 获取所有顶点属性
    v_names <- igraph::V(g)$name
    g_list$vertices <- data.frame(
      name = v_names,
      stringsAsFactors = FALSE
    )
    # 获取所有其他顶点属性
    v_attrs <- igraph::vertex_attr(g)
    if (length(v_attrs) > 1) {
      for (attr_name in names(v_attrs)) {
        if (attr_name != "name") {
          g_list$vertices[[attr_name]] <- v_attrs[[attr_name]]
        }
      }
    }
  }
  
  # 确保边数据框有权重
  if (is.null(g_list$edges$weight)) {
    g_list$edges$weight <- igraph::E(g)$weight
  }
  
  # 在原始图上进行聚类（使用原始图对象）
  if (method == "infomap") {
    comm <- igraph::cluster_infomap(g, e.weights = weights, ...)
  } else if (method == "louvain") {
    comm <- igraph::cluster_louvain(g, weights = weights, ...)
  } else if (method == "walktrap") {
    comm <- igraph::cluster_walktrap(g, weights = weights, ...)
  } else if (method == "leiden") {
    if (!"cluster_leiden" %in% ls("package:igraph")) {
      stop("cluster_leiden is not available in your igraph version.")
    }
    comm <- igraph::cluster_leiden(g, weights = weights, ...)
  }
  
  # 获取节点数量和 membership
  n_nodes <- igraph::vcount(g)
  membership_vec <- as.integer(comm$membership)
  membership_vec[is.na(membership_vec)] <- 1L
  
  # 验证长度匹配
  if (length(membership_vec) != n_nodes) {
    stop(sprintf("Cluster membership length (%d) does not match number of vertices (%d).", 
                 length(membership_vec), n_nodes))
  }
  
  # 关键：确保 membership_vec 的顺序与 g_list$vertices 中的顶点顺序完全一致
  # membership_vec 的顺序是 igraph::V(g) 的顺序，需要确保与 g_list$vertices 的 name 顺序一致
  v_names_from_g <- as.character(igraph::V(g)$name)
  v_names_from_list <- as.character(g_list$vertices$name)
  
  # 如果顶点名称顺序不一致，需要重新排序 membership_vec
  if (!identical(v_names_from_g, v_names_from_list)) {
    # 创建名称到 membership 的映射
    name_to_membership <- stats::setNames(membership_vec, v_names_from_g)
    # 按照 g_list$vertices 的顺序重新排序 membership
    membership_vec <- as.integer(name_to_membership[v_names_from_list])
    membership_vec[is.na(membership_vec)] <- 1L
  }
  
  # 将 community 属性添加到顶点数据框
  g_list$vertices$community <- membership_vec
  
  # 确保所有顶点属性都是简单类型（不是因子）
  for (col_name in colnames(g_list$vertices)) {
    if (col_name != "name") {
      if (is.factor(g_list$vertices[[col_name]])) {
        g_list$vertices[[col_name]] <- as.character(g_list$vertices[[col_name]])
      }
    }
  }
  
  # 重建图，此时所有属性都在构建时设置，完全避免 V<- 问题
  # 这是最底层、最安全的方法
  g_new <- igraph::graph_from_data_frame(
    d = g_list$edges,
    directed = igraph::is_directed(g),
    vertices = g_list$vertices
  )
  
  # 验证重建后的图
  if (igraph::vcount(g_new) != n_nodes) {
    stop(sprintf("Graph reconstruction failed: original has %d vertices, new has %d vertices.", 
                 n_nodes, igraph::vcount(g_new)))
  }
  
  g_new
}


