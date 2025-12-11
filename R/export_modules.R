#' Export network modules as separate files and visualizations
#'
#' This module extracts closely connected node groups from the network and exports
#' each module as a table file and a separate network visualization.
#'
#' @param g igraph object with community/subclade attributes
#' @param module_by character, grouping method: "subclade", "clade", or "community"
#' @param output_dir character, directory to save output files (default: "network_modules")
#' @param min_nodes integer, minimum number of nodes in a module to export (default: 2)
#' @param plot_each logical, whether to generate plots for each module (default: TRUE)
#' @param plot_width numeric, width for individual module plots (default: 8)
#' @param plot_height numeric, height for individual module plots (default: 8)
#' @param plot_dpi numeric, DPI for plots (default: 300)
#' @return list with module_info (summary), module_tables (list of data.frames), module_graphs (list of igraph objects)
#' @export
export_network_modules <- function(g, module_by = c("subclade", "clade", "community"),
                                   output_dir = "network_modules",
                                   min_nodes = 2, plot_each = TRUE,
                                   plot_width = 8, plot_height = 8, plot_dpi = 300) {
  module_by <- match.arg(module_by)
  
  # 获取模块分组信息
  if (module_by == "subclade") {
    module_vec <- igraph::V(g)$subclade
  } else if (module_by == "clade") {
    module_vec <- igraph::V(g)$clade
  } else {
    module_vec <- igraph::V(g)$community
  }
  
  module_vec <- as.character(module_vec)
  module_vec[is.na(module_vec) | module_vec == ""] <- "Unknown"
  
  # 获取所有模块
  unique_modules <- unique(module_vec)
  module_counts <- table(module_vec)
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  module_tables <- list()
  module_graphs <- list()
  module_info <- data.frame(
    module = character(),
    node_count = integer(),
    edge_count = integer(),
    table_file = character(),
    plot_file = character(),
    stringsAsFactors = FALSE
  )
  
  # 处理每个模块
  for (mod in unique_modules) {
    mod_nodes <- which(module_vec == mod)
    
    if (length(mod_nodes) < min_nodes) {
      next
    }
    
    # 提取子图
    subg <- igraph::induced_subgraph(g, mod_nodes)
    
    # 创建节点表格
    node_table <- data.frame(
      id = igraph::V(subg)$name,
      module = mod,
      degree = igraph::degree(subg),
      stringsAsFactors = FALSE
    )
    
    # 添加其他属性（如果存在）
    if ("clade" %in% igraph::vertex_attr_names(subg)) {
      node_table$clade <- igraph::V(subg)$clade
    }
    if ("subclade" %in% igraph::vertex_attr_names(subg)) {
      node_table$subclade <- igraph::V(subg)$subclade
    }
    if ("community" %in% igraph::vertex_attr_names(subg)) {
      node_table$community <- igraph::V(subg)$community
    }
    
    # 创建边表格
    edge_list <- igraph::as_edgelist(subg, names = TRUE)
    edge_table <- data.frame(
      from = edge_list[, 1],
      to = edge_list[, 2],
      weight = igraph::E(subg)$weight,
      stringsAsFactors = FALSE
    )
    
    # 保存表格文件
    safe_mod_name <- gsub("[^A-Za-z0-9_]", "_", mod)
    node_file <- file.path(output_dir, paste0("module_", safe_mod_name, "_nodes.tsv"))
    edge_file <- file.path(output_dir, paste0("module_", safe_mod_name, "_edges.tsv"))
    
    utils::write.table(node_table, file = node_file, sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(edge_table, file = edge_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    module_tables[[mod]] <- list(nodes = node_table, edges = edge_table)
    module_graphs[[mod]] <- subg
    
    # 生成可视化（如果启用）
    plot_file <- NULL
    if (plot_each && length(mod_nodes) >= 2) {
      plot_file <- file.path(output_dir, paste0("module_", safe_mod_name, "_network.png"))
      plot_module_network(subg, mod, plot_file, width = plot_width, height = plot_height, dpi = plot_dpi)
    }
    
    # 记录模块信息
    module_info <- rbind(module_info, data.frame(
      module = mod,
      node_count = length(mod_nodes),
      edge_count = igraph::ecount(subg),
      table_file = paste0("module_", safe_mod_name, "_*.tsv"),
      plot_file = ifelse(is.null(plot_file), "", basename(plot_file)),
      stringsAsFactors = FALSE
    ))
  }
  
  # 保存模块汇总信息
  summary_file <- file.path(output_dir, "module_summary.tsv")
  utils::write.table(module_info, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("Exported %d modules to: %s\n", nrow(module_info), output_dir))
  cat(sprintf("Summary file: %s\n", summary_file))
  
  list(
    module_info = module_info,
    module_tables = module_tables,
    module_graphs = module_graphs,
    output_dir = output_dir
  )
}

#' Plot individual module network
#'
#' @param g_module igraph object for a single module
#' @param module_name character, name of the module
#' @param output_file character, path to save the plot
#' @param width numeric, plot width
#' @param height numeric, plot height
#' @param dpi numeric, plot DPI
#' @return ggplot object (invisibly)
plot_module_network <- function(g_module, module_name, output_file,
                                width = 8, height = 8, dpi = 300) {
  # 计算布局
  if (igraph::vcount(g_module) <= 1) {
    warning(sprintf("Module %s has too few nodes for layout", module_name))
    return(NULL)
  }
  
  # 使用力导向布局
  layout_mat <- igraph::layout_with_fr(g_module, 
                                       weights = igraph::E(g_module)$weight,
                                       niter = 500)
  
  # 准备绘图数据
  node_df <- data.frame(
    x = layout_mat[, 1],
    y = layout_mat[, 2],
    name = igraph::V(g_module)$name,
    degree = igraph::degree(g_module),
    size = (1.5 + sqrt(igraph::degree(g_module)) * 0.8) * 1.5,
    stringsAsFactors = FALSE
  )
  
  edge_list <- igraph::as_edgelist(g_module, names = TRUE)
  edge_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    weight = igraph::E(g_module)$weight,
    stringsAsFactors = FALSE
  )
  
  # 添加坐标
  edge_df$x <- node_df$x[match(edge_df$from, node_df$name)]
  edge_df$y <- node_df$y[match(edge_df$from, node_df$name)]
  edge_df$xend <- node_df$x[match(edge_df$to, node_df$name)]
  edge_df$yend <- node_df$y[match(edge_df$to, node_df$name)]
  
  # 绘制
  gg <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_df,
                         ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                         color = scales::alpha("#808080", 0.3), linewidth = 0.2) +
    ggplot2::geom_point(data = node_df,
                      ggplot2::aes(x = x, y = y, size = size),
                      color = "#E66F73", fill = "#E66F73", alpha = 0.7, shape = 21) +
    ggplot2::scale_size(range = c(2, 8)) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste("Module:", module_name),
                 subtitle = sprintf("Nodes: %d, Edges: %d", 
                                  igraph::vcount(g_module), 
                                  igraph::ecount(g_module))) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                 plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
                 legend.position = "none")
  
  ggplot2::ggsave(filename = output_file, plot = gg, width = width, height = height, dpi = dpi)
  invisible(gg)
}

#' Export modules from fabnetsyn_build_plot result
#'
#' Convenience function to export modules directly from fabnetsyn_build_plot output
#'
#' @param plot_result list, result from fabnetsyn_build_plot
#' @param module_by character, grouping method (default: "subclade")
#' @param output_dir character, output directory
#' @param min_nodes integer, minimum nodes per module
#' @param plot_each logical, whether to plot each module
#' @param ... additional parameters passed to export_network_modules
#' @return list with exported module information
#' @export
export_modules_from_result <- function(plot_result, module_by = c("subclade", "clade", "community"),
                                       output_dir = "network_modules", min_nodes = 2,
                                       plot_each = TRUE, ...) {
  module_by <- match.arg(module_by)
  
  # 需要从原始图中获取模块信息
  # 由于 plot_result 只包含 node.data 和 edge.data，我们需要重建图
  # 或者从 node.data 中提取模块信息
  
  if (!"node.data" %in% names(plot_result)) {
    stop("plot_result must contain 'node.data' from fabnetsyn_build_plot")
  }
  
  # 从 node.data 重建基本的图结构用于导出
  node_data <- plot_result$node.data
  edge_data <- plot_result$edge.data
  
  # 创建简化图用于导出
  g_export <- igraph::graph_from_data_frame(
    d = edge_data[, c("from", "to", "weight")],
    directed = FALSE,
    vertices = data.frame(
      name = node_data$name,
      Module = node_data$Module,
      stringsAsFactors = FALSE
    )
  )
  
  # 根据 module_by 设置顶点属性
  if (module_by == "subclade") {
    igraph::V(g_export)$subclade <- node_data$Module
    igraph::V(g_export)$clade <- node_data$clade
  } else if (module_by == "clade") {
    igraph::V(g_export)$clade <- node_data$Module
  } else {
    # 尝试从 Module 推断 community
    igraph::V(g_export)$community <- as.integer(as.factor(node_data$Module))
  }
  
  # 调用导出函数
  export_network_modules(g_export, module_by = module_by, output_dir = output_dir,
                        min_nodes = min_nodes, plot_each = plot_each, ...)
}

#' Extract and export highly connected clusters from specific node groups
#'
#' This function extracts nodes from a specific group (e.g., NPR1), identifies
#' highly connected clusters within them, and exports the clusters with their
#' connections preserved.
#'
#' @param plot_result list, result from fabnetsyn_build_plot containing graph object
#' @param hit_ids character vector, IDs of nodes to analyze (e.g., NPR1 nodes)
#' @param cluster_method character, clustering method: "infomap", "louvain", "walktrap", "leiden"
#' @param top_n integer, number of top clusters to export (by edge density, default: NULL for all)
#' @param min_edges integer, minimum number of edges per cluster (default: 2)
#' @param output_dir character, directory to save output files (default: "hit_clusters")
#' @param plot_each logical, whether to generate plots for each cluster (default: TRUE)
#' @param plot_width numeric, width for individual cluster plots (default: 10)
#' @param plot_height numeric, height for individual cluster plots (default: 10)
#' @param plot_dpi numeric, DPI for plots (default: 300)
#' @param keep_connections logical, whether to keep connections between clusters (default: TRUE)
#' @return list with cluster_info (summary), cluster_tables (list), cluster_graphs (list)
#' @export
export_hit_clusters <- function(plot_result, hit_ids, 
                                cluster_method = c("infomap", "louvain", "walktrap", "leiden"),
                                top_n = NULL, min_edges = 2,
                                output_dir = "hit_clusters",
                                plot_each = TRUE,
                                plot_width = 10, plot_height = 10, plot_dpi = 300,
                                keep_connections = TRUE) {
  cluster_method <- match.arg(cluster_method)
  
  # 获取图对象
  if (!"graph" %in% names(plot_result)) {
    stop("plot_result must contain 'graph' object from fabnetsyn_build_plot")
  }
  
  g <- plot_result$graph
  
  # ===== 规范化 hit_ids：支持数据框、列表、向量等 =====
  # 如果用户没有提供 hit_ids（缺省或 NULL），则默认把“所有节点”当作 hit，
  # 这样可以导出整个图上所有高度连接的聚类块。
  if (missing(hit_ids) || is.null(hit_ids)) {
    hit_ids <- igraph::V(g)$name
  }
  
  if (is.data.frame(hit_ids)) {
    # 如果是数据框，提取第一列
    if (ncol(hit_ids) >= 1) {
      hit_ids <- hit_ids[[1]]
    } else {
      stop("hit_ids data.frame must have at least one column")
    }
  } else if (is.list(hit_ids) && !is.vector(hit_ids)) {
    # 如果是列表（而不是普通向量），提取第一个元素
    hit_ids <- hit_ids[[1]]
  }
  
  # 转换为字符向量并去重
  hit_ids <- as.character(hit_ids)
  hit_ids <- hit_ids[!is.na(hit_ids) & hit_ids != ""]
  hit_ids <- unique(trimws(hit_ids))  # 去除空格
  
  # 如果归一化后仍然为空，回退为“全图所有节点”，避免直接报错，
  # 便于用户传入 result 等对象时依然能得到聚类结果。
  if (length(hit_ids) == 0) {
    warning("hit_ids is empty after processing; using all nodes in the graph as hit_ids.")
    hit_ids <- igraph::V(g)$name
  }
  
  # 获取图中所有节点名称（去除空格以便匹配）
  all_nodes <- as.character(igraph::V(g)$name)
  all_nodes_trimmed <- trimws(all_nodes)
  
  # 找到匹配的节点索引（支持精确匹配和去除空格后的匹配）
  matched_nodes <- which(all_nodes %in% hit_ids | all_nodes_trimmed %in% hit_ids)
  
  # 如果仍然没有匹配，尝试模糊匹配
  if (length(matched_nodes) == 0) {
    # 尝试大小写不敏感匹配
    hit_ids_lower <- tolower(hit_ids)
    all_nodes_lower <- tolower(all_nodes)
    matched_nodes <- which(all_nodes_lower %in% hit_ids_lower)
  }
  
  if (length(matched_nodes) == 0) {
    # 提供详细的调试信息
    cat("=== Debug Information ===\n")
    cat(sprintf("Input hit_ids count: %d\n", length(hit_ids)))
    cat(sprintf("First 10 hit_ids: %s\n", paste(head(hit_ids, 10), collapse = ", ")))
    cat(sprintf("Graph nodes count: %d\n", length(all_nodes)))
    cat(sprintf("First 10 graph nodes: %s\n", paste(head(all_nodes, 10), collapse = ", ")))
    
    # 检查是否有部分匹配
    partial_matches <- sapply(hit_ids, function(id) {
      any(grepl(id, all_nodes, fixed = TRUE)) || any(grepl(all_nodes, id, fixed = TRUE))
    })
    if (any(partial_matches)) {
      cat(sprintf("Found partial matches for: %s\n", paste(hit_ids[partial_matches], collapse = ", ")))
    }
    
    warning(sprintf("No matching nodes found in the graph. Input %d IDs, graph has %d nodes. Check if IDs match or if nodes were filtered out.", 
                   length(hit_ids), length(all_nodes)))
    return(list(cluster_info = data.frame(), cluster_tables = list(), cluster_graphs = list()))
  }
  
  # 显示匹配信息
  matched_node_names <- all_nodes[matched_nodes]
  unmatched_ids <- setdiff(hit_ids, matched_node_names)
  if (length(unmatched_ids) > 0 && length(unmatched_ids) <= 10) {
    cat(sprintf("Matched %d/%d nodes. Unmatched IDs: %s\n", 
               length(matched_nodes), length(hit_ids), 
               paste(unmatched_ids, collapse = ", ")))
  } else if (length(unmatched_ids) > 10) {
    cat(sprintf("Matched %d/%d nodes. %d unmatched IDs (showing first 10): %s\n", 
               length(matched_nodes), length(hit_ids), length(unmatched_ids),
               paste(head(unmatched_ids, 10), collapse = ", ")))
  } else {
    cat(sprintf("Matched %d/%d nodes successfully.\n", length(matched_nodes), length(hit_ids)))
  }
  
  # 更新hit_ids为实际匹配的节点名称（用于后续标记）
  matched_hit_ids <- matched_node_names
  
  # 提取包含hit_ids的子图（包括它们之间的连接）
  # 如果keep_connections=TRUE，提取包含所有hit_ids及其邻居的扩展子图
  if (keep_connections) {
    # 获取hit_ids节点的邻居
    neighbors_list <- lapply(matched_nodes, function(i) {
      igraph::neighbors(g, i)$name
    })
    all_neighbors <- unique(unlist(neighbors_list))
    
    # 合并hit_ids和邻居节点
    extended_nodes <- unique(c(matched_hit_ids, all_neighbors))
    extended_indices <- which(all_nodes %in% extended_nodes)
    
    # 提取子图
    subg <- igraph::induced_subgraph(g, extended_indices)
    
    # 标记哪些是原始hit_ids（使用匹配的节点名称）
    igraph::V(subg)$is_hit <- as.character(igraph::V(subg)$name) %in% matched_hit_ids
  } else {
    # 只提取hit_ids节点及其之间的连接
    subg <- igraph::induced_subgraph(g, matched_nodes)
    igraph::V(subg)$is_hit <- rep(TRUE, length(matched_nodes))
  }
  
  if (igraph::vcount(subg) < 2) {
    warning("Subgraph has too few nodes for clustering")
    return(list(cluster_info = data.frame(), cluster_tables = list(), cluster_graphs = list()))
  }
  
  # 在子图上进行聚类
  cluster_result <- cluster_graph(subg, method = cluster_method)
  cluster_membership <- igraph::V(cluster_result)$community
  
  # 获取所有聚类
  unique_clusters <- unique(cluster_membership)
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cluster_info_list <- list()
  cluster_tables <- list()
  cluster_graphs <- list()
  
  # 计算每个聚类的统计信息
  cluster_stats <- data.frame(
    cluster_id = integer(),
    node_count = integer(),
    hit_node_count = integer(),
    edge_count = integer(),
    edge_density = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cid in unique_clusters) {
    cluster_nodes <- which(cluster_membership == cid)
    cluster_subg <- igraph::induced_subgraph(cluster_result, cluster_nodes)
    
    edge_count <- igraph::ecount(cluster_subg)
    
    if (edge_count < min_edges) {
      next
    }
    
    # 计算连接密度（边数/可能的边数）
    n_nodes <- igraph::vcount(cluster_subg)
    max_edges <- n_nodes * (n_nodes - 1) / 2
    edge_density <- ifelse(max_edges > 0, edge_count / max_edges, 0)
    
    # 统计hit节点数量
    hit_node_count <- sum(igraph::V(cluster_subg)$is_hit)
    
    cluster_stats <- rbind(cluster_stats, data.frame(
      cluster_id = cid,
      node_count = n_nodes,
      hit_node_count = hit_node_count,
      edge_count = edge_count,
      edge_density = edge_density,
      stringsAsFactors = FALSE
    ))
  }
  
  # 按边数或密度排序（优先边数）
  cluster_stats <- cluster_stats[order(cluster_stats$edge_count, decreasing = TRUE), ]
  
  # 选择top_n个聚类
  if (!is.null(top_n) && top_n > 0) {
    cluster_stats <- cluster_stats[1:min(top_n, nrow(cluster_stats)), ]
  }
  # 为选中的聚类生成连续编号（用于图中编号与表格编号）
  cluster_stats$cluster_index <- seq_len(nrow(cluster_stats))
  
  # 导出选定的聚类
  for (i in 1:nrow(cluster_stats)) {
    cid <- cluster_stats$cluster_id[i]
    cidx <- cluster_stats$cluster_index[i]
    cluster_nodes <- which(cluster_membership == cid)
    cluster_subg <- igraph::induced_subgraph(cluster_result, cluster_nodes)
    
    # 为该聚类内的节点生成连续编号（在该聚类内部从 1,2,... 编号）
    node_names <- igraph::V(cluster_subg)$name
    node_index <- seq_along(node_names)
    igraph::V(cluster_subg)$node_index <- node_index
    
    # 创建节点表格
    node_table <- data.frame(
      id = igraph::V(cluster_subg)$name,
      cluster_id = cid,
      cluster_index = cidx,
      node_index = node_index,
      is_hit = igraph::V(cluster_subg)$is_hit,
      degree = igraph::degree(cluster_subg),
      stringsAsFactors = FALSE
    )
    
    # 添加其他属性（如果存在）
    if ("clade" %in% igraph::vertex_attr_names(cluster_subg)) {
      node_table$clade <- igraph::V(cluster_subg)$clade
    }
    if ("subclade" %in% igraph::vertex_attr_names(cluster_subg)) {
      node_table$subclade <- igraph::V(cluster_subg)$subclade
    }
    
    # 创建边表格（保持连接关系）
    edge_list <- igraph::as_edgelist(cluster_subg, names = TRUE)
    edge_table <- data.frame(
      from = edge_list[, 1],
      to = edge_list[, 2],
      weight = igraph::E(cluster_subg)$weight,
      cluster_id = cid,
      cluster_index = cidx,
      stringsAsFactors = FALSE
    )
    
    # 保存表格文件：使用 cluster_index 统一编号，便于与图形编号对齐
    safe_cluster_name <- paste0("cluster_", cidx)
    node_file <- file.path(output_dir, paste0(safe_cluster_name, "_nodes.tsv"))
    edge_file <- file.path(output_dir, paste0(safe_cluster_name, "_edges.tsv"))
    
    utils::write.table(node_table, file = node_file, sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(edge_table, file = edge_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cluster_tables[[as.character(cid)]] <- list(nodes = node_table, edges = edge_table)
    cluster_graphs[[as.character(cid)]] <- cluster_subg
    
    # 生成可视化（如果启用）
    plot_file <- NULL
    if (plot_each && igraph::vcount(cluster_subg) >= 2) {
      plot_file <- file.path(output_dir, paste0(safe_cluster_name, "_network.png"))
      plot_hit_cluster_network(cluster_subg, cid, plot_file, 
                               cluster_index = cidx,  # 传递 cluster_index 以保持与 grid 版本一致
                               width = plot_width, height = plot_height, dpi = plot_dpi)
    }
    
    # 记录聚类信息
    cluster_info_list[[i]] <- data.frame(
      cluster_id = cid,
      node_count = cluster_stats$node_count[i],
      hit_node_count = cluster_stats$hit_node_count[i],
      edge_count = cluster_stats$edge_count[i],
      edge_density = round(cluster_stats$edge_density[i], 4),
      cluster_index = cidx,
      table_file = paste0(safe_cluster_name, "_*.tsv"),
      plot_file = ifelse(is.null(plot_file), "", basename(plot_file)),
      stringsAsFactors = FALSE
    )
  }
  
  cluster_info <- do.call(rbind, cluster_info_list)
  
  # 保存聚类汇总信息
  summary_file <- file.path(output_dir, "cluster_summary.tsv")
  utils::write.table(cluster_info, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("Exported %d clusters from %d hit nodes to: %s\n", 
              nrow(cluster_info), length(hit_ids), output_dir))
  cat(sprintf("Summary file: %s\n", summary_file))
  
  list(
    cluster_info = cluster_info,
    cluster_tables = cluster_tables,
    cluster_graphs = cluster_graphs,
    output_dir = output_dir
  )
}

#' Plot individual hit cluster network with connection relationships
#'
#' @param g_cluster igraph object for a single cluster
#' @param cluster_id integer, ID of the cluster
#' @param output_file character, path to save the plot
#' @param cluster_index integer, optional, index of the cluster (for consistent labeling with grid plot)
#' @param width numeric, plot width
#' @param height numeric, plot height
#' @param dpi numeric, plot DPI
#' @return ggplot object (invisibly)
plot_hit_cluster_network <- function(g_cluster, cluster_id, output_file,
                                     cluster_index = NULL,
                                     width = 10, height = 10, dpi = 300) {
  if (igraph::vcount(g_cluster) <= 1) {
    warning(sprintf("Cluster %d has too few nodes for layout", cluster_id))
    return(NULL)
  }
  
  # 使用力导向布局，保持连接关系
  layout_mat <- igraph::layout_with_fr(g_cluster, 
                                       weights = igraph::E(g_cluster)$weight,
                                       niter = 1000)
  
  # 准备绘图数据（与 plot_hit_clusters_grid 保持一致）
  node_df <- data.frame(
    x = layout_mat[, 1],
    y = layout_mat[, 2],
    name = igraph::V(g_cluster)$name,
    node_index = if ("node_index" %in% igraph::vertex_attr_names(g_cluster)) {
      igraph::V(g_cluster)$node_index
    } else {
      seq_len(igraph::vcount(g_cluster))
    },
    degree = igraph::degree(g_cluster),
    is_hit = if ("is_hit" %in% igraph::vertex_attr_names(g_cluster)) {
      igraph::V(g_cluster)$is_hit
    } else {
      rep(TRUE, igraph::vcount(g_cluster))
    },
    stringsAsFactors = FALSE
  )
  
  # 获取 subclade 属性（与 plot_hit_clusters_grid 保持一致）
  if ("subclade" %in% igraph::vertex_attr_names(g_cluster)) {
    node_df$subclade <- as.character(igraph::V(g_cluster)$subclade)
  }
  if ("clade" %in% igraph::vertex_attr_names(g_cluster)) {
    node_df$clade <- as.character(igraph::V(g_cluster)$clade)
  }
  
  # 若缺少 subclade 属性，统一归为 Others（与 plot_hit_clusters_grid 保持一致）
  if (!"subclade" %in% colnames(node_df)) {
    node_df$subclade <- "Others"
  } else {
    node_df$subclade[is.na(node_df$subclade) | node_df$subclade == "" |
                       tolower(node_df$subclade) == "unknown"] <- "Others"
  }
  
  # 节点大小计算（与 plot_hit_clusters_grid 保持一致）
  node_df$size <- (1.5 + sqrt(node_df$degree) * 1.0) * 2.0
  
  # 为 subclade 构建调色板（与 plot_hit_clusters_grid 使用相同的默认调色板）
  sub_levels <- sort(unique(node_df$subclade))
  default_palette <- c("#E66F73", "#5c95e0", "#cd7560", "#efbed6", "#5D6193",
                      "#b6d37f", "#A9917E", "#9CB79F", "#ffcead", "#589336", "#CCCCCC")
  if (length(sub_levels) <= length(default_palette)) {
    sub_colors <- default_palette[seq_along(sub_levels)]
  } else {
    sub_colors <- rep(default_palette, length.out = length(sub_levels))
  }
  names(sub_colors) <- sub_levels
  
  edge_list <- igraph::as_edgelist(g_cluster, names = TRUE)
  edge_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    weight = igraph::E(g_cluster)$weight,
    stringsAsFactors = FALSE
  )
  
  # 添加坐标
  edge_df$x <- node_df$x[match(edge_df$from, node_df$name)]
  edge_df$y <- node_df$y[match(edge_df$from, node_df$name)]
  edge_df$xend <- node_df$x[match(edge_df$to, node_df$name)]
  edge_df$yend <- node_df$y[match(edge_df$to, node_df$name)]
  
  # 创建标题（与 plot_hit_clusters_grid 的样式保持一致）
  # 使用数字编号（Cluster 1, Cluster 2, ...）以保持与所有输出一致
  if (!is.null(cluster_index)) {
    title_text <- sprintf("Cluster %d", cluster_index)
  } else {
    title_text <- paste("Cluster", cluster_id)
  }
  
  # 绘制网络图（与 plot_hit_clusters_grid 的样式保持一致）
  gg <- ggplot2::ggplot() +
    # 绘制边（连接线）- 与 grid 版本样式一致
    ggplot2::geom_segment(data = edge_df,
                         ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                         color = scales::alpha("#808080", 0.35), 
                         linewidth = 0.25) +
    # 绘制节点：按 subclade 上色（与 grid 版本一致）
    ggplot2::geom_point(data = node_df,
                      ggplot2::aes(x = x, y = y, size = size, fill = subclade),
                      color = "white",
                      alpha = 0.9,
                      shape = 21,
                      stroke = 0.5) +
    # 在节点中心添加编号标签（与 grid 版本一致）
    ggplot2::geom_text(data = node_df,
                      ggplot2::aes(x = x, y = y, label = node_index),
                      color = "#222222",
                      size = 2.2,
                      vjust = 0.4) +
    ggplot2::scale_fill_manual(values = sub_colors, name = "Subclade") +
    ggplot2::scale_size(range = c(1.5, 6), guide = "none") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = title_text) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 8)),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::coord_fixed()
  
  # 保存图片（与 plot_hit_clusters_grid 的样式保持一致，不再显示额外的 ID 列表）
  ggplot2::ggsave(filename = output_file, plot = gg, width = width, height = height, dpi = dpi)
  invisible(gg)
}

#' Plot multiple hit clusters in a tiled grid layout
#'
#' 将 `export_hit_clusters()` 输出的多个聚类块组合成一个多面板图：
#' - 1 个聚类时占满整张画布
#' - 2–4 个聚类排成 2x2 网格
#' - 5–9 个聚类排成 3x3 网格，以此类推（取最小的正方形网格覆盖所有聚类）
#'
#' @param export_result list, `export_hit_clusters()` 的返回结果
#' @param clusters 可选，指定要绘制的聚类 ID 向量（与 `cluster_info$cluster_id` 对应）
#' @param n_max 整数，最多绘制的聚类数（默认 16，防止面板过多）
#' @param output_file 可选，若提供则保存为该路径的图片文件
#' @param width,height,dpi 保存图片时的尺寸与分辨率
#' @param title 总标题
#' @param show_node_index 逻辑值，是否在每个节点中心显示该聚类内的节点编号（默认 TRUE）
#' @param color_palette 可选，颜色向量（如 c("#E66F73", "#5c95e0", "#cd7560")），函数会按顺序自动分配给不同的 subclade。若为 NULL，使用默认调色板
#' @return ggplot 对象（多面板图）
#' @export
plot_hit_clusters_grid <- function(export_result,
                                   clusters = NULL,
                                   n_max = 16,
                                   output_file = NULL,
                                   width = 12, height = 12, dpi = 300,
                                   title = "Hit clusters (combined)",
                                   show_node_index = TRUE,
                                   color_palette = NULL) {
  # 基本检查
  if (is.null(export_result) ||
      !is.list(export_result) ||
      !all(c("cluster_info", "cluster_graphs") %in% names(export_result))) {
    stop("export_result must be the list returned by export_hit_clusters().")
  }
  cluster_info <- export_result$cluster_info
  cluster_graphs <- export_result$cluster_graphs
  if (is.null(cluster_info) || nrow(cluster_info) == 0) {
    stop("export_result contains no clusters to plot.")
  }
  # 若缺少 cluster_index（兼容旧结果），按当前顺序补一个
  if (!"cluster_index" %in% colnames(cluster_info)) {
    cluster_info$cluster_index <- seq_len(nrow(cluster_info))
  }
  
  # 选择要绘制的聚类 ID（保持 cluster_info 中的顺序，通常已按边数排序）
  all_ids <- as.character(cluster_info$cluster_id)
  if (is.null(clusters)) {
    cluster_ids <- all_ids
  } else {
    clusters_chr <- as.character(clusters)
    keep <- clusters_chr %in% all_ids
    if (!any(keep)) {
      stop("None of the requested clusters are present in export_result$cluster_info$cluster_id.")
    }
    cluster_ids <- clusters_chr[keep]
  }
  
  if (!is.null(n_max) && length(cluster_ids) > n_max) {
    cluster_ids <- cluster_ids[seq_len(n_max)]
  }
  
  k <- length(cluster_ids)
  if (k == 0) {
    stop("No clusters selected for plotting.")
  }
  
  # 面板布局：1 个聚类 → 1x1；否则取最小正方形（2x2, 3x3, ...）
  if (k == 1L) {
    ncol <- 1L
    } else {
    side <- ceiling(sqrt(k))
    ncol <- side
  }
  
  # 为 cluster_id -> cluster_index 构建映射，便于在小图中显示连续编号
  index_map <- stats::setNames(cluster_info$cluster_index, as.character(cluster_info$cluster_id))
  
  # 收集所有聚类的节点和边数据
  node_dfs <- list()
  edge_dfs <- list()
  panel_labels <- character(0)
  
  for (i in seq_along(cluster_ids)) {
    cid <- as.character(cluster_ids[i])
    cidx <- if (!is.null(index_map[[cid]])) index_map[[cid]] else i
    g_cluster <- cluster_graphs[[cid]]
    if (is.null(g_cluster) || !igraph::is_igraph(g_cluster) || igraph::vcount(g_cluster) <= 1) {
      next
    }
    
    # 计算布局（力导向）
    layout_mat <- igraph::layout_with_fr(
      g_cluster,
      weights = igraph::E(g_cluster)$weight,
      niter = 800
    )
    
    # 节点数据
    node_df <- data.frame(
      x = layout_mat[, 1],
      y = layout_mat[, 2],
      name = igraph::V(g_cluster)$name,
      node_index = if ("node_index" %in% igraph::vertex_attr_names(g_cluster)) {
        igraph::V(g_cluster)$node_index
    } else {
        seq_len(igraph::vcount(g_cluster))
      },
      degree = igraph::degree(g_cluster),
      is_hit = if ("is_hit" %in% igraph::vertex_attr_names(g_cluster)) {
        igraph::V(g_cluster)$is_hit
      } else {
        rep(TRUE, igraph::vcount(g_cluster))
      },
      stringsAsFactors = FALSE
    )
    
    # 额外属性（若存在，可用于后续扩展）
    if ("subclade" %in% igraph::vertex_attr_names(g_cluster)) {
      node_df$subclade <- as.character(igraph::V(g_cluster)$subclade)
    }
    if ("clade" %in% igraph::vertex_attr_names(g_cluster)) {
      node_df$clade <- as.character(igraph::V(g_cluster)$clade)
    }
    
    # 若缺少 subclade 属性，统一归为 Others，便于按颜色区分
    if (!"subclade" %in% colnames(node_df)) {
      node_df$subclade <- "Others"
    } else {
      node_df$subclade[is.na(node_df$subclade) | node_df$subclade == "" |
                         tolower(node_df$subclade) == "unknown"] <- "Others"
    }
    
    # 美化：根据度数计算节点大小（高度连接节点更突出）
    node_df$size <- (1.5 + sqrt(node_df$degree) * 1.0) * 2.0
    
    # 聚类标签：用于 facet 标题
    # 使用数字编号（Cluster 1, Cluster 2, ...）以保持与所有输出一致
    panel_label <- sprintf("Cluster %d", cidx)
    node_df$cluster_panel <- panel_label
    panel_labels <- c(panel_labels, panel_label)
    
    # 边数据
    el <- igraph::as_edgelist(g_cluster, names = TRUE)
  edge_df <- data.frame(
      from = el[, 1],
      to = el[, 2],
    weight = igraph::E(g_cluster)$weight,
    stringsAsFactors = FALSE
  )
  edge_df$x <- node_df$x[match(edge_df$from, node_df$name)]
  edge_df$y <- node_df$y[match(edge_df$from, node_df$name)]
  edge_df$xend <- node_df$x[match(edge_df$to, node_df$name)]
  edge_df$yend <- node_df$y[match(edge_df$to, node_df$name)]
    edge_df$cluster_panel <- panel_label
    
    node_dfs[[length(node_dfs) + 1L]] <- node_df
    edge_dfs[[length(edge_dfs) + 1L]] <- edge_df
  }
  
  if (length(node_dfs) == 0L || length(edge_dfs) == 0L) {
    stop("No valid clusters with >=2 nodes to plot.")
  }
  
  node_all <- do.call(rbind, node_dfs)
  edge_all <- do.call(rbind, edge_dfs)
  
  # 统一因子顺序，保持 panel 顺序稳定
  panel_labels <- unique(panel_labels)
  node_all$cluster_panel <- factor(node_all$cluster_panel, levels = panel_labels)
  edge_all$cluster_panel <- factor(edge_all$cluster_panel, levels = panel_labels)
  
  # 为 subclade 构建调色板：如果用户提供了颜色向量，按顺序自动分配；否则使用默认调色板
  sub_levels <- sort(unique(node_all$subclade))
  
  # 默认调色板
  default_palette <- c("#E66F73", "#5c95e0", "#cd7560", "#efbed6", "#5D6193",
                      "#b6d37f", "#A9917E", "#9CB79F", "#ffcead", "#589336", "#CCCCCC")
  
  if (!is.null(color_palette) && is.vector(color_palette) && length(color_palette) > 0) {
    # 用户提供了颜色向量：验证颜色格式并循环使用
    # 移除命名（如果有），只保留颜色值
    if (!is.null(names(color_palette))) {
      color_palette <- unname(color_palette)
    }
    
    # 验证颜色值格式（简单检查是否为有效的十六进制颜色或 R 颜色名称）
    valid_colors <- sapply(color_palette, function(x) {
      is.character(x) && (
        grepl("^#[0-9A-Fa-f]{6}$", x) || 
        grepl("^#[0-9A-Fa-f]{3}$", x) ||
        x %in% colors()  # R 内置颜色名称
      )
    })
    
    if (!all(valid_colors)) {
      warning("Some colors in color_palette are not valid. Invalid entries will be skipped.")
      color_palette <- color_palette[valid_colors]
    }
    
    if (length(color_palette) == 0) {
      warning("No valid colors in color_palette. Using default palette.")
      color_palette <- default_palette
    } else {
      # 使用用户提供的颜色向量，按顺序分配给 sub_levels（循环使用）
      sub_colors <- rep(color_palette, length.out = length(sub_levels))
      names(sub_colors) <- sub_levels
    }
  } else {
    # 使用默认调色板
    if (length(sub_levels) <= length(default_palette)) {
      sub_colors <- default_palette[seq_along(sub_levels)]
    } else {
      sub_colors <- rep(default_palette, length.out = length(sub_levels))
    }
    names(sub_colors) <- sub_levels
  }
  
  gg <- ggplot2::ggplot() +
    # 边
    ggplot2::geom_segment(
      data = edge_all,
                         ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = scales::alpha("#808080", 0.35),
      linewidth = 0.25
    ) +
    # 节点：按 subclade 上色
    ggplot2::geom_point(
      data = node_all,
      ggplot2::aes(x = x, y = y, size = size, fill = subclade),
      color = "white",
      alpha = 0.9,
      shape = 21,
      stroke = 0.5
    ) +
    # 可选：在节点中心添加编号标签（每个聚类内的局部编号）
    { if (isTRUE(show_node_index)) ggplot2::geom_text(
        data = node_all,
        ggplot2::aes(x = x, y = y, label = node_index),
        color = "#222222",
        size = 2.2,
        vjust = 0.4
      ) else NULL } +
    ggplot2::scale_fill_manual(values = sub_colors, name = "Subclade") +
    ggplot2::scale_size(
      range = c(1.5, 6),
      guide = "none"
    ) +
    # 每个子图使用各自的坐标范围，使得布局在各 panel 内充分展开
    ggplot2::facet_wrap(~cluster_panel, scales = "free", ncol = ncol) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      # 保持每个小图接近正方形，布局更均衡
      aspect.ratio = 1,
      strip.text = ggplot2::element_text(size = 9, face = "bold", margin = ggplot2::margin(b = 4)),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 8)),
                 legend.position = "right",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      panel.spacing = grid::unit(6, "mm"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(title = title)
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file, plot = gg, width = width, height = height, dpi = dpi)
  }
  
  gg
}

#' Export species-by-module matrix for heatmap visualization
#'
#' 从 `export_hit_clusters()` 的结果中提取物种（species）× 模块（cluster）的节点数量矩阵，
#' 用于后续绘制热图。
#'
#' @param export_result list, `export_hit_clusters()` 的返回结果
#' @param species_column character, 节点表中用于标识物种的列名（默认 "clade" 或从节点属性推断）
#' @param output_file 可选，若提供则保存矩阵为 TSV 文件
#' @return list with:
#'   - matrix: 物种 × 模块的节点数量矩阵（data.frame，行为物种，列为模块）
#'   - long_data: 长格式数据（data.frame，包含 species, module, count 列）
#'   - species_order: 物种顺序（字符向量）
#'   - module_order: 模块顺序（字符向量）
#' @export
export_heatmap_data <- function(export_result,
                                species_column = NULL,
                                output_file = NULL) {
  # 基本检查
  if (is.null(export_result) ||
      !is.list(export_result) ||
      !all(c("cluster_info", "cluster_tables") %in% names(export_result))) {
    stop("export_result must be the list returned by export_hit_clusters().")
  }
  
  cluster_info <- export_result$cluster_info
  cluster_tables <- export_result$cluster_tables
  
  if (is.null(cluster_info) || nrow(cluster_info) == 0) {
    stop("export_result contains no clusters.")
  }
  
  # 收集所有节点的物种信息
  all_nodes_list <- list()
  
  for (i in seq_len(nrow(cluster_info))) {
    cid <- as.character(cluster_info$cluster_id[i])
    cidx <- if ("cluster_index" %in% colnames(cluster_info)) {
      cluster_info$cluster_index[i]
  } else {
      i
    }
    
    if (!cid %in% names(cluster_tables)) {
      next
    }
    
    node_table <- cluster_tables[[cid]]$nodes
    
    if (is.null(node_table) || nrow(node_table) == 0) {
      next
    }
    
    # 确定物种列
    if (is.null(species_column)) {
      # 自动推断：优先 clade，其次 subclade
      if ("clade" %in% colnames(node_table)) {
        species_col <- "clade"
      } else if ("subclade" %in% colnames(node_table)) {
        species_col <- "subclade"
  } else {
        warning(sprintf("No species column found in cluster %s. Using 'Unknown'.", cid))
        node_table$species <- "Unknown"
        species_col <- "species"
      }
    } else {
      if (!species_column %in% colnames(node_table)) {
        warning(sprintf("Column '%s' not found in cluster %s. Trying to infer.", species_column, cid))
        if ("clade" %in% colnames(node_table)) {
          species_col <- "clade"
        } else if ("subclade" %in% colnames(node_table)) {
          species_col <- "subclade"
        } else {
          node_table$species <- "Unknown"
          species_col <- "species"
        }
      } else {
        species_col <- species_column
      }
    }
    
    # 提取节点信息
    # 使用数字编号格式（Cluster_1, Cluster_2, ...）以保持与所有输出一致
    module_label <- sprintf("Cluster_%d", cidx)
    
    node_subset <- data.frame(
      id = node_table$id,
      species = as.character(node_table[[species_col]]),
      module = module_label,
      cluster_id = cid,
      cluster_index = cidx,
      stringsAsFactors = FALSE
    )
    
    # 处理缺失值
    node_subset$species[is.na(node_subset$species) | node_subset$species == ""] <- "Unknown"
    
    all_nodes_list[[length(all_nodes_list) + 1L]] <- node_subset
  }
  
  if (length(all_nodes_list) == 0) {
    stop("No node data found in cluster_tables.")
  }
  
  # 合并所有节点数据
  all_nodes <- do.call(rbind, all_nodes_list)
  
  # 构建物种 × 模块的计数矩阵
  count_matrix <- stats::xtabs(~ species + module, data = all_nodes)
  count_df <- as.data.frame.matrix(count_matrix)
  
  # 转换为数值矩阵（保持行列名）
  count_matrix_final <- as.matrix(count_df)
  
  # 按 cluster_index 对模块列进行排序（确保与 plot_hit_clusters_grid 的顺序一致）
  # 从模块名称中提取 cluster_index（统一使用数字编号格式）
  extract_cluster_index <- function(module_name) {
    if (grepl("^Cluster_([0-9]+)$", module_name)) {
      return(as.integer(gsub("^Cluster_([0-9]+)$", "\\1", module_name)))
    } else if (grepl("^Cluster ([0-9]+)$", module_name)) {
      return(as.integer(gsub("^Cluster ([0-9]+)$", "\\1", module_name)))
    } else {
      return(Inf)  # 无法识别的格式放在最后
    }
  }
  
  module_indices <- sapply(colnames(count_matrix_final), extract_cluster_index)
  module_order_sorted <- colnames(count_matrix_final)[order(module_indices)]
  count_matrix_final <- count_matrix_final[, module_order_sorted, drop = FALSE]
  count_df <- count_df[, module_order_sorted, drop = FALSE]
  
  # 获取排序后的物种和模块顺序
  species_order <- rownames(count_matrix_final)
  module_order <- colnames(count_matrix_final)
  
  # 创建长格式数据（便于后续 ggplot2 绘图）
  # 使用 base R 方法，避免依赖 reshape2
  long_data <- expand.grid(
    species = rownames(count_matrix_final),
    module = colnames(count_matrix_final),
    stringsAsFactors = FALSE
  )
  long_data$count <- as.vector(count_matrix_final[cbind(
    match(long_data$species, rownames(count_matrix_final)),
    match(long_data$module, colnames(count_matrix_final))
  )])
  
  # 保存矩阵文件（如果指定）
  if (!is.null(output_file)) {
    # 添加物种列（作为第一列）
    count_df_out <- data.frame(
      species = rownames(count_df),
      count_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    utils::write.table(count_df_out, file = output_file, sep = "\t", 
                      quote = FALSE, row.names = FALSE)
    cat(sprintf("Heatmap data matrix saved to: %s\n", output_file))
  }
  
  list(
    matrix = count_df,
    matrix_numeric = count_matrix_final,
    long_data = long_data,
    species_order = species_order,
    module_order = module_order
  )
}

