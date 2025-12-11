#' Prepare plotting data from layout and graph
#'
#' @param g_kept igraph of kept nodes
#' @param layout_kept matrix
#' @return list with node.data, edge.data
#' @export
prepare_plot_data <- function(g_kept, layout_kept, node_groups = NULL, module_by = c("subclade","clade","community")) {
  module_by <- match.arg(module_by)
  rownames(layout_kept) <- igraph::V(g_kept)$name
  edge_list <- igraph::as_edgelist(g_kept, names = TRUE)
  edge.data <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    weight = igraph::E(g_kept)$weight,
    stringsAsFactors = FALSE
  )
  # 模块优先从 node_groups 精确匹配，再回退到图属性
  v_names <- igraph::V(g_kept)$name
  module_vec <- NULL
  clade_vec <- igraph::V(g_kept)$clade
  clade_vec <- as.character(clade_vec)
  if (!is.null(node_groups) && is.data.frame(node_groups)) {
    if ("id" %in% colnames(node_groups)) node_groups$id <- as.character(node_groups$id)
    want_field <- switch(module_by,
                         subclade = if ("subclade" %in% colnames(node_groups)) "subclade" else NULL,
                         clade = if ("clade" %in% colnames(node_groups)) "clade" else NULL,
                         community = NULL)
    if (!is.null(want_field)) {
      tmp <- node_groups[match(v_names, node_groups$id), , drop = FALSE]
      if (nrow(tmp) > 0 && want_field %in% colnames(tmp)) {
        module_vec <- as.character(tmp[[want_field]])
      }
    }
  }
  if (is.null(module_vec)) {
    module_vec <- switch(module_by,
                         subclade = igraph::V(g_kept)$subclade,
                         clade = igraph::V(g_kept)$clade,
                         community = igraph::V(g_kept)$community)
    module_vec <- as.character(module_vec)
  }
  if (module_by == "subclade") {
    sub_vec <- module_vec
    sub_vec[is.na(sub_vec) | sub_vec == ""] <- clade_vec[is.na(sub_vec) | sub_vec == ""]
    sub_vec[is.na(sub_vec) | sub_vec == "" | tolower(sub_vec) == "unknown"] <- "Others"
    module_vec <- sub_vec
  } else if (module_by == "clade") {
    cvec <- module_vec
    cvec[is.na(cvec) | cvec == ""] <- "Unknown"
    cvec[is.na(cvec) | cvec == "" | tolower(cvec) == "unknown"] <- "Others"
    module_vec <- cvec
  } else {
    cmt <- suppressWarnings(as.integer(module_vec))
    cmt[is.na(cmt)] <- 1L
    module_vec <- as.character(cmt)
  }

  node.data <- data.frame(
    pos.x = layout_kept[, 1],
    pos.y = layout_kept[, 2],
    Module = module_vec,
    name = igraph::V(g_kept)$name,
    clade = clade_vec,
    degree = igraph::degree(g_kept),
    stringsAsFactors = FALSE
  )
  # 确保模块中不出现 Unknown，统一显示为 Others
  node.data$Module[node.data$Module == "Unknown" | tolower(node.data$Module) == "unknown" | node.data$Module == ""] <- "Others"
  edge.data$from.x <- node.data$pos.x[match(edge.data$from, node.data$name)]
  edge.data$from.y <- node.data$pos.y[match(edge.data$from, node.data$name)]
  edge.data$to.x <- node.data$pos.x[match(edge.data$to, node.data$name)]
  edge.data$to.y <- node.data$pos.y[match(edge.data$to, node.data$name)]
  # 根据节点的度（连接数）计算节点大小：度越大，节点越大
  # 使用可调参数，确保大小与连接数正相关
  node.data$size <- (1.5 + sqrt(node.data$degree) * 0.5) * 1.5
  list(node.data = node.data, edge.data = edge.data)
}

#' Prepare plotting data with customizable node size scaling
#'
#' @param g_kept igraph of kept nodes
#' @param layout_kept matrix
#' @param node_groups optional node groups data.frame
#' @param module_by module grouping method
#' @param node_size_base base size for nodes with degree 0
#' @param node_size_scale scaling factor for degree contribution
#' @param node_size_multiplier overall size multiplier
#' @param node_size_mode scaling mode for size calculation
#' @return list with node.data, edge.data
#' @export
prepare_plot_data_with_size <- function(g_kept, layout_kept, node_groups = NULL, 
                                        module_by = c("subclade","clade","community"),
                                        node_size_base = 1.5, node_size_scale = 1.2, 
                                        node_size_multiplier = 1.5, node_size_mode = "power075") {
  pd <- prepare_plot_data(g_kept, layout_kept, node_groups = node_groups, module_by = module_by)
  # 根据度重新计算节点大小，使用可调参数
  pd$node.data$size <- calculate_node_size(pd$node.data$degree, 
                                           base_size = node_size_base,
                                           scale_factor = node_size_scale,
                                           size_multiplier = node_size_multiplier,
                                           size_mode = node_size_mode)
  pd
}

#' Calculate node size based on degree with customizable scaling
#'
#' @param degree numeric vector of node degrees
#' @param base_size numeric, base size for nodes with degree 0 (default: 1.5)
#' @param scale_factor numeric, scaling factor for degree contribution (default: 1.2)
#' @param size_multiplier numeric, overall size multiplier (default: 1.5)
#' @param size_mode character, scaling mode: "sqrt" (smooth), "linear" (more sensitive), "power075" (balanced, default)
#' @return numeric vector of node sizes
#' @export
calculate_node_size <- function(degree, base_size = 1.5, scale_factor = 1.2, size_multiplier = 1.5,
                                size_mode = c("power075", "sqrt", "linear")) {
  size_mode <- match.arg(size_mode)
  # 根据模式选择缩放方式，使连接数差异更明显
  if (size_mode == "linear") {
    # 线性缩放：差异最明显，但可能产生极端值
    scaled_degree <- degree * scale_factor
  } else if (size_mode == "sqrt") {
    # 平方根缩放：最平滑，差异较小
    scaled_degree <- sqrt(degree) * scale_factor
  } else {
    # power075 (默认)：度^0.75，平衡差异明显度和平滑度
    scaled_degree <- (degree^0.75) * scale_factor
  }
  (base_size + scaled_degree) * size_multiplier
}

#' Plot Fabaceae network with customizable palette
#'
#' @param node.data data.frame from prepare_plot_data
#' @param edge.data data.frame from prepare_plot_data
#' @param sc_points list of subclade hull points
#' @param palette Optional named vector of colors for Modules; if NULL, uses defaults
#' @param title Plot title
#' @param show_module_boundary logical, whether to show module background boundaries (default: TRUE)
#' @param boundary_alpha numeric, transparency of module boundaries (0-1, default: 0.03, lower = more transparent)
#' @return ggplot object
#' @export
plot_fab_network <- function(node.data, edge.data, sc_points = NULL, palette = NULL,
                             title = "Fabaceae ANK Network",
                             show_module_boundary = TRUE, boundary_alpha = 0.03,
                             labels_map = NULL,
                             id_color_map = NULL, id_color_map_file = NULL,
                             hit_ids = NULL, hit_ids_file = NULL, hit_label = "Selected",
                             hit_color = NULL) {
  # 规范化模块名称中的 Unknown -> Others，确保图例显示一致
  node.data$Module[node.data$Module == "Unknown" | tolower(node.data$Module) == "unknown" | node.data$Module == ""] <- "Others"
  unique_modules <- sort(unique(node.data$Module))
  module_counts <- sort(table(node.data$Module), decreasing = TRUE)
  ordered_modules <- names(module_counts)
  num_mod <- length(ordered_modules)
  default_palette <- c("#E66F73", "#5c95e0", "#cd7560", "#efbed6", "#5D6193",
                       "#b6d37f", "#A9917E", "#9CB79F", "#ffcead", "#589336", "#CCCCCC")
  if (is.null(palette)) {
    if (num_mod <= length(default_palette)) {
      assigned <- c(default_palette[1], default_palette[2:num_mod])
    } else {
      rest <- default_palette[2:length(default_palette)]
      need <- num_mod - 1
      assigned <- c(default_palette[1], rep(rest, length.out = need))
    }
    color_vec <- assigned
    names(color_vec) <- ordered_modules
  } else {
    # palette can be named or unnamed; ensure all modules covered
    if (is.null(names(palette))) {
      if (length(palette) < num_mod) palette <- rep(palette, length.out = num_mod)
      color_vec <- stats::setNames(palette[seq_len(num_mod)], ordered_modules)
    } else {
      color_vec <- palette
      missing <- setdiff(ordered_modules, names(color_vec))
      if (length(missing) > 0) {
        extras <- rep(default_palette, length.out = length(missing))
        names(extras) <- missing
        color_vec <- c(color_vec, extras)
      }
      color_vec <- color_vec[ordered_modules]
    }
  }
  # 处理按 id 强制配色覆盖（仅影响节点颜色，不改变模块/边配色逻辑）
  if (!is.null(id_color_map_file) && is.null(id_color_map)) {
    if (file.exists(id_color_map_file)) {
      tmp_map <- tryCatch({
        read.table(id_color_map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      if (!is.null(tmp_map) && ncol(tmp_map) >= 2) {
        # 期望列名: id, color（若无列名则取前两列）
        colnames(tmp_map)[1:2] <- c("id", "color")
        id_color_map <- stats::setNames(as.character(tmp_map$color), as.character(tmp_map$id))
      }
    }
  }
  override_ids <- character(0)
  if (!is.null(id_color_map)) {
    # 命名向量: names = id, values = color
    id_color_map <- id_color_map[!is.na(names(id_color_map)) & nzchar(names(id_color_map))]
    override_ids <- intersect(names(id_color_map), node.data$name)
  }

  # 边的模块保持基于原始 Module（不受覆盖影响）
  edge.data$from.module <- node.data$Module[match(edge.data$from, node.data$name)]
  edge.data$to.module <- node.data$Module[match(edge.data$to, node.data$name)]
  edge.data$same_module <- edge.data$from.module == edge.data$to.module
  edge.data$edge.color <- ifelse(edge.data$same_module, color_vec[edge.data$from.module],
                                 scales::alpha("#808080", 0.25))
  gg <- ggplot2::ggplot()
  gg <- gg + ggplot2::geom_curve(mapping = ggplot2::aes(x = from.x, y = from.y, xend = to.x, yend = to.y, color = edge.color),
                                  curvature = 0.30, linewidth = 0.216, lineend = "round", data = edge.data, ncp = 5) +
    ggplot2::scale_color_identity(guide = "none")
  if (!is.null(sc_points) && show_module_boundary) {
    for (sc in names(sc_points)) {
      pts <- sc_points[[sc]]
      # 统一成矩阵并清洗无效/重复点
      if (!is.null(pts)) {
        pts <- as.matrix(pts)
        if (ncol(pts) >= 2 && nrow(pts) >= 1) {
          keep_ok <- is.finite(pts[, 1]) & is.finite(pts[, 2])
          pts <- pts[keep_ok, , drop = FALSE]
          if (nrow(pts) > 1) {
            pts <- unique(pts)
          }
        }
      }
      if (!is.null(pts) && is.matrix(pts) && nrow(pts) >= 3 && ncol(pts) >= 2) {
        hull_idx <- grDevices::chull(pts)
        if (length(hull_idx) >= 3 && max(hull_idx) <= nrow(pts)) {
          hull_idx <- c(hull_idx, hull_idx[1])
          this_col <- color_vec[sc]
          if (is.na(this_col)) this_col <- "#CCCCCC"
          col_fill <- scales::alpha(this_col, boundary_alpha)
          # 构造独立数据框，避免在 aes 中引用外部变量导致评估越界
          hull_df <- data.frame(x = pts[hull_idx, 1], y = pts[hull_idx, 2])
          gg <- gg + ggplot2::geom_polygon(
            data = hull_df,
            mapping = ggplot2::aes(x = x, y = y),
            fill = col_fill, color = NA, linewidth = 0
          )
        }
      }
    }
  }
  # 读取命中ID列表（单列），命中的点作为一个类别一起配色（不改布局）
  if (!is.null(hit_ids_file) && is.null(hit_ids)) {
    if (file.exists(hit_ids_file)) {
      tmp_ids <- tryCatch({
        read.table(hit_ids_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      if (!is.null(tmp_ids) && ncol(tmp_ids) >= 1) {
        hit_ids <- as.character(tmp_ids[[1]])
      }
    }
  }
  # 保持原行为：仅接受向量/字符；若传 data.frame 请手动传第一列字符向量
  hit_ids <- unique(as.character(hit_ids))
  hit_ids <- hit_ids[!is.na(hit_ids) & nzchar(hit_ids)]

  # 为节点构造显示用模块（仅用于颜色图例）
  node.data$ModuleDisplay <- node.data$Module
  # 1) 覆盖颜色（每个命中的id独立一个类别，保持之前逻辑）
  if (length(override_ids) > 0) {
    node.data$ModuleDisplay[node.data$name %in% override_ids] <- paste0("__ID_", node.data$name[node.data$name %in% override_ids])
  }
  # 2) 命中ID列表：统一归入一个类别 hit_label（如“Selected”）
  if (length(hit_ids) > 0) {
    node.data$ModuleDisplay[node.data$name %in% hit_ids] <- hit_label
  }

  gg <- gg + ggplot2::geom_point(mapping = ggplot2::aes(x = pos.x, y = pos.y, fill = ModuleDisplay, size = size),
                                  data = node.data, shape = 21, color = "#000000", stroke = 0.5, alpha = 0.6)
  label_df <- stats::aggregate(cbind(pos.x, pos.y) ~ Module, data = node.data, FUN = mean)
  # 支持外部强制标签映射（命名向量：Module -> label）
  if (!is.null(labels_map)) {
    label_df$label <- ifelse(label_df$Module %in% names(labels_map),
                             as.character(labels_map[label_df$Module]),
                             label_df$Module)
  } else {
    label_df$label <- label_df$Module
  }
  # 将默认的 "Unknown" 标签显示为 "Others"（仅影响显示文本，不改变分组/配色）
  label_df$label[label_df$label == "Unknown" | tolower(label_df$label) == "unknown" | label_df$label == ""] <- "Others"
  gg <- gg + ggrepel::geom_label_repel(data = label_df,
                                       mapping = ggplot2::aes(x = pos.x, y = pos.y, label = label),
                                       size = 6.4, label.size = 0, fill = scales::alpha('#FFFFFF', 0.6),
                                       color = '#333333', max.overlaps = Inf, seed = 123,
                                       box.padding = 0.3, point.padding = 0.2)
  gg <- gg + ggplot2::scale_size(range = c(0.8, 5))
  gg <- gg + ggplot2::theme_void()
  gg <- gg + ggplot2::labs(x = "", y = "", title = title)
  # 扩展颜色映射：为命中的 id 指派自定义颜色
  if (length(override_ids) > 0) {
    override_labels <- paste0("__ID_", override_ids)
    extra_cols <- as.character(id_color_map[override_ids])
    names(extra_cols) <- override_labels
    color_vec <- c(color_vec, extra_cols)
  }
  # 若存在命中列表类别，且提供了 hit_color，则为该类别指定颜色
  if (length(hit_ids) > 0 && !is.null(hit_color)) {
    color_vec[hit_label] <- hit_color
  }
  gg <- gg + ggplot2::scale_fill_manual(values = color_vec)
  # 放大右侧图例示例标记大小
  gg <- gg + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 6)))
  gg <- gg + ggplot2::theme(legend.key.size = grid::unit(8, "mm"))
  gg <- gg + ggplot2::coord_fixed()
  gg
}

#' One-shot pipeline to build and plot Fabaceae network
#'
#' @param edge_file edges tsv
#' @param node_group_file nodes tsv
#' @param cluster_method clustering method
#' @param cluster_params list of extra params passed to clustering
#' @param palette color palette mapping or vector
#' @param output_file optional path to save PNG
#' @param width width in inches for saving
#' @param height height in inches for saving
#' @param dpi dpi for saving
#' @param show_module_boundary logical, whether to show module background boundaries (default: TRUE)
#' @param boundary_alpha numeric, transparency of module boundaries (0-1, default: 0.03)
#' @return list with gg (plot), node.data, edge.data
#' @export
fabnetsyn_build_plot <- function(edge_file, node_group_file,
                                 cluster_method = c("infomap","louvain","walktrap","leiden"),
                                 cluster_params = list(),
                                 palette = NULL, output_file = NULL,
                                 width = 12, height = 12, dpi = 300,
                                 show_module_boundary = TRUE, boundary_alpha = 0.03,
                                 labels_map = NULL,
                                 filter_subclade = NULL,   # 新增：按 subclade 过滤
                                 edge_filter = c("any", "both"),  # 新增：边筛选策略
                                 # 命中ID统一上色（单列向量/表格），不改布局
                                 hit_ids = NULL, hit_ids_file = NULL, hit_label = "Selected", hit_color = NULL,
                                 # 单点独立颜色覆盖（命名向量/两列表格）
                                 id_color_map = NULL, id_color_map_file = NULL,
                                 # 节点大小控制：根据连接数（度）调整大小
                                 node_size_base = 1.5,      # 基础大小（度=0时的节点大小）
                                 node_size_scale = 1.2,      # 度对大小的缩放因子（越大，连接数影响越大）
                                 node_size_multiplier = 1.5, # 整体大小倍数
                                 node_size_mode = "power075", # 缩放模式："power075"(平衡，默认), "linear"(差异最明显), "sqrt"(最平滑)
                                 ...) { # 吸收多余参数，避免旧会话报 未用参数
  edge_filter <- match.arg(edge_filter)
  dat <- read_network_tables(edge_file, node_group_file)
  # 若指定了 subclade 过滤，则先筛选边与节点
  if (!is.null(filter_subclade)) {
    ng <- dat$node_groups
    if (!"id" %in% colnames(ng)) {
      stop("node_group_file 缺少 id 列，无法按 subclade 过滤")
    }
    # 标准化类型
    ng$id <- as.character(ng$id)
    if (!"subclade" %in% colnames(ng)) {
      stop("node_group_file 缺少 subclade 列，无法按 subclade 过滤")
    }
    # 选择目标 subclade 的 id
    target_ids <- ng$id[!is.na(ng$subclade) & ng$subclade %in% filter_subclade]
    target_ids <- unique(as.character(target_ids))
    if (length(target_ids) == 0) {
      stop("按给定 subclade 未匹配到任何 id，请检查参数与 node_group_file")
    }
    nd <- dat$network_data
    nd$node1 <- as.character(nd$node1)
    nd$node2 <- as.character(nd$node2)
    if (edge_filter == "any") {
      keep <- (nd$node1 %in% target_ids) | (nd$node2 %in% target_ids)
    } else {
      keep <- (nd$node1 %in% target_ids) & (nd$node2 %in% target_ids)
    }
    nd2 <- nd[keep, , drop = FALSE]
    if (nrow(nd2) == 0) {
      stop("按 subclade 过滤后，edge_file 中没有可用的边记录")
    }
    # 更新节点分组，仅保留出现在边中的节点
    kept_ids <- unique(c(nd2$node1, nd2$node2))
    ng2 <- ng[ng$id %in% kept_ids, , drop = FALSE]
    if (nrow(ng2) == 0) {
      stop("按 subclade 过滤后，node_group_file 中没有与边对应的 id")
    }
    dat$network_data <- nd2
    dat$node_groups <- ng2
  }
  g <- build_fab_graph(dat$network_data, dat$node_groups)
  g <- do.call(cluster_graph, c(list(g = g, method = match.arg(cluster_method)), cluster_params))
  lay <- compute_subclade_layouts(g)
  circ <- assemble_circular_layout(g, lay$sub_layouts, lay$sub_keep_indices, lay$weight_df)
  kept_idx <- which(circ$keep_global)
  if (length(kept_idx) == 0) stop("No nodes left after outlier filtering. Adjust thresholds.")
  g_kept <- igraph::induced_subgraph(g, kept_idx)
  layout_kept <- circ$final_layout[kept_idx, , drop = FALSE]
  rownames(layout_kept) <- igraph::V(g_kept)$name
  adj <- adjust_outer_ring(layout_kept, circ$target_edge_radius)
  layout_kept <- adj$layout_kept
  layout_kept <- force_optimize_layout(g_kept, layout_kept, circ$target_edge_radius,
                                       min_distance = (circ$outer_radius - circ$inner_radius) * 0.015)
  # 准备绘图数据，并根据度计算节点大小
  pd <- prepare_plot_data(g_kept, layout_kept, node_groups = dat$node_groups, module_by = "subclade")
  # 根据连接数（度）重新计算节点大小：连接越多，节点越大
  pd$node.data$size <- calculate_node_size(pd$node.data$degree,
                                           base_size = node_size_base,
                                           scale_factor = node_size_scale,
                                           size_multiplier = node_size_multiplier,
                                           size_mode = node_size_mode)
  # 规范化命中ID：允许传入 data.frame/list，且仅以 node_group_file 的 id 为准
  norm_hit_ids <- hit_ids
  if (!is.null(norm_hit_ids)) {
    if (is.data.frame(norm_hit_ids)) {
      if (ncol(norm_hit_ids) >= 1) norm_hit_ids <- norm_hit_ids[[1]]
    }
    if (is.list(norm_hit_ids) && !is.vector(norm_hit_ids)) {
      norm_hit_ids <- unlist(norm_hit_ids, use.names = FALSE)
    }
    norm_hit_ids <- as.character(norm_hit_ids)
    # 只保留在 node_groups$id 中存在的成员
    if ("id" %in% colnames(dat$node_groups)) {
      norm_hit_ids <- intersect(norm_hit_ids, as.character(dat$node_groups$id))
    }
    norm_hit_ids <- unique(norm_hit_ids)
  }
  gg <- plot_fab_network(pd$node.data, pd$edge.data, sc_points = circ$sc_points, palette = palette,
                         title = "Fabaceae ANK Network - FabnetSyn",
                         show_module_boundary = show_module_boundary, boundary_alpha = boundary_alpha,
                         labels_map = labels_map,
                         id_color_map = id_color_map, id_color_map_file = id_color_map_file,
                         hit_ids = norm_hit_ids, hit_ids_file = hit_ids_file, hit_label = hit_label,
                         hit_color = hit_color)
  if (!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file, plot = gg, width = width, height = height, dpi = dpi)
  }
  # 返回图对象，便于后续导出模块
  list(gg = gg, node.data = pd$node.data, edge.data = pd$edge.data, graph = g_kept)
}


