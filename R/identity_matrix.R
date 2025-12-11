#' Plot identity/co-linearity matrix heatmap for node IDs
#'
#' This function builds an N x N matrix where rows/columns are the IDs from
#' the provided node_group_file. The cell intensity encodes similarity/collinearity:
#' - Same feature (e.g., same subclade/clade) increases the score
#' - Collinearity via edges between the pair increases the score
#' - Diagonal (i==j) is set to the maximum score for strongest color
#'
#' Colors: default cold-to-warm palette; more edges -> warmer color
#'
#' @param node_group_file character, path to node groups (must contain column 'id')
#' @param edge_file optional character, path to edges (two or three columns; weight optional)
#' @param feature_by character, which feature to consider as identity: 'subclade' or 'clade'
#' @param group_by optional character, column name used only for ordering/grouping on axes
#' @param weight_same_feature numeric, weight added when two IDs share the same feature (default 1)
#' @param weight_edge numeric, weight added per edge between two IDs (default 1)
#' @param reorder logical, whether to reorder IDs by feature for block structure (default TRUE)
#' @param color_by character, one of 'edge_count' (按共线条数), 'edge_weight' (按权重和), or 'score' (综合分)
#' @param weight_column integer or 'auto', which column in edge_file holds the weight (default 3; use 'auto' to detect)
#' @param weight_aggregate character, how to aggregate multiple edges between two IDs: 'sum' (default), 'mean', or 'max'
#' @param palette vector of 5 hex colors from cold to warm (default provided)
#' @param verbose logical, print diagnostic information to help debug (default FALSE)
#' @param show_axis_labels logical, whether to show axis labels (default FALSE)
#' @param output_file optional character, if provided saves the heatmap to this file
#' @param width numeric, plot width when saving (inches)
#' @param height numeric, plot height when saving (inches)
#' @param dpi numeric, resolution when saving
#' @return list with gg (plot), matrix (used value matrix), ids (row/col order)
#' @export
plot_identity_matrix <- function(node_group_file,
                                 edge_file = NULL,
                                 feature_by = c("subclade", "clade"),
                                 group_by = NULL,
                                 weight_same_feature = 1,
                                 weight_edge = 1,
                                 reorder = FALSE,
                                 color_by = c("edge_count", "edge_weight", "score"),
                                 weight_column = 3,
                                 weight_aggregate = c("sum","mean","max"),
                                 palette = c("#5e878f", "#dbdb99", "#f1a163", "#dd574b", "#a32440"),
                                 show_axis_labels = FALSE,
                                 verbose = FALSE,
                                 output_file = NULL,
                                 width = 10, height = 10, dpi = 300) {
  feature_by <- match.arg(feature_by)
  color_by <- match.arg(color_by)
  weight_aggregate <- match.arg(weight_aggregate)
  
  if (!file.exists(node_group_file)) {
    stop("node_group_file not found: ", node_group_file)
  }
  
  # Read node groups
  ng <- tryCatch({
    read.table(node_group_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "\"", comment.char = "")
  }, error = function(e) {
    stop("Failed to read node_group_file: ", conditionMessage(e))
  })
  
  if (!"id" %in% colnames(ng)) stop("node_group_file must have column 'id'")
  ids <- unique(as.character(ng$id))
  ids <- trimws(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  n <- length(ids)
  if (n == 0) stop("No IDs found in node_group_file column 'id'")
  if (verbose) {
    cat(sprintf("[identity] IDs loaded: %d\n", n))
  }
  
  # Feature vector (for identity score)
  feat <- rep("Others", n)
  if (feature_by %in% colnames(ng)) {
    m <- ng[match(ids, as.character(ng$id)), feature_by, drop = TRUE]
    m <- as.character(m)
    m[is.na(m) | m == "" | tolower(m) == "unknown"] <- "Others"
    feat <- m
  }

  # Grouping feature vector (for ordering on axes), default to feat if group_by is NULL
  group_feat <- feat
  if (!is.null(group_by) && group_by %in% colnames(ng)) {
    g <- ng[match(ids, as.character(ng$id)), group_by, drop = TRUE]
    g <- as.character(g)
    g[is.na(g) | g == "" | tolower(g) == "unknown"] <- "Others"
    group_feat <- g
  }
  
  # Initialize score and edge count matrices
  score <- matrix(0, nrow = n, ncol = n)
  edge_count <- matrix(0, nrow = n, ncol = n)
  edge_weight_sum <- matrix(0, nrow = n, ncol = n)
  edge_weight_max <- matrix(0, nrow = n, ncol = n)
  edge_weight_n <- matrix(0, nrow = n, ncol = n)
  rownames(score) <- colnames(score) <- ids
  rownames(edge_count) <- colnames(edge_count) <- ids
  rownames(edge_weight_sum) <- colnames(edge_weight_sum) <- ids
  rownames(edge_weight_max) <- colnames(edge_weight_max) <- ids
  rownames(edge_weight_n) <- colnames(edge_weight_n) <- ids
  
  # Same feature contribution (to score)
  if (weight_same_feature != 0) {
    split_idx <- split(seq_len(n), feat)
    for (grp in split_idx) {
      if (length(grp) >= 1) {
        score[grp, grp] <- score[grp, grp] + weight_same_feature
      }
    }
  }
  
  # Edge contribution (collinearity)
  if (!is.null(edge_file)) {
    if (!file.exists(edge_file)) stop("edge_file not found: ", edge_file)
    edges <- tryCatch({
      read.table(edge_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "\"", comment.char = "")
    }, error = function(e) {
      stop("Failed to read edge_file: ", conditionMessage(e))
    })
    if (ncol(edges) < 2) stop("edge_file must have at least two columns")
    if (verbose) {
      cat(sprintf("[identity] edge_file read: %d rows, %d cols | color_by=%s | weight_column=%s | aggregate=%s\n",
                  nrow(edges), ncol(edges), color_by, as.character(weight_column), weight_aggregate))
      print(utils::head(edges, 3))
    }
    # 自动检测或按指定列解析权重与节点列
    parse_edges_auto <- FALSE
    if (is.character(weight_column) && identical(tolower(weight_column), "auto")) {
      parse_edges_auto <- TRUE
    }
    if (parse_edges_auto && ncol(edges) >= 3) {
      # 计算每列与 ids 的重叠数量、以及可数值转换程度
      col_vals <- lapply(edges[1:3], function(x) trimws(as.character(x)))
      overlap_counts <- sapply(col_vals, function(v) sum(v %in% ids, na.rm = TRUE))
      numeric_counts <- sapply(col_vals, function(v) sum(is.finite(suppressWarnings(as.numeric(v))), na.rm = TRUE))
      # 选择两个与 ids 重叠最多的列作为节点列
      node_cols <- order(overlap_counts, decreasing = TRUE)[1:2]
      weight_col <- setdiff(1:3, node_cols)
      # 如果权重列的数值性较差，尝试用数值性最强的列
      if (length(weight_col) != 1 || numeric_counts[weight_col] == 0) {
        weight_col <- order(numeric_counts, decreasing = TRUE)[1]
        node_cols <- setdiff(1:3, weight_col)[1:2]
      }
      edges$node1 <- col_vals[[node_cols[1]]]
      edges$node2 <- col_vals[[node_cols[2]]]
      suppressWarnings({ edges$weight <- as.numeric(col_vals[[weight_col]]) })
      if (verbose) {
        cat(sprintf("[identity] auto-detected columns: node1=V%d, node2=V%d, weight=V%d\n", node_cols[1], node_cols[2], weight_col))
        cat(sprintf("[identity] node overlaps: [%d, %d, %d] | numeric counts: [%d, %d, %d]\n",
                    overlap_counts[1], overlap_counts[2], overlap_counts[3], numeric_counts[1], numeric_counts[2], numeric_counts[3]))
      }
    } else {
      # 兼容权重列位置：当 weight_column==2，格式假定为: node1, weight, node2
      # 当 weight_column==3，格式假定为: node1, node2, weight
      if (weight_column == 2) {
        if (ncol(edges) < 3) stop("edge_file needs >=3 columns when weight_column=2 (node1, weight, node2)")
        edges$node1 <- trimws(as.character(edges[[1]]))
        suppressWarnings({ edges$weight <- as.numeric(edges[[2]]) })
        edges$node2 <- trimws(as.character(edges[[3]]))
      } else {
        # 默认第三列为权重（若不存在，则权重=1）
        edges$node1 <- trimws(as.character(edges[[1]]))
        edges$node2 <- trimws(as.character(edges[[2]]))
        if (ncol(edges) >= 3) {
          suppressWarnings({ edges$weight <- as.numeric(edges[[3]]) })
        } else {
          edges$weight <- NA_real_
        }
      }
    }
    edges$weight[is.na(edges$weight)] <- 1
    edges$node1 <- trimws(as.character(edges$node1))
    edges$node2 <- trimws(as.character(edges$node2))
    
    keep <- edges$node1 %in% ids & edges$node2 %in% ids
    if (any(keep)) {
      e2 <- edges[keep, , drop = FALSE]
      if (verbose) {
        cat(sprintf("[identity] edges kept within ID set: %d (of %d)\n", nrow(e2), nrow(edges)))
        cat(sprintf("[identity] weight summary: min=%.4f, max=%.4f\n", min(e2$weight, na.rm = TRUE), max(e2$weight, na.rm = TRUE)))
      }
      # index map
      id_to_idx <- stats::setNames(seq_len(n), ids)
      i_idx <- id_to_idx[e2$node1]
      j_idx <- id_to_idx[e2$node2]
      for (k in seq_len(nrow(e2))) {
        i <- i_idx[k]; j <- j_idx[k]
        if (is.na(i) || is.na(j)) next
        w <- as.numeric(e2$weight[k])
        if (!is.finite(w)) w <- 1
        # 计数
        edge_count[i, j] <- edge_count[i, j] + 1
        edge_count[j, i] <- edge_count[j, i] + 1
        # 权重累计
        edge_weight_sum[i, j] <- edge_weight_sum[i, j] + w
        edge_weight_sum[j, i] <- edge_weight_sum[j, i] + w
        edge_weight_max[i, j] <- max(edge_weight_max[i, j], w)
        edge_weight_max[j, i] <- max(edge_weight_max[j, i], w)
        edge_weight_n[i, j] <- edge_weight_n[i, j] + 1
        edge_weight_n[j, i] <- edge_weight_n[j, i] + 1
        # 综合分（保持原逻辑）
        score[i, j] <- score[i, j] + weight_edge
        score[j, i] <- score[j, i] + weight_edge
      }
    } else {
      if (verbose) cat("[identity] No edges match the provided IDs.\n")
    }
  }
  
  # 选择用于着色的矩阵
  if (!is.null(edge_file) && color_by == "edge_weight") {
    if (weight_aggregate == "sum") {
      value_mat <- edge_weight_sum
    } else if (weight_aggregate == "mean") {
      value_mat <- ifelse(edge_weight_n > 0, edge_weight_sum / edge_weight_n, 0)
    } else { # max
      value_mat <- edge_weight_max
    }
  } else if (!is.null(edge_file) && color_by == "edge_count") {
    value_mat <- edge_count
  } else {
    value_mat <- score
  }

  # If using edge_file, remove IDs with no collinearity (no edges at all, based on edge_count only)
  if (!is.null(edge_file)) {
    # 使用 edge_count（真实共线边计数）来判断是否有共线关系，
    # 避免被同 feature 的得分或其它度量干扰
    deg_edges <- rowSums(edge_count, na.rm = TRUE) + colSums(edge_count, na.rm = TRUE)
    keep <- deg_edges > 0
    if (any(keep) && sum(keep) < length(keep)) {
      if (verbose) cat(sprintf("[identity] Removing no-edge IDs: %d removed, %d kept\n", sum(!keep), sum(keep)))
      value_mat <- value_mat[keep, keep, drop = FALSE]
      ids <- ids[keep]
      feat <- feat[keep]
      group_feat <- group_feat[keep]
    }
    if (length(ids) == 0 || nrow(value_mat) == 0) {
      stop("All IDs have no collinearity edges after filtering; nothing to plot.")
    }
  }
  
  # Compute off-diagonal range for color scaling
  is_offdiag <- row(value_mat) != col(value_mat)
  off_vals <- if (any(is_offdiag)) value_mat[is_offdiag] else as.numeric(value_mat)
  vmax_off <- suppressWarnings(max(off_vals, na.rm = TRUE))
  vmin_off <- suppressWarnings(min(off_vals, na.rm = TRUE))
  if (!is.finite(vmax_off)) vmax_off <- 0
  if (!is.finite(vmin_off)) vmin_off <- 0
  if (vmax_off == vmin_off) {
    vmin_off <- 0
    vmax_off <- max(1, vmax_off)
  }
  if (verbose) {
    cat(sprintf("[identity] off-diagonal range: [%.4f, %.4f]\n", vmin_off, vmax_off))
  }
  
  # Set diagonal to the top of the scale to show strongest warm color
  diag(value_mat) <- vmax_off
  
  # Optional reorder
  ord <- seq_len(n)
  if (reorder) {
    deg <- rowSums(value_mat, na.rm = TRUE)
    # 使用 group_feat 的“首次出现顺序”作为物种/分组顺序，而不是字母顺序，
    # 这样可以保持 node_group_file 中分组列的自然顺序
    group_factor <- factor(group_feat, levels = unique(group_feat))
    ord <- order(group_factor, -deg, ids)
    value_mat <- value_mat[ord, ord, drop = FALSE]
    ids <- ids[ord]
    feat <- feat[ord]
    group_feat <- group_feat[ord]
    if (verbose) cat("[identity] Reordered matrix by feature and strength.\n")
  }
  
  # Build plotting data
  df <- as.data.frame(value_mat)
  df$y <- factor(ids, levels = ids)
  df <- tidyr::pivot_longer(df, cols = -y, names_to = "x", values_to = "value")
  df$x <- factor(df$x, levels = ids)
  
  # Use off-diagonal range for limits so diagonal does not compress the palette
  vmin <- vmin_off
  vmax <- vmax_off
  
  # Ensure palette length and order
  if (length(palette) < 5) {
    palette <- rep(palette, length.out = 5)
  }
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_gradientn(colors = palette, limits = c(vmin, vmax), name = if (color_by == "edge_count") "Edges" else if (color_by == "edge_weight") "Weight" else "Score") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = if (show_axis_labels) ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6) else ggplot2::element_blank(),
      axis.text.y = if (show_axis_labels) ggplot2::element_text(size = 6) else ggplot2::element_blank(),
      axis.title.x = if (show_axis_labels) ggplot2::element_text() else ggplot2::element_blank(),
      axis.title.y = if (show_axis_labels) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = if (show_axis_labels) "ID" else NULL, y = if (show_axis_labels) "ID" else NULL,
                  title = sprintf("Identity/Collinearity Matrix (%s | %s)", feature_by, color_by))
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(filename = output_file, plot = gg, width = width, height = height, dpi = dpi)
    if (verbose) cat(sprintf("[identity] Saved heatmap to %s\n", output_file))
  }
  
  list(gg = gg, matrix = value_mat, ids = ids)
}
