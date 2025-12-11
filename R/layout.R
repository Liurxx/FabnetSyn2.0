#' Compute subclade internal layouts and weights
#'
#' @param g igraph with vertex attribute `subclade`
#' @param seed integer seed
#' @return list with sub_layouts, sub_keep_indices, weight_df
#' @export
compute_subclade_layouts <- function(g, seed = 123) {
  set.seed(seed)
  subclades <- unique(stats::na.omit(igraph::V(g)$subclade))

  # compute weights
  subclade_weights <- list()
  for (sc in subclades) {
    sc_nodes <- which(igraph::V(g)$subclade == sc)
    if (length(sc_nodes) <= 1) {
      subclade_weights[[sc]] <- 1
      next
    }
    subg <- igraph::induced_subgraph(g, sc_nodes)
    internal_weight <- sum(igraph::E(subg)$weight, na.rm = TRUE)
    node_count <- length(sc_nodes)
    subclade_weights[[sc]] <- internal_weight * sqrt(node_count)
  }

  weight_df <- data.frame(
    subclade = names(subclade_weights),
    weight = unlist(subclade_weights),
    stringsAsFactors = FALSE
  )
  
  # 确保 weight 列存在并且是数值类型
  if (!"weight" %in% colnames(weight_df)) {
    stop("Failed to create weight column in weight_df")
  }
  
  # 确保 weight 是数值类型
  weight_df$weight <- as.numeric(weight_df$weight)
  weight_df$weight[is.na(weight_df$weight)] <- 1
  
  # 按权重降序排列
  weight_df <- weight_df[order(weight_df$weight, decreasing = TRUE), , drop = FALSE]
  
  # 重置行名
  rownames(weight_df) <- NULL

  total_weight <- sum(weight_df$weight)
  weight_df$cumulative_weight <- cumsum(weight_df$weight)
  weight_df$angle_start <- 2 * pi * c(0, head(weight_df$cumulative_weight, -1)) / total_weight
  weight_df$angle_end <- 2 * pi * weight_df$cumulative_weight / total_weight

  sub_layouts <- list()
  sub_keep_indices <- list()

  for (i in seq_len(nrow(weight_df))) {
    sc <- weight_df$subclade[i]
    sc_nodes <- which(igraph::V(g)$subclade == sc)
    if (length(sc_nodes) == 1) {
      sub_layouts[[sc]] <- matrix(c(0, 0), nrow = 1)
      sub_keep_indices[[sc]] <- sc_nodes
      next
    }
    subg <- igraph::induced_subgraph(g, sc_nodes)
    set.seed(abs(sum(utf8ToInt(as.character(sc)))) %% 10000)
    weight_factor <- weight_df$weight[i] / max(weight_df$weight)
    niter <- 800 + 400 * weight_factor
    loc <- igraph::layout_with_fr(subg, weights = igraph::E(subg)$weight, niter = niter)
    if (max(abs(loc)) > 0) loc <- loc / max(abs(loc))
    jitter_sd <- 0.03 * (1 - weight_factor * 0.5)
    jitter_mat <- matrix(stats::rnorm(length(sc_nodes) * 2, sd = jitter_sd), ncol = 2)
    loc <- loc + jitter_mat
    r <- sqrt(rowSums(loc^2))
    quantile_threshold <- 0.85 + 0.1 * weight_factor
    r_thr <- stats::quantile(r, probs = quantile_threshold, na.rm = TRUE)
    keep_mask <- r <= r_thr
    sub_layouts[[sc]] <- loc[keep_mask, , drop = FALSE]
    sub_keep_indices[[sc]] <- sc_nodes[keep_mask]
  }

  list(sub_layouts = sub_layouts, sub_keep_indices = sub_keep_indices, weight_df = weight_df)
}

#' Assemble circular packed layout from sub-layouts
#'
#' @param g igraph
#' @param sub_layouts list from compute_subclade_layouts
#' @param sub_keep_indices list from compute_subclade_layouts
#' @param weight_df data.frame from compute_subclade_layouts
#' @param outer_radius numeric
#' @param inner_radius numeric
#' @param angle_offset numeric
#' @return list with final_layout (matrix), keep_global (logical), sc_points (list)
#' @export
assemble_circular_layout <- function(g, sub_layouts, sub_keep_indices, weight_df,
                                     outer_radius = 1.35, inner_radius = 0.05,
                                     angle_offset = pi/90) {
  ring_thickness <- outer_radius - inner_radius
  target_edge_radius <- outer_radius * 0.80
  final_layout <- matrix(NA_real_, nrow = igraph::vcount(g), ncol = 2)
  rownames(final_layout) <- igraph::V(g)$name
  keep_global <- rep(FALSE, igraph::vcount(g))
  sc_points <- list()

  for (i in seq_len(nrow(weight_df))) {
    sc <- weight_df$subclade[i]
    if (!sc %in% names(sub_layouts)) next
    loc <- sub_layouts[[sc]]
    kept_nodes <- sub_keep_indices[[sc]]
    if (length(kept_nodes) == 0) next
    a1 <- weight_df$angle_start[i]
    a2 <- weight_df$angle_end[i]
    a_mid <- (a1 + a2) / 2
    weight_factor <- weight_df$weight[i] / max(weight_df$weight)
    rot_x <- loc[, 1] * cos(a_mid) - loc[, 2] * sin(a_mid)
    rot_y <- loc[, 1] * sin(a_mid) + loc[, 2] * cos(a_mid)
    tx_min <- min(rot_x); tx_max <- max(rot_x)
    ty_min <- min(rot_y); ty_max <- max(rot_y)
    if (tx_max - tx_min == 0) { tx_max <- tx_min + 1 }
    if (ty_max - ty_min == 0) { ty_max <- ty_min + 1 }
    t_norm <- (rot_x - tx_min) / (tx_max - tx_min)
    r_norm <- (rot_y - ty_min) / (ty_max - ty_min)
    angle_spread <- 0.995
    ang <- a1 + (a2 - a1) * (0.5 - angle_spread / 2 + angle_spread * t_norm) + angle_offset
    radial_spread <- 0.9 + 0.1 * weight_factor
    margin <- ring_thickness * 0.01
    radial <- inner_radius + margin + (ring_thickness - 2 * margin) *
      (0.02 + radial_spread * (pmax(pmin(r_norm, 1), 0))^0.98)
    ang_jitter <- stats::rnorm(length(ang), sd = (a2 - a1) * 0.015)
    radial_jitter <- stats::rnorm(length(radial), sd = ring_thickness * 0.012)
    ang <- ang + ang_jitter
    radial <- pmin(pmax(radial + radial_jitter, inner_radius + margin), outer_radius - margin)
    gx <- radial * cos(ang)
    gy <- radial * sin(ang)
    final_layout[kept_nodes, 1] <- gx
    final_layout[kept_nodes, 2] <- gy
    keep_global[kept_nodes] <- TRUE
    sc_points[[sc]] <- cbind(gx, gy)
  }

  list(final_layout = final_layout, keep_global = keep_global,
       sc_points = sc_points, target_edge_radius = target_edge_radius,
       outer_radius = outer_radius, inner_radius = inner_radius)
}

#' Adjust outer ring nodes to target radius
#'
#' @param layout_kept matrix of coordinates
#' @param target_edge_radius numeric
#' @param adjust_fraction fraction of outer nodes to consider (default 0.15)
#' @return list with layout_kept (matrix), outer_nodes (integer indices)
#' @export
adjust_outer_ring <- function(layout_kept, target_edge_radius, adjust_fraction = 0.15) {
  r_all <- sqrt(layout_kept[, 1]^2 + layout_kept[, 2]^2)
  outer_threshold <- target_edge_radius * (1 - adjust_fraction)
  outer_nodes <- which(r_all > outer_threshold)
  if (length(outer_nodes) > 0) {
    for (i in outer_nodes) {
      current_r <- r_all[i]
      if (current_r > target_edge_radius) {
        scale_factor <- target_edge_radius / current_r
        layout_kept[i, 1] <- layout_kept[i, 1] * scale_factor
        layout_kept[i, 2] <- layout_kept[i, 2] * scale_factor
      }
    }
  }
  list(layout_kept = layout_kept, outer_nodes = outer_nodes)
}

#' Local force optimization to reduce overlaps and tighten modules
#'
#' @param g_kept igraph subgraph of kept nodes
#' @param layout_kept coordinate matrix
#' @param target_edge_radius numeric
#' @param min_distance numeric
#' @param num_iter integer
#' @param step_size numeric
#' @param repel_strength numeric
#' @param attract_base numeric
#' @return optimized layout matrix
#' @export
force_optimize_layout <- function(g_kept, layout_kept, target_edge_radius,
                                  min_distance = 0.02, num_iter = 25,
                                  step_size = 0.015, repel_strength = 0.002,
                                  attract_base = 0.003) {
  edge_list <- igraph::as_edgelist(g_kept, names = TRUE)
  edge.data <- data.frame(from = edge_list[, 1], to = edge_list[, 2],
                          weight = igraph::E(g_kept)$weight, stringsAsFactors = FALSE)
  edge_w <- edge.data$weight
  edge_w <- ifelse(is.finite(edge_w), edge_w, 0)
  edge_w <- edge_w / max(edge_w, na.rm = TRUE)
  edge_w[is.na(edge_w)] <- 0
  name_to_idx <- stats::setNames(seq_len(nrow(layout_kept)), rownames(layout_kept))
  edge_i <- name_to_idx[edge.data$from]
  edge_j <- name_to_idx[edge.data$to]
  node_module_map <- stats::setNames(igraph::V(g_kept)$subclade, igraph::V(g_kept)$name)
  edge_same <- node_module_map[edge.data$from] == node_module_map[edge.data$to]
  edge_wt <- edge_w
  for (it in seq_len(num_iter)) {
    disp <- matrix(0, nrow = nrow(layout_kept), ncol = 2)
    for (i in 1:(nrow(layout_kept) - 1)) {
      xi <- layout_kept[i, 1]; yi <- layout_kept[i, 2]
      for (j in (i + 1):nrow(layout_kept)) {
        dx <- xi - layout_kept[j, 1]
        dy <- yi - layout_kept[j, 2]
        d2 <- dx * dx + dy * dy
        if (d2 <= 0) next
        d <- sqrt(d2)
        if (d < min_distance * 2.5) {
          f <- repel_strength / (d2 + 1e-6)
          fx <- f * dx
          fy <- f * dy
          disp[i, 1] <- disp[i, 1] + fx
          disp[i, 2] <- disp[i, 2] + fy
          disp[j, 1] <- disp[j, 1] - fx
          disp[j, 2] <- disp[j, 2] - fy
        }
      }
    }
    for (e in seq_along(edge_i)) {
      i <- edge_i[e]; j <- edge_j[e]
      if (is.na(i) || is.na(j)) next
      dx <- layout_kept[j, 1] - layout_kept[i, 1]
      dy <- layout_kept[j, 2] - layout_kept[i, 2]
      d <- sqrt(dx * dx + dy * dy) + 1e-6
      same_bonus <- ifelse(edge_same[e], 1.6, 0.9)
      f <- attract_base * edge_wt[e] * same_bonus
      fx <- f * dx
      fy <- f * dy
      disp[i, 1] <- disp[i, 1] + fx
      disp[i, 2] <- disp[i, 2] + fy
      disp[j, 1] <- disp[j, 1] - fx
      disp[j, 2] <- disp[j, 2] - fy
    }
    max_step <- step_size
    mag <- sqrt(disp[, 1]^2 + disp[, 2]^2) + 1e-9
    scale <- pmin(1, max_step / mag)
    disp[, 1] <- disp[, 1] * scale
    disp[, 2] <- disp[, 2] * scale
    layout_kept[, 1] <- layout_kept[, 1] + disp[, 1]
    layout_kept[, 2] <- layout_kept[, 2] + disp[, 2]
    rr <- sqrt(layout_kept[, 1]^2 + layout_kept[, 2]^2)
    overs <- which(rr > target_edge_radius)
    if (length(overs) > 0) {
      scale_b <- target_edge_radius / rr[overs]
      layout_kept[overs, 1] <- layout_kept[overs, 1] * scale_b
      layout_kept[overs, 2] <- layout_kept[overs, 2] * scale_b
    }
  }
  layout_kept
}


