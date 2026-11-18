extract_ddrtree_ordering_compat <- function(cds, root_cell, verbose = TRUE) {
  dp <- monocle::cellPairwiseDistances(cds); dp_mst <- monocle::minSpanningTree(cds)
  states <- rep(1, ncol(dp)); names(states) <- igraph::V(dp_mst)$name
  pseudotimes <- rep(0, ncol(dp)); names(pseudotimes) <- igraph::V(dp_mst)$name
  parents <- rep(NA_character_, ncol(dp)); names(parents) <- igraph::V(dp_mst)$name
  traversal <- igraph::dfs(dp_mst, root = root_cell, mode = "all", unreachable = FALSE, parent = TRUE)
  parent_index <- as.numeric(traversal$parent); current_state <- 1L
  for (i in seq_along(traversal$order)) {
    current_index <- as.numeric(traversal$order[[i]]); current_name <- igraph::V(dp_mst)[current_index]$name
    if (!is.na(parent_index[[current_index]])) { parent_name <- igraph::V(dp_mst)[parent_index[[current_index]]]$name; current_pseudotime <- pseudotimes[[parent_name]] + dp[current_name, parent_name]; if (igraph::degree(dp_mst, v = parent_name) > 2L) current_state <- current_state + 1L } else { parent_name <- NA_character_; current_pseudotime <- 0 }
    states[[current_name]] <- current_state; pseudotimes[[current_name]] <- current_pseudotime; parents[[current_name]] <- parent_name
  }
  result <- data.frame(sample_name = names(states), cell_state = factor(states), pseudo_time = as.vector(pseudotimes), parent = parents, stringsAsFactors = FALSE); rownames(result) <- result$sample_name; result
}
dense_prim_mst_graph <- function(dp) {
  if (!is.matrix(dp) || nrow(dp) != ncol(dp) || nrow(dp) < 2L || any(!is.finite(dp))) stop("dense_prim_distance_matrix_invalid")
  ids <- rownames(dp); if (is.null(ids) || !identical(ids, colnames(dp))) stop("dense_prim_distance_names_invalid")
  n <- nrow(dp); selected <- rep(FALSE, n); key <- rep(Inf, n); parent <- integer(n); key[[1L]] <- 0
  for (step in seq_len(n)) {
    available <- which(!selected); u <- available[[which.min(key[available])]]
    if (!is.finite(key[[u]])) stop("dense_prim_disconnected")
    selected[[u]] <- TRUE; remaining <- which(!selected)
    if (length(remaining)) { candidate <- dp[u, remaining]; update <- candidate < key[remaining]; if (any(update)) { target <- remaining[update]; key[target] <- candidate[update]; parent[target] <- u } }
  }
  child <- seq.int(2L, n); if (any(parent[child] == 0L)) stop("dense_prim_parent_missing")
  edge_index <- cbind(parent[child], child)
  graph <- igraph::graph_from_edgelist(cbind(ids[edge_index[, 1L]], ids[edge_index[, 2L]]), directed = FALSE)
  igraph::E(graph)$weight <- dp[edge_index]
  graph
}
project2MST_compat <- function(cds, Projection_Method) {
  dp_mst <- monocle::minSpanningTree(cds); Z <- monocle::reducedDimS(cds); Y <- monocle::reducedDimK(cds); cds <- monocle:::findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex; closest_names <- colnames(Y)[closest_vertex]; closest_matrix <- as.matrix(closest_vertex); rownames(closest_matrix) <- rownames(closest_vertex)
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1L))
  if (!is.function(Projection_Method)) { P <- Y[, closest_vertex] } else { P <- matrix(0, nrow = nrow(Z), ncol = ncol(Z)); for (i in seq_along(closest_vertex)) { neighbors_i <- igraph::V(dp_mst)[igraph::neighbors(dp_mst, v = closest_names[[i]], mode = "all")]$name; projections <- NULL; distances <- NULL; Z_i <- Z[, i]; for (neighbor in neighbors_i) { candidate <- if (closest_names[[i]] %in% tip_leaves) monocle:::projPointOnLine(Z_i, Y[, c(closest_names[[i]], neighbor)]) else Projection_Method(Z_i, Y[, c(closest_names[[i]], neighbor)]); projections <- rbind(projections, candidate); distances <- c(distances, stats::dist(rbind(Z_i, candidate))) }; if (!is.matrix(projections)) projections <- as.matrix(projections); P[, i] <- projections[which(distances == min(distances))[[1L]], , drop = FALSE] } }
  colnames(P) <- colnames(Z); dp <- as.matrix(stats::dist(t(P))); minimum_nonzero <- min(dp[dp != 0]); dp <- dp + minimum_nonzero; diag(dp) <- 0
  cds <- monocle::`cellPairwiseDistances<-`(cds, value = dp); graph <- dense_prim_mst_graph(dp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- graph; cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P; cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_matrix; cds
}
orderCells_compat <- function(cds, ...) { compat_env <- new.env(parent = environment(monocle::orderCells)); assign("extract_ddrtree_ordering", extract_ddrtree_ordering_compat, envir = compat_env); assign("project2MST", project2MST_compat, envir = compat_env); fn <- monocle::orderCells; environment(fn) <- compat_env; fn(cds, ...) }


# Only the projected-row -> graph-name conversion is changed. No namespace write.
map_projected_rows <- function(cv, rows, principal_ids, root_node, vertex_names) {
  stopifnot(is.matrix(cv), ncol(cv)==1L, !is.null(rownames(cv)),
            !anyDuplicated(rownames(cv)), !anyDuplicated(principal_ids),
            !anyDuplicated(vertex_names), setequal(rownames(cv),vertex_names),
            length(root_node)==1L, root_node %in% principal_ids,
            all(is.finite(cv)), all(cv==as.integer(cv)),
            all(cv>=1L & cv<=length(principal_ids)), length(rows)>0L,
            !anyDuplicated(rows), all(rows>=1L & rows<=nrow(cv)))
  nodes <- principal_ids[as.integer(cv[,1])]
  stopifnot(identical(as.integer(rows),which(nodes==root_node)))
  barcodes <- rownames(cv)[rows]
  stopifnot(!anyNA(match(barcodes,vertex_names)))
  # Returning names preserves Monocle's original tip filtering and root choice.
  barcodes
}


fixed_order_cells <- function(cds, root_state=3) {
  # Resolve original historical helpers locally, preserving installed namespace.
  ce <- new.env(parent=environment(monocle::orderCells))
  ce$extract_ddrtree_ordering <- extract_ddrtree_ordering_compat
  ce$project2MST <- project2MST_compat
  ce$projection_rows_to_names <- function(cds,rows,root_node,K_old) {
    map_projected_rows(cds@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex,rows,colnames(K_old),root_node,igraph::V(monocle::minSpanningTree(cds))$name)
  }
  fn <- monocle::orderCells; environment(fn) <- ce
  old <- quote(V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name)
  replacement <- quote(projection_rows_to_names(cds,cells_mapped_to_graph_root,root_cell,K_old))
  changes <- 0L
  replace_one <- function(x) {
    if (identical(x,old)) { changes <<- changes+1L; return(replacement) }
    if (is.call(x)) for (i in seq_along(x)) if (i>1L && !identical(x[[i]],quote(expr=))) x[[i]] <- replace_one(x[[i]])
    x
  }
  body(fn) <- replace_one(body(fn))
  stopifnot(changes == 1L)
  fn(cds,root_state=root_state)
}
