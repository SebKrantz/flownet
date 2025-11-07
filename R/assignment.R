#' Traffic Assignment Functions
#'
#' Functions for assigning traffic to routes using nested logit models.


#' @title Run Assignment
#' @description Run traffic assignment for the network.
#'
#' @param network Character string specifying the network name.
#' @param scenario Character string specifying the scenario name.
#' @param year Character string specifying the analysis year (e.g., "2040").
#' @param enumeration_run_directory Optional character string path to enumeration results directory.
#'   If not provided, will use the most recent enumeration run.
#' @param cargo_type_map Optional named character vector mapping from OD matrix cargo types
#'   to enumeration cargo types. Defaults to common mappings (Container->Containerized, etc.).
#'
#' @details
#' This function:
#' \itemize{
#'   \item Loads network, zone nodes, OD matrix, and choice parameters
#'   \item Loads alternative routes from enumeration results
#'   \item Calculates tonnage assignment using nested logit model
#'   \item Aggregates results by link and direction
#'   \item Saves results as shapefile
#' }
#'
#' @return Invisibly returns NULL. Results are saved to the assignment results directory.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom sf st_write
run_assignment <- function(graph_df, od_matrix_long, directed = FALSE,
                           cost_col = "cost", mode_col = NULL,
                           method = "PSL", lambda = -1, beta = -1,
                           return.extra = c("graph_df", "od_matrix_long", "paths", "edges", "costs", "weights")) {

  cost <- if(is.character(cost_col) && length(cost_col) == 1L) graph_df[[cost_col]] else
    if(is.numeric(cost_col) && length(cost_col) == fnrow(graph_df)) cost_col else
      stop("cost_col needs to be a column name in graph_df or a numeric vector matching nrow(graph_df)")

  # Create Graph igraph: best so far
  g <- graph_df |> fselect(from, to) |> graph_from_data_frame(directed = FALSE)
  g <- delete_vertex_attr(g, "name")
  iopt <- igraph_options(return.vs.es = FALSE) # sparsematrices = TRUE
  on.exit(igraph_options(iopt))

  # Results object
  res <- list(call = match.call())
  if(length(return.extra) == 1L && return.extra == "all")
    return.extra <- c("graph_df", "od_matrix_long", "dmat", "paths", "edges", "costs", "weights")

  # Distance Matrix
  dmat <- distances(g, mode = "out", weights = cost)
  if(!identical(dim(dmat), rep(fnrow(nodes_df), 2))) stop("Nodes and Distance Matrix Mismatch")
  if(anyv(return.extra, "dmat")) res$dmat <- dmat |> setDimnames(alloc(nodes_df$node, 2))

  # Edge incidence across selected routes
  delta_ks <- integer(length(cost))
  # Final flows vector
  final_flows <- numeric(length(cost))
  if(anyv(return.extra, "graph_df")) {
    graph_df$final_fows <- final_fows
    res$graph_df <- graph_df
  } else res$final_flows <- final_flows

  if(!all(c("from", "to", "flow") %in% names(od_matrix_long))) stop("od_matrix_long needs to have columns 'from', 'to' and 'flow'")
  od_matrix_long %<>% fsubset(is.finite(flow) & flow > 0)
  if(anyv(return.extra, "od_matrix_long")) res$od_matrix_long <- od_matrix_long
  from <- od_matrix_long$from
  to <- od_matrix_long$to
  flow <- od_matrix_long$from
  ckmatch(from, nodes_df$node, e = "Unknown origin nodes:")
  ckmatch(to, nodes_df$node, e = "Unknown destination nodes:")

  retvals <- any(return.extra %in% c("paths", "edges", "costs", "weights"))
  if(retvals) {
    if(anyv(return.extra, "paths")) {
      pathsl <- TRUE
      paths <- vector("list", length(fnrow(od_matrix_long)))
    } else pathsl <- FALSE
  }
  # Now iterating across OD-pairs
  # TODO: could restrict that other nodes must be in the direction of travel and not behind destination node
  for (i in seq_row(od_matrix_long)) {

    d_ij <- dmat[from[i], to[i]] # Shortest path cost
    d_ikj <- dmat[from[i], ] + dmat[, to[i]] # from i to all other nodes k and from these nodes k to j (basically dmat + t(dmat)?)

    short_detour_ij <- d_ikj < 1.1 * d_ij
    short_detour_ij[d_ikj <= d_ij + .Machine$double.eps*1e3] <- FALSE # Exclude nodes k that are on the shortest path
    # which(d_ij == d_ikj) # These are the nodes on the direct path from i to j which yield the shortest distance.
    ks <- which(short_detour_ij)
    cost_ks <- d_ikj[ks]
    # Now Path-Sized Logit: Need to compute overlap between routes
    # Could still optimize calls to shortest_paths(), e.g., go to C directly.
    # We add the shortest path at the end of paths1
    paths1 <- shortest_paths(g, from = from[i], to = c(ks, to[i]), weights = cost, output = "epath")$epath
    paths2 <- shortest_paths(g, from = to[i], to = ks, weights = cost, output = "epath")$epath
    shortest_path <- paths1[[length(paths1)]]

    # # Check
    # cost_ks[k] == sum(cost[paths1[[k]]]) + sum(cost[paths2[[k]]])

    # Get indices of paths that do not contain duplicate edges
    no_dups <- check_path_duplicates(paths1, paths2, delta_ks)

    # # Number of routes in choice set that use link j
    # for (k in no_dups) {
    #   delta_ks[paths1[[k]]] <- delta_ks[paths1[[k]]] + 1L
    #   delta_ks[paths2[[k]]] <- delta_ks[paths2[[k]]] + 1L
    # }
    # delta_ks[shortest_path] <- delta_ks[shortest_path] + 1L
    #
    # # Correction factors for each route k
    # gamma_ks <- sapply(no_dups, function(k) {
    #   path <- c(paths1[[k]], paths2[[k]])
    #   sum(cost[path] / delta_ks[path]) / cost_ks[k]
    # })
    # gamma_1 <- sum(cost[shortest_path] / delta_ks[shortest_path]) / d_ij
    #
    # # Now the PS-MNL
    # prob_ks <- proportions(exp(c(cost_ks[no_dups], d_ij) + beta_PSL * c(gamma_ks, gamma_1)))
    #
    # # Need to reset delta_ks
    # delta_ks[] <- 0L
    #
    # # Assign result to edges:
    # for (k in no_dups) {
    #   final_flows[paths1[[k]]] <- final_flows[paths1[[k]]] + flow[i] * prob_ks[k]
    # }
    # final_flows[shortest_path] <- final_flows[shortest_path] + flow[i] * prob_ks[length(prob_ks)]
    compute_path_sized_logit(paths1, paths2, no_dups, shortest_path,
                             cost, cost_ks, d_ij, beta_PSL, flow[i],
                             delta_ks, final_flows)

    if(retvals) {
      if(pathsl) paths[[i]] <- c(list(shortest_path), lapply(no_dups, function(k) c(paths1[[k]], paths2[[k]])))
    }
  }

  if(retvals) {
    if(pathsl) res$paths <- paths

  }

}


