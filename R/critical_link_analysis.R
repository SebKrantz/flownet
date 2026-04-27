#' @title Critical Link Analysis
#' @description Compute edge-level detour costs and vulnerability ratios for a graph.
#'
#' @param graph_df A data.frame with columns \code{from} and \code{to}, and optionally
#'   a cost column. Each row represents one physical edge.
#' @param directed Logical (default: \code{FALSE}). Whether to treat the graph as directed.
#'   For undirected graphs, each edge can be traversed in both directions, but the returned
#'   result remains aligned with the original rows of \code{graph_df}.
#' @param cost.column Character string (default: \code{"cost"}) or numeric vector. Name of
#'   the cost column in \code{graph_df}, or a numeric vector of edge costs with length equal
#'   to \code{nrow(graph_df)}.
#'
#' @return A data.frame with one row per input edge and columns:
#' \itemize{
#'   \item \code{edge} - Edge identifier (index of the edge in \code{graph_df})
#'   \item \code{from} - Original from-node ID
#'   \item \code{to} - Original to-node ID
#'   \item \code{edge_cost} - Cost of traversing the edge
#'   \item \code{detour_cost} - Shortest path cost between the edge endpoints after removing the edge
#'   \item \code{detour_ratio} - \code{detour_cost / edge_cost}
#' }
#'
#' @details
#' The detour cost for edge \code{e = (u, v)} is the shortest path cost from \code{u}
#' to \code{v} after removing that physical edge row. Parallel edges are kept as distinct
#' alternatives. If no detour exists, \code{detour_cost} and \code{detour_ratio} are \code{Inf}.
#'
#' The implementation uses a compiled multi-label Dijkstra routine per source node. For each
#' target, it tracks the two best distances whose first physical edge differs, which gives the
#' best endpoint detour for each outgoing edge without repeated graph deletion.
#'
#' @seealso \link{distances_from_graph} \link{run_assignment} \link{flownet-package}
#'
#' @examples
#' library(flownet)
#'
#' graph <- data.frame(
#'   from = c(1, 2, 1),
#'   to = c(2, 3, 3),
#'   cost = c(1, 1, 3)
#' )
#'
#' critical_link_analysis(graph, cost.column = "cost")
#'
#' @export
#' @importFrom collapse fmatch fnrow funique.default
critical_link_analysis <- function(graph_df, directed = FALSE, cost.column = "cost") {
  if(!is.data.frame(graph_df)) stop("graph_df must be a data.frame")
  if(!all(c("from", "to") %in% names(graph_df))) {
    stop("graph_df must have columns 'from' and 'to'. Missing: ",
         paste(setdiff(c("from", "to"), names(graph_df)), collapse = ", "))
  }

  n_edges <- fnrow(graph_df)
  cost <- if(is.character(cost.column) && length(cost.column) == 1L) {
    graph_df[[cost.column]]
  } else if(is.numeric(cost.column) && length(cost.column) == n_edges) {
    cost.column
  } else {
    stop("cost.column needs to be a column name in graph_df or a numeric vector matching nrow(graph_df)")
  }

  if(length(cost) != n_edges) {
    stop("cost.column needs to be provided either externally or found in the dataset")
  }
  cost <- as.numeric(cost)
  if(!all(is.finite(cost))) stop("cost.column must contain only finite values")
  if(any(cost <= 0)) stop("cost.column must contain strictly positive values")

  directed <- if(is.logical(directed) && length(directed) == 1L && !is.na(directed)) {
    directed
  } else {
    stop("directed must be TRUE or FALSE")
  }

  nodes <- funique.default(c(graph_df$from, graph_df$to), sort = TRUE)
  from_node <- fmatch(graph_df$from, nodes)
  to_node <- fmatch(graph_df$to, nodes)

  detour_cost <- .Call(
    C_critical_link_detours, # nolint: object_usage_linter
    from_node,
    to_node,
    cost,
    length(nodes),
    directed
  )

  data.frame(
    edge = if(is.numeric(graph_df[["edge"]])) graph_df[["edge"]] else seq_len(n_edges),
    from = graph_df$from,
    to = graph_df$to,
    edge_cost = cost,
    detour_cost = detour_cost,
    detour_ratio = detour_cost / cost
  )
}
