
#' @title Run Traffic Assignment
#' @description Assign traffic flows to network edges using either Path-Sized Logit (PSL)
#'   or All-or-Nothing (AoN) assignment methods.
#'
#' @param graph_df A data.frame with columns \code{from}, \code{to}, and optionally a cost column.
#'   Represents the network graph with edges between nodes.
#' @param od_matrix_long A data.frame with columns \code{from}, \code{to}, and \code{flow}.
#'   Represents the origin-destination matrix in long format with flow values.
#' @param directed Logical (default: FALSE). Whether the graph is directed.
#' @param cost.column Character string (default: "cost") or numeric vector. Name of the cost column
#'   in \code{graph_df}, or a numeric vector of edge costs with length equal to \code{nrow(graph_df)}.
#'   The cost values are used to compute shortest paths and determine route probabilities.
#' @param method Character string (default: "PSL"). Assignment method:
#'   \itemize{
#'     \item \code{"PSL"}: Path-Sized Logit model considering multiple routes with overlap correction
#'     \item \code{"AoN"}: All-or-Nothing assignment, assigns all flow to the shortest path (faster but no route alternatives)
#'   }
#' @param beta Numeric (default: 1). Path-sized logit parameter (beta_PSL). Only used for PSL method.
#' @param \dots Additional arguments (currently ignored).
#' @param detour.max Numeric (default: 1.5). Maximum detour factor for alternative routes (applied to shortest paths cost). Only used for PSL method. This is a key parameter controlling the execution time of the algorithm: considering more routes (higher \code{detour.max}) substantially increases computation time.
#' @param angle.max Numeric (default: 90). Maximum detour angle (in degrees, two sided). Only used for PSL method. I.e., nodes not within this angle measured against a straight line from origin to destination node will not be considered for detours.
#' @param unique.cost Logical (default: TRUE). Deduplicates paths based on the total cost prior to generating them. Only used for PSL method. Since multiple 'intermediate nodes' may be on the same path, there is likely a significant number of duplicate paths which this option removes.
#' @param npaths.max Integer (default: Inf). Maximum number of paths to compute per OD-pair. Only used for PSL method. If the number of paths exceeds this number, a random sample will be taken from all but the shortest path.
#' @param return.extra Character vector specifying additional results to return. Options include:
#'   \code{"graph"}, \code{"dmat"}, \code{"paths"}, \code{"edges"} (PSL only),
#'   \code{"counts"}, \code{"costs"}, and \code{"weights"} (PSL only).
#'   For AoN: \code{"paths"} returns a list of shortest paths (one integer vector per OD pair),
#'   \code{"costs"} returns a numeric vector of shortest path costs, and \code{"counts"} returns a global integer vector of edge traversal counts.
#'   Use \code{"all"} to return all available extra results for the selected method.
#' @param precompute.dmat Logical (default: TRUE). Should distance matrices be precomputed or computed on the fly. Only used for PSL method.
#'   The former is more memory intensive but dramatically speeds up the OD-iterations.
#' @param verbose Logical (default: TRUE). Show progress bar and intermediate steps completion status?
#' @param nthreads Integer (default: 1L). Number of threads (daemons) to use for parallel processing with \code{\link[mirai]{mirai}}. Should not exceed the number of logical processors.
#'
#'
#' @return A list of class \code{"flowr"} containing:
#'   \itemize{
#'     \item \code{call} - The function call
#'     \item \code{final_flows} - Numeric vector of assigned flows for each edge (same length as \code{nrow(graph_df)})
#'     \item \code{od_pairs_used} - Indices of OD pairs with valid flows
#'     \item Additional elements as specified in \code{return.extra}:
#'       \itemize{
#'         \item \code{graph} - The igraph graph object
#'         \item \code{dmat} - Distance matrix with node IDs as dimnames
#'         \item \code{paths} - For PSL: list of lists (multiple routes per OD pair); for AoN: list of integer vectors (one shortest path per OD pair)
#'         \item \code{edges} - List of edge indices used for each OD pair (PSL only)
#'         \item \code{edge_counts} - For PSL: list of edge visit counts per OD pair; for AoN: integer vector of global edge traversal counts
#'         \item \code{path_costs} - For PSL: list of path costs per OD pair; for AoN: numeric vector of shortest path costs
#'         \item \code{path_weights} - List of path weights (probabilities) for each OD pair (PSL only)
#'       }
#'   }
#'
#' @details
#' This function performs traffic assignment using one of two methods:
#'
#' \strong{All-or-Nothing (AoN) Method:}
#' A simple assignment method that assigns all flow from each OD pair to the single shortest path.
#' This is much faster than PSL but does not consider route alternatives or overlaps.
#' Parameters \code{detour.max}, \code{angle.max}, \code{unique.cost}, \code{npaths.max},
#' \code{beta}, and \code{precompute.dmat} are ignored for AoN.
#'
#' \strong{Path-Sized Logit (PSL) Method:}
#' A more sophisticated assignment method that considers multiple alternative routes:
#' \itemize{
#'   \item Creates a graph from \code{graph_df} using igraph, normalizing node IDs internally
#'   \item Computes shortest path distance matrix for all node pairs (if \code{precompute.dmat = TRUE})
#'   \item For each origin-destination pair in \code{od_matrix_long}:
#'     \itemize{
#'       \item Identifies alternative routes (detours) that are within \code{detour.max} of shortest path cost. This is done by considering all other nodes and adding the shortest paths costs from origin node to an intermediate node (any other node) and from intermediate node to destination node. If \code{angle.max} is specified, filters detours to those within the specified angle from origin to destination. This means we only consider intermediate nodes that are roughly 'in the direction' of the destination node, and also not further away in terms of geodesic distance. If also \code{unique.cost = TRUE}, duplicate paths are removed based on the total cost. Thus, using only the shortest-paths-cost matrix and a matrix of the geodesic distances of the nodes (if \code{is.finite(angle.max)}), the algorithm pre-selects unique paths that are plausible both in terms of detour factor (cost) and direction before actually computing them. This speeds up route enumeration and PSL computations considerably.
#'       \item Finds shortest paths from origin to intermediate nodes and from intermediate nodes to destination
#'       \item Filters paths to remove those with duplicate edges, i.e., where the intermediate node is approached and departed from via the same edge(s).
#'       \item Computes path-sized logit probabilities accounting for route overlap
#'       \item Assigns flows to edges based on probabilities weighted by route overlap
#'     }
#'   \item Returns results including final flows and optionally additional information
#' }
#'
#' The path-sized logit model accounts for route overlap by adjusting probabilities based on
#' the number of alternative routes using each edge. Flows are assigned proportionally to
#' the computed probabilities. The model uses the parameter \code{beta} to control the
#' sensitivity to route overlap.
#'
#' When \code{angle.max} is specified and \code{graph_df} contains coordinate columns
#' (\code{FX}, \code{FY}, \code{TX}, \code{TY}), the function uses geographic distance
#' calculations to restrict detours to those within the specified angle, improving
#' computational efficiency and route realism.
#'
#' @seealso \link{flowr-package}
#'
#' @examples
#' library(flowr)
#' library(sf)
#'
#' # Load Graph
#' graph <- linestrings_to_graph(network_gcc)
#' nodes <- nodes_from_graph(graph, sf = TRUE)
#'
#' # Map Zones to Nodes
#' nearest_nodes <- nodes$node[st_nearest_feature(zones_gcc, nodes)]
#'
#' # Process OD Matrix
#' od_matrix_long <- melt_od_matrix(od_matrices_gcc$container, nodes = nearest_nodes)
#'
#' # Run Traffic Assignment
#' result <- run_assignment(graph, od_matrix_long, cost.column = "generalized_cost",
#'                          return.extra = "all")
#' print(result)
#'
#' \dontrun{
#' # Visualize Results
#' network_gcc$final_flows_log10 <- log10(result$final_flows + 1)
#'
#' mapview(network_gcc, zcol = "final_flows_log10",
#'         layer.name = "Container Flows (Log10)") +
#' mapview(zones_gcc, alpha = 0, cex = 3,
#'         col.regions = "red", alpha.regions = 0.8,
#'         layer.name = "Zone Nodes")
#' }
#'
#' @export
#' @importFrom collapse funique.default ss fnrow seq_row ckmatch anyv whichv setDimnames fmatch %+=% gsplit setv any_duplicated fduplicated GRP vlengths
#' @importFrom igraph V graph_from_data_frame delete_vertex_attr igraph_options distances shortest_paths vcount ecount
#' @importFrom geodist geodist_vec
#' @importFrom mirai mirai_map daemons everywhere
#' @importFrom progress progress_bar
run_assignment <- function(graph_df, od_matrix_long,
                           directed = FALSE,
                           cost.column = "cost", # mode_col = NULL,
                           method = c("PSL", "AoN"),
                           beta = 1,
                           ...,
                           detour.max = 1.5,
                           angle.max = 90,
                           unique.cost = TRUE,
                           npaths.max = Inf,
                           return.extra = NULL,
                           precompute.dmat = TRUE,
                           verbose = TRUE,
                           nthreads = 1L) {

  cost <- if(is.character(cost.column) && length(cost.column) == 1L) graph_df[[cost.column]] else
    if(is.numeric(cost.column) && length(cost.column) == fnrow(graph_df)) cost.column else
    stop("cost.column needs to be a column name in graph_df or a numeric vector matching nrow(graph_df)")

  # Validate method
  method <- match.arg(method)
  is_aon <- method == "AoN"

  # Results object
  res <- list(call = match.call())
  if(length(return.extra) == 1L && return.extra == "all") {
    return.extra <- if(is_aon) c("graph", "dmat", "paths", "costs", "counts")
                    else c("graph", "dmat", "paths", "edges", "counts", "costs", "weights")
  }

  # Create Igraph Graph
  from_node <- as.integer(graph_df$from)
  to_node <- as.integer(graph_df$to)
  nodes <- funique.default(c(from_node, to_node), sort = TRUE)
  if(anyv(return.extra, "graph")) {
    res$graph <- graph_from_data_frame(data.frame(from = from_node, to = to_node),
                                       directed = directed,
                                       vertices = data.frame(name = nodes))
  }
  # Internally use normalized graph node indices
  from_node <- fmatch(from_node, nodes)
  to_node <- fmatch(to_node, nodes)
  g <- data.frame(from = from_node, to = to_node) |>
    graph_from_data_frame(directed = directed,
                          vertices = data.frame(name = seq_along(nodes))) |>
    delete_vertex_attr("name")

  if(verbose) cat("Created graph with", vcount(g), "nodes and", ecount(g), "edges...\n")

  # Geolocation and distance matrix are only used for PSL
  if(!is_aon) {
    geol <- is.finite(angle.max) && angle.max > 0 && angle.max < 180
    if(geol) {
      if(!all(c("FX", "FY", "TX", "TY") %in% names(graph_df))) {
        geol <- FALSE
        message("graph_df needs to have columns FX, FY, TX and TY to compute angle-based detour restrictions")
      } else {
        nodes_df <- nodes_from_graph(graph_df, sf = FALSE)
        X <- nodes_df$X
        Y <- nodes_df$Y
      }
    }
  }

  # Distance Matrix
  if((precompute.dmat && !is_aon) || anyv(return.extra, "dmat")) {
    dmat <- distances(g, mode = "out", weights = cost)
    if(nrow(dmat) != ncol(dmat)) stop("Distance matrix must be square")
    if(anyv(return.extra, "dmat")) res$dmat <- setDimnames(dmat, list(nodes, nodes))
    if(!is_aon) {
      dimnames(dmat) <- NULL
      if(geol) dmat_geo <- geodist_vec(X, Y, measure = "haversine")
    }
    if(verbose) cat("Computed distance matrix of dimensions", nrow(dmat), "x", ncol(dmat), "...\n")
  } else if(!is_aon) v <- V(g)

  # Final flows vector (just for placement)
  res$final_flows <- numeric(0)

  # Process/Check OD Matrix
  if(!all(c("from", "to", "flow") %in% names(od_matrix_long))) stop("od_matrix_long needs to have columns 'from', 'to' and 'flow'")
  od_pairs <- which(is.finite(od_matrix_long$flow) & od_matrix_long$flow > 0)
  res$od_pairs_used <- numeric(0) # Just for placement
  if(length(od_pairs) != fnrow(od_matrix_long)) od_matrix_long <- ss(od_matrix_long, od_pairs, check = FALSE)
  from <- ckmatch(od_matrix_long$from, nodes, e = "Unknown origin nodes in od_matrix:")
  to <- ckmatch(od_matrix_long$to, nodes, e = "Unknown destination nodes in od_matrix:")
  flow <- od_matrix_long[["flow"]]
  N <- length(flow)

  # Return block for AoN
  retvals_aon <- FALSE
  if(is_aon) {
    retvals_aon <- any(return.extra %in% c("paths", "costs", "counts"))
    if(retvals_aon) {
      if(anyv(return.extra, "paths")) {
        pathsl <- TRUE
        paths <- vector("list", N)
      } else pathsl <- FALSE
      if(anyv(return.extra, "costs")) {
        costsl <- TRUE
        sp_costs <- numeric(N)
      } else costsl <- FALSE
      if(anyv(return.extra, "counts")) {
        countsl <- TRUE
        edge_counts <- integer(length(cost))
      } else countsl <- FALSE
    }
  }

  # Return block (PSL only)
  retvals <- FALSE
  if(!is_aon) {
    retvals <- any(return.extra %in% c("paths", "edges", "counts", "costs", "weights"))
    if(retvals) {
      if(anyv(return.extra, "paths")) {
        pathsl <- TRUE
        paths <- vector("list", length(flow))
      } else pathsl <- FALSE
      if(anyv(return.extra, "edges")) {
        edgesl <- TRUE
        edges <- vector("list", length(flow))
      } else edgesl <- FALSE
      if(anyv(return.extra, "counts")) {
        countsl <- TRUE
        counts <- vector("list", length(flow))
      } else countsl <- FALSE
      if(anyv(return.extra, "costs")) {
        costsl <- TRUE
        costs <- vector("list", length(flow))
      } else costsl <- FALSE
      if(anyv(return.extra, "weights")) {
        weightsl <- TRUE
        weights <- vector("list", length(flow))
      } else weightsl <- FALSE
    }
  }

  # AoN Core Function - Batched by origin node for efficiency
  run_assignment_core_aon <- function(indices, verbose = FALSE, session = FALSE) {

    if(!session) {
      if(verbose) progress_bar <- progress::progress_bar
      shortest_paths <- igraph::shortest_paths
      igraph_options <- igraph::igraph_options
      GRP <- collapse::GRP
      gsplit <- collapse::gsplit
      flowr <- getNamespace("flowr")
      C_assign_flows_to_paths <- flowr$C_assign_flows_to_paths
      if(retvals_aon) {
        sve <- flowr$sve
        if(countsl) C_mark_edges_traversed <- flowr$C_mark_edges_traversed
        if(costsl) C_sum_path_costs <- flowr$C_sum_path_costs
      }
    }

    iopt <- igraph_options(return.vs.es = FALSE)
    on.exit(igraph_options(iopt))

    final_flows <- numeric(length(cost))

    # Group OD pairs by origin node for batched shortest path computation
    grp <- GRP(from[indices], return.groups = TRUE, call = FALSE, sort = FALSE)
    fromi <- grp$groups[[1L]]  # Unique origin nodes
    toi_idx <- gsplit(indices, grp)  # Indices grouped by origin

    if(verbose) {
      pb <- progress_bar$new(
        format = "Processed :current/:total OD-pairs (:percent) at :tick_rate/sec [Elapsed::elapsed | ETA::eta]",
        total = N, clear = FALSE
      )
      divp <- max(as.integer(nthreads), 1L)
      if(nthreads > 1L && length(indices) * divp > N)
         divp <- as.integer(N / length(indices))
    }

    for (i in seq_along(fromi)) {
      idx <- toi_idx[[i]]  # OD pair indices for this origin

      # Compute all shortest paths from this origin to all its destinations
      sp <- shortest_paths(g, from = fromi[i], to = to[idx], weights = cost,
                           mode = "out", output = "epath", algorithm = "automatic")$epath

      # Check for empty paths and mark as NA
      if(anyv(splen <- vlengths(sp), 0L)) od_pairs[idx[whichv(splen, 0L)]] <- NA_integer_

      # Assign flows to paths (batch operation)
      .Call(C_assign_flows_to_paths, sp, flow[idx], final_flows)

      # Handle return.extra for AoN
      if(retvals_aon) {
        if(pathsl) for(k in seq_along(idx)) sve(paths, idx[k], as.integer(sp[[k]]))
        if(costsl) .Call(C_sum_path_costs, sp, cost, sp_costs, idx)
        if(countsl) .Call(C_mark_edges_traversed, sp, edge_counts)
      }

      if(verbose) pb$tick(divp * length(idx))
    }

    if(verbose && !pb$finished) pb$tick(N - divp*length(indices))

    if(!session) {
      res <- list(final_flows = final_flows, od_pairs = od_pairs)
      if(retvals_aon) {
        if(pathsl) res$paths <- paths
        if(costsl) res$sp_costs <- sp_costs
        if(countsl) res$edge_counts <- edge_counts
      }
      return(res)
    }

    return(final_flows)
  }

  # PSL Core Function for parallel setup with mirai
  run_assignment_core_psl <- function(indices, verbose = FALSE, session = FALSE) {

    if(!session) {
      # Load required functions
      if(verbose) progress_bar <- progress::progress_bar
      distances <- igraph::distances
      shortest_paths <- igraph::shortest_paths
      igraph_options <- igraph::igraph_options
      geodist_vec <- geodist::geodist_vec
      whichv <- collapse::whichv
      setv <- collapse::setv
      `%+=%` <- collapse::`%+=%`
      any_duplicated <- collapse::any_duplicated
      fduplicated <- collapse::fduplicated
      flowr <- getNamespace("flowr")
      sve <- flowr$sve
      C_check_path_duplicates <- flowr$C_check_path_duplicates
      C_compute_path_sized_logit <- flowr$C_compute_path_sized_logit
      C_free_delta_ks <- flowr$C_free_delta_ks
    }

    # Don't return vertex/edge names
    iopt <- igraph_options(return.vs.es = FALSE) # sparsematrices = TRUE
    on.exit(igraph_options(iopt))

    # Final flows vector
    final_flows <- numeric(length(cost))

    # Edge incidence across selected routes
    delta_ks <- integer(length(cost) + 10L)

    if(verbose) {
      pb <- progress_bar$new(
        format = "Processed :current/:total OD-pairs (:percent) at :tick_rate/sec [Elapsed::elapsed | ETA::eta]", # [:bar] :percent eta: :eta",
        total = N, clear = FALSE #, # width = 60
      )
      div <- if(N > 1e4) 100L else 10L
      divp <- div * max(as.integer(nthreads), 1L)
      if(nthreads > 1L && as.integer(length(indices) / div) * divp > N)
         divp <- as.integer(N / as.integer(length(indices) / div))
    }

    # Now iterating across OD-pairs
    j <- 0L
    for (i in indices) {

      fi = from[i]
      ti = to[i]

      if(verbose) {
        j = j + 1L
        if(j %% div == 0L) pb$tick(divp)
      }

      if(precompute.dmat) {
        d_ij = dmat[fi, ti] # Shortest path cost
        d_ikj = dmat[fi, ] + dmat[, ti] # from i to all other nodes k and from these nodes k to j (basically dmat + t(dmat)?)
        if(geol) {
          b = dmat_geo[fi, ]
          a = b[ti]
          theta = acos((a^2 + b^2 - dmat_geo[, ti]^2)/(2*a*b)) * 180 / pi # Angle between a and b
        }
      } else {
        d_ikj = drop(distances(g, fi, v, mode = "out", weights = cost))
        d_ij = d_ikj[ti]
        d_ikj %+=% distances(g, v, ti, mode = "out", weights = cost)
        if(geol) {
          b = drop(geodist_vec(X[fi], Y[fi], X, Y, measure = "haversine"))
          a = b[ti]
          c = drop(geodist_vec(X, Y, X[ti], Y[ti], measure = "haversine"))
          theta = acos((a^2 + b^2 - c^2)/(2*a*b)) * 180 / pi # Angle between a and b
        }
      }
      short_detour_ij = if(geol) d_ikj < detour.max * d_ij & b < a & theta < angle.max else
                                 d_ikj < detour.max * d_ij
      short_detour_ij[d_ikj < d_ij + .Machine$double.eps*1e3] <- FALSE # Exclude nodes k that are on the shortest path
      # which(d_ij == d_ikj) # These are the nodes on the direct path from i to j which yield the shortest distance.
      ks = which(short_detour_ij)
      cost_ks = d_ikj[ks]
      # Remove paths that are likely equivalent
      if(unique.cost && any_duplicated(cost_ks)) {
        ndup = whichv(fduplicated(cost_ks), FALSE)
        cost_ks = cost_ks[ndup]
        ks = ks[ndup]
      }
      # If still too many paths: sample
      if(length(ks) > npaths.max) {
        ks = ks[sample.int(length(ks), npaths.max, useHash = FALSE)]
        cost_ks = d_ikj[ks]
      }

      # We add the shortest path at the end of paths1
      # TODO: Could still optimize calls to shortest_paths(), e.g., go to C directly.
      paths1 = shortest_paths(g, from = fi, to = c(ks, ti), weights = cost,
                              mode = "out", output = "epath", algorithm = "automatic")$epath
      paths2 = shortest_paths(g, from = ti, to = ks, weights = cost,
                              mode = "in", output = "epath", algorithm = "automatic")$epath
      shortest_path = paths1[[length(paths1)]]

      # # Check
      # cost_ks[k] == sum(cost[paths1[[k]]]) + sum(cost[paths2[[k]]])

      # Get indices of paths that do not contain duplicate edges
      no_dups = .Call(C_check_path_duplicates, paths1, paths2, delta_ks)

      # Now Path-Sized Logit: Need to compute overlap between routes
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
      # prob_ks <- proportions(exp(-c(cost_ks[no_dups], d_ij) + beta_PSL * log(c(gamma_ks, gamma_1))))
      #
      # # Need to reset delta_ks
      # delta_ks[] <- 0L
      #
      # # Assign result to edges:
      # for (k in no_dups) {
      #   final_flows[paths1[[k]]] <- final_flows[paths1[[k]]] + flow[i] * prob_ks[k]
      # }
      # final_flows[shortest_path] <- final_flows[shortest_path] + flow[i] * prob_ks[length(prob_ks)]
      wi = .Call(C_compute_path_sized_logit, paths1, paths2, no_dups, shortest_path,
                 cost, cost_ks, d_ij, beta, flow[i], delta_ks, final_flows, !retvals)
      if(is.null(wi)) {
        sve(od_pairs, i, NA_integer_)
        next
      }

      if(retvals) {
        if(pathsl) sve(paths, i, c(list(as.integer(shortest_path)), lapply(no_dups,
                      function(k) c(as.integer(paths1[[k]]), rev.default(as.integer(paths2[[k]]))))))
        if(countsl) {
          ei = whichv(delta_ks, 0L, invert = TRUE)
          if(edgesl) sve(edges, i, ei)
          sve(counts, i, delta_ks[ei])
        } else if(edgesl) sve(edges, i, whichv(delta_ks, 0L, invert = TRUE))
        if(costsl) sve(costs, i, c(d_ij, cost_ks[no_dups]))
        if(weightsl) sve(weights, i, wi)
        .Call(C_free_delta_ks, delta_ks, no_dups, paths1, paths2, shortest_path)
      }
    }

    if(verbose && !pb$finished) pb$tick(N - as.integer(j/div)*divp)
    # pb$terminate()

    if(!session) {
      res <- list(final_flows = final_flows, od_pairs = od_pairs)
      if(retvals) {
        if(pathsl) res$paths <- paths
        if(countsl) res$counts <- counts
        if(edgesl) res$edges <- edges
        if(costsl) res$costs <- costs
        if(weightsl) res$weights <- weights
      }
      return(res)
    }

    return(final_flows)
  }

  # Select core function based on method
  run_assignment_core <- if(is_aon) run_assignment_core_aon else run_assignment_core_psl

  if(!is.finite(nthreads) || nthreads <= 1L) {
    res$final_flows <- run_assignment_core(seq_len(N), verbose, TRUE)
  } else {
    envir <- environment()
    # Split OD matrix in equal parts
    ind_list <- gsplit(g = sample.int(as.integer(nthreads), N, replace = TRUE))
    daemons(n = nthreads - 1L)
    # Pass current environment dynamically
    everywhere({}, envir)
    # Now run the map in the background
    res_other <- mirai_map(ind_list[-1L], run_assignment_core)
    # Runs the first instance in the current session
    final_flows <- run_assignment_core(ind_list[[1L]], verbose, TRUE)
    # Collect other mirai's results
    res_other <- res_other[.stop] # [.stop, .progress] # collect_mirai()
    # Deactivate Daemons
    daemons(0)
    # Combine Results
    # return(environment())
    for (i in seq_along(res_other)) {
      resi <- res_other[[i]]
      ind <- ind_list[[i+1L]]
      final_flows %+=% resi$final_flows
      setv(od_pairs, ind, resi$od_pairs, vind1 = TRUE)
      if(retvals) {
        if(pathsl) paths[ind] <- resi$paths[ind]
        if(countsl) counts[ind] <- resi$counts[ind]
        if(edgesl) edges[ind] <- resi$edges[ind]
        if(costsl) costs[ind] <- resi$costs[ind]
        if(weightsl) weights[ind] <- resi$weights[ind]
      }
      if(retvals_aon) {
        if(pathsl) paths[ind] <- resi$paths[ind]
        if(costsl) sp_costs[ind] <- resi$sp_costs[ind]
        if(countsl) edge_counts %+=% resi$edge_counts
      }
    }
    res$final_flows <- final_flows
    rm(res_other, envir, ind_list, final_flows)
  }

  if(anyNA(od_pairs)) {
    nmiss_od <- whichNA(od_pairs, invert = TRUE)
    if(verbose) cat(length(od_pairs) - length(nmiss_od), "OD-pairs have zero or non-finite flow values and will be skipped...\n")
    res$od_pairs_used <- od_pairs[nmiss_od]
    if(retvals) {
      if(pathsl) res$paths <- paths[nmiss_od]
      if(edgesl) res$edges <- edges[nmiss_od]
      if(countsl) res$edge_counts <- counts[nmiss_od]
      if(costsl) res$path_costs <- costs[nmiss_od]
      if(weightsl) res$path_weights <- weights[nmiss_od]
    }
    if(retvals_aon) {
      if(pathsl) res$paths <- paths[nmiss_od]
      if(costsl) res$path_costs <- sp_costs[nmiss_od]
      if(countsl) res$edge_counts <- edge_counts
    }
  } else {
    res$od_pairs_used <- od_pairs
    if(retvals) {
      if(pathsl) res$paths <- paths
      if(edgesl) res$edges <- edges
      if(countsl) res$edge_counts <- counts
      if(costsl) res$path_costs <- costs
      if(weightsl) res$path_weights <- weights
    }
    if(retvals_aon) {
      if(pathsl) res$paths <- paths
      if(costsl) res$path_costs <- sp_costs
      if(countsl) res$edge_counts <- edge_counts
    }
  }

  class(res) <- "flowr" # , method
  return(res)
}

#' @rdname run_assignment
#'
#' @param x An object of class \code{flowr}, typically returned by \code{\link{run_assignment}}.
#'
#' @export
#' @importFrom collapse fmean fsd vlengths
print.flowr <- function(x, ...) {
  cat("Flowr object\n")
  cat("Call:", deparse(x$call), "\n\n")
  if (!is.null(x$dmat) && is.matrix(x$dmat))
    cat("Number of nodes:", nrow(x$dmat), "\n")
  cat("Number of edges:", length(x$final_flows), "\n")
  if (!is.null(x$od_pairs_used) && length(x$od_pairs_used))
    cat("Number of simulations/OD-pairs:", length(x$od_pairs_used), "\n")
  if (!is.null(x$paths) && length(x$paths)) {
    if (is.null(x$od_pairs_used) || !length(x$od_pairs_used))
      cat("Number of simulations/OD-pairs:", length(x$paths), "\n")
  }
  cat("\n")

  if (!is.null(x$paths) && length(x$paths)) {
    pls <- vlengths(x$paths)
    if(is.numeric(x$paths[[1L]])) {
      # For AoN, vlengths gives path length (number of edges per path)
      cat("Average path length in edges (SD): ", fmean(pls), "  (", fsd(pls, stable.algo = FALSE), ")\n", sep = "")
    } else {
      # For PSL, vlengths gives number of alternative paths per OD pair
      cat("Average number of paths per simulation (SD): ", fmean(pls), "  (", fsd(pls, stable.algo = FALSE), ")\n", sep = "")
    }
  }
  if (!is.null(x$edges) && length(x$edges)) {
    els <- vlengths(x$edges)
    cat("Average number of edges utilized per simulation (SD): ", fmean(els), "  (", fsd(els, stable.algo = FALSE), ")\n", sep = "")
  }
  if (!is.null(x$edge_counts) && length(x$edge_counts)) {
    if(is.numeric(x$edge_counts)) {
      # For AoN, edge_counts is a single integer vector (global counts)
      cat("Average number of visits per edge (SD): ", fmean(x$edge_counts), "  (", fsd(x$edge_counts, stable.algo = FALSE), ")\n", sep = "")
    } else {
      # For PSL, edge_counts is a list of counts per OD pair
      cat("Average number of visits per edge (SD): ", fmean(fmean(x$edge_counts)), "  (", fmean(fsd(x$edge_counts, stable.algo = FALSE)), ")\n", sep = "")
    }
  }
  if (!is.null(x$path_costs) && length(x$path_costs)) {
    if(is.numeric(x$path_costs)) {
      # For AoN, path_costs is a single numeric vector
      cat("Average path cost (SD): ", fmean(x$path_costs), "  (", fsd(x$path_costs, stable.algo = FALSE), ")\n", sep = "")
    } else {
      # For PSL, path_costs is a list of costs per OD pair
      cat("Average path cost (SD): ", fmean(fmean(x$path_costs)), "  (", fmean(fsd(x$path_costs, stable.algo = FALSE)), ")\n", sep = "")
    }
  }
  if (!is.null(x$path_weights) && length(x$path_weights))
    cat("Average path weight (SD): ", fmean(fmean(x$path_weights)), "  (", fmean(fsd(x$path_weights, stable.algo = FALSE)), ")\n", sep = "")
}
