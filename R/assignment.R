
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
#' @param dmat.max.size Integer (default: 1e4^2). Maximum size of distance matrices (both shortest paths and geodesic) to precompute. If smaller than \code{n_nodes^2}, then the full matrix is precomputed. Otherwise, it is computed in chunks as needed, where each chunk has \code{dmat.max.size} elements. Only used for PSL method.
#' @param return.extra Character vector specifying additional results to return.
#'   Use \code{"all"} to return all available extras for the selected method.
#'
#'   \tabular{llll}{
#'     \strong{Option} \tab \strong{PSL} \tab \strong{AoN} \tab \strong{Description} \cr
#'     \code{"graph"} \tab Yes \tab Yes \tab The igraph graph object \cr
#'     \code{"paths"} \tab Yes \tab Yes \tab PSL: list of lists of edge indices (multiple routes per OD); AoN: list of edge index vectors (one path per OD) \cr
#'     \code{"edges"} \tab Yes \tab No \tab List of edge indices used for each OD pair \cr
#'     \code{"counts"} \tab Yes \tab Yes \tab PSL: list of edge visit counts per OD; AoN: integer vector of global edge traversal counts \cr
#'     \code{"costs"} \tab Yes \tab Yes \tab PSL: list of path costs per OD; AoN: numeric vector of shortest path costs \cr
#'     \code{"weights"} \tab Yes \tab No \tab List of path weights (probabilities) for each OD pair \cr
#'   }
#' @param verbose Logical (default: TRUE). Show progress bar and intermediate steps completion status?
#' @param nthreads Integer (default: 1L). Number of threads (daemons) to use for parallel processing with \code{\link[mirai]{mirai}}. Should not exceed the number of logical processors.
#'
#'
#' @return A list of class \code{"flownet"} containing:
#'   \itemize{
#'     \item \code{call} - The function call
#'     \item \code{final_flows} - Numeric vector of assigned flows for each edge (same length as \code{nrow(graph_df)})
#'     \item \code{od_pairs_used} - Indices of OD pairs with valid flows
#'     \item Additional elements as specified in \code{return.extra}:
#'       \itemize{
#'         \item \code{graph} - The igraph graph object
#'         \item \code{paths} - For PSL: list of lists of edge indices (multiple routes per OD pair); for AoN: list of edge index vectors (one shortest path per OD pair)
#'         \item \code{edges} - List of edge indices used for each OD pair (PSL only)
#'         \item \code{edge_counts} - For PSL: list of edge visit counts per OD pair; for AoN: integer vector of global edge traversal counts
#'         \item \code{path_costs} - For PSL: list of path costs per OD pair; for AoN: numeric vector of shortest path costs
#'         \item \code{path_weights} - List of path weights (probabilities) for each OD pair (PSL only)
#'       }
#'   }
#'
#' @details
#' This function performs traffic assignment using one of two methods:
#' \strong{All-or-Nothing (AoN)} is fast but assigns all flow to a single shortest path;
#' \strong{Path-Sized Logit (PSL)} considers multiple routes with overlap correction for
#' more realistic flow distribution.
#'
#' \subsection{All-or-Nothing (AoN) Method}{
#' A simple assignment method that assigns all flow from each OD pair to the single shortest path.
#' This is much faster than PSL but does not consider route alternatives or overlaps.
#' Parameters \code{detour.max}, \code{angle.max}, \code{unique.cost}, \code{npaths.max},
#' \code{beta}, and \code{dmat.max.size} are ignored for AoN.
#' }
#'
#' \subsection{Path-Sized Logit (PSL) Method}{
#' A more sophisticated assignment method that considers multiple alternative routes and
#' accounts for route overlap when assigning flows. The PSL model adjusts choice probabilities
#' based on how much each route overlaps with other alternatives, preventing overestimation
#' of flow on shared segments. The \code{beta} parameter controls the sensitivity to overlap.
#' }
#'
#' \subsection{PSL Model Formulation}{
#' The probability \eqn{P_k} of choosing route \eqn{k} from the set of alternatives \eqn{K} is:
#' \deqn{P_k = \frac{e^{V_k}}{\sum_{j \in K} e^{V_j}}}
#' where the utility \eqn{V_k} is defined as:
#' \deqn{V_k = -C_k + \beta_{PSL} \ln(PS_k)}
#' Here \eqn{C_k} is the generalized cost of route \eqn{k}, \eqn{\beta_{PSL}} is the
#' path-size parameter (the \code{beta} argument), and \eqn{PS_k} is the path-size factor.
#'
#' The path-size factor quantifies route uniqueness:
#' \deqn{PS_k = \frac{1}{C_k} \sum_{a \in \Gamma_k} \frac{c_a}{\delta_a}}
#' where \eqn{\Gamma_k} is the set of edges in path \eqn{k}, \eqn{c_a} is the cost of
#' edge \eqn{a}, and \eqn{\delta_a} is the number of alternative routes using edge \eqn{a}.
#'
#' If a path is unique (\eqn{\delta_a = 1} for all edges), then \eqn{PS_k = 1} and the
#' model reduces to standard MNL. For overlapping routes, \eqn{PS_k < 1} and
#' \eqn{\ln(PS_k) < 0}, so a positive \code{beta} penalizes overlap. Higher \code{beta}
#' values strengthen penalization; \code{beta = 0} gives standard MNL behavior.
#'
#' For more information about the PSL model consult some of the references below.
#' }
#'
#' \subsection{Route Enumeration Algorithm}{
#' For each origin-destination pair, the algorithm identifies alternative routes as follows:
#' \enumerate{
#'   \item Compute the shortest path cost from origin to destination.
#'   \item For each potential intermediate node, calculate the total cost of going
#'         origin -> intermediate -> destination.
#'   \item Keep only routes where total cost is within \code{detour.max} times the
#'         shortest path cost.
#'   \item If \code{angle.max} is specified, filter to intermediate nodes that lie
#'         roughly in the direction of the destination (within the specified angle).
#'   \item If \code{unique.cost = TRUE}, remove duplicate routes based on total cost.
#'   \item Compute the actual paths and filter out those with duplicate edges
#'         (where the intermediate node is approached and departed via the same edge).
#' }
#' This pre-selection using distance matrices speeds up route enumeration considerably
#' by avoiding computation of implausible paths.
#' }
#'
#' \subsection{Coordinate-Based Filtering}{
#' When \code{angle.max} is specified and \code{graph_df} contains coordinate columns
#' (\code{FX}, \code{FY}, \code{TX}, \code{TY}), the function uses geographic distance
#' calculations to restrict detours. Only intermediate nodes that are (a) closer to the
#' origin than the destination is, and (b) within the specified angle from the
#' origin-destination line are considered. This improves both computational efficiency
#' and route realism by excluding geographically implausible detours.
#' }
#'
#' @seealso \link{flownet-package}
#'
#' @references
#' Ben-Akiva, M., & Bierlaire, M. (1999). Discrete choice methods and their
#' applications to short term travel decisions. In R. W. Hall (Ed.), *Handbook
#' of Transportation Science* (pp. 5–33). Springer US. https://doi.org/10.1007/978-1-4615-5203-1_2
#'
#' Cascetta, E. (2001). *Transportation systems engineering: Theory and methods*.
#' Springer.
#'
#' Ben-Akiva, M., & Lerman, S. R. (1985). *Discrete choice analysis: Theory and
#' application to travel demand*. The MIT Press.
#'
#' Ramming, M. S. (2002). *Network knowledge and route choice* (Doctoral
#' dissertation). Massachusetts Institute of Technology.
#'
#' Prato, C. G. (2009). Route choice modeling: Past, present and future research
#' directions. *Journal of Choice Modelling, 2*(1), 65–100. https://doi.org/10.1016/S1755-5345(13)70005-8
#'
#' AequilibiaE Python Documentation: https://www.aequilibrae.com/develop/python/route_choice/path_size_logit.html
#'
#' @examples
#' library(flownet)
#' library(collapse)
#' library(sf)
#'
#' # Load existing network edges (exclude proposed new links)
#' africa_net <- africa_network[!africa_network$add, ]
#'
#' # Convert to graph (use atomic_elem to drop sf geometry, qDF for data.frame)
#' graph <- atomic_elem(africa_net) |> qDF()
#' nodes <- nodes_from_graph(graph, sf = TRUE)
#'
#' # Map cities/ports to nearest network nodes
#' nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]
#'
#' # Simple gravity-based OD matrix
#' od_mat <- outer(africa_cities_ports$population, africa_cities_ports$population) / 1e12
#' dimnames(od_mat) <- list(nearest_nodes, nearest_nodes)
#' od_matrix_long <- melt_od_matrix(od_mat)
#'
#' # Run Traffic Assignment (All-or-Nothing method)
#' result_aon <- run_assignment(graph, od_matrix_long, cost.column = "duration",
#'                              method = "AoN", return.extra = "all")
#' print(result_aon)
#' \donttest{
#' # Run Traffic Assignment (Path-Sized Logit method)
#' # Note: PSL is slower but produces more realistic flow distribution
#' result_psl <- run_assignment(graph, od_matrix_long, cost.column = "duration",
#'                              method = "PSL", nthreads = 1L,
#'                              return.extra = c("edges", "counts", "costs", "weights"))
#' print(result_psl)
#'
#' # Visualize AoN Results
#' africa_net$final_flows_log10 <- log10(result_psl$final_flows + 1)
#' plot(africa_net["final_flows_log10"], main = "PSL Assignment")
#' }
#'
#'
#' # --- Trade Flow Disaggregation Example ---
#' # Disaggregate country-level trade to city-level using population shares
#'
#' # Compute each city's share of its country's population
#' city_pop <- africa_cities_ports |> atomic_elem() |> qDF() |>
#'   fcompute(node = nearest_nodes,
#'            city = qF(city_country),
#'            pop_share = fsum(population, iso3, TRA = "/"),
#'            keep = "iso3")
#'
#' # Aggregate trade to country-country level and disaggregate to cities
#' trade_agg <- africa_trade |> collap(quantity ~ iso3_o + iso3_d, fsum)
#' od_matrix_trade <- trade_agg |>
#'   join(city_pop |> add_stub("_o", FALSE), multiple = TRUE) |>
#'   join(city_pop |> add_stub("_d", FALSE), multiple = TRUE) |>
#'   fmutate(flow = quantity * pop_share_o * pop_share_d) |>
#'   frename(from = node_o, to = node_d) |>
#'   fsubset(flow > 0 & from != to)
#'
#' # Run AoN assignment with trade flows
#' result_trade_aon <- run_assignment(graph, od_matrix_trade, cost.column = "duration",
#'                                    method = "AoN", return.extra = "all")
#' print(result_trade_aon)
#' \donttest{
#' # Visualize trade flow results
#' africa_net$trade_flows_log10 <- log10(result_trade_aon$final_flows + 1)
#' plot(africa_net["trade_flows_log10"], main = "Trade Flow Assignment (AoN)")
#'
#' # Run PSL assignment with trade flows (nthreads can be increased for speed)
#' result_trade_psl <- run_assignment(graph, od_matrix_trade, cost.column = "duration",
#'                                    method = "PSL", nthreads = 1L,
#'                                    return.extra = c("edges", "counts", "costs", "weights"))
#' print(result_trade_psl)
#'
#' # Compare PSL vs AoN: PSL typically shows more distributed flows
#' africa_net$trade_flows_psl_log10 <- log10(result_trade_psl$final_flows + 1)
#' plot(africa_net["trade_flows_psl_log10"], main = "Trade Flow Assignment (PSL)")
#' }
#'
#' @export
#' @importFrom collapse funique.default ss fnrow seq_row ckmatch anyv whichv setDimnames fmatch %+=% gsplit setv any_duplicated fduplicated GRP
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
                           dmat.max.size = 10000^2,
                           return.extra = NULL,
                           verbose = TRUE,
                           nthreads = 1L) {

  cost <- if(is.character(cost.column) && length(cost.column) == 1L) as.numeric(graph_df[[cost.column]]) else
    if(is.numeric(cost.column) && length(cost.column) == fnrow(graph_df)) as.numeric(cost.column) else
    stop("cost.column needs to be a column name in graph_df or a numeric vector matching nrow(graph_df)")

  if(length(cost) != fnrow(graph_df)) stop("cost.column needs to be provided either externally or found in the dataset")

  # Validate method
  method <- match.arg(method)
  is_aon <- method == "AoN"

  # Results object
  res <- list(call = match.call())
  if(length(return.extra) == 1L && return.extra == "all") {
    return.extra <- if(is_aon) c("graph", "paths", "costs", "counts") # "dmat"
                    else c("graph", "paths", "edges", "counts", "costs", "weights") # "dmat"
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
  precompute.dmat <- dmat.max.size >= length(nodes)^2
  if(precompute.dmat && !is_aon) { #  || anyv(return.extra, "dmat")
    dmat <- distances(g, mode = "out", weights = cost)
    if(nrow(dmat) != ncol(dmat)) stop("Distance matrix must be square")
    if(anyv(return.extra, "dmat")) res$dmat <- setDimnames(dmat, list(nodes, nodes))
    # if(!is_aon) {
    dimnames(dmat) <- NULL
    if(geol) dmat_geo <- geodist_vec(X, Y, measure = "haversine")
    # }
    if(verbose) cat("Computed distance matrix of dimensions", nrow(dmat), "x", ncol(dmat), "...\n")
  } else if(!is_aon) {
    dmat_chunk_nrow <- as.integer(dmat.max.size / length(nodes))
    v <- V(g)
  }

  # Final flows vector (just for placement)
  res$final_flows <- numeric(0)

  # Process/Check OD Matrix
  if(!all(c("from", "to", "flow") %in% names(od_matrix_long))) stop("od_matrix_long must have columns 'from', 'to', 'flow'. Missing: ", paste(setdiff(c("from", "to", "flow"), names(od_matrix_long)), collapse = ", "))
  od_pairs <- which(is.finite(od_matrix_long$flow) & od_matrix_long$flow > 0)
  res$od_pairs_used <- numeric(0) # Just for placement
  if(length(od_pairs) != fnrow(od_matrix_long)) od_matrix_long <- ss(od_matrix_long, od_pairs, check = FALSE)
  from <- ckmatch(od_matrix_long$from, nodes, e = "Unknown origin nodes in od_matrix:")
  to <- ckmatch(od_matrix_long$to, nodes, e = "Unknown destination nodes in od_matrix:")
  flow <- as.numeric(od_matrix_long[["flow"]])
  N <- length(flow)

  # Return block
  retvals <- any(return.extra %in% c("paths", "edges", "counts", "costs", "weights"))
  if(retvals) {
    if(anyv(return.extra, "paths")) {
      pathsl <- TRUE
      paths <- vector("list", N)
    } else pathsl <- FALSE
    if(!is_aon && anyv(return.extra, "edges")) {
      edgesl <- TRUE
      edges <- vector("list", N)
    } else edgesl <- FALSE
    if(anyv(return.extra, "counts")) {
      countsl <- TRUE
      counts <- if(is_aon) integer(length(cost)) else vector("list", N)
    } else countsl <- FALSE
    if(anyv(return.extra, "costs")) {
      costsl <- TRUE
      costs <- if(is_aon) numeric(N) else vector("list", N)
    } else costsl <- FALSE
    if(!is_aon && anyv(return.extra, "weights")) {
      weightsl <- TRUE
      weights <- vector("list", N)
    } else weightsl <- FALSE
  }


  # AoN Core Function - Batched by origin node for efficiency
  run_assignment_core_aon <- function(indices, verbose = FALSE, session = FALSE) {

    if(!session) {
      if(verbose) progress_bar <- progress::progress_bar
      shortest_paths <- igraph::shortest_paths
      igraph_options <- igraph::igraph_options
      GRP <- collapse::GRP
      gsplit <- collapse::gsplit
      flownet <- getNamespace("flownet")
      C_assign_flows_to_paths <- flownet$C_assign_flows_to_paths
      if(retvals) {
        sve <- flownet$sve
        if(countsl) C_mark_edges_traversed <- flownet$C_mark_edges_traversed
        if(costsl) C_sum_path_costs <- flownet$C_sum_path_costs
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
      idx = toi_idx[[i]]  # OD pair indices for this origin

      # Compute all shortest paths from this origin to all its destinations
      sp = shortest_paths(g, from = fromi[i], to = to[idx], weights = cost,
                           mode = "out", output = "epath", algorithm = "automatic")$epath

      # Assign flows to paths (batch operation) + Check missing paths
      .Call(C_assign_flows_to_paths, sp, flow, final_flows, idx, od_pairs)

      # Handle return.extra for AoN
      if(retvals) {
        if(pathsl) for(k in seq_along(idx)) sve(paths, idx[k], as.integer(sp[[k]]))
        if(costsl) .Call(C_sum_path_costs, sp, cost, costs, idx)
        if(countsl) .Call(C_mark_edges_traversed, sp, counts)
      }

      if(verbose) pb$tick(divp * length(idx))
    }

    if(verbose && !pb$finished) pb$tick(N - divp*length(indices))

    if(!session) {
      res <- list(final_flows = final_flows, od_pairs = od_pairs)
      if(retvals) {
        if(pathsl) res$paths <- paths
        if(costsl) res$costs <- costs
        if(countsl) res$counts <- counts
      }
      return(res)
    }

    return(final_flows)
  }

  # PSL Core Function for parallel setup with mirai
  run_assignment_core_psl <- function(indices, verbose = FALSE, session = FALSE) {

    n <- length(indices)

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
      flownet <- getNamespace("flownet")
      sve <- flownet$sve
      C_check_path_duplicates <- flownet$C_check_path_duplicates
      C_compute_path_sized_logit <- flownet$C_compute_path_sized_logit
      C_free_delta_ks <- flownet$C_free_delta_ks
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
      if(nthreads > 1L && as.integer(n / div) * divp > N)
         divp <- as.integer(N / as.integer(n / div))
    }

    # Now iterating across OD-pairs
    j <- 0L
    for (i in indices) {

      fi = from[i]
      ti = to[i]

      if(precompute.dmat) {
        if(verbose) {
          j = j + 1L
          if(j %% div == 0L) pb$tick(divp)
        }
        d_ij = dmat[fi, ti] # Shortest path cost
        d_ikj = dmat[fi, ] + dmat[, ti] # from i to all other nodes k and from these nodes k to j (basically dmat + t(dmat)?)
        if(geol) {
          b = dmat_geo[fi, ]
          a = b[ti]
          theta = acos((a^2 + b^2 - dmat_geo[, ti]^2)/(2*a*b)) * 180 / pi # Angle between a and b
        }
      } else {
        if(j %% dmat_chunk_nrow == 0L) {
          k = 0L
          ind = indices[(j + 1L):min(j + dmat_chunk_nrow, n)]
          from_ind = from[ind]
          to_ind = to[ind]
          dmat_total = distances(g, from_ind, mode = "out", weights = cost)
          dimnames(dmat_total) <- NULL
          dmat_sp = dmat_total[cbind(seq_along(ind), to_ind)]
          dmat_total %+=% distances(g, to_ind, mode = "in", weights = cost)
          if(geol) {
            dmat_geo_rows = geodist_vec(X[from_ind], Y[from_ind], X, Y, measure = "haversine")
            dmat_geo_cols = geodist_vec(X, Y, X[to_ind], Y[to_ind], measure = "haversine")
          }
        }
        j = j + 1L
        k = k + 1L
        if(verbose && j %% div == 0L) pb$tick(divp)
        d_ij = dmat_sp[k] # Shortest path cost
        d_ikj = dmat_total[k, ] # from i to all other nodes k and from these nodes k to j
        if(geol) {
          b = dmat_geo_rows[k, ]
          a = b[ti]
          theta = acos((a^2 + b^2 - dmat_geo_cols[, k]^2)/(2*a*b)) * 180 / pi # Angle between a and b
        }
        # if(anyNA(theta)) stop(sprintf("k is %d, j is %d, chunk_size is %d, window size is %d", k, j, dmat_chunk_nrow, length(window)))
        # if(k > dmat_chunk_nrow) stop("k too large")
      }
      # # Skip self-loops (d_ij == 0) and invalid paths
      # if(d_ij <= .Machine$double.eps || !is.finite(d_ij)) {
      #   sve(od_pairs, i, NA_integer_)
      #   next
      # }
      short_detour_ij = if(geol) d_ikj < detour.max * d_ij & b < a & theta < angle.max else
                                 d_ikj < detour.max * d_ij
      short_detour_ij[d_ikj < d_ij + .Machine$double.eps*1e3] <- FALSE # Exclude nodes k that are on the shortest path
      # which(d_ij == d_ikj) # These are the nodes on the direct path from i to j which yield the shortest distance.
      ks = which(short_detour_ij)
      cost_ks = d_ikj[ks]
      # Remove paths that are likely equivalent (round to avoid floating point issues)
      if(unique.cost && any_duplicated(cost_ks_r <- floor(cost_ks * 1e8))) {
        ndup = whichv(fduplicated(cost_ks_r), FALSE)
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
    ind <- sample.int(as.integer(nthreads), N, replace = TRUE)
    ind_list <- gsplit(g = if(is_aon) sort(ind) else ind) # Since AoN should reduce calls to shortest_paths()
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
        if(edgesl) edges[ind] <- resi$edges[ind]
        if(costsl) costs[ind] <- resi$costs[ind]
        if(weightsl) weights[ind] <- resi$weights[ind]
        if(countsl) {
          if(is_aon) counts %+=% resi$counts
          else counts[ind] <- resi$counts[ind]
        }
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
      if(countsl) res$edge_counts <- if(is_aon) counts else counts[nmiss_od]
      if(costsl) res$path_costs <- costs[nmiss_od]
      if(weightsl) res$path_weights <- weights[nmiss_od]
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
  }

  class(res) <- "flownet" # , method
  return(res)
}

#' @rdname run_assignment
#'
#' @param x An object of class \code{flownet}, typically returned by \code{\link{run_assignment}}.
#'
#' @export
#' @importFrom collapse fmean fsd vlengths descr print.qsu
print.flownet <- function(x, ...) {
  cat("FlowNet object\n")
  cat("Call:", deparse(x$call), "\n\n")
  if (!is.null(x$dmat) && is.matrix(x$dmat))
    cat("Number of nodes:", nrow(x$dmat), "\n")
  else if(!is.null(x$graph) && inherits(x$graph, "igraph"))
    cat("Number of nodes:", vcount(x$graph), "\n")
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
  if(length(x$final_flows)) {
    if(length(x$call$return.extra)) cat("\n")
    dff <- descr(x$final_flows)
    cat("Final flows summary statistics:\n")
    print.qsu(dff$final_flows$Stats, digits = 2)
    print.qsu(dff$final_flows$Quant, digits = 2)
  }
}
