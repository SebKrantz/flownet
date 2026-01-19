
utils::globalVariables(c(
  "from", "to", "edge", "FX", "FY", "TX", "TY", "X", "Y", "cost", "flow", ".stop"
  # Add any other variable names that appear in the notes
  # "." # Often needed if you use the data.table or magrittr pipe syntax
))

sve <- function(x, i, elt) .Call(C_set_vector_elt, x, i, elt)

#' @title Convert Linestring to Graph
#'
#' @param lines An sf data frame of LINESTRING geometries.
#' @param digits Numeric rounding applied to coordinates (to ensure that matching points across different linestrings is not impaired by numeric precision issues). Set to \code{NA/Inf/FALSE} to disable.
#' @param keep.cols Character vector of column names to keep from the input data frame.
#' @param compute.length Applies \code{st_length()} to and saves it as an additional column named \code{".length"}.
#'
#' @return A data.frame representing the graph with columns:
#' \itemize{
#'  \item \code{edge} - Edge identifier
#'  \item \code{from} - Starting node ID
#'  \item \code{FX} - Starting node X-coordinate (longitude)
#'  \item \code{FY} - Starting node Y-coordinate (latitude)
#'  \item \code{to} - Ending node ID
#'  \item \code{TX} - Ending node X-coordinate (longitude)
#'  \item \code{TY} - Ending node Y-coordinate (latitude)
#' }
#'
#' @seealso \link{simplify_network} \link{flowr-package}
#'
#' @examples
#' library(flowr)
#'
#' # Convert network LINESTRING geometries to graph
#' graph <- linestrings_to_graph(network_gcc)
#' head(graph)
#'
#' # Graph contains edge, from/to nodes, and coordinates
#' names(graph)
#'
#' @export
#' @importFrom sf st_geometry_type st_coordinates st_length
#' @importFrom collapse qDF GRP get_vars get_vars<- add_vars add_vars<- fselect ffirst flast add_stub fmutate group fmatch %+=% fmax colorder whichNA setv unattrib ss
linestrings_to_graph <- function(lines, digits = 6, keep.cols = is.atomic, compute.length = TRUE) {
  gt <- st_geometry_type(lines, by_geometry = FALSE)
  if(length(gt) != 1L || gt != "LINESTRING") stop("lines needs to be a sf data frame of LINESTRING's")
  graph <- st_coordinates(lines) |> qDF()
  g <- GRP(list(edge = graph$L1), return.order = FALSE)
  graph <- add_vars(fselect(graph, X, Y) |> ffirst(g, na.rm = FALSE, use.g.names = FALSE) |> add_stub("F"),
                    fselect(graph, X, Y) |> flast(g, na.rm = FALSE, use.g.names = FALSE) |> add_stub("T")) |>
           add_vars(g$groups, pos = "front")
  if(compute.length) add_vars(graph) <- list(.length = st_length(lines))
  if(is.numeric(digits) && is.finite(digits)) {
    coords <- c("FX", "FY", "TX", "TY")
    get_vars(graph, coords) <- get_vars(graph, coords) |>
      lapply(round, digits = digits)
  }
  graph <- graph |>
    fmutate(from = unattrib(group(FX, FY)),
            to = from[fmatch(list(TX, TY), list(FX, FY))]) |>
    colorder(edge, from, FX, FY, to, TX, TY)
  if(anyNA(graph$to)) {
    miss <- whichNA(graph$to)
    setv(graph$to, miss, group(ss(graph, miss, c("TX", "TY"), check = FALSE)) %+=% fmax(graph$from), vind1 = TRUE)
  }
  if(!is.null(keep.cols)) add_vars(graph) <- get_vars(unclass(lines), keep.cols)
  graph
}


#' @title Convert Graph to Linestrings
#' @description Convert a graph data frame with node coordinates to an sf object with LINESTRING geometries.
#'
#' @param graph_df A data frame representing a graph with columns:
#'   \code{FX}, \code{FY}, \code{TX}, \code{TY} (starting and ending node coordinates),
#'   and optionally other columns to preserve.
#' @param crs Numeric or character (default: 4326). Coordinate reference system
#'   to assign to the output sf object.
#'
#' @return An sf data frame with LINESTRING geometry, containing all columns from
#'   \code{graph_df} except \code{FX}, \code{FY}, \code{TX}, and \code{TY}. Each row
#'   represents an edge as a LINESTRING connecting the from-node (\code{FX}, \code{FY})
#'   to the to-node (\code{TX}, \code{TY}).
#'
#' @details
#' This function is the inverse operation of \code{\link{linestrings_to_graph}}. It:
#' \itemize{
#'   \item Creates LINESTRING geometries from node coordinates (\code{FX}, \code{FY}, \code{TX}, \code{TY})
#'   \item Removes the coordinate columns from the output
#'   \item Preserves all other columns from the input graph data frame
#'   \item Returns an sf object suitable for spatial operations and visualization
#' }
#'
#' @seealso \link{linestrings_to_graph} \link{flowr-package}
#'
#' @examples
#' library(flowr)
#' library(sf)
#'
#' # Convert network to graph and back to linestrings
#' graph <- linestrings_to_graph(network_gcc)
#' lines <- linestrings_from_graph(graph)
#'
#' # Result is an sf object with LINESTRING geometry
#' class(lines)
#' head(lines)
#'
#' @export
#' @importFrom sf st_linestring st_sfc st_sf
#' @importFrom collapse seq_row fselect add_vars
linestrings_from_graph <- function(graph_df, crs = 4326) {
  if(!is.data.frame(graph_df)) stop("graph_df needs to be a data frame")
  if(inherits(graph_df, "sf")) stop("graph_df should not be a spatial object/data frame")
  if(!all(c("FX", "FY", "TX", "TY") %in% names(graph_df))) stop("graph_df needs to have columns FX, FY, TX and TY")
  # Create Geometries
  lines_list <- with(graph_df, lapply(seq_row(graph_df), function(i) {
    matrix(c(FX[i], FY[i], TX[i], TY[i]), ncol = 2, byrow = TRUE) |>
    st_linestring()
  })) |> st_sfc(crs = crs)
  # Create sf data frame with all columns
  graph_df |>
    fselect(-FX, -FY, -TX, -TY) |>
    add_vars(list(geometry = lines_list)) |>
    st_sf(sf_column_name = "geometry", crs = crs)
}

#' @title Create Undirected Graph
#' @description Convert a directed graph to an undirected graph by normalizing edges and aggregating duplicate connections.
#'
#' @param graph_df A data frame representing a directed graph including columns:
#'   \code{from}, \code{to}, and (optionally) \code{edge}, \code{FX}, \code{FY}, \code{TX}, \code{TY}.
#' @param by Link characteristics to preserve/not aggregate across, passed as a one-sided formula or character vector of column names. Typically this includes attributes like \emph{mode}, \emph{type}, or \emph{capacity} to ensure that only edges with the same characteristics are aggregated.
#' @param \dots Arguments passed to \code{\link[collapse]{collap}()} for aggregation across duplicated (diretional) edges. The defaults are \code{FUN = fmean} for numeric columns and \code{catFUN = fmode} for categorical columns. Select columns using \code{cols} or use argument \code{custom = list(fmean = cols1, fsum = cols2, fmode = cols3)} to map different columns to specific aggregation functions. You can weight the aggregation (using \code{w = ~ weight_col}).
#'
#' @return A data frame representing an undirected graph with:
#'   \itemize{
#'     \item \code{edge} - Edge identifier (first value from duplicates)
#'     \item \code{from} - Starting node ID (normalized to be < \code{to})
#'     \item \code{to} - Ending node ID (normalized to be > \code{from})
#'     \item \code{FX} - Starting node X-coordinate (first value from duplicates)
#'     \item \code{FY} - Starting node Y-coordinate (first value from duplicates)
#'     \item \code{TX} - Ending node X-coordinate (first value from duplicates)
#'     \item \code{TY} - Ending node Y-coordinate (first value from duplicates)
#'     \item Aggregated columns
#'   }
#'
#' @details
#' This function converts a directed graph to an undirected graph by:
#' \itemize{
#'   \item Normalizing edge directions so that \code{from < to} for all edges
#'   \item Collapsing duplicate edges (same \code{from} and \code{to} nodes)
#'   \item For spatial/identifier columns (\code{edge}, \code{FX}, \code{FY}, \code{TX}, \code{TY}),
#'     taking the first value from duplicates
#'   \item For aggregation columns, \code{\link[collapse]{collap}()} will be applied.
#' }
#'
#' @examples
#' library(flowr)
#'
#' # Create graph and convert to undirected
#' graph <- linestrings_to_graph(network_gcc)
#' graph_undir <- create_undirected_graph(graph)
#'
#' # Use 'by' to preserve mode-specific edges
#' graph_undir_mm <- create_undirected_graph(graph, by = ~ mode)
#'
#' # Fewer edges after removing directional duplicates
#' c(directed = nrow(graph), undirected = nrow(graph_undir))
#'
#' @export
#' @importFrom collapse ftransform GRP get_vars add_vars add_vars<- ffirst colorderv %!in% collap
create_undirected_graph <- function(graph_df, by = NULL, ...) {
  if(length(by)) {
    if(is.call(by)) by <- all.vars(by)
    if(!is.character(by)) stop("by needs to be a one-sided formula or a character vector of column names")
  }
  graph_df <- ftransform(graph_df, from = pmin(from, to), to = pmax(from, to))
  g <- GRP(graph_df, c("from", "to", by), sort = FALSE)
  nam <- names(graph_df)
  agg_first <- c("edge", "FX", "FY", "TX", "TY")
  ord <- c("edge", "from", "FX", "FY", "to", "TX", "TY", by)
  res <- g$groups
  if(any(agg_first %in% nam)) {
    add_vars(res) <- ffirst(get_vars(graph_df, nam[nam %in% agg_first]), g, use.g.names = FALSE)
    res <- colorderv(res, ord[ord %in% nam])
  }
  if(any(nam %!in% ord)) {
    add_vars(res) <- collap(get_vars(graph_df, nam[nam %!in% ord]), g, keep.by = FALSE, ...)
  }
  attr(res, "group.starts") <- g$group.starts
  res
}

#' @title Extract Nodes from Graph
#' @description Extract unique nodes with their coordinates from a graph data frame.
#'
#' @param graph_df A data frame representing a graph with columns:
#'   \code{from}, \code{to}, \code{FX}, \code{FY}, \code{TX}, \code{TY}.
#' @param sf Logical. If TRUE, returns result as an \code{sf} POINT object. Default: FALSE.
#' @param crs Coordinate reference system for sf output; default is 4326.
#'
#' @return A data frame (or sf object if \code{sf = TRUE}) with unique nodes and coordinates:
#'   \itemize{
#'     \item \code{node} - Node ID
#'     \item \code{X} - Node X-coordinate (typically longitude)
#'     \item \code{Y} - Node Y-coordinate (typically latitude)
#'   }
#'   Result is sorted by node ID.
#'
#' @details
#' This function extracts all unique nodes from both the \code{from} and \code{to}
#' columns of the graph, along with their corresponding coordinates. Duplicate nodes
#' are removed, keeping only unique node IDs with their coordinates.
#'
#' @examples
#' library(flowr)
#' library(sf)
#'
#' # Extract nodes from graph
#' graph <- linestrings_to_graph(network_gcc)
#' nodes <- nodes_from_graph(graph)
#' head(nodes)
#'
#' # Get nodes as sf POINT object for spatial operations
#' nodes_sf <- nodes_from_graph(graph, sf = TRUE)
#' class(nodes_sf)
#'
#' # Find nearest network nodes to zone centroids
#' nearest <- nodes_sf$node[st_nearest_feature(zones_gcc, nodes_sf)]
#'
#' @export
#' @importFrom collapse rowbind fselect funique
#' @importFrom stats setNames
#' @importFrom sf st_as_sf
nodes_from_graph <- function(graph_df, sf = FALSE, crs = 4326) {
  nodes <- rowbind(graph_df |> fselect(from, FX, FY),
                   graph_df |> fselect(to, TX, TY), use.names = FALSE) |>
    setNames(c("node", "X", "Y")) |>
    funique(cols = "node", sort = TRUE)
  if(sf) return(st_as_sf(nodes, coords = c("X", "Y"), crs = crs))
  nodes
}

#' @title Compute Distance Matrix from Graph
#' @description Compute a distance matrix for all node pairs in a graph using cppRouting.
#'
#' @param graph_df A data frame representing a graph with columns:
#'   \code{from}, \code{to}, and \code{cost}.
#' @param directed Logical (default: FALSE). If TRUE, treats the graph as directed;
#'   if FALSE, treats it as undirected.
#' @param cost.column Character string (optional). Name of the cost column in \code{graph_df}.
#'   Alternatively, a numeric vector of edge costs with length equal to \code{nrow(graph_df)}.
#' @param \dots Additional arguments passed to \code{\link[igraph]{distances}()}, such as \code{v} (from) and \code{to} to compute paths between specific nodes.
#'
#' @return A matrix of distances between all node pairs, where rows and columns
#'   correspond to node IDs. The matrix contains the shortest path distances
#'   (based on the \code{cost} column) between all pairs of nodes.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Converts the graph data frame to a cppRouting graph object
#'   \item Contracts the graph for efficient distance computation
#'   \item Computes the distance matrix for all node pairs using the specified algorithm
#' }
#'
#' @examples
#' library(flowr)
#'
#' # Create a simple graph
#' graph <- data.frame(
#'   from = c(1, 2, 2, 3),
#'   to = c(2, 3, 4, 4),
#'   cost = c(1, 2, 3, 1)
#' )
#'
#' # Compute distance matrix
#' dmat <- distances_from_graph(graph, cost.column = "cost")
#' dmat
#'
#' @export
#' @importFrom collapse fselect fnrow funique.default
#' @importFrom igraph graph_from_data_frame distances
distances_from_graph <- function(graph_df, directed = FALSE, cost.column = "cost", ...) {
  cost <- if(is.character(cost.column) && length(cost.column) == 1L) graph_df[[cost.column]] else
    if(is.numeric(cost.column) && length(cost.column) == fnrow(graph_df)) cost.column else
    stop("cost.column needs to be a column name in graph_df or a numeric vector matching nrow(graph_df)")

  # Create Igraph Graph
  vertices <- data.frame(name = funique.default(c(graph_df$from, graph_df$to), sort = TRUE))
  g <- graph_df |> fselect(from, to) |>
    graph_from_data_frame(directed = directed, vertices = vertices)

  distances(g, mode = "out", weights = cost, ...)
  # graph <- makegraph(graph_df |> fselect(from, to, cost), directed = directed) # directed = FALSE # cpp_simplify()
  # nodes <- graph$dict$ref
  # get_distance_matrix(cpp_contract(graph), from = nodes, to = nodes, algorithm = algorithm, ...)
}

#' @title Normalize Graph Node IDs
#' @description Normalize node IDs in a graph to be consecutive integers starting from 1.
#'   This is useful for ensuring compatibility with graph algorithms that require sequential node IDs.
#'
#' @param graph_df A data frame representing a graph with columns:
#'   \code{from} and \code{to} (node IDs).
#'
#' @return A data frame with the same structure as \code{graph_df}, but with \code{from}
#'   and \code{to} columns remapped to consecutive integer IDs starting from 1.
#'   All other columns are preserved unchanged.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts all unique node IDs from both \code{from} and \code{to} columns
#'   \item Sorts them in ascending order
#'   \item Remaps the original node IDs to sequential integers (1, 2, 3, ...)
#'   \item Updates both \code{from} and \code{to} columns with the normalized IDs
#' }
#'
#' Normalization is useful when:
#' \itemize{
#'   \item Node IDs are non-consecutive (e.g., 1, 5, 10, 20)
#'   \item Node IDs are non-numeric or contain gaps
#'   \item Graph algorithms require sequential integer node IDs starting from 1
#' }
#'
#' Note: This function only normalizes the node IDs; it does not modify the graph structure
#' or any other attributes. The mapping preserves the relative ordering of nodes.
#'
#' @seealso \link{nodes_from_graph} \link{flowr-package}
#'
#' @examples
#' library(flowr)
#'
#' # Create graph with non-consecutive node IDs
#' graph <- data.frame(
#'   from = c(10, 20, 20),
#'   to = c(20, 30, 40),
#'   cost = c(1, 2, 3)
#' )
#'
#' # Normalize to consecutive integers (1, 2, 3, 4)
#' graph_norm <- normalize_graph(graph)
#' graph_norm
#'
#' @export
#' @importFrom collapse funique.default get_vars get_vars<- fmatch
normalize_graph <- function(graph_df) {
  id_cols <- c("from", "to")
  if(!all(id_cols %in% names(graph_df))) stop("graph_df must have columns 'from' and 'to'")
  nodes <- funique.default(c(graph_df$from, graph_df$to), sort = TRUE)
  get_vars(graph_df, id_cols) <- lapply(get_vars(graph_df, id_cols), fmatch, nodes)
  graph_df
}


#' @title Consolidate Graph
#' @description Consolidate a graph by removing intermediate nodes (nodes that occur exactly twice) and optionally dropping loop, duplicate, and singleton edges (leading to dead ends). This simplifies the network topology while preserving connectivity.
#'
#' @param graph_df A data frame representing a graph with columns:
#'   \code{from} and \code{to} (node IDs), and optionally other columns to preserve.
#'   If coordinate columns (\code{FX}, \code{FY}, \code{TX}, \code{TY}) are present, they will be
#'   preserved and updated based on the consolidated node coordinates.
#' @param directed Logical (default: FALSE). Whether the graph is directed.
#' @param drop.edges Character vector (default: \code{c("loop", "duplicate", "single")}). Types of edges to drop:
#'   \code{"loop"} removes self-loops (edges where from == to),
#'   \code{"duplicate"} removes duplicate edges (same from-to pair),
#'   \code{"single"} removes edges leading to singleton nodes (nodes that occur only once).
#'   Set to \code{NULL} to keep all edges.
#' @param consolidate Logical (default: TRUE). If TRUE, consolidates the graph by removing
#'   intermediate nodes (nodes that occur exactly twice) and merging connecting edges.
#'   If FALSE, only drops edges as specified in \code{drop.edges}.
#' @param by Link characteristics to preserve/not consolidate across, passed as a one-sided formula or character vector of column names. Typically this includes attributes like \emph{mode}, \emph{type}, or \emph{capacity} to ensure that only edges with the same characteristics are consolidated.
#' @param keep.nodes Numeric vector (optional). Node IDs to preserve during consolidation,
#'   even if they occur exactly twice. Also used to preserve nodes when dropping singleton edges.
#' @param \dots Arguments passed to \code{\link[collapse]{collap}()} for aggregation across consolidated edges. The defaults are \code{FUN = fmean} for numeric columns and \code{catFUN = fmode} for categorical columns. Select columns using \code{cols} or use argument \code{custom = list(fmean = cols1, fsum = cols2, fmode = cols3)} to map different columns to specific aggregation functions. It is highly recommended to weight the aggregation (using \code{w = ~ weight_col}) by the length/cost of the edges.
#' @param recursive One of \code{"none"/FALSE}, \code{"partial"} (recurse on dropping single edges and consolidation but only aggregate once), or \code{"full"/TRUE} (recursively consolidates and aggregates the graph
#'   until no further consolidation is possible). This ensures that long chains of intermediate
#'   nodes are fully consolidated in a single call.
#' @param verbose Logical (default: TRUE). Whether to print messages about dropped edges
#'   and consolidation progress.
#'
#' @return A data frame representing the consolidated graph with:
#'   \itemize{
#'     \item \code{edge} - Edge identifier (added as first column)
#'     \item All columns from \code{graph_df} (aggregated if consolidation occurred),
#'       excluding \code{from}, \code{to}, and optionally \code{FX}, \code{FY}, \code{TX}, \code{TY}
#'       (which are re-added if present in original)
#'     \item \code{from}, \code{to} - Node IDs (updated after consolidation)
#'     \item Coordinate columns (\code{FX}, \code{FY}, \code{TX}, \code{TY}) if present in original
#'     \item Attribute \code{"keep.edges"} - Indices of original edges that were kept
#'     \item Attribute \code{"gid"} - Edge group IDs mapping consolidated edges to original edges
#'   }
#'
#' @details
#' This function simplifies a graph by:
#' \itemize{
#'   \item \strong{Dropping edges}: Optionally removes self-loops, duplicate edges, and edges
#'     leading to singleton nodes (nodes that appear only once in the graph)
#'   \item \strong{Consolidating nodes}: Removes intermediate nodes (nodes that occur exactly twice)
#'     by merging the two edges connected through them into a single longer edge
#'   \item \strong{Aggregating attributes}: When edges are merged, attributes/columns are aggregated
#'     using \code{\link[collapse]{collap}()}. The default aggregation is mean for numeric columns and mode for categorical columns.
#'   \item \strong{Recursive consolidation}: If \code{recursive = TRUE}, the function continues
#'     consolidating until no more nodes can be consolidated, ensuring complete simplification
#' }
#'
#' Consolidation is useful for simplifying network topology while preserving connectivity.
#' For example, if node B connects A->B and B->C, it will be removed and replaced with A->C.
#' With \code{recursive = TRUE}, long chains (A->B->C->D) are fully consolidated to A->D in
#' a single call.
#'
#' For undirected graphs, the algorithm also handles cases where a node appears twice
#' as either origin or destination (circular cases).
#'
#' If coordinate columns (\code{FX}, \code{FY}, \code{TX}, \code{TY}) are present in the input,
#' they are preserved and updated based on the consolidated node coordinates from the original graph.
#'
#' @seealso \link{create_undirected_graph} \link{simplify_network} \link{flowr-package}
#'
#' @examples
#' library(flowr)
#' library(sf)
#'
#' # Create undirected graph
#' graph <- linestrings_to_graph(network_gcc) |>
#'   create_undirected_graph(by = ~ mode)
#'
#' # Get nodes to preserve (zone locations)
#' nodes <- nodes_from_graph(graph, sf = TRUE)
#' nearest_nodes <- nodes$node[st_nearest_feature(zones_gcc, nodes)]
#'
#' # Consolidate graph, preserving zone nodes
#' graph_cons <- consolidate_graph(graph, by = ~ mode,
#'                                 w = ~ length_km,
#'                                 keep = nearest_nodes)
#'
#' # Consolidated graph has fewer edges
#' c(original = nrow(graph), consolidated = nrow(graph_cons))
#'
#' @export
#' @importFrom collapse fnrow get_vars anyv setv ss seq_row fduplicated fmatch whichv whichNA allNA ffirst GRP collap %iin% %!in% %!iin% join colorderv funique.default %!=% %==% missing_cases qtab flast varying radixorderv groupv na_rm
#' @importFrom kit countOccur
#' @importFrom stats setNames
consolidate_graph <- function(graph_df, directed = FALSE,
                              drop.edges = c("loop", "duplicate", "single"),
                              consolidate = TRUE, by = NULL, keep.nodes = NULL, ...,
                              recursive = "full",
                              verbose = TRUE) {

  if(verbose) namg <- flast(as.character(substitute(graph_df)))

  reci <- switch(as.character(recursive), none =, `FALSE` = 0L, partial = 1L, full =, `TRUE` = 2L,
                 stop("recursive needs to be one of 'none'/FALSE, 'partial', or 'full'/TRUE"))

  if(length(by)) {
    if(is.call(by)) by <- all.vars(by)
    if(!is.character(by)) stop("by needs to be a character vector or a formula of column names")
  }

  if(length(attr(graph_df, "group.starts"))) attr(graph_df, "group.starts") <- NULL

  nam <- names(graph_df)
  nam_rm <- c("from", "to", "FX", "FY", "TX", "TY", "edge", by)
  nam_keep <- nam[nam %!iin% nam_rm]

  if(verbose) {
    cat(sprintf("Consolidating %s graph %s with %d edges using %s recursion\n", if(directed) "directed" else "undirected", namg, fnrow(graph_df), as.character(recursive)))
    print(qtab(countOccur(c(graph_df$from, graph_df$to))$Count, dnn = "Initial node degrees:"))
    cat("\n")
  }

  res <- consolidate_graph_core(graph_df, directed = directed,
                                drop.edges = drop.edges,
                                consolidate = consolidate,
                                by = by,
                                keep.nodes = keep.nodes,
                                reci = reci, nam_keep = nam_keep,
                                verbose = verbose, ...)

  if(length(attr(res, ".early.return"))) {
    attr(res, ".early.return") <- NULL
    return(res)
  }

  if(reci == 2L && fnrow(res)) {
    nam_keep <- nam_keep[nam_keep %iin% names(res)]
    prev_fnrow <- fnrow(graph_df)
    while(prev_fnrow > (nrow_res <- fnrow(res))) {
      prev_fnrow <- nrow_res
      res <- consolidate_graph_core(res, directed = directed,
                                    drop.edges = drop.edges,
                                    consolidate = consolidate,
                                    by = by,
                                    keep.nodes = keep.nodes,
                                    reci = reci, nam_keep = nam_keep,
                                    verbose = verbose, ...)
    }
  }

  if(length(attr(res, ".early.return"))) attr(res, ".early.return") <- NULL

  if(verbose) {
    cat(sprintf("\nConsolidated %s graph %s from %d edges to %d edges (%s%%)\n", if(directed) "directed" else "undirected", namg, fnrow(graph_df), fnrow(res), as.character(signif(100*fnrow(res)/fnrow(graph_df), 3))))
    print(qtab(countOccur(c(res$from, res$to))$Count, dnn = "Final node degrees:"))
  }

  if(any(nam_rm[3:6] %in% nam)) {
    nodes <- nodes_from_graph(graph_df, sf = FALSE)
    if(any(nam_rm[3:4] %in% nam)) res <- join(res, setNames(nodes, c("from", "FX", "FY")), on = "from", verbose = 0L)
    if(any(nam_rm[5:6] %in% nam)) res <- join(res, setNames(nodes, c("to", "TX", "TY")), on = "to", verbose = 0L)
  }
  add_vars(res, pos = "front") <- list(edge = seq_row(res))

  # Reordering columns
  res <- colorderv(res, radixorderv(fmatch(names(res), nam)))
  res
}


# Corec function that can be called recursively
consolidate_graph_core <- function(graph_df, directed = FALSE,
                              drop.edges = c("loop", "duplicate", "single"),
                              consolidate = TRUE, by = NULL, keep.nodes = NULL, ...,
                              reci, nam_keep, verbose = TRUE) {

  keep <- seq_row(graph_df) # Global variable tracking utilized edges
  gft <- get_vars(graph_df, c("from", "to", by)) |> unclass() # Local variable representing the current graph worked on

  if(length(by)) {
    by_id <- groupv(get_vars(graph_df, by))
    # We keep nodes where there are changes (e.g., different mode).
    keep.nodes <- funique.default(c(keep.nodes,
     as.integer(names(which(varying(c(by_id, by_id), c(graph_df$from, graph_df$to), any = FALSE))))))
  }

  if(anyv(drop.edges, "loop") && length(loop <- gft$from %==% gft$to)) {
    keep <- keep[-loop]
    gft <- ss(gft, keep, check = FALSE)
    if(verbose) cat(sprintf("Dropped %d loop edges\n", length(loop)))
  }

  if(anyv(drop.edges, "duplicate") && any(dup <- fduplicated(gft))) {
    if(verbose) cat(sprintf("Dropped %d duplicate edges\n", sum(dup)))
    dup <- whichv(dup, FALSE)
    keep <- keep[dup]
    gft <- ss(gft, dup, check = FALSE)
  }

  if(anyv(drop.edges, "single") && fnrow(gft)) {
    repeat {
      nodes_rm <- unclass(countOccur(c(gft$from, gft$to)))
      if(!anyv(nodes_rm$Count, 1L)) break
      nodes_rm <- nodes_rm$Variable[nodes_rm$Count %==% 1L]
      if(length(keep.nodes)) nodes_rm <- nodes_rm[nodes_rm %!iin% keep.nodes]
      if(length(nodes_rm)) {
        ind <- which(gft$from %!in% nodes_rm & gft$to %!in% nodes_rm)
        if(verbose) cat(sprintf("Dropped %d edges leading to singleton nodes\n", fnrow(gft) - length(ind)))
        keep <- keep[ind]
        gft <- ss(gft, ind, check = FALSE)
      } else break
      if(reci == 0L) break
    }
  }

  if(!consolidate) {
    res <- ss(graph_df, keep, check = FALSE)
    if(reci < 2L) attr(res, "keep.edges") <- keep
    attr(res, ".early.return") <- TRUE
    return(res)
  }
  # TODO: How does not dropping loop or duplicate edges affect the algorithm?

  gid <- seq_row(gft)  # Local variable mapping current edges to groups
  consolidated_any <- FALSE

  merge_linear_nodes <- function(nodes) {
    if(!length(nodes)) return(FALSE)
    from_ind <- fmatch(nodes, gft$from)
    to_ind <- fmatch(nodes, gft$to)
    if(anyNA(from_ind) || anyNA(to_ind)) {
      valid <- whichv(missing_cases(list(from_ind, to_ind)), FALSE)
      if(!length(valid)) return(FALSE)
      from_ind <- from_ind[valid]
      to_ind <- to_ind[valid]
      nodes <- nodes[valid]
    }
    setv(gft$from, from_ind, NA, vind1 = TRUE) # gft$from[from_ind] <<- NA
    to_ind_prev <- integer(0)
    repeat {
      setv(gid, from_ind, gid[to_ind], vind1 = TRUE) # gid[from_ind] <<- gid[to_ind]
      setv(gft$to, to_ind, gft$to[from_ind], vind1 = TRUE) # gft$to[to_ind] <<- gft$to[from_ind]
      to_ind <- fmatch(nodes, gft$to)
      if(length(to_ind) == 0L || identical(to_ind, to_ind_prev) || allNA(to_ind)) break
      valid <- whichNA(to_ind, invert = TRUE)
      from_ind <- from_ind[valid]
      to_ind <- to_ind[valid]
      nodes <- nodes[valid]
      to_ind_prev <- to_ind
    }
    ffirst(gft$from, gid, "fill", set = TRUE)
    TRUE
  }

  orient_undirected_nodes <- function(nodes) {
    if(!length(nodes)) return(FALSE)
    N <- length(gft$from)
    Np <- N + 1L
    Ninv <- N:1 # Not strictly necessary to take second match, but appears faster and catch more nodes...
    idx_from <- Np - na_rm(fmatch(nodes, gft$from[Ninv]))
    if(length(idx_from)) {
      tmp_from <- gft$from[idx_from]
      tmp_from_to <- gft$to[idx_from]
    }
    idx_to <- Np - na_rm(fmatch(nodes, gft$to[Ninv]))
    if(length(idx_to)) {
      tmp_to <- gft$to[idx_to]
      tmp_to_from <- gft$from[idx_to]
    }
    if(length(idx_from)) {
      setv(gft$from, idx_from, tmp_from_to, vind1 = TRUE) # gft$from[idx_from] <<- tmp_from_to
      setv(gft$to, idx_from, tmp_from, vind1 = TRUE) # gft$to[idx_from] <<- tmp_from
    }
    if(length(idx_to)) {
      setv(gft$to, idx_to, tmp_to_from, vind1 = TRUE) # gft$to[idx_to] <<- tmp_to_from
      setv(gft$from, idx_to, tmp_to, vind1 = TRUE)  # gft$from[idx_to] <<- tmp_to
    }
    # # Old slow (iterative) version
    # for(node in nodes) {
    #   if(length(idx <- whichv(gft$from, node))) {
    #     idx <- idx[2L]
    #     tmp <- gft$from[idx]
    #     gft$from[idx] <<- gft$to[idx]
    #     gft$to[idx] <<- tmp
    #   } else if(length(idx <- whichv(gft$to, node))) {
    #     idx <- idx[2L]
    #     tmp <- gft$to[idx]
    #     gft$to[idx] <<- gft$from[idx]
    #     gft$from[idx] <<- tmp
    #   }
    # }
    TRUE
  }

  repeat {

    degree_table <- compute_degrees(gft$from, gft$to)
    if(!fnrow(degree_table)) break

    if(anyv(drop.edges, "single") && anyv(degree_table$deg_total, 1L)) {
      nodes <- degree_table$node[degree_table$deg_total %==% 1L]
      if(length(keep.nodes)) nodes <- nodes[nodes %!iin% keep.nodes]
      if(length(nodes)) {
        ind <- which(gft$from %!in% nodes & gft$to %!in% nodes)
        dropped <- fnrow(gft) - length(ind)
        if(dropped > 0L) {
          if(verbose) cat(sprintf("Dropped %d edges leading to singleton nodes\n", dropped))
          gft <- ss(gft, ind, check = FALSE)
          gid <- gid[ind]
          keep <- keep[ind]
          if(reci > 0L) next
        }
      }
    }

    if(!anyv(degree_table$deg_total, 2L)) break

    if(directed) {
      nodes <- degree_table$node[degree_table$deg_from == 1L & degree_table$deg_to == 1L]
      if(length(keep.nodes)) nodes <- nodes[nodes %!iin% keep.nodes]
      if(!length(nodes)) break
    } else {
      nodes <- degree_table$node[degree_table$deg_total %==% 2L]
      if(length(keep.nodes)) nodes <- nodes[nodes %!iin% keep.nodes]
      if(!length(nodes)) break
      idx <- fmatch(nodes, degree_table$node)
      need_orientation <- nodes[degree_table$deg_from[idx] == 2L | degree_table$deg_to[idx] == 2L]
      if(length(need_orientation)) {
        if(!orient_undirected_nodes(need_orientation)) stop("Failed to orient undirected edges for consolidation; please verify the input graph.")
        if(verbose) cat(sprintf("Oriented %d undirected intermediate edges\n", length(need_orientation)))
      }
    }
    if(!merge_linear_nodes(nodes)) stop("Failed to consolidate oriented undirected edges; please verify the graph topology.")
    consolidated_any <- TRUE
    if(verbose) cat(sprintf("Consolidated %d intermediate nodes\n", length(nodes)))
    if(reci == 0L) break
  }

  if(!consolidated_any) {
    if(verbose) cat("No nodes to consolidate, returning graph\n")
    res <- ss(graph_df, keep, check = FALSE)
    if(reci < 2L) attr(res, "keep.edges") <- keep
    attr(res, ".early.return") <- TRUE
    return(res)
  }

  # Grouping
  if(anyv(drop.edges, "duplicate")) {
    g <- GRP(gft, sort = TRUE)
  } else {
    g <- GRP(c(gft, list(gid = gid)), sort = TRUE)
    g$groups <- g$groups[seq_along(gft)]
    g$group.vars <- g$group.vars[seq_along(gft)]
  }

  # Aggregation
  res <- ss(graph_df, keep, nam_keep, check = FALSE)
  if(verbose) cat("Aggregated", length(keep), "edges down to", g$N.groups, "edges\n")
  res <- collap(res, g, keep.col.order = FALSE, ...)
  if(reci < 2L) {
    attr(res, "keep.edges") <- keep
    attr(res, "group.id") <- g$group.id
  }
  res
}


# Helper for consolidate_graph()
compute_degrees <- function(from_vec, to_vec) {
  nodes <- funique.default(c(from_vec, to_vec))
  deg_from <- integer(length(nodes))
  deg_to <- integer(length(nodes))
  if(length(nodes)) {
    if(length(from_vec)) {
      counts <- unclass(countOccur(from_vec))
      idx <- fmatch(nodes, counts$Variable, nomatch = 0L)
      has <- idx %!=% 0L
      if(length(has)) setv(deg_from, has, counts$Count[idx[has]], vind1 = TRUE)
    }
    if(length(to_vec)) {
      counts <- unclass(countOccur(to_vec))
      idx <- fmatch(nodes, counts$Variable, nomatch = 0L)
      has <- idx %!=% 0L
      if(length(has)) setv(deg_to, has, counts$Count[idx[has]], vind1 = TRUE)
    }
  }
  list(
    node = nodes,
    deg_from = deg_from,
    deg_to = deg_to,
    deg_total = deg_from + deg_to
  )
}

#' @title Simplify Network
#' @description Simplify a network graph using shortest paths or node clustering methods.
#'
#' @param graph A data.frame with columns \code{from} and \code{to} representing the graph edges.
#'   For the cluster method, the graph must also have columns \code{FX}, \code{FY}, \code{TX}, \code{TY}
#'   representing node coordinates.
#' @param nodes For \code{method = "shortest-paths"}: either an atomic vector of node IDs, or a
#'   data.frame with columns \code{from} and \code{to} specifying origin-destination pairs.
#'   For \code{method = "cluster"}: an atomic vector of node IDs to preserve. These nodes will
#'   be kept as cluster centroids, and nearby nodes (within \code{radius_km$nodes}) will be
#'   assigned to their clusters. Remaining nodes are clustered using \code{\link[leaderCluster]{leaderCluster}}.
#' @param method Character string (default: "shortest-paths"). Method to use for simplification:
#'   \code{"shortest-paths"} computes shortest paths between nodes and keeps only traversed edges;
#'   \code{"cluster"} clusters nodes using the \code{\link[leaderCluster]{leaderCluster}} algorithm and contracts the graph.
#' @param directed Logical (default: FALSE). Whether the graph is directed.
#'   For \code{method = "shortest-paths"}: controls path computation direction.
#'   For \code{method = "cluster"}: if TRUE, A->B and B->A remain as separate edges after
#'   contraction; if FALSE, edges are normalized so that \code{from < to} before grouping.
#' @param cost.column Character string (default: "cost"). Name of the cost column in \code{graph}.
#'   Alternatively, a numeric vector of edge costs with length equal to \code{nrow(graph)}.
#'   With \code{method = "cluster"}, a numeric vector of node weights matching \code{nodes_from_graph(graph)} can be provided.
#' @param by Link characteristics to preserve/not simplify across, passed as a one-sided
#'   formula or character vector of column names. Typically includes attributes like
#'   \emph{mode}, \emph{type}, or \emph{capacity}.
#'   For \code{method = "shortest-paths"}: paths are computed separately for each group
#'   defined by \code{by}, with edges not in the current group penalized (cost set to
#'   100x max cost) to compel mode-specific routes.
#'   For \code{method = "cluster"}: edges are grouped by \code{from}, \code{to}, AND
#'   \code{by} columns, preventing consolidation across different modes/types.
#' @param radius_km Named list with elements \code{nodes} (default: 7) and \code{cluster} (default: 20).
#'   Only used for \code{method = "cluster"}.
#'   \code{nodes}: radius in kilometers around preserved nodes. Graph nodes within this radius
#'   will be assigned to the nearest preserved node's cluster.
#'   \code{cluster}: radius in kilometers for clustering remaining nodes using leaderCluster.
#' @param \dots For \code{method = "cluster"}: additional arguments passed to
#'   \code{\link[collapse]{collap}} for edge attribute aggregation.
#'
#' @return A data.frame containing the simplified graph with:
#'   \itemize{
#'     \item For \code{method = "shortest-paths"}:
#'       \itemize{
#'         \item All columns from the input \code{graph} (for edges that were kept)
#'         \item Attribute \code{"edges"}: integer vector of edge indices from the original graph
#'         \item Attribute \code{"edge_counts"}: integer vector indicating how many times each edge was traversed
#'       }
#'     \item For \code{method = "cluster"}:
#'       \itemize{
#'         \item \code{edge} - New edge identifier
#'         \item \code{from}, \code{to} - Cluster centroid node IDs
#'         \item \code{FX}, \code{FY}, \code{TX}, \code{TY} - Coordinates of cluster centroid nodes
#'         \item Aggregated edge attributes from the original graph
#'         \item Attribute \code{"group.id"}: mapping from original edges to simplified edges
#'         \item Attribute \code{"group.starts"}: start indices of each group
#'         \item Attribute \code{"group.sizes"}: number of original edges per simplified edge
#'       }
#'   }
#'
#' @details
#' \strong{Method: "shortest-paths"}
#' \itemize{
#'   \item Validates that all origin and destination nodes exist in the network
#'   \item Computes shortest paths from each origin to all destinations using igraph
#'   \item Marks all edges that are traversed by at least one shortest path
#'   \item Returns only the subset of edges that were traversed
#'   \item If \code{nodes} is a data frame with \code{from} and \code{to} columns, paths are computed
#'     from each unique origin to its specified destinations
#' }
#'
#' \strong{Method: "cluster"}
#' \itemize{
#'   \item Requires the graph to have spatial coordinates (\code{FX}, \code{FY}, \code{TX}, \code{TY})
#'   \item If \code{nodes} is provided, these nodes are preserved as cluster centroids
#'   \item Nearby nodes (within \code{radius_km$nodes} km) are assigned to the nearest preserved node
#'   \item Remaining nodes are clustered using \code{\link[leaderCluster]{leaderCluster}} with
#'     \code{radius_km$cluster} as the clustering radius
#'   \item For each cluster, the node closest to the cluster centroid is selected as representative
#'   \item The graph is contracted by mapping all nodes to their cluster representatives
#'   \item Self-loops (edges where both endpoints map to the same cluster) are dropped
#'   \item For undirected graphs (\code{directed = FALSE}), edges are normalized so \code{from < to},
#'     merging opposite-direction edges; for directed graphs, A->B and B->A remain separate
#'   \item Edge attributes are aggregated using \code{\link[collapse]{collap}} (default: mean for
#'     numeric, mode for categorical); customize via \code{\dots}
#' }
#'
#' @examples
#' # Load example data
#' library(flowr)
#'
#' # Create graph from network
#' graph <- linestrings_to_graph(network_gcc)
#'
#' # Get zone nodes
#' nodes_df <- nodes_from_graph(graph, sf = TRUE)
#' nearest_nodes <- nodes_df$node[sf::st_nearest_feature(zones_gcc, nodes_df)]
#'
#' # Method 1: Shortest-paths simplification (keeps only traversed edges)
#' graph_simple <- simplify_network(graph, nearest_nodes,
#'                                  method = "shortest-paths",
#'                                  cost.column = "generalized_cost")
#' nrow(graph_simple)  # Reduced number of edges
#'
#' # With multimodal support: compute paths separately per mode
#' graph_simple_mm <- simplify_network(graph, nearest_nodes,
#'                                     method = "shortest-paths",
#'                                     cost.column = "generalized_cost",
#'                                     by = ~ mode)
#'
#' # Method 2: Cluster-based simplification (contracts graph spatially)
#' \donttest{
#' graph_cluster <- simplify_network(graph, nearest_nodes,
#'                                   method = "cluster",
#'                                   radius_km = list(nodes = 10, cluster = 30))
#'
#' # Cluster with multimodal support: don't merge edges of different modes
#' graph_cluster_mm <- simplify_network(graph, nearest_nodes,
#'                                      method = "cluster",
#'                                      radius_km = list(nodes = 10, cluster = 30),
#'                                      by = ~ mode)
#' }
#'
#' @export
#' @importFrom collapse fnrow ss ckmatch funique.default fmatch gsplit fmin dapply whichv %+=% GRP add_vars seq_row add_stub colorderv %!in% collap get_vars alloc
#' @importFrom igraph graph_from_data_frame delete_vertex_attr igraph_options shortest_paths
#' @importFrom geodist geodist_vec geodist_min
#' @importFrom leaderCluster leaderCluster
simplify_network <- function(graph, nodes, method = c("shortest-paths", "cluster"),
                             directed = FALSE, cost.column = "cost", by = NULL,
                             radius_km = list(nodes = 7, cluster = 20), ...) {

  method <- match.arg(method)

  # Validate graph input
  if (!is.data.frame(graph) || !all(c("from", "to") %in% names(graph)))
    stop("graph must be a data.frame with 'from' and 'to' columns")

  # Validate by argument
  if(length(by)) {
    if(is.call(by)) by <- all.vars(by)
    if(!is.character(by)) stop("by needs to be a one-sided formula or a character vector of column names")
  }

  if (method == "shortest-paths") {

    # Shortest-paths method
    cost <- if(is.character(cost.column) && length(cost.column) == 1L) graph[[cost.column]] else
      if(is.numeric(cost.column) && length(cost.column) == fnrow(graph)) cost.column else
        stop("cost.column needs to be a column name in graph or a numeric vector matching nrow(graph)")

    from_node <- as.integer(graph$from)
    to_node <- as.integer(graph$to)
    all_nodes <- funique.default(c(from_node, to_node), sort = TRUE)

    # Internally use normalized graph node indices
    from_node <- fmatch(from_node, all_nodes)
    to_node <- fmatch(to_node, all_nodes)
    ig <- data.frame(from = from_node, to = to_node) |>
      graph_from_data_frame(directed = directed,
                            vertices = data.frame(name = seq_along(all_nodes))) |>
      delete_vertex_attr("name")

    # Don't return vertex/edge names
    iopt <- igraph_options(return.vs.es = FALSE)
    on.exit(igraph_options(iopt))

    edges_traversed <- integer(fnrow(graph))

    # Pre-process nodes input ONCE (before any group loop)
    if(is.atomic(nodes)) {
      nodes_matched <- ckmatch(nodes, all_nodes, e = "Unknown nodes:")
      nodes_split <- NULL
    } else if(is.data.frame(nodes)) {
      if(!all(c("from", "to") %in% names(nodes)))
        stop("nodes data.frame must have columns 'from' and 'to'")
      nodes_split <- gsplit(ckmatch(nodes$to, all_nodes, e = "Unknown 'to' nodes:"),
                            ckmatch(nodes$from, all_nodes, e = "Unknown 'from' nodes:"),
                            use.g.names = TRUE)
      nodes_matched <- NULL
    } else stop("nodes must be an atomic vector or a data.frame")

    # Helper to compute paths with given cost vector
    compute_paths <- function(cost_vec) {
      if(length(nodes_matched)) {
        if(directed) {
          for (i in nodes_matched) {
            pathsi <- shortest_paths(ig, from = i, to = nodes_matched,
                                     weights = cost_vec, mode = "out", output = "epath")$epath
            .Call(C_mark_edges_traversed, pathsi, edges_traversed)
          }
        } else {
          n <- length(nodes_matched)
          for (i in 1:(n - 1L)) {
            pathsi <- shortest_paths(ig, from = nodes_matched[i], to = nodes_matched[(i + 1L):n],
                                     weights = cost_vec, mode = "out", output = "epath")$epath
            .Call(C_mark_edges_traversed, pathsi, edges_traversed)
          }
        }
      } else {
        from_nodes <- as.integer(names(nodes_split))
        for (i in seq_along(from_nodes)) {
          pathsi <- shortest_paths(ig, from = from_nodes[i], to = nodes_split[[i]],
                                   weights = cost_vec, mode = "out", output = "epath")$epath
          .Call(C_mark_edges_traversed, pathsi, edges_traversed)
        }
      }
    }

    if(length(by)) {
      # Multimodal: iterate over groups, penalizing edges not in the current group
      by_grp <- GRP(graph, by, return.order = FALSE)
      max_cost <- 100 * max(cost, na.rm = TRUE)

      for(grp_idx in seq_len(by_grp$N.groups)) {
        cost_penalized <- cost
        cost_penalized[by_grp$group.id != grp_idx] <- max_cost
        compute_paths(cost_penalized)
      }
    } else {
      # Single mode: compute paths once
      compute_paths(cost)
    }

    edges <- which(edges_traversed > 0L)
    result <- ss(graph, edges, check = FALSE)
    attr(result, "edges") <- edges
    attr(result, "edge_counts") <- edges_traversed[edges]

  } else {

    if (!all(c("FX", "FY", "TX", "TY") %in% names(graph)))
      stop("graph must have columns 'FX', 'FY', 'TX', 'TY' for method = 'cluster'")

    # Extract nodes with coordinates
    nodes_df <- nodes_from_graph(graph, sf = FALSE)
    # Optional node weights
    if(is.numeric(cost.column) && length(cost.column) == fnrow(nodes_df)) nodes_df$weights <- cost.column
    # Cluster method
    cl <- cluster_nodes(nodes_df, nodes, radius_km$nodes, radius_km$cluster)
    # Graph Contraction to Clusters
    result <- contract_edges(graph, nodes = nodes_df, clusters = cl$clusters,
                             centroids = cl$centroids, directed = directed, by = by, ...)

  }

  return(result)
}

# Helper functions for simplify_network() with method = "cluster"
cluster_nodes <- function(nodes, keep,
                          nodes_radius_km = 7,
                          cluster_radius_km = 20) {
  # Nodes to preserve
  if(length(keep)) {
    clusters <- integer(fnrow(nodes))
    keep <- ckmatch(keep, nodes$node, "Unknown nodes to preserve:")
    clusters[keep] <- seq_along(keep)
    # Cluster nodes close to cities
    dmat <- geodist_vec(nodes$X[keep], nodes$Y[keep],
                        nodes$X[-keep], nodes$Y[-keep],
                        measure = "haversine")
    close <- fmin(dmat) < nodes_radius_km * 1000
    if(any(close)) {
      clusters[-keep][close] <- seq_along(keep)[dapply(dmat[, close, drop = FALSE], which.min)]
    }
    ind <- whichv(clusters, 0L)
    if(length(ind)) {
      mat <- cbind(Y = nodes$Y[ind], X = nodes$X[ind])
      weights <- if(length(nodes$weights)) nodes$weights[ind] else alloc(1, nrow(mat))
      res <- leaderCluster(mat, cluster_radius_km, weights, max_iter = 1000L, distance = "haversine")
      clusters[ind] <- res$cluster_id %+=% length(keep)
      centroids <- integer(length(keep) + res$num_clusters)
      centroids[seq_along(keep)] <- keep
      centroids[-seq_along(keep)] <- suppressMessages(ind[geodist_min(res$cluster_centroids[,2:1], mat[,2:1], measure = "haversine", quiet = TRUE)])
    } else centroids <- keep
  } else {
    mat <- cbind(Y = nodes$Y, X = nodes$X)
    weights <- if(length(nodes$weights)) nodes$weights[ind] else alloc(1, nrow(mat))
    res <- leaderCluster(mat, cluster_radius_km, weights, max_iter = 1000L, distance = "haversine")
    clusters <- res$cluster_id
    centroids <- res$cluster_centroids[,2:1]
    dimnames(centroids)[[2L]] <- c("X", "Y")
    centroids <- suppressMessages(geodist_min(centroids, mat[,2:1], measure = "haversine", quiet = TRUE))
  }
  list(clusters = clusters, centroids = nodes$node[centroids])
}

contract_edges <- function(graph, nodes, clusters, centroids, directed = FALSE, by = NULL, ...) {
    node_centroids <- centroids[clusters]
    graph$from <- node_centroids[ckmatch(graph$from, nodes$node)]
    graph$to <- node_centroids[ckmatch(graph$to, nodes$node)]

    # Drop self-loops (edges where both endpoints map to same cluster)
    self_loops <- graph$from == graph$to
    if(any(self_loops)) {
      graph <- ss(graph, !self_loops, check = FALSE)
    }

    # For undirected graphs, normalize edge direction so from < to
    # This merges A->B and B->A into a single edge
    if(!directed) {
      swap <- graph$from > graph$to
      if(any(swap)) {
        tmp <- graph$from[swap]
        graph$from[swap] <- graph$to[swap]
        graph$to[swap] <- tmp
      }
    }

    # Include 'by' columns in grouping to prevent cross-mode consolidation
    g <- GRP(graph, c("from", "to", by))
    nam <- names(graph)
    res <- g$groups
    add_vars(res, pos = "front") <- list(edge = seq_row(res))
    if(any(nam %in% c("FX", "FY", "TX", "TY"))) {
      add_vars(res) <- ss(nodes, ckmatch(res$from, nodes$node), c("X", "Y")) |> add_stub("F")
      add_vars(res) <- ss(nodes, ckmatch(res$to, nodes$node), c("X", "Y")) |> add_stub("T")
    }
    ord <- c("edge", "from", "FX", "FY", "to", "TX", "TY")
    res <- colorderv(res, ord[ord %in% nam])
    if(any(nam %!in% ord)) {
      add_vars(res) <- collap(get_vars(graph, nam[nam %!in% ord]), g, keep.by = FALSE, ...)
    }
    attr(res, "group.id") <- g$group.id
    attr(res, "group.starts") <- g$group.starts
    attr(res, "group.sizes") <- g$group.sizes
    res
  }

# # Contracting graph
# graph_straight |>
#   join(nodes_clustered |> qDF() |>
#          ftransform(qDF(st_coordinates(geometry))) |>
#          fselect(fx = X, fy = Y, lon, lat)) |>
#   ftransform(fx = lon, fy = lat, lon = NULL, lat = NULL) |>
#   join(nodes_clustered |> qDF() |>
#          ftransform(qDF(st_coordinates(geometry))) |>
#          fselect(tx = X, ty = Y, lon, lat)) |>
#   ftransform(tx = lon, ty = lat, lon = NULL, lat = NULL) |>
#   fgroup_by(1:4) |> fsum() |> na_omit() |>
#   st_as_sf(coords = c("fx", "fy"), crs = 4326, sf_column_name = "from") |> qDF() |>
#   st_as_sf(coords = c("tx", "ty"), crs = 4326, sf_column_name = "to") |> qDF() |>
#   fmutate(row = seq_along(from)) |>
#   pivot(c("row", attrib), c("from", "to")) |>
#   collap(~ row, custom = list(ffirst = attrib, st_combine = c("geometry" = "value")),
#          keep.col.order = FALSE) |>
#   st_as_sf(crs = 4326) |>
#   st_cast("LINESTRING") |>
#   fmutate(row = NULL)


#' @title Melt Origin-Destination Matrix to Long Format
#' @description Convert an origin-destination (OD) matrix to a long-format data frame
#'   with columns \code{from}, \code{to}, and \code{flow}.
#'
#' @param od_matrix A numeric matrix with origin-destination flows. Rows represent origins,
#'   columns represent destinations. The matrix should be square (same number of rows and columns).
#' @param nodes (Optional) Numeric or integer vector of node IDs in the graph matching the rows
#'   and columns of the matrix. If provided, must have length equal to both \code{nrow(od_matrix)}
#'   and \code{ncol(od_matrix)}. When \code{nodes} is provided, these IDs are used directly,
#'   ignoring row and column names. This is particularly useful when mapping zone-based OD matrices
#'   to graph node IDs (e.g., using \code{\link[sf]{st_nearest_feature}} to find nearest graph nodes).
#'   If omitted, row and column names (if present) will be used as node IDs, coerced to integer
#'   if possible. If names are not available or cannot be coerced to integer, sequential integers
#'   will be used.
#' @param sort Sort long OD-matrix in ascending order of from and to columns. This can have computational benefits, e.g., when multithreading with \code{method = "AoN"}.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{from} - Origin node ID (integer)
#'     \item \code{to} - Destination node ID (integer)
#'     \item \code{flow} - Flow value (numeric)
#'   }
#'   Only rows with finite, positive flow values are included.
#'
#' @details
#' This function converts a square OD matrix to long format, which is required by
#' \code{\link[=run_assignment]{run_assignment()}}. The behavior depends on whether \code{nodes}
#' is provided:
#'
#' \strong{When \code{nodes} is provided:}
#' \itemize{
#'   \item The \code{nodes} vector is used directly as node IDs for both origins and destinations
#'   \item Row and column names are ignored (but must match if both are present)
#'   \item This is the recommended approach when working with zone-based OD matrices that need to be
#'     mapped to graph nodes, as it ensures the node IDs match those in the graph
#' }
#'
#' \strong{When \code{nodes} is omitted:}
#' \itemize{
#'   \item Row and column names are extracted from the matrix (if available)
#'   \item Names are coerced to integer if possible; if coercion fails or names are missing,
#'     sequential integers are used
#'   \item This approach works well when the matrix row/column names already correspond to graph node IDs
#' }
#'
#' In both cases, the function:
#' \itemize{
#'   \item Creates a long-format data frame with all origin-destination pairs
#'   \item Filters out non-finite and zero flow values
#' }
#'
#' The function is useful for converting OD matrices (such as those in \code{\link{od_matrices_gcc}})
#' to the long format required by \code{\link[=run_assignment]{run_assignment()}}.
#'
#' @seealso \code{\link{od_matrices_gcc}}, \code{\link{zones_gcc}}, \code{\link[=nodes_from_graph]{nodes_from_graph()}},
#'   \code{\link[=run_assignment]{run_assignment()}}, \link{flowr-package}
#'
#' @examples
#' library(flowr)
#' library(sf)
#'
#' # Example 1: Using nodes argument to map zones to graph nodes
#' # Load Graph
#' graph <- linestrings_to_graph(network_gcc)
#' nodes <- nodes_from_graph(graph, sf = TRUE)
#'
#' # Map Zones to Nodes
#' nearest_nodes <- nodes$node[st_nearest_feature(zones_gcc, nodes)]
#'
#' # Convert OD matrix with node mapping
#' od_matrix_long <- melt_od_matrix(od_matrices_gcc$container, nodes = nearest_nodes)
#' head(od_matrix_long)
#'
#' # Example 2: Using matrix row/column names (when they match graph node IDs)
#' # If matrix names already correspond to graph nodes, nodes argument can be omitted
#' od_matrix_long2 <- melt_od_matrix(od_matrices_gcc$container)
#' head(od_matrix_long2)
#'
#' @export
#' @importFrom collapse vec fsubset seq_row seq_col all_identical roworder
melt_od_matrix <- function(od_matrix, nodes = NULL, sort = TRUE) {
  if(!is.matrix(od_matrix)) stop("od_matrix must be a matrix")
  if(!is.numeric(od_matrix)) stop("od_matrix must be numeric")

  dn <- dimnames(od_matrix)

  if(length(nodes)) {
    if(length(nodes) != nrow(od_matrix)) stop("Length of nodes must match number of rows in od_matrix")
    if(length(nodes) != ncol(od_matrix)) stop("Length of nodes must match number of columns in od_matrix")
    if(length(dn) && !all_identical(dn)) stop("Row- and column-names of od_matrix must match if nodes are provided")
    row_ids <- col_ids <- nodes
  } else {
    # Try to coerce row names to integer
    if(!is.null(dn[[1L]])) {
      row_ids <- as.integer(dn[[1L]])
      if(anyNA(row_ids)) row_ids <- seq_row(od_matrix)
    } else row_ids <- seq_row(od_matrix)
    # Try to coerce column names to integer
    if(!is.null(dn[[2L]])) {
      col_ids <- as.integer(dn[[2L]])
      if(anyNA(col_ids)) col_ids <- seq_len(ncol(od_matrix))
    } else col_ids <- seq_col(od_matrix)
  }

  # Create long format data frame
  od_matrix_long <- data.frame(
    from = rep(row_ids, ncol(od_matrix)),
    to = rep(col_ids, each = nrow(od_matrix)),
    flow = vec(od_matrix)
  ) |>
    fsubset(is.finite(flow) & flow > 0)
  if(sort) roworder(od_matrix_long, from, to) else od_matrix_long
}




