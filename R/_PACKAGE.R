#' Transport Modeling
#'
#' @description
#'
#' \emph{flownet} provides efficient tools for transportation modeling, specifically route
#' enumeration and traffic assignment tasks. The package implements the path-sized logit model for
#' traffic assignment and provides utilities for network processing.
#'
#'
#' \strong{Network Processing}
#'
#'  \code{\link[=linestrings_to_graph]{linestrings_to_graph()}} --- Convert LINESTRING geometries to graph\cr
#'  \code{\link[=create_undirected_graph]{create_undirected_graph()}} --- Convert directed graph to undirected\cr
#'  \code{\link[=consolidate_graph]{consolidate_graph()}} --- Consolidate graph by removing intermediate nodes\cr
#'  \code{\link[=simplify_network]{simplify_network()}} --- Simplify network graph\cr
#'
#'
#' \strong{Traffic Assignment}
#'
#'  \code{\link[=run_assignment]{run_assignment()}} --- Run traffic assignment using path-sized logit model\cr
#'
#'
#' \strong{Graph Utilities}
#'
#'  \code{\link[=normalize_graph]{normalize_graph()}} --- Normalize node IDs to consecutive integers\cr
#'  \code{\link[=nodes_from_graph]{nodes_from_graph()}} --- Extract unique nodes from graph\cr
#'  \code{\link[=linestrings_from_graph]{linestrings_from_graph()}} --- Convert graph to LINESTRING geometries\cr
#'  \code{\link[=distances_from_graph]{distances_from_graph()}} --- Compute distance matrix from graph\cr
#'
#'
#' \strong{OD Matrix Utilities}
#'
#'  \code{\link[=melt_od_matrix]{melt_od_matrix()}} --- Convert origin-destination matrix to long format\cr
#'
#'
#' \strong{Data}
#'
#'  \code{\link{africa_trade}} --- Average BACI HS96 2012-22 trade flows by section between 47 continental African countries\cr
#'  \code{\link{africa_cities_ports}} --- The 453 largest (port-)cities in continental Africa within a 70km radius - from Krantz (2024), \doi{10.1596/1813-9450-10893}\cr
#'  \code{\link{africa_network}} --- African continental road network + extensions to optimally connect the 453 cities - from Krantz (2024), \doi{10.1596/1813-9450-10893}\cr
#'  \code{\link{africa_segments}} --- Primary segments derived from OpenStreetMap routes between the 453 cities - from Krantz (2024), \doi{10.1596/1813-9450-10893}\cr
#'
#'  Replication materials: \url{https://github.com/SebKrantz/OptimalAfricanRoads}
#'
#' @details
#' The package uses efficient C implementations for critical path operations and leverages:
#' \itemize{
#'   \item \code{collapse} - Fast data transformations
#'   \item \code{geodist} - Fast geodesic distance computations
#'   \item \code{igraph} - Graph operations and shortest path algorithms
#'   \item \code{leaderCluster} - Fast spatial clustering for network simplification
#' }
#'
#' @author Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch} and Kamol Roy \email{kamol.roy08@gmail.com}
#' @name flownet-package
#' @aliases flownet
#' @useDynLib flownet, .registration = TRUE
NULL

