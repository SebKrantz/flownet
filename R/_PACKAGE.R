#' Efficient Transport Modeling
#'
#' @description
#'
#' \emph{flownet} provides efficient tools for transportation modeling in R, supporting network processing, route enumeration,
#' and traffic assignment tasks. It implements the path-sized logit (PSL) model for traffic assignment and provides powerful
#' utilities for network processing/preparation.
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
#' @author Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch} and Kamol Roy \email{kamol.roy08@gmail.com}
#' @name flownet-package
#' @aliases flownet
#' @useDynLib flownet, .registration = TRUE
NULL

