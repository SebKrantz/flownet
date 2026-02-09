# Changelog

## flownet 0.2.0

- Release blog post at:
  <https://sebkrantz.github.io/Rblog/2026/02/09/introducing-flownet-efficient-transport-modeling-in-r/>

- Fixed bug in run_assignment with `return.extra = "edges"` where edge
  indices were incorrectly returned. Due to zero indexing in C
  vs. 1-indexing in R they where offset by one, thus, in *flownet*
  versions \<= 0.1.2, the correct edge indices can be retained by
  subtracting 1 from them.

- Improved [Step
  7](https://sebkrantz.github.io/flownet/reference/run_assignment.html#route-enumeration-algorithm)
  (elimination of path with duplicate edges) in the route enumeration
  algorithm to properly handle directed paths, i.e., candidate paths
  where an intermediate node is approached via an edge and departed from
  via the same edge but in a different direction, are now also removed.

- Added option `return.extra = "PSF"` which adds `"path_size_factors"`
  to the results object. This should be useful to calibrate the PSL
  model.

- Added option `return.extra = "eweights"` which adds `"egde_weights"`
  to the results object. Edge weights are the sum of the path
  probabilities across all paths using that edge. These weights are
  computed efficiently at C-level, with minimal overhead, if requested.

- Added verbosity to
  [`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)
  (default `verbose = TRUE`).

- Reordered elements in results object (class ‘flownet’) in a natural
  way: providing first the final flows, followed by path and edge-level
  additional results if requested through `return.extra`.

- Added more unit tests covering `return.extra` options.

- Added variable labels to included `africa_*` datasets - try
  [`collapse::namlab()`](https://fastverse.org/collapse/reference/small-helpers.html)
  on them.

- Small improvements to recursive logic in
  [`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md).

- Bumped *kit* dependency to 0.0.21 (for `fpmin()` and `fpmax()`).

## flownet 0.1.2

CRAN release: 2026-02-03

- In
  [`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md):
  argument `consolidate` was renamed to `contract` for improved clarity,
  while ensuring backwards compatibility.

- Minor improvements to documentation and vignette.

## flownet 0.1.1

CRAN release: 2026-01-27

- Implemented minor CRAN comments

### Documentation Improvements

- Restructured
  [`run_assignment()`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
  `@details` section using `\subsection{}` for better organization:
  - Separate sections for AoN method, PSL method, route enumeration
    algorithm, and coordinate-based filtering
- Added `return.extra` parameter documentation as a clear table showing
  which options are available for each method
- Clarified that `paths` returns edge indices, not node indices
- Added comprehensive examples including:
  - PSL method usage with `nthreads` parameter
  - Trade flow disaggregation workflow (country-to-city level)
- Added academic references for the Path-Sized Logit model (Ben-Akiva &
  Bierlaire, 1999)
- Added link to AequilibraE (Python) documentation for additional PSL
  resources
- Improved error messages to be more informative (show which columns are
  missing, what class was received)
- Increased vignette table of contents depth for better navigation
- Added citation for the PSL model in vignette

### Testing

- Added testthat test suite with 150 tests covering:
  - [`run_assignment()`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
    (AoN and PSL methods)
  - [`normalize_graph()`](https://sebkrantz.github.io/flownet/reference/normalize_graph.md),
    [`nodes_from_graph()`](https://sebkrantz.github.io/flownet/reference/nodes_from_graph.md),
    [`distances_from_graph()`](https://sebkrantz.github.io/flownet/reference/distances_from_graph.md)
  - [`linestrings_to_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_to_graph.md),
    [`linestrings_from_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_from_graph.md),
    [`create_undirected_graph()`](https://sebkrantz.github.io/flownet/reference/create_undirected_graph.md)
  - [`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md),
    [`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)
  - [`melt_od_matrix()`](https://sebkrantz.github.io/flownet/reference/melt_od_matrix.md)
- Added test coverage workflow

### Minor Improvements

- Better penalization strategy for multimodal networks in
  [`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)
  shortest-paths method
- Enforced numeric type conversion for cost columns
- Various minor fixes and improvements

------------------------------------------------------------------------

## flownet 0.1.0 (Initial Release)

### Major Features

#### Traffic Assignment and Route Enumeration

- **[`run_assignment()`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)**:
  Core function for traffic assignment using the path-sized logit (PSL)
  model
  - Efficient route enumeration algorithm that generates sensible
    alternative routes between origin-destination pairs
  - Accounts for route overlap using path-sized logit probabilities
  - Supports both directed and undirected graphs
  - Configurable detour factors and angle constraints for route
    selection
  - Optional return of detailed path information, edge counts, costs,
    and weights
  - High-performance C implementation for critical path operations
  - Multithreading (asynchronous parallelism) using
    [`mirai`](https://github.com/r-lib/mirai).

#### Network Processing

- **[`linestrings_to_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_to_graph.md)**:
  Convert LINESTRING geometries (sf objects) to graph data frames
  - Extracts node coordinates and creates edge representations
  - Optional coordinate rounding for precision handling
  - Preserves additional columns from input data
  - Automatic computation of edge lengths
- **[`create_undirected_graph()`](https://sebkrantz.github.io/flownet/reference/create_undirected_graph.md)**:
  Convert directed graphs to undirected graphs
  - Normalizes edge directions and aggregates duplicate connections
  - Flexible aggregation using
    [`collapse::collap()`](https://fastverse.org/collapse/reference/collap.html)
    with customizable functions
  - Preserves spatial coordinates and line identifiers
- **[`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md)**:
  Simplify network topology by removing intermediate nodes
  - Removes nodes that occur exactly twice, merging connecting edges
  - Optional removal of loops, duplicate edges, and singleton edges
  - Recursive consolidation for complete simplification
  - Weighted aggregation of edge attributes
  - Preserves coordinate information when present
- **[`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)**:
  Simplify networks using shortest-paths or spatial clustering methods
  - **Method “shortest-paths”**: Filters network to edges used in
    shortest paths between specified nodes
  - **Method “cluster”**: Spatially clusters nodes using `leaderCluster`
    algorithm and contracts the graph

#### Graph Utilities

- **[`normalize_graph()`](https://sebkrantz.github.io/flownet/reference/normalize_graph.md)**:
  Normalize node IDs to consecutive integers starting from 1
  - Essential for graph algorithms requiring sequential node IDs
  - Preserves graph structure while remapping identifiers
- **[`nodes_from_graph()`](https://sebkrantz.github.io/flownet/reference/nodes_from_graph.md)**:
  Extract unique nodes with coordinates from graph
  - Returns data frame or sf POINT object with node locations
  - Useful for mapping zones to network nodes
- **[`linestrings_from_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_from_graph.md)**:
  Convert graph data frames back to LINESTRING geometries
  - Inverse operation of
    [`linestrings_to_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_to_graph.md)
  - Preserves all graph attributes in output sf object
- **[`distances_from_graph()`](https://sebkrantz.github.io/flownet/reference/distances_from_graph.md)**:
  Compute distance matrices for all node pairs
  - Uses igraph for efficient shortest path distance computation
  - Supports both directed and undirected graphs

#### OD Matrix Utilities

- **[`melt_od_matrix()`](https://sebkrantz.github.io/flownet/reference/melt_od_matrix.md)**:
  Convert origin-destination matrices to long format
  - Transforms matrix format to edge list suitable for traffic
    assignment
  - Handles missing values and zero flows

#### Example Data

- **`africa_network`**: A road transport network with 2,825 LINESTRING
  features representing existing roads (2,344 edges) and proposed new
  links (481 edges). Each edge includes attributes such as distance,
  travel duration, border crossing costs, terrain ruggedness, and road
  upgrade costs.

- **`africa_cities_ports`**: 453 African cities with population \>
  100,000 and international ports. Includes population data, capital
  status, and port cargo outflows.

- **`africa_segments`**: 14,358 raw network segments representing
  intersected road routes. Useful for demonstrating network
  consolidation and simplification functions.

- **`africa_trade`**: Bilateral trade flows between 47 African countries
  aggregated by HS section (21 product categories). Values represent
  annual averages over 2012-2022.

The `africa_network`, `africa_cities_ports`, and `africa_segments`
datasets are from Krantz, S. (2024). [Optimal Investments in Africa’s
Road Network](https://doi.org/10.1596/1813-9450-10893). Policy Research
Working Paper 10893. World Bank. Replication materials are available at
[github.com/SebKrantz/OptimalAfricanRoads](https://github.com/SebKrantz/OptimalAfricanRoads).

### Technical Details

- High-performance C implementations for path-sized logit computations
- Efficient memory management for large networks
- Integration with `collapse` package for fast data transformations
- Uses `igraph` for graph operations and shortest path algorithms
- Leverages `geodist` for fast geodesic distance computations
- Leverages `leaderCluster` for efficient spatial clustering
- Comprehensive documentation with examples and vignettes

### Documentation

- Complete function documentation with roxygen2
- README with quick start guide and examples
- Introduction vignette demonstrating package workflow
- Package-level documentation in `?flownet-package`

### Dependencies

- **R** (\>= 4.1)
- **collapse** (\>= 2.1.5) - Fast data transformations and memory
  efficient programming
- **kit** (\>= 0.0.5) - Fast tabulation and vectorized switches
- **igraph** (\>= 2.1.4) - Graph operations and shortest path algorithms
- **sf** (\>= 1.0.0) - Spatial data handling
- **geodist** (\>= 0.1.1) - Fast geodesic distance computations
- **leaderCluster** (\>= 1.5.0) - Efficient spatial clustering
  algorithms
- **mirai** (\>= 2.5.2) - Asynchronous parallelism for R
- **progress** (\>= 1.2.3) - Progress bars for long-running operations

### Suggested Packages

- **fastverse** (\>= 0.3.4) - Efficient data manipulation workflow
- **mapview** (\>= 2.11.2) - Interactive visualization of results
- **tmap** (\>= 4.0) - Static visualization of results

### License

GPL-3
