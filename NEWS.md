# flownet 0.3.0

*flownet* 0.3.0 introduces one new function for critical link analysis and substantive performance improvements to `consolidate_graph()` and `simplify_network()`. 

## New Functions

- **`critical_link_analysis()`**: New function that computes edge-level detour costs and vulnerability ratios for any graph. For each edge `(u, v)` it finds the shortest alternative path from `u` to `v` after removing that edge, without repeated graph deletion. The underlying C implementation runs a multi-label Dijkstra once per unique source node, tracking the two best distances that begin via different first edges. Parallel edges are handled as distinct alternatives; edges with no detour receive `Inf`. Supports directed and undirected graphs, and parallelises per-source Dijkstra runs across `nthreads` OpenMP threads.

## `consolidate_graph()` Improvements

- Singleton-edge removal now uses a linear-time peel (`drop_singletons_linear()`), avoiding repeated full-graph scans on large networks.

- Degree-2 chain contraction is now performed at C level (`C_contract_linear_nodes`) by default. The C routine walks each chain once using dense array indexing with no hash table. Set `options(flownet.use_c_contraction = FALSE)` to force the R-only fallback.

- `consolidate_graph_core()` remaps `from` / `to` to dense internal node indices (`funique.default()` + `fmatch()`) before contraction and maps aggregated endpoints back to original node IDs, so results and downstream joins are unchanged.

- `consolidate_graph_core()` now handles topology-only graphs (no coordinate or attribute columns): when no nodes are kept for aggregation, it correctly returns a data frame of contracted groups.

- When `verbose = TRUE`, phase timings are reported for singleton peel, orientation, contraction, and aggregation steps.

- More `testthat` coverage for edge cases: recursion, `keep.nodes`, long chains, C vs. R contraction equivalence, and singleton-peel equivalence.

## `simplify_network()` Improvements

- New `nthreads` parameter (default `1`): parallelises the per-chunk shortest-path computations in `method = "shortest-paths"` using `mirai` daemons, matching the multithreading interface of `run_assignment()`.

- New `cluster.algo.dbscan` parameter (default `FALSE`): when `TRUE`, spatial clustering of remaining nodes uses `dbscan::dbscan()` instead of `leaderCluster`. DBSCAN is substantially faster for larger datasets. Requires the **dbscan** package.

- New `nodes.algo.rann` parameter (default `FALSE`): when `TRUE`, node-to-kept-node assignment uses **RANN**'s fast k-d tree nearest-neighbour search on 3-D Cartesian coordinates (converted from lon/lat via `lonlat_to_xyz()`). Requires the **RANN** package.

- The simplified graph returned by `method = "cluster"` now carries a `keep.edges` attribute containing the integer indices of the edges retained from the original graph.

## Result Attributes

- The grouping/edge-tracking attributes attached to results are now consistent and documented. `consolidate_graph()` attaches `keep.edges` and `group.id` together (aligned, same length) for `recursive = "none"` or `"partial"`; with the default `recursive = "full"` they are omitted (the mapping would span multiple passes), which is now stated in the documentation. `create_undirected_graph()` now also exposes `group.id` and `group.sizes` (previously only `group.starts`, and undocumented), and `consolidate_graph()` strips any stale grouping attributes carried in from upstream operations.

# flownet 0.2.2

- Fixed issue in `consolidate_graph()` which used to modify columns (`from` and `to` in-place). Users in older versions are advised to input a `data.table::copy()` of the graph to retain it. 

- Fixes issue with multithreading for newer versions of *mirai* (or R). Thanks @kent37 (#69).

# flownet 0.2.1

- `angle.max` constraint in `run_assignment()` is now two-sided (angle measured from origin and destination node against the straight line between them), rather than just one-sided (from origin). Also, the implementation is slightly more efficient.

# flownet 0.2.0

- Release blog post at: https://sebkrantz.github.io/Rblog/2026/02/09/introducing-flownet-efficient-transport-modeling-in-r/ 

- Fixed bug in run_assignment with `return.extra = "edges"` where edge indices were incorrectly returned. Due to zero indexing in C vs. 1-indexing in R they where offset by one, thus, in *flownet* versions <= 0.1.2, the correct edge indices can be retained by subtracting 1 from them.

- Improved [Step 7](https://sebkrantz.github.io/flownet/reference/run_assignment.html#route-enumeration-algorithm) (elimination of path with duplicate edges) in the route enumeration algorithm to properly handle directed paths, i.e., candidate paths where an intermediate node is approached via an edge and departed from via the same edge but in a different direction, are now also removed. 

- Added option `return.extra = "PSF"` which adds `"path_size_factors"` to the results object. This should be useful to calibrate the PSL model.

- Added option `return.extra = "eweights"` which adds `"egde_weights"` to the results object. Edge weights are the sum of the path probabilities across all paths using that edge. These weights are computed efficiently at C-level, with minimal overhead, if requested.

- Added verbosity to `simplify_network()` (default `verbose = TRUE`). 

- Reordered elements in results object (class 'flownet') in a natural way: providing first the final flows, followed by path and edge-level additional results if requested through `return.extra`.

- Added more unit tests covering `return.extra` options.

- Added variable labels to included `africa_*` datasets - try `collapse::namlab()` on them.

- Small improvements to recursive logic in `consolidate_graph()`.

- Bumped *kit* dependency to 0.0.21 (for `fpmin()` and `fpmax()`).


# flownet 0.1.2

- In `consolidate_graph()`: argument `consolidate` was renamed to `contract` for improved clarity, while ensuring backwards compatibility.

- Minor improvements to documentation and vignette.

# flownet 0.1.1

- Implemented minor CRAN comments

## Documentation Improvements

- Restructured `run_assignment()` `@details` section using `\subsection{}` for better organization:
  - Separate sections for AoN method, PSL method, route enumeration algorithm, and coordinate-based filtering
- Added `return.extra` parameter documentation as a clear table showing which options are available for each method
- Clarified that `paths` returns edge indices, not node indices
- Added comprehensive examples including:
  - PSL method usage with `nthreads` parameter
  - Trade flow disaggregation workflow (country-to-city level)
- Added academic references for the Path-Sized Logit model (Ben-Akiva & Bierlaire, 1999)
- Added link to AequilibraE (Python) documentation for additional PSL resources
- Improved error messages to be more informative (show which columns are missing, what class was received)
- Increased vignette table of contents depth for better navigation
- Added citation for the PSL model in vignette

## Testing

- Added testthat test suite with 150 tests covering:
  - `run_assignment()` (AoN and PSL methods)
  - `normalize_graph()`, `nodes_from_graph()`, `distances_from_graph()`
  - `linestrings_to_graph()`, `linestrings_from_graph()`, `create_undirected_graph()`
  - `consolidate_graph()`, `simplify_network()`
  - `melt_od_matrix()`
- Added test coverage workflow

## Minor Improvements

- Better penalization strategy for multimodal networks in `simplify_network()` shortest-paths method
- Enforced numeric type conversion for cost columns
- Various minor fixes and improvements

---

# flownet 0.1.0 (Initial Release)

## Major Features

### Traffic Assignment and Route Enumeration
- **`run_assignment()`**: Core function for traffic assignment using the path-sized logit (PSL) model
  - Efficient route enumeration algorithm that generates sensible alternative routes between origin-destination pairs
  - Accounts for route overlap using path-sized logit probabilities
  - Supports both directed and undirected graphs
  - Configurable detour factors and angle constraints for route selection
  - Optional return of detailed path information, edge counts, costs, and weights
  - High-performance C implementation for critical path operations
  - Multithreading (asynchronous parallelism) using [`mirai`](https://github.com/r-lib/mirai). 

### Network Processing
- **`linestrings_to_graph()`**: Convert LINESTRING geometries (sf objects) to graph data frames
  - Extracts node coordinates and creates edge representations
  - Optional coordinate rounding for precision handling
  - Preserves additional columns from input data
  - Automatic computation of edge lengths

- **`create_undirected_graph()`**: Convert directed graphs to undirected graphs
  - Normalizes edge directions and aggregates duplicate connections
  - Flexible aggregation using `collapse::collap()` with customizable functions
  - Preserves spatial coordinates and line identifiers

- **`consolidate_graph()`**: Simplify network topology by removing intermediate nodes
  - Removes nodes that occur exactly twice, merging connecting edges
  - Optional removal of loops, duplicate edges, and singleton edges
  - Recursive consolidation for complete simplification
  - Weighted aggregation of edge attributes
  - Preserves coordinate information when present

- **`simplify_network()`**: Simplify networks using shortest-paths or spatial clustering methods
  - **Method "shortest-paths"**: Filters network to edges used in shortest paths between specified nodes
  - **Method "cluster"**: Spatially clusters nodes using `leaderCluster` algorithm and contracts the graph

### Graph Utilities
- **`normalize_graph()`**: Normalize node IDs to consecutive integers starting from 1
  - Essential for graph algorithms requiring sequential node IDs
  - Preserves graph structure while remapping identifiers

- **`nodes_from_graph()`**: Extract unique nodes with coordinates from graph
  - Returns data frame or sf POINT object with node locations
  - Useful for mapping zones to network nodes

- **`linestrings_from_graph()`**: Convert graph data frames back to LINESTRING geometries
  - Inverse operation of `linestrings_to_graph()`
  - Preserves all graph attributes in output sf object

- **`distances_from_graph()`**: Compute distance matrices for all node pairs
  - Uses igraph for efficient shortest path distance computation
  - Supports both directed and undirected graphs

### OD Matrix Utilities
- **`melt_od_matrix()`**: Convert origin-destination matrices to long format
  - Transforms matrix format to edge list suitable for traffic assignment
  - Handles missing values and zero flows

### Example Data
- **`africa_network`**: A road transport network with 2,825 LINESTRING features representing existing roads (2,344 edges) and proposed new links (481 edges). Each edge includes attributes such as distance, travel duration, border crossing costs, terrain ruggedness, and road upgrade costs.

- **`africa_cities_ports`**: 453 African cities with population > 100,000 and international ports. Includes population data, capital status, and port cargo outflows.

- **`africa_segments`**: 14,358 raw network segments representing intersected road routes. Useful for demonstrating network consolidation and simplification functions.

- **`africa_trade`**: Bilateral trade flows between 47 African countries aggregated by HS section (21 product categories). Values represent annual averages over 2012-2022.

The `africa_network`, `africa_cities_ports`, and `africa_segments` datasets are from Krantz, S. (2024). [Optimal Investments in Africa's Road Network](https://doi.org/10.1596/1813-9450-10893). Policy Research Working Paper 10893. World Bank. Replication materials are available at [github.com/SebKrantz/OptimalAfricanRoads](https://github.com/SebKrantz/OptimalAfricanRoads).

## Technical Details
- High-performance C implementations for path-sized logit computations
- Efficient memory management for large networks
- Integration with `collapse` package for fast data transformations
- Uses `igraph` for graph operations and shortest path algorithms
- Leverages `geodist` for fast geodesic distance computations
- Leverages `leaderCluster` for efficient spatial clustering
- Comprehensive documentation with examples and vignettes

## Documentation
- Complete function documentation with roxygen2
- README with quick start guide and examples
- Introduction vignette demonstrating package workflow
- Package-level documentation in `?flownet-package`

## Dependencies
- **R** (>= 4.1)
- **collapse** (>= 2.1.5) - Fast data transformations and memory efficient programming
- **kit** (>= 0.0.5) - Fast tabulation and vectorized switches
- **igraph** (>= 2.1.4) - Graph operations and shortest path algorithms
- **sf** (>= 1.0.0) - Spatial data handling
- **geodist** (>= 0.1.1) - Fast geodesic distance computations
- **leaderCluster** (>= 1.5.0) - Efficient spatial clustering algorithms
- **mirai** (>= 2.5.2) - Asynchronous parallelism for R
- **progress** (>= 1.2.3) - Progress bars for long-running operations

## Suggested Packages
- **fastverse** (>= 0.3.4) - Efficient data manipulation workflow
- **mapview** (>= 2.11.2) - Interactive visualization of results
- **tmap** (>= 4.0) - Static visualization of results

## License
GPL-3


