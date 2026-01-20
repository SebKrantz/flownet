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
- **R** (>= 3.5)
- **collapse** (>= 2.1.5) - Fast data transformations and memory efficient programming
- **kit** (>= 0.0.5) - Fast tabulation and vectorized switches
- **igraph** (>= 2.1.4) - Graph operations and shortest path algorithms
- **sf** (>= 1.0.0) - Spatial data handling
- **geodist** (>= 0.1.1) - Fast geodesic distance computations
- **leaderCluster** (>= 0.1.5) - Efficient spatial clustering algorithms
- **mirai** (>= 2.5.2) - Asynchronous parallelism for R
- **progress** (>= 1.2.3) - Progress bars for long-running operations

## Suggested Packages
- **fastverse** (>= 0.3.4) - Efficient data manipulation workflow
- **mapview** (>= 2.11.2) - Interactive visualization of results
- **tmap** (>= 4.0) - Static visualization of results

## License
GPL-3


