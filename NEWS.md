# flowr 0.1.0 (Initial Release)

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

- **`simplify_network()`**: Simplify networks by keeping only traversed edges
  - Filters network to edges used in shortest paths between OD pairs
  - Works with both sf LINESTRING objects and graph data frames
  - Efficient shortest path computation using igraph

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
- **`network_gcc`**: Multimodal transport network for the Gulf Cooperation Council (GCC) region
  - Includes road, rail, and maritime connections
  - Pre-processed and ready for analysis

- **`od_matrices_gcc`**: Origin-destination matrices for five cargo types
  - Container, Drybulk, Liquidbulk, General, and HighValue cargo flows
  - Multiple time periods (2019, 2030, 2040)

- **`zones_gcc`**: Zone locations and descriptions for OD-matrix locations
  - Geographic coordinates and zone metadata
  - Compatible with GCC network and OD matrices

## Technical Details
- High-performance C implementations for path-sized logit computations
- Efficient memory management for large networks
- Integration with `collapse` package for fast data transformations
- Uses `igraph` for graph operations and shortest path algorithms
- Leverages `geodist` for fast geodesic distance computations
- Comprehensive documentation with examples and vignettes

## Documentation
- Complete function documentation with roxygen2
- README with quick start guide and examples
- Introduction vignette demonstrating package workflow
- Package-level documentation in `?flowr-package`

## Dependencies
- **R** (>= 3.5)
- **collapse** (>= 2.1.5) - Fast data transformations and memory efficient programming
- **kit** (>= 0.0.5) - Fast tabulation and vectorized switches
- **igraph** (>= 2.1.4) - Graph operations and shortest path algorithms
- **sf** (>= 1.0.0) - Spatial data handling
- **geodist** (>= 0.1.1) - Fast geodesic distance computations
- **mirai** (>= 2.5.2) - Asynchronous parallelism for R
- **progress** (>= 1.2.3) - Progress bars for long-running operations

## Suggested Packages
- **fastverse** (>= 0.3.4) - Enhanced data manipulation workflow
- **mapview** (>= 2.11.2) - Interactive visualization of results
- **testthat** (>= 3.0.0) - Unit testing framework

## License
GPL-3


