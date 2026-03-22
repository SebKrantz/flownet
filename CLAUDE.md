# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Overview

`flownet` is an R package for transport modeling that implements network
processing, route enumeration, and traffic assignment algorithms. The
package is maintained by CPCS transport consultants and provides
high-performance tools through a combination of R, fastverse packages,
and custom C implementations.

## Development Commands

### Package Building and Testing

``` r
# Build and install the package
devtools::install()

# Run tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-assignment.R")

# Check package (like R CMD check)
devtools::check()

# Build documentation
devtools::document()

# Build vignettes
devtools::build_vignettes()
```

### C Code Compilation

The package includes custom C implementations in `src/` for
performance-critical operations: - `path_sized_logit.c` - Core
path-sized logit algorithm - `utils.c` - Utility functions for path
operations - `init.c` - Registration of C functions

After modifying C code:

``` r
devtools::clean_dll()  # Clean compiled objects
devtools::load_all()   # Recompile and reload
```

### Testing Commands

``` bash
# Run R CMD check from terminal
R CMD build . && R CMD check flownet_*.tar.gz

# Run tests with coverage
Rscript -e "covr::package_coverage()"
```

## Architecture

### Core Algorithms

The package centers on two traffic assignment methods:

1.  **All-or-Nothing (AoN)**: Fast assignment that allocates all flow to
    shortest paths. Implementation uses batched shortest path
    computation grouped by origin nodes for efficiency.

2.  **Path-Sized Logit (PSL)**: Sophisticated stochastic assignment
    considering multiple alternative routes with overlap correction. The
    algorithm:

    - Enumerates candidate routes using distance matrices and geographic
      filtering
    - Filters routes by detour factor (`detour.max`) and angle
      constraints (`angle.max`)
    - Computes path-size factors to penalize overlapping routes
    - Assigns flows probabilistically based on generalized costs

### Route Enumeration Strategy

The PSL method uses a two-stage approach to avoid computing implausible
paths:

1.  **Pre-selection**: Uses precomputed distance matrices to identify
    promising intermediate nodes based on total cost
    (origin→intermediate→destination)
2.  **Geographic filtering**: When `angle.max` is specified and
    coordinates are available, filters nodes using the triangle equation
    with geodesic distances
3.  **Path computation**: Only computes actual paths for pre-selected
    candidates, then filters duplicates

This strategy dramatically reduces computational cost compared to
enumerating all possible paths.

### Network Processing Pipeline

Typical workflow for processing spatial networks:

1.  **[`linestrings_to_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_to_graph.md)**:
    Convert sf LINESTRING geometries to graph data frames with node
    coordinates (FX, FY, TX, TY)
2.  **[`create_undirected_graph()`](https://sebkrantz.github.io/flownet/reference/create_undirected_graph.md)**:
    Normalize edge directions and aggregate bidirectional links
3.  **[`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md)**:
    Contract intermediate nodes (degree-2 nodes) recursively, merging
    edges while preserving important nodes
4.  **[`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)**:
    Further reduce network size using either:
    - Shortest-paths method: Keep only edges traversed by shortest paths
      between key nodes
    - Cluster method: Spatially cluster nodes using leaderCluster and
      contract graph

### Parallelization

The package uses `mirai` for asynchronous parallelism: - Work is split
across OD-pairs and distributed to daemon processes - Each daemon
processes a subset independently - Results are aggregated after
collection - The `nthreads` parameter controls the number of parallel
workers

### C Integration

Performance-critical operations are implemented in C and called via
[`.Call()`](https://rdrr.io/r/base/CallExternal.html): -
`C_compute_path_sized_logit`: Core PSL computation including overlap
detection and probability calculation - `C_check_path_duplicates`:
Detect paths with duplicate edges (invalid routes) -
`C_assign_flows_to_paths`: Batch flow assignment for AoN method -
`C_mark_edges_traversed`: Track edge usage for simplify_network -
`C_set_vector_elt`: Efficient list element assignment

## Code Organization

### Main Source Files

- **`R/assignment.R`** (802 lines): Contains
  [`run_assignment()`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
  and both AoN and PSL core functions. The PSL implementation handles
  distance matrix chunking for large networks and coordinates with C
  functions for path overlap calculations.

- **`R/utils.R`** (1273 lines): Network processing functions including:

  - Graph conversion utilities (`linestrings_to_graph`,
    `nodes_from_graph`, etc.)
  - [`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md):
    Recursive node contraction with sophisticated degree tracking
  - [`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md):
    Two methods for network reduction
  - [`melt_od_matrix()`](https://sebkrantz.github.io/flownet/reference/melt_od_matrix.md):
    OD matrix format conversion

- **`R/data.R`**: Documentation for included datasets (Africa network,
  cities, trade flows)

### C Source Files

- **`src/path_sized_logit.c`**: Implements path-size factor computation
  and flow assignment
- **`src/utils.c`**: Helper functions for path operations
- **`src/init.c`**: C function registration for .Call interface

### Tests

Tests are organized by functionality in `tests/testthat/`: -
`test-assignment.R`: Traffic assignment methods - `test-graph-utils.R`:
Graph utility functions - `test-network-processing.R`: Network
conversion and processing - `test-od-matrix.R`: OD matrix operations -
`test-consolidation.R`: Graph consolidation

## Dependencies

### Core Dependencies

- **collapse** (≥ 2.1.5): Fast data transformations, used extensively
  for grouping, aggregation, and memory-efficient operations
- **igraph** (≥ 2.1.4): Shortest path algorithms via Dijkstra
- **sf** (≥ 1.0.0): Spatial data handling for LINESTRING networks
- **geodist** (≥ 0.1.1): Fast haversine distance calculations for
  geographic filtering
- **leaderCluster** (≥ 1.5.0): Efficient spatial clustering for network
  simplification
- **mirai** (≥ 2.5.2): Asynchronous parallelism
- **kit** (≥ 0.0.21): Fast tabulation and vectorized operations

## Key Implementation Details

### Distance Matrix Strategy

The package uses adaptive distance matrix computation: - If network size
≤ `sqrt(dmat.max.size)`, precompute full distance matrix once -
Otherwise, compute in chunks as needed during OD-pair iteration -
Separate geodesic distance matrices are used for angle-based filtering

### Graph Representation

Graphs are represented as data frames with: - `from`, `to`: Node IDs
(integers) - `FX`, `FY`, `TX`, `TY`: Node coordinates (for spatial
operations) - `edge`: Edge identifier (optional, regenerated by many
functions) - Cost/attribute columns (e.g., `duration`, `cost`,
`distance`)

Internally, igraph is used for shortest path computation, but the
primary data structure is a data frame for flexibility and integration
with fastverse tools.

### Node Consolidation Algorithm

The
[`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md)
function uses a sophisticated multi-pass approach: 1. Drop loop edges,
duplicates, and singleton edges (optional) 2. Identify degree-2 nodes
(or nodes with deg_from=1 and deg_to=1 for directed graphs) 3. For
undirected graphs, orient edges so intermediate nodes appear as “from”
in one edge and “to” in another 4. Merge edges through intermediate
nodes, tracking groups via `gid` vector 5. Aggregate edge attributes
using collapse::collap() 6. Repeat recursively if `recursive = "full"`
until no more consolidation possible

The `by` parameter allows preserving mode/type boundaries by preventing
consolidation across different link characteristics.

## Working with Spatial Data

The package integrates with sf for spatial operations:

``` r
# Typical pattern for mapping OD zones to network nodes
nodes <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(od_zones, nodes)]
```

Coordinate columns (FX, FY, TX, TY) are preserved through most
operations and can be used to convert back to sf LINESTRING objects with
[`linestrings_from_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_from_graph.md).

## Performance Considerations

- Use `method = "AoN"` for large networks when route alternatives are
  not needed (much faster)
- Adjust `detour.max` and `angle.max` to control PSL computation time
  (lower values = faster)
- Set `unique.cost = TRUE` to deduplicate routes with same total cost
- Use `dmat.max.size` to control memory usage for large networks
- Enable multithreading with `nthreads` for large OD matrices
- Consider consolidating and simplifying networks before assignment to
  reduce computational burden

## Package Structure

This is a standard R package with: - DESCRIPTION: Package metadata and
dependencies - NAMESPACE: Exported functions (managed by roxygen2) - R/:
R source code - src/: C source code with compiled .so/.dll - man/:
Documentation (auto-generated from roxygen2 comments) - tests/testthat/:
Test suite - vignettes/: Package vignettes - data/: Included datasets
(Africa network, cities, trade)

Use roxygen2 for documentation - add `#'` comments above functions and
run `devtools::document()` to update NAMESPACE and man/ files.
