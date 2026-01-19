# flowr

**Transport Modeling: Route Enumeration and Traffic Assignment with the
Path-Sized Logit**

***NOTE: Package is still under development***

`flowr` provides efficient tools for transportation modeling,
specifically route enumeration and traffic assignment tasks. The package
implements the path-sized logit (PSL) model for traffic assignment and
provides powerful utilities for network processing.

## Features

- **Path-Sized Logit Model**: Efficient traffic assignment accounting
  for route overlap
- **Network Processing**: Convert LINESTRING geometries to graphs,
  consolidate graphs, simplify networks, and handle directed/undirected
  graphs
- **Route Enumeration**: Efficient algorithm for finding alternative
  routes between origin-destination pairs
- **High Performance**: C implementations for critical path operations
- **Multithreading**: Asynchronous parallelism using `mirai` for faster
  processing of large networks

## Installation

``` r
# Install development version from GitHub
remotes::install_github("CPCS-IAU/flowr")

# Alternatively download the repo, extract it, and run
install.packages("path/to/flowr", repos = NULL, type = "source")
```

## Dependencies

- `collapse` (\>= 2.1.5) - Fast data transformations and memory
  efficient programming
- `kit` (\>= 0.0.5) - Fast tabulation and vectorized switches
- `igraph` (\>= 2.1.4) - Graph operations and shortest path algorithms
- `sf` (\>= 1.0.0) - Spatial data handling
- `geodist` (\>= 0.1.1) - Fast geodesic distance computations
- `leaderCluster` (\>= 0.3.0) - Fast spatial clustering for network
  simplification
- `mirai` (\>= 2.5.2) - Asynchronous parallelism for R
- `progress` (\>= 1.2.3) - Progress bars for long-running operations

## Quick Start

### Basic Usage

``` r
library(flowr)

# Create a small graph data frame
graph <- data.frame(from = c(1, 2, 2, 3),
                    to = c(2, 3, 4, 4), cost = c(5, 3, 2, 4))

# Prepare OD matrix with the same node IDs as in graph
od_matrix_long <- data.frame(from = c(1, 2, 3),
                             to = c(4, 4, 4), flow = c(100, 80, 60))

# Run traffic assignment (and route enumeration)
result <- run_assignment(graph, od_matrix_long, angle.max = NA)
#> Created graph with 4 nodes and 4 edges...
#> Computed distance matrix of dimensions 4 x 4 ...

# Access results
result$final_flows
#> [1] 100.00000  16.13649 196.13649  43.86351
```

^(Created on 2025-11-24 with [reprex v2.1.1](https://reprex.tidyverse.org))

### Working with Spatial Networks

``` r
library(flowr)
library(sf)

# Read network from shapefile and create undirected graph (optional)
graph <- st_read("data/network.shp") |> 
  linestrings_to_graph() |>
  create_undirected_graph()

# Read zone centroids and get nearest nodes
od_zones <- st_read("data/od_zones.shp") |> st_centroid()
nodes <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(od_zones, nodes)]

# Consolidate Graph (optional)
graph <- consolidate_graph(graph, keep = nearest_nodes, w = ~ cost)

# Simplify network by keeping only traversed edges along shortest paths (optional)
# Use 'by' for multimodal networks to compute paths separately per mode
graph <- simplify_network(graph, nearest_nodes, cost.column = "cost", by = ~ mode)
```

## Main Functions

### Traffic Assignment and Route Enumeration

- **[`run_assignment()`](https://sebkrantz.github.io/flowr/reference/run_assignment.md)**:
  - Iterates through OD-pairs generating sensible alternative routes
  - Assigns traffic flows to network edges using path-sized logit model
  - Supports directed and undirected graphs
  - Returns flows and optional path/route information
  - **Key Parameters**:
    - **`beta`** (default: 1): Path-sized logit parameter (beta_PSL)
    - **`detour.max`** (default: 1.5): Maximum detour factor for
      alternative routes. Higher values consider more routes but
      increase computation time
    - **`angle.max`** (default: 90): Maximum detour angle in degrees
      (two-sided)
    - **`return.extra`**: Additional results to return from the route
      enumeration stage (`"graph"`, `"dmat"`, `"paths"`, `"edges"`,
      `"counts"`, `"costs"`, `"weights"`, or `"all"`)

### Network Processing

- **[`linestrings_to_graph()`](https://sebkrantz.github.io/flowr/reference/linestrings_to_graph.md)** -
  Convert LINESTRING geometries to graph data frame
- **[`create_undirected_graph()`](https://sebkrantz.github.io/flowr/reference/create_undirected_graph.md)** -
  Convert directed graph to undirected with edge aggregation
- **[`consolidate_graph()`](https://sebkrantz.github.io/flowr/reference/consolidate_graph.md)** -
  Consolidate graph by removing intermediate nodes and merging edges
- **[`simplify_network()`](https://sebkrantz.github.io/flowr/reference/simplify_network.md)** -
  Simplify network using shortest-paths or spatial clustering methods.
  Supports multimodal networks via `by` argument

### Graph Utilities

- **[`nodes_from_graph()`](https://sebkrantz.github.io/flowr/reference/nodes_from_graph.md)** -
  Extract unique nodes with coordinates from graph
- **[`normalize_graph()`](https://sebkrantz.github.io/flowr/reference/normalize_graph.md)** -
  Normalize node IDs to consecutive integers starting from 1
- **[`linestrings_from_graph()`](https://sebkrantz.github.io/flowr/reference/linestrings_from_graph.md)** -
  Convert graph to LINESTRING geometries
- **[`distances_from_graph()`](https://sebkrantz.github.io/flowr/reference/distances_from_graph.md)** -
  Compute distance matrix for all node pairs

### OD Matrix Utilities

- **[`melt_od_matrix()`](https://sebkrantz.github.io/flowr/reference/melt_od_matrix.md)** -
  Convert origin-destination matrices to long format

## Example Workflow

``` r
library(fastverse)
fastverse_extend(flowr, sf, mapview)

# 1. Load network and OD zone nodes
network <- st_read("data/network.shp")
od_zones <- st_read("data/od_zones.shp") |> st_centroid()
od_matrix <- fread("data/od_container_flows.csv") |> qM(1)
if(!all(dim(od_matrix) == nrow(od_zones))) stop("zones and OD matrix must match")

# 2. Convert network to graph
graph <- network |>
  linestrings_to_graph() |>
  create_undirected_graph()

# 3. Map zones to nearest network nodes
nodes <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(od_zones, nodes)]

# 4. Prepare OD matrix
od_matrix_long <- melt_od_matrix(od_matrix, nodes = nearest_nodes)

# 5. Run assignment
result <- run_assignment(graph, od_matrix_long, cost.column = "cost_column")

# 6. Visualize results (optional)
network$final_flows <- NA_real_
network$final_flows[attr(graph, "group.starts")] <- result$final_flows
mapview(network, zcol = "final_flows")
```

## Example Data

The package includes four example datasets for Africa:

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

## Suggested Packages

- **`fastverse`** (\>= 0.3.4) - Enhanced data manipulation workflow
- **`mapview`** (\>= 2.11.2) - Interactive visualization of results
- **`testthat`** (\>= 3.0.0) - Unit testing framework

## Authors

- Sebastian Krantz (<sebastian.krantz@graduateinstitute.ch>)
- Kamol Roy (<kamol.roy08@gmail.com>)

## License

GPL-3

## Citation

If you use `flowr` in your research, please cite:

``` r
citation("flowr")
```
