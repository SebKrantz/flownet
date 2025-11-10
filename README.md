# flowr

**Transport Modeling: Route Enumeration and Traffic Assignment with the Path-Sized Logit**

`flowr` provides efficient tools for transportation modeling, specifically route enumeration and traffic assignment tasks. The package implements the path-sized logit (PSL) model for traffic assignment and provides powerful utilities for network processing.

## Features

- **Path-Sized Logit Model**: Efficient traffic assignment accounting for route overlap
- **Network Processing**: Convert LINESTRING geometries to graphs, simplify networks, and handle directed/undirected graphs
- **Route Enumeration**: Efficient algorithm for finding alternative routes between origin-destination pairs
- **High Performance**: C implementations for critical path operations

## Installation

```r
# Install from source
install.packages("flowr", repos = NULL, type = "source")
```

## Dependencies

- `collapse` (>= 2.0.0) - Fast data transformations
- `igraph` (>= 1.3.0) - Graph operations and shortest path algorithms
- `sf` (>= 1.0.0) - Spatial data handling

## Quick Start

### Basic Usage

```r
library(flowr)

# Convert LINESTRING network to graph
graph_df <- network_linestrings |>
  linestrings_to_graph() |>
  add_vars(cost = network_linestrings$cost) |>
  create_undirected_graph()

# Prepare OD matrix (origin-destination flows)
od_matrix_long <- data.frame(
  from = c(1, 2, 3),
  to = c(10, 11, 12),
  flow = c(100, 200, 150)
)

# Run traffic assignment
result <- run_assignment(
  graph_df = graph_df,
  od_matrix_long = od_matrix_long,
  method = "PSL",
  beta = -1,
  detour.max = 1.5
)

# Access results
final_flows <- result$final_flows
```

### Working with Spatial Networks

```r
# Read network from shapefile
network <- st_read("data/network/base_network.shp")

# Convert to graph
graph_df <- network |>
  linestrings_to_graph() |>
  add_vars(cost = st_length(network)) |>
  create_undirected_graph()

# Simplify network by keeping only traversed edges
simplified <- simplify_network(
  x = network,
  od_matrix_long = od_matrix_long,
  cost.column = NULL  # Uses st_length() by default
)
```

## Main Functions

### Traffic Assignment

- **`run_assignment()`** - Assign traffic flows to network edges using path-sized logit model
  - Supports directed and undirected graphs
  - Configurable detour factors and angle constraints
  - Returns flows and optional path information

### Network Processing

- **`linestrings_to_graph()`** - Convert LINESTRING geometries to graph data frame
- **`create_undirected_graph()`** - Convert directed graph to undirected with edge aggregation
- **`simplify_network()`** - Simplify network by keeping only edges traversed by shortest paths

### Graph Utilities

- **`nodes_from_graph()`** - Extract unique nodes with coordinates from graph
- **`dist_mat_from_graph()`** - Compute distance matrix for all node pairs

## Key Parameters

### `run_assignment()`

- **`beta`** (default: -1): Path-sized logit parameter (beta_PSL)
- **`detour.max`** (default: 1.5): Maximum detour factor for alternative routes. Higher values consider more routes but increase computation time
- **`angle.max`** (default: 90): Maximum detour angle in degrees (two-sided)
- **`return.extra`**: Additional results to return (`"graph"`, `"dmat"`, `"paths"`, `"edges"`, `"costs"`, `"weights"`, or `"all"`)

## Example Workflow

```r
library(flowr)
library(sf)

# 1. Load network and zone nodes
network <- st_read("data/network/base_network.shp")
zone_nodes <- st_read("data/zone_nodes/network_nodes.shp")

# 2. Convert network to graph
graph_df <- network |>
  linestrings_to_graph() |>
  add_vars(cost = network$cost) |>
  create_undirected_graph()

# 3. Map zones to nearest network nodes
nodes_df <- nodes_from_graph(graph_df)
nearest_nodes <- st_nearest_feature(
  zone_nodes,
  st_as_sf(nodes_df, coords = c("X", "Y"), crs = 4326)
)

# 4. Prepare OD matrix
od_matrix_long <- data.frame(
  from = rep.int(nearest_nodes, ncol(od_matrix)),
  to = rep(nearest_nodes, each = nrow(od_matrix)),
  flow = vec(od_matrix)
) |> fsubset(is.finite(flow) & flow > 0)

# 5. Run assignment
result <- run_assignment(
  graph_df = graph_df,
  od_matrix_long = od_matrix_long,
  return.extra = "all"
)

# 6. Visualize results
network$final_flows <- NA_real_
network$final_flows[attr(graph_df, "group.starts")] <- result$final_flows
```

## Authors

- Sebastian Krantz (sebastian.krantz@graduateinstitute.ch)
- Kamol Roy (kamol.roy08@gmail.com)

## License

GPL-3

## Citation

If you use `flowr` in your research, please cite:

```r
citation("flowr")
```

