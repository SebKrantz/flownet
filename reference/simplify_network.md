# Simplify Network

Simplify a network graph using shortest paths or node clustering
methods.

## Usage

``` r
simplify_network(
  graph_df,
  nodes,
  method = c("shortest-paths", "cluster"),
  directed = FALSE,
  cost.column = "cost",
  by = NULL,
  radius_km = list(nodes = 7, cluster = 20),
  ...
)
```

## Arguments

- graph_df:

  A data.frame with columns `from` and `to` representing the graph
  edges. For the cluster method, the graph must also have columns `FX`,
  `FY`, `TX`, `TY` representing node coordinates.

- nodes:

  For `method = "shortest-paths"`: either an atomic vector of node IDs,
  or a data.frame with columns `from` and `to` specifying
  origin-destination pairs. For `method = "cluster"`: an atomic vector
  of node IDs to preserve. These nodes will be kept as cluster
  centroids, and nearby nodes (within `radius_km$nodes`) will be
  assigned to their clusters. Remaining nodes are clustered using
  [`leaderCluster`](https://rdrr.io/pkg/leaderCluster/man/leaderCluster.html).

- method:

  Character string (default: "shortest-paths"). Method to use for
  simplification: `"shortest-paths"` computes shortest paths between
  nodes and keeps only traversed edges; `"cluster"` clusters nodes using
  the
  [`leaderCluster`](https://rdrr.io/pkg/leaderCluster/man/leaderCluster.html)
  algorithm and contracts the graph.

- directed:

  Logical (default: FALSE). Whether the graph is directed. For
  `method = "shortest-paths"`: controls path computation direction. For
  `method = "cluster"`: if TRUE, A-\>B and B-\>A remain as separate
  edges after contraction; if FALSE, edges are normalized so that
  `from < to` before grouping.

- cost.column:

  Character string (default: "cost"). Name of the cost column in
  `graph_df`. Alternatively, a numeric vector of edge costs with length
  equal to `nrow(graph_df)`. With `method = "cluster"`, a numeric vector
  of node weights matching `nodes_from_graph(graph_df)` can be provided.

- by:

  Link characteristics to preserve/not simplify across, passed as a
  one-sided formula or character vector of column names. Typically
  includes attributes like *mode*, *type*, or *capacity*. For
  `method = "shortest-paths"`: paths are computed separately for each
  group defined by `by`, with edges not in the current group penalized
  (cost multiplied by 100) to compel mode-specific routes. For
  `method = "cluster"`: edges are grouped by `from`, `to`, AND `by`
  columns, preventing consolidation across different modes/types.

- radius_km:

  Named list with elements `nodes` (default: 7) and `cluster` (default:
  20). Only used for `method = "cluster"`. `nodes`: radius in kilometers
  around preserved nodes. Graph nodes within this radius will be
  assigned to the nearest preserved node's cluster. `cluster`: radius in
  kilometers for clustering remaining nodes using leaderCluster.

- ...:

  For `method = "cluster"`: additional arguments passed to
  [`collap`](https://fastverse.org/collapse/reference/collap.html) for
  edge attribute aggregation.

## Value

A data.frame containing the simplified graph with:

- For `method = "shortest-paths"`:

  - All columns from the input `graph_df` (for edges that were kept)

  - Attribute `"edges"`: integer vector of edge indices from the
    original graph

  - Attribute `"edge_counts"`: integer vector indicating how many times
    each edge was traversed

- For `method = "cluster"`:

  - `edge` - New edge identifier

  - `from`, `to` - Cluster centroid node IDs

  - `FX`, `FY`, `TX`, `TY` - Coordinates of cluster centroid nodes

  - Aggregated edge attributes from the original graph

  - Attribute `"group.id"`: mapping from original edges to simplified
    edges

  - Attribute `"group.starts"`: start indices of each group

  - Attribute `"group.sizes"`: number of original edges per simplified
    edge

## Details

**Method: "shortest-paths"**

- Validates that all origin and destination nodes exist in the network

- Computes shortest paths from each origin to all destinations using
  igraph

- Marks all edges that are traversed by at least one shortest path

- Returns only the subset of edges that were traversed

- If `nodes` is a data frame with `from` and `to` columns, paths are
  computed from each unique origin to its specified destinations

**Method: "cluster"**

- Requires the graph to have spatial coordinates (`FX`, `FY`, `TX`,
  `TY`)

- If `nodes` is provided, these nodes are preserved as cluster centroids

- Nearby nodes (within `radius_km$nodes` km) are assigned to the nearest
  preserved node

- Remaining nodes are clustered using
  [`leaderCluster`](https://rdrr.io/pkg/leaderCluster/man/leaderCluster.html)
  with `radius_km$cluster` as the clustering radius

- For each cluster, the node closest to the cluster centroid is selected
  as representative

- The graph is contracted by mapping all nodes to their cluster
  representatives

- Self-loops (edges where both endpoints map to the same cluster) are
  dropped

- For undirected graphs (`directed = FALSE`), edges are normalized so
  `from < to`, merging opposite-direction edges; for directed graphs,
  A-\>B and B-\>A remain separate

- Edge attributes are aggregated using
  [`collap`](https://fastverse.org/collapse/reference/collap.html)
  (default: mean for numeric, mode for categorical); customize via `...`

## Examples

``` r
library(flownet)
library(sf)

# Convert segments to undirected graph
graph <- africa_segments |>
  linestrings_from_graph() |>
  linestrings_to_graph() |>
  create_undirected_graph(FUN = "fsum")

# Get city/port nodes to preserve
nodes_df <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes_df$node[st_nearest_feature(africa_cities_ports, nodes_df)]

# Initial consolidation
graph <- consolidate_graph(graph, keep = nearest_nodes, w = ~ passes)
#> Consolidating undirected graph graph with 11385 edges using full recursion
#> Initial node degrees:
#>    1    2    3    4    5    6    7    8    9   10 
#>  125 5293 3240  395  102   31    4    2    1    1 
#> 
#> Dropped 44 loop edges
#> Dropped 11 edges leading to singleton nodes
#> Oriented 3431 undirected intermediate edges
#> Consolidated 5079 intermediate nodes
#> Oriented 10 undirected intermediate edges
#> Consolidated 10 intermediate nodes
#> Aggregated 11330 edges down to 6597 edges
#> Oriented 361 undirected intermediate edges
#> Consolidated 375 intermediate nodes
#> Oriented 8 undirected intermediate edges
#> Consolidated 8 intermediate nodes
#> Aggregated 6597 edges down to 6264 edges
#> Oriented 41 undirected intermediate edges
#> Consolidated 43 intermediate nodes
#> Oriented 1 undirected intermediate edges
#> Consolidated 1 intermediate nodes
#> Aggregated 6264 edges down to 6224 edges
#> Oriented 2 undirected intermediate edges
#> Consolidated 3 intermediate nodes
#> Aggregated 6224 edges down to 6221 edges
#> No nodes to consolidate, returning graph
#> 
#> Consolidated undirected graph graph from 11385 edges to 6221 edges (54.6%)
#> Final node degrees:
#>    1    2    3    4    5    6    7    8 
#>  114  236 3246  392   84   18    2    1 

# Method 1: Shortest-paths simplification (keeps only traversed edges)
graph_simple <- simplify_network(graph, nearest_nodes,
                                 method = "shortest-paths",
                                 cost.column = ".length")
nrow(graph_simple)  # Reduced number of edges
#> [1] 4422

# \donttest{
# Method 2: Cluster-based simplification (contracts graph spatially)
# Compute node weights for clustering
node_weights <- collapse::rowbind(
  collapse::fselect(graph, node = from, gravity_rd),
  collapse::fselect(graph, to, gravity_rd),
  use.names = FALSE) |>
  collapse::collap(~ node, "fsum")

graph_cluster <- simplify_network(graph, nearest_nodes,
                                  method = "cluster",
                                  cost.column = node_weights$gravity_rd,
                                  radius_km = list(nodes = 30, cluster = 27),
                                  w = ~ passes)
nrow(graph_cluster)
#> [1] 2430
# }
```
