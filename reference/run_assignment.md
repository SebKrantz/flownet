# Run Traffic Assignment

Assign traffic flows to network edges using either Path-Sized Logit
(PSL) or All-or-Nothing (AoN) assignment methods.

## Usage

``` r
run_assignment(
  graph_df,
  od_matrix_long,
  directed = FALSE,
  cost.column = "cost",
  method = c("PSL", "AoN"),
  beta = 1,
  ...,
  detour.max = 1.5,
  angle.max = 90,
  unique.cost = TRUE,
  npaths.max = Inf,
  dmat.max.size = 10000^2,
  return.extra = NULL,
  verbose = TRUE,
  nthreads = 1L
)

# S3 method for class 'flownet'
print(x, ...)
```

## Arguments

- graph_df:

  A data.frame with columns `from`, `to`, and optionally a cost column.
  Represents the network graph with edges between nodes.

- od_matrix_long:

  A data.frame with columns `from`, `to`, and `flow`. Represents the
  origin-destination matrix in long format with flow values.

- directed:

  Logical (default: FALSE). Whether the graph is directed.

- cost.column:

  Character string (default: "cost") or numeric vector. Name of the cost
  column in `graph_df`, or a numeric vector of edge costs with length
  equal to `nrow(graph_df)`. The cost values are used to compute
  shortest paths and determine route probabilities.

- method:

  Character string (default: "PSL"). Assignment method:

  - `"PSL"`: Path-Sized Logit model considering multiple routes with
    overlap correction

  - `"AoN"`: All-or-Nothing assignment, assigns all flow to the shortest
    path (faster but no route alternatives)

- beta:

  Numeric (default: 1). Path-sized logit parameter (beta_PSL). Only used
  for PSL method.

- ...:

  Additional arguments (currently ignored).

- detour.max:

  Numeric (default: 1.5). Maximum detour factor for alternative routes
  (applied to shortest paths cost). Only used for PSL method. This is a
  key parameter controlling the execution time of the algorithm:
  considering more routes (higher `detour.max`) substantially increases
  computation time.

- angle.max:

  Numeric (default: 90). Maximum detour angle (in degrees, two sided).
  Only used for PSL method. I.e., nodes not within this angle measured
  against a straight line from origin to destination node will not be
  considered for detours.

- unique.cost:

  Logical (default: TRUE). Deduplicates paths based on the total cost
  prior to generating them. Only used for PSL method. Since multiple
  'intermediate nodes' may be on the same path, there is likely a
  significant number of duplicate paths which this option removes.

- npaths.max:

  Integer (default: Inf). Maximum number of paths to compute per
  OD-pair. Only used for PSL method. If the number of paths exceeds this
  number, a random sample will be taken from all but the shortest path.

- dmat.max.size:

  Integer (default: 1e4^2). Maximum size of distance matrices (both
  shortest paths and geodesic) to precompute. If smaller than
  `n_nodes^2`, then the full matrix is precomputed. Otherwise, it is
  computed in chunks as needed, where each chunk has `dmat.max.size`
  elements. Only used for PSL method.

- return.extra:

  Character vector specifying additional results to return. Options
  include: `"graph"`, `"paths"`, `"edges"` (PSL only), `"counts"`,
  `"costs"`, and `"weights"` (PSL only). For AoN: `"paths"` returns a
  list of shortest paths (one integer vector per OD pair), `"costs"`
  returns a numeric vector of shortest path costs, and `"counts"`
  returns a global integer vector of edge traversal counts. Use `"all"`
  to return all available extra results for the selected method.

- verbose:

  Logical (default: TRUE). Show progress bar and intermediate steps
  completion status?

- nthreads:

  Integer (default: 1L). Number of threads (daemons) to use for parallel
  processing with
  [`mirai`](https://mirai.r-lib.org/reference/mirai.html). Should not
  exceed the number of logical processors.

- x:

  An object of class `flownet`, typically returned by `run_assignment`.

## Value

A list of class `"flownet"` containing:

- `call` - The function call

- `final_flows` - Numeric vector of assigned flows for each edge (same
  length as `nrow(graph_df)`)

- `od_pairs_used` - Indices of OD pairs with valid flows

- Additional elements as specified in `return.extra`:

  - `graph` - The igraph graph object

  - `paths` - For PSL: list of lists (multiple routes per OD pair); for
    AoN: list of integer vectors (one shortest path per OD pair)

  - `edges` - List of edge indices used for each OD pair (PSL only)

  - `edge_counts` - For PSL: list of edge visit counts per OD pair; for
    AoN: integer vector of global edge traversal counts

  - `path_costs` - For PSL: list of path costs per OD pair; for AoN:
    numeric vector of shortest path costs

  - `path_weights` - List of path weights (probabilities) for each OD
    pair (PSL only)

## Details

This function performs traffic assignment using one of two methods:

**All-or-Nothing (AoN) Method:** A simple assignment method that assigns
all flow from each OD pair to the single shortest path. This is much
faster than PSL but does not consider route alternatives or overlaps.
Parameters `detour.max`, `angle.max`, `unique.cost`, `npaths.max`,
`beta`, and `dmat.max.size` are ignored for AoN.

**Path-Sized Logit (PSL) Method:** A more sophisticated assignment
method that considers multiple alternative routes:

- Creates a graph from `graph_df` using igraph, normalizing node IDs
  internally

- Computes shortest path distance matrix for all node pairs (chunkwise
  if `dmat.max.size < n_nodes^2`)

- For each origin-destination pair in `od_matrix_long`:

  - Identifies alternative routes (detours) that are within `detour.max`
    of shortest path cost. This is done by considering all other nodes
    and adding the shortest paths costs from origin node to an
    intermediate node (any other node) and from intermediate node to
    destination node. If `angle.max` is specified, filters detours to
    those within the specified angle from origin to destination. This
    means we only consider intermediate nodes that are roughly 'in the
    direction' of the destination node, and also not further away in
    terms of geodesic distance. If also `unique.cost = TRUE`, duplicate
    paths are removed based on the total cost. Thus, using only the
    shortest-paths-cost matrix and a matrix of the geodesic distances of
    the nodes (if `is.finite(angle.max)`), the algorithm pre-selects
    unique paths that are plausible both in terms of detour factor
    (cost) and direction before actually computing them. This speeds up
    route enumeration and PSL computations considerably.

  - Finds shortest paths from origin to intermediate nodes and from
    intermediate nodes to destination

  - Filters paths to remove those with duplicate edges, i.e., where the
    intermediate node is approached and departed from via the same
    edge(s).

  - Computes path-sized logit probabilities accounting for route overlap

  - Assigns flows to edges based on probabilities weighted by route
    overlap

- Returns results including final flows and optionally additional
  information

The path-sized logit model accounts for route overlap by adjusting
probabilities based on the number of alternative routes using each edge.
Flows are assigned proportionally to the computed probabilities. The
model uses the parameter `beta` to control the sensitivity to route
overlap.

When `angle.max` is specified and `graph_df` contains coordinate columns
(`FX`, `FY`, `TX`, `TY`), the function uses geographic distance
calculations to restrict detours to those within the specified angle,
improving computational efficiency and route realism.

## See also

[flownet-package](https://sebkrantz.github.io/flownet/reference/flownet-package.md)

## Examples

``` r
library(flownet)
library(sf)

# Load existing network edges (exclude proposed new links)
africa_net <- africa_network[!africa_network$add, ]

# Convert to graph (use atomic_elem to drop sf geometry, qDF for data.frame)
graph <- collapse::atomic_elem(africa_net) |> collapse::qDF()
nodes <- nodes_from_graph(graph, sf = TRUE)

# Map cities/ports to nearest network nodes
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

# Simple gravity-based OD matrix
od_mat <- outer(africa_cities_ports$population, africa_cities_ports$population) / 1e12
dimnames(od_mat) <- list(nearest_nodes, nearest_nodes)
od_matrix_long <- melt_od_matrix(od_mat)

if (FALSE) { # \dontrun{
# Run Traffic Assignment (All-or-Nothing method)
result <- run_assignment(graph, od_matrix_long, cost.column = "duration",
                         method = "AoN", return.extra = "all")
print(result)

# Visualize Results
africa_net$final_flows_log10 <- log10(result$final_flows + 1)
plot(africa_net["final_flows_log10"])
} # }
```
