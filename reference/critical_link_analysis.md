# Critical Link Analysis

Compute edge-level detour costs and vulnerability ratios for a graph.

## Usage

``` r
critical_link_analysis(
  graph_df,
  directed = FALSE,
  cost.column = "cost",
  nthreads = 1L
)
```

## Arguments

- graph_df:

  A data.frame with columns `from` and `to`, and optionally a cost
  column. Each row represents one physical edge.

- directed:

  Logical (default: `FALSE`). Whether to treat the graph as directed.
  For undirected graphs, each edge can be traversed in both directions,
  but the returned result remains aligned with the original rows of
  `graph_df`.

- cost.column:

  Character string (default: `"cost"`) or numeric vector. Name of the
  cost column in `graph_df`, or a numeric vector of edge costs with
  length equal to `nrow(graph_df)`.

- nthreads:

  Integer (default: `1L`). Number of OpenMP threads used to parallelize
  the per-source Dijkstra computations. Has no effect when the package
  was compiled without OpenMP support.

## Value

A data.frame with one row per input edge and columns:

- `edge` - Edge identifier (index of the edge in `graph_df`)

- `from` - Original from-node ID

- `to` - Original to-node ID

- `edge_cost` - Cost of traversing the edge

- `detour_cost` - Shortest path cost between the edge endpoints after
  removing the edge

- `detour_ratio` - `detour_cost / edge_cost`

## Details

The detour cost for edge `e = (u, v)` is the shortest path cost from `u`
to `v` after removing that physical edge row. Parallel edges are kept as
distinct alternatives. If no detour exists, `detour_cost` and
`detour_ratio` are `Inf`.

The implementation uses a compiled multi-label Dijkstra routine, run
once per unique `from` node. For each target node, it tracks the two
best distances whose first physical edge differs, which yields the best
endpoint detour for every outgoing edge of that source without repeated
graph deletion. Per-source Dijkstras are independent and are
parallelized across `nthreads` OpenMP threads.

## See also

[distances_from_graph](https://sebkrantz.github.io/flownet/reference/distances_from_graph.md)
[run_assignment](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
[flownet-package](https://sebkrantz.github.io/flownet/reference/flownet-package.md)

## Examples

``` r
library(flownet)

graph <- data.frame(
  from = c(1, 2, 1),
  to = c(2, 3, 3),
  cost = c(1, 1, 3)
)

critical_link_analysis(graph, cost.column = "cost")
#>   edge from to edge_cost detour_cost detour_ratio
#> 1    1    1  2         1           4    4.0000000
#> 2    2    2  3         1           4    4.0000000
#> 3    3    1  3         3           2    0.6666667
```
