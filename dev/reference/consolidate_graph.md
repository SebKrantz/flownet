# Consolidate Graph

Consolidate a graph by contracting/removing intermediate nodes (nodes
that occur exactly twice) and dropping loop, duplicate, and singleton
edges (leading to dead ends). This simplifies the network topology while
preserving connectivity.

## Usage

``` r
consolidate_graph(
  graph_df,
  directed = FALSE,
  drop.edges = c("loop", "duplicate", "single"),
  contract = TRUE,
  by = NULL,
  keep.nodes = NULL,
  ...,
  recursive = "full",
  verbose = TRUE
)
```

## Arguments

- graph_df:

  A data frame representing a graph with columns: `from` and `to` (node
  IDs), and optionally other columns to preserve. If coordinate columns
  (`FX`, `FY`, `TX`, `TY`) are present, they will be preserved and
  updated based on the consolidated node coordinates.

- directed:

  Logical (default: FALSE). Whether the graph is directed.

- drop.edges:

  Character vector (default: `c("loop", "duplicate", "single")`). Types
  of edges to drop: `"loop"` removes self-loops (edges where from ==
  to), `"duplicate"` removes duplicate edges (same from-to pair),
  `"single"` removes edges leading to singleton nodes (nodes that occur
  only once). Set to `NULL` to keep all edges.

- contract:

  Logical (default: TRUE). If TRUE, contracts the graph by removing
  intermediate nodes (nodes that occur exactly twice) and merging
  connecting edges. If FALSE, only drops edges as specified in
  `drop.edges`.

- by:

  Link characteristics to preserve/not contract across, passed as a
  one-sided formula or character vector of column names. Typically this
  includes attributes like *mode*, *type*, or *capacity* to ensure that
  only edges with the same characteristics are contracted.

- keep.nodes:

  Numeric vector (optional). Node IDs to preserve during consolidation,
  even if they occur exactly twice. Also used to preserve nodes when
  dropping singleton edges.

- ...:

  Arguments passed to
  [`collap()`](https://fastverse.org/collapse/reference/collap.html) for
  aggregation across contracted edges. The defaults are `FUN = fmean`
  for numeric columns and `catFUN = fmode` for categorical columns.
  Select columns using `cols` or use argument
  `custom = list(fmean = cols1, fsum = cols2, fmode = cols3)` to map
  different columns to specific aggregation functions. It is highly
  recommended to weight the aggregation (using `w = ~ weight_col`) by
  the length/cost of the edges.

- recursive:

  One of `"none"/FALSE` (drop edges, contract, and aggregate once),
  `"partial"` (recursively drop edges and contract but only aggregate
  once), or `"full"/TRUE` (recursively drop edges, contract, and
  aggregate the graph until no further consolidation is possible).

- verbose:

  Logical (default: TRUE). Whether to print messages about dropped edges
  and consolidation progress.

## Value

A data frame representing the consolidated graph with:

- All columns from `graph_df` (aggregated if consolidation occurred),
  excluding `from`, `to`, and optionally `edge` and `FX`, `FY`, `TX`,
  `TY` (which are re-added if present in original)

- `from`, `to`, `edge` - Node/edge IDs (updated after consolidation)

- Coordinate columns (`FX`, `FY`, `TX`, `TY`) if present in original

- Attribute `"keep.edges"` - Indices of original edges that were kept
  (before aggregation)

- Attribute `"group.id"` - Integer mapping each kept edge to its row in
  the result (after aggregation)

## Details

This function consolidates/simplifies a graph by:

- **Dropping edges**: Optionally removes self-loops, duplicate edges,
  and edges leading to singleton nodes (nodes that appear only once in
  the graph)

- **Contracting nodes**: Removes intermediate nodes (nodes that occur
  exactly twice) by merging the two edges connected through them into a
  single longer edge

- **Aggregating attributes**: When edges are merged, attributes/columns
  are aggregated using
  [`collap()`](https://fastverse.org/collapse/reference/collap.html).
  The default aggregation is mean for numeric columns and mode for
  categorical columns.

- **Recursive consolidation**: If `recursive = TRUE`, the function
  continues consolidating until no more nodes can be dropped or
  contracted, ensuring complete consolidation

Consolidation is useful for simplifying network topology while
preserving connectivity. For example, if node B connects A-\>B and
B-\>C, it will be removed and replaced with A-\>C. With
`recursive = TRUE`, long chains (A-\>B-\>C-\>D) are fully contracted to
A-\>D in a single call.

For undirected graphs, the algorithm also handles cases where a node
appears twice as either origin or destination (circular cases).

If coordinate columns (`FX`, `FY`, `TX`, `TY`) are present in the input,
they are preserved and updated based on the consolidated node
coordinates from the original graph.

## See also

[create_undirected_graph](https://sebkrantz.github.io/flownet/dev/reference/create_undirected_graph.md)
[simplify_network](https://sebkrantz.github.io/flownet/dev/reference/simplify_network.md)
[flownet-package](https://sebkrantz.github.io/flownet/dev/reference/flownet-package.md)

## Examples

``` r
library(flownet)
library(sf)

# Convert segments to undirected graph
graph <- africa_segments |>
  linestrings_from_graph() |>
  linestrings_to_graph() |>
  create_undirected_graph(FUN = "fsum")

# Get nodes to preserve (city/port locations)
nodes <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

# Consolidate graph, preserving city nodes
graph_cons <- consolidate_graph(graph, keep = nearest_nodes, w = ~ passes)
#> Consolidating undirected graph graph with 11385 edges using full recursion
#> Initial node degrees:
#>    1    2    3    4    5    6    7    8    9   10 
#>  125 5293 3240  395  102   31    4    2    1    1 
#> 
#> Dropped 44 loop edges
#> Dropped 11 edges leading to singleton nodes
#> Oriented 3431 undirected intermediate edges
#> Contracted 5079 intermediate nodes
#> Oriented 10 undirected intermediate edges
#> Contracted 10 intermediate nodes
#> Aggregated 11330 edges down to 6597 edges
#> Oriented 361 undirected intermediate edges
#> Contracted 375 intermediate nodes
#> Oriented 8 undirected intermediate edges
#> Contracted 8 intermediate nodes
#> Aggregated 6597 edges down to 6264 edges
#> Oriented 41 undirected intermediate edges
#> Contracted 43 intermediate nodes
#> Oriented 1 undirected intermediate edges
#> Contracted 1 intermediate nodes
#> Aggregated 6264 edges down to 6224 edges
#> Oriented 2 undirected intermediate edges
#> Contracted 3 intermediate nodes
#> Aggregated 6224 edges down to 6221 edges
#> No nodes to contract, returning graph
#> 
#> Consolidated undirected graph graph from 11385 edges to 6221 edges (54.6%)
#> Final node degrees:
#>    1    2    3    4    5    6    7    8 
#>  114  236 3246  392   84   18    2    1 

# Consolidated graph has fewer edges
c(original = nrow(graph), consolidated = nrow(graph_cons))
#>     original consolidated 
#>        11385         6221 
```
