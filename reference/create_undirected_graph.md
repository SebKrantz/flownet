# Create Undirected Graph

Convert a directed graph to an undirected graph by normalizing edges and
aggregating duplicate connections.

## Usage

``` r
create_undirected_graph(graph_df, by = NULL, ...)
```

## Arguments

- graph_df:

  A data frame representing a directed graph including columns: `from`,
  `to`, and (optionally) `edge`, `FX`, `FY`, `TX`, `TY`.

- by:

  Link characteristics to preserve/not aggregate across, passed as a
  one-sided formula or character vector of column names. Typically this
  includes attributes like *mode*, *type*, or *capacity* to ensure that
  only edges with the same characteristics are aggregated.

- ...:

  Arguments passed to
  [`collap()`](https://fastverse.org/collapse/reference/collap.html) for
  aggregation across duplicated (diretional) edges. The defaults are
  `FUN = fmean` for numeric columns and `catFUN = fmode` for categorical
  columns. Select columns using `cols` or use argument
  `custom = list(fmean = cols1, fsum = cols2, fmode = cols3)` to map
  different columns to specific aggregation functions. You can weight
  the aggregation (using `w = ~ weight_col`).

## Value

A data frame representing an undirected graph with:

- `edge` - Edge identifier (first value from duplicates)

- `from` - Starting node ID (normalized to be \< `to`)

- `to` - Ending node ID (normalized to be \> `from`)

- `FX` - Starting node X-coordinate (first value from duplicates)

- `FY` - Starting node Y-coordinate (first value from duplicates)

- `TX` - Ending node X-coordinate (first value from duplicates)

- `TY` - Ending node Y-coordinate (first value from duplicates)

- Aggregated columns

## Details

This function converts a directed graph to an undirected graph by:

- Normalizing edge directions so that `from < to` for all edges

- Collapsing duplicate edges (same `from` and `to` nodes)

- For spatial/identifier columns (`edge`, `FX`, `FY`, `TX`, `TY`),
  taking the first value from duplicates

- For aggregation columns,
  [`collap()`](https://fastverse.org/collapse/reference/collap.html)
  will be applied.

## Examples

``` r
library(flowr)

# Convert segments to graph and make undirected
graph <- africa_segments |>
  linestrings_from_graph() |>
  linestrings_to_graph()
graph_undir <- create_undirected_graph(graph, FUN = "fsum")

# Fewer edges after removing directional duplicates
c(directed = nrow(graph), undirected = nrow(graph_undir))
#>   directed undirected 
#>      14358      11385 
```
