# Melt Origin-Destination Matrix to Long Format

Convert an origin-destination (OD) matrix to a long-format data frame
with columns `from`, `to`, and `flow`.

## Usage

``` r
melt_od_matrix(od_matrix, nodes = NULL, sort = TRUE)
```

## Arguments

- od_matrix:

  A numeric matrix with origin-destination flows. Rows represent
  origins, columns represent destinations. The matrix should be square
  (same number of rows and columns).

- nodes:

  (Optional) Numeric or integer vector of node IDs in the graph matching
  the rows and columns of the matrix. If provided, must have length
  equal to both `nrow(od_matrix)` and `ncol(od_matrix)`. When `nodes` is
  provided, these IDs are used directly, ignoring row and column names.
  This is particularly useful when mapping zone-based OD matrices to
  graph node IDs (e.g., using
  [`st_nearest_feature`](https://r-spatial.github.io/sf/reference/st_nearest_feature.html)
  to find nearest graph nodes). If omitted, row and column names (if
  present) will be used as node IDs, coerced to integer if possible. If
  names are not available or cannot be coerced to integer, sequential
  integers will be used.

- sort:

  Sort long OD-matrix in ascending order of from and to columns. This
  can have computational benefits, e.g., when multithreading with
  `method = "AoN"`.

## Value

A data frame with columns:

- `from` - Origin node ID (integer)

- `to` - Destination node ID (integer)

- `flow` - Flow value (numeric)

Only rows with finite, positive flow values are included.

## Details

This function converts a square OD matrix to long format, which is
required by
[`run_assignment()`](https://sebkrantz.github.io/flowr/reference/run_assignment.md).
The behavior depends on whether `nodes` is provided:

**When `nodes` is provided:**

- The `nodes` vector is used directly as node IDs for both origins and
  destinations

- Row and column names are ignored (but must match if both are present)

- This is the recommended approach when working with zone-based OD
  matrices that need to be mapped to graph nodes, as it ensures the node
  IDs match those in the graph

**When `nodes` is omitted:**

- Row and column names are extracted from the matrix (if available)

- Names are coerced to integer if possible; if coercion fails or names
  are missing, sequential integers are used

- This approach works well when the matrix row/column names already
  correspond to graph node IDs

In both cases, the function:

- Creates a long-format data frame with all origin-destination pairs

- Filters out non-finite and zero flow values

The function is useful for converting OD matrices to the long format
required by
[`run_assignment()`](https://sebkrantz.github.io/flowr/reference/run_assignment.md).

## See also

[`africa_cities_ports`](https://sebkrantz.github.io/flowr/reference/africa_cities_ports.md),
[`africa_network`](https://sebkrantz.github.io/flowr/reference/africa_network.md),
[`nodes_from_graph()`](https://sebkrantz.github.io/flowr/reference/nodes_from_graph.md),
[`run_assignment()`](https://sebkrantz.github.io/flowr/reference/run_assignment.md),
[flowr-package](https://sebkrantz.github.io/flowr/reference/flowr-package.md)

## Examples

``` r
library(flowr)
library(sf)

# Load existing network and convert to graph
africa_net <- africa_network[!africa_network$add, ]
graph <- linestrings_to_graph(africa_net)
nodes <- nodes_from_graph(graph, sf = TRUE)

# Map cities/ports to nearest network nodes
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

# Example 1: Simple gravity-based OD matrix
od_mat <- outer(africa_cities_ports$population, africa_cities_ports$population) / 1e12
dimnames(od_mat) <- list(nearest_nodes, nearest_nodes)
od_long <- melt_od_matrix(od_mat)
head(od_long)
#>   from to       flow
#> 1    1  1 10.3925963
#> 2    1  6  0.3804031
#> 3    1  8  4.6445025
#> 4    1 11  3.1488124
#> 5    1 13  0.6058467
#> 6    1 15  0.4982475

# Example 2: Using nodes argument (when matrix has zone IDs, not node IDs)
# Here zones are 1:n_cities, nodes argument maps them to graph nodes
dimnames(od_mat) <- NULL
od_long2 <- melt_od_matrix(od_mat, nodes = nearest_nodes)
head(od_long2)
#>   from to       flow
#> 1    1  1 10.3925963
#> 2    1  6  0.3804031
#> 3    1  8  4.6445025
#> 4    1 11  3.1488124
#> 5    1 13  0.6058467
#> 6    1 15  0.4982475
```
