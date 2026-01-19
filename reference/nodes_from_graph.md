# Extract Nodes from Graph

Extract unique nodes with their coordinates from a graph data frame.

## Usage

``` r
nodes_from_graph(graph_df, sf = FALSE, crs = 4326)
```

## Arguments

- graph_df:

  A data frame representing a graph with columns: `from`, `to`, `FX`,
  `FY`, `TX`, `TY`.

- sf:

  Logical. If TRUE, returns result as an `sf` POINT object. Default:
  FALSE.

- crs:

  Coordinate reference system for sf output; default is 4326.

## Value

A data frame (or sf object if `sf = TRUE`) with unique nodes and
coordinates:

- `node` - Node ID

- `X` - Node X-coordinate (typically longitude)

- `Y` - Node Y-coordinate (typically latitude)

Result is sorted by node ID.

## Details

This function extracts all unique nodes from both the `from` and `to`
columns of the graph, along with their corresponding coordinates.
Duplicate nodes are removed, keeping only unique node IDs with their
coordinates.

## Examples

``` r
library(flowr)
library(sf)

# Load existing network edges and convert to graph
africa_net <- africa_network[!africa_network$add, ]
graph <- linestrings_to_graph(africa_net)

# Extract nodes from graph
nodes <- nodes_from_graph(graph)
head(nodes)
#>   node         X        Y
#> 1    1 -17.44671 14.69281
#> 2    2 -17.04453 14.70297
#> 3    3 -16.63849 15.11222
#> 4    4 -16.46332 14.75651
#> 5    5 -16.42054 14.34236
#> 6    6 -17.03333 20.93330

# Get nodes as sf POINT object for spatial operations
nodes_sf <- nodes_from_graph(graph, sf = TRUE)
class(nodes_sf)
#> [1] "sf"         "data.frame"

# Find nearest network nodes to cities/ports
nearest_nodes <- nodes_sf$node[st_nearest_feature(africa_cities_ports, nodes_sf)]
head(nearest_nodes)
#> [1]  903  275  609  553  880 1372
```
