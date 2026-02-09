# Package index

- [`flownet-package`](https://sebkrantz.github.io/flownet/reference/flownet-package.md)
  [`flownet`](https://sebkrantz.github.io/flownet/reference/flownet-package.md)
  : Efficient Transport Modeling

## Network Processing

Utility functions to prepare networks graphs (represented using data
frames) for analysis.

- [`linestrings_to_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_to_graph.md)
  : Convert Linestring to Graph
- [`create_undirected_graph()`](https://sebkrantz.github.io/flownet/reference/create_undirected_graph.md)
  : Create Undirected Graph
- [`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md)
  : Consolidate Graph
- [`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)
  : Simplify Network

## Route Choice and Traffic Assignment

Run traffic assignment using path-sized logit model, including efficient
built-in route choice algorithm, or all-or-nothing assignment to the
shortest path.

- [`run_assignment()`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
  [`print(`*`<flownet>`*`)`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
  : Run Traffic Assignment

## Graph Utilities

Utility functions to normalize graphs (represented using data frames)
and extract information from them.

- [`normalize_graph()`](https://sebkrantz.github.io/flownet/reference/normalize_graph.md)
  : Normalize Graph Node IDs
- [`nodes_from_graph()`](https://sebkrantz.github.io/flownet/reference/nodes_from_graph.md)
  : Extract Nodes from Graph
- [`linestrings_from_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_from_graph.md)
  : Convert Graph to Linestrings
- [`distances_from_graph()`](https://sebkrantz.github.io/flownet/reference/distances_from_graph.md)
  : Compute Distance Matrix from Graph

## OD-Matrix Utilities

Utility functions to reshape origin-destination flow/trip matrices.

- [`melt_od_matrix()`](https://sebkrantz.github.io/flownet/reference/melt_od_matrix.md)
  : Melt Origin-Destination Matrix to Long Format

## Data

Trans-African Road Network and Trade Data

- [`africa_cities_ports`](https://sebkrantz.github.io/flownet/reference/africa_cities_ports.md)
  : African Cities and International Ports
- [`africa_network`](https://sebkrantz.github.io/flownet/reference/africa_network.md)
  : Trans-African Road Transport Network
- [`africa_segments`](https://sebkrantz.github.io/flownet/reference/africa_segments.md)
  : Raw Network Segments for Trans-African Transport Network
- [`africa_trade`](https://sebkrantz.github.io/flownet/reference/africa_trade.md)
  : Intra-African Trade Flows by HS Section
