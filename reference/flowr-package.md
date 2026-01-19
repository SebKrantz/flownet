# Transport Modeling

*flowr* provides efficient tools for transportation modeling,
specifically route enumeration and traffic assignment tasks. The package
implements the path-sized logit model for traffic assignment and
provides utilities for network processing.

**Network Processing**

[`linestrings_to_graph()`](https://sebkrantz.github.io/flowr/reference/linestrings_to_graph.md)
— Convert LINESTRING geometries to graph  
[`create_undirected_graph()`](https://sebkrantz.github.io/flowr/reference/create_undirected_graph.md)
— Convert directed graph to undirected  
[`consolidate_graph()`](https://sebkrantz.github.io/flowr/reference/consolidate_graph.md)
— Consolidate graph by removing intermediate nodes  
[`simplify_network()`](https://sebkrantz.github.io/flowr/reference/simplify_network.md)
— Simplify network graph  

**Traffic Assignment**

[`run_assignment()`](https://sebkrantz.github.io/flowr/reference/run_assignment.md)
— Run traffic assignment using path-sized logit model  

**Graph Utilities**

[`normalize_graph()`](https://sebkrantz.github.io/flowr/reference/normalize_graph.md)
— Normalize node IDs to consecutive integers  
[`nodes_from_graph()`](https://sebkrantz.github.io/flowr/reference/nodes_from_graph.md)
— Extract unique nodes from graph  
[`linestrings_from_graph()`](https://sebkrantz.github.io/flowr/reference/linestrings_from_graph.md)
— Convert graph to LINESTRING geometries  
[`distances_from_graph()`](https://sebkrantz.github.io/flowr/reference/distances_from_graph.md)
— Compute distance matrix from graph  

**OD Matrix Utilities**

[`melt_od_matrix()`](https://sebkrantz.github.io/flowr/reference/melt_od_matrix.md)
— Convert origin-destination matrix to long format  

**Data**

[`africa_trade`](https://sebkrantz.github.io/flowr/reference/africa_trade.md)
— Average BACI HS96 2012-22 trade flows by section between 47
continental African countries  
[`africa_cities_ports`](https://sebkrantz.github.io/flowr/reference/africa_cities_ports.md)
— The 453 largest (port-)cities in continental Africa within a 70km
radius - from Krantz (2024),
[doi:10.1596/1813-9450-10893](https://doi.org/10.1596/1813-9450-10893)  
[`africa_network`](https://sebkrantz.github.io/flowr/reference/africa_network.md)
— African continental road network + extensions to optimally connect the
453 cities - from Krantz (2024),
[doi:10.1596/1813-9450-10893](https://doi.org/10.1596/1813-9450-10893)  
[`africa_segments`](https://sebkrantz.github.io/flowr/reference/africa_segments.md)
— Primary segments derived from OpenStreetMap routes between the 453
cities - from Krantz (2024),
[doi:10.1596/1813-9450-10893](https://doi.org/10.1596/1813-9450-10893)  

Replication materials:
<https://github.com/SebKrantz/OptimalAfricanRoads>

## Details

The package uses efficient C implementations for critical path
operations and leverages:

- `collapse` - Fast data transformations

- `geodist` - Fast geodesic distance computations

- `igraph` - Graph operations and shortest path algorithms

- `leaderCluster` - Fast spatial clustering for network simplification

## Author

Sebastian Krantz <sebastian.krantz@graduateinstitute.ch> and Kamol Roy
<kamol.roy08@gmail.com>
