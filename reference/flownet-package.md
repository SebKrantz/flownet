# Efficient Transport Modeling

*flownet* provides efficient tools for transportation modeling in R,
supporting network processing, route enumeration, and traffic assignment
tasks. It implements the path-sized logit (PSL) model for traffic
assignment and provides powerful utilities for network
processing/preparation.

**Network Processing**

[`linestrings_to_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_to_graph.md)
— Convert LINESTRING geometries to graph  
[`create_undirected_graph()`](https://sebkrantz.github.io/flownet/reference/create_undirected_graph.md)
— Convert directed graph to undirected  
[`consolidate_graph()`](https://sebkrantz.github.io/flownet/reference/consolidate_graph.md)
— Consolidate graph by removing intermediate nodes  
[`simplify_network()`](https://sebkrantz.github.io/flownet/reference/simplify_network.md)
— Simplify network graph  

**Traffic Assignment**

[`run_assignment()`](https://sebkrantz.github.io/flownet/reference/run_assignment.md)
— Run traffic assignment using path-sized logit model  

**Graph Utilities**

[`normalize_graph()`](https://sebkrantz.github.io/flownet/reference/normalize_graph.md)
— Normalize node IDs to consecutive integers  
[`nodes_from_graph()`](https://sebkrantz.github.io/flownet/reference/nodes_from_graph.md)
— Extract unique nodes from graph  
[`linestrings_from_graph()`](https://sebkrantz.github.io/flownet/reference/linestrings_from_graph.md)
— Convert graph to LINESTRING geometries  
[`distances_from_graph()`](https://sebkrantz.github.io/flownet/reference/distances_from_graph.md)
— Compute distance matrix from graph  

**OD Matrix Utilities**

[`melt_od_matrix()`](https://sebkrantz.github.io/flownet/reference/melt_od_matrix.md)
— Convert origin-destination matrix to long format  

**Data**

[`africa_trade`](https://sebkrantz.github.io/flownet/reference/africa_trade.md)
— Average BACI HS96 2012-22 trade flows by section between 47
continental African countries  
[`africa_cities_ports`](https://sebkrantz.github.io/flownet/reference/africa_cities_ports.md)
— The 453 largest (port-)cities in continental Africa within a 70km
radius - from Krantz (2024),
[doi:10.1596/1813-9450-10893](https://doi.org/10.1596/1813-9450-10893)  
[`africa_network`](https://sebkrantz.github.io/flownet/reference/africa_network.md)
— African continental road network + extensions to optimally connect the
453 cities - from Krantz (2024),
[doi:10.1596/1813-9450-10893](https://doi.org/10.1596/1813-9450-10893)  
[`africa_segments`](https://sebkrantz.github.io/flownet/reference/africa_segments.md)
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

## References

Ben-Akiva, M., & Bierlaire, M. (1999). Discrete choice methods and their
applications to short term travel decisions. In R. W. Hall (Ed.),
\*Handbook of Transportation Science\* (pp. 5–33). Springer US.
https://doi.org/10.1007/978-1-4615-5203-1_2

Cascetta, E. (2001). \*Transportation systems engineering: Theory and
methods\*. Springer.

Ben-Akiva, M., & Lerman, S. R. (1985). \*Discrete choice analysis:
Theory and application to travel demand\*. The MIT Press.

Ramming, M. S. (2002). \*Network knowledge and route choice\* (Doctoral
dissertation). Massachusetts Institute of Technology.

Prato, C. G. (2009). Route choice modeling: Past, present and future
research directions. \*Journal of Choice Modelling, 2\*(1), 65–100.
https://doi.org/10.1016/S1755-5345(13)70005-8

AequilibiaE Python Documentation:
https://www.aequilibrae.com/develop/python/route_choice/path_size_logit.html

## Author

Sebastian Krantz <sebastian.krantz@graduateinstitute.ch> and Kamol Roy
<kamol.roy08@gmail.com>
