library(fastverse)
fastverse_extend(flownet, sf, mapview)

# --------------------------------------
# Basic Africa Example
# --------------------------------------

# Load Network (Keep only existing)
africa_net <- fsubset(africa_network, !add, -add)
mapview(africa_net, zcol = "speed_kmh", layer.name = "Speed in Kmh")

# Load Graph
graph <- atomic_elem(africa_net) |> qDF()
nodes <- nodes_from_graph(graph, sf = TRUE)

# Test Distance Matrix
anyNA(distances_from_graph(graph, cost.column = "duration"))

# Map Zones to Nodes
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

# Trade Flow Disaggregation (OD Matrix Generation)
city_pop <- africa_cities_ports |> atomic_elem() |> qDF() |>
  fcompute(node = nearest_nodes, city = qF(city_country),
           pop_share = fsum(population, iso3, TRA = "/"),
           keep = "iso3")

od_matrix_long <- africa_trade |>
  collap(quantity ~ iso3_o + iso3_d, fsum) |>
  join(city_pop |> add_stub("_o", FALSE), multiple = TRUE) |>
  join(city_pop |> add_stub("_d", FALSE), multiple = TRUE) |>
  fmutate(flow = quantity * pop_share_o * pop_share_d) |>
  frename(from = node_o, to = node_d) |>
  fsubset(is.finite(flow) & flow > 0)

# Test Trade Flow Disaggregation
all.equal(sum(od_matrix_long$flow), sum(africa_trade$quantity))

# Run Traffic Assignment (using subset for testing - remove [1:5000,] for full analysis)
result <- run_assignment(graph, od_matrix_long[1:5000,],
                         cost.column = "duration", method = "PSL")
print(result)

# Visualize Results
africa_net$final_flows_log10 <- log10(result$final_flows+1)

mapview(africa_net, zcol = "final_flows_log10",
        layer.name = "Trade Flows in Tons (Log10)") +
mapview(africa_cities_ports, alpha = 0, cex = 2,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "African (Port-)Cities")


# --------------------------------------
# Africa Network Consolidation Example
# --------------------------------------

# Load Segments
graph <- africa_segments |>
  linestrings_from_graph() |>
  linestrings_to_graph() |>
  create_undirected_graph(FUN = "fsum")
nodes <- nodes_from_graph(graph, sf = TRUE)

# Map Zones to Nodes
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

# Graph Consolidation
graph <- consolidate_graph(graph, keep = nearest_nodes, w = ~ passes)

# Check
mapview(linestrings_from_graph(graph)) + mapview(africa_cities_ports, col.regions = "red", alpha.regions = 0.8)

# Now Graph Simplification Via Spatial Clustering
node_weights <- rowbind(slt(graph, node = from, gravity_rd), slt(graph, to, gravity_rd),
                        use.names = FALSE) |> collap( ~ node, fsum)
graph <- simplify_network(graph, nearest_nodes, method = "cluster",
                          cost.column = node_weights$gravity_rd,
                          radius_km = list(nodes = 30, cluster = 27),
                          w = ~ passes)
# Final Graph Consolidation
graph <- consolidate_graph(graph, keep = nearest_nodes, w = ~ passes)
# Test Distance Matrix
anyNA(distances_from_graph(graph, cost.column = "duration"))

# Examine consolidation result
mapview(linestrings_from_graph(africa_segments), color = "dodgerblue4",
        layer.name = "Original") +
mapview(linestrings_from_graph(graph), color = "orange",
        layer.name = "Consolidated") +
mapview(africa_cities_ports, alpha = 0, cex = 2,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "African (Port-)Cities")

# Map Zones to Nodes
nodes <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

# Trade Flow Disaggregation (OD Matrix Generation)
city_pop <- africa_cities_ports |> atomic_elem() |> qDF() |>
  fcompute(node = nearest_nodes, city = qF(city_country),
           pop_share = fsum(population, iso3, TRA = "/"),
           keep = "iso3")

od_matrix_long <- africa_trade |>
  collap(quantity ~ iso3_o + iso3_d, fsum) |>
  join(city_pop |> add_stub("_o", FALSE), multiple = TRUE) |>
  join(city_pop |> add_stub("_d", FALSE), multiple = TRUE) |>
  fmutate(flow = quantity * pop_share_o * pop_share_d) |>
  frename(from = node_o, to = node_d) |>
  fsubset(is.finite(flow) & flow > 0)

# Test Trade Flow Disaggregation
all.equal(sum(od_matrix_long$flow), sum(africa_trade$quantity))

# Run Traffic Assignment
result <- run_assignment(graph, od_matrix_long,
                         cost.column = ".length", method = "PSL", nthreads = 8)
print(result)

# Visualize Results
graph$final_flows_log10 <- log10(result$final_flows+1)

mapview(linestrings_from_graph(graph), zcol = "final_flows_log10",
        layer.name = "Trade Flows in Tons (Log10)") +
  mapview(africa_cities_ports, alpha = 0, cex = 1.5,
          col.regions = "red", alpha.regions = 0.8,
          layer.name = "African (Port-)Cities")
