library(fastverse)
fastverse_extend(flowr, sf, mapview)

# --------------------------------------
# Basic GCC Example
# --------------------------------------

mapview(network_gcc, zcol = "mode", layer.name = "Network Mode")

# Load Graph
graph <- linestrings_to_graph(network_gcc)
nodes <- nodes_from_graph(graph, sf = TRUE)

# Map Zones to Nodes
nearest_nodes <- nodes$node[st_nearest_feature(zones_gcc, nodes)]

# Process OD Matrix
od_matrix_long <- melt_od_matrix(od_matrices_gcc$container, nodes = nearest_nodes)

# Run Traffic Assignment
result <- run_assignment(graph, od_matrix_long, cost = "generalized_cost",
                         return.extra = "all")
print(result)

# Visualize Results
network_gcc$final_flows_log10 <- log10(result$final_flows+1)

mapview(network_gcc, zcol = "final_flows_log10",
        layer.name = "Container Flows (Log10)") +
mapview(zones_gcc, alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Zone Nodes")


# --------------------------------------
# GCC Example with Network Consolidation
# --------------------------------------

# Load Graph
graph <- linestrings_to_graph(network_gcc) |>
  create_undirected_graph(by = ~ project + mode)

nodes <- nodes_from_graph(graph, sf = TRUE)

# Map Zones to Nodes
nearest_nodes <- nodes$node[st_nearest_feature(zones_gcc, nodes)]

# Consolidate Graph
graph_cons <- consolidate_graph(graph, by = ~ project + mode, w = ~ length_km, keep = nearest_nodes,
                custom = list(fsum_uw = .c(tariff_cost, tp_cost, tp_time, link_time,
                                           total_cost, total_time, generalized_cost),
                              fmean = .c(speed, tariff),
                              fmode = .c(numlanes, region, country)))

# Examine consolidation result
mapview(linestrings_from_graph(graph), color = "dodgerblue4",
        layer.name = "Original") +
mapview(linestrings_from_graph(graph_cons), color = "orange",
        layer.name = "Consolidated") +
mapview(zones_gcc, alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Zone Nodes")

# Process OD Matrix
od_matrix_long <- melt_od_matrix(od_matrices_gcc$container, nodes = nearest_nodes)

# Run Traffic Assignment
result <- run_assignment(graph_cons, od_matrix_long, cost = "generalized_cost",
                         return.extra = "all")
print(result)

# Visualize Results
graph_cons$final_flows_log10 <- log10(result$final_flows+1)

mapview(linestrings_from_graph(graph_cons), zcol = "final_flows_log10",
        layer.name = "Container Flows (Log10)") +
mapview(zones_gcc, alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Zone Nodes")


# --------------------------------------
# Thailand Example
# --------------------------------------

# Read Network Data
network <- st_read("~/Documents/CPCS/Thailand/CPCSThailandModelling/Modeling_data/Network/TDS_NAM_NET_TDS42.gpkg")
zones <- st_read("~/Documents/CPCS/Thailand/CPCSThailandModelling/Modeling_data/Zones/TDS_NAM_zones.gpkg")

mapview(zones) + mapview(tfm(network, type = qF(linktype)),
        zcol = "type", layer.name = "Linktype",
        color = viridis::turbo)

# Create Undirected Graph
graph <- linestrings_to_graph(network) |>
  create_undirected_graph(by = ~ linktype)

nodes <- nodes_from_graph(graph, sf = TRUE)

# Map Zones to Nodes
nearest_nodes <- nodes$node[st_nearest_feature(st_centroid(zones), nodes)]

# Consolidate Graph (using default aggregation: fmean for numeric, fmode for categorical)
graph_cons <- consolidate_graph(graph, by = ~ linktype, w = ~ distance, keep = nearest_nodes)
all(nearest_nodes %in% c(graph_cons$from, graph_cons$to))

# Examine consolidation result
mapview(linestrings_from_graph(graph), color = "dodgerblue4",
        layer.name = "Original") +
mapview(linestrings_from_graph(graph_cons), color = "orange",
        layer.name = "Consolidated") +
mapview(st_centroid(zones), alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Zone Nodes")

# Process OD Matrix (Trips between zones/districts)
od_person <- fread("~/Documents/CPCS/Thailand/CPCSThailandModelling/Modeling_data/OD/TDS_NAM_total_flows/OD_Person2042_district.csv")
names(od_person) <- c("from", "to", "flow")
od_person$from <- nearest_nodes[ckmatch(od_person$from, zones$zone)]
od_person$to <- nearest_nodes[ckmatch(od_person$to, zones$zone)]

# Run Traffic Assignment (takes 1h without parallelism)
result_trips <- run_assignment(graph_cons, od_person,
                               cost = "distance", nthreads = 10)
print(result_trips)

# Visualize
graph_cons$final_trips_log10 <- log10(result_trips$final_flows+1)
mapview(linestrings_from_graph(graph_cons), zcol = "final_trips_log10",
        layer.name = "Passenger Trips (Log10)") +
mapview(st_centroid(zones), alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Zone Nodes")

# Let's instead use aggregate freight data only available at province level
provinces <- st_centroid(zones) |> tfm(qDF(st_coordinates(geom))) |> atomic_elem() |> qDF() |>
       collap(~ province, w = ~ pop42, custom = list(fsum_uw = .c(area, pop17), fmean = .c(X, Y))) |>
       st_as_sf(coords = .c(X, Y), crs = 4326)
provinces %<>% join(fread("~/Documents/CPCS/Thailand/CPCSThailandModelling/Modeling_data/Zones/provinces.csv"))
mapview(provinces, zcol = "pop42")

nodes_cons <- nodes_from_graph(graph_cons, sf = TRUE)
nearest_nodes <- nodes_cons$node[st_nearest_feature(provinces, nodes_cons)]

# Process OD Matrix (Freight in tons between provinces)
od_flows <- fread("~/Documents/CPCS/Thailand/CPCSThailandModelling/Modeling_data/OD/TDS_NAM_total_flows/OD_Freight2042_province.csv")
names(od_flows) <- c("from", "to", "flow")
od_flows %<>% fsubset(from != 77)
od_flows$from <- nearest_nodes[ckmatch(od_flows$from, provinces$zone)]
od_flows$to <- nearest_nodes[ckmatch(od_flows$to, provinces$zone)]
od_flows$flow <- od_flows$flow / 1e6 # Million tons

# Run Traffic Assignment
result_freight <- run_assignment(graph_cons, od_flows, cost = "distance",
                                 return.extra = "all")
print(result_freight)

# Visualize
graph_cons$final_flows_log10 <- log10(result_freight$final_flows+1)
mapview(linestrings_from_graph(graph_cons), zcol = "final_flows_log10",
        layer.name = "Container Flows (Log10)") +
mapview(provinces, alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Province Population Centroid")

# Simplify Network
graph_simp <- simplify_network(graph_cons, nearest_nodes, cost.column = "distance")$graph
graph_simp <- consolidate_graph(graph_simp, by = ~ linktype, w = ~ distance, keep = nearest_nodes)

result_freight <- run_assignment(graph_simp, od_flows, cost = "distance",
                                 return.extra = "all")
print(result_freight)

# Visualize
graph_simp$final_flows_log10 <- log10(result_freight$final_flows+1)
mapview(linestrings_from_graph(graph_simp), zcol = "final_flows_log10",
        layer.name = "Container Flows (Log10)") +
  mapview(provinces, alpha = 0, cex = 3,
          col.regions = "red", alpha.regions = 0.8,
          layer.name = "Province Population Centroid")

# --------------------------------------
# US FAF5 Network
# --------------------------------------

graph <- st_read("/Users/sebastiankrantz/Documents/CPCS/FAF5/Networks/Geodatabase Format/FAF5Network.gdb",
                layer = "FAF5_Links") |> st_cast("LINESTRING") |> st_transform(4326) |>
                linestrings_to_graph()

zones <- st_read("/Users/sebastiankrantz/Documents/CPCS/FAF5/Networks/Geodatabase Format/FAF5Network.gdb",
                layer = "FAF5_Nodes") |> subset(!is.na(CentroidID)) |> st_transform(4326)

(graph |> subset(STATE == "TX") |>
  linestrings_from_graph() |>
  mapview()) +
mapview(zones |> subset(StateID == 48), alpha = 0, cex = 3,
        col.regions = "red", alpha.regions = 0.8,
        layer.name = "Zone Nodes")

# Map Zones to Nodes
nodes <- nodes_from_graph(graph, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(zones, nodes)]

# Consolidate Graph
graph_cons <- consolidate_graph(graph, w = ~ LENGTH, keep = nearest_nodes)

mapview(linestrings_from_graph(subset(graph, STATE == "TX")), color = "dodgerblue4",
        layer.name = "Original") +
  mapview(linestrings_from_graph(subset(graph_cons, STATE == "TX")), color = "orange",
          layer.name = "Consolidated") +
  mapview(subset(zones, StateID == 48), alpha = 0, cex = 3,
          col.regions = "red", alpha.regions = 0.8,
          layer.name = "Zone Nodes")

nodes <- nodes_from_graph(graph_cons, sf = TRUE)
nearest_nodes <- nodes$node[st_nearest_feature(zones, nodes)]

# Simplify graph
graph_simp <- simplify_network(graph_cons, nearest_nodes, cost = "LENGTH")
graph_simp <- consolidate_graph(graph_simp$graph, w = ~ LENGTH, keep = nearest_nodes)

# Generate OD Matrix
od_matrix <- matrix(abs(rnorm(length(nearest_nodes)^2)), nrow = length(nearest_nodes))
diag(od_matrix) <- 0
od_matrix[od_matrix < 3] <- NA
od_matrix_long <- melt_od_matrix(od_matrix, nodes = nearest_nodes)

# Run Traffic Assignment
result <- run_assignment(graph_simp, od_matrix_long, cost = "LENGTH")
print(result)

# Visualize Results
graph_simp$final_flows_log10 <- log10(result$final_flows+1)

mapview(linestrings_from_graph(graph_simp), zcol = "final_flows_log10",
        layer.name = "Container Flows (Log10)") +
  mapview(zones_gcc, alpha = 0, cex = 3,
          col.regions = "red", alpha.regions = 0.8,
          layer.name = "Zone Nodes")
