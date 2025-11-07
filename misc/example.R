library(mmflowr)
library(sf)

# GCC Scenario
od_matrix <- process_od_matrix(od_matrix_directory, cargo_type, period)
base_network <- st_read("data/network/base_network.shp")
network_gc_cnt <- st_read("data/network/base_network_with_gc_containerized")
zone_nodes <- st_read("data/zone_nodes/network_nodes.shp")

mapview::mapview(network_gc_cnt) + mapview::mapview(zone_nodes, col.regions = "red")

graph_df <- network_gc_cnt |>
  linestring_to_graph() |>
  add_vars(cost = network_gc_cnt$generalize) |>
  create_undirected_graph()

nodes_df <- nodes_from_graph(graph_df)

graph <- makegraph(graph_df |> fselect(from, to, cost), directed = FALSE) # |>
  # cpp_contract() # Could use cpp_simplify() to simplify
nodes <- graph$dict$ref
dmat <- get_distance_matrix(graph, from = nodes, to = nodes, algorithm = "mch")
# dmat <- dist_mat_from_df(graph_df)
dimnames(dmat) <- NULL

get_multi_paths(graph, nodes[1:3], nodes[1:3])

if(!identical(dim(dmat), rep(nrow(nodes_df), 2))) stop("Nodes and Distance Matrix Mismatch")

# Mapping zones to nodes
nearest_nodes <- st_nearest_feature(zone_nodes, st_as_sf(nodes_df, coords = c("X", "Y"), crs = 4326))

# Process OD Matrix
od_matrix <- od_matrix$Container
if(!identical(dim(od_matrix), rep(nrow(zone_nodes), 2))) stop("Zones and OD-Matrix Matrix Mismatch")

od_matrix_long <- data.frame(from = rep.int(nearest_nodes, ncol(od_matrix)),
                             to = rep(nearest_nodes, each = nrow(od_matrix)),
                             flow = vec(od_matrix)) |>
                  fsubset(is.finite(flow) & flow > 0)


# Experiment with dodgr
system.time(dodgr::dodgr_paths(graph_df |> fselect(from, to, dist = cost),
                   from = seq_row(nodes_df), to = seq_row(nodes_df), vertices = FALSE))
res <- dodgr::dodgr_paths(graph_df |> fselect(from, to, d = cost),
                          from = ks, to = ks, vertices = FALSE)


# Experiment with igraph: best so far
library(igraph)
g <- graph_from_data_frame(graph_df |> fselect(from, to), directed = FALSE)
g <- delete_vertex_attr(g, "name")
igraph_options(return.vs.es = FALSE) # sparsematrices = TRUE
# Distance Matrix
D <- distances(g, mode = "out", weights = graph_df$cost)
# Shortest Paths
shortest_paths(g, from = i, to = ks, weights = graph_df$cost, output = "epath")$epath

cost <- graph_df$cost
delta_ks <- integer(nrow(graph_df))
# Final flows vector
final_flows <- numeric(nrow(graph_df))

attach(od_matrix_long)
# Now iterating across OD-pairs
# TODO: could restrict that other nodes must be in the direction of travel and not behind destination node
system.time({
for (i in seq_row(od_matrix_long)) {
  d_ij <- dmat[from[i], to[i]]
  d_ikj <- dmat[from[i], ] + dmat[, to[i]] # from i to all other nodes k and from these nodes k to j (basically dmat + t(dmat)?)
  short_detour_ij <- d_ikj < 1.1 * d_ij
  short_detour_ij[d_ikj <= d_ij + .Machine$double.eps*1e3] <- FALSE # Exclude nodes k that are on the shortest path
  # which(d_ij == d_ikj) # These are the nodes on the direct path from i to j which yield the shortest distance.
  ks <- which(short_detour_ij)
  cost_ks <- d_ikj[ks]
  # Now Path-Sized Logit: Need to compute overlap between routes
  # Could still optimize calls to shortest_paths(), e.g., go to C directly.
  # We add the shortest path at the end of paths1
  paths1 <- shortest_paths(g, from = from[i], to = c(ks, to[i]), weights = cost, output = "epath")$epath
  paths2 <- shortest_paths(g, from = to[i], to = ks, weights = cost, output = "epath")$epath
  shortest_path <- paths1[[length(paths1)]]

  # # Check
  # cost_ks[k] == sum(cost[paths1[[k]]]) + sum(cost[paths2[[k]]])

  # Get indices of paths that do not contain duplicate edges
  no_dups <- mmflowr:::check_path_duplicates(paths1, paths2, delta_ks)

  # # Number of routes in choice set that use link j
  # for (k in no_dups) {
  #   delta_ks[paths1[[k]]] <- delta_ks[paths1[[k]]] + 1L
  #   delta_ks[paths2[[k]]] <- delta_ks[paths2[[k]]] + 1L
  # }
  # delta_ks[shortest_path] <- delta_ks[shortest_path] + 1L
  #
  # # Correction factors for each route k
  # gamma_ks <- sapply(no_dups, function(k) {
  #   path <- c(paths1[[k]], paths2[[k]])
  #   sum(cost[path] / delta_ks[path]) / cost_ks[k]
  # })
  # gamma_1 <- sum(cost[shortest_path] / delta_ks[shortest_path]) / d_ij
  #
  # # Now the PS-MNL
  # prob_ks <- proportions(exp(c(cost_ks[no_dups], d_ij) + beta_PSL * c(gamma_ks, gamma_1)))
  #
  # # Need to reset delta_ks
  # delta_ks[] <- 0L
  #
  # # Assign result to edges:
  # for (k in no_dups) {
  #   final_flows[paths1[[k]]] <- final_flows[paths1[[k]]] + flow[i] * prob_ks[k]
  # }
  # final_flows[shortest_path] <- final_flows[shortest_path] + flow[i] * prob_ks[length(prob_ks)]
  mmflowr:::compute_path_sized_logit(paths1, paths2, no_dups, shortest_path,
                                     cost, cost_ks, d_ij, beta_PSL, flow[i],
                                     delta_ks, final_flows)
}
})
detach(od_matrix_long)

# US Freight Analysis Framework 5.0 Model Network Database
st_layers(dsn = "/Users/sebastiankrantz/Downloads/Networks/Geodatabase Format/FAF5Network.gdb")
nodes <- st_read(dsn = "/Users/sebastiankrantz/Downloads/Networks/Geodatabase Format/FAF5Network.gdb", layer = "FAF5_Nodes")
edges <- st_read(dsn = "/Users/sebastiankrantz/Downloads/Networks/Geodatabase Format/FAF5Network.gdb", layer = "FAF5_Links")


#Choose number of cores used by cppRouting
RcppParallel::setThreadOptions(numThreads = 1)

if(requireNamespace("igraph",quietly = TRUE)){

  #Generate fully connected graph
  gf<- igraph::make_full_graph(400)
  igraph::V(gf)$names<-1:400

  #Convert to data frame and add random weights
  df<-igraph::as_long_data_frame(gf)
  df$dist<-sample(1:100,nrow(df),replace = TRUE)

  #Construct cppRouting graph
  graph<-makegraph(df[,c(1,2,5)],directed = FALSE)

  #Pick up random origin and destination node
  origin<-sample(1:400,1)
  destination<-sample(1:400,1)

  #Compute distance from origin to all nodes
  or_to_all<-get_distance_matrix(graph,from=origin,to=1:400)

  #Compute distance from all nodes to destination
  all_to_dest<-get_distance_matrix(graph,from=1:400,to=destination,)

  #Get all shortest paths from origin to destination, passing by each node of the graph
  total_paths<-rowSums(cbind(t(or_to_all),all_to_dest))

  #Compute shortest path between origin and destination
  distance<-get_distance_pair(graph,from=origin,to=destination)

  #Compute detour with an additional cost of 3
  det<-get_detour(graph,from=origin,to=destination,extra=3)

  #Check result validity
  length(unlist(det))
  length(total_paths[total_paths < distance + 3])
