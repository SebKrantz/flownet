graph <- makegraph(graph_df |> fselect(from, to, cost), directed = FALSE) # |>
# cpp_contract() # Could use cpp_simplify() to simplify
nodes <- graph$dict$ref
dmat <- get_distance_matrix(graph, from = nodes, to = nodes, algorithm = "mch")
# dmat <- dist_mat_from_df(graph_df)
dimnames(dmat) <- NULL

get_multi_paths(graph, nodes[1:3], nodes[1:3])




# Experiment with dodgr
system.time(dodgr::dodgr_paths(graph_df |> fselect(from, to, dist = cost),
                               from = seq_row(nodes_df), to = seq_row(nodes_df), vertices = FALSE))
res <- dodgr::dodgr_paths(graph_df |> fselect(from, to, d = cost),
                          from = ks, to = ks, vertices = FALSE)



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
