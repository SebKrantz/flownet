# Tests for consolidate_graph() and simplify_network()

# --- consolidate_graph() Tests ---

test_that("consolidate_graph reduces edge count for intermediate nodes", {
  # Linear chain: 1 -> 2 -> 3 -> 4 (node 2 and 3 are intermediate)
  graph <- data.frame(
    from = c(1, 2, 3),
    to = c(2, 3, 4),
    cost = c(1, 2, 3)
  )

  result <- consolidate_graph(graph, verbose = FALSE)

  # Should consolidate to single edge 1 -> 4
  expect_lt(nrow(result), nrow(graph))
})

test_that("consolidate_graph keeps specified nodes", {
  # Linear chain with node 2 preserved
  graph <- data.frame(
    from = c(1, 2, 3),
    to = c(2, 3, 4),
    cost = c(1, 2, 3)
  )

  # Keep nodes 1, 2, and 4 (endpoints and one intermediate)
  result <- consolidate_graph(graph, keep.nodes = c(1, 2, 4), verbose = FALSE)

  # Node 2 should still exist in the result
  all_nodes <- unique(c(result$from, result$to))
  expect_true(2 %in% all_nodes)
})

test_that("consolidate_graph recursive='full' completes consolidation", {
  # Chain that needs multiple passes, keep endpoints
  graph <- data.frame(
    from = c(1, 2, 3, 4, 5),
    to = c(2, 3, 4, 5, 6),
    cost = c(1, 1, 1, 1, 1)
  )

  result <- consolidate_graph(graph, keep.nodes = c(1, 6),
                              recursive = "full", verbose = FALSE)

  # Should consolidate to single edge 1->6
  expect_equal(nrow(result), 1)
  expect_equal(result$from, 1)
  expect_equal(result$to, 6)
})

test_that("consolidate_graph removes loops", {
  graph <- data.frame(
    from = c(1, 2, 3),
    to = c(2, 2, 4),  # Edge 2->2 is a loop
    cost = c(1, 5, 2)
  )

  result <- consolidate_graph(graph, drop.edges = "loop",
                              consolidate = FALSE, verbose = FALSE)

  # Loop should be removed
  expect_false(any(result$from == result$to))
  expect_lt(nrow(result), nrow(graph))
})

test_that("consolidate_graph removes duplicates", {
  graph <- data.frame(
    from = c(1, 1, 2),
    to = c(2, 2, 3),  # Two 1->2 edges
    cost = c(1, 2, 3)
  )

  result <- consolidate_graph(graph, drop.edges = "duplicate",
                              consolidate = FALSE, verbose = FALSE)

  # Duplicates should be removed
  expect_lt(nrow(result), nrow(graph))
})

test_that("consolidate_graph removes singleton edges", {
  graph <- data.frame(
    from = c(1, 2, 3, 5),  # Node 5->6 is a dead end
    to = c(2, 3, 4, 6),
    cost = c(1, 1, 1, 1)
  )

  result <- consolidate_graph(graph, drop.edges = "single", verbose = FALSE)

  # Singleton edge should be removed
  all_nodes <- c(result$from, result$to)
  expect_false(5 %in% all_nodes)
  expect_false(6 %in% all_nodes)
})

test_that("consolidate_graph by parameter preserves mode groups", {
  graph <- data.frame(
    from = c(1, 2, 1, 2),
    to = c(2, 3, 2, 3),
    mode = c("road", "road", "rail", "rail"),
    cost = c(1, 2, 3, 4)
  )

  result <- consolidate_graph(graph, by = ~ mode, verbose = FALSE)

  # Should not consolidate across modes
  expect_true("mode" %in% names(result))
})

test_that("consolidate_graph adds edge column", {
  # Graph with branching (not just a chain)
  graph <- data.frame(
    from = c(1, 2, 2, 3, 4),
    to = c(2, 3, 4, 5, 5),
    cost = c(1, 2, 3, 1, 1)
  )

  result <- consolidate_graph(graph, keep.nodes = c(1, 5), verbose = FALSE)

  expect_true("edge" %in% names(result))
  expect_true(nrow(result) > 0)
})

test_that("consolidate_graph aggregates with weights", {
  # Chain with endpoints preserved
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 3),
    cost = c(10, 20),
    weight = c(1, 3)
  )

  result <- consolidate_graph(graph, keep.nodes = c(1, 3),
                              w = ~ weight, verbose = FALSE)

  # Weighted mean of 10 and 20 with weights 1 and 3 = (10*1 + 20*3) / 4 = 17.5
  expect_equal(nrow(result), 1)
  expect_equal(result$cost, 17.5, tolerance = 0.1)
})

# --- simplify_network() shortest-paths Tests ---

test_that("simplify_network shortest-paths returns subset of edges", {
  # Simple graph
  graph <- data.frame(
    from = c(1, 2, 1, 3),
    to = c(2, 3, 3, 4),
    cost = c(1, 1, 10, 1)
  )

  # Simplify keeping paths between nodes 1 and 4
  result <- simplify_network(graph, nodes = c(1, 4),
                             method = "shortest-paths",
                             cost.column = "cost")

  expect_lte(nrow(result), nrow(graph))
})

test_that("simplify_network shortest-paths has edges attribute", {
  graph <- data.frame(
    from = c(1, 2, 1, 3),
    to = c(2, 3, 3, 4),
    cost = c(1, 1, 10, 1)
  )

  result <- simplify_network(graph, nodes = c(1, 4),
                             method = "shortest-paths",
                             cost.column = "cost")

  expect_true(!is.null(attr(result, "edges")))
  expect_true(all(attr(result, "edges") <= nrow(graph)))
})

test_that("simplify_network shortest-paths has edge_counts attribute", {
  graph <- data.frame(
    from = c(1, 2, 1, 3),
    to = c(2, 3, 3, 4),
    cost = c(1, 1, 10, 1)
  )

  result <- simplify_network(graph, nodes = c(1, 4),
                             method = "shortest-paths",
                             cost.column = "cost")

  edge_counts <- attr(result, "edge_counts")
  expect_true(!is.null(edge_counts))
  expect_true(all(edge_counts > 0))
})

test_that("simplify_network shortest-paths keeps shortest path edges", {
  # Graph where 1->2->3 is shorter than 1->3 directly
  graph <- data.frame(
    from = c(1, 2, 1),
    to = c(2, 3, 3),
    cost = c(1, 1, 10)
  )

  result <- simplify_network(graph, nodes = c(1, 3),
                             method = "shortest-paths",
                             cost.column = "cost")

  # Should keep edges 1->2 and 2->3, may or may not keep 1->3
  expect_gte(nrow(result), 2)
})

test_that("simplify_network shortest-paths with OD pairs data.frame", {
  graph <- data.frame(
    from = c(1, 2, 2, 3),
    to = c(2, 3, 4, 4),
    cost = c(1, 2, 3, 1)
  )

  od_pairs <- data.frame(from = c(1, 2), to = c(4, 4))

  result <- simplify_network(graph, nodes = od_pairs,
                             method = "shortest-paths",
                             cost.column = "cost")

  expect_true(nrow(result) > 0)
})

test_that("simplify_network errors on missing columns", {
  graph <- data.frame(from = 1:3, cost = 1:3)

  expect_error(
    simplify_network(graph, nodes = c(1, 3), cost.column = "cost"),
    "to"
  )
})

test_that("simplify_network errors on unknown nodes", {
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 3),
    cost = c(1, 1)
  )

  expect_error(
    simplify_network(graph, nodes = c(1, 99), cost.column = "cost"),
    "Unknown"
  )
})

# --- simplify_network() cluster Tests ---

test_that("simplify_network cluster returns contracted graph", {
  # Convert africa_segments to a proper graph (it only has coordinates)
  graph <- linestrings_from_graph(africa_segments[1:100, ]) |>
    linestrings_to_graph()

  # Get some nodes to preserve
  nodes_df <- nodes_from_graph(graph)
  keep_nodes <- nodes_df$node[1:5]

  result <- simplify_network(graph, nodes = keep_nodes,
                             method = "cluster",
                             cost.column = ".length",
                             radius_km = list(nodes = 50, cluster = 100))

  # Should have fewer edges after clustering
  expect_lt(nrow(result), nrow(graph))
})

test_that("simplify_network cluster has no self-loops", {
  graph <- linestrings_from_graph(africa_segments[1:100, ]) |>
    linestrings_to_graph()
  nodes_df <- nodes_from_graph(graph)
  keep_nodes <- nodes_df$node[1:5]

  result <- simplify_network(graph, nodes = keep_nodes,
                             method = "cluster",
                             cost.column = ".length",
                             radius_km = list(nodes = 50, cluster = 100))

  # No self-loops
  expect_false(any(result$from == result$to))
})

test_that("simplify_network cluster adds group attributes", {
  graph <- linestrings_from_graph(africa_segments[1:100, ]) |>
    linestrings_to_graph()
  nodes_df <- nodes_from_graph(graph)
  keep_nodes <- nodes_df$node[1:5]

  result <- simplify_network(graph, nodes = keep_nodes,
                             method = "cluster",
                             cost.column = ".length",
                             radius_km = list(nodes = 50, cluster = 100))

  expect_true(!is.null(attr(result, "group.id")))
  expect_true(!is.null(attr(result, "group.starts")))
})

test_that("simplify_network cluster errors without coordinate columns", {
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 3),
    cost = c(1, 1)
  )

  expect_error(
    simplify_network(graph, nodes = c(1, 3),
                     method = "cluster",
                     cost.column = "cost"),
    "FX.*FY.*TX.*TY"
  )
})
