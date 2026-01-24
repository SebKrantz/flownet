# Tests for graph utility functions

# Simple test graph
simple_graph <- data.frame(
  from = c(1, 2, 2, 3),
  to = c(2, 3, 4, 4),
  cost = c(1, 2, 3, 1)
)

# Graph with coordinates
graph_with_coords <- data.frame(
  from = c(1, 2, 3),
  to = c(2, 3, 4),
  FX = c(0, 1, 2),
  FY = c(0, 0, 0),
  TX = c(1, 2, 3),
  TY = c(0, 0, 0),
  cost = c(1, 1, 1)
)

# --- normalize_graph() Tests ---

test_that("normalize_graph converts non-consecutive IDs to 1:n", {
  graph <- data.frame(
    from = c(10, 20, 20),
    to = c(20, 30, 40),
    cost = c(1, 2, 3)
  )

  result <- normalize_graph(graph)

  # Should have nodes 1, 2, 3, 4 (from 10, 20, 30, 40)
  expect_equal(sort(unique(c(result$from, result$to))), 1:4)
})

test_that("normalize_graph preserves other columns", {
  graph <- data.frame(
    from = c(10, 20),
    to = c(20, 30),
    cost = c(5, 10),
    name = c("a", "b")
  )

  result <- normalize_graph(graph)

  expect_equal(result$cost, c(5, 10))
  expect_equal(result$name, c("a", "b"))
})

test_that("normalize_graph errors on missing columns", {
  bad_graph <- data.frame(from = c(1, 2), cost = c(1, 2))

  expect_error(normalize_graph(bad_graph), "from.*to")
})

test_that("normalize_graph handles already normalized graph", {
  result <- normalize_graph(simple_graph)

  # Should remain unchanged if already 1:n
  expect_equal(result$from, simple_graph$from)
  expect_equal(result$to, simple_graph$to)
})

# --- nodes_from_graph() Tests ---

test_that("nodes_from_graph returns correct number of unique nodes", {
  result <- nodes_from_graph(graph_with_coords)

  # Graph has nodes 1, 2, 3, 4

  expect_equal(nrow(result), 4)
  expect_equal(sort(result$node), 1:4)
})

test_that("nodes_from_graph returns correct columns", {
  result <- nodes_from_graph(graph_with_coords)

  expect_true(all(c("node", "X", "Y") %in% names(result)))
})

test_that("nodes_from_graph sf=TRUE returns sf object", {
  result <- nodes_from_graph(graph_with_coords, sf = TRUE)

  expect_s3_class(result, "sf")
  expect_equal(sf::st_crs(result)$epsg, 4326)
})

test_that("nodes_from_graph sf=TRUE with custom CRS", {
  result <- nodes_from_graph(graph_with_coords, sf = TRUE, crs = 3857)

  expect_s3_class(result, "sf")
  expect_equal(sf::st_crs(result)$epsg, 3857)
})

test_that("nodes_from_graph result is sorted by node ID", {
  result <- nodes_from_graph(graph_with_coords)

  expect_equal(result$node, sort(result$node))
})

# --- distances_from_graph() Tests ---

test_that("distances_from_graph returns square matrix", {
  result <- distances_from_graph(simple_graph, cost.column = "cost")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), ncol(result))
})

test_that("distances_from_graph has correct dimensions", {
  result <- distances_from_graph(simple_graph, cost.column = "cost")

  # Graph has 4 nodes
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
})

test_that("distances_from_graph diagonal is zero for undirected", {
  result <- distances_from_graph(simple_graph, directed = FALSE, cost.column = "cost")

  expect_equal(unname(diag(result)), rep(0, 4))
})

test_that("distances_from_graph is symmetric for undirected", {
  result <- distances_from_graph(simple_graph, directed = FALSE, cost.column = "cost")

  expect_equal(result, t(result))
})

test_that("distances_from_graph matches expected distances", {
  # Simple graph: 1->2 (cost 1), 2->3 (cost 2), 2->4 (cost 3), 3->4 (cost 1)
  result <- distances_from_graph(simple_graph, directed = FALSE, cost.column = "cost")

  # Distance from 1 to 4: 1->2->4 = 1+3 = 4, or 1->2->3->4 = 1+2+1 = 4
  expect_equal(result["1", "4"], 4)

  # Distance from 1 to 2
  expect_equal(result["1", "2"], 1)

  # Distance from 2 to 3
  expect_equal(result["2", "3"], 2)
})

test_that("distances_from_graph works with numeric cost.column vector", {
  costs <- c(1, 2, 3, 1)
  result <- distances_from_graph(simple_graph, cost.column = costs)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
})

test_that("distances_from_graph errors on invalid cost.column", {
  expect_error(
    distances_from_graph(simple_graph, cost.column = "nonexistent"),
    "cost.column"
  )
})

test_that("distances_from_graph directed=TRUE gives asymmetric result",
{
  result <- distances_from_graph(simple_graph, directed = TRUE, cost.column = "cost")

  # In a directed graph, distance 1->4 may differ from 4->1
 # 4->1 should be Inf in directed case since no path exists
 expect_equal(result["4", "1"], Inf)
})
