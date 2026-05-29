# Tests for critical_link_analysis()

test_that("critical_link_analysis computes finite undirected detours", {
  graph <- data.frame(
    from = c(1, 2, 1),
    to = c(2, 3, 3),
    cost = c(1, 1, 3)
  )

  result <- critical_link_analysis(graph, cost.column = "cost")

  expect_equal(result$edge, 1:3)
  expect_equal(result$edge_cost, graph$cost)
  expect_equal(result$detour_cost, c(4, 4, 2))
  expect_equal(result$detour_ratio, c(4, 4, 2 / 3))
})

test_that("critical_link_analysis returns Inf for bridge detours", {
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 3),
    cost = c(1, 1)
  )

  result <- critical_link_analysis(graph, cost.column = "cost")

  expect_equal(result$detour_cost, c(Inf, Inf))
  expect_equal(result$detour_ratio, c(Inf, Inf))
})

test_that("critical_link_analysis treats parallel edges as distinct alternatives", {
  graph <- data.frame(
    from = c(1, 1),
    to = c(2, 2),
    cost = c(5, 7)
  )

  result <- critical_link_analysis(graph, directed = TRUE, cost.column = "cost")

  expect_equal(result$detour_cost, c(7, 5))
  expect_equal(result$detour_ratio, c(7 / 5, 5 / 7))
})

test_that("critical_link_analysis respects directed and undirected traversal", {
  graph <- data.frame(
    from = c(1, 1, 3),
    to = c(2, 3, 2),
    cost = c(1, 2, 2)
  )

  directed_result <- critical_link_analysis(graph, directed = TRUE, cost.column = "cost")
  undirected_result <- critical_link_analysis(graph, directed = FALSE, cost.column = "cost")

  expect_equal(directed_result$detour_cost[1], 4)
  expect_equal(directed_result$detour_cost[2], Inf)
  expect_equal(undirected_result$detour_cost[2], 3)
})

test_that("critical_link_analysis does not allow removed directed edge after returning to source", {
  graph <- data.frame(
    from = c(1, 1, 3),
    to = c(2, 3, 1),
    cost = c(10, 1, 1)
  )

  result <- critical_link_analysis(graph, directed = TRUE, cost.column = "cost")

  expect_equal(result$detour_cost[1], Inf)
})

test_that("critical_link_analysis preserves row order and original node IDs", {
  graph <- data.frame(
    from = c(10, 10, 20),
    to = c(20, 30, 30),
    cost = c(2, 5, 1)
  )

  result <- critical_link_analysis(graph, cost.column = "cost")

  expect_equal(result$edge, 1:3)
  expect_equal(result$from, graph$from)
  expect_equal(result$to, graph$to)
  expect_equal(result$detour_cost, c(6, 3, 7))
})

test_that("critical_link_analysis accepts numeric cost vectors", {
  graph <- data.frame(
    from = c(1, 2, 1),
    to = c(2, 3, 3)
  )

  result <- critical_link_analysis(graph, cost.column = c(1, 1, 3))

  expect_equal(result$edge_cost, c(1, 1, 3))
  expect_equal(result$detour_cost, c(4, 4, 2))
})

test_that("critical_link_analysis validates inputs", {
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 3),
    cost = c(1, 1)
  )

  expect_error(critical_link_analysis(graph, cost.column = "missing"), "cost.column")
  expect_error(critical_link_analysis(transform(graph, cost = c(1, NA))), "finite")
  expect_error(critical_link_analysis(transform(graph, cost = c(1, 0))), "strictly positive")
  expect_error(critical_link_analysis(graph[, c("from", "cost")]), "from.*to")
  expect_error(critical_link_analysis(graph, directed = NA), "directed")
  expect_error(critical_link_analysis(graph, nthreads = 0L), "nthreads")
  expect_error(critical_link_analysis(graph, nthreads = c(1L, 2L)), "nthreads")
})

test_that("critical_link_analysis matches single-threaded result with nthreads > 1", {
  set.seed(123)
  n_nodes <- 20
  graph <- data.frame(
    from = sample.int(n_nodes, 60, replace = TRUE),
    to = sample.int(n_nodes, 60, replace = TRUE),
    cost = runif(60, 0.1, 5)
  )
  graph <- graph[graph$from != graph$to, ]

  serial <- critical_link_analysis(graph, cost.column = "cost", nthreads = 1L)
  parallel <- critical_link_analysis(graph, cost.column = "cost", nthreads = 4L)

  expect_equal(parallel, serial)
})

test_that("critical_link_analysis only runs Dijkstra for unique from nodes", {
  graph <- data.frame(
    from = c(1, 1, 2, 2),
    to = c(2, 3, 1, 3),
    cost = c(1, 5, 1, 1)
  )

  # Node 3 is never a 'from' node, so it should not trigger a Dijkstra run.
  # e1=1->2: remove e1 -> only outgoing from 1 left is e2 (1->3); 3 has no outgoing -> Inf
  # e2=1->3: remove e2 -> 1->2 (1) + 2->3 (1) = 2
  # e3=2->1: remove e3 -> only outgoing from 2 left is e4 (2->3); 3 has no outgoing -> Inf
  # e4=2->3: remove e4 -> 2->1 (1) + 1->3 (5) = 6
  result <- critical_link_analysis(graph, directed = TRUE, cost.column = "cost")

  expect_equal(result$detour_cost, c(Inf, 2, Inf, 6))
})
