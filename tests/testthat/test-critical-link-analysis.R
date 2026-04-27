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
})
