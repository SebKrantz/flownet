# Tests for consolidate_graph() and simplify_network()

# Strip volatile attrs and sort rows so C vs R contraction paths are comparable.
canonical_consolidate_df <- function(x) {
  x <- as.data.frame(x)
  for (a in c("keep.edges", "group.id", ".early.return")) attr(x, a) <- NULL
  nm <- intersect(names(x), c("from", "to", "mode", "class", "cost", "weight"))
  ord <- if (length(nm)) do.call(order, x[nm]) else seq_len(nrow(x))
  x <- x[ord, , drop = FALSE]
  rownames(x) <- NULL
  x
}

expect_consolidate_equiv_c_r <- function(graph, ..., tol = 1e-9) {
  o <- options()
  on.exit(options(o), add = TRUE)
  options(flownet.use_c_contraction = TRUE)
  res_c <- consolidate_graph(graph, ..., verbose = FALSE)
  options(flownet.use_c_contraction = FALSE)
  res_r <- consolidate_graph(graph, ..., verbose = FALSE)
  expect_equal(
    canonical_consolidate_df(res_c),
    canonical_consolidate_df(res_r),
    tolerance = tol
  )
}

legacy_drop_singletons_recursive <- function(from_vec, to_vec, keep.nodes = NULL, recursive = TRUE) {
  n_edges <- length(from_vec)
  if(!n_edges) return(list(keep = integer(0), dropped = 0L, rounds = 0L))
  keep <- seq_len(n_edges)
  rounds <- 0L
  repeat {
    rounds <- rounds + 1L
    cnt <- unclass(kit::countOccur(c(from_vec, to_vec)))
    deg1 <- cnt$Variable[cnt$Count == 1L]
    if(length(keep.nodes)) deg1 <- setdiff(deg1, keep.nodes)
    if(!length(deg1)) break
    drop <- (from_vec %in% deg1) | (to_vec %in% deg1)
    if(!any(drop)) break
    keep <- keep[!drop]
    from_vec <- from_vec[!drop]
    to_vec <- to_vec[!drop]
    if(!recursive || !length(from_vec)) break
  }
  list(keep = keep, dropped = as.integer(n_edges - length(keep)), rounds = rounds)
}

expect_singleton_peel_equiv <- function(from_vec, to_vec, keep.nodes = NULL, recursive = TRUE) {
  ref <- legacy_drop_singletons_recursive(from_vec, to_vec, keep.nodes, recursive)
  cur <- drop_singletons_linear(from_vec, to_vec, keep.nodes, recursive)
  expect_equal(cur$keep, ref$keep)
  expect_equal(cur$dropped, ref$dropped)
}

test_that("consolidate_graph C vs R contraction are equivalent (simple chain)", {
  graph <- data.frame(
    from = c(1L, 2L, 3L),
    to = c(2L, 3L, 4L),
    cost = c(1, 2, 3)
  )
  expect_consolidate_equiv_c_r(graph)
})

test_that("consolidate_graph C vs R contraction are equivalent (full recursion chain)", {
  graph <- data.frame(
    from = c(1L, 2L, 3L, 4L, 5L),
    to = c(2L, 3L, 4L, 5L, 6L),
    cost = c(1, 1, 1, 1, 1)
  )
  expect_consolidate_equiv_c_r(
    graph,
    keep.nodes = c(1L, 6L),
    recursive = "full"
  )
})

test_that("consolidate_graph C vs R contraction are equivalent (parallel 3-edge bundles)", {
  k <- 150L
  u1 <- seq.int(3L, k + 2L)
  u2 <- seq.int(k + 3L, 2L * k + 2L)
  graph <- data.frame(
    from = as.integer(c(rep.int(1L, k), u1, u2)),
    to = as.integer(c(u1, u2, rep.int(2L, k))),
    mode = "road",
    class = "a",
    cost = 1,
    stringsAsFactors = FALSE
  )
  expect_consolidate_equiv_c_r(
    graph,
    drop.edges = "loop",
    by = ~ mode + class,
    keep.nodes = c(1L, 2L),
    recursive = "partial"
  )
})

test_that("consolidate_graph C vs R contraction are equivalent (by + keep.nodes)", {
  graph <- data.frame(
    from = c(1L, 2L, 3L, 4L),
    to = c(2L, 3L, 4L, 5L),
    mode = c("road", "road", "rail", "rail"),
    class = c("a", "a", "b", "b"),
    cost = c(1, 1, 1, 1)
  )
  expect_consolidate_equiv_c_r(
    graph,
    by = ~ mode + class,
    keep.nodes = c(1L, 5L),
    recursive = "partial"
  )
})

test_that("consolidate_graph C vs R contraction are equivalent (weighted chain)", {
  graph <- data.frame(
    from = c(1L, 2L),
    to = c(2L, 3L),
    cost = c(10, 20),
    weight = c(1, 3)
  )
  expect_consolidate_equiv_c_r(graph, keep.nodes = c(1L, 3L), w = ~ weight)
})

test_that("drop_singletons_linear matches legacy recursive peeling", {
  # Tail attached to triangle core.
  g1 <- data.frame(
    from = c(1L, 2L, 3L, 3L, 4L),
    to = c(2L, 3L, 1L, 4L, 5L)
  )
  expect_singleton_peel_equiv(g1$from, g1$to, recursive = TRUE)
  expect_singleton_peel_equiv(g1$from, g1$to, recursive = FALSE)

  # Keep-node suppression on a long chain.
  g2 <- data.frame(
    from = c(1L, 2L, 3L, 4L, 5L, 6L),
    to = c(2L, 3L, 4L, 5L, 6L, 7L)
  )
  expect_singleton_peel_equiv(g2$from, g2$to, keep.nodes = c(1L, 4L, 7L), recursive = TRUE)
  expect_singleton_peel_equiv(g2$from, g2$to, keep.nodes = c(1L, 4L, 7L), recursive = FALSE)

  # Deterministic random sparse graph.
  set.seed(42)
  n <- 300L
  g3 <- data.frame(
    from = sample.int(80L, n, replace = TRUE),
    to = sample.int(80L, n, replace = TRUE)
  )
  g3 <- g3[g3$from != g3$to, , drop = FALSE]
  expect_singleton_peel_equiv(g3$from, g3$to, recursive = TRUE)
  expect_singleton_peel_equiv(g3$from, g3$to, recursive = FALSE)
  expect_singleton_peel_equiv(g3$from, g3$to, keep.nodes = sample.int(80L, 5L), recursive = TRUE)
  expect_singleton_peel_equiv(g3$from, g3$to, keep.nodes = sample.int(80L, 5L), recursive = FALSE)
})

test_that("drop_singletons_linear matches legacy on varied graph topologies", {
  # Star graph: all edges are singleton incidents in one round.
  star <- data.frame(
    from = rep.int(1L, 12L),
    to = 2L:13L
  )
  expect_singleton_peel_equiv(star$from, star$to, recursive = TRUE)
  expect_singleton_peel_equiv(star$from, star$to, recursive = FALSE)
  expect_singleton_peel_equiv(star$from, star$to, keep.nodes = 1L, recursive = TRUE)
  expect_singleton_peel_equiv(star$from, star$to, keep.nodes = 1L, recursive = FALSE)

  # Disconnected components with one stable cycle and one peelable path.
  mixed <- data.frame(
    from = c(1L, 2L, 3L, 4L, 10L, 11L, 12L, 12L),
    to = c(2L, 3L, 4L, 5L, 11L, 12L, 10L, 13L)
  )
  expect_singleton_peel_equiv(mixed$from, mixed$to, recursive = TRUE)
  expect_singleton_peel_equiv(mixed$from, mixed$to, recursive = FALSE)
  expect_singleton_peel_equiv(mixed$from, mixed$to, keep.nodes = c(5L, 13L), recursive = TRUE)
  expect_singleton_peel_equiv(mixed$from, mixed$to, keep.nodes = c(5L, 13L), recursive = FALSE)

  # Dense-ish random graph with reproducible seeds, multiple replicates.
  for (s in c(7L, 11L, 19L, 23L)) {
    set.seed(s)
    g <- data.frame(
      from = sample.int(120L, 700L, replace = TRUE),
      to = sample.int(120L, 700L, replace = TRUE)
    )
    g <- g[g$from != g$to, , drop = FALSE]
    keep_nodes <- sample.int(120L, 8L)
    expect_singleton_peel_equiv(g$from, g$to, recursive = TRUE)
    expect_singleton_peel_equiv(g$from, g$to, recursive = FALSE)
    expect_singleton_peel_equiv(g$from, g$to, keep.nodes = keep_nodes, recursive = TRUE)
    expect_singleton_peel_equiv(g$from, g$to, keep.nodes = keep_nodes, recursive = FALSE)
  }
})

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

test_that("consolidate_graph handles topology-only graphs", {
  graph <- data.frame(
    from = c(1L, 2L, 3L),
    to = c(2L, 3L, 4L)
  )

  result <- consolidate_graph(graph, keep.nodes = c(1L, 4L), verbose = FALSE)

  expect_named(result, c("from", "to"))
  expect_equal(nrow(result), 1L)
  expect_equal(result$from, 1L)
  expect_equal(result$to, 4L)
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

test_that("consolidate_graph recursive modes differ on singleton tail peeling", {
  # Triangle core plus a two-edge tail from node 3
  graph <- data.frame(
    from = c(1, 2, 3, 3, 4),
    to = c(2, 3, 1, 4, 5),
    cost = 1
  )

  res_none <- consolidate_graph(graph, drop.edges = "single",
                                recursive = "none", contract = FALSE,
                                verbose = FALSE)
  res_partial <- consolidate_graph(graph, drop.edges = "single",
                                   recursive = "partial", contract = FALSE,
                                   verbose = FALSE)
  res_full <- consolidate_graph(graph, drop.edges = "single",
                                recursive = "full", contract = FALSE,
                                verbose = FALSE)

  # One peel for "none" keeps one tail edge, recursive modes remove both.
  expect_equal(nrow(res_none), 4L)
  expect_equal(nrow(res_partial), 3L)
  expect_equal(nrow(res_full), 3L)
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

test_that("consolidate_graph keeps by-boundary nodes under partial recursion", {
  graph <- data.frame(
    from = c(1, 2, 3, 4),
    to = c(2, 3, 4, 5),
    mode = c("road", "road", "rail", "rail"),
    class = c("a", "a", "b", "b"),
    cost = c(1, 1, 1, 1)
  )

  result <- consolidate_graph(graph, by = ~ mode + class,
                              keep.nodes = c(1, 5),
                              recursive = "partial",
                              verbose = FALSE)

  all_nodes <- unique(c(result$from, result$to))
  # Node 3 is a by-boundary (road->rail) and should be preserved.
  expect_true(3 %in% all_nodes)
})

test_that("consolidate_graph handles long singleton tails robustly", {
  core <- data.frame(
    from = c(1, 2, 3),
    to = c(2, 3, 1),
    cost = 1
  )
  # Attach a long dead-end chain to node 1.
  tail_n <- 250
  tail <- data.frame(
    from = c(1, 4:(tail_n + 2)),
    to = c(4, 5:(tail_n + 3)),
    cost = 1
  )
  graph <- rbind(core, tail)

  result <- consolidate_graph(graph, drop.edges = "single",
                              recursive = "partial", contract = FALSE,
                              verbose = FALSE)

  # Recursive singleton peeling should remove the complete tail.
  expect_equal(nrow(result), nrow(core))
  expect_true(all(result$from %in% c(1, 2, 3)))
  expect_true(all(result$to %in% c(1, 2, 3)))
})

test_that("consolidate_graph summarizes singleton peeling logs", {
  graph <- data.frame(
    from = c(1, 2, 3, 4),
    to = c(2, 3, 4, 5),
    cost = 1
  )

  logs <- capture.output(
    consolidate_graph(graph, drop.edges = "single",
                      recursive = "partial", contract = FALSE,
                      verbose = TRUE)
  )
  drop_logs <- logs[grepl("Dropped .*edges leading to singleton nodes", logs)]

  # New implementation should print summarized peel logs, not one line per edge.
  expect_true(length(drop_logs) <= 2L)
  expect_true(any(grepl("peeling steps", drop_logs)))
})

test_that("consolidate_graph adds edge column only if present", {
  # Graph with branching (not just a chain)
  graph <- data.frame(
    from = c(1, 2, 2, 3, 4),
    to = c(2, 3, 4, 5, 5),
    cost = c(1, 2, 3, 1, 1)
  )

  result <- consolidate_graph(graph, keep.nodes = c(1, 5), verbose = FALSE)

  expect_false("edge" %in% names(result))
  expect_true(nrow(result) > 0)

  graph$edge <- seq_len(nrow(graph))
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
                             cost.column = "cost", verbose = FALSE)

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
                             cost.column = "cost", verbose = FALSE)

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
                             cost.column = "cost", verbose = FALSE)

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
                             cost.column = "cost", verbose = FALSE)

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
                             cost.column = "cost", verbose = FALSE)

  expect_true(nrow(result) > 0)
})

test_that("simplify_network shortest-paths nthreads matches single-thread output", {
  graph <- data.frame(
    from = c(1, 1, 2, 2, 3, 4, 4, 5),
    to = c(2, 3, 3, 4, 5, 5, 6, 6),
    cost = c(1, 2, 1, 2, 1, 1, 3, 1),
    mode = c("road", "rail", "road", "road", "rail", "road", "rail", "road")
  )
  od_pairs <- data.frame(
    from = c(1, 1, 2, 3),
    to = c(5, 6, 6, 6)
  )

  res_single <- simplify_network(graph, nodes = od_pairs,
                                 method = "shortest-paths",
                                 cost.column = "cost",
                                 by = ~ mode,
                                 nthreads = 1L, verbose = FALSE)
  res_multi <- simplify_network(graph, nodes = od_pairs,
                                method = "shortest-paths",
                                cost.column = "cost",
                                by = ~ mode,
                                nthreads = 2L, verbose = FALSE)

  expect_identical(res_multi, res_single)
  expect_identical(attr(res_multi, "edges"), attr(res_single, "edges"))
  expect_identical(attr(res_multi, "edge_counts"), attr(res_single, "edge_counts"))
})

test_that("simplify_network errors on missing columns", {
  graph <- data.frame(from = 1:3, cost = 1:3)

  expect_error(
    simplify_network(graph, nodes = c(1, 3), cost.column = "cost", verbose = FALSE),
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
    simplify_network(graph, nodes = c(1, 99), cost.column = "cost", verbose = FALSE),
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
                             radius_km = list(nodes = 50, cluster = 100), verbose = FALSE)

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
                             radius_km = list(nodes = 50, cluster = 100), verbose = FALSE)

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
                             radius_km = list(nodes = 50, cluster = 100), verbose = FALSE)

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
                     cost.column = "cost", verbose = FALSE),
    "FX.*FY.*TX.*TY"
  )
})
