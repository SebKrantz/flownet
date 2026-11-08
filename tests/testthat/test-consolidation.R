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

  expect_true(!is.null(attr(result, "keep.edges")))
  expect_true(all(attr(result, "keep.edges") <= nrow(graph)))
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

  edge_counts <- attr(result, "edge.counts")
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
  expect_identical(attr(res_multi, "keep.edges"), attr(res_single, "keep.edges"))
  expect_identical(attr(res_multi, "edge.counts"), attr(res_single, "edge.counts"))
})

test_that("consolidate_graph leaves pure cycles intact (no crash/hang, C == R)", {
  # All-degree-2 cycles have no contractible chain head: contracting would
  # collapse to a self-loop. Both the C and R contraction paths must return the
  # cycle unchanged rather than crashing collap() or looping forever.
  triangle <- data.frame(from = c(1, 2, 3), to = c(2, 3, 1), cost = c(1, 1, 1))
  cycle4   <- data.frame(from = c(1, 2, 3, 4), to = c(2, 3, 4, 1), cost = c(1, 2, 3, 4))
  two_tri  <- data.frame(from = c(1, 2, 3, 4, 5, 6), to = c(2, 3, 1, 5, 6, 4), cost = 1:6)

  for (g in list(triangle, cycle4, two_tri)) {
    for (dir in c(TRUE, FALSE)) {
      expect_consolidate_equiv_c_r(g, directed = dir)
      res <- consolidate_graph(g, directed = dir, verbose = FALSE)
      expect_equal(nrow(res), nrow(g))  # cycle preserved, nothing contracted
    }
  }
})

test_that("consolidate_graph result attributes match the documented contract", {
  # Junction at node 3 so some edges survive consolidation.
  g <- data.frame(from = c(1, 2, 3, 3, 6, 3, 8), to = c(2, 3, 4, 6, 7, 8, 9),
                  cost = 1, edge = 1:7)

  # none / partial: keep.edges and group.id are attached together and aligned.
  for (rec in c("none", "partial")) {
    r <- consolidate_graph(g, recursive = rec, verbose = FALSE)
    ke <- attr(r, "keep.edges"); gid <- attr(r, "group.id")
    expect_false(is.null(ke))
    expect_false(is.null(gid))
    expect_identical(length(ke), length(gid))
    expect_true(all(gid >= 1L & gid <= nrow(r)))
  }

  # full (default): omitted by design (mapping would span multiple passes).
  rf <- consolidate_graph(g, recursive = "full", verbose = FALSE)
  expect_null(attr(rf, "keep.edges"))
  expect_null(attr(rf, "group.id"))

  # contract = FALSE (non-full) still exposes both, with identity group.id.
  rc <- consolidate_graph(g, contract = FALSE, recursive = "none", verbose = FALSE)
  expect_identical(attr(rc, "group.id"), seq_along(attr(rc, "keep.edges")))
})

test_that("create_undirected_graph attaches documented grouping attributes", {
  g <- data.frame(from = c(1, 2, 2, 1), to = c(2, 3, 1, 2),
                  edge = 1:4, FX = 0, FY = 0, TX = 1, TY = 1, cost = 1)
  r <- create_undirected_graph(g)
  expect_false(is.null(attr(r, "group.id")))
  expect_false(is.null(attr(r, "group.starts")))
  expect_false(is.null(attr(r, "group.sizes")))
  expect_identical(length(attr(r, "group.id")), nrow(g))      # one entry per input edge
  expect_identical(length(attr(r, "group.sizes")), nrow(r))   # one entry per output edge
})

test_that("simplify_network cluster method sets keep.edges (with and without self-loops)", {
  # Two distant clusters -> no self-loops after contraction.
  g_nosl <- data.frame(
    from = c(1, 2, 3, 4), to = c(2, 3, 4, 5),
    FX = c(0, 0, 0, 0), FY = c(0, 1, 2, 3),
    TX = c(0, 0, 0, 0), TY = c(1, 2, 3, 4), cost = 1
  )
  res <- simplify_network(g_nosl, nodes = 1:5, method = "cluster",
                          radius_km = list(nodes = 1e-4, cluster = 1e-4), verbose = FALSE)
  expect_false(is.null(attr(res, "keep.edges")))
  expect_identical(attr(res, "keep.edges"), seq_len(nrow(g_nosl)))
  expect_identical(length(attr(res, "keep.edges")), length(attr(res, "group.id")))

  # Force a self-loop: nodes 1 and 2 collapse into one cluster.
  g_sl <- data.frame(
    from = c(1, 2, 3), to = c(2, 3, 1),
    FX = c(0, 0, 5), FY = c(0, 0, 5),
    TX = c(0, 5, 0), TY = c(0, 5, 0), cost = 1
  )
  res2 <- simplify_network(g_sl, nodes = c(1, 3), method = "cluster",
                           radius_km = list(nodes = 1, cluster = 1), verbose = FALSE)
  expect_false(is.null(attr(res2, "keep.edges")))
  expect_identical(length(attr(res2, "keep.edges")), length(attr(res2, "group.id")))
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

test_that("simplify_network cluster with by returns each by column once", {
  graph <- linestrings_from_graph(africa_segments[1:100, ]) |>
    linestrings_to_graph()
  graph$mode <- rep_len(c("road", "rail"), nrow(graph))
  nodes_df <- nodes_from_graph(graph)
  keep_nodes <- nodes_df$node[1:5]

  result <- simplify_network(graph, nodes = keep_nodes,
                             method = "cluster", by = ~ mode,
                             cost.column = ".length",
                             radius_km = list(nodes = 50, cluster = 100), verbose = FALSE)

  expect_false(anyDuplicated(names(result)) > 0L)
  expect_equal(sum(names(result) == "mode"), 1L)
  # every output edge carries the mode of all the input edges merged into it
  ke <- attr(result, "keep.edges"); gid <- attr(result, "group.id")
  expect_true(all(graph$mode[ke] == result$mode[gid]))
  # the result is a valid graph data frame for downstream functions
  expect_no_error(create_undirected_graph(result, by = ~ mode))
})

test_that("simplify_network cluster adds group attributes", {
  testthat::skip_if_not_installed("RANN")
  graph <- linestrings_from_graph(africa_segments[1:100, ]) |>
    linestrings_to_graph()
  nodes_df <- nodes_from_graph(graph)
  keep_nodes <- nodes_df$node[1:5]

  result <- simplify_network(graph, nodes = keep_nodes,
                             method = "cluster",
                             cost.column = ".length",
                             radius_km = list(nodes = 50, cluster = 100),
                             nodes.algo.rann = TRUE, verbose = FALSE)

  expect_true(!is.null(attr(result, "group.id")))
  expect_true(!is.null(attr(result, "group.starts")))
})

test_that("RANN keep-node assignment matches dense assignment", {
  testthat::skip_if_not_installed("RANN")
  nodes <- data.frame(
    node = seq_len(10),
    X = c(0, 0.005, 0.02, 1, 1.004, 1.03, 10, 10.01, -20, -20.02),
    Y = c(0, 0.004, 0.01, 1, 1.003, 1.02, 0, 0.02, 5, 5.01)
  )
  keep <- c(1L, 4L, 9L)

  dense <- assign_nodes_to_keep(nodes, keep, radius_km = 3, verbose = FALSE)
  rann <- assign_nodes_to_keep_rann(nodes, keep, radius_km = 3,
                                    max_query_size = 3L, verbose = FALSE)

  expect_identical(rann$nodes, dense$nodes)
  expect_identical(rann$clusters, dense$clusters)
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

# Directed consolidation of two-way links given as reciprocal arcs

two_way <- function(from, to, cost = 1, back = cost) {
  data.frame(from = c(from, to), to = c(to, from), cost = c(rep_len(cost, length(from)), rep_len(back, length(from))))
}

drop_attrs <- function(x) {
  attr(x, "keep.edges") <- NULL
  attr(x, "group.id") <- NULL
  x
}

# Grid road network: most links two-way (as reciprocal arcs), some one-way, some removed.
# A share p_dual of the two-way links is a dual carriageway: one-way a -> m -> b and b -> m' -> a.
make_road_grid <- function(n, p_drop = 0.35, p_oneway = 0.25, p_rail = 0, p_dual = 0) {
  id <- matrix(seq_len(n * n), n)
  e <- rbind(cbind(c(id[-n, ]), c(id[-1, ])), cbind(c(id[, -n]), c(id[, -1])))
  e <- e[runif(nrow(e)) > p_drop, , drop = FALSE]
  ow <- runif(nrow(e)) < p_oneway
  rev <- ow & runif(nrow(e)) < 0.5
  e[rev, ] <- e[rev, 2:1]
  mode <- ifelse(runif(nrow(e)) < p_rail, "rail", "road")
  g <- rbind(data.frame(from = e[, 1], to = e[, 2], cost = runif(nrow(e), 1, 10), mode = mode),
             data.frame(from = e[!ow, 2], to = e[!ow, 1], cost = runif(sum(!ow), 1, 10), mode = mode[!ow]))
  if (p_dual <= 0) return(g)
  tw <- which(!ow)
  dl <- tw[runif(length(tw)) < p_dual]
  dual <- c(dl, nrow(e) + match(dl, tw))  # rows of both arcs of each dual link
  mid <- n * n + seq_along(dual)
  rbind(g[-dual, ], transform(g[dual, ], to = mid, cost = cost / 2), transform(g[dual, ], from = mid, cost = cost / 2))
}

keep_distances <- function(g, keep, nodes) {
  ig <- igraph::graph_from_data_frame(g[c("from", "to", "cost")], directed = TRUE,
                                      vertices = data.frame(name = as.character(nodes)))
  igraph::distances(ig, v = as.character(keep), to = as.character(keep), mode = "out",
                    weights = igraph::E(ig)$cost)
}

test_that("directed consolidate_graph contracts reciprocal arcs per direction", {
  # 1 <-> 2 <-> 3, with different costs in each direction
  g <- data.frame(from = c(1, 2, 3, 2), to = c(2, 3, 2, 1), cost = c(1, 2, 10, 20))
  res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 3),
                           custom = list(fsum = "cost"), verbose = FALSE)
  res <- res[order(res$from), ]
  expect_equal(res$from, c(1, 3))
  expect_equal(res$to, c(3, 1))
  expect_equal(res$cost, c(3, 30))

  # 1 <-> 2 <-> 3 <-> 4 <-> 5
  g5 <- two_way(1:4, 2:5, cost = 1:4, back = 5:8)
  res <- consolidate_graph(g5, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 5),
                           custom = list(fsum = "cost"), verbose = FALSE)
  res <- res[order(res$from), ]
  expect_equal(res$from, c(1, 5))
  expect_equal(res$to, c(5, 1))
  expect_equal(res$cost, c(10, 26))

  for (gr in list(g, g5))
    expect_consolidate_equiv_c_r(gr, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 5, 3)[seq_len(2 + identical(gr, g))],
                                 custom = list(fsum = "cost"))
})

test_that("directed consolidate_graph stops at a change between one-way and two-way", {
  # 1 <-> 2 <-> 3 -> 4 -> 5: node 3 can be left towards 2 or 4, so it stays
  g <- rbind(two_way(1:2, 2:3), data.frame(from = 3:4, to = 4:5, cost = 1))
  res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 5),
                           custom = list(fsum = "cost"), verbose = FALSE)
  res <- res[order(res$from, res$to), ]
  expect_equal(res$from, c(1, 3, 3))
  expect_equal(res$to, c(3, 1, 5))
  expect_equal(res$cost, c(2, 2, 2))
  expect_consolidate_equiv_c_r(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 5),
                               custom = list(fsum = "cost"))
})

test_that("directed consolidate_graph does not contract a crossing of one-way roads", {
  # 1 -> 5 -> 2 and 3 -> 5 -> 4: two arcs in and two out, but not reciprocal
  g <- data.frame(from = c(1, 5, 3, 5), to = c(5, 2, 5, 4), cost = 1)
  res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = 1:4, verbose = FALSE)
  expect_equal(nrow(res), 4L)
  # two in from and two out to the same neighbour (parallel arcs) is not a two-way link either
  gp <- data.frame(from = c(1, 1, 2, 2), to = c(2, 2, 1, 3), cost = 1)
  res <- consolidate_graph(gp, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 3), verbose = FALSE)
  expect_equal(nrow(res), 4L)
})

test_that("directed consolidate_graph respects by and keep.nodes at two-way nodes", {
  g <- cbind(two_way(1:2, 2:3), mode = c("road", "rail", "road", "rail"))
  res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, by = ~ mode,
                           keep.nodes = c(1, 3), verbose = FALSE)
  expect_equal(nrow(res), 4L)

  g4 <- two_way(1:3, 2:4)
  res <- consolidate_graph(g4, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 2, 4),
                           custom = list(fsum = "cost"), verbose = FALSE)
  res <- res[order(res$from, res$to), ]
  expect_equal(res$from, c(1, 2, 2, 4))
  expect_equal(res$to, c(2, 1, 4, 2))
  expect_equal(res$cost, c(1, 1, 2, 2))
})

test_that("directed consolidate_graph leaves two-way rings intact", {
  ring <- two_way(1:4, c(2:4, 1), cost = 1:4)
  for (dr in list(NULL, c("loop", "duplicate", "single"))) {
    for (rec in c("none", "partial", "full")) {
      res <- consolidate_graph(ring, directed = TRUE, drop.edges = dr, recursive = rec, verbose = FALSE)
      expect_equal(nrow(res), 8L)
      expect_true(all(c(res$from, res$to) %in% 1:4))
      expect_false(any(res$from == res$to))
    }
  }
  expect_consolidate_equiv_c_r(ring, directed = TRUE, drop.edges = NULL)
})

test_that("directed consolidate_graph keeps or peels two-way dead ends without creating loops", {
  # one-way trunk 1 -> 2 -> 3 with a two-way spur 2 <-> 8 <-> 9
  g <- rbind(data.frame(from = 1:2, to = 2:3, cost = 1), two_way(c(2, 8), c(8, 9)))
  res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 3),
                           custom = list(fsum = "cost"), verbose = FALSE)
  expect_false(any(res$from == res$to))
  res <- res[order(res$from, res$to), ]
  expect_equal(res$from, c(1, 2, 2, 9))
  expect_equal(res$to, c(2, 3, 9, 2))
  expect_equal(res$cost, c(1, 1, 2, 2))
  expect_consolidate_equiv_c_r(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 3),
                               custom = list(fsum = "cost"))

  for (rec in c("none", "partial", "full")) {
    res <- consolidate_graph(g, directed = TRUE, keep.nodes = c(1, 3), recursive = rec,
                             custom = list(fsum = "cost"), verbose = FALSE)
    if (rec == "none") {
      expect_false(any(c(res$from, res$to) == 9))  # outer end peeled once
    } else {
      expect_equal(drop_attrs(res), data.frame(from = 1, to = 3, cost = 2))
    }
  }
})

test_that("directed consolidate_graph partial recursion maps each arc to one output arc", {
  g <- rbind(two_way(c(1, 2, 3, 3), c(2, 3, 4, 6)), data.frame(from = c(4, 5), to = c(5, 3), cost = 1))
  g$cost <- seq_len(nrow(g))
  res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 6),
                           recursive = "partial", custom = list(fsum = "cost"), verbose = FALSE)
  ke <- attr(res, "keep.edges"); gid <- attr(res, "group.id")
  expect_identical(sort(ke), seq_len(nrow(g)))
  expect_true(nrow(res) < nrow(g))
  for (r in seq_len(nrow(res))) {
    m <- ke[gid == r]
    # the members form one directed path from res$from[r] to res$to[r]
    expect_equal(setdiff(g$from[m], g$to[m]), res$from[r])
    expect_equal(setdiff(g$to[m], g$from[m]), res$to[r])
    expect_equal(anyDuplicated(g$from[m]), 0L)
    expect_equal(sum(g$cost[m]), res$cost[r])
  }
  for (rec in c("none", "full"))
    expect_no_error(consolidate_graph(g, directed = TRUE, keep.nodes = c(1, 6), recursive = rec, verbose = FALSE))
  expect_consolidate_equiv_c_r(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 6),
                               recursive = "partial", custom = list(fsum = "cost"))
})

test_that("directed consolidate_graph contracts two-way nodes created by contracting one-way chains", {
  # Node 3 is reached from 1 and 5 by one-way chains and left towards both by one-way chains
  # (where the carriageways of a dual road meet). Contracting those chains turns 3 into a
  # two-way through node, whose arcs are then multi-row chains.
  g <- data.frame(from = c(1, 2, 3, 4, 5, 6, 3, 7),
                  to   = c(2, 3, 4, 1, 6, 3, 7, 5), cost = 1:8)
  for (rec in c("partial", "full")) {
    res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 5), recursive = rec,
                             custom = list(fsum = "cost"), verbose = FALSE)
    res <- drop_attrs(res)
    res <- res[order(res$from), ]
    rownames(res) <- NULL
    expect_equal(res, data.frame(from = c(1, 5), to = c(5, 1), cost = c(1 + 2 + 7 + 8, 5 + 6 + 3 + 4)))
    expect_consolidate_equiv_c_r(g, directed = TRUE, drop.edges = NULL, keep.nodes = c(1, 5), recursive = rec,
                                 custom = list(fsum = "cost"))
  }

  # Excerpt of a German road graph on which contracting such nodes one row per arc never terminated
  g <- data.frame(from = c(6, 7, 5, 14, 16, 13, 4, 8, 9, 3, 12, 15, 1, 13, 2, 11, 15, 7, 11, 6, 16, 1, 10),
                  to = c(7, 11, 11, 15, 6, 16, 13, 1, 1, 5, 14, 10, 4, 8, 3, 12, 2, 6, 7, 16, 13, 9, 15),
                  cost = c(12, 2, 15, 12, 49, 10, 15, 31, 18, 54, 61, 79, 31, 15, 6, 15, 12, 12, 2, 49, 10, 18, 79))
  dist <- function(x, v) keep_distances(x, v, 1:16)
  for (rec in c("partial", "full")) {
    res <- consolidate_graph(g, directed = TRUE, drop.edges = NULL, keep.nodes = 9, recursive = rec,
                             custom = list(fsum = "cost"), verbose = FALSE)
    nd <- sort(unique(c(res$from, res$to)))
    expect_equal(nd, c(9, 10))  # 10 is the end of a two-way dead end, kept without "single"
    expect_equal(dist(res, nd), dist(g, nd))
    expect_consolidate_equiv_c_r(g, directed = TRUE, drop.edges = NULL, keep.nodes = 9, recursive = rec,
                                 custom = list(fsum = "cost"))
  }
})

test_that("directed consolidate_graph preserves directed distances between keep nodes", {
  set.seed(42)
  for (i in 1:9) {
    g <- make_road_grid(12, p_rail = if (i %in% 4:6) 0.2 else 0, p_dual = if (i > 6) 0.3 else 0)
    nodes <- unique(c(g$from, g$to))
    keep <- sample(nodes, 8)
    d0 <- keep_distances(g, keep, nodes)
    for (rec in c("partial", "full")) {
      res <- consolidate_graph(g, directed = TRUE, drop.edges = c("loop", "single"), recursive = rec,
                               by = if (i %in% 4:6) ~ mode, keep.nodes = keep, custom = list(fsum = "cost"), verbose = FALSE)
      expect_true(nrow(res) < nrow(g))
      expect_equal(keep_distances(res, keep, nodes), d0, tolerance = 1e-9)
      # contraction only removes nodes, so distances between all remaining nodes are preserved too
      nd <- unique(c(res$from, res$to))
      expect_equal(keep_distances(res, nd, nodes), keep_distances(g, nd, nodes), tolerance = 1e-9)
    }
    expect_consolidate_equiv_c_r(g, directed = TRUE, drop.edges = c("loop", "single"),
                                 keep.nodes = keep, custom = list(fsum = "cost"))
  }
})

test_that("directed consolidation of an all-two-way graph gives two arcs per undirected edge", {
  set.seed(7)
  for (i in 1:4) {
    g <- make_road_grid(12, p_oneway = 0)
    und <- g[g$from < g$to, ]
    keep <- sample(unique(c(g$from, g$to)), 8)
    dr <- c("loop", "single")
    rd <- consolidate_graph(g, directed = TRUE, drop.edges = dr, keep.nodes = keep, verbose = FALSE)
    ru <- consolidate_graph(und, directed = FALSE, drop.edges = dr, keep.nodes = keep, verbose = FALSE)
    expect_equal(nrow(rd), 2L * nrow(ru))
    # Undirected mode only detects parallel edges with the same orientation as duplicates,
    # so with "duplicate" it is compared at its fixpoint over re-normalised edges.
    dr <- c("loop", "duplicate", "single")
    rd <- consolidate_graph(g, directed = TRUE, drop.edges = dr, keep.nodes = keep, verbose = FALSE)
    ru <- und; n0 <- Inf
    while (nrow(ru) < n0) {
      n0 <- nrow(ru)
      ru <- consolidate_graph(create_undirected_graph(ru), directed = FALSE, drop.edges = dr, keep.nodes = keep, verbose = FALSE)
    }
    expect_equal(nrow(rd), 2L * nrow(ru))
  }
})
