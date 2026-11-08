# Tests for run_assignment()

# Simple test graph from README
simple_graph <- data.frame(
 from = c(1, 2, 2, 3),
 to = c(2, 3, 4, 4),
 cost = c(5, 3, 2, 4)
)

simple_od <- data.frame(
 from = c(1, 2, 3),
 to = c(4, 4, 4),
 flow = c(100, 80, 60)
)

# --- AoN Method Tests ---

test_that("run_assignment AoN returns correct structure", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "AoN",
                          verbose = FALSE)

 expect_s3_class(result, "flownet")
 expect_true("final_flows" %in% names(result))
 expect_true("od_pairs_used" %in% names(result))
 expect_true("call" %in% names(result))
})

test_that("run_assignment AoN final_flows has correct length", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "AoN",
                          verbose = FALSE)

 expect_equal(length(result$final_flows), nrow(simple_graph))
})

test_that("run_assignment AoN return.extra='all' returns expected elements", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "AoN",
                          return.extra = "all",
                          verbose = FALSE)

 expect_true("graph" %in% names(result))
 expect_true("paths" %in% names(result))
 expect_true("path_costs" %in% names(result))
 expect_true("edge_counts" %in% names(result))
})

test_that("run_assignment AoN assigns flow to shortest paths", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "AoN",
                          verbose = FALSE)

 # All flow should be assigned (no flow lost)
 expect_true(sum(result$final_flows) > 0)

 # For AoN, flow goes only through shortest paths
 # Edge 3 (2->4, cost=2) should have high flow as it's on shortest paths
 expect_true(result$final_flows[3] > 0)
})

# --- PSL Method Tests ---

test_that("run_assignment PSL returns correct structure", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          verbose = FALSE)

 expect_s3_class(result, "flownet")
 expect_true("final_flows" %in% names(result))
 expect_equal(length(result$final_flows), nrow(simple_graph))
})

test_that("run_assignment PSL with angle.max=NA matches expected output", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          verbose = FALSE)

 # Expected values from README
 expected <- c(100.00000, 16.13649, 196.13649, 43.86351)
 expect_equal(result$final_flows, expected, tolerance = 1e-4)
})

test_that("run_assignment PSL return.extra='all' returns path weights", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = "all",
                          verbose = FALSE)

 expect_true("path_weights" %in% names(result))
 expect_true("paths" %in% names(result))
 expect_true("edges" %in% names(result))
})

test_that("run_assignment PSL produces more distributed flows than AoN", {
 result_psl <- run_assignment(simple_graph, simple_od,
                              cost.column = "cost",
                              method = "PSL",
                              angle.max = NA,
                              verbose = FALSE)

 result_aon <- run_assignment(simple_graph, simple_od,
                              cost.column = "cost",
                              method = "AoN",
                              verbose = FALSE)

 # PSL should have more non-zero edges or more even distribution
 # Check that PSL uses edge 4 (3->4) which AoN might not use as much
 expect_true(result_psl$final_flows[4] > 0)
})

test_that("run_assignment PSL beta parameter affects distribution", {
 result_beta0 <- run_assignment(simple_graph, simple_od,
                                cost.column = "cost",
                                method = "PSL",
                                beta = 0,
                                angle.max = NA,
                                verbose = FALSE)

 result_beta1 <- run_assignment(simple_graph, simple_od,
                                cost.column = "cost",
                                method = "PSL",
                                beta = 1,
                                angle.max = NA,
                                verbose = FALSE)

 # beta=0 (no overlap penalty) vs beta=1 should differ
 # At minimum, verify both run without error and produce valid flows
 expect_true(sum(result_beta0$final_flows) > 0)
 expect_true(sum(result_beta1$final_flows) > 0)
})

# --- Error Handling Tests ---

test_that("run_assignment errors on missing graph columns", {
 bad_graph <- data.frame(origin = c(1, 2), dest = c(2, 3), cost = c(1, 2))

 expect_error(
   run_assignment(bad_graph, simple_od, cost.column = "cost", verbose = FALSE)
 )
})

test_that("run_assignment errors on missing OD matrix columns", {
 bad_od <- data.frame(from = c(1, 2), to = c(3, 4))

 expect_error(
   run_assignment(simple_graph, bad_od, cost.column = "cost", verbose = FALSE),
   "flow"
 )
})

test_that("run_assignment errors on invalid cost.column", {
 expect_error(
   run_assignment(simple_graph, simple_od, cost.column = "nonexistent", verbose = FALSE),
   "cost.column"
 )
})

test_that("run_assignment errors on unknown nodes in OD matrix", {
 bad_od <- data.frame(from = c(1, 99), to = c(4, 4), flow = c(100, 50))

 expect_error(
   run_assignment(simple_graph, bad_od, cost.column = "cost", verbose = FALSE),
   "Unknown"
 )
})

test_that("run_assignment verbose output reports skipped OD pairs by cause", {
 g <- data.frame(from = c(1, 2, 3), to = c(2, 3, 4), cost = c(1, 2, 3))
 # Directed: pair 2 is unreachable, pair 4 has origin == destination, pairs 5-6 have zero/NA flow
 od <- data.frame(from = c(1, 4, 1, 2, 3, 2), to = c(4, 1, 3, 2, 4, 4), flow = c(10, 5, 2, 1, 0, NA))

 for (method in c("AoN", "PSL")) {
   expect_output(
     suppressWarnings(run_assignment(g, od, directed = TRUE, method = method, angle.max = NA)),
     "Skipped 4 OD-pairs: 2 with zero or non-finite flow and 2 with no path found", fixed = TRUE
   )
 }
 # Nothing skipped: no message
 out <- capture.output(invisible(run_assignment(g, od[c(1, 3), ], directed = TRUE, method = "AoN")))
 expect_false(any(grepl("Skipped", out, fixed = TRUE)))
})

# --- print.flownet Tests ---

test_that("print.flownet works without error", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "AoN",
                          return.extra = "all",
                          verbose = FALSE)

 expect_output(print(result), "FlowNet object")
})

# --- return.extra Options and PSF Tests ---

test_that("run_assignment PSL return.extra='all' returns PSF", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = "all",
                          verbose = FALSE)

 expect_true("path_size_factors" %in% names(result))
 expect_true("edge_weights" %in% names(result))

 # PSF should be list of same length as od_pairs_used
 expect_equal(length(result$path_size_factors), length(result$od_pairs_used))

 # Each PSF vector should have same length as corresponding path_weights
 for (i in seq_along(result$path_size_factors)) {
   expect_equal(length(result$path_size_factors[[i]]),
                length(result$path_weights[[i]]))
 }
})

test_that("PSF values are between 0 and 1", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = "PSF",
                          verbose = FALSE)

 for (psf in result$path_size_factors) {
   expect_true(all(psf > 0 & psf <= 1))
 }
})

test_that("PSF computation is mathematically correct", {
 # Use a slightly larger graph for meaningful overlap
 graph <- data.frame(
   from = c(1, 1, 2, 2, 3, 3),
   to = c(2, 3, 3, 4, 4, 4),
   cost = c(1, 2, 1, 3, 1, 2)
 )
 od <- data.frame(from = 1, to = 4, flow = 100)

 result <- run_assignment(graph, od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = c("paths", "edges", "counts", "costs", "PSF", "weights"),
                          verbose = FALSE)

 # PSF formula: gamma_k = sum(cost_a / delta_a) / cost_k
 # where delta_a is number of paths using edge a
 cost <- graph$cost
 paths <- result$paths[[1]]
 edges <- result$edges[[1]]
 ecounts <- result$edge_counts[[1]]
 pcosts <- result$path_costs[[1]]
 psf <- result$path_size_factors[[1]]

 # Manually compute PSF for each path
 for (k in seq_along(paths)) {
   path_edges <- paths[[k]]
   edge_costs <- cost[path_edges]
   edge_overlaps <- ecounts[match(path_edges, edges)]
   expected_psf <- sum(edge_costs / edge_overlaps) / pcosts[k]
   expect_equal(psf[k], expected_psf, tolerance = 1e-10)
 }
})

test_that("path_weights computation is mathematically correct (PSL formula)", {
 graph <- data.frame(
   from = c(1, 1, 2, 2, 3, 3),
   to = c(2, 3, 3, 4, 4, 4),
   cost = c(1, 2, 1, 3, 1, 2)
 )
 od <- data.frame(from = 1, to = 4, flow = 100)

 beta <- 1
 result <- run_assignment(graph, od,
                          cost.column = "cost",
                          method = "PSL",
                          beta = beta,
                          angle.max = NA,
                          return.extra = c("costs", "PSF", "weights"),
                          verbose = FALSE)

 pcosts <- result$path_costs[[1]]
 psf <- result$path_size_factors[[1]]
 pweights <- result$path_weights[[1]]

 # PSL formula: prob_k = exp(-cost_k + beta * log(PSF_k)) / sum(...)
 expected <- proportions(exp(-pcosts + beta * log(psf)))
 expect_equal(pweights, expected, tolerance = 1e-10)
})

test_that("edge_weights are sum of path probabilities per edge", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = c("paths", "edges", "weights", "eweights"),
                          verbose = FALSE)

 # For each OD pair, verify edge_weights
 for (i in seq_along(result$od_pairs_used)) {
   paths <- result$paths[[i]]
   edges <- result$edges[[i]]
   pweights <- result$path_weights[[i]]
   eweights <- result$edge_weights[[i]]

   # Manually compute edge weights
   expected_eweights <- numeric(length(edges))
   for (k in seq_along(paths)) {
     path_edges <- paths[[k]]
     for (e in path_edges) {
       idx <- match(e, edges)
       expected_eweights[idx] <- expected_eweights[idx] + pweights[k]
     }
   }

   expect_equal(eweights, expected_eweights, tolerance = 1e-10)
 }
})

test_that("path_costs are sum of edge costs along path", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = c("paths", "costs"),
                          verbose = FALSE)

 cost <- simple_graph$cost

 for (i in seq_along(result$od_pairs_used)) {
   paths <- result$paths[[i]]
   pcosts <- result$path_costs[[i]]

   for (k in seq_along(paths)) {
     expected_cost <- sum(cost[paths[[k]]])
     expect_equal(pcosts[k], expected_cost, tolerance = 1e-10)
   }
 }
})

test_that("edge_counts match number of paths using each edge", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = c("paths", "edges", "counts"),
                          verbose = FALSE)

 for (i in seq_along(result$od_pairs_used)) {
   paths <- result$paths[[i]]
   edges <- result$edges[[i]]
   ecounts <- result$edge_counts[[i]]

   # Manually count edge usage
   expected_counts <- integer(length(edges))
   for (k in seq_along(paths)) {
     for (e in paths[[k]]) {
       idx <- match(e, edges)
       expected_counts[idx] <- expected_counts[idx] + 1L
     }
   }

   expect_equal(ecounts, expected_counts)
 }
})

test_that("final_flows equals sum of flow * path_weight * edge membership", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = c("paths", "weights"),
                          verbose = FALSE)

 # Manually compute final_flows
 n_edges <- nrow(simple_graph)
 expected_flows <- numeric(n_edges)

 for (i in seq_along(result$od_pairs_used)) {
   flow_i <- simple_od$flow[result$od_pairs_used[i]]
   paths <- result$paths[[i]]
   pweights <- result$path_weights[[i]]

   for (k in seq_along(paths)) {
     for (e in paths[[k]]) {
       expected_flows[e] <- expected_flows[e] + flow_i * pweights[k]
     }
   }
 }

 expect_equal(result$final_flows, expected_flows, tolerance = 1e-10)
})

test_that("PSL probabilities sum to 1 for each OD pair", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "PSL",
                          angle.max = NA,
                          return.extra = "weights",
                          verbose = FALSE)

 for (pweights in result$path_weights) {
   expect_equal(sum(pweights), 1, tolerance = 1e-10)
 }
})

test_that("beta parameter affects PSF penalty correctly", {
 graph <- data.frame(
   from = c(1, 1, 2, 2, 3),
   to = c(2, 3, 3, 4, 4),
   cost = c(1, 2, 1, 2, 1)
 )
 od <- data.frame(from = 1, to = 4, flow = 100)

 result_b0 <- run_assignment(graph, od,
                             cost.column = "cost",
                             method = "PSL",
                             beta = 0,
                             angle.max = NA,
                             return.extra = c("costs", "PSF", "weights"),
                             verbose = FALSE)

 result_b5 <- run_assignment(graph, od,
                             cost.column = "cost",
                             method = "PSL",
                             beta = 5,
                             angle.max = NA,
                             return.extra = c("costs", "PSF", "weights"),
                             verbose = FALSE)

 # With beta=0, PSF should have no effect - pure cost-based logit
 pcosts <- result_b0$path_costs[[1]]
 expected_b0 <- proportions(exp(-pcosts))
 expect_equal(result_b0$path_weights[[1]], expected_b0, tolerance = 1e-10)

 # With higher beta, overlapping paths should be penalized more
 # PSF values should be the same, but weights different
 expect_equal(result_b0$path_size_factors[[1]],
              result_b5$path_size_factors[[1]], tolerance = 1e-10)

 # Verify beta=5 formula
 psf <- result_b5$path_size_factors[[1]]
 expected_b5 <- proportions(exp(-pcosts + 5 * log(psf)))
 expect_equal(result_b5$path_weights[[1]], expected_b5, tolerance = 1e-10)
})

test_that("PSL arguments passed as variables work with nthreads > 1", {
  # Regression test: beta, detour.max, unique.cost and npaths.max are only
  # referenced inside the parallel worker closure. Unless they are forced in
  # run_assignment(), mirai serializes them as unevaluated promises and the
  # daemons fail with "object '<name>' not found".
  graph <- data.frame(
    from = c(1, 1, 2, 2, 3, 4, 4, 5),
    to = c(2, 3, 3, 4, 5, 5, 6, 6),
    cost = c(1, 2, 1, 2, 1, 1, 3, 1)
  )
  od <- data.frame(from = c(1, 1, 2, 3), to = c(5, 6, 6, 6), flow = c(10, 20, 30, 40))

  beta_var <- 2
  detour_var <- 1.5
  unique_var <- TRUE
  npaths_var <- 50

  expect_no_error(
    res_multi <- run_assignment(graph, od, cost.column = "cost", method = "PSL",
                                beta = beta_var, detour.max = detour_var,
                                unique.cost = unique_var, npaths.max = npaths_var,
                                angle.max = NA, nthreads = 2L, verbose = FALSE)
  )

  res_single <- run_assignment(graph, od, cost.column = "cost", method = "PSL",
                               beta = 2, detour.max = 1.5,
                               unique.cost = TRUE, npaths.max = 50,
                               angle.max = NA, nthreads = 1L, verbose = FALSE)

  expect_equal(res_multi$final_flows, res_single$final_flows, tolerance = 1e-10)
  expect_identical(res_multi$od_pairs_used, res_single$od_pairs_used)
})

test_that("AoN paths match igraph::shortest_paths(), ties included, for any nthreads", {
  # The C search mirrors igraph's Dijkstra (heap, relaxation, incident-edge order), so ties
  # between equal-cost paths, zero-cost edges and parallel edges resolve as in igraph.
  igraph_epaths <- function(g, od, directed) {
    nodes <- sort(unique(c(g$from, g$to)))
    ig <- igraph::graph_from_data_frame(data.frame(from = match(g$from, nodes), to = match(g$to, nodes)),
                                        directed = directed, vertices = data.frame(name = seq_along(nodes)))
    o <- match(od$from, nodes); d <- match(od$to, nodes)
    out <- vector("list", nrow(od))
    for (s in unique(o)) {
      i <- which(o == s)
      sp <- suppressWarnings(igraph::shortest_paths(ig, from = s, to = d[i], weights = g$cost,
                                                    mode = "out", output = "epath")$epath)
      for (k in seq_along(i)) out[[i[k]]] <- as.integer(sp[[k]])
    }
    out
  }
  set.seed(101)
  for (it in 1:40) {
    nn <- sample(5:40, 1); E <- sample(nn:(3 * nn), 1)
    f <- sample.int(nn, E, TRUE); t <- sample.int(nn, E, TRUE); dup <- sample.int(E, E %/% 4, TRUE)
    g <- data.frame(from = c(f, f[dup]) * 3L, to = c(t, t[dup]) * 3L)
    g$cost <- sample(c(0, 1, 1, 2), nrow(g), TRUE)
    nds <- unique(c(g$from, g$to))
    od <- data.frame(from = sample(nds, 60, TRUE), to = sample(nds, 60, TRUE), flow = 1)
    for (directed in c(TRUE, FALSE)) {
      ref <- igraph_epaths(g, od, directed)
      for (nt in c(1L, 2L)) {
        r <- suppressWarnings(run_assignment(g, od, directed = directed, method = "AoN",
                                             return.extra = c("paths", "costs"), verbose = FALSE, nthreads = nt))
        used <- which(lengths(ref) > 0)
        expect_identical(as.integer(r$od_pairs_used), used)
        expect_identical(r$paths, ref[used])
        expect_identical(r$final_flows, tabulate(unlist(ref), nrow(g)) * 1)
        expect_equal(r$path_costs, vapply(ref[used], function(p) sum(g$cost[p]), 0))
      }
    }
  }
})
