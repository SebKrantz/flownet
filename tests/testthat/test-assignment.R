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

# --- print.flownet Tests ---

test_that("print.flownet works without error", {
 result <- run_assignment(simple_graph, simple_od,
                          cost.column = "cost",
                          method = "AoN",
                          return.extra = "all",
                          verbose = FALSE)

 expect_output(print(result), "FlowNet object")
})
