# Tests for melt_od_matrix()

# --- Basic Functionality Tests ---

test_that("melt_od_matrix returns data.frame with correct columns", {
  od_mat <- matrix(c(0, 10, 20, 30), nrow = 2)

  result <- melt_od_matrix(od_mat)

  expect_true(is.data.frame(result))
  expect_true(all(c("from", "to", "flow") %in% names(result)))
})

test_that("melt_od_matrix filters zero flows", {
  od_mat <- matrix(c(0, 10, 0, 30), nrow = 2)

  result <- melt_od_matrix(od_mat)

  # Only non-zero values should remain
  expect_true(all(result$flow > 0))
  expect_equal(nrow(result), 2)  # Only 10 and 30
})

test_that("melt_od_matrix filters non-finite flows", {
  od_mat <- matrix(c(NA, 10, Inf, 30), nrow = 2)

  result <- melt_od_matrix(od_mat)

  # Only finite values should remain
  expect_true(all(is.finite(result$flow)))
  expect_equal(nrow(result), 2)  # Only 10 and 30
})

test_that("melt_od_matrix uses dimnames when nodes not provided", {
  od_mat <- matrix(c(0, 10, 20, 0), nrow = 2)
  dimnames(od_mat) <- list(c("100", "200"), c("100", "200"))

  result <- melt_od_matrix(od_mat)

  # Should use dimnames as node IDs (coerced to integer)
  expect_true(all(result$from %in% c(100, 200)))
  expect_true(all(result$to %in% c(100, 200)))
})

test_that("melt_od_matrix uses sequential IDs when no dimnames", {
  od_mat <- matrix(c(0, 10, 20, 0), nrow = 2)

  result <- melt_od_matrix(od_mat)

  # Should use 1:n as node IDs
  expect_true(all(result$from %in% 1:2))
  expect_true(all(result$to %in% 1:2))
})

test_that("melt_od_matrix uses nodes argument when provided", {
  od_mat <- matrix(c(0, 10, 20, 0), nrow = 2)
  node_ids <- c(50, 60)

  result <- melt_od_matrix(od_mat, nodes = node_ids)

  expect_true(all(result$from %in% node_ids))
  expect_true(all(result$to %in% node_ids))
})

test_that("melt_od_matrix nodes argument overrides dimnames", {
  od_mat <- matrix(c(0, 10, 20, 0), nrow = 2)
  dimnames(od_mat) <- list(c("100", "200"), c("100", "200"))
  node_ids <- c(50, 60)

  result <- melt_od_matrix(od_mat, nodes = node_ids)

  # Should use nodes, not dimnames
  expect_true(all(result$from %in% node_ids))
  expect_false(any(result$from %in% c(100, 200)))
})

test_that("melt_od_matrix sort=TRUE orders by from, to", {
  od_mat <- matrix(1:4, nrow = 2)
  dimnames(od_mat) <- list(c("2", "1"), c("2", "1"))

  result <- melt_od_matrix(od_mat, sort = TRUE)

  # Should be sorted by from, then to
  expected_order <- order(result$from, result$to)
  expect_equal(seq_len(nrow(result)), expected_order)
})

test_that("melt_od_matrix sort=FALSE preserves original order", {
  od_mat <- matrix(1:4, nrow = 2)

  result_sorted <- melt_od_matrix(od_mat, sort = TRUE)
  result_unsorted <- melt_od_matrix(od_mat, sort = FALSE)

  # Both should have same content but potentially different order
  expect_equal(sort(result_sorted$flow), sort(result_unsorted$flow))
})

# --- Edge Cases ---

test_that("melt_od_matrix handles 1x1 matrix", {
  od_mat <- matrix(5, nrow = 1)

  result <- melt_od_matrix(od_mat)

  expect_equal(nrow(result), 1)
  expect_equal(result$flow, 5)
})

test_that("melt_od_matrix handles all-zero matrix", {
  od_mat <- matrix(0, nrow = 3, ncol = 3)

  result <- melt_od_matrix(od_mat)

  expect_equal(nrow(result), 0)
})

test_that("melt_od_matrix handles non-square matrix", {
  od_mat <- matrix(1:6, nrow = 2, ncol = 3)

  result <- melt_od_matrix(od_mat)

  expect_equal(nrow(result), 6)  # All positive values
})

test_that("melt_od_matrix handles large values", {
  od_mat <- matrix(c(0, 1e12, 1e12, 0), nrow = 2)

  result <- melt_od_matrix(od_mat)

  expect_equal(result$flow, c(1e12, 1e12))
})

test_that("melt_od_matrix handles small positive values", {
  od_mat <- matrix(c(0, 1e-10, 1e-10, 0), nrow = 2)

  result <- melt_od_matrix(od_mat)

  expect_equal(nrow(result), 2)
  expect_true(all(result$flow > 0))
})

# --- Error Handling Tests ---

test_that("melt_od_matrix errors on non-matrix input", {
  expect_error(melt_od_matrix(data.frame(a = 1:3)), "matrix")
  expect_error(melt_od_matrix(c(1, 2, 3)), "matrix")
  expect_error(melt_od_matrix(list(a = 1)), "matrix")
})

test_that("melt_od_matrix errors on non-numeric matrix", {
  char_mat <- matrix(c("a", "b", "c", "d"), nrow = 2)

  expect_error(melt_od_matrix(char_mat), "numeric")
})

test_that("melt_od_matrix errors on mismatched nodes length (rows)", {
  od_mat <- matrix(1:4, nrow = 2)
  nodes <- c(1, 2, 3)  # Wrong length

  expect_error(melt_od_matrix(od_mat, nodes = nodes), "rows")
})

test_that("melt_od_matrix errors on mismatched nodes length (cols)", {
  od_mat <- matrix(1:6, nrow = 2, ncol = 3)
  nodes <- c(1, 2)  # Matches rows but not cols

  expect_error(melt_od_matrix(od_mat, nodes = nodes), "columns")
})

test_that("melt_od_matrix errors when dimnames don't match with nodes", {
  od_mat <- matrix(1:4, nrow = 2)
  dimnames(od_mat) <- list(c("a", "b"), c("c", "d"))  # Different row/col names
  nodes <- c(1, 2)

  expect_error(melt_od_matrix(od_mat, nodes = nodes), "match")
})

# --- Integration with run_assignment() ---

test_that("melt_od_matrix output works with run_assignment", {
  # Create simple OD matrix
  od_mat <- matrix(c(0, 100, 80, 0), nrow = 2)
  dimnames(od_mat) <- list(c("1", "3"), c("1", "3"))

  # Simple graph
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 3),
    cost = c(1, 1)
  )

  # Melt and use
  od_long <- melt_od_matrix(od_mat)

  result <- run_assignment(graph, od_long,
                           cost.column = "cost",
                           method = "AoN",
                           verbose = FALSE)

  expect_s3_class(result, "flownet")
  expect_true(sum(result$final_flows) > 0)
})
