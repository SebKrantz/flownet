# Tests for network processing functions

# --- linestrings_to_graph() Tests ---

test_that("linestrings_to_graph returns data.frame with required columns", {
  # Use a small subset of africa_network
  lines <- africa_network[1:10, ]

  result <- linestrings_to_graph(lines)

  expect_true(is.data.frame(result))
  expect_true(all(c("edge", "from", "to", "FX", "FY", "TX", "TY") %in% names(result)))
})
  
test_that("linestrings_to_graph computes .length when compute.length=TRUE", {
  lines <- africa_network[1:10, ]

  result <- linestrings_to_graph(lines, compute.length = TRUE)

  expect_true(".length" %in% names(result))
  # .length is a units object, convert to numeric for comparison
  expect_true(all(as.numeric(result$.length) > 0))
})

test_that("linestrings_to_graph omits .length when compute.length=FALSE", {
  lines <- africa_network[1:10, ]

  result <- linestrings_to_graph(lines, compute.length = FALSE)

  expect_false(".length" %in% names(result))
})

test_that("linestrings_to_graph preserves columns with keep.cols", {
  lines <- africa_network[1:10, ]

  result <- linestrings_to_graph(lines, keep.cols = is.atomic)

  # Should preserve atomic columns from original
  expect_true("distance" %in% names(result) || "duration" %in% names(result))
})

test_that("linestrings_to_graph creates consistent node IDs", {
  lines <- africa_network[1:10, ]

  result <- linestrings_to_graph(lines)

  # All from and to should be positive integers
  expect_true(all(result$from > 0))
  expect_true(all(result$to > 0))
  expect_true(all(result$from == as.integer(result$from)))
  expect_true(all(result$to == as.integer(result$to)))
})

test_that("linestrings_to_graph errors on non-LINESTRING input", {
  # Create a non-LINESTRING sf object (POINT)
  points <- sf::st_as_sf(data.frame(x = 1:3, y = 1:3), coords = c("x", "y"))

  expect_error(linestrings_to_graph(points), "LINESTRING")
})

# --- linestrings_from_graph() Tests ---

test_that("linestrings_from_graph returns sf object with LINESTRING", {
  graph_df <- data.frame(
    FX = c(0, 1),
    FY = c(0, 0),
    TX = c(1, 2),
    TY = c(0, 0),
    cost = c(1, 2)
  )

  result <- linestrings_from_graph(graph_df)

  expect_s3_class(result, "sf")
  expect_equal(as.character(sf::st_geometry_type(result, by_geometry = FALSE)), "LINESTRING")
})

test_that("linestrings_from_graph uses correct CRS", {
  graph_df <- data.frame(
    FX = c(0, 1),
    FY = c(0, 0),
    TX = c(1, 2),
    TY = c(0, 0),
    id = c(1, 2)  # Need a non-coordinate column
  )

  result <- linestrings_from_graph(graph_df, crs = 4326)
  expect_equal(sf::st_crs(result)$epsg, 4326)

  result2 <- linestrings_from_graph(graph_df, crs = 3857)
  expect_equal(sf::st_crs(result2)$epsg, 3857)
})

test_that("linestrings_from_graph preserves other columns", {
  graph_df <- data.frame(
    FX = c(0, 1),
    FY = c(0, 0),
    TX = c(1, 2),
    TY = c(0, 0),
    cost = c(1, 2),
    name = c("a", "b")
  )

  result <- linestrings_from_graph(graph_df)

  expect_true("cost" %in% names(result))
  expect_true("name" %in% names(result))
  expect_equal(result$cost, c(1, 2))
})

test_that("linestrings_from_graph removes coordinate columns", {
  graph_df <- data.frame(
    FX = c(0, 1),
    FY = c(0, 0),
    TX = c(1, 2),
    TY = c(0, 0),
    id = c(1, 2)  # Need a non-coordinate column
  )

  result <- linestrings_from_graph(graph_df)

  expect_false("FX" %in% names(result))
  expect_false("FY" %in% names(result))
  expect_false("TX" %in% names(result))
  expect_false("TY" %in% names(result))
  expect_true("id" %in% names(result))  # Other columns preserved
})

test_that("linestrings_from_graph errors on non-data.frame", {
  expect_error(linestrings_from_graph("not a data frame"), "data frame")
})

test_that("linestrings_from_graph errors on sf input", {
  lines <- africa_network[1:5, ]
  expect_error(linestrings_from_graph(lines), "spatial object")
})

test_that("linestrings_from_graph errors on missing coordinate columns", {
  graph_df <- data.frame(from = 1:2, to = 2:3)
  expect_error(linestrings_from_graph(graph_df), "FX.*FY.*TX.*TY")
})

# --- Roundtrip Test ---

test_that("linestrings roundtrip preserves structure", {
  # Start with segments data (has coordinates)
  graph_df <- africa_segments[1:20, ]

  # Convert to sf and back
  sf_obj <- linestrings_from_graph(graph_df)
  result <- linestrings_to_graph(sf_obj, keep.cols = is.atomic)

  # Should have same number of rows
 expect_equal(nrow(result), nrow(graph_df))

  # Should have from/to columns
  expect_true(all(c("from", "to") %in% names(result)))
})

# --- create_undirected_graph() Tests ---

test_that("create_undirected_graph reduces edge count", {
  # Create a directed graph with bidirectional edges
  graph <- data.frame(
    from = c(1, 2, 2, 3),
    to = c(2, 1, 3, 2),
    cost = c(1, 1, 2, 2)
  )

  result <- create_undirected_graph(graph)

  # Should have fewer edges after removing directional duplicates
  expect_lt(nrow(result), nrow(graph))
})

test_that("create_undirected_graph ensures from < to", {
  graph <- data.frame(
    from = c(5, 3, 4),
    to = c(3, 5, 2),
    cost = c(1, 1, 2)
  )

  result <- create_undirected_graph(graph)

  expect_true(all(result$from < result$to))
})

test_that("create_undirected_graph aggregates with FUN", {
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 1),
    cost = c(10, 20)
  )

  result_mean <- create_undirected_graph(graph, FUN = "fmean")
  result_sum <- create_undirected_graph(graph, FUN = "fsum")

  # fmean should give 15, fsum should give 30
  expect_equal(result_mean$cost, 15)
  expect_equal(result_sum$cost, 30)
})

test_that("create_undirected_graph preserves coordinate columns", {
  graph <- data.frame(
    from = c(1, 2),
    to = c(2, 1),
    FX = c(0, 1),
    FY = c(0, 0),
    TX = c(1, 0),
    TY = c(0, 0),
    cost = c(1, 1)
  )

  result <- create_undirected_graph(graph)

  expect_true(all(c("FX", "FY", "TX", "TY") %in% names(result)))
})

test_that("create_undirected_graph with by preserves groups", {
  graph <- data.frame(
    from = c(1, 2, 1, 2),
    to = c(2, 1, 2, 1),
    mode = c("road", "road", "rail", "rail"),
    cost = c(1, 1, 2, 2)
  )

  result <- create_undirected_graph(graph, by = ~ mode)

  # Should have 2 edges (one per mode)
  expect_equal(nrow(result), 2)
  expect_true(all(c("road", "rail") %in% result$mode))
})

test_that("create_undirected_graph sets group.starts attribute", {
  graph <- data.frame(
    from = c(1, 2, 2, 3),
    to = c(2, 1, 3, 2),
    cost = c(1, 1, 2, 2)
  )

  result <- create_undirected_graph(graph)

  expect_true(!is.null(attr(result, "group.starts")))
})
