context("Cluster manipulation")

# Preparation -------------------------------------------------------------
library(beadplexr)
data("lplex")
.data <- lplex[[1]]
.nrow_data <- nrow(.data)
.parameters <- c("FSC-A", "SSC-A")

# Calculate distance ------------------------------------------------------
test_that("Distance to center is correct", {
  expect_equal(beadplexr:::calc_dist_to_centre(.x = 5, .c = 9), 4)
  expect_equal(beadplexr:::calc_dist_to_centre(.x = 9, .c = 9), 0)

  expect_equal(beadplexr:::calc_dist_to_centre(.x = c(5, 5), .c = c(9, 9)), 5.66, tolerance = .004)
  expect_equal(beadplexr:::calc_dist_to_centre(.x = c(5,9), .c = c(5,9)), 0)

  expect_equal(beadplexr:::calc_dist_to_centre(.x = c(5, 4, 6), .c = c(9, 8, 7)), 5.74, tolerance = .01)
  expect_equal(beadplexr:::calc_dist_to_centre(.x = c(5,9, 3), .c = c(5,9, 3)), 0)
  })

test_that("Calculation fails", {
  expect_error(beadplexr:::calc_dist_to_centre(.x = "A", .c = 9))
  expect_warning(beadplexr:::calc_dist_to_centre(.x = c(5, 4, 6), .c = c(9, 8)))
})

# Trim clusters -----------------------------------------------------------
.data$population <- 1L
test_that("Clusters are trimmed", {
  expect_is(trim_population(.data, .parameter = .parameters, .column_name = "population", .trim = 0.1), "data.frame")
  expect_equal(nrow(trim_population(.data, .parameter = .parameters, .column_name = "population", .trim = 0.1)), .nrow_data)

  expect_is(trim_population(.data, .parameter = .parameters, .column_name = "population", .trim = 0), "data.frame")
  expect_equal(nrow(trim_population(.data, .parameter = .parameters, .column_name = "population", .trim = 0)), .nrow_data)

  expect_is(trim_population(.data, .parameter = .parameters[1], .column_name = "population", .trim = 0), "data.frame")
  expect_equal(nrow(trim_population(.data, .parameter = .parameters[2], .column_name = "population", .trim = 0)), .nrow_data)
})

test_that("Clusters fail", {
  expect_error(trim_population(.data, .parameter = .parameters, .column_name = "xxx", .trim = 0.1))
  expect_error(trim_population(.data, .parameter = "xxx", .column_name = "population", .trim = 0.1))
})
