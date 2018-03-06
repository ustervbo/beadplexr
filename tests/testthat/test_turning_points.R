context("Turning points")

# Preparation -------------------------------------------------------------
set.seed(1234)
.x2 <- c(rnorm(100, 2, 1), rnorm(100, 9, 1))
.x3 <- c(rnorm(100, 2, 1), rnorm(100, 3, 1), rnorm(100, 4, 1))
.x4 <- c(rnorm(100, 2, 1), rnorm(100, 3, 1), rnorm(100, 4, 1), rnorm(100, 6, 1))

# Test do_find_turning_points ---------------------------------------------
test_that("do_find_turning_points", {
  expect_is(beadplexr:::do_find_turning_points(.x2), "list")
  expect_is(beadplexr:::do_find_turning_points(.x2, .return = "index"), "list")

  expect_equal(length(beadplexr:::do_find_turning_points(.x2)[[1]]), 2)
  expect_equal(length(beadplexr:::do_find_turning_points(.x2)[[2]]), 1)
  expect_equal(length(beadplexr:::do_find_turning_points(.x2, adjust = 0.2)[[1]]), 5)
  expect_equal(length(beadplexr:::do_find_turning_points(.x2, adjust = 0.2)[[2]]), 4)
})

# Test approx_adjust ------------------------------------------------------

test_that("approx_adjust works", {
  # Gives 0.5
  expect_equal(approx_adjust(.x2, .k = 2), 0.4)
  # Equals 0.42
  expect_equal(approx_adjust(.x = .x2,
                             .k = 3, .lower = 0.2,
                             .upper = 1,
                             .step = 0.001), 0.281)
  # Gives warning and return NA
  expect_warning(ret_val <- approx_adjust(.x2,
                             .k = 3,
                             .lower = 0.1,
                             .upper = 1,
                             .step = 0.1))
  expect_true(is.na(ret_val))

  expect_warning(ret_val <- approx_adjust(.x2, .k = 3, .lower = -1))
  expect_true(is.na(ret_val))
})

# Test turning_point ------------------------------------------------------
test_that("turning_point works", {
  expect_is(turning_point(.x2, .which = "min", .return = "ind"), "list")
  expect_is(turning_point(.x2, .which = "max", .return = "ind"), "list")
  expect_is(turning_point(.x2, .which = "both", .return = "ind"), "list")

  expect_is(turning_point(.x2, .which = "min", .return = "value"), "list")
  expect_is(turning_point(.x2, .which = "max", .return = "value"), "list")
  expect_is(turning_point(.x2, .which = "both", .return = "value"), "list")

  expect_is(turning_point(.x2, .which = "min", .return = "value")[[1]], "data.frame")
  expect_is(turning_point(.x2, .which = "max", .return = "value")[[1]], "data.frame")
  expect_is(turning_point(.x2, .which = "both", .return = "value")[[1]], "data.frame")

  expect_is(turning_point(list(A = .x2, B = .x2), .which = "both", .return = "value"), "list")
  expect_is(turning_point(list(A = .x2, B = .x2), .which = "both", .return = "value")[[1]], "data.frame")
})


test_that("turning_points warn because of different maxima in the two parameters", {
  # warning comes from cbind
  expect_warning(turning_point(.x = list(A = .x3, B = .x4), .which = "min", .return = "ind", .adjust = 0.8))
  expect_warning(turning_point(.x = list(A = .x3, B = .x4), .which = "max", .return = "ind", .adjust = 0.8))
  expect_warning(turning_point(.x = list(A = .x3, B = .x4), .which = "both", .return = "ind", .adjust = 0.8))
})

test_that("Adjust can be estimated", {
  expect_message(turning_point(.x2, .which = "both", .return = "ind", .k = 2))
  expect_message(turning_point(list(A = .x2, B = .x2), .which = "both", .return = "ind", .k = 2))

  expect_error(suppressWarnings(turning_point(.x = .x2, .which = "both", .return = "value", .k = 3)))
  expect_is(ret_val <- suppressMessages(
    turning_point(list(.x3), .which = "both", .return = "value", .k = 4)
  ), "list")
})
