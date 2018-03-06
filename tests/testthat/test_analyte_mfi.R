context("Calculate analyte MFI")

# Preparation -------------------------------------------------------------
library(beadplexr)

data("lplex")
.data <- lplex[[1]]
.parameters <- c("FSC-A", "SSC-A")

.data <- subset(lplex[[1]], `FSC-A` > 4e5L & `FSC-A` < 6.3e5L)
.data <- identify_analyte(.data, .parameter = "FL6-H",
                   .analyte_id = as.character(c(1:7)))
.nrow_data <- nrow(.data)

# Analyte MFI -------------------------------------------------------------
test_that("Analyte MFI - 1 parameter", {
  expect_is(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "geometric"), "data.frame")
  expect_is(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "harmonic"), "data.frame")
  expect_is(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "arithmetic"), "data.frame")

  expect_lt(nrow(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "geometric")), .nrow_data)
  expect_lt(nrow(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "harmonic")), .nrow_data)
  expect_lt(nrow(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "arithmetic")), .nrow_data)

  expect_equal(ncol(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "geometric")), 2)
  expect_equal(ncol(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "harmonic")), 2)
  expect_equal(ncol(calc_analyte_mfi(.data, .parameter = "FL2-H", .mean_fun = "arithmetic")), 2)
})

test_that("Analyte MFI - 2 parameter", {
  expect_is(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "geometric"), "data.frame")
  expect_is(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "harmonic"), "data.frame")
  expect_is(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "arithmetic"), "data.frame")

  expect_lt(nrow(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "geometric")), .nrow_data)
  expect_lt(nrow(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "harmonic")), .nrow_data)
  expect_lt(nrow(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "arithmetic")), .nrow_data)

  expect_equal(ncol(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "geometric")), 3)
  expect_equal(ncol(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "harmonic")), 3)
  expect_equal(ncol(calc_analyte_mfi(.data, .parameter = c("FL2-H", "FL6-H"), .mean_fun = "arithmetic")), 3)
})

# Harmonic mean -----------------------------------------------------------
set.seed(12345)
x <- runif(10)

test_that("Harmonic mean is correct", {
  expect_equal(beadplexr:::harmonic_mean(x), 0.4855018, tolerance = 3e-8)
  expect_equal(beadplexr:::harmonic_mean(c(x, NA)), 0.4855018, tolerance = 3e-8)
  expect_equal(beadplexr:::harmonic_mean(c(x, Inf)), 0.534052, tolerance = 3e-8)
})

test_that("Harmonic mean fails", {
  expect_error(beadplexr:::harmonic_mean(c(x, "D")))
})

# Geometric mean ----------------------------------------------------------
set.seed(12345)
x <- runif(10)

test_that("Geometric mean is correct", {
  expect_equal(beadplexr:::geometric_mean(x), 0.5735927, tolerance = 3e-8)
  expect_equal(beadplexr:::geometric_mean(c(x, NA)), 0.5735927, tolerance = 3e-8)
  expect_equal(beadplexr:::geometric_mean(c(x, Inf)), Inf)
  expect_equal(beadplexr:::geometric_mean(c(x, -Inf)), 0.6033214, tolerance = 3e-8)
})

test_that("Harmonic mean fails", {
  expect_error(beadplexr:::geometric_mean(c(x, "D")))
})


