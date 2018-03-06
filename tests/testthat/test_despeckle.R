context("Despeckle")

# Preparation -------------------------------------------------------------
library(beadplexr)
data("lplex")
.data <- lplex[[1]]
.nrow_data <- nrow(.data)

# Perform test ------------------------------------------------------------
test_that("Despeckle works", {
  # The call to raster::clump results in the following warnings:
  # no non-missing arguments to min; returning Inf
  # no non-missing arguments to max; returning -Inf
  # I don't know where exactly this originates from, but it related to the
  # raster::clump call. In the despecklefunction the warnings from the clump
  # call are suppressed
  expect_lt(nrow(despeckle(.data, .parameters = c("FL6-H", "FL2-H"))), .nrow_data)
  expect_lt(nrow(despeckle(.data, .parameters = c("FL6-H", "FL2-H"), .bins = 128L)), .nrow_data)
  expect_lt(nrow(despeckle(.data, .parameters = c("FL6-H", "FL2-H"), .neighbours = 8L)), .nrow_data)
  # If we have no points to remove, we should stay equal (the filter is neighbours < criteria)
  expect_equal(despeckle(.data, .parameters = c("FL6-H", "FL2-H"), .neighbours = 1L), .data)
})

test_that("Errors are thrown",{
  expect_error(despeckle(.data, .parameters = c("FL6-H")))
  expect_error(despeckle(.data, .parameters = c("xxxFL6-H", "FL2-H")))
})
