context("Cluster functions")

exclude_dbscan <- TRUE

# Preparation -------------------------------------------------------------
library(testthat)
library(beadplexr)
data("lplex")
.data <- lplex[[1]]
.parameters <- c("FSC-A", "SSC-A")

# Test bp_kmeans ----------------------------------------------------------
test_that("bp_kmeans works", {
  # One cluster
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "clu ster", .k = 1), "data.frame")
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "clu ster", .k = 1, .trim = 0.2), "data.frame")
  # Two clusters
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2), "data.frame")
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2, .trim = 0.2), "data.frame")
  # Pass redudant parameter to kmeans
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 1,  centers = 3), "data.frame")
  # Additional parameters to kmeans
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  nstart = 3), "data.frame")
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  algorithm = "Lloyd"), "data.frame")
  expect_output(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  trace = TRUE))
  # One dimentional
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters[1], .column_name = "cluster", .k = 1, .trim = 0.2), "data.frame")
  expect_is(bp_kmeans(.data = .data, .parameter = .parameters[2], .column_name = "cluster", .k = 2, .trim = 0.2), "data.frame")
})

test_that("bp_kmeans fails", {
  expect_error(bp_kmeans(.data = .data, .parameter = "xxx", .column_name = "cluster", .k = 2, .trim = 0.2))
  expect_error(bp_kmeans(.data = .data, .parameter = .parameters, .column_name = "cluster", centers = 2, .trim = 0.2))
})

# Test bp_clara -----------------------------------------------------------
test_that("bp_clara works", {
  # One cluster
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "clus ter", .k = 1), "data.frame")
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "clu ster", .k = 1, .trim = 0.2), "data.frame")
  # Two clusters
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2), "data.frame")
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2, .trim = 0.2), "data.frame")
  # Pass redudant parameter to clara
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 1,  k = 3), "data.frame")
  # Additional parameters to clara
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  samples = 50), "data.frame")
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  metric = "manhattan"), "data.frame")
  expect_is(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  pamLike = FALSE), "data.frame")
  # One dimentional
  expect_is(bp_clara(.data = .data, .parameter = .parameters[1], .column_name = "cluster", .k = 1, .trim = 0.2), "data.frame")
  expect_is(bp_clara(.data = .data, .parameter = .parameters[2], .column_name = "cluster", .k = 2, .trim = 0.2), "data.frame")
})

test_that("bp_clara fails", {
  expect_error(bp_clara(.data = .data, .parameter = "xxx", .column_name = "cluster", .k = 2, .trim = 0.2))
  expect_error(bp_clara(.data = .data, .parameter = .parameters, .column_name = "cluster", k = 2, .trim = 0.2))
})

# Test bp_dbscan ----------------------------------------------------------
test_that("bp_dbscan works", {
  # Defaults
  expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "clu ster"), "data.frame")
  if(exclude_dbscan){
    skip("Indepth test of bp_dbscan")
    # Set eps and MinPts
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", .eps = 0.1), "data.frame")
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", .MinPts = 100), "data.frame")
    # Pass redudant parameter to dbscan
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", .eps = 0.2, eps = 0.1), "data.frame")
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", .MinPts = 50, MinPts = 50), "data.frame")
    # Additional parameters to dbscan
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", method = "raw"), "data.frame")
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", seeds = FALSE), "data.frame")
    expect_is(bp_dbscan(.data = .data, .parameter = .parameters, .column_name = "cluster", scale = TRUE), "data.frame")
  }
})

test_that("bp_dbscan fails", {
  expect_error(bp_dbscan(.data = .data, .parameter = "xxx", .column_name = "cluster"))
})


# Test bp_mclust ----------------------------------------------------------
test_that("bp_mclust works", {
  # One cluster
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "clus ter", .k = 1), "data.frame")
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "clu ster", .k = 1, .trim = 0.2), "data.frame")
  # Two clusters
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2), "data.frame")
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2, .trim = 0.2), "data.frame")
  # Pass redudant parameter to clara
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 1,  k = 3), "data.frame")
  # Additional parameters to clara
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  samples = 50), "data.frame")
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  metric = "manhattan"), "data.frame")
  expect_is(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", .k = 2,  pamLike = FALSE), "data.frame")
  # One dimentional
  expect_is(bp_mclust(.data = .data, .parameter = .parameters[1], .column_name = "cluster", .k = 1, .trim = 0.2), "data.frame")
  expect_is(bp_mclust(.data = .data, .parameter = .parameters[2], .column_name = "cluster", .k = 2, .trim = 0.2), "data.frame")
})

test_that("bp_mclust fails", {
  expect_error(bp_mclust(.data = .data, .parameter = "xxx", .column_name = "cluster", .k = 2, .trim = 0.2))
  expect_error(bp_mclust(.data = .data, .parameter = .parameters, .column_name = "cluster", k = 2, .trim = 0.2))
})
