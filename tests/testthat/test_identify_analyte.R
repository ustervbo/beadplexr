context("Identify analytes")

# Preparation -------------------------------------------------------------
library(beadplexr)
data("lplex")
.data <- lplex[[1]]
.parameters <- c("FSC-A", "SSC-A")

# identify_analyte() ---------------------------------------------------
test_that("Identify analyte works", {
  expect_is(identify_analyte(.data = .data,
                                .parameter = c("FSC-A", "SSC-A"),
                                .analyte_id = c("A", "B"),
                                .column_name = "ana lyte",
                                .method = "clara", .trim = 0.02), "data.frame")
  expect_is(identify_analyte(.data = .data,
                                .parameter = c("FSC-A", "SSC-A"),
                                .analyte_id = c("A", "B"),
                                .column_name = "analyte",
                                .method = "kmeans", .trim = 0.02), "data.frame")
  expect_is(identify_analyte(.data = .data,
                                .parameter = c("FSC-A", "SSC-A"),
                                .analyte_id = c("A", "B"),
                                .column_name = "analyte",
                                .method = "dbscan"), "data.frame")
})

test_that("Identify analyte give warnings", {
  expect_warning(identify_analyte(.data = .data,
                      .parameter = c("FSC-A", "SSC-A"),
                      .analyte_id = c("A", "B"),
                      .column_name = "analyte",
                      .method = "clara", .k = 4))
})


# assign_analyte_id() -----------------------------------------------------

.data <- bp_clara(.data, .parameter = c("FSC-A", "SSC-A"), .column_name = "analyte", .k = 2)

test_that("AnalyteIDs are assigned", {
  expect_is(beadplexr:::assign_analyte_id(.data = .data,
                                          .parameter = c("FSC-A", "SSC-A"),
                                          .analyte_id = c("A", "B"),
                                          .column_name = "pop name",
                                          .cluster_column_name = "analyte"), "data.frame")
  expect_is(beadplexr:::assign_analyte_id(.data = .data,
                                          .parameter = c("FSC-A", "SSC-A"),
                                          .analyte_id = c("A", "B"),
                                          .column_name = "pop name",
                                          .cluster_column_name = "analyte", .desc = TRUE), "data.frame")
  expect_true("pop name" %in% names(beadplexr:::assign_analyte_id(.data = .data,
                                                      .parameter = c("FSC-A", "SSC-A"),
                                                      .analyte_id = c("A", "B"),
                                                      .column_name = "pop name",
                                                      .cluster_column_name = "analyte", .desc = TRUE)))
  expect_false("analyte" %in% names(beadplexr:::assign_analyte_id(.data = .data,
                                                                  .parameter = c("FSC-A", "SSC-A"),
                                                                  .analyte_id = c("A", "B"),
                                                                  .column_name = "pop name",
                                                                  .cluster_column_name = "analyte", .desc = TRUE)))
})

test_that("AnalyteID column is overwritten", {
  .data2 <- beadplexr:::assign_analyte_id(.data = .data,
                                          .parameter = c("FSC-A", "SSC-A"),
                                          .analyte_id = c("A", "B"),
                                          .column_name = "pop name",
                                          .cluster_column_name = "analyte")

  .data3 <- bp_clara(.data2, .parameter = c("FSC-A", "SSC-A"), .column_name = "analyte", .k = 2)

  expect_is(beadplexr:::assign_analyte_id(.data = .data3,
                                          .parameter = c("FSC-A", "SSC-A"),
                                          .analyte_id = c("A", "B"),
                                          .column_name = "pop name",
                                          .cluster_column_name = "analyte"), "data.frame")
  expect_identical(beadplexr:::assign_analyte_id(.data = .data3,
                                                 .parameter = c("FSC-A", "SSC-A"),
                                                 .analyte_id = c("A", "B"),
                                                 .column_name = "pop name",
                                                 .cluster_column_name = "analyte"), .data2)
})
