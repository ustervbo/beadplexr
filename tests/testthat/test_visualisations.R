context("Visualise data")

# Preparation -------------------------------------------------------------
.file_name <- system.file("extdata", "K2-C07-A7.fcs", package = "beadplexr")
.data <- read_fcs(.file_name = .file_name,
                  .filter = list("FSC-A" = c(2e5L, 6.3e5L),
                                 "SSC-A" = c(2e5, 1e6L)))
.data$bead_group <- ifelse(.data$`FSC-A` < 4e5L, "A", "B")

# FACS plots --------------------------------------------------------------

test_that("FACS data can be visualised", {
  # Test main wrapper function
  expect_is(facs_plot(.data, .type = "scatter"), "ggplot")
  expect_is(facs_plot(.data, .type = "density1d"), "ggplot")
  expect_is(facs_plot(.data, .type = "density2d"), "ggplot")
  expect_is(facs_plot(.data, .type = "hexbin"), "ggplot")
  expect_is(facs_plot(.data, .type = "scatter", .beads = "bead_group"), "ggplot")
  expect_is(facs_plot(.data, .type = "density1d", .beads = "bead_group"), "ggplot")
  expect_is(facs_plot(.data, .type = "density2d", .beads = "bead_group"), "ggplot")
  expect_is(facs_plot(.data, .type = "hexbin", .bins = 50), "ggplot")
  expect_is(facs_plot(.data, .x = "FL2-H", .type = "scatter", .beads = "bead_group"), "ggplot")
  # Test individual functions
  expect_is(facs_scatter(.data), "ggplot")
  expect_is(facs_scatter(.data, .beads = "bead_group", .plot_distinct = FALSE), "ggplot")
  expect_is(facs_scatter(.data, .beads = "bead_group"), "ggplot")
  expect_is(facs_scatter(.data, .x = "FL2-H", .y = "FL6-H", .beads = "bead_group"), "ggplot")
  expect_is(facs_density1d(.data), "ggplot")
  expect_is(facs_density1d(.data, .beads = "bead_group"), "ggplot")
  expect_is(facs_density2d(.data), "ggplot")
  expect_is(facs_density2d(.data, .beads = "bead_group"), "ggplot")
  expect_is(facs_hexbin(.data), "ggplot")
  expect_is(facs_hexbin(.data, .bins = 30), "ggplot")
})
