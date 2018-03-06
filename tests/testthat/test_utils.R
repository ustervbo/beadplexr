context("Panel info list to data.frame")

# Preparation -------------------------------------------------------------
library(beadplexr)
panel_info <- load_panel(.panel_name = "Human Growth Factor Panel (13-plex)")

# Perform test ------------------------------------------------------------
test_that("Conversion of list", {
  .analytes <-  list(A = list(
    A1 = list(name = "name_a1", concentration = 500),
    A2 = list(name = "name_a2", concentration = 50000)))

  expect_is(as_data_frame_analyte(.analytes), "data.frame")
  expect_is(as_data_frame_analyte(panel_info$analytes), "data.frame")
})

