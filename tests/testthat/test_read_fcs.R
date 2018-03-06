context("Read and transform fcs file")

# Preparation -------------------------------------------------------------
.file_name <- system.file("extdata", "K2-C07-A7.fcs", package = "beadplexr")
.flow_frame <-flowCore::read.FCS(filename = .file_name, transformation = FALSE)


# Boundary events can be removed ------------------------------------------
test_that("Boundary events can be removed", {
  expect_is(remove_boundary_events(.flow_frame = .flow_frame, .channels = c("FSC-A", "SSC-A")), "flowFrame")
  expect_is(remove_boundary_events(.flow_frame = .flow_frame, .channels = c("FSC-A")), "flowFrame")
  # Must have parameters
  # expect_error(remove_boundary_events(.flow_frame = .flow_frame, .channels = NULL))
})


# Channels can be filtered ------------------------------------------------
test_that("Channels can be filtered", {
  expect_is(subset_channels(.flow_frame = .flow_frame, .filter = list("FSC-A" = c(1e5L, 8e5L))), "flowFrame")
  expect_is(subset_channels(.flow_frame = .flow_frame, .filter = list("FSC-A" = c(-Inf, 8e5L))), "flowFrame")
  expect_is(subset_channels(.flow_frame = .flow_frame, .filter = NULL), "flowFrame")
  # Each list element must be a numeric vector of two
  expect_error(subset_channels(.flow_frame = .flow_frame, .filter = list("FSC-A" = c(8e5L))))
  expect_error(subset_channels(.flow_frame = .flow_frame, .filter = list("FSC-A" = c("A", "B"))))
})

# Transformation can be applied -------------------------------------------
test_that("Transformation can be applied", {
  expect_is(transform_bead_channels(.flow_frame = .flow_frame, .bead_channels = c("FL6-H", "FL2-H")), "flowFrame")
  # The channel must exist
  expect_error(transform_bead_channels(.flow_frame = .flow_frame, .bead_channels = c("FL6-H", "FL23-H")))
  # Beads are identified by intensity in one channel, and the relative level of the product is given by another channel
  expect_error(transform_bead_channels(.flow_frame = .flow_frame, .bead_channels = c("FL6-H")))
})

# Expression set is exportet as data.frame --------------------------------
test_that("Expression set is exportet as data.frame", {
  expect_is(as_data_frame_flow_frame(.flow_frame = .flow_frame), "data.frame")
  expect_is(as_data_frame_flow_frame(.flow_frame = .flow_frame), "tbl")
  expect_is(as_data_frame_flow_frame(.flow_frame = .flow_frame, .channels = c("FSC-A", "SSC-A")), "tbl")
  # Parameters must exist
  expect_error(as_data_frame_flow_frame(.flow_frame = .flow_frame, .channels = c(" ")))
})


# Compensation application ------------------------------------------------
odd_ff <- .flow_frame
flowCore::keyword(odd_ff) <- list("SPILL" = flowCore::keyword(odd_ff, "$SPILLOVER")[[1]])

comp_matrix <- diag(1, nrow = 6, ncol = 6)
colnames(comp_matrix) <- c("FL1-H", "FL2-H", "FL6-H", "FL1-A", "FL2-A", "FL6-A")

test_that("Compensation can be applied", {
  expect_is(apply_compensation(.flow_frame = .flow_frame, .compensation = NULL), "flowFrame")
  expect_is(apply_compensation(.flow_frame = .flow_frame, .compensation = "guess"), "flowFrame")
  expect_is(apply_compensation(.flow_frame = .flow_frame, .compensation = "SPILL"), "flowFrame")
  expect_is(apply_compensation(.flow_frame = .flow_frame, .compensation = "spillover"), "flowFrame")
  expect_is(apply_compensation(.flow_frame = .flow_frame, .compensation = comp_matrix), "flowFrame")

  # Warnings
  expect_warning(apply_compensation(.flow_frame = .flow_frame, .compensation = c("guess", "xxx")))
  expect_warning(apply_compensation(.flow_frame = .flow_frame, .compensation = "xxx"))
  expect_warning(apply_compensation(.flow_frame = odd_ff, .compensation = "guess"))
  expect_warning(apply_compensation(.flow_frame = odd_ff, .compensation = "spill"))

  # Expect error
  expect_error(apply_compensation(.flow_frame = odd_ff, .compensation = list(list("guess", "XXX"))))
  expect_error(apply_compensation(.flow_frame = odd_ff, .compensation = data.frame("guess", "XXX")))
})
# Clean
rm(odd_ff, comp_matrix)

# FACS files can be read and processed ------------------------------------
test_that("FACS files can be read and processed", {
  expect_is(read_fcs(.file_name = .file_name), "data.frame")
  expect_is(read_fcs(.file_name = .file_name, .filter = list("FSC-A" =  c(Inf, 5))), "data.frame")
  expect_is(read_fcs(.file_name = .file_name, .fsc_ssc = c("FSC-H", "SSC-H")), "data.frame")
  expect_is(read_fcs(.file_name = .file_name, .bead_channels = c("FSC-H", "SSC-H")), "data.frame")
  # Parameters must exist
  expect_error(read_fcs(.file_name = .file_name, .fsc_ssc = c("FSC-H", " ")))
  # Two forward side scatter channels
  expect_error(read_fcs(.file_name = .file_name, .fsc_ssc = c("FSC-H")))
  # filter must be list
  expect_error(read_fcs(.file_name = .file_name, .filter = c("FSC-H")))
})
