context("Structure of panel information files")

# Preparation -------------------------------------------------------------
library(beadplexr)

info_files <- system.file(package = "beadplexr")
info_files <- list.files(file.path(info_files, "resources"), pattern = "*.yml$", full.names = TRUE)

test_file <- function(.file_path){
  # File can be loaded
  test_that(paste(basename(.file_path), "can be loaded"), {
    expect_is(panel_info <- load_panel(.file_name = .file_path), "list")
  })

  # Stil here?
  panel_info <- load_panel(.file_name = .file_path)

  # Slots are present
  test_that(paste(basename(.file_path), "slots present"), {
    slot_names <- names(panel_info)
    expect_equal(TRUE %in% grepl("panel_name", slot_names), TRUE)
    expect_equal(TRUE %in% grepl("cytokine_unit", slot_names), TRUE)
    expect_equal(TRUE %in% grepl("std_dilution", slot_names), TRUE)
    expect_equal(TRUE %in% grepl("analytes", slot_names), TRUE)
  })

  # Analyte can be converted to data frame
  test_that(paste(basename(.file_path), "as_data_frame_analyte"), {
    expect_is(as_data_frame_analyte(panel_info$analytes), "data.frame")
  })
}

# Run test ----------------------------------------------------------------
for(.file_path in info_files){
  test_file(.file_path)
}


