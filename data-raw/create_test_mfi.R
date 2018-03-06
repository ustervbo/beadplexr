library(beadplexr)
data("lplex")

panel_info <- load_panel(.panel_name = "Human Growth Factor Panel (13-plex)")

args_ident_analyte <- list(fs = list(.parameter = c("FSC-A", "SSC-A"),
                                     .column_name = "Bead group",
                                     .trim = 0.03),
                           analytes = list(.parameter = "FL6-H",
                                           .column_name = "Analyte ID"))


find_and_trim <- function(.data){
  identify_legendplex_analyte(.data, .analytes = panel_info$analytes,
                              .method_args = args_ident_analyte) %>%
    group_by(`Analyte ID`) %>%
    do(trim_population(., .parameter = c("FL6-H", "FL2-H"), .column_name = "Analyte ID", .trim = 0.1))
}

analyte_mfi <- lplex %>% lapply(find_and_trim) %>%
  map_df(calculate_analyte_mfi,
         .parameter = "FL2-H",
         .column_name = "Analyte ID",
         .mean_fun = "geometric",
         .id = "Sample") %>%
  filter(!is.na(`Analyte ID`))

saveRDS(analyte_mfi,file = "./tests/testthat/testdata/analyte_mfi.rds")
