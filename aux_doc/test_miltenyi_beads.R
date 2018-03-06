# Test Miltenyi data. Kindly provided by Miltenyi
#
# The fcs files are in data-raw/_miltenyi_test_fcs


# Libraries ---------------------------------------------------------------
library(beadplexr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# Functions ---------------------------------------------------------------
plot_side_by_side <- function(plot_data, cur_sample){

  plot_all_beads <- plot_data %>%
    facs_plot(.x = "FSC-A", .y = "SSC-A", .type = "hex") +
    ggtitle("All events")

  plot_beads <- plot_data %>%
    facs_plot(.x = "FL2-A", .y = "FL3-A", .beads = "analyte")

  plot_pe <- plot_data %>%
    filter(!is.na(analyte)) %>%
    facs_density1d(.x = "FL3-A", .beads = "analyte") +
    ggtitle("APC")

  arrangeGrob(plot_all_beads, plot_beads, plot_pe,
              nrow = 1, ncol = 3, top = cur_sample)
}

# Load data ---------------------------------------------------------------
fcs_files <- list.files(path = "./data-raw/_miltenyi_test_fcs", pattern = "*.fcs", full.names = TRUE)
names(fcs_files) <-  basename(fcs_files)

comp_matrix <- matrix(c(1, 0.5, 0.24, 0.24, 1, 0, 0, 0, 1),
       nrow = 3, ncol = 3, byrow = TRUE,
       dimnames = list(c("FL2-A", "FL3-A", "FL6-A"), c("FL2-A", "FL3-A", "FL6-A")))

bead_data <- fcs_files %>% map(read_fcs, .fsc_ssc = c("FSC-A", "SSC-A"),
                                          .bead_channels = c("FL2-A", "FL3-A", "FL6-A"),
                                          .filter = list("FSC-A" = c(0L, 50L), "SSC-A"= c(0L, 100L),
                                                         "FL2-A" = c(1.75, Inf),
                                                         # "FL3-A" = c(2, Inf),
                                                         "FL6-A" = c(0, Inf)),
                                          .compensation = comp_matrix,
                                          emptyValue = FALSE)


bead_analyte <- bead_data %>%
  lapply(identify_analyte, .parameter = c("FL2-A", "FL3-A"), .analyte_id = LETTERS[1:11], .trim = 0.1)

bead_analyte_mfi <- bead_analyte %>%
  map_df(calc_analyte_mfi,
         .parameter = "FL6-A",
         .column_name = "analyte",
         .mean_fun = "geometric",
         .id = "Sample") %>%
  filter(!is.na(`analyte`))

# Plots -------------------------------------------------------------------
all_plots <- bead_analyte %>%
  map2(names(.), plot_side_by_side) %>%
  marrangeGrob(ncol = 1, nrow = 4, top = NA)

ggsave(filename = "miltinyi_dot_plot.pdf", plot = all_plots, width = 8.27, height = 11.69, path = "aux_doc/")

bead_analyte_mfi %>%
  mutate(Sample = str_replace(Sample, ".*_", "")) %>%
  mutate(Sample = str_replace_all(Sample, "\\.fcs|00", "")) %>%
  # mutate(Sample = as.integer(Sample)) %>%
  mutate(Sample = factor(Sample, levels = c(1:8))) %>%
  ggplot() +
  aes(x = Sample, y = `FL6-A`, colour = analyte, group = analyte) +
  geom_point() +
  geom_line() +
  labs(x = "Sample number", y = "APC")
