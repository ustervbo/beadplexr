# Test CBA data. Kindly provided by Colin McGuckin, CTI Biotech
#
# The fcs files are in data-raw/_cba_test_fcs


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
    facs_plot(.x = "CBA Red-H", .y = "CBA NIR-H", .beads = "analyte")

  plot_pe <- plot_data %>%
    filter(!is.na(analyte)) %>%
    facs_density1d(.x = "PE-H", .beads = "analyte") +
    ggtitle("PE")

  arrangeGrob(plot_all_beads, plot_beads, plot_pe,
              nrow = 1, ncol = 3, top = cur_sample)
}

# Load data ---------------------------------------------------------------
fcs_files <- list.files(path = "./data-raw/_cba_test_fcs", pattern = "*.fcs", full.names = TRUE)
names(fcs_files) <-  basename(fcs_files)

cba_data <- fcs_files %>% map(read_fcs, .fsc_ssc = c("FSC-A", "SSC-A"),
                              .bead_channels = c("CBA Red-H", "CBA NIR-H", "PE-H"),
                              .filter = list("FSC-A" = c(2e3L, 5e4L), "SSC-A"= c(2e3L, 5e4L),
                                             "CBA Red-H" = c(8.5, Inf), "CBA NIR-H" = c(7.5, Inf)))


cba_analyte <- cba_data %>%
  lapply(identify_analyte, .parameter = c("CBA Red-H", "CBA NIR-H"), .analyte_id = LETTERS[1:5], .trim = 0.1)

cba_analyte_mfi <- cba_analyte %>%
  map_df(calc_analyte_mfi,
         .parameter = "PE-H",
         .column_name = "analyte",
         .mean_fun = "geometric",
         .id = "Sample") %>%
  filter(!is.na(`analyte`))

# Plots -------------------------------------------------------------------
all_plots <- cba_analyte %>%
  map2(names(.), plot_side_by_side) %>%
  marrangeGrob(ncol = 1, nrow = 4, top = NA)

ggsave(filename = "cba_dot_plot.pdf", plot = all_plots, width = 8.27, height = 11.69, path = "./aux_doc/")

cba_analyte_mfi %>%
  mutate(Sample = str_replace(Sample, "_.*", "")) %>%
  mutate(Sample = str_replace(Sample, ".*\\s", "")) %>%
  # mutate(Sample = as.integer(Sample)) %>%
  mutate(Sample = factor(Sample, levels = c(1:10))) %>%
  ggplot() +
  aes(x = Sample, y = `PE-H`, colour = analyte, group = analyte) +
  geom_point() +
  geom_line() +
  labs(x = "Sample number", y = "PE")
