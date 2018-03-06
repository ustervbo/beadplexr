library(dplyr)
library(purrr)

# Use our own functions to read the fcs-files
source("./R/read_fcs.R")

# The fcs files
fcs_files <- list.files(path = "data-raw/raw_fcs/", pattern = "*.fcs", full.names = TRUE)
names(fcs_files) <-  basename(fcs_files)

# Load the data
lplex <- fcs_files %>% map(read_fcs)

devtools::use_data(lplex, compress = "xz", overwrite = TRUE)
