#' Generate an artificial bead population
#'
#' @param .range Named list of one or two elements (for one or two dimentions,
#'   respectively), giving the smallest and largest value.
#' @param .n_pop Numeric vector of 1 giving the number of populations to create
#'   within the range. Can only be \code{NULL} if \code{.pop_mean} is not.
#' @param .pop_mean List of means of each population. Each element is a
#'   population giving the coordinates for the mean. Can only be \code{NULL} if
#'   \code{.n_pop} is not.
#' @param .pop_sd Numeric vector giving the standard deviation of the populations. Is recycled if needed.
#' @param .pop_size Numeric vector giving the size of the populations. Is recycled if needed.
#' @param .pop_skew Numeric vector of the length of 1 giving the variance of
#'   skew to be applied. The actual skew is randomly shoen from a normal
#'   distribution with the given variance. \code{NA} is no skew.
#' @param noise Numeric vector giving the number of random events. Defaults to
#'   \code{.n_pop * .pop_size / 10} or \code{length(.pop_mean) * .pop_size / 10}
#' @param .min_dist Numeric giving the minimum distance between the points. Defaults to \code{.pop_sd * 2.5}.
#' @param .max_dist Numeric giving the maximum distance between the points. Defaults to \code{Inf}.
#'
#' @details
#'
#' \code{.range} gives the range of possible obervations. In case the
#' \code{.pop_mean} is not given, the means of the \code{.n_pop} are picked at
#' random within the range.
#'
#' @section Reason
#'
#' # CBA and MACSPlex data are different from LEGENXplex data in that they
#' utilize two dimentions to identify the beads. The package should be able to
#' handle all sorts of variations, so in lack of experimental data, we simulate
#' some.
#'
#' @return A \code{data.frame}
#' @export
#'

library(fGarch)
library(purrr)
library(dplyr)
library(ggplot2)
generate_bead_population <- function(.range,
           .n_pop = NULL,
           .pop_mean = NULL,
           .pop_sd = 0.5,
           .pop_size = 300,
           .pop_skew = NA,
           .noise = NULL,
           .min_dist = NULL,
           .max_dist = Inf) {

    # Argument check
    if(is.null(.n_pop) & is.null(.pop_mean)){
      stop("Both .n_pop and .pop_mean cannot be NULL")
    }
    if(!is.null(.n_pop) & !is.null(.pop_mean)){
      warning("Neiter .n_pop nor .pop_mean was NULL, which is confusing. I am just using .pop_mean")
      .n_pop <- NULL
    }

    if(length(.range) > 2){
      stop("I work within 1 or 2 dimensions")
    }

    if(is.null(.min_dist)){
      .min_dist <- .pop_sd * 5
    }

    # Adust the range
    .adjust_for_sd <- function(.range, .pop_sd){
      .range - .min_dist
    }
    .range %>% map(.adjust_for_sd, .pop_sd = .pop_sd)

    # Generate population centers
    if(!is.null(.n_pop)){
      .select_pop_mean <- function(.range){
        .range %>% map_dbl(~ runif(n = 1, min = .x[1], max = .x[2]))
      }

      # We need to find centers that are far enough from each other
      .pop_mean <- matrix(ncol = length(.range), nrow = .n_pop)

      # While we still have NA in the matrix we create random values and remove
      # them if they are too close to existing values
      while(TRUE %in% is.na(.pop_mean)){
        .free_row <- .pop_mean[, 1] %>% is.na() %>% which() %>% min
        .pop_mean[.free_row, ] <- .range %>% .select_pop_mean()

        .pop_dist <- dist(.pop_mean)
        if(FALSE %in% ((.pop_dist > .min_dist) & (.pop_dist < .max_dist))){
          .pop_mean[.free_row, ] <- rep(NA, ncol(.pop_mean))
        }
      }

      .pop_mean <- .pop_mean %>%
        t %>% data.frame %>% as.list %>%
        setNames(seq_len(.n_pop))

    }else{
      .n_pop <- length(.pop_mean)
    }

    # Calculate number noisy points
    if(is.null(.noise)){
      .noise <- (.n_pop * .pop_size / 10) %>% mean()
    }

    # Prepare to generate data
    .pop_sd <- rep_len(.pop_sd, .n_pop)
    .pop_size <- rep_len(.pop_size, .n_pop)

    # Generate the data. Requires fGarch to work.
    .generate_data <- function(.pop_mean, .pop_size, .pop_sd, .pop_skew){
      if(is.na(.pop_skew)){
        .pop_skew <- 1
      }else{
        # .pop_skew <- rnorm(1, mean = 0, sd = .pop_skew)
        .pop_skew <- .pop_skew
      }

      .pop_mean %>%
        map(~ rsnorm(n = .pop_size, mean = .x, sd = .pop_sd, xi = .pop_skew)) %>%
        do.call(cbind, .) %>% as.data.frame()
    }

    .generate_noise <- function(.range, .noise){
      runif(n = .noise, min = .range[1], max = .range[2])
    }

    population <- pmap_df(list(.pop_mean, .pop_size, .pop_sd, .pop_skew), .generate_data) %>%
      setNames(names(.range))

    noise <- .range %>%
      map(.generate_noise, .noise = .noise) %>%
      do.call(cbind, .) %>% as.data.frame() %>%
      setNames(names(.range))
    bind_rows(population, noise)
}

set.seed(12345)
# Create simulated LegendPlex set -----------------------------------------
.n_pop <- 6
.pop_size <- 300
bead_intesity_a <- generate_bead_population(
    .range = list(PE = c(5, 40)),
    .n_pop = .n_pop,
    .pop_sd = 0.15,
    .pop_skew = 1.6,
    .pop_size = .pop_size
  )

bead_intesity_b <- generate_bead_population(
  .range = list(PE = c(5, 20)),
  .n_pop = .n_pop,
  .pop_sd = 0.15,
  .pop_skew = 1,
  .pop_size = .pop_size
)

bead_intesity <- bind_rows(bead_intesity_a, bead_intesity_b)

bead_id_a <- generate_bead_population(
  .range = list(APC = c(5, 30)),
  .n_pop = .n_pop,
  .pop_sd = 0.5,
  .pop_skew = 1.6,
  .pop_size = .pop_size
)

bead_id_b <- generate_bead_population(
  .range = list(APC = c(5, 25)),
  .n_pop = .n_pop,
  .pop_sd = 0.5,
  .pop_skew = 1,
  .pop_size = .pop_size
)

bead_id <- bind_rows(bead_id_a, bead_id_b)

bead_events <- bind_cols(bead_intesity, bead_id)

fsc_ssc <- generate_bead_population(
  .range = list(FSC = c(5, 25), SSC = c(5, 25)),
  .n_pop = 2,
  .pop_sd = 2,
  .pop_skew = 1.6,
  .pop_size = nrow(bead_events)/2,
  .min_dist = 1, .max_dist = 9.5
) %>%
  sample_n(nrow(bead_events))

legendplex_example <- bead_events %>% bind_cols(fsc_ssc)

legendplex_example %>%
  ggplot() +
  aes(x = FSC, y = SSC) +
  geom_point()

bead_events %>%
  ggplot() +
  aes(x = PE, y = APC) +
  geom_point()

# saveRDS(legendplex_example, "data-raw/sim_lplex.rds")

# Create simulated CBA set -------------------------------------------------
.n_pop <- 30
.pop_size <- 300
bead_intesity <- generate_bead_population(
    .range = list(PE = c(5, 40)),
    .n_pop = .n_pop,
    .pop_sd = 0.5,
    .pop_skew = 1.6,
    .pop_size = .pop_size,
    .min_dist = 0
  )

bead_id <- generate_bead_population(
  .range = list(APC = c(5, 30), `APC-Cy7` = c(5, 25)),
  .n_pop = .n_pop,
  .pop_sd = 0.25,
  .pop_skew = 1.6,
  .pop_size = .pop_size,
  .min_dist = 3
)

bead_events <- bind_cols(bead_intesity, bead_id)

fsc_ssc <- generate_bead_population(
  .range = list(FSC = c(5, 25), SSC = c(5, 25)),
  .pop_mean = list(c(15, 15)),
  .pop_sd = 2,
  .pop_skew = 1.6,
  .pop_size = nrow(bead_events)
) %>%
  sample_n(nrow(bead_events))

cba_example <- bead_events %>% bind_cols(fsc_ssc)

cba_example %>%
  ggplot() +
  aes(x = FSC, y = SSC) +
  geom_point()

cba_example %>%
  ggplot() +
  aes(x = APC, y = `APC-Cy7`) +
  geom_point()

cba_example %>%
  ggplot() +
  aes(x = PE) +
  geom_density()

# saveRDS(cba_example, "data-raw/sim_cba.rds")

# Create simulated MACSPlex set -------------------------------------------------
.n_pop <- 10
.pop_size <- 300
bead_intesity <- generate_bead_population(
  .range = list(APC = c(5, 40)),
  .n_pop = .n_pop,
  .pop_sd = 0.5,
  .pop_skew = 1.6,
  .pop_size = .pop_size,
  .min_dist = 0
)

bead_id <- generate_bead_population(
  .range = list(PE = c(5, 15), FITC = c(5, 15)),
  .n_pop = .n_pop,
  .pop_sd = 0.25,
  .pop_skew = 1.6,
  .pop_size = .pop_size,
  .min_dist = 2
)

bead_events <- bind_cols(bead_intesity, bead_id)

fsc_ssc <- generate_bead_population(
  .range = list(FSC = c(5, 25), SSC = c(5, 25)),
  .pop_mean = list(c(15, 15)),
  .pop_sd = 2,
  .pop_skew = 1.6,
  .pop_size = nrow(bead_events)
) %>%
  sample_n(nrow(bead_events))

mplex_example <- bead_events %>% bind_cols(fsc_ssc)

mplex_example %>%
  ggplot() +
  aes(x = FSC, y = SSC) +
  geom_point()

mplex_example %>%
  ggplot() +
  aes(x = PE, y = FITC) +
  geom_point()

mplex_example %>%
  ggplot() +
  aes(x = APC) +
  geom_density()

# saveRDS(mplex_example, "data-raw/sim_mplex.rds")


# Combine and save simulated example --------------------------------------
# lplex <- readRDS("data-raw/sim_lplex.rds")
# mplex <- readRDS("data-raw/sim_mplex.rds")
# cba <- readRDS("data-raw/sim_cba.rds")

simplex <- list(lplex = legendplex_example, mplex = mplex_example, cba = cba_example)
devtools::use_data(simplex, compress = "xz", overwrite = TRUE)
# saveRDS(sim_data, "data-raw/sim_data_example.rds")
