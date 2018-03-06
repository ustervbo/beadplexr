#' Calculate the analyte intensity
#'
#' The mean fluorescence intensity (MFI) of the analyte
#'
#' @inheritParams cluster_events
#' @param .mean_fun A character giving the mean function to use.
#'
#' @section Mean calculation:
#'
#' @inheritSection harmonic_mean Mean calculation
#'
#' @inheritSection geometric_mean Mean calculation
#'
#' @return A summarised `data.frame`
#' @export
#'
#' @examples
#' library(beadplexr)
#' library(dplyr)
#'
#' data("lplex")
#'
#' .data <- lplex[[1]] %>%
#'   filter(`FSC-A` > 4e5L, `FSC-A` < 6.3e5L) %>%
#'   identify_analyte(.parameter = "FL6-H",
#'                    .analyte_id = as.character(c(1:7)))
#'
#' .data %>%
#'   calc_analyte_mfi(.parameter = "FL2-H")
#'
#' .data %>%
#'   calc_analyte_mfi(.parameter = "FL2-H",
#'               .mean_fun = "harmonic")
calc_analyte_mfi <- function(.data,
                              .parameter,
                              .column_name = "analyte",
                              .mean_fun = c("geometric", "harmonic", "arithmetic")){
  .mean_fun <- match.arg(.mean_fun)

  .mean_fun <- switch(.mean_fun,
    geometric = geometric_mean,
    harmonic = harmonic_mean,
    arithmetic = mean
  )

  .create_summarise_fun <- function(.parameter){
    lazyeval::interp(~.mean_fun(.the_parameter), .the_parameter = as.name(.parameter))
  }

  .summarise_fun <- lapply(.parameter, .create_summarise_fun) %>% stats::setNames(.parameter)

  .data %>%
    dplyr::group_by_(as.name(.column_name)) %>%
    dplyr::summarise_(.dots = .summarise_fun)
}

#' Calculate harmonic mean
#'
#' @param .x A numeric to take the harmonic mean of
#'
#' @section Mean calculation:
#' The calculation of the harmonic mean is as follows:
#' \deqn{
#'  n\frac{1}{\sum_{i=1}^n \frac{1}{x_i}}
#' }{length(x)/sum(1/x)}
#'
#' `NA`s are removed before calculation
#'
#' @return A single numeric value
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' harmonic_mean(runif(10))
#' }
harmonic_mean <- function(.x){
  .x <- .x[!is.na(.x)]
  length(.x)/sum(1/.x)
}

#' Calculate geometric mean
#'
#' @param .x A numeric to take the geometric mean of
#'
#' @section Mean calculation:
#' The geometric mean is given by:
#'
#' \deqn{
#'   \left(\prod_{i=1}^n x_i \right)^\frac{1}{n}
#' }{prod(x)^(length(x))}
#'
#'
#' but implemented as:
#'
#' \deqn{
#'   \frac{1}{n}\exp{\sum_{i=1}^n\log(x_i)}
#' }{exp(sum(log(x[x > 0])) / length(x))}
#'
#' `NA`s are removed before calculation
#'
#'
#' @return A single numeric value
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' geometric_mean(runif(10))
#' }
geometric_mean <- function(.x){
  .x <- .x[!is.na(.x)]
  exp(sum(log(.x[.x > 0])) / length(.x))
}
