#' Calculate standard concentration
#'
#' Given a start concentration and dilution factor, the concentration of the given standard samples is calculated
#'
#' @param .standard_sample a vector giving the standard samples. The sample with
#'   the highest value is given the start concentration, and a
#'   `.standard_sample` with the value of 0, is set to 0 if it exists. See
#'   Details for details on how order is assessed.
#' @param .start_concentration a numeric vector giving the initial standard
#'   concentration. If longer than one the maximum value is taken as start
#'   concentration.
#' @param .dilution_factor a numeric vector giving the dilution factor. If a
#'   single element is passed, this is applied to all standard samples as a
#'   dilution series. If more then one value is given, it must be of equal
#'   length as the `.standard_sample`, and each element is taken as the dilution
#'   factor to the previous element, using 1 for the first element. The order of
#'   dilution factors must match that of the ordered `.standard_sample`.
#'
#' @details
#'
#' In the manuals to the LEGENDplex system, standard are labeled 0 to 8, where 8
#' indicate the highest concentration and 0 the background (no analyte). The standard is diluted at 1:4 so that
#'
#' ```
#' [s7] = [start]
#' [s6] = [s7]/4
#' [s5] = [s6]/4
#' [s4] = [s5]/4
#' [s3] = [s4]/4
#' [s2] = [s3]/4
#' [s1] = [s2]/4
#' [s0] = 0
#' ```
#' It might happen, that a dilution step is missing in which case the dilution
#' is corrected to accommodate the missing step. However, since it is inspired
#' guess work and out of the ordinary, a warning is thrown, see Examples.
#'
#' @section Standard sample order:
#' If the vector is numeric, the values are ordered numerically from high to low.
#'
#' If the vector is not numeric, things become a little more difficult, because
#' sorting a vector like `c("a", "c", "0", "b")` by default results in `c("0",
#' "a", "b", "c")`, which means that '0' is the highest value and will be
#' assigned the start concentration and the sample 'a' is then the first
#' dilution.
#'
#' To avoid this problem, the vector is split into two: one containing numerical
#' values and one containing alphabetical. Each vector is then sorted
#' appropriately and combined, see Examples.
#'
#' @return A numeric vector
#' @export
#'
#' @examples
#'
#' calc_std_conc(.standard_sample = c(7:0),
#'                                  .start_concentration = 5000)
#' # Sample 5 is missing
#' calc_std_conc(.standard_sample = c(7, 6, 4, 3, 2, 1, 0),
#'                                  .start_concentration = 5000)
#' calc_std_conc(.standard_sample = rep(c(7:0), 2),
#'                                  .start_concentration = 5000)
#' calc_std_conc(.standard_sample = c(9:0),
#'                                 .start_concentration = 5000)
#'
#' calc_std_conc(.standard_sample = c(letters[1:7], 0),
#'                                  .start_concentration = 5000)
#' calc_std_conc(.standard_sample = c(letters[1:7], 0, 1),
#'                                  .start_concentration = 5000)
#'
#' calc_std_conc(.standard_sample = c(7:1, 0),
#'                                  .start_concentration = 5000,
#'                                  .dilution_factor = c(1, 2, 2, 2, 4, 6, 6, 0))
#'
#' # If 0 exists it is always set to 0
#' calc_std_conc(.standard_sample = c(7:1, 0),
#'                                  .start_concentration = 5000,
#'                                  .dilution_factor = c(1, 2, 2, 2, 4, 6, 6, 100000))
#' calc_std_conc(.standard_sample = c(8:1),
#'                                  .start_concentration = 5000,
#'                                  .dilution_factor = c(1, 2, 2, 2, 4, 6, 6, 100000))
calc_std_conc <- function(.standard_sample, .start_concentration, .dilution_factor = 4L){
  if(length(.dilution_factor) > 1 & length(.dilution_factor) < length(.standard_sample)){
    stop(".dilution_factor must be either equal for all samples, or given for each sample")
  }
  .start_concentration <- max(.start_concentration)

  # Order the standard samples from high to low, but remember the original position
  .ss_org <-  .standard_sample
  .standard_sample <- .standard_sample %>% unique()

  if(length(.standard_sample) < 8){
    warning(paste(8 - length(.standard_sample), "standard samples are missing. The result might not be correct"))
  }

  if("numeric" %in% class(".standard_sample")){
    .standard_sample <- .standard_sample %>% sort(decreasing = TRUE)

    # Get the difference between each sample.
    std_dilution <- .standard_sample %>% diff() %>% abs
  }else{
    # We cannot calculate the difference between the standards because they are
    # not numeric. We must then assume equal difference while preserving the
    # logic of order
    .numeric_values <- grep("^[0-9]*$", .standard_sample)
    .alpha_values <- .standard_sample[!.numeric_values] %>% sort
    .numeric_values <- .standard_sample[.numeric_values] %>% sort(decreasing = TRUE)

    .standard_sample <- c(.alpha_values, .numeric_values)

    std_dilution <- seq_along(1:length(.standard_sample)) %>% diff()
  }

  if(length(.dilution_factor) == 1){
    # It is possible that a sample wasn't collected, and by taking the difference
    # (which is 1 when everything is well) we can calculate the correct dilution if
    # the samples
    std_dilution[std_dilution != 1] <- .dilution_factor^(std_dilution[std_dilution != 1] - 1)
    # The first element is missing, so we insert is. The standards are diluted in a
    # 1:4 serial dilution and by taking the cumulative product we can calculate the
    # actual dilution from the first sample
    std_dilution <-  c(1, (std_dilution * .dilution_factor)) %>% cumprod()
  }else{
    std_dilution <- .dilution_factor %>% cumprod()
  }

  # Calculate the standard concentration for each standard
  std_concentration <- .start_concentration[1]/std_dilution

  if(0 %in% .standard_sample){
    std_concentration[.standard_sample == 0] <- 0
  }

  names(std_concentration) <- .standard_sample
  std_concentration[as.character(.ss_org)] %>% unname()
}

#' Fit a standard curve
#'
#' Fit a logistic function to the standard concentrations.
#'
#' @param .data A tidy data.frame.
#' @param .parameter A character giving the name of column(s) where populations
#'   are identified.
#' @param .concentration A character giving the name of the column with the
#'   standard concentration.
#' @param .fct A character giving the name of the logistic function to use in
#'   the fit, see [drc::drm()] for details.
#' @param ... Other arguments to [drc::drm()]
#'
#' @return An object of class `drc`
#' @export
#'
#' @examples
#'
#' library(beadplexr)
#' library(drc)
#' data(ryegrass)
#'
#' ryegrass_m <-
#'   fit_standard_curve(.data = ryegrass,
#'                      .parameter = "rootl",
#'                      .concentration = "conc")
#'
#' summary(ryegrass_m)
#'
fit_standard_curve <- function(.data,
                               .parameter = "FL2.H",
                               .concentration = "Concentration",
                               .fct = "LL.5", ...){

  .filter_criteria <- lazyeval::interp(~ !is.infinite(conc_column) &
                                         !is.infinite(param_column),
                                       conc_column = as.name(.concentration),
                                       param_column = as.name(.parameter))

  .data <- .data %>%
    dplyr::filter_(.filter_criteria)

  fct <- switch (.fct,
                 LL.2 = drc::LL.2(),
                 l2 = drc::l2(),
                 LL2.2 = drc::LL2.2(),
                 LL.3 = drc::LL.3(),
                 LL.3u = drc::LL.3u(),
                 l3 = drc::l3(),
                 l3u = drc::l3u(),
                 LL2.3 = drc::LL2.3(),
                 LL2.3u= drc::LL2.3u(),
                 LL.4 = drc::LL.4(),
                 l4 = drc::l4(),
                 LL2.4 = drc::LL2.4(),
                 LL.5 = drc::LL.5(),
                 l5 = drc::l5(),
                 LL2.5 = drc::LL2.5()
  )

  # Give function arguments useful values, if they are not set by the user
  .ellipsis <- list(...)

  .fit_formula <- stats::as.formula(paste(.parameter, .concentration, sep = "~"))

  .all_args <- list(formula = .fit_formula, data = .data, fct = fct)

  if("fct" %in% names(.ellipsis)){
    .all_args[["fct"]] <- NULL
  }

  .all_args <- c(.all_args, .ellipsis)

  do.call(drc::drm, .all_args)
}

#' Calculate concentration.
#'
#' Calculate the concentration in a sample
#'
#' @param .data A tidy data.frame.
#' @param .model An object of class `drc` with the fitted dose-response model.
#' @param .parameter A character giving the name of column(s) where populations
#'   are identified.
#' @param .value A character giving the name of the column to store the calculated concentration
#'
#' @return The `.data` with the calculated concentration and error added in two columns.
#' @export
#'
#' @examples
#'
#' library(beadplexr)
#' library(drc)
#' data(ryegrass)
#'
#' ryegrass_m <-
#'   fit_standard_curve(.data = ryegrass,
#'                      .parameter = "rootl",
#'                      .concentration = "conc")
#'
#' sample_data <-
#'   calculate_concentration(.data = ryegrass[sample(1:nrow(ryegrass), 5),],
#'                           .model = ryegrass_m,
#'                           .parameter = "rootl")
calculate_concentration <- function(.data, .model,
                                    .parameter = "FL2.H",
                                    .value = "Calc.conc"){

  .response_values <- .data %>%
    dplyr::select(dplyr::one_of(.parameter)) %>%
    purrr::flatten_dbl()

  ed_res <- suppressWarnings(drc::ED(.model, .response_values,
                                     type = "absolute",
                                     od = TRUE,
                                     display = FALSE)) %>%
    as.data.frame %>%
    stats::setNames(c(.value, paste(.value, "error")))

  .data %>%
    dplyr::bind_cols(ed_res)
}

