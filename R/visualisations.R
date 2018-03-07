#' Plot FACS data.
#'
#' Wrappers around building a ggplot with \code{geom_point},
#' \code{geom_density_2d}, and \code{geom_hex}.
#'
#' @param .data The data to be plotted in a \code{data.frame}.
#' @param .x,.y Character vector with the column name for the variable to plot
#'   on the x or y-axis.
#' @param .beads Character vector to with the column name with identification of
#'   beads. If used it will show up with the aesthetic 'colour'. Defaults to not
#'   being used.
#' @param .plot_distinct Boolean to decide if only distinct events should be
#'   plotted. If used, the number of data points might be greatly reduced which
#'   could make for faster plotting. Defaults to TRUE.
#' @param .bins Numeric vector giving number of bins in both vertical and
#'   horizontal directions. Set to 75 by default.
#' @param .type Character vector giving the type of plot being used. Options are
#'   one of \code{"scatter", "density", "hexbin"}.
#' @param ... Arguments passed to the individual functions.
#'
#' @details
#' These plot functions are meant to provide a quick way of viewing the FACS
#' data. For more control, use \code{ggplot2} directly.
#'
#' @return A \code{ggplot}
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' library(beadplexr)
#'
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                           package = "beadplexr")
#'
#' .data <- read_fcs(.file_name = .file_name,
#'                   .filter = list("FSC-A" = c(2e5L, 6.3e5L),
#'                                  "SSC-A" = c(2e5, 1e6L)))
#' .data$bead_group <- ifelse(.data$`FSC-A` < 4e5L, "A", "B")
#'
#' # Using facs_plot
#' facs_plot(.data, .type = "scatter")
#' facs_plot(.data, .type = "density1d")
#' facs_plot(.data, .type = "density2d")
#' facs_plot(.data, .type = "hexbin")
#'
#' facs_plot(.data, .type = "scatter", .beads = "bead_group")
#' facs_plot(.data, .type = "density1d", .beads = "bead_group")
#' facs_plot(.data, .type = "hexbin", .bins = 50)
#'
#' facs_plot(.data, .x = "FL2-H", .type = "scatter", .beads = "bead_group")
#'
#' # Individual functions
#' facs_scatter(.data)
#'
#' facs_scatter(.data, .beads = "bead_group", .plot_distinct = FALSE)
#' facs_scatter(.data, .beads = "bead_group")
#'
#' facs_scatter(.data, .x = "FL2-H", .y = "FL6-H", .beads = "bead_group")
#'
#' facs_density1d(.data)
#' facs_density1d(.data, .beads = "bead_group")
#'
#' facs_density2d(.data)
#' facs_density2d(.data, .beads = "bead_group")
#'
#' facs_hexbin(.data)
#' facs_hexbin(.data, .bins = 30)
#' }
facs_plot <- function(.data, .x = "FSC-A", .y = "SSC-A", .type = c("scatter", "density1d", "density2d", "hexbin"), ...){
  .type <- match.arg(.type)
  switch(.type,
         scatter = facs_scatter(.data = .data, .x = .x, .y = .y, ...),
         density1d = facs_density1d(.data = .data, .x = .x, ...),
         density2d = facs_density2d(.data = .data, .x = .x, .y = .y, ...),
         hexbin = facs_hexbin(.data = .data, .x = .x, .y = .y, ...))
}

#' @rdname facs_plot
#' @export
facs_scatter <- function(.data, .x = "FSC-A", .y = "SSC-A", .beads = NULL, .plot_distinct = TRUE){

  .data <- .data %>% dplyr::ungroup()

  if(.plot_distinct){
    if(!is.null(.beads)){
      .data <- .data %>%
        dplyr::distinct_(as.name(.x), as.name(.y), as.name(.beads))
    }else{
      .data <- .data %>%
        dplyr::distinct_(as.name(.x), as.name(.y))
    }

  }

  tmp_plot <- .data %>%
    ggplot2::ggplot() +
    ggplot2::aes_(x = as.name(.x), y = as.name(.y)) +
    ggplot2::geom_point(size = 0.5)
    # geom_density2d() +

  if(!is.null(.beads)){
    tmp_plot <- tmp_plot + ggplot2::aes_(colour = as.name(.beads))
  }
  tmp_plot
}

#' @rdname facs_plot
#' @export
facs_density2d <- function(.data, .x = "FSC-A", .y = "SSC-A", .beads = NULL){

  .data <- .data %>% dplyr::ungroup()

  tmp_plot <- .data %>%
    ggplot2::ggplot() +
    ggplot2::aes_(x = as.name(.x), y = as.name(.y)) +
    ggplot2::geom_density_2d()

    if(!is.null(.beads)){
      tmp_plot <- tmp_plot + ggplot2::aes_(colour = as.name(.beads))
    }
  tmp_plot
}

#' @rdname facs_plot
#' @export
facs_density1d <- function(.data, .x = "FSC-A", .beads = NULL){

  .data <- .data %>% dplyr::ungroup()

  tmp_plot <- .data %>%
    ggplot2::ggplot() +
    ggplot2::aes_(x = as.name(.x)) +
    ggplot2::geom_density()

  if(!is.null(.beads)){
    tmp_plot <- tmp_plot + ggplot2::aes_(fill = as.name(.beads))
  }
  tmp_plot
}

#' @rdname facs_plot
#' @export
facs_hexbin <- function(.data, .x = "FSC-A", .y = "SSC-A", .bins = 75){
  .data <- .data %>% dplyr::ungroup()
  .data %>%
    ggplot2::ggplot() +
    ggplot2::aes_(x = as.name(.x), y = as.name(.y)) +
    ggplot2::geom_hex(bins = .bins) +
    ggplot2::scale_fill_gradientn(colours = c("blue", "green", "yellow", "red")) +
    ggplot2:: theme(legend.position = "none")
}

#' Plot concentrations
#'
#' @param .data A `data.frame` with the data to be plotted.
#' @param .sample_data A `data.frame` with the calculated sample concentrations.
#' @param .standard_data A `data.frame` with the calculated standard
#'   concentrations.
#' @param .model An object of class `drc` with the fitted dose-response model.
#' @param .title A character giving the title of the plot.
#' @param .parameter A character giving the name of the column with the MFI
#' @param .concentration  A character giving the name of the column with the with the calculated concentrations.
#' @param .std_concentration A character giving the name of the column with the standard concentration.
#'
#' @return A \code{ggplot}
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @name plot_concentrations
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
#' recalc_std <-
#'   calculate_concentration(.data = ryegrass,
#'                           .model = ryegrass_m,
#'                           .parameter = "rootl")
#' sample_data <-
#'   calculate_concentration(.data = ryegrass[sample(1:nrow(ryegrass), 5),],
#'                           .model = ryegrass_m,
#'                           .parameter = "rootl")
#'
#' plot_std_curve(ryegrass,
#'                ryegrass_m,
#'                .parameter = "rootl",
#'                .concentration = "conc")
#'
#' plot_target_est_conc(.data = recalc_std,
#'                      .concentration = "Calc.conc",
#'                      .std_concentration = "conc")
#' plot_estimate(
#'   .sample_data = sample_data,
#'   .standard_data = ryegrass,
#'   .model = ryegrass_m,
#'   .parameter = "rootl",
#'   .concentration = "conc")
#'
plot_std_curve <- function(.data, .model, .title = NULL,
                           .parameter = "FL2.H",
                           .concentration = "Concentration"){

  # Create date for displaying fit line and interval
  .filter_criteria <- lazyeval::interp(~ !is.infinite(conc_column) &
                                         !is.infinite(param_column),
                                       conc_column = as.name(.concentration),
                                       param_column = as.name(.parameter))
  fit_line_range <- .data %>%
    dplyr::filter_(.filter_criteria) %>%
    dplyr::select(dplyr::one_of(.concentration)) %>%
    purrr::flatten_dbl() %>%
    range()

  fit_line_range <- data.frame(x = seq(min(fit_line_range), max(fit_line_range), length = 100))
  fit_line_range <- dplyr::bind_cols(fit_line_range,
                                     suppressWarnings(stats::predict(.model,
                                                    newdata = fit_line_range,
                                                    interval = "confidence")) %>%
                                       as.data.frame) %>%
    stats::setNames(c(.concentration, .parameter, "Lower", "Upper"))

  .data %>%
    ggplot2::ggplot() +
    ggplot2::aes_string(x = .concentration, y = .parameter) +
    ggplot2::geom_point()  +
    ggplot2::geom_ribbon(data = fit_line_range,
                         ggplot2::aes_string(ymin = "Lower", ymax = "Upper"), alpha=0.2) +
    ggplot2::geom_line(data = fit_line_range, color = "blue") +
    ggplot2::labs(title = .title)
}

#' @rdname plot_concentrations
#' @export
plot_target_est_conc <- function(.data, .title = NULL,
                                 .concentration = "Calc.conc",
                                 .std_concentration = "Concentration"){
  filter_fun <- function(.p){
    lazyeval::interp(~ !is.infinite(.col) &
                       !is.na(.col) &
                       !is.nan(.col), .col = as.name(.p))
  }

  .data <- .data %>%
    dplyr::filter_(filter_fun(.std_concentration)) %>%
    dplyr::filter_(filter_fun(.concentration))

  .fit_formula <- stats::as.formula(paste(.concentration, .std_concentration, sep = "~"))

  lm_res <- stats::lm(formula = .fit_formula, data = .data)
  lm_res_summary <- summary(lm_res)
  r_squared <- lm_res_summary$r.squared
  p_value <- 1 - stats::pf(lm_res_summary$fstatistic[1],
                    lm_res_summary$fstatistic[2],
                    lm_res_summary$fstatistic[3])
  slope <- stats::coef(lm_res_summary)[2]

  .filter_criteria <- lazyeval::interp(~ !is.infinite(est_col) &
                                         !is.infinite(target_col),
                                       est_col = as.name(.concentration),
                                       target_col = as.name(.std_concentration))
  .get_min_x <- lazyeval::interp(~ min(x, na.rm = TRUE), x = as.name(.std_concentration))
  .get_max_y <- lazyeval::interp(~ max(y, na.rm = TRUE), y = as.name(.concentration))

  fit_text <- .data %>%
    dplyr::filter_(.filter_criteria) %>%
    dplyr::summarise_(.dots = stats::setNames(list(.get_min_x, .get_max_y),
                                c(.std_concentration, .concentration))) %>%
    dplyr::mutate(label = paste("R2:", format(r_squared, digits = 3), "\n",
                         "Slope:", format.pval(slope, digits = 3), "\n",
                         "p:", format.pval(p_value, digits = 3)))

  .data %>%
    ggplot2::ggplot() +
    ggplot2::aes_string(x = .std_concentration, y = .concentration) +
    ggplot2::geom_point() +
    ggplot2::geom_text(data = fit_text,
                       ggplot2::aes_string(label = "label"),
                       vjust = 1, hjust = 0, size = 3) +
    ggplot2::stat_smooth(method = "lm") +
    ggplot2::labs(title = .title)
}

#' @rdname plot_concentrations
#' @export
plot_estimate <- function(.sample_data, .standard_data, .model, .title = NULL,
                          .parameter = "FL2.H",
                          .concentration = "Concentration"){


  # Create date for displaying fit line and interval
  .filter_criteria <- lazyeval::interp(~ !is.infinite(conc_column) &
                                         !is.infinite(param_column),
                                       conc_column = as.name(.concentration),
                                       param_column = as.name(.parameter))

  fit_line_range <- .standard_data %>%
    dplyr::filter_(.filter_criteria) %>%
    dplyr::select(dplyr::one_of(.concentration)) %>%
    purrr::flatten_dbl() %>%
    range()

  fit_line_range <- data.frame(x = seq(min(fit_line_range),
                                       max(fit_line_range),
                                       length = 100))
  fit_line_range <-  dplyr::bind_cols(fit_line_range,
                                      suppressWarnings(stats::predict(.model,
                                                     newdata = fit_line_range,
                                                     interval = "confidence")) %>%
                                        as.data.frame) %>%
    stats::setNames(c(.concentration, .parameter, "Lower", "Upper"))

  .standard_data %>%
    ggplot2::ggplot() +
    ggplot2::aes_string(x = .concentration, y = .parameter) +
    ggplot2::geom_ribbon(data = fit_line_range,
                         ggplot2::aes_string(ymin = "Lower", ymax = "Upper"),
                         alpha=0.2) +
    ggplot2::geom_line(data = fit_line_range, color = "blue") +
    ggplot2::geom_hline(data = .sample_data,
                        ggplot2::aes_string(yintercept = .parameter),
                        linetype = 2) +
    ggplot2::labs(title = .title)
}

