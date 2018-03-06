#' Despeckle parameters
#'
#' Remove lonely, noisy data points in a 2D scatter matrix
#'
#' @param .data A tidy data.frame
#' @param .parameters A character of the length of two giving the parameters to despeckle.
#' @param .bins A numeric giving the resolution of the raster matrix.
#' @param .neighbours A numeric giving the minimum number of neighbours. Points with fewer neighbours are removed.
#' @param ... Additional parameters passed to [raster::clump()]
#'
#' @details
#' The values of the two parameters are binned into the given number of bins.
#' They are then cast into a 2D matrix, with the bins of the first of the
#' parameters ending up as rows, the bins of the second parameter as
#' columns, and combinations are marked by `1`.
#'
#' This matrix is turned into a `RasterLayer` by [raster::raster()] and the
#' number of neighbours are calculated by [raster::clump()].
#'
#' The rows of the `.data` where lonely points are found in `.parameters` are removed.
#'
#' @note
#' This function requires that the `igraph` package is available.
#'
#' @return A `data.frame` with noisy points removed.
#'
#' @importFrom magrittr "%>%"
#' @importFrom raster freq "%in%"
#'
#' @export
#'
#' @examples
#' library(beadplexr)
#' library(magrittr)
#' library(ggplot2)
#'
#' data("lplex")
#'
#' lplex[[1]] %>%
#'   ggplot() +
#'   aes(x = `FL6-H`, y = `FL2-H`) +
#'   geom_point()
#'
#' lplex[[1]] %>%
#'   despeckle(.parameters = c("FL6-H", "FL2-H"), .neighbours = 8) %>%
#'   ggplot() +
#'   aes(x = `FL6-H`, y = `FL2-H`) +
#'   geom_point()
#'
#' lplex[[1]] %>%
#'   despeckle(.parameters = c("FL6-H", "FL2-H"), .bin = 128, direction = 4) %>%
#'   ggplot() +
#'   aes(x = `FL6-H`, y = `FL2-H`) +
#'   geom_point()
#'
despeckle <- function(.data, .parameters, .bins = 256L, .neighbours = 4L, ...){
  if(length(.parameters) != 2){
    stop(".parameters must be of length two. No more, and no less")
  }
  if(FALSE %in% (.parameters %in% names(.data))){
    stop(".parameters not found in data")
  }

  .col_names_bin <- c("x_bin", "y_bin")
  names(.parameters) <- .col_names_bin
  # We use the raster package to calculate the number of neighbours each point
  # have when the two paramters are plotted against eachother. This works best
  # if we have integers, but we must also be able to  identify the points we
  # have discharge in the original data. By breaking the values of the two
  # parameters into numbered bins, and by keeping this information, we
  # have what we need
  .set_bins <- function(.parameters, .bins){
    lazyeval::interp(~base::cut(.parameters, breaks = .bins, labels = seq_along(1:.bins)),
                     .bins = .bins,
                     .parameters = as.name(.parameters))
  }

  .binned_data <- .data %>%
    dplyr::mutate_(.dots = lapply(.parameters, .set_bins, .bins = .bins)) %>%
    dplyr::mutate_at(.col_names_bin, as.character)

  # Raster only works with a matrix so this we need to have.
  # The columns of interest in have 'bins' in their names, and since the same
  # x_bin - y_bin combination can be found more than once we remove all
  # duplicates before spreading, setting rownames and cast to matrix
  .present_value <-  stats::setNames(list(lazyeval::interp(~1L)), "present")

  .binned_matrix <- .binned_data %>%
    dplyr::select(dplyr::one_of(.col_names_bin)) %>%
    dplyr::distinct() %>%
    dplyr::mutate_(.dots = .present_value) %>%
    tidyr::spread_(key_col = .col_names_bin[2], value_col = "present", fill = 0L) %>%
    # Rownames are not allowed on tibbles so we cast to a regular data.frame
    as.data.frame() %>%
    # For some reason, the rownames is set already, so it must be removed
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = .col_names_bin[1]) %>% as.matrix

  # It a little bit of work getting to the values with low number of neighbours.
  # First we have to find the number of neighboutrs ofcourse and the ones below
  # the cutoff so we can the number of neightbours to nothing
  .raster_res <- raster::raster(.binned_matrix)
  # This function makes test_that give warnings. I don't know it is stemms from
  # igraph or where, but annoying it is
  .clump_res <- suppressWarnings(raster::clump(.raster_res, ...))

  # clump assigns each data point a number and a count of neighbours
  .neighbour_filter <-lazyeval::interp(
    ~ .the_column < .the_neighbours,
    .the_column = as.name("count"),
    .the_neighbours = .neighbours)

  .speckles <- .clump_res %>% freq %>%
    as.data.frame() %>%
    dplyr::filter_(.dots = .neighbour_filter)

  .clump_res[.clump_res %in% .speckles$value] <- NA
  # The number of neighbours are stored in a single vector which we must turn
  # into a tidy data.frame. The column and row names of the binned matrix
  # corrospond to the bin number, which is what we need to filter the original
  # data
  .neighbour_filter <- lazyeval::interp(~!is.na(.the_column), .the_column = as.name("raster_number"))

  .set_x_bin <- lazyeval::interp(".the_rownames", .the_rownames = rownames(.binned_matrix))
  .set_x_bin <- stats::setNames(list(.set_x_bin), .col_names_bin[1])

  .despeckled_data <- .clump_res@data@values %>%
    as.integer() %>%
    matrix(ncol = ncol(.binned_matrix), byrow = TRUE) %>%
    dplyr::as_data_frame() %>%
    stats::setNames(colnames(.binned_matrix)) %>%
    dplyr::mutate_(.dots = .set_x_bin) %>%
    # The row names are converted to double, though they are stored as character
    dplyr::mutate_at(names(.set_x_bin), as.character) %>%
    tidyr::gather_(key_col = .col_names_bin[2], value_col = "raster_number", colnames(.binned_matrix)) %>%
    dplyr::filter_(.dots = .neighbour_filter)

  # Having the bin numbers of all data points with more than the minimum number
  # of neighbours we can now filtermclean and return
  .despeckled_data %>%
    dplyr::left_join(.binned_data, by = .col_names_bin) %>%
    dplyr::select_("-x_bin", "-y_bin", "-raster_number")
}

#' Chebyshev distance
#'
#' @param x a numeric matrix or data frame.
#' @param diag logical value indicating whether the diagonal of the distance
#'   matrix should be printed by `print.dist`.
#' @param upper logical value indicating whether the upper triangle of the
#'   distance matrix should be printed by `print.dist` .
#'
#' @return Chebyshev distance returns an object of class "`dist`".
#' @export
#'
#' @examples
#'
#' x <- matrix(rnorm(100), nrow = 5)
#' dist_chebyshev(x)
#'
dist_chebyshev <- function(x, diag = FALSE, upper =FALSE) {
  x <- x %>% as.matrix()
  N <- nrow(x)
  .pairs <- utils::combn(N, 2)

  .dist <- apply(.pairs, MARGIN = 2, function(i){
    max(
      abs(
        x[i[1], ] - x[i[2], ]
      )
    )
  })

  attr(.dist, "Size") <- N
  attr(.dist, "Labels") <- dimnames(x)[[1]]
  attr(.dist, "Diag") <- diag
  attr(.dist, "Upper") <- upper
  attr(.dist, "method") <- "Chebyshev"
  attr(.dist, "call") <- match.call()
  class(.dist) <- "dist"
  .dist
}

#' Cast list of analytes to `data.frame`
#'
#' A well structured list, such at those loaded by [load_panel()], is cast to a data.frame.
#'
#' @param .analytes The named list to be cast. It usually is loaded using
#'   [load_panel()]. See Details for expected structure.
#' @param .id_bead,.id_analyte The name of the column to hold the bead group and the analyte ID, respectively.
#'
#' @details
#'
#' Each analyte in the `list` passed to the function is expected to be a named
#' list with named elements `name` and `concentration`. The name of the `list`
#' with the analyte specific information is the analyte ID.
#'
#' Because of the particular setup of the LEGENDplex assay with two bead groups,
#' the analytes are expected to be wrapped in another `list`.
#'
#' @return A `data-frame`
#' @export
#'
#' @examples
#' .analytes <-  list(A = list(
#'                 A1 = list(name = "name_a1", concentration = 500),
#'                 A2 = list(name = "name_a2", concentration = 50000)))
#'
#' as_data_frame_analyte(.analytes)
as_data_frame_analyte <- function(.analytes, .id_bead = "Bead group", .id_analyte = "Analyte ID"){
  if(! "list" %in% class(.analytes)){
    stop("I expect .analytes to be a list")
  }

  .analytes %>%
    purrr::map_df(function(.l){
      .l %>% purrr::map_df(function(.x){
        tibble::data_frame(name = .x$name, concentration = .x$concentration)
      }, .id = .id_analyte)
    }, .id = .id_bead)

}
