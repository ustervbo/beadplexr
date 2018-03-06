#' Find turning points
#'
#' Find turning points (minima and maxima) in a vector.
#'
#' @inheritParams turning_point
#' @param .x A numeric vector
#' @param ... Other parameters passed to [stats::density()]
#'
#' @return A list with the two elements `maxima` and `minima`
#' @keywords internal
#' @export
#'
#' @examples
#' set.seed(1234)
#' .x <- c(rnorm(100, 2, 1), rnorm(100, 9, 1))
#' tpi <- do_find_turning_points(.x, .return = "index", .adjust = 1)
#'
#' dx <- density(.x, adjust = 1, n = length(.x))
#' plot(dx)
#' points(dx$x[tpi$maxima], dx$y[tpi$maxima], pch = 19, col = "red")
#' points(dx$x[tpi$minima], dx$y[tpi$minima], pch = 19, col = "blue")
#'
#' do_find_turning_points(.x, .return = "value", .adjust = 1)
do_find_turning_points <- function(.x,
                                   .return = c("value", "index"),
                                   .adjust = 1.5, ...){
  .return <- match.arg(.return)

  if(length(.adjust) > 1){
    .adjust <- .adjust[1]
    warning(".adjust has length > 1 and only the first element will be used", call. = FALSE)
  }

  # Give function arguments useful values, if they are not set by the user
  .ellipsis <- list(...)
  .all_args <- list(x = .x,
                    n = length(.x),
                    adjust = .adjust)

  if("n" %in% names(.ellipsis)){
    .all_args[["n"]] <- NULL
  }
  if("x" %in% names(.ellipsis)){
    .all_args[["x"]] <- NULL
  }
  if("adjust" %in% names(.ellipsis)){
    .all_args[["adjust"]] <- NULL
  }

  .all_args <- c(.all_args, .ellipsis)

  .data_dens <- do.call(stats::density, .all_args)
  .data_dens_y <- .data_dens$y

  r <- rle(.data_dens_y)
  # These functions ignore the extremes if they are the first or last point
  maxima_index <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2,  times = r$lengths))
  minima_index <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == 2,  times = r$lengths))

  # The functions also return several indicies if the maxima or minima consost
  # of several neighbouring points with identical value
  dup_val_max <- duplicated(.data_dens_y[maxima_index])
  maxima_index <- maxima_index[!dup_val_max]

  dup_val_min <- duplicated(.data_dens_y[minima_index])
  minima_index <- minima_index[!dup_val_min]

  if(.return == "value"){
    maxima_index <- .data_dens$x[maxima_index]
    minima_index <- .data_dens$x[minima_index]
  }

  list(maxima = maxima_index, minima = minima_index)
}

#' Turning points
#'
#' Find turning points (minima and maxima) in a vector.
#'
#' @inheritParams approx_adjust
#' @param .x A numeric vector or a list of numeric vectors. If the list is
#'   named, the names become column names in the returned data.frames
#' @param .which A character indicating the values of interest.
#' @param .return A character giving the desired return type.
#' @param .adjust A numeric giving the adjustment to the `adjust` argument of [stats::density()].
#' @inheritDotParams approx_adjust .lower .upper .step
#'
#' @return A list with the two elements `maxima` and `minima`. each element
#'   consist of a single `data.frame`.
#' @export
#'
#' @examples
#'
#' set.seed(1234)
#' .x <- c(rnorm(100, 2, 1), rnorm(100, 9, 1))
#'
#' turning_point(.x = .x, .adjust = 1)
#' turning_point(.x = .x, .k = 2)
#'
#' turning_point(.x = .x, .which = "minima")
#' turning_point(.x = .x, .which = "maxima")
#'
#' turning_point(.x = .x, .return = "index")
turning_point <- function(.x,
                          .which = c("both", "minima", "maxima"),
                          .return = c("value", "index"),
                          .adjust = 1.5,
                          .k = NULL, ...){

  .which <- match.arg(.which)
  .return <- match.arg(.return)

  if(!"list" %in% class(.x)){
    .x <- list(.x)
  }

  if(!is.null(.k)){
    find_approx_adjust <- function(.p){
      approx_adjust(.x = .x, .k = .k)
    }

    .adjust <- .x %>%
      purrr::map(find_approx_adjust) %>%
      purrr::flatten_dbl()
    .adjust <- .adjust[!is.na(.adjust)]

    if(length(.adjust) == 0){
      stop("I found not a single useful adjust value, and give up!")
    }else{
      message(paste("Using the band width adjustment of: ", paste(.adjust, collapse = ", ")))
    }
  }

  map_turning_point <- function(.x, .a){
    do_find_turning_points(.x = .x, .return = .return, .adjust = .a)
  }

  revert_list <- function(.l) { # @Josh O'Brien
    # get sub-elements in same order
    x <- lapply(.l, `[`, names(.l[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list)
  }

  rtn_lst <- .x %>% purrr::map2(.adjust, map_turning_point) %>% stats::setNames(names(.x)) %>%
  revert_list() %>% purrr::map(function(x){do.call(cbind, x) %>% as.data.frame()})

  if(.which == "maxima"){
    rtn_lst$minima <-NULL
  }else  if(.which == "minima"){
    rtn_lst$maxima <-NULL
  }

  rtn_lst
}

#' Approximate bandwidth adjustment.
#'
#' Approximates the adjust argument to `stats::density()` needed to find the required number of clusters.
#'
#' @param .x A numeric vector.
#' @param .k Numeric giving the number of expected clusters.
#' @param .lower,.upper The interval for possible value of adjust.
#' @param .step A numeric giving the increment to adjust. Sometimes low values
#'   are needed to find a proper adjust value.
#'
#' @details
#' This function finds the first value of the [stats::density()] `adjust`
#' argument which gives the `.k` number of clusters. It it quite crude in that
#' every value of adjust from `.lower` to `.upper` is tested until the desired
#' number of clusters is found. A cluster is defined by a peak, and should no
#' suitable `adjust` value be found, NA is returned.
#'
#'
#' @return A numeric.
#' @export
#'
#' @examples
#'
#' set.seed(1234)
#' .x <- c(rnorm(100, 2, 1), rnorm(100, 9, 1))
#' approx_adjust(.x, 2)
approx_adjust <- function(.x, .k, .lower = 0.4, .upper = 2, .step = 0.2){
  .all_args <- list(.x = .x,
                    .which = c("maxima"),
                    .return = "index",
                    .k = NULL)

  counter <- 0
  max_it <- ceiling((.upper - .lower)/.step)

  cur_adjust <- .lower

  while(cur_adjust <= .upper){
    if(cur_adjust <= 0){
      cur_adjust <- NA
      warning("Found no positive adjustment")
      break
    }

    if(counter >= max_it){
      warning(paste("Did not converge after all possible steps. Last value I tried was",
                    cur_adjust,
                    "- maybe increasing the range and decreasing the step size helps"))
      cur_adjust <- NA
      break
    }
    counter <- counter + 1

    .all_args$.adjust <- cur_adjust

    turn_points <- do.call(turning_point, .all_args)
    num_maxima <- turn_points$maxima %>% nrow

    if(num_maxima < .k){
      cur_adjust <- cur_adjust - .step
    }else if(num_maxima > .k){
      cur_adjust <- cur_adjust + .step
    }else{
      break
    }
  }
  cur_adjust
}
