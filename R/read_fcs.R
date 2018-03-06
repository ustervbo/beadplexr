#' Read a fcs file.
#'
#' This is a wrapper around reading a FACS data file using \code{flowCore},
#' transforming the channels with bead signal events using
#' \code{arcsinhTransform()} from \code{flowCore}, automatic removing boundary
#' events of the forward and side scatter channels, subsetting channels (if
#' needed), and cast to a \code{data.frame}.
#'
#' @param .file_name The path and name of the file to be read.
#' @param .fsc_ssc The names of the forward and side scatter channels. A
#'   character vector of length of two.
#' @param .bead_channels The names of the channels with bead events. A character
#'   vector of length of at least two.
#' @param .filter Optional list of upper and lower cutoff for individual
#'   channels. Use \code{.filter = NULL} for no filtering at all.
#' @param .compensation A character vector, a compensation matrix, or
#'   \code{NULL}. See 'Details' for extended information of the argument.
#'
#' @param ... additional arguments passed to \code{flowCore::read.FCS}
#'
#' @details
#'
#' The \code{.compensation} argument takes a character vector or an actual
#' compensation matrix. In case of the latter, it must be a object of class
#' \code{matrix}. If an object of class \code{character} is given to the
#' \code{.compensation} argument, it can be the keyword 'guess' or a search
#' pattern matching the keyword in the fcs file that contains the compensation
#' matrix. If the keyword 'guess' is passed, the function looks for a matrix at
#' the keywords "SPILL" and "SPILLOVER", and if none are found it takes the
#' first matrix which contains the parameter names as column names. Finally, the
#' values of \code{.compensation} can also be \code{NULL}, in which case, no
#' compensation is performed.
#'
#' To summarise, the argument \code{.compensation} can be:
#'
#' \describe{
#'   \item{A matrix}{The compensation matrix to apply}
#'   \item{A character}{
#'     \describe{
#'       \item{The word 'guess'}{Looks for a matrix at the keywords "SPILL"
#'       and "SPILLOVER". If none found the first matrix with parameter names is
#'       applied}
#'       \item{A search string}{The first matrix found within the keywords
#'       matching the given search string is applied. The search string can be a
#'       regular expression}
#'     }
#'   }
#'   \item{NULL}{Nothing is done}
#' }
#'
#' @return A \code{data.frame} with the columns given in \code{.fsc_ssc} and \code{.bead_channels}.
#' @export
#'
#' @seealso \code{\link[flowCore]{rectangleGate}} for instructions to the
#'   \code{.filter} list argument, \code{\link[flowCore]{boundaryFilter}}
#'   for automatic removal of boundary events, and
#'   \code{\link[flowCore]{arcsinhTransform}} for the transformation.
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' library(beadplexr)
#'
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                           package = "beadplexr")
#'
#' # Load the fcs file, with no filter
#' .data <- read_fcs(.file_name = .file_name, .filter = NULL)
#'
#' plot(.data[, c("FSC-A", "SSC-A")])
#' plot(.data[, c("FL2-H", "FL6-H")])
#'
#' # Load the fcs file, with default filter
#' .data <- read_fcs(.file_name = .file_name, .fsc_ssc = c("FSC-H", "SSC-H"))
#'
#' plot(.data[, c("FSC-H", "SSC-H")])
#' plot(.data[, c("FL2-H", "FL6-H")])
#'
#' # Load the fcs file, with custom filter
#' .data <- read_fcs(.file_name = .file_name,
#'                   .filter = list("FSC-A" = c(2e5L, 3.3e5L),
#'                                  "SSC-A" = c(2e5, 1e6L),
#'                                  "FL2-H" = c(11L, 14L),
#'                                  "FL6-H" = c(9L, Inf)))
#'
#' plot(.data[, c("FSC-A", "SSC-A")])
#' plot(.data[, c("FL2-H", "FL6-H")])
#'
#' # Specify three bead channels
#'
#' .data <- read_fcs(.file_name = .file_name,
#'                   .bead_channels = c("FL6-H", "FL2-H", "FL1-H"))
#'
#' plot(.data[, c("FSC-A", "SSC-A")])
#' plot(.data[, c("FL2-H", "FL6-H")])
#' plot(.data[, c("FL1-H", "FL2-H")])
#' plot(.data[, c("FL1-H", "FL6-H")])
#'
read_fcs <- function(.file_name, .fsc_ssc = c("FSC-A", "SSC-A"), .bead_channels = c("FL6-H", "FL2-H"),
                     .filter = list("FSC-A" = c(2e5L, 8e5L), "SSC-A"= c(2e5L, 1e6L), "FL6-H" = c(7.3, Inf)),
                     .compensation = "guess", ...){

  if(!is.null(.fsc_ssc)){
    if(length(.fsc_ssc) != 2 | !is.character(.fsc_ssc)){
      stop(".fsc_ssc must be a character vector of length 2")
    }
  }

  if(!is.null(.bead_channels)){
    if(length(.bead_channels) < 2 | !is.character(.bead_channels)){
      stop(".bead_channels must be a character vector of at least 2")
    }
  }


  flowCore::read.FCS(filename = .file_name, transformation = FALSE, ...) %>%
    transform_bead_channels(.bead_channels = .bead_channels) %>%
    apply_compensation(.compensation = .compensation) %>%
    remove_boundary_events(.channels = c(.fsc_ssc)) %>%
    subset_channels(.filter = .filter) %>%
    as_data_frame_flow_frame(.channels = c(.fsc_ssc, .bead_channels))
}

#' Remove boundary events.
#'
#' A wrapper around application of \code{flowCore}'s \code{boundaryFilter}.
#'
#' @param .flow_frame A \code{flowFrame}. Usually the result of \code{read.FCS} from \code{flowCore}.
#' @param .channels A character vector with the channels to apply the filter to (mostly just forward and side scatter).
#'
#' @return A \code{FlowFrame} with border events removed.
#'
#' @seealso \code{\link[flowCore]{boundaryFilter}} for comments on the filter.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' library(beadplexr)
#' library(flowCore)
#'
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                           package = "beadplexr")
#'  # Load the fcs file
#'  .flow_frame <- read.FCS(filename = .file_name,
#'                          transformation = FALSE)
#'
#' # Plot with all events
#' w_bdr <- .flow_frame@exprs
#' plot(w_bdr[, c("FSC-A", "SSC-A")])
#'
#' # Remove boundary events
#' .flow_frame <- remove_boundary_events(.flow_frame,
#'                           .channels = c("FSC-A", "SSC-A"))
#'
#' wo_bdr <- as.data.frame(.flow_frame@exprs)
#' plot(wo_bdr[, c("FSC-A", "SSC-A")])
#' }
remove_boundary_events <- function(.flow_frame, .channels){
   if(is.null(.channels)){
     return(.flow_frame)
   }

  bf <- flowCore::boundaryFilter(x = .channels, side = "both", filterId = "Automatic Boundary")
  bf <- flowCore::filter(x = .flow_frame, filter = bf)
  flowCore::Subset(.flow_frame, bf)
}

#' Subset channels
#'
#' A wrapper around the application of \code{flowCore}'s \code{rectangleGate}.
#'
#' @param .flow_frame A \code{flowFrame}. Usually the result of \code{read.FCS}
#'   from \code{flowCore}.
#' @param .filter An Optional list of upper and lower cutoff for individual
#'   channels. Default is no filtering at all.
#'
#' @return A \code{FlowFrame} with border events removed.
#'
#' @seealso \code{\link[flowCore]{rectangleGate}} for instructions to the
#'   \code{.filter} list argument
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' library(beadplexr)
#' library(flowCore)
#'
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                             package = "beadplexr")
#' # Load the fcs file
#'  .flow_frame <- read.FCS(filename = .file_name,
#'                           transformation = FALSE)
#'
#'  # Plot with all events
#'  all_events <- .flow_frame@exprs
#'  plot(all_events[, c("FSC-A", "SSC-A")])
#'  # Events are untransformed
#'  plot(all_events[, c("FL2-H", "FL6-H")])
#'
#'  # A silly thing that does nothing
#'  .flow_frame <- subset_channels(.flow_frame)
#'  filter_1 <- .flow_frame@exprs
#'  plot(filter_1[, c("FSC-A", "SSC-A")])
#'
#'  # Filter on forward scatter
#'  .flow_frame <- subset_channels(.flow_frame,
#'                                 .filter = list("FSC-A" = c(1e3L, 2e6)))
#'  filter_2 <- .flow_frame@exprs
#'  plot(filter_2[, c("FSC-A", "SSC-A")])
#'
#'  # Filter on forward and side scatter
#'  .flow_frame <- subset_channels(.flow_frame,
#'                                 .filter = list("FSC-A" = c(2e5L, 6.5e5L),
#'                                                 "SSC-A" = c(2e5, 1e6L)))
#'  filter_3 <- .flow_frame@exprs
#'  plot(filter_3[, c("FSC-A", "SSC-A")])
#'
#'  # Filter on the unfiltered bead channels (please transform first)
#'  .flow_frame <- subset_channels(.flow_frame,
#'                                 .filter = list("FL2-H" = c(5e4L, 6e5L),
#'                                                "FL6-H" = c(0, 3e5L)))
#'  filter_4 <- .flow_frame@exprs
#'  plot(filter_4[, c("FL2-H", "FL6-H")])
#'
#'  # And everything at once
#'  .flow_frame <- read.FCS(filename = .file_name,
#'                         transformation = FALSE)
#' .flow_frame <- subset_channels(.flow_frame,
#'                                .filter = list("FSC-A" = c(2e5L, 6.5e5L),
#'                                               "SSC-A" = c(2e5, 1e6L),
#'                                               "FL2-H" = c(5e4L, 6e5L),
#'                                               "FL6-H" = c(0, 3e5L)))
#' filter_5 <- .flow_frame@exprs
#' plot(filter_5[, c("FSC-A", "SSC-A")])
#' plot(filter_5[, c("FL2-H", "FL6-H")])
#' }
subset_channels <- function(.flow_frame, .filter = NULL){
  if(is.null(.filter)){
    return(.flow_frame)
  }

  gf <- flowCore::rectangleGate(.filter, filterId = "Manual Boundary")
  gf <- flowCore::filter(x = .flow_frame, filter = gf)
  flowCore::Subset(.flow_frame, gf)
}

#' Cast \code{flowFrame} to \code{data_frame}.
#'
#' Extracts the flow events from the \code{flowFrame} as a \code{data_frame}
#'
#' @param .flow_frame A \code{flowFrame}. Usually the result of \code{read.FCS} from \code{flowCore}.
#' @param .channels A character vector with the event channels to extract. Default is all channels.
#'
#' @return A \code{data_frame}
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' library(beadplexr)
#' library(flowCore)
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                           package = "beadplexr")
#' # Load the fcs file
#' .flow_frame <- read.FCS(filename = .file_name,
#'                         transformation = FALSE)
#'
#' # Get all channels
#' as_data_frame_flow_frame(.flow_frame)
#'
#' # Just interesting
#' as_data_frame_flow_frame(.flow_frame,
#'                          .channels = c("FSC-A", "SSC-A",
#'                                        "FL2-H", "FL6-H"))
#' }
as_data_frame_flow_frame <- function(.flow_frame, .channels = NULL){
  if(is.null(.channels)){
    .channels <- flowCore::featureNames(.flow_frame)
  }

  .flow_frame@exprs[, .channels] %>%
    as.data.frame.matrix() %>% tibble::as_data_frame()
}

#' Transform parameters with bead events.
#'
#' Channels with bead events are transformed using \code{flowCore}'s \code{arcsinhTransform()}.
#'
#' @param .flow_frame A \code{flowFrame}. Usually the result of \code{read.FCS} from \code{flowCore}.
#' @param .bead_channels A character vector of the length of 2 with the names of the channels to transform.
#'
#' @return A \code{flowFrame} with the beads channels transformed
#'
#' @keywords internal
#'
#' @seealso \code{\link[flowCore]{arcsinhTransform}} for the transformation.
#'
#' @examples
#' \dontrun{
#' library(beadplexr)
#' library(flowCore)
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                           package = "beadplexr")
#' # Load the fcs file
#' .flow_frame <- read.FCS(filename = .file_name,
#'                                  transformation = FALSE)
#' # Transform channels
#' .flow_frame <- transform_bead_channels(.flow_frame = .flow_frame,
#'                             .bead_channels =
#'                                 c("FL6-H", "FL2-H"))
#' }
transform_bead_channels <- function(.flow_frame, .bead_channels){
  if(is.null(.bead_channels)){
    return(.flow_frame)
  }
  if(length(.bead_channels) < 2){
    stop(".bead_channels must at least be of length 2")
  }
  # If any of the values in bead_channels are not found we cannot transform or
  # do anything
  in_fn <- .bead_channels %in% flowCore::featureNames(.flow_frame)

  if(FALSE %in% in_fn){
    stop(paste(paste(.bead_channels[!in_fn], collapse = " and "), "not found in the flow frame", sep = " "))
  }

  # Transform channels and return
  transformtion_list <- flowCore::transformList(.bead_channels, flowCore::arcsinhTransform(), transformationId = "defaultArcsinhTransform")
  flowCore::transform(.flow_frame, transformtion_list)
}

#' Apply compensation
#'
#' A compensation matrix is applied
#'
#' @param .flow_frame A \code{flowFrame}. Usually the result of \code{read.FCS}
#'   from \code{flowCore}.
#' @param .compensation A character vector, a compensation matrix, or
#'   \code{NULL}. See 'Details' for extended information of the argument.
#'
#' @details
#'
#' The \code{.compensation} argument takes a character vector or an actual
#' compensation matrix. In case of the latter, it must be a object of class
#' \code{matrix}. If an object of class \code{character} is given to the
#' \code{.compensation} argument, it can be the keyword 'guess' or a search
#' pattern matching the keyword in the fcs file that contains the compensation
#' matrix. If the keyword 'guess' is passed, the function looks for a matrix at
#' the keywords "SPILL" and "SPILLOVER", and if none are found it takes the
#' first matrix which contains the parameter names as column names. Finally, the
#' values of \code{.compensation} can also be \code{NULL}, in which case, no
#' compensation is performed.
#'
#' To summarise, the argument \code{.compensation} can be:
#'
#' \describe{
#'   \item{A numerical matrix}{The compensation matrix to apply}
#'   \item{A character}{
#'     \describe{
#'       \item{The word 'guess'}{Looks for a matrix at the keywords "SPILL"
#'       and "SPILLOVER". If none found the first matrix with parameter names is
#'       applied}
#'       \item{A search string}{The first matrix found within the keywords
#'       matching the given search string is applied. The search string can be a
#'       regular expression}
#'     }
#'   }
#'   \item{NULL}{Nothing is done}
#' }
#'
#' @return A compensated \code{flowFrame}
#'
#' @keywords internal
#'
#' @importFrom magrittr "%>%"
#'
#' @seealso \code{\link[flowCore]{compensation-class}} for the compensation.
#'
#' @examples
#' \dontrun{
#' library(beadplexr)
#' library(flowCore)
#' .file_name <- system.file("extdata", "K2-C07-A7.fcs",
#'                           package = "beadplexr")
#' # Load the fcs file
#' .flow_frame <- read.FCS(filename = .file_name,
#'                                  transformation = FALSE)
#' # Apply compensation by guessing the matrix
#' .flow_frame <- apply_compensation(.flow_frame = .flow_frame,
#'                             .compensation = "guess")
#' }
#'
apply_compensation <- function(.flow_frame, .compensation){
  # Got nothing, do nothing
  if(is.null(.compensation)){
    return(.flow_frame)
  }

  if(purrr::is_character(.compensation) & length(.compensation) > 1){
    warning(".compensation can only have the length of 1. I am using just the first element")
    .compensation <- .compensation[1]
  }

  # Test if .x is a numerical matrix with the parameters
  .get_potential_comp_matrices <- function(.x, .p){
    if((class(.x) == "matrix") & purrr::is_numeric(.x)){
      TRUE %in% (.p %in% colnames(.x))
    }else{
      FALSE
    }
  }

  # Search keywords and unify warning
  .search_for_keywords <- function(.keywords, .pattern){
    .matches <- grepl(.pattern, .keywords, ignore.case = TRUE) %>% which
    if(length(.matches) > 1){
      warning("Found more than one potential compensation matrix. Applying the first.")
      .matches <- .matches[1]
    }
    .matches
  }

  # We get either a character or something which can be turned into a numerical matrix
  if(purrr::is_character(.compensation)){
    ff_keywords <- flowCore::keyword(.flow_frame)
    ff_parameters <- flowCore::colnames(.flow_frame)

    # No matter if we guess or were given a search word, we need a numerical
    # matrix where the parameters are in the column names. Otherwise, the actual
    # compensation application will fail.
    num_matrix <- ff_keywords %>%
      purrr::map_lgl(.get_potential_comp_matrices, .p = ff_parameters)
    num_matrix <- ff_keywords[num_matrix]

    # If we found nothing we cannot do nothing - except give a warning
    if(length(num_matrix) == 0){
      warning("No compensation applied. I could not find an appropriate matrix in the flow frame.")
      return(.flow_frame)
    }

    if(.compensation == "guess"){
      # Do the numerical matricies have names that we could expect?
      match_names <- names(num_matrix) %>% .search_for_keywords(.pattern = "SPILL")

      if(length(match_names) > 0){
        comp_matrix <- num_matrix[[match_names]]
      }else{
        # We found no standard keywords. If we have more than one potential
        # matrix, we warn
        if(length(num_matrix) > 1){
          warning("Found more than one potential compensation matrix. Applying the first.")
        }
        comp_matrix <- num_matrix[[1]]
      }
    }else{
      # A search we have
      match_names <- names(num_matrix) %>% .search_for_keywords(.pattern = .compensation)
      if(length(match_names) > 0){
        comp_matrix <- num_matrix[[match_names]]
      }else{
        warning("Your search word gave no match. No compensation applied")
        return(.flow_frame)
      }

    }
  }
  else{
    if(purrr::is_numeric(.compensation)){
      comp_matrix <- .compensation
    }else{
      stop("I don't know how to deal with what you gave me. Please check the documentation.")
    }
  }

  flowCore::compensate(.flow_frame, comp_matrix)
}

