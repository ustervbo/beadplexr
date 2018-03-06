#' LEGENDplex example data
#'
#' Data from a "Human Growth Factor Panel (13-plex)" LEGENDplex experiment, with
#' 8 controls and 4 human serum samples, all in duplicates. The beads were
#' measured on a CytoFLEX cytometer, and the fcs-files were processed using
#' [read_fcs()], with default settings.
#'
#' @docType data
#'
#' @usage data("lplex")
#'
#' @format A list with 18 elements. Each element is a data.frame about 5000 rows
#'   and 4 columns (the exact number varies a little due to the data
#'   acquisition):
#' \describe{
#'   \item{FSC-A}{The forward scatter parameter}
#'   \item{SSC-A}{The side scatter parameter}
#'   \item{FL6-H}{Intensity in the FL6 channel}
#'   \item{FL2-H}{Intensity in the FL2 channel}
#' }
#'
#' The list contains 8 standard samples in duplicates, and one serum sample,
#' also in duplicate. The names of each element have the format K3 (internal
#' panel shorthand), C:num: for standards and S:num: for serum sample, and a
#' number indicating the replicate (1 or 2).
#'
#' @source Ulrik Stervbo, 2016, Unpublished
#'
"lplex"

#' Simulated beadplex data
#'
#' Very simple, simulated multiplex data to demonstrate the clustering
#' functionality of the **beadplexr** package on CBA and MACSPlex assays.
#'
#' @docType data
#'
#' @usage data(simplex)
#'
#' @format A list with three elements. Each element is a `data.frame` of 3000 to 9000 rows. The exact format depends on the assay simulated:
#' \describe{
#'   \item{lplex}{
#'     Simulated LEGENDplex data. A single `data.frame` with the columns:
#'     \describe{
#'       \item{FSC}{The forward scatter parameter}
#'       \item{SSC}{The side scatter parameter}
#'       \item{APC}{Intensity in the APC channel}
#'       \item{PE}{Intensity in the PE channel}
#'     }
#'   }
#'   \item{mplex}{
#'     Simulated MACSPlex data. A single `data.frame` with the columns:
#'     \describe{
#'       \item{FSC}{The forward scatter parameter}
#'       \item{SSC}{The side scatter parameter}
#'       \item{FITC}{Intensity in the FITC channel}
#'       \item{PE}{Intensity in the PE channel}
#'       \item{APC}{Intensity in the PE channel}
#'     }
#'   }
#'   \item{cba}{
#'     Simulated CBA data. A single `data.frame` with the columns:
#'     \describe{
#'       \item{FSC}{The forward scatter parameter}
#'       \item{SSC}{The side scatter parameter}
#'       \item{APC}{Intensity in the APC channel}
#'       \item{APC-Cy7}{Intensity in the APC-Cy7 channel}
#'       \item{PE}{Intensity in the PE channel}
#'     }
#'   }
#' }
#'
#' @source Artificial
"simplex"
