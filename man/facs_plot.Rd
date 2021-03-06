% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisations.R
\name{facs_plot}
\alias{facs_plot}
\alias{facs_scatter}
\alias{facs_density2d}
\alias{facs_density1d}
\alias{facs_hexbin}
\title{Plot FACS data.}
\usage{
facs_plot(.data, .x = "FSC-A", .y = "SSC-A", .type = c("scatter",
  "density1d", "density2d", "hexbin"), ...)

facs_scatter(.data, .x = "FSC-A", .y = "SSC-A", .beads = NULL,
  .plot_distinct = TRUE)

facs_density2d(.data, .x = "FSC-A", .y = "SSC-A", .beads = NULL)

facs_density1d(.data, .x = "FSC-A", .beads = NULL)

facs_hexbin(.data, .x = "FSC-A", .y = "SSC-A", .bins = 75)
}
\arguments{
\item{.data}{The data to be plotted in a \code{data.frame}.}

\item{.x, .y}{Character vector with the column name for the variable to plot
on the x or y-axis.}

\item{.type}{Character vector giving the type of plot being used. Options are
one of \code{"scatter", "density", "hexbin"}.}

\item{...}{Arguments passed to the individual functions.}

\item{.beads}{Character vector to with the column name with identification of
beads. If used it will show up with the aesthetic 'colour'. Defaults to not
being used.}

\item{.plot_distinct}{Boolean to decide if only distinct events should be
plotted. If used, the number of data points might be greatly reduced which
could make for faster plotting. Defaults to TRUE.}

\item{.bins}{Numeric vector giving number of bins in both vertical and
horizontal directions. Set to 75 by default.}
}
\value{
A \code{ggplot}
}
\description{
Wrappers around building a ggplot with \code{geom_point},
\code{geom_density_2d}, and \code{geom_hex}.
}
\details{
These plot functions are meant to provide a quick way of viewing the FACS
data. For more control, use \code{ggplot2} directly.
}
\examples{
\dontrun{
library(beadplexr)

.file_name <- system.file("extdata", "K2-C07-A7.fcs",
                          package = "beadplexr")

.data <- read_fcs(.file_name = .file_name,
                  .filter = list("FSC-A" = c(2e5L, 6.3e5L),
                                 "SSC-A" = c(2e5, 1e6L)))
.data$bead_group <- ifelse(.data$`FSC-A` < 4e5L, "A", "B")

# Using facs_plot
facs_plot(.data, .type = "scatter")
facs_plot(.data, .type = "density1d")
facs_plot(.data, .type = "density2d")
facs_plot(.data, .type = "hexbin")

facs_plot(.data, .type = "scatter", .beads = "bead_group")
facs_plot(.data, .type = "density1d", .beads = "bead_group")
facs_plot(.data, .type = "hexbin", .bins = 50)

facs_plot(.data, .x = "FL2-H", .type = "scatter", .beads = "bead_group")

# Individual functions
facs_scatter(.data)

facs_scatter(.data, .beads = "bead_group", .plot_distinct = FALSE)
facs_scatter(.data, .beads = "bead_group")

facs_scatter(.data, .x = "FL2-H", .y = "FL6-H", .beads = "bead_group")

facs_density1d(.data)
facs_density1d(.data, .beads = "bead_group")

facs_density2d(.data)
facs_density2d(.data, .beads = "bead_group")

facs_hexbin(.data)
facs_hexbin(.data, .bins = 30)
}
}
