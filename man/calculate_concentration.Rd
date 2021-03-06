% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_conc.R
\name{calculate_concentration}
\alias{calculate_concentration}
\title{Calculate concentration.}
\usage{
calculate_concentration(.data, .model, .parameter = "FL2.H",
  .value = "Calc.conc")
}
\arguments{
\item{.data}{A tidy data.frame.}

\item{.model}{An object of class \code{drc} with the fitted dose-response model.}

\item{.parameter}{A character giving the name of column(s) where populations
are identified.}

\item{.value}{A character giving the name of the column to store the calculated concentration}
}
\value{
The \code{.data} with the calculated concentration and error added in two columns.
}
\description{
Calculate the concentration in a sample
}
\examples{

library(beadplexr)
library(drc)
data(ryegrass)

ryegrass_m <-
  fit_standard_curve(.data = ryegrass,
                     .parameter = "rootl",
                     .concentration = "conc")

sample_data <-
  calculate_concentration(.data = ryegrass[sample(1:nrow(ryegrass), 5),],
                          .model = ryegrass_m,
                          .parameter = "rootl")
}
