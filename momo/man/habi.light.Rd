% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/habi.R
\name{habi.light}
\alias{habi.light}
\title{Light habitat function}
\usage{
habi.light(liv, XR, YR, ienv, time.cont)
}
\arguments{
\item{XR}{Limits of the 2-D fields in x direction.}

\item{YR}{Limits of the 2-D fields in y direction.}

\item{ienv}{Indicator matrix mapping each model time step to the time steps
of each field. The first dimension corresponds to the number of fields
and the second to the number of model time steps.}

\item{time.cont}{Vector with continuous model time steps.}

\item{FIELDS}{A list with 3-D array fields that are used for interpolation,
where the first 2 dimensions span the x and y direction of the spatial
field and the third dimension indicates the time dimension.}
}
\value{
List with functions.
}
\description{
Creates a light habitat class that is being used within
\emph{momo}.
}
