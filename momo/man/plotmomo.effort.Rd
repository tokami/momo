% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotmomo.effort}
\alias{plotmomo.effort}
\title{Plot effort fields}
\usage{
plotmomo.effort(
  x,
  select = NULL,
  main = "Effort fields",
  labels = TRUE,
  plot.land = FALSE,
  keep.gpar = FALSE,
  xlab = "x",
  ylab = "y",
  bg = NULL,
  plot.contour = TRUE,
  plot.legend = TRUE,
  ...
)
}
\arguments{
\item{x}{a grid or a list of class \code{momo.data}, \code{momo.sim}, or \code{momo.fit}.}

\item{select}{optional; allows to select specific time steps. Default:
\code{NULL}.}

\item{main}{a main title for the plot. Default: \verb{Environmental fields}.}

\item{labels}{logical; If \code{TRUE} (default), plot numbers in cells.}

\item{plot.land}{logical; If \code{TRUE}, plot land masses using the function
\link[maps:map]{maps::map}. Default: \code{FALSE}.}

\item{keep.gpar}{logical; If \code{TRUE}, do not overwrite the graphical
parameters. Default: \code{FALSE}.}

\item{xlab}{a label for the x axis. Default: "x".}

\item{ylab}{a label for the y axis. Default: "y".}

\item{bg}{background color of the graphical window. By default (\code{NULL}),
background is transparent.}

\item{plot.contour}{logical; If \code{TRUE} (default), add contours to plot.}

\item{plot.legend}{logical; If \code{TRUE} (default), add legend to plot.}

\item{...}{additional arguments for the function \link{plot}.}
}
\value{
Nothing.
}
\description{
Plot the effort fields of any momo object that contains
effort data.
}
