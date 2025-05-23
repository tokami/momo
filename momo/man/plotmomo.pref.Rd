% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotmomo.pref}
\alias{plotmomo.pref}
\title{Plot habitat preference function}
\usage{
plotmomo.pref(
  x,
  type = "taxis",
  select = NULL,
  main = NULL,
  cols = momo.cols(10),
  lwd = 1,
  ci = 0.95,
  keep.gpar = FALSE,
  add = FALSE,
  xlab = "Covariate",
  ylab = "Preference",
  bg = NULL,
  ylim = NULL,
  xlim = NULL,
  return.limits = FALSE,
  data.range = FALSE,
  asp = 2,
  ...
)
}
\arguments{
\item{x}{a list of class \code{momo.sim} or \code{momo.fit}.}

\item{type}{the process that should be plotted, either "taxis" (default) or
"diffusion".}

\item{main}{a main title for the plot. Default: "Preference".}

\item{lwd}{line width. Default: \code{1}.}

\item{ci}{level for confidence intervals. Default: \code{0.95}.}

\item{keep.gpar}{logical; If \code{TRUE}, do not overwrite the graphical
parameters. Default: \code{FALSE}.}

\item{add}{logical; if \code{TRUE}, no new graphical window is created and tags
are added to existing graphical window. Default: \code{FALSE}.}

\item{xlab}{a label for the x axis. Default: "Covariate".}

\item{ylab}{a label for the y axis. Default: "Preference".}

\item{bg}{background color of the graphical window. By default (\code{NULL}),
background is transparent.}

\item{ylim}{y limits}

\item{xlim}{x limits}

\item{return.limits}{logical; If \code{TRUE}, no graph is created but x and y
limits are returned. Default: \code{FALSE}.}

\item{...}{additional arguments for the function \link{plot}.}

\item{col}{color of line. Default: "black".}
}
\value{
Nothing.
}
\description{
Plot the habitat preference function as a function of
any covariate.
}
