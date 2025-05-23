% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{get.peclet}
\alias{get.peclet}
\title{Get peclet number}
\usage{
get.peclet(grid, env, par, knots.tax = NULL, knots.dif = NULL)
}
\arguments{
\item{grid}{a grid object of class \code{momo.grid} as returned by the function
\link{create.grid}.}

\item{env}{a list with environmental covariates of class \code{momo.env} as
returned by the function \link{prep.env}.}

\item{par}{parameter list with initial values as produced by the function
\link{def.par}.}

\item{knots.tax}{knots for the taxis component. Default: \code{NULL}.}

\item{knots.dif}{knots for the diffusion component. Default: \code{NULL}.}
}
\value{
The peclet number.
}
\description{
\code{get.peclet} calculates the peclet number.
}
