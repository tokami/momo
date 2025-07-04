## ----ReaddataLoadLibraries, message=FALSE, include=FALSE, echo=FALSE----------
knitr::opts_chunk$set(echo = TRUE,
                      cache = FALSE,
                      warning = FALSE,
                      eval = TRUE,
                      error = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      include = TRUE,
                      collapse = TRUE,
                      comment = "#>",
                      fig.show = "hold",
                      fig.width=8, fig.height=7)

## ----eval = FALSE-------------------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("tokami/momo")

## ----eval = FALSE-------------------------------------------------------------
# remotes::install_github("tokami/momo", ref = "dev")

## ----eval = TRUE, echo = FALSE, fig.alt = ""----------------------------------
library(momo)

## ----eval = TRUE, echo = FALSE, fig.alt = ""----------------------------------
# Simulate a data set with tagging data
sim <- sim.momo(n.ctags = 200, n.atags = 20)

## ----eval = TRUE, echo = FALSE, fig.alt = ""----------------------------------
# Fit the movement model
fit <- fit.momo(sim)

## ----eval = TRUE, echo = FALSE, fig.alt = "Simulated data and model predictions"----
# Plot the results
plotmomo.compare(sim = sim, fit = fit, bg = "white")

