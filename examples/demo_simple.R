## Simple demo of momo with simulated data
########################

## Install required packages (requires R (>= 4.4)):
## install.packages(c("Rcpp","RcppEigen","remotes"))
## remotes::install_github("kaskr/adcomp/TMB")
## remotes::install_github("kaskr/RTMB/RTMB")
## remotes::install_github("calbertsen/RTMBconvenience/RTMBconvenience")
## remotes::install_github("tokami/momo/momo",
##                         auth_token = "...") ## HERE: enter your personal auth_token to install private repo from github


## Load momo
require(momo)


## Simulate a simple grid
grid <- create.grid()


## Simulate an environmental field
set.seed(410)
env <- sim.env(grid, nt = 1, rho_s = 0.6, nu = 4)

## Plot env field
plotmomo.env(env)


## Define some parameters for simulation
par.true <- list(
    alpha = matrix(c(0, 0.03, 0.05), 3, 1),
    beta = matrix(log(0.001), 1, 1))


## Simulate some conventional tags
ctags <- sim.ctags(grid, par = par.true, env = env, n = 100)

## Plot tags
plotmomo.ctags(ctags)


## Setup data
dat <- setup.momo.data(grid = grid,
                       env = env,
                       ctags = ctags)


## Get default configurations
conf <- def.conf(dat)


## Get default parameters and inital values
par <- def.par(dat, conf)


## Fit momo (takes around 30 seconds)
fit1 <- fit.momo(dat, conf, par)


## Plot results
plotmomo.pref(fit1, par = par.true)



## Run matrix exponential approach on same data
conf$use.expm <- TRUE


## Get default parameters and inital values
par <- def.par(dat, conf)


## Fit momo (takes around 1 minute)
fit2 <- fit.momo(dat, conf, par)


## Plot results
plotmomo.pref(fit2, par = par.true)
