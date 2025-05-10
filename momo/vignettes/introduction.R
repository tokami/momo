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

## ----echo = TRUE, include = TRUE, eval = FALSE--------------------------------
# install.packages("remotes")
# remotes::install_github("tokami/momo/momo")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# install.packages("remotes")
# remotes::install_github("tokami/momo/momo", ref = "dev")

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
library(momo)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
data(skjepo)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
## Load ctags
ctags <- prep.ctags(skjepo.ctags,
                    names = c("date_time","date_caught",
                              "rel_lon","recap_lon",
                              "rel_lat","recap_lat"),
                    origin = "1899-12-30",
                    speed.limit = 200)
plotmomo.ctags(ctags, plot.land = TRUE)

## Load atags
atags <- prep.atags(skjepo.atags,
                    names = c("time","mptlon","mptlat"),
                    origin = "1899-12-30")
plotmomo.atags(atags, plot.land = TRUE)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
# ## Create a grid based on tagging data
# grid <- create.grid(c(-150, -70), c(-30, 35),
#                     dxdy = c(10,10),
#                     select = 2, plot.land = TRUE)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
grid <- create.grid(grid = skjepo.grid, dxdy = c(10,10))
plotmomo.grid(grid, plot.land = TRUE)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
## Env data
env <- prep.env(skjepo.env)
plotmomo.env(env[,,1:4], plot.land = TRUE,
             xlab = "lon", ylab = "lat")

## ----eval=TRUE, echo=FALSE----------------------------------------------------
ctags <- ctags[sample(1:nrow(ctags),100),]

## ----eval=TRUE, echo=FALSE----------------------------------------------------
## Combine and check data
dat <- setup.momo.data(grid = grid,
                       env = env,
                       ctags = ctags,
                       atags = atags)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
## Default configurations
conf <- def.conf(dat)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
## Default parameters
par <- def.par(dat, conf)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
## Fitting movement model
fit <- fit.momo(dat, conf, par,
                do.sdreport = FALSE, ## only for vignette
                verbose = TRUE)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
fit$opt$par

## ----eval=TRUE, echo=FALSE----------------------------------------------------
plotmomo.compare(sim = skjepo, fit = fit)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
tmp <- get.release.events(dat,
                          grid = create.grid(dat$xrange, dat$yrange, c(1,1)),
                          time.cont = seq(dat$trange[1], dat$trange[2],
                                          1/(52*diff(dat$trange))))
dat$rel.events <- tmp$rel.events
dat$ctags$rel.event <- tmp$idx

## ----eval=FALSE, echo=FALSE---------------------------------------------------
# conf <- def.conf(dat)
# conf$use.expm <- TRUE

## ----eval=FALSE, echo=FALSE---------------------------------------------------
# map <- def.map(dat, conf, par)

