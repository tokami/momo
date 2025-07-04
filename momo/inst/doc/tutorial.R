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

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
library(momo)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
data(skjepo)

## ----eval=TRUE, echo=TRUE, fig.alt = "Mark-recapture and archival tags"-------
## Load ctags
ctags <- prep.ctags(skjepo.ctags,
                    names = c("date_time","date_caught",
                              "rel_lon","recap_lon",
                              "rel_lat","recap_lat"),
                    origin = "1899-12-30",
                    speed.limit = 200)
plotmomo.ctags(ctags, plot.land = TRUE, bg = "white")

## Load atags
atags <- prep.atags(skjepo.atags,
                    names = c("time","mptlon","mptlat"),
                    origin = "1899-12-30")
plotmomo.atags(atags, plot.land = TRUE, bg = "white")

## ----eval=FALSE, echo=TRUE, fig.alt = "Spatial grid"--------------------------
# ## Create a grid based on tagging data
# grid <- create.grid(c(-150, -70), c(-30, 35),
#                     dxdy = c(10,10),
#                     select = 2, plot.land = TRUE)

## ----eval=TRUE, echo=TRUE, fig.alt = "Spatial grid"---------------------------
grid <- create.grid(grid = skjepo.grid, dxdy = c(10,10))
plotmomo.grid(grid, plot.land = TRUE, bg = "white")

## ----eval=TRUE, echo=TRUE, fig.alt = "First four environmental fields"--------
## Env data
env <- prep.env(skjepo.env)
plotmomo.env(env[,,1:4], plot.land = TRUE,
             xlab = "lon", ylab = "lat", bg = "white")

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
ctags <- ctags[sample(1:nrow(ctags),100),]

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
## Combine and check data
dat <- setup.momo.data(grid = grid,
                       env = env,
                       ctags = ctags,
                       atags = atags)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
## Default configurations
conf <- def.conf(dat)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
## Default parameters
par <- def.par(dat, conf)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
## Fitting movement model
fit <- fit.momo(dat, conf, par,
                verbose = TRUE)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
fit$opt$par

## ----eval=TRUE, echo=TRUE, fig.alt = "Simulated data and model predictions"----
plotmomo.compare(sim = skjepo, fit = fit,
                 plot.land = TRUE,
                 cor.dif = 0.5,
                 bg = "white")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# tmp <- get.release.events(dat,
#                           grid = create.grid(dat$xrange, dat$yrange, c(1,1)),
#                           time.cont = seq(dat$trange[1], dat$trange[2],
#                                           1/(52*diff(dat$trange))))
# dat$rel.events <- tmp$rel.events
# dat$ctags$rel.event <- tmp$idx

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# conf <- def.conf(dat)
# conf$use.expm <- TRUE

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# map <- def.map(dat, conf, par)

