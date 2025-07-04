##' Rescale an env field
##'
##' @description Rescale an env field
##'
##' @param env field
##' @param zrange range
##' @param i NULL
##'
##' @return Rescaled env field
rescale.env <- function(env, zrange, i = NULL){
    env.range <- range(unlist(env), na.rm = TRUE)
    envi <- if(!is.null(i)) env[[i]] else env
    res <- (envi - env.range[1]) /
        (env.range[2] - env.range[1]) * (zrange[2] - zrange[1]) + zrange[1]
    return(res)
}


##' Simulate a simple data set
##'
##' @description `sim.momo` allows to quickly simulate all data sets and
##'     components required to fit the movement model, [fit.momo].
##'
##' @param fit optional; allows to provide a list of class `momo.fit` to
##'     simulate based on the fitted object. By default (`NULL`), no fitted
##'     object is used.
##' @param xrange range of the x-dimension of spatial domain. Default: `c(0,1)`.
##' @param yrange range of the y-dimension of spatial domain. Default: `c(0,1)`.
##' @param dxdy resolution of grid in x and y direction. Default: `c(0.1,0.1)`.
##' @param select logical; if `TRUE`, allows to select cells in spatial grid.
##'     Default: `FALSE`.
##' @param plot.land logical; If `TRUE`, plot land masses using the function
##'     [maps::map]. Default: `FALSE`.
##' @param keep.gpar logical; If `TRUE`, do not overwrite the graphical
##'     parameters. Default: `FALSE`.
##' @param trange range of model time. Default: `c(0,1)`.
##' @param dt resolution of model time steps. Default: `0.1`.
##' @param nt number of environmental fields. Default: `1`.
##' @param rho_t autocorrelation coefficient of temporal autocorrelation for
##'     environmental fields. Default: `0.85`.
##' @param sd standard deviation for simulation of environmental fields.
##'     Default: `2`.
##' @param h parameter of the matern covariance structure. Default: `0.2`.
##' @param nu parameter of the matern covariance structure. Default: `2`.
##' @param rho_s parameter of the matern covariance structure. Default: `0.8`.
##' @param delta parameter of the matern covariance structure. Default: `0.1`.
##' @param zrange range of the environmental covariate. Default: `c(20, 28)`.
##' @param matern logical; if `TRUE` (default), matern covariance structure is
##'     used for simulation of environmental fields.
##' @param diagonal logical; if `TRUE`, diagonal neighbours are considered in
##'     neighbouring structure. Default: `FALSE`.
##' @param par optional; allows to specify all or some of the parameters used
##'     for simulation. Default: `NULL`.
##' @param n.alpha number of knots for taxis. Default: `3`.
##' @param const.dif logical; if `TRUE` (default), constant diffusion is
##'     assumed.
##' @param knots.tax optional; allows to specify knots for the taxis component.
##'     Default: `NULL`.
##' @param knots.dif optional; allows to specify knots for the diffusion
##'     component. Default: `NULL`.
##' @param correct.peclet logical; if `TRUE`, the diffusion is increased until
##'     the peclet number is fulfilled. Default: `FALSE`.
##' @param use.ctags logical; if `TRUE` (default), mark-recapture tags are
##'     simulated.
##' @param n.ctags number of mark-recapture tags. Default: `500`.
##' @param use.atags logical; if `TRUE` (default), archival tags are simulated.
##' @param n.atags number of archival tags. Default: `50`.
##' @param trange.rel optional; allows to specify a specific part of the time
##'     series in which tags are released. Default: `NULL`.
##' @param xrange.rel optional; allows to specify a specific part of the x range
##'     of the spatial domain in which tags are released. Default: `c(0.4,0.6)`.
##' @param yrange.rel optional; allows to specify a specific part of the y range
##'     of the spatial domain in which tags are released. Default: `c(0.4,0.6)`.
##' @param verbose If `TRUE`, print information to console. Default: `TRUE`.
##'
##' @return A list of class `momo.sim`.
##'
##' @examples
##' sim <- sim.momo()
##'
##' @export
sim.momo <- function(fit = NULL,
                     xrange = c(0,1),
                     yrange = c(0,1),
                     dxdy = c(0.1, 0.1),
                     select = FALSE,
                     plot.land = FALSE,
                     keep.gpar = FALSE,
                     trange = c(0,1),
                     dt = 0.1,
                     nt = 1,
                     rho_t = 0.85, sd = 2, h = 0.2, nu = 2,
                     rho_s = 0.8, delta = 0.1,
                     zrange = c(20, 28), matern = TRUE,
                     diagonal = FALSE,
                     par = NULL,
                     n.alpha = 3,
                     const.dif = TRUE,
                     knots.tax = NULL,
                     knots.dif = NULL,
                     correct.peclet = FALSE,
                     use.ctags = TRUE,
                     n.ctags = 500,
                     use.atags = TRUE,
                     n.atags = 50,
                     trange.rel = c(0,0.2),
                     xrange.rel = c(0.4,0.6),
                     yrange.rel = c(0.4,0.6),
                     verbose = TRUE){


    grid <- create.grid(xrange = xrange,
                        yrange = yrange,
                        dxdy = dxdy,
                        select = select,
                        plot.land = plot.land,
                        keep.gpar = keep.gpar,
                        verbose = verbose,
                        fit = fit)

    if(is.null(xrange.rel)){
        xrange.rel <- attr(grid, "xrange")
    }
    if(is.null(yrange.rel)){
        yrange.rel <- attr(grid, "yrange")
    }

    env <- sim.env(grid, nt = nt,
                   rho_t = rho_t, sd = sd, h = h, nu = nu,
                   rho_s = rho_s, delta = delta,
                   zrange = zrange,
                   matern = matern,
                   diagonal = diagonal,
                   sim.buffer = TRUE)

    attributes(env)$dimnames[[3]] <- as.character(seq(trange[1], trange[2], length.out = nt))

    par.out <- list(alpha = array(sim.alpha(0.01, 0.2, n.alpha),  ## 0.01-0.04
                              dim = c(n.alpha,1,1)))
    if(const.dif){
        par.out$beta <- array(log(runif(1, 0.01, 0.1)), dim = c(1,1,1)) ## 0.02 - 0.1
    }else{
        par.out$beta <- array(log(runif(3, 0.01, 0.1)), dim = c(3,1,1))
    }
    par.out$logSdObsATS <- log(runif(1, 0.01, 0.1))
    par.out$logSdObsSTS <- log(runif(1, 0.01, 0.1))  ## 0.005 - 0.08
    if(!is.null(par)){
        for(i in 1:length(par)){
            par.out[names(par)[i]] <- par[names(par)[i]]
        }
    }
    par <- par.out

    if(is.null(knots.tax)){
        env.vals <- as.numeric(env[as.integer(cut(0.2,attr(grid,"xgr"))):
                                   as.integer(cut(0.8,attr(grid,"xgr"))),
                                   as.integer(cut(0.2,attr(grid,"ygr"))):
                                   as.integer(cut(0.8,attr(grid,"ygr"))),1])
        ## knots.tax <- matrix(quantile(env.vals, probs = c(0.1,0.5,0.9)),3,1)
        knots.tax <- matrix(seq(min(env.vals, na.rm = TRUE),
                                max(env.vals, na.rm = TRUE),
                                length.out = 7)[c(2,4,6)],3,1)
    }
    if(is.null(knots.dif)){
        if(const.dif){
            knots.dif <- NULL
        }else{
            knots.dif <- knots.tax
        }
    }

    ## Peclet number criterion
    if(correct.peclet){
        peclet <- get.peclet(grid, env, par,
                                 knots.tax = knots.tax,
                                 knots.dif = knots.dif)
        ## quant <- quantile(abs(unlist(peclet)), 0.9)
        quant <- median(abs(unlist(peclet)))
        maxi <- 1
        while(quant > 2 && maxi < 20){
            par$beta <- par$beta + log(1.5)
            peclet <- get.peclet(grid, env, par,
                                 knots.tax = knots.tax,
                                 knots.dif = knots.dif)
            ## quant <- quantile(abs(unlist(peclet)), 0.9)
            quant <- median(abs(unlist(peclet)))
            maxi <- maxi + 1
        }
    }

    if(use.ctags){
        ctags <- sim.ctags(grid, par, env, n.ctags,
                           trange = trange,
                           dt = dt,
                           trange.rel = trange.rel,
                           xrange.rel = xrange.rel,
                           yrange.rel = yrange.rel,
                           knots.tax = knots.tax,
                           knots.dif = knots.dif
                           )
    }else{
        ctags <- NULL
    }

    if(use.atags){
        atags <- sim.atags(grid, par, env, n.atags,
                           trange = trange,
                           dt = dt,
                           trange.rel = trange.rel,
                           xrange.rel = xrange.rel,
                           yrange.rel = yrange.rel,
                           knots.tax = knots.tax,
                           knots.dif = knots.dif
                           )
    }else{
        atags <- NULL
    }

    effort <- NULL

    dat <- setup.momo.data(grid = grid,
                           env = env,
                           ctags = ctags,
                           atags = atags,
                           effort = effort,
                           knots.tax = knots.tax,
                           knots.dif = knots.dif,
                           trange = trange,
                           dt = dt)

    res <- list()
    res$grid <- grid
    res$env <- env
    res$par.sim <- par
    res$ctags <- ctags
    res$atags <- atags
    res$dat <- dat

    res$conf <- def.conf(dat)
    res$par <- def.par(dat, res$conf)
    res$map <- def.map(dat, res$conf, res$par)

    res <- add.class(res, "momo.sim")

    return(res)
}



##' Simulate environmental fields
##'
##' @description `sim.env` allows to simulate environmental fields.
##'
##' @param grid a grid.
##' @param nt number of environmental fields. Default: `1`.
##' @param rho_t autocorrelation coefficient of temporal autocorrelation for
##'     environmental fields. Default: `0.85`.
##' @param sd standard deviation for simulation of environmental fields.
##'     Default: `2`.
##' @param h parameter of the matern covariance structure. Default: `0.2`.
##' @param nu parameter of the matern covariance structure. Default: `2`.
##' @param rho_s parameter of the matern covariance structure. Default: `0.8`.
##' @param delta parameter of the matern covariance structure. Default: `0.1`.
##' @param zrange range of the environmental covariate. Default: `c(20, 28)`.
##' @param matern logical; if `TRUE` (default), matern covariance structure is
##'     used for simulation of environmental fields.
##' @param diagonal logical; if `TRUE`, diagonal neighbours are considered in
##'     neighbouring structure. Default: `FALSE`.
##' @param sim.buffer logical; if `TRUE` extends the provided grid so that the
##'     2D interpolation function does not return `NaN` close to the boundary.
##' @param verbose If `TRUE`, print information to console. Default: `TRUE`.
##'
##' @return A list with simulated environmental fields.
##'
##' @examples
##' env <- sim.env()
##'
##' @export
sim.env <- function(grid,
                    nt = 9,
                    rho_t = 0.85,
                    sd = 2,
                    h = 0.2,
                    nu = 2,
                    rho_s = NULL, ## 0.8,
                    delta = 0.1,
                    zrange = c(20,28),
                    matern = TRUE,
                    diagonal = FALSE,
                    sim.buffer = TRUE,
                    verbose = TRUE){

    if(sim.buffer){
        grid <- add.buffer(grid)
    }

    if(is.null(rho_s)){
        rho_s <- mean(attr(grid,"dxdy")) / 0.125
    }

    get.env.one <- function(grid, sd, h, nu, rho, delta, matern, diagonal){

        ## Generate a random field
        rf <- rnorm(nrow(grid$xygrid), 0, sd = sd)

        ## GMRF
        Q <- get.precision.matrix(grid, h, nu, rho, delta, matern, diagonal)
        L <- chol(Q)
        S <- solve(L, rf)
        rf.smooth <- matrix(NA,
                            nrow = attr(grid,"nx"),
                            ncol = attr(grid,"ny"))
        rf.smooth[!is.na(grid$celltable)] <- S

        return(rf.smooth)
    }

    env0 <- env <- vector("list", nt)
    env0[[1]] <- get.env.one(grid, sd, h, nu, rho_s, delta, matern, diagonal)
    if(nt > 1){
        for(i in 2:nt){
            env0[[i]] <- rho_t * env0[[i-1]] + sqrt(1 - rho_t^2) *
                get.env.one(grid, sd, h, nu, rho_s, delta, matern, diagonal)
        }
    }

    ## Rescale
    for(i in 1:nt){
        env[[i]] <- rescale.env(env0, zrange, i)
        rownames(env[[i]]) <- attr(grid, "xcen")
        colnames(env[[i]]) <- attr(grid, "ycen")
    }

    env <- prep.env(env)

    return(env)
}


##' Simulate mark-recapture tags
##'
##' @description `sim.ctags` allows to simulate mark-recapture tags.
##'
##' @param grid a grid.
##' @param par optional; allows to specify all or some of the parameters used
##'     for simulation. Default: `NULL`.
##' @param env optional; allows to specify environmental fields that are used
##'     for simulation. Default: `NULL`.
##' @param n number of tags. Default: `300`.
##' @param effort optional; allows to specify effort fields that are used for
##'     simulation. Default: `NULL`.
##' @param trange range model time. Default: `NULL`.
##' @param dt optional; resolution of model time steps. If `NULL` (default), the
##'     default time step of [setup.momo.data] is used (`0.1`).
##' @param trange.rel optional; allows to specify a specific part of the time
##'     series in which tags are released. Default: `NULL`.
##' @param xrange.rel optional; allows to specify a specific part of the x range
##'     of the spatial domain in which tags are released. Default: `c(0.4,0.6)`.
##' @param yrange.rel optional; allows to specify a specific part of the y range
##'     of the spatial domain in which tags are released. Default: `c(0.4,0.6)`.
##' @param trange.rec optional; allows to specify a specific part of the time
##'     series in which tags are recaptured. Default: `NULL`.
##' @param by time step of the tag simulation model. Default: `0.05`.
##' @param knots.tax optional; allows to specify knots for the taxis component.
##'     Default: `NULL`.
##' @param knots.dif optional; allows to specify knots for the diffusion
##'     component. Default: `NULL`.
##' @param funcs optional; allows to specify movement functions. If `NULL`
##'     (default), the function [get.sim.funcs] is used.
##' @param save.tracks logical; if `TRUE`, full tracks with all
##'     intermediate steps are saved. Default: `FALSE`.
##' @param verbose If `TRUE`, print information to console. Default: `TRUE`.
##'
##' @return A data frame with simulate mark-recaptured tags.
##'
##' @examples
##' data(skjepo)
##'
##' ctags <- sim.ctags(skjepo$grid)
##'
##' @export
sim.ctags <- function(grid,
                      par = NULL,
                      env = NULL,
                      n = 300,
                      effort = NULL,
                      trange = NULL,
                      dt = NULL,
                      trange.rel = NULL,
                      xrange.rel = NULL,
                      yrange.rel = NULL,
                      trange.rec = NULL,
                      by = 0.05,
                      knots.tax = NULL,
                      knots.dif = NULL,
                      funcs = NULL,
                      save.tracks = FALSE,
                      verbose = TRUE){

    flag.const.dif <- ifelse(nrow(par$beta) == 1, TRUE, FALSE)
    flag.effort <- ifelse(!is.null(effort), TRUE, FALSE)

    env <- check.that.list(env)

    ## Dimensions
    if(is.null(trange) && is.null(env)){
        trange <- c(0,1)
    }else if(is.null(trange) && !is.null(env)){
        if(!is.null(attributes(env[[1]])$dimnames[[3]])){
            ## stop("Implement to simulate tags from dates in env, but dates could have any format.")
            trange <- range(as.numeric(attributes(env[[1]])$dimnames[[3]]))
        }else{
            trange <- c(0, max(sapply(env, function(x) dim(x)[3])))
            if(verbose) warning("Time range extracted from number of environmental fields.")
        }
    }

    if(is.null(trange.rel)){
        trange.rel <- trange
    }
    if(trange.rel[2] > trange[2]) stop("trange.rel[2] > trange[2]!")
    if(is.null(trange.rec)){
        trange.rec <- trange
    }

    xrange <- attributes(grid)$xrange
    yrange <- attributes(grid)$yrange

    if(is.null(xrange.rel)){
        xrange.rel <- xrange
    }
    if(is.null(yrange.rel)){
        yrange.rel <- yrange
    }

    ## Checks
    trange.rel <- sort(trange.rel)
    trange.rec <- sort(trange.rec)
    xrange.rel <- sort(xrange.rel)
    yrange.rel <- sort(yrange.rel)

    ## Setup default data and conf
    dat <- setup.momo.data(grid, env,
                           effort = effort,
                           knots.tax = knots.tax,
                           knots.dif = knots.dif,
                           trange = trange,
                           dt = dt)

    conf <- def.conf(dat)

    ## Parameters
    par <- get.sim.par(par)

    ## Functions
    funcs <- get.sim.funcs(funcs, dat, conf, env, par)

    doone <- function(by = 0.01){
        t0 <- round(runif(1, trange.rel[1],
                          ifelse(trange.rel[2]-by < trange.rel[1],
                                 trange.rel[1], trange.rel[2]-by))/by) * by
        if(!flag.effort){
            trec1 <- ifelse(trange.rec[1] < t0, t0+by, trange.rec[1])
            t1 <- round(runif(1, trec1,
                              ifelse(trange.rec[2]-by < trec1,
                                 trec1, trange.rec[2]-by))/by) * by
        }else{
            t1 <- trange[2]-by
        }
        x0 <- runif(1, xrange.rel[1], xrange.rel[2])
        y0 <- runif(1, yrange.rel[1], yrange.rel[2])
        xy <- matrix(c(x0, y0),1,2)
        t <- t0
        state <- ifelse(flag.effort, 1, 2)
        if(save.tracks){
            track <- data.frame(t = t, x = xy[1], y = xy[2], state = state)
        }else{
            track <- NULL
        }
        while((t + by) < t1){
            move <- c(0,0)
            if(conf$use.taxis){
                move <- move + funcs$tax(xy, t) * by
            }
            if(conf$use.advection){
                move <- move + funcs$adv(xy, t) * by
            }
            ## Diffusion
            move <- move + rnorm(2, mean = 0,
                                 sd = sqrt(2 * exp(funcs$dif(xy, t)) * by))
            if(conf$use.boundaries){
                move <- move * funcs$bound(xy, t)
            }
            ## Move
            xy <- xy + move
            if(flag.effort){
                FF <- funcs$fish.mort(xy, t)
                M <- funcs$nat.mort(xy, t)
                psurv <- exp(-(FF + M) * by)
                state <- sample(1:3, 1, prob=c(psurv,
                                               FF/(FF+M)*(1-psurv),
                                               M/(FF+M)*(1-psurv)))
                if(state != 1){
                    track$state[nrow(track)] <- state
                    break
                }
            }
            t <- t + by
            if(save.tracks){
                track <- rbind(track,
                               c(t = t, x = xy[1], y = xy[2], state = state))
            }
        }
        if(!flag.effort || state == 2){
            t1 <- t
            x1 <- xy[1]
            y1 <- xy[2]
        }else if(state == 3){
            t1 <- t
            x1 <- NA
            y1 <- NA
        }else{
            t1 <- NA
            x1 <- NA
            y1 <- NA
        }
        res <- c(t0 = t0, t1 = t1, x0 = x0, x1 = x1, y0 = y0, y1 = y1)
        ret <- list(res = res,
                    track = track)
        return(ret)
    }

    res.list <- track.list <- vector("list", n)
    for(i in 1:n){
        tmp <- doone(by = by)
        res.list[[i]] <- tmp$res
        track.list[[i]] <- tmp$track
    }
    res <- as.data.frame(do.call(rbind, res.list))

    if(save.tracks){
        ret <- list(ctags = res,
                    tracks = track.list)
    }else{
        ret <- res
    }

    ## res <- as.data.frame(t(replicate(n, doone(by = by))))
    ## res <- res[is.na(res$t1) |
    ##            res$t0 != res$t1,]

    return(ret)
}


##' Simulate archival tags
##'
##' @description `sim.atags` allows to simulate archival tags.
##'
##' @param grid a grid.
##' @param par optional; allows to specify all or some of the parameters used
##'     for simulation. Default: `NULL`.
##' @param env optional; allows to specify environmental fields that are used
##'     for simulation. Default: `NULL`.
##' @param n number of tags. Default: `30`.
##' @param effort optional; allows to specify effort fields that are used for
##'     simulation. Default: `NULL`.
##' @param trange range model time. Default: `NULL`.
##' @param dt optional; resolution of model time steps. If `NULL` (default), the
##'     default time step of [setup.momo.data] is used (`0.1`).
##' @param trange.rel optional; allows to specify a specific part of the time
##'     series in which tags are released. Default: `NULL`.
##' @param xrange.rel optional; allows to specify a specific part of the x range
##'     of the spatial domain in which tags are released. Default: `c(0.4,0.6)`.
##' @param yrange.rel optional; allows to specify a specific part of the y range
##'     of the spatial domain in which tags are released. Default: `c(0.4,0.6)`.
##' @param trange.rec optional; allows to specify a specific part of the time
##'     series in which tags are recaptured. Default: `NULL`.
##' @param by time step of the tag simulation model. Default: `0.01`.
##' @param knots.tax optional; allows to specify knots for the taxis component.
##'     Default: `NULL`.
##' @param knots.dif optional; allows to specify knots for the diffusion
##'     component. Default: `NULL`.
##' @param funcs optional; allows to specify movement functions. If `NULL`
##'     (default), the function [get.sim.funcs] is used.
##' @param boundary.grid optional; allows to provide a list of class momo.grid
##'     with boundary effects. Default: `NULL`.
##' @param verbose If `TRUE`, print information to console. Default: `TRUE`.
##'
##' @return A data frame with simulate archival tags.
##'
##' @examples
##' data(skjepo)
##'
##' atags <- sim.atags(skjepo$grid)
##'
##' @export
sim.atags <- function(grid,
                      par = NULL,
                      env = NULL,
                      n = 30,
                      effort = NULL,
                      trange = NULL,
                      dt = NULL,
                      trange.rel = NULL,
                      xrange.rel = NULL,
                      yrange.rel = NULL,
                      trange.rec = NULL,
                      by = 0.01,
                      knots.tax = NULL,
                      knots.dif = NULL,
                      funcs = NULL,
                      boundary.grid = NULL,
                      verbose = TRUE){

    flag.env <- ifelse(is.null(env), FALSE, TRUE)
    flag.const.dif <- ifelse(nrow(par$beta) == 1, TRUE, FALSE)
    flag.effort <- ifelse(!is.null(effort), TRUE, FALSE)

    env <- check.that.list(env)

    ## Dimensions
    if(is.null(trange) && is.null(env)){
        trange <- c(0,1)
    }else if(is.null(trange) && !is.null(env)){
        if(!is.null(attributes(env[[1]])$dimnames[[3]])){
            ## stop("Implement to simulate tags from dates in env, but dates could have any format.")
            trange <- range(as.numeric(attributes(env[[1]])$dimnames[[3]]))
        }else{
            trange <- c(0, max(sapply(env, function(x) dim(x)[3])))
        }
    }

    if(is.null(trange.rel)){
        trange.rel <- trange
    }
    if(trange.rel[2] > trange[2]) stop("trange.rel[2] > trange[2]!")
    if(is.null(trange.rec)){
        trange.rec <- trange
    }

    xrange <- attributes(grid)$xrange
    yrange <- attributes(grid)$yrange

    if(is.null(xrange.rel)){
        xrange.rel <- xrange
    }
    if(is.null(yrange.rel)){
        yrange.rel <- yrange
    }

    ## Setup default data and conf
    dat <- setup.momo.data(grid, env,
                           knots.tax = knots.tax,
                           knots.dif = knots.dif,
                           trange = trange,
                           boundary.grid = boundary.grid,
                           dt = dt)

    conf <- def.conf(dat)


    ## Where is this coming from?
    ## conf$ienvS$tax <- matrix(c(1,2), byrow = TRUE, 3, ncol(conf$ienvS$tax))

    ## needed for montagus' harrier
    ## conf$ienvS$tax <- matrix(rep(c(1,2), length.out = ncol(conf$ienvS$tax)),
    ##                          byrow = TRUE, 1, ncol(conf$ienvS$tax))

    ## Parameters
    par <- get.sim.par(par)

    ## Functions
    funcs <- get.sim.funcs(funcs, dat, conf, env, par)

    dooneAT <- function(by = 0.01, sdObs = 0.05){
        t0 <- round(runif(1, trange.rel[1],
                          ifelse(trange.rel[2]-by < trange.rel[1],
                                 trange.rel[1], trange.rel[2]-by))/by) * by
        if(!flag.effort){
            trec1 <- ifelse(trange.rec[1] < t0, t0+by, trange.rec[1])
            t1 <- round(runif(1, trec1,
                              ifelse(trange.rec[2]-by < trec1,
                                 trec1, trange.rec[2]-by))/by) * by
        }else{
            t1 <- trange[2]-by
        }
        x0 <- runif(1, xrange.rel[1], xrange.rel[2])
        y0 <- runif(1, yrange.rel[1], yrange.rel[2])
        xy <- matrix(c(x0, y0),1,2)
        ret <- matrix(NA, nrow = 1, ncol = 3)
        ret[1,] <- c(t0, x0, y0)
        t <- t0
        state <- ifelse(flag.effort, 1, 2)
        while((t + by) < t1){
            move <- c(0,0)
            if(conf$use.taxis){
                move <- move + funcs$tax(xy, t) * by
            }
            if(conf$use.advection){
                move <- move + funcs$adv(xy, t) * by
            }
            ## Diffusion
            move <- move + rnorm(2, mean = 0,
                                 sd = sqrt(2 * exp(funcs$dif(xy, t)) * by))
            if(conf$use.boundaries){
                move <- move * funcs$bound(xy, t)
            }
            ## Move
            xy.new <- xy + move
            if(conf$use.boundaries){
                if(all(!is.na(xy.new)) &&
                   !is.na(funcs$bound(xy.new, t)) &&
                   funcs$bound(xy.new, t) <= 0){
                    xy.new <- xy + (-1 * move)
                    ## browser()
                    ## xy
                    ## points(xy[1,1], xy[1,2], col = 3)
                    ## points(xy[1,1]+move[1,2], xy[1,2]+move[1,2], col = 2)
                    ## points(xy.new[1,1], xy.new[1,2], col = 6)
                }
            }
            xy <- xy.new
            if(flag.effort){
                FF <- funcs$fish.mort(xy, t)
                M <- funcs$nat.mort(xy, t)
                psurv <- exp(-(FF + M) * by)
                state <- sample(1:3, 1, prob = c(psurv,
                                                 FF/(FF+M)*(1-psurv),
                                                 M/(FF+M)*(1-psurv)))
            }
            t <- t + by
            ret <- rbind(ret, c(t, xy[1], xy[2]))
            if(flag.effort && state != 1) break
        }
        if(nrow(ret) >= 3){
            ret[-c(1,nrow(ret)),c(2,3)] <- ret[-c(1,nrow(ret)),c(2,3)] +
                rnorm(2*(nrow(ret)-2), 0, sdObs)  ## noise in pos obs
        }
        if(flag.effort && state != 2){
            ret <- ret[1,,drop=FALSE]
        }
        colnames(ret) <- c("t", "x", "y")
        return(ret)
    }

    res <- replicate(n, dooneAT(by = by, sdObs = exp(par$logSdObsATS)))
##    res <- res[lapply(res, nrow) > 3]

    return(res)
}



##' Simulate alpha
##'
##' @description Simulate alpha parameter specifying taxis component.
##'
##' @param ll lower limit of [runif] to simulate parameters.
##' @param ul upper limit of [runif] to simulate parameters.
##' @param n number of alpha parameters. Default: `3`.
##' @param digits number of significant digits for the parameters. Default: `2`.
##'
##' @return A vector with simulate alpha values.
##'
sim.alpha <- function(ll, ul, n = 3, digits = 2){
    alphas <- rep(NA, n)
    alphas[1] <- 0
    if(n > 1){
        for(i in 2:n){
            alphas[i] <- alphas[i-1] + sample(c(-1,1), 1) * runif(1, ll, ul)
        }
    }
    return(signif(alphas, digits))
}
