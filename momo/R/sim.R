##' Simulate alpha
##'
##' @param ll ll
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



##' Simulate all data needed for fit.momo
##'
##' @param fit optional. momo fit
##'
##' @export
sim.momo <- function(fit = NULL,
                     xrange = c(0, 1),
                     yrange = c(0, 1),
                     dxdy = c(0.1, 0.1),
                     select = FALSE,
                     plot.land = FALSE,
                     keep.par = FALSE,
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
                     trange.rel = NULL,
                     xrange.rel = c(0.4,0.6),
                     yrange.rel = c(0.4,0.6),
                     verbose = TRUE){


    grid <- create.grid(xrange = xrange,
                        yrange = yrange,
                        dxdy = dxdy,
                        select = FALSE,
                        plot.land = FALSE,
                        keep.par = FALSE,
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
                   diagonal = diagonal)


    par.out <- list(alpha = array(sim.alpha(0.01, 0.04, n.alpha),
                              dim = c(n.alpha,1,1)))
    if(const.dif){
        par.out$beta <- array(log(runif(1, 0.01, 0.04)), dim = c(1,1,1))
    }else{
        par.out$beta <- array(log(runif(3, 0.01, 0.04)), dim = c(3,1,1))
    }
    par.out$logSdObsATS <- log(runif(1, 0.01, 0.04))
    if(!is.null(par)){
        for(i in 1:length(par)){
            par.out[names(par)[i]] <- par[names(par)[i]]
        }
    }
    par <- par.out

    if(is.null(knots.tax)){
        env.vals <- as.numeric(env[as.integer(cut(0.3,attr(grid,"xgr"))):
                                   as.integer(cut(0.7,attr(grid,"xgr"))),
                                   as.integer(cut(0.3,attr(grid,"ygr"))):
                                   as.integer(cut(0.7,attr(grid,"ygr"))),1])
        ## knots.tax <- matrix(quantile(env.vals, probs = c(0.1,0.5,0.9)),3,1)
        knots.tax <- matrix(seq(min(env.vals, na.rm = TRUE),
                                max(env.vals, na.rm = TRUE),
                                length.out = 5)[-c(1,5)],3,1)
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
                           knots.dif = knots.dif)

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
                    diagonal = FALSE){

    nx <- nrow(grid$celltable)
    ny <- ncol(grid$celltable)

    if(is.null(rho_s)){
        rho_s <- mean(attr(grid,"dxdy")) / 0.125
    }

    rescale.env <- function(env, zrange, i = NULL){
        env.range <- range(unlist(env), na.rm = TRUE)
        envi <- if(!is.null(i)) env[[i]] else env
        res <- (envi - env.range[1]) /
            (env.range[2] - env.range[1]) * (zrange[2] - zrange[1]) + zrange[1]
        return(res)
    }

    get.env <- function(nx, ny, sd, h, nu, rho, delta, matern, diagonal){

        ## Generate a random field
        rf <- rnorm(nrow(grid$xygrid), 0, sd = sd)

        ## GMRF
        Q <- get.precision.matrix(grid, h, nu, rho, delta, matern, diagonal)
        L <- chol(Q)
        S <- solve(L, rf)
        rf.smooth <- matrix(NA, nrow = nx, ncol = ny)
        rf.smooth[!is.na(grid$celltable)] <- S

        return(rf.smooth)
    }

    env0 <- env <- vector("list", nt)
    env0[[1]] <- get.env(nx, ny, sd, h, nu, rho_s, delta, matern, diagonal)
    if(nt > 1){
        for(i in 2:nt){
            env0[[i]] <- rho_t * env0[[i-1]] + sqrt(1 - rho_t^2) * get.env(nx, ny, sd, h, nu, rho_s, delta, matern, diagonal)
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

##' Simulate effort data
##' @export
sim.effort <- function(grid,
                       nt = 9,
                       rho_t = 0.85,
                       sd = 2,
                       h = 0.2,
                       nu = 2,
                       rho_s = 0.8,
                       delta = 0.1,
                       zrange = c(20,28),
                       matern = TRUE,
                       diagonal = FALSE){

    nx <- nrow(grid$celltable)
    ny <- ncol(grid$celltable)

    rescale.effort <- function(eff, zrange, i = NULL){
        eff.range <- range(unlist(eff))
        effi <- if(!is.null(i)) eff[[i]] else eff
        res <- (effi - eff.range[1]) /
            (eff.range[2] - eff.range[1]) * (zrange[2] - zrange[1]) + zrange[1]
        return(res)
    }

    get.effort <- function(nx, ny, sd, h, nu, rho, delta, matern, diagonal){

        ## Generate a random field
        rf <- rnorm(nx * ny, 0, sd = sd)

        ## GMRF
        Q <- get.precision.matrix(grid, h, nu, rho, delta, matern, diagonal)
        L <- chol(Q)
        S <- solve(L, rf)
        rf.smooth <- matrix(S, nrow = nx, ncol = ny)

        return(rf.smooth)
    }

    eff0 <- eff <- vector("list", nt)
    eff0[[1]] <- get.effort(nx, ny, sd, h, nu, rho_s, delta, matern, diagonal)
    if(nt > 1){
        for(i in 2:nt){
            eff0[[i]] <- rho_t * eff0[[i-1]] + sqrt(1 - rho_t^2) * get.effort(nx, ny, sd, h, nu, rho_s, delta, matern, diagonal)
        }
    }

    ## Rescale
    for(i in 1:nt){
        eff[[i]] <- rescale.effort(eff0, zrange, i)
    }

    eff <- prep.effort(eff)

    return(eff)
}


##' Simulate conventional tags
##' @export
sim.ctags <- function(grid,
                      par = NULL,
                      env = NULL,
                      n = 300,
                      effort = NULL,
                      trange = NULL,
                      trange.rel = NULL,
                      xrange.rel = NULL,
                      yrange.rel = NULL,
                      trange.rec = NULL,
                      by = 0.01,
                      knots.tax = NULL,
                      knots.dif = NULL,
                      funcs = NULL){

    flag.const.dif <- ifelse(nrow(par$beta) == 1, TRUE, FALSE)
    flag.effort <- ifelse(!is.null(effort), TRUE, FALSE)

    env <- check.that.list(env)

    ## Dimensions
    if(is.null(trange) && is.null(env)){
        trange <- c(0,1)
    }
    if(!is.null(env)){
        trange <- c(0, max(sapply(env, function(x) dim(x)[3])))
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
                           knots.tax = knots.tax,
                           knots.dif = knots.dif,
                           trange = trange)

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
        xy <- c(x0, y0)
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
            xy <- xy + move
            if(flag.effort){
                FF <- funcs$fish.mort(xy, t)
                M <- funcs$nat.mort(xy, t)
                psurv <- exp(-(FF + M) * by)
                state <- sample(1:3, 1, prob=c(psurv,
                                               FF/(FF+M)*(1-psurv),
                                               M/(FF+M)*(1-psurv)))
                if(state != 1) break
            }
            t <- t + by
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
        return(res)
    }

    res <- as.data.frame(t(replicate(n, doone(by = by))))
    ## res <- res[is.na(res$t1) |
    ##            res$t0 != res$t1,]

    return(res)
}


##' Simulate archival tags
##' @export
sim.atags <- function(grid,
                      par = NULL,
                      env = NULL,
                      n = 30,
                      effort = NULL,
                      trange = NULL,
                      trange.rel = NULL,
                      xrange.rel = NULL,
                      yrange.rel = NULL,
                      trange.rec = NULL,
                      by = 0.01,
                      knots.tax = NULL,
                      knots.dif = NULL,
                      funcs = NULL){

    flag.env <- ifelse(is.null(env), FALSE, TRUE)
    flag.const.dif <- ifelse(nrow(par$beta) == 1, TRUE, FALSE)
    flag.effort <- ifelse(!is.null(effort), TRUE, FALSE)

    env <- check.that.list(env)

    ## Dimensions
    if(is.null(trange) && is.null(env)){
        trange <- c(0,1)
    }
    if(!is.null(env)){
        trange <- c(0, max(sapply(env, function(x) dim(x)[3])))
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
                           trange = trange)
    conf <- def.conf(dat)


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
        xy <- c(x0, y0)
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
            xy <- xy + move
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
