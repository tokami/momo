
##' Simulate environmental fields
##' @export
sim.env <- function(grid,
                    nt = 9,
                    rho_t = 0.85,
                    sd = 2,
                    h = 0.2,
                    nu = 2,
                    rho_s = 0.8,
                    delta = 0.1,
                    zrange = c(20,28),
                    matern = TRUE,
                    diagonal = TRUE){

    nx <- nrow(grid$celltable)
    ny <- ncol(grid$celltable)

    rescale.env <- function(env, zrange, i = NULL){
        env.range <- range(unlist(env))
        envi <- if(!is.null(i)) env[[i]] else env
        res <- (envi - env.range[1]) /
            (env.range[2] - env.range[1]) * (zrange[2] - zrange[1]) + zrange[1]
        return(res)
    }

    get.env <- function(nx, ny, sd, h, nu, rho, delta, matern, diagonal){

        ## Generate a random field
        rf <- rnorm(nx * ny, 0, sd = sd)

        ## GMRF
        Q <- get.precision.matrix(grid, h, nu, rho, delta, matern, diagonal)
        L <- chol(Q)
        S <- solve(L, rf)
        rf.smooth <- matrix(S, nrow = nx, ncol = ny)

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
                       diagonal = TRUE){

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
                      n = 100,
                      effort = NULL,
                      trange = NULL,
                      trange.rel = NULL,
                      xrange.rel = NULL,
                      yrange.rel = NULL,
                      by = 0.01,
                      knots.tax = NULL,
                      knots.dif = NULL,
                      funcs = NULL){

    flag.const.dif <- ifelse(nrow(par$beta) == 1, TRUE, FALSE)
    flag.effort <- ifelse(!is.null(effort), TRUE, FALSE)

    env <- check.that.list(env)

    ## Dimensions
    if(is.null(trange) && is.null(env)){
        trange <- c(0, 1)
    }else if(!is.null(env)){
        trange <- c(0, max(sapply(env, function(x) dim(x)[3])))
    }
    if(is.null(trange.rel)){
        trange.rel <- trange
    }

    if(trange.rel[2] > trange[2]) stop("trange.rel[2] > trange[2]!")

    xrange <- attributes(grid)$xrange
    yrange <- attributes(grid)$yrange

    if(is.null(xrange.rel)){
        xrange.rel <- xrange
    }
    if(is.null(yrange.rel)){
        yrange.rel <- yrange
    }

    ## Setup default data and conf
    dat <- setup.momo.data(grid, env, trange = trange)
    conf <- def.conf(dat)
    if(!is.null(knots.tax)){
        conf$knots.tax <- knots.tax
    }
    if(!is.null(knots.dif)){
        conf$knots.dif <- knots.dif
    }

    ## Parameters
    par <- get.sim.par(par)

    ## Functions
    funcs <- get.sim.funcs(funcs, dat, conf, env, par)

    doone <- function(by = 0.01){
        t0 <- round(runif(1, trange.rel[1], trange.rel[2])/by) * by
        if(!flag.effort){
            t1 <- round(runif(1, t0, trange[2])/by) * by
        }else{
            t1 <- trange[2]
        }
        x0 <- runif(1, xrange.rel[1], xrange.rel[2])
        y0 <- runif(1, yrange.rel[1], yrange.rel[2])
        xy <- c(x0, y0)
        t <- t0
        state <- 1
        while(t < t1){
            xy <- xy +
                funcs$tax(xy, t) * by +
                funcs$adv(xy, t) * by +
                rnorm(2, mean = 0,
                      sd = sqrt(2 * exp(funcs$dif(xy, t)) * by))
            FF <- funcs$fish.mort(xy, t)
            M <- funcs$nat.mort(xy, t)
            psurv <- exp(-(FF + M) * by)
            state <- sample(1:3, 1, prob=c(psurv,
                                           FF/(FF+M)*(1-psurv),
                                           M/(FF+M)*(1-psurv)))
            if(flag.effort && state != 1) break
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
                      by = 0.01,
                      knots.tax = NULL,
                      knots.dif = NULL,
                      funcs = NULL){

    flag.env <- ifelse(is.null(env), FALSE, TRUE)
    flag.const.dif <- ifelse(nrow(par$beta) == 1, TRUE, FALSE)
    flag.effort <- ifelse(!is.null(effort), TRUE, FALSE)

    env <- check.that.list(env)

    if(is.null(trange) && is.null(env)){
        trange <- c(0,1)
    }else if(!is.null(env)){
        trange <- c(0,max(sapply(env, function(x) dim(x)[3])))
    }
    if(is.null(trange.rel)){
        trange.rel <- trange
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
    dat <- setup.momo.data(grid, env, trange = trange)
    conf <- def.conf(dat)
    if(!is.null(knots.tax)){
        conf$knots.tax <- knots.tax
    }
    if(!is.null(knots.dif)){
        conf$knots.dif <- knots.dif
    }

    ## Parameters
    par <- get.sim.par(par)

    ## Functions
    funcs <- get.sim.funcs(funcs, dat, conf, env, par)

    dooneAT <- function(by = 0.01, sdObs = 0.05){
        t0 <- round(runif(1, trange.rel[1], trange.rel[2])/by) * by
        if(!flag.effort){
            t1 <- round(runif(1, t0, trange[2])/by) * by
        }else{
            t1 <- trange[2]
        }
        x0 <- runif(1, xrange.rel[1], xrange.rel[2])
        y0 <- runif(1, yrange.rel[1], yrange.rel[2])
        xy <- c(x0, y0)
        ret <- matrix(NA, nrow = 1, ncol = 3)
        ret[1,] <- c(t0, x0, y0)
        t <- t0 + by
        state <- 1
        while(t < t1){
            xy <- xy +
                funcs$tax(xy, t) * by +
                funcs$adv(xy, t) * by +
                rnorm(2, mean = 0,
                      sd = sqrt(2 * exp(funcs$dif(xy, t)) * by))
            FF <- funcs$fish.mort(xy, t)
            M <- funcs$nat.mort(xy, t)
            psurv <- exp(-(FF + M) * by)
            state <- sample(1:3, 1, prob = c(psurv,
                                           FF/(FF+M)*(1-psurv),
                                           M/(FF+M)*(1-psurv)))
            ret <- rbind(ret, c(t, xy[1], xy[2]))
            if(flag.effort && state != 1) break
            t <- t + by
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

    res <- replicate(n, dooneAT())
    res <- res[lapply(res, nrow) > 3]

    return(res)
}
