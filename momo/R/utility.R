##' Create a grid
##'
##' @param xrange Range x
##' @param yrange Range y
##' @param dxdy dxdy
##' @param select.cells Select cells
##' @param plot.land Optional plot land (only when select.cells = TRUE)
##' @param verbose Print messages?
##'
##' @export
create.grid <- function(xrange = c(0,1),
                        yrange = c(0,1),
                        dxdy = c(0.1,0.1),
                        select = FALSE,
                        plot.land = FALSE,
                        keep.par = FALSE,
                        verbose = TRUE){

    if(length(dxdy) == 1) dxdy <- rep(dxdy, 2)

    xrange <- sort(xrange)
    yrange <- sort(yrange)

    xgr <- seq(xrange[1], xrange[2], by = dxdy[1])
    xcen <- xgr[-1] - 0.5 * dxdy[1]
    ygr <- seq(yrange[1], yrange[2], by = dxdy[2])
    ycen <- ygr[-1] - 0.5 * dxdy[2]

    xygrid <- expand.grid(xcen, ycen)
    igrid <- expand.grid(idx = 1:length(xcen), idy = 1:length(ycen))
    celltable <- matrix(rep(NA, ((length(xgr)-1)*(length(ygr)-1))),
                        ncol = (length(xgr)-1))

    if(as.integer(select) != 0){
        opts <- options()
        on.exit(options(opts))
        options(locatorBell = FALSE)

        if(!keep.par){
            opar <- par(no.readonly = TRUE)
            on.exit(par(opar))
            par(mfrow = c(1,1))
        }

        if(verbose) writeLines("Please refer to the graphical device and select the cells to include/exclude in the grid by clicking onto the map. Right-click to exit the selection function.")

        if(as.integer(select) == 1){
            idx <- NULL
        }else if(as.integer(select) == 2){
            idx <- 1:(length(xcen) * length(ycen))
        }

        repeat{
            plot(xygrid[,1], xygrid[,2],
                 type = "n",
                 main = "Click on points to select/de-select,\nright-click to exit",
                 xlim = xrange,
                 ylim = yrange,
                 xaxs = "i", yaxs = "i",
                 xlab = "x", ylab = "y")
            c0 <- t(matrix(rep(NA, length(xcen) * length(ycen)),
                           ncol = length(xcen)))
            c0[idx] <- 1
            c0 <- t(c0)
            image(xgr, ygr,
                  t(c0),
                  col = adjustcolor("dodgerblue2",0.2),
                 xlim = xrange,
                 ylim = yrange,
                 add = TRUE)
            if(plot.land){
                maps::map("world",
                          xlim = xrange,
                          ylim = yrange,
                          fill = TRUE, plot = TRUE, add = TRUE,
                          col = grey(0.95),
                          border = grey(0.5))
            }
            abline(v = xgr)
            abline(h = ygr)
            points(xygrid[,1], xygrid[,2], col = "dodgerblue2")
            points(xygrid[idx,1], xygrid[idx,2], pch = 16,
                   col = "dodgerblue2")
            idx.sel <- identify(xygrid[,1], xygrid[,2], n = 1)
            idx <- if(any(idx == idx.sel)) idx[idx != idx.sel] else c(idx, idx.sel)

            if (length(idx.sel) == 0) {
                cat("No more points selected, exiting...\n")
                break
            }
        }

        idx <- sort(unique(idx))
        xygrid <- xygrid[idx,]
        igrid <- igrid[idx,]

    }

    celltable[cbind(igrid$idx,igrid$idy)] <- 1:nrow(igrid)

    grid <- list(xygrid = xygrid,
                 igrid = igrid,
                 celltable = celltable)
    attributes(grid) <- list(names = attributes(grid)$names,
                             xrange = xrange,
                             yrange = yrange,
                             dxdy = dxdy,
                             xgr = xgr,
                             ygr = ygr,
                             xcen = xcen,
                             ycen = ycen,
                             nx = length(xcen),
                             ny = length(ycen))

    return(grid)
}

##' Get dimensions from tagging data
##'
##' @param ctags Conventional tags
##' @param atags Archival tags
##'
##' @export
get.dim <- function(ctags = NULL, atags = NULL){

    trange.c <- trange.a <- xrange.c <- xrange.a <- yrange.c <- yrange.a <- NULL
    if(!is.null(ctags)){
        var <- c("t0","t1")
        trange.c <- range(ctags[,var], na.rm = TRUE)
        var <- c("x0","x1")
        xrange.c <- range(ctags[,var], na.rm = TRUE)
        var <- c("y0","y1")
        yrange.c <- range(ctags[,var], na.rm = TRUE)
    }
    if(!is.null(atags)){
        atags <- do.call(rbind, atags)
        var <- c("t")
        trange.a <- range(atags[,var], na.rm = TRUE)
        var <- c("x")
        xrange.a <- range(atags[,var], na.rm = TRUE)
        var <- c("y")
        yrange.a <- range(atags[,var], na.rm = TRUE)
    }

    res <- list(trange = range(trange.c, trange.a, na.rm = TRUE),
                xrange = range(xrange.c, xrange.a, na.rm = TRUE),
                yrange = range(yrange.c, yrange.a, na.rm = TRUE))
    return(res)
}


##' Get neighbours for each cell of a grid
##'
##' @param celltable celltable
##'
##' @export
get.neighbours <- function(celltable, diagonal = FALSE){
    if(diagonal){
        AT <- array(c(celltable,
                      ## top
                      cbind(celltable[,-1], NA),
                      ## down
                      cbind(NA, celltable[,-ncol(celltable)]),
                      ## left
                      rbind(celltable[-1,], NA),
                      ## right
                      rbind(NA, celltable[-nrow(celltable),]),
                      ## top-right
                      cbind(rbind(celltable[-1, -1], NA), NA),
                      ## top-left
                      cbind(NA, rbind(celltable[-1, -ncol(celltable)], NA)),
                      ## bottom-right
                      cbind(rbind(NA, celltable[-nrow(celltable), -1]), NA),
                      ## bottom-left
                      cbind(NA,rbind(NA, celltable[-nrow(celltable), -ncol(celltable)]))),
                    dim = c(nrow(celltable), ncol(celltable), 9))
        labs <- c("c", "right", "left", "top", "down",
                  "top-right", "top-left",
                  "bottom-right", "bottom-left")
    }else{
        AT <- array(c(celltable,
                      ## top
                      cbind(celltable[,-1], NA),
                      ## down
                      cbind(NA, celltable[,-ncol(celltable)]),
                      ## left
                      rbind(NA, celltable[-nrow(celltable),]),
                      ## right
                      rbind(celltable[-1,], NA)
                      ),
                    dim = c(nrow(celltable), ncol(celltable), 5))
        labs <- c("c", "top", "down", "left", "right")
    }
    nextTo <- apply(AT, 3, function(x) x)
    nextTo <- nextTo[!is.na(nextTo[,1]),]
    nextTo <- nextTo[order(nextTo[,1]),]
    colnames(nextTo) <- labs
    return(nextTo)
}


##' Get precision matrix
##'
##' @param grid Grid
##' @param delta Parameter
##'
##' @export
get.precision.matrix <- function(grid, h, nu, rho, delta,
                                 matern = TRUE, diagonal = FALSE){
    n <- nrow(grid$xygrid)

    if(matern){

        ## Compute the pairwise distance matrix between all points in the grid
        dist_matrix <- as.matrix(dist(grid$xygrid))

        ## Define the Matérn covariance function for distances
        matern_covariance <- function(h, nu, rho) {
            if(h == 0) return(1)  ## Variance at distance 0
            scale_factor <- (2^(1 - nu)) / gamma(nu)
            distance_factor <- (sqrt(2 * nu) * h / rho)^nu
            bessel_part <- besselK(sqrt(2 * nu) * h / rho, nu)
            return(scale_factor * distance_factor * bessel_part)
        }

        ## Construct the covariance matrix using the Matérn function
        C <- matrix(0, n, n)
        for (i in 1:n) {
            for (j in 1:n) {
                C[i, j] <- matern_covariance(dist_matrix[i, j], nu, rho)
            }
        }

        ## Invert the covariance matrix to get the precision matrix Q
        Q0 <- solve(C)

    }else{

        nextTo <- get.neighbours(grid$celltable, diagonal)
        Q0 <- Matrix::sparseMatrix(1:n, 1:n, x=0, dims = c(n, n))
        diag(Q0) <- rowSums(!is.na(nextTo[,-1]), na.rm=TRUE)
        for(i in 1:nrow(nextTo)){
            Q0[i,nextTo[,-1][i,!is.na(nextTo[,-1][i,])]] <- -1
        }

    }

    ## Add delta * I to ensure positive definiteness
    I <- Matrix::sparseMatrix(1:n, 1:n, x=0, dims = c(n, n))
    diag(I) <- 1
    Q <- Q0 + delta * I

    return(Q)
}


inv.logit <- function(x){
  1 / (1 + exp(-x))
}


get.momo.cols <- function(n = 1, alpha = 1, type = NULL){

    if(is.null(type) || is.na(type)){
        adjustcolor(c("dodgerblue3","goldenrod2",
                      "darkgreen","purple",
                      "hotpink")[1:n], alpha)
    }else if(type == "true"){
        adjustcolor("darkorange", alpha)
    }else if(type == "pos"){
        adjustcolor("dodgerblue3", alpha)
    }else if(type == "neg"){
        adjustcolor("firebrick3", alpha)
    }else if(type == "notsig"){
        adjustcolor("chartreuse4", alpha)
    }else if(type == "sig"){
        adjustcolor("firebrick3", alpha)
    }else{
        adjustcolor(c("dodgerblue3","goldenrod2",
                      "darkgreen","purple",
                      "hotpink")[1:n], alpha)
    }
}


is.empty <- function(x){
    length(x) == 0
}

check.that.list <- function(x){
    if(!is.null(x) && !inherits(x, "list")){
        x <- list(x)
    }
    return(x)
}

get.map <- function(x){
    x[!is.na(x)] <- 1:sum(!is.na(x))
    return(factor(x))
}

pfx <- function(x, eps=0.01){
  eps * log(exp(x/eps) + 1)
}

smooth.identity <- function(x, from=0, to=1){
  z <- (x - from) / (to - from)
  y <- 1 - pfx(1 - pfx(z))
  y * (to - from) + from
}


matexpo <- function(mat, log2steps){
    if((log2steps <= 0) || (log2steps > 100)){
        return(Matrix::expm(mat))
    }else{
        ng <- nrow(mat)
        Imat <- RTMB::matrix(0, ng, ng)
        Imat[cbind(1:ng,1:ng)] <- 1
        mat <- Imat + mat / 2^log2steps
        for(i in 1:log2steps){
            mat <- mat %*% mat
        }
        return(mat)
    }
}


poly.fun <- function(xp, yp, deriv = FALSE, adv = FALSE){

    if(adv){

        f <- local({
            val <- yp[1]
            function(x){
                RTMB::sapply(x, function(x) x * val)
            }
        })

        df <- local({
            val <- yp[1]
            function(x){
                RTMB::sapply(x, function(x) x * val)
            }
        })

    }else{

        n <- length(xp)

        if(n == 1){

            f <- local({
                val <- yp[1]
                function(x){
                    RTMB::sapply(x, function(x) val)
                }
            })

            df <- local({
                val <- yp[1]
                function(x){
                    RTMB::sapply(x, function(x) val)
                }
            })


        }else{

            A <- outer(xp, 0:(n-1), "^")
            alpha <- RTMB::solve(A, yp)

            f <- local({
                vals <- alpha
                function(x){
                    RTMB::sapply(x,
                                 function(x) vals[1] +
                                             sum(vals[-1] * x^(1:(n-1))))
                }
            })

            df <- local({
                vals <- alpha
                function(x){
                    RTMB::sapply(x,
                                 function(x) vals[2] +
                                             sum(vals[-c(1:2)] *
                                                 (2:(n-1)) * x^(1:(n-2))))
                }
            })

        }
    }

    if(!deriv) return(f) else return(df)
}



add.class <- function(x, class){
    if(!inherits(x, class)){
        class(x) <- c(class(x), class)
    }
    return(x)
}


get.env.funcs <- function(dat, conf, par){

    nenv <- length(dat$env)

    ## Setup environmental functions --------------------------
    env.func.tax <- env.func.dif <-
        env.func.adv.x <- env.func.adv.y <-
            env.dfunc.tax <- env.dfunc.dif <-
                env.dfunc.adv <- vector("list", nenv)
    for(i in 1:nenv){
        env.func.tax[[i]] <- poly.fun(conf$knots.tax[,i], par$alpha[,i])
        env.dfunc.tax[[i]] <- poly.fun(conf$knots.tax[,i], par$alpha[,i],
                                       deriv = TRUE)
        env.func.dif[[i]] <- poly.fun(conf$knots.dif[,i], par$beta[,i])
        env.func.adv.x[[i]] <- poly.fun(NULL, par$gamma[1,i], adv = TRUE)
        env.func.adv.y[[i]] <- poly.fun(NULL, par$gamma[2,i], adv = TRUE)
    }

    ## Setup habitat objects ---------------------------------
    habitat.tax <- habi.full(dat$env, dat$xranges, dat$yranges,
                             conf$ienv$tax, dat$time.cont,
                             env.func.tax, env.dfunc.tax)
    habitat.dif <- habi.full(dat$env, dat$xranges, dat$yranges,
                             conf$ienv$dif, dat$time.cont,
                             env.func.dif, env.dfunc.dif)
    habitat.adv.x <- habi.full(dat$env, dat$xranges, dat$yranges,
                               conf$ienv$adv.x, dat$time.cont,
                               env.func.adv.x, env.dfunc.adv)
    habitat.adv.y <- habi.full(dat$env, dat$xranges, dat$yranges,
                               conf$ienv$adv.y, dat$time.cont,
                               env.func.adv.y, env.dfunc.adv)

    if(conf$use.effort){
        effort <- habi.light(dat$effort, dat$xranges.eff, dat$yranges.eff,
                             dat$ieff, dat$time.cont)
    }

    diffusion.fun <- function(xy, t){
        habitat.dif$val(xy, t)
    }

    taxis.fun <- function(xy, t){
        habitat.tax$grad(xy, t)
    }

    ## TODO: advection + flags to turn stuff on and off


    res <- list(dif = diffusion.fun,
                tax = taxis.fun)

    if(conf$use.boundaries){
        bound <- habi.light(conf$boundaries,
                            dat$xranges,## + c(-1,1) * dat$dxdy[1],
                            dat$yranges,## + c(-1,1) * dat$dxdy[2],
                            conf$ibound, dat$time.cont)
        bound.fun <- function(xy, t){
            bound$val(xy, t)
        }
        res$bound <- bound.fun

    }

    return(res)
}


get.sim.par <- function(par = NULL){

    par.out <- list(alpha = matrix(c(0, 0.1, 0.001), 3, 1),
                    beta = matrix(log(0.02), 1, 1),
                    gamma = matrix(0, 2, 1),
                    logLambda = matrix(log(0.5), 1, 1),
                    logM = log(0.3),
                    logSdObsATS = log(0.05))

    if(!is.null(par)){
        for(i in 1:length(par)){
            par.out[names(par)[i]] <- par[names(par)[i]]
        }
    }

    return(par.out)
}


get.sim.funcs <- function(funcs = NULL, dat, conf, env, par){

    if(!is.null(env)){
        funcs <- c(funcs, get.env.funcs(dat, conf, par))
    }

    diffusion.fun <- function(xy, t){ par$beta[1,1] }
    taxis.fun <- function(xy, t){c(-0.5 * xy[1] + 0.25,
                                   -0.5 * xy[2] + 0.25)}
    ## taxis.fun <- function(xy, t){
    ##     atr1 <- c(0.15, 0.15)
    ##     atr2 <- c(0.85, 0.85)
    ##     vec1 <- (atr1 - xy) ## / sqrt(sum((atr1 - xy)^2))
    ##     vec2 <- (atr2 - xy) ## / sqrt(sum((atr2 - xy)^2))
    ##     p <- if(xy[1] < (1 - xy[2])) 1 else 0
    ##     c(p * vec1[1] + (1-p) * vec2[1],
    ##       p * vec1[2] + (1-p) * vec2[2])
    ## }
## taxis.fun <- function(xy, t){c(ifelse(xy[1] > (1 - xy[2]),
##                                       -0.5 * xy[1] + 0.4,
##                                       -0.5 * xy[1] + 0.15),
##                                ifelse(xy[1] > (1 - xy[2]),
##                                       -0.5 * xy[2] + 0.4,
##                                       -0.5 * xy[2] + 0.15))}
    adv.fun <- function(xy, t){ c(0,0) }

    bound.fun <- function(xy, t){ 1 }

    fish.mort.fun <- function(xy, t){ exp(par$logLambda) }
    ## 0.5 * dnorm(xy[1],.2,.3)*dnorm(xy[2],.2,.3)/(dnorm(0,sd=.3)^2)
    nat.mort.fun <- function(xy, t){ exp(par$logM) }


    funcs.out <- list(
        tax = taxis.fun,
        dif = diffusion.fun,
        adv = adv.fun,
        bound = bound.fun,
        fish.mort = fish.mort.fun,
        nat.mort = nat.mort.fun)

    if(!is.null(funcs)){
        for(i in 1:length(funcs)){
            funcs.out[names(funcs)[i]] <- funcs[names(funcs)[i]]
        }
    }

    return(funcs.out)
}



##' Get neighbours for each cell of a grid
##'
##' @param fit fit
##'
##' @export
get.diag <- function(fit){

    if(!inherits(fit, "momo.fit")) stop("fit must be a fitted momo object ('momo.fit')!")

    get.diag.single <- function(fit, var,
                                is.atags = FALSE,
                                select = NULL, do.box = TRUE){

        if(is.atags){
            browser()
            tmp <- data.frame(unlist(fit$dat$atags),
                              resid = fit$rep[[var]][,select])
        }else{
            tmp <- data.frame(fit$dat$ctags,
                              resid = fit$rep[[var]][,select])
            tmp <- tmp[fit$conf$excl.ctags == 0,]
        }
        tmp <- tmp[!is.na(tmp$resid),]
        tmp <- tmp[order(tmp$t1),]

        res <- list(mean = mean(tmp$resid, na.rm = TRUE),
                    ttest = t.test(tmp$resid))
        moransI <- try(spdep::moran.test(
                                  tmp$resid,
                                  spdep::nb2listw(
                                             spdep::dnearneigh(
                                                        tmp[,c("x1","y1")],
                                                        d1 = 0,
                                                        d2 =  mean(fit$dat$dxdy)))))
        if(inherits(moransI, "try-error")){
            moransI <- list(p.value = NA)
        }
        res$moransI <- moransI

        if(do.box){
            res$box <- Box.test(tmp$resid, lag = 4, type = "Ljung-Box", fitdf = 1)
        }
        res$shapiro <- shapiro.test(tmp$resid)
        return(res)
    }

    ret <- list()
    if(fit$conf$use.ctags){
        ret$ctags.x <- get.diag.single(fit, "resid.ctags", select = 1)
        ret$ctags.y <- get.diag.single(fit, "resid.ctags", select = 2)
        if(fit$conf$use.effort){
            ret$ctags.t <- get.diag.single(fit, "resid.ctags", select = 4)
        }
    }
    if(fit$conf$use.atags){
        ret$atags.x <- get.diag.single(fit, "resid.atags", select = 1)
        ret$atags.y <- get.diag.single(fit, "resid.atags", select = 2)
        if(fit$conf$use.effort){
            ret$atags.t <- get.diag.single(fit, "resid.atags", select = 4)
        }
    }

    if(length(ret) > 0){
        fit$diag <- ret
    }

    return(fit)
}



get.adv <- function(dat, par, conf, funcs = NULL){
    par <- get.sim.par(par)
    env <- check.that.list(dat$env)
    dat$env <- env
    dat$env.pred <- NULL
    conf <- def.conf(dat)
    funcs <- get.sim.funcs(funcs, dat, conf, env, par)
    hTdx.true <- sapply(dat$time.cont.pred,
                        function(t) apply(dat$xygrid.pred, 1, function(x) funcs$tax(x,t)[1]))
    hTdy.true <- sapply(dat$time.cont.pred,
                        function(t) apply(dat$xygrid.pred, 1, function(x) funcs$tax(x,t)[2]))
    uv.true <- matrix(NA, nrow(dat$xygrid), 2)
    uv.true[,1] <- rowMeans(hTdx.true)
    uv.true[,2] <- rowMeans(hTdy.true)

    return(uv.true)
}


##' Get peclet number
##'
##' @param grid grid
##'
##' @export
get.peclet <- function(grid, env, par){

    dat <- setup.momo.data(grid = grid,
                           env = env)
    conf <- def.conf(dat)

    uv <- get.adv(dat, par, conf)

    ## Péclet number
    peclet <- data.frame(u = uv[,1] * dat$dxdy[1]/mean(exp(par$beta)),
                         v = uv[,2] * dat$dxdy[2]/mean(exp(par$beta)))
    ## TODO: make similar get.dif function!

    ## For numerical stability: PE <= 2
    ## or: u/D <= 2/\{Delta}x

    return(peclet)
}


##' Get Courant–Friedrichs–Lewy condition
##'
##' @param grid grid
##'
##' @export
get.cfl <- function(grid, env, par){

    dat <- setup.momo.data(grid = grid,
                           env = env)
    conf <- def.conf(dat)

    uv <- get.adv(dat, par, conf)

    u <- momo:::get.adv(dat, true1, conf)[,1]
    cfl <- min(dat$dxdy[1] / abs(u), prod(dat$dxdy) / (2 * mean(exp(true1$beta))))
    ## improve beta component! get.dif()


    ## Courant–Friedrichs–Lewy (CFL) condition: delta t <= min(Delta x/ u, Delta x^2 / 2D)

    return(cfl)
}


## TODO: improve by allowing to use for diffusion, advection,
## TODO: option to select rel and/or rec!
get.env <- function(dat, conf){

    if(is.null(dat$ctags) && is.null(dat$atags)){
        return(NULL)
    }

    nenv <- length(dat$env)

    xind <- yind <- tind <- NULL
    if(!is.null(dat$ctags)){
        xind <- c(xind, as.integer(cut(c(dat$ctags$x0,dat$ctags$x1), dat$xgr,
                               include.lowest = TRUE)))
        yind <- c(yind, as.integer(cut(c(dat$ctags$y0,dat$ctags$y1), dat$ygr,
                               include.lowest = TRUE)))
        tind <- c(tind, as.integer(cut(c(dat$ctags$t0,dat$ctags$t1), dat$time.cont, include.lowest = TRUE)))
    }

    if(!is.null(dat$atags)){
        xind <- c(xind, as.integer(cut(unlist(sapply(dat$atags, function(x) x$x)),
                                       dat$xgr,
                               include.lowest = TRUE)))
        yind <- c(yind, as.integer(cut(unlist(sapply(dat$atags, function(x) x$y)),
                                       dat$ygr,
                               include.lowest = TRUE)))
        tind <- c(tind, as.integer(cut(unlist(sapply(dat$atags, function(x) x$t)),
                                       dat$time.cont, include.lowest = TRUE)))
    }

    env.obs <- vector("list", nenv)
    for(i in 1:nenv){
        env.obs[[i]] <- dat$env[[i]][cbind(xind,yind,conf$ienv$tax[i,tind])]
    }

    return(env.obs)
}
