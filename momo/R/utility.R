##' Create a grid
##'
##' @description `create.grid` allows to create a new grid or modify an existing
##'     grid that is used by [fit.momo] to predict movement rates and required
##'     for the matrix exponential approach.
##'
##' @param xrange range of the x-dimension of spatial domain. Default: `c(0,1)`.
##' @param yrange range of the y-dimension of spatial domain. Default: `c(0,1)`.
##' @param dxdy resolution of grid in x and y direction. Default: `c(0.1,0.1)`.
##' @param select logical; if `TRUE`, allows to select cells in spatial grid.
##'     Default: `FALSE`. If a vector of length > 1 is provided, the vector is
##'     interpreted as an index of the selected grid cells of the full grid.
##' @param plot.land plot.land logical; If `TRUE`, plot land masses using the
##'     function [maps::map]. Default: `FALSE`.
##' @param keep.gpar logical; If `TRUE`, do not overwrite the graphical
##'     parameters. Default: `FALSE`.
##' @param verbose if `TRUE`, print information to console. Default: `TRUE`.
##' @param fit optional; allows to extract and update grid from fitted object of
##'     class `momo.fit` returned by function [fit.momo]. Default: `NULL`. NULL.
##' @param grid optional; allows to update a grid of class `momo.grid` as
##'     returned by the function [create.grid]. Default: `NULL`.
##'
##' @return A list of class `momo.grid` with a information about the spatial
##'     domain and grid.
##'
##' @examples
##' grid <- create.grid()
##'
##' @export
create.grid <- function(xrange = c(0,1),
                        yrange = c(0,1),
                        dxdy = c(0.1,0.1),
                        select = FALSE,
                        plot.land = FALSE,
                        keep.gpar = FALSE,
                        verbose = TRUE,
                        fit = NULL,
                        grid = NULL){

    if(!is.null(fit)){

        grid <- list(xygrid = fit$dat$xygrid,
                     igrid = fit$dat$igrid,
                     celltable = fit$dat$celltable)
        attributes(grid) <- list(names = attributes(grid)$names,
                                 xrange = range(grid$xygrid[,1]) +
                                     c(-1,1) * fit$dat$dxdy[1]/2,
                                 yrange = range(grid$xygrid[,2]) +
                                     c(-1,1) * fit$dat$dxdy[2]/2,
                                 dxdy = fit$dat$dxdy,
                                 xgr = fit$dat$xgr,
                                 ygr = fit$dat$ygr,
                                 xcen = unique(grid$xygrid[,1]),
                                 ycen = unique(grid$xygrid[,2]),
                                 nx = fit$dat$nx,
                                 ny = fit$dat$ny)

        return(grid)
    }

    if(!is.null(grid)){
        xrange <- attr(grid, "xrange")
        yrange <- attr(grid, "yrange")
        idx <- which(is.na(grid$celltable), arr.ind = TRUE)
    }else{
        xrange <- sort(xrange)
        yrange <- sort(yrange)
    }

    if(length(dxdy) == 1) dxdy <- rep(dxdy, 2)


    xgr <- seq(xrange[1], xrange[2], by = dxdy[1])
    xcen <- xgr[-1] - 0.5 * dxdy[1]
    ygr <- seq(yrange[1], yrange[2], by = dxdy[2])
    ycen <- ygr[-1] - 0.5 * dxdy[2]

    xygrid <- expand.grid(x = xcen, y = ycen)
    igrid <- expand.grid(idx = 1:length(xcen), idy = 1:length(ycen))
    celltable <- matrix(rep(NA, ((length(xgr)-1)*(length(ygr)-1))),
                        nrow = (length(xgr)-1))

    if(length(select) > 1){

        idx <- sort(unique(select))
        xygrid <- xygrid[idx,]
        igrid <- igrid[idx,]

    }else if(as.integer(select) != 0){
        opts <- options()
        on.exit(options(opts))
        options(locatorBell = FALSE)

        if(!keep.gpar){
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
                plot.land(xrange, yrange,
                          shift = ifelse(max(xrange) > 180, TRUE, FALSE))
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

    if(!is.null(grid)){
        idx1 <- grid$celltable[cbind(cut(xygrid$x, attr(grid, "xgr")),
                                    cut(xygrid$y, attr(grid, "ygr")))]
        idx2 <- grid$celltable[cbind(cut(xygrid$x, attr(grid, "xgr"),
                                         right = FALSE),
                                     cut(xygrid$y, attr(grid, "ygr"),
                                         right = FALSE))]
        idx3 <- grid$celltable[cbind(cut(xygrid$x, attr(grid, "xgr")),
                                     cut(xygrid$y, attr(grid, "ygr"),
                                         right = FALSE))]
        idx4 <- grid$celltable[cbind(cut(xygrid$x, attr(grid, "xgr"),
                                         right = FALSE),
                                     cut(xygrid$y, attr(grid, "ygr")))]
        idx <- 1:nrow(xygrid)
        idx[apply(cbind(idx1,idx2,idx3,idx4),1,function(x) sum(is.na(x)) > 2)] <- NA
        idx <- which(!is.na(idx))
        if(length(idx) > 0){
            xygrid <- xygrid[idx,]
            igrid <- igrid[idx,]
        }
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

    ## Return
    grid <- add.class(grid, "momo.grid")
    return(grid)
}



##' Get dimensions from tagging data
##'
##' @description `get.dim` extracts the time and spatial dimensions from tagging
##'     data.
##'
##' @param ctags a data frame with mark-recapture tags of class `momo.ctags` as
##'     returned by the function [prep.ctags]. Default: `NULL`.
##' @param atags a list with archival tags of class `momo.atags`. as returned by
##'     the function [prep.atags]. Default: `NULL`.
##' @param stags a list with mark-resight tags of class `momo.stags`. as
##'     returned by the function [prep.stags]. Default: `NULL`.
##'
##' @return A list with the temporal and spatial (in x and y direction)
##'     dimensions.
##'
##' @export
get.dim <- function(ctags = NULL, atags = NULL, stags = NULL){

    trange.c <- trange.a <- trange.s <-
        xrange.c <- xrange.a <- xrange.s <-
            yrange.c <- yrange.a <- yrange.s <-NULL
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
    if(!is.null(stags)){
        stags <- do.call(rbind, stags)
        var <- c("t")
        trange.s <- range(stags[,var], na.rm = TRUE)
        var <- c("x")
        xrange.s <- range(stags[,var], na.rm = TRUE)
        var <- c("y")
        yrange.s <- range(stags[,var], na.rm = TRUE)
    }

    res <- list(trange = range(trange.c, trange.a, trange.s, na.rm = TRUE),
                xrange = range(xrange.c, xrange.a, xrange.s, na.rm = TRUE),
                yrange = range(yrange.c, yrange.a, yrange.s, na.rm = TRUE))
    return(res)
}


##' Get neighbours for each cell of a grid
##'
##' @description `get.neighbours` allows to get a
##'
##' @param celltable data frame with numbered cells as returned by the function
##'     [create.grid].
##' @param diagonal logical; if `TRUE`, diagonal neighbours are included.
##'     Default: `FALSE`.
##'
##' @return A matix with each cell and corresponding neighbouring cells.
##'
##' @examples
##' grid <- create.grid()
##'
##' neighbours <- get.neighbours(grid$celltable)
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
##' @description `get.precision.matrix` allows to create a precision matrix
##'     based on a grid.
##'
##' @param grid a grid object of class `momo.grid` as returned by the function
##'     [create.grid].
##' @param h parameter of the matern covariance structure.
##' @param nu parameter of the matern covariance structure.
##' @param rho parameter of the matern covariance structure.
##' @param delta parameter of the matern covariance structure.
##' @param matern logical; if `TRUE` (default), matern covariance structure is
##'     used for simulation of environmental fields.
##' @param diagonal logical; if `TRUE`, diagonal neighbours are considered in
##'     neighbouring structure. Default: `FALSE`.
##'
##' @return A precision matrix.
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


##' Momo colors
##'
##' @description `momo.cols` returns a vector with colors for different
##'     scenarios.
##'
##' @param n number of colors. Default: `1`.
##' @param alpha transparency value. Default: `1`.
##' @param type optional; allows to extract colors for specific occasions.
##'    Default: `NULL`.
##'
##' @return A vector with colors.
##'
##' @export
momo.cols <- function(n = 1, alpha = 1, type = NULL){

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
    is.null(x) || length(x) == 0
}

check.that.list <- function(x){
    x0 <- x
    if(!is.null(x) && !inherits(x, "list")){
        x <- list(x)
        attributes(x[[1]]) <- attributes(x0)
    }
    return(x)
}

get.map <- function(x){
    x[!is.na(x)] <- 1:sum(!is.na(x))
    return(factor(x))
}

pfx <- function(x, eps=0.01){
  return(eps * log(exp(x/eps) + 1))
}

smooth.identity <- function(x, from=0, to=1){
  z <- (x - from) / (to - from)
  y <- 1 - pfx(1 - pfx(z))
  return(y * (to - from) + from)
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


poly.fun.old <- function(xp, yp, deriv = FALSE, adv = FALSE){

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


poly.fun <- function(xp, yp, deriv = FALSE, adv = FALSE){

    if(!adv && length(xp) > 1 && all(diff(xp) == 0)) return(NULL)

    if(adv){
        ## Simple linear case
        val <- yp[1]
        f <- function(x) x * val
        df <- function(x) rep(val, length(x))
    }else{
        if(length(xp) == 1){
            ## For one knot return parameter (e.g. constant diffusion)
            val <- yp[1]
            f <- function(x) val
            df <- function(x) rep(val, length(x))
        }else{
            ## Solve for polynomial coefficients
            n <- length(xp)
            A <- outer(xp, 0:(n-1), "^")
            alpha <- RTMB::solve(A, yp)

            f <- function(x){
                ## Evaluate polynomial: sum(alpha[j+1] * x^j)
                v <- outer(x, 0:(n-1), "^")
                as.vector(v %*% alpha)
                ## as.vector(alpha[1] + sum(alpha[-1] * x^(1:(n-1))))
            }

            df <- function(x){
                ## Evaluate derivative: sum(j * alpha[j+1] * x^(j-1))
                if (n == 2) {
                  rep(alpha[2], length(x))  # Linear case
                }else{
                  j <- 1:(n-1)
                  v <- outer(x, j - 1, "^")
                    as.vector(v %*% (j * alpha[-1]))
                }
                ## as.vector(alpha[2] +
                ##           sum(alpha[-c(1:2)] *
                ##               (2:(n-1)) * x^(1:(n-2))))
            }
        }
  }

  if (!deriv) return(f) else return(df)
}




add.class <- function(x, class){
    if(!inherits(x, class)){
        class(x) <- c(class(x), class)
    }
    return(x)
}

check.class <- function(x, class){
    if(!inherits(x, class)){
        stop(paste0("The object ",deparse(substitute(x)),
                       " does not inherit class ", class,
                       ". Please check your code."))
    }
    return(invisible(NULL))
}


get.env.funcs <- function(dat, conf, par){

    nenv <- length(dat$env)

    ## Setup environmental functions --------------------------
    ## env.func.tax <- env.func.dif <-
    ##     env.func.adv.x <- env.func.adv.y <-
    ##         env.dfunc.tax <- env.dfunc.dif <-
    ##             env.dfunc.adv <- vector("list", nenv)
    ## for(i in 1:nenv){
    ##     env.func.tax[[i]] <- poly.fun(dat$knots.tax[,i], par$alpha[,i])
    ##     env.dfunc.tax[[i]] <- poly.fun(dat$knots.tax[,i], par$alpha[,i],
    ##                                    deriv = TRUE)
    ##     env.func.dif[[i]] <- poly.fun(dat$knots.dif[,i], par$beta[,i])
    ##     env.func.adv.x[[i]] <- poly.fun(NULL, par$gamma[1,i], adv = TRUE)
    ##     env.func.adv.y[[i]] <- poly.fun(NULL, par$gamma[2,i], adv = TRUE)
    ## }

    ## TRY:
    env.func.tax <- env.func.dif <-
        env.func.adv.x <- env.func.adv.y <-
            env.dfunc.tax <- env.dfunc.dif <-
                env.dfunc.adv <- vector("list", nenv)
    for(i in 1:nenv){
        env.func.tax[[i]] <- env.dfunc.tax[[i]] <- vector("list", dim(par$alpha)[3])
        for(j in 1:dim(par$alpha)[3]){
            env.func.tax[[i]][[j]] <- poly.fun(dat$knots.tax[,i], par$alpha[,i,j])
            env.dfunc.tax[[i]][[j]] <- poly.fun(dat$knots.tax[,i], par$alpha[,i,j], deriv = TRUE)
        }
        env.func.dif[[i]] <- env.dfunc.dif[[i]] <- vector("list", dim(par$beta)[3])
        for(j in 1:dim(par$beta)[3]){
            env.func.dif[[i]][[j]] <- poly.fun(dat$knots.dif[,i], par$beta[,i,j])
        }
        env.func.adv.x[[i]] <- env.func.adv.y[[i]] <- env.dfunc.adv[[i]] <- vector("list", dim(par$gamma)[3])
        for(j in 1:dim(par$gamma)[3]){
            env.func.adv.x[[i]][[j]] <- poly.fun(NULL, par$gamma[1,i,j], adv = TRUE)
            env.func.adv.y[[i]][[j]] <- poly.fun(NULL, par$gamma[2,i,j], adv = TRUE)
        }
    }


    ## Setup habitat objects ---------------------------------

    ## TODO: build in check that envs are equally spaced grid (maybe put in check function)
    dxdy <- cbind(do.call(rbind, lapply(dat$env, function(x)
        diff(as.numeric(attributes(x)$dimnames[[1]]))[1])),
        do.call(rbind, lapply(dat$env, function(x)
        diff(as.numeric(attributes(x)$dimnames[[2]]))[1])))

    liv.all <- get.liv(dat$env)
                       ## dat$xranges + c(1,-1) * dxdy[,1]/2,
                       ## dat$yranges + c(1,-1) * dxdy[,2]/2)
    liv.tax <- liv.all
    liv.dif <- liv.all
    liv.adv.x <- liv.all
    liv.adv.y <- liv.all

    if(conf$use.effort){
        liv.effort <- get.liv(dat$effort) ## , dat$xranges.eff, dat$yranges.eff)
    }else{
        liv.effort <- NULL
    }

    if(conf$use.boundaries){
        browser()
        liv.bound <- get.liv(list(array(dat$boundaries,
                                   dim = c(nrow(dat$boundaries),
                                           ncol(dat$boundaries),1))))
                             ## matrix(dat$boundary.xrange + c(1,-1) *
                             ##        dat$boundary.dxdy[1]/2,1,2),
                             ## matrix(dat$boundary.yrange+ c(1,-1) *
                             ##        dat$boundary.dxdy[2]/2,1,2))
        liv.bound
    }else{
        liv.bound <- NULL
    }

    habitat.tax <- habi.full(liv.tax, dat$xranges, dat$yranges,
                             conf$ienv$tax, dat$time.cont,
                             env.func.tax, env.dfunc.tax, conf$ienvS$tax)
    habitat.dif <- habi.full(liv.dif, dat$xranges, dat$yranges,
                             conf$ienv$dif, dat$time.cont,
                             env.func.dif, env.dfunc.dif, conf$ienvS$dif)
    habitat.adv.x <- habi.full(liv.adv.x, dat$xranges, dat$yranges,
                               conf$ienv$adv.x, dat$time.cont,
                               env.func.adv.x, env.dfunc.adv, confienvS$adv)
    habitat.adv.y <- habi.full(liv.adv.y, dat$xranges, dat$yranges,
                               conf$ienv$adv.y, dat$time.cont,
                               env.func.adv.y, env.dfunc.adv, conf$ienvS$adv)

    diffusion.fun <- function(xy, t){
        habitat.dif$val(xy, t)
    }

    taxis.fun <- function(xy, t){
        habitat.tax$grad(xy, t)
    }

    ## TODO: advection + flags to turn stuff on and off


    res <- list(dif = diffusion.fun,
                tax = taxis.fun)

    if(conf$use.effort){
        effort <- habi.light(liv.effort, dat$xranges.eff,
                             dat$yranges.eff,
                             dat$ieff, dat$time.cont)

        fish.mort.fun <- function(xy, t){
            exp(par$logLambda[1,1]) * effort$val(xy, t)
        }

        res$fish.mort <- fish.mort.fun
    }


    if(conf$use.boundaries){
        bound <- habi.light(liv.bound,
                            dat$boundary.xrange + c(1,-1) *
                                    dat$boundary.dxdy[1]/2,
                            dat$boundary.yrange + c(1,-1) *
                                    dat$boundary.dxdy[2]/2,
                            conf$ibound, dat$time.cont)
        bound.fun <- function(xy, t){
            bound$val(xy, t)
        }
        res$bound <- bound.fun

    }

    return(res)
}


get.sim.par <- function(par = NULL){

    par.out <- list(alpha = array(c(0, 0.1, 0.001), dim = c(3,1,1)),
                    beta = array(log(0.02), dim = c(1,1,1)),
                    gamma = array(0, dim = c(2,1,1)),
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

##' Simulate functions
##'
##' @description Simulate movement functions.
##'
##'
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



##' Get diagnostics for residuals
##'
##' @description `get.diag` applies diagnostics tests to the residuals.
##'
##' @param fit a list of class `momo.fit` as returned by the function [fit.momo].
##'
##' @return Fitted object with the list element `diag` with the results of the
##'     diagnostic tests.
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
##' @description `get.peclet` calculates the peclet number.
##'
##' @param grid a grid object of class `momo.grid` as returned by the function
##'     [create.grid].
##' @param env a list with environmental covariates of class `momo.env` as
##'     returned by the function [prep.env].
##' @param par parameter list with initial values as produced by the function
##'     [def.par].
##' @param knots.tax knots for the taxis component. Default: `NULL`.
##' @param knots.dif knots for the diffusion component. Default: `NULL`.
##'
##' @return The peclet number.
##'
##' @export
get.peclet <- function(grid, env, par, knots.tax = NULL, knots.dif = NULL){

    dat <- setup.momo.data(grid = grid,
                           env = env,
                           knots.tax = knots.tax,
                           knots.dif = knots.dif)
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
##' @description `get.cfl` allows to calculate the Courant-Friedrichs-Lewy
##'     condition.
##'
##' @param grid a grid object of class `momo.grid` as returned by the function
##'     [create.grid].
##' @param env a list with environmental covariates of class `momo.env` as
##'     returned by the function [prep.env].
##' @param par parameter list with initial values as produced by the function
##'     [def.par].
##'
##' @return The Courant-Friedrichs-Lewy condition.
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

    env.obs <- vector("list", nenv)
    for(i in 1:nenv){

        envi <- dat$env[[i]]
        xgr <- as.numeric(attributes(envi)$dimnames[[1]])
        ygr <- as.numeric(attributes(envi)$dimnames[[2]])

        xind <- yind <- tind <- NULL
        if(!is.null(dat$ctags)){
            xind <- c(xind, as.integer(cut(c(dat$ctags$x0,dat$ctags$x1), xgr,
                                           include.lowest = TRUE)))
            yind <- c(yind, as.integer(cut(c(dat$ctags$y0,dat$ctags$y1), ygr,
                                           include.lowest = TRUE)))
            tind <- c(tind, as.integer(cut(c(dat$ctags$t0,dat$ctags$t1),
                                           dat$time.cont, include.lowest = TRUE)))
        }

        if(!is.null(dat$atags)){
            xind <- c(xind, as.integer(cut(unlist(sapply(dat$atags, function(x) x$x)),
                                           xgr,
                                           include.lowest = TRUE)))
            yind <- c(yind, as.integer(cut(unlist(sapply(dat$atags, function(x) x$y)),
                                           ygr,
                                           include.lowest = TRUE)))
            tind <- c(tind, as.integer(cut(unlist(sapply(dat$atags, function(x) x$t)),
                                           dat$time.cont, include.lowest = TRUE)))
        }

        env.obs[[i]] <- dat$env[[i]][cbind(xind,yind,conf$ienv$tax[i,tind])]
    }

    return(env.obs)
}


##' Add a label to a plot
##'
##' @param lab label
##'
##' @export
add.lab <- function(lab){
    legend("topleft", legend = lab,
           bg = "white", x.intersp = -0.4,
           cex = 1.8, text.font = 2)
}


##' Get release events
##'
##' @description `get.release.events` allows to group mark-recapture tags to
##'     release events that are close in space and time.
##'
##' @param dat data frame with input data as produced by the function
##'     [check.momo.data].
##' @param grid optional; allows to provide an extra spatial grid (e.g. with a
##'     finer resolution) that is used to define the release events. If `NULL`
##'     (default), the grid in `dat` is used.
##' @param time.cont optional; allows to provide an extra time vector (e.g. with
##'     a finer resolution) that is used to define the release events. If `NULL`
##'     (default), the time vector in `dat` is used.
##' @param age.max optional; allows to specify a maximum age of the species
##'     which is used as the maximum time of each release event. If `NULL`
##'     (default), the maximum time difference of the release event and each
##'     associated tag that was recovered is used.
##'
##' @return A list with a data frame containing the release events and and index
##'     vector matching each mark-recapture tag to a release event.
##'
##' @examples
##' data(skjepo)
##'
##' dat <- skjepo$dat
##'
##' tmp <- get.release.events(dat,
##'                           grid = create.grid(dat$xrange, dat$yrange, c(1,1)),
##'                           time.cont = seq(dat$trange[1], dat$trange[2],
##'                                           1/(52*diff(dat$trange))))
##' dat$rel.events <- tmp$rel.events
##' dat$ctags$rel.event <- tmp$idx
##'
##' @export
get.release.events <- function(dat,
                               grid = NULL,
                               time.cont = NULL,
                               age.max = NULL){

    ctags <- dat$ctags

    ## Grid for aggregation
    if(is.null(grid)){
        celltable <- dat$celltable
        xygrid <- dat$xygrid
        xgr <- dat$xgr
        ygr <- dat$ygr
    }else{
        celltable <- grid$celltable
        xygrid <- grid$xygrid
        xgr <- attr(grid, "xgr")
        ygr <- attr(grid, "ygr")
    }

    ## TODO: make check if provided grid is not equal to model grid => problems for expm! but here we don't have expm=TRUE yet...

    ## Time for aggregation
    if(is.null(time.cont)){
        time.cont <- dat$time.cont
    }

    idx.space <- celltable[cbind(cut(ctags$x0, xgr),
                                 cut(ctags$y0, ygr))]

    idx.time <- as.integer(cut(ctags$t0, time.cont, include.lowest = TRUE))

    rel.all <- data.frame(xygrid[idx.space,],
                          time.cont[idx.time])
    colnames(rel.all) <- c("x0","y0","t0")

    ## Unique release events
    rel.events <- rel.all[!duplicated(rel.all),]
    rownames(rel.events) <- NULL

    ## Index to match tag to release event
    idx.all <- apply(rel.all, 1, paste, collapse = ":")
    idx.uni <- apply(rel.events, 1, paste, collapse = ":")
    idx <- match(idx.all, idx.uni)

    ## Maximum time for each release event
    if(is.null(age.max)){
        age.max <- ceiling(max(ctags$t1[!is.na(ctags$x0)] - ctags$t0[!is.na(ctags$x0)],
                               na.rm = TRUE))
    }
    rel.events$t1 <- rel.events$t0 + age.max
    rel.events$t1[rel.events$t1 > dat$trange[2]] <- dat$trange[2]

    ## Discretisation for expm
    rel.events$itrel <- as.integer(cut(rel.events$t0, time.cont,
                                       include.lowest = TRUE))
    rel.events$itrec <- as.integer(cut(rel.events$t1, time.cont,
                                       include.lowest = TRUE))
    rel.events$icrel <- dat$celltable[cbind(as.integer(cut(rel.events$x0, dat$xgr)),
                                            as.integer(cut(rel.events$y0, dat$ygr)))]

    res <- list(rel.events = rel.events,
                idx = idx)

    return(res)
}



get.par.names <- function(fit){
    tab <- table(names(fit$opt$par))
    res <- unlist(sapply(seq_along(tab), function(x) if(tab[x] > 1) paste0(names(tab)[x],1:tab[x]) else names(tab)[x]))
    return(res)
}



unscented.transform <- function(x, Px){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
    L <- length(x)
    W0 <- 0
    chi <- array(0, dim = c(L, 2*L+1))
    chi[,1] <- x
    if(!RTMB:::ad_context()){
        sqrtMat <- chol((L/(1-W0)) * Px)
    }else{
        sqrtMat <- RTMB:::chol.advector((L/(1-W0)) * Px)
    }
    for(i in 1:L){
        chi[, i + 1] <- x + sqrtMat[, i]
        chi[, i + 1 + L] <- x - sqrtMat[, i]
    }
    Wm <- Wc <- RTMB::matrix(c(W0, rep((1-W0)/(2*L), 2*L)), 1, 2*L+1)
    return(list(chi = chi, Wm = Wm, Wc = Wc))
}



date.2.decimal.year <- function(dates){
    leap_year <- function(year) {
        (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
    }
    year <- as.numeric(format(dates, "%Y"))
    doy <- as.numeric(format(dates, "%j"))
    is_leap <- leap_year(year)
    days_in_year <- ifelse(is_leap, 366, 365)
    dec_year <- year + (doy - 1) / days_in_year
    return(dec_year)
}



shift.map <- function(database = "world"){
  map_data <- maps::map(database, plot = FALSE, fill = TRUE)

  x <- map_data$x
  y <- map_data$y
  names <- map_data$names

  na_idx <- which(is.na(x))
  start_idx <- c(1, na_idx + 1)
  end_idx <- c(na_idx - 1, length(x))

  x_new <- c()
  y_new <- c()
  names_new <- c()

  for (i in seq_along(start_idx)) {
    xi <- x[start_idx[i]:end_idx[i]]
    yi <- y[start_idx[i]:end_idx[i]]

    lon_range <- range(xi, na.rm = TRUE)
    if (lon_range[1] < 0 && lon_range[2] > 0) next

    xi_shifted <- ifelse(xi < 0, xi + 360, xi)

    x_new <- c(x_new, xi_shifted, NA)
    y_new <- c(y_new, yi, NA)
    names_new <- c(names_new, names[i])
  }

  list(x = x_new, y = y_new, names = names_new, range = range(x_new, na.rm = TRUE))
}



group.consecutive.ranges <- function(x){
  x <- sort(unique(x))
  breaks <- c(0, which(diff(x) != 1), length(x))

  ranges <- mapply(function(i, j) {
    start <- x[i + 1]
    end <- x[j]
    if (start == end) {
      as.character(start)
    } else {
      paste0(start, ":", end)
    }
  }, breaks[-length(breaks)], breaks[-1])

  res <- paste0("c(",paste(ranges,collapse = ","), ")")

  return(res)
}


t2index <- function(t, time.cont){
    findInterval(t, time.cont, rightmost.closed = TRUE)
}
